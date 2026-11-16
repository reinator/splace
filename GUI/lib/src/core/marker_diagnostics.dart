/// Why a marker is missing from a genome, and how the copies of a
/// duplicated marker relate: the explanations of step 4's "Missing
/// markers" and "Duplicate markers" panels.
library;

import 'genbank.dart';
import 'genbank_markers.dart';
import 'syngenes.dart';

enum MissingReason {
  /// The gene is annotated under a feature type not selected in step 3.
  typeNotSelected,

  /// Only a `gene` feature, without the CDS/rRNA feature SPLACE reads.
  geneFeatureOnly,

  /// The record has no CDS (or rRNA) annotated at all.
  noAnnotation,

  /// The record has CDS/rRNA names SPLACE did not recognize; one of them
  /// may be this gene.
  unknownName,

  /// The record is much shorter than the others (partial genome).
  partialGenome,

  /// Annotated genes do not include it (absent, or not annotated).
  notAnnotated,

  /// FASTA input: no sequence of this locus for the sample.
  notInFiles,
}

class MissingCause {
  const MissingCause(this.reason, [this.detail]);
  final MissingReason reason;

  /// Feature type, unrecognized names, location or lengths, by reason.
  final String? detail;
}

/// Median length of the records (the "typical genome" of the dataset).
int typicalRecordLength(List<GbRecord> records) {
  if (records.isEmpty) return 0;
  final l = [for (final r in records) r.sequence.length]..sort();
  return l[l.length ~/ 2];
}

/// Why [gene] is missing from record [index] of [records].
MissingCause diagnoseMissing({
  required List<GbRecord> records,
  required int index,
  required String gene,
  required GeneCatalog catalog,
  required Set<String> types,
  Map<String, String> overrides = const {},
  int? typical,
}) {
  final r = records[index];
  // Feature types the gene comes from in the other genomes.
  final geneTypes = {
    for (final t in markerFeatureTypes)
      if (catalog.byType[t]?[gene] != null) t,
  };
  for (final t in markerFeatureTypes) {
    if (types.contains(t)) continue;
    if (catalog.byType[t]?[gene]?.byRecord.containsKey(index) ?? false) {
      return MissingCause(MissingReason.typeNotSelected, featureTypeLabels[t] ?? t);
    }
  }
  for (final f in r.features) {
    if (f.type != 'gene') continue;
    if (standardizeGeneName(f.q('gene'), f.q('product'), catalog.organelle, overrides) == gene) {
      return MissingCause(MissingReason.geneFeatureOnly, f.location);
    }
  }
  final wanted = geneTypes.isEmpty ? const {'CDS', 'rRNA'} : geneTypes;
  if (!r.features.any((f) => wanted.contains(f.type))) {
    return MissingCause(MissingReason.noAnnotation, wanted.map((t) => featureTypeLabels[t] ?? t).join('/'));
  }
  final unknown = [
    for (final u in catalog.unknown)
      if (u.records.contains(index) && wanted.contains(u.type)) u.name,
  ];
  if (unknown.isNotEmpty) return MissingCause(MissingReason.unknownName, unknown.join(', '));
  final t = typical ?? typicalRecordLength(records);
  if (t > 0 && r.sequence.length < t * 0.8) {
    return MissingCause(MissingReason.partialGenome, '${r.sequence.length}|$t');
  }
  return const MissingCause(MissingReason.notAnnotated);
}

enum CopyRelation {
  /// Same sequence (e.g. inverted-repeat copies of a plastome).
  identical,

  /// The shorter copies are part of the longest (a fragment or a
  /// truncated duplicate annotation).
  contained,

  /// The copies differ (paralogues, NUMT-like copies or annotation
  /// errors): the choice changes the alignment.
  different,
}

/// How the copies [seqs] of one gene in one genome relate.
CopyRelation copyRelation(List<String> seqs) {
  final s = [for (final x in seqs) x.toUpperCase().replaceAll(RegExp('[^A-Z]'), '')]..sort((a, b) => b.length.compareTo(a.length));
  if (s.length < 2) return CopyRelation.identical;
  if (s.every((x) => x == s.first)) return CopyRelation.identical;
  final longest = s.first;
  final rc = reverseComplement(longest).toUpperCase();
  if (s.skip(1).every((x) => x.isNotEmpty && (longest.contains(x) || rc.contains(x)))) return CopyRelation.contained;
  return CopyRelation.different;
}
