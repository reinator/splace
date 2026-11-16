/// Markers (genes) found in GenBank records: feature types, standardized
/// names, duplicates, FASTA headers and per-gene FASTA extraction — the
/// steps 3–5a of the HTML app.
library;

import 'fasta.dart';
import 'genbank.dart';
import 'genetic_code.dart';
import 'syngenes.dart';

const markerFeatureTypes = ['CDS', 'rRNA', 'tRNA'];
const featureTypeLabels = {'CDS': 'PCGs', 'rRNA': 'rRNAs', 'tRNA': 'tRNAs'};

class GeneHit {
  GeneHit(this.record, this.feature, this.type, this.length);
  final int record;
  final int feature;
  final String type;
  final int length;
}

class GeneSummary {
  GeneSummary(this.name, this.product);
  final String name;
  final String product;
  final Set<String> types = {};

  /// record index -> hits (more than one = duplicated annotation).
  final Map<int, List<GeneHit>> byRecord = {};

  int get count => byRecord.values.fold(0, (a, l) => a + l.length);
  Iterable<int> get lengths => byRecord.values.expand((l) => l).map((h) => h.length).where((l) => l > 0);
}

class UnknownGene {
  UnknownGene(this.name, this.type);
  final String name;
  final String type;
  final List<int> records = [];
}

class GeneCatalog {
  GeneCatalog(this.organelle, this.featureTypeCounts, this.byType, this.unknown, this.recordGenes);

  final OrganelleType organelle;

  /// CDS/rRNA/tRNA → number of features.
  final Map<String, int> featureTypeCounts;

  /// feature type → gene name → summary.
  final Map<String, Map<String, GeneSummary>> byType;

  /// CDS/rRNA names not found in SynGenes.
  final List<UnknownGene> unknown;

  /// record index → feature index → standardized gene name.
  final List<Map<int, String>> recordGenes;

  /// Genes of the selected feature types, merged across types.
  Map<String, GeneSummary> visible(Set<String> types) {
    final out = <String, GeneSummary>{};
    for (final t in markerFeatureTypes) {
      if (!types.contains(t)) continue;
      for (final g in byType[t]?.values ?? const <GeneSummary>[]) {
        final merged = out.putIfAbsent(g.name, () => GeneSummary(g.name, g.product));
        merged.types.addAll(g.types);
        for (final e in g.byRecord.entries) {
          merged.byRecord.putIfAbsent(e.key, () => []).addAll(e.value);
        }
      }
    }
    return Map.fromEntries(out.entries.toList()..sort((a, b) => a.key.compareTo(b.key)));
  }
}

/// Name of a feature as used in the marker list, or null.
String? featureGeneName(GbFeature f, OrganelleType type, Map<String, String> overrides) {
  if (f.type == 'tRNA') return standardizeTrnaName(f);
  if (f.q('gene') == null && f.q('product') == null) {
    return f.type == 'CDS' ? resolveCodonGeneName(f, type, overrides) : null;
  }
  return standardizeGeneName(f.q('gene'), f.q('product'), type, overrides);
}

/// Key of one feature in [buildGeneCatalog]'s `featureOverrides`.
String featureKey(GbRecord r, int feature) => '${r.accession}#$feature';

/// [overrides] map raw names to standard ones in every record;
/// [featureOverrides] (by [featureKey]) fix single features, e.g. one
/// genome whose /gene and /product disagree.
GeneCatalog buildGeneCatalog(
  List<GbRecord> records, {
  Map<String, String> overrides = const {},
  Map<String, String> featureOverrides = const {},
  OrganelleType? organelle,
}) {
  final type = organelle ?? detectOrganelleType(records);
  final counts = <String, int>{};
  final byType = <String, Map<String, GeneSummary>>{};
  final unknown = <String, UnknownGene>{};
  final recordGenes = <Map<int, String>>[];
  for (var ri = 0; ri < records.length; ri++) {
    final r = records[ri];
    final names = <int, String>{};
    final trnas = <String>{};
    for (var fi = 0; fi < r.features.length; fi++) {
      final f = r.features[fi];
      if (!markerFeatureTypes.contains(f.type)) continue;
      counts.update(f.type, (v) => v + 1, ifAbsent: () => 1);
      var name = featureOverrides[featureKey(r, fi)] ?? featureGeneName(f, type, overrides);
      if (name == null || name.isEmpty) continue;
      if (f.type == 'tRNA') {
        // Two copies with the same name in one genome: the second is the
        // other Ser/Leu isoacceptor.
        if (trnas.contains(name)) {
          if (name == 'tRNA-Ser1') name = 'tRNA-Ser2';
          if (name == 'tRNA-Leu1') name = 'tRNA-Leu2';
        }
        trnas.add(name);
      } else if (!isKnownGeneName(name) && !overrides.containsValue(name)) {
        unknown.putIfAbsent(name, () => UnknownGene(name!, f.type)).records.add(ri);
      }
      names[fi] = name;
      final summary = byType.putIfAbsent(f.type, () => {}).putIfAbsent(name, () => GeneSummary(name!, f.q('product') ?? name));
      summary.types.add(f.type);
      summary.byRecord.putIfAbsent(ri, () => []).add(GeneHit(ri, fi, f.type, locationLength(f.location)));
    }
    recordGenes.add(names);
  }
  return GeneCatalog(type, counts, byType, unknown.values.toList()..sort((a, b) => a.name.compareTo(b.name)), recordGenes);
}

// ---------------------------------------------------------------------------
// FASTA headers (5a "Configure FASTA Header")
// ---------------------------------------------------------------------------

const headerFieldKeys = ['species', 'accession', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'authorship'];

String headerFieldValue(GbRecord r, String key) {
  final t = r.taxonomy;
  String na(String v) => v.trim().isEmpty ? 'NA' : v.trim();
  return switch (key) {
    'species' => (r.organismName.isEmpty ? 'Unknown' : r.organismName).replaceAll(RegExp(r'\s+'), '_'),
    'accession' => na(r.accession),
    'kingdom' => na(t.kingdom),
    'phylum' => na(t.phylum),
    'class' => na(t.klass),
    'order' => na(t.order),
    'family' => na(t.family),
    'genus' => na(t.genus),
    'authorship' => na(t.authorship),
    _ => key,
  };
}

class HeaderTemplate {
  const HeaderTemplate({this.fields = const ['accession', 'family', 'genus', 'species'], this.separator = '_'});
  final List<String> fields;
  final String separator;

  String build(GbRecord r) => sanitizeTaxonName(fields.map((k) => headerFieldValue(r, k)).join(separator));

  /// Headers shared by more than one record.
  Map<String, List<String>> duplicates(List<GbRecord> records) {
    final seen = <String, List<String>>{};
    for (final r in records) {
      seen.putIfAbsent(build(r), () => []).add(r.accession);
    }
    return {
      for (final e in seen.entries)
        if (e.value.length > 1) e.key: e.value,
    };
  }

  HeaderTemplate copyWith({List<String>? fields, String? separator}) =>
      HeaderTemplate(fields: fields ?? this.fields, separator: separator ?? this.separator);
}

String duplicateChoiceKey(GbRecord r, String gene) => '${r.accession}::$gene';

/// Index of the hit used for a duplicated gene: the user's choice, else
/// the longest copy.
int chosenHit(List<GeneHit> hits, int? choice) {
  if (choice != null && choice >= 0 && choice < hits.length) return choice;
  var best = 0;
  for (var i = 1; i < hits.length; i++) {
    if (hits[i].length > hits[best].length) best = i;
  }
  return best;
}

/// One FASTA per selected gene (nucleotide), one sequence per record.
Map<String, List<SeqRecord>> extractGeneFastas(
  List<GbRecord> records,
  GeneCatalog catalog,
  Set<String> types,
  Set<String> genes,
  HeaderTemplate header,
  Map<String, int> duplicateChoices,
) {
  final visible = catalog.visible(types);
  final out = <String, List<SeqRecord>>{};
  for (final gene in genes) {
    final summary = visible[gene];
    if (summary == null) continue;
    for (final e in summary.byRecord.entries) {
      final r = records[e.key];
      final hit = e.value[chosenHit(e.value, duplicateChoices[duplicateChoiceKey(r, gene)])];
      final seq = extractLocation(r.sequence, r.features[hit.feature].location);
      if (seq.isEmpty) continue;
      out.putIfAbsent(gene, () => []).add(SeqRecord(header.build(r), seq));
    }
  }
  return out;
}

/// A feature that GenBank marks as incomplete: its location starts with
/// `<` or ends with `>` (a fragment such as a COI barcode, "partial cds").
bool isPartialLocation(String location) => location.contains('<') || location.contains('>');

/// Records of each gene whose chosen feature is partial ([isPartialLocation]).
Map<String, int> partialCounts(
  List<GbRecord> records,
  GeneCatalog catalog,
  Set<String> types,
  Set<String> genes,
  Map<String, int> duplicateChoices,
) {
  final visible = catalog.visible(types);
  final out = <String, int>{};
  for (final gene in genes) {
    final summary = visible[gene];
    if (summary == null) continue;
    for (final e in summary.byRecord.entries) {
      final r = records[e.key];
      final hit = e.value[chosenHit(e.value, duplicateChoices[duplicateChoiceKey(r, gene)])];
      if (isPartialLocation(r.features[hit.feature].location)) out.update(gene, (v) => v + 1, ifAbsent: () => 1);
    }
  }
  return out;
}

/// A CDS prepared for codon-aware alignment.
class CodonCandidate {
  CodonCandidate(this.header, this.dna, this.code, this.codeFromFeature);
  final String header;

  /// CDS shifted by /codon_start.
  final String dna;
  final int code;
  final bool codeFromFeature;
}

/// CDS of the selected genes (codon mode). Pseudo CDS are skipped unless
/// [allowPseudo].
Map<String, List<CodonCandidate>> extractCodonCandidates(
  List<GbRecord> records,
  GeneCatalog catalog,
  Set<String> genes,
  HeaderTemplate header,
  Map<String, int> duplicateChoices, {
  required int fallbackCode,
  bool allowPseudo = false,
  void Function(String)? log,
}) {
  final cds = catalog.byType['CDS'] ?? const {};
  final out = <String, List<CodonCandidate>>{};
  for (final gene in genes) {
    final summary = cds[gene];
    if (summary == null) continue;
    for (final e in summary.byRecord.entries) {
      final r = records[e.key];
      final name = header.build(r);
      final usable = [
        for (final h in e.value)
          if (allowPseudo || !r.features[h.feature].isPseudo) h,
      ];
      if (usable.length < e.value.length) log?.call('[CDS extraction] [$gene] SKIP "$name": /pseudo CDS excluded by policy');
      if (usable.isEmpty) continue;
      final hit = usable[chosenHit(usable, duplicateChoices[duplicateChoiceKey(r, gene)])];
      final f = r.features[hit.feature];
      final raw = extractLocation(r.sequence, f.location);
      if (raw.isEmpty) {
        log?.call('[CDS extraction] [$gene] SKIP "$name": could not extract the CDS from its location');
        continue;
      }
      final codonStart = int.tryParse(f.q('codon_start') ?? '') ?? 1;
      final offset = (codonStart - 1).clamp(0, 2);
      if (offset > 0) log?.call('[CDS extraction] [$gene] "$name": /codon_start=$codonStart, trimmed $offset bp');
      final table = int.tryParse(f.q('transl_table') ?? '');
      out.putIfAbsent(gene, () => []).add(CodonCandidate(
            name,
            raw.substring(offset),
            table != null && isKnownGeneticCode(table) ? table : fallbackCode,
            table != null,
          ));
    }
  }
  return out;
}

/// Taxonomic tree node (Kingdom → … → Species) for the records panel.
class TaxonNode {
  TaxonNode(this.rank, this.name);
  final String rank;
  final String name;
  final List<int> records = [];
  final Map<String, TaxonNode> children = {};

  List<TaxonNode> get sortedChildren => children.values.toList()..sort((a, b) => a.name.compareTo(b.name));
}

TaxonNode buildTaxonomyTree(List<GbRecord> records) {
  final root = TaxonNode('root', '');
  for (var i = 0; i < records.length; i++) {
    final r = records[i];
    if (!r.hasFetchedTaxonomy) continue;
    final t = r.taxonomy;
    var node = root;
    for (final (rank, name) in [
      ('Kingdom', t.kingdom),
      ('Phylum', t.phylum),
      ('Class', t.klass),
      ('Order', t.order),
      ('Family', t.family),
      ('Genus', t.genus),
      ('Species', t.binomial.isEmpty ? r.accession : t.binomial),
    ]) {
      if (name.trim().isEmpty) continue;
      node = node.children.putIfAbsent('$rank:$name', () => TaxonNode(rank, name.trim()));
      node.records.add(i);
    }
  }
  return root;
}

/// Cell statistics for the PCG/rRNA/tRNA count columns: the mode is the
/// expected count; lower values hint at gene loss, higher at duplication.
class MarkerCountStats {
  MarkerCountStats(this.min, this.max, this.mode);
  final int min, max, mode;

  factory MarkerCountStats.of(Iterable<int> values) {
    final list = values.toList();
    if (list.isEmpty) return MarkerCountStats(0, 0, 0);
    final counts = <int, int>{};
    for (final v in list) {
      counts.update(v, (c) => c + 1, ifAbsent: () => 1);
    }
    final mode = (counts.entries.toList()..sort((a, b) => b.value != a.value ? b.value - a.value : a.key - b.key)).first.key;
    return MarkerCountStats(list.reduce((a, b) => a < b ? a : b), list.reduce((a, b) => a > b ? a : b), mode);
  }
}

/// Records of single markers (a COI from one voucher, a CYTB from another)
/// joined by taxon: the FASTA header is the taxon, and when a taxon has
/// several sequences of one marker the longest is kept. Returns the kept
/// sequences and how many were dropped per marker.
({Map<String, List<SeqRecord>> fastas, Map<String, int> dropped}) onePerTaxon(Map<String, List<SeqRecord>> fastas) {
  final out = <String, List<SeqRecord>>{};
  final dropped = <String, int>{};
  int bases(SeqRecord r) => r.sequence.replaceAll(RegExp(r'[-?Nn]'), '').length;
  for (final e in fastas.entries) {
    final best = <String, SeqRecord>{};
    for (final r in e.value) {
      final kept = best[r.header];
      if (kept == null || bases(r) > bases(kept)) best[r.header] = r;
    }
    out[e.key] = best.values.toList();
    if (e.value.length > best.length) dropped[e.key] = e.value.length - best.length;
  }
  return (fastas: out, dropped: dropped);
}

/// The same for CDS of the codon-aware mode.
List<CodonCandidate> oneCdsPerTaxon(List<CodonCandidate> cds) {
  final best = <String, CodonCandidate>{};
  for (final c in cds) {
    final kept = best[c.header];
    if (kept == null || c.dna.length > kept.dna.length) best[c.header] = c;
  }
  return best.values.toList();
}

/// FASTA header that names the taxon only: records of the same species
/// share it, so their markers end up in one row of the matrix.
const taxonHeader = HeaderTemplate(fields: ['family', 'species']);
