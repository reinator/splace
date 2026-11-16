/// Suggestions for gene names that SynGenes could not standardize, and
/// checks of features whose /gene and /product name different genes.
///
/// Two kinds of evidence:
/// - the qualifiers: /gene and /product looked up on their own, also after
///   normalizing spelling ("NADH-dehydrogenase subunit 4 L" → "nadhdehydrogenase4l");
/// - the sequence: the feature is compared with the copies of every
///   standardized marker found in the other records (protein 3-mers for
///   CDS, nucleotide 8-mers for rRNA), like a quick BLAST against the
///   dataset itself.
library;

import 'genbank.dart';
import 'genbank_markers.dart';
import 'genetic_code.dart';
import 'syngenes.dart';
import 'syngenes_data.dart';

enum SuggestionSource { product, gene, sequence }

/// Best match of a sequence against the standardized markers.
class SequenceMatch {
  const SequenceMatch(this.gene, this.score, this.references, {this.ranking = const []});
  final String gene;

  /// The closest markers, best first (gene, share of k-mers 0–1).
  final List<(String, double)> ranking;

  /// Share of the query's k-mers found in the closest copy (0–1).
  final double score;

  /// Copies of [gene] compared.
  final int references;
}

class GeneSuggestion {
  GeneSuggestion({
    required this.name,
    required this.type,
    this.suggested,
    this.source,
    this.match,
    this.genes = const {},
    this.products = const {},
    this.byGene,
    this.byProduct,
    this.ranking = const [],
    this.copies = 0,
  });

  /// Closest markers to this name's sequences (gene, similarity 0–1), also
  /// when none is clearly ahead.
  final List<(String, double)> ranking;

  /// Features carrying this name.
  final int copies;

  /// Name that SynGenes did not recognize.
  final String name;
  final String type;

  /// Standard name to use, or null when nothing is convincing.
  final String? suggested;
  final SuggestionSource? source;

  /// Sequence evidence, also when the suggestion came from a qualifier.
  final SequenceMatch? match;

  /// /gene and /product values seen for this name.
  final Set<String> genes;
  final Set<String> products;

  /// Standard names the /gene and /product values point to, if any.
  final String? byGene, byProduct;

  /// Both qualifiers are recognized but name different genes.
  bool get qualifiersDisagree => byGene != null && byProduct != null && byGene != byProduct;
}

/// A feature whose /gene and /product resolve to two different genes.
class QualifierConflict {
  QualifierConflict({
    required this.record,
    required this.feature,
    required this.type,
    required this.gene,
    required this.product,
    required this.geneName,
    required this.productName,
    required this.used,
    this.match,
    this.geneScore,
    this.productScore,
    this.scores = const {},
  });

  /// Similarity of the sequence to every marker (0–1), best first.
  final Map<String, double> scores;

  /// Similarity of the sequence to the genes named by /gene and /product
  /// in the other genomes (0–1), null when there is no reference.
  final double? geneScore, productScore;

  final int record;
  final int feature;
  final String type;

  /// Raw qualifiers.
  final String gene, product;

  /// Standard names each qualifier points to.
  final String geneName, productName;

  /// Name SPLACE uses (the product wins, as in the HTML app).
  final String used;

  /// What the sequence looks like.
  final SequenceMatch? match;

  /// The sequence agrees with the /gene qualifier rather than the name used.
  bool get sequenceFavorsGene => match != null && match!.gene == geneName && geneName != used;
}

class GeneReview {
  GeneReview(this.suggestions, this.conflicts);
  final Map<String, GeneSuggestion> suggestions;
  final List<QualifierConflict> conflicts;
}

/// Minimum k-mer containment to trust a sequence match, and the lead it
/// needs over the second-best marker.
const _minScore = 0.35;
const _minLead = 0.1;
const _refsPerGene = 25;

String _norm(String s) => s
    .toLowerCase()
    .replaceAll(RegExp(r'\b(subunit|protein|gene|putative|partial)\b'), '')
    .replaceAll(RegExp(r'[^a-z0-9]'), '');

final Map<String, String> _normalized = () {
  final out = <String, String>{};
  for (final table in [mitochondrialGenes, chloroplastGenes]) {
    for (final e in table.entries) {
      out.putIfAbsent(_norm(e.key), () => e.key);
      for (final s in e.value) {
        out.putIfAbsent(_norm(s), () => e.key);
      }
    }
  }
  return out;
}();

/// Standard name of one qualifier value: exact SynGenes lookup first, then
/// the spelling-normalized tables. Null when unknown.
String? lookupGeneName(String? value, OrganelleType type, {bool asProduct = false}) {
  if (value == null || value.trim().isEmpty) return null;
  final exact = asProduct ? standardizeGeneName(null, value, type) : standardizeGeneName(value, null, type);
  if (isKnownGeneName(exact)) return exact;
  final n = _norm(value);
  return n.isEmpty ? null : _normalized[n];
}

Set<String> _kmers(String s, int k) {
  final out = <String>{};
  for (var i = 0; i + k <= s.length; i++) {
    final w = s.substring(i, i + k);
    if (!w.contains('X') && !w.contains('*') && !w.contains('N')) out.add(w);
  }
  return out;
}

/// Protein (CDS) or nucleotide (rRNA) sequence of a feature.
String featureSequence(GbRecord r, GbFeature f, int geneticCode) {
  if (f.type == 'CDS') {
    final given = f.q('translation')?.replaceAll(RegExp(r'\s'), '');
    if (given != null && given.isNotEmpty) return given.toUpperCase();
    final dna = extractLocation(r.sequence, f.location);
    final start = (int.tryParse(f.q('codon_start') ?? '1') ?? 1) - 1;
    final code = int.tryParse(f.q('transl_table') ?? '') ?? geneticCode;
    return dna.length > start ? translate(dna.substring(start), code) : '';
  }
  return extractLocation(r.sequence, f.location);
}

class _Reference {
  _Reference(this.record, this.kmers);
  final int record;
  final Set<String> kmers;
}

/// Compares sequences with the standardized markers of the dataset.
class MarkerMatcher {
  MarkerMatcher(List<GbRecord> records, GeneCatalog catalog, this.geneticCode) {
    for (final type in const ['CDS', 'rRNA']) {
      for (final g in catalog.byType[type]?.values ?? const <GeneSummary>[]) {
        if (!isKnownGeneName(g.name)) continue;
        final refs = _refs.putIfAbsent(type, () => {}).putIfAbsent(g.name, () => []);
        for (final hits in g.byRecord.values) {
          if (refs.length >= _refsPerGene) break;
          final h = hits.first;
          final r = records[h.record];
          final seq = featureSequence(r, r.features[h.feature], geneticCode);
          final k = _kmers(seq, _k(type));
          if (k.length >= 10) refs.add(_Reference(h.record, k));
        }
      }
    }
  }

  final int geneticCode;
  final _refs = <String, Map<String, List<_Reference>>>{};

  static int _k(String type) => type == 'CDS' ? 3 : 8;

  /// Similarity of [seq] to every marker of [type] (best copy, 0–1) and
  /// how many copies were compared, ignoring [excludeRecord]; best first.
  List<(String, double, int)> scores(String seq, String type, {int? excludeRecord}) {
    final q = _kmers(seq, _k(type));
    if (q.length < 10) return const [];
    final out = <(String, double, int)>[];
    for (final e in (_refs[type] ?? const <String, List<_Reference>>{}).entries) {
      var best = 0.0;
      var n = 0;
      for (final ref in e.value) {
        if (ref.record == excludeRecord) continue;
        n++;
        final s = q.where(ref.kmers.contains).length / q.length;
        if (s > best) best = s;
      }
      if (n > 0) out.add((e.key, best, n));
    }
    return out..sort((a, b) => b.$2.compareTo(a.$2));
  }

  /// Best marker for [seq] of feature [type], ignoring copies from
  /// [excludeRecord]; null when no marker is clearly ahead.
  SequenceMatch? match(String seq, String type, {int? excludeRecord}) {
    final all = scores(seq, type, excludeRecord: excludeRecord);
    if (all.isEmpty) return null;
    final top = all.first;
    final second = all.length > 1 ? all[1].$2 : 0.0;
    if (top.$2 < _minScore || top.$2 - second < _minLead) return null;
    return SequenceMatch(top.$1, top.$2, top.$3, ranking: [for (final s in all.take(3)) (s.$1, s.$2)]);
  }
}

/// Suggestions for every unknown name of [catalog] and the /gene vs
/// /product conflicts of all CDS and rRNA features.
GeneReview reviewGeneNames(
  List<GbRecord> records,
  GeneCatalog catalog, {
  required int geneticCode,
  Map<String, String> overrides = const {},
  Map<String, String> featureOverrides = const {},
}) {
  final type = catalog.organelle;
  final matcher = MarkerMatcher(records, catalog, geneticCode);
  final unknown = {for (final u in catalog.unknown) u.name: u};

  // Features carrying each unknown name.
  final features = <String, List<(int, int)>>{};
  final conflicts = <QualifierConflict>[];
  for (var ri = 0; ri < records.length; ri++) {
    final r = records[ri];
    for (final e in catalog.recordGenes[ri].entries) {
      final f = r.features[e.key];
      if (f.type == 'tRNA') continue;
      if (unknown.containsKey(e.value)) features.putIfAbsent(e.value, () => []).add((ri, e.key));

      // Both qualifiers known but naming different genes.
      final gene = f.q('gene'), product = f.q('product');
      if (gene == null || product == null || overrides.containsKey(gene) || overrides.containsKey(product)) continue;
      if (featureOverrides.containsKey(featureKey(r, e.key))) continue;
      final byGene = lookupGeneName(gene, type);
      final byProduct = lookupGeneName(product, type, asProduct: true);
      if (byGene == null || byProduct == null || byGene == byProduct) continue;
      final seq = featureSequence(r, f, geneticCode);
      final all = {for (final s in matcher.scores(seq, f.type, excludeRecord: ri)) s.$1: s.$2};
      conflicts.add(QualifierConflict(
        record: ri,
        feature: e.key,
        type: f.type,
        gene: gene,
        product: product,
        geneName: byGene,
        productName: byProduct,
        used: e.value,
        match: matcher.match(seq, f.type, excludeRecord: ri),
        geneScore: all[byGene],
        productScore: all[byProduct],
        scores: all,
      ));
    }
  }

  final suggestions = <String, GeneSuggestion>{};
  for (final u in catalog.unknown) {
    final hits = features[u.name] ?? const [];
    final genes = <String>{}, products = <String>{};
    // One vote per copy; the marker with most votes wins, keeping its best
    // score.
    final votes = <String, int>{};
    final bestOf = <String, SequenceMatch>{};
    var ranking = const <(String, double)>[];
    for (final (ri, fi) in hits) {
      final f = records[ri].features[fi];
      if (ranking.isEmpty) {
        ranking = [
          for (final s in matcher.scores(featureSequence(records[ri], f, geneticCode), f.type, excludeRecord: ri).take(3)) (s.$1, s.$2),
        ];
      }
      if (f.q('gene') case final g?) genes.add(g);
      if (f.q('product') case final p?) products.add(p);
      final m = matcher.match(featureSequence(records[ri], f, geneticCode), f.type, excludeRecord: ri);
      if (m == null) continue;
      votes.update(m.gene, (v) => v + 1, ifAbsent: () => 1);
      if ((bestOf[m.gene]?.score ?? -1) < m.score) bestOf[m.gene] = m;
    }
    final winner = votes.isEmpty
        ? null
        : (votes.entries.toList()..sort((a, b) => b.value != a.value ? b.value - a.value : bestOf[b.key]!.score.compareTo(bestOf[a.key]!.score))).first.key;
    final match = winner == null ? null : bestOf[winner];

    String? byProduct, byGene;
    for (final p in products) {
      byProduct ??= lookupGeneName(p, type, asProduct: true);
    }
    for (final g in genes) {
      byGene ??= lookupGeneName(g, type);
    }
    final (suggested, source) = byProduct != null
        ? (byProduct, SuggestionSource.product)
        : byGene != null
            ? (byGene, SuggestionSource.gene)
            : match != null
                ? (match.gene, SuggestionSource.sequence)
                : (null, null);
    suggestions[u.name] = GeneSuggestion(
      name: u.name,
      type: u.type,
      suggested: suggested,
      source: source,
      match: match,
      genes: genes,
      products: products,
      byGene: byGene,
      byProduct: byProduct,
      ranking: match?.ranking.isNotEmpty == true ? match!.ranking : ranking,
      copies: hits.length,
    );
  }
  return GeneReview(suggestions, conflicts);
}
