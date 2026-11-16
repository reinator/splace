import 'dart:math';

import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/core/gene_suggestions.dart';
import 'package:splace/src/core/genbank.dart';
import 'package:splace/src/core/genbank_markers.dart';

const _aa = 'ACDEFGHIKLMNPQRSTVWY';

String _protein(Random rng, int n) => List.generate(n, (_) => _aa[rng.nextInt(_aa.length)]).join();

/// [p] with about [rate] of its residues replaced.
String _mutate(Random rng, String p, double rate) =>
    p.split('').map((c) => rng.nextDouble() < rate ? _aa[rng.nextInt(_aa.length)] : c).join();

GbFeature _cds(String? gene, String? product, String translation) => GbFeature('CDS', '1..600', {
      'gene': ?gene,
      'product': ?product,
      'translation': translation,
    });

GbRecord _record(String acc, List<GbFeature> features) => GbRecord(
      accession: acc,
      organism: 'Caranx hippos',
      lineage: '',
      definition: 'Caranx hippos voucher MZUSP 1234 mitochondrion, complete genome.',
      features: features,
      sequence: '',
    );

void main() {
  final rng = Random(7);
  final cox1 = _protein(rng, 220);
  final nd4 = _protein(rng, 220);

  final records = [
    _record('NC_1', [
      _cds('COX1', 'cytochrome c oxidase subunit 1', cox1),
      _cds('ND4', 'NADH dehydrogenase subunit 4', nd4),
    ]),
    _record('NC_2', [
      _cds('COX1', 'cytochrome c oxidase subunit 1', _mutate(rng, cox1, 0.05)),
      _cds('ND4', 'NADH dehydrogenase subunit 4', _mutate(rng, nd4, 0.05)),
    ]),
    _record('NC_3', [
      // Unknown name whose sequence is COI.
      _cds('orf1', 'hypothetical protein', _mutate(rng, cox1, 0.08)),
      // /gene says ND4, /product says COI; the sequence is ND4.
      _cds('ND4', 'cytochrome c oxidase subunit 1', _mutate(rng, nd4, 0.08)),
    ]),
    _record('NC_4', [
      // Misspelled product that only the normalized lookup recognizes.
      _cds('xyz', 'NADH-dehydrogenase  subunit-4', _mutate(rng, nd4, 0.1)),
    ]),
  ];
  final catalog = buildGeneCatalog(records);
  final review = reviewGeneNames(records, catalog, geneticCode: 2);

  test('unknown name is suggested from its sequence', () {
    final s = review.suggestions['orf1']!;
    expect(s.suggested, 'COI');
    expect(s.source, SuggestionSource.sequence);
    expect(s.match!.score, greaterThan(0.5));
  });

  test('misspelled product is matched after normalizing', () {
    final s = review.suggestions['xyz']!;
    expect(s.suggested, 'ND4');
    expect(s.source, SuggestionSource.product);
    expect(s.match?.gene, 'ND4');
  });

  test('/gene vs /product conflict is reported with the sequence verdict', () {
    final k = review.conflicts.single;
    expect(k.geneName, 'ND4');
    expect(k.productName, 'COI');
    expect(k.used, 'COI');
    expect(k.match?.gene, 'ND4');
    expect(k.sequenceFavorsGene, isTrue);
  });

  test('a per-feature fix renames only that genome and clears the conflict', () {
    final k = review.conflicts.single;
    final fix = {featureKey(records[k.record], k.feature): k.geneName};
    final fixed = buildGeneCatalog(records, featureOverrides: fix);
    expect(fixed.recordGenes[k.record][k.feature], 'ND4');
    // NC_1's COI is untouched.
    expect(fixed.recordGenes[0][0], 'COI');
    expect(reviewGeneNames(records, fixed, geneticCode: 2, featureOverrides: fix).conflicts, isEmpty);
  });

  test('voucher comes from the genome title', () {
    expect(records.first.voucher, 'voucher MZUSP 1234');
    expect(titleVoucher('Bufo bufo isolate BB-2, complete genome'), 'isolate BB-2');
    expect(titleVoucher('Bufo bufo mitochondrion, complete genome'), isNull);
  });
}
