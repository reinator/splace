import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/core/genbank.dart';
import 'package:splace/src/core/genbank_markers.dart';
import 'package:splace/src/core/marker_diagnostics.dart';
import 'package:splace/src/services/ncbi.dart';

GbRecord _rec(String acc, List<GbFeature> features, int length) => GbRecord(
      accession: acc,
      organism: 'Caranx hippos',
      lineage: '',
      definition: '',
      features: features,
      sequence: 'a' * length,
    );

GbFeature _f(String type, String gene, String loc) => GbFeature(type, loc, {'gene': gene});

void main() {
  test('why a gene is missing', () {
    final full = [_f('CDS', 'COX1', '1..1500'), _f('CDS', 'ND1', '1600..2500'), _f('rRNA', 'rrnL', '3000..4600')];
    final records = [
      _rec('A', full, 16500),
      _rec('B', full, 16600),
      // ND1 only as a 'gene' feature.
      _rec('C', [_f('CDS', 'COX1', '1..1500'), _f('gene', 'ND1', '1600..2500')], 16500),
      // Raw submission: nothing annotated.
      _rec('D', const [], 16400),
      // Partial genome.
      _rec('E', [_f('CDS', 'COX1', '1..1500')], 4000),
      // An unrecognized CDS name.
      _rec('F', [_f('CDS', 'COX1', '1..1500'), _f('CDS', 'NADH-like orf', '1600..2500')], 16500),
      // Plain absence.
      _rec('G', [_f('CDS', 'COX1', '1..1500')], 16500),
    ];
    final cat = buildGeneCatalog(records);
    MissingReason why(int i, String gene, {Set<String> types = const {'CDS'}}) =>
        diagnoseMissing(records: records, index: i, gene: gene, catalog: cat, types: types).reason;
    expect(why(2, 'ND1'), MissingReason.geneFeatureOnly);
    expect(why(3, 'COI'), MissingReason.noAnnotation);
    expect(why(4, 'ND1'), MissingReason.partialGenome);
    expect(why(5, 'ND1'), MissingReason.unknownName);
    expect(why(6, 'ND1'), MissingReason.notAnnotated);
    // 16S (rrnL) is in A as an rRNA, a type not selected.
    expect(why(0, '16S'), MissingReason.typeNotSelected);
  });

  test('copies of a duplicated gene', () {
    expect(copyRelation(['ACGTACGTAA', 'acgtacgtaa']), CopyRelation.identical);
    expect(copyRelation(['ACGTACGTAAGG', 'GTACGT']), CopyRelation.contained);
    // Reverse-complemented fragment.
    expect(copyRelation(['AAACCCGGGTTTA', reverseComplement('CCCGGG').toUpperCase()]), CopyRelation.contained);
    expect(copyRelation(['ACGTACGTAA', 'TTTTGGGGCC']), CopyRelation.different);
  });

  test('download estimates grow with the number of genomes and the key helps', () {
    Duration est(int n, bool key) => NcbiClient.estimateFetch(n, withKey: key);
    for (final key in [true, false]) {
      expect(est(50, key) < est(200, key), isTrue);
      expect(est(200, key) < est(1000, key), isTrue);
    }
    for (final n in [50, 200, 1000]) {
      expect(est(n, true) < est(n, false), isTrue);
    }
    // Measured on Nov 8 2026 with Testudines: 200 genomes in 4.4 s and 700
    // in 10.1 s with a key, 200 in 16.9 s and 700 in 64.6 s without.
    expect(est(200, true).inMilliseconds, inInclusiveRange(3500, 5500));
    expect(est(700, true).inMilliseconds, inInclusiveRange(8500, 12000));
    expect(est(200, false).inMilliseconds, inInclusiveRange(15000, 21000));
    expect(est(700, false).inMilliseconds, inInclusiveRange(58000, 70000));
  });
}
