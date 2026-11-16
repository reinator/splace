import 'dart:math' as math;

import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/core/fasta.dart';
import 'package:splace/src/core/saturation.dart';

const _bases = ['A', 'C', 'G', 'T'];

/// Simulates codon-structured sequences under JC69 on a balanced binary
/// tree with [depth] levels and equal branch lengths, using a different
/// rate multiplier for each codon position.
List<SeqRecord> simulateCodons({
  required int codons,
  required List<double> rates,
  int depth = 4,
  double branch = 0.1,
  int seed = 1,
}) {
  final rnd = math.Random(seed);
  String mutate(String seq, double t) {
    final out = StringBuffer();
    for (var i = 0; i < seq.length; i++) {
      final rate = rates[i % 3];
      final pChange = 0.75 * (1 - math.exp(-4 / 3 * rate * t));
      if (rnd.nextDouble() < pChange) {
        final others = _bases.where((b) => b != seq[i]).toList();
        out.write(others[rnd.nextInt(3)]);
      } else {
        out.write(seq[i]);
      }
    }
    return out.toString();
  }

  var level = [List.generate(codons * 3, (_) => _bases[rnd.nextInt(4)]).join()];
  for (var d = 0; d < depth; d++) {
    level = [for (final s in level) ...[mutate(s, branch), mutate(s, branch)]];
  }
  return [for (var i = 0; i < level.length; i++) SeqRecord('taxon$i', level[i])];
}

void main() {
  group('distances', () {
    test('K2P equals JC when there are no transversion/transition biases', () {
      // JC: d = -3/4 ln(1 - 4p/3). With P = 2p/3·½... use p = 0.3 split
      // as Ts = 0.1, Tv = 0.2 (JC expectation: Ts = p/3, Tv = 2p/3).
      final d = k2pDistance(0.1, 0.2)!;
      final jc = -0.75 * math.log(1 - 4 * 0.3 / 3);
      expect(d, closeTo(jc, 1e-9));
    });

    test('TN93 reduces to K2P with equal base frequencies', () {
      final tn = tn93Distance(0.05, 0.05, 0.2, [0.25, 0.25, 0.25, 0.25])!;
      expect(tn, closeTo(k2pDistance(0.1, 0.2)!, 1e-9));
    });

    test('undefined beyond the correctable limit', () {
      expect(k2pDistance(0.5, 0.3), isNull);
    });
  });

  group('frame', () {
    test('codon position mapping honours the frame offset', () {
      expect([for (var c = 0; c < 6; c++) codonPositionOf(c, 0)], [1, 2, 3, 1, 2, 3]);
      expect([for (var c = 0; c < 6; c++) codonPositionOf(c, 1)], [3, 1, 2, 3, 1, 2]);
    });

    test('stop codons pick the right frame', () {
      // Frame 0 reads ATC AAA TAG ... (internal stop); frame 1 reads
      // TCA AAT AGC ... (none).
      final aln = [SeqRecord('a', 'ATCAAATAGCCCGGGTTT'), SeqRecord('b', 'ATCAAATAGCCCGGGTTT')];
      expect(countInternalStops(aln, 0, 1), 2);
      expect(countInternalStops(aln, 1, 1), 0);
    });
  });

  group('saturation test on simulated data', () {
    final alignment = simulateCodons(codons: 400, rates: [0.3, 0.1, 6.0]);
    final report = testSaturation('sim', alignment, const SaturationOptions(frame: 0));

    test('3rd positions are flagged saturated, 1st and 2nd are not', () {
      expect(report.subsets[CodonSubset.pos3]!.verdict, SaturationVerdict.saturated);
      expect(report.subsets[CodonSubset.pos2]!.verdict, SaturationVerdict.ok);
      expect(report.subsets[CodonSubset.pos1]!.verdict, isNot(SaturationVerdict.saturated));
    });

    test('Iss orders positions by rate', () {
      final iss = [CodonSubset.pos2, CodonSubset.pos1, CodonSubset.pos3]
          .map((s) => report.subsets[s]!.iss!)
          .toList();
      expect(iss[0], lessThan(iss[1]));
      expect(iss[1], lessThan(iss[2]));
      expect(iss[2], greaterThan(0.85), reason: 'fast 3rd positions look nearly random');
    });

    test('recommendation removes 3rd positions only', () {
      final rec = report.recommendation;
      expect(rec.kind, RecommendationKind.excludePositions);
      expect(rec.exclude, {3});
    });

    test('slowly evolving data keeps all positions', () {
      final slow = simulateCodons(codons: 400, rates: [0.2, 0.1, 0.5], seed: 3);
      final r = testSaturation('slow', slow, const SaturationOptions(frame: 0));
      expect(r.recommendation.kind, RecommendationKind.keepAll);
    });

    test('everything saturated suggests amino acids', () {
      final fast = simulateCodons(codons: 300, rates: [8, 8, 8], seed: 5);
      final r = testSaturation('fast', fast, const SaturationOptions(frame: 0));
      expect(r.recommendation.kind, RecommendationKind.useAminoAcids);
    });

    test('plot points stay within bounds', () {
      for (final s in report.subsets.values) {
        expect(s.points.length, lessThanOrEqualTo(3000));
        for (final p in s.points) {
          expect(p.p, inInclusiveRange(0, 0.75));
          expect(p.distance, greaterThanOrEqualTo(p.p - 1e-9));
        }
      }
    });
  });

  group('complementary analyses', () {
    test('chi-square upper tail matches tables', () {
      expect(chiSquareP(3.841, 1), closeTo(0.05, 0.001));
      expect(chiSquareP(7.815, 3), closeTo(0.05, 0.001));
      expect(chiSquareP(0, 3), 1);
      expect(chiSquareP(30, 3), lessThan(1e-5));
    });

    test('per-position p and composition on simulated data', () {
      final aln = simulateCodons(codons: 300, rates: [0.3, 0.1, 6.0], seed: 3);
      final r = testSaturation('sim', aln, const SaturationOptions(frame: 0));
      expect(r.positionPairs, isNotEmpty);
      final mean3 = r.positionPairs.map((p) => p.p3).reduce((a, b) => a + b) / r.positionPairs.length;
      final mean2 = r.positionPairs.map((p) => p.p2).reduce((a, b) => a + b) / r.positionPairs.length;
      expect(mean3, greaterThan(mean2));
      final s3 = r.subsets[CodonSubset.pos3]!;
      expect(s3.variableSites, greaterThan(r.subsets[CodonSubset.pos2]!.variableSites));
      expect(s3.informativeSites, lessThanOrEqualTo(s3.variableSites));
      // JC simulation: homogeneous composition.
      expect(s3.composition!.p, greaterThan(0.01));
    });

    test('composition heterogeneity is detected', () {
      final aln = [
        for (var i = 0; i < 4; i++) SeqRecord('at$i', 'ATATATAATTATATTAAT' * 4),
        for (var i = 0; i < 4; i++) SeqRecord('gc$i', 'GCGCGCGGCCGCGCCGGC' * 4),
      ];
      final r = testSaturation('comp', aln, const SaturationOptions(frame: 0));
      final c = r.subsets[CodonSubset.all]!.composition!;
      expect(c.p, lessThan(0.001));
      expect(c.failing, 8);
    });
  });

  test('rejects unaligned input', () {
    expect(
      () => testSaturation('x', [SeqRecord('a', 'ACGT'), SeqRecord('b', 'ACG')]),
      throwsArgumentError,
    );
  });
}
