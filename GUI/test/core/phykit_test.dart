import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/core/phykit.dart';
import 'package:splace/src/core/saturation.dart';

void main() {
  test('PhyKIT saturation JSON: value from --json, pairs from -v --json', () {
    // Verbatim PhyKIT 2.8.0 output (shortened).
    const summary = '{"exclude_gaps": false, "summary": {"one_minus_slope_abs": 0.7223, "slope": 0.2777}, "verbose": false}';
    const verbose = '{"exclude_gaps": false, "pairs": [], "rows": ['
        '{"patristic_distance": 0.3359, "taxon_a": "S1", "taxon_b": "S2", "uncorrected_distance": 0.2128}, '
        '{"patristic_distance": 1.4476, "taxon_a": "S1", "taxon_b": "S3", "uncorrected_distance": 0.4275}, '
        '{"patristic_distance": 1.338, "taxon_a": "S2", "taxon_b": "S3", "uncorrected_distance": 0.4081}], "verbose": true}';
    final r = parsePhykitSaturation(summary, verbose);
    expect(r.slope, 0.2777);
    expect(r.oneMinus, closeTo(0.7223, 1e-9));
    expect(r.pairs.length, 3);
    expect(r.pairs.last.patristic, 1.338);
    final s = phykitSubset(CodonSubset.pos3, 200, r, phykitThresholds);
    expect(s.verdict, SaturationVerdict.saturated);
    expect(s.meanDistance, closeTo((0.3359 + 1.4476 + 1.338) / 3, 1e-9));
    // Too few pairs to judge.
    final few = phykitSubset(CodonSubset.pos3, 200, r, const SaturationThresholds(slopeModerate: 0.5, slopeSaturated: 0.3, minPairs: 10));
    expect(few.verdict, SaturationVerdict.insufficient);
    expect(subsetColumns(CodonSubset.pos3, 9, 1), [0, 3, 6]);
  });
}
