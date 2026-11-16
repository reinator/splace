/// Saturation by PhyKIT (Steenwyk et al. 2021): `phykit saturation -a
/// ALIGNMENT -t TREE -v --json` regresses the uncorrected p-distance of
/// every pair of taxa on their patristic distance in the tree (Philippe et
/// al. 2011). A slope ("saturation") near 1 means no saturation; the
/// closer to 0, the more saturated.
library;

import 'dart:convert';
import 'dart:math' as math;

import 'saturation.dart';

/// One run of `phykit saturation`.
class PhykitSaturation {
  const PhykitSaturation(this.slope, this.pairs);

  /// PhyKIT's saturation value (slope of p on patristic distance).
  final double? slope;

  /// Every pair of taxa: patristic distance and p-distance.
  final List<({String a, String b, double patristic, double p})> pairs;

  /// |1 − saturation|, PhyKIT's second value.
  double? get oneMinus => slope == null ? null : (1 - slope!).abs();
}

Map<String, dynamic> _json(String text) {
  final start = text.indexOf('{');
  if (start < 0) throw FormatException('PhyKIT printed no JSON: ${text.length > 200 ? text.substring(0, 200) : text}');
  return jsonDecode(text.substring(start)) as Map<String, dynamic>;
}

/// Reads PhyKIT's JSON: the saturation value from [summary] (`--json`)
/// and the pairs from [verbose] (`-v --json`, under "pairs" or "rows").
PhykitSaturation parsePhykitSaturation(String summary, [String? verbose]) {
  final d = _json(summary);
  final sum = d['summary'] as Map<String, dynamic>?;
  final slope = (sum?['slope'] ?? d['slope']) as num?;
  final v = verbose == null ? d : _json(verbose);
  final a = v['pairs'] as List? ?? const [], b = v['rows'] as List? ?? const [];
  final list = a.length >= b.length ? a : b;
  final pairs = [
    for (final x in list.whereType<Map<String, dynamic>>())
      (
        a: '${x['taxon_a']}',
        b: '${x['taxon_b']}',
        patristic: (x['patristic_distance'] as num?)?.toDouble() ?? 0,
        p: (x['uncorrected_distance'] as num?)?.toDouble() ?? 0,
      ),
  ];
  return PhykitSaturation(slope?.toDouble(), pairs);
}

/// Verdict of a PhyKIT [slope]: too little data, too little divergence to
/// matter, then the slope thresholds.
SaturationVerdict phykitVerdict(SaturationThresholds th, {required int sites, required int pairs, required double meanP, required double? slope}) {
  if (sites < th.minSites || pairs < th.minPairs || slope == null) return SaturationVerdict.insufficient;
  if (meanP < th.lowDivergence) return SaturationVerdict.ok;
  if (slope < th.slopeSaturated) return SaturationVerdict.saturated;
  if (slope < th.slopeModerate) return SaturationVerdict.moderate;
  return SaturationVerdict.ok;
}

/// One codon subset of a marker as the saturation report shows it, from
/// the PhyKIT run [r] on its [sites] columns ([model]: the model of the
/// IQ-TREE branch lengths).
SubsetSaturation phykitSubset(CodonSubset subset, int sites, PhykitSaturation r, SaturationThresholds th,
    {String model = 'GTR+G', int maxPoints = 3000, int variableSites = 0, int informativeSites = 0}) {
  final n = r.pairs.length;
  final meanP = n == 0 ? 0.0 : r.pairs.fold(0.0, (a, x) => a + x.p) / n;
  final meanD = n == 0 ? 0.0 : r.pairs.fold(0.0, (a, x) => a + x.patristic) / n;
  var points = [for (final x in r.pairs) SaturationPoint(x.patristic, x.p, 0, 0)];
  if (points.length > maxPoints) {
    final step = points.length / maxPoints;
    points = [for (var i = 0; i < maxPoints; i++) points[(i * step).floor()]];
  }
  return SubsetSaturation(
    subset: subset,
    sites: sites,
    pairs: n,
    undefinedPairs: 0,
    meanP: meanP,
    meanDistance: meanD,
    slope: r.slope,
    r2: null,
    tsTvRatio: null,
    iss: null,
    baseFrequencies: const [0.25, 0.25, 0.25, 0.25],
    model: model,
    points: points,
    verdict: phykitVerdict(th, sites: sites, pairs: n, meanP: meanP, slope: r.slope),
    variableSites: variableSites,
    informativeSites: informativeSites,
  );
}

/// p of the 1st, 2nd and 3rd positions of each pair against its patristic
/// distance in the tree of all sites (for the positions plot).
List<PositionPair> phykitPositionPairs(PhykitSaturation all, PhykitSaturation p1, PhykitSaturation p2, PhykitSaturation p3) {
  String key(String a, String b) => a.compareTo(b) < 0 ? '$a\t$b' : '$b\t$a';
  Map<String, double> of(PhykitSaturation r) => {for (final x in r.pairs) key(x.a, x.b): x.p};
  final m1 = of(p1), m2 = of(p2), m3 = of(p3);
  return [
    for (final x in all.pairs)
      if (m1[key(x.a, x.b)] case final a?)
        if (m2[key(x.a, x.b)] case final b?)
          if (m3[key(x.a, x.b)] case final c?) PositionPair(x.patristic, a, b, c),
  ];
}

/// Default thresholds for PhyKIT's saturation value (heuristic; the 2nd
/// positions of coding genes sit near 1 and saturated 3rd positions well
/// below 0.5).
const phykitThresholds = SaturationThresholds(slopeModerate: 0.5, slopeSaturated: 0.3);

/// Columns of [length] in [subset] for reading frame [frame].
List<int> subsetColumns(CodonSubset subset, int length, int frame) => [
      for (var c = 0; c < length; c++)
        if (subset == CodonSubset.all || subset.positions.contains(codonPositionOf(c, frame))) c,
    ];

/// Variable and parsimony-informative columns of [seqs] among [cols]
/// (bases only; gaps and ambiguities ignored).
(int, int) variableInformative(List<String> seqs, List<int> cols) {
  var variable = 0, informative = 0;
  for (final c in cols) {
    final counts = <int, int>{};
    for (final s in seqs) {
      if (c >= s.length) continue;
      final ch = s.codeUnitAt(c) | 0x20;
      if (ch == 0x61 || ch == 0x63 || ch == 0x67 || ch == 0x74 || ch == 0x75) counts.update(ch, (v) => v + 1, ifAbsent: () => 1);
    }
    if (counts.length > 1) variable++;
    if (counts.values.where((v) => v >= 2).length >= 2) informative++;
  }
  return (variable, math.min(variable, informative));
}
