/// Substitution saturation tests per codon position.
///
/// For each codon-position subset (1st, 2nd, 3rd, 1st+2nd, all) this module
/// computes, over all pairs of sequences, the uncorrected p-distance, the
/// proportions of transitions (Ts) and transversions (Tv), and the TN93
/// distance (K2P when base frequencies are degenerate). It summarizes them
/// as:
///
/// * `slope`: least-squares slope through the origin of p on the corrected
///   distance. ~1 means observed differences capture nearly all inferred
///   substitutions; values well below 1 mean multiple hits are hidden.
/// * `undefinedFraction`: pairs whose corrected distance cannot be computed
///   (log of a non-positive number), i.e. divergence beyond what the model
///   can correct — a strong saturation signal.
/// * `iss`: Xia et al. (2003) index of substitution saturation,
///   Iss = mean site entropy / expected entropy under full saturation. The
///   expectation is computed for the actual number of unambiguous sequences
///   at each site (multinomial sampling), so small data sets are not biased
///   towards "unsaturated". The critical value Iss.c from Xia's simulations
///   is NOT computed here; compare with DAMBE if a formal test is needed.
///
/// Complementary analyses, per subset:
///
/// * Composition homogeneity: chi-square test of the base frequencies of
///   each sequence against the whole subset (as IQ-TREE reports it, df = 3)
///   and the contingency test across all taxa (PAUP*). Saturated third
///   positions often drift in composition (e.g. GC3), which also biases
///   trees.
/// * Variable and parsimony-informative sites, and GC content.
/// * Per pair of sequences, the p-distance of the 1st, 2nd and 3rd positions
///   against the corrected distance of all positions ([PositionPair]): a
///   position that levels off while the others still grow is saturated.
///
/// Verdicts use configurable heuristic thresholds ([SaturationThresholds]).
library;

import 'dart:math' as math;
import 'dart:typed_data';

import 'fasta.dart';
import 'genetic_code.dart';

enum CodonSubset { pos1, pos2, pos3, pos12, all }

extension CodonSubsetX on CodonSubset {
  /// Codon positions (1-based) included in this subset.
  Set<int> get positions => switch (this) {
        CodonSubset.pos1 => {1},
        CodonSubset.pos2 => {2},
        CodonSubset.pos3 => {3},
        CodonSubset.pos12 => {1, 2},
        CodonSubset.all => {1, 2, 3},
      };

  String get label => switch (this) {
        CodonSubset.pos1 => '1st',
        CodonSubset.pos2 => '2nd',
        CodonSubset.pos3 => '3rd',
        CodonSubset.pos12 => '1st+2nd',
        CodonSubset.all => 'All',
      };
}

/// Codon position (1, 2 or 3) of alignment column [column] when the first
/// complete codon starts at [frame] (0, 1 or 2).
int codonPositionOf(int column, int frame) => (column - frame) % 3 + 1;

const _a = 0, _c = 1, _g = 2, _t = 3, _other = 4;

Uint8List encodeNucleotides(String seq) {
  final out = Uint8List(seq.length);
  for (var i = 0; i < seq.length; i++) {
    out[i] = switch (seq.codeUnitAt(i) | 0x20) {
      0x61 => _a, // a
      0x63 => _c, // c
      0x67 => _g, // g
      0x74 || 0x75 => _t, // t, u
      _ => _other,
    };
  }
  return out;
}

bool _isTransition(int x, int y) =>
    (x == _a && y == _g) || (x == _g && y == _a) || (x == _c && y == _t) || (x == _t && y == _c);

class SaturationPoint {
  const SaturationPoint(this.distance, this.p, this.ts, this.tv);

  /// Corrected (TN93/K2P) distance.
  final double distance;

  /// Uncorrected p-distance.
  final double p;

  /// Proportion of sites differing by a transition / transversion.
  final double ts;
  final double tv;
}

enum SaturationVerdict { ok, moderate, saturated, insufficient }

class SaturationThresholds {
  const SaturationThresholds({
    this.slopeModerate = 0.75,
    this.slopeSaturated = 0.5,
    this.undefinedModerate = 0.02,
    this.undefinedSaturated = 0.10,
    this.issSaturated = 0.80,
    this.minPairs = 3,
    this.minSites = 30,
    this.lowDivergence = 0.03,
  });

  /// Slope of p on corrected distance below which a subset is flagged.
  /// Under JC, slope ≈ 0.75 corresponds to p ≈ 0.4 and 0.5 to p ≈ 0.6.
  final double slopeModerate;
  final double slopeSaturated;
  final double undefinedModerate;
  final double undefinedSaturated;

  /// Iss at or above which sites look close to random.
  final double issSaturated;
  final int minPairs;
  final int minSites;

  /// Mean p below this is too little divergence for saturation to matter.
  final double lowDivergence;
}

class SubsetSaturation {
  SubsetSaturation({
    required this.subset,
    required this.sites,
    required this.pairs,
    required this.undefinedPairs,
    required this.meanP,
    required this.meanDistance,
    required this.slope,
    required this.r2,
    required this.tsTvRatio,
    required this.iss,
    required this.baseFrequencies,
    required this.model,
    required this.points,
    required this.verdict,
    this.composition,
    this.variableSites = 0,
    this.informativeSites = 0,
  });

  final CodonSubset subset;
  final int sites;

  /// Pairs with enough comparable sites.
  final int pairs;
  final int undefinedPairs;
  final double meanP;
  final double meanDistance;
  final double? slope;
  final double? r2;
  final double? tsTvRatio;
  final double? iss;

  /// A, C, G, T frequencies.
  final List<double> baseFrequencies;
  final String model;

  /// Possibly subsampled points for plotting; only pairs with a defined
  /// corrected distance.
  final List<SaturationPoint> points;
  final SaturationVerdict verdict;

  /// Base-composition homogeneity across taxa.
  final CompositionTest? composition;
  final int variableSites;
  final int informativeSites;

  /// G + C frequency.
  double get gc => baseFrequencies[_c] + baseFrequencies[_g];

  double get undefinedFraction => pairs == 0 ? 0 : undefinedPairs / pairs;

  /// 0 (none) .. 1 (complete): share of inferred substitutions not visible
  /// as differences, used for heatmap colouring.
  double get saturationIndex {
    if (pairs == 0) return 0;
    final fromSlope = slope == null ? 0.0 : (1 - slope!).clamp(0.0, 1.0);
    final fromUndefined = (undefinedFraction / 0.2).clamp(0.0, 1.0);
    return math.max(fromSlope, fromUndefined);
  }
}

/// Chi-square tests of base-composition homogeneity.
class CompositionTest {
  const CompositionTest({required this.chi2, required this.df, required this.p, required this.failing, required this.taxa});

  /// Contingency chi-square across all taxa (bases × taxa).
  final double chi2;
  final int df;
  final double p;

  /// Taxa whose composition differs from the subset's (per-sequence test,
  /// df = 3, p < 0.05).
  final int failing;
  final int taxa;
}

/// One pair of sequences: corrected distance over all codon positions
/// and the p-distance of each position.
class PositionPair {
  const PositionPair(this.distance, this.p1, this.p2, this.p3);
  final double distance;
  final double p1, p2, p3;
}

class SaturationReport {
  SaturationReport(this.name, this.taxa, this.length, this.frame, this.subsets, this.frameStops, [this.positionPairs = const []]);
  final String name;
  final int taxa;
  final int length;
  final int frame;
  final Map<CodonSubset, SubsetSaturation> subsets;

  /// p of each codon position against the all-positions distance.
  final List<PositionPair> positionPairs;

  /// Internal stop codons counted for each candidate frame (when a genetic
  /// code was supplied for frame detection).
  final List<int>? frameStops;

  CodonRecommendation get recommendation => recommendCodonPositions(this);
}

class SaturationOptions {
  const SaturationOptions({
    this.frame,
    this.geneticCode = 1,
    this.thresholds = const SaturationThresholds(),
    this.maxPlotPoints = 3000,
    this.subsets = CodonSubset.values,
    this.issReplicates = 400,
  });

  /// Random sequences simulated per site size for the expected entropy of
  /// Xia's Iss under full saturation.
  final int issReplicates;

  /// 0, 1, 2 or null to detect the frame from stop codons.
  final int? frame;
  final int geneticCode;
  final SaturationThresholds thresholds;
  final int maxPlotPoints;
  final List<CodonSubset> subsets;
}

/// Runs the saturation analysis on an aligned set of nucleotide sequences.
///
/// [onProgress] gets the fraction done after each stage (frame, each codon
/// subset, the position pairs of the plot).
SaturationReport testSaturation(String name, List<SeqRecord> alignment,
    [SaturationOptions options = const SaturationOptions(), void Function(double fraction)? onProgress]) {
  if (alignment.isEmpty) {
    return SaturationReport(name, 0, 0, options.frame ?? 0, {}, null);
  }
  if (!isAligned(alignment)) {
    throw ArgumentError('Sequences in "$name" are not aligned (different lengths).');
  }
  final length = alignment.first.length;
  List<int>? stops;
  var frame = options.frame;
  if (frame == null) {
    stops = [for (var f = 0; f < 3; f++) countInternalStops(alignment, f, options.geneticCode)];
    frame = 0;
    for (var f = 1; f < 3; f++) {
      if (stops[f] < stops[frame!]) frame = f;
    }
  }
  final encoded = [for (final r in alignment) encodeNucleotides(r.sequence)];
  // Codon positions against the whole gene only make sense with codons.
  final codons = options.subsets.any((s) => s != CodonSubset.all);
  final stages = options.subsets.length + (codons ? 1 : 0) + 1;
  var stage = 1;
  onProgress?.call(stage / stages);
  final subsets = <CodonSubset, SubsetSaturation>{};
  for (final subset in options.subsets) {
    final columns = [
      for (var c = 0; c < length; c++)
        if (subset.positions.contains(codonPositionOf(c, frame!))) c,
    ];
    subsets[subset] = _analyzeSubset(subset, encoded, columns, options);
    onProgress?.call(++stage / stages);
  }
  final pairs = codons ? _positionPairs(encoded, length, frame!, options.maxPlotPoints) : const <PositionPair>[];
  onProgress?.call(1);
  return SaturationReport(name, alignment.length, length, frame!, subsets, stops, pairs);
}

/// Per pair: TN93/K2P distance of all positions and p of each position.
List<PositionPair> _positionPairs(List<Uint8List> seqs, int length, int frame, int maxPoints) {
  final n = seqs.length;
  final pos = Uint8List(length);
  for (var c = 0; c < length; c++) {
    pos[c] = codonPositionOf(c, frame) - 1;
  }
  final counts = List<int>.filled(5, 0);
  for (final s in seqs) {
    for (final b in s) {
      counts[b]++;
    }
  }
  final total = counts[_a] + counts[_c] + counts[_g] + counts[_t];
  final freqs = total == 0 ? [0.25, 0.25, 0.25, 0.25] : [for (var b = 0; b < 4; b++) counts[b] / total];
  final useTn93 = freqs.every((f) => f > 0.001);
  final out = <PositionPair>[];
  final valid = List<int>.filled(3, 0), diff = List<int>.filled(3, 0);
  for (var i = 0; i < n; i++) {
    for (var j = i + 1; j < n; j++) {
      final si = seqs[i], sj = seqs[j];
      valid.fillRange(0, 3, 0);
      diff.fillRange(0, 3, 0);
      var all = 0, ag = 0, ct = 0, tv = 0;
      for (var c = 0; c < length; c++) {
        final x = si[c], y = sj[c];
        if (x == _other || y == _other) continue;
        final k = pos[c];
        valid[k]++;
        all++;
        if (x == y) continue;
        diff[k]++;
        if (_isTransition(x, y)) {
          if (x == _a || x == _g) {
            ag++;
          } else {
            ct++;
          }
        } else {
          tv++;
        }
      }
      if (all < 30 || valid.any((v) => v < 10)) continue;
      final d = useTn93 ? tn93Distance(ag / all, ct / all, tv / all, freqs) : k2pDistance((ag + ct) / all, tv / all);
      if (d == null) continue;
      out.add(PositionPair(d, diff[0] / valid[0], diff[1] / valid[1], diff[2] / valid[2]));
    }
  }
  if (out.length <= maxPoints) return out;
  return (List.of(out)..shuffle(math.Random(7))).sublist(0, maxPoints);
}

/// Chi-square tests of base composition over [cols] (see [CompositionTest]).
CompositionTest? compositionTest(List<Uint8List> seqs, Int32List cols) {
  final perTaxon = <List<int>>[];
  final totals = List<int>.filled(4, 0);
  for (final s in seqs) {
    final c = List<int>.filled(4, 0);
    for (var k = 0; k < cols.length; k++) {
      final b = s[cols[k]];
      if (b != _other) c[b]++;
    }
    if (c.reduce((a, b) => a + b) == 0) continue;
    perTaxon.add(c);
    for (var b = 0; b < 4; b++) {
      totals[b] += c[b];
    }
  }
  final grand = totals.reduce((a, b) => a + b);
  if (perTaxon.length < 2 || grand == 0) return null;
  final present = [for (var b = 0; b < 4; b++) if (totals[b] > 0) b];
  if (present.length < 2) return null;
  var chi2 = 0.0;
  var failing = 0;
  for (final c in perTaxon) {
    final n = c.reduce((a, b) => a + b);
    var own = 0.0;
    for (final b in present) {
      final e = n * totals[b] / grand;
      if (e <= 0) continue;
      final d = c[b] - e;
      own += d * d / e;
    }
    chi2 += own;
    if (chiSquareP(own, present.length - 1) < 0.05) failing++;
  }
  final df = (present.length - 1) * (perTaxon.length - 1);
  return CompositionTest(chi2: chi2, df: df, p: chiSquareP(chi2, df), failing: failing, taxa: perTaxon.length);
}

/// Upper tail of the chi-square distribution: P(X ≥ [x]) with [df]
/// degrees of freedom (regularized upper incomplete gamma Q(df/2, x/2)).
double chiSquareP(double x, int df) {
  if (df <= 0) return 1;
  if (x <= 0) return 1;
  return _gammaQ(df / 2, x / 2);
}

double _gammaQ(double a, double x) {
  final gln = _lnGamma(a);
  if (x < a + 1) {
    // Series for P, then Q = 1 - P.
    var ap = a, sum = 1 / a, del = sum;
    for (var n = 0; n < 500; n++) {
      ap++;
      del *= x / ap;
      sum += del;
      if (del.abs() < sum.abs() * 1e-12) break;
    }
    return (1 - sum * math.exp(-x + a * math.log(x) - gln)).clamp(0.0, 1.0);
  }
  // Continued fraction for Q (Lentz).
  const tiny = 1e-300;
  var b = x + 1 - a, c = 1 / tiny, d = 1 / b, h = d;
  for (var i = 1; i < 500; i++) {
    final an = -i * (i - a);
    b += 2;
    d = an * d + b;
    if (d.abs() < tiny) d = tiny;
    c = b + an / c;
    if (c.abs() < tiny) c = tiny;
    d = 1 / d;
    final del = d * c;
    h *= del;
    if ((del - 1).abs() < 1e-12) break;
  }
  return (math.exp(-x + a * math.log(x) - gln) * h).clamp(0.0, 1.0);
}

double _lnGamma(double x) {
  const cof = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
  var y = x;
  final tmp = x + 5.5 - (x + 0.5) * math.log(x + 5.5);
  var ser = 1.000000000190015;
  for (final c in cof) {
    ser += c / ++y;
  }
  return -tmp + math.log(2.5066282746310005 * ser / x);
}

/// Variable and parsimony-informative columns among [cols].
(int, int) variableSites(List<Uint8List> seqs, Int32List cols) {
  var variable = 0, informative = 0;
  final c = List<int>.filled(4, 0);
  for (var k = 0; k < cols.length; k++) {
    c.fillRange(0, 4, 0);
    for (final s in seqs) {
      final b = s[cols[k]];
      if (b != _other) c[b]++;
    }
    final states = c.where((v) => v > 0).length;
    if (states >= 2) variable++;
    if (c.where((v) => v >= 2).length >= 2) informative++;
  }
  return (variable, informative);
}

SubsetSaturation _analyzeSubset(
  CodonSubset subset,
  List<Uint8List> seqs,
  List<int> columns,
  SaturationOptions options,
) {
  final th = options.thresholds;
  final n = seqs.length;
  final m = columns.length;
  final cols = Int32List.fromList(columns);

  // Base frequencies over the subset.
  final counts = List<int>.filled(5, 0);
  for (final s in seqs) {
    for (var k = 0; k < m; k++) {
      counts[s[cols[k]]]++;
    }
  }
  final totalBases = counts[_a] + counts[_c] + counts[_g] + counts[_t];
  final freqs = totalBases == 0
      ? [0.25, 0.25, 0.25, 0.25]
      : [for (var b = 0; b < 4; b++) counts[b] / totalBases];
  final useTn93 = freqs.every((f) => f > 0.001);

  final allPoints = <SaturationPoint>[];
  var pairs = 0, undefined = 0;
  var sumP = 0.0, sumD = 0.0, sumPD = 0.0, sumDD = 0.0, sumPP = 0.0;
  var sumTs = 0.0, sumTv = 0.0;

  // Need enough comparable sites to estimate a distance.
  final minComparable = math.max(10, (m * 0.2).floor());

  for (var i = 0; i < n; i++) {
    final si = seqs[i];
    for (var j = i + 1; j < n; j++) {
      final sj = seqs[j];
      var valid = 0, ag = 0, ct = 0, tv = 0;
      for (var k = 0; k < m; k++) {
        final x = si[cols[k]], y = sj[cols[k]];
        if (x == _other || y == _other) continue;
        valid++;
        if (x == y) continue;
        if (_isTransition(x, y)) {
          if (x == _a || x == _g) {
            ag++;
          } else {
            ct++;
          }
        } else {
          tv++;
        }
      }
      if (valid < minComparable) continue;
      pairs++;
      final p1 = ag / valid, p2 = ct / valid, q = tv / valid;
      final p = p1 + p2 + q;
      final d = useTn93 ? tn93Distance(p1, p2, q, freqs) : k2pDistance(p1 + p2, q);
      sumTs += p1 + p2;
      sumTv += q;
      if (d == null) {
        undefined++;
        continue;
      }
      allPoints.add(SaturationPoint(d, p, p1 + p2, q));
      sumP += p;
      sumD += d;
      sumPD += p * d;
      sumDD += d * d;
      sumPP += p * p;
    }
  }

  final defined = allPoints.length;
  double? slope, r2;
  if (defined >= th.minPairs && sumDD > 0) {
    slope = sumPD / sumDD;
    final cov = sumPD / defined - (sumP / defined) * (sumD / defined);
    final varP = sumPP / defined - math.pow(sumP / defined, 2);
    final varD = sumDD / defined - math.pow(sumD / defined, 2);
    if (varP > 0 && varD > 0) r2 = (cov * cov) / (varP * varD);
  }
  final meanP = defined == 0 ? 0.0 : sumP / defined;
  final iss = xiaIss(seqs, cols, freqs, replicates: options.issReplicates);
  final (variable, informative) = variableSites(seqs, cols);

  final verdict = _verdict(
    th,
    sites: m,
    pairs: pairs,
    undefinedFraction: pairs == 0 ? 0 : undefined / pairs,
    meanP: meanP,
    slope: slope,
    iss: iss,
  );

  return SubsetSaturation(
    subset: subset,
    sites: m,
    pairs: pairs,
    undefinedPairs: undefined,
    meanP: meanP,
    meanDistance: defined == 0 ? 0 : sumD / defined,
    slope: slope,
    r2: r2,
    tsTvRatio: sumTv > 0 ? sumTs / sumTv : null,
    iss: iss,
    baseFrequencies: freqs,
    model: useTn93 ? 'TN93' : 'K2P',
    points: _subsample(allPoints, options.maxPlotPoints),
    verdict: verdict,
    composition: compositionTest(seqs, cols),
    variableSites: variable,
    informativeSites: informative,
  );
}

SaturationVerdict _verdict(
  SaturationThresholds th, {
  required int sites,
  required int pairs,
  required double undefinedFraction,
  required double meanP,
  required double? slope,
  required double? iss,
}) {
  if (sites < th.minSites || pairs < th.minPairs) return SaturationVerdict.insufficient;
  if (undefinedFraction >= th.undefinedSaturated) return SaturationVerdict.saturated;
  if (slope == null) return SaturationVerdict.insufficient;
  if (meanP < th.lowDivergence) return SaturationVerdict.ok;
  if (slope < th.slopeSaturated) return SaturationVerdict.saturated;
  if (iss != null && iss >= th.issSaturated) return SaturationVerdict.saturated;
  if (slope < th.slopeModerate || undefinedFraction >= th.undefinedModerate) {
    return SaturationVerdict.moderate;
  }
  return SaturationVerdict.ok;
}

List<SaturationPoint> _subsample(List<SaturationPoint> points, int max) {
  if (points.length <= max) return points;
  final rnd = math.Random(7);
  final copy = List.of(points)..shuffle(rnd);
  return copy.sublist(0, max);
}

/// Kimura 2-parameter distance from transition (P) and transversion (Q)
/// proportions; null when undefined.
double? k2pDistance(double p, double q) {
  final w1 = 1 - 2 * p - q;
  final w2 = 1 - 2 * q;
  if (w1 <= 0 || w2 <= 0) return null;
  return -0.5 * math.log(w1) - 0.25 * math.log(w2);
}

/// Tamura–Nei (1993) distance. [p1]: A<->G proportion, [p2]: C<->T,
/// [q]: transversions, [pi]: A, C, G, T frequencies. Null when undefined.
double? tn93Distance(double p1, double p2, double q, List<double> pi) {
  final piA = pi[_a], piC = pi[_c], piG = pi[_g], piT = pi[_t];
  final gR = piA + piG, gY = piC + piT;
  final w1 = 1 - gR / (2 * piA * piG) * p1 - q / (2 * gR);
  final w2 = 1 - gY / (2 * piC * piT) * p2 - q / (2 * gY);
  final w3 = 1 - q / (2 * gR * gY);
  if (w1 <= 0 || w2 <= 0 || w3 <= 0) return null;
  return -2 * piA * piG / gR * math.log(w1) -
      2 * piC * piT / gY * math.log(w2) -
      2 * (gR * gY - piA * piG * gY / gR - piC * piT * gR / gY) * math.log(w3);
}

double _entropy(List<int> counts, int total) {
  var h = 0.0;
  for (final c in counts) {
    if (c == 0) continue;
    final f = c / total;
    h -= f * math.log(f) / math.ln2;
  }
  return h;
}

/// Xia's Iss. Gaps and ambiguous characters are excluded per site; sites
/// with fewer than 4 unambiguous characters are skipped. The full-saturation
/// entropy for a site with N characters is the expected entropy of a
/// multinomial sample of size N from [freqs] (Monte Carlo, fixed seed).
double? xiaIss(List<Uint8List> seqs, Int32List cols, List<double> freqs, {int replicates = 400}) {
  final expected = <int, double>{};
  final rnd = math.Random(2003);
  final cumulative = [freqs[0], freqs[0] + freqs[1], freqs[0] + freqs[1] + freqs[2]];

  double expectedEntropy(int n) => expected.putIfAbsent(n, () {
        var sum = 0.0;
        final c = [0, 0, 0, 0];
        for (var r = 0; r < replicates; r++) {
          c[0] = c[1] = c[2] = c[3] = 0;
          for (var i = 0; i < n; i++) {
            final u = rnd.nextDouble();
            c[u < cumulative[0] ? 0 : u < cumulative[1] ? 1 : u < cumulative[2] ? 2 : 3]++;
          }
          sum += _entropy(c, n);
        }
        return sum / replicates;
      });

  var sumH = 0.0, sumHfss = 0.0;
  var sites = 0;
  final counts = [0, 0, 0, 0];
  for (var k = 0; k < cols.length; k++) {
    counts[0] = counts[1] = counts[2] = counts[3] = 0;
    var n = 0;
    for (final s in seqs) {
      final b = s[cols[k]];
      if (b == _other) continue;
      counts[b]++;
      n++;
    }
    if (n < 4) continue;
    sumH += _entropy(counts, n);
    sumHfss += expectedEntropy(n);
    sites++;
  }
  if (sites == 0 || sumHfss == 0) return null;
  return sumH / sumHfss;
}

// ---------------------------------------------------------------------------
// Reading frame
// ---------------------------------------------------------------------------

/// Stop codons for an NCBI genetic code.
Set<String> stopCodons(int geneticCode) => stopCodonsFor(geneticCode);

/// Counts in-frame stop codons over all sequences, ignoring the last codon
/// of each sequence and codons containing gaps or ambiguities.
int countInternalStops(List<SeqRecord> alignment, int frame, int geneticCode) {
  final stops = stopCodons(geneticCode);
  var total = 0;
  for (final r in alignment) {
    final seq = r.sequence.toUpperCase().replaceAll('U', 'T');
    final ungappedEnd = seq.lastIndexOf(RegExp('[ACGT]'));
    for (var i = frame; i + 3 <= seq.length; i += 3) {
      if (i + 2 >= ungappedEnd) break; // terminal codon may be a stop
      final codon = seq.substring(i, i + 3);
      if (stops.contains(codon)) total++;
    }
  }
  return total;
}

// ---------------------------------------------------------------------------
// Recommendation
// ---------------------------------------------------------------------------

enum RecommendationKind {
  /// No position is saturated: keep everything, partition by codon position.
  keepAll,

  /// Remove the positions in [CodonRecommendation.exclude].
  excludePositions,

  /// All three positions saturated: nucleotides unreliable at this depth.
  useAminoAcids,

  /// Not enough data to decide.
  undetermined,
}

class CodonRecommendation {
  CodonRecommendation(this.kind, this.exclude, this.caution);
  final RecommendationKind kind;

  /// Codon positions (1..3) to exclude.
  final Set<int> exclude;

  /// Positions with moderate saturation: keep, but partition separately and
  /// test sensitivity (e.g. with and without them, or RY-coding).
  final Set<int> caution;
}

CodonRecommendation recommendCodonPositions(SaturationReport report) => recommendFromVerdicts({
      1: report.subsets[CodonSubset.pos1]?.verdict,
      2: report.subsets[CodonSubset.pos2]?.verdict,
      3: report.subsets[CodonSubset.pos3]?.verdict,
    });

/// What to do with codon positions 1–3 given their verdicts.
CodonRecommendation recommendFromVerdicts(Map<int, SaturationVerdict?> verdicts) {
  if (verdicts.values.every((v) => v == null || v == SaturationVerdict.insufficient)) {
    return CodonRecommendation(RecommendationKind.undetermined, {}, {});
  }
  final saturated = {for (final e in verdicts.entries) if (e.value == SaturationVerdict.saturated) e.key};
  final moderate = {for (final e in verdicts.entries) if (e.value == SaturationVerdict.moderate) e.key};
  if (saturated.length == 3) {
    return CodonRecommendation(RecommendationKind.useAminoAcids, {}, {});
  }
  if (saturated.isEmpty) {
    return CodonRecommendation(RecommendationKind.keepAll, {}, moderate);
  }
  return CodonRecommendation(RecommendationKind.excludePositions, saturated, moderate);
}

/// The coding genes tested one by one, summed up per codon position (the
/// genes are never concatenated for the test). A position counts as
/// saturated when it is saturated in more than half of the genes where it
/// could be assessed, as moderate when it is saturated in some of them or
/// saturated or moderate in more than half, and as unsaturated otherwise.
class PooledSaturation {
  PooledSaturation._(this.reports, this.genes);

  /// The coding genes' reports.
  final List<SaturationReport> reports;

  /// For each position (1–3): the genes of each verdict.
  final Map<int, Map<SaturationVerdict, List<String>>> genes;

  static PooledSaturation? of(Iterable<SaturationReport> all) {
    final coding = [for (final r in all) if (r.subsets.containsKey(CodonSubset.pos1)) r];
    if (coding.isEmpty) return null;
    final genes = <int, Map<SaturationVerdict, List<String>>>{};
    for (final (p, subset) in [(1, CodonSubset.pos1), (2, CodonSubset.pos2), (3, CodonSubset.pos3)]) {
      final m = genes[p] = <SaturationVerdict, List<String>>{};
      for (final r in coding) {
        final v = r.subsets[subset]?.verdict ?? SaturationVerdict.insufficient;
        m.putIfAbsent(v, () => []).add(r.name);
      }
    }
    return PooledSaturation._(coding, genes);
  }

  /// Genes where position [p] could be assessed.
  int assessed(int p) => genes[p]!.entries.where((e) => e.key != SaturationVerdict.insufficient).fold(0, (a, e) => a + e.value.length);

  int count(int p, SaturationVerdict v) => genes[p]![v]?.length ?? 0;

  SaturationVerdict verdict(int p) {
    final n = assessed(p);
    if (n == 0) return SaturationVerdict.insufficient;
    final sat = count(p, SaturationVerdict.saturated), mod = count(p, SaturationVerdict.moderate);
    if (2 * sat > n) return SaturationVerdict.saturated;
    // Saturated in some genes (not most), or moderate in most: kept, with
    // the caution of a moderate position.
    if (sat > 0 || 2 * (sat + mod) > n) return SaturationVerdict.moderate;
    return SaturationVerdict.ok;
  }

  CodonRecommendation get recommendation => recommendFromVerdicts({for (final p in [1, 2, 3]) p: verdict(p)});
}
