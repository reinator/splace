/// Outgroup suggestion from the alignments: pairwise p-distances over the
/// sites two taxa share in every locus, an average-linkage (UPGMA) tree of
/// those distances, and the smaller of the two clusters joined last, which
/// is the group farthest from the rest.
library;

import 'dart:math' as math;
import 'dart:typed_data';

import 'supermatrix.dart';

class OutgroupSuggestion {
  const OutgroupSuggestion({
    required this.taxa,
    required this.meanDistance,
    required this.between,
    required this.within,
    required this.confident,
    required this.loci,
    required this.taxonCount,
  });

  /// Suggested outgroup, farthest first.
  final List<String> taxa;

  /// Mean p-distance of each taxon to all the others.
  final Map<String, double> meanDistance;

  /// Mean distance between the suggested group and the rest, and within
  /// the rest (the ingroup).
  final double between, within;

  /// The deepest split leaves a small group apart (false: the taxa are
  /// about equally distant and the most distant taxon is only a guess).
  final bool confident;
  final int loci, taxonCount;

  /// Taxa ranked by their mean distance to the others, farthest first.
  List<String> get ranked => meanDistance.keys.toList()..sort((a, b) => meanDistance[b]!.compareTo(meanDistance[a]!));
}

/// Residue code: 0 for gaps and unknown states, so they are not compared.
Uint8List _encode(String s, bool dna) {
  final out = Uint8List(s.length);
  for (var i = 0; i < s.length; i++) {
    var c = s.codeUnitAt(i);
    if (c >= 0x61 && c <= 0x7a) c -= 0x20;
    if (dna) {
      out[i] = switch (c) { 0x41 => 1, 0x43 => 2, 0x47 => 3, 0x54 || 0x55 => 4, _ => 0 };
    } else {
      out[i] = c >= 0x41 && c <= 0x5a && c != 0x58 && c != 0x42 && c != 0x5a && c != 0x4a ? c : 0;
    }
  }
  return out;
}

/// [maxSize]: a group joined last larger than this (or than a third of the
/// taxa) is not an outgroup; the most distant taxon is offered instead.
OutgroupSuggestion? suggestOutgroup(List<LocusAlignment> loci, {int maxSize = 6, int minSites = 20}) {
  final names = <String>[];
  final index = <String, int>{};
  for (final l in loci) {
    for (final r in l.records) {
      index.putIfAbsent(r.header, () {
        names.add(r.header);
        return names.length - 1;
      });
    }
  }
  final n = names.length;
  if (n < 3) return null;
  final diff = Float64List(n * n), comp = Float64List(n * n);
  var used = 0;
  for (final l in loci) {
    if (l.records.length < 2) continue;
    final sample = l.records.first.sequence.toUpperCase().replaceAll(RegExp('[^A-Z]'), '');
    final acgt = sample.replaceAll(RegExp('[^ACGTUN]'), '').length;
    final dna = sample.isEmpty || acgt / sample.length > 0.9;
    final ids = [for (final r in l.records) index[r.header]!];
    final seqs = [for (final r in l.records) _encode(r.sequence, dna)];
    used++;
    for (var a = 0; a < seqs.length; a++) {
      final sa = seqs[a];
      for (var b = a + 1; b < seqs.length; b++) {
        final sb = seqs[b];
        final len = math.min(sa.length, sb.length);
        var d = 0, c = 0;
        for (var k = 0; k < len; k++) {
          final x = sa[k], y = sb[k];
          if (x == 0 || y == 0) continue;
          c++;
          if (x != y) d++;
        }
        final i = ids[a], j = ids[b];
        diff[i * n + j] += d;
        comp[i * n + j] += c;
        diff[j * n + i] += d;
        comp[j * n + i] += c;
      }
    }
  }
  // Distances (NaN: the two taxa share too few sites).
  final dist = Float64List(n * n);
  for (var i = 0; i < n * n; i++) {
    dist[i] = comp[i] >= minSites ? diff[i] / comp[i] : double.nan;
  }
  final mean = <String, double>{};
  for (var i = 0; i < n; i++) {
    var s = 0.0, k = 0;
    for (var j = 0; j < n; j++) {
      if (i == j || dist[i * n + j].isNaN) continue;
      s += dist[i * n + j];
      k++;
    }
    mean[names[i]] = k == 0 ? 0 : s / k;
  }

  // UPGMA with the sums and counts of known distances between clusters.
  final clusters = [for (var i = 0; i < n; i++) <int>[i]];
  final sum = <int, Map<int, double>>{}, cnt = <int, Map<int, int>>{};
  for (var i = 0; i < n; i++) {
    sum[i] = {};
    cnt[i] = {};
    for (var j = 0; j < n; j++) {
      if (i == j) continue;
      final v = dist[i * n + j];
      sum[i]![j] = v.isNaN ? 0 : v;
      cnt[i]![j] = v.isNaN ? 0 : 1;
    }
  }
  final alive = {for (var i = 0; i < n; i++) i};
  double avg(int a, int b) => cnt[a]![b]! == 0 ? double.infinity : sum[a]![b]! / cnt[a]![b]!;
  while (alive.length > 2) {
    var best = double.infinity;
    int? ba, bb;
    final list = alive.toList();
    for (var x = 0; x < list.length; x++) {
      for (var y = x + 1; y < list.length; y++) {
        final v = avg(list[x], list[y]);
        if (v < best) {
          best = v;
          ba = list[x];
          bb = list[y];
        }
      }
    }
    if (ba == null) {
      // Nothing comparable left: join the first two.
      ba = list[0];
      bb = list[1];
    }
    final a = ba, b = bb!;
    clusters[a].addAll(clusters[b]);
    alive.remove(b);
    for (final k in alive) {
      if (k == a) continue;
      final s = sum[a]![k]! + sum[b]![k]!;
      final c = cnt[a]![k]! + cnt[b]![k]!;
      sum[a]![k] = s;
      sum[k]![a] = s;
      cnt[a]![k] = c;
      cnt[k]![a] = c;
    }
  }
  final last = alive.toList();
  final g1 = clusters[last[0]], g2 = clusters[last[1]];
  final small = g1.length <= g2.length ? g1 : g2;
  final confident = n >= 4 && small.length <= maxSize && small.length <= (n / 3).ceil();
  final ranked = mean.keys.toList()..sort((a, b) => mean[b]!.compareTo(mean[a]!));
  final pick = confident ? ([for (final i in small) names[i]]..sort((a, b) => mean[b]!.compareTo(mean[a]!))) : [ranked.first];
  final pickIdx = {for (final t in pick) index[t]!};
  double meanOf(Iterable<int> xs, Iterable<int> ys, {bool same = false}) {
    var s = 0.0, k = 0;
    for (final i in xs) {
      for (final j in ys) {
        if (i == j || (same && j < i)) continue;
        final v = dist[i * n + j];
        if (v.isNaN) continue;
        s += v;
        k++;
      }
    }
    return k == 0 ? 0 : s / k;
  }

  final ingroup = [for (var i = 0; i < n; i++) if (!pickIdx.contains(i)) i];
  return OutgroupSuggestion(
    taxa: pick,
    meanDistance: mean,
    between: meanOf(pickIdx, ingroup),
    within: meanOf(ingroup, ingroup, same: true),
    confident: confident,
    loci: used,
    taxonCount: n,
  );
}
