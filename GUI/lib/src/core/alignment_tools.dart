/// Alignment post-processing that must stay codon-aware.
library;

import 'fasta.dart';

/// Parses the `#ColumnsMap\t0, 1, 2, ...` line that trimAl prints with
/// `-colnumbering` (0-based indices of the kept columns).
List<int>? parseTrimalColumnMap(String stdout) {
  final m = RegExp(r'#ColumnsMap\s+([\d,\s]*)').firstMatch(stdout);
  if (m == null) return null;
  return m.group(1)!
      .split(',')
      .map((s) => int.tryParse(s.trim()))
      .whereType<int>()
      .toList();
}

/// Columns of [full] kept in [trimmed] (the same sequences in the same
/// order), for tools that do not report them (ClipKIT): each column of
/// [trimmed] is matched to the next identical column of [full]. Trimmers
/// decide from the column content alone, so identical columns share their
/// fate and the first match is the kept one. Null when [trimmed] is not a
/// column subset of [full].
List<int>? keptColumns(List<SeqRecord> full, List<SeqRecord> trimmed) {
  if (full.isEmpty || full.length != trimmed.length) return null;
  final a = [for (final r in full) r.sequence.toUpperCase()];
  final b = [for (final r in trimmed) r.sequence.toUpperCase()];
  final n = a.first.length, m = b.first.length;
  if (a.any((s) => s.length != n) || b.any((s) => s.length != m)) return null;
  bool same(int i, int j) {
    for (var r = 0; r < a.length; r++) {
      if (a[r].codeUnitAt(i) != b[r].codeUnitAt(j)) return false;
    }
    return true;
  }

  final out = <int>[];
  var i = 0;
  for (var j = 0; j < m; j++) {
    while (i < n && !same(i, j)) {
      i++;
    }
    if (i >= n) return null;
    out.add(i++);
  }
  return out;
}

/// 0-based columns kept by ClipKIT, from its `-l` log (one line per
/// position: `<1-based position> <keep|trim> <site class> <gappyness>`).
List<int>? parseClipkitLog(String text) {
  final out = <int>[];
  var lines = 0;
  for (final line in text.split('\n')) {
    final f = line.trim().split(RegExp(r'\s+'));
    if (f.length < 2) continue;
    final pos = int.tryParse(f[0]);
    if (pos == null || (f[1] != 'keep' && f[1] != 'trim')) continue;
    lines++;
    if (f[1] == 'keep') out.add(pos - 1);
  }
  return lines == 0 ? null : (out..sort());
}

/// Sequences of a DNA alignment whose gaps are not multiples of three
/// inside the sequence: from such a gap on, the bases sit in the wrong
/// codon position of the column (a frame shift). For each sequence the
/// offset between column and base index (mod 3) should stay constant; a
/// sequence is reported when at least one codon's worth of bases (3) is off
/// the dominant offset.
List<String> frameShiftedSequences(List<SeqRecord> alignment) {
  final out = <String>[];
  for (final r in alignment) {
    final counts = [0, 0, 0];
    var i = 0;
    final s = r.sequence;
    for (var c = 0; c < s.length; c++) {
      final ch = s.codeUnitAt(c);
      if (ch == 0x2d || ch == 0x2e || ch == 0x3f) continue; // - . ?
      counts[(c - i) % 3]++;
      i++;
    }
    final total = counts[0] + counts[1] + counts[2];
    final dominant = [counts[0], counts[1], counts[2]].reduce((a, b) => a > b ? a : b);
    if (total - dominant >= 3) out.add(r.header);
  }
  return out;
}

/// MAFFT `--adjustdirection` prefixes reversed sequences with `_R_`.
({List<SeqRecord> records, List<String> reversed}) stripMafftReversedPrefix(List<SeqRecord> records) {
  final reversed = <String>[];
  final out = <SeqRecord>[];
  for (final r in records) {
    if (r.header.startsWith('_R_')) {
      final name = r.header.substring(3);
      reversed.add(name);
      out.add(r.copyWith(header: name));
    } else {
      out.add(r);
    }
  }
  return (records: out, reversed: reversed);
}

class AlignmentStats {
  AlignmentStats({
    required this.taxa,
    required this.length,
    required this.gapFraction,
    required this.variableSites,
    required this.informativeSites,
  });

  final int taxa;
  final int length;

  /// Fraction of `-`, `?`, `N` cells.
  final double gapFraction;
  final int variableSites;

  /// Parsimony-informative: at least two states each present in ≥ 2 taxa.
  final int informativeSites;
}

AlignmentStats alignmentStats(List<SeqRecord> records) {
  if (records.isEmpty) {
    return AlignmentStats(taxa: 0, length: 0, gapFraction: 0, variableSites: 0, informativeSites: 0);
  }
  final length = records.first.length;
  var missing = 0, variable = 0, informative = 0;
  final counts = <int, int>{};
  for (var c = 0; c < length; c++) {
    counts.clear();
    for (final r in records) {
      if (c >= r.sequence.length) continue;
      final ch = r.sequence.codeUnitAt(c) | 0x20; // lower-case
      if (ch == 0x61 || ch == 0x63 || ch == 0x67 || ch == 0x74 || ch == 0x75) {
        counts.update(ch == 0x75 ? 0x74 : ch, (v) => v + 1, ifAbsent: () => 1);
      } else {
        missing++;
      }
    }
    if (counts.length > 1) {
      variable++;
      if (counts.values.where((v) => v >= 2).length >= 2) informative++;
    }
  }
  return AlignmentStats(
    taxa: records.length,
    length: length,
    gapFraction: missing / (records.length * length),
    variableSites: variable,
    informativeSites: informative,
  );
}
