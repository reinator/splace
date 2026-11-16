import 'dart:async';
import 'dart:math' as math;

import 'package:flutter/material.dart';

import '../../state/project_controller.dart' show SeqMode;
import '../common.dart';
import '../fa.dart';
import '../theme.dart';
import 'alignment_viewer.dart' show aaColors, ntColors, residueInk;

/// Small drawn examples that explain an option at a glance (like the
/// partition pictures of step 7): a toy alignment, codon positions or a
/// saturation curve. They use plain widgets, so they also work inside
/// IntrinsicHeight rows.

/// Codon position colours (the same as the partition pictures).
const codonPositionColors = [SplaceBlue.c300, SplacePalette.orange, SplacePalette.red700];

class _Caption extends StatelessWidget {
  const _Caption(this.text);
  final String text;

  @override
  Widget build(BuildContext context) =>
      Text(text, style: TextStyle(fontSize: 11.5, color: context.tokens.textMuted, height: 1.35));
}

/// A frame with the "Example" label.
class ExampleFrame extends StatelessWidget {
  const ExampleFrame({super.key, required this.child, this.caption, this.label, this.trailing});
  final Widget child;
  final String? caption;
  final String? label;

  /// Small control at the right of the label (e.g. replay).
  final Widget? trailing;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    // Spans the card it is in (also when its content is narrower).
    return Container(
      width: double.infinity,
      padding: const EdgeInsets.fromLTRB(10, 8, 10, 10),
      decoration: BoxDecoration(color: t.surface2, borderRadius: BorderRadius.circular(8), border: Border.all(color: t.border)),
      child: Column(crossAxisAlignment: CrossAxisAlignment.start, mainAxisSize: MainAxisSize.min, children: [
        Row(mainAxisSize: MainAxisSize.min, children: [
          Text((label ?? context.l10n.exampleLabel).toUpperCase(),
              style: TextStyle(fontSize: 9.5, letterSpacing: 0.6, fontWeight: FontWeight.w700, color: t.textFaint)),
          if (trailing != null) ...[const SizedBox(width: 8), trailing!],
        ]),
        const SizedBox(height: 6),
        child,
        if (caption != null) ...[const SizedBox(height: 6), _Caption(caption!)],
      ]),
    );
  }
}

const _cellH = 15.0;

/// One letter of a toy alignment, coloured like the alignment viewer
/// (nucleotides on [ntColors], amino acids on [aaColors]) unless [color]
/// is given.
Widget _cell(BuildContext context, String ch, {Color? color, double w = 12, bool aa = false, bool faded = false}) {
  final t = context.tokens;
  final gap = ch == '-';
  final c = gap ? t.surface3 : (color ?? (aa ? aaColors[ch] : ntColors[ch]) ?? const Color(0xFFCBD5E1));
  final ink = gap ? t.textFaint : residueInk(c);
  return Opacity(
    opacity: faded ? 0.28 : 1,
    child: Container(
      width: w,
      height: _cellH,
      margin: const EdgeInsets.only(right: 1),
      alignment: Alignment.center,
      decoration: BoxDecoration(color: c, borderRadius: BorderRadius.circular(2)),
      child: Text(ch, style: TextStyle(fontFamily: monoFont, fontSize: 9.5, height: 1, fontWeight: FontWeight.w700, color: ink)),
    ),
  );
}

/// Species names of a toy alignment: one column of the same width, so
/// every sequence starts at the same place.
class _Names extends StatelessWidget {
  const _Names(this.count, {this.top = 0});
  final int count;

  /// Space above the first name (a ruler or markers row).
  final double top;

  static const width = 64.0;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    return SizedBox(
      width: width,
      child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
        SizedBox(height: top),
        for (var i = 0; i < count; i++)
          Container(
            height: _cellH,
            margin: const EdgeInsets.only(bottom: 2),
            alignment: Alignment.centerLeft,
            child: Text(
              context.l10n.exampleSpecies(String.fromCharCode(65 + i)),
              maxLines: 1,
              overflow: TextOverflow.ellipsis,
              style: TextStyle(fontSize: 10.5, height: 1, fontStyle: FontStyle.italic, color: t.textMuted),
            ),
          ),
      ]),
    );
  }
}

/// The same three genes aligned in each mode of step 5b, as one alignment
/// (same columns, same colours as the viewer): nucleotides (the gap of
/// species B cuts codons 3 and 4), amino acids (one letter per codon) or
/// codons (the gap is the whole codon 3).
class SeqModeExample extends StatelessWidget {
  const SeqModeExample({super.key, required this.mode});
  final SeqMode mode;

  static const _w = 12.0, _codonGap = 4.0;

  @override
  Widget build(BuildContext context) {
    final l10n = context.l10n;
    final t = context.tokens;
    final rows = switch (mode) {
      SeqMode.nucleotide => const ['ATG GCT TAT TCA', 'ATG GCT T-- -CA', 'ATG ACT AAT TCA'],
      SeqMode.aminoAcid => const ['M A Y S', 'M A - S', 'M T N S'],
      SeqMode.codon => const ['ATG GCT TAT TCA', 'ATG GCT --- TCA', 'ATG ACT AAT TCA'],
    };
    final aa = mode == SeqMode.aminoAcid;
    // One amino acid takes the width of its codon, so the three pictures
    // share the codon ruler.
    const codonW = 3 * (_w + 1) - 1;
    Widget group(String g) => Row(mainAxisSize: MainAxisSize.min, children: [
          for (final ch in g.split('')) _cell(context, ch, w: aa ? codonW : _w, aa: aa),
        ]);
    final ruler = Row(mainAxisSize: MainAxisSize.min, children: [
      for (var i = 0; i < 4; i++) ...[
        if (i > 0) const SizedBox(width: _codonGap),
        Container(
          width: codonW + 1,
          height: 12,
          alignment: Alignment.center,
          margin: const EdgeInsets.only(bottom: 3),
          decoration: BoxDecoration(border: Border(bottom: BorderSide(color: t.textFaint, width: 1))),
          child: Text('${i + 1}', style: TextStyle(fontSize: 9, height: 1, fontWeight: FontWeight.w700, color: t.textFaint)),
        ),
      ],
    ]);
    return ExampleFrame(
      caption: switch (mode) {
        SeqMode.nucleotide => l10n.exampleNt,
        SeqMode.aminoAcid => l10n.exampleAa,
        SeqMode.codon => l10n.exampleCodon,
      },
      child: Row(mainAxisSize: MainAxisSize.min, crossAxisAlignment: CrossAxisAlignment.start, children: [
        const _Names(3, top: 15),
        Column(crossAxisAlignment: CrossAxisAlignment.start, mainAxisSize: MainAxisSize.min, children: [
          ruler,
          for (final r in rows)
            Padding(
              padding: const EdgeInsets.only(bottom: 2),
              child: Row(mainAxisSize: MainAxisSize.min, children: [
                for (final (i, g) in r.split(' ').indexed) ...[
                  if (i > 0) const SizedBox(width: _codonGap),
                  group(g),
                ],
              ]),
            ),
        ]),
      ]),
    );
  }
}

/// A toy MAFFT alignment of 10 species × 24 columns, animated: the
/// alignment, the columns the chosen method removes, then the trimmed
/// alignment. Like a real alignment it has no column made only of gaps:
/// the gappy columns are insertions of one or a few species, a ragged
/// region and short deletions. It fills the width it is given, plays on
/// its own, rests on the result and starts again.
class TrimExample extends StatefulWidget {
  const TrimExample({super.key, required this.method, this.gapThreshold});

  /// trimAl option (`-gappyout`, `-nogaps`, `manual`…) or ClipKIT mode
  /// (`gappy`, `kpic`, `kpic-gappy`, `smart-gap`).
  final String method;

  /// trimAl `-gt` in the manual mode (fraction of sequences without a gap).
  final double? gapThreshold;

  static const rows = [
    'ATGT--CTAGGCACATAGCATCGA',
    'ATGT--CTAGGCACATAGCATCGA',
    'ATGC--CTAGGTACCTAGCA--GA',
    'ATGC--CTAGG-ACGTAGCA--GA',
    'ATGC--CTGGG-ACTTAGCA--GA',
    'ATGC--CTAGA-ACTTA-CA--GA',
    'ATGC--CTAGA-ACCTA-CA--GA',
    'ATGC--CTAGA-ACATAGCA--GA',
    'ATGC--CTAGA-ACGTAGCT--GA',
    'ATGCTACTAGA-ACTTAGCA--GA',
  ];

  /// Fraction of the sequences with a gap in column [c].
  static double gapFraction(int c) => rows.where((r) => r[c] == '-').length / rows.length;

  /// Columns removed by each method in this example.
  static Set<int> removedColumns(String method, double? gapThreshold) {
    final cols = rows.first.length;
    Set<int> where(bool Function(double gaps) f) => {for (var c = 0; c < cols; c++) if (f(gapFraction(c))) c};
    // Insertion of species J (4–5) and of species A–B (20–21).
    final gappy = where((g) => g >= 0.75);
    const lowSimilarity = 14; // four different bases
    const singletons = {8, 19}; // one species differs: no phylogenetic signal
    return switch (method) {
      '-nogaps' => where((g) => g > 0),
      // A MAFFT alignment has no column made only of gaps.
      '-noallgaps' => where((g) => g == 1),
      '-strict' => {...gappy, lowSimilarity},
      '-strictplus' => {...gappy, lowSimilarity, 11},
      'manual' => where((g) => g > 1 - (gapThreshold ?? 0.9) + 1e-9),
      'gappy' => where((g) => g >= 0.9),
      'kpic' => singletons,
      'kpic-gappy' => {...singletons, ...where((g) => g >= 0.9)},
      _ => gappy, // -automated1, -gappyout, smart-gap
    };
  }

  @override
  State<TrimExample> createState() => _TrimExampleState();
}

class _TrimExampleState extends State<TrimExample> with SingleTickerProviderStateMixin {
  late final _anim = AnimationController(vsync: this, duration: const Duration(milliseconds: 5200));
  Timer? _rest;

  static const _markH = 14.0, _barH = 16.0;
  static const _red = SplacePalette.red700;

  // Phases of one cycle (fractions of the controller).
  static const _mark = Interval(0.16, 0.30, curve: Curves.easeOut);
  static const _fade = Interval(0.42, 0.58, curve: Curves.easeIn);
  static const _close = Interval(0.62, 0.84, curve: Curves.easeInOut);

  @override
  void initState() {
    super.initState();
    _play();
  }

  @override
  void didUpdateWidget(TrimExample old) {
    super.didUpdateWidget(old);
    if (old.method != widget.method || old.gapThreshold != widget.gapThreshold) _play();
  }

  void _play() {
    _rest?.cancel();
    _anim.forward(from: 0).then((_) {
      if (!mounted) return;
      // Rest on the trimmed alignment, then show it again.
      _rest = Timer(const Duration(milliseconds: 2600), () {
        if (mounted) _play();
      });
    });
  }

  @override
  void dispose() {
    _rest?.cancel();
    _anim.dispose();
    super.dispose();
  }

  @override
  Widget build(BuildContext context) {
    final l10n = context.l10n;
    final t = context.tokens;
    final removed = TrimExample.removedColumns(widget.method, widget.gapThreshold);
    final rows = TrimExample.rows;
    final cols = rows.first.length;
    return ExampleFrame(
      caption: removed.isEmpty ? l10n.exampleTrimNone('$cols') : l10n.exampleTrim('${removed.length}', '$cols'),
      trailing: SizedBox(
        width: 18,
        height: 14,
        child: IconButton(
          tooltip: l10n.examplePlayAgain,
          padding: EdgeInsets.zero,
          iconSize: 10,
          color: t.textMuted,
          onPressed: _play,
          icon: const Icon(FA.rotateLeft),
        ),
      ),
      child: LayoutBuilder(builder: (context, constraints) {
        // The columns share the width left by the names.
        final room = constraints.maxWidth.isFinite ? constraints.maxWidth : 520.0;
        final w = ((room - _Names.width) / cols - 1).clamp(10.0, 34.0);
        return AnimatedBuilder(
          animation: _anim,
          builder: (context, _) {
            final v = _anim.value;
            final mark = _mark.transform(v), fade = _fade.transform(v), close = _close.transform(v);
            final phase = v < _mark.begin ? 0 : (v < _close.begin ? 1 : 2);
            Widget step(int i, String text) {
              final on = phase == i;
              return AnimatedContainer(
                duration: const Duration(milliseconds: 200),
                padding: const EdgeInsets.symmetric(horizontal: 7, vertical: 3),
                decoration: BoxDecoration(
                  color: on ? (i == 1 ? _red.withValues(alpha: 0.14) : t.accentSoft) : Colors.transparent,
                  borderRadius: BorderRadius.circular(8),
                  border: Border.all(color: on ? (i == 1 ? _red : t.accent) : t.border),
                ),
                child: Text('${i + 1}. $text',
                    style: TextStyle(fontSize: 10.5, fontWeight: on ? FontWeight.w700 : FontWeight.w500, color: on ? (i == 1 ? _red : t.accentText) : t.textFaint)),
              );
            }

            // One column: the × mark, the bases and how many gaps it has.
            Widget column(int c, {bool marked = false}) {
              final g = TrimExample.gapFraction(c);
              return Column(mainAxisSize: MainAxisSize.min, children: [
                SizedBox(
                  width: w + 1,
                  height: _markH + 2,
                  child: marked
                      ? Opacity(
                          opacity: mark,
                          child: const Center(child: Text('×', style: TextStyle(fontSize: 16, height: 1, fontWeight: FontWeight.w800, color: _red))),
                        )
                      : null,
                ),
                for (final r in rows) Padding(padding: const EdgeInsets.only(bottom: 2), child: _cell(context, r[c], w: w)),
                SizedBox(
                  width: w + 1,
                  height: _barH,
                  child: Align(
                    alignment: Alignment.bottomCenter,
                    child: Container(
                      width: w - 2,
                      height: math.max(1, (_barH - 2) * g),
                      decoration: BoxDecoration(
                        color: g == 0 ? t.border : (marked ? _red : t.textFaint),
                        borderRadius: BorderRadius.circular(1),
                      ),
                    ),
                  ),
                ),
              ]);
            }

            return Column(crossAxisAlignment: CrossAxisAlignment.start, mainAxisSize: MainAxisSize.min, children: [
              Wrap(spacing: 4, runSpacing: 4, children: [
                step(0, l10n.trimStepOriginal),
                step(1, l10n.trimStepMarked),
                step(2, l10n.trimStepResult),
              ]),
              const SizedBox(height: 8),
              Row(mainAxisSize: MainAxisSize.min, crossAxisAlignment: CrossAxisAlignment.start, children: [
                Column(crossAxisAlignment: CrossAxisAlignment.start, mainAxisSize: MainAxisSize.min, children: [
                  _Names(rows.length, top: _markH + 2),
                  SizedBox(
                    width: _Names.width,
                    height: _barH,
                    child: Align(
                      alignment: Alignment.bottomLeft,
                      child: Text(l10n.exampleGapsRow, style: TextStyle(fontSize: 10, height: 1, color: t.textFaint)),
                    ),
                  ),
                ]),
                for (var c = 0; c < cols; c++)
                  if (!removed.contains(c))
                    column(c)
                  else
                    ClipRect(
                      child: Align(
                        alignment: Alignment.centerLeft,
                        widthFactor: 1 - close,
                        child: Opacity(
                          opacity: 1 - fade,
                          child: DecoratedBox(
                            decoration: BoxDecoration(color: _red.withValues(alpha: 0.22 * mark), borderRadius: BorderRadius.circular(3)),
                            child: column(c, marked: true),
                          ),
                        ),
                      ),
                    ),
              ]),
            ]);
          },
        );
      }),
    );
  }
}

/// Trimming a coding gene, both cases in one animation over the toy
/// alignment of [TrimExample]: "column by column" (trimAl [method]) removes
/// the columns it marks (×), so the codon positions (the coloured ruler)
/// fall out of order; "keep the reading frame" (ClipKIT [codonMethod]
/// `--codon`) removes the whole codon of every column it marks (+), so
/// 1-2-3 stays in order. The lane of the option chosen ([keepFrame]) is
/// outlined. Plays on its own.
class TrimFrameExample extends StatefulWidget {
  const TrimFrameExample(
      {super.key, required this.method, this.gapThreshold, this.codonMethod = 'smart-gap', this.codonTool = 'trimAl -backtrans', required this.keepFrame});
  final String method;
  final double? gapThreshold;

  /// Method of the whole-codon trimming (trimAl option or ClipKIT mode)
  /// and the tool's name ("trimAl -backtrans", "ClipKIT --codon").
  final String codonMethod, codonTool;

  /// "Keep the reading frame" is on.
  final bool keepFrame;

  /// The columns of every codon that has a removed column.
  static Set<int> wholeCodons(Set<int> removed) => {for (final c in removed) ...[c - c % 3, c - c % 3 + 1, c - c % 3 + 2]};

  /// Kept columns whose codon position no longer matches their place in
  /// the trimmed alignment.
  static Set<int> outOfPhase(Set<int> removed, int cols) {
    final out = <int>{};
    var i = 0;
    for (var c = 0; c < cols; c++) {
      if (removed.contains(c)) continue;
      if (c % 3 != i % 3) out.add(c);
      i++;
    }
    return out;
  }

  @override
  State<TrimFrameExample> createState() => _TrimFrameExampleState();
}

class _TrimFrameExampleState extends State<TrimFrameExample> with SingleTickerProviderStateMixin {
  late final _anim = AnimationController(vsync: this, duration: const Duration(milliseconds: 7600));
  Timer? _rest;

  static const _red = SplacePalette.red700;
  static const _codonGap = 4.0, _rulerH = 14.0, _markH = 15.0, _barH = 12.0;

  // Phases of one cycle (fractions of the controller).
  static const _mark = Interval(0.07, 0.18, curve: Curves.easeOut);
  static const _grow = Interval(0.24, 0.36, curve: Curves.easeOut);
  static const _fade = Interval(0.44, 0.54, curve: Curves.easeIn);
  static const _close = Interval(0.56, 0.76, curve: Curves.easeInOut);
  static const _verdict = Interval(0.78, 0.90, curve: Curves.easeOut);

  @override
  void initState() {
    super.initState();
    _play();
  }

  @override
  void didUpdateWidget(TrimFrameExample old) {
    super.didUpdateWidget(old);
    if (old.method != widget.method || old.gapThreshold != widget.gapThreshold || old.codonMethod != widget.codonMethod) _play();
  }

  void _play() {
    _rest?.cancel();
    _anim.forward(from: 0).then((_) {
      if (!mounted) return;
      // Rest on the result, then show it again.
      _rest = Timer(const Duration(milliseconds: 3000), () {
        if (mounted) _play();
      });
    });
  }

  @override
  void dispose() {
    _rest?.cancel();
    _anim.dispose();
    super.dispose();
  }

  @override
  Widget build(BuildContext context) {
    final l10n = context.l10n;
    final t = context.tokens;
    const rows = TrimExample.rows;
    final cols = rows.first.length;
    // Columns each tool marks: trimAl (column by column), ClipKIT (codons).
    final trimalMarks = TrimExample.removedColumns(widget.method, widget.gapThreshold);
    final clipMarks = TrimExample.removedColumns(widget.codonMethod, widget.gapThreshold);
    final codons = TrimFrameExample.wholeCodons(clipMarks);
    final marked = {...trimalMarks, ...clipMarks};
    return ExampleFrame(
      label: l10n.trimCasesLabel,
      caption: marked.isEmpty ? l10n.exampleTrimNone('$cols') : l10n.trimCasesCaptionTool(widget.codonTool),
      trailing: SizedBox(
        width: 18,
        height: 14,
        child: IconButton(
          tooltip: l10n.examplePlayAgain,
          padding: EdgeInsets.zero,
          iconSize: 10,
          color: t.textMuted,
          onPressed: _play,
          icon: const Icon(FA.rotateLeft),
        ),
      ),
      child: LayoutBuilder(builder: (context, constraints) {
        final room = constraints.maxWidth.isFinite ? constraints.maxWidth : 720.0;
        // Side by side when both fit, one above the other otherwise.
        final wide = room >= 820;
        final lane = wide ? (room - 12) / 2 : room;
        final w = ((lane - 20 - _Names.width - (cols ~/ 3 - 1) * _codonGap) / cols - 1).clamp(8.0, 26.0);
        return AnimatedBuilder(
          animation: _anim,
          builder: (context, _) {
            final v = _anim.value;
            final mark = _mark.transform(v), grow = _grow.transform(v), fade = _fade.transform(v);
            final close = _close.transform(v), verdict = _verdict.transform(v);
            final phase = v < _mark.begin ? 0 : (v < _grow.begin ? 1 : (v < _close.begin ? 2 : 3));

            Widget step(int i, String text, {bool red = false}) {
              final on = phase == i;
              return AnimatedContainer(
                duration: const Duration(milliseconds: 200),
                padding: const EdgeInsets.symmetric(horizontal: 7, vertical: 3),
                decoration: BoxDecoration(
                  color: on ? (red ? _red.withValues(alpha: 0.14) : t.accentSoft) : Colors.transparent,
                  borderRadius: BorderRadius.circular(8),
                  border: Border.all(color: on ? (red ? _red : t.accent) : t.border),
                ),
                child: Text('${i + 1}. $text',
                    style: TextStyle(
                        fontSize: 10.5, fontWeight: on ? FontWeight.w700 : FontWeight.w500, color: on ? (red ? _red : t.accentText) : t.textFaint)),
              );
            }

            Widget laneView({required bool frame}) {
              final marked = frame ? clipMarks : trimalMarks;
              final removed = frame ? codons : marked;
              final shifted = TrimFrameExample.outOfPhase(removed, cols);
              final selected = frame == widget.keepFrame;

              // One column: codon position, the × / + mark, bases, gaps.
              Widget column(int c) {
                final first = marked.contains(c);
                final extra = frame && !first && codons.contains(c);
                final red = first ? mark : (extra ? grow : 0.0);
                final pos = c % 3;
                final wrong = shifted.contains(c) ? verdict : 0.0;
                final g = TrimExample.gapFraction(c);
                return DecoratedBox(
                  decoration: BoxDecoration(color: _red.withValues(alpha: 0.2 * red), borderRadius: BorderRadius.circular(3)),
                  child: Column(mainAxisSize: MainAxisSize.min, children: [
                    Container(
                      width: w,
                      height: _rulerH - 2,
                      margin: const EdgeInsets.only(right: 1, bottom: 1),
                      alignment: Alignment.center,
                      decoration: BoxDecoration(
                        color: codonPositionColors[pos],
                        borderRadius: BorderRadius.circular(2),
                        border: wrong > 0 ? Border.all(color: Color.lerp(codonPositionColors[pos], _red, wrong)!, width: 2) : null,
                      ),
                      child: Text('${pos + 1}', style: TextStyle(fontSize: 8.5, height: 1, fontWeight: FontWeight.w800, color: residueInk(codonPositionColors[pos]))),
                    ),
                    SizedBox(
                      width: w + 1,
                      height: _markH,
                      child: first || extra
                          ? Opacity(
                              opacity: red,
                              child: Center(
                                child: Text(first ? '×' : '+',
                                    style: const TextStyle(fontSize: 14, height: 1, fontWeight: FontWeight.w800, color: _red)),
                              ),
                            )
                          : wrong > 0
                              ? Opacity(
                                  opacity: wrong,
                                  child: const Center(child: Icon(FA.triangleExclamation, size: 8, color: _red)),
                                )
                              : null,
                    ),
                    for (final r in rows) Padding(padding: const EdgeInsets.only(bottom: 2), child: _cell(context, r[c], w: w)),
                    SizedBox(
                      width: w + 1,
                      height: _barH,
                      child: Align(
                        alignment: Alignment.bottomCenter,
                        child: Container(
                          width: math.max(2, w - 2),
                          height: math.max(1, (_barH - 2) * g),
                          decoration: BoxDecoration(
                            color: g == 0 ? t.border : (first ? _red : t.textFaint),
                            borderRadius: BorderRadius.circular(1),
                          ),
                        ),
                      ),
                    ),
                  ]),
                );
              }

              Widget shrink(Widget child) => ClipRect(
                    child: Align(
                      alignment: Alignment.centerLeft,
                      widthFactor: 1 - close,
                      child: Opacity(opacity: 1 - fade, child: child),
                    ),
                  );

              final children = <Widget>[];
              for (var k = 0; k < cols ~/ 3; k++) {
                final group = [for (var c = 3 * k; c < 3 * k + 3; c++) removed.contains(c) ? shrink(column(c)) : column(c)];
                children.add(Row(mainAxisSize: MainAxisSize.min, crossAxisAlignment: CrossAxisAlignment.start, children: group));
                final gone = removed.containsAll([3 * k, 3 * k + 1, 3 * k + 2]);
                if (k < cols ~/ 3 - 1) children.add(SizedBox(width: gone ? _codonGap * (1 - close) : _codonGap));
              }

              final keptCodons = cols ~/ 3 - codons.length ~/ 3;
              final ok = shifted.isEmpty;
              final verdictText = frame
                  ? l10n.trimCaseFrameKept('$keptCodons', '${cols ~/ 3}')
                  : ok
                      ? l10n.trimCaseFrameKeptCols
                      : l10n.trimCaseFrameLost('${shifted.length}');
              return Container(
                padding: const EdgeInsets.fromLTRB(9, 8, 9, 8),
                decoration: BoxDecoration(
                  color: t.surface,
                  borderRadius: BorderRadius.circular(8),
                  border: Border.all(color: selected ? t.accent : t.border, width: selected ? 1.5 : 1),
                ),
                child: Column(crossAxisAlignment: CrossAxisAlignment.start, mainAxisSize: MainAxisSize.min, children: [
                  Row(children: [
                    Icon(frame ? FA.dna : FA.scissors, size: 11, color: selected ? t.accent : t.textMuted),
                    const SizedBox(width: 6),
                    Flexible(
                      child: Text(frame ? l10n.trimCaseCodonsTool(widget.codonTool) : l10n.trimCaseColumns,
                          maxLines: 1,
                          overflow: TextOverflow.ellipsis,
                          style: TextStyle(fontSize: 11.5, fontWeight: FontWeight.w700, color: selected ? t.accentText : t.text)),
                    ),
                    if (selected) ...[
                      const SizedBox(width: 6),
                      Container(
                        padding: const EdgeInsets.symmetric(horizontal: 6, vertical: 1),
                        decoration: BoxDecoration(color: t.accentSoft, borderRadius: BorderRadius.circular(8)),
                        child: Text(l10n.trimCaseSelected, style: TextStyle(fontSize: 9.5, fontWeight: FontWeight.w700, color: t.accentText)),
                      ),
                    ],
                  ]),
                  const SizedBox(height: 6),
                  FittedBox(
                    fit: BoxFit.scaleDown,
                    alignment: Alignment.centerLeft,
                    child: Row(mainAxisSize: MainAxisSize.min, crossAxisAlignment: CrossAxisAlignment.start, children: [
                      Column(crossAxisAlignment: CrossAxisAlignment.start, mainAxisSize: MainAxisSize.min, children: [
                        _Names(rows.length, top: _rulerH + _markH - 1),
                        SizedBox(
                          width: _Names.width,
                          height: _barH,
                          child: Align(
                            alignment: Alignment.bottomLeft,
                            child: Text(l10n.exampleGapsRow, style: TextStyle(fontSize: 9.5, height: 1, color: t.textFaint)),
                          ),
                        ),
                      ]),
                      ...children,
                    ]),
                  ),
                  const SizedBox(height: 6),
                  Opacity(
                    opacity: marked.isEmpty ? 1 : verdict,
                    child: Row(children: [
                      Icon(ok ? FA.circleCheck : FA.triangleExclamation, size: 11, color: ok ? SplaceTokens.success : _red),
                      const SizedBox(width: 6),
                      Flexible(
                        child: Text(verdictText,
                            style: TextStyle(fontSize: 11, fontWeight: FontWeight.w600, color: ok ? SplaceTokens.success : _red)),
                      ),
                    ]),
                  ),
                ]),
              );
            }

            final a = laneView(frame: false), b = laneView(frame: true);
            return Column(crossAxisAlignment: CrossAxisAlignment.start, mainAxisSize: MainAxisSize.min, children: [
              Wrap(spacing: 4, runSpacing: 4, children: [
                step(0, l10n.trimStepOriginal),
                step(1, l10n.trimCasesStepMark, red: true),
                step(2, l10n.frameStepCodonTool(widget.codonTool), red: true),
                step(3, l10n.trimStepResult),
              ]),
              const SizedBox(height: 8),
              if (wide)
                IntrinsicHeight(
                  child: Row(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
                    Expanded(child: a),
                    const SizedBox(width: 12),
                    Expanded(child: b),
                  ]),
                )
              else ...[
                a,
                const SizedBox(height: 10),
                b,
              ],
            ]);
          },
        );
      }),
    );
  }
}

/// Codons of a coding gene with the excluded positions struck out, and
/// what is left for the tree.
class CodonExclusionExample extends StatelessWidget {
  const CodonExclusionExample({super.key, required this.excluded});
  final Set<int> excluded;

  static const _codons = ['ATG', 'GCT', 'TAT', 'TCA', 'AAG', 'CTG', 'GGA', 'ACC', 'TTC', 'GAT', 'CCA', 'AGC'];
  static const _codonGap = 6.0;

  @override
  Widget build(BuildContext context) {
    final l10n = context.l10n;
    final t = context.tokens;
    final kept = [for (final p in [1, 2, 3]) if (!excluded.contains(p)) p];
    Widget label(String text) => SizedBox(
          width: _Names.width + 16,
          child: Text(text, maxLines: 1, overflow: TextOverflow.ellipsis, style: TextStyle(fontSize: 10.5, height: 1, color: t.textMuted)),
        );
    return ExampleFrame(
      caption: excluded.isEmpty
          ? l10n.exampleCodonKeepAll
          : l10n.exampleCodonKept(kept.map((p) => '$pª').join(' + '), '${(100 * kept.length / 3).round()}'),
      // Fills the width: the cells grow with the window.
      child: LayoutBuilder(builder: (context, box) {
        final room = (box.maxWidth.isFinite ? box.maxWidth : 600.0) - _Names.width - 16;
        final n = _codons.length;
        final w = ((room - _codonGap * (n - 1)) / (3 * n) - 1).clamp(10.0, 30.0);
        Widget codonRow(bool after) => Row(mainAxisSize: MainAxisSize.min, children: [
              for (final (i, codon) in _codons.indexed) ...[
                if (i > 0) const SizedBox(width: _codonGap),
                for (final (p, ch) in codon.split('').indexed)
                  if (!after)
                    Stack(alignment: Alignment.center, children: [
                      _cell(context, ch, color: codonPositionColors[p], w: w, faded: excluded.contains(p + 1)),
                      if (excluded.contains(p + 1)) Container(width: w, height: 1.6, color: SplacePalette.red700),
                    ])
                  else if (!excluded.contains(p + 1))
                    _cell(context, ch, color: codonPositionColors[p], w: w),
              ],
            ]);
        return Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
          Row(children: [
            SizedBox(width: _Names.width + 16),
            for (final (i, _) in _codons.indexed) ...[
              if (i > 0) const SizedBox(width: _codonGap),
              for (var p = 1; p <= 3; p++)
                SizedBox(
                  width: w + 1,
                  child: Text('$p', textAlign: TextAlign.center, style: TextStyle(fontSize: 9, fontWeight: FontWeight.w700, color: codonPositionColors[p - 1])),
                ),
            ],
          ]),
          const SizedBox(height: 3),
          Row(children: [label(l10n.exampleCodonGene), codonRow(false)]),
          const SizedBox(height: 6),
          Row(children: [label(l10n.exampleCodonMatrix), codonRow(true)]),
        ]);
      }),
    );
  }
}

/// How saturation looks in the "p × corrected distance" plot: points on a
/// rising line (no saturation) against points that level off (saturated).
class SaturationReadingExample extends StatelessWidget {
  const SaturationReadingExample({super.key});

  @override
  Widget build(BuildContext context) {
    final l10n = context.l10n;
    final t = context.tokens;
    Widget plot(bool saturated, String title, String caption, Color color) => Expanded(
          child: Row(crossAxisAlignment: CrossAxisAlignment.start, children: [
            CustomPaint(
              size: const Size(150, 104),
              painter: _CurvePainter(saturated: saturated, color: color, axis: t.textMuted, grid: t.border),
            ),
            const SizedBox(width: 12),
            Expanded(
              child: Column(crossAxisAlignment: CrossAxisAlignment.start, mainAxisSize: MainAxisSize.min, children: [
                Row(children: [
                  Container(width: 8, height: 8, decoration: BoxDecoration(color: color, shape: BoxShape.circle)),
                  const SizedBox(width: 6),
                  Expanded(child: Text(title, style: TextStyle(fontSize: 12.5, fontWeight: FontWeight.w700, color: t.text))),
                ]),
                const SizedBox(height: 4),
                _Caption(caption),
              ]),
            ),
          ]),
        );
    return ExampleFrame(
      label: l10n.exampleSatLabel,
      child: Row(crossAxisAlignment: CrossAxisAlignment.start, children: [
        plot(false, l10n.exampleSatOk, l10n.exampleSatOkBody, SplaceBlue.c500),
        const SizedBox(width: 24),
        plot(true, l10n.exampleSatBad, l10n.exampleSatBadBody, SplaceRed.c500),
      ]),
    );
  }
}

class _CurvePainter extends CustomPainter {
  _CurvePainter({required this.saturated, required this.color, required this.axis, required this.grid});
  final bool saturated;
  final Color color, axis, grid;

  @override
  void paint(Canvas canvas, Size size) {
    const l = 18.0, b = 14.0;
    final plot = Rect.fromLTRB(l, 4, size.width - 4, size.height - b);
    final axisPaint = Paint()
      ..color = axis
      ..strokeWidth = 1;
    canvas.drawLine(plot.bottomLeft, plot.topLeft, axisPaint);
    canvas.drawLine(plot.bottomLeft, plot.bottomRight, axisPaint);
    // y = x, dashed.
    final dash = Paint()
      ..color = axis.withValues(alpha: 0.6)
      ..strokeWidth = 1;
    for (var f = 0.0; f < 1; f += 0.06) {
      canvas.drawLine(Offset(plot.left + plot.width * f, plot.bottom - plot.height * f),
          Offset(plot.left + plot.width * math.min(1, f + 0.03), plot.bottom - plot.height * math.min(1, f + 0.03)), dash);
    }
    final rnd = math.Random(saturated ? 3 : 7);
    final dot = Paint()..color = color.withValues(alpha: 0.75);
    for (var i = 0; i < 26; i++) {
      final x = 0.04 + 0.92 * i / 25;
      final y = (saturated ? 0.62 * (1 - math.exp(-4.2 * x)) : 0.82 * x) + (rnd.nextDouble() - 0.5) * 0.05;
      canvas.drawCircle(Offset(plot.left + plot.width * x, plot.bottom - plot.height * y.clamp(0.0, 1.0)), 2.2, dot);
    }
    void label(String s, Offset at, {bool vertical = false}) {
      final tp = TextPainter(text: TextSpan(text: s, style: TextStyle(fontFamily: 'Montserrat', fontSize: 9.5, color: axis, fontWeight: FontWeight.w600)), textDirection: TextDirection.ltr)
        ..layout();
      canvas.save();
      canvas.translate(at.dx, at.dy);
      if (vertical) canvas.rotate(-math.pi / 2);
      tp.paint(canvas, Offset(-tp.width / 2, -tp.height / 2));
      canvas.restore();
    }

    label('p', Offset(7, plot.center.dy), vertical: true);
    label('d', Offset(plot.center.dx, size.height - 6));
  }

  @override
  bool shouldRepaint(_CurvePainter old) => old.color != color || old.axis != axis;
}
