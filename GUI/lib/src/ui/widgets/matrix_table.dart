import 'dart:math' as math;

import 'package:flutter/gestures.dart';
import 'package:flutter/material.dart';
import 'package:flutter/services.dart';

import '../../core/supermatrix.dart';
import '../common.dart';
import '../theme.dart';
import 'alignment_viewer.dart' show ScrollTrack;
import 'kit.dart' show Seg;
import 'splace_tooltip.dart';

/// The supermatrix as SequenceMatrix shows it: one row per taxon, one
/// column per locus. Each cell holds the bases the taxon has in that locus
/// (darker = more complete) or "—" when the locus is missing; the header
/// gives each alignment's length, how many taxa have it and its gaps.
///
/// The table has a fixed height and scrolls inside itself (wheel,
/// shift+wheel, the tracks): only the visible rows and columns are drawn,
/// so thousands of taxa and hundreds of loci stay smooth and the page
/// does not grow with the number of taxa.
class MatrixOverviewTable extends StatefulWidget {
  const MatrixOverviewTable({super.key, required this.overview});
  final MatrixOverview overview;

  @override
  State<MatrixOverviewTable> createState() => _MatrixOverviewTableState();
}

class _MatrixOverviewTableState extends State<MatrixOverviewTable> {
  final _x = ValueNotifier<double>(0);
  final _y = ValueNotifier<double>(0);

  /// (row, column): row -1 is the header; column -1 the fixed columns.
  final _hover = ValueNotifier<(int, int)?>(null);

  static const rowH = 22.0, headerH = 50.0;
  static const statW = 74.0;
  static const maxBodyH = 420.0;

  /// Width of the taxon names (drag the edge; double-click fits the
  /// longest name). Shared by every table.
  static double _nameW = 190;
  double get _fixedW => _nameW + 3 * statW;

  double _viewW = 0, _maxX = 0, _maxY = 0;

  /// Column width now: the loci fill the table when they fit.
  double _colWNow = 64;

  /// Cells show bases (false) or % of the marker's length; rows in the
  /// matrix order or the most complete taxa first.
  bool _percent = false, _byCompleteness = false;
  List<int> _order = const [];

  List<int> _orderOf(MatrixOverview o) {
    final idx = List<int>.generate(o.taxa.length, (i) => i);
    if (_byCompleteness) idx.sort((a, b) => o.basesOf(b).compareTo(o.basesOf(a)));
    return idx;
  }

  void _fitNames() {
    var w = 0.0;
    for (final n in widget.overview.taxa.take(3000)) {
      final tp = TextPainter(text: TextSpan(text: n, style: const TextStyle(fontSize: 12, fontFamily: _sans)), textDirection: TextDirection.ltr)..layout();
      w = math.max(w, tp.width);
    }
    setState(() => _nameW = (w + 24).clamp(90.0, 900.0));
  }

  @override
  void didUpdateWidget(MatrixOverviewTable old) {
    super.didUpdateWidget(old);
    if (old.overview != widget.overview) {
      _x.value = 0;
      _y.value = 0;
      _hover.value = null;
    }
  }

  @override
  void dispose() {
    _x.dispose();
    _y.dispose();
    _hover.dispose();
    super.dispose();
  }

  void _setX(double v) => _x.value = v.clamp(0.0, _maxX);
  void _setY(double v) => _y.value = v.clamp(0.0, _maxY);

  void _onSignal(PointerSignalEvent e) {
    if (e is! PointerScrollEvent) return;
    var dx = e.scrollDelta.dx, dy = e.scrollDelta.dy;
    if (dx == 0 && HardwareKeyboard.instance.isShiftPressed) {
      dx = dy;
      dy = 0;
    }
    final canH = dx != 0 && (dx < 0 ? _x.value > 0 : _x.value < _maxX);
    final canV = dy != 0 && (dy < 0 ? _y.value > 0 : _y.value < _maxY);
    // At the ends the wheel scrolls the page instead.
    if (!canH && !canV) return;
    GestureBinding.instance.pointerSignalResolver.register(e, (_) {
      if (canH) _setX(_x.value + dx);
      if (canV) _setY(_y.value + dy);
    });
  }

  /// Row under a point of the body (null in the header or past the end).
  int? _row(double dy) {
    if (dy < headerH) return -1;
    final r = ((dy - headerH + _y.value) / rowH).floor();
    return r >= 0 && r < widget.overview.taxa.length ? r : null;
  }

  @override
  Widget build(BuildContext context) {
    final l10n = context.l10n;
    final t = context.tokens;
    final o = widget.overview;
    final rowsH = rowH * o.taxa.length;
    _order = _orderOf(o);
    final order = _order;
    final bodyH = headerH + math.min(rowsH, maxBodyH);
    _maxY = math.max(0, rowsH - (bodyH - headerH));

    Widget swatch(Color a, Color b, [Color? c]) => Container(
          width: 64,
          height: 10,
          decoration: BoxDecoration(
            borderRadius: BorderRadius.circular(3),
            gradient: LinearGradient(colors: [a, ?c, b]),
          ),
        );
    Widget legend(Widget sw, String lo, String text, String hi) => Row(mainAxisSize: MainAxisSize.min, children: [
          Text(text, style: TextStyle(fontSize: 12, fontWeight: FontWeight.w600, color: t.text)),
          const SizedBox(width: 8),
          Text(lo, style: TextStyle(fontSize: 11, color: t.textMuted)),
          const SizedBox(width: 4),
          sw,
          const SizedBox(width: 4),
          Text(hi, style: TextStyle(fontSize: 11, color: t.textMuted)),
        ]);
    return Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
      // How the cells read, and the colour scales.
      Wrap(spacing: 18, runSpacing: 8, crossAxisAlignment: WrapCrossAlignment.center, children: [
        Seg<bool>(value: _percent, options: [(false, l10n.mtShowBases), (true, l10n.mtShowPercent)], onChanged: (v) => setState(() => _percent = v)),
        Seg<bool>(value: _byCompleteness, options: [(false, l10n.mtOrderMatrix), (true, l10n.mtOrderComplete)], onChanged: (v) => setState(() => _byCompleteness = v)),
        legend(swatch(_fullness(0), _fullness(1)), '0%', l10n.mtScaleComplete, '100%'),
        legend(swatch(_missing(0), _missing(1), _missing(0.5)), '0%', l10n.mtScaleMissing, '100%'),
        Row(mainAxisSize: MainAxisSize.min, children: [
          Container(width: 14, height: 10, decoration: BoxDecoration(color: t.surface3, borderRadius: BorderRadius.circular(3))),
          const SizedBox(width: 6),
          Text(l10n.mtMissingCell, style: TextStyle(fontSize: 11, color: t.textMuted)),
        ]),
      ]),
      const SizedBox(height: 8),
      RepaintBoundary(
        child: Container(
          decoration: BoxDecoration(border: Border.all(color: t.border), borderRadius: BorderRadius.circular(10), color: t.surface),
          clipBehavior: Clip.antiAlias,
          child: LayoutBuilder(builder: (context, box) {
            final vBar = _maxY > 0;
            final fixedW = _fixedW;
            _viewW = math.max(0, box.maxWidth - fixedW - _handleW - (vBar ? ScrollTrack.size : 0));
            // The loci share the whole width when they fit; otherwise they
            // keep their width and the table scrolls sideways.
            final colW = o.loci.isEmpty ? _colW(o) : math.max(_colW(o), _viewW / o.loci.length);
            _colWNow = colW;
            final gridW = colW * o.loci.length;
            _maxX = math.max(0, gridW - _viewW);
            if (_x.value > _maxX) _x.value = _maxX;
            if (_y.value > _maxY) _y.value = _maxY;
            final hBar = _maxX > 0;
            return Listener(
              onPointerSignal: _onSignal,
              child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
                SizedBox(
                  height: bodyH,
                  child: Row(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
                    MouseRegion(
                      onHover: (e) {
                        final r = _row(e.localPosition.dy);
                        _hover.value = r == null || r < 0 ? null : (r, -1);
                      },
                      onExit: (_) => _hover.value = null,
                      child: CustomPaint(
                        size: Size(fixedW, double.infinity),
                        painter: _FixedPainter(o, _y, _hover, t, nameW: _nameW, order: order,
                            heads: [l10n.mtTaxon, l10n.mtLoci, l10n.mtBases, l10n.mtMissing]),
                      ),
                    ),
                    // Drag to widen the names; double-click to fit them.
                    SplaceTooltip(
                      message: l10n.mtResizeNames,
                      child: MouseRegion(
                        cursor: SystemMouseCursors.resizeColumn,
                        child: GestureDetector(
                          behavior: HitTestBehavior.opaque,
                          onHorizontalDragUpdate: (d) => setState(() => _nameW = (_nameW + d.delta.dx).clamp(90.0, 900.0)),
                          onDoubleTap: _fitNames,
                          child: SizedBox(width: _handleW, child: Center(child: Container(width: 2, color: t.border))),
                        ),
                      ),
                    ),
                    Expanded(
                      child: MouseRegion(
                        onHover: (e) {
                          final c = ((e.localPosition.dx + _x.value) / _colWNow).floor();
                          final r = _row(e.localPosition.dy);
                          _hover.value = r == null || c < 0 || c >= o.loci.length ? null : (r, c);
                        },
                        onExit: (_) => _hover.value = null,
                        child: ClipRect(
                          child: CustomPaint(size: Size.infinite, painter: _GridPainter(o, _x, _y, _hover, t, colW, order, _percent)),
                        ),
                      ),
                    ),
                    if (vBar)
                      Padding(
                        padding: const EdgeInsets.only(top: headerH),
                        child: ScrollTrack(axis: Axis.vertical, offset: _y, viewport: bodyH - headerH, content: rowsH, onChanged: _setY),
                      ),
                  ]),
                ),
                if (hBar)
                  Padding(
                    padding: EdgeInsets.only(left: fixedW + _handleW, right: vBar ? ScrollTrack.size : 0),
                    child: ScrollTrack(axis: Axis.horizontal, offset: _x, viewport: _viewW, content: gridW, onChanged: _setX),
                  ),
              ]),
            );
          }),
        ),
      ),
      const SizedBox(height: 6),
      // What the pointer is on, or the legend.
      ValueListenableBuilder<(int, int)?>(
        valueListenable: _hover,
        builder: (context, h, _) {
          String text;
          if (h == null || h.$2 >= o.loci.length) {
            text = l10n.mtLegend;
          } else if (h.$2 < 0) {
            final r = order[h.$1];
            text = '${o.taxa[r]}: ${o.lociOf(r)}/${o.loci.length} · ${o.basesOf(r)} bp';
          } else {
            final l = h.$2, len = o.lengths[l];
            final locus = l10n.mtLocusInfo(o.loci[l], '$len', '${o.taxaWith(l)}', '${o.taxa.length}', '${(100 * o.gapFraction(l)).round()}');
            if (h.$1 < 0) {
              text = locus;
            } else {
              final r = order[h.$1];
              final b = o.bases[r][l], g = o.gaps[r][l];
              text = b == 0
                  ? '${o.taxa[r]} × ${o.loci[l]}: ${l10n.mtMissingCell}'
                  : '${o.taxa[r]} × ${o.loci[l]}: ${l10n.mtCellInfo('$b', '$len', '$g')}';
            }
          }
          return Text(text, maxLines: 2, overflow: TextOverflow.ellipsis, style: TextStyle(fontSize: 12, color: t.textMuted));
        },
      ),
    ]);
  }
}

/// Narrowest column: "766/766 · 12%" needs more room than "20/20 · 3%".
double _colW(MatrixOverview o) => o.taxa.length >= 100 ? 88.0 : 64.0;

/// The draggable edge between the names and the loci.
const _handleW = 7.0;

// A painter does not inherit the theme: name the app's font.
const _sans = 'Montserrat';

void _text(Canvas canvas, String s, Offset at, double width, TextStyle style, {TextAlign align = TextAlign.center}) {
  final tp = TextPainter(text: TextSpan(text: s, style: style), textDirection: TextDirection.ltr, maxLines: 1, ellipsis: '…', textAlign: align)
    ..layout(minWidth: width, maxWidth: width);
  tp.paint(canvas, at);
}

/// Heat colours of the table (app palette): completeness in the app blues,
/// missing data and gaps from pale blue to crimson.
Color _fullness(double f) => Color.lerp(SplaceBlue.c50, SplaceBlue.c600, 0.12 + 0.88 * f.clamp(0.0, 1.0))!;
Color _missing(double f) {
  final v = f.clamp(0.0, 1.0);
  return v <= 0.5
      ? Color.lerp(SplacePalette.sky100, SplacePalette.orange, v / 0.5)!
      : Color.lerp(SplacePalette.orange, SplacePalette.red900, (v - 0.5) / 0.5)!;
}
Color _ink(Color c) => c.computeLuminance() > 0.4 ? const Color(0xFF1F2937) : Colors.white;

/// A value on its heat colour, in a rounded cell.
void _heatCell(Canvas canvas, Rect r, String text, Color fill) {
  canvas.drawRRect(RRect.fromRectAndRadius(r, const Radius.circular(3)), Paint()..color = fill);
  final style = TextStyle(fontSize: 10.5, fontWeight: FontWeight.w600, color: _ink(fill), fontFamily: _sans);
  _text(canvas, text, Offset(r.left + 2, r.top + (r.height - 13) / 2), r.width - 4, style);
}

/// Visible rows for a body of [height] scrolled by [y].
(int, int) _rows(int count, double y, double height) {
  const rowH = _MatrixOverviewTableState.rowH, headerH = _MatrixOverviewTableState.headerH;
  final first = math.max(0, (y / rowH).floor());
  final last = math.min(count - 1, ((y + height - headerH) / rowH).ceil());
  return (first, last);
}

/// Taxon, loci present, bases and missing, with the header kept on top.
class _FixedPainter extends CustomPainter {
  _FixedPainter(this.o, this.y, this.hover, this.t, {required this.nameW, required this.heads, required this.order})
      : super(repaint: Listenable.merge([y, hover]));
  final MatrixOverview o;

  /// Display row → taxon.
  final List<int> order;
  final ValueNotifier<double> y;
  final ValueNotifier<(int, int)?> hover;
  final SplaceTokens t;
  final double nameW;
  final List<String> heads;

  static const rowH = _MatrixOverviewTableState.rowH, headerH = _MatrixOverviewTableState.headerH;
  static const statW = _MatrixOverviewTableState.statW;

  @override
  void paint(Canvas canvas, Size size) {
    canvas.save();
    canvas.clipRect(Rect.fromLTWH(0, headerH, size.width, size.height - headerH));
    final (first, last) = _rows(o.taxa.length, y.value, size.height);
    final h = hover.value;
    final name = TextStyle(fontSize: 12, color: t.text, fontFamily: _sans);
    final maxBases = [for (var i = 0; i < o.taxa.length; i++) o.basesOf(i)].fold(1, math.max);
    for (var r = first; r <= last; r++) {
      final top = headerH + r * rowH - y.value;
      if (h != null && h.$1 == r) {
        canvas.drawRect(Rect.fromLTWH(0, top, size.width, rowH), Paint()..color = t.accentSoft);
      } else if (r.isOdd) {
        canvas.drawRect(Rect.fromLTWH(0, top, size.width, rowH), Paint()..color = t.surface2);
      }
      final i = order[r];
      _text(canvas, o.taxa[i], Offset(10, top + 4), nameW - 16, name, align: TextAlign.left);
      // Markers, bases and missing data of the taxon as heat cells.
      final missing = 1 - o.basesOf(i) / math.max(1, o.length);
      final stats = [
        ('${o.lociOf(i)}/${o.loci.length}', _fullness(o.lociOf(i) / math.max(1, o.loci.length))),
        ('${o.basesOf(i)}', _fullness(o.basesOf(i) / math.max(1, maxBases))),
        ('${(100 * missing).toStringAsFixed(1)}%', _missing(missing)),
      ];
      for (final (k, (v, c)) in stats.indexed) {
        _heatCell(canvas, Rect.fromLTWH(nameW + k * statW + 2, top + 2, statW - 6, rowH - 4), v, c);
      }
    }
    canvas.restore();
    final head = TextStyle(fontSize: 11, fontWeight: FontWeight.w700, color: t.text, fontFamily: _sans);
    _text(canvas, heads[0], const Offset(10, headerH - 20), nameW - 16, head, align: TextAlign.left);
    for (var i = 1; i < heads.length; i++) {
      _text(canvas, heads[i], Offset(nameW + (i - 1) * statW, headerH - 20), statW - 8, head, align: TextAlign.right);
    }
    canvas.drawLine(Offset(0, headerH - 0.5), Offset(size.width, headerH - 0.5), Paint()..color = t.border);
  }

  @override
  bool shouldRepaint(_FixedPainter old) =>
      old.o != o || old.t != t || old.nameW != nameW || old.heads.join() != heads.join() || !identical(old.order, order);
}

/// The loci: only the visible cells are drawn; the header stays on top.
class _GridPainter extends CustomPainter {
  _GridPainter(this.o, this.x, this.y, this.hover, this.t, this.colW, this.order, this.percent) : super(repaint: Listenable.merge([x, y, hover]));
  final MatrixOverview o;
  final List<int> order;

  /// Cells show the % of the marker's length instead of the bases.
  final bool percent;
  final ValueNotifier<double> x, y;
  final ValueNotifier<(int, int)?> hover;
  final SplaceTokens t;
  final double colW;

  static const rowH = _MatrixOverviewTableState.rowH, headerH = _MatrixOverviewTableState.headerH;

  @override
  void paint(Canvas canvas, Size size) {
    if (o.loci.isEmpty) return;
    final x0 = x.value, y0 = y.value;
    final first = math.max(0, (x0 / colW).floor());
    final last = math.min(o.loci.length - 1, ((x0 + size.width) / colW).ceil());
    final (r0, r1) = _rows(o.taxa.length, y0, size.height);
    final h = hover.value;
    final headStyle = TextStyle(fontSize: 11, fontWeight: FontWeight.w700, color: t.text, fontFamily: monoFont);
    final subStyle = TextStyle(fontSize: 10, color: t.textMuted, fontFamily: _sans);
    final dash = TextStyle(fontSize: 11, color: t.textFaint, fontFamily: _sans);
    final missing = Paint()..color = t.surface3;
    final stripe = Paint()..color = t.surface2;
    final hot = Paint()..color = t.accentSoft;

    // Body, clipped below the header.
    canvas.save();
    canvas.clipRect(Rect.fromLTWH(0, headerH, size.width, size.height - headerH));
    for (var r = r0; r <= r1; r++) {
      final top = headerH + r * rowH - y0;
      if (h != null && h.$1 == r && h.$2 < 0) {
        canvas.drawRect(Rect.fromLTWH(0, top, size.width, rowH), hot);
      } else if (r.isOdd) {
        canvas.drawRect(Rect.fromLTWH(0, top, size.width, rowH), stripe);
      }
    }
    for (var l = first; l <= last; l++) {
      final left = l * colW - x0;
      if (h != null && h.$2 == l) canvas.drawRect(Rect.fromLTWH(left, headerH, colW, size.height - headerH), hot);
      for (var r = r0; r <= r1; r++) {
        final top = headerH + r * rowH - y0;
        final b = o.bases[order[r]][l];
        final cell = RRect.fromRectAndRadius(Rect.fromLTWH(left + 2, top + 2, colW - 4, rowH - 4), const Radius.circular(3));
        if (b == 0) {
          canvas.drawRRect(cell, missing);
          _text(canvas, '—', Offset(left + 2, top + 4), colW - 4, dash);
          continue;
        }
        final f = (b / math.max(1, o.lengths[l])).clamp(0.0, 1.0);
        final fill = _fullness(f);
        canvas.drawRRect(cell, Paint()..color = fill);
        _text(canvas, percent ? '${(100 * f).round()}%' : '$b', Offset(left + 2, top + 4), colW - 4,
            TextStyle(fontSize: 11, fontWeight: FontWeight.w600, color: _ink(fill), fontFamily: monoFont));
        if (h != null && h.$1 == r && h.$2 == l) {
          canvas.drawRRect(
              cell,
              Paint()
                ..style = PaintingStyle.stroke
                ..strokeWidth = 2
                ..color = t.text);
        }
      }
    }
    canvas.restore();

    // Header: name, length, taxa with the locus and gaps.
    canvas.drawRect(Rect.fromLTWH(0, 0, size.width, headerH), Paint()..color = t.surface);
    for (var l = first; l <= last; l++) {
      final left = l * colW - x0;
      if (h != null && h.$2 == l) canvas.drawRect(Rect.fromLTWH(left, 0, colW, headerH), hot);
      final gap = o.gapFraction(l);
      _text(canvas, o.loci[l], Offset(left + 2, 4), colW - 4, headStyle);
      _text(canvas, '${o.lengths[l]} bp', Offset(left + 2, 19), colW - 4, subStyle);
      // Taxa with the marker (blue) and its gaps (pale to crimson).
      final half = (colW - 6) / 2;
      _heatCell(canvas, Rect.fromLTWH(left + 2, 31, half, 15), '${o.taxaWith(l)}/${o.taxa.length}', _fullness(o.taxaWith(l) / math.max(1, o.taxa.length)));
      _heatCell(canvas, Rect.fromLTWH(left + 4 + half, 31, half, 15), '${(100 * gap).round()}%', _missing(gap));
    }
    canvas.drawLine(Offset(0, headerH - 0.5), Offset(size.width, headerH - 0.5), Paint()..color = t.border);
  }

  @override
  bool shouldRepaint(_GridPainter old) => old.o != o || old.t != t || old.colW != colW || old.percent != percent || !identical(old.order, order);
}
