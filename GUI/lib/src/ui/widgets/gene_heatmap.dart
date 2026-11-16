import 'dart:math' as math;

import 'package:flutter/gestures.dart';
import 'package:flutter/material.dart';

import '../theme.dart';

/// Rows × columns of a marker heatmap. `values[row][col]` is the length
/// (bp) of the marker in that sample, or null when absent.
class HeatmapData {
  HeatmapData({
    required this.rows,
    required this.rowItalic,
    required this.rowDetail,
    required this.cols,
    required this.values,
    required this.colActive,
  });

  final List<String> rows;
  final List<bool> rowItalic;
  final List<String?> rowDetail;
  final List<String> cols;
  final List<List<int?>> values;
  final List<bool> colActive;

  /// Identity for repaint checks.
  late final String signature = '${rows.length}x${cols.length}:${colActive.map((a) => a ? 1 : 0).join()}:${cols.join('|')}';

  /// Samples carrying each marker.
  late final List<int> colPresent = [
    for (var c = 0; c < cols.length; c++) values.where((r) => r[c] != null).length,
  ];

  /// Markers present in each sample.
  late final List<int> rowPresent = [for (final r in values) r.where((v) => v != null).length];
}

enum HeatmapMode {
  /// Gene Presence Heatmap: blue present, orange absent.
  presence,

  /// Marker Length Heatmap: each marker scaled by its own min–max length.
  length,
}

/// Heatmap with frozen headers: marker names stay on top and sample names
/// on the left while the grid scrolls, so large datasets keep their
/// labels. The hovered row and column are highlighted and described in the
/// bar above; clicking a marker name toggles it.
class GeneHeatmap extends StatefulWidget {
  const GeneHeatmap({super.key, required this.data, this.mode = HeatmapMode.presence, this.onToggleColumn, this.hint = ' '});

  final HeatmapData data;

  /// Shown in the info bar until something is hovered.
  final String hint;
  final HeatmapMode mode;
  final ValueChanged<String>? onToggleColumn;

  static const present = SplacePalette.blue600;
  static const absent = SplacePalette.red500;
  static const lengthLow = Color(0xFFEEF0FB);
  static const lengthHigh = Color(0xFF323795);

  @override
  State<GeneHeatmap> createState() => _GeneHeatmapState();
}

class _GeneHeatmapState extends State<GeneHeatmap> {
  static const headerHeight = 104.0;
  final _gridH = ScrollController(), _headH = ScrollController();
  final _gridV = ScrollController(), _labelV = ScrollController();

  /// Hovered (row, column); either may be null (header or label hover).
  final _hover = ValueNotifier<(int?, int?)>((null, null));

  @override
  void initState() {
    super.initState();
    // Headers follow the grid.
    _gridH.addListener(() => _follow(_gridH, _headH));
    _gridV.addListener(() => _follow(_gridV, _labelV));
  }

  static void _follow(ScrollController from, ScrollController to) {
    if (to.hasClients && to.offset != from.offset) to.jumpTo(from.offset.clamp(0, to.position.maxScrollExtent));
  }

  @override
  void dispose() {
    for (final c in [_gridH, _headH, _gridV, _labelV]) {
      c.dispose();
    }
    _hover.dispose();
    super.dispose();
  }

  double _labelWidth(HeatmapData d) {
    var longest = 0;
    for (var i = 0; i < d.rows.length; i++) {
      longest = math.max(longest, d.rows[i].length + (d.rowDetail[i]?.length ?? 0) + 3);
    }
    return (longest * 7.2 + 24).clamp(160, 380);
  }

  /// Mouse wheel over the frozen labels scrolls the grid.
  void _wheel(PointerSignalEvent e, ScrollController c, {required bool vertical}) {
    if (e is! PointerScrollEvent || !c.hasClients) return;
    final delta = vertical ? e.scrollDelta.dy : (e.scrollDelta.dx != 0 ? e.scrollDelta.dx : e.scrollDelta.dy);
    c.jumpTo((c.offset + delta).clamp(0, c.position.maxScrollExtent));
  }

  @override
  Widget build(BuildContext context) => LayoutBuilder(builder: (context, box) => _build(context, box.maxWidth));

  Widget _build(BuildContext context, double width) {
    final d = widget.data;
    final t = context.tokens;
    final lw = _labelWidth(d);
    final length = widget.mode == HeatmapMode.length;
    // Minimum cell size: square in the presence map, wide enough for
    // "1545" in the length map. Cells then grow to fill the whole width
    // (any screen, ultrawide included) and rows grow a little with them.
    final small = d.cols.length > 60 ? 20.0 : (d.cols.length > 30 ? 24.0 : 30.0);
    final minW = length ? 46.0 : small;
    final room = width - lw - 1 - 14;
    final cw = d.cols.isEmpty ? minW : math.max(minW, room / d.cols.length);
    final ch = length ? (cw * 0.4).clamp(28.0, 40.0) : (cw * 0.6).clamp(small, 40.0);
    final gridSize = Size(cw * d.cols.length, ch * d.rows.length + 4);
    final maxViewport = math.max(320.0, MediaQuery.sizeOf(context).height * 0.62);
    final viewport = math.min(maxViewport, gridSize.height + 14);

    String describe((int?, int?) h) {
      final (r, c) = h;
      if (c != null && r == null) return '${d.cols[c]}: ${d.colPresent[c]} / ${d.rows.length}';
      if (r != null && c == null) {
        return '${d.rows[r]}${d.rowDetail[r] == null ? '' : ' (${d.rowDetail[r]})'}: ${d.rowPresent[r]} / ${d.cols.length}';
      }
      if (r != null && c != null) {
        final v = d.values[r][c];
        return '${d.rows[r]}${d.rowDetail[r] == null ? '' : ' (${d.rowDetail[r]})'} × ${d.cols[c]}: ${v == null ? '—' : '$v bp'}'
            '   ·   ${d.cols[c]} ${d.colPresent[c]} / ${d.rows.length}';
      }
      return widget.hint;
    }

    final colors = (text: t.text, muted: t.textFaint, accent: t.accentText, band: t.accentSoft, empty: t.surface3);
    int? colAt(double x, ScrollController c) {
      final i = (x + (c.hasClients ? c.offset : 0)) ~/ cw;
      return i >= 0 && i < d.cols.length ? i : null;
    }

    return Container(
      decoration: BoxDecoration(borderRadius: BorderRadius.circular(12), border: Border.all(color: t.border)),
      clipBehavior: Clip.antiAlias,
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        ValueListenableBuilder<(int?, int?)>(
          valueListenable: _hover,
          builder: (context, h, _) => Container(
            color: t.surface2,
            padding: const EdgeInsets.symmetric(horizontal: 12, vertical: 7),
            child: Text(describe(h),
                maxLines: 1, overflow: TextOverflow.ellipsis, style: TextStyle(fontSize: 12.5, color: h == (null, null) ? t.textMuted : t.text)),
          ),
        ),
        // Frozen header row: corner + marker names.
        SizedBox(
          height: headerHeight,
          child: Row(children: [
            Container(
              width: lw,
              color: t.surface,
              alignment: Alignment.bottomLeft,
              padding: const EdgeInsets.fromLTRB(12, 0, 8, 8),
              child: Text('${d.rows.length} × ${d.cols.length}', style: TextStyle(fontSize: 12, color: t.textMuted)),
            ),
            const SizedBox(width: 1),
            Expanded(
              child: Listener(
                onPointerSignal: (e) => _wheel(e, _gridH, vertical: false),
                child: MouseRegion(
                  onHover: (e) => _hover.value = (null, colAt(e.localPosition.dx, _headH)),
                  onExit: (_) => _hover.value = (null, null),
                  cursor: widget.onToggleColumn == null ? MouseCursor.defer : SystemMouseCursors.click,
                  child: GestureDetector(
                    onTapUp: (e) {
                      final c = colAt(e.localPosition.dx, _headH);
                      if (widget.onToggleColumn != null && c != null) widget.onToggleColumn!(d.cols[c]);
                    },
                    child: SingleChildScrollView(
                      controller: _headH,
                      scrollDirection: Axis.horizontal,
                      physics: const NeverScrollableScrollPhysics(),
                      child: ValueListenableBuilder<(int?, int?)>(
                        valueListenable: _hover,
                        builder: (context, h, _) => CustomPaint(
                          size: Size(gridSize.width + 14, headerHeight),
                          painter: _HeaderPainter(data: d, cw: cw, hoverCol: h.$2, colors: colors),
                        ),
                      ),
                    ),
                  ),
                ),
              ),
            ),
          ]),
        ),
        Divider(height: 1, color: t.border),
        // Frozen sample names + scrolling grid.
        SizedBox(
          height: viewport,
          child: Row(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
            SizedBox(
              width: lw,
              child: Listener(
                onPointerSignal: (e) => _wheel(e, _gridV, vertical: true),
                child: MouseRegion(
                  onHover: (e) {
                    final r = ((e.localPosition.dy + (_labelV.hasClients ? _labelV.offset : 0)) ~/ ch);
                    _hover.value = (r >= 0 && r < d.rows.length ? r : null, null);
                  },
                  onExit: (_) => _hover.value = (null, null),
                  child: SingleChildScrollView(
                    controller: _labelV,
                    physics: const NeverScrollableScrollPhysics(),
                    child: ValueListenableBuilder<(int?, int?)>(
                      valueListenable: _hover,
                      builder: (context, h, _) => CustomPaint(
                        size: Size(lw, gridSize.height + 14),
                        painter: _LabelPainter(data: d, ch: ch, width: lw, hoverRow: h.$1, colors: colors),
                      ),
                    ),
                  ),
                ),
              ),
            ),
            Container(width: 1, color: t.border),
            Expanded(
              child: Scrollbar(
                controller: _gridV,
                thumbVisibility: gridSize.height + 14 > viewport + 0.5,
                child: SingleChildScrollView(
                  controller: _gridV,
                  child: Scrollbar(
                    controller: _gridH,
                    thumbVisibility: gridSize.width > room + 0.5,
                    notificationPredicate: (n) => n.depth == 0,
                    child: SingleChildScrollView(
                      controller: _gridH,
                      scrollDirection: Axis.horizontal,
                      child: Padding(
                        padding: const EdgeInsets.only(right: 14, bottom: 14),
                        child: MouseRegion(
                          onHover: (e) {
                            final r = e.localPosition.dy ~/ ch, c = e.localPosition.dx ~/ cw;
                            final hit = r >= 0 && c >= 0 && r < d.rows.length && c < d.cols.length ? (r, c) : (null, null);
                            if (hit != _hover.value) _hover.value = hit;
                          },
                          onExit: (_) => _hover.value = (null, null),
                          child: Stack(children: [
                            RepaintBoundary(
                              child: CustomPaint(size: gridSize, painter: _CellPainter(data: d, mode: widget.mode, cw: cw, ch: ch, empty: t.surface3)),
                            ),
                            IgnorePointer(
                              child: ValueListenableBuilder<(int?, int?)>(
                                valueListenable: _hover,
                                builder: (context, h, _) => CustomPaint(size: gridSize, painter: _HoverPainter(hover: h, cw: cw, ch: ch, color: t.text, band: t.accentSoft)),
                              ),
                            ),
                          ]),
                        ),
                      ),
                    ),
                  ),
                ),
              ),
            ),
          ]),
        ),
      ]),
    );
  }
}

typedef _Colors = ({Color text, Color muted, Color accent, Color band, Color empty});

TextPainter _tp(InlineSpan span, double maxWidth) =>
    TextPainter(text: span, textDirection: TextDirection.ltr, maxLines: 1, ellipsis: '…')..layout(maxWidth: maxWidth);

/// Rotated marker names; the hovered one is bold on a band.
class _HeaderPainter extends CustomPainter {
  _HeaderPainter({required this.data, required this.cw, required this.hoverCol, required this.colors});
  final HeatmapData data;
  final double cw;
  final int? hoverCol;
  final _Colors colors;

  @override
  void paint(Canvas canvas, Size size) {
    final d = data;
    const hh = _GeneHeatmapState.headerHeight;
    if (hoverCol != null) canvas.drawRect(Rect.fromLTWH(hoverCol! * cw, 0, cw, hh), Paint()..color = colors.band);
    for (var c = 0; c < d.cols.length; c++) {
      final hover = c == hoverCol;
      final tp = _tp(
        TextSpan(
          text: d.cols[c],
          style: TextStyle(
            fontFamily: 'Montserrat',
            fontSize: cw >= 28 ? 13.5 : 12,
            fontWeight: hover ? FontWeight.w800 : FontWeight.w500,
            color: d.colActive[c] ? colors.accent : colors.muted,
          ),
        ),
        hh - 12,
      );
      canvas.save();
      canvas.translate(c * cw + cw / 2 - tp.height / 2, hh - 6);
      canvas.rotate(-math.pi / 2);
      tp.paint(canvas, Offset.zero);
      canvas.restore();
    }
  }

  @override
  bool shouldRepaint(_HeaderPainter old) => old.hoverCol != hoverCol || !identical(old.data, data) || old.cw != cw || old.colors != colors;
}

/// Sample names (italic species and the accession); hovered one on a band.
class _LabelPainter extends CustomPainter {
  _LabelPainter({required this.data, required this.ch, required this.width, required this.hoverRow, required this.colors});
  final HeatmapData data;
  final double ch, width;
  final int? hoverRow;
  final _Colors colors;

  @override
  void paint(Canvas canvas, Size size) {
    final d = data;
    if (hoverRow != null) canvas.drawRect(Rect.fromLTWH(0, hoverRow! * ch, width, ch), Paint()..color = colors.band);
    for (var r = 0; r < d.rows.length; r++) {
      final tp = _tp(
        TextSpan(children: [
          TextSpan(
            text: d.rows[r].replaceAll('_', ' '),
            style: TextStyle(
              fontFamily: 'Montserrat',
              fontSize: ch >= 28 ? 13 : 12,
              fontWeight: r == hoverRow ? FontWeight.w700 : FontWeight.w400,
              fontStyle: d.rowItalic[r] ? FontStyle.italic : FontStyle.normal,
              color: colors.text,
            ),
          ),
          if (d.rowDetail[r] != null) TextSpan(text: '  ${d.rowDetail[r]}', style: TextStyle(fontFamily: 'Montserrat', fontSize: 11, color: colors.muted)),
        ]),
        width - 20,
      );
      tp.paint(canvas, Offset(12, r * ch + (ch - tp.height) / 2));
    }
  }

  @override
  bool shouldRepaint(_LabelPainter old) => old.hoverRow != hoverRow || !identical(old.data, data) || old.ch != ch || old.colors != colors;
}

/// The coloured cells (drawn once into their own layer); the length map
/// writes each length (bp) inside its cell.
class _CellPainter extends CustomPainter {
  _CellPainter({required this.data, required this.mode, required this.cw, required this.ch, required this.empty});
  final HeatmapData data;
  final HeatmapMode mode;
  final double cw, ch;
  final Color empty;

  @override
  void paint(Canvas canvas, Size size) {
    final d = data;
    // Per-column ranges for the length mode.
    final lo = List<int>.filled(d.cols.length, 1 << 30), hi = List<int>.filled(d.cols.length, 0);
    for (final row in d.values) {
      for (var c = 0; c < d.cols.length; c++) {
        final v = row[c];
        if (v == null) continue;
        lo[c] = math.min(lo[c], v);
        hi[c] = math.max(hi[c], v);
      }
    }
    final gap = math.min(cw, ch) >= 28 ? 2.0 : 1.5;
    for (var r = 0; r < d.rows.length; r++) {
      for (var c = 0; c < d.cols.length; c++) {
        final v = d.values[r][c];
        Color color;
        var k = 0.0;
        if (mode == HeatmapMode.presence) {
          color = v != null ? GeneHeatmap.present : GeneHeatmap.absent;
        } else if (v == null) {
          color = empty;
        } else {
          k = hi[c] == lo[c] ? 1.0 : (v - lo[c]) / (hi[c] - lo[c]);
          color = Color.lerp(GeneHeatmap.lengthLow, GeneHeatmap.lengthHigh, 0.15 + 0.85 * k)!;
        }
        if (!d.colActive[c]) color = color.withValues(alpha: 0.25);
        final rect = Rect.fromLTWH(c * cw + gap, r * ch + gap, cw - 2 * gap, ch - 2 * gap);
        canvas.drawRRect(RRect.fromRectAndRadius(rect, Radius.circular(math.min(cw, ch) / 7)), Paint()..color = color);
        if (mode == HeatmapMode.length && v != null) {
          final dark = 0.15 + 0.85 * k > 0.5 && d.colActive[c];
          final tp = _tp(
            TextSpan(
              text: '$v',
              style: TextStyle(fontFamily: 'Montserrat', fontSize: (ch * 0.36).clamp(10.5, 13.5), fontWeight: FontWeight.w600, color: dark ? Colors.white : GeneHeatmap.lengthHigh),
            ),
            rect.width - 2,
          );
          tp.paint(canvas, Offset(rect.center.dx - tp.width / 2, rect.center.dy - tp.height / 2));
        }
      }
    }
  }

  @override
  bool shouldRepaint(_CellPainter old) =>
      old.data.signature != data.signature || old.mode != mode || !identical(old.data, data) || old.cw != cw || old.ch != ch;
}

/// Row/column bands and the outline of the hovered cell.
class _HoverPainter extends CustomPainter {
  _HoverPainter({required this.hover, required this.cw, required this.ch, required this.color, required this.band});
  final (int?, int?) hover;
  final double cw, ch;
  final Color color, band;

  @override
  void paint(Canvas canvas, Size size) {
    final (r, c) = hover;
    final p = Paint()..color = band.withValues(alpha: 0.35);
    if (r != null) canvas.drawRect(Rect.fromLTWH(0, r * ch, size.width, ch), p);
    if (c != null) canvas.drawRect(Rect.fromLTWH(c * cw, 0, cw, size.height), p);
    if (r != null && c != null) {
      canvas.drawRRect(
        RRect.fromRectAndRadius(Rect.fromLTWH(c * cw + 1, r * ch + 1, cw - 2, ch - 2), const Radius.circular(5)),
        Paint()
          ..color = color
          ..style = PaintingStyle.stroke
          ..strokeWidth = 2,
      );
    }
  }

  @override
  bool shouldRepaint(_HoverPainter old) => old.hover != hover || old.cw != cw || old.ch != ch;
}
