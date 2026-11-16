import 'dart:math' as math;

import 'package:flutter/material.dart';

import '../../core/saturation.dart';
import '../common.dart';
import '../theme.dart';
import 'examples.dart' show codonPositionColors;

/// Blue and red of the app palette, used by the saturation plots.
const plotBlue = SplaceBlue.c500;
const plotBlueLight = SplaceBlue.c300;
const plotRed = SplaceRed.c500;

class PlotSeries {
  const PlotSeries(this.label, this.color, this.points);
  final String label;
  final Color color;

  /// (x, y) pairs.
  final List<(double, double)> points;
}

/// Scatter plot with optional y = x reference line, drawn on one canvas so
/// thousands of pairwise points render without per-point widgets. Only the
/// left and bottom axes are drawn; the grid marks every tick.
class ScatterPlot extends StatelessWidget {
  const ScatterPlot({
    super.key,
    required this.title,
    required this.xLabel,
    required this.yLabel,
    required this.series,
    this.identityLine = false,
    this.trend = true,
    this.subtitle,
  });

  final String title;
  final String? subtitle;
  final String xLabel;
  final String yLabel;
  final List<PlotSeries> series;
  final bool identityLine;

  /// Least-squares line of each series (its slope is the "Slope" of the
  /// table).
  final bool trend;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    // Title and explanation on top, legend and plot at the bottom: in a row
    // of plots of equal height (EqualHeightRow) the plots line up and no
    // text is cut.
    return Column(crossAxisAlignment: CrossAxisAlignment.stretch, mainAxisAlignment: MainAxisAlignment.spaceBetween, children: [
      Column(crossAxisAlignment: CrossAxisAlignment.stretch, mainAxisSize: MainAxisSize.min, children: [
        Text(title, style: TextStyle(fontSize: 13.5, fontWeight: FontWeight.w700, color: t.text)),
        if (subtitle != null) ...[
          const SizedBox(height: 2),
          Text(subtitle!, style: TextStyle(fontSize: 12, color: t.textMuted, height: 1.35)),
        ],
      ]),
      Column(crossAxisAlignment: CrossAxisAlignment.stretch, mainAxisSize: MainAxisSize.min, children: [
      const SizedBox(height: 6),
      Wrap(spacing: 14, runSpacing: 4, children: [
        for (final s in series)
          Row(mainAxisSize: MainAxisSize.min, children: [
            Container(width: 10, height: 10, decoration: BoxDecoration(color: s.color, shape: BoxShape.circle)),
            const SizedBox(width: 5),
            Text(s.label, style: TextStyle(fontSize: 12, color: t.textMuted)),
          ]),
        if (trend)
          Row(mainAxisSize: MainAxisSize.min, children: [
            Container(width: 16, height: 1.6, color: t.textMuted),
            const SizedBox(width: 5),
            Text(context.l10n.plotTrend, style: TextStyle(fontSize: 12, color: t.textMuted)),
          ]),
        if (identityLine)
          Row(mainAxisSize: MainAxisSize.min, children: [
            SizedBox(width: 16, child: CustomPaint(size: const Size(16, 10), painter: _DashSample(t.textMuted))),
            const SizedBox(width: 5),
            Text('y = x', style: TextStyle(fontSize: 12, color: t.textMuted)),
          ]),
      ]),
      const SizedBox(height: 4),
      AspectRatio(
        aspectRatio: 1.3,
        child: CustomPaint(
          painter: _ScatterPainter(
            series: series,
            identityLine: identityLine,
            trend: trend,
            xLabel: xLabel,
            yLabel: yLabel,
            axis: t.textMuted,
            axisLine: t.text.withValues(alpha: 0.7),
            grid: t.textFaint.withValues(alpha: 0.32),
          ),
        ),
      ),
      ]),
    ]);
  }
}

class _DashSample extends CustomPainter {
  _DashSample(this.color);
  final Color color;

  @override
  void paint(Canvas canvas, Size size) {
    final p = Paint()
      ..color = color
      ..strokeWidth = 1.2;
    for (var x = 0.0; x < size.width; x += 6) {
      canvas.drawLine(Offset(x, size.height / 2), Offset(math.min(size.width, x + 3.5), size.height / 2), p);
    }
  }

  @override
  bool shouldRepaint(_DashSample old) => old.color != color;
}

class _ScatterPainter extends CustomPainter {
  _ScatterPainter({
    required this.series,
    required this.identityLine,
    required this.trend,
    required this.xLabel,
    required this.yLabel,
    required this.axis,
    required this.axisLine,
    required this.grid,
  });

  final List<PlotSeries> series;
  final bool identityLine, trend;
  final String xLabel, yLabel;
  final Color axis, axisLine, grid;

  static double _niceMax(double v) {
    if (v <= 0) return 1;
    final exp = math.pow(10, (math.log(v) / math.ln10).floor()).toDouble();
    for (final m in [1, 1.25, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10]) {
      if (m * exp >= v) return m * exp;
    }
    return 10 * exp;
  }

  void _text(Canvas canvas, String s, Offset at, {TextAlign align = TextAlign.center, bool vertical = false, bool bold = false}) {
    final tp = TextPainter(
      text: TextSpan(
        text: s,
        style: TextStyle(fontFamily: 'Montserrat', fontSize: 10.5, color: axis, fontWeight: bold ? FontWeight.w600 : FontWeight.w400),
      ),
      textDirection: TextDirection.ltr,
    )..layout();
    canvas.save();
    canvas.translate(at.dx, at.dy);
    if (vertical) canvas.rotate(-math.pi / 2);
    final dx = align == TextAlign.center ? -tp.width / 2 : (align == TextAlign.right ? -tp.width : 0.0);
    tp.paint(canvas, Offset(dx, -tp.height / 2));
    canvas.restore();
  }

  @override
  void paint(Canvas canvas, Size size) {
    const left = 48.0, bottom = 36.0, top = 8.0, right = 10.0, ticks = 5;
    final plot = Rect.fromLTRB(left, top, size.width - right, size.height - bottom);
    final all = [for (final s in series) ...s.points];
    final xMax = _niceMax(all.fold<double>(0, (a, p) => math.max(a, p.$1)));
    final yMax = _niceMax(all.fold<double>(0, (a, p) => math.max(a, p.$2)));
    String label(double v, double max) => v.toStringAsFixed(max < 0.1 ? 3 : (max < 1 ? 2 : 1));

    // Grid at every tick, between the axes.
    final gridPaint = Paint()
      ..color = grid
      ..strokeWidth = 0.8;
    for (var i = 0; i <= ticks; i++) {
      final x = plot.left + plot.width * i / ticks, y = plot.bottom - plot.height * i / ticks;
      if (i > 0) {
        canvas.drawLine(Offset(x, plot.top), Offset(x, plot.bottom), gridPaint);
        canvas.drawLine(Offset(plot.left, y), Offset(plot.right, y), gridPaint);
      }
      _text(canvas, label(xMax * i / ticks, xMax), Offset(x, plot.bottom + 10));
      _text(canvas, label(yMax * i / ticks, yMax), Offset(plot.left - 6, y), align: TextAlign.right);
    }
    _text(canvas, xLabel, Offset(plot.center.dx, size.height - 8), bold: true);
    _text(canvas, yLabel, Offset(10, plot.center.dy), vertical: true, bold: true);

    Offset map(double x, double y) => Offset(plot.left + plot.width * (x / xMax), plot.bottom - plot.height * (y / yMax));

    canvas.save();
    canvas.clipRect(plot);
    if (identityLine) {
      final m = math.min(xMax, yMax);
      final a = map(0, 0), b = map(m, m);
      final dash = Paint()
        ..color = axis
        ..strokeWidth = 1.2;
      final len = (b - a).distance;
      final dir = (b - a) / len;
      for (var d = 0.0; d < len; d += 7) {
        canvas.drawLine(a + dir * d, a + dir * math.min(len, d + 4), dash);
      }
    }
    for (final s in series) {
      final paint = Paint()..color = s.color.withValues(alpha: 0.5);
      for (final (x, y) in s.points) {
        canvas.drawCircle(map(x, y), 2.4, paint);
      }
    }
    if (trend) {
      for (final s in series) {
        final n = s.points.length;
        if (n < 3) continue;
        var sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0, x0 = double.infinity, x1 = 0.0;
        for (final (x, y) in s.points) {
          sx += x;
          sy += y;
          sxx += x * x;
          sxy += x * y;
          x0 = math.min(x0, x);
          x1 = math.max(x1, x);
        }
        final d = n * sxx - sx * sx;
        if (d.abs() < 1e-12) continue;
        final b = (n * sxy - sx * sy) / d, a = (sy - b * sx) / n;
        canvas.drawLine(
          map(x0, a + b * x0),
          map(x1, a + b * x1),
          Paint()
            ..color = s.color
            ..strokeWidth = 1.6,
        );
      }
    }
    canvas.restore();
    // Left and bottom axes only.
    final axisPaint = Paint()
      ..color = axisLine
      ..strokeWidth = 1.2;
    canvas.drawLine(plot.bottomLeft, plot.topLeft, axisPaint);
    canvas.drawLine(plot.bottomLeft, plot.bottomRight, axisPaint);
    for (var i = 0; i <= ticks; i++) {
      final x = plot.left + plot.width * i / ticks, y = plot.bottom - plot.height * i / ticks;
      canvas.drawLine(Offset(x, plot.bottom), Offset(x, plot.bottom + 4), axisPaint);
      canvas.drawLine(Offset(plot.left - 4, y), Offset(plot.left, y), axisPaint);
    }
  }

  @override
  bool shouldRepaint(_ScatterPainter old) => old.series != series || old.axis != axis || old.grid != grid || old.trend != trend;
}

/// The two classic saturation plots for one codon-position subset.
List<PlotSeries> observedVsCorrected(SubsetSaturation s, [Color color = plotBlue]) =>
    [PlotSeries('p', color, [for (final p in s.points) (p.distance, p.p)])];

List<PlotSeries> transitionsTransversions(SubsetSaturation s, String ts, String tv) => [
      PlotSeries(ts, plotRed, [for (final p in s.points) (p.distance, p.ts)]),
      PlotSeries(tv, plotBlue, [for (final p in s.points) (p.distance, p.tv)]),
    ];

/// p of each codon position against the all-positions distance.
List<PlotSeries> positionsVsDistance(SaturationReport r, String first, String second, String third) => [
      PlotSeries(first, codonPositionColors[0], [for (final p in r.positionPairs) (p.distance, p.p1)]),
      PlotSeries(second, codonPositionColors[1], [for (final p in r.positionPairs) (p.distance, p.p2)]),
      PlotSeries(third, codonPositionColors[2], [for (final p in r.positionPairs) (p.distance, p.p3)]),
    ];
