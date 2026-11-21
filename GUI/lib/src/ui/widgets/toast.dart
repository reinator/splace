import 'dart:async';
import 'dart:math' as math;

import 'package:flutter/material.dart';

import '../fa.dart';
import '../theme.dart';

/// A notice in the app's colours at the bottom right of the window: an
/// animated badge (a check drawn inside a ring that fills, red for an
/// error), the [title], a [detail] line (a file name, in mono) and a
/// [path] under it; a bar in the progress colours runs out until it goes.
/// [action] opens what was saved (a folder). Replaces the one shown before.
void showSplaceToast(
  BuildContext context, {
  required String title,
  String? detail,
  String? path,
  bool error = false,
  ({String label, VoidCallback onPressed})? action,
  Duration duration = const Duration(milliseconds: 4200),
}) {
  final overlay = Overlay.maybeOf(context, rootOverlay: true);
  if (overlay == null) return;
  _current?.remove();
  late final OverlayEntry entry;
  entry = OverlayEntry(
    builder: (_) => _Toast(
      title: title,
      detail: detail,
      path: path,
      error: error,
      action: action,
      duration: duration,
      onGone: () {
        if (identical(_current, entry)) _current = null;
        if (entry.mounted) entry.remove();
      },
    ),
  );
  _current = entry;
  overlay.insert(entry);
}

OverlayEntry? _current;

class _Toast extends StatefulWidget {
  const _Toast({
    required this.title,
    required this.detail,
    required this.path,
    required this.error,
    required this.action,
    required this.duration,
    required this.onGone,
  });
  final String title;
  final String? detail, path;
  final bool error;
  final ({String label, VoidCallback onPressed})? action;
  final Duration duration;
  final VoidCallback onGone;

  @override
  State<_Toast> createState() => _ToastState();
}

class _ToastState extends State<_Toast> with TickerProviderStateMixin {
  // In and out (slide + fade), the badge (ring, then check), the time left.
  late final _show = AnimationController(vsync: this, duration: const Duration(milliseconds: 320));
  late final _badge = AnimationController(vsync: this, duration: const Duration(milliseconds: 900));
  late final _life = AnimationController(vsync: this, duration: widget.duration);
  bool _closing = false;

  @override
  void initState() {
    super.initState();
    _show.forward();
    _badge.forward();
    _life.forward().then((_) => _close());
  }

  Future<void> _close() async {
    if (_closing || !mounted) return;
    _closing = true;
    await _show.reverse();
    widget.onGone();
  }

  @override
  void dispose() {
    _show.dispose();
    _badge.dispose();
    _life.dispose();
    super.dispose();
  }

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final tone = widget.error ? SplaceTokens.danger : t.accent;
    final slide = CurvedAnimation(parent: _show, curve: Curves.easeOutCubic, reverseCurve: Curves.easeInCubic);
    return Positioned(
      right: 24,
      bottom: 24,
      child: FadeTransition(
        opacity: slide,
        child: SlideTransition(
          position: Tween(begin: const Offset(0, 0.35), end: Offset.zero).animate(slide),
          child: MouseRegion(
            // The time stops while the pointer is on it.
            onEnter: (_) => _life.stop(),
            onExit: (_) {
              if (!_closing) _life.forward().then((_) => _close());
            },
            child: Material(
              type: MaterialType.transparency,
              child: Container(
                width: 400,
                clipBehavior: Clip.antiAlias,
                decoration: BoxDecoration(
                  color: t.surface,
                  borderRadius: BorderRadius.circular(14),
                  border: Border.all(color: widget.error ? SplaceTokens.danger.withValues(alpha: 0.45) : t.accentSoftBorder, width: 1.5),
                  boxShadow: shadowMd,
                ),
                child: Column(mainAxisSize: MainAxisSize.min, children: [
                  Padding(
                    padding: const EdgeInsets.fromLTRB(14, 14, 8, 12),
                    child: Row(crossAxisAlignment: CrossAxisAlignment.start, children: [
                      AnimatedBuilder(
                        animation: _badge,
                        builder: (context, _) => CustomPaint(
                          size: const Size.square(40),
                          painter: _BadgePainter(_badge.value, tone: tone, soft: widget.error ? SplaceRed.c50 : t.accentSoft, error: widget.error),
                        ),
                      ),
                      const SizedBox(width: 12),
                      Expanded(
                        child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
                          const SizedBox(height: 1),
                          Text(widget.title, style: TextStyle(fontSize: 14, fontWeight: FontWeight.w700, color: t.text)),
                          if (widget.detail != null) ...[
                            const SizedBox(height: 3),
                            Text(widget.detail!,
                                maxLines: 1,
                                overflow: TextOverflow.ellipsis,
                                style: TextStyle(fontFamily: monoFont, fontSize: 12, fontWeight: FontWeight.w600, color: widget.error ? SplaceTokens.danger : t.accentText)),
                          ],
                          if (widget.path != null) ...[
                            const SizedBox(height: 2),
                            Text(widget.path!, maxLines: 2, overflow: TextOverflow.ellipsis, style: TextStyle(fontFamily: monoFont, fontSize: 11, color: t.textFaint)),
                          ],
                          if (widget.action != null) ...[
                            const SizedBox(height: 8),
                            InkWell(
                              borderRadius: BorderRadius.circular(6),
                              onTap: () {
                                widget.action!.onPressed();
                                _close();
                              },
                              child: Padding(
                                padding: const EdgeInsets.symmetric(vertical: 2),
                                child: Row(mainAxisSize: MainAxisSize.min, children: [
                                  Icon(FA.folderOpen, size: 11, color: t.accent),
                                  const SizedBox(width: 6),
                                  Text(widget.action!.label, style: TextStyle(fontSize: 12, fontWeight: FontWeight.w600, color: t.accent)),
                                ]),
                              ),
                            ),
                          ],
                        ]),
                      ),
                      IconButton(
                        tooltip: MaterialLocalizations.of(context).closeButtonTooltip,
                        visualDensity: VisualDensity.compact,
                        iconSize: 12,
                        color: t.textFaint,
                        onPressed: _close,
                        icon: const Icon(FA.xmark),
                      ),
                    ]),
                  ),
                  // Time left, in the progress colours (red for an error).
                  AnimatedBuilder(
                    animation: _life,
                    builder: (context, _) => Align(
                      alignment: Alignment.centerLeft,
                      child: FractionallySizedBox(
                        widthFactor: 1 - _life.value,
                        child: Container(
                          height: 3,
                          decoration: BoxDecoration(
                            gradient: widget.error ? null : splaceProgress,
                            color: widget.error ? SplaceTokens.danger : null,
                          ),
                        ),
                      ),
                    ),
                  ),
                ]),
              ),
            ),
          ),
        ),
      ),
    );
  }
}

/// The badge: a soft disc, a ring drawn around it, then the check (or the
/// cross of an error) traced inside, with a small pop at the end.
class _BadgePainter extends CustomPainter {
  _BadgePainter(this.v, {required this.tone, required this.soft, required this.error});
  final double v;
  final Color tone, soft;
  final bool error;

  @override
  void paint(Canvas canvas, Size size) {
    final c = size.center(Offset.zero);
    final ring = Curves.easeInOut.transform((v / 0.55).clamp(0.0, 1.0));
    final mark = Curves.easeOut.transform(((v - 0.45) / 0.4).clamp(0.0, 1.0));
    final pop = v < 0.85 ? 1.0 : 1 + 0.08 * math.sin((v - 0.85) / 0.15 * math.pi);
    final r = size.width / 2 * pop;
    canvas.drawCircle(c, r, Paint()..color = soft);
    final stroke = Paint()
      ..color = tone
      ..style = PaintingStyle.stroke
      ..strokeWidth = 2.4
      ..strokeCap = StrokeCap.round;
    canvas.drawArc(Rect.fromCircle(center: c, radius: r - 2), -math.pi / 2, 2 * math.pi * ring, false, stroke);
    if (mark == 0) return;
    final s = size.width / 40 * pop;
    final path = Path();
    const segs = [
      [Offset(-6, -6), Offset(6, 6)],
      [Offset(6, -6), Offset(-6, 6)],
    ];
    const check = [
      [Offset(-7, 0.5), Offset(-2, 5.5), Offset(7.5, -5)],
    ];
    final lines = error ? segs : check;
    final total = lines.length;
    for (var i = 0; i < total; i++) {
      final part = ((mark * total) - i).clamp(0.0, 1.0);
      if (part == 0) continue;
      final pts = [for (final o in lines[i]) c + o * s];
      final lens = [for (var k = 1; k < pts.length; k++) (pts[k] - pts[k - 1]).distance];
      var left = lens.fold<double>(0, (a, b) => a + b) * part;
      path.moveTo(pts[0].dx, pts[0].dy);
      for (var k = 1; k < pts.length && left > 0; k++) {
        final f = math.min(1.0, left / lens[k - 1]);
        final p = Offset.lerp(pts[k - 1], pts[k], f)!;
        path.lineTo(p.dx, p.dy);
        left -= lens[k - 1];
      }
    }
    canvas.drawPath(path, stroke..strokeWidth = 2.8);
  }

  @override
  bool shouldRepaint(_BadgePainter old) => old.v != v || old.tone != tone || old.soft != soft;
}
