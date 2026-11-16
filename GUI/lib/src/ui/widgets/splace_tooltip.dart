import 'dart:async';

import 'package:flutter/material.dart';

/// Tooltip with an arrow, styled like the HTML app's Tippy.js "splace"
/// theme: SPLACE blue bubble (dark slate in dark mode), placed above the
/// target and flipped below when there is no room.
class SplaceTooltip extends StatefulWidget {
  const SplaceTooltip({super.key, this.message = '', this.richMessage, required this.child, this.waitDuration = const Duration(milliseconds: 250)});

  final String message;
  final InlineSpan? richMessage;
  final Widget child;
  final Duration waitDuration;

  @override
  State<SplaceTooltip> createState() => _SplaceTooltipState();
}

class _SplaceTooltipState extends State<SplaceTooltip> {
  final _portal = OverlayPortalController();
  Timer? _timer;

  bool get _empty => widget.richMessage == null && widget.message.trim().isEmpty;

  void _show() {
    if (_empty) return;
    _timer?.cancel();
    _timer = Timer(widget.waitDuration, () {
      if (mounted && !_portal.isShowing) _portal.show();
    });
  }

  void _hide() {
    _timer?.cancel();
    if (_portal.isShowing) _portal.hide();
  }

  @override
  void dispose() {
    _timer?.cancel();
    super.dispose();
  }

  @override
  Widget build(BuildContext context) {
    if (_empty) return widget.child;
    return OverlayPortal(
      controller: _portal,
      overlayChildBuilder: (overlayContext) {
        final box = context.findRenderObject() as RenderBox?;
        if (box == null || !box.attached) return const SizedBox.shrink();
        final target = box.localToGlobal(Offset.zero) & box.size;
        return Positioned.fill(child: IgnorePointer(child: _TooltipBubble(target: target, message: widget.message, richMessage: widget.richMessage)));
      },
      child: MouseRegion(
        onEnter: (_) => _show(),
        onExit: (_) => _hide(),
        child: Listener(
          behavior: HitTestBehavior.translucent,
          onPointerDown: (_) => _hide(),
          child: GestureDetector(behavior: HitTestBehavior.translucent, onLongPress: () => _portal.show(), child: widget.child),
        ),
      ),
    );
  }
}

class _TooltipBubble extends StatelessWidget {
  const _TooltipBubble({required this.target, required this.message, this.richMessage});
  final Rect target;
  final String message;
  final InlineSpan? richMessage;

  @override
  Widget build(BuildContext context) {
    final dark = Theme.of(context).brightness == Brightness.dark;
    final bg = dark ? const Color(0xFF1B2234) : const Color(0xFF323795);
    final fg = dark ? const Color(0xFFEEF2FF) : Colors.white;
    final text = DefaultTextStyle(
      style: TextStyle(fontFamily: 'Montserrat', fontSize: 13, height: 1.5, color: fg, decoration: TextDecoration.none, fontWeight: FontWeight.w500),
      // Justified like the HTML tooltips; the last line stays left-aligned.
      child: richMessage != null
          ? Text.rich(richMessage!, textAlign: TextAlign.justify)
          : Text(message, textAlign: TextAlign.justify),
    );
    return TweenAnimationBuilder<double>(
      tween: Tween(begin: 0, end: 1),
      duration: const Duration(milliseconds: 140),
      builder: (context, v, child) => Opacity(opacity: v, child: child),
      child: CustomMultiChildLayout(
        delegate: _TooltipLayout(target),
        children: [
          LayoutId(
            id: _Part.bubble,
            child: ConstrainedBox(
              // Wider than Tippy's 320 so justified lines have fewer gaps.
              constraints: const BoxConstraints(maxWidth: 360),
              child: Container(
                padding: const EdgeInsets.symmetric(horizontal: 12, vertical: 7),
                decoration: BoxDecoration(
                  color: bg,
                  borderRadius: BorderRadius.circular(6),
                  border: Border.all(color: dark ? const Color(0x38949FB8) : const Color(0x1FFFFFFF)),
                  boxShadow: [BoxShadow(color: dark ? const Color(0x7A020617) : const Color(0x59323795), blurRadius: dark ? 30 : 16, offset: const Offset(0, 4))],
                ),
                child: text,
              ),
            ),
          ),
          LayoutId(id: _Part.arrowDown, child: CustomPaint(size: _arrow, painter: _ArrowPainter(bg, down: true))),
          LayoutId(id: _Part.arrowUp, child: CustomPaint(size: _arrow, painter: _ArrowPainter(bg, down: false))),
        ],
      ),
    );
  }
}

const _arrow = Size(14, 7);

enum _Part { bubble, arrowDown, arrowUp }

/// Bubble above the target (arrow pointing down), or below it (arrow
/// pointing up) when there is no room; kept inside the window.
class _TooltipLayout extends MultiChildLayoutDelegate {
  _TooltipLayout(this.target);
  final Rect target;

  @override
  void performLayout(Size size) {
    const gap = 4.0, margin = 8.0;
    final bubble = layoutChild(_Part.bubble, BoxConstraints.loose(Size(size.width - 2 * margin, size.height)));
    layoutChild(_Part.arrowDown, BoxConstraints.tight(_arrow));
    layoutChild(_Part.arrowUp, BoxConstraints.tight(_arrow));
    final above = target.top - gap - _arrow.height - bubble.height >= margin;
    final x = (target.center.dx - bubble.width / 2).clamp(margin, size.width - bubble.width - margin);
    final y = above ? target.top - gap - _arrow.height - bubble.height : target.bottom + gap + _arrow.height;
    positionChild(_Part.bubble, Offset(x, y));
    final ax = (target.center.dx - _arrow.width / 2).clamp(x + 6, x + bubble.width - _arrow.width - 6);
    const hidden = Offset(-1000, -1000);
    positionChild(_Part.arrowDown, above ? Offset(ax, y + bubble.height - 0.5) : hidden);
    positionChild(_Part.arrowUp, above ? hidden : Offset(ax, y - _arrow.height + 0.5));
  }

  @override
  bool shouldRelayout(_TooltipLayout old) => old.target != target;
}

class _ArrowPainter extends CustomPainter {
  _ArrowPainter(this.color, {required this.down});
  final Color color;
  final bool down;

  @override
  void paint(Canvas canvas, Size size) {
    final path = down
        ? (Path()
          ..moveTo(0, 0)
          ..lineTo(size.width, 0)
          ..lineTo(size.width / 2, size.height)
          ..close())
        : (Path()
          ..moveTo(0, size.height)
          ..lineTo(size.width, size.height)
          ..lineTo(size.width / 2, 0)
          ..close());
    canvas.drawPath(path, Paint()..color = color);
  }

  @override
  bool shouldRepaint(_ArrowPainter old) => old.color != color || old.down != down;
}
