import 'dart:math' as math;

import 'package:flutter/material.dart';

import '../common.dart';
import '../fa.dart';
import '../theme.dart';
import 'kit.dart';

/// Sizes of the app's modals: every content modal is [standard] (one
/// width and one height for all), confirmations are [compact]. [full] is
/// kept as an alias of [standard].
enum DialogSize { standard, compact, full }

/// `.modal-panel` of the HTML app: gradient header with icon, title and
/// subtitle, scrollable body and a footer with actions.
class SplaceDialog extends StatelessWidget {
  const SplaceDialog({
    super.key,
    required this.icon,
    required this.title,
    this.subtitle,
    required this.child,
    this.actions = const [],
    this.footerLeading,
    this.size = DialogSize.standard,
    this.danger = false,
    this.scroll = true,
    this.fill = false,
  });

  static const standardWidth = 1400.0;
  static const maxHeight = 1000.0;
  static const compactWidth = 560.0;

  /// Width of a standard modal in a window of [screen] size.
  static double widthFor(Size screen) => math.min(standardWidth, screen.width - 48);

  /// Height of a standard modal in a window of [screen] size.
  static double heightFor(Size screen) => math.min(maxHeight, screen.height - 48);

  final IconData icon;
  final String title;
  final String? subtitle;
  final Widget child;
  final List<Widget> actions;

  /// Shown at the left of the footer (counts, status).
  final Widget? footerLeading;
  final DialogSize size;
  final bool danger;

  /// False when [child] scrolls by itself (tables).
  final bool scroll;

  /// Always use the full height (lists and tables); otherwise the modal
  /// hugs its content.
  final bool fill;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final screen = MediaQuery.sizeOf(context);
    // Content modals all share one size; only confirmations hug their text.
    final standard = size != DialogSize.compact;
    final full = standard;
    final height = standard ? heightFor(screen) : math.min(maxHeight, screen.height - 64);
    const pad = EdgeInsets.all(24);
    final body = scroll ? SingleChildScrollView(padding: pad, child: child) : Padding(padding: pad, child: child);
    final column = Column(mainAxisSize: fill || full ? MainAxisSize.max : MainAxisSize.min, crossAxisAlignment: CrossAxisAlignment.stretch, children: [
      DialogHeader(
        leading: DialogHeader.iconTile(context, icon, danger: danger),
        title: title,
        subtitle: subtitle,
        onClose: () => Navigator.of(context).pop(),
      ),
      if (fill || full) Expanded(child: body) else Flexible(child: body),
      if (actions.isNotEmpty || footerLeading != null)
        Container(
          padding: const EdgeInsets.symmetric(horizontal: 24, vertical: 14),
          decoration: BoxDecoration(color: t.surface2, border: Border(top: BorderSide(color: t.border))),
          child: Row(children: [
            Expanded(child: footerLeading ?? const SizedBox.shrink()),
            for (final (i, a) in actions.indexed) ...[if (i > 0) const SizedBox(width: 10), a],
          ]),
        ),
    ]);
    return Dialog(
      clipBehavior: Clip.antiAlias,
      insetPadding: const EdgeInsets.all(24),
      child: SizedBox(
        width: standard ? widthFor(screen) : math.min(compactWidth, screen.width - 64),
        height: fill || full ? height : null,
        child: ConstrainedBox(constraints: BoxConstraints(maxHeight: height), child: column),
      ),
    );
  }
}

/// Header of every modal: icon tile, title and subtitle on the modal
/// surface, separated from the body by the splash's animated rainbow line
/// instead of a coloured background.
class DialogHeader extends StatelessWidget {
  const DialogHeader({super.key, required this.leading, required this.title, this.subtitle, this.onClose});
  final Widget leading;
  final String title;
  final String? subtitle;

  /// Shows the close button when set.
  final VoidCallback? onClose;

  /// 40 px rounded tile with [icon]: accent tint, red for [danger], or
  /// tinted with [color].
  static Widget iconTile(BuildContext context, IconData icon, {bool danger = false, Color? color}) {
    final t = context.tokens;
    final tint = danger ? SplaceTokens.danger : color;
    final fg = tint ?? t.accent;
    return Container(
      width: 40,
      height: 40,
      decoration: BoxDecoration(
        color: tint?.withValues(alpha: 0.1) ?? t.accentSoft,
        borderRadius: BorderRadius.circular(10),
        border: Border.all(color: tint?.withValues(alpha: 0.25) ?? t.accentSoftBorder),
      ),
      child: Icon(icon, size: 16, color: fg),
    );
  }

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    return Column(mainAxisSize: MainAxisSize.min, crossAxisAlignment: CrossAxisAlignment.stretch, children: [
      Container(
        color: t.surface,
        padding: const EdgeInsets.fromLTRB(24, 18, 20, 18),
        child: Row(children: [
          leading,
          const SizedBox(width: 14),
          Expanded(
            child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
              Text(title, style: TextStyle(color: t.text, fontSize: 16, fontWeight: FontWeight.w700)),
              if (subtitle != null) ...[
                const SizedBox(height: 2),
                Text(subtitle!, style: TextStyle(color: t.textMuted, fontSize: 13)),
              ],
            ]),
          ),
          if (onClose != null)
            // Bordered square with a dark cross, so it stands out on the
            // white header.
            Tooltip(
              message: context.l10n.close,
              child: Material(
                color: t.surface2,
                shape: RoundedRectangleBorder(borderRadius: BorderRadius.circular(9), side: BorderSide(color: t.border)),
                child: InkWell(
                  borderRadius: BorderRadius.circular(9),
                  hoverColor: SplaceTokens.danger.withValues(alpha: 0.12),
                  onTap: onClose,
                  child: SizedBox(width: 34, height: 34, child: Icon(FA.xmark, size: 15, color: t.text)),
                ),
              ),
            ),
        ]),
      ),
      const AnimatedRainbowBar(height: 3),
    ]);
  }
}

/// The one "Cancel" button of every modal: SPLACE red, same height and
/// text as the filled action next to it (as the HTML app's red Cancel).
class DialogCancelButton extends StatelessWidget {
  const DialogCancelButton({super.key, required this.onPressed, this.label});
  final VoidCallback? onPressed;
  final String? label;

  /// Text and size of the theme's filled button, so both sit level.
  static ButtonStyle style(BuildContext context) {
    final filled = Theme.of(context).filledButtonTheme.style;
    return FilledButton.styleFrom(
      backgroundColor: SplaceRed.c600,
      foregroundColor: Colors.white,
      disabledBackgroundColor: SplaceRed.c600.withValues(alpha: 0.4),
      minimumSize: const Size(112, 44),
      padding: const EdgeInsets.symmetric(horizontal: 20, vertical: 12),
      textStyle: filled?.textStyle?.resolve(const {}),
    );
  }

  @override
  Widget build(BuildContext context) =>
      FilledButton(style: style(context), onPressed: onPressed, child: Text(label ?? context.l10n.cancel));
}

/// Confirmation dialog (clear records, remove record, skip taxonomy…).
Future<bool> confirmDialog(
  BuildContext context, {
  required String title,
  required String body,
  String? body2,
  required String confirm,
  String? cancel,
  bool danger = true,
  IconData icon = FA.triangleExclamation,
}) async {
  final ok = await showDialog<bool>(
    context: context,
    builder: (ctx) => SplaceDialog(
      icon: icon,
      title: title,
      danger: danger,
      size: DialogSize.compact,
      actions: [
        DialogCancelButton(onPressed: () => Navigator.pop(ctx, false), label: cancel),
        FilledButton(
          style: danger ? FilledButton.styleFrom(backgroundColor: SplaceRed.c600) : null,
          onPressed: () => Navigator.pop(ctx, true),
          child: Text(confirm),
        ),
      ],
      child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
        Text(body, style: const TextStyle(fontSize: 14)),
        if (body2 != null) ...[const SizedBox(height: 10), Text(body2, style: TextStyle(fontSize: 13, color: ctx.tokens.textMuted))],
      ]),
    ),
  );
  return ok ?? false;
}

/// Small labelled value used in dialog summaries.
class InfoLine extends StatelessWidget {
  const InfoLine(this.label, this.value, {super.key});
  final String label;
  final String value;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    // Label and value side by side in the same weight; the label is muted.
    return Padding(
      padding: const EdgeInsets.only(bottom: 6),
      child: Row(crossAxisAlignment: CrossAxisAlignment.start, children: [
        SizedBox(width: 120, child: Text(label, style: TextStyle(fontSize: 13, color: t.textMuted))),
        Expanded(child: Text(value, style: TextStyle(fontSize: 13, color: t.text))),
      ]),
    );
  }
}

/// Kicker + child block used inside dialogs.
class DialogBlock extends StatelessWidget {
  const DialogBlock({super.key, required this.label, required this.child});
  final String label;
  final Widget child;

  @override
  Widget build(BuildContext context) => Padding(
        padding: const EdgeInsets.only(bottom: 14),
        child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [Kicker(label), const SizedBox(height: 6), child]),
      );
}
