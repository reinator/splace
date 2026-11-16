import 'dart:math' as math;

import 'package:flutter/material.dart';
import 'package:flutter/rendering.dart';

import '../../services/ncbi.dart' show foldForSearch;
import '../fa.dart';
import '../theme.dart';
import 'splace_tooltip.dart';

/// Flutter versions of the HTML interface's building blocks. Class names in
/// the docs refer to the CSS in docs/theme.css and docs/css/steps/*.css.

/// Chevron used by every select, like the native `<select>` arrow in the
/// HTML app.
const dropdownChevron = Padding(padding: EdgeInsets.only(left: 6), child: Icon(FA.chevronDown, size: 11));

/// `.step-badge`: filled accent circle with the step number.
class StepBadge extends StatelessWidget {
  const StepBadge(this.label, {super.key, this.active = false, this.pending = false, this.size = 28});
  final String label;
  final bool active;
  final bool pending;
  final double size;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    return AnimatedContainer(
      duration: const Duration(milliseconds: 250),
      width: size,
      height: size,
      alignment: Alignment.center,
      decoration: BoxDecoration(
        color: pending ? t.textFaint : t.accent,
        shape: BoxShape.circle,
        boxShadow: active ? [BoxShadow(color: t.accentSoft, spreadRadius: 4)] : null,
      ),
      child: Text(label,
          style: TextStyle(color: Colors.white, fontWeight: FontWeight.w700, fontSize: label.length > 1 ? 11 : 13)),
    );
  }
}

/// `section.step-section`: white rounded card with badge, icon and title.
class StepSection extends StatelessWidget {
  const StepSection({
    super.key,
    this.badge,
    required this.icon,
    required this.title,
    this.description,
    this.count,
    this.actions = const [],
    this.active = false,
    required this.child,
  });

  final String? badge;
  final IconData icon;
  final String title;
  final String? description;
  final int? count;
  final List<Widget> actions;
  final bool active;
  final Widget child;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    return Container(
      margin: EdgeInsets.only(bottom: t.gapSection),
      padding: EdgeInsets.all(t.padCard),
      decoration: BoxDecoration(
        color: t.surface,
        borderRadius: BorderRadius.circular(14),
        border: Border.all(color: t.border),
      ),
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        Wrap(
          alignment: WrapAlignment.spaceBetween,
          crossAxisAlignment: WrapCrossAlignment.center,
          runSpacing: 8,
          spacing: 8,
          children: [
            Row(mainAxisSize: MainAxisSize.min, children: [
              if (badge != null) ...[StepBadge(badge!, active: active), const SizedBox(width: 8)],
              Icon(icon, size: 16, color: t.accent),
              const SizedBox(width: 8),
              Text(title, style: TextStyle(fontSize: 16, fontWeight: FontWeight.w600, color: t.text)),
              if (count != null) ...[const SizedBox(width: 8), CountPill('$count')],
            ]),
            if (actions.isNotEmpty) Wrap(spacing: 8, runSpacing: 8, children: actions),
          ],
        ),
        if (description != null) ...[
          const SizedBox(height: 10),
          Text(description!, style: TextStyle(fontSize: 14, color: t.textMuted)),
        ],
        const SizedBox(height: 14),
        child,
      ]),
    );
  }
}

/// `#panel6a`-style sub panel inside a step ("5a", "5b"...).
class SubPanel extends StatelessWidget {
  const SubPanel({super.key, required this.badge, required this.title, this.body, required this.child, this.optional = false});
  final String badge;
  final String title;
  final String? body;
  final bool optional;
  final Widget child;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    return Container(
      margin: const EdgeInsets.only(bottom: 16),
      padding: EdgeInsets.all(t.padCard),
      decoration: BoxDecoration(color: t.surface, borderRadius: BorderRadius.circular(12), border: Border.all(color: t.border)),
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        Row(children: [
          StepBadge(badge, size: 28),
          const SizedBox(width: 8),
          Expanded(child: Text(title, style: TextStyle(fontSize: 16, fontWeight: FontWeight.w600, color: t.text))),
        ]),
        if (body != null) ...[
          const SizedBox(height: 10),
          Text(body!, style: TextStyle(fontSize: 14, color: t.textMuted)),
        ],
        const SizedBox(height: 12),
        child,
      ]),
    );
  }
}

/// Keys of the lettered parts of the steps ("2a", "7b"…), so the step
/// navigator can scroll to them and follow which one is in view.
class SubAnchors extends InheritedWidget {
  const SubAnchors({super.key, required this.keys, required super.child});
  final Map<String, GlobalKey> keys;

  GlobalKey keyFor(String id) => keys.putIfAbsent(id, GlobalKey.new);

  static SubAnchors? maybeOf(BuildContext context) => context.getInheritedWidgetOfExactType<SubAnchors>();

  @override
  bool updateShouldNotify(SubAnchors old) => !identical(old.keys, keys);
}

/// Marks [child] as the part [id] of a step (see [SubAnchors]).
class SubAnchor extends StatelessWidget {
  const SubAnchor({super.key, required this.id, required this.child});
  final String id;
  final Widget child;

  @override
  Widget build(BuildContext context) {
    final anchors = SubAnchors.maybeOf(context);
    return anchors == null ? child : KeyedSubtree(key: anchors.keyFor(id), child: child);
  }
}

/// Title of a lettered part of a step ("7a Which alignments"), with what
/// it decides below.
class SubHeading extends StatelessWidget {
  const SubHeading({super.key, required this.id, required this.title, this.body, this.trailing});
  final String id, title;
  final String? body;
  final Widget? trailing;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    return Padding(
      padding: const EdgeInsets.only(bottom: 10),
      child: Row(crossAxisAlignment: CrossAxisAlignment.start, children: [
        StepBadge(id, size: 26),
        const SizedBox(width: 10),
        Expanded(
          child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
            Padding(
              padding: const EdgeInsets.only(top: 3),
              child: Text(title, style: TextStyle(fontSize: 15.5, fontWeight: FontWeight.w700, color: t.text)),
            ),
            if (body != null) ...[
              const SizedBox(height: 2),
              Text(body!, style: TextStyle(fontSize: 13, color: t.textMuted, height: 1.4)),
            ],
          ]),
        ),
        if (trailing != null) ...[const SizedBox(width: 8), trailing!],
      ]),
    );
  }
}

/// `.step1-card`: option card with kicker ("OPTION 1"), title and info.
class OptionCard extends StatelessWidget {
  const OptionCard({
    super.key,
    this.badge,
    this.kicker,
    required this.title,
    this.info,
    this.text,
    this.action,
    this.footer,
    required this.child,
  });

  /// Letter of the part of the step ("1a"), before the title.
  final String? badge;
  final String? kicker;
  final String title;
  final String? info;
  final String? text;
  final Widget? action;
  final Widget child;

  /// Kept at the bottom of the card (run buttons), also when the card is
  /// stretched to the height of its neighbour ([EqualHeightRow]).
  final Widget? footer;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    return Container(
      padding: const EdgeInsets.all(16),
      decoration: BoxDecoration(
        color: t.surface,
        borderRadius: BorderRadius.circular(16),
        border: Border.all(color: t.border),
        boxShadow: shadowSm,
      ),
      child: Column(
        crossAxisAlignment: CrossAxisAlignment.stretch,
        mainAxisAlignment: MainAxisAlignment.spaceBetween,
        children: [
          Column(crossAxisAlignment: CrossAxisAlignment.stretch, mainAxisSize: MainAxisSize.min, children: [
            Row(crossAxisAlignment: CrossAxisAlignment.center, children: [
              if (badge != null) ...[StepBadge(badge!, size: 26), const SizedBox(width: 10)],
              Expanded(
                child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
                  if (kicker != null) Kicker(kicker!),
                  Row(children: [
                    Flexible(child: Text(title, style: TextStyle(fontSize: 16, fontWeight: FontWeight.w600, color: t.text))),
                    if (info != null) ...[const SizedBox(width: 7), InfoDot(info!)],
                  ]),
                ]),
              ),
              ?action,
            ]),
            if (text != null) ...[
              const SizedBox(height: 8),
              Text(text!, style: TextStyle(fontSize: 14, color: t.textMuted)),
            ],
            const SizedBox(height: 14),
            child,
          ]),
          ?footer,
        ],
      ),
    );
  }
}

/// Uppercase small label (`.step1-card-kicker`, `.mm-title`, `QUERY PREVIEW`).
class Kicker extends StatelessWidget {
  const Kicker(this.text, {super.key, this.color, this.size = 11, this.singleLine = false});
  final String text;
  final Color? color;
  final double size;
  final bool singleLine;

  @override
  Widget build(BuildContext context) => Text(
        text.toUpperCase(),
        maxLines: singleLine ? 1 : null,
        overflow: singleLine ? TextOverflow.ellipsis : null,
        style: TextStyle(
          fontSize: size,
          fontWeight: FontWeight.w700,
          letterSpacing: size * 0.08,
          color: color ?? context.tokens.textMuted,
        ),
      );
}

/// `.step1-help-btn`: small accent info icon with a tooltip.
class InfoDot extends StatelessWidget {
  const InfoDot(this.message, {super.key});
  final String message;

  @override
  Widget build(BuildContext context) => SplaceTooltip(
        message: message,
        child: Icon(FA.circleInfo, size: 15, color: context.tokens.accent),
      );
}

/// `.step1-tips-accordion`.
class TipsAccordion extends StatefulWidget {
  const TipsAccordion({super.key, required this.kicker, required this.title, required this.children, this.icon = FA.lightbulb});
  final String kicker;
  final String title;
  final IconData icon;
  final List<Widget> children;

  @override
  State<TipsAccordion> createState() => _TipsAccordionState();
}

class _TipsAccordionState extends State<TipsAccordion> {
  bool _open = false;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final dark = Theme.of(context).brightness == Brightness.dark;
    return Container(
      margin: const EdgeInsets.only(bottom: 20),
      decoration: BoxDecoration(
        color: t.accentSoft,
        borderRadius: BorderRadius.circular(14),
        border: Border.all(color: dark ? t.accentSoftBorder : SplaceBlue.c200),
      ),
      clipBehavior: Clip.antiAlias,
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        InkWell(
          onTap: () => setState(() => _open = !_open),
          hoverColor: dark ? t.accentSoft : SplaceBlue.c100,
          child: Padding(
            padding: const EdgeInsets.symmetric(horizontal: 18, vertical: 16),
            child: Row(children: [
              // `.step1-tips-icon`: blue-500 → blue-600 gradient tile.
              Container(
                width: 36,
                height: 36,
                decoration: BoxDecoration(
                  gradient: const LinearGradient(begin: Alignment.topLeft, end: Alignment.bottomRight, colors: [SplaceBlue.c500, SplaceBlue.c600]),
                  borderRadius: BorderRadius.circular(11),
                ),
                child: Icon(widget.icon, size: 17, color: Colors.white),
              ),
              const SizedBox(width: 12),
              Expanded(
                child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
                  Kicker(widget.kicker, color: dark ? SplaceBlue.c300 : t.accent),
                  const SizedBox(height: 2),
                  Text(widget.title,
                      style: TextStyle(fontSize: 15, fontWeight: FontWeight.w700, color: dark ? SplaceBlue.c100 : SplaceBlue.c900)),
                ]),
              ),
              AnimatedRotation(
                turns: _open ? 0.5 : 0,
                duration: const Duration(milliseconds: 200),
                child: Icon(FA.chevronDown, size: 14, color: dark ? SplaceBlue.c300 : t.accent),
              ),
            ]),
          ),
        ),
        if (_open)
          Padding(
            padding: const EdgeInsets.fromLTRB(18, 0, 18, 18),
            child: DefaultTextStyle.merge(
              style: TextStyle(fontSize: 14, color: t.accentText, height: 1.6),
              child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: widget.children),
            ),
          ),
      ]),
    );
  }
}

/// Bullet line for tips lists. With [icon], the marker is that icon in
/// [iconColor] (`.step1-tips-list li::before`: red warning triangle).
class Bullet extends StatelessWidget {
  const Bullet(this.child, {super.key, this.icon, this.iconColor = SplaceRed.c600, this.gap = 7});
  final Widget child;
  final IconData? icon;
  final Color iconColor;
  final double gap;

  @override
  Widget build(BuildContext context) => Padding(
        padding: EdgeInsets.only(top: gap),
        child: Row(crossAxisAlignment: CrossAxisAlignment.start, children: [
          if (icon == null)
            const Padding(padding: EdgeInsets.only(left: 4, right: 10), child: Text('•'))
          else
            Padding(
              padding: const EdgeInsets.only(top: 3, right: 10),
              child: SizedBox(width: 22, child: Icon(icon, size: 15, color: iconColor)),
            ),
          Expanded(child: child),
        ]),
      );
}

/// Text with `**bold**` keywords, for translated strings that highlight
/// their key terms (tips, notes). Bold runs take [strongColor], or the
/// main text colour.
class MarkupText extends StatelessWidget {
  const MarkupText(this.text, {super.key, this.style, this.strongColor});
  final String text;
  final TextStyle? style;
  final Color? strongColor;

  /// The spans of [text]: odd `**`-delimited runs are bold.
  static List<TextSpan> spans(String text, {TextStyle? strong}) {
    final parts = text.split('**');
    return [
      for (final (i, part) in parts.indexed)
        if (part.isNotEmpty) TextSpan(text: part, style: i.isOdd ? strong : null),
    ];
  }

  @override
  Widget build(BuildContext context) {
    final dark = Theme.of(context).brightness == Brightness.dark;
    final strong = TextStyle(
      fontWeight: FontWeight.w700,
      color: strongColor ?? (dark ? context.tokens.text : SplaceBlue.c900),
    );
    return Text.rich(TextSpan(style: style, children: spans(text, strong: strong)));
  }
}

/// `.step1-primary-btn` (accent, up to 15rem wide).
class PrimaryButton extends StatelessWidget {
  const PrimaryButton({super.key, required this.label, this.icon, this.onPressed, this.expand = false});
  final String label;
  final IconData? icon;
  final VoidCallback? onPressed;
  final bool expand;

  @override
  Widget build(BuildContext context) {
    final button = FilledButton(
      onPressed: onPressed,
      // Buttons show text only; [icon] is kept for call sites but not drawn.
      child: Text(label, overflow: TextOverflow.ellipsis),
    );
    if (expand) return SizedBox(width: double.infinity, child: button);
    return ConstrainedBox(constraints: const BoxConstraints(minWidth: 240), child: button);
  }
}

enum ChipButtonStyle {
  /// `.step1-example-btn` (red).
  danger,

  /// `.step1-inline-btn` / "Select All" (grey).
  neutral,

  /// "Shared Only" (blue-100).
  soft,

  /// `.step2-toolbar-btn` (accent).
  accent,
}

class ChipButton extends StatelessWidget {
  const ChipButton({super.key, required this.label, this.icon, this.onPressed, this.style = ChipButtonStyle.neutral, this.tooltip, this.keepIcon = false});
  final String label;

  /// Drawn only with [keepIcon] (table column toggles); other buttons show
  /// text only.
  final IconData? icon;
  final bool keepIcon;
  final VoidCallback? onPressed;
  final ChipButtonStyle style;
  final String? tooltip;

  /// Height of every small button.
  static const height = 32.0;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final (bg, fg, border) = switch (style) {
      ChipButtonStyle.danger => (SplaceRed.c600, Colors.white, SplaceRed.c600),
      ChipButtonStyle.neutral => (t.surface2, t.textMuted, t.border),
      ChipButtonStyle.soft => (t.accentSoftBorder, t.accentText, t.accentSoftBorder),
      ChipButtonStyle.accent => (t.accent, Colors.white, t.accent),
    };
    final enabled = onPressed != null;
    Widget w = Material(
      color: enabled ? bg : bg.withValues(alpha: 0.5),
      shape: RoundedRectangleBorder(borderRadius: BorderRadius.circular(8), side: BorderSide(color: border)),
      child: InkWell(
        borderRadius: BorderRadius.circular(8),
        onTap: onPressed,
        // Small buttons are all [height] tall (icon buttons beside them too).
        child: ConstrainedBox(
          constraints: const BoxConstraints(minHeight: height),
          child: Padding(
            padding: const EdgeInsets.symmetric(horizontal: 13, vertical: 6),
            child: Row(mainAxisSize: MainAxisSize.min, children: [
              if (icon != null && keepIcon) ...[Icon(icon, size: 13, color: fg), const SizedBox(width: 7)],
              Text(label, style: TextStyle(color: fg, fontSize: 13, height: 1.2, fontWeight: FontWeight.w600)),
            ]),
          ),
        ),
      ),
    );
    if (tooltip != null) w = SplaceTooltip(message: tooltip!, child: w);
    return w;
  }
}

/// `.icon-btn` of the header and log toolbar.
class IconBtn extends StatefulWidget {
  const IconBtn({super.key, required this.icon, required this.tooltip, this.onPressed, this.active = false, this.size = 36, this.bordered = false});
  final IconData icon;
  final String tooltip;
  final VoidCallback? onPressed;
  final bool active;
  final double size;

  /// `.logs-toolbar button` look (always bordered).
  final bool bordered;

  @override
  State<IconBtn> createState() => _IconBtnState();
}

class _IconBtnState extends State<IconBtn> {
  bool _hover = false;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final active = widget.active;
    final bg = active ? t.accentSoft : (_hover ? t.surface2 : (widget.bordered ? t.surface : Colors.transparent));
    final fg = active || _hover ? t.accent : t.textMuted;
    final border = active
        ? t.accentSoftBorder
        : (_hover ? (widget.bordered ? t.accentSoftBorder : t.border) : (widget.bordered ? t.border : Colors.transparent));
    return SplaceTooltip(
      message: widget.tooltip,
      child: MouseRegion(
        cursor: SystemMouseCursors.click,
        onEnter: (_) => setState(() => _hover = true),
        onExit: (_) => setState(() => _hover = false),
        child: GestureDetector(
          onTap: widget.onPressed,
          child: AnimatedContainer(
            duration: const Duration(milliseconds: 150),
            width: widget.size,
            height: widget.size,
            decoration: BoxDecoration(
              color: bg,
              borderRadius: BorderRadius.circular(widget.bordered ? 6 : 10),
              border: Border.all(color: border),
            ),
            child: Icon(widget.icon, size: widget.bordered ? 12 : 17, color: fg),
          ),
        ),
      ),
    );
  }
}

/// `.lang-toggle`.
class LangToggle extends StatelessWidget {
  const LangToggle({super.key, required this.current, required this.onChanged});
  final String current;
  final ValueChanged<String> onChanged;

  static const _langs = [('en', 'us', 'EN'), ('es', 'es', 'ES'), ('pt', 'br', 'PT')];

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    return Container(
      padding: const EdgeInsets.all(2),
      decoration: BoxDecoration(color: t.surface2, borderRadius: BorderRadius.circular(10), border: Border.all(color: t.border)),
      child: Row(mainAxisSize: MainAxisSize.min, children: [
        for (final (code, flag, label) in _langs)
          GestureDetector(
            onTap: () => onChanged(code),
            child: MouseRegion(
              cursor: SystemMouseCursors.click,
              child: AnimatedContainer(
                duration: const Duration(milliseconds: 150),
                padding: const EdgeInsets.symmetric(horizontal: 10, vertical: 4),
                decoration: BoxDecoration(
                  color: current == code ? t.surface : Colors.transparent,
                  borderRadius: BorderRadius.circular(6),
                  boxShadow: current == code ? shadowSm : null,
                ),
                child: Row(mainAxisSize: MainAxisSize.min, children: [
                  ClipRRect(
                    borderRadius: BorderRadius.circular(2),
                    child: Image.asset('assets/images/flag_$flag.png', width: 16),
                  ),
                  const SizedBox(width: 6),
                  Text(label,
                      style: TextStyle(
                        fontSize: 12,
                        fontWeight: FontWeight.w600,
                        color: current == code ? t.accent : t.textMuted,
                      )),
                ]),
              ),
            ),
          ),
      ]),
    );
  }
}

/// `.badge` count next to a step title.
class CountPill extends StatelessWidget {
  const CountPill(this.text, {super.key});
  final String text;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    return Container(
      padding: const EdgeInsets.symmetric(horizontal: 8, vertical: 2),
      decoration: BoxDecoration(color: t.accentSoftBorder, borderRadius: BorderRadius.circular(8)),
      child: Text(text, style: TextStyle(fontSize: 12, fontWeight: FontWeight.w700, color: t.accentText)),
    );
  }
}

/// Pipeline/method badge (`Codon-aware`, `MAFFT --amino`…), drawn like the
/// soft buttons of the step 4 card ([ChipButtonStyle.soft]): rounded
/// rectangle, soft blue fill, bold blue text. A [color] tints it instead
/// (warnings, frames).
class TagChip extends StatelessWidget {
  const TagChip(this.text, {super.key, this.icon, this.color, this.leading});
  final String text;
  final IconData? icon;
  final Color? color;

  /// Drawn before the icon and text (the step number of a pipeline).
  final Widget? leading;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final c = color;
    final (bg, fg, border) = c == null ? (t.accentSoftBorder, t.accentText, t.accentSoftBorder) : (c.withValues(alpha: 0.12), c, c.withValues(alpha: 0.32));
    return Container(
      padding: const EdgeInsets.symmetric(horizontal: 10, vertical: 4),
      decoration: BoxDecoration(color: bg, borderRadius: BorderRadius.circular(8), border: Border.all(color: border)),
      child: Row(mainAxisSize: MainAxisSize.min, children: [
        if (leading != null) ...[leading!, const SizedBox(width: 6)],
        if (icon != null) ...[Icon(icon, size: 11, color: fg), const SizedBox(width: 6)],
        Flexible(
          child: Text(text, maxLines: 1, overflow: TextOverflow.ellipsis, style: TextStyle(fontSize: 12, fontWeight: FontWeight.w600, color: fg)),
        ),
      ]),
    );
  }
}

/// Labelled form field (`.step1-field`).
class Field extends StatelessWidget {
  const Field({super.key, required this.label, required this.child, this.info});
  final String label;
  final String? info;
  final Widget child;

  @override
  Widget build(BuildContext context) => Column(crossAxisAlignment: CrossAxisAlignment.stretch, mainAxisSize: MainAxisSize.min, children: [
        Row(children: [
          Flexible(child: Text(label, style: TextStyle(fontSize: 14, fontWeight: FontWeight.w600, color: context.tokens.text))),
          if (info != null) ...[const SizedBox(width: 6), InfoDot(info!)],
        ]),
        const SizedBox(height: 6),
        child,
      ]);
}

/// Responsive grid like `grid-template-columns: repeat(auto-fit, minmax(w, 1fr))`.
class AutoGrid extends StatelessWidget {
  const AutoGrid({super.key, required this.children, this.minItemWidth = 180, this.spacing = 12});
  final List<Widget> children;
  final double minItemWidth;
  final double spacing;

  @override
  Widget build(BuildContext context) => LayoutBuilder(builder: (context, c) {
        // Never more columns than children, so a short row fills the width
        // instead of leaving empty slots on the right.
        final cols = math.max(1, math.min(children.length, ((c.maxWidth + spacing) / (minItemWidth + spacing)).floor()));
        final w = (c.maxWidth - spacing * (cols - 1)) / cols;
        return Wrap(spacing: spacing, runSpacing: spacing, children: [
          for (final child in children) SizedBox(width: w, child: child),
        ]);
      });
}

/// Toggle row inside a bordered box ("Complete Only" switch in step 1). The
/// whole box toggles and takes the accent tint when on.
class SwitchBox extends StatelessWidget {
  const SwitchBox({super.key, required this.label, required this.value, required this.onChanged, this.detail, this.mono = false});
  final String label;
  final bool value;
  final ValueChanged<bool>? onChanged;

  /// Second line (what the option does), muted.
  final String? detail;

  /// Command-line flag as the label: monospace, never wrapped.
  final bool mono;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final on = value && onChanged != null;
    return GestureDetector(
      onTap: onChanged == null ? null : () => onChanged!(!value),
      child: MouseRegion(
        cursor: onChanged == null ? SystemMouseCursors.basic : SystemMouseCursors.click,
        child: AnimatedContainer(
          duration: const Duration(milliseconds: 150),
          // At least the 44 px control height; long labels wrap.
          constraints: const BoxConstraints(minHeight: 44),
          alignment: Alignment.centerLeft,
          padding: const EdgeInsets.symmetric(horizontal: 10, vertical: 6),
          decoration: BoxDecoration(
            borderRadius: BorderRadius.circular(8),
            border: Border.all(color: on ? t.accentSoftBorder : t.border),
            color: on ? t.accentSoft : t.surface,
          ),
          child: Row(children: [
            // M3 switch at 26 px (≈ 42×26), inside the 44 px control height.
            SizedBox(height: 26, child: FittedBox(child: Switch(value: value, onChanged: onChanged))),
            const SizedBox(width: 8),
            Expanded(
              child: Column(crossAxisAlignment: CrossAxisAlignment.start, mainAxisSize: MainAxisSize.min, children: [
                Text(
                  label,
                  maxLines: mono ? 1 : null,
                  softWrap: !mono,
                  overflow: mono ? TextOverflow.fade : null,
                  style: TextStyle(
                    fontSize: mono ? 13 : 14,
                    fontFamily: mono ? monoFont : null,
                    fontWeight: FontWeight.w600,
                    color: on ? t.accentText : t.text,
                  ),
                ),
                if (detail != null) ...[
                  const SizedBox(height: 2),
                  Text(detail!, style: TextStyle(fontSize: 12.5, color: t.textMuted, height: 1.35)),
                ],
              ]),
            ),
          ]),
        ),
      ),
    );
  }
}

/// One option of a [SplaceSelect].
class SelectOption<T> {
  const SelectOption(this.value, this.label, {this.icon, this.detail});
  final T value;
  final String label;
  final IconData? icon;
  final String? detail;
}

/// Material 3 select (`DropdownMenu`): read-only field showing the value in
/// body text, and a rounded menu with icons, details and a check on the
/// current value.
class SplaceSelect<T> extends StatelessWidget {
  const SplaceSelect({super.key, required this.value, required this.options, required this.onChanged});
  final T value;
  final List<SelectOption<T>> options;
  final ValueChanged<T>? onChanged;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    return DropdownMenu<T>(
      // initialSelection is read once; a new key follows outside changes.
      key: ValueKey(value),
      initialSelection: value,
      enabled: onChanged != null,
      expandedInsets: EdgeInsets.zero,
      requestFocusOnTap: false,
      enableFilter: false,
      enableSearch: false,
      // Value in body text, like the text fields next to it.
      textStyle: TextStyle(fontSize: 14, color: t.text),
      // A plain chevron instead of the default 48 px IconButton, so the
      // field keeps the 44 px height and the icon stays centred.
      decorationBuilder: (context, menu) => InputDecoration(
        suffixIcon: Icon(menu.isOpen ? FA.chevronUp : FA.chevronDown, size: 11, color: menu.isOpen ? t.accent : t.textMuted),
        suffixIconConstraints: const BoxConstraints(minWidth: 36),
      ),
      menuStyle: MenuStyle(
        backgroundColor: WidgetStatePropertyAll(t.surface),
        surfaceTintColor: const WidgetStatePropertyAll(Colors.transparent),
        elevation: const WidgetStatePropertyAll(8),
        shadowColor: const WidgetStatePropertyAll(Color(0x40000000)),
        padding: const WidgetStatePropertyAll(EdgeInsets.all(6)),
        shape: WidgetStatePropertyAll(RoundedRectangleBorder(borderRadius: BorderRadius.circular(12), side: BorderSide(color: t.border))),
      ),
      onSelected: (v) {
        if (v != null) onChanged?.call(v);
      },
      dropdownMenuEntries: [
        for (final o in options)
          DropdownMenuEntry<T>(
            value: o.value,
            label: o.label,
            labelWidget: Column(crossAxisAlignment: CrossAxisAlignment.start, mainAxisSize: MainAxisSize.min, children: [
              Text(o.label, style: TextStyle(fontSize: 14, color: o.value == value ? t.accentText : t.text)),
              if (o.detail != null) Text(o.detail!, style: TextStyle(fontSize: 12, color: t.textMuted)),
            ]),
            leadingIcon: o.icon == null ? null : Icon(o.icon, size: 14, color: o.value == value ? t.accent : t.textMuted),
            trailingIcon: o.value == value ? Icon(FA.check, size: 12, color: t.accent) : null,
            style: ButtonStyle(
              padding: const WidgetStatePropertyAll(EdgeInsets.symmetric(horizontal: 12, vertical: 10)),
              shape: WidgetStatePropertyAll(RoundedRectangleBorder(borderRadius: BorderRadius.circular(8))),
              backgroundColor: WidgetStateProperty.resolveWith(
                (s) => o.value == value ? t.accentSoft : (s.contains(WidgetState.hovered) ? t.surface2 : null),
              ),
            ),
          ),
      ],
    );
  }
}

/// Radio card (`border-2 border-splace-blue-500 bg-splace-blue-50` when selected).
class RadioCard<T> extends StatelessWidget {
  const RadioCard({super.key, required this.value, required this.group, required this.onChanged, required this.title, this.subtitle, this.trailing});
  final T value;
  final T group;
  final ValueChanged<T> onChanged;
  final String title;
  final String? subtitle;
  final Widget? trailing;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final selected = value == group;
    return GestureDetector(
      onTap: () => onChanged(value),
      child: MouseRegion(
        cursor: SystemMouseCursors.click,
        child: Container(
          padding: const EdgeInsets.all(12),
          decoration: BoxDecoration(
            color: selected ? t.accentSoft : t.surface,
            borderRadius: BorderRadius.circular(12),
            border: Border.all(color: selected ? SplaceBlue.c500 : t.border, width: 2),
          ),
          child: Row(crossAxisAlignment: CrossAxisAlignment.start, children: [
            Padding(
              padding: const EdgeInsets.only(top: 2),
              child: Container(
                width: 16,
                height: 16,
                decoration: BoxDecoration(
                  shape: BoxShape.circle,
                  border: Border.all(color: selected ? t.accent : t.textFaint, width: selected ? 5 : 1.5),
                ),
              ),
            ),
            const SizedBox(width: 12),
            Expanded(
              child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
                Text(title, style: TextStyle(fontSize: 14, fontWeight: FontWeight.w600, color: t.text)),
                if (subtitle != null) ...[
                  const SizedBox(height: 2),
                  Text(subtitle!, style: TextStyle(fontSize: 13, color: t.textMuted)),
                ],
              ]),
            ),
            ?trailing,
          ]),
        ),
      ),
    );
  }
}

/// Checkbox card (`rounded-xl border bg-gray-50`).
class CheckCard extends StatelessWidget {
  const CheckCard({super.key, required this.value, required this.onChanged, required this.title, this.subtitle});
  final bool value;
  final ValueChanged<bool> onChanged;
  final String title;
  final String? subtitle;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    return GestureDetector(
      onTap: () => onChanged(!value),
      child: MouseRegion(
        cursor: SystemMouseCursors.click,
        child: Container(
          padding: const EdgeInsets.all(12),
          decoration: BoxDecoration(color: t.surface2, borderRadius: BorderRadius.circular(12), border: Border.all(color: t.border)),
          child: Row(children: [
            SizedBox(width: 22, height: 22, child: Checkbox(value: value, onChanged: (v) => onChanged(v!))),
            const SizedBox(width: 12),
            Expanded(
              child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
                Text(title, style: TextStyle(fontSize: 14, fontWeight: FontWeight.w600, color: t.text)),
                if (subtitle != null) Text(subtitle!, style: TextStyle(fontSize: 13, color: t.textMuted)),
              ]),
            ),
          ]),
        ),
      ),
    );
  }
}

enum PanelTone { info, warn, duplicate, success, danger }

/// Coloured panels: blue info box (step 4 legend), amber "missing genes",
/// orange "duplicate genes", green/red results.
class TonePanel extends StatelessWidget {
  const TonePanel({super.key, required this.tone, this.title, this.icon, required this.child});
  final PanelTone tone;
  final String? title;
  final IconData? icon;
  final Widget child;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final dark = Theme.of(context).brightness == Brightness.dark;
    final (bg, border, fg) = switch (tone) {
      PanelTone.info => (t.accentSoft, t.accentSoftBorder, t.accentText),
      // Warnings and duplicates in the yellow/orange of the app scale.
      PanelTone.warn => dark
          ? (const Color(0x1FFEE090), const Color(0x4DFEE090), SplacePalette.yellow)
          : (const Color(0x40FEE090), SplacePalette.yellow, SplacePalette.red900),
      PanelTone.duplicate => dark
          ? (const Color(0x1FFDAE61), const Color(0x4DFDAE61), SplacePalette.orange)
          : (const Color(0x33FDAE61), SplacePalette.orange, SplacePalette.red900),
      // Success panels share the app blue with info panels.
      PanelTone.success => (t.accentSoft, t.accentSoftBorder, t.accentText),
      PanelTone.danger => dark
          ? (const Color(0x1FD73027), const Color(0x4DD73027), SplacePalette.red500)
          : (SplaceRed.c50, const Color(0x59D73027), SplacePalette.red900),
    };
    return Container(
      padding: const EdgeInsets.all(14),
      decoration: BoxDecoration(color: bg, borderRadius: BorderRadius.circular(12), border: Border.all(color: border)),
      child: DefaultTextStyle.merge(
        style: TextStyle(color: fg, fontSize: 14),
        child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
          if (title != null) ...[
            Row(children: [
              if (icon != null) ...[Icon(icon, size: 13, color: fg), const SizedBox(width: 6)],
              Expanded(child: Kicker(title!, color: fg, size: 12)),
            ]),
            const SizedBox(height: 10),
          ],
          child,
        ]),
      ),
    );
  }
}

/// `.gene-chip` with switch, name, count and mean length.
class MarkerChip extends StatelessWidget {
  const MarkerChip({super.key, required this.name, required this.active, required this.onChanged, this.count, this.meta, this.enabled = true});
  final String name;
  final bool active;
  final bool enabled;
  final ValueChanged<bool> onChanged;
  final String? count;
  final String? meta;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    return Opacity(
      opacity: enabled ? 1 : 0.5,
      child: GestureDetector(
        onTap: enabled ? () => onChanged(!active) : null,
        child: MouseRegion(
          cursor: enabled ? SystemMouseCursors.click : SystemMouseCursors.basic,
          child: AnimatedContainer(
            duration: const Duration(milliseconds: 150),
            padding: const EdgeInsets.symmetric(horizontal: 10, vertical: 4),
            decoration: BoxDecoration(
              color: active ? t.accentSoft : t.surface,
              borderRadius: BorderRadius.circular(8),
              border: Border.all(color: active ? t.accentSoftBorder : t.border),
            ),
            child: Row(children: [
              SizedBox(
                height: 28,
                child: FittedBox(child: Switch(value: active, onChanged: enabled ? onChanged : null)),
              ),
              const SizedBox(width: 6),
              Flexible(
                child: Text(name,
                    overflow: TextOverflow.ellipsis,
                    style: TextStyle(fontSize: 14, fontWeight: FontWeight.w500, color: active ? t.accentText : t.text)),
              ),
              if (count != null) ...[
                const SizedBox(width: 6),
                Text(count!, style: TextStyle(fontSize: 12, color: t.textMuted)),
              ],
              if (meta != null) ...[
                const SizedBox(width: 6),
                Flexible(
                  child: Text(meta!,
                      overflow: TextOverflow.ellipsis,
                      style: TextStyle(fontSize: 12, fontStyle: FontStyle.italic, color: t.textFaint)),
                ),
              ],
            ]),
          ),
        ),
      ),
    );
  }
}

/// `#dropZone`: dashed area to drop or browse files.
class DropZone extends StatelessWidget {
  const DropZone({super.key, required this.title, required this.subtitle, required this.onTap, this.highlight = false});
  final String title;
  final String subtitle;
  final VoidCallback onTap;
  final bool highlight;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    return GestureDetector(
      onTap: onTap,
      child: MouseRegion(
        cursor: SystemMouseCursors.click,
        child: CustomPaint(
          painter: _DashedBorder(color: highlight ? t.accent : const Color(0xFFD1D5DB), radius: 12),
          child: AnimatedContainer(
            duration: const Duration(milliseconds: 150),
            padding: const EdgeInsets.all(24),
            decoration: BoxDecoration(
              color: highlight ? t.accentSoft : Colors.transparent,
              borderRadius: BorderRadius.circular(12),
            ),
            child: Column(children: [
              Icon(FA.cloudArrowUp, size: 30, color: highlight ? t.accent : t.textFaint),
              const SizedBox(height: 8),
              Text(title, textAlign: TextAlign.center, style: TextStyle(fontWeight: FontWeight.w500, color: t.textMuted, fontSize: 14)),
              const SizedBox(height: 4),
              Text(subtitle, textAlign: TextAlign.center, style: TextStyle(fontSize: 14, color: t.textFaint)),
            ]),
          ),
        ),
      ),
    );
  }
}

class _DashedBorder extends CustomPainter {
  _DashedBorder({required this.color, required this.radius});
  final Color color;
  final double radius;

  @override
  void paint(Canvas canvas, Size size) {
    final paint = Paint()
      ..color = color
      ..style = PaintingStyle.stroke
      ..strokeWidth = 2;
    final path = Path()..addRRect(RRect.fromRectAndRadius(Offset.zero & size, Radius.circular(radius)));
    for (final metric in path.computeMetrics()) {
      for (var d = 0.0; d < metric.length; d += 10) {
        canvas.drawPath(metric.extractPath(d, d + 6), paint);
      }
    }
  }

  @override
  bool shouldRepaint(_DashedBorder old) => old.color != color;
}

/// `.step1-query-preview`: grey box with an uppercase label and mono text.
class PreviewBox extends StatelessWidget {
  const PreviewBox({super.key, required this.label, required this.text, this.maxHeight = 220});
  final String label;
  final String text;
  final double maxHeight;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    return Container(
      padding: const EdgeInsets.symmetric(horizontal: 14, vertical: 12),
      decoration: BoxDecoration(color: t.surface2, borderRadius: BorderRadius.circular(14), border: Border.all(color: t.border)),
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        Kicker(label),
        const SizedBox(height: 6),
        ConstrainedBox(
          constraints: BoxConstraints(maxHeight: maxHeight),
          child: SingleChildScrollView(
            child: SelectableText(text, style: TextStyle(fontFamily: monoFont, fontSize: 12, color: t.text)),
          ),
        ),
      ]),
    );
  }
}

/// Header + rows table styled like `.records-table` (uppercase header,
/// zebra rows) inside a bordered, height-limited, virtualized box. With
/// [resizable], each column can be widened or narrowed by dragging the
/// edge of its header (double-click restores it), and the table scrolls
/// sideways when the columns are wider than the box.
class RecordsTable extends StatefulWidget {
  const RecordsTable({
    super.key,
    required this.columns,
    required this.rowCount,
    required this.cellBuilder,
    this.flex,
    this.height = 420,
    this.rowHeight = 34,
    this.fill = false,
    this.headerCells,
    this.rowColor,
    this.onRowTap,
    this.cellColor,
    this.resizable = false,
    this.minWidths,
  });

  /// Highlight for some rows (selected genomes…).
  final Color? Function(int row)? rowColor;

  /// Background of single cells (gene counts coloured against the mode).
  final Color? Function(int row, int col)? cellColor;
  final void Function(int row)? onRowTap;

  /// Expands to the available height (dialogs) instead of [height].
  final bool fill;

  /// Optional custom header widgets (e.g. a select-all checkbox).
  final List<Widget?>? headerCells;
  final List<String> columns;
  final List<int>? flex;
  final int rowCount;
  final Widget Function(BuildContext, int row, int col) cellBuilder;
  final double height;
  final double rowHeight;

  /// Columns resized by dragging their header edge.
  final bool resizable;

  /// Starting width of each column at least (resizable tables).
  final List<double>? minWidths;

  static const minColumnWidth = 48.0;

  @override
  State<RecordsTable> createState() => _RecordsTableState();
}

class _RecordsTableState extends State<RecordsTable> {
  /// Widths set by the user, by column name (they survive hiding others).
  final Map<String, double> _widths = {};
  final _h = ScrollController();

  @override
  void dispose() {
    _h.dispose();
    super.dispose();
  }

  RecordsTable get w => widget;

  Widget _header(BuildContext context, int c) {
    final t = context.tokens;
    // Column names as written (PCGs, rRNAs, tRNAs read better than in
    // capitals).
    return w.headerCells?[c] ??
        SplaceTooltip(
          message: w.columns[c],
          child: Text(w.columns[c], maxLines: 1, overflow: TextOverflow.ellipsis, style: TextStyle(fontSize: 12.5, fontWeight: FontWeight.w700, color: t.textMuted)),
        );
  }

  Widget _cell(BuildContext context, int r, int c) {
    final t = context.tokens;
    final color = w.cellColor?.call(r, c);
    return Container(
      color: color,
      margin: color == null ? null : const EdgeInsets.only(right: 4),
      padding: color == null ? null : const EdgeInsets.symmetric(horizontal: 8),
      // Centred vertically, full cell width (right-aligned actions keep
      // their place).
      child: Column(mainAxisAlignment: MainAxisAlignment.center, crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        DefaultTextStyle.merge(
          style: TextStyle(fontSize: 13, color: t.text),
          maxLines: 1,
          overflow: TextOverflow.ellipsis,
          child: w.cellBuilder(context, r, c),
        ),
      ]),
    );
  }

  Widget _row(BuildContext context, int r, Widget Function(int c, Widget child) sized) {
    final t = context.tokens;
    return _tappable(
      r,
      Container(
        padding: const EdgeInsets.symmetric(horizontal: 12),
        decoration: BoxDecoration(
          color: w.rowColor?.call(r) ?? (r.isEven ? t.surface : t.surface2.withValues(alpha: 0.5)),
          border: Border(bottom: BorderSide(color: t.borderSoft)),
        ),
        child: Row(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
          for (var c = 0; c < w.columns.length; c++) sized(c, _cell(context, r, c)),
        ]),
      ),
    );
  }

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final f = w.flex ?? List.filled(w.columns.length, 1);
    final box = BoxDecoration(borderRadius: BorderRadius.circular(12), border: Border.all(color: t.border));
    if (!w.resizable) {
      Widget sized(int c, Widget child) => Expanded(flex: f[c], child: child);
      return Container(
        decoration: box,
        clipBehavior: Clip.antiAlias,
        child: Column(children: [
          Container(
            color: t.surface2,
            padding: const EdgeInsets.symmetric(horizontal: 12, vertical: 10),
            child: Row(children: [for (var c = 0; c < w.columns.length; c++) sized(c, _header(context, c))]),
          ),
          _sized(ListView.builder(itemCount: w.rowCount, itemExtent: w.rowHeight, itemBuilder: (context, r) => _row(context, r, sized))),
        ]),
      );
    }
    return Container(
      decoration: box,
      clipBehavior: Clip.antiAlias,
      child: LayoutBuilder(builder: (context, constraints) {
        final room = constraints.maxWidth - 24;
        final total = f.fold<int>(0, (a, b) => a + b);
        final widths = [
          for (var c = 0; c < w.columns.length; c++)
            _widths[w.columns[c]] ??
                math.max(w.minWidths != null && c < w.minWidths!.length ? w.minWidths![c] : RecordsTable.minColumnWidth, room * f[c] / math.max(1, total)),
        ];
        final content = widths.fold<double>(0, (a, b) => a + b);
        final width = math.max(content, room) + 24;
        Widget sized(int c, Widget child) => SizedBox(width: widths[c], child: child);

        // Header cell with the drag handle on its right edge.
        Widget head(int c) => SizedBox(
              width: widths[c],
              child: Row(children: [
                Expanded(child: _header(context, c)),
                MouseRegion(
                  cursor: SystemMouseCursors.resizeColumn,
                  child: GestureDetector(
                    behavior: HitTestBehavior.opaque,
                    onHorizontalDragUpdate: (d) =>
                        setState(() => _widths[w.columns[c]] = math.max(RecordsTable.minColumnWidth, widths[c] + d.delta.dx)),
                    onDoubleTap: () => setState(() => _widths.remove(w.columns[c])),
                    child: SizedBox(
                      width: 12,
                      height: 22,
                      child: Center(child: Container(width: 2, height: 16, decoration: BoxDecoration(color: t.border, borderRadius: BorderRadius.circular(1)))),
                    ),
                  ),
                ),
              ]),
            );

        final table = SizedBox(
          width: width,
          child: Column(children: [
            Container(
              color: t.surface2,
              padding: const EdgeInsets.symmetric(horizontal: 12, vertical: 8),
              child: Row(children: [for (var c = 0; c < w.columns.length; c++) head(c)]),
            ),
            _sized(ListView.builder(itemCount: w.rowCount, itemExtent: w.rowHeight, itemBuilder: (context, r) => _row(context, r, sized))),
          ]),
        );
        final wide = content > room + 0.5;
        return Scrollbar(
          controller: _h,
          thumbVisibility: wide,
          child: SingleChildScrollView(controller: _h, scrollDirection: Axis.horizontal, child: table),
        );
      }),
    );
  }

  Widget _tappable(int r, Widget row) => w.onRowTap == null
      ? row
      : MouseRegion(cursor: SystemMouseCursors.click, child: GestureDetector(onTap: () => w.onRowTap!(r), child: row));

  Widget _sized(Widget list) =>
      w.fill ? Expanded(child: list) : SizedBox(height: math.min(w.height, w.rowCount * w.rowHeight + 2), child: list);
}

/// Status line with a coloured icon (done / running / error / idle).
class StatusDot extends StatelessWidget {
  const StatusDot({super.key, required this.state, this.tooltip});

  /// 0 idle, 1 running, 2 done, 3 error, 4 skipped.
  final int state;
  final String? tooltip;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final (icon, color) = switch (state) {
      1 => (FA.spinner, t.accent),
      2 => (FA.circleCheck, SplaceTokens.success),
      3 => (FA.circleXmark, SplaceTokens.danger),
      4 => (FA.chevronRight, t.textFaint),
      _ => (FA.circleInfo, t.surface3),
    };
    return SplaceTooltip(message: tooltip ?? '', child: Icon(icon, size: 14, color: color));
  }
}

/// `.seg`: segmented control used in the settings popover.
class Seg<T> extends StatelessWidget {
  const Seg({super.key, required this.value, required this.options, required this.onChanged});
  final T value;
  final List<(T, String)> options;
  final ValueChanged<T> onChanged;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    return Container(
      padding: const EdgeInsets.all(2),
      decoration: BoxDecoration(color: t.surface2, borderRadius: BorderRadius.circular(8), border: Border.all(color: t.border)),
      child: Row(mainAxisSize: MainAxisSize.min, children: [
        for (final (v, label) in options)
          GestureDetector(
            onTap: () => onChanged(v),
            child: MouseRegion(
              cursor: SystemMouseCursors.click,
              child: AnimatedContainer(
                duration: const Duration(milliseconds: 150),
                padding: const EdgeInsets.symmetric(horizontal: 10, vertical: 5),
                decoration: BoxDecoration(
                  color: v == value ? t.surface : Colors.transparent,
                  borderRadius: BorderRadius.circular(6),
                  boxShadow: v == value ? shadowSm : null,
                ),
                child: Text(label,
                    style: TextStyle(fontSize: 12, fontWeight: FontWeight.w600, color: v == value ? t.accent : t.textMuted)),
              ),
            ),
          ),
      ]),
    );
  }
}

/// Heading inside a step that groups GenBank and FASTA content.
class SectionHeading extends StatelessWidget {
  const SectionHeading({super.key, required this.icon, required this.title, this.trailing});
  final IconData icon;
  final String title;
  final Widget? trailing;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    return Padding(
      padding: const EdgeInsets.only(bottom: 12),
      child: Row(children: [
        Icon(icon, size: 14, color: t.accent),
        const SizedBox(width: 8),
        Expanded(child: Text(title, style: TextStyle(fontSize: 14, fontWeight: FontWeight.w600, color: t.text))),
        ?trailing,
      ]),
    );
  }
}

/// Count cell of the records table (`.marker-heatmap-cell`): red when below
/// the dataset's most common count (gene loss / missing annotation), blue
/// when above (duplication).
class CountCell extends StatelessWidget {
  const CountCell({super.key, required this.value, required this.mode, required this.min, required this.max, required this.label});
  final int value, mode, min, max;
  final String label;

  /// Cell background and text for a gene count: red below the mode, blue
  /// above it, none at the mode.
  static (Color?, Color) colors(BuildContext context, int value, int mode, int min, int max) {
    final t = context.tokens;
    if (value < mode) {
      final k = math.max(0.28, (mode - value) / math.max(1, mode - min));
      return (Color.lerp(const Color(0xFFFEF2F2), const Color(0xFFBF0030), k)!, k > 0.62 ? Colors.white : const Color(0xFF8A001F));
    }
    if (value > mode) {
      final k = math.max(0.28, (value - mode) / math.max(1, max - mode));
      return (Color.lerp(const Color(0xFFEEF0FB), const Color(0xFF3A47AB), k)!, k > 0.62 ? Colors.white : const Color(0xFF191B4D));
    }
    return (null, t.text);
  }

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    var bg = t.surface2;
    var fg = t.text;
    var border = t.border;
    if (value < mode) {
      final k = math.max(0.28, (mode - value) / math.max(1, mode - min));
      bg = Color.lerp(const Color(0xFFFEF2F2), const Color(0xFFBF0030), k)!;
      fg = k > 0.62 ? Colors.white : const Color(0xFF8A001F);
      border = const Color(0x8CBF0030);
    } else if (value > mode) {
      final k = math.max(0.28, (value - mode) / math.max(1, max - mode));
      bg = Color.lerp(const Color(0xFFEEF0FB), const Color(0xFF3A47AB), k)!;
      fg = k > 0.62 ? Colors.white : const Color(0xFF191B4D);
      border = const Color(0x8C3A47AB);
    }
    return SplaceTooltip(
      message: '$label: $value (mode $mode)',
      child: Container(
        constraints: const BoxConstraints(minWidth: 38),
        padding: const EdgeInsets.symmetric(horizontal: 8, vertical: 2),
        decoration: BoxDecoration(color: bg, borderRadius: BorderRadius.circular(8), border: Border.all(color: border)),
        child: Text('$value', textAlign: TextAlign.center, style: TextStyle(fontSize: 12, fontWeight: FontWeight.w700, color: fg)),
      ),
    );
  }
}

/// `.splash-hero::before`: the heatmap-palette line, its gradient sliding
/// back and forth (`splashGradientShift`, 4 s, alternate).
class AnimatedRainbowBar extends StatefulWidget {
  const AnimatedRainbowBar({super.key, this.height = 8});
  final double height;

  @override
  State<AnimatedRainbowBar> createState() => _AnimatedRainbowBarState();
}

class _AnimatedRainbowBarState extends State<AnimatedRainbowBar> with SingleTickerProviderStateMixin {
  late final _ctl = AnimationController(vsync: this, duration: const Duration(seconds: 4))..repeat(reverse: true);
  late final _shift = CurvedAnimation(parent: _ctl, curve: Curves.easeInOut);

  @override
  void dispose() {
    _shift.dispose();
    _ctl.dispose();
    super.dispose();
  }

  @override
  Widget build(BuildContext context) => SizedBox(
        height: widget.height,
        child: ClipRect(
          child: LayoutBuilder(builder: (context, box) {
            // background-size: 160%; background-position 0% → 100%.
            final w = box.maxWidth;
            return AnimatedBuilder(
              animation: _shift,
              builder: (context, child) => OverflowBox(
                alignment: Alignment.centerLeft,
                minWidth: w * 1.6,
                maxWidth: w * 1.6,
                child: Transform.translate(offset: Offset(-0.6 * w * _shift.value, 0), child: child),
              ),
              child: const DecoratedBox(decoration: BoxDecoration(gradient: splaceRainbow), child: SizedBox.expand()),
            );
          }),
        ),
      );
}

/// Progress bar: the red → blue gradient spans the whole track and fills it
/// from left to right as the work advances. It only ever shows real progress: without a [value] the
/// track stays empty. [color] replaces the gradient (errors).
class RainbowProgressBar extends StatelessWidget {
  const RainbowProgressBar({super.key, this.value, this.color, this.height = 8});
  final double? value;
  final Color? color;
  final double height;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    return ClipRRect(
      borderRadius: BorderRadius.circular(height),
      child: Container(
        height: height,
        color: t.surface3,
        alignment: Alignment.centerLeft,
        child: LayoutBuilder(
          builder: (context, box) => TweenAnimationBuilder<double>(
            tween: Tween(end: (value ?? 0).clamp(0.0, 1.0)),
            duration: const Duration(milliseconds: 300),
            curve: Curves.easeOut,
            builder: (context, f, _) => SizedBox(
              width: box.maxWidth * f,
              child: ClipRect(
                child: OverflowBox(
                  alignment: Alignment.centerLeft,
                  minWidth: box.maxWidth,
                  maxWidth: box.maxWidth,
                  child: DecoratedBox(
                    decoration: color == null ? const BoxDecoration(gradient: splaceProgress) : BoxDecoration(color: color),
                    child: const SizedBox.expand(),
                  ),
                ),
              ),
            ),
          ),
        ),
      ),
    );
  }
}

/// The SPLACE logo card filling its box (no white margin around it).
class SplaceLogo extends StatelessWidget {
  const SplaceLogo({super.key, required this.height, this.radius = 16, this.shadow = true});
  final double height, radius;
  final bool shadow;

  /// splace_logo_card.png is 720×600.
  static const aspect = 1.2;

  @override
  Widget build(BuildContext context) => Container(
        width: height * aspect,
        height: height,
        decoration: BoxDecoration(
          color: const Color(0xFFEEEEEE),
          borderRadius: BorderRadius.circular(radius),
          boxShadow: shadow ? const [BoxShadow(color: Color(0x14000000), blurRadius: 12, offset: Offset(0, 4))] : null,
        ),
        clipBehavior: Clip.antiAlias,
        child: Image.asset('assets/images/splace_logo_card.png', fit: BoxFit.cover, filterQuality: FilterQuality.medium),
      );
}

/// [text] split into spans where every occurrence of [query] (case and
/// accent insensitive) takes [mark].
List<TextSpan> highlightSpans(String text, String query, TextStyle mark) {
  final q = foldForSearch(query.trim());
  if (q.isEmpty) return [TextSpan(text: text)];
  final folded = foldForSearch(text);
  // Rare characters change length when lower-cased; skip the highlight.
  if (folded.length != text.length) return [TextSpan(text: text)];
  final spans = <TextSpan>[];
  var start = 0;
  for (var i = folded.indexOf(q); i >= 0; i = folded.indexOf(q, i + q.length)) {
    if (i > start) spans.add(TextSpan(text: text.substring(start, i)));
    spans.add(TextSpan(text: text.substring(i, i + q.length), style: mark));
    start = i + q.length;
  }
  if (start < text.length) spans.add(TextSpan(text: text.substring(start)));
  return spans;
}

/// Style of the search highlight: bold on a light accent background.
TextStyle searchMark(BuildContext context) {
  final t = context.tokens;
  return TextStyle(
    backgroundColor: Theme.of(context).brightness == Brightness.dark ? t.accentSoftBorder : SplaceBlue.c100,
    color: t.accentText,
    fontWeight: FontWeight.w700,
  );
}

/// One line of [text] with [query] highlighted.
class HighlightText extends StatelessWidget {
  const HighlightText(this.text, this.query, {super.key, this.style, this.maxLines = 1});
  final String text, query;
  final TextStyle? style;
  final int maxLines;

  @override
  Widget build(BuildContext context) => Text.rich(
        TextSpan(style: style, children: highlightSpans(text, query, searchMark(context))),
        maxLines: maxLines,
        overflow: TextOverflow.ellipsis,
      );
}

/// Number input with − / + buttons, clamped to [min]..[max]. Typing a
/// value applies it when it parses; the buttons move by [step].
class NumberField extends StatefulWidget {
  const NumberField({
    super.key,
    required this.value,
    required this.onChanged,
    this.min = 0,
    this.max = 100,
    this.step = 1,
    this.decimals = 0,
    this.width = 120,
    this.suffix,
  });

  final double value, min, max, step;
  final int decimals;
  final double width;
  final String? suffix;
  final ValueChanged<double> onChanged;

  @override
  State<NumberField> createState() => _NumberFieldState();
}

class _NumberFieldState extends State<NumberField> {
  late final _ctrl = TextEditingController(text: _fmt(widget.value));
  final _focus = FocusNode();

  String _fmt(double v) => v.toStringAsFixed(widget.decimals);

  @override
  void didUpdateWidget(NumberField old) {
    super.didUpdateWidget(old);
    // Follow outside changes, but never rewrite what is being typed.
    if (!_focus.hasFocus && double.tryParse(_ctrl.text.replaceAll(',', '.')) != widget.value) _ctrl.text = _fmt(widget.value);
  }

  @override
  void dispose() {
    _ctrl.dispose();
    _focus.dispose();
    super.dispose();
  }

  void _set(double v) {
    final c = v.clamp(widget.min, widget.max).toDouble();
    _ctrl.text = _fmt(c);
    widget.onChanged(c);
  }

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    Widget btn(IconData icon, double delta) => InkWell(
          borderRadius: BorderRadius.circular(6),
          onTap: () => _set((double.tryParse(_ctrl.text.replaceAll(',', '.')) ?? widget.value) + delta),
          child: SizedBox(width: 26, height: 34, child: Icon(icon, size: 10, color: t.textMuted)),
        );
    return Container(
      width: widget.width,
      height: 36,
      decoration: BoxDecoration(color: t.surface, borderRadius: BorderRadius.circular(8), border: Border.all(color: t.border)),
      child: Row(children: [
        btn(FA.minus, -widget.step),
        Expanded(
          child: TextField(
            controller: _ctrl,
            focusNode: _focus,
            textAlign: TextAlign.center,
            keyboardType: const TextInputType.numberWithOptions(decimal: true),
            style: TextStyle(fontSize: 13, color: t.text),
            decoration: InputDecoration(
              isDense: true,
              border: InputBorder.none,
              enabledBorder: InputBorder.none,
              focusedBorder: InputBorder.none,
              filled: false,
              contentPadding: const EdgeInsets.symmetric(vertical: 9),
              suffixText: widget.suffix,
              suffixStyle: TextStyle(fontSize: 12, color: t.textMuted),
            ),
            onChanged: (s) {
              final v = double.tryParse(s.replaceAll(',', '.'));
              if (v != null && v >= widget.min && v <= widget.max) widget.onChanged(v);
            },
            onSubmitted: (s) => _set(double.tryParse(s.replaceAll(',', '.')) ?? widget.value),
            onTapOutside: (_) {
              _focus.unfocus();
              _ctrl.text = _fmt(widget.value);
            },
          ),
        ),
        btn(FA.plus, widget.step),
      ]),
    );
  }
}


/// Children side by side with equal widths and the height of the tallest
/// (each is laid out once to measure, then again at that height), so cards
/// with buttons at the bottom line up. Works with LayoutBuilder inside the
/// children, unlike IntrinsicHeight.
class EqualHeightRow extends MultiChildRenderObjectWidget {
  const EqualHeightRow({super.key, required super.children, this.gap = 16});
  final double gap;

  @override
  RenderObject createRenderObject(BuildContext context) => RenderEqualHeight(gap);

  @override
  void updateRenderObject(BuildContext context, RenderEqualHeight renderObject) => renderObject.gap = gap;
}

class _EqualHeightParentData extends ContainerBoxParentData<RenderBox> {}

class RenderEqualHeight extends RenderBox
    with ContainerRenderObjectMixin<RenderBox, _EqualHeightParentData>, RenderBoxContainerDefaultsMixin<RenderBox, _EqualHeightParentData> {
  RenderEqualHeight(this._gap);
  double _gap;
  set gap(double v) {
    if (v == _gap) return;
    _gap = v;
    markNeedsLayout();
  }

  @override
  void setupParentData(RenderBox child) {
    if (child.parentData is! _EqualHeightParentData) child.parentData = _EqualHeightParentData();
  }

  @override
  void performLayout() {
    final n = childCount;
    if (n == 0) {
      size = constraints.smallest;
      return;
    }
    final w = (constraints.maxWidth - _gap * (n - 1)) / n;
    var h = 0.0;
    var child = firstChild;
    while (child != null) {
      child.layout(BoxConstraints(minWidth: w, maxWidth: w), parentUsesSize: true);
      h = h > child.size.height ? h : child.size.height;
      child = childAfter(child);
    }
    // A minimum height, not a tight one: a tight child would become a
    // relayout boundary and could later grow alone past the row.
    var x = 0.0, tallest = h;
    child = firstChild;
    while (child != null) {
      child.layout(BoxConstraints(minWidth: w, maxWidth: w, minHeight: h), parentUsesSize: true);
      (child.parentData! as _EqualHeightParentData).offset = Offset(x, 0);
      tallest = tallest > child.size.height ? tallest : child.size.height;
      x += w + _gap;
      child = childAfter(child);
    }
    size = constraints.constrain(Size(constraints.maxWidth, tallest));
  }

  @override
  void paint(PaintingContext context, Offset offset) => defaultPaint(context, offset);

  @override
  bool hitTestChildren(BoxHitTestResult result, {required Offset position}) => defaultHitTestChildren(result, position: position);
}

/// Toggle pill styled like the selected/unselected gene chips.
class MarkerChipLike extends StatelessWidget {
  const MarkerChipLike({super.key, required this.label, required this.selected, required this.onTap});
  final String label;
  final bool selected;
  final VoidCallback onTap;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    return GestureDetector(
      onTap: onTap,
      child: MouseRegion(
        cursor: SystemMouseCursors.click,
        child: AnimatedContainer(
          duration: const Duration(milliseconds: 150),
          padding: const EdgeInsets.symmetric(horizontal: 14, vertical: 8),
          decoration: BoxDecoration(
            color: selected ? SplaceRed.c50 : t.surface,
            borderRadius: BorderRadius.circular(8),
            border: Border.all(color: selected ? SplaceRed.c600 : t.border),
          ),
          child: Row(mainAxisSize: MainAxisSize.min, children: [
            Icon(selected ? FA.scissors : FA.check, size: 12, color: selected ? SplaceRed.c600 : t.textFaint),
            const SizedBox(width: 8),
            Text(label, style: TextStyle(fontSize: 14, fontWeight: FontWeight.w600, color: selected ? SplaceRed.c600 : t.text)),
          ]),
        ),
      ),
    );
  }
}
