import 'package:flutter/material.dart';

import '../common.dart';
import '../fa.dart';
import 'project_actions.dart';
import '../theme.dart';
import '../widgets/citations.dart';
import '../widgets/kit.dart';
import 'tutorial_dialog.dart';

/// Welcome screen (`.splash-overlay`), shown at start and from the help
/// button.
class SplashOverlay extends StatelessWidget {
  const SplashOverlay({super.key});

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final t = context.tokens;
    final l10n = context.l10n;
    final dark = Theme.of(context).brightness == Brightness.dark;
    final lang = c.locale?.languageCode ?? Localizations.localeOf(context).languageCode;
    final heroBold = TextStyle(fontWeight: FontWeight.w700, color: t.text);
    // Each step in a full sentence; the tools and methods it uses
    // ([[Tool]] in the strings) show as pills with the citation.
    final steps = [
      (l10n.step1, l10n.splashStep1),
      (l10n.step2, l10n.splashStep2),
      (l10n.stepMarkers, l10n.splashStepMarkers),
      (l10n.step5, l10n.splashStep5),
      (l10n.step6, l10n.splashStep6),
      (l10n.step7, l10n.splashStep7),
      (l10n.step8, l10n.splashStep8),
    ];

    // Short windows (720 px) get a smaller hero so the two-line steps fit.
    final compact = MediaQuery.sizeOf(context).height < 820;

    return Material(
      color: t.surface,
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        // Hero
        Container(
          color: dark ? t.surface2 : const Color(0xFFF1F2F6),
          padding: EdgeInsets.fromLTRB(40, compact ? 16 : 26, 40, compact ? 16 : 26),
          child: Stack(children: [
            Row(children: [
              SplaceLogo(height: compact ? 112 : 136, radius: 18),
              SizedBox(width: compact ? 22 : 28),
              Expanded(
                child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
                  Text('SPLACE', style: TextStyle(fontSize: 44, fontWeight: FontWeight.w800, height: 1, letterSpacing: -0.9, color: t.text)),
                  const SizedBox(height: 8),
                  Text.rich(TextSpan(style: TextStyle(fontSize: 18, color: t.textMuted), children: [
                    TextSpan(text: 'SPL', style: heroBold),
                    const TextSpan(text: 'it, '),
                    TextSpan(text: 'A', style: heroBold),
                    const TextSpan(text: 'lign and '),
                    TextSpan(text: 'C', style: heroBold),
                    const TextSpan(text: 'oncatenat'),
                    TextSpan(text: 'E', style: heroBold),
                  ])),
                  const SizedBox(height: 8),
                  Text.rich(TextSpan(children: authorSpans(TextStyle(fontSize: 16, color: t.text)))),
                  const SizedBox(height: 8),
                  Text('v$appVersion · GPL-3.0 · ${githubUrl.replaceFirst('https://', '')}', style: TextStyle(fontSize: 13, color: t.textMuted, letterSpacing: 0.2)),
                ]),
              ),
            ]),
            Positioned(
              top: 0,
              right: 0,
              child: LangToggle(current: lang, onChanged: (code) => c.setLocale(Locale(code))),
            ),
          ]),
        ),
        const AnimatedRainbowBar(height: 9),
        Padding(
          padding: EdgeInsets.fromLTRB(40, compact ? 12 : 18, 40, 6),
          child: Text(l10n.splashSummary, style: TextStyle(fontSize: 14, color: t.textMuted)),
        ),
        // The pipeline as a timeline: one line per step, the title followed
        // by what the step does.
        Expanded(
          child: LayoutBuilder(builder: (context, box) {
            // Two lines at 1.5 line spacing fit in each row.
            final rowHeight = (box.maxHeight / steps.length).clamp(60.0, 72.0);
            return SingleChildScrollView(
              padding: const EdgeInsets.symmetric(horizontal: 40),
              child: Column(children: [
                for (var i = 0; i < steps.length; i++)
                  _TimelineStep(
                    number: i + 1,
                    title: steps[i].$1,
                    body: steps[i].$2,
                    height: rowHeight,
                    first: i == 0,
                    last: i == steps.length - 1,
                  ),
              ]),
            );
          }),
        ),
        Padding(
          padding: EdgeInsets.fromLTRB(40, 10, 40, compact ? 16 : 22),
          child: Row(children: [
            // The tutorial on its own, without starting.
            // In the logo's red, like the step numbers.
            FilledButton(
              style: FilledButton.styleFrom(backgroundColor: SplaceRed.c700, foregroundColor: Colors.white),
              onPressed: () => showTutorial(context),
              child: Text(l10n.tutorialButton),
            ),
            const SizedBox(width: 12),
            // A saved project (.splc) continues where it stopped.
            OutlinedButton.icon(
              onPressed: c.busy ? null : () => openProjectAction(context),
              icon: const Icon(FA.folderOpen, size: 14),
              label: Text(l10n.splashOpenProject),
            ),
            const SizedBox(width: 12),
            Expanded(
              child: FilledButton(
                onPressed: () {
                  // The navigator outlives the welcome screen, which goes
                  // away now; the tutorial opens over the app.
                  final nav = Navigator.of(context);
                  c.setSplashVisible(false);
                  if (c.tutorialAtStart) showTutorial(nav.context);
                },
                child: Text(l10n.splashStart),
              ),
            ),
          ]),
        ),
      ]),
    );
  }
}

/// One step of the welcome timeline: numbered node on a vertical rail and a
/// single-line card with "Title · description".
class _TimelineStep extends StatelessWidget {
  const _TimelineStep({
    required this.number,
    required this.title,
    required this.body,
    required this.height,
    required this.first,
    required this.last,
  });

  final int number;
  final String title, body;
  final double height;
  final bool first, last;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    const node = 36.0, inset = 14.0;
    // The rail runs through the cards at the centre of the numbers and
    // across the gaps between them, linking the steps.
    final rail = Positioned(
      left: inset + node / 2 - 1,
      width: 2,
      top: first ? height / 2 : 0,
      bottom: last ? height / 2 : 0,
      child: ColoredBox(color: t.accentSoftBorder),
    );
    return SizedBox(
      height: height,
      child: Stack(children: [
        Positioned.fill(
          child: Container(
            margin: const EdgeInsets.symmetric(vertical: 3),
            decoration: BoxDecoration(
              color: t.surface2,
              borderRadius: BorderRadius.circular(10),
              border: Border.all(color: t.border),
            ),
          ),
        ),
        rail,
        Positioned.fill(
          child: Padding(
            padding: const EdgeInsets.only(left: inset, right: 16),
            child: Row(children: [
              // Step number: filled circle in the logo's red (the blue and
              // red of the palette).
              Container(
                width: node,
                height: node,
                alignment: Alignment.center,
                decoration: BoxDecoration(
                  color: SplaceRed.c700,
                  shape: BoxShape.circle,
                  border: Border.all(color: t.surface2, width: 2),
                  boxShadow: shadowSm,
                ),
                child: Text('$number', style: const TextStyle(fontWeight: FontWeight.w700, fontSize: 15, color: Colors.white)),
              ),
              const SizedBox(width: 14),
              Expanded(
                child: Text.rich(
                  TextSpan(style: const TextStyle(height: 1.5), children: [
                    TextSpan(text: title, style: TextStyle(fontSize: 14, fontWeight: FontWeight.w600, color: t.text)),
                    TextSpan(text: '  —  ', style: TextStyle(fontSize: 13, color: t.textMuted)),
                    TextSpan(style: TextStyle(fontSize: 13, color: t.textMuted, height: 1.5), children: citationSpans(body, fontSize: 11.5)),
                  ]),
                  maxLines: 2,
                  overflow: TextOverflow.ellipsis,
                ),
              ),
            ]),
          ),
        ),
      ]),
    );
  }
}
