import 'package:flutter/material.dart';

import '../../l10n/app_localizations.dart';
import '../core/saturation.dart';
import '../state/project_controller.dart';

const appVersion = '4.2.2';
const githubUrl = 'https://github.com/luanrabelo/splace';

extension L10nX on BuildContext {
  AppLocalizations get l10n => AppLocalizations.of(this);
}

/// Provides the [ProjectController] to the widget tree and rebuilds
/// dependents when it notifies.
class ProjectScope extends InheritedNotifier<ProjectController> {
  const ProjectScope({super.key, required ProjectController controller, required super.child})
      : super(notifier: controller);

  static ProjectController of(BuildContext context) =>
      context.dependOnInheritedWidgetOfExactType<ProjectScope>()!.notifier!;

  /// Reads without subscribing (for callbacks).
  static ProjectController read(BuildContext context) =>
      context.getInheritedWidgetOfExactType<ProjectScope>()!.notifier!;
}

/// Letter of the FASTA part of step 2: after the GenBank parts (2a–2c)
/// when there are GenBank records too.
String fastaHeadersSub(ProjectController c) => c.records.isNotEmpty ? '2d' : '2a';

/// Letter of the FASTA part of step 3 (after the feature types, 3a).
String fastaFieldsSub(ProjectController c) => c.records.isNotEmpty && c.taxonomyResolved ? '3b' : '3a';

int stepState(StepStatus s) => switch (s) {
      StepStatus.idle => 0,
      StepStatus.running => 1,
      StepStatus.done => 2,
      StepStatus.error => 3,
      StepStatus.skipped => 4,
    };

/// Blue (no saturation) → red (saturated) scale of the app palette, used by
/// the saturation heatmap and its legend.
const saturationScaleStops = [
  (0.0, Color(0xFF313695)),
  (0.2, Color(0xFF4575B4)),
  (0.4, Color(0xFFABD9E9)),
  (0.5, Color(0xFFFEE090)),
  (0.62, Color(0xFFFDAE61)),
  (0.78, Color(0xFFF46D43)),
  (1.0, Color(0xFFA50026)),
];

Color saturationScale(double x) {
  final v = x.clamp(0.0, 1.0);
  for (var i = 1; i < saturationScaleStops.length; i++) {
    final (b, cb) = saturationScaleStops[i];
    if (v <= b) {
      final (a, ca) = saturationScaleStops[i - 1];
      return Color.lerp(ca, cb, (v - a) / (b - a))!;
    }
  }
  return saturationScaleStops.last.$2;
}

/// Position of a subset on [saturationScale]: OK cells stay in the blues,
/// moderate ones in the oranges and saturated ones in the reds; inside each
/// band the saturation index (slope and uncorrectable pairs) sets the shade.
/// Null when there is too little data.
double? saturationScore(SubsetSaturation s) {
  final i = s.saturationIndex;
  double within(double lo, double hi) => ((i - lo) / (hi - lo)).clamp(0.0, 1.0);
  return switch (s.verdict) {
    SaturationVerdict.ok => 0.04 + 0.34 * within(0, 0.25),
    SaturationVerdict.moderate => 0.56 + 0.2 * within(0.25, 0.5),
    SaturationVerdict.saturated => 0.8 + 0.2 * within(0.5, 1),
    SaturationVerdict.insufficient => null,
  };
}

/// Colour of a verdict in chips, dots and plots (one point of the scale).
Color verdictColor(SaturationVerdict v) => switch (v) {
      SaturationVerdict.ok => const Color(0xFF4575B4),
      SaturationVerdict.moderate => const Color(0xFFF46D43),
      SaturationVerdict.saturated => const Color(0xFFA50026),
      SaturationVerdict.insufficient => const Color(0xFF9CA3AF),
    };

/// Readable text colour on [background].
Color onColor(Color background) => background.computeLuminance() > 0.42 ? const Color(0xFF1F2937) : Colors.white;

/// Why [s] got its verdict (PhyKIT's saturation value against the
/// thresholds), in the order the verdict is decided.
String verdictReason(BuildContext context, SubsetSaturation s, [SaturationThresholds th = const SaturationThresholds()]) {
  final l10n = context.l10n;
  final slope = s.slope;
  if (s.verdict == SaturationVerdict.insufficient) return l10n.reasonInsufficient('${s.sites}', '${s.pairs}');
  if (s.meanP < th.lowDivergence) return l10n.reasonLowDivergence(fmt(s.meanP, 3));
  if (slope != null && slope < th.slopeSaturated) return l10n.reasonSlopeSaturated(fmt(slope, 3), fmt(th.slopeSaturated));
  if (slope != null && slope < th.slopeModerate) return l10n.reasonSlopeModerate(fmt(slope, 3), fmt(th.slopeModerate));
  return l10n.reasonSlopeOk(fmt(slope, 3));
}

/// "3rd" or "1st and 3rd".
String codonPositions(BuildContext context, Iterable<int> positions) {
  final l10n = context.l10n;
  final names = [
    for (final p in positions.toList()..sort())
      switch (p) {
        1 => l10n.ord1,
        2 => l10n.ord2,
        _ => l10n.ord3,
      }
  ];
  if (names.length < 2) return names.join();
  return '${names.sublist(0, names.length - 1).join(', ')} ${l10n.listAnd} ${names.last}';
}

String verdictLabel(BuildContext context, SaturationVerdict v) => switch (v) {
      SaturationVerdict.ok => context.l10n.verdictOk,
      SaturationVerdict.moderate => context.l10n.verdictModerate,
      SaturationVerdict.saturated => context.l10n.verdictSaturated,
      SaturationVerdict.insufficient => context.l10n.verdictInsufficient,
    };

String fmt(double? v, [int digits = 2]) => v == null ? '—' : v.toStringAsFixed(digits);
String pct(double v, [int digits = 1]) => (v * 100).toStringAsFixed(digits);
