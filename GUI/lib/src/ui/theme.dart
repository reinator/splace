import 'package:flutter/material.dart';

import 'fa.dart';

/// Design tokens copied from `docs/theme.css` so the Flutter app looks like
/// the HTML interface (light and dark themes, comfortable/compact density).
class SplaceBlue {
  static const c50 = Color(0xFFEEF0FB);
  static const c100 = Color(0xFFD5D9F3);
  static const c200 = Color(0xFFA8AFE5);
  static const c300 = Color(0xFF7B85D7);
  static const c400 = Color(0xFF4E5BC9);
  static const c500 = Color(0xFF3A47AB);
  static const c600 = Color(0xFF323795);
  static const c700 = Color(0xFF2A2E7D);
  static const c800 = Color(0xFF212465);
  static const c900 = Color(0xFF191B4D);
}

class SplaceRed {
  static const c50 = Color(0xFFFEF2F2);
  static const c500 = Color(0xFFBF0030);
  static const c600 = Color(0xFFA50026);
  static const c700 = Color(0xFF8A001F);
}

/// The app's colours: the red-to-blue scale of the splash rule and the
/// progress bars (ColorBrewer RdYlBu) and the app blues ([SplaceBlue]).
/// Examples, animations, figures and categories take their colours from
/// here (greys and white aside).
class SplacePalette {
  static const red900 = Color(0xFFA50026);
  static const red700 = Color(0xFFD73027);
  static const red500 = Color(0xFFF46D43);
  static const orange = Color(0xFFFDAE61);
  static const yellow = Color(0xFFFEE090);
  static const sky100 = Color(0xFFE0F3F8);
  static const sky300 = Color(0xFFABD9E9);
  static const sky500 = Color(0xFF74ADD1);
  static const blue600 = Color(0xFF4575B4);
  static const blue900 = Color(0xFF313695);

  /// Categories (taxa, loci, programs, biomes…): blues and reds/oranges in
  /// turn, so neighbours always differ.
  static const categorical = [
    blue600, red500, SplaceBlue.c300, orange, blue900, red900, sky500, red700, SplaceBlue.c500, yellow, sky300, SplaceBlue.c700,
  ];

  /// Colour [i] of a categorical series of any length, never repeated:
  /// the [categorical] colours, then new hues (golden angle, so neighbours
  /// are far apart on the wheel) at three lightness levels.
  static Color categoricalAt(int i) {
    if (i < categorical.length) return categorical[i];
    final k = i - categorical.length;
    const light = [0.40, 0.60, 0.30];
    const sat = [0.62, 0.70, 0.55];
    final level = (k ~/ 7) % 3;
    return HSLColor.fromAHSL(1, (k * 137.508 + 20) % 360, sat[level], light[level]).toColor();
  }
}

/// Gradient of the splash hero rule (`.splash-hero::before`).
const splaceRainbow = LinearGradient(colors: [
  Color(0xFF313695), Color(0xFF4E5BC9), Color(0xFFABD9E9), Color(0xFFFDAE61),
  Color(0xFFF46D43), Color(0xFFA50026), Color(0xFF313695),
], stops: [0, .12, .28, .52, .74, .88, 1]);

/// Progress bars: red at the start of the track to blue at the end, so the
/// bar turns blue as the work finishes.
const splaceProgress = LinearGradient(colors: [
  Color(0xFFA50026), Color(0xFFF46D43), Color(0xFFFDAE61), Color(0xFFABD9E9), Color(0xFF4E5BC9), Color(0xFF313695),
], stops: [0, .2, .4, .6, .82, 1]);

@immutable
class SplaceTokens extends ThemeExtension<SplaceTokens> {
  const SplaceTokens({
    required this.bg,
    required this.surface,
    required this.surface2,
    required this.surface3,
    required this.border,
    required this.borderSoft,
    required this.text,
    required this.textMuted,
    required this.textFaint,
    required this.accent,
    required this.accentSoft,
    required this.accentSoftBorder,
    required this.accentText,
    required this.compact,
  });

  final Color bg, surface, surface2, surface3, border, borderSoft;
  final Color text, textMuted, textFaint;
  final Color accent, accentSoft, accentSoftBorder, accentText;
  final bool compact;

  /// Errors and warnings in the brand reds (the RdYlBu oranges are for
  /// charts and categories only).
  static const danger = SplaceRed.c600;
  static const warn = SplaceRed.c500;
  /// "Done" states use the app blue, not green (readable in both themes).
  static const success = SplaceBlue.c400;

  /// `--pad-card` and `--gap-section`.
  double get padCard => compact ? 14 : 20;
  double get gapSection => compact ? 14 : 20;

  static const light = SplaceTokens(
    bg: Color(0xFFF7F8FB),
    surface: Colors.white,
    surface2: Color(0xFFF3F4F6),
    surface3: Color(0xFFE5E7EB),
    border: Color(0xFFE5E7EB),
    borderSoft: Color(0xFFF1F2F6),
    text: Color(0xFF1F2937),
    textMuted: Color(0xFF6B7280),
    textFaint: Color(0xFF9CA3AF),
    accent: SplaceBlue.c600,
    accentSoft: SplaceBlue.c50,
    accentSoftBorder: SplaceBlue.c100,
    accentText: SplaceBlue.c700,
    compact: false,
  );

  static const dark = SplaceTokens(
    bg: Color(0xFF0F1116),
    surface: Color(0xFF181B24),
    surface2: Color(0xFF1F2230),
    surface3: Color(0xFF2A2E3E),
    border: Color(0xFF2A2E3E),
    borderSoft: Color(0xFF232734),
    text: Color(0xFFE5E7EB),
    textMuted: Color(0xFF9CA3AF),
    textFaint: Color(0xFF6B7280),
    accent: Color(0xFF7B85D7),
    accentSoft: Color(0x1F7B85D7),
    accentSoftBorder: Color(0x477B85D7),
    accentText: Color(0xFFC7CCEF),
    compact: false,
  );

  @override
  SplaceTokens copyWith({bool? compact}) => SplaceTokens(
        bg: bg,
        surface: surface,
        surface2: surface2,
        surface3: surface3,
        border: border,
        borderSoft: borderSoft,
        text: text,
        textMuted: textMuted,
        textFaint: textFaint,
        accent: accent,
        accentSoft: accentSoft,
        accentSoftBorder: accentSoftBorder,
        accentText: accentText,
        compact: compact ?? this.compact,
      );

  @override
  SplaceTokens lerp(SplaceTokens? other, double t) => t < 0.5 ? this : (other ?? this);
}

extension TokensX on BuildContext {
  SplaceTokens get tokens => Theme.of(this).extension<SplaceTokens>()!;
}

const shadowSm = [BoxShadow(color: Color(0x0F0F172A), blurRadius: 2, offset: Offset(0, 1))];
const shadowMd = [BoxShadow(color: Color(0x140F172A), blurRadius: 12, offset: Offset(0, 4))];
const shadowLg = [BoxShadow(color: Color(0x1F0F172A), blurRadius: 32, offset: Offset(0, 12))];

const monoFont = 'JetBrainsMono';

ThemeData buildSplaceTheme(Brightness brightness, {bool compact = false}) {
  final t = (brightness == Brightness.dark ? SplaceTokens.dark : SplaceTokens.light).copyWith(compact: compact);
  final radius8 = BorderRadius.circular(8);
  final scheme = ColorScheme.fromSeed(seedColor: SplaceBlue.c600, brightness: brightness).copyWith(
    primary: t.accent,
    onPrimary: Colors.white,
    secondary: t.accent,
    surface: t.surface,
    onSurface: t.text,
    onSurfaceVariant: t.textMuted,
    outline: t.border,
    outlineVariant: t.borderSoft,
    error: SplaceTokens.danger,
    surfaceContainerLowest: t.surface,
    surfaceContainerLow: t.surface,
    surfaceContainer: t.surface2,
    surfaceContainerHigh: t.surface2,
    surfaceContainerHighest: t.surface3,
    primaryContainer: t.accentSoft,
    onPrimaryContainer: t.accentText,
  );
  TextStyle s(double size, [FontWeight w = FontWeight.w400, Color? c]) =>
      TextStyle(fontFamily: 'Montserrat', fontSize: size, fontWeight: w, color: c ?? t.text, height: 1.4);
  final input = OutlineInputBorder(borderRadius: radius8, borderSide: BorderSide(color: t.border));

  return ThemeData(
    useMaterial3: true,
    brightness: brightness,
    colorScheme: scheme,
    fontFamily: 'Montserrat',
    scaffoldBackgroundColor: t.bg,
    canvasColor: t.surface,
    dividerColor: t.border,
    visualDensity: compact ? VisualDensity.compact : VisualDensity.standard,
    extensions: [t],
    textTheme: TextTheme(
      headlineSmall: s(22, FontWeight.w700),
      titleLarge: s(20, FontWeight.w700),
      titleMedium: s(14),
      titleSmall: s(14, FontWeight.w600),
      // Field values: TextField uses bodyLarge and DropdownButton titleMedium;
      // both are body text so inputs and selects look the same.
      bodyLarge: s(14),
      bodyMedium: s(14),
      bodySmall: s(12.5, FontWeight.w400, t.textMuted),
      labelLarge: s(14, FontWeight.w600),
      labelMedium: s(12.5, FontWeight.w600),
      labelSmall: s(11, FontWeight.w600, t.textMuted),
    ),
    iconTheme: IconThemeData(color: t.textMuted, size: 16),
    dividerTheme: DividerThemeData(color: t.border, space: 1, thickness: 1),
    cardTheme: CardThemeData(
      color: t.surface,
      elevation: 0,
      margin: EdgeInsets.zero,
      shape: RoundedRectangleBorder(borderRadius: BorderRadius.circular(14), side: BorderSide(color: t.border)),
    ),
    inputDecorationTheme: InputDecorationTheme(
      isDense: true,
      filled: true,
      fillColor: brightness == Brightness.dark ? t.surface2 : t.surface,
      contentPadding: const EdgeInsets.symmetric(horizontal: 12, vertical: 12),
      border: input,
      enabledBorder: input,
      disabledBorder: input,
      focusedBorder: OutlineInputBorder(borderRadius: radius8, borderSide: const BorderSide(color: SplaceBlue.c400, width: 2)),
      labelStyle: s(14, FontWeight.w500, t.textMuted),
      floatingLabelStyle: s(14, FontWeight.w600, t.accent),
      hintStyle: s(14, FontWeight.w400, t.textFaint),
      helperStyle: s(12, FontWeight.w400, t.textFaint),
    ),
    filledButtonTheme: FilledButtonThemeData(
      style: FilledButton.styleFrom(
        backgroundColor: t.accent,
        foregroundColor: Colors.white,
        disabledBackgroundColor: SplaceBlue.c600.withValues(alpha: 0.4),
        disabledForegroundColor: Colors.white,
        minimumSize: const Size(0, 44),
        padding: const EdgeInsets.symmetric(horizontal: 20, vertical: 12),
        shape: RoundedRectangleBorder(borderRadius: radius8),
        textStyle: s(15, FontWeight.w600),
      ),
    ),
    outlinedButtonTheme: OutlinedButtonThemeData(
      style: OutlinedButton.styleFrom(
        foregroundColor: t.text,
        side: BorderSide(color: t.border),
        backgroundColor: t.surface,
        // Every action button has the size of the filled one.
        minimumSize: const Size(0, 44),
        padding: const EdgeInsets.symmetric(horizontal: 20, vertical: 12),
        shape: RoundedRectangleBorder(borderRadius: radius8),
        textStyle: s(15, FontWeight.w600),
      ),
    ),
    textButtonTheme: TextButtonThemeData(
      style: TextButton.styleFrom(
        foregroundColor: t.accent,
        minimumSize: const Size(0, 44),
        padding: const EdgeInsets.symmetric(horizontal: 20, vertical: 12),
        shape: RoundedRectangleBorder(borderRadius: radius8),
        textStyle: s(15, FontWeight.w600),
      ),
    ),
    tooltipTheme: TooltipThemeData(
      decoration: BoxDecoration(
        color: brightness == Brightness.dark ? const Color(0xFF1B2234) : SplaceBlue.c600,
        borderRadius: BorderRadius.circular(6),
        boxShadow: shadowMd,
      ),
      textStyle: s(12.5, FontWeight.w500, Colors.white),
      padding: const EdgeInsets.symmetric(horizontal: 12, vertical: 8),
      waitDuration: const Duration(milliseconds: 350),
    ),
    // Material 3 switch: outlined track and small thumb when off, filled
    // track and a check icon when on.
    switchTheme: SwitchThemeData(
      thumbColor: WidgetStateProperty.resolveWith((s) => s.contains(WidgetState.selected) ? Colors.white : t.textFaint),
      trackColor: WidgetStateProperty.resolveWith((s) => s.contains(WidgetState.selected) ? t.accent : t.surface2),
      trackOutlineColor: WidgetStateProperty.resolveWith((s) => s.contains(WidgetState.selected) ? Colors.transparent : t.textFaint),
      thumbIcon: WidgetStateProperty.resolveWith(
        (s) => s.contains(WidgetState.selected) ? Icon(FA.check, size: 10, color: t.accent) : null,
      ),
      materialTapTargetSize: MaterialTapTargetSize.shrinkWrap,
    ),
    chipTheme: ChipThemeData(
      backgroundColor: t.surface,
      selectedColor: t.accentSoft,
      checkmarkColor: t.accent,
      side: BorderSide(color: t.border),
      labelStyle: s(13, FontWeight.w500),
      shape: RoundedRectangleBorder(borderRadius: radius8),
    ),
    badgeTheme: BadgeThemeData(backgroundColor: t.accent, textColor: Colors.white),
    checkboxTheme: CheckboxThemeData(
      fillColor: WidgetStateProperty.resolveWith((s) => s.contains(WidgetState.selected) ? t.accent : Colors.transparent),
      side: BorderSide(color: t.textFaint, width: 1.5),
      shape: RoundedRectangleBorder(borderRadius: BorderRadius.circular(4)),
    ),
    radioTheme: RadioThemeData(fillColor: WidgetStateProperty.resolveWith((s) => s.contains(WidgetState.selected) ? t.accent : t.textFaint)),
    sliderTheme: SliderThemeData(
      activeTrackColor: t.accent,
      inactiveTrackColor: t.surface3,
      thumbColor: t.accent,
      overlayColor: t.accentSoft,
      trackHeight: 4,
      valueIndicatorColor: t.accent,
    ),
    progressIndicatorTheme: ProgressIndicatorThemeData(color: t.accent, linearTrackColor: t.surface3),
    scrollbarTheme: ScrollbarThemeData(
      thickness: const WidgetStatePropertyAll(8),
      radius: const Radius.circular(8),
      thumbColor: WidgetStatePropertyAll(t.textFaint.withValues(alpha: 0.5)),
    ),
    popupMenuTheme: PopupMenuThemeData(color: t.surface, textStyle: s(14)),
    dialogTheme: DialogThemeData(
      backgroundColor: t.surface,
      shape: RoundedRectangleBorder(borderRadius: BorderRadius.circular(18)),
    ),
    segmentedButtonTheme: SegmentedButtonThemeData(
      style: ButtonStyle(
        textStyle: WidgetStatePropertyAll(s(12.5, FontWeight.w600)),
        backgroundColor: WidgetStateProperty.resolveWith((st) => st.contains(WidgetState.selected) ? t.surface : t.surface2),
        foregroundColor: WidgetStateProperty.resolveWith((st) => st.contains(WidgetState.selected) ? t.accent : t.textMuted),
        side: WidgetStatePropertyAll(BorderSide(color: t.border)),
        shape: WidgetStatePropertyAll(RoundedRectangleBorder(borderRadius: radius8)),
      ),
    ),
  );
}
