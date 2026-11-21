import 'dart:async';
import 'dart:convert';
import 'dart:io';
import 'dart:math' as math;
import 'dart:ui' as ui;

import 'package:file_selector/file_selector.dart';
import 'package:flutter/gestures.dart';
import 'package:flutter/material.dart';
import 'package:flutter/services.dart';

import '../../../l10n/app_localizations.dart';
import '../../core/newick.dart';
import '../../core/tree_compare.dart';
import '../common.dart';
import '../fa.dart';
import '../theme.dart';
import 'kit.dart';
import 'alignment_viewer.dart' show residueInk;
import 'salve_widgets.dart' show iucnColors;
import 'splace_tooltip.dart';

part 'tree_circular.dart';
part 'tree_figure_state.dart';
part 'tree_tangle.dart';

enum TreeRooting { asIs, outgroup, midpoint, manual }

enum TreeShape { phylogram, cladogram }

/// How branch support is drawn: numbers (with a coloured dot), only the
/// coloured dots (values in the info bar), or nothing.
enum TreeSupport { numbers, dots, none }

/// Parts a tip label can be built from.
enum TipField { label, species, genus, family, order, klass, phylum, accession, voucher, country, locality, conservation, commonName, endemic, biomes, trend, officialList }

/// Heatmap columns beside the tip labels, from the SALVE sheets.
enum HeatField { status, trend, officialList, threats, endemic, biomes }

/// Colours of the Brazilian biomes (SALVE's names).
const biomeColors = {
  'Amazônia': SplacePalette.blue900,
  'Caatinga': SplacePalette.orange,
  'Cerrado': SplacePalette.red500,
  'Mata Atlântica': SplacePalette.blue600,
  'Pampa': SplacePalette.yellow,
  'Pantanal': SplacePalette.sky500,
  'Marinho': SplaceBlue.c300,
};

/// Population trend (SALVE) on the app scale: declining red, stable
/// yellow, increasing blue; null when unknown.
Color? trendColor(String trend) {
  final v = trend.toLowerCase();
  if (v.startsWith('declin') || v.startsWith('decreas')) return SplacePalette.red700;
  if (v.startsWith('estáv') || v.startsWith('estav') || v.startsWith('stable')) return SplacePalette.yellow;
  if (v.startsWith('aument') || v.startsWith('cresc') || v.startsWith('increas')) return SplacePalette.blue600;
  if (v.startsWith('descon') || v.startsWith('unknown')) return const Color(0xFF9CA3AF);
  return null;
}

/// "Sim"/"Não" as a blue / grey cell.
_HeatCell _yesNo(String v) => switch (v.toLowerCase()) {
      'sim' || 'yes' || 's' => const _HeatCell(Color(0xFF3A47AB)),
      'não' || 'nao' || 'no' || 'n' => const _HeatCell(Color(0xFFD1D5DB)),
      _ => const _HeatCell(null),
    };

/// One cell of the heatmap: its colour (null: empty) and a short text.
class _HeatCell {
  const _HeatCell(this.fill, [this.text = '', this.ink = Colors.white, this.group]);
  final Color? fill;
  final String text;
  final Color ink;

  /// Joined columns: consecutive rows with the same group form one bar.
  final String? group;
}

/// The heatmap: column titles and the cells of each tip (by tip name).
class _Heat {
  const _Heat(this.titles, this.cells, this.cell, this.headerH, [this.gap = 2, this.joined = const {}]);
  final List<String> titles;

  /// Columns drawn as one bar per run of rows of the same group (species
  /// delimitation) instead of a square per row.
  final Set<int> joined;
  final Map<String, List<_HeatCell>> cells;

  /// Side of a cell and height of the turned titles above the columns.
  final double cell, headerH;
  /// Space between the columns.
  final double gap;
  double get width => titles.isEmpty ? 0 : titles.length * (cell + gap) + 6;
}

/// A species delimitation that can be drawn beside the tips: [label] on
/// its button, [title] above its column, the tips of each species.
class TreeDelimitation {
  const TreeDelimitation({required this.id, required this.label, required this.title, required this.species});
  final String id, label, title;
  final List<List<String>> species;
}

/// Colours of the delimitation bars: consecutive species alternate.
const delimitationColors = [SplaceBlue.c700, SplaceBlue.c300];

/// Ranks the rows can be coloured by.
enum TaxRank { species, genus, family, order, klass, phylum, kingdom }

/// One support class: values from [from] up to the next class.
typedef SupportClass = ({double from, Color color});

/// Legend ranges of support classes starting at [froms] (highest first; the
/// last may be 0, "below the previous one"), as (class index, label). Only
/// values at or above [min] ("show support from") get a mark on the tree,
/// so classes below it are left out and the lowest one starts there.
List<(int, String)> supportClassRanges(List<int> froms, int min) {
  String range(int lo, int hi) => lo >= hi ? '$lo' : '$lo–$hi';
  final out = <(int, String)>[];
  for (var i = 0; i < froms.length; i++) {
    final lo = math.max(froms[i], min);
    if (i == 0) {
      out.add((0, '≥ $lo'));
      continue;
    }
    final hi = froms[i - 1] - 1;
    if (lo > hi) continue;
    out.add((i, froms[i] == 0 && min == 0 ? '< ${froms[i - 1]}' : range(lo, hi)));
  }
  return out;
}

/// The three classic classes (strong, moderate, weak) of [supportClassRanges].
List<(int, String)> supportClassLabels(int strong, int moderate, int min) => supportClassRanges([strong, moderate, 0], min);

/// Colours of support classes, strongest first (the app's red-to-blue
/// scale, read from blue).
const supportPalette = [SplacePalette.blue900, SplaceBlue.c400, SplacePalette.sky500, SplacePalette.orange, SplacePalette.red500, SplacePalette.red900];

/// [k] colours spread over [supportPalette].
List<Color> supportColors(int k) => [
      for (var i = 0; i < k; i++) supportPalette[k == 1 ? 0 : (i * (supportPalette.length - 1) / (k - 1)).round()],
    ];

/// Ready-made classes: thresholds for UFBoot/SH-aLRT (0–100, high values
/// mean more) and for the standard bootstrap.
enum SupportPreset { auto, two, three, four, five, quantiles }

/// Where the figure legend goes.
enum LegendPosition { below, right }

/// Shape of the support marks on the nodes.
enum MarkShape { circle, square, diamond, triangle }

/// How a branch joins its parent: right angles, right angles with rounded
/// corners, or a straight diagonal line.
enum BranchStyle { square, rounded, slanted }

/// Pixels (1/96 in) per typographic point: sizes are set in points and the
/// figure width in centimetres, so the exported figure has its real size.
const _pt = 96 / 72;
const _pxPerCm = 96 / 2.54;

/// Usual figure widths of journals (cm): one column, 1.5 columns, full page.
const _widthPresets = [8.5, 12.0, 17.5];

/// Same width and height for every field of the editing panel.
const _fieldW = 176.0, _fieldH = 36.0;

/// Width of the editor's option panel: two fields side by side.
const _panelW = 2 * _fieldW + 14;

/// Text given to a tip in place of its built label.
typedef _Rename = ({String text, bool italic, bool bold});

/// Label drawn as a bar beside a clade (manual, or one per taxon of a rank).
typedef _CladeMark = ({String text, Color color, bool italic});

/// One part of a tip label and its style.
class LabelPart {
  const LabelPart(this.field, {this.italic = false, this.bold = false, this.parens = false});
  final TipField field;
  final bool italic, bold, parens;

  LabelPart copyWith({bool? italic, bool? bold, bool? parens}) =>
      LabelPart(field, italic: italic ?? this.italic, bold: bold ?? this.bold, parens: parens ?? this.parens);
}

/// Ready-made label formats.
const _presets = <String, List<LabelPart>>{
  'original': [LabelPart(TipField.label, italic: true)],
  'species': [LabelPart(TipField.species, italic: true)],
  'speciesVoucher': [LabelPart(TipField.species, italic: true), LabelPart(TipField.voucher, bold: true, parens: true)],
  'speciesAccession': [LabelPart(TipField.species, italic: true), LabelPart(TipField.accession, bold: true, parens: true)],
  // Tips without a locality show the species alone.
  'speciesLocality': [LabelPart(TipField.species, italic: true), LabelPart(TipField.locality, parens: true)],
  'familySpecies': [LabelPart(TipField.family, bold: true), LabelPart(TipField.species, italic: true), LabelPart(TipField.accession, parens: true)],
};

/// Fonts of the figure: the app's font and the usual figure fonts of the
/// system.
const _fonts = ['Montserrat', 'Arial', 'Calibri', 'Times New Roman', 'Georgia', 'Courier New'];

/// Colours of the row groups (one per taxon, in order down the tree).

/// Swatches of the colour picker.
const _swatches = [
  Color(0xFF111827), Color(0xFF374151), Color(0xFF6B7280), Color(0xFF9CA3AF), Color(0xFFD1D5DB), Color(0xFFFFFFFF),
  SplaceBlue.c900, SplaceBlue.c700, SplaceBlue.c500, SplaceBlue.c400, SplaceBlue.c300, SplaceBlue.c100,
  SplacePalette.blue900, SplacePalette.blue600, SplacePalette.sky500, SplacePalette.sky300, SplacePalette.sky100, SplacePalette.yellow,
  SplacePalette.orange, SplacePalette.red500, SplacePalette.red700, SplacePalette.red900, SplaceRed.c700, SplaceRed.c50,
];

enum _Tab { tree, branches, labels, colours, legend, export, compare }

/// Copy of the editable state, for undo and redo.
class _Snap {
  _Snap(this.tree, this.collapsed, this.truncated, this.cladeColors, this.highlights, this.cladeLabels, this.renamed, this.rooting);
  final TreeNode tree;
  final Set<TreeNode> collapsed, truncated;
  final Map<TreeNode, Color> cladeColors, highlights;
  final Map<TreeNode, _CladeMark> cladeLabels;
  final Map<String, _Rename> renamed;
  final TreeRooting rooting;
}

/// Look of the figure: what the painter and the SVG writer need besides the
/// layout. The figure has its own colours (white background by default),
/// so the drawing on screen is what gets exported.
class _Style {
  const _Style({
    required this.family,
    required this.font,
    required this.supportFont,
    required this.legendFont,
    required this.lineW,
    required this.dot,
    required this.branch,
    required this.label,
    required this.background,
    required this.classes,
    required this.support,
    required this.supportMin,
    required this.supportBelow,
    this.supportArrow = false,
    this.supportCurved = false,
    this.markShape = MarkShape.circle,
    this.barWidth = 0,
    this.barFont = 0,
    this.heatGap = 2,
    required this.proportional,
    required this.tipDots,
    required this.lengths,
    required this.showLabels,
    required this.align,
    required this.labelGap,
    required this.groupRows,
    required this.groupLabels,
    required this.legendBox,
    required this.legendTitle,
    required this.title,
    required this.titleFont,
    required this.margin,
    required this.branchStyle,
    required this.markOutgroup,
    required this.highlightAlpha,
  });

  final String family, legendTitle, title;

  /// Sizes in pixels (1/96 in).
  final double font, supportFont, legendFont, titleFont, lineW, dot, supportMin, labelGap, margin, highlightAlpha;
  final Color branch, label, background;

  /// Support classes, highest first.
  final List<SupportClass> classes;

  /// Shape of the support marks.
  final MarkShape markShape;

  /// Clade bars: line width and label size in px (0: from the branches and
  /// the tip labels); gap between heatmap columns in px.
  final double barWidth, barFont, heatGap;

  double get barStroke => barWidth > 0 ? barWidth : math.max(2, lineW * 1.6);
  double get barText => barFont > 0 ? barFont : font;
  final TreeSupport support;
  final BranchStyle branchStyle;
  final bool supportBelow, supportArrow, supportCurved, proportional, tipDots, lengths, showLabels, align, groupRows, groupLabels, legendBox, markOutgroup;

  static const muted = Color(0xFF6B7280);
  static const guide = Color(0xFF9CA3AF);
  static const outgroup = Color(0xFFE5E7EB);

  Color supportColor(double? v) {
    if (v == null || classes.isEmpty) return muted;
    for (final c in classes) {
      if (v >= c.from) return c.color;
    }
    return classes.last.color;
  }

  double radius(double? v) {
    if (proportional) return math.max(1.0, dot * (0.35 + 0.65 * ((v ?? 0) / 100).clamp(0.0, 1.0)));
    return support == TreeSupport.dots ? dot : dot * 0.7;
  }

  /// Title height on top of the figure.
  double get titleH => title.trim().isEmpty ? 0 : titleFont * 1.35 + 8;
}

/// Rectangular tree drawing of a Newick tree: phylogram or cladogram,
/// rooted as inferred, on the outgroup, at the midpoint or on a branch the
/// user picks. Clades can be rotated, collapsed, coloured and used as root;
/// long branches can be shortened (marked `//`); every edit can be undone.
/// Tip labels are built from the taxonomy of the records (species in
/// italic, voucher in bold…), rows, branches and labels can be coloured by
/// any rank, and the legend is drawn in the figure. Exports Newick, SVG
/// and PNG at the chosen DPI.
class TreeViewer extends StatefulWidget {
  const TreeViewer({
    super.key,
    required this.newick,
    this.outgroup = const [],
    this.fileName = 'tree',
    this.bootstrapScale = false,
    this.tipInfo = const {},
    this.title = '',
    this.compareWith = const [],
    this.delimitations = const [],
    this.figure,
    this.onFigureChanged,
    this.onIdentify,
  });

  final String newick;

  /// Figure saved before (options and edits, see [treeFigureVersion]): the
  /// viewer opens as it was left.
  final Map<String, dynamic>? figure;

  /// Called (shortly after the edits stop) with the figure as JSON, so the
  /// project keeps it and saves it in the project file.
  final ValueChanged<Map<String, dynamic>>? onFigureChanged;

  /// "Identify by BLAST" of a tip or clade: gets its tips (names as in
  /// the tree) and gives back the tips to rename in the figure (tip → new
  /// name, in italic), or null.
  final Future<Map<String, String>?> Function(BuildContext context, List<String> tips)? onIdentify;

  /// Species delimitations (mPTP…) that can be shown as bars beside the
  /// tips, each with its own button; tips are matched by name.
  final List<TreeDelimitation> delimitations;

  /// Name of the tree (editor header, comparison heading).
  final String title;

  /// Other trees it can be compared with, face to face.
  final List<TreeOption> compareWith;
  final List<String> outgroup;
  final String fileName;

  /// Standard bootstrap (strong ≥ 70) instead of UFBoot (strong ≥ 95).
  final bool bootstrapScale;

  /// Taxonomy and identifiers by tip name.
  final Map<String, TipInfo> tipInfo;

  @override
  State<TreeViewer> createState() => _TreeViewerState();
}

class _TreeViewerState extends State<TreeViewer> {
  TreeNode? _tree;
  String? _error;
  late TreeRooting _rooting = widget.outgroup.isNotEmpty ? TreeRooting.outgroup : TreeRooting.midpoint;
  TreeShape _shape = TreeShape.phylogram;

  /// Circular drawing: the opening of the circle and where it starts
  /// (degrees, 0 = the first row at the top).
  bool _circular = false;
  double _arcDeg = 350, _rotateDeg = 0;
  TreeSupport _support = TreeSupport.numbers;
  _Tab _tab = _Tab.tree;

  // Figure.
  String _family = 'Montserrat', _title = '', _legendTitle = '';
  bool _align = false, _tipDots = false, _showLabels = true, _lengths = false, _scaleBar = true, _figLegend = true, _proportional = false;
  bool _supportBelow = false, _supportArrow = false, _supportCurved = false, _legendBox = true, _transparent = false, _markOutgroup = true;
  bool _groupRows = true, _groupBranches = false, _groupLabels = false;
  BranchStyle _branchStyle = BranchStyle.square;

  /// Sizes in points (converted to pixels for the drawing).
  double _font = 9, _supportFont = 7, _legendFont = 8, _titleFont = 12, _lineW = 1, _dot = 2.5, _labelGap = 4, _margin = 10;
  double _xScale = 1, _yScale = 1, _dpi = 300, _highlightOpacity = 18;

  /// Width of the figure in centimetres (the same on any screen); its
  /// height follows from the number of rows.
  double _widthCm = 17.5;

  /// Width of the branches in cm (0: automatic, at least 45 % of the
  /// figure) and side of the heatmap cells in pt (0: from the row height).
  double _treeWidthCm = 0, _heatCellPt = 0;

  /// Support mark shape; clade bar width and label size and heatmap column
  /// gap in pt (0: automatic).
  MarkShape _markShape = MarkShape.circle;
  double _barWidthPt = 0, _barFontPt = 0, _heatGapPt = 1.5;

  Color _branchColor = const Color(0xFF111827), _labelColor = const Color(0xFF111827), _background = const Color(0xFFFFFFFF);
  LegendPosition _legendPos = LegendPosition.below;

  /// Length of the scale bar in substitutions per site (0 = automatic).
  double _barLength = 0;

  /// Support below this is not drawn (0 = all).
  double _supportMin = 0;

  /// Support classes: a preset (automatic: 95/80 for UFBoot and SH-aLRT,
  /// 70/50 for the bootstrap) until the user edits them; then their own.
  SupportPreset _supportPreset = SupportPreset.auto;
  List<SupportClass>? _classesCustom;
  bool _supportInLegend = true;
  String _supportLegendTitle = '';
  List<LabelPart> _format = _presets['original']!;
  TaxRank? _colorBy;

  /// SALVE heatmap columns beside the labels.
  final Set<HeatField> _heat = {};

  /// Delimitations shown as bars beside the labels (ids, in the order
  /// they were turned on).
  final List<String> _delims = [];

  List<TreeDelimitation> get _shownDelims => [
        for (final id in _delims) ?widget.delimitations.where((d) => d.id == id).firstOrNull,
      ];

  void _toggleDelim(String id) => setState(() => _delims.contains(id) ? _delims.remove(id) : _delims.add(id));
  String _query = '';
  TreeNode? _selected, _hover;
  final Set<TreeNode> _collapsed = {};

  /// Branches drawn shorter than they are (marked `//`).
  final Set<TreeNode> _truncated = {};

  /// Colours given to clades (their branches).
  final Map<TreeNode, Color> _cladeColors = {};

  /// Coloured boxes behind clades.
  final Map<TreeNode, Color> _highlights = {};

  /// Bars with a name beside clades.
  final Map<TreeNode, _CladeMark> _cladeLabels = {};

  /// Tips given their own text (by tip name).
  final Map<String, _Rename> _renamed = {};

  /// Rank whose taxa get a bar beside their clades (null = none).
  TaxRank? _barRank;

  /// The full-window editor is a route of its own, rebuilt on every change.
  final _changes = ValueNotifier(0);
  Route<void>? _editor;
  final _editorFocus = FocusNode();

  /// Editor canvas: where the figure's corner is on screen and its scale
  /// (screen only; the figure keeps its size in cm).
  Offset _pan = Offset.zero;
  double _scale = 1;
  bool _fitPending = true, _fitPage = false;
  Size _viewport = Size.zero;

  // Comparison with another tree (tanglegram).
  bool _comparing = false;
  int _otherIndex = 0;
  int? _otherLoaded;
  TreeNode? _otherTree;
  String? _otherError;
  bool _untangle = true, _markConflicts = true, _groupConnectors = true;

  /// The tree being edited on the right, the other one on the left.
  bool _swapSides = false;
  Color _conflictColor = SplacePalette.red500;
  double _gapPct = 14;
  String? _tangleHover;
  final List<_Snap> _undo = [], _redo = [];
  final _scroll = ScrollController();
  final _tipCache = <String, TipInfo>{};
  late Map<String, TipInfo> _byKey = _keyed(widget.tipInfo);

  /// Width of the figure in pixels at 96 dpi.
  double get _width => _widthCm * _pxPerCm;

  /// Background of the searched text inside the labels.
  static const _markColor = SplacePalette.yellow;

  /// Row height: 1.7 times the label size, times the vertical scale.
  double get _rowH => _font * _pt * 1.7 * _yScale;
  /// Thresholds of a preset for this tree's kind of support.
  List<double> _presetFroms(SupportPreset p) {
    final boot = widget.bootstrapScale;
    switch (p) {
      case SupportPreset.auto:
      case SupportPreset.three:
        return boot ? [70, 50, 0] : [95, 80, 0];
      case SupportPreset.two:
        return boot ? [70, 0] : [95, 0];
      case SupportPreset.four:
        return boot ? [90, 70, 50, 0] : [98, 95, 80, 0];
      case SupportPreset.five:
        return boot ? [90, 80, 70, 50, 0] : [99, 95, 90, 80, 0];
      case SupportPreset.quantiles:
        // Three classes cut at the tertiles of this tree's own values.
        final values = [for (final n in _tree?.nodes ?? const <TreeNode>[]) if (!n.isTip && n.support != null) n.support!]..sort();
        if (values.length < 3) return _presetFroms(SupportPreset.three);
        final hi = values[(values.length * 2 / 3).floor().clamp(0, values.length - 1)].roundToDouble();
        final lo = values[(values.length / 3).floor().clamp(0, values.length - 1)].roundToDouble();
        return [hi, if (lo < hi) lo, 0];
    }
  }

  List<SupportClass> _classesOf(List<double> froms) {
    final colors = supportColors(froms.length);
    return [for (final (i, f) in froms.indexed) (from: f, color: colors[i])];
  }

  /// The classes in use, in the editor's order.
  List<SupportClass> get _classes => _classesCustom ?? _classesOf(_presetFroms(_supportPreset));

  /// The classes in use, highest first (the figure and its legend).
  List<SupportClass> get _sortedClasses => [..._classes]..sort((a, b) => b.from.compareTo(a.from));

  /// The support classes as the figure shows them (see [supportClassRanges]).
  List<_LegendItem> _supportClasses() {
    final cls = _sortedClasses;
    return [
      for (final (i, label) in supportClassRanges([for (final c in cls) c.from.round()], _supportMin.round())) _LegendItem(label, color: cls[i].color),
    ];
  }

  /// Edits the classes (switching from the preset to the user's own). The
  /// rows keep their order while being typed in; the figure sorts them.
  void _setClasses(List<SupportClass> v) => _classesCustom = v;

  @override
  void initState() {
    super.initState();
    // Records give species and vouchers: start with "Species (voucher)".
    if (widget.tipInfo.values.any((i) => i.voucher.isNotEmpty)) {
      _format = _presets['speciesVoucher']!;
    } else if (widget.tipInfo.isNotEmpty) {
      _format = _presets['speciesAccession']!;
    }
    _build();
    final saved = widget.figure;
    if (saved != null) {
      try {
        _applyFigure(saved);
        _lastFigure = jsonEncode(saved);
      } catch (_) {
        // A figure that cannot be read: the defaults stay.
      }
    }
    _changes.addListener(_figureChanged);
  }

  /// The figure as last reported to [TreeViewer.onFigureChanged].
  String? _lastFigure;
  Timer? _figureTimer;

  void _figureChanged() {
    if (widget.onFigureChanged == null) return;
    _figureTimer?.cancel();
    _figureTimer = Timer(const Duration(milliseconds: 600), _reportFigure);
  }

  void _reportFigure() {
    _figureTimer?.cancel();
    _figureTimer = null;
    final report = widget.onFigureChanged;
    if (report == null || _tree == null) return;
    final j = _figureJson();
    final text = jsonEncode(j);
    if (text == _lastFigure) return;
    _lastFigure = text;
    report(j);
  }

  @override
  void didUpdateWidget(TreeViewer old) {
    super.didUpdateWidget(old);
    if (!identical(old.tipInfo, widget.tipInfo)) {
      _byKey = _keyed(widget.tipInfo);
      _tipCache.clear();
    }
    if (old.newick != widget.newick) {
      _undo.clear();
      _redo.clear();
      _build();
    }
  }

  @override
  void setState(VoidCallback fn) {
    super.setState(fn);
    _changes.value++;
  }

  @override
  void dispose() {
    // Edits made in the last moment are kept too.
    if (_figureTimer != null) _reportFigure();
    _changes.removeListener(_figureChanged);
    final editor = _editor;
    final changes = _changes;
    if (editor != null && editor.isActive) {
      WidgetsBinding.instance.addPostFrameCallback((_) {
        if (editor.isActive) editor.navigator?.removeRoute(editor);
        changes.dispose();
      });
    } else {
      changes.dispose();
    }
    _scroll.dispose();
    _editorFocus.dispose();
    super.dispose();
  }

  // ---------------------------------------------------------------------
  // Tree and edits
  // ---------------------------------------------------------------------

  /// Parses the Newick and roots it as chosen (edits start over).
  void _build() {
    try {
      var t = parseNewick(widget.newick);
      t = switch (_rooting) {
        TreeRooting.outgroup => rootOnOutgroup(t, widget.outgroup),
        TreeRooting.midpoint => midpointRoot(t),
        _ => t,
      };
      ladderize(t);
      _tree = t;
      _error = null;
    } catch (e) {
      _tree = null;
      _error = '$e';
    }
    _selected = null;
    _hover = null;
    _collapsed.clear();
    _truncated.clear();
    _cladeColors.clear();
    _highlights.clear();
    _cladeLabels.clear();
  }

  _Snap _snap() {
    final map = <TreeNode, TreeNode>{};
    TreeNode copy(TreeNode n) {
      final c = TreeNode(name: n.name, length: n.length, label: n.label);
      map[n] = c;
      for (final k in n.children) {
        c.add(copy(k));
      }
      return c;
    }

    final t = copy(_tree!);
    return _Snap(
      t,
      {for (final n in _collapsed) ?map[n]},
      {for (final n in _truncated) ?map[n]},
      {for (final e in _cladeColors.entries) if (map[e.key] != null) map[e.key]!: e.value},
      {for (final e in _highlights.entries) if (map[e.key] != null) map[e.key]!: e.value},
      {for (final e in _cladeLabels.entries) if (map[e.key] != null) map[e.key]!: e.value},
      Map.of(_renamed),
      _rooting,
    );
  }

  void _restore(_Snap s) {
    _tree = s.tree;
    _collapsed
      ..clear()
      ..addAll(s.collapsed);
    _truncated
      ..clear()
      ..addAll(s.truncated);
    _cladeColors
      ..clear()
      ..addAll(s.cladeColors);
    _highlights
      ..clear()
      ..addAll(s.highlights);
    _cladeLabels
      ..clear()
      ..addAll(s.cladeLabels);
    _renamed
      ..clear()
      ..addAll(s.renamed);
    _rooting = s.rooting;
    _selected = null;
    _hover = null;
  }

  /// Applies [change] to the tree, keeping the state before it for undo.
  void _edit(VoidCallback change) {
    if (_tree == null) return;
    setState(() {
      _undo.add(_snap());
      if (_undo.length > 60) _undo.removeAt(0);
      _redo.clear();
      change();
    });
  }

  void _undoEdit() {
    if (_undo.isEmpty) return;
    setState(() {
      _redo.add(_snap());
      _restore(_undo.removeLast());
    });
  }

  void _redoEdit() {
    if (_redo.isEmpty) return;
    setState(() {
      _undo.add(_snap());
      _restore(_redo.removeLast());
    });
  }

  void _setRooting(TreeRooting v) => _edit(() {
        _rooting = v;
        _build();
      });

  void _rootOnSelected() {
    final s = _selected;
    if (s == null) return;
    _edit(() {
      final t = rerootAbove(_tree!, s);
      ladderize(t);
      _tree = t;
      _rooting = TreeRooting.manual;
      _selected = null;
      _collapsed.clear();
    });
  }

  static void _reverse(TreeNode n) {
    final r = n.children.reversed.toList();
    n.children
      ..clear()
      ..addAll(r);
  }

  /// Swaps the order of the selected clade's children.
  void _rotateSelected() {
    final s = _selected;
    if (s == null || s.isTip) return;
    _edit(() => _reverse(s));
  }

  /// Mirrors the whole tree top to bottom.
  void _flip() => _edit(() {
        for (final n in _tree!.nodes) {
          _reverse(n);
        }
      });

  void _ladderize(bool up) => _edit(() => ladderize(_tree!, up: up));

  void _toggleCollapse(TreeNode n) => _edit(() => _collapsed.contains(n) ? _collapsed.remove(n) : _collapsed.add(n));

  void _toggleTruncate(TreeNode n) => _edit(() => _truncated.contains(n) ? _truncated.remove(n) : _truncated.add(n));

  Future<void> _colorClade(TreeNode n) async {
    final c = await _pickColor(context, _cladeColors[n] ?? SplacePalette.categoricalAt(_cladeColors.length));
    if (c != null) _edit(() => _cladeColors[n] = c);
  }

  Future<void> _highlightClade(TreeNode n) async {
    final c = await _pickColor(context, _highlights[n] ?? SplacePalette.categoricalAt(_highlights.length + 2));
    if (c != null) _edit(() => _highlights[n] = c);
  }

  /// Name and colour of the bar beside clade [n].
  Future<void> _labelClade(TreeNode n) async {
    final current = _cladeLabels[n];
    final guess = current?.text ?? _commonTaxon(n) ?? '';
    final r = await _textDialog(
      context,
      title: context.l10n.treeCladeLabel,
      initial: guess,
      italic: current?.italic ?? false,
      color: current?.color ?? const Color(0xFF111827),
    );
    if (r == null) return;
    _edit(() => r.text.trim().isEmpty ? _cladeLabels.remove(n) : _cladeLabels[n] = (text: r.text.trim(), color: r.color!, italic: r.italic));
  }

  /// Own text for tip [n] (empty restores the built label).
  Future<void> _renameTip(TreeNode n) async {
    final current = _renamed[n.name];
    final r = await _textDialog(
      context,
      title: context.l10n.treeRenameTip,
      initial: current?.text ?? _plain(n.name),
      italic: current?.italic ?? false,
      bold: current?.bold ?? false,
    );
    if (r == null) return;
    _edit(() => r.text.trim().isEmpty ? _renamed.remove(n.name) : _renamed[n.name] = (text: r.text.trim(), italic: r.italic, bold: r.bold));
  }

  /// Sends the tips of [n] to [TreeViewer.onIdentify]; the names it gives
  /// back rename the tips in the figure (one undoable edit).
  Future<void> _identify(TreeNode n) async {
    final f = widget.onIdentify;
    if (f == null) return;
    final renames = await f(context, [for (final tip in n.tips) tip.name]);
    if (renames == null || renames.isEmpty || !mounted) return;
    _edit(() {
      for (final e in renames.entries) {
        _renamed[e.key] = (text: e.value, italic: true, bold: false);
      }
    });
  }

  /// Removes tip or clade [n] from the figure (undoable).
  void _prune(TreeNode n) => _edit(() {
        for (final k in n.nodes) {
          _collapsed.remove(k);
          _truncated.remove(k);
          _cladeColors.remove(k);
          _highlights.remove(k);
          _cladeLabels.remove(k);
        }
        _tree = pruneNode(_tree!, n);
        _selected = null;
        _hover = null;
      });

  /// Lowest rank all tips of [n] share (genus, family, order…), for the
  /// default name of its bar.
  String? _commonTaxon(TreeNode n) {
    for (final r in [TaxRank.genus, TaxRank.family, TaxRank.order, TaxRank.klass, TaxRank.phylum]) {
      final v = {for (final t in n.tips) _rank(r, t.name)};
      if (v.length == 1 && v.first.isNotEmpty) return v.first;
    }
    return null;
  }

  /// Distance from the root to the farthest tip.
  double _height(TreeNode root) {
    double h(TreeNode n) => n.isTip ? 0 : n.children.map((c) => (c.length ?? 0) + h(c)).reduce(math.max);
    return h(root);
  }

  double _median(Iterable<double> v) {
    final s = v.where((x) => x > 0).toList()..sort();
    return s.isEmpty ? 0 : s[s.length ~/ 2];
  }

  /// Length a shortened branch is drawn with: three times the median branch
  /// length, at least 5% of the tree height.
  double get _cap {
    final t = _tree!;
    return math.max(3 * _median(t.nodes.map((n) => n.length ?? 0)), 0.05 * _height(t));
  }

  /// Branches much longer than the rest: over 5× the median and over 30% of
  /// the tree height.
  List<TreeNode> _longBranches() {
    final t = _tree!;
    final m = _median(t.nodes.map((n) => n.length ?? 0)), h = _height(t);
    return [
      for (final n in t.nodes)
        if (n.parent != null && (n.length ?? 0) > 5 * m && (n.length ?? 0) > 0.3 * h) n,
    ];
  }

  double _drawnLength(TreeNode n) {
    final l = n.length ?? 0;
    return _truncated.contains(n) ? math.min(l, _cap) : l;
  }

  /// Rows of the drawing: tips and collapsed clades, top to bottom.
  List<TreeNode> _rows(TreeNode tree) {
    final rows = <TreeNode>[];
    void collect(TreeNode n) {
      if (n.isTip || _collapsed.contains(n)) {
        rows.add(n);
      } else {
        n.children.forEach(collect);
      }
    }

    collect(tree);
    return rows;
  }

  // ---------------------------------------------------------------------
  // Labels and groups
  // ---------------------------------------------------------------------

  static String _key(String s) => s.toLowerCase().replaceAll(RegExp(r'[^a-z0-9]+'), '_');

  static Map<String, TipInfo> _keyed(Map<String, TipInfo> m) => {for (final e in m.entries) _key(e.key): e.value};

  TipInfo _info(String name) => _tipCache.putIfAbsent(name, () => widget.tipInfo[name] ?? _byKey[_key(name)] ?? TipInfo.guess(name));

  String _display(String name) => name.replaceAll('_', ' ');

  String _field(TipField f, String name) {
    final i = _info(name);
    return switch (f) {
      TipField.label => _display(name),
      TipField.species => i.species,
      TipField.genus => i.genus,
      TipField.family => i.family,
      TipField.order => i.order,
      TipField.klass => i.klass,
      TipField.phylum => i.phylum,
      TipField.accession => i.accession,
      TipField.voucher => i.voucher,
      TipField.country => i.country,
      TipField.locality => i.locality,
      TipField.conservation => i.status.isEmpty ? '' : (i.statusName.isEmpty ? i.status : '${i.status} · ${i.statusName}'),
      TipField.commonName => i.commonName,
      TipField.endemic => i.endemic.isEmpty ? '' : '${context.l10n.colEndemic}: ${i.endemic}',
      TipField.biomes => i.biomes,
      TipField.trend => i.trend,
      TipField.officialList => i.officialList.isEmpty ? '' : '${context.l10n.colOfficialList}: ${i.officialList}',
    };
  }

  /// Label parts of tip [name] with their text (empty parts are left out;
  /// with none left, the tip name).
  List<(String, LabelPart)> _parts(String name) {
    if (_renamed[name] case final r?) return [(r.text, LabelPart(TipField.label, italic: r.italic, bold: r.bold))];
    final out = <(String, LabelPart)>[
      for (final p in _format)
        if (_field(p.field, name) case final v when v.isNotEmpty) (p.parens ? '($v)' : v, p),
    ];
    return out.isEmpty ? [(_display(name), const LabelPart(TipField.label))] : out;
  }

  String _plain(String name) => _parts(name).map((p) => p.$1).join(' ');

  /// Text of a row: the tip label, or the first tip and the count of a
  /// collapsed clade.
  String _rowText(TreeNode n) {
    if (!_collapsed.contains(n)) return _plain(n.name);
    final name = _cladeLabels[n]?.text;
    return name != null ? '$name (${n.tipCount})' : '${_plain(n.tips.first.name)} … (${n.tipCount})';
  }

  /// Styled label of a row; the searched text is highlighted where it
  /// occurs (on screen only: [search] is false for the exported figure).
  TextSpan _rowSpan(TreeNode n, TextStyle base, {bool search = true}) {
    final q = search ? _query : '';
    const mark = TextStyle(backgroundColor: _markColor, fontWeight: FontWeight.w800);
    if (_collapsed.contains(n)) return TextSpan(style: base, children: highlightSpans(_rowText(n), q, mark));
    final parts = _parts(n.name);
    return TextSpan(style: base, children: [
      for (final (i, (text, p)) in parts.indexed) ...[
        if (i > 0) const TextSpan(text: ' '),
        TextSpan(
          style: TextStyle(fontStyle: p.italic ? FontStyle.italic : FontStyle.normal, fontWeight: p.bold ? FontWeight.w700 : null),
          children: highlightSpans(text, q, mark),
        ),
      ],
    ]);
  }

  String _rank(TaxRank r, String name) {
    final i = _info(name);
    return switch (r) {
      TaxRank.species => i.species,
      TaxRank.genus => i.genus,
      TaxRank.family => i.family,
      TaxRank.order => i.order,
      TaxRank.klass => i.klass,
      TaxRank.phylum => i.phylum,
      TaxRank.kingdom => i.kingdom,
    };
  }

  /// Ranks with a value for at least one tip.
  List<TaxRank> _ranks(TreeNode tree) => [
        for (final r in TaxRank.values)
          if (tree.tips.any((t) => _rank(r, t.name).isNotEmpty)) r,
      ];

  /// Taxon of each row for [_colorBy] (a collapsed clade has one when all
  /// its tips share it) and the colour of each taxon.
  (Map<TreeNode, String>, Map<String, Color>) _groups(List<TreeNode> rows) {
    final by = _colorBy;
    if (by == null) return (const {}, const {});
    final of = <TreeNode, String>{};
    final colors = <String, Color>{};
    for (final n in rows) {
      final values = {for (final t in n.tips) _rank(by, t.name)};
      if (values.length != 1 || values.first.isEmpty) continue;
      final v = values.first;
      of[n] = v;
      // One colour per taxon, however many there are (no repeats).
      colors.putIfAbsent(v, () => SplacePalette.categoricalAt(colors.length));
    }
    return (of, colors);
  }

  /// Colour of the branch above each node: the clade colour of the nearest
  /// coloured ancestor, else (when branches follow the groups) the group all
  /// its tips share, else the branch colour.
  Map<TreeNode, Color> _branchColors(TreeNode tree, Map<TreeNode, String> groupOf, Map<String, Color> groupColors) {
    final out = <TreeNode, Color>{};
    String? group(TreeNode n) {
      if (groupOf.containsKey(n)) return groupOf[n];
      if (n.isTip || _collapsed.contains(n)) return null;
      final gs = {for (final c in n.children) group(c)};
      return gs.length == 1 ? gs.first : null;
    }

    final groupOfNode = <TreeNode, String?>{};
    if (_groupBranches && groupOf.isNotEmpty) {
      for (final n in tree.nodes) {
        groupOfNode[n] = group(n);
      }
    }
    void walk(TreeNode n, Color? inherited) {
      final own = _cladeColors[n] ?? inherited;
      final g = groupOfNode[n];
      out[n] = own ?? (g != null ? groupColors[g]! : _branchColor);
      if (!_collapsed.contains(n)) {
        for (final c in n.children) {
          walk(c, own);
        }
      }
    }

    walk(tree, null);
    return out;
  }

  String _rankLabel(AppLocalizations l10n, TaxRank r) => switch (r) {
        TaxRank.species => l10n.rankSpecies,
        TaxRank.genus => l10n.rankGenus,
        TaxRank.family => l10n.rankFamily,
        TaxRank.order => l10n.rankOrder,
        TaxRank.klass => l10n.rankClass,
        TaxRank.phylum => l10n.rankPhylum,
        TaxRank.kingdom => l10n.rankKingdom,
      };

  String _fieldLabel(AppLocalizations l10n, TipField f) => switch (f) {
        TipField.label => l10n.tipFieldOriginal,
        TipField.species => l10n.rankSpecies,
        TipField.genus => l10n.rankGenus,
        TipField.family => l10n.rankFamily,
        TipField.order => l10n.rankOrder,
        TipField.klass => l10n.rankClass,
        TipField.phylum => l10n.rankPhylum,
        TipField.accession => l10n.tipFieldAccession,
        TipField.voucher => l10n.tipFieldVoucher,
        TipField.country => l10n.tipFieldCountry,
        TipField.locality => l10n.tipFieldLocality,
        TipField.conservation => l10n.tipFieldConservation,
        TipField.commonName => l10n.colCommonNames,
        TipField.endemic => l10n.colEndemic,
        TipField.biomes => l10n.colBiomes,
        TipField.trend => l10n.colTrend,
        TipField.officialList => l10n.colOfficialList,
      };

  // ---------------------------------------------------------------------
  // Figure
  // ---------------------------------------------------------------------

  _Style get _style => _Style(
        family: _family,
        font: _font * _pt,
        supportFont: _supportFont * _pt,
        legendFont: _legendFont * _pt,
        titleFont: _titleFont * _pt,
        lineW: _lineW * _pt,
        dot: _dot * _pt,
        branch: _branchColor,
        label: _labelColor,
        background: _background,
        classes: _sortedClasses,
        support: _support,
        supportMin: _supportMin,
        supportBelow: _supportBelow,
        supportArrow: _supportArrow,
        supportCurved: _supportCurved,
        markShape: _markShape,
        barWidth: _barWidthPt * _pt,
        barFont: _barFontPt * _pt,
        heatGap: _heatGapPt * _pt,
        proportional: _proportional,
        tipDots: _tipDots,
        lengths: _lengths,
        showLabels: _showLabels,
        align: _align,
        labelGap: _labelGap * _pt,
        margin: _margin * _pt,
        branchStyle: _branchStyle,
        markOutgroup: _markOutgroup,
        highlightAlpha: _highlightOpacity / 100,
        groupRows: _groupRows,
        groupLabels: _groupLabels,
        legendBox: _legendBox,
        legendTitle: _legendTitle,
        title: _title,
      );

  /// The figure legend: support classes, outgroup, shortened branches,
  /// clade colours and the taxa of the colour groups, one group per line.
  List<List<_LegendItem>> _legendLines(Map<String, Color> groupColors) {
    final l10n = context.l10n;
    return [
      if (_support != TreeSupport.none && _supportInLegend)
        [
          _LegendItem(_supportLegendTitle.trim().isEmpty ? l10n.treeSupportLegend : _supportLegendTitle.trim(), heading: true),
          ..._supportClasses(),
          if (_supportMin > 0) _LegendItem(l10n.treeLegendUnmarked('${_supportMin.round()}')),
        ],
      if (widget.outgroup.isNotEmpty || _truncated.isNotEmpty)
        [
          _LegendItem(l10n.treeLegendMarks, heading: true),
          if (widget.outgroup.isNotEmpty) _LegendItem(l10n.treeOutgroupLegend, color: _Style.outgroup, box: true),
          if (_truncated.isNotEmpty) _LegendItem(l10n.treeShortenedLegend, symbol: '//'),
        ],
      if (groupColors.isNotEmpty && _colorBy != null)
        [
          _LegendItem(_rankLabel(l10n, _colorBy!), heading: true),
          for (final e in groupColors.entries)
            _LegendItem(e.key, color: e.value.withValues(alpha: _groupRows ? 0.45 : 1), box: _groupRows, italic: _colorBy == TaxRank.species || _colorBy == TaxRank.genus),
        ],
      if (_tree != null) ..._heatLegend(_tree!),
    ];
  }

  /// Some tip has a SALVE sheet (category, endemism, biomes).
  bool get _hasSalve => widget.tipInfo.values.any((i) => i.status.isNotEmpty);

  /// Biomes of the tips, in the order of [biomeColors] (others after).
  List<String> _biomesOf(TreeNode tree) {
    final found = <String>{
      for (final t in tree.tips)
        for (final b in _info(t.name).biomes.split(','))
          if (b.trim().isNotEmpty) b.trim(),
    };
    return [...biomeColors.keys.where(found.contains), ...(found.difference(biomeColors.keys.toSet()).toList()..sort())];
  }

  /// The heatmap columns chosen in the editor (rectangular layout).
  _Heat? _heatOf(TreeNode tree, _Style st) {
    final delims = _shownDelims;
    if ((_heat.isEmpty && delims.isEmpty) || _circular) return null;
    final l10n = context.l10n;
    final biomes = _heat.contains(HeatField.biomes) ? _biomesOf(tree) : const <String>[];
    final maxThreats = tree.tips.map((t) => _info(t.name).threats).fold(0, math.max);
    final titles = [
      if (_heat.contains(HeatField.status)) l10n.treeHeatStatus,
      if (_heat.contains(HeatField.trend)) l10n.colTrend,
      if (_heat.contains(HeatField.officialList)) l10n.colOfficialList,
      if (_heat.contains(HeatField.threats)) l10n.colThreats,
      if (_heat.contains(HeatField.endemic)) l10n.colEndemic,
      ...biomes,
    ];
    final salve = titles.length;
    titles.addAll([for (final d in delims) d.title]);
    if (titles.isEmpty) return null;
    // Delimitations: species number of each tip (by name, also with
    // spaces for underscores), and the species coloured alternately in the
    // order they appear down the tree.
    String key(String n) => n.toLowerCase().replaceAll('_', ' ').trim();
    final delimCells = <String, List<_HeatCell>>{};
    for (final (j, d) in delims.indexed) {
      final of = <String, int>{};
      for (final (i, sp) in d.species.indexed) {
        for (final tip in sp) {
          of[key(tip)] = i;
        }
      }
      final order = <int, int>{};
      for (final t in _rows(tree)) {
        for (final tip in t.tips) {
          final sp = of[key(tip.name)];
          if (sp != null) order.putIfAbsent(sp, () => order.length);
        }
      }
      for (final t in tree.tips) {
        final sp = of[key(t.name)];
        final list = delimCells.putIfAbsent(t.name, () => List.filled(delims.length, const _HeatCell(null)));
        if (sp != null) list[j] = _HeatCell(delimitationColors[(order[sp] ?? sp) % delimitationColors.length], '', Colors.white, '${d.id}#$sp');
      }
    }
    final cells = <String, List<_HeatCell>>{};
    for (final t in tree.tips) {
      final i = _info(t.name);
      final own = {for (final b in i.biomes.split(',')) b.trim()};
      cells[t.name] = [
        if (_heat.contains(HeatField.status))
          i.status.isEmpty ? const _HeatCell(null) : () {
            final (bg, fg) = iucnColors(i.status);
            return _HeatCell(bg, i.status, fg);
          }(),
        if (_heat.contains(HeatField.trend)) _HeatCell(trendColor(i.trend)),
        if (_heat.contains(HeatField.officialList)) _yesNo(i.officialList),
        if (_heat.contains(HeatField.threats))
          i.status.isEmpty
              ? const _HeatCell(null)
              : () {
                  // More threats, deeper red.
                  final c = Color.lerp(SplacePalette.yellow, SplacePalette.red900, maxThreats == 0 ? 0 : i.threats / maxThreats)!;
                  return _HeatCell(c, '${i.threats}', residueInk(c));
                }(),
        if (_heat.contains(HeatField.endemic))
          switch (i.endemic.toLowerCase()) {
            'sim' || 'yes' || 's' => const _HeatCell(Color(0xFF3A47AB)),
            'não' || 'nao' || 'no' || 'n' => const _HeatCell(Color(0xFFD1D5DB)),
            _ => const _HeatCell(null),
          },
        for (final b in biomes) _HeatCell(own.contains(b) ? (biomeColors[b] ?? const Color(0xFF6B7280)) : null),
        ...delimCells[t.name] ?? List.filled(delims.length, const _HeatCell(null)),
      ];
    }
    final cell = _heatCellPt > 0 ? _heatCellPt * _pt : _rowH * 0.78;
    final headerH = titles.map((x) => _Pen.width([_Run(x)], st.family, st.legendFont)).fold(0.0, math.max) + 6;
    return _Heat(titles, cells, cell, headerH, _heatGapPt * _pt, {for (var j = salve; j < titles.length; j++) j});
  }

  /// Legend of the heatmap: the categories present and the biome colours.
  List<List<_LegendItem>> _heatLegend(TreeNode tree) {
    if (_heat.isEmpty || _circular) return const [];
    final l10n = context.l10n;
    final codes = <String>{for (final t in tree.tips) _info(t.name).status}..remove('');
    const order = ['EX', 'EW', 'RE', 'CR', 'EN', 'VU', 'NT', 'LC', 'DD', 'NA', 'NE'];
    return [
      if (_heat.contains(HeatField.status) && codes.isNotEmpty)
        [
          _LegendItem(l10n.treeHeatStatusLegend, heading: true),
          for (final c in [...order.where(codes.contains), ...codes.where((c) => !order.contains(c))]) _LegendItem(c, color: iucnColors(c).$1, box: true),
        ],
      if (_heat.contains(HeatField.trend))
        [
          _LegendItem(l10n.colTrend, heading: true),
          for (final (v, label) in [('declinando', l10n.trendDeclining), ('estável', l10n.trendStable), ('aumentando', l10n.trendIncreasing), ('desconhecida', l10n.trendUnknown)])
            _LegendItem(label, color: trendColor(v) ?? _Style.outgroup, box: true),
        ],
      if (_heat.contains(HeatField.officialList))
        [
          _LegendItem(l10n.colOfficialList, heading: true),
          _LegendItem(l10n.treeHeatYes, color: const Color(0xFF3A47AB), box: true),
          _LegendItem(l10n.treeHeatNo, color: const Color(0xFFD1D5DB), box: true),
        ],
      if (_heat.contains(HeatField.threats))
        [
          _LegendItem(l10n.colThreats, heading: true),
          _LegendItem(l10n.treeHeatFew, color: SplacePalette.yellow, box: true),
          _LegendItem(l10n.treeHeatMany, color: SplacePalette.red900, box: true),
        ],
      if (_heat.contains(HeatField.endemic))
        [
          _LegendItem(l10n.colEndemic, heading: true),
          _LegendItem(l10n.treeHeatYes, color: const Color(0xFF3A47AB), box: true),
          _LegendItem(l10n.treeHeatNo, color: const Color(0xFFD1D5DB), box: true),
        ],
    ];
  }

  /// Bars beside clades: the ones named by hand and, with [_barRank], one
  /// per taxon of that rank (on the largest clade whose tips all share it).
  List<(TreeNode, _CladeMark)> _bars(TreeNode tree, Map<String, Color> groupColors) {
    final out = <(TreeNode, _CladeMark)>[for (final e in _cladeLabels.entries) (e.key, e.value)];
    final rank = _barRank;
    if (rank != null) {
      final italic = rank == TaxRank.genus || rank == TaxRank.species;
      final shared = <TreeNode, String?>{};
      String? value(TreeNode n) => shared.putIfAbsent(n, () {
            final v = {for (final t in n.tips) _rank(rank, t.name)};
            return v.length == 1 && v.first.isNotEmpty ? v.first : null;
          });
      void walk(TreeNode n) {
        final v = value(n);
        if (v != null && !_cladeLabels.containsKey(n)) {
          final color = rank == _colorBy && groupColors[v] != null ? _darker(groupColors[v]!) : _labelColor;
          out.add((n, (text: v, color: color, italic: italic)));
          return;
        }
        if (v == null && !_collapsed.contains(n)) n.children.forEach(walk);
      }

      walk(tree);
    }
    return out;
  }

  _Layout _layout(TreeNode tree, double width, Map<TreeNode, String> groupOf, Map<String, Color> groupColors) {
    final style = _style;
    double measure(TextSpan span) {
      final tp = TextPainter(text: span, textDirection: TextDirection.ltr)..layout();
      return tp.width;
    }

    final base = TextStyle(fontFamily: _family, fontSize: style.font);
    final bold = TextStyle(fontFamily: _family, fontSize: style.font, fontWeight: FontWeight.w700);
    final groupW = groupOf.isEmpty ? 0.0 : {for (final v in groupOf.values) v}.map((v) => measure(TextSpan(text: v, style: bold))).fold(0.0, math.max) + 22;
    final legend = _figLegend ? _legendLines(groupColors) : const <List<_LegendItem>>[];
    return _Layout.compute(
      tree,
      _collapsed,
      _shape,
      (n) => _showLabels ? measure(_rowSpan(n, base, search: false)) + (_tipDots && n.isTip ? style.dot : 0) + style.labelGap : 0,
      width,
      _rowH,
      _drawnLength,
      math.min(groupW, width * 0.25),
      xScale: _xScale,
      barLength: _scaleBar ? _barLength : -1,
      legend: legend,
      legendPos: _legendPos,
      style: style,
      bars: _bars(tree, groupColors),
      heat: _heatOf(tree, style),
      treeW: _treeWidthCm * _pxPerCm,
    );
  }

  _TreePainter _painter(_Layout layout, Map<TreeNode, String> groupOf, Map<String, Color> groupColors, {bool export = false, double zoom = 1}) => _TreePainter(
        layout: layout,
        zoom: zoom,
        highlights: _highlights,
        style: _style,
        accent: SplaceBlue.c500,
        query: export ? '' : _query,
        selected: export ? null : _selected,
        hover: export ? null : _hover,
        span: (n, base) => _rowSpan(n, base, search: !export),
        text: _rowText,
        outgroup: {for (final o in widget.outgroup) o.toLowerCase().replaceAll('_', ' ').trim()},
        truncated: _truncated,
        groupOf: groupOf,
        groupColors: groupColors,
        branchColors: _branchColors(layout.root, groupOf, groupColors),
        transparent: export && _transparent,
      );

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final l10n = context.l10n;
    if (_error != null || _tree == null) {
      return TonePanel(tone: PanelTone.danger, child: Text(l10n.treeParseError(_error ?? ''), style: const TextStyle(fontSize: 13)));
    }
    final tree = _tree!;
    final rows = _rows(tree);
    final (groupOf, groupColors) = _groups(rows);
    final circle = _circular ? _circle(tree, groupOf, groupColors) : null;
    final layout = circle == null ? _layout(tree, _width, groupOf, groupColors) : null;
    final figure = circle?.size ?? layout!.size;

    // On the page: the figure as it will be exported, and the way into the
    // full-window editor.
    return Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
      Wrap(spacing: 8, runSpacing: 8, crossAxisAlignment: WrapCrossAlignment.center, children: [
        PrimaryButton(label: l10n.treeOpenEditor, icon: FA.penToSquare, onPressed: _openEditor),
        if (widget.compareWith.isNotEmpty) ChipButton(label: l10n.treeCompare, icon: FA.clone, keepIcon: true, onPressed: () => _openEditor(compare: true)),
        const SizedBox(width: 8),
        ChipButton(label: 'Newick', onPressed: () => _save('nwk', toNewick(tree))),
        ChipButton(label: 'SVG', onPressed: () => _save('svg', _svg(tree))),
        ChipButton(label: 'PNG (${_dpi.round()} dpi)', onPressed: _savePng),
      ]),
      if (widget.delimitations.isNotEmpty) ...[
        const SizedBox(height: 10),
        _delimButtons(context),
      ],
      const SizedBox(height: 12),
      Container(
        decoration: BoxDecoration(color: t.surface2, borderRadius: BorderRadius.circular(12), border: Border.all(color: t.border)),
        clipBehavior: Clip.antiAlias,
        child: LayoutBuilder(builder: (context, box) {
          // Fits the width of the page (the figure keeps its own size).
          final zoom = ((box.maxWidth - 32) / figure.width).clamp(0.25, 2.0);
          final shown = figure * zoom;
          final height = math.min(shown.height + 32, math.max(320.0, MediaQuery.sizeOf(context).height * 0.6));
          return SizedBox(
            height: height,
            child: Scrollbar(
              controller: _scroll,
              thumbVisibility: shown.height + 32 > height,
              child: SingleChildScrollView(
                controller: _scroll,
                padding: const EdgeInsets.all(16),
                child: Center(
                  child: MouseRegion(
                    cursor: SystemMouseCursors.click,
                    child: GestureDetector(
                      onTap: _openEditor,
                      child: DecoratedBox(
                        decoration: const BoxDecoration(boxShadow: shadowMd),
                        child: CustomPaint(
                          size: shown,
                          painter: circle != null
                              ? _FnPainter((canvas, _) {
                                  canvas.scale(zoom);
                                  _drawCircle(_CanvasPen(canvas), circle, groupOf, groupColors, export: true);
                                })
                              : _painter(layout!, groupOf, groupColors, export: true, zoom: zoom),
                        ),
                      ),
                    ),
                  ),
                ),
              ),
            ),
          );
        }),
      ),
      const SizedBox(height: 6),
      Text(l10n.treeEditorHint, style: TextStyle(fontSize: 12, color: t.textMuted)),
      const SizedBox(height: 4),
      _statsRow(context, tree, figure),
    ]);
  }

  /// One button per species delimitation (tool · tree): each adds or
  /// removes its column of bars beside the labels.
  Widget _delimButtons(BuildContext context, {bool heading = true}) {
    final l10n = context.l10n;
    final t = context.tokens;
    return Wrap(spacing: 8, runSpacing: 6, crossAxisAlignment: WrapCrossAlignment.center, children: [
      if (heading)
      SplaceTooltip(
        message: l10n.treeDelimitationInfo,
        child: Row(mainAxisSize: MainAxisSize.min, children: [
          Kicker(l10n.treeDelimitation, size: 10.5),
          const SizedBox(width: 4),
          Icon(FA.circleInfo, size: 11, color: t.textFaint),
        ]),
      ),
      for (final d in widget.delimitations)
        ChipButton(
          label: d.label,
          icon: _delims.contains(d.id) ? FA.check : null,
          keepIcon: _delims.contains(d.id),
          style: _delims.contains(d.id) ? ChipButtonStyle.accent : ChipButtonStyle.neutral,
          tooltip: _circular ? l10n.treeDelimitationCircular : null,
          onPressed: _circular ? null : () => _toggleDelim(d.id),
        ),
    ]);
  }

  // ---------------------------------------------------------------------
  // Full-window editor
  // ---------------------------------------------------------------------

  void _openEditor({bool compare = false}) {
    if (_editor != null) return;
    setState(() {
      _comparing = compare && widget.compareWith.isNotEmpty;
      if (_comparing) {
        _tab = _Tab.compare;
        _loadOther();
      } else if (_tab == _Tab.compare) {
        _tab = _Tab.tree;
      }
      _fitPending = true;
      _fitPage = false;
      _selected = null;
      _hover = null;
    });
    final route = PageRouteBuilder<void>(
      transitionDuration: const Duration(milliseconds: 160),
      reverseTransitionDuration: const Duration(milliseconds: 120),
      pageBuilder: (context, _, _) => ValueListenableBuilder<int>(valueListenable: _changes, builder: (context, _, _) => _workspace(context)),
      transitionsBuilder: (context, a, _, child) => FadeTransition(opacity: a, child: child),
    );
    _editor = route;
    Navigator.of(context).push(route).whenComplete(() {
      _editor = null;
      if (mounted) setState(() => _hover = null);
    });
  }

  void _closeEditor() {
    final r = _editor;
    if (r != null && r.isActive) r.navigator?.removeRoute(r);
  }

  /// Reads and roots the tree chosen for the comparison.
  void _loadOther() {
    if (widget.compareWith.isEmpty) return;
    final i = _otherIndex.clamp(0, widget.compareWith.length - 1);
    if (_otherLoaded == i) return;
    _otherLoaded = i;
    final o = widget.compareWith[i];
    try {
      var t = parseNewick(o.load());
      try {
        t = o.outgroup.isNotEmpty ? rootOnOutgroup(t, o.outgroup) : midpointRoot(t);
      } catch (_) {
        t = midpointRoot(t);
      }
      ladderize(t);
      _otherTree = t;
      _otherError = null;
    } catch (e) {
      _otherTree = null;
      _otherError = '$e';
    }
  }

  /// Puts the figure in the middle of the canvas: as wide as it fits, or
  /// (with [page]) wholly visible.
  void _fit(Size figure, Size view, {bool page = false}) {
    const pad = 36.0, status = 40.0;
    final sw = (view.width - 2 * pad) / figure.width;
    final sh = (view.height - 2 * pad - status) / figure.height;
    _scale = (page ? math.min(sw, sh) : sw).clamp(0.1, 4.0);
    final w = figure.width * _scale, h = figure.height * _scale;
    _pan = Offset((view.width - w) / 2, page ? math.max(pad, (view.height - status - h) / 2) : pad);
  }

  void _zoomAt(Offset at, double scale) {
    final s = scale.clamp(0.1, 8.0);
    final fig = (at - _pan) / _scale;
    _scale = s;
    _pan = at - fig * s;
  }

  void _wheel(PointerScrollEvent e) {
    final keys = HardwareKeyboard.instance;
    setState(() {
      if (keys.isControlPressed || keys.isMetaPressed) {
        _zoomAt(e.localPosition, _scale * math.pow(0.998, e.scrollDelta.dy).toDouble());
      } else if (keys.isShiftPressed) {
        _pan -= Offset(e.scrollDelta.dy, 0);
      } else {
        _pan -= e.scrollDelta;
      }
    });
  }

  /// Smallest clade holding [a] and [b].
  TreeNode _mrca(TreeNode a, TreeNode b) {
    final up = <TreeNode>{};
    for (TreeNode? n = a; n != null; n = n.parent) {
      up.add(n);
    }
    for (TreeNode? n = b; n != null; n = n.parent) {
      if (up.contains(n)) return n;
    }
    return _tree!;
  }

  void _select(TreeNode? hit) {
    final cur = _selected;
    setState(() {
      _selected = hit != null && cur != null && HardwareKeyboard.instance.isShiftPressed ? _mrca(cur, hit) : hit;
    });
  }

  /// Right-click menu of node [n].
  Future<void> _nodeMenu(Offset global, TreeNode n) async {
    final l10n = context.l10n;
    final tree = _tree!;
    setState(() => _selected = n);
    final truncated = _truncated.contains(n);
    final items = <(String, VoidCallback?, bool)>[
      (l10n.treeRootHere, identical(n, tree) ? null : _rootOnSelected, false),
      if (n.parent != null) (l10n.treeSelectParent, () => setState(() => _selected = n.parent), false),
      if (!n.isTip) (l10n.treeRotateClade, _rotateSelected, false),
      if (!n.isTip) (_collapsed.contains(n) ? l10n.treeExpand : l10n.treeCollapse, () => _toggleCollapse(n), false),
      (l10n.treeColorClade, () => _colorClade(n), false),
      if (_cladeColors.containsKey(n)) (l10n.treeRemoveColor, () => _edit(() => _cladeColors.remove(n)), false),
      (l10n.treeHighlightClade, () => _highlightClade(n), false),
      if (_highlights.containsKey(n)) (l10n.treeRemoveHighlight, () => _edit(() => _highlights.remove(n)), false),
      (l10n.treeCladeLabel, () => _labelClade(n), false),
      if (n.isTip) (l10n.treeRenameTip, () => _renameTip(n), false),
      if (widget.onIdentify != null) (l10n.treeIdentify, () => _identify(n), false),
      if (n.parent != null && (n.length ?? 0) > 0) (truncated ? l10n.treeRestoreBranch : l10n.treeShortenBranch, () => _toggleTruncate(n), false),
      if (n.parent != null && tree.tipCount - n.tipCount >= 2) (l10n.treeRemoveFromFigure, () => _prune(n), true),
    ];
    final t = context.tokens;
    final v = await showMenu<int>(
      context: context,
      position: RelativeRect.fromLTRB(global.dx, global.dy, global.dx + 1, global.dy + 1),
      items: [
        for (final (i, (label, f, danger)) in items.indexed)
          PopupMenuItem<int>(
            value: i,
            enabled: f != null,
            height: 34,
            child: Text(label, style: TextStyle(fontSize: 13, color: danger ? SplaceTokens.danger : t.text, fontWeight: danger ? FontWeight.w600 : null)),
          ),
      ],
    );
    if (v != null) items[v].$2?.call();
  }

  Widget _workspace(BuildContext context) {
    final t = context.tokens;
    final tree = _tree;
    if (tree == null) return const SizedBox.shrink();
    String norm(String s) => s.toLowerCase().replaceAll('_', ' ').trim();
    final matches = _query.isEmpty ? 0 : tree.tips.where((n) => _plain(n.name).toLowerCase().contains(_query) || norm(n.name).contains(_query)).length;
    bool canvasFocused() => _editorFocus.hasPrimaryFocus;
    return Material(
      color: t.bg,
      child: CallbackShortcuts(
        bindings: {
          const SingleActivator(LogicalKeyboardKey.keyZ, control: true): _undoEdit,
          const SingleActivator(LogicalKeyboardKey.keyZ, meta: true): _undoEdit,
          const SingleActivator(LogicalKeyboardKey.keyY, control: true): _redoEdit,
          const SingleActivator(LogicalKeyboardKey.keyZ, control: true, shift: true): _redoEdit,
          const SingleActivator(LogicalKeyboardKey.keyZ, meta: true, shift: true): _redoEdit,
          const SingleActivator(LogicalKeyboardKey.delete): () {
            final s = _selected;
            if (canvasFocused() && s != null && s.parent != null && tree.tipCount - s.tipCount >= 2) _prune(s);
          },
          const SingleActivator(LogicalKeyboardKey.escape): () {
            if (_selected != null) {
              setState(() => _selected = null);
            } else {
              _closeEditor();
            }
          },
        },
        child: Focus(
          focusNode: _editorFocus,
          autofocus: true,
          child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
            _topBar(context, tree),
            Expanded(
              child: Row(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
                _sidePanel(context, tree, matches),
                Expanded(child: _canvas(context, tree)),
              ]),
            ),
          ]),
        ),
      ),
    );
  }

  Widget _topBar(BuildContext context, TreeNode tree) {
    final t = context.tokens;
    final l10n = context.l10n;
    Widget icon(IconData i, String tip, VoidCallback? onPressed) =>
        IconButton(tooltip: tip, onPressed: onPressed, visualDensity: VisualDensity.compact, icon: Icon(i, size: 14));
    Widget divider() => Container(width: 1, height: 26, margin: const EdgeInsets.symmetric(horizontal: 10), color: t.border);
    final pngLabel = 'PNG (${_dpi.round()} dpi)';
    return Container(
      height: 60,
      padding: const EdgeInsets.symmetric(horizontal: 12),
      decoration: BoxDecoration(color: t.surface, border: Border(bottom: BorderSide(color: t.border))),
      child: Row(children: [
        icon(FA.xmark, l10n.close, _closeEditor),
        const SizedBox(width: 6),
        Container(
          width: 34,
          height: 34,
          alignment: Alignment.center,
          decoration: BoxDecoration(color: t.accentSoft, borderRadius: BorderRadius.circular(9), border: Border.all(color: t.accentSoftBorder)),
          child: Icon(FA.diagramProject, size: 14, color: t.accent),
        ),
        const SizedBox(width: 10),
        Flexible(
          child: Column(mainAxisAlignment: MainAxisAlignment.center, crossAxisAlignment: CrossAxisAlignment.start, children: [
            Text(l10n.treeEditorTitle, maxLines: 1, overflow: TextOverflow.ellipsis, style: TextStyle(fontSize: 15, fontWeight: FontWeight.w700, color: t.text)),
            if (widget.title.isNotEmpty)
              Text(widget.title, maxLines: 1, overflow: TextOverflow.ellipsis, style: TextStyle(fontSize: 12, color: t.textMuted)),
          ]),
        ),
        const SizedBox(width: 16),
        if (widget.compareWith.isNotEmpty)
          Seg<bool>(
            value: _comparing,
            options: [(false, l10n.treeModeEdit), (true, l10n.treeCompare)],
            onChanged: (v) => setState(() {
              _comparing = v;
              _selected = null;
              _tangleHover = null;
              if (v) {
                _loadOther();
                _tab = _Tab.compare;
              } else if (_tab == _Tab.compare) {
                _tab = _Tab.tree;
              }
              _fitPending = true;
            }),
          ),
        if (_comparing) ...[
          const SizedBox(width: 8),
          IconBtn(
            icon: FA.arrowRightArrowLeft,
            tooltip: '${l10n.treeSwapSides} · ${l10n.treeSwapSidesInfo}',
            active: _swapSides,
            bordered: true,
            onPressed: () => setState(() {
              _swapSides = !_swapSides;
              _tangleHover = null;
            }),
          ),
        ],
        const Spacer(),
        icon(FA.rotateLeft, '${l10n.treeUndo} (Ctrl+Z)', _undo.isEmpty ? null : _undoEdit),
        icon(FA.rotateRight, '${l10n.treeRedo} (Ctrl+Y)', _redo.isEmpty ? null : _redoEdit),
        const SizedBox(width: 4),
        ChipButton(
          label: l10n.treeReset,
          onPressed: _undo.isEmpty
              ? null
              : () => _edit(() {
                    _rooting = widget.outgroup.isNotEmpty ? TreeRooting.outgroup : TreeRooting.midpoint;
                    _build();
                  }),
        ),
        divider(),
        icon(FA.magnifyingGlassMinus, l10n.treeZoomOut, () => setState(() => _zoomAt(_viewport.center(Offset.zero), _scale / 1.25))),
        SizedBox(
          width: 48,
          child: Text('${(_scale * 100).round()}%',
              textAlign: TextAlign.center, style: TextStyle(fontSize: 12.5, fontWeight: FontWeight.w600, color: t.textMuted)),
        ),
        icon(FA.magnifyingGlassPlus, l10n.treeZoomIn, () => setState(() => _zoomAt(_viewport.center(Offset.zero), _scale * 1.25))),
        const SizedBox(width: 4),
        ChipButton(label: l10n.treeFitWidth, onPressed: () => setState(() => _fitPending = true)),
        const SizedBox(width: 6),
        ChipButton(
          label: l10n.treeFitPage,
          onPressed: () => setState(() {
            _fitPending = true;
            _fitPage = true;
          }),
        ),
        divider(),
        PopupMenuButton<String>(
          tooltip: l10n.treeExport,
          position: PopupMenuPosition.under,
          onSelected: (v) => switch (v) {
            'nwk' => _save('nwk', toNewick(tree)),
            'svg' => _comparing ? _saveTangle(svg: true) : _save('svg', _svg(tree)),
            _ => _comparing ? _saveTangle(svg: false) : _savePng(),
          },
          itemBuilder: (context) => [
            const PopupMenuItem(value: 'nwk', height: 36, child: Text('Newick', style: TextStyle(fontSize: 13))),
            const PopupMenuItem(value: 'svg', height: 36, child: Text('SVG', style: TextStyle(fontSize: 13))),
            PopupMenuItem(value: 'png', height: 36, child: Text(pngLabel, style: const TextStyle(fontSize: 13))),
          ],
          child: IgnorePointer(child: ChipButton(label: l10n.treeExportMenu, icon: FA.fileExport, keepIcon: true, style: ChipButtonStyle.soft, onPressed: () {})),
        ),
        const SizedBox(width: 8),
        FilledButton(onPressed: _closeEditor, child: Text(l10n.treeDone)),
      ]),
    );
  }

  Widget _sidePanel(BuildContext context, TreeNode tree, int matches) {
    final t = context.tokens;
    final l10n = context.l10n;
    final tabs = [
      (_Tab.tree, FA.diagramProject, l10n.treeTabTree),
      (_Tab.branches, FA.codeBranch, l10n.treeTabBranches),
      (_Tab.labels, FA.font, l10n.treeTabLabels),
      (_Tab.colours, FA.palette, l10n.treeTabColours),
      (_Tab.legend, FA.listCheck, l10n.treeFigureLegend),
      (_Tab.export, FA.fileExport, l10n.treeTabFigure),
      if (widget.compareWith.isNotEmpty) (_Tab.compare, FA.clone, l10n.treeTabCompare),
    ];
    Widget railItem(_Tab tab, IconData icon, String label) {
      final on = tab == _tab;
      return Padding(
        padding: const EdgeInsets.only(bottom: 4),
        child: Material(
          color: on ? t.accentSoft : Colors.transparent,
          borderRadius: BorderRadius.circular(10),
          child: InkWell(
            borderRadius: BorderRadius.circular(10),
            onTap: () => setState(() {
              _tab = tab;
              if (tab == _Tab.compare && !_comparing) {
                _comparing = true;
                _loadOther();
                _fitPending = true;
              }
            }),
            child: Container(
              padding: const EdgeInsets.symmetric(vertical: 10, horizontal: 4),
              decoration: BoxDecoration(borderRadius: BorderRadius.circular(10), border: Border.all(color: on ? t.accentSoftBorder : Colors.transparent)),
              child: Column(children: [
                Icon(icon, size: 15, color: on ? t.accent : t.textMuted),
                const SizedBox(height: 6),
                Text(label,
                    textAlign: TextAlign.center,
                    maxLines: 3,
                    style: TextStyle(fontSize: 11, height: 1.2, fontWeight: FontWeight.w600, color: on ? t.accentText : t.textMuted)),
              ]),
            ),
          ),
        ),
      );
    }

    final title = tabs.where((x) => x.$1 == _tab).map((x) => x.$3).firstOrNull ?? '';
    return Container(
      width: 88 + _panelW + 34,
      decoration: BoxDecoration(color: t.surface, border: Border(right: BorderSide(color: t.border))),
      child: Row(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        Container(
          width: 88,
          decoration: BoxDecoration(color: t.surface2, border: Border(right: BorderSide(color: t.border))),
          child: ListView(padding: const EdgeInsets.all(8), children: [for (final (tab, icon, label) in tabs) railItem(tab, icon, label)]),
        ),
        Expanded(
          child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
            Padding(
              padding: const EdgeInsets.fromLTRB(16, 14, 16, 12),
              child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
                TextField(
                  decoration: InputDecoration(isDense: true, hintText: l10n.treeSearch, prefixIcon: Icon(FA.magnifyingGlass, size: 13, color: t.textFaint)),
                  onChanged: (v) => setState(() => _query = v.trim().toLowerCase().replaceAll('_', ' ')),
                ),
                if (_query.isNotEmpty) ...[
                  const SizedBox(height: 6),
                  Text(l10n.treeMatches('$matches'),
                      style: TextStyle(fontSize: 12.5, fontWeight: FontWeight.w600, color: matches > 0 ? t.accentText : t.textFaint)),
                ],
              ]),
            ),
            Divider(height: 1, color: t.border),
            Expanded(
              child: SingleChildScrollView(
                padding: const EdgeInsets.fromLTRB(16, 14, 16, 24),
                child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
                  Text(title, style: TextStyle(fontSize: 15, fontWeight: FontWeight.w700, color: t.text)),
                  const SizedBox(height: 14),
                  Wrap(spacing: 14, runSpacing: 14, crossAxisAlignment: WrapCrossAlignment.end, children: _fields(context, tree, _tab)),
                ]),
              ),
            ),
          ]),
        ),
      ]),
    );
  }

  Widget _canvas(BuildContext context, TreeNode tree) {
    final t = context.tokens;
    final l10n = context.l10n;
    final other = _comparing ? _otherTree : null;
    final rows = _rows(tree);
    final (groupOf, groupColors) = _groups(rows);
    final tangle = other == null ? null : _tangle(tree, _facing(other));
    final circle = tangle == null && _circular ? _circle(tree, groupOf, groupColors) : null;
    final layout = tangle == null && circle == null ? _layout(tree, _width, groupOf, groupColors) : null;
    final figure = tangle?.size ?? circle?.size ?? layout!.size;
    TreeNode? hitAt(Offset p) => circle?.hit(p) ?? layout?.hit(p);
    final rightTitle = _comparing && widget.compareWith.isNotEmpty ? widget.compareWith[_otherIndex.clamp(0, widget.compareWith.length - 1)].label : '';
    return LayoutBuilder(builder: (context, box) {
      final view = box.biggest;
      _viewport = view;
      if (_fitPending) {
        _fit(figure, view, page: _fitPage);
        _fitPending = false;
        _fitPage = false;
        // The zoom shown in the top bar was built before the fit.
        WidgetsBinding.instance.addPostFrameCallback((_) {
          if (mounted) _changes.value++;
        });
      }
      Offset toFigure(Offset p) => (p - _pan) / _scale;
      void paint(Canvas canvas, Size size) {
        final page = _pan & (figure * _scale);
        canvas.drawRect(
          page.shift(const Offset(0, 3)).inflate(2),
          Paint()
            ..color = const Color(0x2A0F172A)
            ..maskFilter = const MaskFilter.blur(BlurStyle.normal, 8),
        );
        canvas.save();
        canvas.clipRect(page);
        canvas.translate(_pan.dx, _pan.dy);
        canvas.scale(_scale);
        if (tangle != null) {
          _drawTangle(_CanvasPen(canvas), tangle, rightTitle, export: false);
        } else if (circle != null) {
          _drawCircle(_CanvasPen(canvas), circle, groupOf, groupColors, export: false);
        } else {
          _painter(layout!, groupOf, groupColors).paint(canvas, figure);
        }
        canvas.restore();
      }

      final sel = _selected;
      final status = tangle != null
          ? (_tangleHover ?? l10n.treeCanvasHint)
          : (_hover != null ? _describe(_hover, l10n) : l10n.treeCanvasHint);
      final (pw, ph) = _pngSize(figure);
      return Stack(children: [
        Positioned.fill(
          child: Listener(
            onPointerSignal: (e) {
              if (e is PointerScrollEvent) GestureBinding.instance.pointerSignalResolver.register(e, (e) => _wheel(e as PointerScrollEvent));
            },
            child: MouseRegion(
              cursor: (tangle == null && _hover != null) ? SystemMouseCursors.click : SystemMouseCursors.grab,
              onHover: (e) {
                final p = toFigure(e.localPosition);
                if (tangle != null) {
                  final k = tangle.hit(p);
                  if (k != _tangleHover) setState(() => _tangleHover = k);
                } else {
                  final n = hitAt(p);
                  if (!identical(n, _hover)) setState(() => _hover = n);
                }
              },
              onExit: (_) => setState(() {
                _hover = null;
                _tangleHover = null;
              }),
              child: GestureDetector(
                behavior: HitTestBehavior.opaque,
                onTapDown: (_) => _editorFocus.requestFocus(),
                onTapUp: (d) {
                  if (tangle == null) _select(hitAt(toFigure(d.localPosition)));
                },
                onSecondaryTapUp: (d) {
                  final n = tangle == null ? hitAt(toFigure(d.localPosition)) : null;
                  if (n != null) _nodeMenu(d.globalPosition, n);
                },
                onPanUpdate: (d) => setState(() => _pan += d.delta),
                child: CustomPaint(size: view, painter: _FnPainter(paint)),
              ),
            ),
          ),
        ),
        if (tangle == null && sel != null) Positioned(left: 14, right: 14, bottom: 46, child: _selectionBar(context, tree, sel)),
        if (_comparing && other == null)
          Positioned(
            left: 14,
            right: 14,
            top: 14,
            child: TonePanel(tone: PanelTone.danger, child: Text(l10n.treeCompareError(_otherError ?? ''), style: const TextStyle(fontSize: 13))),
          ),
        Positioned(
          left: 0,
          right: 0,
          bottom: 0,
          child: Container(
            height: 34,
            padding: const EdgeInsets.symmetric(horizontal: 14),
            decoration: BoxDecoration(color: t.surface.withValues(alpha: 0.94), border: Border(top: BorderSide(color: t.border))),
            child: Row(children: [
              Expanded(
                child: Text(status, maxLines: 1, overflow: TextOverflow.ellipsis, style: TextStyle(fontSize: 12, color: t.textMuted)),
              ),
              const SizedBox(width: 12),
              Text(
                l10n.treeFigureSize((figure.width / _pxPerCm).toStringAsFixed(1), (figure.height / _pxPerCm).toStringAsFixed(1), '$pw', '$ph', '${_dpi.round()}'),
                style: TextStyle(fontSize: 12, fontWeight: FontWeight.w600, color: t.textMuted),
              ),
            ]),
          ),
        ),
      ]);
    });
  }

  /// Options of the comparison: the other tree, untangling, marks and the
  /// distance between the two.
  List<Widget> _compareFields(
    BuildContext context,
    Widget Function(String, Widget, {String? info, bool wide}) group,
    Widget Function<T>(String, T, List<(T, String)>, ValueChanged<T>, {String? info, String? display, bool fonts}) select,
    Widget Function(String, String, bool, VoidCallback, {String? info}) check,
    Widget Function(String, double, double, double, ValueChanged<double>, {double step, int decimals, String? suffix, String? info}) number,
    Widget Function(String, Color, ValueChanged<Color>) color,
  ) {
    final t = context.tokens;
    final l10n = context.l10n;
    final other = _otherTree;
    // In the order drawn, so "left" and "right" in the result match.
    final cmp = other == null ? null : (_swapSides ? compareTrees(other, _tree!) : compareTrees(_tree!, other));
    return [
      group(l10n.treeCompare, Text(l10n.treeCompareIntro, style: TextStyle(fontSize: 12.5, height: 1.45, color: t.textMuted)), wide: true),
      group(
        l10n.treeCompareWith,
        wide: true,
        DropdownButtonFormField<int>(
          key: ValueKey(_otherIndex),
          initialValue: _otherIndex.clamp(0, widget.compareWith.length - 1),
          isExpanded: true,
          icon: dropdownChevron,
          decoration: const InputDecoration(isDense: true),
          items: [
            for (final (i, o) in widget.compareWith.indexed)
              DropdownMenuItem(value: i, child: Text(o.label, overflow: TextOverflow.ellipsis, style: const TextStyle(fontSize: 13))),
          ],
          onChanged: (v) => setState(() {
            _otherIndex = v ?? 0;
            _loadOther();
            _fitPending = true;
          }),
        ),
      ),
      check(l10n.treeSwapSides, l10n.treeSwapSidesText, _swapSides, () => _swapSides = !_swapSides, info: l10n.treeSwapSidesInfo),
      check(l10n.treeUntangle, l10n.treeUntangleText, _untangle, () => _untangle = !_untangle, info: l10n.treeUntangleInfo),
      check(l10n.treeConflicts, l10n.treeConflictsText, _markConflicts, () => _markConflicts = !_markConflicts, info: l10n.treeConflictsInfo),
      color(l10n.treeConflictColor, _conflictColor, (c) => _conflictColor = c),
      check(l10n.treeConnectors, l10n.treeConnectorsGroup, _groupConnectors, () => _groupConnectors = !_groupConnectors, info: l10n.treeConnectorsInfo),
      number(l10n.treeCompareGap, _gapPct, 4, 40, (v) => _gapPct = v, suffix: '%'),
      if (cmp != null)
        group(
          l10n.treeCompareResult,
          info: l10n.treeCompareRfInfo,
          wide: true,
          cmp.common.length < 4
              ? Text(l10n.treeCompareTooFew, style: TextStyle(fontSize: 12.5, color: SplaceTokens.warn))
              : Container(
                  width: _panelW,
                  padding: const EdgeInsets.symmetric(horizontal: 12, vertical: 10),
                  decoration: BoxDecoration(color: t.accentSoft, borderRadius: BorderRadius.circular(8), border: Border.all(color: t.accentSoftBorder)),
                  child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
                    Text(l10n.treeCompareRf('${cmp.rf}', cmp.normalizedRf.toStringAsFixed(2)),
                        style: TextStyle(fontSize: 13, fontWeight: FontWeight.w700, color: cmp.rf == 0 ? SplaceTokens.success : t.accentText)),
                    const SizedBox(height: 4),
                    for (final line in [
                      l10n.treeCompareShared('${cmp.shared}', '${cmp.splitsA.length - cmp.shared}', '${cmp.splitsB.length - cmp.shared}'),
                      l10n.treeCompareCommon('${cmp.common.length}'),
                      if (cmp.onlyA.isNotEmpty || cmp.onlyB.isNotEmpty) l10n.treeCompareOnly('${cmp.onlyA.length}', '${cmp.onlyB.length}'),
                    ])
                      Text(line, style: TextStyle(fontSize: 12.5, height: 1.5, color: t.text)),
                  ]),
                ),
        ),
    ];
  }

  /// Options of editor tab [tab]. Every field (number, colour, text,
  /// choice, on/off) has the same size, so the options line up in two
  /// columns; groups of buttons take the whole width.
  List<Widget> _fields(BuildContext context, TreeNode tree, _Tab tab) {
    final t = context.tokens;
    final l10n = context.l10n;
    final box = BoxDecoration(color: t.surface, borderRadius: BorderRadius.circular(8), border: Border.all(color: t.border));
    Widget group(String label, Widget child, {String? info, bool wide = false}) {
      final w = Column(crossAxisAlignment: CrossAxisAlignment.start, mainAxisSize: MainAxisSize.min, children: [
        // Boxes that continue the group before them have no heading.
        if (label.trim().isNotEmpty) ...[
          Row(children: [
            Flexible(child: Kicker(label, size: 10, singleLine: true)),
            if (info != null) ...[const SizedBox(width: 5), InfoDot(info)],
          ]),
          const SizedBox(height: 6),
        ],
        child,
      ]);
      return SizedBox(width: wide ? _panelW : _fieldW, child: w);
    }

    Widget select<T>(String label, T value, List<(T, String)> options, ValueChanged<T> onChanged, {String? info, String? display, bool fonts = false}) {
      final current = display ?? options.where((o) => o.$1 == value).map((o) => o.$2).firstOrNull ?? '';
      return group(
        label,
        info: info,
        PopupMenuButton<T>(
          tooltip: label,
          position: PopupMenuPosition.under,
          constraints: const BoxConstraints(minWidth: _fieldW),
          onSelected: (v) => setState(() => onChanged(v)),
          itemBuilder: (context) => [
            for (final (v, text) in options)
              PopupMenuItem<T>(
                value: v,
                height: 36,
                child: Row(children: [
                  Expanded(
                    child: Text(text,
                        style: TextStyle(fontSize: 13, fontFamily: fonts ? text : null, color: v == value ? t.accentText : t.text, fontWeight: v == value ? FontWeight.w600 : null)),
                  ),
                  if (v == value) Icon(FA.check, size: 11, color: t.accent),
                ]),
              ),
          ],
          child: Container(
            width: _fieldW,
            height: _fieldH,
            padding: const EdgeInsets.symmetric(horizontal: 10),
            decoration: box,
            child: Row(children: [
              Expanded(
                child: Text(current,
                    maxLines: 1, overflow: TextOverflow.ellipsis, style: TextStyle(fontSize: 13, color: t.text, fontFamily: fonts ? current : null)),
              ),
              Icon(FA.chevronDown, size: 10, color: t.textMuted),
            ]),
          ),
        ),
      );
    }

    Widget check(String label, String text, bool on, VoidCallback onTap, {String? info}) => group(
          label,
          info: info,
          InkWell(
            borderRadius: BorderRadius.circular(8),
            onTap: () => setState(onTap),
            child: Container(
              width: _fieldW,
              height: _fieldH,
              padding: const EdgeInsets.symmetric(horizontal: 8),
              decoration: BoxDecoration(
                color: on ? t.accentSoft : t.surface,
                borderRadius: BorderRadius.circular(8),
                border: Border.all(color: on ? t.accentSoftBorder : t.border),
              ),
              child: Row(children: [
                Container(
                  width: 16,
                  height: 16,
                  alignment: Alignment.center,
                  decoration: BoxDecoration(
                    color: on ? t.accent : Colors.transparent,
                    borderRadius: BorderRadius.circular(4),
                    border: Border.all(color: on ? t.accent : t.textFaint, width: 1.5),
                  ),
                  child: on ? const Icon(FA.check, size: 9, color: Colors.white) : null,
                ),
                const SizedBox(width: 8),
                Expanded(
                  child: Text(text,
                      maxLines: 1,
                      overflow: TextOverflow.ellipsis,
                      style: TextStyle(fontSize: 13, fontWeight: FontWeight.w500, color: on ? t.accentText : t.text)),
                ),
              ]),
            ),
          ),
        );
    Widget number(String label, double value, double min, double max, ValueChanged<double> onChanged,
            {double step = 1, int decimals = 0, String? suffix, String? info}) =>
        group(
          label,
          info: info,
          NumberField(
            value: value,
            min: min,
            max: max,
            step: step,
            decimals: decimals,
            suffix: suffix,
            width: _fieldW,
            height: _fieldH,
            onChanged: (v) => setState(() => onChanged(v)),
          ),
        );
    Widget color(String label, Color value, ValueChanged<Color> onChanged) => group(
          label,
          _ColorButton(color: value, onPick: () async {
            final c = await _pickColor(context, value);
            if (c != null) setState(() => onChanged(c));
          }),
        );
    Widget text(String label, String value, ValueChanged<String> onChanged, {String? hint}) => group(
          label,
          SizedBox(
            width: _fieldW,
            height: _fieldH,
            child: TextFormField(
              initialValue: value,
              style: TextStyle(fontSize: 13, color: t.text),
              decoration: InputDecoration(
                isDense: true,
                hintText: hint,
                hintStyle: TextStyle(fontSize: 12.5, color: t.textFaint),
                contentPadding: const EdgeInsets.symmetric(horizontal: 10, vertical: 10),
              ),
              onChanged: (v) => setState(() => onChanged(v)),
            ),
          ),
        );
    Widget buttons(String label, List<Widget> children, {String? info}) =>
        group(label, info: info, wide: true, Wrap(spacing: 6, runSpacing: 6, crossAxisAlignment: WrapCrossAlignment.center, children: children));

    final long = _longBranches();
    final ranks = _ranks(tree);
    final canRotate = _selected != null && !_selected!.isTip;
    // Ranks by index (-1 = none): a menu item cannot carry null.
    final rankOptions = [(-1, l10n.treeColorNone), for (final r in ranks) (r.index, _rankLabel(l10n, r))];
    final preset = _presets.entries.where((e) => _sameFormat(_format, e.value)).map((e) => e.key).firstOrNull;
    final widthPreset = _widthPresets.where((w) => (w - _widthCm).abs() < 0.01).firstOrNull;
    return switch (tab) {
      _Tab.tree => [
          select<TreeRooting>(
            l10n.treeRootLabel,
            _rooting == TreeRooting.manual ? TreeRooting.asIs : _rooting,
            [
              (TreeRooting.asIs, l10n.treeRootAsIs),
              if (widget.outgroup.isNotEmpty) (TreeRooting.outgroup, l10n.treeRootOutgroup),
              (TreeRooting.midpoint, l10n.treeRootMidpoint),
            ],
            _setRooting,
            display: _rooting == TreeRooting.manual ? l10n.treeRootManual : null,
          ),
          select<bool>(l10n.treeLayoutLabel, _circular, [(false, l10n.treeLayoutRect), (true, l10n.treeLayoutCircle)], (v) {
            _circular = v;
            _fitPending = true;
          }, info: l10n.treeLayoutInfo),
          select<TreeShape>(l10n.treeShapeLabel, _shape, [(TreeShape.phylogram, l10n.treePhylogram), (TreeShape.cladogram, l10n.treeCladogram)], (v) => _shape = v),
          if (_circular) ...[
            number(l10n.treeArc, _arcDeg, 30, 360, (v) => _arcDeg = v, step: 10, suffix: '°', info: l10n.treeArcInfo),
            number(l10n.treeRotation, _rotateDeg, -180, 180, (v) => _rotateDeg = v, step: 15, suffix: '°'),
          ],
          select<BranchStyle>(
            l10n.treeBranchStyle,
            _branchStyle,
            [(BranchStyle.square, l10n.treeBranchSquare), (BranchStyle.rounded, l10n.treeBranchRounded), (BranchStyle.slanted, l10n.treeBranchSlanted)],
            (v) => _branchStyle = v,
          ),
          if (!_circular) ...[
            number(l10n.treeBranchWidth, _treeWidthCm, 0, 80, (v) => _treeWidthCm = v, step: 0.5, decimals: 1, suffix: 'cm', info: l10n.treeBranchWidthInfo),
            number(l10n.treeXScale, _xScale, 0.2, 10, (v) => _xScale = v, step: 0.1, decimals: 1, suffix: '×', info: l10n.treeXScaleInfo),
            number(l10n.treeYScale, _yScale, 0.5, 4, (v) => _yScale = v, step: 0.1, decimals: 1, suffix: '×', info: l10n.treeYScaleInfo),
          ],
          buttons(l10n.treeOrderLabel, info: l10n.treeRotateInfo, [
            ChipButton(label: l10n.treeLadderUp, onPressed: () => _ladderize(true)),
            ChipButton(label: l10n.treeLadderDown, onPressed: () => _ladderize(false)),
            ChipButton(label: l10n.treeFlip, onPressed: _flip),
            ChipButton(label: l10n.treeRotateClade, onPressed: canRotate ? _rotateSelected : null),
          ]),
        ],
      _Tab.branches => [
          number(l10n.treeLineWidth, _lineW, 0.25, 6, (v) => _lineW = v, step: 0.25, decimals: 2, suffix: 'pt'),
          color(l10n.treeBranchColor, _branchColor, (c) => _branchColor = c),
          check(l10n.treeBranchLengths, l10n.treeShowLengths, _lengths, () => _lengths = !_lengths),
          select<TreeSupport>(
            l10n.treeShowSupport,
            _support,
            [(TreeSupport.numbers, l10n.treeSupportValues), (TreeSupport.dots, l10n.treeSupportDots), (TreeSupport.none, l10n.treeSupportHidden)],
            (v) => _support = v,
          ),
          number(l10n.treeSupportMin, _supportMin, 0, 100, (v) => _supportMin = v, step: 5, info: l10n.treeSupportMinInfo),
          select<bool>(l10n.treeSupportPosition, _supportBelow, [(false, l10n.treeAboveBranch), (true, l10n.treeBelowBranch)], (v) => _supportBelow = v),
          if (_support == TreeSupport.numbers) ...[
            check(l10n.treeSupportArrowLabel, l10n.treeSupportArrow, _supportArrow, () => _supportArrow = !_supportArrow, info: l10n.treeSupportArrowInfo),
            select<bool>(l10n.treeArrowShape, _supportCurved, [(false, l10n.treeArrowStraight), (true, l10n.treeArrowCurved)], (v) => _supportCurved = v,
                info: l10n.treeArrowShapeInfo),
          ],
          number(l10n.treeSupportFont, _supportFont, 4, 24, (v) => _supportFont = v, step: 0.5, decimals: 1, suffix: 'pt'),
          number(l10n.treePointSize, _dot, 0.5, 12, (v) => _dot = v, step: 0.5, decimals: 1, suffix: 'pt'),
          select<MarkShape>(
            l10n.treeMarkShape,
            _markShape,
            [
              (MarkShape.circle, l10n.treeShapeCircle),
              (MarkShape.square, l10n.treeShapeSquare),
              (MarkShape.diamond, l10n.treeShapeDiamond),
              (MarkShape.triangle, l10n.treeShapeTriangle),
            ],
            (v) => _markShape = v,
          ),
          check(l10n.treePoints, l10n.treeProportional, _proportional, () => _proportional = !_proportional),
          check(l10n.treeTipDotsLabel, l10n.treeTipDots, _tipDots, () => _tipDots = !_tipDots),
          group(l10n.treeSupportClasses, _supportClassesEditor(context), info: l10n.treeSupportClassesInfo, wide: true),
          buttons(l10n.treeLongBranches, info: l10n.treeLongBranchesInfo, [
            ChipButton(
              label: l10n.treeShortenLong('${long.where((n) => !_truncated.contains(n)).length}'),
              onPressed: long.every(_truncated.contains) ? null : () => _edit(() => _truncated.addAll(long)),
            ),
            ChipButton(label: l10n.treeRestoreLengths, onPressed: _truncated.isEmpty ? null : () => _edit(_truncated.clear)),
          ]),
        ],
      _Tab.labels => [
          check(l10n.treeDisplay, l10n.treeShowLabels, _showLabels, () => _showLabels = !_showLabels),
          select<String>(
            l10n.treeLabelFormat,
            preset ?? '',
            [for (final k in _presets.keys) (k, _presetLabel(l10n, k))],
            (k) => _format = _presets[k]!,
            display: preset == null ? l10n.treeLabelCustom : null,
          ),
          select<String>(l10n.treeFontFamily, _family, [for (final f in _fonts) (f, f)], (v) => _family = v, fonts: true),
          number(l10n.treeFont, _font, 4, 32, (v) => _font = v, step: 0.5, decimals: 1, suffix: 'pt'),
          color(l10n.treeLabelColor, _labelColor, (c) => _labelColor = c),
          select<bool>(l10n.treeLabelPosition, _align, [(false, l10n.treeAtTips), (true, l10n.treeAlignLabels)], (v) => _align = v),
          number(l10n.treeLabelGap, _labelGap, 0, 40, (v) => _labelGap = v, step: 0.5, decimals: 1, suffix: 'pt', info: l10n.treeLabelGapInfo),
          if (widget.outgroup.isNotEmpty) check(l10n.treeOutgroupLegend, l10n.treeMarkOutgroup, _markOutgroup, () => _markOutgroup = !_markOutgroup),
          group(l10n.treeLabelParts, _partsEditor(context), info: l10n.treeLabelPartsInfo, wide: true),
          // SALVE (Brazilian fauna): category, endemism and biomes as
          // coloured columns beside the labels; as text, they are label
          // parts above.
          if (_hasSalve)
            buttons(l10n.treeHeatmap, info: l10n.treeHeatmapInfo, [
              for (final f in HeatField.values)
                ChipButton(
                  label: switch (f) {
                    HeatField.status => l10n.treeHeatStatus,
                    HeatField.trend => l10n.colTrend,
                    HeatField.officialList => l10n.colOfficialList,
                    HeatField.threats => l10n.colThreats,
                    HeatField.endemic => l10n.colEndemic,
                    HeatField.biomes => l10n.colBiomes,
                  },
                  icon: _heat.contains(f) ? FA.check : null,
                  keepIcon: _heat.contains(f),
                  style: _heat.contains(f) ? ChipButtonStyle.accent : ChipButtonStyle.neutral,
                  onPressed: _circular ? null : () => setState(() => _heat.contains(f) ? _heat.remove(f) : _heat.add(f)),
                ),
            ]),
          if (widget.delimitations.isNotEmpty) group(l10n.treeDelimitation, _delimButtons(context, heading: false), info: l10n.treeDelimitationInfo, wide: true),
          if ((_hasSalve && _heat.isNotEmpty) || _delims.isNotEmpty) ...[
            number(l10n.treeHeatCell, _heatCellPt, 0, 40, (v) => _heatCellPt = v, step: 0.5, decimals: 1, suffix: 'pt', info: l10n.treeHeatCellInfo),
            number(l10n.treeHeatGap, _heatGapPt, 0, 20, (v) => _heatGapPt = v, step: 0.5, decimals: 1, suffix: 'pt'),
            number(l10n.treeBranchWidth, _treeWidthCm, 0, 80, (v) => _treeWidthCm = v, step: 0.5, decimals: 1, suffix: 'cm', info: l10n.treeBranchWidthInfo),
          ],
          buttons(l10n.treeRenamed, info: l10n.treeRenamedInfo, [
            ChipButton(label: l10n.treeRenameTip, onPressed: _selected != null && _selected!.isTip ? () => _renameTip(_selected!) : null),
            ChipButton(label: l10n.treeRestoreNames('${_renamed.length}'), onPressed: _renamed.isEmpty ? null : () => _edit(_renamed.clear)),
          ]),
        ],
      _Tab.colours => [
          if (ranks.isEmpty)
            group(l10n.treeColorRows, Text(l10n.treeNoTaxonomy, style: TextStyle(fontSize: 12.5, color: t.textMuted)), wide: true)
          else ...[
            select<int>(l10n.treeColorRows, _colorBy?.index ?? -1, rankOptions, (v) => _colorBy = v < 0 ? null : TaxRank.values[v], info: l10n.treeColorRowsInfo),
            if (_colorBy != null) ...[
              check(l10n.treeColorApply, l10n.treeColorApplyRows, _groupRows, () => _groupRows = !_groupRows),
              // One heading for the three boxes (blank above the others).
              check(' ', l10n.treeColorApplyBranches, _groupBranches, () => _groupBranches = !_groupBranches),
              check(' ', l10n.treeColorApplyLabels, _groupLabels, () => _groupLabels = !_groupLabels),
            ],
            select<int>(l10n.treeBarsByRank, _barRank?.index ?? -1, rankOptions, (v) => _barRank = v < 0 ? null : TaxRank.values[v], info: l10n.treeBarsByRankInfo),
          ],
          color(l10n.treeBackground, _background, (c) => _background = c),
          number(l10n.treeHighlightOpacity, _highlightOpacity, 5, 60, (v) => _highlightOpacity = v, step: 5, suffix: '%'),
          number(l10n.treeBarWidth, _barWidthPt, 0, 20, (v) => _barWidthPt = v, step: 0.5, decimals: 1, suffix: 'pt', info: l10n.treeBarWidthInfo),
          number(l10n.treeBarFont, _barFontPt, 0, 36, (v) => _barFontPt = v, step: 0.5, decimals: 1, suffix: 'pt', info: l10n.treeBarFontInfo),
          buttons(l10n.treeCladeColors, info: l10n.treeCladeColorsInfo, [
            ChipButton(label: l10n.treeColorClade, onPressed: _selected == null ? null : () => _colorClade(_selected!)),
            ChipButton(label: l10n.treeHighlightClade, onPressed: _selected == null ? null : () => _highlightClade(_selected!)),
            ChipButton(label: l10n.treeCladeLabel, onPressed: _selected == null ? null : () => _labelClade(_selected!)),
            ChipButton(label: l10n.treeClearCladeColors('${_cladeColors.length}'), onPressed: _cladeColors.isEmpty ? null : () => _edit(_cladeColors.clear)),
            ChipButton(label: l10n.treeClearHighlights('${_highlights.length}'), onPressed: _highlights.isEmpty ? null : () => _edit(_highlights.clear)),
            ChipButton(label: l10n.treeClearCladeLabels('${_cladeLabels.length}'), onPressed: _cladeLabels.isEmpty ? null : () => _edit(_cladeLabels.clear)),
          ]),
        ],
      _Tab.legend => [
          check(l10n.treeFigureLegend, l10n.treeLegendInFigure, _figLegend, () => _figLegend = !_figLegend, info: l10n.treeFigureLegendInfo),
          if (_figLegend) ...[
            select<LegendPosition>(
              l10n.treeLegendPosition,
              _legendPos,
              [(LegendPosition.below, l10n.treeLegendBelow), (LegendPosition.right, l10n.treeLegendRight)],
              (v) => _legendPos = v,
            ),
            check(l10n.treeLegendStyle, l10n.treeLegendBox, _legendBox, () => _legendBox = !_legendBox),
            text(l10n.treeLegendTitle, _legendTitle, (v) => _legendTitle = v, hint: l10n.treeLegendTitleHint),
            number(l10n.treeLegendFont, _legendFont, 4, 24, (v) => _legendFont = v, step: 0.5, decimals: 1, suffix: 'pt'),
            // The support part follows what the tree shows (classes and
            // "show support from"); its limits are set here too.
            if (_support != TreeSupport.none) ...[
              check(l10n.treeLegendSupport, l10n.treeLegendSupportShow, _supportInLegend, () => _supportInLegend = !_supportInLegend),
              if (_supportInLegend) ...[
                text(l10n.treeLegendSupportTitle, _supportLegendTitle, (v) => _supportLegendTitle = v, hint: l10n.treeSupportLegend),
                number(l10n.treeSupportMin, _supportMin, 0, 100, (v) => _supportMin = v, step: 5, info: l10n.treeSupportMinInfo),
                group(l10n.treeSupportClasses, _supportClassesEditor(context), info: l10n.treeSupportClassesInfo, wide: true),
              ],
            ],
          ],
        ],
      _Tab.export => [
          number(l10n.treeFigureWidth, _widthCm, 4, 60, (v) => _widthCm = v, step: 0.5, decimals: 1, suffix: 'cm', info: l10n.treeFigureWidthInfo),
          if (!_circular)
            number(l10n.treeBranchWidth, _treeWidthCm, 0, 80, (v) => _treeWidthCm = v, step: 0.5, decimals: 1, suffix: 'cm', info: l10n.treeBranchWidthInfo),
          select<double>(
            l10n.treeWidthPreset,
            widthPreset ?? -1,
            [(8.5, l10n.treeWidthOneCol), (12.0, l10n.treeWidthOneHalf), (17.5, l10n.treeWidthTwoCol)],
            (v) => _widthCm = v,
            display: widthPreset == null ? l10n.treeWidthCustom : null,
          ),
          number(l10n.treeMargin, _margin, 0, 72, (v) => _margin = v, step: 2, suffix: 'pt'),
          text(l10n.treeFigureTitle, _title, (v) => _title = v, hint: l10n.treeFigureTitleHint),
          number(l10n.treeTitleFont, _titleFont, 6, 36, (v) => _titleFont = v, step: 0.5, decimals: 1, suffix: 'pt'),
          check(l10n.treeScaleBar, l10n.treeShowScaleBar, _scaleBar, () => _scaleBar = !_scaleBar),
          if (_scaleBar)
            number(l10n.treeScaleBarLength, _barLength, 0, 100, (v) => _barLength = v, step: 0.01, decimals: 3, info: l10n.treeScaleBarLengthInfo),
          number(l10n.treeDpi, _dpi, 72, 1200, (v) => _dpi = v, step: 50, suffix: 'dpi', info: l10n.treeDpiInfo),
          check(l10n.treePngBackground, l10n.treeTransparent, _transparent, () => _transparent = !_transparent),
          buttons(l10n.treeExport, info: l10n.treeExportInfo, [
            ChipButton(label: 'Newick', onPressed: () => _save('nwk', toNewick(tree))),
            ChipButton(label: 'SVG', onPressed: () => _save('svg', _svg(tree))),
            ChipButton(label: 'PNG (${_dpi.round()} dpi)', style: ChipButtonStyle.accent, onPressed: _savePng),
          ]),
        ],
      _Tab.compare => _compareFields(context, group, select, check, number, color),
    };

  }

  static bool _sameFormat(List<LabelPart> a, List<LabelPart> b) {
    if (a.length != b.length) return false;
    for (var i = 0; i < a.length; i++) {
      if (a[i].field != b[i].field || a[i].italic != b[i].italic || a[i].bold != b[i].bold || a[i].parens != b[i].parens) return false;
    }
    return true;
  }

  String _presetLabel(AppLocalizations l10n, String key) => switch (key) {
        'original' => l10n.tipFieldOriginal,
        'species' => l10n.rankSpecies,
        'speciesVoucher' => '${l10n.rankSpecies} (${l10n.tipFieldVoucher})',
        'speciesAccession' => '${l10n.rankSpecies} (${l10n.tipFieldAccession})',
        'speciesLocality' => '${l10n.rankSpecies} (${l10n.tipFieldLocality})',
        _ => '${l10n.rankFamily} · ${l10n.rankSpecies}',
      };

  /// The parts of the label in order, each with italic / bold / brackets,
  /// and a menu to add more.
  /// Support classes: a preset (automatic for the kind of support, two to
  /// five classes, or cut at this tree's own tertiles) and one row per
  /// class with where it starts and its colour; classes can be added and
  /// removed.
  Widget _supportClassesEditor(BuildContext context) {
    final t = context.tokens;
    final l10n = context.l10n;
    final cls = _classes;
    final sorted = _sortedClasses;
    // Range of each threshold as the legend shows it.
    final ranges = {
      for (final (i, label) in supportClassRanges([for (final c in sorted) c.from.round()], _supportMin.round())) sorted[i].from.round(): label,
    };
    String presetName(SupportPreset p) => switch (p) {
          SupportPreset.auto => l10n.treeClassesAuto(widget.bootstrapScale ? 'bootstrap' : 'UFBoot / SH-aLRT'),
          SupportPreset.two => l10n.treeClassesN('2'),
          SupportPreset.three => l10n.treeClassesN('3'),
          SupportPreset.four => l10n.treeClassesN('4'),
          SupportPreset.five => l10n.treeClassesN('5'),
          SupportPreset.quantiles => l10n.treeClassesTertiles,
        };
    void edit(List<SupportClass> v) => setState(() => _setClasses(v));
    return Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
      Wrap(spacing: 6, runSpacing: 6, children: [
        for (final p in SupportPreset.values)
          ChipButton(
            label: presetName(p),
            style: _classesCustom == null && _supportPreset == p ? ChipButtonStyle.accent : ChipButtonStyle.neutral,
            onPressed: () => setState(() {
              _supportPreset = p;
              _classesCustom = null;
            }),
          ),
      ]),
      const SizedBox(height: 8),
      for (final (i, c) in cls.indexed)
        Padding(
          padding: const EdgeInsets.only(bottom: 6),
          child: Row(children: [
            _ColorButton(color: c.color, onPick: () async {
              final picked = await _pickColor(context, c.color);
              if (picked != null) edit([...cls]..[i] = (from: c.from, color: picked));
            }),
            const SizedBox(width: 8),
            Text('≥', style: TextStyle(fontSize: 13, color: t.textMuted)),
            const SizedBox(width: 6),
            NumberField(
              value: c.from,
              min: 0,
              max: 100,
              step: 1,
              width: 92,
              height: _fieldH,
              onChanged: (v) => edit([...cls]..[i] = (from: v, color: c.color)),
            ),
            const SizedBox(width: 10),
            Expanded(
              child: Text(ranges[c.from.round()] ?? l10n.treeClassHidden,
                  maxLines: 1,
                  overflow: TextOverflow.ellipsis,
                  style: TextStyle(fontSize: 12.5, color: ranges[c.from.round()] == null ? t.textFaint : t.text)),
            ),
            IconBtn(
              icon: FA.xmark,
              tooltip: l10n.treeClassRemove,
              size: 26,
              onPressed: cls.length <= 1 ? null : () => edit([...cls]..removeAt(i)),
            ),
          ]),
        ),
      Align(
        alignment: Alignment.centerLeft,
        child: ChipButton(
          label: l10n.treeClassAdd,
          icon: FA.plus,
          onPressed: () {
            // Halfway into the widest class, with the colour between its
            // neighbours.
            var best = 0, width = -1.0;
            for (var i = 0; i < sorted.length; i++) {
              final hi = i == 0 ? 100.0 : sorted[i - 1].from;
              if (hi - sorted[i].from > width) {
                width = hi - sorted[i].from;
                best = i;
              }
            }
            final hi = best == 0 ? 100.0 : sorted[best - 1].from;
            final from = ((hi + sorted[best].from) / 2).roundToDouble();
            final color = Color.lerp(best == 0 ? sorted[0].color : sorted[best - 1].color, sorted[best].color, 0.5)!;
            edit([...sorted, (from: from, color: color)]);
          },
        ),
      ),
    ]);
  }

  Widget _partsEditor(BuildContext context) {
    final t = context.tokens;
    final l10n = context.l10n;
    Widget style(String label, bool on, String tip, TextStyle textStyle, VoidCallback onTap) => SplaceTooltip(
          message: tip,
          child: InkWell(
            borderRadius: BorderRadius.circular(5),
            onTap: onTap,
            child: Container(
              width: 26,
              height: 22,
              alignment: Alignment.center,
              decoration: BoxDecoration(
                color: on ? t.accent : Colors.transparent,
                borderRadius: BorderRadius.circular(5),
                border: Border.all(color: on ? t.accent : t.border),
              ),
              child: Text(label, style: textStyle.copyWith(fontSize: 12, color: on ? Colors.white : t.textMuted)),
            ),
          ),
        );
    void update(int i, LabelPart p) => setState(() => _format = [..._format]..[i] = p);
    void moveLeft(int i) => setState(() {
          final f = [..._format];
          final p = f.removeAt(i);
          f.insert(i - 1, p);
          _format = f;
        });
    final used = {for (final p in _format) p.field};
    return Wrap(spacing: 8, runSpacing: 8, crossAxisAlignment: WrapCrossAlignment.center, children: [
      for (final (i, p) in _format.indexed)
        Container(
          padding: const EdgeInsets.fromLTRB(8, 3, 2, 3),
          decoration: BoxDecoration(color: t.surface, borderRadius: BorderRadius.circular(8), border: Border.all(color: t.border)),
          child: Row(mainAxisSize: MainAxisSize.min, children: [
            if (i > 0)
              SplaceTooltip(
                message: l10n.treeMovePart,
                child: InkWell(
                  onTap: () => moveLeft(i),
                  child: Padding(padding: const EdgeInsets.only(right: 6), child: Text('‹', style: TextStyle(fontSize: 15, color: t.textMuted))),
                ),
              ),
            Text(_fieldLabel(l10n, p.field), style: TextStyle(fontSize: 12.5, fontWeight: FontWeight.w600, color: t.text)),
            const SizedBox(width: 8),
            style('I', p.italic, l10n.treeItalic, const TextStyle(fontStyle: FontStyle.italic, fontWeight: FontWeight.w600), () => update(i, p.copyWith(italic: !p.italic))),
            const SizedBox(width: 3),
            style('B', p.bold, l10n.treeBold, const TextStyle(fontWeight: FontWeight.w800), () => update(i, p.copyWith(bold: !p.bold))),
            const SizedBox(width: 3),
            style('( )', p.parens, l10n.treeParens, const TextStyle(fontWeight: FontWeight.w600), () => update(i, p.copyWith(parens: !p.parens))),
            IconButton(
              tooltip: l10n.treeRemovePart,
              visualDensity: VisualDensity.compact,
              iconSize: 11,
              onPressed: _format.length == 1 ? null : () => setState(() => _format = [..._format]..removeAt(i)),
              icon: const Icon(FA.xmark),
            ),
          ]),
        ),
      if (used.length < TipField.values.length)
        PopupMenuButton<TipField>(
          tooltip: l10n.treeAddPart,
          onSelected: (f) => setState(() => _format = [..._format, LabelPart(f, italic: f == TipField.species || f == TipField.genus)]),
          itemBuilder: (context) => [
            for (final f in TipField.values)
              if (!used.contains(f)) PopupMenuItem(value: f, child: Text(_fieldLabel(l10n, f))),
          ],
          child: Container(
            padding: const EdgeInsets.symmetric(horizontal: 10, vertical: 6),
            decoration: BoxDecoration(borderRadius: BorderRadius.circular(8), border: Border.all(color: t.accentSoftBorder), color: t.accentSoft),
            child: Row(mainAxisSize: MainAxisSize.min, children: [
              Icon(FA.plus, size: 10, color: t.accentText),
              const SizedBox(width: 6),
              Text(l10n.treeAddPart, style: TextStyle(fontSize: 12.5, fontWeight: FontWeight.w600, color: t.accentText)),
            ]),
          ),
        ),
    ]);
  }

  Widget _selectionBar(BuildContext context, TreeNode tree, TreeNode sel, {bool scroll = false}) {
    final t = context.tokens;
    final l10n = context.l10n;
    final truncated = _truncated.contains(sel);
    return Container(
      padding: const EdgeInsets.symmetric(horizontal: 12, vertical: 8),
      decoration: BoxDecoration(color: t.accentSoft, borderRadius: BorderRadius.circular(10), border: Border.all(color: t.accentSoftBorder), boxShadow: shadowMd),
      child: _bar(scroll, [
        Row(mainAxisSize: MainAxisSize.min, children: [
          Icon(FA.circleCheck, size: 13, color: t.accent),
          const SizedBox(width: 8),
          ConstrainedBox(
            constraints: const BoxConstraints(maxWidth: 420),
            child: Text(
              sel.isTip ? _plain(sel.name) : '${l10n.treeClade(sel.tipCount.toString())}${sel.label == null ? '' : ' · ${l10n.treeSupportValue(sel.label!)}'}',
              overflow: TextOverflow.ellipsis,
              style: TextStyle(fontSize: 13, fontWeight: FontWeight.w600, color: t.accentText),
            ),
          ),
          const SizedBox(width: 8),
        ]),
        ChipButton(label: l10n.treeRootHere, style: ChipButtonStyle.accent, onPressed: identical(sel, tree) ? null : _rootOnSelected),
        if (sel.parent != null) ChipButton(label: l10n.treeSelectParent, onPressed: () => setState(() => _selected = sel.parent)),
        if (!sel.isTip) ChipButton(label: l10n.treeRotateClade, onPressed: _rotateSelected),
        if (!sel.isTip) ChipButton(label: _collapsed.contains(sel) ? l10n.treeExpand : l10n.treeCollapse, onPressed: () => _toggleCollapse(sel)),
        ChipButton(label: l10n.treeColorClade, onPressed: () => _colorClade(sel)),
        if (_cladeColors.containsKey(sel)) ChipButton(label: l10n.treeRemoveColor, onPressed: () => _edit(() => _cladeColors.remove(sel))),
        ChipButton(label: l10n.treeHighlightClade, onPressed: () => _highlightClade(sel)),
        if (_highlights.containsKey(sel)) ChipButton(label: l10n.treeRemoveHighlight, onPressed: () => _edit(() => _highlights.remove(sel))),
        ChipButton(label: l10n.treeCladeLabel, onPressed: () => _labelClade(sel)),
        if (sel.isTip) ChipButton(label: l10n.treeRenameTip, onPressed: () => _renameTip(sel)),
        if (widget.onIdentify != null) ChipButton(label: l10n.treeIdentify, icon: FA.magnifyingGlass, onPressed: () => _identify(sel)),
        if (sel.parent != null && (sel.length ?? 0) > 0)
          ChipButton(label: truncated ? l10n.treeRestoreBranch : l10n.treeShortenBranch, onPressed: () => _toggleTruncate(sel)),
        if (sel.parent != null && tree.tipCount - sel.tipCount >= 2)
          ChipButton(label: l10n.treeRemoveFromFigure, style: ChipButtonStyle.danger, onPressed: () => _prune(sel)),
        ChipButton(label: l10n.treeClearSelection, onPressed: () => setState(() => _selected = null)),
      ]),
    );
  }

  /// Buttons in one scrolling row ([scroll]) or wrapped in lines.
  static Widget _bar(bool scroll, List<Widget> children) => scroll
      ? SingleChildScrollView(
          scrollDirection: Axis.horizontal,
          child: Row(children: [for (final (i, w) in children.indexed) ...[if (i > 0) const SizedBox(width: 8), w]]),
        )
      : Wrap(spacing: 8, runSpacing: 8, crossAxisAlignment: WrapCrossAlignment.center, children: children);

  Widget _statsRow(BuildContext context, TreeNode tree, Size figure) {
    final t = context.tokens;
    final l10n = context.l10n;
    final total = tree.nodes.fold<double>(0, (a, n) => a + (n.length ?? 0));
    final (pw, ph) = _pngSize(figure);
    return Wrap(spacing: 14, runSpacing: 6, crossAxisAlignment: WrapCrossAlignment.center, children: [
      Text(l10n.treeStats('${tree.tipCount}', total.toStringAsFixed(3)), style: TextStyle(fontSize: 12, color: t.textMuted)),
      Text(
        l10n.treeFigureSize((figure.width / _pxPerCm).toStringAsFixed(1), (figure.height / _pxPerCm).toStringAsFixed(1), '$pw', '$ph', '${_dpi.round()}'),
        style: TextStyle(fontSize: 12, color: t.textMuted),
      ),
      if (_rooting == TreeRooting.manual) TagChip(l10n.treeRootManual, color: t.accent),
      SplaceTooltip(
        message: widget.bootstrapScale ? l10n.treeSupportInfoBootstrap : l10n.treeSupportInfoUfboot,
        child: Row(mainAxisSize: MainAxisSize.min, children: [
          Text(l10n.treeSupportLegend, style: TextStyle(fontSize: 12, color: t.textMuted)),
          const SizedBox(width: 4),
          Icon(FA.circleInfo, size: 12, color: t.textFaint),
        ]),
      ),
    ]);
  }

  /// Info bar text: the hovered tip or clade, or the hint.
  String _describe(TreeNode? n, AppLocalizations l10n) {
    if (n == null) return l10n.treeHint;
    var length = n.length == null ? '—' : n.length!.toStringAsFixed(5);
    if (_truncated.contains(n)) length = '$length (${l10n.treeShortenedLegend})';
    if (n.isTip) {
      final i = _info(n.name);
      final extra = [
        if (i.family.isNotEmpty) i.family,
        if (i.accession.isNotEmpty) i.accession,
        if (i.voucher.isNotEmpty) i.voucher,
        if (i.locality.isNotEmpty) i.locality,
        if (i.latLon.isNotEmpty) i.latLon,
      ].join(' · ');
      return '${l10n.treeHoverTip(_display(n.name), length)}${extra.isEmpty ? '' : '  ·  $extra'}';
    }
    return l10n.treeHoverClade('${n.tipCount}', n.label ?? '—', length);
  }

  Future<void> _save(String ext, String content, {String? name}) async {
    final loc = await getSaveLocation(suggestedName: '${name ?? widget.fileName}.$ext');
    if (loc == null) return;
    await File(loc.path).writeAsString(content);
  }

  /// Pixels per figure pixel (1/96 in) of the PNG: the chosen DPI, with the
  /// longest side kept within 16 000 px.
  double _pngRatio(Size figure) {
    final ratio = _dpi / 96;
    final longest = math.max(figure.width, figure.height);
    return longest * ratio > 16000 ? 16000 / longest : ratio;
  }

  (int, int) _pngSize(Size figure) {
    final r = _pngRatio(figure);
    return ((figure.width * r).ceil(), (figure.height * r).ceil());
  }

  /// PNG of the figure at the chosen DPI, drawn off screen (no selection
  /// or search marks), with the DPI written in the file. Its size in
  /// centimetres is the figure's, whatever the screen.
  Future<void> _savePng() async {
    final tree = _tree;
    if (tree == null) return;
    final (groupOf, groupColors) = _groups(_rows(tree));
    final circle = _circular ? _circle(tree, groupOf, groupColors) : null;
    final layout = circle == null ? _layout(tree, _width, groupOf, groupColors) : null;
    final size = circle?.size ?? layout!.size;
    final ratio = _pngRatio(size);
    final recorder = ui.PictureRecorder();
    final canvas = Canvas(recorder)..scale(ratio);
    if (circle != null) {
      _drawCircle(_CanvasPen(canvas), circle, groupOf, groupColors, export: true);
    } else {
      _painter(layout!, groupOf, groupColors, export: true).paint(canvas, size);
    }
    final image = await recorder.endRecording().toImage((size.width * ratio).ceil(), (size.height * ratio).ceil());
    final bytes = await image.toByteData(format: ui.ImageByteFormat.png);
    image.dispose();
    if (bytes == null) return;
    final loc = await getSaveLocation(suggestedName: '${widget.fileName}.png');
    if (loc == null) return;
    await File(loc.path).writeAsBytes(pngWithDpi(bytes.buffer.asUint8List(), _dpi * ratio / (_dpi / 96)));
  }

  /// SVG file of [size] (its real size in centimetres) around [body].
  String _svgDocument(Size size, String body) {
    final family = '${_SvgPen.esc(_family)}, Arial, sans-serif';
    return '<svg xmlns="http://www.w3.org/2000/svg" width="${(size.width / _pxPerCm).toStringAsFixed(2)}cm" '
        'height="${(size.height / _pxPerCm).toStringAsFixed(2)}cm" viewBox="0 0 ${_SvgPen.f(size.width)} ${_SvgPen.f(size.height)}" '
        'font-family="$family">\n$body</svg>\n';
  }

  /// The two trees face to face, as SVG or as PNG at the chosen DPI.
  Future<void> _saveTangle({required bool svg}) async {
    final tree = _tree, other = _otherTree;
    if (tree == null || other == null || widget.compareWith.isEmpty) return;
    final l = _tangle(tree, _facing(other));
    final right = widget.compareWith[_otherIndex.clamp(0, widget.compareWith.length - 1)].label;
    final name = '${widget.fileName}_compare';
    if (svg) {
      final pen = _SvgPen();
      _drawTangle(pen, l, right, export: true);
      await _save('svg', _svgDocument(l.size, pen.b.toString()), name: name);
      return;
    }
    final ratio = _pngRatio(l.size);
    final recorder = ui.PictureRecorder();
    final canvas = Canvas(recorder)..scale(ratio);
    _drawTangle(_CanvasPen(canvas), l, right, export: true);
    final image = await recorder.endRecording().toImage((l.size.width * ratio).ceil(), (l.size.height * ratio).ceil());
    final bytes = await image.toByteData(format: ui.ImageByteFormat.png);
    image.dispose();
    if (bytes == null) return;
    final loc = await getSaveLocation(suggestedName: '$name.png');
    if (loc == null) return;
    await File(loc.path).writeAsBytes(pngWithDpi(bytes.buffer.asUint8List(), _dpi * ratio / (_dpi / 96)));
  }

  /// The same figure as SVG (vector, for journals and editing).
  String _svg(TreeNode tree) {
    final (groupOf, groupColors) = _groups(_rows(tree));
    if (_circular) {
      final c = _circle(tree, groupOf, groupColors);
      final pen = _SvgPen();
      _drawCircle(pen, c, groupOf, groupColors, export: true);
      return _svgDocument(c.size, pen.b.toString());
    }
    final layout = _layout(tree, _width, groupOf, groupColors);
    final st = _style;
    final colors = _branchColors(tree, groupOf, groupColors);
    String hex(Color c) => '#${(c.toARGB32() & 0xFFFFFF).toRadixString(16).padLeft(6, '0')}';
    String op(Color c) => c.a.toStringAsFixed(2);
    String esc(String s) => s.replaceAll('&', '&amp;').replaceAll('<', '&lt;').replaceAll('>', '&gt;');
    String f(double v) => v.toStringAsFixed(1);
    final family = '${esc(st.family)}, Arial, sans-serif';
    // Size in centimetres (the figure's real size), drawn in 1/96 in units.
    final cm = (layout.size.width / _pxPerCm, layout.size.height / _pxPerCm);
    final b = StringBuffer()
      ..writeln('<svg xmlns="http://www.w3.org/2000/svg" width="${cm.$1.toStringAsFixed(2)}cm" height="${cm.$2.toStringAsFixed(2)}cm" '
          'viewBox="0 0 ${f(layout.size.width)} ${f(layout.size.height)}" font-family="$family">');
    if (!_transparent) b.writeln('<rect width="100%" height="100%" fill="${hex(st.background)}"/>');
    // Title.
    if (st.title.trim().isNotEmpty) {
      b.writeln('<text x="${f(layout.left)}" y="${f(st.margin + st.titleFont)}" font-size="${f(st.titleFont)}" font-weight="bold" fill="${hex(st.label)}">${esc(st.title)}</text>');
    }
    // Row backgrounds by taxon, named on the right.
    final groups = layout.groups(groupOf);
    if (st.groupRows) {
      for (final (name, y0, y1) in groups) {
        final c = groupColors[name]!;
        b
          ..writeln('<rect x="0" y="${f(y0)}" width="${f(layout.groupRight)}" height="${f(y1 - y0)}" fill="${hex(c)}" fill-opacity="0.2"/>')
          ..writeln('<text x="${f(layout.groupRight - 8)}" y="${f((y0 + y1) / 2 + st.font * 0.35)}" font-size="${st.font}" font-weight="bold" '
              'text-anchor="end" fill="${hex(c)}">${esc(name)}</text>');
      }
    }
    // Coloured boxes behind clades.
    for (final e in _highlights.entries.toList()..sort((a, b) => b.key.tipCount.compareTo(a.key.tipCount))) {
      final sp = layout.span(e.key);
      if (sp == null) continue;
      final x0 = (e.key.parent == null ? layout.point(e.key).dx : (layout.point(e.key).dx + layout.point(e.key.parent!).dx) / 2) - 2;
      b.writeln('<rect x="${f(x0)}" y="${f(sp.$1 - _rowH / 2 + 1)}" width="${f(layout.groupRight - 4 - x0)}" height="${f(sp.$2 - sp.$1 + _rowH - 2)}" rx="5" '
          'fill="${hex(e.value)}" fill-opacity="${st.highlightAlpha.toStringAsFixed(2)}"/>');
    }
    // Branches, each in its colour.
    final cap = st.branchStyle == BranchStyle.square ? 'square' : 'round';
    b.writeln('<g stroke-width="${st.lineW.toStringAsFixed(2)}" fill="none" stroke-linecap="$cap" stroke-linejoin="round">');
    for (final (cmds, n) in layout.branchPaths(st.branchStyle)) {
      final d = [for (final (op, pts) in cmds) '$op${pts.map((p) => '${f(p.dx)},${f(p.dy)}').join(' ')}'].join(' ');
      b.writeln('<path d="$d" stroke="${hex(colors[n] ?? st.branch)}"/>');
    }
    b.writeln('</g>');
    // Shortened branches: // on the branch.
    if (_shape == TreeShape.phylogram) {
      for (final n in _truncated) {
        if (layout.hidden(n) || n.parent == null) continue;
        final m = layout.breakPoint(n, st.branchStyle);
        var gap = st.background;
        if (st.groupRows) {
          for (final (name, y0, y1) in groups) {
            if (m.dy >= y0 && m.dy < y1) gap = Color.alphaBlend(groupColors[name]!.withValues(alpha: 0.2), gap);
          }
        }
        b.writeln('<rect x="${f(m.dx - 4)}" y="${f(m.dy - 6)}" width="8" height="12" fill="${hex(gap)}"/>');
        for (final dx in [-2.5, 2.5]) {
          b.writeln('<line x1="${f(m.dx + dx - 2)}" y1="${f(m.dy + 5)}" x2="${f(m.dx + dx + 2)}" y2="${f(m.dy - 5)}" stroke="${hex(st.branch)}" stroke-width="1.3"/>');
        }
      }
    }
    // Branch lengths.
    if (st.lengths) {
      for (final n in tree.nodes) {
        if (n.parent == null || n.length == null || layout.hidden(n)) continue;
        final m = layout.breakPoint(n, st.branchStyle);
        b.writeln('<text x="${f(m.dx)}" y="${f(m.dy - st.lineW / 2 - 2)}" font-size="${math.max(7.0, st.supportFont - 1).toStringAsFixed(1)}" '
            'text-anchor="middle" fill="${hex(_Style.muted)}">${_lengthText(n.length!)}</text>');
      }
    }
    // Tip labels.
    final out = {for (final o in widget.outgroup) o.toLowerCase().replaceAll('_', ' ').trim()};
    final end = layout.left + layout.maxX * layout.scale;
    for (final n in layout.rows) {
      final o = layout.point(n);
      final collapsed = _collapsed.contains(n);
      if (collapsed) {
        b.writeln('<path d="M${f(o.dx)},${f(o.dy)} L${f(end)},${f(o.dy - _rowH * 0.42)} L${f(end)},${f(o.dy + _rowH * 0.42)} Z" '
            'fill="${hex(colors[n] ?? st.branch)}" fill-opacity="0.25" stroke="${hex(colors[n] ?? st.branch)}"/>');
      }
      if (st.tipDots && n.isTip) b.writeln('<circle cx="${f(o.dx)}" cy="${f(o.dy)}" r="${f(st.dot)}" fill="${hex(colors[n] ?? st.branch)}"/>');
      if (!st.showLabels) continue;
      var x = (collapsed || st.align ? end : o.dx) + st.labelGap;
      if (st.align && !collapsed && end - o.dx > 6) {
        b.writeln('<line x1="${f(o.dx + 4)}" y1="${f(o.dy)}" x2="${f(end + st.labelGap - 3)}" y2="${f(o.dy)}" stroke="${hex(_Style.guide)}" stroke-width="1" stroke-dasharray="1.5,2.5"/>');
      }
      if (st.tipDots && n.isTip && !st.align) x += st.dot;
      final isOut = st.markOutgroup && n.isTip && out.contains(n.name.toLowerCase().replaceAll('_', ' ').trim());
      final labelColor = st.groupLabels && groupOf[n] != null ? _darker(groupColors[groupOf[n]]!) : st.label;
      final spans = collapsed
          ? esc(_rowText(n))
          : [
              for (final (text, p) in _parts(n.name))
                '<tspan${p.italic ? ' font-style="italic"' : ''}${p.bold || isOut ? ' font-weight="bold"' : ''}>${esc(text)}</tspan>',
            ].join(' ');
      if (isOut) {
        final w = TextPainter(text: _rowSpan(n, TextStyle(fontFamily: st.family, fontSize: st.font), search: false), textDirection: TextDirection.ltr)..layout();
        b.writeln('<rect x="${f(x - 3)}" y="${f(o.dy - w.height / 2 - 1)}" width="${f(w.width + 6)}" height="${f(w.height + 2)}" rx="4" fill="${hex(_Style.outgroup)}"/>');
      }
      b.writeln('<text x="${f(x)}" y="${f(o.dy + st.font * 0.35)}" font-size="${st.font}" fill="${hex(labelColor)}" xml:space="preserve">$spans</text>');
    }
    // Bars with names beside clades.
    for (final bar in layout.bars) {
      final sp = layout.span(bar.node);
      if (sp == null) continue;
      final y0 = sp.$1 - _rowH * 0.38, y1 = sp.$2 + _rowH * 0.38;
      final c = hex(bar.mark.color);
      b
        ..writeln('<line x1="${f(bar.x)}" y1="${f(y0)}" x2="${f(bar.x)}" y2="${f(y1)}" stroke="$c" stroke-width="${f(st.barStroke)}"/>')
        ..writeln('<text x="${f(bar.x + 6)}" y="${f((y0 + y1) / 2 + st.barText * 0.35)}" font-size="${f(st.barText)}" font-weight="bold"'
            '${bar.mark.italic ? ' font-style="italic"' : ''} fill="$c">${esc(bar.mark.text)}</text>');
    }
    {
      final pen = _SvgPen();
      _drawHeat(pen, layout, st);
      b.write(pen.b);
    }
    // Support.
    if (st.support != TreeSupport.none) {
      final items = <_SupportItem>[];
      for (final n in tree.nodes) {
        if (n.isTip || n.parent == null || n.label == null || layout.hidden(n)) continue;
        final v = n.support ?? 0;
        if (v < st.supportMin) continue;
        items.add((o: layout.point(n), text: n.label!, r: st.radius(n.support), color: st.supportColor(n.support)));
      }
      if (st.support == TreeSupport.numbers) {
        final pen = _SvgPen();
        _drawSupports(pen, st, items);
        b.write(pen.b);
      }
      for (final it in items) {
        final pen = _SvgPen();
        _drawMark(pen, it.o, it.r, it.color, st.markShape);
        b.write(pen.b);
      }
    }
    // Scale bar.
    final bar = layout.scaleBar();
    if (bar != null) {
      final c = hex(st.branch);
      final w = f(st.lineW);
      b
        ..writeln('<line x1="${f(bar.$1.dx)}" y1="${f(bar.$1.dy)}" x2="${f(bar.$2.dx)}" y2="${f(bar.$2.dy)}" stroke="$c" stroke-width="$w"/>')
        ..writeln('<line x1="${f(bar.$1.dx)}" y1="${f(bar.$1.dy - 4)}" x2="${f(bar.$1.dx)}" y2="${f(bar.$1.dy + 4)}" stroke="$c" stroke-width="$w"/>')
        ..writeln('<line x1="${f(bar.$2.dx)}" y1="${f(bar.$2.dy - 4)}" x2="${f(bar.$2.dx)}" y2="${f(bar.$2.dy + 4)}" stroke="$c" stroke-width="$w"/>')
        ..writeln('<text x="${f((bar.$1.dx + bar.$2.dx) / 2)}" y="${f(bar.$1.dy + 6 + st.legendFont)}" font-size="${f(st.legendFont)}" text-anchor="middle" '
            'fill="${hex(_Style.muted)}">${bar.$3}</text>');
    }
    // Legend.
    final lg = layout.legendBox;
    if (lg != null) {
      final o = layout.legendOrigin;
      if (st.legendBox) {
        b.writeln('<rect x="${f(o.dx)}" y="${f(o.dy)}" width="${f(lg.size.width)}" height="${f(lg.size.height)}" rx="6" '
            'fill="${hex(st.background)}" stroke="#d1d5db"/>');
      }
      if (lg.title.isNotEmpty) {
        b.writeln('<text x="${f(o.dx + lg.pad)}" y="${f(o.dy + lg.pad + st.legendFont)}" font-size="${f(st.legendFont + 1)}" font-weight="bold" '
            'fill="${hex(st.label)}">${esc(lg.title)}</text>');
      }
      for (final (it, at) in lg.items) {
        final p = o + at;
        var x = p.dx;
        final mid = p.dy + st.legendFont * 0.65;
        if (it.symbol != null) {
          b.writeln('<text x="${f(x)}" y="${f(p.dy + st.legendFont)}" font-size="${st.legendFont}" font-weight="bold" fill="${hex(st.branch)}">${esc(it.symbol!)}</text>');
          x += st.legendFont * 1.5;
        } else if (it.color != null) {
          final c = it.color!;
          b.writeln(it.box
              ? '<rect x="${f(x)}" y="${f(mid - 5)}" width="13" height="10" rx="2" fill="${hex(c)}" fill-opacity="${op(c)}"/>'
              : '<circle cx="${f(x + 5)}" cy="${f(mid)}" r="4.5" fill="${hex(c)}"/>');
          x += 17;
        }
        b.writeln('<text x="${f(x)}" y="${f(p.dy + st.legendFont)}" font-size="${st.legendFont}"${it.heading ? ' font-weight="bold"' : ''}'
            '${it.italic ? ' font-style="italic"' : ''} fill="${hex(it.heading ? st.label : const Color(0xFF4B5563))}">${esc(it.label)}</text>');
      }
    }
    b.writeln('</svg>');
    return b.toString();
  }
}

Color _darker(Color c) => Color.lerp(c, const Color(0xFF000000), 0.25)!;

// ---------------------------------------------------------------------------
// Colour picker
// ---------------------------------------------------------------------------

class _ColorButton extends StatelessWidget {
  const _ColorButton({required this.color, required this.onPick});
  final Color color;
  final VoidCallback onPick;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    return InkWell(
      borderRadius: BorderRadius.circular(8),
      onTap: onPick,
      child: Container(
        width: _fieldW,
        height: _fieldH,
        padding: const EdgeInsets.symmetric(horizontal: 8),
        decoration: BoxDecoration(color: t.surface, borderRadius: BorderRadius.circular(8), border: Border.all(color: t.border)),
        child: Row(children: [
          Container(
            width: 28,
            height: 20,
            decoration: BoxDecoration(color: color, borderRadius: BorderRadius.circular(4), border: Border.all(color: t.border)),
          ),
          const SizedBox(width: 10),
          Expanded(child: Text(_hex(color), style: TextStyle(fontFamily: monoFont, fontSize: 12.5, color: t.text))),
          Icon(FA.palette, size: 11, color: t.textFaint),
        ]),
      ),
    );
  }
}

String _hex(Color c) => '#${(c.toARGB32() & 0xFFFFFF).toRadixString(16).padLeft(6, '0').toUpperCase()}';

/// Swatches of the app palette and a hex field.
Future<Color?> _pickColor(BuildContext context, Color current) {
  final hex = TextEditingController(text: _hex(current));
  return showDialog<Color>(
    context: context,
    builder: (context) {
      final t = context.tokens;
      final l10n = context.l10n;
      return AlertDialog(
        title: Text(l10n.treePickColor, style: const TextStyle(fontSize: 16)),
        content: SizedBox(
          width: 300,
          child: Column(mainAxisSize: MainAxisSize.min, crossAxisAlignment: CrossAxisAlignment.stretch, children: [
            Wrap(spacing: 8, runSpacing: 8, children: [
              for (final c in _swatches)
                InkWell(
                  borderRadius: BorderRadius.circular(6),
                  onTap: () => Navigator.of(context).pop(c),
                  child: Container(
                    width: 38,
                    height: 30,
                    decoration: BoxDecoration(
                      color: c,
                      borderRadius: BorderRadius.circular(6),
                      border: Border.all(color: c.toARGB32() == current.toARGB32() ? t.accent : t.border, width: c.toARGB32() == current.toARGB32() ? 2.5 : 1),
                    ),
                  ),
                ),
            ]),
            const SizedBox(height: 14),
            TextField(
              controller: hex,
              style: const TextStyle(fontFamily: monoFont, fontSize: 13),
              decoration: const InputDecoration(isDense: true, labelText: 'HEX', hintText: '#3A47AB'),
              onSubmitted: (v) {
                final c = _parseHex(v);
                if (c != null) Navigator.of(context).pop(c);
              },
            ),
          ]),
        ),
        actions: [
          TextButton(onPressed: () => Navigator.of(context).pop(), child: Text(l10n.cancel)),
          FilledButton(
            onPressed: () {
              final c = _parseHex(hex.text);
              if (c != null) Navigator.of(context).pop(c);
            },
            child: const Text('OK'),
          ),
        ],
      );
    },
  );
}

/// Text with italic / bold and, when [color] is given, a colour: the name
/// of a clade bar or a tip's own label. An empty text removes it.
Future<({String text, bool italic, bool bold, Color? color})?> _textDialog(
  BuildContext context, {
  required String title,
  required String initial,
  bool italic = false,
  bool bold = false,
  Color? color,
}) {
  final ctrl = TextEditingController(text: initial);
  var it = italic, bd = bold;
  var col = color;
  return showDialog(
    context: context,
    builder: (context) => StatefulBuilder(builder: (context, setLocal) {
      final t = context.tokens;
      final l10n = context.l10n;
      Widget flag(String label, bool on, TextStyle style, VoidCallback onTap) => InkWell(
            borderRadius: BorderRadius.circular(8),
            onTap: () => setLocal(onTap),
            child: Container(
              width: 36,
              height: 36,
              alignment: Alignment.center,
              decoration: BoxDecoration(
                color: on ? t.accent : t.surface,
                borderRadius: BorderRadius.circular(8),
                border: Border.all(color: on ? t.accent : t.border),
              ),
              child: Text(label, style: style.copyWith(fontSize: 14, color: on ? Colors.white : t.textMuted)),
            ),
          );
      void done() => Navigator.of(context).pop((text: ctrl.text, italic: it, bold: bd, color: col));
      return AlertDialog(
        title: Text(title, style: const TextStyle(fontSize: 16)),
        content: SizedBox(
          width: 360,
          child: Column(mainAxisSize: MainAxisSize.min, crossAxisAlignment: CrossAxisAlignment.stretch, children: [
            Row(children: [
              Expanded(
                child: TextField(
                  controller: ctrl,
                  autofocus: true,
                  style: TextStyle(fontSize: 14, fontStyle: it ? FontStyle.italic : null, fontWeight: bd ? FontWeight.w700 : null),
                  decoration: InputDecoration(isDense: true, hintText: l10n.treeTextHint),
                  onSubmitted: (_) => done(),
                ),
              ),
              const SizedBox(width: 8),
              SplaceTooltip(message: l10n.treeItalic, child: flag('I', it, const TextStyle(fontStyle: FontStyle.italic, fontWeight: FontWeight.w600), () => it = !it)),
              if (color == null) ...[
                const SizedBox(width: 6),
                SplaceTooltip(message: l10n.treeBold, child: flag('B', bd, const TextStyle(fontWeight: FontWeight.w800), () => bd = !bd)),
              ],
            ]),
            if (col != null) ...[
              const SizedBox(height: 12),
              _ColorButton(
                color: col!,
                onPick: () async {
                  final c = await _pickColor(context, col!);
                  if (c != null) setLocal(() => col = c);
                },
              ),
            ],
            const SizedBox(height: 8),
            Text(l10n.treeTextEmptyHint, style: TextStyle(fontSize: 12, color: t.textMuted)),
          ]),
        ),
        actions: [
          TextButton(onPressed: () => Navigator.of(context).pop(), child: Text(l10n.cancel)),
          FilledButton(onPressed: done, child: const Text('OK')),
        ],
      );
    }),
  );
}

Color? _parseHex(String s) {
  final h = s.trim().replaceFirst('#', '');
  if (!RegExp(r'^[0-9a-fA-F]{6}$').hasMatch(h)) return null;
  return Color(0xFF000000 | int.parse(h, radix: 16));
}

// ---------------------------------------------------------------------------
// PNG resolution
// ---------------------------------------------------------------------------

/// [png] with a pHYs chunk saying [dpi] (so figures open at the right size
/// in editors and journals' systems).
Uint8List pngWithDpi(Uint8List png, double dpi) {
  // Signature (8) + IHDR chunk (4 length + 4 type + 13 data + 4 CRC).
  const ihdrEnd = 8 + 25;
  if (png.length < ihdrEnd) return png;
  final ppm = (dpi / 0.0254).round();
  final data = ByteData(9)
    ..setUint32(0, ppm)
    ..setUint32(4, ppm)
    ..setUint8(8, 1);
  final type = 'pHYs'.codeUnits;
  final body = Uint8List.fromList([...type, ...data.buffer.asUint8List()]);
  final chunk = ByteData(4 + body.length + 4)..setUint32(0, 9);
  final bytes = chunk.buffer.asUint8List()..setRange(4, 4 + body.length, body);
  chunk.setUint32(4 + body.length, _crc32(body));
  return Uint8List.fromList([...png.sublist(0, ihdrEnd), ...bytes, ...png.sublist(ihdrEnd)]);
}

final _crcTable = List<int>.generate(256, (n) {
  var c = n;
  for (var k = 0; k < 8; k++) {
    c = (c & 1) != 0 ? 0xEDB88320 ^ (c >> 1) : c >> 1;
  }
  return c;
});

int _crc32(List<int> bytes) {
  var c = 0xFFFFFFFF;
  for (final b in bytes) {
    c = _crcTable[(c ^ b) & 0xFF] ^ (c >> 8);
  }
  return (c ^ 0xFFFFFFFF) & 0xFFFFFFFF;
}

// ---------------------------------------------------------------------------
// Layout
// ---------------------------------------------------------------------------

/// Positions of the nodes: x from the drawn branch lengths (phylogram) or
/// depth (cladogram), y one row per visible tip or collapsed clade; the
/// title above, the row-group names, the scale bar and the legend.
class _Layout {
  _Layout._(this.root, this.rows, this._x, this._y, this.size, this.left, this.top, this.scale, this.collapsed, this.maxX, this.shape, this.rowH,
      this.groupW, this.groupRight, this.barLength, this.legendBox, this.legendOrigin, this.bars, this.heat, this.heatX);

  /// SALVE heatmap beside the labels, from [heatX].
  final _Heat? heat;
  final double heatX;

  final TreeNode root;
  final List<TreeNode> rows;
  final Map<TreeNode, double> _x, _y;
  final Size size;
  final double left, top, scale, maxX, rowH, groupW;

  /// Right edge of the row bands and their names (left of a side legend).
  final double groupRight;
  final Set<TreeNode> collapsed;
  final TreeShape shape;

  /// Scale bar length (0 = automatic, negative = hidden).
  final double barLength;
  final _LegendBox? legendBox;
  final Offset legendOrigin;

  /// Bars beside clades with their x (nested ones in columns further out).
  final List<({TreeNode node, _CladeMark mark, double x})> bars;

  /// Bottom of the last row.
  double get rowsBottom => top + rows.length * rowH;

  static _Layout compute(TreeNode root, Set<TreeNode> collapsed, TreeShape shape, double Function(TreeNode) labelWidth, double width,
      double rowH, double Function(TreeNode) length, double groupW,
      {required _Style style,
      double xScale = 1,
      double barLength = 0,
      List<List<_LegendItem>> legend = const [],
      LegendPosition legendPos = LegendPosition.below,
      List<(TreeNode, _CladeMark)> bars = const [],
      _Heat? heat,
      double treeW = 0}) {
    final rows = <TreeNode>[];
    final x = <TreeNode, double>{}, y = <TreeNode, double>{};
    int height(TreeNode n) => n.isTip || collapsed.contains(n) ? 0 : 1 + n.children.map(height).reduce(math.max);
    final depth = height(root);
    double place(TreeNode n, double px) {
      final nx = shape == TreeShape.phylogram ? px + (n.parent == null ? 0 : length(n)) : (depth - height(n)).toDouble();
      x[n] = nx;
      if (n.isTip || collapsed.contains(n)) {
        rows.add(n);
        y[n] = (rows.length - 1).toDouble();
      } else {
        final ys = [for (final c in n.children) place(c, nx)];
        y[n] = (ys.first + ys.last) / 2;
      }
      return y[n]!;
    }

    place(root, 0);
    var maxX = x.values.fold<double>(0, math.max);
    // Collapsed clades reach their deepest tip.
    for (final n in collapsed) {
      if (!x.containsKey(n)) continue;
      double deepest(TreeNode k, double acc) =>
          k.isTip ? acc : k.children.map((c) => deepest(c, acc + (shape == TreeShape.phylogram ? length(c) : 1))).reduce(math.max);
      maxX = math.max(maxX, x[n]! + deepest(n, 0));
    }
    if (maxX <= 0) maxX = 1;
    var labelW = 0.0;
    for (final n in rows) {
      labelW = math.max(labelW, labelWidth(n));
    }
    labelW = math.min(labelW + 8, width * 0.5);
    final left = style.margin + 4;
    final heatW = heat?.width ?? 0;
    final top = style.margin + 6 + style.titleH + (heat?.headerH ?? 0);
    final side = legend.isNotEmpty && legendPos == LegendPosition.right;
    final sideBox = side ? _LegendBox.compute(legend, style, vertical: true, maxWidth: width * 0.35) : null;
    final sideW = sideBox == null ? 0.0 : sideBox.size.width + 24;
    // Bars beside clades: a clade holding other labelled clades goes one
    // column further out; each column is as wide as its longest name.
    final shown = [for (final b in bars) if (x.containsKey(b.$1)) b];
    final level = <TreeNode, int>{};
    for (final b in [...shown]..sort((a, b) => a.$1.tipCount.compareTo(b.$1.tipCount))) {
      var lv = 0;
      for (final o in shown) {
        if (!identical(o.$1, b.$1) && level.containsKey(o.$1) && b.$1.nodes.contains(o.$1)) lv = math.max(lv, level[o.$1]! + 1);
      }
      level[b.$1] = lv;
    }
    final colW = <int, double>{};
    for (final (n, m) in shown) {
      final tp = TextPainter(
        text: TextSpan(text: m.text, style: TextStyle(fontFamily: style.family, fontSize: style.barText, fontWeight: FontWeight.w700, fontStyle: m.italic ? FontStyle.italic : null)),
        textDirection: TextDirection.ltr,
      )..layout();
      colW[level[n]!] = math.max(colW[level[n]!] ?? 0, tp.width + 16);
    }
    final barW = colW.values.fold(0.0, (a, b) => a + b);
    // The horizontal scale stretches (or shrinks) the branches only; the
    // drawing grows sideways past the set width when it is over 1.
    // Branches: the width set by the user ([treeW]) or what is left beside
    // the labels, heatmap, bars and legend, but never under 45 % of the
    // figure: with many extras the figure grows instead of squeezing the
    // tree.
    final room = width - left - labelW - heatW - barW - groupW - sideW - style.margin;
    final branchW = treeW > 0 ? treeW : math.max(room, width * 0.45);
    final scale = math.max(10.0, branchW) * xScale / maxX;
    final w = math.max(width, left + maxX * scale + labelW + heatW + barW + groupW + sideW + style.margin);
    final labelsRight = left + maxX * scale + labelW;
    double colX(int lv) => labelsRight + heatW + [for (var i = 0; i < lv; i++) colW[i] ?? 0].fold(0.0, (a, b) => a + b) + 4;
    final barList = [for (final (n, m) in shown) (node: n, mark: m, x: colX(level[n]!))];
    final barH = barLength < 0 || shape != TreeShape.phylogram ? 0.0 : 40.0;
    final rowsBottom = top + rows.length * rowH;
    _LegendBox? box;
    Offset origin;
    double h;
    final bottom = style.margin + 4;
    if (legend.isEmpty) {
      origin = Offset.zero;
      h = rowsBottom + barH + bottom;
    } else if (side) {
      box = sideBox;
      origin = Offset(w - sideW + 12, top);
      h = math.max(rowsBottom + barH, top + box!.size.height) + bottom;
    } else {
      box = _LegendBox.compute(legend, style, vertical: false, maxWidth: w - 2 * left);
      origin = Offset(left, rowsBottom + barH + 8);
      h = origin.dy + box.size.height + bottom;
    }
    return _Layout._(root, rows, x, {for (final e in y.entries) e.key: top + e.value * rowH + rowH / 2}, Size(w, h), left, top, scale, collapsed, maxX,
        shape, rowH, groupW, w - sideW, barLength, box, origin, barList, heat, labelsRight + 4);
  }

  bool hidden(TreeNode n) => !_x.containsKey(n);

  Offset point(TreeNode n) => Offset(left + _x[n]! * scale, _y[n]!);

  /// Middle of the branch above [n] (where a shortened branch shows `//`).
  Offset breakPoint(TreeNode n, [BranchStyle style = BranchStyle.square]) {
    final o = point(n);
    if (n.parent == null) return Offset(o.dx - 6, o.dy);
    final p = point(n.parent!);
    return style == BranchStyle.slanted ? Offset.lerp(p, o, 0.5)! : Offset((p.dx + o.dx) / 2, o.dy);
  }

  /// Centre of the first and of the last row of clade [n] (null when it is
  /// hidden inside a collapsed clade).
  (double, double)? span(TreeNode n) {
    if (hidden(n)) return null;
    if (n.isTip || collapsed.contains(n)) return (_y[n]!, _y[n]!);
    double? a, b;
    for (final k in n.nodes) {
      if ((k.isTip || collapsed.contains(k)) && _y.containsKey(k)) {
        final v = _y[k]!;
        a = a == null ? v : math.min(a, v);
        b = b == null ? v : math.max(b, v);
      }
    }
    return a == null ? null : (a, b!);
  }

  /// Branches as drawing commands (`M`, `L`, `Q` with their points) and the
  /// node whose colour each takes. Square: the vertical bar of a node takes
  /// the node's colour, the horizontal branch to a child the child's.
  /// Rounded and slanted: one line from the node to each child, in the
  /// child's colour.
  Iterable<(List<(String, List<Offset>)>, TreeNode)> branchPaths(BranchStyle style) sync* {
    for (final n in _x.keys) {
      if (n.isTip || collapsed.contains(n)) continue;
      final o = point(n);
      if (style == BranchStyle.square) {
        final ys = n.children.map((c) => _y[c]!).toList();
        yield ([('M', [Offset(o.dx, ys.reduce(math.min))]), ('L', [Offset(o.dx, ys.reduce(math.max))])], n);
      }
      for (final c in n.children) {
        final e = point(c);
        switch (style) {
          case BranchStyle.square:
            yield ([('M', [Offset(o.dx, e.dy)]), ('L', [e])], c);
          case BranchStyle.slanted:
            yield ([('M', [o]), ('L', [e])], c);
          case BranchStyle.rounded:
            final dy = e.dy - o.dy;
            final r = math.min(rowH * 0.5, math.min(dy.abs(), e.dx - o.dx));
            if (r <= 0.5) {
              yield ([('M', [o]), ('L', [Offset(o.dx, e.dy)]), ('L', [e])], c);
            } else {
              final sign = dy.sign;
              yield ([
                ('M', [o]),
                ('L', [Offset(o.dx, e.dy - sign * r)]),
                ('Q', [Offset(o.dx, e.dy), Offset(o.dx + r, e.dy)]),
                ('L', [e]),
              ], c);
            }
        }
      }
    }
  }

  /// Runs of consecutive rows of one group: (name, top, bottom).
  List<(String, double, double)> groups(Map<TreeNode, String> of) {
    final out = <(String, double, double)>[];
    String? current;
    double start = 0;
    for (final (i, n) in rows.indexed) {
      final v = of[n];
      final y0 = top + i * rowH;
      if (v != current) {
        if (current != null) out.add((current, start, y0));
        current = v;
        start = y0;
      }
    }
    if (current != null) out.add((current, start, top + rows.length * rowH));
    return out;
  }

  /// Node under [p]: its horizontal branch (from the parent's x to its own)
  /// within half a row, or the tip label.
  TreeNode? hit(Offset p) {
    TreeNode? best;
    var bestD = double.infinity;
    final reach = math.max(7.0, rowH / 2 - 1);
    for (final n in _x.keys) {
      final o = point(n);
      final x0 = n.parent == null ? o.dx - 12 : point(n.parent!).dx;
      final x1 = n.isTip || collapsed.contains(n) ? groupRight : o.dx + 4;
      final dy = (p.dy - o.dy).abs();
      if (p.dx < x0 - 4 || p.dx > x1 || dy > reach) continue;
      if (dy < bestD) {
        bestD = dy;
        best = n;
      }
    }
    return best;
  }

  /// Scale bar of a round length about a fifth of the tree (phylogram).
  (Offset, Offset, String)? scaleBar() {
    if (shape != TreeShape.phylogram || barLength < 0) return null;
    var nice = barLength;
    if (nice <= 0) {
      final target = maxX / 5;
      final mag = math.pow(10, (math.log(target) / math.ln10).floor()).toDouble();
      nice = [1, 2, 5, 10].map((k) => k * mag).lastWhere((v) => v <= target, orElse: () => mag);
    }
    final y = rowsBottom + 14;
    return (Offset(left, y), Offset(left + nice * scale, y), _fmt(nice));
  }

  static String _fmt(double v) =>
      v >= 1 ? (v == v.roundToDouble() ? v.toStringAsFixed(0) : v.toStringAsFixed(2)) : double.parse(v.toStringAsPrecision(2)).toString();
}

// ---------------------------------------------------------------------------
// Legend
// ---------------------------------------------------------------------------

/// One entry of the figure legend: a coloured dot or box, a symbol (`//`)
/// or a heading, with its text.
class _LegendItem {
  const _LegendItem(this.label, {this.color, this.box = false, this.heading = false, this.italic = false, this.symbol});
  final String label;
  final Color? color;
  final bool box, heading, italic;
  final String? symbol;
}

/// Legend items placed in a box: below the tree each group is a line with
/// its heading in a left column (the items wrap under each other); on the
/// right each item has its own line under its heading.
class _LegendBox {
  _LegendBox(this.items, this.size, this.title, this.pad);
  final List<(_LegendItem, Offset)> items;
  final Size size;
  final String title;
  final double pad;

  static TextStyle textStyle(_LegendItem it, _Style s, Color color) => TextStyle(
        fontFamily: s.family,
        fontSize: s.legendFont,
        fontWeight: it.heading ? FontWeight.w700 : FontWeight.w400,
        fontStyle: it.italic ? FontStyle.italic : null,
        color: color,
      );

  static double _textW(String text, TextStyle style) => (TextPainter(text: TextSpan(text: text, style: style), textDirection: TextDirection.ltr)..layout()).width;

  static double itemWidth(_LegendItem it, _Style s) =>
      _textW(it.label, textStyle(it, s, const Color(0xFF000000))) + (it.symbol != null ? s.legendFont * 1.5 : (it.color != null ? 17 : 0));

  static _LegendBox compute(List<List<_LegendItem>> lines, _Style s, {required bool vertical, required double maxWidth}) {
    final pad = s.legendBox ? 10.0 : 0.0;
    final lineH = s.legendFont * 1.7;
    final title = s.legendTitle.trim();
    var y = pad + (title.isEmpty ? 0 : lineH + 2);
    final out = <(_LegendItem, Offset)>[];
    var right = 0.0;
    if (vertical) {
      for (final (i, line) in lines.indexed) {
        if (i > 0) y += lineH * 0.4;
        for (final it in line) {
          out.add((it, Offset(pad + (it.heading ? 0 : 8), y)));
          right = math.max(right, pad + (it.heading ? 0 : 8) + itemWidth(it, s));
          y += lineH;
        }
      }
    } else {
      // A grid: headings in the first column, the entries of every line in
      // shared columns, so each entry sits exactly below the one above it.
      final headingW = lines.map((l) => l.where((it) => it.heading).map((it) => itemWidth(it, s)).fold(0.0, math.max)).fold(0.0, math.max);
      final x0 = pad + (headingW > 0 ? headingW + 14 : 0);
      const gap = 16.0;
      final entries = [for (final l in lines) l.where((it) => !it.heading).toList()];
      final most = entries.map((e) => e.length).fold(0, math.max);
      List<double> widths(int cols) {
        final w = List<double>.filled(math.max(1, cols), 0);
        for (final e in entries) {
          for (final (j, it) in e.indexed) {
            w[j % cols] = math.max(w[j % cols], itemWidth(it, s));
          }
        }
        return w;
      }

      var cols = math.max(1, most);
      var colW = widths(cols);
      while (cols > 1 && x0 + colW.fold(0.0, (a, b) => a + b) + gap * (cols - 1) > maxWidth - pad) {
        cols--;
        colW = widths(cols);
      }
      final colX = <double>[];
      var cx = x0;
      for (final w in colW) {
        colX.add(cx);
        cx += w + gap;
      }
      for (final (i, line) in lines.indexed) {
        for (final it in line.where((it) => it.heading)) {
          out.add((it, Offset(pad, y)));
        }
        final e = entries[i];
        for (final (j, it) in e.indexed) {
          if (j > 0 && j % cols == 0) y += lineH;
          out.add((it, Offset(colX[j % cols], y)));
          right = math.max(right, colX[j % cols] + colW[j % cols]);
        }
        y += lineH;
      }
    }
    if (title.isNotEmpty) right = math.max(right, pad + _textW(title, TextStyle(fontFamily: s.family, fontSize: s.legendFont + 1, fontWeight: FontWeight.w700)));
    return _LegendBox(out, Size(right + pad, y + pad - lineH * 0.3), title, pad);
  }
}

/// Branch length as a label: three significant digits.
String _lengthText(double v) => v == 0 ? '0' : (v >= 1 ? v.toStringAsFixed(2) : v.toStringAsPrecision(2));

// ---------------------------------------------------------------------------
// Painter
// ---------------------------------------------------------------------------

class _TreePainter extends CustomPainter {
  _TreePainter({
    required this.layout,
    required this.zoom,
    required this.highlights,
    required this.style,
    required this.accent,
    required this.query,
    required this.selected,
    required this.hover,
    required this.span,
    required this.text,
    required this.outgroup,
    required this.truncated,
    required this.groupOf,
    required this.groupColors,
    required this.branchColors,
    required this.transparent,
  });

  final _Layout layout;
  final _Style style;

  /// Scale of the drawing on screen (1 for the PNG).
  final double zoom;

  /// Coloured boxes behind clades.
  final Map<TreeNode, Color> highlights;

  /// Selection and hover marks (screen only).
  final Color accent;
  final String query;
  final TreeNode? selected, hover;
  final TextSpan Function(TreeNode, TextStyle) span;
  final String Function(TreeNode) text;

  /// Outgroup names (lower case, spaces).
  final Set<String> outgroup;
  final Set<TreeNode> truncated;
  final Map<TreeNode, String> groupOf;
  final Map<String, Color> groupColors;
  final Map<TreeNode, Color> branchColors;
  final bool transparent;

  @override
  void paint(Canvas canvas, Size size) {
    canvas.save();
    canvas.scale(zoom);
    _paint(canvas);
    canvas.restore();
  }

  void _paint(Canvas canvas) {
    final l = layout;
    final s = style;
    final size = l.size;
    final rowH = l.rowH;
    final inSel = <TreeNode>{if (selected != null) ...selected!.nodes};
    if (!transparent) canvas.drawRect(Offset.zero & size, Paint()..color = s.background);
    // Title.
    if (s.title.trim().isNotEmpty) {
      final tp = TextPainter(
        text: TextSpan(text: s.title, style: TextStyle(fontFamily: s.family, fontSize: s.titleFont, fontWeight: FontWeight.w700, color: s.label)),
        textDirection: TextDirection.ltr,
      )..layout(maxWidth: size.width - 2 * l.left);
      tp.paint(canvas, Offset(l.left, s.margin));
    }
    // Row backgrounds by taxon, with the taxon's name on the right.
    final groups = s.groupRows ? l.groups(groupOf) : const <(String, double, double)>[];
    Color background(double y) {
      for (final (name, y0, y1) in groups) {
        if (y >= y0 && y < y1) return Color.alphaBlend(groupColors[name]!.withValues(alpha: 0.2), s.background);
      }
      return s.background;
    }

    for (final (name, y0, y1) in groups) {
      final c = groupColors[name]!;
      canvas.drawRect(Rect.fromLTRB(0, y0, l.groupRight, y1), Paint()..color = c.withValues(alpha: 0.2));
      final tp = TextPainter(
        text: TextSpan(text: name, style: TextStyle(fontFamily: s.family, fontSize: s.font, fontWeight: FontWeight.w700, color: _darker(c))),
        textDirection: TextDirection.ltr,
        maxLines: 1,
        ellipsis: '…',
      )..layout(maxWidth: math.max(10, l.groupW - 12));
      tp.paint(canvas, Offset(l.groupRight - 8 - tp.width, (y0 + y1) / 2 - tp.height / 2));
    }
    // Coloured boxes behind clades (larger clades first, so nested boxes
    // stay visible).
    for (final e in highlights.entries.toList()..sort((a, b) => b.key.tipCount.compareTo(a.key.tipCount))) {
      final sp = l.span(e.key);
      if (sp == null) continue;
      final x0 = (e.key.parent == null ? l.point(e.key).dx : (l.point(e.key).dx + l.point(e.key.parent!).dx) / 2) - 2;
      canvas.drawRRect(
        RRect.fromRectAndRadius(Rect.fromLTRB(x0, sp.$1 - rowH / 2 + 1, l.groupRight - 4, sp.$2 + rowH / 2 - 1), const Radius.circular(5)),
        Paint()..color = e.value.withValues(alpha: s.highlightAlpha),
      );
    }
    // Band behind the selected clade's rows.
    if (selected != null) {
      final ys = [for (final n in l.rows) if (inSel.contains(n)) l.point(n).dy];
      if (ys.isNotEmpty) {
        canvas.drawRRect(
          RRect.fromRectAndRadius(Rect.fromLTRB(l.point(selected!).dx - 6, ys.first - rowH / 2, l.groupRight - 4, ys.last + rowH / 2), const Radius.circular(6)),
          Paint()..color = accent.withValues(alpha: 0.10),
        );
      }
    }
    // Row of the hovered tip, to follow it to its name.
    final hv = hover;
    if (hv != null && (hv.isTip || l.collapsed.contains(hv)) && !l.hidden(hv)) {
      final o = l.point(hv);
      canvas.drawRect(Rect.fromLTRB(0, o.dy - rowH / 2, l.groupRight, o.dy + rowH / 2), Paint()..color = accent.withValues(alpha: 0.10));
    }
    // Branches, each in its colour (the selected clade in the accent).
    final cap = s.branchStyle == BranchStyle.square ? StrokeCap.square : StrokeCap.round;
    final paints = <Color, Paint>{};
    Paint line(Color c) => paints.putIfAbsent(
        c,
        () => Paint()
          ..color = c
          ..strokeWidth = s.lineW
          ..strokeCap = cap
          ..strokeJoin = StrokeJoin.round
          ..style = PaintingStyle.stroke);
    final selLine = Paint()
      ..color = accent
      ..strokeWidth = s.lineW + 1.2
      ..strokeCap = cap
      ..strokeJoin = StrokeJoin.round
      ..style = PaintingStyle.stroke;
    for (final (cmds, n) in l.branchPaths(s.branchStyle)) {
      final path = Path();
      for (final (op, pts) in cmds) {
        switch (op) {
          case 'M':
            path.moveTo(pts[0].dx, pts[0].dy);
          case 'L':
            path.lineTo(pts[0].dx, pts[0].dy);
          case 'Q':
            path.quadraticBezierTo(pts[0].dx, pts[0].dy, pts[1].dx, pts[1].dy);
        }
      }
      canvas.drawPath(path, inSel.contains(n) ? selLine : line(branchColors[n] ?? s.branch));
    }
    // Collapsed clades: a triangle reaching the deepest tip.
    final end = l.left + l.maxX * l.scale;
    for (final n in l.collapsed) {
      if (l.hidden(n)) continue;
      final o = l.point(n);
      final c = inSel.contains(n) ? accent : (branchColors[n] ?? s.branch);
      final path = Path()
        ..moveTo(o.dx, o.dy)
        ..lineTo(end, o.dy - rowH * 0.42)
        ..lineTo(end, o.dy + rowH * 0.42)
        ..close();
      canvas.drawPath(path, Paint()..color = c.withValues(alpha: 0.25));
      canvas.drawPath(path, Paint()
        ..color = c
        ..style = PaintingStyle.stroke);
    }
    // Shortened branches: a gap with // in the middle of the branch.
    if (l.shape == TreeShape.phylogram) {
      for (final n in truncated) {
        if (l.hidden(n) || n.parent == null) continue;
        final m = l.breakPoint(n, s.branchStyle);
        canvas.drawRect(Rect.fromCenter(center: m, width: 8, height: s.lineW + 10), Paint()..color = background(m.dy));
        final mark = Paint()
          ..color = s.branch
          ..strokeWidth = math.max(1.2, s.lineW * 0.8);
        for (final dx in [-2.5, 2.5]) {
          canvas.drawLine(Offset(m.dx + dx - 2, m.dy + 5), Offset(m.dx + dx + 2, m.dy - 5), mark);
        }
      }
    }
    // Tip labels.
    for (final n in s.showLabels ? l.rows : const <TreeNode>[]) {
      final o = l.point(n);
      final collapsed = l.collapsed.contains(n);
      var x = (collapsed || s.align ? end : o.dx) + s.labelGap;
      final out = s.markOutgroup && n.isTip && outgroup.contains(n.name.toLowerCase().replaceAll('_', ' ').trim());
      // Dotted guide from the tip to its aligned name.
      if (s.align && !collapsed && end - o.dx > 6) {
        final guide = Paint()
          ..color = _Style.guide
          ..strokeWidth = 1;
        for (var gx = o.dx + 4; gx < end + s.labelGap - 3; gx += 4) {
          canvas.drawLine(Offset(gx, o.dy), Offset(math.min(gx + 1.5, end + s.labelGap - 3), o.dy), guide);
        }
      }
      if (s.tipDots && n.isTip && !s.align) x += s.dot;
      final color = inSel.contains(n) ? accent : (s.groupLabels && groupOf[n] != null ? _darker(groupColors[groupOf[n]]!) : s.label);
      final tp = TextPainter(
        text: span(
          n,
          TextStyle(
            fontFamily: s.family,
            fontSize: s.font,
            fontWeight: identical(n, hover) || inSel.contains(n) ? FontWeight.w700 : FontWeight.w400,
            color: color,
          ),
        ),
        textDirection: TextDirection.ltr,
      )..layout();
      if (out) {
        canvas.drawRRect(
          RRect.fromRectAndRadius(Rect.fromLTWH(x - 3, o.dy - tp.height / 2 - 1, tp.width + 6, tp.height + 2), const Radius.circular(4)),
          Paint()..color = _Style.outgroup,
        );
      }
      tp.paint(canvas, Offset(x, o.dy - tp.height / 2));
    }
    // Bars with names beside clades.
    for (final b in l.bars) {
      final sp = l.span(b.node);
      if (sp == null) continue;
      final y0 = sp.$1 - rowH * 0.38, y1 = sp.$2 + rowH * 0.38;
      canvas.drawLine(Offset(b.x, y0), Offset(b.x, y1), Paint()
        ..color = b.mark.color
        ..strokeWidth = s.barStroke
        ..strokeCap = StrokeCap.butt);
      final tp = TextPainter(
        text: TextSpan(
          text: b.mark.text,
          style: TextStyle(fontFamily: s.family, fontSize: s.barText, fontWeight: FontWeight.w700, fontStyle: b.mark.italic ? FontStyle.italic : null, color: b.mark.color),
        ),
        textDirection: TextDirection.ltr,
      )..layout();
      tp.paint(canvas, Offset(b.x + 6, (y0 + y1) / 2 - tp.height / 2));
    }
    _drawHeat(_CanvasPen(canvas), l, s);
    // Points on the tips.
    if (s.tipDots) {
      for (final n in l.rows) {
        if (n.isTip) canvas.drawCircle(l.point(n), s.dot, Paint()..color = inSel.contains(n) ? accent : (branchColors[n] ?? s.branch));
      }
    }
    // Support values on internal branches (from the chosen minimum).
    if (s.support != TreeSupport.none) {
      final items = <_SupportItem>[];
      for (final n in l._x.keys) {
        if (n.isTip || n.parent == null || n.label == null) continue;
        final v = n.support;
        if ((v ?? 0) < s.supportMin) continue;
        items.add((o: l.point(n), text: n.label!, r: s.radius(v), color: s.supportColor(v)));
      }
      // On a background so a value never mixes with the branches; placed
      // so that no two values overlap.
      if (s.support == TreeSupport.numbers) {
        _drawSupports(_CanvasPen(canvas), s, items, background: (o) => background(o.dy).withValues(alpha: 0.85));
      }
      for (final it in items) {
        _drawMark(_CanvasPen(canvas), it.o, it.r, it.color, s.markShape);
      }
    }
    // Branch lengths above the branches.
    if (s.lengths) {
      for (final n in l._x.keys) {
        final len = n.length;
        if (n.parent == null || len == null) continue;
        final m = l.breakPoint(n, s.branchStyle);
        final tp = TextPainter(
          text: TextSpan(text: _lengthText(len), style: TextStyle(fontFamily: s.family, fontSize: math.max(7, s.supportFont - 1), color: _Style.muted)),
          textDirection: TextDirection.ltr,
        )..layout();
        tp.paint(canvas, Offset(m.dx - tp.width / 2, m.dy - tp.height - s.lineW / 2));
      }
    }
    // Hovered node.
    if (hover != null && !l.hidden(hover!)) {
      canvas.drawCircle(l.point(hover!), s.dot + 2.5, Paint()
        ..color = accent
        ..style = PaintingStyle.stroke
        ..strokeWidth = 2);
    }
    // Scale bar.
    final bar = l.scaleBar();
    if (bar != null) {
      final barLine = Paint()
        ..color = s.branch
        ..strokeWidth = s.lineW;
      canvas.drawLine(bar.$1, bar.$2, barLine);
      for (final x in [bar.$1.dx, bar.$2.dx]) {
        canvas.drawLine(Offset(x, bar.$1.dy - 4), Offset(x, bar.$1.dy + 4), barLine);
      }
      final tp = TextPainter(
        text: TextSpan(text: bar.$3, style: TextStyle(fontFamily: s.family, fontSize: s.legendFont, color: _Style.muted)),
        textDirection: TextDirection.ltr,
      )..layout();
      tp.paint(canvas, Offset((bar.$1.dx + bar.$2.dx) / 2 - tp.width / 2, bar.$1.dy + 6));
    }
    // Legend.
    final lg = l.legendBox;
    if (lg != null) {
      final o = l.legendOrigin;
      if (s.legendBox) {
        final r = RRect.fromRectAndRadius(o & lg.size, const Radius.circular(6));
        canvas.drawRRect(r, Paint()..color = s.background);
        canvas.drawRRect(r, Paint()
          ..color = const Color(0xFFD1D5DB)
          ..style = PaintingStyle.stroke);
      }
      if (lg.title.isNotEmpty) {
        final tp = TextPainter(
          text: TextSpan(text: lg.title, style: TextStyle(fontFamily: s.family, fontSize: s.legendFont + 1, fontWeight: FontWeight.w700, color: s.label)),
          textDirection: TextDirection.ltr,
        )..layout();
        tp.paint(canvas, o + Offset(lg.pad, lg.pad));
      }
      for (final (it, at) in lg.items) {
        final p = o + at;
        var x = p.dx;
        final mid = p.dy + s.legendFont * 0.65;
        if (it.symbol != null) {
          final tp = TextPainter(
            text: TextSpan(text: it.symbol, style: TextStyle(fontFamily: s.family, fontSize: s.legendFont, fontWeight: FontWeight.w700, color: s.branch)),
            textDirection: TextDirection.ltr,
          )..layout();
          tp.paint(canvas, Offset(x, p.dy));
          x += s.legendFont * 1.5;
        } else if (it.color != null) {
          if (it.box) {
            canvas.drawRRect(RRect.fromRectAndRadius(Rect.fromLTWH(x, mid - 5, 13, 10), const Radius.circular(2)), Paint()..color = it.color!);
          } else {
            canvas.drawCircle(Offset(x + 5, mid), 4.5, Paint()..color = it.color!);
          }
          x += 17;
        }
        final tp = TextPainter(text: TextSpan(text: it.label, style: _LegendBox.textStyle(it, s, it.heading ? s.label : const Color(0xFF4B5563))), textDirection: TextDirection.ltr)
          ..layout();
        tp.paint(canvas, Offset(x, p.dy));
      }
    }
  }

  @override
  bool shouldRepaint(_TreePainter old) => true;
}

/// The SALVE heatmap: one column per field, a cell per tip row, the column
/// titles turned above the first row.
void _drawHeat(_Pen pen, _Layout l, _Style st) {
  final h = l.heat;
  if (h == null) return;
  final step = h.cell + h.gap;
  for (final (j, title) in h.titles.indexed) {
    final cx = l.heatX + j * step + h.cell / 2;
    pen.text([_Run(title)], Offset(cx, l.top - 4), family: st.family, size: st.legendFont, color: _Style.muted, angle: -math.pi / 2);
  }
  // Joined columns (species delimitation): one bar per run of rows of the
  // same species, a small gap between species.
  for (final j in h.joined) {
    final x = l.heatX + j * step + h.cell * 0.2;
    final w = h.cell * 0.6;
    String? group;
    Color? fill;
    double top = 0, bottom = 0;
    void flush() {
      final f = fill;
      if (group != null && f != null) pen.rect(Rect.fromLTRB(x, top + 1, x + w, bottom - 1), f, radius: math.min(2, w / 2));
    }

    for (final (i, n) in l.rows.indexed) {
      final cells = h.cells[n.name];
      final c = n.isTip && cells != null && j < cells.length ? cells[j] : null;
      final y0 = l.top + i * l.rowH, y1 = y0 + l.rowH;
      if (c?.group != null && c!.group == group) {
        bottom = y1;
        continue;
      }
      flush();
      group = c?.group;
      fill = c?.fill;
      top = y0;
      bottom = y1;
    }
    flush();
  }
  for (final n in l.rows) {
    final cells = h.cells[n.name];
    if (cells == null) continue;
    final y = l.point(n).dy - h.cell / 2;
    for (final (j, c) in cells.indexed) {
      if (h.joined.contains(j)) continue;
      final r = Rect.fromLTWH(l.heatX + j * step, y, h.cell, h.cell);
      if (c.fill == null) {
        pen.rect(r, _Style.outgroup.withValues(alpha: 0.6), radius: 2);
        continue;
      }
      pen.rect(r, c.fill!, radius: 2);
      if (c.text.isNotEmpty) {
        final size = math.min(st.legendFont, h.cell * 0.55);
        final tw = _Pen.width([_Run(c.text)], st.family, size, bold: true);
        pen.text([_Run(c.text)], Offset(r.center.dx - tw / 2, r.center.dy), family: st.family, size: size, color: c.ink, bold: true);
      }
    }
  }
}

/// A support mark of radius [r] at [o] in [shape].
void _drawMark(_Pen pen, Offset o, double r, Color color, MarkShape shape) {
  switch (shape) {
    case MarkShape.circle:
      pen.circle(o, r, color);
    case MarkShape.square:
      pen.rect(Rect.fromCircle(center: o, radius: r * 0.9), color);
    case MarkShape.diamond:
      final d = r * 1.25;
      pen.fill([('M', [o.translate(0, -d)]), ('L', [o.translate(d, 0)]), ('L', [o.translate(0, d)]), ('L', [o.translate(-d, 0)])], color);
    case MarkShape.triangle:
      final d = r * 1.3;
      pen.fill([('M', [o.translate(0, -d)]), ('L', [o.translate(d * 0.95, d * 0.65)]), ('L', [o.translate(-d * 0.95, d * 0.65)])], color);
  }
}
