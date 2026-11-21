part of 'tree_viewer.dart';

/// The figure of a tree as JSON, saved in the project file (`.splc`) and
/// kept by the project while the viewer is closed: every option of the
/// editor (style, labels, colours, legend, export, comparison) and the
/// edits of the tree (rooting, rotations, pruned tips, collapsed and
/// shortened clades, clade colours, boxes and bars, renamed tips).
///
/// The tree is stored as it is drawn (nested nodes, exactly), and the
/// edited clades by their place in it (pre-order index), so they come back
/// on the same nodes.
const treeFigureVersion = 1;

Map<String, dynamic> _nodeJson(TreeNode n) => {
      if (n.name.isNotEmpty) 'n': n.name,
      if (n.length != null) 'l': n.length,
      if (n.label != null) 's': n.label,
      if (n.children.isNotEmpty) 'c': [for (final k in n.children) _nodeJson(k)],
    };

TreeNode _nodeFrom(Map<String, dynamic> j) {
  final n = TreeNode(name: j['n'] as String? ?? '', length: (j['l'] as num?)?.toDouble(), label: j['s'] as String?);
  for (final k in (j['c'] as List? ?? const [])) {
    n.add(_nodeFrom((k as Map).cast<String, dynamic>()));
  }
  return n;
}

int _argb(Color c) => c.toARGB32();
Color _color(Object? v, Color fallback) => v is int ? Color(v) : fallback;

T _enum<T extends Enum>(List<T> values, Object? name, T fallback) => values.where((e) => e.name == name).firstOrNull ?? fallback;

extension _TreeFigureJson on _TreeViewerState {
  Map<String, dynamic> _figureJson() {
    final tree = _tree;
    final index = <TreeNode, int>{};
    if (tree != null) {
      for (final (i, n) in tree.nodes.indexed) {
        index[n] = i;
      }
    }
    List<int> nodes(Iterable<TreeNode> s) => [for (final n in s) ?index[n]]..sort();
    Map<String, int> colors(Map<TreeNode, Color> m) => {for (final e in m.entries) if (index[e.key] != null) '${index[e.key]}': _argb(e.value)};
    return {
      'version': treeFigureVersion,
      // Tree and its edits.
      if (tree != null) 'tree': _nodeJson(tree),
      'rooting': _rooting.name,
      'collapsed': nodes(_collapsed),
      'truncated': nodes(_truncated),
      'cladeColors': colors(_cladeColors),
      'highlights': colors(_highlights),
      'cladeLabels': {
        for (final e in _cladeLabels.entries)
          if (index[e.key] != null) '${index[e.key]}': {'text': e.value.text, 'color': _argb(e.value.color), 'italic': e.value.italic},
      },
      'renamed': {for (final e in _renamed.entries) e.key: {'text': e.value.text, 'italic': e.value.italic, 'bold': e.value.bold}},
      // Layout.
      'shape': _shape.name,
      'circular': _circular,
      'arcDeg': _arcDeg,
      'rotateDeg': _rotateDeg,
      'branchStyle': _branchStyle.name,
      'widthCm': _widthCm,
      'treeWidthCm': _treeWidthCm,
      'xScale': _xScale,
      'yScale': _yScale,
      'margin': _margin,
      'lineW': _lineW,
      'branchColor': _argb(_branchColor),
      'background': _argb(_background),
      'transparent': _transparent,
      'lengths': _lengths,
      'scaleBar': _scaleBar,
      'barLength': _barLength,
      'markOutgroup': _markOutgroup,
      // Support.
      'support': _support.name,
      'supportFont': _supportFont,
      'supportMin': _supportMin,
      'supportBelow': _supportBelow,
      'supportArrow': _supportArrow,
      'supportCurved': _supportCurved,
      'proportional': _proportional,
      'dot': _dot,
      'markShape': _markShape.name,
      'supportPreset': _supportPreset.name,
      if (_classesCustom != null) 'supportClasses': [for (final c in _classesCustom!) {'from': c.from, 'color': _argb(c.color)}],
      'supportInLegend': _supportInLegend,
      'supportLegendTitle': _supportLegendTitle,
      // Labels.
      'family': _family,
      'font': _font,
      'labelColor': _argb(_labelColor),
      'labelGap': _labelGap,
      'align': _align,
      'tipDots': _tipDots,
      'showLabels': _showLabels,
      'format': [for (final f in _format) {'field': f.field.name, 'italic': f.italic, 'bold': f.bold, 'parens': f.parens}],
      // Colours, bars and heatmap.
      if (_colorBy != null) 'colorBy': _colorBy!.name,
      'groupRows': _groupRows,
      'groupBranches': _groupBranches,
      'groupLabels': _groupLabels,
      'highlightOpacity': _highlightOpacity,
      if (_barRank != null) 'barRank': _barRank!.name,
      'barWidthPt': _barWidthPt,
      'barFontPt': _barFontPt,
      'heat': [for (final h in _heat) h.name],
      'heatCellPt': _heatCellPt,
      'heatGapPt': _heatGapPt,
      'delimitations': [..._delims],
      // Title and legend.
      'title': _title,
      'titleFont': _titleFont,
      'figLegend': _figLegend,
      'legendTitle': _legendTitle,
      'legendFont': _legendFont,
      'legendBox': _legendBox,
      'legendPos': _legendPos.name,
      // Export.
      'dpi': _dpi,
      // Comparison with another tree.
      'untangle': _untangle,
      'markConflicts': _markConflicts,
      'groupConnectors': _groupConnectors,
      'swapSides': _swapSides,
      'conflictColor': _argb(_conflictColor),
      'gapPct': _gapPct,
    };
  }

  /// Applies a saved figure. Its tree is used only when it has the same tips
  /// as the tree file (an analysis run again may have changed the tree);
  /// otherwise only the options come back.
  void _applyFigure(Map<String, dynamic> j) {
    double d(String k, double v) => (j[k] as num?)?.toDouble() ?? v;
    bool b(String k, bool v) => j[k] as bool? ?? v;
    String s(String k, String v) => j[k] as String? ?? v;

    final saved = j['tree'];
    if (saved is Map && _tree != null) {
      final t = _nodeFrom(saved.cast<String, dynamic>());
      List<String> tips(TreeNode n) => [for (final x in n.tips) x.name]..sort();
      if (tips(t).join('\u0000') == tips(_tree!).join('\u0000')) {
        _tree = t;
        _rooting = _enum(TreeRooting.values, j['rooting'], _rooting);
        final nodes = t.nodes.toList();
        TreeNode? at(Object? i) {
          final k = i is int ? i : int.tryParse('$i');
          return k != null && k >= 0 && k < nodes.length ? nodes[k] : null;
        }

        _collapsed
          ..clear()
          ..addAll([for (final i in j['collapsed'] as List? ?? const []) ?at(i)]);
        _truncated
          ..clear()
          ..addAll([for (final i in j['truncated'] as List? ?? const []) ?at(i)]);
        void colors(Map<TreeNode, Color> into, Object? m) {
          into.clear();
          for (final e in (m as Map? ?? const {}).entries) {
            final n = at(e.key);
            if (n != null && e.value is int) into[n] = Color(e.value as int);
          }
        }

        colors(_cladeColors, j['cladeColors']);
        colors(_highlights, j['highlights']);
        _cladeLabels.clear();
        for (final e in (j['cladeLabels'] as Map? ?? const {}).entries) {
          final n = at(e.key);
          final v = e.value as Map;
          if (n != null) _cladeLabels[n] = (text: v['text'] as String? ?? '', color: _color(v['color'], _Style.muted), italic: v['italic'] as bool? ?? false);
        }
      }
    }
    _renamed.clear();
    for (final e in (j['renamed'] as Map? ?? const {}).entries) {
      final v = e.value as Map;
      _renamed[e.key as String] = (text: v['text'] as String? ?? '', italic: v['italic'] as bool? ?? false, bold: v['bold'] as bool? ?? false);
    }

    _shape = _enum(TreeShape.values, j['shape'], _shape);
    _circular = b('circular', _circular);
    _arcDeg = d('arcDeg', _arcDeg);
    _rotateDeg = d('rotateDeg', _rotateDeg);
    _branchStyle = _enum(BranchStyle.values, j['branchStyle'], _branchStyle);
    _widthCm = d('widthCm', _widthCm);
    _treeWidthCm = d('treeWidthCm', _treeWidthCm);
    _xScale = d('xScale', _xScale);
    _yScale = d('yScale', _yScale);
    _margin = d('margin', _margin);
    _lineW = d('lineW', _lineW);
    _branchColor = _color(j['branchColor'], _branchColor);
    _background = _color(j['background'], _background);
    _transparent = b('transparent', _transparent);
    _lengths = b('lengths', _lengths);
    _scaleBar = b('scaleBar', _scaleBar);
    _barLength = d('barLength', _barLength);
    _markOutgroup = b('markOutgroup', _markOutgroup);

    _support = _enum(TreeSupport.values, j['support'], _support);
    _supportFont = d('supportFont', _supportFont);
    _supportMin = d('supportMin', _supportMin);
    _supportBelow = b('supportBelow', _supportBelow);
    _supportArrow = b('supportArrow', _supportArrow);
    _supportCurved = b('supportCurved', _supportCurved);
    _proportional = b('proportional', _proportional);
    _dot = d('dot', _dot);
    _markShape = _enum(MarkShape.values, j['markShape'], _markShape);
    _supportPreset = _enum(SupportPreset.values, j['supportPreset'], _supportPreset);
    final classes = j['supportClasses'];
    _classesCustom = classes is List
        ? [
            for (final c in classes.cast<Map>()) (from: (c['from'] as num?)?.toDouble() ?? 0, color: _color(c['color'], _Style.muted)),
          ]
        : null;
    _supportInLegend = b('supportInLegend', _supportInLegend);
    _supportLegendTitle = s('supportLegendTitle', _supportLegendTitle);

    _family = s('family', _family);
    _font = d('font', _font);
    _labelColor = _color(j['labelColor'], _labelColor);
    _labelGap = d('labelGap', _labelGap);
    _align = b('align', _align);
    _tipDots = b('tipDots', _tipDots);
    _showLabels = b('showLabels', _showLabels);
    final format = j['format'];
    if (format is List) {
      final parts = [
        for (final f in format.cast<Map>())
          if (TipField.values.where((x) => x.name == f['field']).firstOrNull case final field?)
            LabelPart(field, italic: f['italic'] as bool? ?? false, bold: f['bold'] as bool? ?? false, parens: f['parens'] as bool? ?? false),
      ];
      if (parts.isNotEmpty) _format = parts;
    }

    _colorBy = TaxRank.values.where((x) => x.name == j['colorBy']).firstOrNull;
    _groupRows = b('groupRows', _groupRows);
    _groupBranches = b('groupBranches', _groupBranches);
    _groupLabels = b('groupLabels', _groupLabels);
    _highlightOpacity = d('highlightOpacity', _highlightOpacity);
    _barRank = TaxRank.values.where((x) => x.name == j['barRank']).firstOrNull;
    _barWidthPt = d('barWidthPt', _barWidthPt);
    _barFontPt = d('barFontPt', _barFontPt);
    _heat
      ..clear()
      ..addAll([for (final h in j['heat'] as List? ?? const []) ?HeatField.values.where((x) => x.name == h).firstOrNull]);
    _heatCellPt = d('heatCellPt', _heatCellPt);
    _heatGapPt = d('heatGapPt', _heatGapPt);
    // Only the delimitations this project still has.
    _delims
      ..clear()
      ..addAll([for (final id in j['delimitations'] as List? ?? const []) if (widget.delimitations.any((x) => x.id == id)) id as String]);

    _title = s('title', _title);
    _titleFont = d('titleFont', _titleFont);
    _figLegend = b('figLegend', _figLegend);
    _legendTitle = s('legendTitle', _legendTitle);
    _legendFont = d('legendFont', _legendFont);
    _legendBox = b('legendBox', _legendBox);
    _legendPos = _enum(LegendPosition.values, j['legendPos'], _legendPos);

    _dpi = d('dpi', _dpi);

    _untangle = b('untangle', _untangle);
    _markConflicts = b('markConflicts', _markConflicts);
    _groupConnectors = b('groupConnectors', _groupConnectors);
    _swapSides = b('swapSides', _swapSides);
    _conflictColor = _color(j['conflictColor'], _conflictColor);
    _gapPct = d('gapPct', _gapPct);
  }
}
