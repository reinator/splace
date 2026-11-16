import 'dart:math' as math;

import 'package:flutter/gestures.dart';
import 'package:flutter/material.dart';
import 'package:flutter/services.dart';

import '../../core/fasta.dart';
import '../../core/genetic_code.dart';
import '../common.dart';
import '../fa.dart';
import '../theme.dart';
import 'kit.dart';
import 'splace_tooltip.dart';

/// Nucleotide colours, from the app palette: four that never mix up
/// (navy, light blue, orange, crimson).
const ntColors = {
  'A': SplacePalette.blue900, 'C': SplacePalette.sky500, 'G': SplacePalette.orange, 'T': SplacePalette.red900, 'U': SplacePalette.red900, //
};

/// Amino-acid colours by property, from the app palette: aliphatic and
/// aromatic (blues), polar (lilac), acidic (orange-red), basic (deep
/// blue), glycine, proline and cysteine (pale / yellow / orange), stop
/// (crimson).
const aaColors = {
  'A': SplacePalette.sky300, 'V': SplacePalette.sky300, 'L': SplacePalette.sky300, 'I': SplacePalette.sky300, 'M': SplacePalette.sky300, //
  'F': SplacePalette.sky500, 'W': SplacePalette.sky500, 'Y': SplacePalette.sky500,
  'S': SplaceBlue.c200, 'T': SplaceBlue.c200, 'N': SplaceBlue.c200, 'Q': SplaceBlue.c200,
  'D': SplacePalette.red500, 'E': SplacePalette.red500,
  'K': SplacePalette.blue600, 'R': SplacePalette.blue600, 'H': SplaceBlue.c300,
  'G': SplacePalette.sky100, 'P': SplacePalette.yellow, 'C': SplacePalette.orange,
  'X': Color(0xFFCBD5E1), '*': SplacePalette.red900, '?': Color(0xFFCBD5E1),
};

/// Readable letter on a residue colour.
Color residueInk(Color background) => background.computeLuminance() > 0.36 ? const Color(0xFF111827) : Colors.white;

/// Residues removed by trimming (grey, readable in both themes).
const _removedResidue = Color(0xFFB8BCC6);

/// Codon view: codons that vary without changing the amino acid, and
/// codons whose amino acid changes.
const _synonymous = SplaceBlue.c300;
const _nonSynonymous = SplacePalette.red700;

/// Column space of a view, so the raw / aligned / trimmed tabs of a marker
/// open at the same place.
enum ViewSpace {
  /// Unaligned input: positions are residues of each sequence.
  raw,

  /// Columns of the MAFFT alignment (aligned and trimmed views).
  mafft,

  /// Protein alignment of the codon-aware mode: one column is three MAFFT
  /// (DNA) columns.
  protein,
}

/// Where a view is scrolled to: the sequence at the top, how many of its
/// residues lie left of the view, and the MAFFT column at the left edge.
class AlignmentAnchor {
  const AlignmentAnchor({this.header, this.residue = 0, this.column, this.cell = 20, this.translate = false, this.codons = false, this.code});

  /// Sequence in the top row.
  final String? header;

  /// Nucleotides (amino acids × 3) of [header] before the left edge.
  final int residue;

  /// MAFFT column at the left edge (null in the raw view).
  final int? column;
  final double cell;
  final bool translate;

  /// The codon view (protein-coding DNA, one block per codon).
  final bool codons;

  /// Genetic code chosen in the viewer for the amino-acid view.
  final int? code;
}

/// Colour-coded, virtualized alignment viewer: only the visible window is
/// painted, so alignments with thousands of columns stay smooth. It has
/// its own scrollbars (outside the grid), scrolls with the wheel
/// (Shift + wheel: sideways), by dragging or with the scrollbars, and
/// reports where it is ([onAnchor]) so another view can open there.
class AlignmentViewer extends StatefulWidget {
  const AlignmentViewer({
    super.key,
    required this.records,
    this.protein = false,
    this.translatable = false,
    this.frame = 0,
    this.geneticCode = 1,
    this.height = 380,
    this.coding = false,
    this.columnMap,
    this.mapFrame = 0,
    this.removed,
    this.space = ViewSpace.mafft,
    this.anchor,
    this.onAnchor,
  });

  final List<SeqRecord> records;
  final bool protein;

  /// Protein-coding DNA alignment: offers the amino-acid view.
  final bool translatable;

  /// Protein-coding DNA: the hover info shows the codon position.
  final bool coding;

  /// Column where the first complete codon starts (translation and codon
  /// positions).
  final int frame;
  final int geneticCode;
  final double height;

  /// Compact trimmed alignments: the MAFFT column of each column, so the
  /// hover shows it and the codon position comes from the MAFFT alignment
  /// (in frame [mapFrame]).
  final List<int>? columnMap;
  final int mapFrame;

  /// Columns removed by trimming, drawn in grey (the trimmed view shows
  /// the whole MAFFT alignment).
  final List<bool>? removed;
  final ViewSpace space;

  /// Position to open at, and where the view moves to.
  final AlignmentAnchor? anchor;
  final ValueChanged<AlignmentAnchor>? onAnchor;

  @override
  State<AlignmentViewer> createState() => _AlignmentViewerState();
}

class _AlignmentViewerState extends State<AlignmentViewer> {
  /// Width of the names column, shared by every viewer (resizable).
  static double _nameW = 190;
  static const rowH = 20.0, consH = 24.0, rulerH = 20.0, track = 12.0;

  /// Scroll offsets (px) of the grid.
  final _x = ValueNotifier<double>(0), _y = ValueNotifier<double>(0);

  /// Hovered (row, column); the column is null over the names.
  final _hover = ValueNotifier<(int, int?)?>(null);
  /// Width of one column: fixed at the largest size (letters readable).
  static const _cell = 20.0;
  late bool _translate = widget.translatable && !widget.protein && (widget.anchor?.translate ?? false);

  /// Codon view: the nucleotides grouped by codon (one block per codon, in
  /// the DNA colours, not translated), with how much each codon varies
  /// above the alignment.
  late bool _codons = widget.translatable && !widget.protein && !_translate && (widget.anchor?.codons ?? false);

  /// Dots instead of the bases in columns where every sequence has the same
  /// base (the first sequence keeps its letters). Shared by every viewer.
  static bool _dots = false;
  bool _table = false;

  /// Genetic code of the amino-acid view, chosen in the toolbar (starts at
  /// the project code; never guessed from the data).
  late int _code = widget.anchor?.code ?? widget.geneticCode;
  int _sourceCode = 0;
  AlignmentAnchor? _pending;
  Size _viewport = Size.zero;
  int _cols = 0;

  /// The columns do not fit: the horizontal scrollbar is shown.
  bool _needH = true;

  /// Grid viewport the horizontal scrollbar was built with (rebuilt when
  /// the layout finds another one, so the thumb spans the real range).
  Size _trackViewport = Size.zero;

  List<SeqRecord>? _shown;
  List<double>? _conservation;
  List<bool>? _removed;
  Object? _source;
  bool _sourceTranslate = false;
  bool _sourceCodons = false;
  int _stops = 0;

  /// Display columns where every sequence has the same residue.
  List<bool> _identical = const [];

  /// Codon view: the amino acid of each codon of each row ('-' for a gap
  /// codon), and per codon: how much it varies
  /// (1 − share of the commonest codon) and how (0 none, 1 synonymous only,
  /// 2 the amino acid changes).
  List<String> _codonAa = const [];
  List<double> _codonVar = const [];
  List<int> _codonKind = const [];
  List<int> _codonCount = const [], _codonAaCount = const [];
  Map<String, int> _stopCodons = const {};
  Map<String, int>? _rowOf;

  bool get _translated => _translate && !widget.protein;
  bool get _codonView => _codons && !_translated && !widget.protein && widget.translatable;

  @override
  void initState() {
    super.initState();
    _pending = widget.anchor;
    _x.addListener(_emit);
    _y.addListener(_emit);
  }

  @override
  void dispose() {
    _x.dispose();
    _y.dispose();
    _hover.dispose();
    super.dispose();
  }

  // ---------------------------------------------------------------------
  // Data
  // ---------------------------------------------------------------------

  void _prepare() {
    if (identical(_source, widget.records) && _sourceTranslate == _translate && _sourceCodons == _codons && _sourceCode == _code && _shown != null) return;
    _source = widget.records;
    _sourceTranslate = _translate;
    _sourceCodons = _codons;
    _sourceCode = _code;
    _rowOf = null;
    final translated = _translated;
    final f = widget.frame;
    _shown = translated
        ? [for (final r in widget.records) SeqRecord(r.header, translate(r.sequence.substring(math.min(f, r.sequence.length)), _code))]
        : widget.records;
    final seqs = _shown!;
    _cols = seqs.isEmpty ? 0 : seqs.map((r) => r.length).reduce(math.max);

    // Removed columns in the shown space (a codon is grey when trimming
    // removed any of its bases).
    final removed = widget.removed;
    _removed = removed == null
        ? null
        : translated
            ? [for (var k = 0; k < _cols; k++) [for (var i = 0; i < 3; i++) f + 3 * k + i].any((c) => c < removed.length && removed[c])]
            : removed;

    // Internal stops (the last residue of each sequence may be a stop).
    _stops = 0;
    final codons = <String, int>{};
    if (translated) {
      for (var r = 0; r < seqs.length; r++) {
        final s = seqs[r].sequence;
        final last = s.lastIndexOf(RegExp('[^-]'));
        for (var k = 0; k < last; k++) {
          if (s.codeUnitAt(k) != 0x2a) continue; // *
          _stops++;
          final dna = widget.records[r].sequence;
          final at = f + 3 * k;
          if (at + 3 <= dna.length) codons.update(dna.substring(at, at + 3).toUpperCase().replaceAll('U', 'T'), (v) => v + 1, ifAbsent: () => 1);
        }
      }
    }
    _stopCodons = codons;

    final aa = widget.protein || translated;
    final cons = List<double>.filled(_cols, 0);
    final same = List<bool>.filled(_cols, false);
    for (var c = 0; c < _cols; c++) {
      final counts = <int, int>{};
      for (final r in seqs) {
        if (c >= r.length) continue;
        final ch = r.sequence.codeUnitAt(c) | 0x20; // lower case
        if (ch == 0x2d || ch == 0x2e || ch == 0x3f || ch == 0x6e && !aa || ch == 0x78 && aa) continue; // - . ? n x
        counts.update(ch, (v) => v + 1, ifAbsent: () => 1);
      }
      cons[c] = counts.isEmpty ? 0 : counts.values.reduce(math.max) / seqs.length;
      same[c] = counts.length == 1 && counts.values.first > 1;
    }
    _conservation = cons;
    _identical = same;
    _prepareCodons();
  }

  /// Amino acids and variation of every codon (codon view only).
  void _prepareCodons() {
    if (!_codonView) {
      _codonAa = const [];
      _codonKind = _codonCount = _codonAaCount = const [];
      _codonVar = const [];
      return;
    }
    final f = widget.frame;
    final seqs = widget.records;
    final n = math.max(0, (_cols - f) ~/ 3);
    final table = codonTable(_code);
    _codonAa = [
      for (final r in seqs)
        String.fromCharCodes([
          for (var k = 0; k < n; k++)
            () {
              final at = f + 3 * k;
              if (at + 3 > r.length) return 0x2d;
              final codon = r.sequence.substring(at, at + 3).toUpperCase().replaceAll('U', 'T');
              if (codon == '---') return 0x2d;
              return (table[codon] ?? 'X').codeUnitAt(0);
            }(),
        ]),
    ];
    final varL = List<double>.filled(n, 0), kindL = List<int>.filled(n, 0);
    final countL = List<int>.filled(n, 0), aaCountL = List<int>.filled(n, 0);
    for (var k = 0; k < n; k++) {
      final codons = <String, int>{};
      final aas = <int>{};
      var counted = 0;
      for (final (i, r) in seqs.indexed) {
        final at = f + 3 * k;
        if (at + 3 > r.length) continue;
        final codon = r.sequence.substring(at, at + 3).toUpperCase();
        // Gaps and ambiguous bases say nothing about variation.
        if (!RegExp(r'^[ACGTU]{3}$').hasMatch(codon)) continue;
        counted++;
        codons.update(codon, (v) => v + 1, ifAbsent: () => 1);
        aas.add(_codonAa[i].codeUnitAt(k));
      }
      countL[k] = codons.length;
      aaCountL[k] = aas.length;
      if (counted == 0) continue;
      varL[k] = 1 - codons.values.reduce(math.max) / counted;
      kindL[k] = codons.length <= 1 ? 0 : (aas.length <= 1 ? 1 : 2);
    }
    _codonVar = varL;
    _codonKind = kindL;
    _codonCount = countL;
    _codonAaCount = aaCountL;
  }

  // ---------------------------------------------------------------------
  // Scrolling and position
  // ---------------------------------------------------------------------

  /// Column positions: codons spaced apart in the codon view.
  _Geo get _geo => _codonView ? _Geo(_cell, widget.frame, (_cell * 0.4).clamp(3.0, 10.0)) : _Geo(_cell);

  double get _maxX => math.max(0, _geo.width(_cols) - _viewport.width);
  double get _maxY => math.max(0, (_shown?.length ?? 0) * rowH - _viewport.height);

  void _setX(double v) => _x.value = v.clamp(0, _maxX).toDouble();
  void _setY(double v) => _y.value = v.clamp(0, _maxY).toDouble();

  bool _isGap(int ch) => ch == 0x2d || ch == 0x2e; // - .

  /// Display column of DNA/protein column [i] of [widget.records].
  int _toDisplay(int i) => _translated ? math.max(0, (i - widget.frame) ~/ 3) : i;

  /// Column of [widget.records] at display column [col].
  int _fromDisplay(int col) => _translated ? widget.frame + 3 * col : col;

  /// MAFFT column (DNA, or amino acid in the amino-acid mode) of display
  /// column [col]; null in the raw view.
  int? _mafftColumn(int col) {
    switch (widget.space) {
      case ViewSpace.raw:
        return null;
      case ViewSpace.protein:
        return col * 3;
      case ViewSpace.mafft:
        final i = _fromDisplay(col);
        final map = widget.columnMap;
        if (map == null) return i;
        if (map.isEmpty) return 0;
        return i < map.length ? map[i] : map.last + 1;
    }
  }

  int _displayOfMafft(int column) {
    switch (widget.space) {
      case ViewSpace.raw:
        return 0;
      case ViewSpace.protein:
        return column ~/ 3;
      case ViewSpace.mafft:
        final map = widget.columnMap;
        if (map == null) return _toDisplay(column);
        // First kept column at or after [column] (the map is increasing).
        var lo = 0, hi = map.length;
        while (lo < hi) {
          final mid = (lo + hi) >> 1;
          if (map[mid] < column) {
            lo = mid + 1;
          } else {
            hi = mid;
          }
        }
        return _toDisplay(lo);
    }
  }

  /// Nucleotides of row [r] before display column [col].
  int _residuesBefore(int r, int col) {
    final s = widget.records[r].sequence;
    final end = math.min(s.length, _fromDisplay(col));
    var n = 0;
    for (var i = 0; i < end; i++) {
      if (!_isGap(s.codeUnitAt(i))) n++;
    }
    return widget.protein ? n * 3 : n;
  }

  /// Display column where row [r] reaches [residue] nucleotides.
  int _columnOfResidue(int r, int residue) {
    final s = widget.records[r].sequence;
    final target = widget.protein ? residue ~/ 3 : residue;
    var n = 0;
    for (var i = 0; i < s.length; i++) {
      if (_isGap(s.codeUnitAt(i))) continue;
      if (n == target) return _toDisplay(i);
      n++;
    }
    return _toDisplay(s.length);
  }

  int _rowOfHeader(String? header) {
    if (header == null) return 0;
    _rowOf ??= {for (final (i, r) in widget.records.indexed) r.header: i};
    return _rowOf![header] ?? 0;
  }

  AlignmentAnchor _anchorNow() {
    final seqs = _shown!;
    final col = _geo.col(_x.value);
    final row = seqs.isEmpty ? 0 : (_y.value / rowH).round().clamp(0, seqs.length - 1);
    return AlignmentAnchor(
      header: seqs.isEmpty ? null : seqs[row].header,
      residue: seqs.isEmpty ? 0 : _residuesBefore(row, col),
      column: _mafftColumn(col),
      cell: _cell,
      translate: _translate,
      codons: _codons,
      code: _code,
    );
  }

  void _emit() {
    if (_pending != null || _shown == null) return;
    widget.onAnchor?.call(_anchorNow());
  }

  void _restore(AlignmentAnchor a) {
    final row = _rowOfHeader(a.header);
    final col = a.column != null && widget.space != ViewSpace.raw ? _displayOfMafft(a.column!) : _columnOfResidue(row, a.residue);
    _setX(_geo.x(col));
    _setY(row * rowH);
  }

  /// 0 DNA, 1 amino acids, 2 codons.
  int get _view => _translate ? 1 : (_codons ? 2 : 0);

  void _setView(int v) {
    final a = _anchorNow();
    setState(() {
      _translate = v == 1;
      _codons = v == 2;
      _pending = AlignmentAnchor(header: a.header, residue: a.residue, column: a.column, cell: _cell, translate: _translate, codons: _codons, code: _code);
      _hover.value = null;
    });
  }

  /// Wheel: vertical, Shift + wheel or a sideways wheel: horizontal. The
  /// page scrolls instead when the grid is already at that edge.
  /// Dragging with the right (or middle) button scrolls the grid, also
  /// over the names; the left button keeps its own drag.
  bool _panning = false;

  static bool _panButtons(int buttons) => buttons & (kSecondaryMouseButton | kMiddleMouseButton) != 0;

  void _onPanDown(PointerDownEvent e) {
    if (e.kind == PointerDeviceKind.mouse && _panButtons(e.buttons)) setState(() => _panning = true);
  }

  void _onPanMove(PointerMoveEvent e, {bool vertical = false}) {
    if (!_panning || !_panButtons(e.buttons)) return;
    if (!vertical) _setX(_x.value - e.delta.dx);
    _setY(_y.value - e.delta.dy);
  }

  void _onPanEnd(PointerEvent e) {
    if (_panning) setState(() => _panning = false);
  }

  void _onSignal(PointerSignalEvent e) {
    if (e is! PointerScrollEvent) return;
    var dx = e.scrollDelta.dx, dy = e.scrollDelta.dy;
    if (dx == 0 && HardwareKeyboard.instance.isShiftPressed) {
      dx = dy;
      dy = 0;
    }
    final canH = dx != 0 && (dx < 0 ? _x.value > 0 : _x.value < _maxX);
    final canV = dy != 0 && (dy < 0 ? _y.value > 0 : _y.value < _maxY);
    if (!canH && !canV) return;
    GestureBinding.instance.pointerSignalResolver.register(e, (_) {
      if (canH) _setX(_x.value + dx);
      if (canV) _setY(_y.value + dy);
    });
  }

  /// Next internal stop codon after the hovered (or first visible) cell.
  void _nextStop() {
    final seqs = _shown!;
    final h = _hover.value;
    final (int, int) start = h?.$2 != null ? (h!.$2!, h.$1) : (_geo.col(_x.value), -1);
    (int, int)? found;
    bool isStop(int r, int c) {
      final s = seqs[r].sequence;
      if (c >= s.length || s.codeUnitAt(c) != 0x2a) return false;
      return s.substring(c + 1).contains(RegExp('[^-]'));
    }

    for (var step = 0; step <= _cols && found == null; step++) {
      final int c = (start.$1 + step) % math.max(1, _cols);
      for (var r = step == 0 ? start.$2 + 1 : 0; r < seqs.length; r++) {
        if (isStop(r, c)) {
          found = (r, c);
          break;
        }
      }
    }
    if (found == null) return;
    final (r, c) = found;
    _setX(_geo.x(c) - _viewport.width / 3);
    _setY(r * rowH - _viewport.height / 2);
    _hover.value = (r, c);
  }

  /// Fits the names column to the longest name.
  void _fitNames(List<SeqRecord> seqs) {
    var w = 0.0;
    for (final r in seqs.take(2000)) {
      final tp = TextPainter(
        text: TextSpan(text: r.header, style: const TextStyle(fontFamily: monoFont, fontSize: 12.5)),
        textDirection: TextDirection.ltr,
      )..layout();
      w = math.max(w, tp.width);
    }
    setState(() => _nameW = (w + 22).clamp(80.0, 900.0));
  }

  // ---------------------------------------------------------------------
  // Build
  // ---------------------------------------------------------------------

  @override
  Widget build(BuildContext context) {
    _prepare();
    final t = context.tokens;
    final l10n = context.l10n;
    final seqs = _shown!;
    final aa = widget.protein || _translated;
    final cols = _cols;
    final bodyH = math.min(widget.height, consH + rulerH + seqs.length * rowH + 2);
    final gridH = bodyH - consH - rulerH;
    final dark = Theme.of(context).brightness == Brightness.dark;
    final removedCount = _removed?.where((r) => r).length ?? 0;

    return Container(
      decoration: BoxDecoration(borderRadius: BorderRadius.circular(10), border: Border.all(color: t.border), color: t.surface),
      clipBehavior: Clip.antiAlias,
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        // Toolbar.
        Container(
          color: t.surface2,
          padding: const EdgeInsets.symmetric(horizontal: 10, vertical: 6),
          child: Wrap(spacing: 10, runSpacing: 6, crossAxisAlignment: WrapCrossAlignment.center, children: [
            if (widget.translatable && !widget.protein)
              Seg<int>(value: _view, options: [(0, l10n.viewDna), (1, l10n.viewTranslate), (2, l10n.viewCodons)], onChanged: _setView),
            SplaceTooltip(
              message: l10n.viewDotsInfo,
              child: ChipButton(
                label: l10n.viewDots,
                icon: FA.ellipsis,
                keepIcon: true,
                style: _dots ? ChipButtonStyle.accent : ChipButtonStyle.neutral,
                onPressed: () => setState(() => _dots = !_dots),
              ),
            ),
            Text('${seqs.length} × $cols ${aa ? 'aa' : 'bp'}', style: TextStyle(fontSize: 12, color: t.textMuted)),
            if (_translated || _codonView) ...[
              SizedBox(
                width: 250,
                child: DropdownButtonFormField<int>(
                  icon: dropdownChevron,
                  initialValue: _code,
                  isExpanded: true,
                  isDense: true,
                  style: Theme.of(context).textTheme.bodyMedium?.copyWith(fontSize: 12.5, color: t.text),
                  decoration: InputDecoration(isDense: true, labelText: l10n.geneticCode, contentPadding: const EdgeInsets.symmetric(horizontal: 10, vertical: 8)),
                  items: [for (final e in geneticCodeNames.entries) DropdownMenuItem(value: e.key, child: Text(e.value, overflow: TextOverflow.ellipsis))],
                  onChanged: (v) => setState(() {
                    _code = v!;
                    _hover.value = null;
                    widget.onAnchor?.call(_anchorNow());
                  }),
                ),
              ),
              if (_translated) TagChip(l10n.possibleStops('$_stops'), color: _stops > 0 ? SplaceTokens.danger : SplaceTokens.success),
              if (_translated && _stops > 0) ChipButton(label: l10n.nextStop, onPressed: _nextStop),
              ChipButton(
                label: l10n.translationTable,
                style: _table ? ChipButtonStyle.accent : ChipButtonStyle.neutral,
                onPressed: () => setState(() => _table = !_table),
              ),
            ],
            if (_codonView)
              TagChip(
                l10n.codonsVarying('${_codonKind.where((k) => k > 0).length}', '${_codonKind.length}', '${_codonKind.where((k) => k == 2).length}'),
                color: t.textMuted,
              ),
            if (removedCount > 0)
              SplaceTooltip(
                message: l10n.removedColumnsInfo,
                child: TagChip(l10n.removedColumns('$removedCount', '$cols'), color: t.textMuted),
              ),
            Row(mainAxisSize: MainAxisSize.min, children: _legend(aa, t, removedCount > 0)),
          ]),
        ),
        if ((_translated || _codonView) && _table) _CodeTable(code: _code, stops: _stopCodons),
        SizedBox(
          height: bodyH,
          child: Row(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
            // Names (follow the vertical offset).
            SizedBox(
              width: _nameW,
              child: Column(children: [
                Container(
                  height: consH + rulerH,
                  alignment: Alignment.bottomLeft,
                  padding: const EdgeInsets.only(left: 8, bottom: 4),
                  child: Column(
                    mainAxisAlignment: MainAxisAlignment.end,
                    crossAxisAlignment: CrossAxisAlignment.start,
                    children: [Kicker(_codonView ? l10n.codonVariation : l10n.conservation, size: 9.5), const SizedBox(height: 6), Kicker(l10n.viewerNames, size: 9.5)],
                  ),
                ),
                Expanded(
                  child: Listener(
                    onPointerSignal: _onSignal,
                    onPointerDown: _onPanDown,
                    onPointerMove: (e) => _onPanMove(e, vertical: true),
                    onPointerUp: _onPanEnd,
                    onPointerCancel: _onPanEnd,
                    child: MouseRegion(
                      cursor: _panning ? SystemMouseCursors.grabbing : MouseCursor.defer,
                      onHover: (e) {
                        final r = ((e.localPosition.dy + _y.value) / rowH).floor();
                        _hover.value = r >= 0 && r < seqs.length ? (r, null) : null;
                      },
                      onExit: (_) => _hover.value = null,
                      child: ClipRect(
                        child: CustomPaint(
                          size: Size.infinite,
                          painter: _NamesPainter(seqs: seqs, rowH: rowH, y: _y, hover: _hover, text: t.text, band: t.accent, bandText: t.accentText),
                        ),
                      ),
                    ),
                  ),
                ),
              ]),
            ),
            // Drag to resize the names column; double-click to fit.
            MouseRegion(
              cursor: SystemMouseCursors.resizeColumn,
              child: GestureDetector(
                behavior: HitTestBehavior.opaque,
                onHorizontalDragUpdate: (d) => setState(() => _nameW = (_nameW + d.delta.dx).clamp(80.0, 900.0)),
                onDoubleTap: () => _fitNames(seqs),
                child: SizedBox(
                  width: 8,
                  child: Center(child: Container(width: 2, color: t.border)),
                ),
              ),
            ),
            Expanded(
              child: LayoutBuilder(builder: (context, box) {
                final needV = seqs.length * rowH > gridH + 0.5;
                final width = box.maxWidth - (needV ? track : 0);
                _viewport = Size(width, gridH);
                // Clamp after a resize or zoom, then open at the anchor.
                WidgetsBinding.instance.addPostFrameCallback((_) {
                  if (!mounted) return;
                  final a = _pending;
                  if (a != null) {
                    _pending = null;
                    _restore(a);
                    _emit();
                  } else {
                    _setX(_x.value);
                    _setY(_y.value);
                  }
                  final need = _geo.width(_cols) > _viewport.width + 0.5;
                  if (need != _needH || _viewport != _trackViewport) {
                    setState(() {
                      _needH = need;
                      _trackViewport = _viewport;
                    });
                  }
                });
                return Column(children: [
                  SizedBox(
                    height: consH,
                    child: ClipRect(
                      child: CustomPaint(
                        size: Size(box.maxWidth, consH),
                        painter: _codonView
                            ? _CodonVariationPainter(_codonVar, _codonKind, widget.frame, _geo, _x, width, t.surface3)
                            : _ConservationPainter(_conservation!, _removed, _cell, _x, width, t.accent, t.surface3, t.textFaint),
                      ),
                    ),
                  ),
                  SizedBox(
                    height: rulerH,
                    child: ClipRect(
                      child: CustomPaint(
                        size: Size(box.maxWidth, rulerH),
                        painter: _RulerPainter(cols, _geo, _x, width, _hover, t.textMuted, t.border, t.accent),
                      ),
                    ),
                  ),
                  Expanded(
                    child: Row(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
                      Expanded(
                        child: Listener(
                          onPointerSignal: _onSignal,
                          onPointerDown: _onPanDown,
                          onPointerMove: _onPanMove,
                          onPointerUp: _onPanEnd,
                          onPointerCancel: _onPanEnd,
                          child: MouseRegion(
                            cursor: _panning ? SystemMouseCursors.grabbing : MouseCursor.defer,
                            onHover: (e) {
                              final r = ((e.localPosition.dy + _y.value) / rowH).floor();
                              final c = _geo.col(e.localPosition.dx + _x.value);
                              _hover.value = r >= 0 && r < seqs.length && c >= 0 && c < cols ? (r, c) : null;
                            },
                            onExit: (_) => _hover.value = null,
                            child: GestureDetector(
                              behavior: HitTestBehavior.opaque,
                              onPanUpdate: (d) {
                                _setX(_x.value - d.delta.dx);
                                _setY(_y.value - d.delta.dy);
                              },
                              child: ClipRect(
                                child: CustomPaint(
                                  size: Size.infinite,
                                  painter: _GridPainter(
                                    seqs: seqs,
                                    removed: _removed,
                                    geo: _geo,
                                    rowH: rowH,
                                    x: _x,
                                    y: _y,
                                    hover: _hover,
                                    aa: aa,
                                    gap: t.surface3,
                                    mark: t.text,
                                    accent: t.accent,
                                    text: dark ? Colors.black : const Color(0xFF111827),
                                    removedGap: t.border,
                                    // Same rule in every view: a dot where the column is the
                                    // same in all sequences (codons included).
                                    dots: _dots ? _identical : null,
                                    dot: t.textFaint,
                                    codonFrame: _codonView ? widget.frame : null,
                                    codonAa: _codonAa,
                                  ),
                                ),
                              ),
                            ),
                          ),
                        ),
                      ),
                      if (needV)
                        ScrollTrack(
                          axis: Axis.vertical,
                          offset: _y,
                          viewport: gridH,
                          content: seqs.length * rowH,
                          onChanged: _setY,
                        ),
                    ]),
                  ),
                ]);
              }),
            ),
          ]),
        ),
        // Horizontal scrollbar below the grid, with the removed regions.
        if (_needH)
          Container(
            color: t.surface,
            padding: const EdgeInsets.only(top: 2, bottom: 2),
            child: Row(children: [
              SizedBox(width: _nameW + 8),
              Expanded(
                child: ScrollTrack(
                  axis: Axis.horizontal,
                  offset: _x,
                  viewport: _trackViewport.width,
                  content: _geo.width(cols),
                  onChanged: _setX,
                  marks: _removed,
                ),
              ),
              if (seqs.length * rowH > gridH + 0.5) const SizedBox(width: track),
            ]),
          ),
        // Position of the hovered base.
        Container(
          padding: const EdgeInsets.symmetric(horizontal: 10, vertical: 6),
          decoration: BoxDecoration(color: t.surface2, border: Border(top: BorderSide(color: t.border))),
          child: ValueListenableBuilder<(int, int?)?>(valueListenable: _hover, builder: (context, hover, _) => _hoverInfo(context, seqs, hover, aa)),
        ),
      ]),
    );
  }

  Widget _hoverInfo(BuildContext context, List<SeqRecord> seqs, (int, int?)? hover, bool aa) {
    final t = context.tokens;
    final l10n = context.l10n;
    final style = TextStyle(fontSize: 12, color: t.textMuted);
    if (hover == null || hover.$1 >= seqs.length) return Text(l10n.viewerHint, maxLines: 1, overflow: TextOverflow.ellipsis, style: style);
    final (row, col0) = hover;
    final name = TextSpan(text: seqs[row].header, style: TextStyle(fontWeight: FontWeight.w700, color: t.text));
    if (col0 == null) {
      final residues = seqs[row].sequence.replaceAll(RegExp('[-.]'), '').length;
      return Text.rich(TextSpan(style: style, children: [name, TextSpan(text: '  ·  $residues ${aa ? 'aa' : 'bp'}')]),
          maxLines: 1, overflow: TextOverflow.ellipsis);
    }
    final col = col0;
    final s = seqs[row].sequence;
    final ch = col < s.length ? s[col].toUpperCase() : '-';
    final gap = ch == '-' || ch == '.';
    var pos = 0;
    for (var i = 0; i <= col && i < s.length; i++) {
      final c = s[i];
      if (c != '-' && c != '.') pos++;
    }
    final map = widget.columnMap;
    final source = map != null && col < map.length ? map[col] : null;
    final codonPos = widget.coding && !aa
        ? source != null
            ? (source - widget.mapFrame) % 3 + 1
            : (col - widget.frame) % 3 + 1
        : null;
    String? codon;
    if (_translated && !gap) {
      final dna = widget.records[row].sequence;
      final at = widget.frame + 3 * col;
      if (at + 3 <= dna.length) codon = dna.substring(at, at + 3).toUpperCase();
    }
    final removed = _removed != null && col < _removed!.length && _removed![col];
    final same = gap ? 0 : seqs.where((r) => col < r.length && r.sequence[col].toUpperCase() == ch).length;
    return Text.rich(
      TextSpan(style: style, children: [
        name,
        TextSpan(text: '  ·  ${l10n.viewerColumn('${col + 1}')}${source == null ? '' : ' (${l10n.viewerMafftColumn('${source + 1}')})'}  ·  '),
        TextSpan(
          text: gap ? l10n.viewerGap : ch,
          style: TextStyle(fontFamily: monoFont, fontWeight: FontWeight.w700, color: gap ? t.textFaint : ((aa ? aaColors[ch] : ntColors[ch]) ?? t.text)),
        ),
        if (codon != null) TextSpan(text: '  ·  ${l10n.viewerCodon(codon)}'),
        ..._codonHover(context, row, col, gap),
        if (!gap) TextSpan(text: '  ·  ${l10n.viewerSeqPos('$pos')}'),
        if (codonPos != null && !gap) TextSpan(text: '  ·  ${l10n.viewerCodonPos('$codonPos')}'),
        if (!gap) TextSpan(text: '  ·  ${l10n.viewerIdentical('${(100 * same / seqs.length).round()}')}'),
        if (removed) TextSpan(text: '  ·  ${l10n.viewerRemoved}', style: TextStyle(fontWeight: FontWeight.w700, color: t.text)),
      ]),
      maxLines: 1,
      overflow: TextOverflow.ellipsis,
    );
  }

  /// Codon view: the codon under the pointer, its amino acid and how much
  /// that codon varies in the alignment.
  List<InlineSpan> _codonHover(BuildContext context, int row, int col, bool gap) {
    if (!_codonView || col < widget.frame) return const [];
    final k = (col - widget.frame) ~/ 3;
    if (k >= _codonKind.length) return const [];
    final l10n = context.l10n;
    final at = widget.frame + 3 * k;
    final s = widget.records[row].sequence;
    final codon = at + 3 <= s.length ? s.substring(at, at + 3).toUpperCase() : '';
    final aa = _codonAa[row][k];
    return [
      TextSpan(text: '  ·  ${l10n.viewerCodonN('${k + 1}')}: '),
      TextSpan(text: codon, style: const TextStyle(fontFamily: monoFont, fontWeight: FontWeight.w700)),
      if (aa != '-') TextSpan(text: ' → ${_aaNames[aa] ?? aa}'),
      TextSpan(text: '  ·  ${l10n.codonColumnInfo('${_codonCount[k]}', '${_codonAaCount[k]}')}'),
    ];
  }

  List<Widget> _legend(bool aa, SplaceTokens t, bool removed) {
    Widget dot(Color c, String s) => Padding(
          padding: const EdgeInsets.only(left: 8),
          child: Row(mainAxisSize: MainAxisSize.min, children: [
            Container(width: 10, height: 10, decoration: BoxDecoration(color: c, borderRadius: BorderRadius.circular(2))),
            const SizedBox(width: 3),
            Text(s, style: TextStyle(fontSize: 11, color: t.textMuted)),
          ]),
        );
    final l10n = context.l10n;
    if (_codonView) {
      return [
        for (final b in ['A', 'C', 'G', 'T']) dot(ntColors[b]!, b),
        dot(_synonymous, l10n.codonLegendSyn),
        dot(_nonSynonymous, l10n.codonLegendNonSyn),
        dot(t.surface3, l10n.gapLegend),
        if (removed) dot(_removedResidue, l10n.removedLegend),
      ];
    }
    return [
      if (aa) ...[dot(aaColors['L']!, 'AA'), dot(aaColors['*']!, l10n.stopLegend)] else ...[
        dot(ntColors['A']!, 'A'),
        dot(ntColors['T']!, 'T/U'),
        dot(ntColors['G']!, 'G'),
        dot(ntColors['C']!, 'C'),
      ],
      dot(t.surface3, l10n.gapLegend),
      if (removed) dot(_removedResidue, l10n.removedLegend),
    ];
  }
}

// ---------------------------------------------------------------------------
// Translation table
// ---------------------------------------------------------------------------

const _aaNames = {
  'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys', 'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile', //
  'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro', 'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val', '*': 'Stop',
};

/// The genetic code used to translate (64 codons, the classic 4 × 16
/// layout): stop codons in red with how many internal stops of the
/// alignment use them, codons that differ from the standard code outlined.
class _CodeTable extends StatelessWidget {
  const _CodeTable({required this.code, required this.stops});
  final int code;
  final Map<String, int> stops;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final l10n = context.l10n;
    final table = codonTable(code), standard = codonTable(1);
    const bases = ['T', 'C', 'A', 'G'];
    Widget cell(String codon) {
      final aa = table[codon]!;
      final stop = aa == '*';
      final differs = standard[codon] != aa;
      final n = stops[codon] ?? 0;
      final fg = residueInk(aaColors[aa] ?? t.surface);
      return SplaceTooltip(
        message: differs ? l10n.codonDiffers(codon, _aaNames[aa]!, _aaNames[standard[codon]]!) : '$codon → ${_aaNames[aa]}',
        child: Container(
          margin: const EdgeInsets.all(1.5),
          padding: const EdgeInsets.symmetric(horizontal: 6, vertical: 3),
          decoration: BoxDecoration(
            color: stop ? SplaceTokens.danger : (aaColors[aa] ?? t.surface3).withValues(alpha: 0.35),
            borderRadius: BorderRadius.circular(5),
            border: Border.all(color: differs ? t.accent : Colors.transparent, width: 1.5),
          ),
          child: Row(children: [
            Text(codon, style: TextStyle(fontFamily: monoFont, fontSize: 11.5, fontWeight: FontWeight.w600, color: fg)),
            const SizedBox(width: 6),
            Expanded(child: Text(_aaNames[aa]!, style: TextStyle(fontSize: 11.5, color: fg))),
            if (n > 0)
              Container(
                padding: const EdgeInsets.symmetric(horizontal: 5),
                decoration: BoxDecoration(color: Colors.white, borderRadius: BorderRadius.circular(99)),
                child: Text('×$n', style: const TextStyle(fontSize: 10.5, fontWeight: FontWeight.w700, color: SplaceTokens.danger)),
              ),
          ]),
        ),
      );
    }

    return Container(
      padding: const EdgeInsets.fromLTRB(12, 10, 12, 12),
      decoration: BoxDecoration(color: t.surface, border: Border(bottom: BorderSide(color: t.border))),
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        Text(l10n.translationTableTitle('$code', geneticCodeName(code)), style: TextStyle(fontSize: 13, fontWeight: FontWeight.w700, color: t.text)),
        const SizedBox(height: 2),
        Text(l10n.translationTableBody, style: TextStyle(fontSize: 12, color: t.textMuted)),
        const SizedBox(height: 8),
        // Rows: first base × third base; columns: second base.
        Row(children: [
          const SizedBox(width: 22),
          for (final b in bases) Expanded(child: Center(child: Kicker(l10n.codonSecondBase(b), size: 10))),
        ]),
        for (final first in bases)
          Padding(
            padding: const EdgeInsets.only(top: 4),
            child: Row(crossAxisAlignment: CrossAxisAlignment.center, children: [
              SizedBox(width: 22, child: Text(first, style: TextStyle(fontFamily: monoFont, fontWeight: FontWeight.w700, color: t.textMuted))),
              for (final second in bases)
                Expanded(child: Column(children: [for (final third in bases) cell('$first$second$third')])),
            ]),
          ),
        const SizedBox(height: 8),
        Wrap(spacing: 14, runSpacing: 4, children: [
          _key(context, SplaceTokens.danger, l10n.stopLegend, filled: true),
          _key(context, t.accent, l10n.codonDiffersLegend),
          if (stops.isNotEmpty) Text(l10n.stopCodonsFound(stops.entries.map((e) => '${e.key} ×${e.value}').join(', ')), style: TextStyle(fontSize: 11.5, color: t.textMuted)),
        ]),
      ]),
    );
  }

  Widget _key(BuildContext context, Color c, String label, {bool filled = false}) => Row(mainAxisSize: MainAxisSize.min, children: [
        Container(
          width: 12,
          height: 12,
          decoration: BoxDecoration(color: filled ? c : null, borderRadius: BorderRadius.circular(3), border: Border.all(color: c, width: 1.5)),
        ),
        const SizedBox(width: 5),
        Text(label, style: TextStyle(fontSize: 11.5, color: context.tokens.textMuted)),
      ]);
}

// ---------------------------------------------------------------------------
// Scrollbar
// ---------------------------------------------------------------------------

/// Scrollbar outside a painted grid (the alignment viewer, the matrix
/// table): a rounded track and a thumb that can be
/// dragged; clicking the track jumps there. The horizontal one shows the
/// columns removed by trimming as grey marks.
class ScrollTrack extends StatefulWidget {
  const ScrollTrack({super.key, required this.axis, required this.offset, required this.viewport, required this.content, required this.onChanged, this.marks});
  final Axis axis;
  final ValueNotifier<double> offset;
  final double viewport, content;
  final ValueChanged<double> onChanged;
  final List<bool>? marks;

  static const size = 12.0;

  @override
  State<ScrollTrack> createState() => _ScrollTrackState();
}

class _ScrollTrackState extends State<ScrollTrack> {
  bool _hot = false;
  double _grab = 0;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final horizontal = widget.axis == Axis.horizontal;
    return LayoutBuilder(builder: (context, box) {
      final len = horizontal ? box.maxWidth : box.maxHeight;
      final content = math.max(widget.content, 1.0);
      final thumb = math.min(len, math.max(32.0, len * widget.viewport / content));
      final travel = math.max(1.0, len - thumb);
      final maxOffset = math.max(0.0, content - widget.viewport);
      double pos(Offset p) => horizontal ? p.dx : p.dy;
      double thumbStart() => maxOffset == 0 ? 0 : widget.offset.value / maxOffset * travel;
      void moveTo(double p) => widget.onChanged(((p - _grab) / travel).clamp(0.0, 1.0) * maxOffset);
      void grab(double p) {
        final s = thumbStart();
        _grab = p >= s && p <= s + thumb ? p - s : thumb / 2;
      }

      return MouseRegion(
        cursor: SystemMouseCursors.click,
        onEnter: (_) => setState(() => _hot = true),
        onExit: (_) => setState(() => _hot = false),
        child: GestureDetector(
          behavior: HitTestBehavior.opaque,
          onTapDown: (d) {
            grab(pos(d.localPosition));
            moveTo(pos(d.localPosition));
          },
          onPanStart: (d) => grab(pos(d.localPosition)),
          onPanUpdate: (d) => moveTo(pos(d.localPosition)),
          child: SizedBox(
            width: horizontal ? null : ScrollTrack.size,
            height: horizontal ? ScrollTrack.size : null,
            child: CustomPaint(
              size: Size.infinite,
              painter: _TrackPainter(
                horizontal: horizontal,
                offset: widget.offset,
                maxOffset: maxOffset,
                thumb: thumb,
                marks: widget.marks,
                track: t.surface2,
                border: t.border,
                thumbColor: t.accent.withValues(alpha: _hot ? 0.9 : 0.6),
                mark: _removedResidue,
              ),
            ),
          ),
        ),
      );
    });
  }
}

class _TrackPainter extends CustomPainter {
  _TrackPainter({
    required this.horizontal,
    required this.offset,
    required this.maxOffset,
    required this.thumb,
    required this.marks,
    required this.track,
    required this.border,
    required this.thumbColor,
    required this.mark,
  }) : super(repaint: offset);
  final bool horizontal;
  final ValueNotifier<double> offset;
  final double maxOffset, thumb;
  final List<bool>? marks;
  final Color track, border, thumbColor, mark;

  @override
  void paint(Canvas canvas, Size size) {
    final len = horizontal ? size.width : size.height;
    final cross = horizontal ? size.height : size.width;
    Rect along(double a, double b, double inset) =>
        horizontal ? Rect.fromLTRB(a, inset, b, cross - inset) : Rect.fromLTRB(inset, a, cross - inset, b);
    canvas.drawRRect(RRect.fromRectAndRadius(along(0, len, 2), const Radius.circular(4)), Paint()..color = track);
    canvas.drawRRect(
      RRect.fromRectAndRadius(along(0, len, 2), const Radius.circular(4)),
      Paint()
        ..color = border
        ..style = PaintingStyle.stroke,
    );
    final m = marks;
    if (m != null && m.isNotEmpty) {
      final p = Paint()..color = mark;
      var start = -1;
      for (var i = 0; i <= m.length; i++) {
        final on = i < m.length && m[i];
        if (on && start < 0) start = i;
        if (!on && start >= 0) {
          canvas.drawRect(along(start / m.length * len, math.max(start / m.length * len + 1, i / m.length * len), 4), p);
          start = -1;
        }
      }
    }
    final s = maxOffset == 0 ? 0.0 : offset.value / maxOffset * (len - thumb);
    canvas.drawRRect(RRect.fromRectAndRadius(along(s, s + thumb, 1), const Radius.circular(5)), Paint()..color = thumbColor);
  }

  @override
  bool shouldRepaint(_TrackPainter old) => old.thumb != thumb || old.maxOffset != maxOffset || old.thumbColor != thumbColor || old.marks != marks;
}

// ---------------------------------------------------------------------------
// Painters
// ---------------------------------------------------------------------------

class _NamesPainter extends CustomPainter {
  _NamesPainter({required this.seqs, required this.rowH, required this.y, required this.hover, required this.text, required this.band, required this.bandText})
      : super(repaint: Listenable.merge([y, hover]));
  final List<SeqRecord> seqs;
  final double rowH;
  final ValueNotifier<double> y;
  final ValueNotifier<(int, int?)?> hover;
  final Color text, band, bandText;

  @override
  void paint(Canvas canvas, Size size) {
    final y0 = y.value;
    final first = (y0 / rowH).floor().clamp(0, seqs.length);
    final last = ((y0 + size.height) / rowH).ceil().clamp(0, seqs.length);
    final h = hover.value?.$1;
    for (var r = first; r < last; r++) {
      final top = r * rowH - y0;
      final on = r == h;
      if (on) {
        canvas.drawRect(Rect.fromLTWH(0, top, size.width, rowH), Paint()..color = band.withValues(alpha: 0.16));
        canvas.drawRect(Rect.fromLTWH(0, top, 3, rowH), Paint()..color = band);
      }
      final tp = TextPainter(
        text: TextSpan(
          text: seqs[r].header,
          style: TextStyle(fontFamily: monoFont, fontSize: 12.5, color: on ? bandText : text, fontWeight: on ? FontWeight.w700 : FontWeight.w400),
        ),
        textDirection: TextDirection.ltr,
        maxLines: 1,
        ellipsis: '…',
      )..layout(maxWidth: math.max(0, size.width - 14));
      tp.paint(canvas, Offset(8, top + (rowH - tp.height) / 2));
    }
  }

  @override
  bool shouldRepaint(_NamesPainter old) => !identical(old.seqs, seqs) || old.text != text;
}

class _ConservationPainter extends CustomPainter {
  _ConservationPainter(this.values, this.removed, this.cell, this.x, this.viewportWidth, this.color, this.base, this.removedColor) : super(repaint: x);
  final List<double> values;
  final List<bool>? removed;
  final double cell;
  final ValueNotifier<double> x;
  final double viewportWidth;
  final Color color, base, removedColor;

  @override
  void paint(Canvas canvas, Size size) {
    final x0 = x.value;
    final first = (x0 / cell).floor().clamp(0, values.length);
    final last = ((x0 + viewportWidth) / cell).ceil().clamp(0, values.length);
    canvas.drawRect(Rect.fromLTWH(0, size.height - 2, viewportWidth, 1), Paint()..color = base);
    final paint = Paint();
    for (var c = first; c < last; c++) {
      final v = values[c];
      final gone = removed != null && c < removed!.length && removed![c];
      paint.color = gone ? removedColor.withValues(alpha: 0.5) : Color.lerp(base, color, v)!;
      final hgt = (size.height - 4) * v;
      canvas.drawRect(Rect.fromLTWH(c * cell - x0, size.height - 2 - hgt, math.max(1, cell - (cell > 4 ? 1 : 0)), hgt), paint);
    }
  }

  @override
  bool shouldRepaint(_ConservationPainter old) =>
      old.values != values || old.cell != cell || old.viewportWidth != viewportWidth || old.removed != removed || old.color != color;
}

class _GridPainter extends CustomPainter {
  _GridPainter({
    required this.seqs,
    required this.removed,
    required this.geo,
    required this.rowH,
    required this.x,
    required this.y,
    required this.aa,
    required this.gap,
    required this.text,
    required this.hover,
    required this.mark,
    required this.accent,
    required this.removedGap,
    this.dots,
    this.dot = const Color(0xFF9CA3AF),
    this.codonFrame,
    this.codonAa = const [],
  }) : super(repaint: Listenable.merge([x, y, hover]));

  final ValueNotifier<(int, int?)?> hover;
  final Color mark, accent, removedGap, dot;

  /// Columns drawn as dots below the first row (every base the same).
  final List<bool>? dots;

  /// Codon view: codons start at this column; each one is a block of its
  /// three nucleotides. [codonAa] (one string per row) is kept for the
  /// synonymous / non-synonymous variation track.
  final int? codonFrame;
  final List<String> codonAa;

  final List<SeqRecord> seqs;
  final List<bool>? removed;
  final _Geo geo;
  final double rowH;
  double get cell => geo.cell;
  final ValueNotifier<double> x, y;
  final bool aa;
  final Color gap, text;
  static final _glyphs = <String, TextPainter>{};

  TextPainter _glyph(String ch, Color color) => _glyphs.putIfAbsent('$ch${color.toARGB32()}', () {
        return TextPainter(
          text: TextSpan(text: ch, style: TextStyle(fontFamily: monoFont, fontSize: 13, color: color, fontWeight: FontWeight.w500)),
          textDirection: TextDirection.ltr,
        )..layout();
      });

  @override
  void paint(Canvas canvas, Size size) {
    final x0 = x.value, y0 = y.value;
    final firstCol = math.max(0, geo.col(x0));
    final lastCol = geo.col(x0 + size.width) + 1;
    final firstRow = (y0 / rowH).floor().clamp(0, seqs.length);
    final lastRow = ((y0 + size.height) / rowH).ceil().clamp(0, seqs.length);
    final paint = Paint();
    final gapPaint = Paint()..color = gap;
    final removedGapPaint = Paint()..color = removedGap;
    final letters = cell >= 9;
    final w = cell > 4 ? cell - 1 : cell;
    final cut = removed;
    final f = codonFrame;
    final codons = f == null || codonAa.isEmpty ? 0 : codonAa.first.length;
    final dotGlyph = _glyph('.', dot);
    bool isGone(int c) => cut != null && c < cut.length && cut[c];
    for (var r = firstRow; r < lastRow; r++) {
      final s = seqs[r].sequence;
      final top = r * rowH - y0 + 1;
      final end = math.min(lastCol, s.length);
      for (var c = firstCol; c < end; c++) {
        // Codon view: drawn like the DNA view (dots included); the codons
        // are set apart by the gaps of [geo].
        final left = geo.x(c) - x0;
        final gone = isGone(c);
        final ch = s[c].toUpperCase();
        if (ch == '-' || ch == '.') {
          canvas.drawRect(Rect.fromLTWH(left, top + rowH / 2 - 1.5, w, 2), gone ? removedGapPaint : gapPaint);
          continue;
        }
        // Same base in every sequence: a dot below the first row.
        if (r > 0 && dots != null && c < dots!.length && dots![c]) {
          if (letters) dotGlyph.paint(canvas, Offset(left + (cell - dotGlyph.width) / 2, top + (rowH - 2 - dotGlyph.height) / 2));
          continue;
        }
        paint.color = gone ? _removedResidue : ((aa ? aaColors[ch] : ntColors[ch]) ?? const Color(0xFFCBD5E1));
        canvas.drawRect(Rect.fromLTWH(left, top, w, rowH - 2), paint);
        if (letters) {
          final g = _glyph(ch, gone ? const Color(0xFF6B7280) : residueInk((aa ? aaColors[ch] : ntColors[ch]) ?? const Color(0xFFCBD5E1)));
          g.paint(canvas, Offset(left + (cell - g.width) / 2, top + (rowH - 2 - g.height) / 2));
        }
      }
    }
    // Hovered column and cell.
    final hv = hover.value;
    if (hv != null) {
      final (r, c) = hv;
      final band = Paint()..color = accent.withValues(alpha: 0.16);
      canvas.drawRect(Rect.fromLTWH(0, r * rowH - y0, size.width, rowH), band);
      if (c != null) {
        // Codon view: the whole codon is the unit.
        final inCodon = f != null && c >= f && (c - f) ~/ 3 < codons;
        final c0 = inCodon ? f + 3 * ((c - f) ~/ 3) : c;
        final span = inCodon ? 3 * cell : cell;
        canvas.drawRect(Rect.fromLTWH(geo.x(c0) - x0, 0, span, size.height), band);
        canvas.drawRect(
          Rect.fromLTWH(geo.x(c0) - x0 - 0.5, r * rowH - y0 + 0.5, span, rowH - 1),
          Paint()
            ..color = mark
            ..style = PaintingStyle.stroke
            ..strokeWidth = 1.5,
        );
      }
    }
  }

  @override
  bool shouldRepaint(_GridPainter old) =>
      !identical(old.seqs, seqs) ||
      old.geo != geo ||
      old.aa != aa ||
      old.gap != gap ||
      old.removed != removed ||
      old.text != text ||
      !identical(old.dots, dots) ||
      old.codonFrame != codonFrame ||
      !identical(old.codonAa, codonAa);
}

/// Codon view: above each codon, how much it varies (1 − share of the
/// commonest codon), orange when only synonymous codons differ, red when
/// the amino acid changes.
class _CodonVariationPainter extends CustomPainter {
  _CodonVariationPainter(this.values, this.kinds, this.frame, this.geo, this.x, this.viewportWidth, this.base) : super(repaint: x);
  final List<double> values;
  final List<int> kinds;
  final int frame;
  final _Geo geo;
  double get cell => geo.cell;
  final ValueNotifier<double> x;
  final double viewportWidth;
  final Color base;

  @override
  void paint(Canvas canvas, Size size) {
    final x0 = x.value;
    canvas.drawRect(Rect.fromLTWH(0, size.height - 2, viewportWidth, 1), Paint()..color = base);
    final first = math.max(0, (geo.col(x0) - frame) ~/ 3);
    final last = math.min(values.length, (geo.col(x0 + viewportWidth) - frame) ~/ 3 + 1);
    final paint = Paint();
    for (var k = first; k < last; k++) {
      if (kinds[k] == 0) continue;
      paint.color = kinds[k] == 1 ? _synonymous : _nonSynonymous;
      final h = math.max(2.0, (size.height - 4) * values[k]);
      canvas.drawRRect(
        RRect.fromRectAndRadius(Rect.fromLTWH(geo.x(frame + 3 * k) - x0 + 1, size.height - 2 - h, 3 * cell - 3, h), const Radius.circular(2)),
        paint,
      );
    }
  }

  @override
  bool shouldRepaint(_CodonVariationPainter old) =>
      !identical(old.values, values) || old.frame != frame || old.geo != geo || old.viewportWidth != viewportWidth;
}

/// Column numbers above the alignment: a tick per column when zoomed in and
/// a number every 10, 20, 50… columns so that labels never overlap. The
/// hovered column is marked.
class _RulerPainter extends CustomPainter {
  _RulerPainter(this.cols, this.geo, this.x, this.viewportWidth, this.hover, this.color, this.line, this.accent)
      : super(repaint: Listenable.merge([x, hover]));
  final int cols;
  final _Geo geo;
  final double viewportWidth;
  double get cell => geo.cell;
  final ValueNotifier<double> x;
  final ValueNotifier<(int, int?)?> hover;
  final Color color, line, accent;

  TextPainter _label(String s, Color c, {bool bold = false}) => TextPainter(
        text: TextSpan(text: s, style: TextStyle(fontFamily: monoFont, fontSize: 10, color: c, fontWeight: bold ? FontWeight.w700 : FontWeight.w500)),
        textDirection: TextDirection.ltr,
      )..layout();

  @override
  void paint(Canvas canvas, Size size) {
    final x0 = x.value;
    final first = geo.col(x0).clamp(0, cols);
    final last = (geo.col(x0 + viewportWidth) + 1).clamp(0, cols);
    final paint = Paint()
      ..color = line
      ..strokeWidth = 1;
    canvas.drawLine(Offset(0, size.height - 0.5), Offset(viewportWidth, size.height - 0.5), paint);
    // Smallest step whose labels are at least ~48 px apart.
    const steps = [5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000];
    final step = steps.firstWhere((s) => s * cell >= 48, orElse: () => 10000);
    final minor = step >= 10 && cell * step / 10 >= 4 ? step ~/ (step % 10 == 0 ? 10 : 5) : step;
    final hv = hover.value?.$2;
    for (var c = first; c < last; c++) {
      final n = c + 1;
      final px = geo.x(c) + cell / 2 - x0;
      if (n % step == 0 || n == 1) {
        canvas.drawLine(Offset(px, size.height - 7), Offset(px, size.height), paint..color = color);
        if (hv == null || (geo.x(hv) - geo.x(c)).abs() > 30) {
          final tp = _label('$n', color);
          tp.paint(canvas, Offset(px - tp.width / 2, 1));
        }
      } else if (minor > 0 && n % minor == 0) {
        canvas.drawLine(Offset(px, size.height - 3.5), Offset(px, size.height), paint..color = line);
      }
    }
    if (hv != null && hv >= first && hv < last) {
      final px = geo.x(hv) + cell / 2 - x0;
      final tp = _label('${hv + 1}', Colors.white, bold: true);
      final r = RRect.fromRectAndRadius(Rect.fromCenter(center: Offset(px, 7.5), width: tp.width + 8, height: 14), const Radius.circular(4));
      canvas.drawRRect(r, Paint()..color = accent);
      tp.paint(canvas, Offset(px - tp.width / 2, 7.5 - tp.height / 2));
      canvas.drawLine(
        Offset(px, 14.5),
        Offset(px, size.height),
        Paint()
          ..color = accent
          ..strokeWidth = 1.5,
      );
    }
  }

  @override
  bool shouldRepaint(_RulerPainter old) => old.cols != cols || old.geo != geo || old.viewportWidth != viewportWidth || old.color != color;
}

/// Column geometry of the grid. In the codon view each codon is followed by
/// a gap ([gap] px) so the codons read as groups of three; elsewhere every
/// column is [cell] px wide.
@immutable
class _Geo {
  const _Geo(this.cell, [this.frame, this.gap = 0]);
  final double cell;

  /// Column where the first codon starts (codon view), null otherwise.
  final int? frame;
  final double gap;

  bool get spaced => frame != null && gap > 0;

  /// Left edge of column [c].
  double x(int c) => !spaced || c < frame! ? c * cell : c * cell + ((c - frame!) ~/ 3) * gap;

  /// Column at [px] (a gap belongs to the codon before it).
  int col(double px) {
    if (!spaced || px < frame! * cell) return (px / cell).floor();
    final unit = 3 * cell + gap;
    final rel = px - frame! * cell;
    final k = (rel / unit).floor();
    final i = math.min(2, ((rel - k * unit) / cell).floor());
    return frame! + 3 * k + i;
  }

  /// Width of [cols] columns.
  double width(int cols) => cols <= 0 ? 0 : x(cols - 1) + cell;

  @override
  bool operator ==(Object other) => other is _Geo && other.cell == cell && other.frame == frame && other.gap == gap;

  @override
  int get hashCode => Object.hash(cell, frame, gap);
}
