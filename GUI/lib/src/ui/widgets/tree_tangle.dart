part of 'tree_viewer.dart';

// ---------------------------------------------------------------------------
// Two trees face to face (tanglegram)
// ---------------------------------------------------------------------------

/// A tree to compare with the one in the viewer.
class TreeOption {
  const TreeOption({required this.label, required this.load, this.outgroup = const []});

  /// Name shown in the list and above the tree in the figure.
  final String label;

  /// Reads the Newick text.
  final String Function() load;

  /// Tips to root it on (when present; otherwise it is midpoint rooted).
  final List<String> outgroup;
}

/// A piece of text with its own style (tip labels are made of several).
class _Run {
  const _Run(this.text, {this.italic = false, this.bold = false});
  final String text;
  final bool italic, bold;
}

/// Where the drawing goes: the screen (or PNG) canvas, or SVG text. The
/// comparison figure is drawn once for both.
abstract class _Pen {
  void rect(Rect r, Color c, {double radius = 0});
  void path(List<(String, List<Offset>)> cmds, Color c, double width, {bool round = false});

  /// Closed shape filled with [c].
  void fill(List<(String, List<Offset>)> cmds, Color c);
  void circle(Offset c, double r, Color color);

  /// [runs] with their left edge (right edge when [end]) at [at].dx and
  /// their middle (top when [top]) at [at].dy, turned by [angle] (radians)
  /// around [at], on a rounded [background] when given.
  void text(List<_Run> runs, Offset at,
      {required String family,
      required double size,
      required Color color,
      bool end = false,
      bool bold = false,
      bool top = false,
      double angle = 0,
      Color? background});

  static TextSpan span(List<_Run> runs, String family, double size, Color color, bool bold) => TextSpan(
        style: TextStyle(fontFamily: family, fontSize: size, color: color, fontWeight: bold ? FontWeight.w700 : FontWeight.w400),
        children: [
          for (final r in runs)
            TextSpan(text: r.text, style: TextStyle(fontStyle: r.italic ? FontStyle.italic : FontStyle.normal, fontWeight: r.bold ? FontWeight.w700 : null)),
        ],
      );

  static double width(List<_Run> runs, String family, double size, {bool bold = false}) =>
      (TextPainter(text: span(runs, family, size, const Color(0xFF000000), bold), textDirection: TextDirection.ltr)..layout()).width;

  /// [runs] cut at the end with "…" so they are at most [maxW] wide.
  static List<_Run> fit(List<_Run> runs, String family, double size, double maxW, {bool bold = false}) {
    if (width(runs, family, size, bold: bold) <= maxW) return runs;
    final out = [...runs];
    while (out.isNotEmpty) {
      final last = out.removeLast();
      for (var n = last.text.length - 1; n > 0; n--) {
        final cut = [...out, _Run('${last.text.substring(0, n).trimRight()}…', italic: last.italic, bold: last.bold)];
        if (width(cut, family, size, bold: bold) <= maxW) return cut;
      }
    }
    return const [_Run('…')];
  }
}

class _CanvasPen implements _Pen {
  _CanvasPen(this.canvas);
  final Canvas canvas;

  @override
  void rect(Rect r, Color c, {double radius = 0}) =>
      radius > 0 ? canvas.drawRRect(RRect.fromRectAndRadius(r, Radius.circular(radius)), Paint()..color = c) : canvas.drawRect(r, Paint()..color = c);

  @override
  void path(List<(String, List<Offset>)> cmds, Color c, double width, {bool round = false}) {
    canvas.drawPath(
        _path(cmds),
        Paint()
          ..color = c
          ..strokeWidth = width
          ..style = PaintingStyle.stroke
          ..strokeCap = round ? StrokeCap.round : StrokeCap.square
          ..strokeJoin = StrokeJoin.round);
  }

  static Path _path(List<(String, List<Offset>)> cmds) {
    final p = Path();
    for (final (op, pts) in cmds) {
      switch (op) {
        case 'M':
          p.moveTo(pts[0].dx, pts[0].dy);
        case 'L':
          p.lineTo(pts[0].dx, pts[0].dy);
        case 'Q':
          p.quadraticBezierTo(pts[0].dx, pts[0].dy, pts[1].dx, pts[1].dy);
      }
    }
    return p;
  }

  @override
  void fill(List<(String, List<Offset>)> cmds, Color c) => canvas.drawPath(_path(cmds)..close(), Paint()..color = c);

  @override
  void circle(Offset c, double r, Color color) => canvas.drawCircle(c, r, Paint()..color = color);

  @override
  void text(List<_Run> runs, Offset at,
      {required String family,
      required double size,
      required Color color,
      bool end = false,
      bool bold = false,
      bool top = false,
      double angle = 0,
      Color? background}) {
    final tp = TextPainter(text: _Pen.span(runs, family, size, color, bold), textDirection: TextDirection.ltr)..layout();
    final local = Offset(end ? -tp.width : 0, top ? 0 : -tp.height / 2);
    if (angle == 0 && background == null) {
      tp.paint(canvas, at + local);
      return;
    }
    canvas.save();
    canvas.translate(at.dx, at.dy);
    if (angle != 0) canvas.rotate(angle);
    if (background != null) {
      canvas.drawRRect(RRect.fromRectAndRadius(Rect.fromLTWH(local.dx - 2, local.dy, tp.width + 4, tp.height), const Radius.circular(3)), Paint()..color = background);
    }
    tp.paint(canvas, local);
    canvas.restore();
  }
}

class _SvgPen implements _Pen {
  final b = StringBuffer();
  static String f(double v) => v.toStringAsFixed(1);
  static String hex(Color c) => '#${(c.toARGB32() & 0xFFFFFF).toRadixString(16).padLeft(6, '0')}';
  static String esc(String s) => s.replaceAll('&', '&amp;').replaceAll('<', '&lt;').replaceAll('>', '&gt;');
  static String fillAttr(Color c) => 'fill="${hex(c)}"${c.a < 1 ? ' fill-opacity="${c.a.toStringAsFixed(2)}"' : ''}';
  static String stroke(Color c) => 'stroke="${hex(c)}"${c.a < 1 ? ' stroke-opacity="${c.a.toStringAsFixed(2)}"' : ''}';

  @override
  void rect(Rect r, Color c, {double radius = 0}) => b.writeln(
      '<rect x="${f(r.left)}" y="${f(r.top)}" width="${f(r.width)}" height="${f(r.height)}"${radius > 0 ? ' rx="${f(radius)}"' : ''} ${fillAttr(c)}/>');

  @override
  void path(List<(String, List<Offset>)> cmds, Color c, double width, {bool round = false}) {
    b.writeln('<path d="${_d(cmds)}" fill="none" ${stroke(c)} stroke-width="${width.toStringAsFixed(2)}" '
        'stroke-linecap="${round ? 'round' : 'square'}" stroke-linejoin="round"/>');
  }

  static String _d(List<(String, List<Offset>)> cmds) => [for (final (op, pts) in cmds) '$op${pts.map((p) => '${f(p.dx)},${f(p.dy)}').join(' ')}'].join(' ');

  @override
  void fill(List<(String, List<Offset>)> cmds, Color c) => b.writeln('<path d="${_d(cmds)} Z" ${fillAttr(c)}/>');

  @override
  void circle(Offset c, double r, Color color) => b.writeln('<circle cx="${f(c.dx)}" cy="${f(c.dy)}" r="${f(r)}" ${fillAttr(color)}/>');

  @override
  void text(List<_Run> runs, Offset at,
      {required String family,
      required double size,
      required Color color,
      bool end = false,
      bool bold = false,
      bool top = false,
      double angle = 0,
      Color? background}) {
    final turned = angle != 0;
    // Turned text is written around (0, 0) inside a rotated group.
    final o = turned ? Offset.zero : at;
    final y = top ? o.dy + size : o.dy + size * 0.35;
    final spans = [
      for (final r in runs)
        '<tspan${r.italic ? ' font-style="italic"' : ''}${r.bold ? ' font-weight="bold"' : ''}>${esc(r.text)}</tspan>',
    ].join();
    if (turned) b.writeln('<g transform="translate(${f(at.dx)} ${f(at.dy)}) rotate(${(angle * 180 / math.pi).toStringAsFixed(2)})">');
    if (background != null) {
      final w = _Pen.width(runs, family, size, bold: bold);
      final h = size * 1.2;
      final x0 = end ? o.dx - w : o.dx;
      final y0 = top ? o.dy : o.dy - h / 2;
      b.writeln('<rect x="${f(x0 - 2)}" y="${f(y0)}" width="${f(w + 4)}" height="${f(h)}" rx="3" ${fillAttr(background)}/>');
    }
    b.writeln('<text x="${f(o.dx)}" y="${f(y)}" font-size="${size.toStringAsFixed(2)}"${bold ? ' font-weight="bold"' : ''}'
        '${end ? ' text-anchor="end"' : ''} ${fillAttr(color)} xml:space="preserve">$spans</text>');
    if (turned) b.writeln('</g>');
  }
}

/// One support value to place: the node's point [o] (radius [r]).
typedef _SupportItem = ({Offset o, String text, double r, Color color});

/// Draws the support values of [items] so that none overlaps another value
/// or another node's point. Each value starts beside the end of its branch
/// (above or below it, see `supportBelow`) or, with `supportArrow`, moved
/// off the branches; when that place is taken it moves further away
/// (up/down, then to the other side, then back along the branch) and an
/// arrow, straight or curved (`supportCurved`), points from the value to
/// its node. [dir] is -1 on a mirrored tree (values on the right).
void _drawSupports(_Pen pen, _Style st, List<_SupportItem> items, {double dir = 1, Color Function(Offset o)? background}) {
  final h = st.supportFont * 1.2;
  final placed = <Rect>[];
  final points = [for (final it in items) Rect.fromCircle(center: it.o, radius: it.r + 0.5)];
  final order = [for (var i = 0; i < items.length; i++) i]..sort((a, b) {
      final c = items[a].o.dy.compareTo(items[b].o.dy);
      return c != 0 ? c : items[a].o.dx.compareTo(items[b].o.dx);
    });
  for (final i in order) {
    final it = items[i];
    final w = _Pen.width([_Run(it.text)], st.family, st.supportFont, bold: true);
    final vy = st.supportBelow ? 1.0 : -1.0;
    // Box of the value with its corner nearest the node at (cx, cy).
    Rect box(double cx, double cy, double side) {
      final left = dir > 0 ? cx - w : cx;
      final top = side > 0 ? cy : cy - h;
      return Rect.fromLTWH(left, top, w, h);
    }

    final gap = math.max(7.0, st.supportFont * 1.1);
    // Nearest places first: beside the node on the chosen side, on the
    // other side, then further away (each with an arrow).
    final candidates = <(Rect, bool)>[
      if (!st.supportArrow)
        for (final side in [vy, -vy]) (box(it.o.dx - dir * (it.r + 3), it.o.dy + side * (st.lineW / 2 + 1), side), false),
      for (var k = 0; k <= 6; k++)
        for (final back in [0.0, 0.6, 1.2])
          for (final side in [vy, -vy])
            (box(it.o.dx - dir * (it.r + gap * (k == 0 ? 1 : 0.6) + back * w), it.o.dy + side * (it.r + gap * 0.8 + k * h * 0.95), side), true),
    ];
    var chosen = candidates.first;
    for (final c in candidates) {
      final r = c.$1.deflate(0.5);
      final clash = placed.any((p) => p.overlaps(r)) || points.indexed.any((e) => e.$1 != i && e.$2.overlaps(r));
      if (!clash) {
        chosen = c;
        break;
      }
    }
    final rect = chosen.$1;
    placed.add(rect.inflate(1));
    final bg = background?.call(it.o);
    pen.text([_Run(it.text)], Offset(dir > 0 ? rect.right : rect.left, rect.top),
        family: st.family, size: st.supportFont, color: it.color, bold: true, top: true, end: dir > 0, background: bg);
    if (chosen.$2) {
      // From the corner of the value nearest the node.
      final from = Offset(dir > 0 ? rect.right : rect.left, it.o.dy < rect.top ? rect.top : (it.o.dy > rect.bottom ? rect.bottom : rect.center.dy));
      _drawArrow(pen, from, it.o, it.r, it.color, st.lineW, curved: st.supportCurved);
    }
  }
}

/// Thin arrow from [from] to the edge of the point of radius [r] at [to],
/// straight or curved (a quarter bend that leaves the value sideways).
void _drawArrow(_Pen pen, Offset from, Offset to, double r, Color color, double lineW, {bool curved = false}) {
  final v = to - from;
  final len = v.distance;
  if (len < r + 3) return;
  final head = math.max(3.0, lineW * 2.6);
  final width = math.max(0.6, lineW * 0.75);
  if (!curved) {
    final d = v / len;
    final n = Offset(-d.dy, d.dx);
    final tip = to - d * (r + 0.8);
    pen.path([('M', [from + d * 1.5]), ('L', [tip - d * head * 0.6])], color, width, round: true);
    pen.fill([('M', [tip]), ('L', [tip - d * head + n * head * 0.5]), ('L', [tip - d * head - n * head * 0.5])], color);
    return;
  }
  // Control point level with the value and above/below the node: the
  // line leaves the value sideways and comes down (or up) onto the node.
  final ctrl = Offset(to.dx, from.dy);
  final end = to - Offset(0, (to.dy - from.dy).sign) * (r + 0.8);
  final tangent = end - ctrl;
  final tl = tangent.distance;
  if (tl < 1) {
    _drawArrow(pen, from, to, r, color, lineW);
    return;
  }
  final d = tangent / tl;
  final n = Offset(-d.dy, d.dx);
  pen.path([('M', [from]), ('Q', [ctrl, end - d * head * 0.6])], color, width, round: true);
  pen.fill([('M', [end]), ('L', [end - d * head + n * head * 0.5]), ('L', [end - d * head - n * head * 0.5])], color);
}

/// One tree of the comparison: tips top to bottom, growing to the right
/// (left tree) or to the left (right tree, [mirror]).
class _TangleSide {
  _TangleSide(this.root, this.mirror);
  final TreeNode root;
  final bool mirror;
  final rows = <TreeNode>[];
  final x = <TreeNode, double>{}, y = <TreeNode, double>{};
  final labelW = <TreeNode, double>{};
  double origin = 0, scale = 1, maxX = 1, rowH = 10;

  /// Widest a tip name may be (longer ones are cut).
  double labelMax = double.infinity;

  double get dir => mirror ? -1 : 1;
  Offset point(TreeNode n) => Offset(origin + dir * x[n]! * scale, y[n]!);

  /// x of the farthest tip.
  double get end => origin + dir * maxX * scale;

  void place(double Function(TreeNode) length, bool cladogram, double top, double rowHeight) {
    rowH = rowHeight;
    final heights = <TreeNode, int>{};
    int height(TreeNode n) => heights[n] ??= n.isTip ? 0 : 1 + n.children.map(height).reduce(math.max);
    final depth = height(root);
    double walk(TreeNode n, double px) {
      final nx = cladogram ? (depth - height(n)).toDouble() : px + (n.parent == null ? 0 : length(n));
      x[n] = nx;
      if (n.isTip) {
        rows.add(n);
        y[n] = top + (rows.length - 0.5) * rowH;
      } else {
        final ys = [for (final c in n.children) walk(c, nx)];
        y[n] = (ys.first + ys.last) / 2;
      }
      return y[n]!;
    }

    walk(root, 0);
    maxX = x.values.fold(0.0, math.max);
    if (maxX <= 0) maxX = 1;
  }

  /// Branches as drawing commands, with the node whose colour each takes
  /// (as [_Layout.branchPaths], mirrored for the right tree).
  Iterable<(List<(String, List<Offset>)>, TreeNode)> paths(BranchStyle style) sync* {
    for (final n in x.keys) {
      if (n.isTip) continue;
      final o = point(n);
      if (style == BranchStyle.square) {
        final ys = n.children.map((c) => y[c]!).toList();
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
            final r = math.min(rowH * 0.5, math.min(dy.abs(), (e.dx - o.dx).abs()));
            if (r <= 0.5) {
              yield ([('M', [o]), ('L', [Offset(o.dx, e.dy)]), ('L', [e])], c);
            } else {
              yield ([
                ('M', [o]),
                ('L', [Offset(o.dx, e.dy - dy.sign * r)]),
                ('Q', [Offset(o.dx, e.dy), Offset(o.dx + dir * r, e.dy)]),
                ('L', [e]),
              ], c);
            }
        }
      }
    }
  }

  /// Scale bar of a round length about a fifth of the tree, under it.
  (double, double, String) scaleBar() {
    final target = maxX / 5;
    final mag = math.pow(10, (math.log(target) / math.ln10).floor()).toDouble();
    final nice = [1, 2, 5, 10].map((k) => k * mag).lastWhere((v) => v <= target, orElse: () => mag);
    final a = origin, b = origin + dir * nice * scale;
    return (math.min(a, b), math.max(a, b), _Layout._fmt(nice));
  }
}

/// Both trees placed in the figure, with the strip between them where the
/// lines join each tip to its pair.
class _Tangle {
  _Tangle(this.left, this.right, this.size, this.zoneL, this.zoneR, this.headY, this.barY, this.legendY, this.comparison);
  final _TangleSide left, right;
  final Size size;
  final double zoneL, zoneR, headY;
  final double? barY, legendY;
  final TreeComparison comparison;

  /// Key of the tip in the row under [p] (on the side of the figure [p] is on).
  String? hit(Offset p) {
    final side = p.dx < size.width / 2 ? left : right;
    for (final n in side.rows) {
      if ((side.y[n]! - p.dy).abs() <= side.rowH / 2) return tipKey(n.name);
    }
    return null;
  }
}

extension _TangleDrawing on _TreeViewerState {
  List<_Run> _tipRuns(String name) {
    final parts = _parts(name);
    return [
      for (final (i, (text, p)) in parts.indexed) ...[
        if (i > 0) const _Run(' '),
        _Run(text, italic: p.italic, bold: p.bold),
      ],
    ];
  }

  /// The other tree rotated to face the one being edited.
  TreeNode _facing(TreeNode other) {
    final copy = copyTree(other);
    if (_untangle) untangle(copy, {for (final (i, t) in _tree!.tips.indexed) tipKey(t.name): i});
    return copy;
  }

  /// Layout of the comparison: the tree being edited on the left (on the
  /// right with [_swapSides]) and [other] facing it.
  _Tangle _tangle(TreeNode edited, TreeNode other) {
    final st = _style;
    final a = _swapSides ? other : edited, b = _swapSides ? edited : other;
    final cmp = compareTrees(a, b);
    final w = _width;
    final m = st.margin + 4;
    final headH = st.font * 2.2;
    final top = st.margin + 6 + st.titleH + headH;
    final cladogram = _shape == TreeShape.cladogram;
    double raw(TreeNode n) => n.length ?? 0;
    final left = _TangleSide(a, false)..place(_swapSides ? raw : _drawnLength, cladogram, top, _rowH);
    final right = _TangleSide(b, true)..place(_swapSides ? _drawnLength : raw, cladogram, top, _rowH);
    final gap = math.max(24.0, w * _gapPct / 100);
    // Room outside each root for the support values drawn beside it.
    final pad = st.support == TreeSupport.numbers
        ? _Pen.width(const [_Run('100')], st.family, st.supportFont, bold: true) + st.dot + 6 + (st.supportArrow ? math.max(7.0, st.supportFont * 1.1) : 0)
        : 0.0;
    final avail = (w - 2 * m - 2 * pad - gap) / 2;
    // Names take what they need, leaving the tree at least a third of its
    // half; longer names are cut ("…") so they never reach the lines
    // between the trees or the other tree.
    final extra = st.labelGap + (_tipDots ? st.dot : 0) + 4;
    final maxLabel = math.max(20.0, avail * 0.67 - extra);
    for (final s in [left, right]) {
      var lw = 0.0;
      for (final n in s.rows) {
        final v = _showLabels ? math.min(maxLabel, _Pen.width(_tipRuns(n.name), st.family, st.font)) : 0.0;
        s.labelW[n] = v;
        lw = math.max(lw, v);
      }
      s.labelMax = maxLabel;
      s.scale = math.max(10.0, avail - lw - extra) / s.maxX;
    }
    left.origin = m + pad;
    right.origin = w - m - pad;
    final rowsBottom = top + math.max(left.rows.length, right.rows.length) * _rowH;
    final bar = _scaleBar && !cladogram;
    final barY = bar ? rowsBottom + 14 : null;
    var h = rowsBottom + (bar ? 40 : 8);
    double? legendY;
    if (_figLegend && _markConflicts && cmp.common.length >= 4) {
      legendY = h;
      h += st.legendFont * 1.7 + 6;
    }
    return _Tangle(left, right, Size(w, h + m), m + pad + avail, w - m - pad - avail, st.margin + st.titleH + headH * 0.35, barY, legendY, cmp);
  }

  /// Draws the comparison on [pen]; [export] leaves out the hover marks.
  /// [rightTitle] names the other tree (on the left with [_swapSides]).
  void _drawTangle(_Pen pen, _Tangle l, String rightTitle, {required bool export}) {
    final st = _style;
    final family = st.family;
    if (!(export && _transparent)) pen.rect(Offset.zero & l.size, st.background);
    if (st.title.trim().isNotEmpty) {
      pen.text([_Run(st.title)], Offset(l.left.origin, st.margin), family: family, size: st.titleFont, color: st.label, bold: true, top: true);
    }
    final (leftTitle, otherTitle) = _swapSides ? (rightTitle, widget.title) : (widget.title, rightTitle);
    // Each title within its half, so the two never run into each other.
    final titleW = (l.right.origin - l.left.origin) / 2 - 12;
    pen.text(_Pen.fit([_Run(leftTitle)], family, st.font * 1.05, titleW, bold: true), Offset(l.left.origin, l.headY),
        family: family, size: st.font * 1.05, color: st.label, bold: true, top: true);
    pen.text(_Pen.fit([_Run(otherTitle)], family, st.font * 1.05, titleW, bold: true), Offset(l.right.origin, l.headY),
        family: family, size: st.font * 1.05, color: st.label, bold: true, top: true, end: true);

    // Colours of the taxa (same taxon, same colour on both sides).
    final groupColors = <String, Color>{};
    String? groupOf(TreeNode tip) {
      final by = _colorBy;
      if (by == null) return null;
      final v = _rank(by, tip.name);
      if (v.isEmpty) return null;
      groupColors.putIfAbsent(v, () => _groupColors[groupColors.length % _groupColors.length]);
      return v;
    }

    for (final n in [...l.left.rows, ...l.right.rows]) {
      groupOf(n);
    }
    final cmp = l.comparison;
    final hoverKey = export ? null : _tangleHover;
    final accent = SplaceBlue.c500;

    // Lines between the pairs: from the end of each label to the strip,
    // then across it.
    final rightByKey = {for (final n in l.right.rows) tipKey(n.name): n};
    double labelStart(_TangleSide s, TreeNode n) => (_align ? s.end : s.point(n).dx) + s.dir * (st.labelGap + (_tipDots ? st.dot : 0));
    for (final n in l.left.rows) {
      final k = tipKey(n.name);
      final r = rightByKey[k];
      if (r == null) continue;
      final y1 = l.left.y[n]!, y2 = l.right.y[r]!;
      final x1 = labelStart(l.left, n) + l.left.labelW[n]! + 3;
      final x2 = labelStart(l.right, r) - l.right.labelW[r]! - 3;
      final g = groupOf(n);
      final hot = k == hoverKey;
      final c = hot ? accent : (_groupConnectors && g != null ? groupColors[g]! : const Color(0xFF9CA3AF));
      pen.path([
        ('M', [Offset(x1, y1)]),
        ('L', [Offset(l.zoneL, y1)]),
        ('L', [Offset(l.zoneR, y2)]),
        ('L', [Offset(x2, y2)]),
      ], c, hot ? st.lineW + 1.2 : math.max(0.6, st.lineW * 0.8), round: true);
    }

    for (final (side, conflicts, own) in [(l.left, cmp.conflictA, !_swapSides), (l.right, cmp.conflictB, _swapSides)]) {
      // Branches: clade colours (left tree), taxon colours, and the clades
      // the other tree lacks.
      final colors = own ? _branchColors(side.root, const {}, const {}) : const <TreeNode, Color>{};
      final groupBranch = <TreeNode, Color>{};
      if (_groupBranches && _colorBy != null) {
        String? shared(TreeNode n) {
          final v = {for (final t in n.tips) groupOf(t)};
          return v.length == 1 ? v.first : null;
        }

        for (final n in side.x.keys) {
          final g = shared(n);
          if (g != null) groupBranch[n] = groupColors[g]!;
        }
      }
      final round = st.branchStyle != BranchStyle.square;
      for (final (cmds, n) in side.paths(st.branchStyle)) {
        final conflict = _markConflicts && conflicts.contains(n);
        final clade = colors[n] != null && colors[n] != _branchColor ? colors[n] : null;
        final c = conflict ? _conflictColor : (clade ?? groupBranch[n] ?? st.branch);
        pen.path(cmds, c, conflict ? st.lineW * 1.8 + 0.6 : st.lineW, round: round);
      }
      // Labels.
      if (_showLabels) {
        for (final n in side.rows) {
          final k = tipKey(n.name);
          final paired = cmp.common.contains(k);
          final g = groupOf(n);
          final c = k == hoverKey
              ? accent
              : !paired
                  ? _Style.muted
                  : (_groupLabels && g != null ? _darker(groupColors[g]!) : st.label);
          pen.text(_Pen.fit(_tipRuns(n.name), family, st.font, side.labelMax, bold: k == hoverKey), Offset(labelStart(side, n), side.y[n]!),
              family: family, size: st.font, color: c, end: side.mirror, bold: k == hoverKey);
        }
      }
      if (_tipDots) {
        for (final n in side.rows) {
          pen.circle(side.point(n), st.dot, groupBranch[n] ?? st.branch);
        }
      }
      // Support.
      if (st.support != TreeSupport.none) {
        final items = <_SupportItem>[];
        for (final n in side.x.keys) {
          if (n.isTip || n.parent == null || n.label == null) continue;
          final v = n.support;
          if ((v ?? 0) < st.supportMin) continue;
          items.add((o: side.point(n), text: n.label!, r: st.radius(v), color: st.supportColor(v)));
        }
        // On the figure background, so a value never mixes with the
        // branches or the lines between the trees.
        if (st.support == TreeSupport.numbers) {
          _drawSupports(pen, st, items, dir: side.dir, background: (_) => st.background.withValues(alpha: 0.85));
        }
        for (final it in items) {
          _drawMark(pen, it.o, it.r, it.color, st.markShape);
        }
      }
      // Scale bar.
      if (l.barY != null) {
        final (a, b, text) = side.scaleBar();
        final y = l.barY!;
        pen.path([('M', [Offset(a, y)]), ('L', [Offset(b, y)])], st.branch, st.lineW);
        for (final x in [a, b]) {
          pen.path([('M', [Offset(x, y - 4)]), ('L', [Offset(x, y + 4)])], st.branch, st.lineW);
        }
        pen.text([_Run(text)], Offset((a + b) / 2 - _Pen.width([_Run(text)], family, st.legendFont) / 2, y + 6),
            family: family, size: st.legendFont, color: _Style.muted, top: true);
      }
    }

    // Legend: the marked clades and the distance.
    if (l.legendY != null) {
      final y = l.legendY! + st.legendFont * 0.85;
      final x = l.left.origin;
      pen.path([('M', [Offset(x, y)]), ('L', [Offset(x + 18, y)])], _conflictColor, st.lineW * 1.8 + 0.6);
      final l10n = context.l10n;
      final text = '${l10n.treeCompareLegend}   ·   ${l10n.treeCompareRf('${cmp.rf}', cmp.normalizedRf.toStringAsFixed(2))}';
      pen.text([_Run(text)], Offset(x + 24, y), family: family, size: st.legendFont, color: const Color(0xFF4B5563));
    }
  }
}

/// Paints with a function (the editor's canvas: page, figure and marks).
class _FnPainter extends CustomPainter {
  _FnPainter(this.draw);
  final void Function(Canvas canvas, Size size) draw;

  @override
  void paint(Canvas canvas, Size size) => draw(canvas, size);

  @override
  bool shouldRepaint(_FnPainter old) => true;
}
