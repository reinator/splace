part of 'tree_viewer.dart';

// ---------------------------------------------------------------------------
// Circular tree
// ---------------------------------------------------------------------------

/// A tree drawn around a centre: the root in the middle, each row (tip or
/// collapsed clade) at its own angle, branches growing outwards. Labels,
/// taxon ring and clade bars sit outside the tips.
class _Circle {
  _Circle(this.root, this.rows, this.angle, this.radius, this.center, this.scale, this.step, this.tipR, this.labelR, this.ringR, this.barR,
      this.barW, this.size, this.top, this.left, this.maxR, this.collapsed, this.bars, this.legendBox, this.legendOrigin, this.barY);

  final TreeNode root;
  final List<TreeNode> rows;

  /// Angle (radians) and distance from the root (tree units) of each node.
  final Map<TreeNode, double> angle, radius;
  final Offset center;

  /// Pixels per tree unit; angle between two rows.
  final double scale, step;

  /// Radius of the farthest tip, where aligned labels start, of the taxon
  /// ring and of the first column of clade bars (pixels).
  final double tipR, labelR, ringR, barR;

  /// Width of the ring (0 without taxon colours).
  final double barW;
  final Size size;
  final double top, left, maxR;
  final Set<TreeNode> collapsed;
  final List<({TreeNode node, _CladeMark mark, double r})> bars;
  final _LegendBox? legendBox;
  final Offset legendOrigin;
  final double? barY;

  bool hidden(TreeNode n) => !angle.containsKey(n);

  Offset polar(double a, double r) => center + Offset(math.cos(a), math.sin(a)) * r;

  /// Pixels from the centre to node [n].
  double r(TreeNode n) => radius[n]! * scale;

  Offset point(TreeNode n) => polar(angle[n]!, r(n));

  /// Angles of the first and last row of clade [n].
  (double, double)? span(TreeNode n) {
    if (hidden(n)) return null;
    if (n.isTip || collapsed.contains(n)) return (angle[n]!, angle[n]!);
    double? a, b;
    for (final k in n.nodes) {
      if ((k.isTip || collapsed.contains(k)) && angle.containsKey(k)) {
        final v = angle[k]!;
        a = a == null ? v : math.min(a, v);
        b = b == null ? v : math.max(b, v);
      }
    }
    return a == null ? null : (a, b!);
  }

  /// Points along the arc of radius [r] from [a0] to [a1].
  List<(String, List<Offset>)> arc(double r, double a0, double a1, {bool move = true}) {
    final n = math.max(1, ((a1 - a0).abs() / (math.pi / 90)).ceil());
    return [
      for (var i = 0; i <= n; i++) (i == 0 && move ? 'M' : 'L', [polar(a0 + (a1 - a0) * i / n, r)]),
    ];
  }

  /// Ring sector between radii [r0] and [r1] and angles [a0] and [a1].
  List<(String, List<Offset>)> sector(double r0, double r1, double a0, double a1) => [
        ...arc(r1, a0, a1),
        ...arc(r0, a1, a0, move: false),
      ];

  /// Text turned along the radius at [a]: on the left half it is flipped
  /// so it still reads from left to right, anchored at its end.
  static (double, bool) turn(double a) {
    final c = math.cos(a);
    return c >= 0 ? (a, false) : (a + math.pi, true);
  }

  /// Node under [p]: the radial branch (or the arc of a node) within a few
  /// pixels, or the row whose label is under it.
  TreeNode? hit(Offset p) {
    final d = p - center;
    final pr = d.distance;
    var pa = math.atan2(d.dy, d.dx);
    final reach = math.max(5.0, math.min(9.0, tipR * step / 2));
    TreeNode? best;
    var bestD = double.infinity;
    double near(double a) {
      // Same angle as [a], up to whole turns.
      var x = pa;
      while (x - a > math.pi) {
        x -= 2 * math.pi;
      }
      while (a - x > math.pi) {
        x += 2 * math.pi;
      }
      return x;
    }

    for (final n in angle.keys) {
      final a = angle[n]!;
      final x = near(a);
      final r1 = n.isTip || collapsed.contains(n) ? ringR : r(n);
      final r0 = n.parent == null ? 0.0 : r(n.parent!);
      // Distance to the radial branch.
      if (pr >= r0 - 2 && pr <= r1 + 2) {
        final dist = (pr * math.sin(x - a)).abs();
        final within = math.cos(x - a) > 0 && (n.isTip || collapsed.contains(n) ? (x - a).abs() <= step / 2 : true);
        if (within && dist < reach && dist < bestD) {
          bestD = dist;
          best = n;
        }
      }
      // Distance to the arc joining the children.
      if (!n.isTip && !collapsed.contains(n)) {
        final sp = span(n);
        if (sp != null) {
          final xa = near((sp.$1 + sp.$2) / 2);
          if (xa >= sp.$1 - 0.01 && xa <= sp.$2 + 0.01) {
            final dist = (pr - r(n)).abs();
            if (dist < reach && dist < bestD) {
              bestD = dist;
              best = n;
            }
          }
        }
      }
    }
    return best;
  }
}

extension _CircleDrawing on _TreeViewerState {
  List<_Run> _rowRuns(TreeNode n) => _collapsed.contains(n) ? [_Run(_rowText(n))] : _tipRuns(n.name);

  _Circle _circle(TreeNode tree, Map<TreeNode, String> groupOf, Map<String, Color> groupColors) {
    final st = _style;
    final w = _width;
    final m = st.margin + 4;
    final top = st.margin + 6 + st.titleH;
    final rows = <TreeNode>[];
    final angle = <TreeNode, double>{}, radius = <TreeNode, double>{};
    final heights = <TreeNode, int>{};
    int height(TreeNode n) => heights[n] ??= n.isTip || _collapsed.contains(n) ? 0 : 1 + n.children.map(height).reduce(math.max);
    final depth = height(tree);
    final cladogram = _shape == TreeShape.cladogram;
    final n = _rows(tree).length;
    final arc = _arcDeg.clamp(30, 360) * math.pi / 180;
    // A full circle leaves one row of space between the last and the first.
    final step = arc / (arc >= 2 * math.pi - 1e-6 ? n : math.max(1, n - 1));
    final start = (_rotateDeg - 90) * math.pi / 180;
    double place(TreeNode k, double pr) {
      final r = cladogram ? (depth - height(k)).toDouble() : pr + (k.parent == null ? 0 : _drawnLength(k));
      radius[k] = r;
      if (k.isTip || _collapsed.contains(k)) {
        angle[k] = start + rows.length * step;
        rows.add(k);
      } else {
        final as = [for (final c in k.children) place(c, r)];
        angle[k] = (as.first + as.last) / 2;
      }
      return angle[k]!;
    }

    place(tree, 0);
    var maxR = radius.values.fold<double>(0, math.max);
    for (final k in _collapsed) {
      if (!radius.containsKey(k)) continue;
      double deepest(TreeNode x, double acc) => x.isTip ? acc : x.children.map((c) => deepest(c, acc + (cladogram ? 1 : _drawnLength(c)))).reduce(math.max);
      maxR = math.max(maxR, radius[k]! + deepest(k, 0));
    }
    if (maxR <= 0) maxR = 1;

    // Room outside the tips: labels, taxon ring, clade bars.
    var labelW = 0.0;
    if (_showLabels) {
      for (final r in rows) {
        labelW = math.max(labelW, _Pen.width(_rowRuns(r), st.family, st.font));
      }
    }
    labelW = math.min(labelW + st.labelGap + (_tipDots ? st.dot : 0) + 4, w * 0.22);
    final ringW = st.groupRows && groupOf.isNotEmpty ? math.max(6.0, st.font * 0.9) : 0.0;
    final shown = [for (final b in _bars(tree, groupColors)) if (angle.containsKey(b.$1)) b];
    final level = <TreeNode, int>{};
    for (final b in [...shown]..sort((a, b) => a.$1.tipCount.compareTo(b.$1.tipCount))) {
      var lv = 0;
      for (final o in shown) {
        if (!identical(o.$1, b.$1) && level.containsKey(o.$1) && b.$1.nodes.contains(o.$1)) lv = math.max(lv, level[o.$1]! + 1);
      }
      level[b.$1] = lv;
    }
    final colW = <int, double>{};
    for (final (k, mk) in shown) {
      final tw = _Pen.width([_Run(mk.text, italic: mk.italic)], st.family, st.font, bold: true);
      colW[level[k]!] = math.max(colW[level[k]!] ?? 0, tw + 16);
    }
    final barsW = colW.values.fold(0.0, (a, b) => a + b);
    final outer = (w - 2 * m) / 2;
    final tipR = math.max(30.0, outer - labelW - (ringW > 0 ? ringW + 4 : 0) - barsW);
    final scale = tipR / maxR;
    final labelR = tipR;
    final ringR = tipR + labelW + (ringW > 0 ? 2 : 0);
    final barR = ringR + (ringW > 0 ? ringW + 4 : 0) + 4;
    double colR(int lv) => barR + [for (var i = 0; i < lv; i++) colW[i] ?? 0].fold(0.0, (a, b) => a + b);
    final center = Offset(w / 2, top + outer);
    final circleBottom = top + 2 * outer;
    final bar = _scaleBar && _barLength >= 0 && !cladogram;
    final barY = bar ? circleBottom + 14 : null;
    var h = circleBottom + (bar ? 40 : 8);
    final legend = _figLegend ? _legendLines(groupColors) : const <List<_LegendItem>>[];
    _LegendBox? box;
    var origin = Offset.zero;
    if (legend.isNotEmpty) {
      box = _LegendBox.compute(legend, st, vertical: false, maxWidth: w - 2 * m);
      origin = Offset(m, h);
      h += box.size.height + 8;
    }
    return _Circle(tree, rows, angle, radius, center, scale, step, tipR, labelR, ringR, barR, ringW, Size(w, h + st.margin), top, m, maxR,
        _collapsed, [for (final (k, mk) in shown) (node: k, mark: mk, r: colR(level[k]!))], box, origin, barY);
  }

  /// Draws the circular figure on [pen]; [export] leaves out the
  /// selection, hover and search marks.
  void _drawCircle(_Pen pen, _Circle l, Map<TreeNode, String> groupOf, Map<String, Color> groupColors, {required bool export}) {
    final st = _style;
    final family = st.family;
    final accent = SplaceBlue.c500;
    final selected = export ? null : _selected;
    final hover = export ? null : _hover;
    final query = export ? '' : _query.trim().toLowerCase();
    final inSel = <TreeNode>{if (selected != null) ...selected.nodes};
    final colors = _branchColors(l.root, groupOf, groupColors);
    final outgroup = {for (final o in widget.outgroup) o.toLowerCase().replaceAll('_', ' ').trim()};
    final half = l.step / 2;
    final labelEnd = l.ringR;

    if (!(export && _transparent)) pen.rect(Offset.zero & l.size, st.background);
    if (st.title.trim().isNotEmpty) {
      pen.text([_Run(st.title)], Offset(l.left, st.margin), family: family, size: st.titleFont, color: st.label, bold: true, top: true);
    }
    // Coloured sectors behind clades (larger first).
    for (final e in _highlights.entries.toList()..sort((a, b) => b.key.tipCount.compareTo(a.key.tipCount))) {
      final sp = l.span(e.key);
      if (sp == null) continue;
      final r0 = e.key.parent == null ? 0.0 : (l.r(e.key) + l.r(e.key.parent!)) / 2;
      pen.fill(l.sector(r0, labelEnd, sp.$1 - half * 0.9, sp.$2 + half * 0.9), e.value.withValues(alpha: st.highlightAlpha));
    }
    if (selected != null) {
      final sp = l.span(selected);
      if (sp != null) pen.fill(l.sector(math.max(0, l.r(selected) - 6), labelEnd, sp.$1 - half, sp.$2 + half), accent.withValues(alpha: 0.10));
    }
    if (hover != null && (hover.isTip || l.collapsed.contains(hover)) && !l.hidden(hover)) {
      final a = l.angle[hover]!;
      pen.fill(l.sector(l.r(hover), labelEnd, a - half, a + half), accent.withValues(alpha: 0.10));
    }
    // Taxon ring outside the labels.
    if (l.barW > 0) {
      String? current;
      double? a0;
      double prev = 0;
      void close() {
        if (current != null) pen.fill(l.sector(l.ringR, l.ringR + l.barW, a0! - half, prev + half), groupColors[current]!);
      }

      for (final n in l.rows) {
        final g = groupOf[n];
        final a = l.angle[n]!;
        if (g != current) {
          close();
          current = g;
          a0 = a;
        }
        prev = a;
      }
      close();
    }
    // Branches: an arc joining the children of each node, a radial line to
    // each child (or a straight line, slanted style).
    final round = st.branchStyle != BranchStyle.square;
    for (final n in l.angle.keys) {
      if (n.isTip || l.collapsed.contains(n)) continue;
      final rn = l.r(n);
      final as = [for (final c in n.children) l.angle[c]!];
      if (st.branchStyle != BranchStyle.slanted && as.length > 1) {
        pen.path(l.arc(rn, as.reduce(math.min), as.reduce(math.max)), inSel.contains(n) ? accent : (colors[n] ?? st.branch),
            inSel.contains(n) ? st.lineW + 1.2 : st.lineW, round: true);
      }
      for (final c in n.children) {
        final from = st.branchStyle == BranchStyle.slanted ? l.point(n) : l.polar(l.angle[c]!, rn);
        pen.path([('M', [from]), ('L', [l.point(c)])], inSel.contains(c) ? accent : (colors[c] ?? st.branch), inSel.contains(c) ? st.lineW + 1.2 : st.lineW,
            round: round);
      }
    }
    // Collapsed clades: a wedge out to the deepest tip.
    final end = l.maxR * l.scale;
    for (final n in l.collapsed) {
      if (l.hidden(n)) continue;
      final a = l.angle[n]!;
      final c = inSel.contains(n) ? accent : (colors[n] ?? st.branch);
      final shape = [('M', [l.point(n)]), ('L', [l.polar(a - half * 0.84, end)]), ...l.arc(end, a - half * 0.84, a + half * 0.84, move: false)];
      pen.fill(shape, c.withValues(alpha: 0.25));
      pen.path([...shape, ('L', [l.point(n)])], c, 1);
    }
    // Shortened branches: // across the branch.
    if (_shape == TreeShape.phylogram) {
      for (final n in _truncated) {
        if (l.hidden(n) || n.parent == null) continue;
        final a = l.angle[n]!;
        final mid = (l.r(n) + l.r(n.parent!)) / 2;
        final along = Offset(math.cos(a), math.sin(a)), across = Offset(-math.sin(a), math.cos(a));
        for (final k in [-2.5, 2.5]) {
          final o = l.polar(a, mid + k);
          pen.path([('M', [o - across * 5 - along * 2]), ('L', [o + across * 5 + along * 2])], st.branch, math.max(1.2, st.lineW * 0.8));
        }
      }
    }
    // Labels, turned along their radius.
    if (_showLabels) {
      for (final n in l.rows) {
        final collapsed = l.collapsed.contains(n);
        final a = l.angle[n]!;
        final from = (collapsed || _align ? l.labelR : l.r(n)) + st.labelGap + (_tipDots && n.isTip && !_align ? st.dot : 0);
        if (_align && !collapsed && l.labelR - l.r(n) > 6) {
          for (var r = l.r(n) + 4; r < l.labelR + st.labelGap - 3; r += 4) {
            pen.path([('M', [l.polar(a, r)]), ('L', [l.polar(a, math.min(r + 1.5, l.labelR + st.labelGap - 3))])], _Style.guide, 1);
          }
        }
        final (turned, flip) = _Circle.turn(a);
        final hit = query.isNotEmpty && _rowText(n).toLowerCase().contains(query);
        final color = inSel.contains(n) || hit
            ? accent
            : (_groupLabels && groupOf[n] != null ? _darker(groupColors[groupOf[n]]!) : st.label);
        final out = _markOutgroup && n.isTip && outgroup.contains(n.name.toLowerCase().replaceAll('_', ' ').trim());
        // Long names are cut ("…") at the room the layout keeps for them,
        // so they never leave the figure.
        final bold = identical(n, hover) || inSel.contains(n) || hit;
        pen.text(_Pen.fit(_rowRuns(n), family, st.font, math.max(12.0, l.ringR - from - 2), bold: bold), l.polar(a, from),
            family: family,
            size: st.font,
            color: color,
            end: flip,
            angle: turned,
            bold: bold,
            background: out ? _Style.outgroup : null);
      }
    }
    // Bars with names beside clades: an arc and the name along the radius.
    for (final b in l.bars) {
      final sp = l.span(b.node);
      if (sp == null) continue;
      final a0 = sp.$1 - half * 0.76, a1 = sp.$2 + half * 0.76;
      pen.path(l.arc(b.r, a0, a1), b.mark.color, st.barStroke);
      final mid = (a0 + a1) / 2;
      final (turned, flip) = _Circle.turn(mid);
      pen.text([_Run(b.mark.text, italic: b.mark.italic)], l.polar(mid, b.r + 6),
          family: family, size: st.font, color: b.mark.color, bold: true, end: flip, angle: turned);
    }
    // Points on the tips.
    if (_tipDots) {
      for (final n in l.rows) {
        if (n.isTip) pen.circle(l.point(n), st.dot, inSel.contains(n) ? accent : (colors[n] ?? st.branch));
      }
    }
    // Support values beside the branch of each node.
    if (st.support != TreeSupport.none) {
      for (final n in l.angle.keys) {
        if (n.isTip || n.parent == null || n.label == null) continue;
        final v = n.support;
        if ((v ?? 0) < st.supportMin) continue;
        final color = st.supportColor(v);
        final o = l.point(n);
        if (st.support == TreeSupport.numbers) {
          final a = l.angle[n]!;
          final (turned, flip) = _Circle.turn(a);
          final unit = Offset(-math.sin(a), math.cos(a)) * ((st.supportBelow ? 1 : -1) * (flip ? -1 : 1));
          // With an arrow the value moves farther off the branch and back
          // towards the centre; the arrow joins its end to the node.
          final gap = st.supportArrow ? math.max(7.0, st.supportFont * 1.1) : 0.0;
          final at = l.polar(a, l.r(n) - 4 - st.radius(v) - gap) + unit * (st.supportFont * 0.65 + st.lineW / 2 + gap * 0.8);
          pen.text([_Run(n.label!)], at,
              family: family, size: st.supportFont, color: color, bold: true, end: !flip, angle: turned, background: st.background.withValues(alpha: 0.85));
          if (st.supportArrow) _drawArrow(pen, at, o, st.radius(v), color, st.lineW, curved: st.supportCurved);
        }
        _drawMark(pen, o, st.radius(v), color, st.markShape);
      }
    }
    // Branch lengths along the branches.
    if (_lengths) {
      for (final n in l.angle.keys) {
        final len = n.length;
        if (n.parent == null || len == null) continue;
        final a = l.angle[n]!;
        final (turned, flip) = _Circle.turn(a);
        final mid = (l.r(n) + l.r(n.parent!)) / 2;
        final w = _Pen.width([_Run(_lengthText(len))], family, math.max(7, st.supportFont - 1));
        final side = Offset(-math.sin(a), math.cos(a)) * ((flip ? 1 : -1) * (st.supportFont * 0.6));
        pen.text([_Run(_lengthText(len))], l.polar(a, mid + (flip ? w / 2 : -w / 2)) + side,
            family: family, size: math.max(7, st.supportFont - 1), color: _Style.muted, angle: turned, end: false);
      }
    }
    // Hovered node.
    if (hover != null && !l.hidden(hover)) {
      final o = l.point(hover);
      final r = st.dot + 2.5;
      pen.path([for (var i = 0; i <= 24; i++) (i == 0 ? 'M' : 'L', [o + Offset(math.cos(i * math.pi / 12), math.sin(i * math.pi / 12)) * r])], accent, 2);
    }
    // Scale bar under the circle.
    if (l.barY != null) {
      var nice = _barLength;
      if (nice <= 0) {
        final target = l.maxR / 5;
        final mag = math.pow(10, (math.log(target) / math.ln10).floor()).toDouble();
        nice = [1, 2, 5, 10].map((k) => k * mag).lastWhere((v) => v <= target, orElse: () => mag);
      }
      final y = l.barY!;
      final a = l.left, b = l.left + nice * l.scale;
      pen.path([('M', [Offset(a, y)]), ('L', [Offset(b, y)])], st.branch, st.lineW);
      for (final x in [a, b]) {
        pen.path([('M', [Offset(x, y - 4)]), ('L', [Offset(x, y + 4)])], st.branch, st.lineW);
      }
      final text = _Layout._fmt(nice);
      pen.text([_Run(text)], Offset((a + b) / 2 - _Pen.width([_Run(text)], family, st.legendFont) / 2, y + 6),
          family: family, size: st.legendFont, color: _Style.muted, top: true);
    }
    // Legend.
    final lg = l.legendBox;
    if (lg != null) {
      final o = l.legendOrigin;
      if (st.legendBox) {
        final r = o & lg.size;
        pen.rect(r, const Color(0xFFD1D5DB), radius: 6);
        pen.rect(r.deflate(1), st.background, radius: 5);
      }
      if (lg.title.isNotEmpty) {
        pen.text([_Run(lg.title)], o + Offset(lg.pad, lg.pad), family: family, size: st.legendFont + 1, color: st.label, bold: true, top: true);
      }
      for (final (it, at) in lg.items) {
        final p = o + at;
        var x = p.dx;
        final mid = p.dy + st.legendFont * 0.65;
        if (it.symbol != null) {
          pen.text([_Run(it.symbol!)], Offset(x, p.dy), family: family, size: st.legendFont, color: st.branch, bold: true, top: true);
          x += st.legendFont * 1.5;
        } else if (it.color != null) {
          if (it.box) {
            pen.rect(Rect.fromLTWH(x, mid - 5, 13, 10), it.color!, radius: 2);
          } else {
            pen.circle(Offset(x + 5, mid), 4.5, it.color!);
          }
          x += 17;
        }
        pen.text([_Run(it.label, italic: it.italic)], Offset(x, p.dy),
            family: family, size: st.legendFont, color: it.heading ? st.label : const Color(0xFF4B5563), bold: it.heading, top: true);
      }
    }
  }
}
