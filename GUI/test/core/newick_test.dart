import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/core/newick.dart';

double _dist(TreeNode root, String a, String b) {
  Map<TreeNode, double> up(TreeNode t) {
    final out = <TreeNode, double>{};
    var d = 0.0;
    for (TreeNode? n = t; n != null; n = n.parent) {
      out[n] = d;
      d += n.length ?? 0;
    }
    return out;
  }

  final ta = root.tips.firstWhere((t) => t.name == a), tb = root.tips.firstWhere((t) => t.name == b);
  final pa = up(ta);
  for (TreeNode? n = tb; n != null; n = n.parent) {
    if (pa.containsKey(n)) return pa[n]! + up(tb)[n]!;
  }
  throw StateError('no path');
}

Set<String> _splits(TreeNode root) {
  final all = root.tips.map((t) => t.name).toSet();
  return {
    for (final n in root.nodes)
      if (!n.isTip && n.parent != null && n.label != null)
        () {
          final side = n.tips.map((t) => t.name).toSet();
          final key = side.contains(all.reduce((x, y) => x.compareTo(y) < 0 ? x : y)) ? all.difference(side) : side;
          return '${(key.toList()..sort()).join(',')}=${n.label}';
        }(),
  };
}

void main() {
  const iq = "(A:0.1,(B:0.2,(C:0.3,'D d':0.4)95/100:0.05)80/90:0.15,E:0.5);";

  test('parses IQ-TREE output with quoted names and SH-aLRT/UFBoot labels', () {
    final t = parseNewick(iq);
    expect(t.tips.map((x) => x.name), ['A', 'B', 'C', 'D d', 'E']);
    final inner = t.nodes.where((n) => !n.isTip && n.parent != null).toList();
    expect(inner.map((n) => n.label), ['80/90', '95/100']);
    expect(inner.map((n) => n.support), [90, 100]);
    expect(toNewick(parseNewick(toNewick(t))), toNewick(t));
  });

  test('rerooting keeps distances and moves supports with their splits', () {
    final t = parseNewick(iq);
    final before = {for (final p in [('A', 'C'), ('B', 'E'), ('D d', 'E'), ('A', 'B')]) p: _dist(t, p.$1, p.$2)};
    final splits = _splits(t);
    final r = rootOnOutgroup(t, ['C', 'D d']);
    expect(r.children, hasLength(2));
    expect(r.children.any((c) => c.tips.map((x) => x.name).toSet().containsAll({'C', 'D d'}) && c.tipCount == 2), isTrue);
    for (final e in before.entries) {
      expect(_dist(r, e.key.$1, e.key.$2), closeTo(e.value, 1e-9), reason: '${e.key}');
    }
    // Same splits with the same supports (the one on the new root's
    // branch is shown on both sides of it).
    expect(_splits(r).containsAll(splits), isTrue);
  });

  test('midpoint root halves the longest path', () {
    final t = parseNewick('((A:1,B:1):1,(C:1,D:7):1);');
    final r = midpointRoot(t);
    final d = _dist(r, 'A', 'D');
    expect(d, closeTo(10, 1e-9));
    double toRoot(String name) {
      var x = 0.0;
      var n = r.tips.firstWhere((t) => t.name == name);
      while (n.parent != null) {
        x += n.length ?? 0;
        n = n.parent!;
      }
      return x;
    }

    expect(toRoot('D'), closeTo(5, 1e-9));
    expect(toRoot('A'), closeTo(5, 1e-9));
  });

  test('ladderize orders by clade size', () {
    final t = parseNewick('((A,(B,C)),D);');
    ladderize(t);
    expect(t.tips.map((x) => x.name), ['D', 'A', 'B', 'C']);
  });
}
