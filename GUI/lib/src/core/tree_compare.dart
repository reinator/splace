import 'dart:math' as math;

import 'newick.dart';

/// Tip name used to match the tips of two trees: lower case, `_` as space.
String tipKey(String name) => name.toLowerCase().replaceAll('_', ' ').trim();

/// What two trees share: their common tips, their clades (bipartitions,
/// unrooted, over the common tips) and the Robinson–Foulds distance.
class TreeComparison {
  TreeComparison._({
    required this.common,
    required this.onlyA,
    required this.onlyB,
    required this.splitsA,
    required this.splitsB,
    required this.conflictA,
    required this.conflictB,
  });

  /// Tip keys ([tipKey]) in both trees, and in one only.
  final Set<String> common, onlyA, onlyB;

  /// Bipartitions of each tree over the common tips.
  final Set<String> splitsA, splitsB;

  /// Internal nodes whose clade is not in the other tree.
  final Set<TreeNode> conflictA, conflictB;

  int get shared => splitsA.intersection(splitsB).length;

  /// Robinson–Foulds distance: clades in one tree and not in the other.
  int get rf => splitsA.difference(splitsB).length + splitsB.difference(splitsA).length;

  /// Largest possible distance for fully resolved trees of the common tips.
  int get maxRf => math.max(0, 2 * (common.length - 3));

  /// [rf] over [maxRf] (0 = same topology, 1 = no clade in common).
  double get normalizedRf => maxRf == 0 ? 0 : rf / maxRf;
}

/// Key of each internal node's bipartition over [tips] (sorted names of the
/// side without the first tip), or nothing when the split is trivial.
Map<TreeNode, String> _splits(TreeNode root, Set<String> tips) {
  final ref = (tips.toList()..sort()).first;
  final out = <TreeNode, String>{};
  Set<String> walk(TreeNode n) {
    if (n.isTip) {
      final k = tipKey(n.name);
      return tips.contains(k) ? {k} : <String>{};
    }
    final s = <String>{};
    for (final c in n.children) {
      s.addAll(walk(c));
    }
    if (n.parent != null && s.length >= 2 && s.length <= tips.length - 2) {
      final side = s.contains(ref) ? tips.difference(s) : s;
      out[n] = (side.toList()..sort()).join('\u0001');
    }
    return s;
  }

  walk(root);
  return out;
}

/// Compares the topologies of [a] and [b] over their common tips.
TreeComparison compareTrees(TreeNode a, TreeNode b) {
  final ta = {for (final t in a.tips) tipKey(t.name)};
  final tb = {for (final t in b.tips) tipKey(t.name)};
  final common = ta.intersection(tb);
  if (common.length < 4) {
    return TreeComparison._(
      common: common,
      onlyA: ta.difference(common),
      onlyB: tb.difference(common),
      splitsA: {},
      splitsB: {},
      conflictA: {},
      conflictB: {},
    );
  }
  final sa = _splits(a, common), sb = _splits(b, common);
  final setA = sa.values.toSet(), setB = sb.values.toSet();
  return TreeComparison._(
    common: common,
    onlyA: ta.difference(common),
    onlyB: tb.difference(common),
    splitsA: setA,
    splitsB: setB,
    conflictA: {for (final e in sa.entries) if (!setB.contains(e.value)) e.key},
    conflictB: {for (final e in sb.entries) if (!setA.contains(e.value)) e.key},
  );
}

/// Rotates the clades of [tree] (without changing its topology) so its tips
/// come as close as possible to [order] (tip key → row in the other tree):
/// the children of every node are sorted by the mean row of their tips.
/// Tips missing from [order] keep their place among their sisters.
void untangle(TreeNode tree, Map<String, int> order) {
  double? mean(TreeNode n) {
    var sum = 0.0, count = 0;
    for (final t in n.tips) {
      final r = order[tipKey(t.name)];
      if (r != null) {
        sum += r;
        count++;
      }
    }
    return count == 0 ? null : sum / count;
  }

  for (final n in tree.nodes.toList()) {
    if (n.isTip) continue;
    final keyed = [for (final (i, c) in n.children.indexed) (c, mean(c) ?? double.infinity, i)];
    keyed.sort((x, y) {
      final d = x.$2.compareTo(y.$2);
      return d != 0 ? d : x.$3.compareTo(y.$3);
    });
    n.children
      ..clear()
      ..addAll([for (final k in keyed) k.$1]);
  }
}

/// Pairs of tips drawn with crossing connectors: tips of the left order
/// whose rows in the right order are not increasing.
int crossings(List<String> left, List<String> right) {
  final pos = {for (final (i, k) in right.indexed) k: i};
  final seq = [for (final k in left) ?pos[k]];
  var n = 0;
  for (var i = 0; i < seq.length; i++) {
    for (var j = i + 1; j < seq.length; j++) {
      if (seq[i] > seq[j]) n++;
    }
  }
  return n;
}
