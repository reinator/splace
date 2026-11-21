/// Newick trees for the tree viewer: parsing, rerooting (on a branch, on
/// outgroups or at the midpoint), ladderizing and writing back.
library;


class TreeNode {
  TreeNode({this.name = '', this.length, this.label});

  /// Tip name (empty for internal nodes).
  String name;

  /// Length of the branch above this node.
  double? length;

  /// Internal-node label: support values such as `100`, `85.2/97` or
  /// `0.98`.
  String? label;
  TreeNode? parent;
  final List<TreeNode> children = [];

  bool get isTip => children.isEmpty;

  void add(TreeNode c) {
    c.parent = this;
    children.add(c);
  }

  Iterable<TreeNode> get tips sync* {
    if (isTip) {
      yield this;
    } else {
      for (final c in children) {
        yield* c.tips;
      }
    }
  }

  Iterable<TreeNode> get nodes sync* {
    yield this;
    for (final c in children) {
      yield* c.nodes;
    }
  }

  int get tipCount => isTip ? 1 : children.fold(0, (a, c) => a + c.tipCount);

  /// Support of an internal node: the last number of its label (UFBoot in
  /// IQ-TREE's `SH-aLRT/UFBoot`), on a 0–100 scale.
  double? get support {
    final l = label;
    if (l == null || l.isEmpty) return null;
    final parts = l.split('/');
    final v = double.tryParse(parts.last);
    if (v == null) return null;
    // Posterior-like values (0–1) are shown in percent.
    return v <= 1 && !l.contains('.') ? v : (v <= 1.0 ? v * 100 : v);
  }
}

class NewickException implements Exception {
  NewickException(this.message);
  final String message;
  @override
  String toString() => 'Newick: $message';
}

/// Parses one Newick tree (quoted names, comments in brackets, branch
/// lengths and internal labels).
TreeNode parseNewick(String text) {
  var i = 0;
  final s = text.trim();

  void skipSpace() {
    while (i < s.length) {
      if (s[i] == '[') {
        final end = s.indexOf(']', i);
        i = end < 0 ? s.length : end + 1;
      } else if (' \t\r\n'.contains(s[i])) {
        i++;
      } else {
        break;
      }
    }
  }

  String readName() {
    skipSpace();
    if (i < s.length && s[i] == "'") {
      final b = StringBuffer();
      i++;
      while (i < s.length) {
        if (s[i] == "'") {
          if (i + 1 < s.length && s[i + 1] == "'") {
            b.write("'");
            i += 2;
            continue;
          }
          i++;
          break;
        }
        b.write(s[i++]);
      }
      return b.toString();
    }
    final start = i;
    while (i < s.length && !'(),:;['.contains(s[i])) {
      i++;
    }
    return s.substring(start, i).trim();
  }

  double? readLength() {
    skipSpace();
    if (i < s.length && s[i] == ':') {
      i++;
      skipSpace();
      final start = i;
      while (i < s.length && RegExp(r'[0-9eE.+-]').hasMatch(s[i])) {
        i++;
      }
      return double.tryParse(s.substring(start, i));
    }
    return null;
  }

  TreeNode node() {
    skipSpace();
    final n = TreeNode();
    if (i < s.length && s[i] == '(') {
      i++;
      while (true) {
        n.add(node());
        skipSpace();
        if (i >= s.length) throw NewickException('unexpected end');
        if (s[i] == ',') {
          i++;
          continue;
        }
        if (s[i] == ')') {
          i++;
          break;
        }
        throw NewickException('unexpected "${s[i]}" at $i');
      }
      final label = readName();
      if (label.isNotEmpty) n.label = label;
    } else {
      n.name = readName();
    }
    n.length = readLength();
    return n;
  }

  final root = node();
  skipSpace();
  if (i < s.length && s[i] != ';') throw NewickException('text after the tree at $i');
  return root;
}

String _quote(String name) => RegExp(r"[\s(),:;\[\]']").hasMatch(name) ? "'${name.replaceAll("'", "''")}'" : name;

String _len(double v) => v == 0 ? '0' : v.toStringAsPrecision(10).replaceFirst(RegExp(r'\.?0+$'), '');

/// Writes [root] as Newick.
String toNewick(TreeNode root) {
  String w(TreeNode n) {
    final b = StringBuffer();
    if (n.isTip) {
      b.write(_quote(n.name));
    } else {
      b
        ..write('(')
        ..write(n.children.map(w).join(','))
        ..write(')');
      if (n.label != null) b.write(n.label);
    }
    if (n.length != null && n.parent != null) b.write(':${_len(n.length!)}');
    return b.toString();
  }

  return '${w(root)};';
}

/// Roots the tree on the branch above [target]: the new root has [target]
/// on one side and the rest of the tree on the other, splitting the branch
/// length in half. Lengths and support labels belong to branches (splits),
/// so they move with them when the branches along the path to the old root
/// change direction.
TreeNode rerootAbove(TreeNode root, TreeNode target) {
  if (identical(target, root)) return root;
  final path = <TreeNode>[];
  for (TreeNode? n = target; n != null; n = n.parent) {
    path.add(n);
  }
  // Branch j joins path[j] and path[j + 1].
  final lengths = [for (var j = 0; j < path.length - 1; j++) path[j].length];
  final labels = [for (var j = 0; j < path.length - 1; j++) path[j].label];
  for (var j = 0; j < path.length - 1; j++) {
    path[j + 1].children.remove(path[j]);
  }
  final newRoot = TreeNode();
  final half = (lengths[0] ?? 0) / 2;
  target.length = half;
  newRoot.add(target);
  final p1 = path[1]
    ..length = half
    ..label = labels[0];
  newRoot.add(p1);
  for (var j = 1; j < path.length - 1; j++) {
    final child = path[j + 1]
      ..length = lengths[j]
      ..label = labels[j];
    path[j].add(child);
  }
  _collapseUnary(newRoot);
  return newRoot;
}

/// Removes nodes with a single child (left by rerooting), joining lengths.
void _collapseUnary(TreeNode root) {
  for (final n in root.nodes.toList()) {
    if (n.children.length == 1 && n.parent != null) {
      final c = n.children.single;
      final p = n.parent!;
      final idx = p.children.indexOf(n);
      c.length = (c.length ?? 0) + (n.length ?? 0);
      c.label ??= n.label;
      p.children[idx] = c;
      c.parent = p;
    }
  }
}

/// Roots on the smallest clade holding all [outgroup] tips (or on the
/// branch above the single outgroup tip). When the outgroup tips include
/// the current root's split, the tree is first rooted on another tip.
TreeNode rootOnOutgroup(TreeNode root, Iterable<String> outgroup) {
  final names = outgroup.toSet();
  final tips = root.tips.where((t) => names.contains(t.name)).toList();
  if (tips.isEmpty) return root;
  var r = root;
  TreeNode? mrca(TreeNode tree) {
    final want = tree.tips.where((t) => names.contains(t.name)).toSet();
    TreeNode? best;
    for (final n in tree.nodes) {
      final ts = n.tips.toSet();
      if (ts.containsAll(want) && (best == null || ts.length < best.tipCount)) best = n;
    }
    return best;
  }

  var m = mrca(r);
  if (m == null) return r;
  if (identical(m, r)) {
    // The outgroup spans the root: root first on an ingroup tip.
    final ingroup = r.tips.firstWhere((t) => !names.contains(t.name), orElse: () => r.tips.first);
    r = rerootAbove(r, ingroup);
    m = mrca(r);
    if (m == null || identical(m, r)) return r;
  }
  return rerootAbove(r, m);
}

/// Midpoint rooting: roots halfway along the longest tip-to-tip path.
TreeNode midpointRoot(TreeNode root) {
  final tips = root.tips.toList();
  if (tips.length < 3) return root;
  // Distances between every pair of tips via their ancestors.
  Map<TreeNode, double> toRoot(TreeNode t) {
    final out = <TreeNode, double>{};
    var d = 0.0;
    for (TreeNode? n = t; n != null; n = n.parent) {
      out[n] = d;
      d += n.length ?? 0;
    }
    return out;
  }

  var best = 0.0;
  TreeNode? a, b, lca;
  final paths = {for (final t in tips) t: toRoot(t)};
  for (var i = 0; i < tips.length; i++) {
    for (var j = i + 1; j < tips.length; j++) {
      final pi = paths[tips[i]]!, pj = paths[tips[j]]!;
      for (TreeNode? n = tips[j]; n != null; n = n.parent) {
        if (pi.containsKey(n)) {
          final d = pi[n]! + pj[n]!;
          if (d > best) {
            best = d;
            a = tips[i];
            b = tips[j];
            lca = n;
          }
          break;
        }
      }
    }
  }
  if (a == null || b == null || lca == null) return root;
  // The midpoint lies between the farther tip and the common ancestor.
  final far = paths[a]![lca]! >= paths[b]![lca]! ? a : b;
  var remaining = best / 2;
  TreeNode n = far;
  while (n.parent != null && !identical(n, lca) && remaining > (n.length ?? 0)) {
    remaining -= n.length ?? 0;
    n = n.parent!;
  }
  final r = rerootAbove(root, n);
  // Split the branch at the midpoint instead of in half.
  if (r.children.length == 2) {
    final total = (r.children[0].length ?? 0) + (r.children[1].length ?? 0);
    final side = r.children.firstWhere((c) => c.tips.contains(far), orElse: () => r.children.first);
    final other = r.children.firstWhere((c) => !identical(c, side));
    side.length = remaining.clamp(0, total).toDouble();
    other.length = total - side.length!;
  }
  return r;
}

/// Orders children by clade size (smaller first when [up], else larger).
void ladderize(TreeNode root, {bool up = true}) {
  for (final n in root.nodes) {
    n.children.sort((x, y) => up ? x.tipCount.compareTo(y.tipCount) : y.tipCount.compareTo(x.tipCount));
  }
}

/// Removes [node] (a tip or a whole clade) from the tree. Its parent, left
/// with one child, is replaced by that child (lengths joined); when the
/// root is left with one child, that child becomes the root. Returns the
/// root, or [root] unchanged when [node] is the root or the tree would be
/// left with fewer than two tips.
TreeNode pruneNode(TreeNode root, TreeNode node) {
  final parent = node.parent;
  if (parent == null || root.tipCount - node.tipCount < 2) return root;
  parent.children.remove(node);
  node.parent = null;
  if (parent.children.length != 1) return root;
  final only = parent.children.single;
  final grand = parent.parent;
  if (grand == null) {
    // The root kept one child: that child is the new root.
    only
      ..parent = null
      ..length = null;
    return only;
  }
  only.length = (only.length ?? 0) + (parent.length ?? 0);
  only.label ??= parent.label;
  grand.children[grand.children.indexOf(parent)] = only;
  only.parent = grand;
  return root;
}

/// Deep copy.
TreeNode copyTree(TreeNode n) {
  final c = TreeNode(name: n.name, length: n.length, label: n.label);
  for (final k in n.children) {
    c.add(copyTree(k));
  }
  return c;
}

/// Taxonomy and identifiers of one tip (from the GenBank records, or
/// guessed from the tip name for FASTA input).
class TipInfo {
  const TipInfo({
    this.species = '',
    this.genus = '',
    this.family = '',
    this.order = '',
    this.klass = '',
    this.phylum = '',
    this.kingdom = '',
    this.accession = '',
    this.voucher = '',
    this.locality = '',
    this.latLon = '',
    this.status = '',
    this.statusName = '',
    this.commonName = '',
    this.endemic = '',
    this.biomes = '',
    this.trend = '',
    this.officialList = '',
    this.threats = 0,
  });

  /// Binomial ("Genus species").
  final String species;
  final String genus, family, order, klass, phylum, kingdom, accession, voucher;

  /// Where the specimen was collected ("Colombia: Meta, San Juan de Arama")
  /// and its coordinates, when the record says (optional).
  final String locality, latLon;

  /// The country of [locality].
  String get country => locality.split(':').first.trim();

  /// From the species' SALVE sheet (Brazilian fauna): IUCN category code
  /// and name, first common name, endemic to Brazil ("Sim"/"Não") and
  /// biomes (comma separated).
  final String status, statusName, commonName, endemic, biomes;

  /// Population trend, in an official list ("Sim"/"Não"), number of threats.
  final String trend, officialList;
  final int threats;

  /// What the tip name tells: the first "Genus_species" pair, a family
  /// ending in -idae/-aceae and a GenBank accession.
  factory TipInfo.guess(String name) {
    final tokens = name.split(RegExp(r'[_\s|]+')).where((t) => t.isNotEmpty).toList();
    var genus = '', species = '';
    for (var i = 0; i + 1 < tokens.length; i++) {
      if (RegExp(r'^[A-Z][a-z]+$').hasMatch(tokens[i]) && RegExp(r'^[a-z][a-z-]+$').hasMatch(tokens[i + 1])) {
        genus = tokens[i];
        species = '${tokens[i]} ${tokens[i + 1]}';
        break;
      }
    }
    final family = tokens.firstWhere((t) => RegExp(r'^[A-Z][a-z]+(idae|aceae)$').hasMatch(t), orElse: () => '');
    final accession = RegExp(r'(?<![A-Za-z0-9])([A-Z]{2}_\d{6,9}|[A-Z]{1,2}\d{5,8})(\.\d+)?').firstMatch(name)?.group(0) ?? '';
    return TipInfo(species: species, genus: genus, family: family, accession: accession);
  }
}
