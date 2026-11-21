/// Species delimitation, by published programs:
/// - on the inferred trees, GMYC (Pons et al. 2006; pGMYC, Zhang et al.
///   2013) on the ML tree made ultrametric by LSD2 (To et al. 2016, relative
///   dates: root 0, tips 1);
/// - on the inferred trees, mPTP (Kapli et al. 2017) with the single-rate
///   PTP model (Zhang et al. 2013, `--single`) and the multi-rate mPTP model
///   (`--multi`), by maximum likelihood. IQ-TREE and RAxML-NG do not delimit
///   species themselves; their ML trees are what PTP and mPTP read;
/// - on the alignments (pairwise distances), ASAP (Puillandre, Brouillet &
///   Achaz 2021) and ABGD (Puillandre et al. 2012), which write SPART files
///   (Miralles et al. 2022).
///
/// Pure Dart: the arguments of each call, the files written for the
/// programs and the parsing of their output.
library;

import 'dart:convert';
import 'dart:math' as math;

import 'newick.dart';

/// Model of the Poisson tree processes.
enum PtpModel {
  /// One coalescent rate for every species (PTP).
  single,

  /// One coalescent rate per species (mPTP).
  multi;

  String get label => this == single ? 'PTP' : 'mPTP';
  String get flag => this == single ? '--single' : '--multi';
}

/// Arguments of `mptp --minbr_auto`: the smallest branch length that
/// still means a substitution, from the p-distances of [fasta]. The
/// answer is written to `<prefix>.txt`. Never with `--outgroup_crop`: the
/// cropped tree has fewer tips than [fasta] and mPTP stops ("FASTA file
/// contains more sequences than expected"); the distances do not depend
/// on the root.
List<String> mptpMinbrArgs({required String tree, required String fasta, required String prefix, List<String> outgroup = const []}) => [
      '--minbr_auto', fasta,
      '--tree_file', tree,
      '--output_file', prefix,
      if (outgroup.isNotEmpty) ...['--outgroup', outgroup.join(',')],
    ];

/// mPTP's significance level of the likelihood-ratio test (`--pvalue`).
const mptpDefaultPvalue = 0.001;

/// GMYC's significance level (as in the splits R package, 0.01 here too).
const gmycDefaultPvalue = 0.01;

/// Arguments of an ML delimitation: the tree is rooted on [outgroup] (left
/// out of the delimitation with [cropOutgroup]); branches up to [minbr]
/// are ignored. Writes `<prefix>.txt` and `<prefix>.svg`.
List<String> mptpArgs({
  required String tree,
  required String prefix,
  required PtpModel model,
  List<String> outgroup = const [],
  bool cropOutgroup = false,
  double? minbr,
  double? pvalue,
  int? seed,
}) =>
    [
      '--ml',
      model.flag,
      '--tree_file', tree,
      '--output_file', prefix,
      if (outgroup.isNotEmpty) ...['--outgroup', outgroup.join(',')],
      if (outgroup.isNotEmpty && cropOutgroup) '--outgroup_crop',
      if (minbr != null) ...['--minbr', _num(minbr)],
      if (pvalue != null) ...['--pvalue', _num(pvalue)],
      if (seed != null) ...['--seed', '$seed'],
    ];

String _num(double v) => v == 0 ? '0' : v.toStringAsPrecision(6).replaceFirst(RegExp(r'\.?0+$'), '');

/// `Minimum branch length (--minbr) should be set to X` of
/// `--minbr_auto`, or null.
double? parseMinbr(String text) {
  final m = RegExp(r'should be set to\s+([0-9.eE+-]+)').firstMatch(text);
  return m == null ? null : double.tryParse(m.group(1)!);
}

/// The tree as mPTP reads it: binary (polytomies resolved with
/// zero-length branches), without support labels, every branch with a
/// length. mPTP refuses trees with other node degrees, and labels such as
/// IQ-TREE's `80.3/95` are of no use to it.
String mptpNewick(TreeNode root) {
  String len(TreeNode n) => ':${_num((n.length ?? 0).isNaN ? 0 : (n.length ?? 0).clamp(0, double.infinity).toDouble())}';
  String quote(String name) => RegExp(r"[\s(),:;\[\]'\x22]").hasMatch(name) ? "'${name.replaceAll("'", '_')}'" : name;

  // Children joined two at a time: (a,(b,(c,d))) with zero-length inner
  // branches for a node with more than [keep] children.
  String join(List<String> parts, int keep) {
    if (parts.length <= keep) return parts.join(',');
    final rest = parts.sublist(keep - 1);
    var inner = rest.last;
    for (var i = rest.length - 2; i >= 0; i--) {
      inner = '(${rest[i]},$inner):0';
    }
    return [...parts.sublist(0, keep - 1), inner].join(',');
  }

  String w(TreeNode n) {
    if (n.isTip) return '${quote(n.name)}${len(n)}';
    // A node with one child is just its child, with both branches.
    if (n.children.length == 1) {
      final c = n.children.single;
      final merged = TreeNode(name: c.name, length: (c.length ?? 0) + (n.length ?? 0));
      for (final k in [...c.children]) {
        merged.add(k);
      }
      return w(merged);
    }
    return '(${join([for (final c in n.children) w(c)], 2)})${len(n)}';
  }

  // The root keeps two or three children (rooted or unrooted).
  final kids = [for (final c in root.children) w(c)];
  if (kids.length == 1) return '${kids.single.replaceFirst(RegExp(r':[^:()]*$'), '')};';
  return '(${join(kids, 3)});';
}

/// The delimitation methods, by program.
enum DelimMethod {
  ptp,
  mptp,
  gmyc,
  asap,
  abgd;

  String get label => switch (this) { ptp => 'PTP', mptp => 'mPTP', gmyc => 'GMYC', asap => 'ASAP', abgd => 'ABGD' };

  /// Reads a tree (PTP, mPTP, GMYC) rather than an alignment (ASAP, ABGD).
  bool get onTree => this == ptp || this == mptp || this == gmyc;

  PtpModel? get ptpModel => switch (this) { ptp => PtpModel.single, mptp => PtpModel.multi, _ => null };

  static DelimMethod ofModel(PtpModel m) => m == PtpModel.single ? ptp : mptp;
}

/// Pairwise distance of ASAP and ABGD (`-d`; both default to Jukes-Cantor).
enum DelimDistance {
  jc69(1),
  k2p(0),
  tn93(2),
  p(3);

  const DelimDistance(this.code);
  final int code;

  String get label => switch (this) { jc69 => 'Jukes-Cantor (JC69)', k2p => 'Kimura 2-parameter (K2P)', tn93 => 'Tamura-Nei (TN93)', p => 'p-distance' };

  /// ABGD offers only JC69 and K2P.
  bool get abgd => this == jc69 || this == k2p;
}

/// Arguments of ASAP on [fasta] (a file name in the working folder),
/// writing to the folder [outDir] (relative too: ASAP and ABGD take the
/// last `/` of a path as the start of the file name).
List<String> asapArgs({required String fasta, required String outDir, DelimDistance distance = DelimDistance.jc69}) =>
    ['-d', '${distance.code}', '-o', outDir, fasta];

/// Arguments of ABGD on [fasta]: priors of intraspecific divergence from
/// [pMin] to [pMax] in [steps] steps (log scale), relative gap width
/// [slope] (`-X`); `-a` writes the partitions as SPART files.
List<String> abgdArgs({
  required String fasta,
  required String outDir,
  DelimDistance distance = DelimDistance.jc69,
  double pMin = 0.001,
  double pMax = 0.1,
  int steps = 10,
  double slope = 1.5,
}) =>
    ['-a', '-d', '${distance.code}', '-p', _num(pMin), '-P', _num(pMax), '-n', '$steps', '-X', _num(slope), '-o', outDir, fasta];

/// The priors ABGD tries (Compute_myDist in its main_abgd.c): [steps]
/// values from [pMin] to [pMax] on a log scale.
List<double> abgdPriors(double pMin, double pMax, int steps) {
  if (steps <= 1) return [pMin];
  final scale = (math.log(pMax / pMin) / math.ln10) / (steps - 1);
  final init = math.log(pMin) / math.ln10;
  return [for (var i = 0; i < steps - 1; i++) i == 0 ? pMin : math.pow(10, init + scale * i).toDouble(), pMax];
}

/// The sequences of [records] as ASAP and ABGD read them: short ids
/// (`s0001`…, so long or unusual names are neither cut at 100 characters
/// nor rewritten), upper case, unknown bases (`?`) as `N`. Returns the
/// FASTA text and the id → name map.
(String, Map<String, String>) delimitationFasta(List<(String, String)> records) {
  final names = <String, String>{};
  final b = StringBuffer();
  for (final (i, (name, seq)) in records.indexed) {
    final id = 's${(i + 1).toString().padLeft(4, '0')}';
    names[id] = name;
    b
      ..writeln('>$id')
      ..writeln(seq.toUpperCase().replaceAll(RegExp(r'[?.]'), 'N').replaceAll(RegExp(r'\s'), ''));
  }
  return (b.toString(), names);
}

/// One partition of a SPART file: its name, the number of subsets
/// (species), its score as written ([score]: ASAP score, ABGD barcode gap
/// distance) and the members of each subset.
class Spartition {
  const Spartition(this.name, this.species, [this.score]);
  final String name;
  final List<List<String>> species;
  final String? score;
}

/// Reads the partitions of a SPART file (ASAP, ABGD), in their order;
/// [names] turns the ids of the file back into the sample names.
List<Spartition> parseSpart(String text, [Map<String, String> names = const {}]) {
  final lines = text.split(RegExp(r'\r?\n'));
  final spartNames = RegExp(r'N_spartitions\s*=\s*\d+\s*:\s*([^;]*);').firstMatch(text)?.group(1)?.split('/').map((s) => s.trim()).toList() ?? <String>[];
  // Scores: the bracketed line after "Asap scores" or "Barcode gap
  // distance".
  List<String> scores = const [];
  for (var i = 0; i < lines.length - 1; i++) {
    if (RegExp(r'Asap scores|Barcode gap distance', caseSensitive: false).hasMatch(lines[i])) {
      scores = lines[i + 1].replaceAll(RegExp(r'[\[\]]'), '').split('/').map((s) => s.trim()).toList();
      break;
    }
  }
  // Individual_assignment: "name : 1 / 2 / 1", the last line ends with ";".
  final assigned = <(String, List<int>)>[];
  final start = lines.indexWhere((l) => l.trim().startsWith('Individual_assignment'));
  if (start >= 0) {
    for (final raw in lines.skip(start + 1)) {
      final line = raw.trim();
      if (line.isEmpty || line.startsWith('[')) continue;
      if (line.toLowerCase().startsWith('end')) break;
      final m = RegExp(r'^(\S+)\s*:\s*([0-9/\s]+);?$').firstMatch(line);
      if (m != null) {
        final ids = [for (final v in m.group(2)!.split('/')) int.tryParse(v.trim()) ?? 0];
        assigned.add((names[m.group(1)!] ?? m.group(1)!, ids));
      }
      if (line.endsWith(';')) break;
    }
  }
  final count = assigned.isEmpty ? 0 : assigned.map((a) => a.$2.length).reduce(math.max);
  return [
    for (var k = 0; k < count; k++)
      () {
        final groups = <int, List<String>>{};
        for (final (name, ids) in assigned) {
          if (k < ids.length) groups.putIfAbsent(ids[k], () => []).add(name);
        }
        final keys = groups.keys.toList()..sort();
        return Spartition(k < spartNames.length ? spartNames[k] : '${k + 1}', [for (final g in keys) groups[g]!], k < scores.length ? scores[k] : null);
      }(),
  ];
}

/// Arguments of LSD2 making [tree] ultrametric for GMYC: relative dates
/// (root 0, every tip 1, `-a 0 -z 1`), [seqLength] sites; rooted on the
/// outgroups listed in [outgroupFile] and without them (`-g … -G`), or with
/// the root estimated by LSD2 (`-r a`). Writes `<prefix>.date.nexus`.
List<String> lsd2Args({required String tree, required String prefix, required int seqLength, String? outgroupFile}) => [
      '-i', tree,
      '-s', '$seqLength',
      '-a', '0',
      '-z', '1',
      if (outgroupFile != null) ...['-g', outgroupFile, '-G'] else ...['-r', 'a'],
      '-o', prefix,
    ];

/// LSD2's outgroup file: the number of outgroups, then one per line.
String lsd2OutgroupFile(List<String> outgroup) => '${outgroup.length}\n${outgroup.join('\n')}\n';

/// The dated tree of LSD2's `.date.nexus` (branch lengths in time) as a
/// plain Newick tree for GMYC: the first `tree … =` statement, without the
/// `[&date=…]` annotations and internal labels. Null when there is none.
String? ultrametricFromNexus(String nexus) {
  final m = RegExp(r'tree\s+[^=]*=\s*(.*?;)', caseSensitive: false, dotAll: true).firstMatch(nexus);
  if (m == null) return null;
  final text = m.group(1)!.replaceAll(RegExp(r'\[[^\]]*\]'), '').replaceAll(RegExp(r'\s+'), '');
  final tree = parseNewick(text);
  for (final n in tree.nodes) {
    if (!n.isTip) {
      n
        ..label = null
        ..name = '';
    }
  }
  return toNewick(tree);
}

/// Arguments of the GMYC command line (tool/gmyc/splace_gmyc.py).
List<String> gmycArgs({required String tree, required String out, double pvalue = 0.01}) => [tree, '--pvalue', _num(pvalue), '--out', out];

/// Reads the JSON the GMYC command line writes.
Delimitation parseGmyc(String json, {double? alpha}) {
  final d = jsonDecode(json) as Map<String, dynamic>;
  double? n(String k) => (d[k] as num?)?.toDouble();
  return Delimitation(
    species: [for (final s in d['species'] as List) [for (final t in s as List) '$t']],
    nullLogL: n('null_logl'),
    logL: n('best_logl'),
    pvalue: n('pvalue'),
    lrtPassed: d['lrt_passed'] as bool?,
    bestSpecies: (d['threshold_species'] as num?)?.toInt(),
    alpha: alpha,
  );
}

/// Result of one delimitation run.
class Delimitation {
  Delimitation({
    this.model,
    required this.species,
    this.nullLogL,
    this.logL,
    this.pvalue,
    this.lrtPassed,
    this.minbr,
    this.bestSpecies,
    this.alpha,
  });

  /// PTP or mPTP (null for ASAP and ABGD).
  final PtpModel? model;

  /// Tips of each delimited species, in the order the program lists them.
  final List<List<String>> species;
  final double? nullLogL, logL, pvalue;

  /// Likelihood-ratio test against the null model (all tips one species).
  final bool? lrtPassed;

  /// Species of the best model (mPTP "Number of delimited species", the
  /// GMYC threshold): when the test fails the program still gives one
  /// species for all tips, this is what it had found.
  final int? bestSpecies;

  /// Significance level of the test (`--pvalue`).
  final double? alpha;

  /// The test failed: every tip is one species (the null model).
  bool get collapsed => lrtPassed == false;

  /// Minimum branch length used (`--minbr`).
  double? minbr;

  /// Species with one tip each.
  int get singletons => species.where((s) => s.length == 1).length;

  /// Species number of each tip.
  Map<String, int> get speciesOf => {
        for (final (i, s) in species.indexed)
          for (final tip in s) tip: i,
      };
}

/// Reads `<prefix>.txt` of an ML run.
Delimitation parseMptp(String text, PtpModel model, {double? alpha}) {
  double? value(String key) {
    final m = RegExp('${RegExp.escape(key)}[^:]*:\\s*(-?[0-9.eE+-]+|nan|-?inf)', caseSensitive: false).firstMatch(text);
    return m == null ? null : double.tryParse(m.group(1)!);
  }

  final lrt = RegExp(r'^LRT:\s*(\w+)', multiLine: true).firstMatch(text)?.group(1)?.toLowerCase();
  final species = <List<String>>[];
  List<String>? current;
  for (final raw in text.split(RegExp(r'\r?\n'))) {
    final line = raw.trim();
    if (RegExp(r'^Species\s+\d+:$').hasMatch(line)) {
      current = <String>[];
      species.add(current);
    } else if (current != null && line.isNotEmpty) {
      current.add(line);
    } else if (line.isEmpty) {
      current = null;
    }
  }
  return Delimitation(
    model: model,
    species: species.where((s) => s.isNotEmpty).toList(),
    nullLogL: value('Null-model score'),
    logL: value('Best score for'),
    pvalue: value('LRT computed p-value'),
    lrtPassed: lrt == null ? null : lrt == 'passed',
    bestSpecies: value('Number of delimited species')?.round(),
    alpha: alpha,
  );
}

/// Delimitation as a table: one line per tip with its species.
String delimitationTsv(Delimitation d) => [
      'species\ttip',
      for (final (i, s) in d.species.indexed)
        for (final tip in s) '${i + 1}\t$tip',
    ].join('\n');
