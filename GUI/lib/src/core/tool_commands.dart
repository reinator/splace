/// Command-line builders and output parsers for the external tools.
/// Pure Dart: paths are plain strings, nothing is executed here.
library;

/// [exe] and [args] as a shell command with one option per line, each line
/// ended by ` \` (a flag and its values stay together). Arguments with
/// spaces or shell characters are quoted.
String multilineCommand(String exe, List<String> args) {
  String quote(String a) => a.isEmpty || RegExp(r'''[\s'"$`*?;&|<>()]''').hasMatch(a) ? "'${a.replaceAll("'", "'\\''")}'" : a;
  bool isFlag(String a) => a.startsWith('-') && a.length > 1 && double.tryParse(a) == null;
  final lines = <String>[exe];
  for (final a in args) {
    if (isFlag(a) || lines.length == 1) {
      lines.add(quote(a));
    } else {
      lines[lines.length - 1] = '${lines.last} ${quote(a)}';
    }
  }
  return [for (final (i, l) in lines.indexed) '${i == 0 ? '' : '  '}$l${i < lines.length - 1 ? ' \\' : ''}'].join('\n');
}

class MafftOptions {
  const MafftOptions({
    this.strategy = '--auto',
    this.threads = 1,
    this.adjustDirection = false,
    this.maxIterate = 0,
    this.gapOpen = 1.53,
    this.gapExtend = 0.0,
    this.reorder = false,
    this.preserveCase = false,
    this.extra = '',
    this.amino = false,
  });

  /// e.g. `--auto`, `--localpair` (L-INS-i), `--retree 2 --reorder`.
  final String strategy;
  final int threads;

  /// `--adjustdirection`: reverse-complement sequences given in the wrong
  /// orientation (common in data assembled from several sources).
  final bool adjustDirection;

  /// `--maxiterate` (0 = no refinement).
  final int maxIterate;

  /// `--op` / `--ep`; written only when they differ from MAFFT's defaults.
  final double gapOpen;
  final double gapExtend;
  final bool reorder;
  final bool preserveCase;
  final String extra;

  /// `--amino` (protein input: amino-acid and codon-aware modes).
  final bool amino;

  MafftOptions copyWith({
    String? strategy,
    int? threads,
    bool? adjustDirection,
    int? maxIterate,
    double? gapOpen,
    double? gapExtend,
    bool? reorder,
    bool? preserveCase,
    String? extra,
    bool? amino,
  }) =>
      MafftOptions(
        strategy: strategy ?? this.strategy,
        threads: threads ?? this.threads,
        adjustDirection: adjustDirection ?? this.adjustDirection,
        maxIterate: maxIterate ?? this.maxIterate,
        gapOpen: gapOpen ?? this.gapOpen,
        gapExtend: gapExtend ?? this.gapExtend,
        reorder: reorder ?? this.reorder,
        preserveCase: preserveCase ?? this.preserveCase,
        extra: extra ?? this.extra,
        amino: amino ?? this.amino,
      );

  List<String> args(String input) {
    final strat = strategy.split(RegExp(r'\s+')).where((a) => a.isNotEmpty).toList();
    return [
      ...strat,
      if (maxIterate > 0 && !strat.contains('--maxiterate')) ...['--maxiterate', '$maxIterate'],
      if (gapOpen != 1.53) ...['--op', '$gapOpen'],
      if (gapExtend != 0.0) ...['--ep', '$gapExtend'],
      if (reorder && !strat.contains('--reorder')) '--reorder',
      if (preserveCase) '--preservecase',
      if (amino) '--amino',
      '--thread', '$threads',
      if (adjustDirection && !amino) '--adjustdirection',
      ...extra.split(RegExp(r'\s+')).where((a) => a.isNotEmpty),
      // No --quiet: SPLACE reads the progress counters MAFFT prints.
      input,
    ];
  }
}

class TrimalOptions {
  const TrimalOptions({this.method = '-automated1', this.gapThreshold = 0.9, this.similarity, this.consistency, this.window, this.extra = ''});

  /// e.g. `-automated1`, `-gappyout`, `-strict`, or `manual` (thresholds).
  final String method;
  final double gapThreshold;
  final double? similarity;
  final int? consistency;
  final int? window;
  final String extra;

  bool get manual => method == 'manual';

  /// Trimming flags without input/output.
  List<String> get flags => [
        if (manual) ...[
          '-gt', '$gapThreshold',
          if (similarity != null) ...['-st', '$similarity'],
          if (consistency != null) ...['-cons', '$consistency'],
          if (window != null) ...['-w', '$window'],
        ] else
          ...method.split(RegExp(r'\s+')).where((a) => a.isNotEmpty),
        ...extra.split(RegExp(r'\s+')).where((a) => a.isNotEmpty),
      ];

  String get label => manual ? flags.join(' ') : method;

  List<String> args(String input, String output, {String? backtranslateFrom}) => [
        '-in', input,
        '-out', output,
        ...flags,
        if (backtranslateFrom != null) ...['-backtrans', backtranslateFrom],
      ];
}

class IqtreeOptions {
  const IqtreeOptions({
    this.model = 'MFP+MERGE',
    this.ufboot = 1000,
    this.alrt = 1000,
    this.classicBootstrap = 0,
    this.threads = 'AUTO',
    this.maxThreads,
    this.seed,
    this.outgroup = const [],
    this.raxmlModelsOnly = false,
    this.bnni = false,
    this.nstop,
    this.perturbation,
    this.runs = 1,
    this.extra = '',
  });

  /// `MFP` (ModelFinder), `MFP+MERGE` (merge partitions), `MF+MERGE`
  /// (only choose the scheme, no tree), or a fixed model like `GTR+F+G4`.
  final String model;

  /// Ultrafast bootstrap replicates (`-B`), 0 disables.
  final int ufboot;

  /// SH-aLRT replicates (`--alrt`), 0 disables.
  final int alrt;

  /// Standard nonparametric bootstrap (`-b`), 0 disables.
  final int classicBootstrap;
  final String threads;
  final int? maxThreads;
  final int? seed;
  final List<String> outgroup;

  /// `--mset raxml`: only consider models RAxML-NG can also use, so the
  /// chosen scheme can be passed on to RAxML-NG.
  final bool raxmlModelsOnly;

  /// `--bnni`: optimises each UFBoot tree by NNI, which reduces the risk of
  /// overestimated support under model violations.
  final bool bnni;

  /// `-nstop`: unsuccessful iterations before the search stops (IQ-TREE
  /// default 100); null keeps the default.
  final int? nstop;

  /// `-pers`: perturbation strength of the search (IQ-TREE default 0.5).
  final double? perturbation;

  /// `--runs`: independent tree searches; the best tree is kept.
  final int runs;
  final String extra;

  IqtreeOptions copyWith({String? threads, int? maxThreads}) => IqtreeOptions(
        model: model,
        ufboot: ufboot,
        alrt: alrt,
        classicBootstrap: classicBootstrap,
        threads: threads ?? this.threads,
        maxThreads: maxThreads ?? this.maxThreads,
        seed: seed,
        outgroup: outgroup,
        raxmlModelsOnly: raxmlModelsOnly,
        bnni: bnni,
        nstop: nstop,
        perturbation: perturbation,
        runs: runs,
        extra: extra,
      );

  bool get treeInference => !model.startsWith('MF+') && model != 'MF';

  /// Partitions above which MERGE uses the relaxed clustering algorithm
  /// (`--rcluster 10`): trying every pair of hundreds of partitions would
  /// take days, and the IQ-TREE manual recommends it for large data sets.
  static const rclusterAbove = 50;

  /// [partitionCount] partitions in the partition file (0 when none).
  List<String> args({required String alignment, String? partitions, required String prefix, int partitionCount = 0}) => [
        '-s', alignment,
        if (partitions != null) ...['-p', partitions],
        '-m', model,
        if (partitions != null && model.contains('MERGE') && partitionCount > rclusterAbove && !extra.contains('rcluster')) ...['--rcluster', '10'],
        if (raxmlModelsOnly) ...['--mset', 'raxml'],
        if (treeInference && ufboot > 0) ...['-B', '$ufboot'],
        if (treeInference && classicBootstrap > 0) ...['-b', '$classicBootstrap'],
        if (treeInference && alrt > 0) ...['--alrt', '$alrt'],
        if (treeInference && bnni && ufboot > 0) '--bnni',
        if (treeInference && nstop != null) ...['-nstop', '$nstop'],
        if (treeInference && perturbation != null) ...['-pers', '$perturbation'],
        if (treeInference && runs > 1) ...['--runs', '$runs'],
        '-T', threads,
        if (maxThreads != null) ...['--threads-max', '$maxThreads'],
        if (seed != null) ...['--seed', '$seed'],
        if (outgroup.isNotEmpty) ...['-o', outgroup.join(',')],
        '--prefix', prefix,
        '--redo',
        ...extra.split(RegExp(r'\s+')).where((a) => a.isNotEmpty),
      ];
}

/// Where RAxML-NG gets its substitution models from.
enum RaxmlModels {
  /// MOOSE, the model selection built into RAxML-NG 2.0 (`--model DNA` or
  /// `--model auto=partitions.txt`).
  moose,

  /// IQ-TREE's ModelFinder restricted to RAxML-NG models (`--mset raxml`).
  modelFinder,

  /// One model for every partition.
  fixed,
}

/// Branch support measures of RAxML-NG (`--bs-metric`). [replicates] need
/// bootstrap trees; the others are computed on the ML tree only.
enum RaxmlSupport {
  fbp('fbp', 'FBP', 'supportFBP', true),
  tbe('tbe', 'TBE', 'supportTBE', true),
  rbs('rbs', 'RBS', 'supportRBS', true),
  sh('sh', 'SH-aLRT', 'supportSH', false),
  ebg('ebg', 'EBG', 'supportEBG', false);

  const RaxmlSupport(this.flag, this.label, this.suffix, this.replicates);
  final String flag;
  final String label;

  /// Output file: `<prefix>.raxml.<suffix>`.
  final String suffix;
  final bool replicates;

  /// Output key in [TreeRun.outputs] (`treeFBP`, `treeSH`…).
  String get outputKey => 'tree${suffix.substring(7)}';

  static List<RaxmlSupport> parse(String metrics) => [
        for (final m in metrics.split(',').map((x) => x.trim().toLowerCase()))
          ...RaxmlSupport.values.where((v) => v.flag == m),
      ];
}

class RaxmlOptions {
  const RaxmlOptions({
    this.bootstraps = 'autoMRE{1000}',
    this.startTrees = '',
    this.maxThreads = 4,
    this.seed = 12345,
    this.outgroup = const [],
    this.supportMetrics = 'fbp,tbe',
    this.fast = false,
    this.criterion = 'BIC',
    this.extra = '',
  });

  /// Number of bootstrap trees or `autoMRE{N}` (bootstopping).
  final String bootstraps;

  /// `--tree`: empty for RAxML-NG 2's adaptive choice of starting trees,
  /// or e.g. `pars{10},rand{10}`.
  final String startTrees;

  /// Used as `--threads auto{N}` so RAxML-NG picks an efficient number up
  /// to N based on the alignment size.
  final int maxThreads;
  final int seed;
  final List<String> outgroup;

  /// Comma-separated [RaxmlSupport] flags (`fbp,tbe`, `sh,ebg`…). FBP and
  /// RBS cannot be combined.
  final String supportMetrics;

  /// `--fast`: simplified search with early stopping (RAxML-NG 2).
  final bool fast;

  /// MOOSE information criterion: BIC, AIC or AICc.
  final String criterion;
  final String extra;

  List<RaxmlSupport> get metrics => RaxmlSupport.parse(supportMetrics);

  RaxmlOptions copyWith({int? maxThreads}) => RaxmlOptions(
        bootstraps: bootstraps,
        startTrees: startTrees,
        maxThreads: maxThreads ?? this.maxThreads,
        seed: seed,
        outgroup: outgroup,
        supportMetrics: supportMetrics,
        fast: fast,
        criterion: criterion,
        extra: extra,
      );

  /// True when a measure needs bootstrap replicate trees.
  bool get needsReplicates => metrics.any((m) => m.replicates);

  /// [dataType] (`DNA`/`AA`) is needed by model selection (`--model DNA`
  /// or `auto=`), which does not guess it from the alignment.
  List<String> parseArgs({required String msa, required String model, required String prefix, String? dataType}) => [
        '--parse', '--msa', msa, '--model', model, if (dataType != null) ...['--data-type', dataType], '--prefix', prefix, '--redo',
      ];

  /// Full analysis: ML search + support values. [msa] may be the
  /// `.raxml.rba` produced by `--parse`, in which case the model is read
  /// from it and [model] can be null. With [moose] the model is `DNA`, `AA`
  /// or `auto=<partition file>` and MOOSE selects the models first (the
  /// alignment must then be the PHYLIP file: a `.rba` keeps its models).
  List<String> allArgs({required String msa, String? model, required String prefix, bool moose = false, String dataType = 'DNA'}) => [
        '--all',
        if (fast) '--fast',
        '--msa', msa,
        if (model != null) ...['--model', model],
        if (moose) ...['--data-type', dataType, '--moose-options', 'criterion=${criterion.toLowerCase()}'],
        '--prefix', prefix,
        if (startTrees.trim().isNotEmpty) ...['--tree', startTrees.trim()],
        if (needsReplicates) ...['--bs-trees', bootstraps],
        '--bs-metric', supportMetrics,
        '--threads', 'auto{$maxThreads}',
        '--seed', '$seed',
        if (outgroup.isNotEmpty) ...['--outgroup', outgroup.join(',')],
        '--redo',
        ...extra.split(RegExp(r'\s+')).where((a) => a.isNotEmpty),
      ];
}

/// Models chosen by MOOSE, from the RAxML-NG output:
/// `Partition #0: HKY+FO+I+G4m (LogLH = -1782.1  BIC = 3672.7)`.
List<String> parseMooseModels(String log) => [
      for (final m in RegExp(r'Partition #\d+: (\S+) \(LogLH').allMatches(log)) m.group(1)!,
    ];

/// Facts read from a RAxML-NG log.
class RaxmlLogSummary {
  RaxmlLogSummary({this.recommendedThreads, this.finalLogLikelihood, this.bootstrapTrees, this.memory});
  final int? recommendedThreads;
  final double? finalLogLikelihood;
  final int? bootstrapTrees;
  final String? memory;
}

RaxmlLogSummary parseRaxmlLog(String log) {
  int? intOf(RegExp re) => int.tryParse(re.firstMatch(log)?.group(1) ?? '');
  return RaxmlLogSummary(
    recommendedThreads: intOf(RegExp(r'Recommended number of threads[^:]*:\s*(\d+)')),
    finalLogLikelihood: double.tryParse(RegExp(r'Final LogLikelihood:\s*(-?[\d.]+)').firstMatch(log)?.group(1) ?? ''),
    bootstrapTrees: intOf(RegExp(r'Bootstrapping converged after (\d+) replicates')),
    memory: RegExp(r'Estimated memory requirements\s*:\s*(.+)').firstMatch(log)?.group(1)?.trim(),
  );
}

/// IQ-TREE base model names → RAxML-NG names. IQ-TREE uses suffix `e` for
/// equal frequencies (TIM2e) and `u` for unequal (TPM2u); RAxML-NG uses
/// `ef` / `uf` the other way round for some families.
const _iqToRaxmlBase = {
  'JC': 'JC', 'JC69': 'JC', 'F81': 'F81',
  'K2P': 'K80', 'K80': 'K80', 'HKY': 'HKY', 'HKY85': 'HKY',
  'TNE': 'TN93ef', 'TN': 'TN93', 'TRN': 'TN93', 'TN93': 'TN93',
  'K3P': 'K81', 'K81': 'K81', 'K3PU': 'K81uf', 'K81U': 'K81uf',
  'TPM2': 'TPM2', 'TPM2U': 'TPM2uf', 'TPM3': 'TPM3', 'TPM3U': 'TPM3uf',
  'TIME': 'TIM1', 'TIM': 'TIM1uf', 'TIM2E': 'TIM2', 'TIM2': 'TIM2uf',
  'TIM3E': 'TIM3', 'TIM3': 'TIM3uf', 'TVME': 'TVMef', 'TVM': 'TVM',
  'SYM': 'SYM', 'GTR': 'GTR',
};

/// Converts an IQ-TREE DNA model string (e.g. `TIM2+F+I+G4`) to RAxML-NG
/// syntax. Returns null when the base model has no RAxML-NG equivalent.
String? iqtreeModelToRaxml(String model) {
  final parts = model.trim().split('+');
  final base = _iqToRaxmlBase[parts.first.toUpperCase()];
  if (base == null) return null;
  final out = [base];
  for (final m in parts.skip(1)) {
    final u = m.toUpperCase();
    if (u == 'F' || u == 'FO' || u == 'FQ') {
      out.add(u == 'FQ' ? 'FE' : u); // RAxML-NG: +FE equal, +FO ML, +F/+FC empirical
    } else if (u == 'I' || u == 'IO') {
      out.add('I');
    } else if (RegExp(r'^G\d*$').hasMatch(u)) {
      out.add(u == 'G' ? 'G4' : u);
    } else if (RegExp(r'^R\d+$').hasMatch(u)) {
      out.add(u);
    } else if (u == 'ASC') {
      // Dropped: RAxML-NG's ASC_LEWIS aborts on partitions with invariant
      // sites, which IQ-TREE's +ASC choice does not guarantee.
      continue;
    } else {
      return null;
    }
  }
  return out.join('+');
}

class SchemeEntry {
  SchemeEntry(this.model, this.name, this.ranges);
  final String model;
  final String name;
  final String ranges;
}

/// Parses an IQ-TREE `.best_scheme` (RAxML-style) partition file:
/// `TIM2+F+G4, COI_pos1 = 1-600\3, 601-900\3`.
List<SchemeEntry> parseRaxmlStyleScheme(String text) {
  final entries = <SchemeEntry>[];
  final line = RegExp(r'^\s*([^,\s]+)\s*,\s*([^=\s]+)\s*=\s*(.+?)\s*;?\s*$');
  for (final l in text.split('\n')) {
    final m = line.firstMatch(l);
    if (m != null) entries.add(SchemeEntry(m.group(1)!, m.group(2)!, m.group(3)!));
  }
  return entries;
}

/// Parses the `charset`/`charpartition` blocks of IQ-TREE's
/// `.best_scheme.nex` or `.best_model.nex`. IQ-TREE 3 writes only `DNA` in
/// the RAxML-style `.best_scheme`, so the chosen models must be read here.
/// Parameter values in braces (`GTR{0.8,...}+F{...}`) are stripped.
List<SchemeEntry> parseIqtreeNexusScheme(String text) {
  final charsets = <String, String>{};
  for (final m in RegExp(r'charset\s+(\S+)\s*=\s*([^;]+);', caseSensitive: false).allMatches(text)) {
    charsets[m.group(1)!] = m.group(2)!.trim().split(RegExp(r'\s+')).join(', ');
  }
  final cp = RegExp(r'charpartition\s+\S+\s*=\s*([^;]+);', caseSensitive: false).firstMatch(text);
  if (cp == null) return const [];
  final entries = <SchemeEntry>[];
  for (final item in cp.group(1)!.replaceAll(RegExp(r'\{[^}]*\}'), '').split(',')) {
    final colon = item.lastIndexOf(':');
    if (colon < 0) continue;
    final model = item.substring(0, colon).trim();
    final name = item.substring(colon + 1).trim();
    final ranges = charsets[name];
    if (ranges != null) entries.add(SchemeEntry(model, name, ranges));
  }
  return entries;
}

/// Turns an IQ-TREE scheme into a RAxML-NG partition file. Models without a
/// RAxML-NG equivalent fall back to [fallback]; every change is reported in
/// `replaced`.
({String partitions, List<String> replaced}) schemeToRaxml(List<SchemeEntry> scheme, {String fallback = 'GTR+G4'}) {
  final b = StringBuffer();
  final replaced = <String>[];
  for (final e in scheme) {
    var model = iqtreeModelToRaxml(e.model);
    if (model == null) {
      replaced.add('${e.name}: ${e.model} -> $fallback');
      model = fallback;
    } else if (RegExp(r'\+ASC', caseSensitive: false).hasMatch(e.model)) {
      replaced.add('${e.name}: ${e.model} -> $model (ascertainment correction dropped)');
    }
    b.writeln('$model, ${e.name} = ${e.ranges}');
  }
  return (partitions: b.toString(), replaced: replaced);
}
