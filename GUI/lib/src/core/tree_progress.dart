import 'dart:math' as math;

/// Stages of a tree inference shown in the progress modal.
enum TreeStage { prepare, modelFinder, search, support, finish }

/// Reads IQ-TREE and RAxML-NG output line by line and tracks the stage and
/// the numbers worth showing (models tested, iterations, bootstrap
/// replicates…). Lines of both programs can be fed to one instance (the
/// RAxML-NG run is preceded by IQ-TREE's ModelFinder).
///
/// IQ-TREE prints, in order:
///
///     ModelFinder will test up to 968 DNA models (…)
///       139  HKY+R4        18949.087    96  38090.174 …
///     Best-fit model: TIM2+F+R4 chosen according to BIC
///     Generating 1000 samples for ultrafast bootstrap (seed: …)
///     Iteration 60 / LogL: -18692.587 / Time: 0h:0m:32s (0h:0m:8s left)
///     BEST SCORE FOUND : -18692.166
///     Creating bootstrap support values...
///
/// RAxML-NG prints:
///
///     Starting model selection with parsimony starting tree (MOOSE)
///     Partition #0: HKY+FO+I+G4m (LogLH = -1782.169337  BIC = 3672.710540)
///     [00:00:02] Starting ML tree search with 20 distinct starting trees
///     [00:00:05] [worker #0] ML tree search #1, logLikelihood: -18700.1
///     [00:00:30] Starting bootstrapping analysis with 1000 replicates.
///     [00:00:33] [worker #1] Bootstrap tree #12, logLikelihood: …
///     [00:05:00] Bootstrapping converged after 300 replicates.
class TreeProgress {
  TreeProgress({required this.stages, this.raxml = false});

  /// Stages of this run, in order (e.g. without [TreeStage.modelFinder]
  /// when the model is fixed).
  final List<TreeStage> stages;

  /// RAxML-NG search: the support stage is its bootstrap.
  final bool raxml;

  TreeStage stage = TreeStage.prepare;
  int? modelsDone, modelsTotal;

  /// RAxML-NG: partitions of the alignment (MOOSE picks one model each).
  int? partitions;

  /// Models selected by RAxML-NG's MOOSE (counts are partitions).
  bool moose = false;
  String? bestModel;
  int? iteration;
  double? logL, bestLogL;
  int? secondsLeft, secondsElapsed;
  int? startsDone, startsTotal;
  int? bootDone, bootTotal;
  bool bootConverged = false;
  double _fraction = 0;

  /// Elapsed time (IQ-TREE's clock) at the first iteration: the search
  /// share is measured from there, without ModelFinder's time.
  int? _searchStart;

  static final _mfStart = RegExp(r'ModelFinder will test up to (\d+)');
  static final _mfLine = RegExp(r'^\s*(\d+)\s+([A-Za-z0-9+{}.,_-]+)\s+(-?\d+\.\d+)\s+\d+\s');
  static final _bestFit = RegExp(r'Best-fit model:\s*(\S+)');
  static final _ufboot = RegExp(r'Generating (\d+) samples for');
  static final _iteration = RegExp(r'Iteration (\d+) / LogL: (-?[\d.]+) / Time: (\S+)(?: \((\S+) left\))?');
  static final _bestScore = RegExp(r'BEST SCORE FOUND : (-?[\d.]+)');
  static final _bootRep = RegExp(r'START BOOTSTRAP REPLICATE NUMBER (\d+)');
  static final _rxStarts = RegExp(r'Starting ML tree search with (\d+) distinct');
  static final _rxSearch = RegExp(r'ML tree search #(\d+), logLikelihood: (-?[\d.]+)');
  static final _rxSearchDone = RegExp(r'ML tree search completed, best tree logLH: (-?[\d.]+)');
  static final _rxBootStart = RegExp(r'Starting bootstrapping analysis with (\d+) replicates');
  static final _rxBoot = RegExp(r'Bootstrap tree #(\d+), logLikelihood');
  static final _rxConverged = RegExp(r'Bootstrapping converged after (\d+) replicates');
  static final _rxPartitions = RegExp(r'Alignment comprises (\d+) partitions');
  static final _mooseModel = RegExp(r'Partition #(\d+): (\S+) \(LogLH');

  /// Overall fraction (0–1), never going back.
  double get fraction => _fraction;

  /// "0h:1m:5s" → 65.
  static int? seconds(String? s) {
    if (s == null) return null;
    final m = RegExp(r'(\d+)h:(\d+)m:(\d+)s').firstMatch(s);
    if (m == null) return null;
    return int.parse(m.group(1)!) * 3600 + int.parse(m.group(2)!) * 60 + int.parse(m.group(3)!);
  }

  void _enter(TreeStage s) {
    if (stages.indexOf(s) > stages.indexOf(stage)) stage = s;
  }

  /// Updates the state with one output line; returns true when something
  /// shown changed.
  bool add(String line) {
    final before = (stage, modelsDone, bestModel, iteration, startsDone, bootDone, bootTotal, _fraction);
    Match? m;
    if ((m = _mfStart.firstMatch(line)) != null) {
      _enter(TreeStage.modelFinder);
      modelsTotal = int.parse(m!.group(1)!);
    } else if (stage == TreeStage.modelFinder && (m = _mfLine.firstMatch(line)) != null) {
      modelsDone = math.max(modelsDone ?? 0, int.parse(m!.group(1)!));
    } else if ((m = _bestFit.firstMatch(line)) != null) {
      bestModel = m!.group(1);
      if (raxml) {
        // ModelFinder before RAxML-NG: the search comes next.
        modelsDone = modelsTotal;
      } else {
        _enter(TreeStage.search);
      }
    } else if ((m = _rxPartitions.firstMatch(line)) != null) {
      partitions = int.parse(m!.group(1)!);
    } else if (line.contains('Starting model selection')) {
      // MOOSE (RAxML-NG 2): one line per partition at the end.
      _enter(TreeStage.modelFinder);
      moose = true;
      modelsTotal = partitions;
      modelsDone = 0;
    } else if (stage == TreeStage.modelFinder && (m = _mooseModel.firstMatch(line)) != null) {
      modelsDone = int.parse(m!.group(1)!) + 1;
      bestModel = modelsTotal == 1 ? m.group(2) : null;
    } else if (line.contains('Model selection time')) {
      modelsDone = modelsTotal;
    } else if ((m = _ufboot.firstMatch(line)) != null) {
      bootTotal = int.parse(m!.group(1)!);
      _enter(TreeStage.search);
    } else if ((m = _iteration.firstMatch(line)) != null) {
      _enter(TreeStage.search);
      iteration = int.parse(m!.group(1)!);
      logL = double.tryParse(m.group(2)!);
      secondsElapsed = seconds(m.group(3));
      secondsLeft = seconds(m.group(4));
      _searchStart ??= secondsElapsed;
    } else if ((m = _bestScore.firstMatch(line)) != null) {
      bestLogL = double.tryParse(m!.group(1)!);
    } else if (line.contains('Creating bootstrap support values') || line.contains('Computing bootstrap consensus tree')) {
      _enter(TreeStage.support);
    } else if ((m = _bootRep.firstMatch(line)) != null) {
      _enter(TreeStage.support);
      bootDone = int.parse(m!.group(1)!);
    } else if ((m = _rxStarts.firstMatch(line)) != null) {
      _enter(TreeStage.search);
      startsTotal = int.parse(m!.group(1)!);
    } else if ((m = _rxSearch.firstMatch(line)) != null) {
      _enter(TreeStage.search);
      startsDone = math.max(startsDone ?? 0, int.parse(m!.group(1)!));
      logL = double.tryParse(m.group(2)!);
    } else if ((m = _rxSearchDone.firstMatch(line)) != null) {
      bestLogL = double.tryParse(m!.group(1)!);
      startsDone = startsTotal;
    } else if ((m = _rxBootStart.firstMatch(line)) != null) {
      _enter(TreeStage.support);
      bootTotal = int.parse(m!.group(1)!);
      bootDone = 0;
    } else if ((m = _rxBoot.firstMatch(line)) != null) {
      _enter(TreeStage.support);
      bootDone = math.max(bootDone ?? 0, int.parse(m!.group(1)!));
    } else if ((m = _rxConverged.firstMatch(line)) != null) {
      bootDone = int.parse(m!.group(1)!);
      bootConverged = true;
    } else if (line.contains('Analysis results written to') || line.contains('Best ML tree saved to') || line.contains('Elapsed time:')) {
      _enter(TreeStage.finish);
    }
    _update();
    return before != (stage, modelsDone, bestModel, iteration, startsDone, bootDone, bootTotal, _fraction);
  }

  /// Share of the run of each stage.
  double _weight(TreeStage s) => switch (s) {
        TreeStage.prepare => 0.03,
        TreeStage.modelFinder => 0.35,
        TreeStage.search => raxml ? 0.25 : 0.52,
        TreeStage.support => raxml ? 0.45 : 0.08,
        TreeStage.finish => 0.02,
      };

  /// Progress inside the current stage (0–1).
  double get stageFraction => switch (stage) {
        TreeStage.prepare => 0.5,
        TreeStage.modelFinder => modelsTotal == null || modelsDone == null ? 0.05 : math.min(1, modelsDone! / modelsTotal!),
        TreeStage.search => raxml
            ? (startsTotal == null ? 0.05 : math.min(1, (startsDone ?? 0) / startsTotal!))
            : secondsLeft != null && secondsElapsed != null && _searchStart != null && secondsElapsed! - _searchStart! + secondsLeft! > 0
                // IQ-TREE estimates the time left; before that, iterations
                // (it stops after 100 without improvement).
                ? math.max((iteration ?? 0) / 400, (secondsElapsed! - _searchStart!) / (secondsElapsed! - _searchStart! + secondsLeft!))
                : math.min(0.25, (iteration ?? 0) / 160),
        TreeStage.support => bootTotal == null || bootDone == null ? 0.3 : math.min(1, bootDone! / bootTotal!),
        TreeStage.finish => 1,
      };

  void _update() {
    final total = stages.fold<double>(0, (a, s) => a + _weight(s));
    var done = 0.0;
    for (final s in stages) {
      if (s == stage) break;
      done += _weight(s);
    }
    final f = (done + _weight(stage) * stageFraction) / total;
    _fraction = math.max(_fraction, math.min(stage == TreeStage.finish ? 1 : 0.99, f));
  }

  void finish() {
    stage = TreeStage.finish;
    _fraction = 1;
  }
}
