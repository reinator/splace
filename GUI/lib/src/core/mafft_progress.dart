import 'dart:math' as math;

/// Reads MAFFT's stderr and estimates how far one alignment is (0–1).
///
/// MAFFT prints counters for each stage, e.g. (L-INS-i, `--auto`):
///
///     All-to-all alignment.        10 / 46 (by thread 1)
///     Constructing a UPGMA tree    20 / 46
///     Progressive alignment        STEP 12 /45 (thread 0)
///     (iterative refinement)       001-0012-0 (thread 2) identical
///
/// FFT-NS-2 prints the distance/tree/progressive block twice instead of the
/// refinement. The stages before the progressive alignment take 0–35 %,
/// the first progressive alignment 35–70 % and the second one or the
/// refinement 70–98 %. The refinement stops when it converges, so its share
/// grows by half of what is left at each iteration. The value never goes
/// back.
class MafftProgress {
  MafftProgress(this.sequences);

  /// Number of input sequences (for the refinement's branch count).
  final int sequences;

  double _value = 0;
  int _counterRuns = 0, _lastCount = -1, _progressiveRuns = 0, _lastStep = -1;
  bool _inProgressive = false;

  double get value => _value;

  static final _counter = RegExp(r'^\s*(\d+)\s*/\s*(\d+)(?:\s|$|\()');
  static final _step = RegExp(r'^\s*STEP\s+(\d+)\s*/\s*(\d+)');
  static final _refine = RegExp(r'^\s*(\d{3})-(\d{4})-\d');

  /// Whether [line] is only a progress counter (not worth logging).
  static bool isProgressLine(String line) => _counter.hasMatch(line) || _step.hasMatch(line) || _refine.hasMatch(line);

  /// Updates the estimate with one line; returns true when it changed.
  bool add(String line) {
    final before = _value;
    final step = _step.firstMatch(line);
    if (step != null) {
      final n = int.parse(step.group(1)!), total = math.max(1, int.parse(step.group(2)!));
      // A new STEP series starts a new progressive alignment.
      if (!_inProgressive || n < _lastStep - total ~/ 2) {
        _progressiveRuns++;
        _inProgressive = true;
      }
      _lastStep = n;
      final f = n / total;
      _set(_progressiveRuns <= 1 ? 0.35 + 0.35 * f : 0.70 + 0.28 * f);
      return _value != before;
    }
    final refine = _refine.firstMatch(line);
    if (refine != null) {
      final iteration = int.parse(refine.group(1)!), branch = int.parse(refine.group(2)!);
      final branches = math.max(1, 2 * (sequences - 1));
      final done = (iteration - 1) + math.min(1.0, branch / branches);
      _set(0.70 + 0.28 * (1 - math.pow(0.5, done)));
      return _value != before;
    }
    final counter = _counter.firstMatch(line);
    if (counter != null) {
      final n = int.parse(counter.group(1)!), total = math.max(1, int.parse(counter.group(2)!));
      if (n <= _lastCount || _lastCount < 0) _counterRuns++;
      _lastCount = n;
      _inProgressive = false;
      // Distance matrix and tree before the first progressive alignment.
      if (_progressiveRuns == 0) _set(0.35 * math.min(0.95, (_counterRuns - 1 + n / total) / 3));
      return _value != before;
    }
    if (line.contains('Progressive alignment')) _inProgressive = false;
    return false;
  }

  void _set(double v) => _value = math.max(_value, math.min(0.98, v));
}
