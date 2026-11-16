import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/core/tree_progress.dart';

void main() {
  test('IQ-TREE output moves through ModelFinder, search and support', () {
    final p = TreeProgress(stages: TreeStage.values);
    final values = <double>[];
    final lines = [
      'Reading alignment file COI_aligned.fasta ... Fasta format detected',
      'ModelFinder will test up to 968 DNA models (sample size: 1551 epsilon: 0.100) ...',
      ' No. Model         -LnL         df  AIC          AICc         BIC',
      '  1  JC            25444.200    89  51066.399    51077.364    51542.252',
      '139  HKY+R4        18949.087    96  38090.174    38102.983    38603.453',
      'Best-fit model: TIM2+F+R4 chosen according to BIC',
      'Generating 1000 samples for ultrafast bootstrap (seed: 915241)...',
      'Iteration 10 / LogL: -18692.167 / Time: 0h:0m:25s',
      'Iteration 30 / LogL: -18692.203 / Time: 0h:0m:28s (0h:0m:15s left)',
      'Iteration 100 / LogL: -18692.242 / Time: 0h:0m:38s (0h:0m:1s left)',
      'BEST SCORE FOUND : -18692.166',
      'Creating bootstrap support values...',
      'Analysis results written to: ',
    ];
    for (final l in lines) {
      p.add(l);
      values.add(p.fraction);
      if (l.startsWith('139')) {
        expect(p.stage, TreeStage.modelFinder);
        expect(p.modelsDone, 139);
        expect(p.modelsTotal, 968);
      }
      if (l.startsWith('Iteration 30')) {
        expect(p.stage, TreeStage.search);
        expect(p.secondsLeft, 15);
        expect(p.bootTotal, 1000);
      }
    }
    for (var i = 1; i < values.length; i++) {
      expect(values[i], greaterThanOrEqualTo(values[i - 1]));
    }
    expect(p.bestModel, 'TIM2+F+R4');
    expect(p.bestLogL, -18692.166);
    expect(p.stage, TreeStage.finish);
    expect(p.fraction, 1);
  });

  test('RAxML-NG output counts starting trees and bootstrap replicates', () {
    final p = TreeProgress(stages: [TreeStage.prepare, TreeStage.search, TreeStage.support, TreeStage.finish], raxml: true);
    p.add('[00:00:02] Starting ML tree search with 20 distinct starting trees');
    p.add('[00:00:05] [worker #0] ML tree search #4, logLikelihood: -18700.100000');
    expect((p.stage, p.startsDone, p.startsTotal), (TreeStage.search, 4, 20));
    p.add('[00:00:30] ML tree search completed, best tree logLH: -18692.200000');
    p.add('[00:00:30] Starting bootstrapping analysis with 1000 replicates.');
    p.add('[00:00:33] [worker #1] Bootstrap tree #12, logLikelihood: -18800.3');
    expect((p.stage, p.bootDone, p.bootTotal), (TreeStage.support, 12, 1000));
    final mid = p.fraction;
    p.add('[00:05:00] Bootstrapping converged after 300 replicates.');
    expect(p.bootConverged, isTrue);
    expect(p.fraction, greaterThan(mid));
  });

  test('RAxML-NG 2 MOOSE: model selection per partition, then the search', () {
    final p = TreeProgress(stages: [TreeStage.prepare, TreeStage.modelFinder, TreeStage.search, TreeStage.finish], raxml: true);
    p.add('Alignment comprises 2 partitions and 332 patterns');
    p.add('Starting model selection with parsimony starting tree using 2 threads');
    expect((p.stage, p.moose, p.modelsDone, p.modelsTotal), (TreeStage.modelFinder, true, 0, 2));
    p.add('Partition #0: HKY+FO+I+G4m (LogLH = -1782.169337  BIC = 3672.710540)');
    p.add('Partition #1: TIM2+FO+I+G4m (LogLH = -2264.475185  BIC = 4655.031998)');
    expect(p.modelsDone, 2);
    expect(p.bestModel, isNull, reason: 'several partitions: no single best model');
    p.add('[00:00:00] Starting ML tree search with 1 distinct starting trees');
    expect(p.stage, TreeStage.search);
  });
}
