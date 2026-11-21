@Tags(['integration'])
@Timeout(Duration(minutes: 10))
library;

import 'dart:io';
import 'dart:math' as math;

import 'package:flutter_test/flutter_test.dart';
import 'package:path/path.dart' as p;
import 'package:splace/src/core/saturation.dart';
import 'package:splace/src/core/tool_commands.dart';
import 'package:splace/src/services/tool_locator.dart';
import 'package:splace/src/state/project_controller.dart';
import 'package:splace/src/ui/methods_text.dart';

/// End-to-end run with the real external tools (bundled in GUI/tools or
/// on PATH). Skipped when MAFFT, trimAl or IQ-TREE are missing. RAxML-NG is
/// exercised when found or when SPLACE_RAXML_NG points to it (a native path
/// or `wsl:/path/to/raxml-ng`).
///
///   flutter test test/integration
void main() {
  late Map<ExternalTool, ResolvedTool?> tools;
  final raxmlOverride = Platform.environment['SPLACE_RAXML_NG'];
  ToolLocator locator() =>
      ToolLocator(overrides: {ExternalTool.raxmlng: ?raxmlOverride});

  setUpAll(() async {
    tools = await locator().resolveAll();
  });

  test('sample-per-file FASTAs → loci → MAFFT → trimAl → saturation → matrix → IQ-TREE', () async {
    for (final t in [ExternalTool.mafft, ExternalTool.trimal, ExternalTool.iqtree]) {
      if (tools[t] == null) {
        markTestSkipped('${t.label} not available');
        return;
      }
    }
    final tmp = Directory.systemTemp.createTempSync('splace_it_');
    addTearDown(() => tmp.deleteSync(recursive: true));

    // 16 taxa, two coding loci with a fast 3rd codon position. One file per
    // sample, header = locus, with a few bases chopped so lengths differ
    // and MAFFT has real work to do.
    final rnd = math.Random(11);
    const bases = 'ACGT';
    String evolve(String s, double t) => String.fromCharCodes(s.codeUnits.indexed.map((e) {
          final rate = [0.3, 0.1, 6.0][e.$1 % 3];
          final change = rnd.nextDouble() < 0.75 * (1 - math.exp(-4 / 3 * rate * t));
          if (!change) return e.$2;
          final others = bases.codeUnits.where((b) => b != e.$2).toList();
          return others[rnd.nextInt(3)];
        }));
    List<String> simulate(int codons) {
      final seed = List.generate(codons * 3, (_) => bases[rnd.nextInt(4)]).join();
      // Remove internal stops so the locus is recognized as coding.
      var level = [seed.replaceAllMapped(RegExp(r'.{3}'), (m) => {'TAA', 'TAG', 'TGA'}.contains(m[0]) ? 'CTG' : m[0]!)];
      for (var d = 0; d < 4; d++) {
        level = [for (final s in level) ...[evolve(s, 0.1), evolve(s, 0.1)]];
      }
      return [
        for (final s in level)
          s.replaceAllMapped(RegExp(r'.{3}'), (m) => {'TAA', 'TAG', 'TGA'}.contains(m[0]) ? 'CTG' : m[0]!),
      ];
    }

    final coi = simulate(200), nd1 = simulate(150);
    final inputs = <String>[];
    for (var i = 0; i < coi.length; i++) {
      final f = File(p.join(tmp.path, 'VOUCHER${i.toString().padLeft(2, '0')}_Genus_species$i.fasta'));
      final trimEnd = i.isEven ? 3 : 0;
      f.writeAsStringSync('>COX1\n${coi[i].substring(0, coi[i].length - trimEnd)}\n>NAD1\n${nd1[i]}\n');
      inputs.add(f.path);
    }

    final c = ProjectController(locator: locator(), autoDetectTools: false)..workDir = p.join(tmp.path, 'project');
    addTearDown(c.dispose);
    c.setThreads(math.min(4, Platform.numberOfProcessors));
    await c.detectTools();

    await c.addFiles(inputs);
    await Future<void>.delayed(const Duration(milliseconds: 50));
    expect(c.organized, isNotNull);
    expect(c.organized!.loci.toSet(), {'COI', 'ND1'}, reason: 'COX1/NAD1 normalized by synonyms');
    expect(c.organized!.samples, hasLength(16));

    await c.useOrganizedLoci();
    expect(File(p.join(c.workDir, '01_loci', 'COI.fasta')).existsSync(), isTrue);

    await c.alignAll();
    for (final l in c.loci.values) {
      expect(l.alignStatus, StepStatus.done, reason: '${l.name}: ${l.message}');
      expect(l.coding, isTrue);
    }

    await c.trimAll();
    for (final l in c.loci.values) {
      expect(l.trimStatus, StepStatus.done, reason: l.name);
      expect(l.trimmed!.first.length % 3, 0, reason: 'codon-preserving trimming keeps whole codons');
      // The kept columns are known, so the viewer can grey the removed ones.
      expect(l.trimMethod, contains('-backtrans'), reason: 'coding genes are trimmed by trimAl in whole codons');
      expect(l.trimmedColumns, isNotNull, reason: '${l.name}: ${l.trimMethod}');
      expect(l.trimmedColumns!.length, l.trimmed!.first.length);
    }

    await c.runSaturation();
    // Each gene is tested on its own (IQ-TREE tree + PhyKIT per codon
    // subset); nothing is concatenated for the test.
    final pooled = c.pooledSaturation!;
    expect(pooled.reports.map((r) => r.name), unorderedEquals(['COI', 'ND1']));
    expect(pooled.verdict(3), SaturationVerdict.saturated);
    expect(pooled.recommendation.exclude, {3});
    expect(File(p.join(c.workDir, '04_saturation', 'saturation_summary.tsv')).existsSync(), isTrue);

    c.applySaturationRecommendation(perLocus: false);
    await c.buildMatrix();
    final m = c.matrix!;
    expect(m.taxa, hasLength(16));
    expect(m.partitions.map((x) => x.name), containsAll(['COI_pos1', 'COI_pos2', 'ND1_pos1', 'ND1_pos2']));
    expect(m.partitions.any((x) => x.name.endsWith('pos3')), isFalse);

    await c.runIqtree(IqtreeOptions(model: 'GTR+F+G4', ufboot: 1000, alrt: 0, threads: '2', seed: 1, outgroup: [m.taxa.first]));
    final run = c.treeRuns.first;
    expect(run.status, StepStatus.done, reason: c.logs.reversed.take(15).map((e) => e.message).join('\n'));
    expect(File(run.treeFile!).readAsStringSync(), contains(m.taxa.last));

    if (tools[ExternalTool.raxmlng] != null) {
      await c.runRaxml(RaxmlOptions(bootstraps: '50', maxThreads: 2, outgroup: [m.taxa.first]), models: RaxmlModels.modelFinder);
      final rax = c.treeRuns.first;
      expect(rax.status, StepStatus.done, reason: c.logs.reversed.take(15).map((e) => e.message).join('\n'));
      expect(rax.outputs.keys, containsAll(['treeFBP', 'treeTBE', 'bestTree']));
      expect(File(p.join(rax.directory, 'partitions_modelfinder.txt')).readAsStringSync(), contains('COI_pos1'));

      // RAxML-NG 2: MOOSE picks one model per partition, fast search,
      // supports computed on the ML tree only.
      await c.runRaxml(RaxmlOptions(maxThreads: 2, fast: true, supportMetrics: 'sh,ebg', outgroup: [m.taxa.first]));
      final moose = c.treeRuns.first;
      expect(moose.status, StepStatus.done, reason: c.logs.reversed.take(15).map((e) => e.message).join('\n'));
      expect(moose.models, RaxmlModels.moose);
      expect(moose.mooseModels, hasLength(m.partitions.length), reason: 'one MOOSE model per partition (partition file found)');
      expect(moose.outputs.keys, containsAll(['treeSH', 'treeEBG']));
      expect(moose.command, contains('auto='));
      expect(moose.command, isNot(contains('--bs-trees')));
      final methods = buildMethodsText(c, 'en');
      final rxText = methods.paragraphs.firstWhere((p) => p.contains('RAxML-NG'));
      expect(rxText, allOf(contains('MOOSE'), contains('--fast'), contains('SH-like aLRT'), contains('Educated Bootstrap Guesser')));
      expect(methods.cited, containsAll(['RAxML-NG', 'RAxML-NG fast', 'SH-aLRT', 'EBG']));
      // ignore: avoid_print
      print(rxText);

      // Both programs at once, threads split between them.
      final before = c.treeRuns.length;
      c.threads = 4;
      await c.runBoth(
        IqtreeOptions(model: 'GTR+F+G4', ufboot: 1000, alrt: 1000, bnni: true, outgroup: [m.taxa.first]),
        RaxmlOptions(fast: true, supportMetrics: 'sh', outgroup: [m.taxa.first]),
      );
      final both = c.treeRuns.take(c.treeRuns.length - before).toList();
      expect(both.map((r) => r.engine).toSet(), {'IQ-TREE', 'RAxML-NG'});
      expect(both.every((r) => r.status == StepStatus.done), isTrue, reason: c.logs.reversed.take(20).map((e) => e.message).join('\n'));
      expect(both.firstWhere((r) => r.engine == 'IQ-TREE').command, allOf(contains('-T 2'), contains('--bnni')));
      expect(both.firstWhere((r) => r.engine == 'RAxML-NG').command, contains('auto{2}'));
    }
  });
}
