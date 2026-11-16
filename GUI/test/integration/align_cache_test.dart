@Tags(['integration'])
@Timeout(Duration(minutes: 5))
library;

import 'dart:io';
import 'dart:math' as math;

import 'package:flutter_test/flutter_test.dart';
import 'package:path/path.dart' as p;
import 'package:splace/src/services/tool_locator.dart';
import 'package:splace/src/state/project_controller.dart';

/// MAFFT/trimAl with the real tools: each marker reports its own progress
/// read from MAFFT, and a second run reuses unchanged markers.
void main() {
  test('per-marker progress and alignment/trimming cache', () async {
    final tools = await ToolLocator().resolveAll();
    if (tools[ExternalTool.mafft] == null || tools[ExternalTool.trimal] == null) {
      markTestSkipped('MAFFT/trimAl not available');
      return;
    }
    final tmp = Directory.systemTemp.createTempSync('splace_cache_');
    addTearDown(() => tmp.deleteSync(recursive: true));
    final rnd = math.Random(3);
    String seq(int n) => List.generate(n, (_) => 'ACGT'[rnd.nextInt(4)]).join().replaceAllMapped(
        RegExp(r'.{3}'), (m) => {'TAA', 'TAG', 'TGA'}.contains(m[0]) ? 'CTG' : m[0]!);
    final coi = seq(600), nd1 = seq(450);
    String mutate(String s) => String.fromCharCodes(s.codeUnits.map((c) => rnd.nextDouble() < 0.05 ? 'ACGT'.codeUnitAt(rnd.nextInt(4)) : c));
    final inputs = <String>[];
    for (var i = 0; i < 30; i++) {
      final f = File(p.join(tmp.path, 'S${i.toString().padLeft(2, '0')}_Genus_species$i.fasta'));
      f.writeAsStringSync('>COX1\n${mutate(coi).substring(0, 600 - 3 * (i % 3))}\n>NAD1\n${mutate(nd1)}\n');
      inputs.add(f.path);
    }
    final c = ProjectController(autoDetectTools: false)..workDir = p.join(tmp.path, 'project');
    addTearDown(c.dispose);
    await c.detectTools();
    await c.addFiles(inputs);
    await Future<void>.delayed(const Duration(milliseconds: 50));
    await c.useOrganizedLoci();

    // Progress values seen for each marker while aligning.
    final seen = <String, Set<double>>{};
    void listen() {
      for (final e in c.modal?.items.entries ?? const <MapEntry<String, ModalItem>>[]) {
        seen.putIfAbsent(e.key, () => {}).add(e.value.progress);
      }
    }

    c.addListener(listen);
    await c.alignAll();
    c.removeListener(listen);
    for (final l in c.loci.values) {
      expect(l.alignStatus, StepStatus.done, reason: '${l.name}: ${l.message}');
    }
    expect(c.modal!.items.values.every((i) => i.state == ItemState.done), isTrue);
    for (final e in seen.entries) {
      final mid = e.value.where((v) => v > 0.05 && v < 0.95).length;
      expect(mid, greaterThan(3), reason: 'MAFFT counters move the bar of ${e.key}');
    }

    // Second run: nothing changed, everything from the cache.
    c.closeModal();
    final watch = Stopwatch()..start();
    await c.alignAll();
    expect(c.modal!.items.values.every((i) => i.state == ItemState.cached), isTrue);
    expect(watch.elapsedMilliseconds, lessThan(3000));
    expect(c.loci.values.every((l) => l.alignStatus == StepStatus.done && l.aligned != null), isTrue);

    await c.trimAll();
    expect(c.loci.values.every((l) => l.trimStatus == StepStatus.done), isTrue);
    c.closeModal();
    await c.trimAll();
    expect(c.modal!.items.values.every((i) => i.state == ItemState.cached), isTrue);

    // Re-aligning keeps the cached trimming of unchanged markers.
    c.closeModal();
    await c.alignAll();
    expect(c.loci.values.every((l) => l.trimmed != null), isTrue);
  });
}
