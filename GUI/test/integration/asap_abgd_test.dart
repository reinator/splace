import 'dart:io';
import 'dart:math' as math;

import 'package:flutter_test/flutter_test.dart';
import 'package:path/path.dart' as p;
import 'package:splace/src/core/delimitation.dart';
import 'package:splace/src/services/process_runner.dart';
import 'package:splace/src/services/tool_locator.dart';

/// The bundled ASAP and ABGD (`GUI/tools/<platform>`, built by
/// tool/delimitation/build.sh) on 12 sequences of three clearly separated
/// species: the SPART files they write and how SPLACE reads them.
void main() {
  // Three species, 60 differences between them, 2–3 within.
  final rnd = math.Random(1);
  String mutate(String s, int n) {
    final b = s.split('');
    for (var k = 0; k < n; k++) {
      final i = rnd.nextInt(b.length);
      b[i] = 'ACGT'.replaceAll(b[i], '')[rnd.nextInt(3)];
    }
    return b.join();
  }

  final base = List.generate(600, (_) => 'ACGT'[rnd.nextInt(4)]).join();
  final records = <(String, String)>[
    for (final sp in ['Alpha_one', 'Beta_two', 'Gamma_three'])
      for (final anc in [mutate(base, 60)])
        for (var k = 1; k <= 4; k++) ('NC_00$k.1_$sp', mutate(anc, 2)),
  ];

  for (final (tool, method) in [(ExternalTool.asap, DelimMethod.asap), (ExternalTool.abgd, DelimMethod.abgd)]) {
    test('${tool.label} finds the three species', () async {
      final resolved = await ToolLocator().resolve(tool);
      if (resolved == null) {
        markTestSkipped('${tool.label} not found');
        return;
      }
      final tmp = Directory.systemTemp.createTempSync('splace_${method.name}_');
      addTearDown(() => tmp.deleteSync(recursive: true));
      final (fasta, names) = delimitationFasta(records);
      File(p.join(tmp.path, 'input.fasta')).writeAsStringSync(fasta);
      final args = method == DelimMethod.asap
          ? asapArgs(fasta: 'input.fasta', outDir: 'out', distance: DelimDistance.k2p)
          : abgdArgs(fasta: 'input.fasta', outDir: 'out', distance: DelimDistance.k2p);
      final res = await runTool(resolved, args, workingDirectory: tmp.path);
      expect(res.exitCode, 0);
      final sparts = Directory(p.join(tmp.path, 'out')).listSync().whereType<File>().where((f) => f.path.endsWith('.spart')).toList();
      expect(sparts, isNotEmpty, reason: 'no SPART file');
      final parts = parseSpart(sparts.first.readAsStringSync(), names);
      expect(parts, isNotEmpty);
      // Some partition has exactly the three species, with the real names.
      final three = parts.where((x) => x.species.length == 3).firstOrNull;
      expect(three, isNotNull, reason: [for (final x in parts) x.species.length].join(','));
      expect({for (final s in three!.species) s.map((n) => n.split('_').skip(2).join('_')).toSet().single}, {'Alpha_one', 'Beta_two', 'Gamma_three'});
    }, tags: ['integration'], timeout: const Timeout(Duration(minutes: 3)));
  }
}
