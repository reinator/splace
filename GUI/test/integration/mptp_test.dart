import 'dart:io';

import 'package:flutter_test/flutter_test.dart';
import 'package:path/path.dart' as p;
import 'package:splace/src/core/delimitation.dart';
import 'package:splace/src/core/newick.dart';
import 'package:splace/src/services/process_runner.dart';
import 'package:splace/src/services/tool_locator.dart';

/// The bundled mPTP (`GUI/tools/<platform>`, or one on PATH) on a small
/// tree with three clear species and a two-tip outgroup: the PTP and mPTP
/// models, `--minbr_auto` and the parsing of the output.
void main() {
  test('mPTP delimits the species of a small tree', () async {
    final tool = await ToolLocator().resolve(ExternalTool.mptp);
    if (tool == null) {
      markTestSkipped('mPTP not found');
      return;
    }
    final tmp = Directory.systemTemp.createTempSync('splace_mptp_');
    addTearDown(() => tmp.deleteSync(recursive: true));
    // Short branches within species, long ones between them; IQ-TREE-style
    // support labels, which mptpNewick removes.
    final tree = parseNewick('((((A1:0.001,A2:0.002)100/100:0.08,(B1:0.001,B2:0.0015,B3:0.001)99/100:0.09)95/98:0.1,'
        '(C1:0.002,C2:0.001)100/100:0.12)80/90:0.05,D1:0.3,D2:0.28);');
    final input = p.join(tmp.path, 'tree.nwk');
    File(input).writeAsStringSync('${mptpNewick(tree)}\n');
    String seq(String core) => '$core${'ACGT' * 10}';
    File(p.join(tmp.path, 'm.fasta')).writeAsStringSync([
      for (final (n, s) in [
        ('A1', seq('AAAAAAAAAA')), ('A2', seq('AAAAAAAAAC')), ('B1', seq('CCCCCCCCCC')), ('B2', seq('CCCCCCCCCA')), ('B3', seq('CCCCCCCCCC')),
        ('C1', seq('GGGGGGGGGG')), ('C2', seq('GGGGGGGGGA')), ('D1', seq('TTTTTTTTTT')), ('D2', seq('TTTTTTTTTA')),
      ])
        '>$n\n$s',
    ].join('\n'));

    final minbr = await runTool(tool, mptpMinbrArgs(tree: input, fasta: p.join(tmp.path, 'm.fasta'), prefix: p.join(tmp.path, 'minbr')),
        workingDirectory: tmp.path);
    expect(minbr.exitCode, 0);
    expect(parseMinbr(File(p.join(tmp.path, 'minbr.txt')).readAsStringSync()), isNotNull);

    for (final model in PtpModel.values) {
      final prefix = p.join(tmp.path, model.name);
      final out = await runTool(tool, mptpArgs(tree: input, prefix: prefix, model: model, outgroup: ['D1', 'D2'], cropOutgroup: true),
          workingDirectory: tmp.path);
      expect(out.exitCode, 0, reason: model.label);
      final d = parseMptp(File('$prefix.txt').readAsStringSync(), model);
      // The outgroup is cropped; every ingroup tip is in one species.
      expect({for (final s in d.species) ...s}, {'A1', 'A2', 'B1', 'B2', 'B3', 'C1', 'C2'}, reason: model.label);
      expect(d.species.length, greaterThanOrEqualTo(1));
    }
  }, tags: ['integration'], timeout: const Timeout(Duration(minutes: 2)));
}
