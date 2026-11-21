import 'dart:io';

import 'package:flutter_test/flutter_test.dart';
import 'package:path/path.dart' as p;
import 'package:splace/src/core/delimitation.dart';
import 'package:splace/src/core/newick.dart';
import 'package:splace/src/services/process_runner.dart';
import 'package:splace/src/services/tool_locator.dart';

/// LSD2 (bundled in GUI/tools) dates a small tree with three clear species
/// and an outgroup; GMYC (tools/gmyc, or Python with ete3 and scipy and
/// tool/gmyc/splace_gmyc.py) delimits them.
void main() {
  test('LSD2 + GMYC delimit the species of a small tree', () async {
    final lsd2 = await ToolLocator().resolve(ExternalTool.lsd2);
    final gmyc = await ToolLocator().resolve(ExternalTool.gmyc);
    if (lsd2 == null || gmyc == null) {
      markTestSkipped('LSD2 or GMYC not found');
      return;
    }
    final tmp = Directory.systemTemp.createTempSync('splace_gmyc_');
    addTearDown(() => tmp.deleteSync(recursive: true));
    // Four tips per species, short branches within, long between.
    String sp(String n) => '((${n}1:0.002,${n}2:0.003):0.001,(${n}3:0.002,${n}4:0.001):0.002)';
    final tree = parseNewick('(((${sp('A')}:0.08,${sp('B')}:0.07):0.05,${sp('C')}:0.1):0.05,(O1:0.01,O2:0.01):0.3);');
    final input = p.join(tmp.path, 'tree.nwk');
    File(input).writeAsStringSync('${mptpNewick(tree)}\n');
    File(p.join(tmp.path, 'og.txt')).writeAsStringSync(lsd2OutgroupFile(['O1', 'O2']));
    final prefix = p.join(tmp.path, 'lsd2');
    final res = await runTool(lsd2, lsd2Args(tree: input, prefix: prefix, seqLength: 1000, outgroupFile: p.join(tmp.path, 'og.txt')), workingDirectory: tmp.path);
    expect(res.exitCode, 0);
    final um = ultrametricFromNexus(File('$prefix.date.nexus').readAsStringSync());
    expect(um, isNotNull);
    final umFile = p.join(tmp.path, 'um.nwk');
    File(umFile).writeAsStringSync('$um\n');
    final json = p.join(tmp.path, 'gmyc.json');
    final g = await runTool(gmyc, gmycArgs(tree: umFile, out: json), workingDirectory: tmp.path, captureStdout: true);
    expect(g.exitCode, 0, reason: g.stdout);
    final d = parseGmyc(File(json).readAsStringSync());
    expect({for (final s in d.species) ...s}, hasLength(12));
  }, tags: ['integration'], timeout: const Timeout(Duration(minutes: 5)));
}
