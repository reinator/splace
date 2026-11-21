import 'dart:io';

import 'package:flutter_test/flutter_test.dart';
import 'package:path/path.dart' as p;
import 'package:splace/src/services/process_runner.dart';
import 'package:splace/src/services/tool_locator.dart';

/// MAFFT for Windows (mafft.bat + MSYS2) in a folder whose path has spaces,
/// accents and parentheses, as in many user names ("C:\Users\Maria Ação (Lab)
/// \AppData\…"): cmd /c broke at the first space and MSYS cannot open
/// non-ASCII paths. Windows only, with the repository's GUI/tools/win.
void main() {
  final win = p.normalize(p.join(Directory.current.path, 'tools', 'win'));
  final available = Platform.isWindows && File(p.join(win, 'mafft.bat')).existsSync();

  test('bundled MAFFT runs from a folder with spaces, accents and parentheses', () async {
    final tmp = Directory.systemTemp.createTempSync('splace Ação (teste) ');
    addTearDown(() => tmp.deleteSync(recursive: true));
    final tools = Directory(p.join(tmp.path, 'tools'))..createSync();
    File(p.join(win, 'mafft.bat')).copySync(p.join(tools.path, 'mafft.bat'));
    final r = await Process.run('robocopy', [p.join(win, 'usr'), p.join(tools.path, 'usr'), '/E', '/NFL', '/NDL', '/NJH', '/NJS']);
    expect(r.exitCode, lessThan(8), reason: 'robocopy');
    final input = File(p.join(tmp.path, 'entrada ção.fasta'))..writeAsStringSync('>a\nATGGCTTATTCACTGGGA\n>b\nATGGCCTATTCACTGGGA\n>c\nATGGCTAATTCACTGGGA\n');
    final tool = ResolvedTool(ExternalTool.mafft, p.join(tools.path, 'mafft.bat'), ToolRuntime.native);
    expect(tool.msysMafft, isTrue);
    final out = await runTool(tool, ['--auto', input.path], workingDirectory: tmp.path, captureStdout: true);
    expect(out.exitCode, 0);
    expect(out.stdout, contains('>a'));
    expect(out.stdout.toLowerCase(), contains('atggcttattcactggga'));
  }, skip: available ? false : 'needs Windows and tools/win/mafft.bat', timeout: const Timeout(Duration(minutes: 3)));
}
