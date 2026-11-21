import 'dart:io';

import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/services/file_association.dart';

/// Windows: SPLACE registers `.splc` for the current user (here under a
/// test key of its own, removed at the end).
void main() {
  test('project files are registered to open in SPLACE', () async {
    final key = 'HKCU\\Software\\SPLACE_test_$pid';
    addTearDown(() => Process.runSync('reg', ['delete', key, '/f']));
    const exe = r'C:\Users\Maria Silva\AppData\Local\SPLACE\portable\4.2.2\splace.exe';
    expect(await registerProjectFileType(executable: exe, classesKey: key), isTrue);
    String value(String sub) => '${Process.runSync('reg', ['query', '$key\\$sub', '/ve']).stdout}';
    expect(value('.splc'), contains('SPLACE.Project'));
    expect(value(r'SPLACE.Project\shell\open\command'), contains('"$exe" "%1"'));
    expect(value(r'SPLACE.Project\DefaultIcon'), contains('"$exe",0'));
    // Already registered: nothing to write.
    expect(await registerProjectFileType(executable: exe, classesKey: key), isTrue);
  }, skip: !Platform.isWindows);
}
