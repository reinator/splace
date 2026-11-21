import 'dart:io';

/// Windows: project files (`.splc`) open in this SPLACE on a double-click.
///
/// The portable program has no installer, so SPLACE registers the type
/// itself for the current user (HKCU, no administrator rights), pointing to
/// the executable that is running. It only writes when the registration is
/// missing or points to another executable (a new version unpacks to a new
/// folder). [classesKey] is the registry key under which the types go
/// (tests use a key of their own).
Future<bool> registerProjectFileType({
  String extension = '.splc',
  String? executable,
  String classesKey = r'HKCU\Software\Classes',
  String progId = 'SPLACE.Project',
}) async {
  if (!Platform.isWindows) return false;
  final exe = executable ?? Platform.resolvedExecutable;
  final command = '"$exe" "%1"';
  try {
    final current = await Process.run('reg', ['query', '$classesKey\\$progId\\shell\\open\\command', '/ve']);
    final ext = await Process.run('reg', ['query', '$classesKey\\$extension', '/ve']);
    if (current.exitCode == 0 && '${current.stdout}'.contains(command) && ext.exitCode == 0 && '${ext.stdout}'.contains(progId)) return true;
    for (final (key, value) in [
      ('$classesKey\\$extension', progId),
      ('$classesKey\\$progId', 'SPLACE project'),
      ('$classesKey\\$progId\\DefaultIcon', '"$exe",0'),
      ('$classesKey\\$progId\\shell\\open\\command', command),
    ]) {
      final r = await Process.run('reg', ['add', key, '/ve', '/d', value, '/f']);
      if (r.exitCode != 0) return false;
    }
    return true;
  } catch (_) {
    return false;
  }
}
