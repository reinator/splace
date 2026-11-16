import 'dart:async';
import 'dart:convert';
import 'dart:io';

import 'package:path/path.dart' as p;

import 'tool_locator.dart';

class ProcessOutcome {
  ProcessOutcome(this.exitCode, this.stdout, {this.cancelled = false});
  final int exitCode;

  /// Captured stdout (MAFFT writes the alignment there).
  final String stdout;
  final bool cancelled;
  bool get ok => exitCode == 0 && !cancelled;

  /// "exit code N" with what it means when the program crashed.
  String get describe => exitCodeMeaning(exitCode);
}

/// Exit code in words: Windows reports crashes as NTSTATUS codes
/// (0xC0000005 is a memory access violation), POSIX as signals.
String exitCodeMeaning(int code) {
  final u = code & 0xFFFFFFFF;
  final hex = '0x${u.toRadixString(16).toUpperCase().padLeft(8, '0')}';
  final meaning = switch (u) {
    0xC0000005 => 'the program crashed (memory access violation)',
    0xC00000FD => 'the program crashed (stack overflow)',
    0xC0000409 => 'the program crashed (stack buffer overrun / fail-fast)',
    0xC0000374 => 'the program crashed (heap corruption)',
    0xC0000017 || 0xC000012D => 'out of memory',
    0xC000001D => 'illegal instruction: the CPU lacks an instruction set the program was built for',
    0xC0000135 => 'a DLL the program needs was not found',
    0xC000013A => 'stopped by Ctrl+C',
    _ when code == 137 || code == -9 => 'killed (often out of memory)',
    _ when code == 139 || code == -11 => 'the program crashed (segmentation fault)',
    _ when code == 134 || code == -6 => 'the program aborted',
    _ => null,
  };
  if (meaning == null) return 'exit code $code';
  return 'exit code $code ($hex: $meaning)';
}

/// Token used to cancel a running process from the UI.
class CancelToken {
  final _listeners = <void Function()>[];
  bool _cancelled = false;
  bool get isCancelled => _cancelled;

  void cancel() {
    if (_cancelled) return;
    _cancelled = true;
    for (final l in List.of(_listeners)) {
      l();
    }
  }

  void _onCancel(void Function() f) => _listeners.add(f);
}

/// Runs an external tool asynchronously, streaming stderr (and stdout,
/// unless captured) line by line to [onLog]. Never blocks the UI isolate.
Future<ProcessOutcome> runTool(
  ResolvedTool tool,
  List<String> args, {
  required String workingDirectory,
  void Function(String line)? onLog,
  bool captureStdout = false,
  CancelToken? cancel,
}) async {
  onLog?.call('\$ ${tool.description} ${args.join(' ')}');
  // The MSYS of MAFFT for Windows cannot open paths with accents or other
  // non-ASCII letters (a user name such as "João" is in most paths): the
  // input files go to a plain temporary folder, where MAFFT runs; the
  // alignment comes back on stdout.
  Directory? scratch;
  var runArgs = args, runDir = workingDirectory;
  if (tool.msysMafft && (args.any((a) => File(a).existsSync() && !ResolvedTool.plainPath(a)) || !ResolvedTool.plainPath(workingDirectory))) {
    final base = ResolvedTool.msysTempDir(p.dirname(tool.executable));
    scratch = Directory(base).createTempSync('splace_mafft_');
    var n = 0;
    runArgs = [
      for (final a in args)
        if (File(a).existsSync() && !ResolvedTool.plainPath(a))
          (File(a).copySync(p.join(scratch.path, 'in${n++}${p.extension(a)}'))).path
        else
          a,
    ];
    runDir = scratch.path;
  }
  final (exe, finalArgs) = tool.command(runArgs);
  final Process process;
  try {
    process = await Process.start(
      exe,
      finalArgs,
      workingDirectory: runDir,
      environment: tool.environment,
    );
  } catch (_) {
    _clean(scratch);
    rethrow;
  }
  cancel?._onCancel(() => process.kill(ProcessSignal.sigterm));

  final stdoutBuffer = StringBuffer();
  final stdoutDone = captureStdout
      ? process.stdout.transform(const Utf8Decoder(allowMalformed: true)).forEach(stdoutBuffer.write)
      : process.stdout
          .transform(const Utf8Decoder(allowMalformed: true))
          .transform(const LineSplitter())
          .forEach((l) => onLog?.call(l));
  final stderrDone = process.stderr
      .transform(const Utf8Decoder(allowMalformed: true))
      .transform(const LineSplitter())
      .where((l) => l.trim().isNotEmpty && !ResolvedTool.launcherNoise.hasMatch(l))
      .forEach((l) => onLog?.call(l));

  final code = await process.exitCode;
  await Future.wait([stdoutDone, stderrDone]);
  _clean(scratch);
  return ProcessOutcome(code, stdoutBuffer.toString(), cancelled: cancel?.isCancelled ?? false);
}

void _clean(Directory? dir) {
  try {
    if (dir != null && dir.existsSync()) dir.deleteSync(recursive: true);
  } catch (_) {}
}
