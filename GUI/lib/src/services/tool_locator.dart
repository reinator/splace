import 'dart:convert';
import 'dart:io';

import 'package:path/path.dart' as p;

enum ExternalTool { mafft, trimal, iqtree, raxmlng, clipkit, phykit, mptp, asap, abgd, lsd2, gmyc }

extension ExternalToolInfo on ExternalTool {
  String get label => switch (this) {
        ExternalTool.mafft => 'MAFFT',
        ExternalTool.trimal => 'trimAl',
        ExternalTool.iqtree => 'IQ-TREE',
        ExternalTool.raxmlng => 'RAxML-NG',
        ExternalTool.clipkit => 'ClipKIT',
        ExternalTool.phykit => 'PhyKIT',
        ExternalTool.mptp => 'mPTP',
        ExternalTool.asap => 'ASAP',
        ExternalTool.abgd => 'ABGD',
        ExternalTool.lsd2 => 'LSD2',
        ExternalTool.gmyc => 'GMYC',
      };

  /// Python package of the tools run as `python -m <module>`.
  String? get pythonModule => switch (this) {
        ExternalTool.clipkit => 'clipkit',
        ExternalTool.phykit => 'phykit',
        _ => null,
      };

  /// Executable names in order of preference.
  List<String> get executableNames => switch (this) {
        ExternalTool.mafft => ['mafft'],
        ExternalTool.trimal => ['trimal'],
        ExternalTool.iqtree => ['iqtree3', 'iqtree2', 'iqtree'],
        ExternalTool.raxmlng => ['raxml-ng'],
        ExternalTool.clipkit => ['clipkit'],
        ExternalTool.phykit => ['phykit'],
        ExternalTool.mptp => ['mptp'],
        ExternalTool.asap => ['asap'],
        ExternalTool.abgd => ['abgd'],
        ExternalTool.lsd2 => ['lsd2'],
        ExternalTool.gmyc => ['gmyc'],
      };

  List<String> get versionArgs => ['--version'];

  /// ASAP and ABGD print no version (their usage has numbers that would be
  /// read as one).
  bool get reportsVersion => this != ExternalTool.asap && this != ExternalTool.abgd && this != ExternalTool.lsd2;

  /// Python script run when no packaged program is found (GMYC: the pGMYC
  /// command line in tool/gmyc), with the Python modules it needs.
  (String, List<String>)? get pythonScript => switch (this) {
        ExternalTool.gmyc => ('gmyc/splace_gmyc.py', ['ete3', 'scipy']),
        _ => null,
      };

  String get installHint => switch (this) {
        ExternalTool.raxmlng =>
          'conda install -c bioconda raxml-ng  (Windows: SPLACE ships raxml-ng.exe in tools/; a raxml-ng inside WSL also works)',
        ExternalTool.iqtree => 'conda install -c bioconda iqtree',
        ExternalTool.mafft => 'conda install -c bioconda mafft',
        ExternalTool.trimal => 'conda install -c bioconda trimal',
        ExternalTool.clipkit => 'pip install clipkit  (or the "Install ClipKIT" button in step 4)',
        ExternalTool.phykit => 'pip install phykit  (or the "Install PhyKIT" button in step 5)',
        ExternalTool.mptp => 'conda install -c bioconda mptp  (SPLACE ships mptp in tools/; binaries at github.com/Pas-Kapli/mptp/releases)',
        ExternalTool.asap || ExternalTool.abgd => 'SPLACE ships it in tools/ (built by GUI/tool/delimitation/build.sh from the iTaxoTools sources)',
        ExternalTool.lsd2 => 'SPLACE ships it in tools/ (official binaries at github.com/tothuhien/lsd2)',
        ExternalTool.gmyc => 'packaged in tools/gmyc, or Python 3 with: pip install ete3 scipy',
      };

  /// Optional tools do not produce warnings when missing.
  bool get optional => this == ExternalTool.clipkit || this == ExternalTool.phykit || this == ExternalTool.mptp || this == ExternalTool.asap || this == ExternalTool.abgd || this == ExternalTool.lsd2 || this == ExternalTool.gmyc;
}

/// How a tool is executed.
enum ToolRuntime { native, wsl }

class ResolvedTool {
  ResolvedTool(this.tool, this.executable, this.runtime, {this.version, this.prefixArgs = const []});
  final ExternalTool tool;

  /// Native path, or the Linux path inside WSL for [ToolRuntime.wsl].
  final String executable;
  final ToolRuntime runtime;
  final String? version;

  /// Arguments placed before the tool's own (`-m clipkit` for Python).
  final List<String> prefixArgs;

  bool get _isBatch => Platform.isWindows && RegExp(r'\.(bat|cmd)$', caseSensitive: false).hasMatch(executable);

  /// Executable + arguments ready for [Process.start]. Windows `.bat`
  /// launchers (bundled MAFFT) go through `cmd /c`; for WSL, Windows paths
  /// in [args] are translated to `/mnt/<drive>/...`.
  (String, List<String>) command(List<String> args) {
    final all = [...prefixArgs, ...args];
    if (runtime == ToolRuntime.wsl) return ('wsl.exe', ['-e', executable, ...all.map(windowsToWslPath)]);
    // The bundled MAFFT: its own bash, straight (see [_msysBash]); Windows
    // paths with forward slashes, which MSYS reads as they are.
    if (_msysBash case final bash?) {
      return (bash, ['/usr/bin/mafft', for (final a in all) RegExp(r'^[A-Za-z]:[\\/]').hasMatch(a) ? a.replaceAll('\\', '/') : a]);
    }
    if (_isBatch) return ('cmd', ['/c', executable, ...all]);
    return (executable, all);
  }

  /// MAFFT for Windows is `mafft.bat` plus an MSYS2 system in `usr\`. The
  /// .bat goes through `cmd /c`, which breaks when the folder has spaces or
  /// other special characters (cmd drops the first and last quotes of the
  /// line: "C:\Users\Maria Silva\…" becomes "C:\Users\Maria"), and it
  /// switches the console code page. So its bash is run directly, with the
  /// variables the .bat would set ([environment]).
  /// The bundled MAFFT for Windows (run through its MSYS bash).
  bool get msysMafft => _msysBash != null;

  String? get _msysBash {
    if (!Platform.isWindows || runtime != ToolRuntime.native || p.basename(executable).toLowerCase() != 'mafft.bat') return null;
    final bash = p.join(p.dirname(executable), 'usr', 'bin', 'bash.exe');
    return File(bash).existsSync() ? bash : null;
  }

  /// Variables for the process (null: inherit as they are).
  Map<String, String>? get environment {
    if (_msysBash == null) return null;
    final root = p.dirname(executable);
    final tmp = msysTempDir(root);
    return {
      'PATH': '${p.join(root, 'usr', 'bin')};${Platform.environment['PATH'] ?? ''}',
      'MAFFT_BINARIES': '/usr/lib/mafft',
      'TMPDIR': tmp,
      'MAFFT_TMPDIR': tmp,
    };
  }

  /// A temporary folder MAFFT's scripts can use: the system one when its
  /// path is plain (ASCII, no spaces), else one under the tools folder or
  /// ProgramData, whichever is plain and writable. Forward slashes.
  /// ASCII and no spaces: a path MSYS reads as it is.
  static bool plainPath(String path) => RegExp(r'^[A-Za-z]:[\\/][\x21-\x7e]*$').hasMatch(path) && !path.contains(' ');

  static String msysTempDir(String toolsRoot) {
    const plain = plainPath;
    final drive = Platform.environment['SystemDrive'] ?? 'C:';
    for (final dir in [Directory.systemTemp.path, p.join(toolsRoot, 'tmp'), p.join('$drive\\', 'ProgramData', 'SPLACE', 'tmp'), p.join('$drive\\', 'SPLACE_tmp')]) {
      if (!plain(dir)) continue;
      try {
        Directory(dir).createSync(recursive: true);
        final probe = File(p.join(dir, '.splace_probe'))..writeAsStringSync('x');
        probe.deleteSync();
        return dir.replaceAll('\\', '/');
      } catch (_) {
        continue;
      }
    }
    return Directory.systemTemp.path.replaceAll('\\', '/');
  }

  /// Lines printed by launchers rather than the tool (bundled mafft.bat).
  static final launcherNoise = RegExp(r'65001|Preparing environment|anti.virus|could not find /tmp', caseSensitive: false);

  String get description =>
      runtime == ToolRuntime.wsl ? '$executable (WSL)' : [executable, ...prefixArgs].join(' ');

  bool get bundled => ToolLocator.bundledBinDirs().any((d) => p.isWithin(p.normalize(d), p.normalize(runtime == ToolRuntime.wsl ? '' : executable)));
}

/// `C:\Users\x\run` -> `/mnt/c/Users/x/run`, also after a `key=` prefix
/// (RAxML-NG's `auto=C:\…\partitions.txt`); other strings unchanged.
String windowsToWslPath(String value) {
  final m = RegExp(r'^([A-Za-z][\w-]*=)?([A-Za-z]):[\\/](.*)$').firstMatch(value);
  if (m == null) return value;
  return '${m.group(1) ?? ''}/mnt/${m.group(2)!.toLowerCase()}/${m.group(3)!.replaceAll('\\', '/')}';
}

String get platformBinFolder => Platform.isWindows ? 'win' : (Platform.isMacOS ? 'mac' : 'linux');

/// Finds external tools: user override → bundled `bin/<platform>` →
/// PATH → (Windows only) WSL.
class ToolLocator {
  ToolLocator({Map<ExternalTool, String>? overrides}) : overrides = overrides ?? {};

  final Map<ExternalTool, String> overrides;
  final Map<ExternalTool, ResolvedTool?> _cache = {};

  void clearCache() => _cache.clear();

  Future<ResolvedTool?> resolve(ExternalTool tool) async {
    if (_cache.containsKey(tool)) return _cache[tool];
    final found = await _find(tool);
    _cache[tool] = found;
    return found;
  }

  Future<Map<ExternalTool, ResolvedTool?>> resolveAll() async =>
      {for (final t in ExternalTool.values) t: await resolve(t)};

  Future<ResolvedTool?> _find(ExternalTool tool) async {
    final override = overrides[tool];
    if (override != null && override.isNotEmpty) {
      if (override.startsWith('wsl:')) {
        return _withVersion(ResolvedTool(tool, override.substring(4), ToolRuntime.wsl));
      }
      if (File(override).existsSync()) return _withVersion(ResolvedTool(tool, override, ToolRuntime.native));
    }
    for (final dir in bundledBinDirs()) {
      for (final name in tool.executableNames) {
        // tools/<name>, or tools/<name>/<name> (programs packed with their
        // own folder, as ClipKIT and PhyKIT).
        for (final candidate in [..._withExtensions(p.join(dir, name)), ..._withExtensions(p.join(dir, name, name))]) {
          if (File(candidate).existsSync()) {
            return _withVersion(ResolvedTool(tool, candidate, ToolRuntime.native));
          }
        }
      }
    }
    // Linux builds shipped for WSL: a fallback when a Windows tool is
    // missing from tools/ (SPLACE ships its own raxml-ng.exe).
    if (Platform.isWindows) {
      for (final dir in bundledWslDirs()) {
        for (final name in tool.executableNames) {
          final candidate = p.join(dir, name);
          if (File(candidate).existsSync() && await wslAvailable()) {
            return _withVersion(ResolvedTool(tool, windowsToWslPath(p.absolute(candidate)), ToolRuntime.wsl));
          }
        }
      }
    }
    for (final name in tool.executableNames) {
      final onPath = _searchPath(name);
      if (onPath != null) return _withVersion(ResolvedTool(tool, onPath, ToolRuntime.native));
    }
    if (tool.pythonScript case (final script, final modules)?) {
      final path = _scriptPath(script);
      final py = path == null ? null : await findPython();
      if (path != null && py != null && await _pythonHas(py, modules.join(', '))) {
        return _withVersion(ResolvedTool(tool, py.executable, ToolRuntime.native, prefixArgs: [...py.args, path]));
      }
      return null;
    }
    if (tool.pythonModule case final module?) {
      final py = await findPython();
      if (py != null && await _pythonHas(py, module)) {
        return _withVersion(ResolvedTool(tool, py.executable, ToolRuntime.native, prefixArgs: [...py.args, '-m', module]));
      }
      return null;
    }
    if (Platform.isWindows) {
      for (final name in tool.executableNames) {
        final wslPath = await _searchWsl(name);
        if (wslPath != null) return _withVersion(ResolvedTool(tool, wslPath, ToolRuntime.wsl));
      }
    }
    return null;
  }

  /// Folder with the external tools shipped with the app: `tools/` next to
  /// the executable (Windows/Linux bundles) or in the macOS Resources.
  static String get bundledToolsDir {
    final exeDir = p.dirname(Platform.resolvedExecutable);
    return Platform.isMacOS ? p.normalize(p.join(exeDir, '..', 'Resources', 'tools')) : p.join(exeDir, 'tools');
  }

  /// Candidate folders with bundled binaries: the app's `tools/` folder for
  /// installed builds, and the repository's `GUI/tools/<platform>` during
  /// development.
  static List<String> bundledBinDirs() {
    final dirs = <String>[];
    final env = Platform.environment['SPLACE_BIN_DIR'];
    if (env != null) dirs.add(env);
    dirs.add(bundledToolsDir);
    // Development: walk up from the working directory to the repository.
    var dir = Directory.current.absolute;
    for (var i = 0; i < 6; i++) {
      dirs.add(p.join(dir.path, 'tools', platformBinFolder));
      dirs.add(p.join(dir.path, 'GUI', 'tools', platformBinFolder));
      final parent = dir.parent;
      if (parent.path == dir.path) break;
      dir = parent;
    }
    return dirs;
  }

  /// Folders with Linux binaries run through WSL on Windows.
  static List<String> bundledWslDirs() => [
        p.join(bundledToolsDir, 'wsl'),
        for (final d in bundledBinDirs())
          if (d.endsWith(p.join('tools', 'win'))) p.join(p.dirname(d), 'linux'),
      ];

  static bool? _wsl;

  static Future<bool> wslAvailable() async {
    if (_wsl != null) return _wsl!;
    try {
      final r = await Process.run('wsl.exe', ['-e', 'true']).timeout(const Duration(seconds: 20));
      _wsl = r.exitCode == 0;
    } catch (_) {
      _wsl = false;
    }
    return _wsl!;
  }

  /// [relative] (under `tool/`) in the app's tools folder or in the
  /// repository's GUI/tool (development).
  static String? _scriptPath(String relative) {
    final candidates = [p.join(bundledToolsDir, relative)];
    var dir = Directory.current.absolute;
    for (var i = 0; i < 6; i++) {
      candidates
        ..add(p.join(dir.path, 'tool', relative))
        ..add(p.join(dir.path, 'GUI', 'tool', relative));
      final parent = dir.parent;
      if (parent.path == dir.path) break;
      dir = parent;
    }
    return candidates.where((c) => File(c).existsSync()).firstOrNull;
  }

  static Future<bool> _pythonHas(PythonCommand py, String module) async {
    try {
      final r = await Process.run(py.executable, [...py.args, '-c', 'import $module'])
          .timeout(const Duration(seconds: 20));
      return r.exitCode == 0;
    } catch (_) {
      return false;
    }
  }

  static List<String> _withExtensions(String base) =>
      Platform.isWindows ? ['$base.exe', '$base.bat', '$base.cmd'] : [base];

  static String? _searchPath(String name) {
    final pathVar = Platform.environment[Platform.isWindows ? 'Path' : 'PATH'] ??
        Platform.environment['PATH'] ??
        '';
    for (final dir in pathVar.split(Platform.isWindows ? ';' : ':')) {
      if (dir.isEmpty) continue;
      for (final candidate in _withExtensions(p.join(dir, name))) {
        if (File(candidate).existsSync()) return candidate;
      }
    }
    return null;
  }

  /// Uses an interactive login shell so conda's PATH from `.bashrc` is
  /// visible.
  static Future<String?> _searchWsl(String name) async {
    try {
      final r = await Process.run('wsl.exe', ['-e', 'bash', '-lic', 'command -v $name'],
          stdoutEncoding: utf8, stderrEncoding: utf8)
          .timeout(const Duration(seconds: 20));
      if (r.exitCode != 0) return null;
      final line = (r.stdout as String)
          .split('\n')
          .map((l) => l.trim())
          .lastWhere((l) => l.startsWith('/'), orElse: () => '');
      return line.isEmpty ? null : line;
    } catch (_) {
      return null;
    }
  }

  static Future<ResolvedTool> _withVersion(ResolvedTool t) async {
    if (!t.tool.reportsVersion) return t;
    if (t.tool.pythonModule != null && t.prefixArgs.isNotEmpty) {
      final code = t.tool == ExternalTool.phykit
          ? 'from phykit.version import __version__ as v; print(v)'
          : 'import clipkit; print(getattr(clipkit, "__version__", ""))';
      try {
        final r = await Process.run(t.executable, [...t.prefixArgs.take(t.prefixArgs.length - 2), '-c', code])
            .timeout(const Duration(seconds: 20));
        final v = (r.stdout as String).trim();
        return ResolvedTool(t.tool, t.executable, t.runtime, prefixArgs: t.prefixArgs, version: v.isEmpty ? null : v);
      } catch (_) {
        return t;
      }
    }
    try {
      final (exe, args) = t.command(t.tool.versionArgs);
      final r = await Process.run(exe, args, environment: t.environment, stdoutEncoding: utf8, stderrEncoding: utf8)
          .timeout(const Duration(seconds: 20));
      final text = '${r.stdout}\n${r.stderr}';
      final version = RegExp(r'v?(\d+\.\d+(?:\.\d+)?(?:[.\-]?rev\d+)?)').firstMatch(text)?.group(1);
      return ResolvedTool(t.tool, t.executable, t.runtime, version: version, prefixArgs: t.prefixArgs);
    } catch (_) {
      return t;
    }
  }
}

class PythonCommand {
  const PythonCommand(this.executable, this.args, this.version);
  final String executable;
  final List<String> args;
  final String version;
  String get display => [executable, ...args].join(' ');
}

/// First working Python 3 (`py -3` on Windows, then python/python3).
Future<PythonCommand?> findPython() async {
  final candidates = Platform.isWindows
      ? [('py', ['-3']), ('python', <String>[]), ('python3', <String>[])]
      : [('python3', <String>[]), ('python', <String>[])];
  for (final (exe, args) in candidates) {
    try {
      final r = await Process.run(exe, [...args, '--version'], stdoutEncoding: utf8, stderrEncoding: utf8)
          .timeout(const Duration(seconds: 15));
      final text = '${r.stdout}${r.stderr}'.trim();
      if (r.exitCode == 0 && RegExp(r'Python\s+3', caseSensitive: false).hasMatch(text)) {
        return PythonCommand(exe, args, text.replaceFirst(RegExp(r'^Python\s*', caseSensitive: false), ''));
      }
    } catch (_) {}
  }
  return null;
}
