import 'dart:async';
import 'dart:convert';
import 'dart:io';

import 'package:flutter/foundation.dart';

/// Bytes the app sent and received over HTTP (NCBI, GBIF, BLAST…), counted
/// by the services as they read responses.
class NetMeter {
  static int received = 0, sent = 0;

  static void add({int received = 0, int sent = 0}) {
    NetMeter.received += received;
    NetMeter.sent += sent;
  }
}

/// One reading of the machine and the app.
class ResourceSnapshot {
  const ResourceSnapshot({
    this.cpu,
    this.appMemory = 0,
    this.memTotal,
    this.memFree,
    this.diskFree,
    this.diskTotal,
    this.netRate = 0,
    this.netTotal = 0,
  });

  /// Machine CPU load (0–100), including MAFFT, IQ-TREE and the others.
  final double? cpu;

  /// Resident memory of the app (bytes).
  final int appMemory;

  /// Physical memory of the machine (bytes).
  final int? memTotal, memFree;

  /// Space on the drive of the project folder (bytes).
  final int? diskFree, diskTotal;

  /// App network traffic: bytes per second over the last interval, and
  /// total since the app started.
  final double netRate;
  final int netTotal;
}

/// Reads CPU, memory, disk and network every few seconds. On Windows one
/// PowerShell process stays open and prints a line per reading (starting a
/// process per reading would cost more than it measures); on Linux the
/// /proc files and `df` are read; on macOS `df` and the app's own numbers.
class ResourceMonitor extends ChangeNotifier {
  ResourceMonitor({required this.folder, this.interval = const Duration(seconds: 4)});

  /// Folder whose drive is measured (the project folder).
  final String Function() folder;
  final Duration interval;

  ResourceSnapshot value = const ResourceSnapshot();
  Process? _ps;
  String? _psDrive;
  Timer? _timer;
  int _lastNet = 0;
  DateTime _lastAt = DateTime.now();
  (int, int)? _lastCpu;
  bool _disposed = false;

  void start() {
    _timer = Timer.periodic(interval, (_) => _tick());
    _tick();
  }

  String _drive() {
    final f = folder();
    final m = RegExp(r'^([A-Za-z]):').firstMatch(f);
    return (m?.group(1) ?? 'C').toUpperCase();
  }

  Future<void> _tick() async {
    if (_disposed) return;
    final now = DateTime.now();
    final net = NetMeter.received + NetMeter.sent;
    final secs = now.difference(_lastAt).inMilliseconds / 1000;
    final rate = secs <= 0 ? 0.0 : (net - _lastNet) / secs;
    _lastNet = net;
    _lastAt = now;
    double? cpu = value.cpu;
    int? memTotal = value.memTotal, memFree = value.memFree, diskFree = value.diskFree, diskTotal = value.diskTotal;
    try {
      if (Platform.isWindows) {
        _ensurePowerShell();
      } else if (Platform.isLinux) {
        cpu = _linuxCpu();
        final mem = File('/proc/meminfo').readAsLinesSync();
        int kb(String key) => int.parse(RegExp(r'(\d+)').firstMatch(mem.firstWhere((l) => l.startsWith(key)))!.group(1)!) * 1024;
        memTotal = kb('MemTotal');
        memFree = kb('MemAvailable');
      }
      if (!Platform.isWindows) {
        final df = await Process.run('df', ['-kP', folder()]);
        final parts = (df.stdout as String).trim().split('\n').last.split(RegExp(r'\s+'));
        if (parts.length >= 4) {
          diskTotal = int.parse(parts[1]) * 1024;
          diskFree = int.parse(parts[3]) * 1024;
        }
      }
    } catch (_) {
      // A reading that fails keeps the previous one.
    }
    if (_disposed) return;
    value = ResourceSnapshot(
      cpu: cpu,
      appMemory: ProcessInfo.currentRss,
      memTotal: memTotal,
      memFree: memFree,
      diskFree: diskFree,
      diskTotal: diskTotal,
      netRate: rate,
      netTotal: net,
    );
    notifyListeners();
  }

  double? _linuxCpu() {
    final f = File('/proc/stat').readAsLinesSync().first.split(RegExp(r'\s+')).skip(1).map(int.parse).toList();
    final idle = f[3] + (f.length > 4 ? f[4] : 0);
    final total = f.fold(0, (a, b) => a + b);
    final last = _lastCpu;
    _lastCpu = (idle, total);
    if (last == null || total == last.$2) return null;
    return 100 * (1 - (idle - last.$1) / (total - last.$2));
  }

  /// Starts (or, when the project moved to another drive, restarts) the
  /// PowerShell loop.
  void _ensurePowerShell() {
    final drive = _drive();
    if (_ps != null && _psDrive == drive) return;
    _ps?.kill();
    _psDrive = drive;
    final seconds = interval.inSeconds.clamp(2, 60);
    final script = r"$ErrorActionPreference='SilentlyContinue'; "
        "\$d = New-Object System.IO.DriveInfo('$drive'); "
        // Ends on its own when SPLACE is gone (closed or killed).
        r'while ($true) { '
        "if (-not (Get-Process -Id $pid -ErrorAction SilentlyContinue)) { exit }; "
        r'$c = (Get-CimInstance Win32_Processor | Measure-Object -Property LoadPercentage -Average).Average; '
        r'$o = Get-CimInstance Win32_OperatingSystem; '
        r'Write-Output ("{0}|{1}|{2}|{3}|{4}" -f $c, $o.TotalVisibleMemorySize, $o.FreePhysicalMemory, $d.AvailableFreeSpace, $d.TotalSize); '
        "Start-Sleep -Seconds $seconds }";
    Process.start('powershell', ['-NoProfile', '-NonInteractive', '-Command', script], runInShell: false).then((p) {
      if (_disposed) {
        p.kill();
        return;
      }
      _ps = p;
      p.stdout.transform(utf8.decoder).transform(const LineSplitter()).listen(_onLine);
      p.stderr.drain<void>();
    }).catchError((_) {});
  }

  void _onLine(String line) {
    final f = line.trim().split('|');
    if (f.length < 5 || _disposed) return;
    final cpu = double.tryParse(f[0].replaceAll(',', '.'));
    int? n(String s) => int.tryParse(s.trim());
    final memTotal = n(f[1]), memFree = n(f[2]);
    value = ResourceSnapshot(
      cpu: cpu ?? value.cpu,
      appMemory: ProcessInfo.currentRss,
      memTotal: memTotal == null ? value.memTotal : memTotal * 1024,
      memFree: memFree == null ? value.memFree : memFree * 1024,
      diskFree: n(f[3]) ?? value.diskFree,
      diskTotal: n(f[4]) ?? value.diskTotal,
      netRate: value.netRate,
      netTotal: value.netTotal,
    );
    notifyListeners();
  }

  @override
  void dispose() {
    _disposed = true;
    _timer?.cancel();
    _ps?.kill();
    super.dispose();
  }
}

/// "1.2 GB", "340 MB", "12 KB".
String formatBytes(num b) {
  const units = ['B', 'KB', 'MB', 'GB', 'TB'];
  var v = b.toDouble();
  var i = 0;
  while (v >= 1024 && i < units.length - 1) {
    v /= 1024;
    i++;
  }
  return '${v.toStringAsFixed(v >= 100 || i == 0 ? 0 : 1)} ${units[i]}';
}
