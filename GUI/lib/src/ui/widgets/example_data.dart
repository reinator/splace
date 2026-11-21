import 'dart:convert';
import 'dart:io';

import 'package:flutter/services.dart';
import 'package:path/path.dart' as p;

import '../../state/project_controller.dart';

/// Example data bundled with the app (`assets/examples/<kind>/`): files the
/// maintainers drop in those folders, loaded by the buttons of "Before you
/// start" (step 1).
enum ExampleKind {
  genbank('genbank', ['.gb', '.gbk', '.genbank', '.gbff']),
  loci('loci', ['.fasta', '.fa', '.fas', '.fna', '.ffn']),
  barcoding('barcoding', ['.fasta', '.fa', '.fas', '.fna']);

  const ExampleKind(this.folder, this.extensions);
  final String folder;
  final List<String> extensions;

  /// Folder of the files in the repository.
  String get source => 'GUI/assets/examples/$folder';
}

/// Example files of each kind in the app's assets (empty folders give none).
Future<Map<ExampleKind, List<String>>> listExampleAssets() async {
  try {
    final manifest = await AssetManifest.loadFromAssetBundle(rootBundle);
    final all = manifest.listAssets();
    return {
      for (final k in ExampleKind.values)
        k: [
          for (final a in all)
            if (a.startsWith('assets/examples/${k.folder}/') && k.extensions.contains(p.extension(a).toLowerCase())) a,
        ]..sort(),
    };
  } catch (_) {
    return {for (final k in ExampleKind.values) k: const []};
  }
}

/// Copies the example files of [kind] to `Documents/SPLACE/examples/<kind>`
/// and loads them: GenBank records, FASTA loci (the Sample and marker dialog
/// opens) or FASTA barcodes (with the barcoding switch on).
Future<int> loadExample(ProjectController c, ExampleKind kind, List<String> assets) async {
  if (assets.isEmpty) return 0;
  final dir = Directory(p.join(c.documentsDir, 'examples', kind.folder))..createSync(recursive: true);
  final paths = <String>[];
  for (final a in assets) {
    final data = await rootBundle.load(a);
    final file = File(p.join(dir.path, p.basename(a)));
    await file.writeAsBytes(data.buffer.asUint8List(data.offsetInBytes, data.lengthInBytes));
    paths.add(file.path);
  }
  c.log('Example data (${kind.folder}): ${paths.length} file(s) from ${kind.source}', source: 'SPLACE');
  switch (kind) {
    case ExampleKind.genbank:
      await c.addGenBankFiles(paths);
    case ExampleKind.loci:
      await c.addFiles(paths);
    case ExampleKind.barcoding:
      await c.addFiles(paths);
      c.setBarcodingMode(true);
  }
  return paths.length;
}

/// Tooltips of the example buttons, edited in `assets/examples/tooltips.json`
/// (HTML: bold, italic, line breaks): key (`1b`, `1c`, `1d_loci`,
/// `1d_barcoding`) → language → text.
Future<Map<String, Map<String, String>>> loadExampleTooltips() async {
  try {
    final json = jsonDecode(await rootBundle.loadString('assets/examples/tooltips.json')) as Map<String, dynamic>;
    return {
      for (final e in json.entries)
        if (e.value is Map) e.key: {for (final l in (e.value as Map).entries) '${l.key}': '${l.value}'},
    };
  } catch (_) {
    return const {};
  }
}

/// The tooltip of [key] in [lang] (English when that language is missing).
String? exampleTooltip(Map<String, Map<String, String>> tips, String key, String lang) {
  final t = tips[key];
  final v = t?[lang] ?? t?['en'];
  return v == null || v.trim().isEmpty ? null : v;
}
