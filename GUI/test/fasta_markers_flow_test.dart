@Tags(['screenshots'])
library;

import 'dart:io';
import 'dart:ui' as ui;

import 'package:flutter/material.dart';
import 'package:flutter/rendering.dart';
import 'package:flutter/services.dart';
import 'package:flutter_test/flutter_test.dart';
import 'package:path/path.dart' as p;
import 'package:splace/src/app.dart';
import 'package:splace/src/core/saturation.dart';
import 'package:splace/src/core/tool_commands.dart';
import 'package:splace/src/state/project_controller.dart';

/// Many markers from per-sample FASTA files (target capture: `>ITV01908-L12`):
/// the "Sample and marker" dialog, the FASTA per marker, then every step
/// with the real tools, timing how long each screen takes to draw. Opt-in:
///
///   SPLACE_MARKERS=teste/amanda SPLACE_SCREENSHOTS=out_dir flutter test test/fasta_markers_flow_test.dart
///
/// SPLACE_MARKERS_QUICK=1 stops after step 4 (no external tools).
Future<void> _loadFonts() async {
  final families = {
    'Montserrat': ['Regular', 'Italic', 'Medium', 'SemiBold', 'Bold', 'ExtraBold'].map((w) => 'assets/fonts/Montserrat-$w.ttf'),
    'JetBrainsMono': ['assets/fonts/JetBrainsMono-Regular.ttf', 'assets/fonts/JetBrainsMono-Medium.ttf'],
    'FASolid': ['assets/fonts/fa-solid-900.ttf'],
    'FABrands': ['assets/fonts/fa-brands-400.ttf'],
  };
  for (final e in families.entries) {
    final loader = FontLoader(e.key);
    for (final f in e.value) {
      loader.addFont(rootBundle.load(f));
    }
    await loader.load();
  }
}

void main() {
  final out = Platform.environment['SPLACE_SCREENSHOTS'];
  final data = Platform.environment['SPLACE_MARKERS'];
  // Only the dialog and the markers (no tools).
  final quick = Platform.environment['SPLACE_MARKERS_QUICK'] == '1';

  testWidgets('many markers from per-sample FASTA files', (tester) async {
    await _loadFonts();
    tester.view.physicalSize = const Size(1280, 720);
    tester.view.devicePixelRatio = 1;
    addTearDown(tester.view.reset);
    final work = Directory.systemTemp.createTempSync('splace_markers_');
    final key = GlobalKey();
    final c = ProjectController(autoDetectTools: false)..workDir = p.join(work.path, 'project');
    addTearDown(c.dispose);
    await tester.pumpWidget(RepaintBoundary(key: key, child: SplaceApp(controller: c)));
    await tester.pump();
    c.setSplashVisible(false);
    await tester.pump();

    final timings = <String, int>{};
    Future<void> shot(String name) async {
      await tester.pump(const Duration(milliseconds: 400));
      if (out == null) return;
      await tester.runAsync(() async {
        final image = await (key.currentContext!.findRenderObject()! as RenderRepaintBoundary).toImage();
        final bytes = await image.toByteData(format: ui.ImageByteFormat.png);
        File('$out/$name.png').writeAsBytesSync(bytes!.buffer.asUint8List());
      });
    }

    /// Opens a step and records how long the frames took (wall clock).
    Future<void> goTo(String step, String name) async {
      final sw = Stopwatch()..start();
      await tester.tap(find.text(step).first);
      await tester.pumpAndSettle();
      timings[name] = sw.elapsedMilliseconds;
      await shot(name);
    }

    final files = Directory(data!).listSync().whereType<File>().where((f) => f.path.endsWith('.fasta')).map((f) => f.path).toList()..sort();
    await tester.runAsync(() async {
      await c.detectTools();
      await c.addFiles(files);
      await Future<void>.delayed(const Duration(milliseconds: 600));
    });
    // The dialog opens by itself after FASTA files are read.
    await tester.pump();
    await tester.pump(const Duration(milliseconds: 300));
    await tester.runAsync(() => Future<void>.delayed(const Duration(milliseconds: 800)));
    await tester.pump(const Duration(milliseconds: 300));
    expect(find.byType(Dialog), findsOneWidget);
    await shot('m1_dialog');
    await tester.dragFrom(const Offset(640, 450), const Offset(0, -480));
    await tester.pump(const Duration(milliseconds: 300));
    await shot('m2_dialog_samples');
    // ITV01054.fasta has its headers written ITV1054: offered as a rename.
    expect(find.textContaining('ITV01054.fasta  →  ITV1054'), findsOneWidget);
    // …and as a sample name with another structure (4 digits, not 5).
    expect(find.text('Standardize to ITV#####'), findsOneWidget);
    await tester.tap(find.text('Standardize to ITV#####'));
    await tester.pump(const Duration(milliseconds: 400));
    await tester.runAsync(() => Future<void>.delayed(const Duration(milliseconds: 300)));
    await tester.pump(const Duration(milliseconds: 300));
    expect(find.text('Standardize to ITV#####'), findsNothing, reason: 'every name now has 5 digits');
    await shot('m2b_dialog_standardized');

    final apply = find.textContaining('marker FASTA files');
    expect(apply, findsOneWidget);
    await tester.tap(apply);
    for (var i = 0; i < 100 && c.lociExportDir == null; i++) {
      await tester.pump();
      await tester.runAsync(() => Future<void>.delayed(const Duration(milliseconds: 100)));
    }
    await tester.pump();
    final org = c.organized!;
    // ignore: avoid_print
    print('organized: ${org.samples.length} samples × ${org.loci.length} markers');
    expect(org.loci.length, greaterThan(400));
    expect(org.samples.length, 20);
    expect(org.samples, contains('ITV01054'));
    expect(c.lociMode, isTrue);
    expect(c.codonAwareTrimming, isFalse, reason: 'loci: the reading frame is not kept by default');
    final written = Directory(p.join(c.workDir, '01_loci')).listSync().where((f) => f.path.endsWith('.fasta')).length;
    expect(written, org.loci.length, reason: 'one FASTA per marker');

    await goTo('Sample and Locus Fields', 'm3_fields');
    await goTo('Select Markers', 'm4_markers');
    if (quick) return;

    final sw = Stopwatch()..start();
    await tester.runAsync(() async {
      c.setThreads(Platform.numberOfProcessors);
      await c.prepareLoci();
      await c.alignAll();
    });
    // ignore: avoid_print
    print('align ${org.loci.length} markers: ${sw.elapsedMilliseconds} ms');
    expect(c.loci.values.where((l) => l.coding), isEmpty, reason: 'unannotated loci are not coding');
    c.closeModal();
    await tester.pump();
    await goTo('Sequence Alignment', 'm5_alignment');
    sw.reset();
    await tester.runAsync(() async {
      await c.trimAll();
      // Optional for loci; without coding markers, all sites are tested.
      await c.runSaturation();
      await c.buildMatrix();
    });
    // ignore: avoid_print
    print('trim + saturation + matrix: ${sw.elapsedMilliseconds} ms');
    c.closeModal();
    await tester.pump();
    expect(c.pooledSaturation, isNull, reason: 'non-coding loci have no codon positions');
    expect(c.loci.values.where((l) => l.saturation != null).length, org.loci.length);
    expect(c.loci.values.first.saturation!.subsets.keys, [CodonSubset.all]);
    expect(c.overview!.loci.length, greaterThan(400));
    await goTo('Saturation Test', 'm6_saturation');
    await goTo('Concatenate and NEXUS', 'm7_concat');
    await tester.ensureVisible(find.text('TAXA × MARKERS'));
    await tester.pumpAndSettle();
    await shot('m7b_matrix_table');
    sw.reset();
    await tester.runAsync(() async {
      await c.runIqtree(IqtreeOptions(model: 'GTR+F+G4', ufboot: 1000, alrt: 0, threads: '${Platform.numberOfProcessors}', seed: 1));
    });
    // ignore: avoid_print
    print('IQ-TREE: ${sw.elapsedMilliseconds} ms, status ${c.treeRuns.first.status}');
    c.closeModal();
    await tester.pump();
    await goTo('Phylogenetic Inference', 'm8_inference');
    // ignore: avoid_print
    print('frame times (ms): $timings');
  }, skip: data == null, timeout: const Timeout(Duration(minutes: 40)));
}
