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
import 'package:splace/src/core/genbank_markers.dart';
import 'package:splace/src/core/marker_search.dart';
import 'package:splace/src/core/supermatrix.dart';
import 'package:splace/src/state/project_controller.dart';
import 'package:splace/src/ui/widgets/matrix_table.dart';

/// A marker search with hundreds of records (Caranx COI: ~800 sequences,
/// many of them partial), through every step with the real tools, timing
/// how long each screen takes to draw. Opt-in, needs the network:
///
///   SPLACE_NETWORK=1 SPLACE_SCREENSHOTS=out_dir flutter test test/marker_search_flow_test.dart
///
/// SPLACE_MARKER_TAXON (default Caranx) and SPLACE_MARKER_GENES (COI).
/// SPLACE_MARKER_ALL=1 keeps every record as its own sample (one per
/// species off), as when the accession stays in the header.
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
  final env = Platform.environment;
  final enabled = env['SPLACE_NETWORK'] == '1';
  final out = env['SPLACE_SCREENSHOTS'];
  final taxon = env['SPLACE_MARKER_TAXON'] ?? 'Caranx';
  final genes = (env['SPLACE_MARKER_GENES'] ?? 'COI').split(',');
  final all = env['SPLACE_MARKER_ALL'] == '1';

  testWidgets('marker search with hundreds of records', (tester) async {
    await _loadFonts();
    // The test binding answers every HTTP request with 400: use the network.
    HttpOverrides.global = null;
    tester.view.physicalSize = const Size(1280, 720);
    tester.view.devicePixelRatio = 1;
    addTearDown(tester.view.reset);
    final work = Directory.systemTemp.createTempSync('splace_msearch_');
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

    // A lettered part from the step navigator ("5d"…).
    Future<void> goToSub(String id, String name) async {
      await tester.ensureVisible(find.text(id).first);
      await tester.tap(find.text(id).first);
      await _settle(tester);
      await shot(name);
    }

    Future<void> goTo(String step, String name) async {
      final sw = Stopwatch()..start();
      await tester.ensureVisible(find.text(step).first);
      await tester.tap(find.text(step).first);
      await _settle(tester);
      timings[name] = sw.elapsedMilliseconds;
      // ignore: avoid_print
      print('$name: ${timings[name]} ms');
      await shot(name);
    }

    final search = MarkerSearch(taxon: taxon, genes: genes);
    final sw = Stopwatch()..start();
    await tester.runAsync(() async {
      await c.detectTools();
      final found = await c.searchMarkers(search);
      expect(found, isNotEmpty);
      final selected = found!.take(900).toList();
      // ignore: avoid_print
      print('found ${found.length}, importing ${selected.length}');
      await c.importMarkers(selected, search);
    });
    // ignore: avoid_print
    print('search + import: ${sw.elapsedMilliseconds} ms, ${c.records.length} records');
    c.closeModal();
    c.skipTaxonomy();
    if (all) c.setOnePerTaxon(false);
    await tester.pump();
    await goTo('Loaded Records', 's2_records');
    await goTo('Features and Markers', 's4_markers');
    if (all) {
      // family_species repeats inside COI: the dialog cannot confirm it and
      // the alignment stays locked until the accession is added.
      c.setHeader(const HeaderTemplate(fields: ['family', 'species']));
      expect(c.headerBlocked, isTrue);
      await goTo('Sequence Alignment', 's5a_blocked');
      await tester.tap(find.text('Open Header Builder'));
      await tester.pump(const Duration(milliseconds: 500));
      await shot('s5a_dialog_blocked');
      final confirm = find.widgetWithText(FilledButton, 'Confirm header');
      expect(tester.widget<FilledButton>(confirm).onPressed, isNull);
      final add = find.descendant(of: find.byType(Dialog), matching: find.text('Add Accession'));
      await tester.ensureVisible(add);
      await tester.tap(add);
      await tester.pump(const Duration(milliseconds: 300));
      await shot('s5a_dialog_fixed');
      await tester.tap(confirm);
      await tester.pump(const Duration(milliseconds: 500));
      expect(c.headerBlocked, isFalse);
    }
    // Records named only to genus (Caranx sp.) cannot be joined as one
    // species: the alignment waits until they are left out.
    if (c.headerBlocked) {
      await goTo('Sequence Alignment', 's5a_undetermined');
      expect(c.undeterminedConflicts, isNotEmpty);
      final drop = find.textContaining('Leave out the');
      await tester.ensureVisible(drop.first);
      await tester.tap(drop.first);
      await tester.pump(const Duration(milliseconds: 300));
      expect(c.headerBlocked, isFalse);
    }
    c.confirmHeader();
    sw.reset();
    await tester.runAsync(() async {
      c.setThreads(Platform.numberOfProcessors);
      await c.prepareLoci();
      await c.alignAll();
    });
    expect(c.loci['COI']!.mostlyPartial, isTrue, reason: 'barcodes are partial CDS');
    expect(c.saturationRequired, isFalse);
    // ignore: avoid_print
    print('align: ${sw.elapsedMilliseconds} ms; ${c.loci.values.map((l) => '${l.name} ${l.aligned?.length} seq coding=${l.coding}').join(', ')}');
    c.closeModal();
    await tester.pump();
    await goTo('Sequence Alignment', 's5_alignment');
    await goToSub('4b', 's5b_seqtype');
    await goToSub('4c', 's5d_trim_prompt');
    c.openTrim();
    await tester.pump();
    await goToSub('4d', 's5e_trim');
    sw.reset();
    await tester.runAsync(() async {
      await c.trimAll();
    });
    c.closeModal();
    await tester.pump();
    await goTo('Sequence Alignment', 's5b_trimmed');
    await tester.runAsync(() async {
      await c.runSaturation();
    });
    c.closeModal();
    await tester.pump();
    // ignore: avoid_print
    print('trim + saturation: ${sw.elapsedMilliseconds} ms');
    await goTo('Saturation Test', 's6_saturation');
    expect(c.pooledSaturation?.reports.map((r) => r.name), ['COI'], reason: 'one gene, tested on its own');
    await goToSub('5c', 's6c_suggestion');
    c.applySaturationRecommendation(perLocus: false);
    await tester.pump();
    sw.reset();
    await tester.runAsync(() async {
      await c.buildMatrix();
    });
    c.closeModal();
    await tester.pump();
    // ignore: avoid_print
    print('matrix: ${sw.elapsedMilliseconds} ms, ${c.matrix!.taxa.length} taxa');
    // One gene: "by locus" would be one partition, so only by codon
    // position (frame known) or none are offered.
    final coi = c.loci['COI']!;
    // ignore: avoid_print
    print('COI coding=${coi.coding} frame=${coi.frame} trimmed by codons=${c.codonTrimmed(coi)} schemes=${c.partitionSchemes} used=${c.effectiveScheme}');
    expect(c.partitionSchemes, c.codonPositionsUsable(coi) ? [PartitionScheme.byCodonPosition, PartitionScheme.none] : [PartitionScheme.none]);
    expect(c.partitionSchemes, contains(c.effectiveScheme));
    await goTo('Concatenate and NEXUS', 's7_concat');
    await goToSub('6b', 's7b_partitions');
    await tester.ensureVisible(find.text('TAXA × MARKERS'));
    await _settle(tester);
    await shot('s7b_matrix_table');
    // The table scrolls inside itself: the page does not grow with the taxa.
    expect(tester.getSize(find.byType(MatrixOverviewTable)).height, lessThan(560));
    await goTo('Phylogenetic Inference', 's8_inference');
    // A rebuild (every controller change) and a scroll frame.
    for (var i = 0; i < 3; i++) {
      sw.reset();
      c.notifyListeners();
      await tester.pump();
      final rebuild = sw.elapsedMilliseconds;
      sw.reset();
      await tester.drag(find.byType(SingleChildScrollView).first, Offset(0, i.isEven ? -300 : 300));
      await tester.pump();
      // ignore: avoid_print
      print('rebuild: $rebuild ms, scroll frame: ${sw.elapsedMilliseconds} ms');
    }
    // ignore: avoid_print
    print('frame times (ms): $timings');
  }, skip: !enabled, timeout: const Timeout(Duration(minutes: 40)));
}

/// Frames for 2 s: some screens keep a spinner running, so waiting for
/// the animations to end would never finish.
Future<void> _settle(WidgetTester tester) async {
  for (var i = 0; i < 20; i++) {
    await tester.pump(const Duration(milliseconds: 100));
  }
}
