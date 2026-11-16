@Tags(['screenshots'])
library;

import 'dart:io';
import 'dart:ui' as ui;

import 'package:flutter/material.dart';
import 'package:flutter/rendering.dart';
import 'package:flutter/services.dart';
import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/app.dart';
import 'package:flutter/gestures.dart';
import 'package:splace/src/core/genbank_markers.dart';
import 'package:splace/src/core/taxonomy_sources.dart';
import 'package:splace/src/services/ncbi.dart';
import 'package:splace/src/ui/chrome/genome_review_dialog.dart';
import 'package:splace/src/ui/fa.dart';
import 'package:splace/src/ui/widgets/splace_tooltip.dart';
import 'package:splace/src/core/tool_commands.dart';
import 'package:splace/src/core/tree_progress.dart';
import 'package:splace/src/state/project_controller.dart';

/// Renders the real UI (bundled fonts, logo, flags) to PNG files so the
/// layout can be compared with the HTML app. Opt-in:
///
///   SPLACE_SCREENSHOTS=out_dir [SPLACE_NETWORK=1 | SPLACE_DATASET=datasets/inputs] flutter test test/screenshots_test.dart
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
  final dataset = Platform.environment['SPLACE_DATASET'];
  final network = Platform.environment['SPLACE_NETWORK'] == '1';

  testWidgets('screenshots', (tester) async {
    // flutter_test replaces HttpClient with a stub that answers 400; the
    // GenBank flow needs the real NCBI/GBIF.
    if (network) HttpOverrides.global = null;
    await _loadFonts();
    // The smallest window the app allows (1280×720).
    tester.view.physicalSize = const Size(1280, 720);
    tester.view.devicePixelRatio = 1;
    addTearDown(tester.view.reset);

    final key = GlobalKey();
    final c = ProjectController(autoDetectTools: false);
    addTearDown(c.dispose);
    await tester.pumpWidget(RepaintBoundary(key: key, child: SplaceApp(controller: c)));
    await tester.runAsync(() async {
      for (final a in ['assets/images/splace_logo.png', 'assets/images/splace_logo_card.png', 'assets/images/flag_us.png', 'assets/images/flag_es.png', 'assets/images/flag_br.png']) {
        await precacheImage(AssetImage(a), tester.element(find.byType(Scaffold)));
      }
    });
    await tester.pump();

    Future<void> shot(String name) async {
      await tester.pump(const Duration(milliseconds: 400));
      await tester.runAsync(() async {
        final boundary = key.currentContext!.findRenderObject()! as RenderRepaintBoundary;
        final image = await boundary.toImage();
        final bytes = await image.toByteData(format: ui.ImageByteFormat.png);
        File('$out/$name.png').writeAsBytesSync(bytes!.buffer.asUint8List());
      });
    }

    // Some screens keep an animation running (rainbow modal headers,
    // spinners), so a settle that never ends is not an error here.
    Future<void> settle() async {
      try {
        await tester.pumpAndSettle(const Duration(milliseconds: 100), EnginePhase.sendSemanticsUpdate, const Duration(seconds: 5));
      } on FlutterError {
        await tester.pump(const Duration(milliseconds: 300));
      }
    }

    await shot('01_splash');
    c.setSplashVisible(false);
    await tester.pump();
    await shot('02_load');
    // Option 1 for individual markers (SynGenes names), then back.
    await tester.tap(find.text('Individual markers'));
    await tester.pump();
    await tester.enterText(find.byType(TextField).first, 'Astyanax');
    await tester.pump();
    await tester.drag(find.byType(Scrollable).first, const Offset(0, -220));
    await tester.pump(const Duration(milliseconds: 300));
    await shot('02f_marker_search');
    expect(find.textContaining('COX1[Gene]'), findsOneWidget);
    await tester.enterText(find.byType(TextField).first, '');
    await tester.tap(find.text('Genomes'));
    await tester.pump();
    await tester.drag(find.byType(Scrollable).first, const Offset(0, 220));
    await tester.pump(const Duration(milliseconds: 300));

    // NCBI API key modal (key button in the top bar).
    // Modal headers keep the rainbow line moving, so open them with a fixed
    // pump instead of pumpAndSettle (which would never settle).
    await tester.tap(find.byIcon(FA.key).first);
    await tester.pump();
    await tester.pump(const Duration(milliseconds: 300));
    await shot('02a_api_key');
    await tester.tap(find.text('Close').last);
    await settle();

    // Genome type select, open.
    await tester.tap(find.text('Mitochondrion').first);
    await settle();
    await shot('02e_select');
    await tester.tapAt(const Offset(640, 700));
    await settle();

    Future<void> goTo(String step, String name) async {
      if (find.text(step).evaluate().isEmpty) {
        // ignore: avoid_print
        print('skip $name: no "$step"');
        return;
      }
      await tester.tap(find.text(step).first);
      await settle();
      await shot(name);
    }

    // Scrolls to [f] and takes the shot; skipped when the screen no longer
    // shows that text.
    Future<void> view(Finder f, String name) async {
      if (f.evaluate().isEmpty) {
        // ignore: avoid_print
        print('skip $name: $f');
        return;
      }
      await tester.ensureVisible(f.first);
      await settle();
      await shot(name);
    }

    if (network) {
      // Arrow tooltip on an info icon of step 1.
      final mouse = await tester.createGesture(kind: PointerDeviceKind.mouse);
      await mouse.addPointer(location: Offset.zero);
      await mouse.moveTo(tester.getCenter(find.byIcon(FA.circleInfo).first));
      await tester.pump(const Duration(milliseconds: 400));
      await shot('02b_tooltip');
      await mouse.moveTo(Offset.zero);
      await tester.pump();

      // "Review genomes before importing" with a live search.
      const search = GenomeSearch(taxon: 'Bufonidae', refseqOnly: false);
      List<GenomeCandidate> candidates = const [];
      await tester.runAsync(() async {
        final ids = await c.ncbi.searchGenomes(search, retmax: 40);
        candidates = await c.ncbi.summaries(ids);
      });
      showGenomeReviewDialog(tester.element(find.byType(Scaffold)), search, candidates);
      await tester.pump();
      await tester.pump(const Duration(milliseconds: 300));
      await shot('02c_review');
      Navigator.of(tester.element(find.byType(Scaffold))).pop();
      await settle();

      // Progress modal (`#progressModal`).
      c.modal = TaskModal(ModalKind.ncbiFetch, '5')
        ..done = 2
        ..total = 5
        ..current = 'NC_047224.1 Anaxyrus americanus';
      c.setThemeMode(c.themeMode);
      await tester.pump();
      await shot('02d_progress_modal');
      c.closeModal();
      await tester.pump();

      // The HTML app's GenBank workflow with live NCBI/GBIF data.
      await tester.runAsync(() async {
        await c.detectTools();
        await c.fetchAccessions('2216740843 2301278816 1915405715 1834962159 1834962145');
        await c.fetchTaxonomy();
      });
      c.closeModal();
      await tester.pump();
      await shot('03_loaded');
      // External tools dialog (wrench), display settings and activity log.
      await tester.tap(find.byIcon(FA.screwdriverWrench).first);
      await tester.pump();
      await tester.pump(const Duration(milliseconds: 300));
      await shot('03a_tools');
      await tester.tap(find.text('Done').last);
      await settle();
      await tester.tap(find.byIcon(FA.sliders).first);
      await settle();
      await shot('03b_settings');
      await tester.tap(find.byIcon(FA.sliders).first);
      await settle();
      await tester.tap(find.byIcon(FA.rectangleList).first);
      await settle();
      await shot('03c_log');
      await tester.tap(find.byIcon(FA.rectangleList).first);
      await settle();
      await goTo('Loaded Records', '04_records');
      c.setTaxonomySource(TaxonomySource.salve);
      await tester.pump();
      await shot('04b_taxonomy_salve');
      c.setTaxonomySource(TaxonomySource.ncbi);
      await tester.runAsync(c.fetchTaxonomy);
      c.closeModal();
      await tester.pump();
      await shot('04c_taxonomy_ncbi');
      final table = find.text('Records table');
      if (table.evaluate().isNotEmpty) {
        await tester.ensureVisible(table.last);
        await settle();
        await shot('04d_records_table');
      }
      c.toggleFeatureType('CDS', true);
      c.toggleFeatureType('rRNA', true);
      await tester.pump();
      await goTo('Feature Types', '05_feature_types');
      c.selectDefaultMarkers();
      await tester.pump();
      await goTo('Select Markers', '06_markers');
      await tester.ensureVisible(find.text('Gene Presence Heatmap').last);
      await settle();
      // Both heatmaps open: they fill the width of the page.
      await tester.tap(find.text('Show heatmap').first);
      await settle();
      await view(find.text('Gene Presence Heatmap').last, '06b_heatmap');
      await tester.ensureVisible(find.text('Marker Length Heatmap').last);
      await settle();
      await tester.tap(find.text('Show heatmap').first);
      await settle();
      await view(find.text('Marker Length Heatmap').last, '06c_length_heatmap');
      for (final m in c.activeMarkers.where((m) => !{'COI', 'CYTB', 'ND2', 'ATP6'}.contains(m)).toList()) {
        c.toggleMarker(m, false);
      }
      c.setHeader(const HeaderTemplate(fields: ['species', 'accession']));
      c.confirmHeader();
      c.setSeqMode(SeqMode.codon);
      await tester.pump();
      await goTo('Sequence Alignment', '07_alignment');
      await goTo('Sequence Type', '07b_seqtype');
      await tester.drag(find.byType(Scrollable).first, const Offset(0, -260));
      await settle();
      await shot('07e_seqtype_examples');
      await goTo('MAFFT Settings', '07c_mafft');
      await tester.runAsync(() async {
        await c.alignAll();
      });
      await tester.pump();
      await shot('07d_alignment_modal');
      c.closeModal();
      await tester.pump();
      await goTo('Run Alignment', '08_trim_prompt');
      if (find.text('Yes, set up trimming').evaluate().isNotEmpty) {
        await tester.tap(find.text('Yes, set up trimming').first);
        await settle();
      }
      if (find.textContaining('Columns removed').evaluate().isNotEmpty) {
        await tester.ensureVisible(find.textContaining('Columns removed').first);
        await settle();
        await view(find.textContaining('ClipKIT mode'), '08e_trim_example');
        // The animation half-way: the columns to remove are marked.
        await tester.tap(find.byTooltip('Play again').first);
        await tester.pump();
        await tester.pump(const Duration(milliseconds: 1700));
        await shot('08f_trim_marked');
        await settle();
      }
      await tester.runAsync(() async {
        await c.trimAll();
        await c.runSaturation();
        c.applySaturationRecommendation(perLocus: false);
        c.setConcat(c.concat.copyWith(outgroups: [(c.loci.values.first.trimmed ?? c.loci.values.first.aligned)!.first.header]));
        await c.buildMatrix();
        await c.runIqtree(const IqtreeOptions(model: 'GTR+F+G4', ufboot: 1000, alrt: 0, threads: '2', seed: 1), perGene: true);
        // RAxML-NG 2: MOOSE models, fast search, SH-aLRT and EBG.
        await c.runRaxml(const RaxmlOptions(maxThreads: 2, fast: true, supportMetrics: 'sh,ebg', seed: 1));
      });
      c.closeModal();
      await tester.pump();
      await tester.ensureVisible(find.text('Alignment Results'));
      await settle();
      await tester.tap(find.byWidgetPredicate((w) => w is SplaceTooltip && w.message == 'Show alignment').first);
      await settle();
      await shot('08b_results_viewer');
      final trimmedTab = find.text('Trimmed');
      if (trimmedTab.evaluate().length > 1) {
        await tester.ensureVisible(trimmedTab.last);
        await settle();
        await tester.tap(trimmedTab.last);
        await settle();
        await shot('08c_trimmed');
      }
      final translate = find.text('Translate to amino acids');
      if (translate.evaluate().isNotEmpty) {
        await tester.ensureVisible(translate.first);
        await settle();
        await tester.tap(translate.first);
        await settle();
        await tester.tap(find.text('Translation table').first);
        await settle();
        await view(find.textContaining('Genetic code '), '08d_translated');
      }
      await goTo('Saturation Test', '09_saturation');
      await view(find.textContaining('The shade inside each colour'), '09b_saturation_heatmap');
      await view(find.text('What to do with the codon positions'), '09e_saturation_recommendation');
      await view(find.textContaining('of each codon position'), '09c_saturation_plots');
      await view(find.text('HOW TO READ THE PLOTS'), '09f_saturation_reading');
      await view(find.text('COMPOSITION AND VARIABILITY'), '09d_saturation_composition');
      await goTo('Concatenate and NEXUS', '11_concat');
      await view(find.text('How to split the matrix into partitions'), '11b_concat_partitions');
      await view(find.textContaining('codon positions are kept').evaluate().isNotEmpty ? find.textContaining('codon positions are kept') : find.textContaining('Kept:'), '11d_codon_example');
      await view(find.text('WHAT WILL BE BUILT'), '11c_concat_summary');
      await goTo('Phylogenetic Inference', '12_inference');
      await tester.ensureVisible(find.text('Run RAxML-NG').first);
      await settle();
      await tester.ensureVisible(find.textContaining('MOOSE (RAxML-NG').first);
      await settle();
      await tester.drag(find.byType(Scrollable).first, const Offset(0, -140));
      await settle();
      await shot('12b_raxml');
      await tester.ensureVisible(find.text('Run both').first);
      await settle();
      await tester.drag(find.byType(Scrollable).first, const Offset(0, 260));
      await settle();
      await shot('12c_run_both');
      // A wide window: IQ-TREE and RAxML-NG side by side, same height.
      tester.view.physicalSize = const Size(1920, 1080);
      await settle();
      await view(find.text('ENGINE 1'), '12d_inference_wide');
      tester.view.physicalSize = const Size(1280, 720);
      await settle();
      // Both programs at once: one row each, with stages and last lines;
      // RAxML-NG failed (the Windows crash code is explained).
      final pi = TreeProgress(stages: [TreeStage.prepare, TreeStage.modelFinder, TreeStage.search, TreeStage.support, TreeStage.finish]);
      final pr = TreeProgress(raxml: true, stages: [TreeStage.prepare, TreeStage.modelFinder, TreeStage.search, TreeStage.support, TreeStage.finish]);
      for (final l in ['Reading alignment file concatenated.phy ...', 'ModelFinder will test up to 484 DNA models (sample size: 4395) ...',
          ' No. Model         -LnL         df  AIC          AICc         BIC', ' 1  GTR+F         11966.214    15  23962.428    23962.537    24058.188']) {
        pi.add(l);
      }
      for (final l in ['[00:00:01] Starting model selection with parsimony starting tree', '[00:05:37] ML tree search completed, best tree logLH: -304948.126482',
          '[00:05:37] Starting bootstrapping analysis with 1000 replicates.', '[00:43:32] Bootstrap tree #34, logLikelihood: -305381.560860']) {
        pr.add(l);
      }
      final both = TaskModal(ModalKind.bothTrees, '7 + 8')
        ..trees['IQ-TREE'] = pi
        ..trees['RAxML-NG'] = pr
        ..items['IQ-TREE'] = (ModalItem()
          ..state = ItemState.running
          ..progress = 0.3)
        ..items['RAxML-NG'] = (ModalItem()
          ..state = ItemState.failed
          ..note = 'RAxML-NG exit code -1073741819 (0xC0000005: the program crashed (memory access violation))');
      both.tails['IQ-TREE'] = ['ModelFinder will test up to 484 DNA models (sample size: 4395) ...', ' 1  GTR+F         11966.214    15  23962.428'];
      both.tails['RAxML-NG'] = ['[00:42:02] Bootstrap tree #33, logLikelihood: -304214.270608', '[00:43:32] Bootstrap tree #34, logLikelihood: -305381.560860'];
      c.modal = both;
      c.setThemeMode(c.themeMode);
      await tester.pump(const Duration(milliseconds: 300));
      await shot('12e_both_modal');
      c.closeModal();
      await tester.pump();
      await view(find.text('Open the tree editor'), '14_tree');
      Future<void> tapText(Finder f) async {
        if (f.evaluate().isEmpty) {
          // ignore: avoid_print
          print('skip tap: $f');
          return;
        }
        await tester.ensureVisible(f);
        await settle();
        await tester.tap(f);
        await settle();
      }

      // The full-window editor.
      await tapText(find.text('Open the tree editor').first);
      await shot('14a_tree_editor');
      // Rounded branches.
      await tapText(find.text('Right angles').first);
      await tapText(find.text('Rounded corners').last);
      await tapText(find.text('Labels').first);
      await shot('14b_tree_labels');
      await tapText(find.text('Colours').first);
      // Rows coloured by genus, bars named by family.
      if (find.text('None').evaluate().length >= 2) {
        await tapText(find.text('None').first);
        await tapText(find.text('Genus').last);
        await tapText(find.text('None').first);
        await tapText(find.text('Family').last);
      }
      await shot('14c_tree_colours');
      await tapText(find.text('Branches and support').first);
      final shorten = find.textContaining('Shorten long branches');
      if (shorten.evaluate().isNotEmpty) {
        await tester.ensureVisible(shorten.first);
        await settle();
        await tester.tap(shorten.first, warnIfMissed: false);
        await settle();
      }
      await shot('14d_tree_branches');
      await tapText(find.text('Legend').first);
      await shot('14f_tree_legend');
      await tapText(find.text('Figure and export').first);
      await shot('14e_tree_figure');
      // A click on a branch selects it (the bar of actions shows up).
      await tester.tapAt(const Offset(760, 200));
      await settle();
      await shot('14g_tree_selected');
      // Circular layout.
      await tapText(find.text('Tree').first);
      await tapText(find.text('Rectangular').first);
      await tapText(find.text('Circular').last);
      await shot('14i_tree_circular');
      await tester.tapAt(const Offset(900, 300));
      await settle();
      await shot('14j_tree_circular_selected');
      // Two trees face to face, then swapped.
      await tapText(find.text('Compare trees').first);
      await shot('14h_tree_compare');
      await tapText(find.text('This tree on the right').first);
      await shot('14k_tree_compare_swapped');
      await tapText(find.text('Done').first);
      await view(find.text('Save everything'), '15a_zip');
      await view(find.text('How to Cite'), '15_cite_after_tree');
      await view(find.text('REFERENCES CITED IN THE TEXT'), '15b_cite_references');
      c.setThemeMode(ThemeMode.dark);
      c.setLocale(const Locale('pt'));
      await settle();
      await goTo('Registros Carregados', '13_dark_pt');
    } else if (dataset != null) {
      final files = Directory(dataset).listSync().whereType<File>().where((f) => f.path.endsWith('.fasta')).take(12).map((f) => f.path).toList();
      await tester.runAsync(() async {
        await c.detectTools();
        await c.addFiles(files);
        await Future<void>.delayed(const Duration(milliseconds: 500));
      });
      await tester.pump();
      await shot('03_loaded');
      await goTo('Loaded Records', '04_records');
      await goTo('Sample and Locus Fields', '05_fields');
      await goTo('Select Markers', '06_markers');

      // A handful of plastid genes, then the real pipeline with the
      // bundled MAFFT / trimAl / IQ-TREE.
      c.selectNoMarkers();
      for (final l in ['psbA', 'rbcL', 'matK', 'atpB', 'ndhF', 'rpoC2', 'psaA', 'atpA']) {
        if (c.markers.containsKey(l)) c.toggleMarker(l, true);
      }
      await tester.runAsync(() async {
        await c.prepareLoci();
        await c.alignAll();
        await c.trimAll();
        await c.runSaturation();
        c.applySaturationRecommendation(perLocus: false);
        await c.buildMatrix();
        await c.runIqtree(const IqtreeOptions(model: 'GTR+F+G4', ufboot: 1000, alrt: 0, threads: '2', seed: 1));
      });
      await tester.pump();
      await goTo('Sequence Alignment', '08_alignment');
      await goTo('Saturation Test', '09_saturation');
      await goTo('Concatenate and NEXUS', '11_concat');
      await view(find.text('How to split the matrix into partitions'), '11b_concat_partitions');
      await view(find.text('WHAT WILL BE BUILT'), '11c_concat_summary');
      await goTo('Phylogenetic Inference', '12_inference');
    }
  }, skip: out == null);
}
