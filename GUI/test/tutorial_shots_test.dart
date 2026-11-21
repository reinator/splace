@Tags(['screenshots'])
library;

import 'dart:convert';
import 'dart:io';
import 'dart:ui' as ui;

import 'package:flutter/gestures.dart';
import 'package:flutter/material.dart';
import 'package:flutter/rendering.dart';
import 'package:flutter/services.dart';
import 'package:flutter_test/flutter_test.dart';
import 'package:path/path.dart' as p;
import 'package:splace/src/app.dart';
import 'package:splace/src/state/project_controller.dart';
import 'package:splace/src/ui/fa.dart';

/// Screens of the real app for the tutorial (screenshots/tut_*.png, then
/// tool/make_tutorial_gifs.py): the top bar and "New project", the tips
/// with the example data, a species name with its GBIF synonyms, the
/// species identification, "Export corrected GenBank" and the barcoding
/// card with a real identification of COI barcodes (GenBank BLAST + BOLD).
///
///   SPLACE_TUTORIAL_SHOTS=screenshots SPLACE_NETWORK=1 flutter test test/tutorial_shots_test.dart
void main() {
  final out = Platform.environment['SPLACE_TUTORIAL_SHOTS'];
  final network = Platform.environment['SPLACE_NETWORK'] == '1';

  testWidgets('tutorial screens', (tester) async {
    if (out == null) {
      markTestSkipped('set SPLACE_TUTORIAL_SHOTS');
      return;
    }
    if (network) HttpOverrides.global = null;
    for (final e in {
      'Montserrat': ['Regular', 'Italic', 'Medium', 'SemiBold', 'Bold', 'ExtraBold'].map((w) => 'assets/fonts/Montserrat-$w.ttf'),
      'JetBrainsMono': ['assets/fonts/JetBrainsMono-Regular.ttf', 'assets/fonts/JetBrainsMono-Medium.ttf'],
      'FASolid': ['assets/fonts/fa-solid-900.ttf'],
      'FABrands': ['assets/fonts/fa-brands-400.ttf'],
    }.entries) {
      final loader = FontLoader(e.key);
      for (final f in e.value) {
        loader.addFont(rootBundle.load(f));
      }
      await loader.load();
    }
    tester.view.physicalSize = const Size(1280, 720);
    tester.view.devicePixelRatio = 1;
    addTearDown(tester.view.reset);

    final key = GlobalKey();
    final c = ProjectController(autoDetectTools: false)
      ..splashVisible = false
      ..tutorialAtStart = false;
    addTearDown(c.dispose);
    await tester.pumpWidget(RepaintBoundary(key: key, child: SplaceApp(controller: c)));
    await tester.runAsync(() async {
      for (final a in ['assets/images/splace_logo.png', 'assets/images/flag_us.png', 'assets/images/flag_es.png', 'assets/images/flag_br.png']) {
        await precacheImage(AssetImage(a), tester.element(find.byType(Scaffold).first));
      }
    });
    await tester.pump();

    Future<void> settle([int ms = 300]) async {
      for (var i = 0; i < 4; i++) {
        await tester.runAsync(() => Future.delayed(Duration(milliseconds: ms ~/ 4)));
        await tester.pump(const Duration(milliseconds: 100));
      }
    }

    Future<void> shot(String name) async {
      await settle();
      await tester.runAsync(() async {
        final image = await (key.currentContext!.findRenderObject()! as RenderRepaintBoundary).toImage();
        final bytes = await image.toByteData(format: ui.ImageByteFormat.png);
        File(p.join(out, '$name.png')).writeAsBytesSync(bytes!.buffer.asUint8List());
      });
    }

    /// Scrolls [f] to about a fifth of the page from the top.
    Future<void> show(Finder f, [double alignment = 0.2]) async {
      final ctx = tester.element(f.first);
      await tester.runAsync(() => Scrollable.ensureVisible(ctx, alignment: alignment));
      await settle();
    }

    // The top bar (with New project) and its confirmation.
    await shot('tut_top_bar');
    await tester.tap(find.byIcon(FA.plus).first);
    await settle();
    await shot('tut_new_project');
    await tester.tap(find.text('Cancel').last);
    await settle();

    // The example buttons of 1c and 1d, one tooltip shown (tooltips.json).
    await show(find.text('Use local GenBank files').last, 0.15);
    final mouse = await tester.createGesture(kind: PointerDeviceKind.mouse);
    await mouse.addPointer(location: Offset.zero);
    addTearDown(mouse.removePointer);
    await mouse.moveTo(tester.getCenter(find.text('Barcoding example')));
    await settle(600);
    await shot('tut_examples');
    await mouse.moveTo(Offset.zero);
    await settle();

    // A species: italics and its synonyms from GBIF.
    final taxon = find.byWidgetPredicate((w) => w is TextField && w.decoration?.hintText == 'Carangidae');
    await show(taxon);
    await tester.enterText(taxon, 'Rhinella marina');
    await settle();
    // The GBIF request must start in the real zone (runAsync) to finish.
    await tester.runAsync(() async {
      await tester.tap(find.text('Also search its synonyms (GBIF)'));
      bool done() =>
          find.textContaining('names searched together').evaluate().isNotEmpty ||
          find.textContaining('GBIF has no synonyms').evaluate().isNotEmpty ||
          find.textContaining('GBIF could not be reached').evaluate().isNotEmpty;
      for (var i = 0; i < 80 && !done(); i++) {
        await Future.delayed(const Duration(milliseconds: 500));
        await tester.pump();
      }
    });
    await settle();
    await show(find.textContaining('Species recognized', findRichText: true));
    await shot('tut_species_search');

    // Individual markers (COI, 12S, CYTB) of the species, NCBI and BOLD.
    await show(find.text('Individual markers'), 0.4);
    await tester.tap(find.text('Individual markers'));
    await settle();
    await show(find.text('Also search BOLD Systems'), 0.3);
    await shot('tut_markers_bold');
    if (network) {
      await tester.runAsync(() async {
        await tester.tap(find.text('Search').last);
        for (var i = 0; i < 240 && find.text('Review genomes before importing').evaluate().isEmpty; i++) {
          await Future.delayed(const Duration(milliseconds: 500));
          await tester.pump();
        }
      });
      await settle(600);
      // The BOLD records first, to see them.
      final filter = find.descendant(of: find.byType(Dialog), matching: find.byType(TextField)).first;
      await tester.enterText(filter, 'BOLD');
      await settle();
      await shot('tut_review_bold');
      await tester.sendKeyEvent(LogicalKeyboardKey.escape);
      await settle();
    }
    await show(find.text('Genomes'), 0.4);
    await tester.tap(find.text('Genomes').first);
    await settle();
    await tester.enterText(taxon, '');
    await settle();

    // A GenBank record: species identification and corrected GenBank.
    final dir = Directory.systemTemp.createTempSync('splace_tut_');
    addTearDown(() => dir.deleteSync(recursive: true));
    final gb = File(p.join(dir.path, 'record.gb'))..writeAsStringSync(File('test/core/export_test.dart').readAsStringSync().split("'''")[1]);
    await tester.runAsync(() => c.addGenBankFiles([gb.path]));
    await settle();
    await show(find.text('Identify species (GenBank + BOLD)'));
    await shot('tut_species_id_card');
    await tester.tap(find.text('Identify species (GenBank + BOLD)').first);
    await settle();
    await shot('tut_species_id_setup');
    await tester.sendKeyEvent(LogicalKeyboardKey.escape);
    await settle();
    c.taxonomyResolved = true;
    c.notifyListeners();
    await settle();
    await show(find.text('Export corrected GenBank'), 0.75);
    await shot('tut_export_genbank');
    c.clearRecords();
    await settle();

    // Barcodes: real COI sequences from NCBI, identified for real.
    if (!network || Platform.environment['SPLACE_TUT_SKIP_BARCODING'] == '1') return;
    final http = HttpClient();
    Future<String> get(String path, Map<String, String> q) async {
      final req = await http.getUrl(Uri.https('eutils.ncbi.nlm.nih.gov', '/entrez/eutils/$path', q));
      return (await req.close()).transform(utf8.decoder).join();
    }

    final fasta = StringBuffer();
    await tester.runAsync(() async {
      var n = 0;
      for (final (sp, code) in [('Rhinella marina', 'RMAR'), ('Anaxyrus americanus', 'AAME')]) {
        final ids = (jsonDecode(await get('esearch.fcgi', {
          'db': 'nucleotide',
          'term': '"$sp"[Organism] AND (COI[Gene] OR COX1[Gene]) AND 550:700[SLEN]',
          'retmax': '2',
          'retmode': 'json',
        }))['esearchresult']['idlist'] as List).cast<String>();
        await Future.delayed(const Duration(milliseconds: 400));
        final text = await get('efetch.fcgi', {'db': 'nucleotide', 'id': ids.join(','), 'rettype': 'fasta', 'retmode': 'text'});
        for (final rec in text.split('>').where((r) => r.trim().isNotEmpty)) {
          final seq = rec.split('\n').skip(1).join().trim();
          fasta.writeln('>$code${(++n).toString().padLeft(2, '0')}_COI\n$seq');
        }
        await Future.delayed(const Duration(milliseconds: 400));
      }
    });
    final fa = File(p.join(dir.path, 'barcodes.fasta'))..writeAsStringSync(fasta.toString());
    await tester.runAsync(() => c.addFiles([fa.path]));
    await settle(600);
    // The Sample and marker dialog: apply its proposal.
    final apply = find.textContaining('marker FASTA file');
    if (apply.evaluate().isNotEmpty) {
      await tester.tap(apply.first);
      await settle(600);
    }
    c.setBarcodingMode(true);
    await settle();
    await tester.runAsync(() => c.identifyBarcodes(), );
    c.closeModal();
    await settle();
    await show(find.text('These sequences are DNA barcodes'));
    await shot('tut_barcoding');
    // The alignments of the first specimen (GenBank from BLAST; BOLD fetched
    // and aligned here).
    await tester.runAsync(() async {
      await tester.tap(find.byIcon(FA.alignLeft).first);
      for (var i = 0; i < 120 && find.textContaining('fetching its sequence').evaluate().isNotEmpty || i < 6; i++) {
        await Future.delayed(const Duration(milliseconds: 500));
        await tester.pump();
      }
    });
    await settle(600);
    await shot('tut_alignments');
    // A BOLD hit: its sequence fetched from the portal and aligned here.
    await tester.runAsync(() async {
      await tester.tap(find.text('BOLD').last);
      for (var i = 0; i < 120 && (find.textContaining('fetching its sequence').evaluate().isNotEmpty || i < 6); i++) {
        await Future.delayed(const Duration(milliseconds: 500));
        await tester.pump();
      }
    });
    await settle(600);
    await shot('tut_alignments_bold');
  }, timeout: const Timeout(Duration(minutes: 20)));
}
