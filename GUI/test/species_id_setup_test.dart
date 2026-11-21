import 'dart:io';
import 'dart:ui' as ui;

import 'package:flutter/material.dart';
import 'package:flutter/rendering.dart';
import 'package:flutter/services.dart';
import 'package:flutter_test/flutter_test.dart';
import 'package:path/path.dart' as p;
import 'package:splace/src/app.dart';
import 'package:splace/src/state/project_controller.dart';

/// The identification dialog at the smallest window (1280×720): genes and
/// both parameter cards in view, nothing to scroll, the cards the same
/// height.
void main() {
  testWidgets('species identification setup fits without scrolling', (tester) async {
    for (final e in {
      'Montserrat': ['Regular', 'Italic', 'Medium', 'SemiBold', 'Bold'].map((w) => 'assets/fonts/Montserrat-$w.ttf'),
      'FASolid': ['assets/fonts/fa-solid-900.ttf'],
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
    final dir = Directory.systemTemp.createTempSync('splace_setup_');
    addTearDown(() => dir.deleteSync(recursive: true));
    final key = GlobalKey();
    final c = ProjectController(autoDetectTools: false)
      ..splashVisible = false
      ..tutorialAtStart = false
      ..workDir = p.join(dir.path, 'project');
    addTearDown(c.dispose);
    await tester.pumpWidget(RepaintBoundary(key: key, child: SplaceApp(controller: c)));
    final gb = File(p.join(dir.path, 'record.gb'))..writeAsStringSync(File('test/core/export_test.dart').readAsStringSync().split("'''")[1]);
    await tester.runAsync(() => c.addGenBankFiles([gb.path]));
    await tester.pump();
    await tester.ensureVisible(find.text('Identify species (GenBank + BOLD)').first);
    await tester.pump(const Duration(milliseconds: 300));
    await tester.tap(find.text('Identify species (GenBank + BOLD)').first);
    for (var i = 0; i < 6; i++) {
      await tester.pump(const Duration(milliseconds: 100));
    }
    final out = Platform.environment['SPLACE_SCREENSHOTS'];
    if (out != null) {
      await tester.runAsync(() async {
        final image = await (key.currentContext!.findRenderObject()! as RenderRepaintBoundary).toImage();
        final bytes = await image.toByteData(format: ui.ImageByteFormat.png);
        File('$out/species_id_setup.png').writeAsBytesSync(bytes!.buffer.asUint8List());
      });
    }
    expect(find.text('GenBank (NCBI BLAST)'), findsOneWidget);
    expect(find.text('BOLD Systems (ID engine)'), findsOneWidget);
    // Nothing to scroll inside the dialog.
    for (final e in find.descendant(of: find.byType(Dialog), matching: find.byType(Scrollable)).evaluate()) {
      final pos = (e as StatefulElement).state as ScrollableState;
      expect(pos.position.maxScrollExtent, 0, reason: 'the dialog scrolls');
    }
    // The two parameter cards are the same height.
    Size card(String title) => tester.getSize(find.ancestor(of: find.text(title), matching: find.byType(Container)).first);
    expect(card('GenBank (NCBI BLAST)').height, card('BOLD Systems (ID engine)').height);
  });
}
