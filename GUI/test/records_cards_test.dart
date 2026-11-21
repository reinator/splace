import 'dart:io';
import 'dart:ui' as ui;

import 'package:flutter/material.dart';
import 'package:flutter/rendering.dart';
import 'package:flutter/services.dart';
import 'package:flutter_test/flutter_test.dart';
import 'package:path/path.dart' as p;
import 'package:splace/l10n/app_localizations.dart';
import 'package:splace/src/state/project_controller.dart';
import 'package:splace/src/ui/common.dart';
import 'package:splace/src/ui/sections/alignment_section.dart';
import 'package:splace/src/ui/sections/records_section.dart';
import 'package:splace/src/ui/theme.dart';

/// Step 2 with one record: the identification and taxonomy cards (plain
/// background, no icon tile, buttons of one size; "Export corrected
/// GenBank" only in the taxonomy card). Picture with SPLACE_SCREENSHOTS.
void main() {
  testWidgets('step 2 cards', (tester) async {
    for (final (family, files) in [
      ('Montserrat', ['Regular', 'Medium', 'SemiBold', 'Bold'].map((w) => 'assets/fonts/Montserrat-$w.ttf')),
      ('FASolid', ['assets/fonts/fa-solid-900.ttf']),
    ]) {
      final loader = FontLoader(family);
      for (final f in files) {
        loader.addFont(rootBundle.load(f));
      }
      await loader.load();
    }
    tester.view.physicalSize = const Size(1100, 1100);
    tester.view.devicePixelRatio = 1;
    addTearDown(tester.view.reset);
    final c = ProjectController(autoDetectTools: false);
    final dir = Directory.systemTemp.createTempSync('splace_cards_');
    addTearDown(() => dir.deleteSync(recursive: true));
    final gb = File(p.join(dir.path, 'r.gb'))..writeAsStringSync(File('test/core/export_test.dart').readAsStringSync().split("'''")[1]);
    await tester.runAsync(() => c.addGenBankFiles([gb.path]));
    c.taxonomyResolved = true;
    final key = GlobalKey();
    await tester.pumpWidget(RepaintBoundary(
      key: key,
      child: ProjectScope(
        controller: c,
        child: MaterialApp(
          theme: buildSplaceTheme(Brightness.light),
          localizationsDelegates: AppLocalizations.localizationsDelegates,
          supportedLocales: AppLocalizations.supportedLocales,
          home: const Scaffold(body: SingleChildScrollView(padding: EdgeInsets.all(16), child: RecordsSection(active: true))),
        ),
      ),
    ));
    await tester.pump();
    expect(find.text('Export corrected GenBank'), findsOneWidget);
    final out = Platform.environment['SPLACE_SCREENSHOTS'];
    if (out != null) {
      await tester.runAsync(() async {
        final image = await (key.currentContext!.findRenderObject()! as RenderRepaintBoundary).toImage();
        final bytes = await image.toByteData(format: ui.ImageByteFormat.png);
        File('$out/step2_cards.png').writeAsBytesSync(bytes!.buffer.asUint8List());
      });
    }
  });

  testWidgets('step 5: MAFFT settings in three columns, run at the end of 5c', (tester) async {
    tester.view.physicalSize = const Size(1100, 1500);
    tester.view.devicePixelRatio = 1;
    addTearDown(tester.view.reset);
    final c = ProjectController(autoDetectTools: false);
    final key = GlobalKey();
    await tester.pumpWidget(RepaintBoundary(
      key: key,
      child: ProjectScope(
        controller: c,
        child: MaterialApp(
          theme: buildSplaceTheme(Brightness.light),
          localizationsDelegates: AppLocalizations.localizationsDelegates,
          supportedLocales: AppLocalizations.supportedLocales,
          home: const Scaffold(body: SingleChildScrollView(padding: EdgeInsets.all(16), child: AlignmentSection(active: true))),
        ),
      ),
    ));
    await tester.pump();
    // The run button is part of 5c; there is no 5d before aligning.
    expect(find.textContaining('Run MAFFT'), findsOneWidget);
    expect(find.text('5d'), findsNothing);
    final out = Platform.environment['SPLACE_SCREENSHOTS'];
    if (out != null) {
      await tester.runAsync(() async {
        final image = await (key.currentContext!.findRenderObject()! as RenderRepaintBoundary).toImage();
        final bytes = await image.toByteData(format: ui.ImageByteFormat.png);
        File('$out/step5.png').writeAsBytesSync(bytes!.buffer.asUint8List());
      });
    }
  });
}
