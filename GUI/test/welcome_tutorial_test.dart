import 'dart:io';
import 'dart:ui' as ui;

import 'package:flutter/material.dart';
import 'package:flutter/rendering.dart';
import 'package:flutter/services.dart';
import 'package:flutter_test/flutter_test.dart';
import 'package:splace/l10n/app_localizations.dart';
import 'package:splace/src/state/project_controller.dart';
import 'package:splace/src/ui/chrome/splash.dart';
import 'package:splace/src/ui/chrome/tutorial_dialog.dart';
import 'package:splace/src/ui/common.dart';
import 'package:splace/src/ui/theme.dart';

/// The welcome screen and the tutorial carousel (pictures with
/// SPLACE_SCREENSHOTS=out_dir).
void main() {
  final out = Platform.environment['SPLACE_SCREENSHOTS'];

  Future<void> fonts() async {
    for (final (family, files) in [
      ('Montserrat', ['Regular', 'Italic', 'Medium', 'SemiBold', 'Bold', 'ExtraBold'].map((w) => 'assets/fonts/Montserrat-$w.ttf')),
      ('FASolid', ['assets/fonts/fa-solid-900.ttf']),
    ]) {
      final loader = FontLoader(family);
      for (final f in files) {
        if (File(f).existsSync()) loader.addFont(rootBundle.load(f));
      }
      await loader.load();
    }
  }

  Future<void> shot(WidgetTester tester, GlobalKey key, String name) async {
    await tester.pump(const Duration(milliseconds: 50));
    if (out == null) return;
    await tester.runAsync(() async {
      final image = await (key.currentContext!.findRenderObject()! as RenderRepaintBoundary).toImage();
      final bytes = await image.toByteData(format: ui.ImageByteFormat.png);
      File('$out/$name.png').writeAsBytesSync(bytes!.buffer.asUint8List());
    });
  }

  testWidgets('welcome screen, then the tutorial pages', (tester) async {
    await fonts();
    tester.view.physicalSize = const Size(1280, 720);
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
          locale: const Locale('pt'),
          home: const Scaffold(body: SplashOverlay()),
        ),
      ),
    ));
    await tester.pump();
    expect(find.textContaining('itvgenomics/splace'), findsOneWidget);
    await shot(tester, key, 'welcome');
    final ctx = tester.element(find.byType(SplashOverlay));
    showTutorial(ctx, start: 1);
    // Scenes (JSON) and screenshots load asynchronously.
    Future<void> settle() async {
      for (var i = 0; i < 3; i++) {
        await tester.runAsync(() => Future.delayed(const Duration(milliseconds: 300)));
        await tester.pump();
      }
    }

    await tester.pump();
    await tester.pump(const Duration(milliseconds: 300));
    await settle();
    expect(find.text('Carregar Dados'), findsWidgets);
    await shot(tester, key, 'tutorial_step1_item1');
    // The third item: its scene, highlight and arrow.
    await tester.tap(find.textContaining('(duas palavras)', findRichText: true));
    await settle();
    await shot(tester, key, 'tutorial_step1_item3');
    await tester.sendKeyEvent(LogicalKeyboardKey.arrowRight);
    await settle();
    await shot(tester, key, 'tutorial_step2');
    // Step 2 at the minimum window: every item visible, no scrolling.
    expect(find.byType(Scrollable).evaluate().where((e) => e.findAncestorWidgetOfExactType<Dialog>() != null), isEmpty);
    expect(find.textContaining('busque, edite ou remova', findRichText: true), findsOneWidget);
    // The last item shows the next step on the progress bar.
    await tester.tap(find.textContaining('busque, edite ou remova', findRichText: true));
    await settle();
    expect(find.textContaining('A seguir: 3'), findsOneWidget);
    await shot(tester, key, 'tutorial_step2_last');
    expect(tester.takeException(), isNull);
    // Closing the dialog stops its timer.
    await tester.sendKeyEvent(LogicalKeyboardKey.escape);
    await tester.pump(const Duration(milliseconds: 400));
  });
}
