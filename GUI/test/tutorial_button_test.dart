import 'package:flutter/material.dart';
import 'package:flutter/services.dart';
import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/app.dart';
import 'package:splace/src/state/project_controller.dart';

/// The "Tutorial" button of the welcome screen opens the tutorial over it,
/// and every page renders (scenes and pictures from the assets).
void main() {
  testWidgets('Tutorial button of the welcome screen', (tester) async {
    tester.view.physicalSize = const Size(1280, 720);
    tester.view.devicePixelRatio = 1;
    addTearDown(tester.view.reset);
    final c = ProjectController(autoDetectTools: false);
    addTearDown(c.dispose);
    await tester.pumpWidget(SplaceApp(controller: c));
    await tester.pump();
    await tester.tap(find.widgetWithText(FilledButton, 'Tutorial'));
    await tester.pump();
    for (var i = 0; i < 3; i++) {
      await tester.runAsync(() => Future.delayed(const Duration(milliseconds: 200)));
      await tester.pump(const Duration(milliseconds: 200));
    }
    expect(find.text('How to use SPLACE'), findsOneWidget);
    for (var page = 0; page < 9; page++) {
      await tester.sendKeyEvent(LogicalKeyboardKey.arrowRight);
      await tester.runAsync(() => Future.delayed(const Duration(milliseconds: 100)));
      await tester.pump(const Duration(milliseconds: 100));
    }
    expect(tester.takeException(), isNull);
    expect(find.text('8 · Tree Viewer'), findsOneWidget);
    await tester.sendKeyEvent(LogicalKeyboardKey.escape);
    for (var i = 0; i < 4; i++) {
      await tester.pump(const Duration(milliseconds: 300));
    }
    expect(find.text('How to use SPLACE'), findsNothing);
  });
}
