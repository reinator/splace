import 'package:flutter/material.dart';
import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/app.dart';
import 'package:splace/src/state/project_controller.dart';
import 'package:splace/src/ui/common.dart';
import 'package:splace/src/ui/widgets/dialogs.dart';

void main() {
  test('species names inside a sentence are set in italics', () {
    final span = italicNames('The closest match is Bufo marinus (99.8%), while the current name Rhinella marina reaches 98%.',
        ['Rhinella marina', 'Bufo marinus']);
    final italic = [for (final c in span.children!.cast<TextSpan>()) if (c.style?.fontStyle == FontStyle.italic) c.text];
    expect(italic, ['Bufo marinus', 'Rhinella marina']);
  });

  testWidgets('a destructive confirmation has one red button: Cancel is neutral', (tester) async {
    tester.view.physicalSize = const Size(1280, 720);
    tester.view.devicePixelRatio = 1;
    addTearDown(tester.view.reset);
    final c = ProjectController(autoDetectTools: false)
      ..splashVisible = false
      ..tutorialAtStart = false;
    addTearDown(c.dispose);
    await tester.pumpWidget(SplaceApp(controller: c));
    await tester.pump();
    confirmDialog(tester.element(find.byType(Scaffold).first),
        title: 'Remove record?', body: 'Rhinella marina (MN1.1) — it leaves the analysis.', confirm: 'Remove', italic: ['Rhinella marina']);
    await tester.pump(const Duration(milliseconds: 400));
    final dialog = find.byType(Dialog);
    expect(find.descendant(of: dialog, matching: find.widgetWithText(OutlinedButton, 'Cancel')), findsOneWidget);
    expect(find.descendant(of: dialog, matching: find.byType(FilledButton)), findsOneWidget, reason: 'only Remove is filled (red)');
    expect(find.textContaining('Rhinella marina', findRichText: true), findsOneWidget);
  });
}
