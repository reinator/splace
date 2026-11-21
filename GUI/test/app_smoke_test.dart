import 'package:flutter/services.dart';
import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/app.dart';
import 'package:splace/src/state/project_controller.dart';
import 'package:splace/src/ui/common.dart';
import 'package:splace/src/ui/fa.dart';
import 'package:splace/src/ui/shell.dart';

void main() {
  testWidgets('shell mirrors the HTML layout: header, step navigator, log', (tester) async {
    tester.view.physicalSize = const Size(1600, 1000);
    tester.view.devicePixelRatio = 1;
    addTearDown(tester.view.reset);

    final controller = ProjectController(autoDetectTools: false);
    addTearDown(controller.dispose);
    await tester.pumpWidget(SplaceApp(controller: controller));
    await tester.pump();

    // Welcome screen first, as in the HTML app.
    expect(find.text('Start App'), findsOneWidget);
    await tester.tap(find.text('Start App'));
    await tester.pump();
    await tester.pump(const Duration(milliseconds: 400));
    // The tutorial opens over the app; it can be closed at any time.
    expect(find.text('How to use SPLACE'), findsOneWidget);
    await tester.sendKeyEvent(LogicalKeyboardKey.escape);
    await tester.pump();
    await tester.pump(const Duration(milliseconds: 400));
    expect(find.text('How to use SPLACE'), findsNothing);

    expect(find.text('SPLACE'), findsWidgets);
    expect(find.text('Pipeline steps'.toUpperCase()), findsOneWidget);
    for (final step in [
      'Load Data', 'Loaded Records', 'Samples and Markers',
      'Sequence Alignment', 'Saturation Test', 'Concatenate and NEXUS', 'Phylogenetic Inference',
    ]) {
      expect(find.text(step), findsWidgets, reason: step);
    }
    // The activity log stays closed until its button is clicked.
    expect(find.text('Activity log'), findsNothing);
    await tester.tap(find.byIcon(FA.rectangleList).first);
    await tester.pump();
    expect(find.text('Activity log'), findsOneWidget);
    // The option card and its entry (1d) in the step navigator.
    expect(find.text('Use local FASTA files'), findsNWidgets(2));
    expect(find.text('1d'), findsNWidgets(2));
    // Later steps stay hidden until data is loaded (locked in the navigator).
    expect(find.text('Run saturation test'), findsNothing);

    await tester.tap(find.text('PT'));
    await tester.pump();
    expect(find.text('Carregar Dados'), findsWidgets);
  });

  testWidgets('New project: a fresh controller, tools kept, in a new folder', (tester) async {
    tester.view.physicalSize = const Size(1600, 1000);
    tester.view.devicePixelRatio = 1;
    addTearDown(tester.view.reset);
    final first = ProjectController(autoDetectTools: false)
      ..splashVisible = false
      ..toolsChecked = true
      ..workDir = 'first_project';
    addTearDown(first.dispose);
    await tester.pumpWidget(SplaceApp(controller: first));
    await tester.pump();
    first.requestNewProject();
    await tester.pump();
    final next = ProjectScope.read(tester.element(find.byType(SplaceShell)));
    expect(identical(next, first), isFalse);
    expect(next.toolsChecked, isTrue);
    expect(next.records, isEmpty);
    expect(next.workDir, isNot('first_project'));
  });
}
