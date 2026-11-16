import 'package:flutter/material.dart';
import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/app.dart';
import 'package:splace/src/state/project_controller.dart';
import 'package:splace/src/ui/fa.dart';

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

    expect(find.text('SPLACE'), findsWidgets);
    expect(find.text('Pipeline steps'.toUpperCase()), findsOneWidget);
    for (final step in [
      'Load Data', 'Loaded Records', 'Sample and Locus Fields', 'Select Markers',
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
}
