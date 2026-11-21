import 'dart:io';
import 'dart:ui' as ui;

import 'package:flutter/rendering.dart';
import 'package:flutter/services.dart';
import 'package:flutter/material.dart';
import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/app.dart';
import 'package:splace/src/services/species_id.dart';
import 'package:splace/src/state/project_controller.dart';
import 'package:splace/src/ui/chrome/id_alignment_dialog.dart';

/// The alignment viewer opens with the query and a GenBank hit drawn.
void main() {
  testWidgets('alignment viewer shows the query and the hit', (tester) async {
    for (final e in {
      'Montserrat': ['Regular', 'Italic', 'SemiBold', 'Bold'].map((w) => 'assets/fonts/Montserrat-$w.ttf'),
      'JetBrainsMono': ['assets/fonts/JetBrainsMono-Regular.ttf', 'assets/fonts/JetBrainsMono-Medium.ttf'],
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
    final c = ProjectController(autoDetectTools: false)
      ..splashVisible = false
      ..tutorialAtStart = false;
    addTearDown(c.dispose);
    final key = GlobalKey();
    await tester.pumpWidget(RepaintBoundary(key: key, child: SplaceApp(controller: c)));
    await tester.pump();
    final q = 'ACGTACGTAC' * 7;
    final r = RecordIdentification('S1', 'Rhinella marina')..queries['COI'] = q;
    r.addEvidence(IdSource.genbank, 'COI', 'Rhinella marina', 98.6, acc: 'MN1.1',
        hit: IdHit(
          source: IdSource.genbank,
          species: 'Rhinella marina',
          accession: 'MN1.1',
          identity: 98.6,
          alignLength: 70,
          qFrom: 1,
          qTo: 70,
          sFrom: 11,
          sTo: 80,
          mismatches: 1,
          qseq: q,
          hseq: '${q.substring(0, 30)}T${q.substring(31)}',
        ));
    showIdAlignments(tester.element(find.byType(Scaffold).first), r);
    for (var i = 0; i < 5; i++) {
      await tester.pump(const Duration(milliseconds: 100));
    }
    expect(tester.takeException(), isNull);
    final out = Platform.environment['SPLACE_SCREENSHOTS'];
    if (out != null) {
      await tester.runAsync(() async {
        final image = await (key.currentContext!.findRenderObject()! as RenderRepaintBoundary).toImage();
        final bytes = await image.toByteData(format: ui.ImageByteFormat.png);
        File('$out/id_alignment.png').writeAsBytesSync(bytes!.buffer.asUint8List());
      });
    }
    // The match line starts under the first base of the query.
    final mid = tester.widgetList<Text>(find.byType(Text)).map((t) => t.data).whereType<String>().firstWhere((d) => d.trimLeft().startsWith('|'));
    expect(mid.length - mid.trimLeft().length, 'Query     1  '.length);
    expect(find.text('MN1.1'), findsOneWidget);
    expect(find.textContaining('Query     1', findRichText: true), findsWidgets);
    expect(find.textContaining('Hit      11', findRichText: true), findsOneWidget);
  });
}
