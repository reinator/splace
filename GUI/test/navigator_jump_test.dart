import 'dart:io';
import 'dart:ui' as ui;

import 'package:flutter/material.dart';
import 'package:flutter/rendering.dart';
import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/app.dart';
import 'package:splace/src/state/project_controller.dart';

/// Clicking a part of the navigator near the end of the page (1d, which
/// never reaches the top) keeps that part selected after the scroll.
void main() {
  testWidgets('1d stays selected after the jump', (tester) async {
    tester.view.physicalSize = const Size(1280, 720);
    tester.view.devicePixelRatio = 1;
    addTearDown(tester.view.reset);
    final key = GlobalKey();
    final c = ProjectController(autoDetectTools: false)
      ..splashVisible = false
      ..tutorialAtStart = false;
    addTearDown(c.dispose);
    await tester.pumpWidget(RepaintBoundary(key: key, child: SplaceApp(controller: c)));
    await tester.pump();
    // The navigator entry (the first match; the section title is the second).
    await tester.tap(find.text('Use local FASTA files').first);
    for (var i = 0; i < 10; i++) {
      await tester.pump(const Duration(milliseconds: 100));
    }
    final out = Platform.environment['SPLACE_SCREENSHOTS'];
    if (out != null) {
      await tester.runAsync(() async {
        final image = await (key.currentContext!.findRenderObject()! as RenderRepaintBoundary).toImage();
        final bytes = await image.toByteData(format: ui.ImageByteFormat.png);
        File('$out/nav_1d.png').writeAsBytesSync(bytes!.buffer.asUint8List());
      });
    }
    // The selected part is drawn in bold in the navigator.
    final label = tester.widget<Text>(find.text('Use local FASTA files').first);
    expect(label.style?.fontWeight, FontWeight.w600);
    final above = tester.widget<Text>(find.text('Use local GenBank files').first);
    expect(above.style?.fontWeight, isNot(FontWeight.w600));
  });
}
