@Tags(['screenshots'])
library;

import 'dart:io';
import 'dart:math' as math;
import 'dart:ui' as ui;

import 'package:flutter/material.dart';
import 'package:flutter/rendering.dart';
import 'package:flutter/services.dart';
import 'package:flutter_test/flutter_test.dart';
import 'package:splace/l10n/app_localizations.dart';
import 'package:splace/src/ui/theme.dart';
import 'package:splace/src/ui/widgets/tree_viewer.dart';

/// Tree editor on a larger random tree (60 tips), to look at the circular
/// layout and the comparison with many rows. Opt-in like the app shots:
///
///   SPLACE_SCREENSHOTS=out_dir flutter test test/tree_viewer_shots_test.dart
String _randomTree(int n, math.Random rnd, {String Function(int)? name}) {
  final nodes = [for (var i = 0; i < n; i++) '${(name ?? (i) => 'Genus${i ~/ 8}_species_$i')(i)}:${(0.01 + rnd.nextDouble() * 0.2).toStringAsFixed(3)}'];
  while (nodes.length > 2) {
    final a = nodes.removeAt(rnd.nextInt(nodes.length));
    final b = nodes.removeAt(rnd.nextInt(nodes.length));
    nodes.add('($a,$b)${50 + rnd.nextInt(51)}:${(0.005 + rnd.nextDouble() * 0.08).toStringAsFixed(3)}');
  }
  return '(${nodes.join(',')});';
}

void main() {
  final out = Platform.environment['SPLACE_SCREENSHOTS'];

  testWidgets('circular tree and comparison shots', (tester) async {
    for (final (family, files) in [
      ('Montserrat', ['Regular', 'Italic', 'Medium', 'SemiBold', 'Bold', 'BoldItalic'].map((w) => 'assets/fonts/Montserrat-$w.ttf')),
      ('FASolid', ['assets/fonts/fa-solid-900.ttf']),
    ]) {
      final loader = FontLoader(family);
      for (final f in files) {
        if (File(f).existsSync()) loader.addFont(rootBundle.load(f));
      }
      await loader.load();
    }
    tester.view.physicalSize = const Size(1440, 900);
    tester.view.devicePixelRatio = 1;
    addTearDown(tester.view.reset);
    // A few long names and titles: they must stay in their half.
    String name(int i) => i % 9 == 0 ? 'Genus${i ~/ 8}_species_with_a_long_epithet_voucher_ITV0${1400 + i}_25_Castanheira' : 'Genus${i ~/ 8}_species_$i';
    final a = _randomTree(60, math.Random(1), name: name);
    final b = _randomTree(60, math.Random(2), name: name);
    final key = GlobalKey();
    await tester.pumpWidget(RepaintBoundary(
      key: key,
      child: MaterialApp(
        theme: buildSplaceTheme(Brightness.light),
        localizationsDelegates: AppLocalizations.localizationsDelegates,
        supportedLocales: AppLocalizations.supportedLocales,
        home: Scaffold(
          body: SingleChildScrollView(
            child: TreeViewer(
              newick: a,
              title: 'Tree A · IQ-TREE 3 · 23:41 · ML tree, SH-aLRT / UFBoot with a long title',
              compareWith: [TreeOption(label: 'Tree B · RAxML-NG · 23:56 · ML tree, predicted bootstrap (EBG) and more', load: () => b)],
            ),
          ),
        ),
      ),
    ));
    await tester.pumpAndSettle();

    Future<void> shot(String name) async {
      await tester.pumpAndSettle();
      if (out == null) return;
      await tester.runAsync(() async {
        final image = await (key.currentContext!.findRenderObject()! as RenderRepaintBoundary).toImage();
        final bytes = await image.toByteData(format: ui.ImageByteFormat.png);
        File('$out/$name.png').writeAsBytesSync(bytes!.buffer.asUint8List());
      });
    }

    Future<void> tap(Finder f) async {
      await tester.ensureVisible(f);
      await tester.pumpAndSettle();
      await tester.tap(f);
      await tester.pumpAndSettle();
    }

    await tap(find.text('Open the tree editor').first);
    await tap(find.text('Rectangular').first);
    await tap(find.text('Circular').last);
    await tap(find.text('Whole figure').first);
    await shot('t1_circular_60');
    // A click on a branch selects its clade.
    await tester.tapAt(const Offset(1000, 420));
    await shot('t2_circular_selected');
    // Half circle.
    await tester.tap(find.text('Colours').first);
    await tester.pumpAndSettle();
    await shot('t3_circular_colours');
    await tap(find.text('Compare trees').first);
    await shot('t4_compare');
    await tap(find.text('This tree on the right').first);
    await shot('t5_compare_swapped');
    // Support values with an arrow to the branch.
    await tap(find.text('Branches and support').first);
    await tap(find.text('Value with an arrow to the branch').first);
    await shot('t6_compare_arrows');
    await tap(find.text('Edit figure').first);
    await shot('t7_circular_arrows');
    await tap(find.text('Tree').first);
    await tap(find.text('Circular').first);
    await tap(find.text('Rectangular').last);
    await tap(find.text('Whole figure').first);
    await shot('t8_rect_arrows');
  });
}
