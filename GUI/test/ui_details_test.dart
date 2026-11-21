import 'package:flutter/material.dart';
import 'package:flutter/services.dart';
import 'package:flutter_test/flutter_test.dart';
import 'package:splace/l10n/app_localizations.dart';
import 'package:splace/src/app.dart';
import 'package:splace/src/services/species_id.dart';
import 'package:splace/src/state/project_controller.dart';
import 'package:splace/src/ui/chrome/tree_identify_dialog.dart';
import 'package:splace/src/ui/theme.dart';
import 'package:splace/src/ui/widgets/examples.dart';
import 'package:splace/src/ui/widgets/toast.dart';

Widget _app(Widget child) => MaterialApp(
      theme: buildSplaceTheme(Brightness.light),
      localizationsDelegates: AppLocalizations.localizationsDelegates,
      supportedLocales: AppLocalizations.supportedLocales,
      home: Scaffold(body: child),
    );

void main() {
  test('categorical colours never repeat, however many groups', () {
    final colors = [for (var i = 0; i < 80; i++) SplacePalette.categoricalAt(i)];
    expect(colors.toSet(), hasLength(80));
    expect(colors.take(SplacePalette.categorical.length), SplacePalette.categorical);
  });

  test('verdict of an identified tip', () {
    RecordIdentification id(String name, List<(String, double)> hits) {
      final r = RecordIdentification('tip', name);
      for (final (sp, pct) in hits) {
        r.addEvidence(IdSource.genbank, 'COI', sp, pct, acc: 'X$sp');
      }
      return r;
    }

    expect(tipVerdict(null), TipVerdict.none);
    expect(tipVerdict(id('Caranx hippos', [])), TipVerdict.none);
    expect(tipVerdict(id('Caranx hippos', [('Caranx hippos', 99.8)])), TipVerdict.confirmed);
    expect(tipVerdict(id('Caranx hippos', [('Caranx latus', 99.5)])), TipVerdict.otherSpecies);
    expect(tipVerdict(id('Caranx hippos', [('Seriola dumerili', 99.1), ('Caranx hippos', 90)])), TipVerdict.otherGenus);
    expect(tipVerdict(id('Caranx hippos', [('Seriola dumerili', 91)])), TipVerdict.noMatch);
    expect(tipVerdict(id('Caranx hippos', [('Caranx latus', 99.5), ('Caranx hippos', 99.5)])), isNot(TipVerdict.otherSpecies));
  });

  testWidgets('the trimming animation fills the card at any width', (tester) async {
    tester.view.physicalSize = const Size(1900, 1000);
    tester.view.devicePixelRatio = 1;
    addTearDown(tester.view.reset);
    await tester.pumpWidget(_app(const Align(alignment: Alignment.topLeft, child: SizedBox(width: 1800, child: TrimExample(method: '-gappyout')))));
    await tester.pump();
    final frame = tester.getRect(find.byType(TrimExample));
    // The rightmost base of the alignment, before any column is removed.
    final letters = find.descendant(of: find.byType(TrimExample), matching: find.byWidgetPredicate((w) => w is Text && RegExp(r'^[ACGT-]$').hasMatch(w.data ?? '')));
    final right = tester.getRect(letters.last).right;
    expect(frame.right - right, lessThan(60), reason: 'the columns reach the right side of an 1800 px card');
    // Let the cycle and its pause finish (timers).
    await tester.pumpWidget(const SizedBox());
  });

  testWidgets('the saved notice shows the file and goes away by itself', (tester) async {
    await tester.pumpWidget(_app(Builder(
      builder: (context) => TextButton(
        onPressed: () => showSplaceToast(context, title: 'Project saved', detail: 'study.splc', path: r'D:\data'),
        child: const Text('save'),
      ),
    )));
    await tester.tap(find.text('save'));
    await tester.pump();
    await tester.pump(const Duration(milliseconds: 400));
    expect(find.text('Project saved'), findsOneWidget);
    expect(find.text('study.splc'), findsOneWidget);
    expect(find.text(r'D:\data'), findsOneWidget);
    await tester.pump(const Duration(seconds: 5));
    await tester.pump(const Duration(milliseconds: 400));
    expect(find.text('Project saved'), findsNothing);
  });

  testWidgets('the window title shows the project file and whether it is saved', (tester) async {
    tester.view.physicalSize = const Size(1600, 1000);
    tester.view.devicePixelRatio = 1;
    addTearDown(tester.view.reset);
    final titles = <String>[];
    tester.binding.defaultBinaryMessenger.setMockMethodCallHandler(const MethodChannel('splace/window'), (call) async {
      if (call.method == 'setTitle') titles.add((call.arguments as Map)['title'] as String);
      return null;
    });
    addTearDown(() => tester.binding.defaultBinaryMessenger.setMockMethodCallHandler(const MethodChannel('splace/window'), null));
    final c = ProjectController(autoDetectTools: false)..splashVisible = false;
    addTearDown(c.dispose);
    await tester.pumpWidget(SplaceApp(controller: c));
    await tester.pump();
    expect(titles.last, 'SPLACE');
    c.projectFile = r'D:\data\fish.splc';
    c.unsaved.value = false;
    c.setThemeMode(ThemeMode.light);
    await tester.pump();
    expect(titles.last, r'D:\data\fish.splc · saved — SPLACE');
    c.unsaved.value = true;
    await tester.pump();
    expect(titles.last, r'● D:\data\fish.splc · unsaved changes — SPLACE');
    // The check of unsaved changes runs a moment after each change.
    await tester.pump(const Duration(seconds: 1));
  });
}
