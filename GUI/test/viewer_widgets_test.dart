import 'dart:convert';
import 'dart:io';
import 'dart:ui' as ui;

import 'package:flutter/material.dart';
import 'package:flutter/rendering.dart';
import 'package:flutter/services.dart';
import 'package:flutter_test/flutter_test.dart';
import 'package:splace/l10n/app_localizations.dart';
import 'package:splace/src/core/fasta.dart';
import 'package:splace/src/core/newick.dart';
import 'package:splace/src/core/saturation.dart';
import 'package:splace/src/core/supermatrix.dart';
import 'package:splace/src/ui/theme.dart';
import 'package:splace/src/ui/widgets/alignment_viewer.dart';
import 'package:splace/src/ui/widgets/examples.dart';
import 'package:splace/src/ui/widgets/kit.dart';
import 'package:splace/src/ui/widgets/salve_widgets.dart';
import 'package:splace/src/ui/widgets/matrix_table.dart';
import 'package:splace/src/ui/widgets/tree_viewer.dart';

/// The alignment viewer's codon and dot views, the "keep the reading frame"
/// animation, the full-width matrix table and the support legend. Pictures
/// with SPLACE_SCREENSHOTS=out_dir.
void main() {
  final out = Platform.environment['SPLACE_SCREENSHOTS'];

  Future<void> fonts() async {
    for (final (family, files) in [
      ('Montserrat', ['Regular', 'Italic', 'Medium', 'SemiBold', 'Bold'].map((w) => 'assets/fonts/Montserrat-$w.ttf')),
      ('FASolid', ['assets/fonts/fa-solid-900.ttf']),
      ('JetBrainsMono', ['assets/fonts/JetBrainsMono-Regular.ttf', 'assets/fonts/JetBrainsMono-Medium.ttf']),
    ]) {
      final loader = FontLoader(family);
      for (final f in files) {
        if (File(f).existsSync()) loader.addFont(rootBundle.load(f));
      }
      await loader.load();
    }
  }

  Future<GlobalKey> pump(WidgetTester tester, Widget child) async {
    tester.view.physicalSize = const Size(1300, 900);
    tester.view.devicePixelRatio = 1;
    addTearDown(tester.view.reset);
    final key = GlobalKey();
    await tester.pumpWidget(RepaintBoundary(
      key: key,
      child: MaterialApp(
        theme: buildSplaceTheme(Brightness.light),
        localizationsDelegates: AppLocalizations.localizationsDelegates,
        supportedLocales: AppLocalizations.supportedLocales,
        home: Scaffold(body: Padding(padding: const EdgeInsets.all(16), child: SingleChildScrollView(child: child))),
      ),
    ));
    await tester.pump();
    return key;
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

  testWidgets('codon view and identical-column dots', (tester) async {
    await fonts();
    // Codon 2 varies synonymously (GCT/GCC), codon 3 changes the amino acid
    // (TAT/AAT), the rest is the same everywhere.
    final records = [
      SeqRecord('Species_A', 'ATGGCTTATTCACTGGGA'),
      SeqRecord('Species_B', 'ATGGCCTATTCACTGGGA'),
      SeqRecord('Species_C', 'ATGGCTAATTCACTGGGA'),
      SeqRecord('Species_D', 'ATGGCT---TCACTGGGA'),
    ];
    final key = await pump(tester, AlignmentViewer(records: records, translatable: true, coding: true, height: 200));
    expect(find.text('Codons'), findsOneWidget);
    await tester.tap(find.text('Codons'));
    await tester.pump();
    // 6 codons: two vary, one of them changes the amino acid.
    expect(find.text('2 of 6 codons vary · 1 change the amino acid'), findsOneWidget);
    await shot(tester, key, 'viewer_codons');
    await tester.tap(find.text('Dots'));
    await tester.pump();
    await shot(tester, key, 'viewer_codons_dots');
    await tester.tap(find.text("DNA view"));
    await tester.pump();
    await shot(tester, key, 'viewer_dna_dots');
    // Back to the default for the other viewers.
    await tester.tap(find.text('Dots'));
    await tester.pump();
  });

  testWidgets('columns removed by trimming are grey in the full view', (tester) async {
    await fonts();
    final records = [
      SeqRecord('Species_A', 'ATGGCTTATTCACTGGGA'),
      SeqRecord('Species_B', 'ATGGCCTA-TCACTGGGA'),
      SeqRecord('Species_C', 'ATGGCTAATTCA---GGA'),
    ];
    final removed = [for (var i = 0; i < 18; i++) i >= 6 && i < 9 || i >= 12 && i < 15];
    final key = await pump(tester, AlignmentViewer(records: records, removed: removed, coding: true, translatable: true, height: 160));
    expect(find.textContaining('6 of 18'), findsWidgets);
    await shot(tester, key, 'viewer_removed_grey');
  });

  test('column-by-column trimming shifts the codon positions; whole codons keep them', () {
    final marked = TrimExample.removedColumns('-gappyout', null);
    expect(TrimFrameExample.outOfPhase(marked, 24), isNotEmpty);
    expect(TrimFrameExample.outOfPhase(TrimFrameExample.wholeCodons(marked), 24), isEmpty);
  });

  testWidgets('trimming animation shows both cases', (tester) async {
    await fonts();
    final key = await pump(
        tester, const SizedBox(width: 1000, child: TrimFrameExample(method: '-gappyout', keepFrame: true)));
    for (final (ms, name) in [(300, 'frame_0'), (1200, 'frame_1'), (1600, 'frame_2'), (2600, 'frame_3'), (1600, 'frame_4')]) {
      await tester.pump(Duration(milliseconds: ms));
      await shot(tester, key, name);
    }
    // Let the timers finish.
    await tester.pumpWidget(const SizedBox());
  });

  testWidgets('a resizable table widens a column by dragging its header edge', (tester) async {
    await fonts();
    const codes = ['EX', 'EW', 'RE', 'CR', 'EN', 'VU', 'NT', 'LC', 'DD', 'NE'];
    final key = await pump(
        tester,
        SizedBox(
          width: 900,
          child: RecordsTable(
            columns: const ['Species', 'Category', 'Common names'],
            flex: const [2, 3, 3],
            rowCount: codes.length,
            resizable: true,
            cellBuilder: (context, r, c) => switch (c) {
              0 => Text('Species ${codes[r]}'),
              1 => Align(alignment: Alignment.centerLeft, child: IucnBadge(code: codes[r], name: 'Category ${codes[r]}', possiblyExtinct: r == 3)),
              _ => const Text('A very long list of common names that does not fit in its column at first'),
            },
          ),
        ));
    await shot(tester, key, 'table_resizable_0');
    final before = tester.getSize(find.text('Species EX')).width;
    final handles = find.byWidgetPredicate((w) => w is MouseRegion && w.cursor == SystemMouseCursors.resizeColumn);
    expect(handles, findsNWidgets(3));
    await tester.drag(handles.at(2), const Offset(400, 0));
    await tester.pump();
    await shot(tester, key, 'table_resizable_1');
    // The last column grew past the box: the table now scrolls sideways.
    expect(find.byType(Scrollbar), findsOneWidget);
    expect(tester.getSize(find.text('Species EX')).width, before);
  });

  testWidgets('the matrix table fills its width and the names can be widened', (tester) async {
    await fonts();
    final m = buildSupermatrix([
      LocusAlignment('COI', [SeqRecord('Caranx_hippos_MK000001', 'ATGCC'), SeqRecord('Caranx_latus', 'ATGCA')]),
      LocusAlignment('CYTB', [SeqRecord('Caranx_hippos_MK000001', 'ATG--'), SeqRecord('Caranx_crysos', 'ATGAA')]),
    ]);
    final key = await pump(tester, MatrixOverviewTable(overview: matrixOverview(m)));
    await shot(tester, key, 'matrix_full_width');
    final table = tester.getRect(find.byType(MatrixOverviewTable));
    final handle = find.byWidgetPredicate((w) => w is MouseRegion && w.cursor == SystemMouseCursors.resizeColumn);
    expect(handle, findsOneWidget);
    final before = tester.getCenter(handle).dx;
    await tester.drag(handle, const Offset(120, 0));
    await tester.pump();
    // The drag starts after the gesture slop (~20 px).
    expect(tester.getCenter(handle).dx, greaterThan(before + 90));
    expect(tester.getRect(find.byType(MatrixOverviewTable)).width, table.width, reason: 'the table keeps the whole width');
    await shot(tester, key, 'matrix_names_wider');
  });

  testWidgets('the legend tab sets the support classes', (tester) async {
    await fonts();
    const newick = '((A_alpha:0.1,B_beta:0.12)99:0.05,((C_gamma:0.2,D_delta:0.1)85:0.03,(E_eps:0.1,F_phi:0.2)60:0.04)97:0.02,G_out:0.3);';
    final key = await pump(tester, const TreeViewer(newick: newick, outgroup: ['G_out']));
    await tester.pumpAndSettle();
    await tester.tap(find.text('Open the tree editor').first);
    await tester.pumpAndSettle();
    await tester.tap(find.text('Legend').first);
    await tester.pumpAndSettle();
    expect(find.text('In the legend'), findsOneWidget);
    expect(find.textContaining(RegExp('support classes', caseSensitive: false)), findsWidgets);
    // "Show support from" 80: the legend drops the weak class.
    final field = find.descendant(of: find.ancestor(of: find.textContaining(RegExp('show support from', caseSensitive: false)), matching: find.byType(Column)).first, matching: find.byType(TextField));
    await tester.enterText(field.first, '80');
    await tester.testTextInput.receiveAction(TextInputAction.done);
    await tester.pumpAndSettle();
    await shot(tester, key, 'tree_legend_min80');
  });

  testWidgets('SALVE heatmap beside the labels and support values that never overlap', (tester) async {
    await fonts();
    // Short internal branches with support: values crowd together.
    const newick = '((((A_alpha:0.1,B_beta:0.1)99:0.002,C_gamma:0.1)97:0.002,(D_delta:0.1,E_eps:0.1)85:0.002)92:0.002,'
        '((F_phi:0.1,G_gam:0.1)60:0.002,H_eta:0.1)75:0.002,I_out:0.3);';
    const info = {
      'A_alpha': TipInfo(species: 'A alpha', status: 'CR', statusName: 'Criticamente em Perigo', endemic: 'Sim', biomes: 'Amazônia, Cerrado'),
      'B_beta': TipInfo(species: 'B beta', status: 'EN', endemic: 'Sim', biomes: 'Mata Atlântica'),
      'C_gamma': TipInfo(species: 'C gamma', status: 'VU', endemic: 'Não', biomes: 'Cerrado, Caatinga'),
      'D_delta': TipInfo(species: 'D delta', status: 'NT', endemic: 'Não', biomes: 'Pantanal'),
      'E_eps': TipInfo(species: 'E eps', status: 'LC', endemic: 'Não', biomes: 'Amazônia'),
      'F_phi': TipInfo(species: 'F phi', status: 'DD', endemic: 'Sim', biomes: 'Pampa'),
    };
    final key = await pump(tester, const TreeViewer(newick: newick, outgroup: ['I_out'], tipInfo: info));
    await tester.pumpAndSettle();
    await tester.tap(find.text('Open the tree editor').first);
    await tester.pumpAndSettle();
    await tester.tap(find.text('Labels').first);
    await tester.pumpAndSettle();
    for (final label in ['Category', 'Endemic to Brazil', 'Biomes']) {
      await tester.tap(find.widgetWithText(ChipButton, label).first);
      await tester.pumpAndSettle();
    }
    await shot(tester, key, 'tree_salve_heatmap');
    // Values with curved arrows: still no two values on top of each other.
    await tester.tap(find.text('Branches and support').first);
    await tester.pumpAndSettle();
    await tester.tap(find.text('Value with an arrow to the branch').first);
    await tester.pumpAndSettle();
    await tester.tap(find.text('Straight').first);
    await tester.pumpAndSettle();
    await tester.tap(find.text('Curved').last);
    await tester.pumpAndSettle();
    await shot(tester, key, 'tree_support_curved');
  });

  testWidgets('species delimitation bars: one button per tool and tree', (tester) async {
    await fonts();
    const newick = '((((A1:0.001,A2:0.002)100:0.08,(B1:0.001,B2:0.0015,B3:0.001)99:0.09)95:0.1,(C1:0.002,C2:0.001)100:0.12)80:0.05,D1:0.3,D2:0.28);';
    const mptp = TreeDelimitation(id: 'iq|multi', label: 'mPTP · IQ-TREE', title: 'mPTP', species: [
      ['A1', 'A2'],
      ['B1', 'B2', 'B3'],
      ['C1', 'C2'],
    ]);
    const ptp = TreeDelimitation(id: 'rx|single', label: 'PTP · RAxML-NG', title: 'PTP', species: [
      ['A1', 'A2', 'B1', 'B2', 'B3'],
      ['C1'],
      ['C2'],
    ]);
    final key = await pump(tester, const TreeViewer(newick: newick, outgroup: ['D1', 'D2'], delimitations: [mptp, ptp]));
    await tester.pumpAndSettle();
    expect(find.textContaining(RegExp('species delimitation', caseSensitive: false)), findsWidgets);
    for (final label in ['mPTP · IQ-TREE', 'PTP · RAxML-NG']) {
      await tester.tap(find.widgetWithText(ChipButton, label).first);
      await tester.pumpAndSettle();
    }
    expect(tester.takeException(), isNull);
    await shot(tester, key, 'tree_delimitation');
    // Off again.
    await tester.tap(find.widgetWithText(ChipButton, 'PTP · RAxML-NG').first);
    await tester.pumpAndSettle();
    expect(tester.takeException(), isNull);
  });

  testWidgets('the figure of a tree is reported and comes back as it was left', (tester) async {
    await fonts();
    const newick = '((((A1:0.001,A2:0.002)100:0.08,(B1:0.001,B2:0.0015,B3:0.001)99:0.09)95:0.1,(C1:0.002,C2:0.001)100:0.12)80:0.05,D1:0.3,D2:0.28);';
    const mptp = TreeDelimitation(id: 'iq|multi', label: 'mPTP · IQ-TREE', title: 'mPTP', species: [
      ['A1', 'A2'],
      ['B1', 'B2', 'B3'],
      ['C1', 'C2'],
    ]);
    const ptp = TreeDelimitation(id: 'rx|single', label: 'PTP · RAxML-NG', title: 'PTP', species: [
      ['A1', 'A2', 'B1', 'B2', 'B3'],
      ['C1'],
      ['C2'],
    ]);
    Map<String, dynamic>? reported;
    await pump(tester, TreeViewer(key: const ValueKey(1), newick: newick, outgroup: const ['D1', 'D2'], delimitations: const [mptp, ptp], onFigureChanged: (f) => reported = f));
    await tester.pumpAndSettle();
    await tester.tap(find.widgetWithText(ChipButton, 'mPTP · IQ-TREE').first);
    await tester.pump(const Duration(seconds: 1));
    expect(reported, isNotNull);
    expect(reported!['delimitations'], ['iq|multi']);
    // An edited figure: another title, size and colour, a collapsed clade
    // and the tree rotated by hand (its children in another order).
    final saved = Map<String, dynamic>.of(reported!);
    final tree = parseNewick(newick);
    final rotated = tree.children.first.children.first;
    rotated.children.setAll(0, rotated.children.reversed.toList());
    saved
      ..['tree'] = jsonDecode(jsonEncode(_treeJson(tree)))
      ..['rooting'] = 'asIs'
      ..['collapsed'] = [3]
      ..['title'] = 'Figure 2'
      ..['font'] = 11.5
      ..['branchColor'] = 0xFF1D4ED8;
    reported = null;
    await pump(tester, TreeViewer(key: const ValueKey(2), newick: newick, outgroup: const ['D1', 'D2'], delimitations: const [mptp, ptp], figure: saved, onFigureChanged: (f) => reported = f));
    await tester.pumpAndSettle();
    await tester.tap(find.widgetWithText(ChipButton, 'PTP · RAxML-NG').first);
    await tester.pump(const Duration(seconds: 1));
    expect(reported, isNotNull);
    expect(reported!['delimitations'], ['iq|multi', 'rx|single']);
    for (final k in ['tree', 'rooting', 'collapsed', 'title', 'font', 'branchColor']) {
      expect(reported![k], saved[k], reason: k);
    }
    // A saved tree with other tips (the analysis was run again): only the
    // options come back.
    final other = Map<String, dynamic>.of(saved)..['tree'] = _treeJson(parseNewick('(X:1,Y:1,Z:1);'));
    reported = null;
    await pump(tester, TreeViewer(key: const ValueKey(3), newick: newick, outgroup: const ['D1', 'D2'], delimitations: const [mptp, ptp], figure: other, onFigureChanged: (f) => reported = f));
    await tester.pumpAndSettle();
    await tester.tap(find.widgetWithText(ChipButton, 'PTP · RAxML-NG').first);
    await tester.pump(const Duration(seconds: 1));
    expect(reported!['title'], 'Figure 2');
    expect(reported!['collapsed'], isEmpty);
    expect(tester.takeException(), isNull);
  });

  test('support legend follows "show support from"', () {
    List<String> labels(int s, int m, int min) => [for (final (_, l) in supportClassLabels(s, m, min)) l];
    expect(labels(95, 80, 0), ['≥ 95', '80–94', '< 80']);
    expect(labels(95, 80, 80), ['≥ 95', '80–94']);
    expect(labels(95, 80, 85), ['≥ 95', '85–94']);
    expect(labels(95, 80, 95), ['≥ 95']);
    expect(labels(95, 80, 50), ['≥ 95', '80–94', '50–79']);
    expect(labels(70, 50, 0), ['≥ 70', '50–69', '< 50']);
    // Any number of classes.
    List<String> ranges(List<int> f, int min) => [for (final (_, l) in supportClassRanges(f, min)) l];
    expect(ranges([98, 95, 80, 0], 0), ['≥ 98', '95–97', '80–94', '< 80']);
    expect(ranges([95, 0], 50), ['≥ 95', '50–94']);
  });

  test('saturation is summed up gene by gene, never concatenated', () {
    SaturationReport report(String name, SaturationVerdict p3) {
      SubsetSaturation sub(CodonSubset s, SaturationVerdict v) => SubsetSaturation(
            subset: s,
            sites: 200,
            pairs: 45,
            undefinedPairs: 0,
            meanP: 0.1,
            meanDistance: 0.12,
            slope: 0.8,
            r2: 0.9,
            tsTvRatio: 3,
            iss: 0.2,
            baseFrequencies: const [0.25, 0.25, 0.25, 0.25],
            model: 'TN93',
            points: const [],
            verdict: v,
          );
      return SaturationReport(name, 10, 600, 0, {
        CodonSubset.pos1: sub(CodonSubset.pos1, SaturationVerdict.ok),
        CodonSubset.pos2: sub(CodonSubset.pos2, SaturationVerdict.ok),
        CodonSubset.pos3: sub(CodonSubset.pos3, p3),
      }, null);
    }

    final pooled = PooledSaturation.of([
      report('COI', SaturationVerdict.saturated),
      report('CYTB', SaturationVerdict.saturated),
      report('ND2', SaturationVerdict.ok),
    ])!;
    expect(pooled.count(3, SaturationVerdict.saturated), 2);
    expect(pooled.verdict(3), SaturationVerdict.saturated, reason: '2 of 3 genes');
    expect(pooled.recommendation.exclude, {3});
    final half = PooledSaturation.of([report('COI', SaturationVerdict.saturated), report('ND2', SaturationVerdict.ok)])!;
    expect(half.verdict(3), SaturationVerdict.moderate, reason: 'saturated in half the genes: kept, with caution');
    expect(half.recommendation.kind, RecommendationKind.keepAll);
    expect(half.recommendation.caution, {3});
  });
}

/// A tree as the viewer saves it (nested nodes).
Map<String, dynamic> _treeJson(TreeNode n) => {
      if (n.name.isNotEmpty) 'n': n.name,
      if (n.length != null) 'l': n.length,
      if (n.label != null) 's': n.label,
      if (n.children.isNotEmpty) 'c': [for (final k in n.children) _treeJson(k)],
    };
