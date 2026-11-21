import 'dart:async';
import 'dart:convert';
import 'dart:math' as math;

import 'package:flutter/material.dart';
import 'package:flutter/rendering.dart';
import 'package:flutter/services.dart';

import '../common.dart';
import '../fa.dart';
import '../theme.dart';
import '../widgets/dialogs.dart';
import '../widgets/kit.dart';

/// The tutorial after "Start App": one page per step (and one on the
/// interface) with what the step is for, what to do there and a tip. Each
/// "what to do" item has a scene (`assets/tutorial/scenes.json`): a screen
/// of the real app with the place to fill in highlighted and a curved
/// arrow from the item to it. The scenes play one after the other, like a
/// GIF; clicking an item shows its scene. A carousel: arrows, the step
/// pills or the keyboard (← →); it can be closed at any time.
Future<void> showTutorial(BuildContext context, {int start = 0}) => showDialog(
      context: context,
      builder: (_) => ProjectScope(controller: ProjectScope.read(context), child: _TutorialDialog(start: start)),
    );

/// One scene: a screenshot (1280×720) and the area to point at.
class _Scene {
  const _Scene(this.shot, this.rect);
  final String shot;
  final Rect? rect;
}

/// Width and height of the screenshots of the scenes.
const _shotSize = Size(1280, 720);

/// Where a picture of [_shotSize] is drawn inside [box] (BoxFit.contain).
Rect _fitted(Size box) {
  final s = math.min(box.width / _shotSize.width, box.height / _shotSize.height);
  final w = _shotSize.width * s, h = _shotSize.height * s;
  return Rect.fromLTWH((box.width - w) / 2, (box.height - h) / 2, w, h);
}

class _Page {
  const _Page(this.title, this.intro, this.todo, this.tip);
  final String title, intro, tip;
  final List<String> todo;
}

class _TutorialDialog extends StatefulWidget {
  const _TutorialDialog({required this.start});
  final int start;

  @override
  State<_TutorialDialog> createState() => _TutorialDialogState();
}

class _TutorialDialogState extends State<_TutorialDialog> with SingleTickerProviderStateMixin {
  late int _page = widget.start;

  /// The item (and scene) shown.
  int _item = 0;
  List<List<_Scene>> _scenes = const [];
  bool _paused = false;
  final _focus = FocusNode();
  final _stackKey = GlobalKey(), _imageKey = GlobalKey();
  final _itemKeys = List.generate(10, (_) => GlobalKey());

  /// Time of the item shown: the bar of the item fills, then the next item
  /// (or, after the last one, the next step) comes.
  late final _clock = AnimationController(vsync: this, duration: const Duration(milliseconds: 4800))..addStatusListener(_onTick);

  /// The arrow from the item to its area (stack coordinates), when both
  /// are visible: start, the x of its vertical part, end.
  (Offset, double, Offset)? _arrow;

  /// Number of pages (the interface and the eight steps).
  static const _pageCount = 9;

  @override
  void initState() {
    super.initState();
    // Without the scenes (an app built before they existed) the pages show
    // their text only.
    rootBundle.loadString('assets/tutorial/scenes.json').catchError((Object e) {
      debugPrint('Tutorial scenes not found: $e');
      return '{"pages": []}';
    }).then((text) {
      final j = jsonDecode(text) as Map<String, dynamic>;
      final pages = [
        for (final page in j['pages'] as List)
          [
            for (final sc in page as List)
              () {
                final r = (sc as Map)['rect'] as List?;
                return _Scene(sc['shot'] as String,
                    r == null ? null : Rect.fromLTWH((r[0] as num).toDouble(), (r[1] as num).toDouble(), (r[2] as num).toDouble(), (r[3] as num).toDouble()));
              }(),
          ],
      ];
      if (mounted) setState(() => _scenes = pages);
    });
    _clock.forward();
  }

  @override
  void dispose() {
    _clock.dispose();
    _focus.dispose();
    super.dispose();
  }

  List<_Scene> get _pageScenes => _page < _scenes.length ? _scenes[_page] : const [];

  /// Items of the page shown (scenes, or at least one).
  int get _itemCount => math.max(1, _pageScenes.length);

  void _onTick(AnimationStatus status) {
    if (status != AnimationStatus.completed || !mounted) return;
    if (_item < _itemCount - 1) {
      setState(() => _item++);
      _clock.forward(from: 0);
    } else if (_page < _pageCount - 1) {
      _go(_page + 1);
    }
  }

  void _go(int i) {
    setState(() {
      _page = i.clamp(0, _pageCount - 1);
      _item = 0;
      _arrow = null;
    });
    _clock.value = 0;
    if (!_paused) _clock.forward();
  }

  void _show(int item) {
    setState(() => _item = item);
    _clock.value = 0;
    if (!_paused) _clock.forward();
  }

  void _togglePause() {
    setState(() => _paused = !_paused);
    _paused ? _clock.stop() : _clock.forward();
  }

  /// Arrow from the left of the active item to the right side of its area
  /// in the picture, in straight lines with right angles: left to the gap
  /// between picture and text, up or down to the area, left to it.
  void _updateArrow() {
    if (!mounted) return;
    (Offset, double, Offset)? next;
    final stack = _stackKey.currentContext?.findRenderObject() as RenderBox?;
    final image = _imageKey.currentContext?.findRenderObject() as RenderBox?;
    final item = _item < _itemKeys.length ? _itemKeys[_item].currentContext?.findRenderObject() as RenderBox? : null;
    final scenes = _pageScenes;
    final rect = _item < scenes.length ? scenes[_item].rect : null;
    if (stack != null && image != null && item != null && rect != null && stack.attached && image.attached && item.attached) {
      final fit = _fitted(image.size);
      final s = fit.width / _shotSize.width;
      final origin = image.localToGlobal(Offset.zero, ancestor: stack);
      final target = Rect.fromLTWH(fit.left + rect.left * s, fit.top + rect.top * s, rect.width * s, rect.height * s).shift(origin);
      final at = item.localToGlobal(Offset.zero, ancestor: stack);
      final itemSize = MatrixUtils.transformRect(item.getTransformTo(stack), Offset.zero & item.size);
      final start = Offset(at.dx - 6, at.dy + math.min(11.0, itemSize.height / 2));
      final gap = origin.dx + image.size.width + (start.dx - origin.dx - image.size.width) / 2;
      final end = Offset(target.right + 4, start.dy.clamp(target.top + 6, target.bottom - 6).toDouble());
      next = (start, gap, end);
    }
    if (next != _arrow) setState(() => _arrow = next);
  }

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    final lang = Localizations.localeOf(context).languageCode;
    final pages = _pages(lang, [
      l10n.tutorialInterface,
      l10n.step1, l10n.step2, l10n.stepMarkers, l10n.step5, l10n.step6, l10n.step7, l10n.step8, l10n.step9,
    ]);
    final p = pages[_page];
    final last = _page == pages.length - 1;
    final scenes = _pageScenes;
    final scene = _item < scenes.length ? scenes[_item] : (scenes.isEmpty ? null : scenes.first);
    WidgetsBinding.instance.addPostFrameCallback((_) => _updateArrow());

    final bold = TextStyle(fontWeight: FontWeight.w700, color: t.accentText);
    // Text of the page: it always fits (scaled down when the window is
    // small), so every item is seen without scrolling.
    final text = Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
      Row(children: [
        Container(
          width: 32,
          height: 32,
          alignment: Alignment.center,
          decoration: const BoxDecoration(color: SplaceRed.c700, shape: BoxShape.circle),
          child: _page == 0
              ? const Icon(FA.circleInfo, size: 14, color: Colors.white)
              : Text('$_page', style: const TextStyle(fontWeight: FontWeight.w700, color: Colors.white)),
        ),
        const SizedBox(width: 12),
        Expanded(child: Text(p.title, style: TextStyle(fontSize: 19, fontWeight: FontWeight.w700, color: t.text))),
      ]),
      const SizedBox(height: 10),
      Text.rich(TextSpan(children: _rich(p.intro, bold, t.accentText, 13)),
          textAlign: TextAlign.justify, style: TextStyle(fontSize: 14, color: t.text, height: 1.5)),
      const SizedBox(height: 12),
      Kicker(l10n.tutorialWhatToDo, size: 10.5),
      const SizedBox(height: 6),
      for (final (i, item) in p.todo.indexed)
        Padding(
          key: i < _itemKeys.length ? _itemKeys[i] : null,
          padding: const EdgeInsets.only(bottom: 4),
          child: Material(
            color: i == _item ? t.accentSoft : Colors.transparent,
            borderRadius: BorderRadius.circular(8),
            child: InkWell(
              borderRadius: BorderRadius.circular(8),
              onTap: () => _show(i),
              child: Padding(
                padding: const EdgeInsets.symmetric(horizontal: 6, vertical: 4),
                child: Row(crossAxisAlignment: CrossAxisAlignment.start, children: [
                  Container(
                    width: 20,
                    height: 20,
                    margin: const EdgeInsets.only(top: 1),
                    alignment: Alignment.center,
                    decoration: BoxDecoration(color: i == _item ? SplaceRed.c700 : t.accentSoft, shape: BoxShape.circle),
                    child: Text('${i + 1}', style: TextStyle(fontSize: 11, fontWeight: FontWeight.w700, color: i == _item ? Colors.white : t.accentText)),
                  ),
                  const SizedBox(width: 10),
                  Expanded(
                    child: Text.rich(TextSpan(children: _rich(item, bold, t.accentText, 12.5)),
                        textAlign: TextAlign.justify, style: TextStyle(fontSize: 13.5, color: t.text, height: 1.45)),
                  ),
                ]),
              ),
            ),
          ),
        ),
      const SizedBox(height: 8),
      TonePanel(
        tone: PanelTone.info,
        icon: FA.lightbulb,
        child: Text.rich(TextSpan(children: _rich(p.tip, bold, t.accentText, 12)),
            textAlign: TextAlign.justify, style: const TextStyle(fontSize: 13, height: 1.45)),
      ),
    ]);

    return Focus(
      focusNode: _focus,
      autofocus: true,
      onKeyEvent: (_, e) {
        if (e is! KeyDownEvent) return KeyEventResult.ignored;
        if (e.logicalKey == LogicalKeyboardKey.arrowRight) {
          _go(_page + 1);
          return KeyEventResult.handled;
        }
        if (e.logicalKey == LogicalKeyboardKey.arrowLeft) {
          _go(_page - 1);
          return KeyEventResult.handled;
        }
        return KeyEventResult.ignored;
      },
      child: SplaceDialog(
        icon: FA.circleQuestion,
        title: l10n.tutorialTitle,
        subtitle: l10n.tutorialSubtitle,
        scroll: false,
        headerTrailing: LangToggle(current: lang, onChanged: (code) => c.setLocale(Locale(code))),
        footerLeading: Row(children: [
          SizedBox(
            width: 22,
            height: 22,
            child: Checkbox(value: c.tutorialAtStart, onChanged: (v) => c.setTutorialAtStart(v ?? true)),
          ),
          const SizedBox(width: 8),
          Flexible(child: Text(l10n.tutorialShowAtStart, style: TextStyle(fontSize: 13, color: t.textMuted))),
        ]),
        actions: [
          OutlinedButton(onPressed: _page == 0 ? null : () => _go(_page - 1), child: Text(l10n.tutorialPrev)),
          FilledButton(
            onPressed: last ? () => Navigator.of(context).pop() : () => _go(_page + 1),
            child: Text(last ? l10n.tutorialFinish : l10n.tutorialNext),
          ),
        ],
        child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
          // The steps in one row: the current one with its name, the others
          // by number (name in the tooltip).
          Row(children: [
            for (var i = 0; i < pages.length; i++) ...[
              if (i > 0) const SizedBox(width: 6),
              i == _page
                  ? ChipButton(label: i == 0 ? pages[0].title : '$i · ${pages[i].title}', style: ChipButtonStyle.accent, onPressed: () => _go(i))
                  : Tooltip(
                      message: pages[i].title,
                      child: ChipButton(label: i == 0 ? 'i' : '$i', onPressed: () => _go(i)),
                    ),
            ],
          ]),
          const SizedBox(height: 14),
          Expanded(
            child: Stack(key: _stackKey, children: [
              Row(crossAxisAlignment: CrossAxisAlignment.start, children: [
                Expanded(
                  flex: 11,
                  child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
                    Expanded(
                      child: Container(
                        key: _imageKey,
                        decoration: BoxDecoration(
                          color: t.surface2,
                          borderRadius: BorderRadius.circular(12),
                          border: Border.all(color: t.border),
                        ),
                        clipBehavior: Clip.antiAlias,
                        child: scene == null
                            ? const SizedBox.expand()
                            : Stack(fit: StackFit.expand, children: [
                                Image.asset(
                                  'assets/tutorial/shots/${scene.shot}.png',
                                  fit: BoxFit.contain,
                                  gaplessPlayback: true,
                                  errorBuilder: (_, _, _) => Center(child: Icon(FA.image, size: 40, color: t.textFaint)),
                                ),
                                if (scene.rect != null)
                                  IgnorePointer(child: CustomPaint(painter: _HighlightPainter(scene.rect!, _item + 1))),
                              ]),
                      ),
                    ),
                    const SizedBox(height: 10),
                    _Progress(
                      clock: _clock,
                      item: _item,
                      count: _itemCount,
                      paused: _paused,
                      onPause: _togglePause,
                      next: last ? null : (_page + 1 == 0 ? pages[0].title : '${_page + 1} · ${pages[_page + 1].title}'),
                    ),
                  ]),
                ),
                const SizedBox(width: 44),
                Expanded(
                  flex: 9,
                  child: ClipRect(child: _FitText(child: text)),
                ),
              ]),
              if (_arrow != null)
                Positioned.fill(child: IgnorePointer(child: CustomPaint(painter: _ArrowPainter(_arrow!.$1, _arrow!.$2, _arrow!.$3)))),
            ]),
          ),
        ]),
      ),
    );
  }
}

/// Progress of the page: one segment per item (done, filling, to come),
/// the item number, the step that comes after the last item, and pause.
class _Progress extends StatelessWidget {
  const _Progress({required this.clock, required this.item, required this.count, required this.paused, required this.onPause, this.next});
  final Animation<double> clock;
  final int item, count;
  final bool paused;
  final VoidCallback onPause;

  /// The next step (null on the last page).
  final String? next;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final l10n = context.l10n;
    return Row(children: [
      IconBtn(icon: paused ? FA.play : FA.pause, tooltip: paused ? l10n.tutorialPlay : l10n.tutorialPause, size: 28, bordered: true, onPressed: onPause),
      const SizedBox(width: 10),
      Expanded(
        child: AnimatedBuilder(
          animation: clock,
          builder: (context, _) => Row(children: [
            for (var i = 0; i < count; i++) ...[
              if (i > 0) const SizedBox(width: 4),
              Expanded(
                child: ClipRRect(
                  borderRadius: BorderRadius.circular(3),
                  child: LinearProgressIndicator(
                    value: i < item ? 1 : (i == item ? clock.value : 0),
                    minHeight: 6,
                    color: SplaceRed.c700,
                    backgroundColor: t.surface3,
                  ),
                ),
              ),
            ],
          ]),
        ),
      ),
      const SizedBox(width: 10),
      Text(l10n.tutorialItemOf('${item + 1}', '$count'), style: TextStyle(fontSize: 12, color: t.textMuted)),
      if (next != null && item == count - 1) ...[
        const SizedBox(width: 10),
        Icon(FA.arrowRight, size: 11, color: SplaceRed.c700),
        const SizedBox(width: 6),
        Flexible(
          child: Text(l10n.tutorialNextStep(next!),
              overflow: TextOverflow.ellipsis, style: const TextStyle(fontSize: 12, fontWeight: FontWeight.w600, color: SplaceRed.c700)),
        ),
      ],
    ]);
  }
}

/// The area of a scene: the rest of the picture dimmed, a red frame and the
/// item's number at its corner.
class _HighlightPainter extends CustomPainter {
  _HighlightPainter(this.rect, this.number);
  final Rect rect;
  final int number;

  @override
  void paint(Canvas canvas, Size size) {
    final fit = _fitted(size);
    final s = fit.width / _shotSize.width;
    final r = Rect.fromLTWH(fit.left + rect.left * s, fit.top + rect.top * s, rect.width * s, rect.height * s).inflate(3);
    final frame = RRect.fromRectAndRadius(r, const Radius.circular(8));
    canvas.drawPath(
      Path.combine(PathOperation.difference, Path()..addRect(fit), Path()..addRRect(frame)),
      Paint()..color = const Color(0x5A0F111E),
    );
    canvas.drawRRect(frame, Paint()
      ..color = SplaceRed.c700
      ..style = PaintingStyle.stroke
      ..strokeWidth = 3);
    final badge = Offset(r.left, r.top);
    canvas.drawCircle(badge, 11, Paint()..color = SplaceRed.c700);
    final tp = TextPainter(
      text: TextSpan(text: '$number', style: const TextStyle(fontFamily: 'Montserrat', fontSize: 12, fontWeight: FontWeight.w700, color: Colors.white)),
      textDirection: TextDirection.ltr,
    )..layout();
    tp.paint(canvas, badge - Offset(tp.width / 2, tp.height / 2));
  }

  @override
  bool shouldRepaint(_HighlightPainter old) => old.rect != rect || old.number != number;
}

/// A red arrow from the active item ([from]) to its area ([to]) in straight
/// lines with right angles: left to [gapX], up or down, left to the area.
class _ArrowPainter extends CustomPainter {
  _ArrowPainter(this.from, this.gapX, this.to);
  final Offset from, to;
  final double gapX;

  @override
  void paint(Canvas canvas, Size size) {
    const r = 8.0;
    final dy = to.dy - from.dy;
    final path = Path()..moveTo(from.dx, from.dy);
    if (dy.abs() < 2 * r) {
      path.lineTo(to.dx, from.dy);
    } else {
      final sy = dy.sign;
      path
        ..lineTo(gapX + r, from.dy)
        ..quadraticBezierTo(gapX, from.dy, gapX, from.dy + sy * r)
        ..lineTo(gapX, to.dy - sy * r)
        ..quadraticBezierTo(gapX, to.dy, gapX - r, to.dy)
        ..lineTo(to.dx, to.dy);
    }
    canvas.drawPath(path, Paint()
      ..color = SplaceRed.c700
      ..style = PaintingStyle.stroke
      ..strokeWidth = 3
      ..strokeCap = StrokeCap.round
      ..strokeJoin = StrokeJoin.round);
    canvas.drawCircle(from, 4, Paint()..color = SplaceRed.c700);
    // Head pointing left, into the area.
    final end = Offset(to.dx, dy.abs() < 2 * r ? from.dy : to.dy);
    final head = Path()
      ..moveTo(end.dx, end.dy)
      ..lineTo(end.dx + 13, end.dy - 7)
      ..lineTo(end.dx + 13, end.dy + 7)
      ..close();
    canvas.drawPath(head, Paint()..color = SplaceRed.c700);
  }

  @override
  bool shouldRepaint(_ArrowPainter old) => old.from != from || old.to != to || old.gapX != gapX;
}

/// The pages in [lang]; [titles]: the interface and the eight steps.
List<_Page> _pages(String lang, List<String> titles) {
  final text = switch (lang) { 'pt' => _pt, 'es' => _es, _ => _en };
  return [for (var i = 0; i < text.length; i++) _Page(titles[i], text[i].$1, text[i].$2, text[i].$3)];
}

typedef _Text = (String, List<String>, String);

// Markup of the texts: **name** in bold (the function), {icon} drawn as the
// icon of the app ([_icons]).

const _en = <_Text>[
  (
    'SPLACE runs the whole analysis on one page, step by step, with published programs bundled in the app. The status bar at the bottom shows the version and the CPU, memory, disk and network in use.',
    [
      '**Step navigator**: lists the nine steps and their parts (1a, 1b…); click one to go there. {check} marks what is done; locked steps wait for the data of the previous ones.',
      '**{key} NCBI API key**: a free key from your NCBI account makes downloads about three times faster (10 instead of 3 requests per second). Recommended whenever you download many records.',
      '**{tools} Tools**, **{sliders} Display** and **{log} Log**: the programs found and their versions (set a path if one is missing), theme and text size, and the record of everything SPLACE did.',
      '**{plus} New project**: starts a clean analysis in a new folder without closing the app — use it to analyse another group.',
      '**Confirm the new project**: the current results stay in their folder (shown here); settings, the API key and the tools are kept.',
    ],
    'Use ← and → or the step buttons to move through the tutorial; the {help} button of the top bar opens the welcome screen again.',
  ),
  (
    'Bring the sequences into the project. Sources can be combined: for example, genomes from NCBI plus your own outgroup in a GenBank file.',
    [
      '**Example buttons** of 1b, 1c and 1d: fill in accessions or load ready-made GenBank, loci or barcoding files — the quickest way to try every step; point at a button to see what its example holds.',
      '**1a Search by taxonomic name**: type a genus or a family and choose **Genomes** (complete or partial mitogenomes/plastomes) or **Individual markers** (COI, CYTB… when the group has no genomes; COI and the plant barcodes rbcL and matK also come from **BOLD Systems**, each record showing its source).',
      '**A species name** (two words) is recognized and shown in italics; turn on **synonyms (GBIF)** to also find records filed under older names.',
      '**Review before importing**: complete, partial and UNVERIFIED records are counted separately; leave UNVERIFIED out and click **Import**.',
      '**1b Accessions**, **1c GenBank files**, **1d FASTA files**: for known records, your own annotated genomes, and loci or barcodes from your lab.',
      '**{key} API key**: paste it once; it is stored only on this computer.',
    ],
    'Individual markers are searched with every synonym of the gene in SynGenes and joined by species, so a missing marker becomes missing data.',
  ),
  (
    'Check the records before the analysis: the species identification, the taxonomy and the table. For FASTA barcodes, this is where each specimen is identified.',
    [
      '**2a Species identification**: compares a barcode gene of each record with GenBank (BLAST) and BOLD — use it to catch misidentified records or to name a newly assembled genome.',
      '**Genes and {sliders} search parameters**: every gene goes to GenBank; COI to BOLD (rbcL and matK with its plant database). The parameters set the BLAST database, program, e-value, word size and number of results, and the BOLD database, minimum identity, overlap and hits; answers already received are reused.',
      '**2b Taxonomy**: GBIF (global names), NCBI Taxonomy (by the taxon id of the record) or ICMBio SALVE (Brazilian fauna, with the extinction-risk category).',
      '**Export corrected GenBank**: saves each record with the identified name, the taxonomy and the gene names standardized by SynGenes — ideal before submitting a new genome.',
      '**DNA barcodes**: with FASTA input, switch on "These sequences are DNA barcodes", identify every specimen, keep or change each name and export the **BOLD submission spreadsheet** (.xlsx).',
      '**{align} Alignments**: your sequence beside each GenBank and BOLD hit, with identity, coverage, mismatches, gaps, e-value and bit score — see what a name is based on before keeping it.',
      '**2c Records table**: search, edit or remove records; gene counts that differ from the rest are highlighted.',
    ],
    'Identification needs the internet and can take minutes; results are cached, so running again is fast.',
  ),
  (
    'Choose which annotations go into the analysis. SynGenes gives the same gene the same name in every record; choose the markers of the matrix after seeing which samples have each one and how long they are.',
    [
      '**PCGs, rRNAs, tRNAs**: switch on the kinds of genes to use — PCGs and rRNAs are the usual choice for phylogenies; tRNAs are short and add little signal.',
      '**SynGenes standardization**: names such as COX1, COI and cox1 become one; names it does not know are listed for you to assign.',
      '**FASTA input**: this step shows how the sample and the locus were read from each header.',
      '**Shortcuts**: Select All, None, **Shared Only** (markers present in every sample) and **Default to Phylogeny** (the usual set).',
      '**Marker switches**: each marker with the samples that have it and its mean length; turn on the ones to analyse.',
      '**Heatmaps**: presence shows the gaps of the dataset; length reveals truncated or unusually long sequences.',
    ],
    'A gene under two names would become two markers full of missing data: check the names before going on. Markers missing from a few samples can stay: they become missing data (?) in the supermatrix.',
  ),
  (
    'Align each marker with MAFFT and remove poorly aligned regions with trimAl or ClipKIT.',
    [
      '**4a Header Builder**: chooses the sequence names (accession, family, genus, species…), the same in every file.',
      '**4b Sequence type**: **Nucleotide** for rRNAs and non-coding loci, **Amino acid** for divergent proteins, **Codon-aware** for coding genes (proteins aligned, DNA back-translated, reading frame kept).',
      '**4c MAFFT**: the method (--auto by default), threads and options; **Run MAFFT Alignment** at the end of 4c aligns every marker in parallel.',
      '**4d Trimming**: trimAl or ClipKIT remove gappy, doubtful columns; **Keep the reading frame** removes whole codons.',
      '**Alignment results**: the viewer of each marker — raw, aligned, trimmed (removed columns in grey) and protein.',
    ],
    'Codon-aware alignment plus whole-codon trimming keeps the 1st, 2nd and 3rd positions in order for partitions and the saturation test.',
  ),
  (
    'Find out whether a codon position has accumulated so many changes that it no longer carries phylogenetic signal.',
    [
      '**5a Run the test**: PhyKIT measures saturation per gene and codon position on the IQ-TREE tree of each gene; the thresholds are editable.',
      '**5b Plots**: points that flatten below the diagonal mean saturation — usually the 3rd position.',
      '**5c Recommendation**: pooled over all genes (for example, remove the 3rd position), applied to the matrix with one click.',
    ],
    'A position is called saturated for the dataset when it is saturated in more than half of the genes.',
  ),
  (
    'Join the alignments into one supermatrix and decide how it is partitioned.',
    [
      '**6a Which alignments**: the trimmed ones (recommended) or the MAFFT output; gaps can become missing data.',
      '**6b Partitions**: by locus and codon position (recommended), by locus, by codon position or none; codon positions can be removed here.',
      '**6c Outgroup**: the taxa that root the tree.',
      '**6d Build the matrix**: writes NEXUS, PHYLIP and FASTA with the partition files of IQ-TREE and RAxML-NG.',
    ],
    'Samples missing a marker are filled with missing data, so they still enter the tree.',
  ),
  (
    'Infer the maximum-likelihood tree and, if you wish, delimit species.',
    [
      '**7a IQ-TREE 3**: model selection (MFP+MERGE), UFBoot2 and SH-aLRT support — the standard choice.',
      '**7b RAxML-NG**: a second, independent ML tree with bootstrap support; compare it with the IQ-TREE tree.',
      '**7c Run both**: both programs at once, the threads split between them.',
      '**7d Runs**: every run with its models, likelihood, command and files; **View tree** opens step 8.',
      '**7e Species delimitation**: mPTP and GMYC on the trees, ASAP and ABGD on the alignments — use them to propose species limits.',
    ],
    'Delimitations are hypotheses: compare the methods side by side on the tree and confirm them with other evidence.',
  ),
  (
    'Edit the figure, save everything and get the methods text with the references.',
    [
      '**Tree editor**: root, ladderize, colour clades, change labels and the legend — in the full window, every change can be undone.',
      '**Species delimitation on the tree**: one button per method and tree adds its bars beside the tips.',
      '**Compare and export**: two trees face to face; Newick, SVG or PNG at 300 dpi, sized in cm and pt.',
      '**8b Save everything**: the whole project in one ZIP.',
      '**8c How to cite**: the methods paragraph and the references, each saying what the tool did in your project.',
    ],
    'The figure keeps its real printed size on any monitor.',
  ),
];

const _pt = <_Text>[
  (
    'O SPLACE faz toda a análise em uma página, etapa por etapa, com programas publicados que vêm com o app. A barra de status, embaixo, mostra a versão e o uso de CPU, memória, disco e rede.',
    [
      '**Navegador de etapas**: lista as nove etapas e suas partes (1a, 1b…); clique para ir até elas. {check} marca o que foi feito; etapas bloqueadas esperam os dados das anteriores.',
      '**{key} Chave de API do NCBI**: uma chave gratuita da sua conta no NCBI deixa os downloads cerca de três vezes mais rápidos (10 em vez de 3 pedidos por segundo). Recomendada sempre que baixar muitos registros.',
      '**{tools} Ferramentas**, **{sliders} Exibição** e **{log} Registro**: os programas encontrados e suas versões (indique o caminho se faltar algum), tema e tamanho do texto, e o histórico de tudo o que o SPLACE fez.',
      '**{plus} Novo projeto**: começa uma análise limpa em outra pasta, sem fechar o app — use para analisar outro grupo.',
      '**Confirme o novo projeto**: os resultados atuais continuam na pasta deles (mostrada aqui); preferências, chave de API e ferramentas são mantidas.',
    ],
    'Use ← e → ou os botões das etapas para navegar pelo tutorial; o botão {help} da barra superior abre de novo a tela de boas-vindas.',
  ),
  (
    'Traga as sequências para o projeto. As fontes podem ser combinadas: por exemplo, genomas do NCBI mais o seu grupo externo num arquivo GenBank.',
    [
      '**Botões de exemplo** de 1b, 1c e 1d: preenchem acessos ou carregam arquivos prontos de GenBank, loci ou barcoding — o jeito mais rápido de testar todas as etapas; aponte para um botão para ver o que o exemplo contém.',
      '**1a Busca por nome taxonômico**: digite um gênero ou uma família e escolha **Genomas** (mitogenomas/plastomas completos ou parciais) ou **Marcadores individuais** (COI, CYTB… quando o grupo não tem genomas; o COI e os barcodes de plantas rbcL e matK vêm também do **BOLD Systems**, cada registro mostrando sua origem).',
      '**Um nome de espécie** (duas palavras) é reconhecido e aparece em itálico; ligue **sinônimos (GBIF)** para encontrar também registros com nomes antigos.',
      '**Revisão antes de importar**: registros completos, parciais e UNVERIFIED são contados à parte; deixe os UNVERIFIED de fora e clique em **Importar**.',
      '**1b Acessos**, **1c Arquivos GenBank**, **1d Arquivos FASTA**: para registros conhecidos, seus próprios genomas anotados, e loci ou barcodes do seu laboratório.',
      '**{key} Chave de API**: cole uma vez; ela fica guardada só neste computador.',
    ],
    'Os marcadores individuais são buscados com todos os sinônimos do gene no SynGenes e reunidos por espécie; um marcador ausente vira dado faltante.',
  ),
  (
    'Confira os registros antes da análise: a identificação das espécies, a taxonomia e a tabela. Para barcodes em FASTA, é aqui que cada espécime é identificado.',
    [
      '**2a Identificação de espécies**: compara um gene de barcode de cada registro com o GenBank (BLAST) e o BOLD — use para achar registros mal identificados ou dar nome a um genoma recém-montado.',
      '**Genes e {sliders} parâmetros da busca**: todo gene vai para o GenBank; o COI para o BOLD (rbcL e matK com o banco de plantas). Os parâmetros definem banco, programa, e-value, tamanho da palavra e número de resultados do BLAST, e banco, identidade mínima, sobreposição e resultados do BOLD; respostas já recebidas são reaproveitadas.',
      '**2b Taxonomia**: GBIF (nomes globais), NCBI Taxonomy (pelo taxon id do registro) ou SALVE/ICMBio (fauna brasileira, com a categoria de ameaça).',
      '**Exportar GenBank corrigido**: salva cada registro com o nome identificado, a taxonomia e os nomes dos genes padronizados pelo SynGenes — ideal antes de submeter um genoma novo.',
      '**DNA barcodes**: com entrada FASTA, ligue "Estas sequências são DNA barcodes", identifique cada espécime, mantenha ou troque cada nome e exporte a **planilha de submissão do BOLD** (.xlsx).',
      '**{align} Alinhamentos**: a sua sequência ao lado de cada resultado do GenBank e do BOLD, com identidade, cobertura, diferenças, gaps, e-value e bit score — veja em que um nome se apoia antes de mantê-lo.',
      '**2c Tabela de registros**: busque, edite ou remova registros; contagens de genes diferentes do resto ficam destacadas.',
    ],
    'A identificação precisa de internet e pode levar minutos; os resultados ficam em cache, então rodar de novo é rápido.',
  ),
  (
    'Escolha quais anotações entram na análise. O SynGenes dá o mesmo nome ao mesmo gene em todos os registros; escolha os marcadores da matriz depois de ver quais amostras têm cada um e qual o seu tamanho.',
    [
      '**PCGs, rRNAs, tRNAs**: ligue os tipos de genes a usar — PCGs e rRNAs são a escolha usual para filogenias; os tRNAs são curtos e trazem pouco sinal.',
      '**Padronização pelo SynGenes**: nomes como COX1, COI e cox1 viram um só; os nomes que ele não conhece aparecem para você associar.',
      '**Entrada FASTA**: esta etapa mostra como a amostra e o locus foram lidos de cada cabeçalho.',
      '**Atalhos**: Selecionar todos, Nenhum, **Só compartilhados** (presentes em todas as amostras) e **Padrão para filogenia** (o conjunto usual).',
      '**Chaves dos marcadores**: cada marcador com as amostras que o têm e o tamanho médio; ligue os que vai analisar.',
      '**Heatmaps**: o de presença mostra as lacunas do conjunto; o de tamanho revela sequências truncadas ou longas demais.',
    ],
    'Um gene com dois nomes viraria dois marcadores cheios de dados faltantes: confira os nomes antes de seguir. Marcadores ausentes em poucas amostras podem ficar: viram dados faltantes (?) na supermatriz.',
  ),
  (
    'Alinhe cada marcador com o MAFFT e remova as regiões mal alinhadas com o trimAl ou o ClipKIT.',
    [
      '**4a Construtor de cabeçalho**: define os nomes das sequências (acesso, família, gênero, espécie…), iguais em todos os arquivos.',
      '**4b Tipo de sequência**: **Nucleotídeo** para rRNAs e loci não codificantes, **Aminoácido** para proteínas divergentes, **Preservando códons** para genes codificadores (proteínas alinhadas, DNA retrotraduzido, fase de leitura mantida).',
      '**4c MAFFT**: o método (--auto por padrão), threads e opções; **Rodar alinhamento MAFFT**, no fim da 4c, alinha todos os marcadores em paralelo.',
      '**4d Trimagem**: o trimAl ou o ClipKIT removem colunas duvidosas e cheias de gaps; **Manter a fase de leitura** remove códons inteiros.',
      '**Resultados do alinhamento**: o visualizador de cada marcador — bruto, alinhado, trimado (colunas removidas em cinza) e proteína.',
    ],
    'Alinhamento preservando códons mais trimagem de códons inteiros mantém as posições 1, 2 e 3 em ordem para as partições e o teste de saturação.',
  ),
  (
    'Descubra se uma posição do códon acumulou tantas mudanças que já não traz sinal filogenético.',
    [
      '**5a Rodar o teste**: o PhyKIT mede a saturação por gene e posição do códon na árvore do IQ-TREE de cada gene; os limiares podem ser editados.',
      '**5b Gráficos**: pontos que achatam abaixo da diagonal indicam saturação — em geral na 3ª posição.',
      '**5c Recomendação**: reunida de todos os genes (por exemplo, remover a 3ª posição), aplicada à matriz com um clique.',
    ],
    'Uma posição é considerada saturada no conjunto quando está saturada em mais da metade dos genes.',
  ),
  (
    'Una os alinhamentos numa supermatriz e decida como ela é particionada.',
    [
      '**6a Quais alinhamentos**: os trimados (recomendado) ou a saída do MAFFT; os gaps podem virar dados faltantes.',
      '**6b Partições**: por locus e posição do códon (recomendado), por locus, por posição do códon ou nenhuma; posições do códon podem ser removidas aqui.',
      '**6c Grupo externo**: os táxons que enraízam a árvore.',
      '**6d Construir a matriz**: grava NEXUS, PHYLIP e FASTA com os arquivos de partição do IQ-TREE e do RAxML-NG.',
    ],
    'Amostras sem um marcador recebem dados faltantes, então continuam entrando na árvore.',
  ),
  (
    'Infira a árvore de máxima verossimilhança e, se quiser, delimite espécies.',
    [
      '**7a IQ-TREE 3**: escolha de modelo (MFP+MERGE), suporte UFBoot2 e SH-aLRT — a escolha padrão.',
      '**7b RAxML-NG**: uma segunda árvore de ML, independente, com bootstrap; compare-a com a do IQ-TREE.',
      '**7c Rodar os dois**: os dois programas ao mesmo tempo, com as threads divididas.',
      '**7d Execuções**: cada execução com modelos, verossimilhança, comando e arquivos; **Ver árvore** abre a etapa 8.',
      '**7e Delimitação de espécies**: mPTP e GMYC nas árvores, ASAP e ABGD nos alinhamentos — use para propor limites de espécies.',
    ],
    'Delimitações são hipóteses: compare os métodos lado a lado na árvore e confirme-as com outras evidências.',
  ),
  (
    'Edite a figura, salve tudo e obtenha o texto de métodos com as referências.',
    [
      '**Editor da árvore**: enraíze, ordene, colora clados, mude rótulos e a legenda — em tela cheia, e toda mudança pode ser desfeita.',
      '**Delimitação de espécies na árvore**: um botão por método e árvore acrescenta as barras ao lado dos terminais.',
      '**Comparar e exportar**: duas árvores frente a frente; Newick, SVG ou PNG a 300 dpi, com tamanhos em cm e pt.',
      '**8b Salvar tudo**: o projeto inteiro em um ZIP.',
      '**8c Como citar**: o parágrafo de métodos e as referências, cada uma dizendo o que a ferramenta fez no seu projeto.',
    ],
    'A figura mantém o tamanho real de impressão em qualquer monitor.',
  ),
];

const _es = <_Text>[
  (
    'SPLACE hace todo el análisis en una página, paso a paso, con programas publicados que vienen con la app. La barra de estado, abajo, muestra la versión y el uso de CPU, memoria, disco y red.',
    [
      '**Navegador de pasos**: lista los nueve pasos y sus partes (1a, 1b…); haga clic para ir a ellos. {check} marca lo hecho; los pasos bloqueados esperan los datos de los anteriores.',
      '**{key} Clave de API del NCBI**: una clave gratuita de su cuenta del NCBI hace las descargas unas tres veces más rápidas (10 en vez de 3 pedidos por segundo). Recomendada siempre que descargue muchos registros.',
      '**{tools} Herramientas**, **{sliders} Visualización** y **{log} Registro**: los programas encontrados y sus versiones (indique la ruta si falta alguno), tema y tamaño del texto, y el historial de todo lo que hizo SPLACE.',
      '**{plus} Proyecto nuevo**: empieza un análisis limpio en otra carpeta, sin cerrar la app — úselo para analizar otro grupo.',
      '**Confirme el proyecto nuevo**: los resultados actuales quedan en su carpeta (mostrada aquí); se conservan las preferencias, la clave de API y las herramientas.',
    ],
    'Use ← y → o los botones de los pasos para recorrer el tutorial; el botón {help} de la barra superior abre de nuevo la pantalla de bienvenida.',
  ),
  (
    'Traiga las secuencias al proyecto. Las fuentes se pueden combinar: por ejemplo, genomas del NCBI más su grupo externo en un archivo GenBank.',
    [
      '**Botones de ejemplo** de 1b, 1c y 1d: rellenan accesos o cargan archivos listos de GenBank, loci o barcoding — la forma más rápida de probar todos los pasos; apunte a un botón para ver qué contiene el ejemplo.',
      '**1a Búsqueda por nombre taxonómico**: escriba un género o una familia y elija **Genomas** (mitogenomas/plastomas completos o parciales) o **Marcadores individuales** (COI, CYTB… cuando el grupo no tiene genomas; el COI y los barcodes de plantas rbcL y matK vienen también de **BOLD Systems**, cada registro mostrando su origen).',
      '**Un nombre de especie** (dos palabras) se reconoce y aparece en cursiva; active **sinónimos (GBIF)** para encontrar también registros con nombres antiguos.',
      '**Revisión antes de importar**: los registros completos, parciales y UNVERIFIED se cuentan aparte; deje fuera los UNVERIFIED y haga clic en **Importar**.',
      '**1b Accesos**, **1c Archivos GenBank**, **1d Archivos FASTA**: para registros conocidos, sus propios genomas anotados, y loci o barcodes de su laboratorio.',
      '**{key} Clave de API**: péguela una vez; queda guardada solo en este equipo.',
    ],
    'Los marcadores individuales se buscan con todos los sinónimos del gen en SynGenes y se reúnen por especie; un marcador ausente pasa a ser dato faltante.',
  ),
  (
    'Revise los registros antes del análisis: la identificación de las especies, la taxonomía y la tabla. Para barcodes en FASTA, aquí se identifica cada espécimen.',
    [
      '**2a Identificación de especies**: compara un gen de barcode de cada registro con GenBank (BLAST) y BOLD — úsela para encontrar registros mal identificados o nombrar un genoma recién ensamblado.',
      '**Genes y {sliders} parámetros de la búsqueda**: todo gen va a GenBank; el COI a BOLD (rbcL y matK con la base de plantas). Los parámetros definen base, programa, e-value, tamaño de palabra y número de resultados de BLAST, y base, identidad mínima, solapamiento y resultados de BOLD; las respuestas ya recibidas se reutilizan.',
      '**2b Taxonomía**: GBIF (nombres globales), NCBI Taxonomy (por el taxon id del registro) o SALVE/ICMBio (fauna brasileña, con la categoría de amenaza).',
      '**Exportar GenBank corregido**: guarda cada registro con el nombre identificado, la taxonomía y los nombres de genes estandarizados por SynGenes — ideal antes de enviar un genoma nuevo.',
      '**DNA barcodes**: con entrada FASTA, active "Estas secuencias son DNA barcodes", identifique cada espécimen, conserve o cambie cada nombre y exporte la **planilla de envío de BOLD** (.xlsx).',
      '**{align} Alineamientos**: su secuencia junto a cada resultado de GenBank y BOLD, con identidad, cobertura, diferencias, gaps, e-value y bit score — vea en qué se apoya un nombre antes de conservarlo.',
      '**2c Tabla de registros**: busque, edite o elimine registros; los conteos de genes distintos del resto se destacan.',
    ],
    'La identificación necesita internet y puede tardar minutos; los resultados quedan en caché, así que repetirla es rápido.',
  ),
  (
    'Elija qué anotaciones entran en el análisis. SynGenes da el mismo nombre al mismo gen en todos los registros; elija los marcadores de la matriz después de ver qué muestras tienen cada uno y su longitud.',
    [
      '**PCGs, rRNAs, tRNAs**: active los tipos de genes a usar — PCGs y rRNAs son la elección habitual para filogenias; los tRNAs son cortos y aportan poca señal.',
      '**Estandarización por SynGenes**: nombres como COX1, COI y cox1 pasan a ser uno; los nombres que no conoce aparecen para que los asigne.',
      '**Entrada FASTA**: este paso muestra cómo se leyeron la muestra y el locus de cada encabezado.',
      '**Atajos**: Seleccionar todos, Ninguno, **Solo compartidos** (presentes en todas las muestras) y **Predeterminado para filogenia** (el conjunto habitual).',
      '**Interruptores de marcadores**: cada marcador con las muestras que lo tienen y su longitud media; active los que va a analizar.',
      '**Heatmaps**: el de presencia muestra los huecos del conjunto; el de longitud revela secuencias truncadas o demasiado largas.',
    ],
    'Un gen con dos nombres se convertiría en dos marcadores llenos de datos faltantes: revise los nombres antes de seguir. Los marcadores ausentes en pocas muestras pueden quedarse: pasan a ser datos faltantes (?) en la supermatriz.',
  ),
  (
    'Alinee cada marcador con MAFFT y elimine las regiones mal alineadas con trimAl o ClipKIT.',
    [
      '**4a Constructor de encabezados**: define los nombres de las secuencias (acceso, familia, género, especie…), iguales en todos los archivos.',
      '**4b Tipo de secuencia**: **Nucleótido** para rRNAs y loci no codificantes, **Aminoácido** para proteínas divergentes, **Preservando codones** para genes codificantes (proteínas alineadas, ADN retrotraducido, marco de lectura conservado).',
      '**4c MAFFT**: el método (--auto por defecto), hilos y opciones; **Ejecutar alineamiento MAFFT**, al final de 4c, alinea todos los marcadores en paralelo.',
      '**4d Recorte**: trimAl o ClipKIT eliminan columnas dudosas y llenas de gaps; **Mantener el marco de lectura** elimina codones enteros.',
      '**Resultados del alineamiento**: el visor de cada marcador — bruto, alineado, recortado (columnas eliminadas en gris) y proteína.',
    ],
    'El alineamiento preservando codones más el recorte de codones enteros mantiene las posiciones 1, 2 y 3 en orden para las particiones y la prueba de saturación.',
  ),
  (
    'Averigüe si una posición del codón acumuló tantos cambios que ya no aporta señal filogenética.',
    [
      '**5a Ejecutar la prueba**: PhyKIT mide la saturación por gen y posición del codón en el árbol de IQ-TREE de cada gen; los umbrales se pueden editar.',
      '**5b Gráficos**: los puntos que se aplanan bajo la diagonal indican saturación — en general en la 3ª posición.',
      '**5c Recomendación**: reunida de todos los genes (por ejemplo, quitar la 3ª posición), aplicada a la matriz con un clic.',
    ],
    'Una posición se considera saturada en el conjunto cuando lo está en más de la mitad de los genes.',
  ),
  (
    'Una los alineamientos en una supermatriz y decida cómo se particiona.',
    [
      '**6a Qué alineamientos**: los recortados (recomendado) o la salida de MAFFT; los gaps pueden pasar a datos faltantes.',
      '**6b Particiones**: por locus y posición del codón (recomendado), por locus, por posición del codón o ninguna; aquí se pueden quitar posiciones del codón.',
      '**6c Grupo externo**: los taxones que enraízan el árbol.',
      '**6d Construir la matriz**: guarda NEXUS, PHYLIP y FASTA con los archivos de particiones de IQ-TREE y RAxML-NG.',
    ],
    'Las muestras sin un marcador reciben datos faltantes, así que siguen entrando en el árbol.',
  ),
  (
    'Infiera el árbol de máxima verosimilitud y, si quiere, delimite especies.',
    [
      '**7a IQ-TREE 3**: elección de modelo (MFP+MERGE), soporte UFBoot2 y SH-aLRT — la elección estándar.',
      '**7b RAxML-NG**: un segundo árbol de ML, independiente, con bootstrap; compárelo con el de IQ-TREE.',
      '**7c Ejecutar ambos**: los dos programas a la vez, con los hilos repartidos.',
      '**7d Ejecuciones**: cada ejecución con modelos, verosimilitud, comando y archivos; **Ver árbol** abre el paso 8.',
      '**7e Delimitación de especies**: mPTP y GMYC en los árboles, ASAP y ABGD en los alineamientos — úselos para proponer límites de especies.',
    ],
    'Las delimitaciones son hipótesis: compare los métodos lado a lado en el árbol y confírmelas con otras evidencias.',
  ),
  (
    'Edite la figura, guarde todo y obtenga el texto de métodos con las referencias.',
    [
      '**Editor del árbol**: enraíce, ordene, coloree clados, cambie etiquetas y la leyenda — a pantalla completa, y todo cambio se puede deshacer.',
      '**Delimitación de especies en el árbol**: un botón por método y árbol añade sus barras junto a los terminales.',
      '**Comparar y exportar**: dos árboles frente a frente; Newick, SVG o PNG a 300 dpi, con tamaños en cm y pt.',
      '**8b Guardar todo**: todo el proyecto en un ZIP.',
      '**8c Cómo citar**: el párrafo de métodos y las referencias, cada una diciendo lo que hizo la herramienta en su proyecto.',
    ],
    'La figura conserva su tamaño real de impresión en cualquier monitor.',
  ),
];

/// Icons of the app named in the texts as {name}.
const _icons = <String, IconData>{
  'key': FA.key,
  'tools': FA.screwdriverWrench,
  'sliders': FA.sliders,
  'log': FA.rectangleList,
  'help': FA.circleQuestion,
  'plus': FA.plus,
  'check': FA.check,
  'github': FA.github,
  'menu': FA.bars,
  'align': FA.alignLeft,
};

/// [text] as spans: **bold** in [bold], {icon} as the app's icon.
List<InlineSpan> _rich(String text, TextStyle bold, Color iconColor, double size) {
  final out = <InlineSpan>[];
  var start = 0;
  for (final m in RegExp(r'\*\*(.+?)\*\*|\{(\w+)\}').allMatches(text)) {
    if (m.start > start) out.add(TextSpan(text: text.substring(start, m.start)));
    if (m.group(1) != null) {
      out.add(TextSpan(style: bold, children: _rich(m.group(1)!, bold, iconColor, size)));
    } else {
      final icon = _icons[m.group(2)];
      out.add(icon == null
          ? TextSpan(text: m.group(0))
          : WidgetSpan(
              alignment: PlaceholderAlignment.middle,
              child: Padding(padding: const EdgeInsets.symmetric(horizontal: 2), child: Icon(icon, size: size, color: iconColor)),
            ));
    }
    start = m.end;
  }
  if (start < text.length) out.add(TextSpan(text: text.substring(start)));
  return out;
}

/// Fits [child] (the text of a page) in the space it has without scrolling:
/// when it is too tall, it is laid out wider and scaled down, so the lines
/// still fill the whole width instead of leaving an empty band.
class _FitText extends SingleChildRenderObjectWidget {
  const _FitText({required Widget super.child});

  @override
  RenderObject createRenderObject(BuildContext context) => _RenderFitText();
}

class _RenderFitText extends RenderProxyBox {
  double _scale = 1;

  @override
  void performLayout() {
    final w = constraints.maxWidth, h = constraints.maxHeight;
    final c = child!;
    _scale = 1;
    // Down to 55 %, in 2.5 % steps: the first scale whose height fits.
    for (var s = 1.0; s >= 0.55; s -= 0.025) {
      _scale = s;
      c.layout(BoxConstraints.tightFor(width: w / s), parentUsesSize: true);
      if (c.size.height * s <= h) break;
    }
    size = constraints.constrain(Size(w, math.min(h, c.size.height * _scale)));
  }

  @override
  void applyPaintTransform(RenderBox child, Matrix4 transform) => transform.scaleByDouble(_scale, _scale, 1, 1);

  @override
  void paint(PaintingContext context, Offset offset) {
    layer = context.pushTransform(needsCompositing, offset, Matrix4.diagonal3Values(_scale, _scale, 1), (ctx, o) => ctx.paintChild(child!, o),
        oldLayer: layer as TransformLayer?);
  }

  @override
  bool hitTestChildren(BoxHitTestResult result, {required Offset position}) => result.addWithPaintTransform(
        transform: Matrix4.diagonal3Values(_scale, _scale, 1),
        position: position,
        hitTest: (r, p) => child!.hitTest(r, position: p),
      );

  @override
  double computeMinIntrinsicHeight(double width) => 0;
  @override
  double computeMaxIntrinsicHeight(double width) => 0;
}
