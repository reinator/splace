import 'dart:async';

import 'package:flutter/material.dart';

import '../../../l10n/app_localizations.dart';
import '../../core/saturation.dart';
import '../../core/tree_progress.dart';
import '../../state/project_controller.dart';
import '../common.dart';
import '../fa.dart';
import '../theme.dart';
import '../widgets/dialogs.dart';
import '../widgets/kit.dart';
import '../widgets/splace_tooltip.dart';

/// `#progressModal` of the HTML app: shown over the page while SPLACE talks
/// to NCBI or the taxonomy services, with the request rate, the API key
/// status, the current item and a progress bar.
class ProgressModal extends StatefulWidget {
  const ProgressModal({super.key, required this.controller});
  final ProjectController controller;

  @override
  State<ProgressModal> createState() => _ProgressModalState();
}

class _ProgressModalState extends State<ProgressModal> {
  /// Refreshes the elapsed time of long runs.
  Timer? _clock;

  @override
  void initState() {
    super.initState();
    _clock = Timer.periodic(const Duration(seconds: 1), (_) {
      final m = widget.controller.modal;
      if (mounted && m != null && !m.finished) setState(() {});
    });
  }

  @override
  void dispose() {
    _clock?.cancel();
    super.dispose();
  }

  @override
  Widget build(BuildContext context) {
    final c = widget.controller;
    final m = c.modal!;
    final t = context.tokens;
    final l10n = context.l10n;
    final ncbi = m.kind == ModalKind.ncbiSearch || m.kind == ModalKind.ncbiFetch || m.subject == 'NCBI';
    final tree = m.tree;
    final fraction = tree != null
        ? tree.fraction
        : m.trees.isNotEmpty
            ? m.trees.values.map((p) => p.fraction).fold(0.0, (a, b) => a + b) / m.trees.length
            // One bar for the whole run: every marker's own progress
            // (MAFFT's stages, validation, back-translation) adds to it.
            : m.items.isNotEmpty
                ? _itemsFraction(m.items)
                : (m.total == 0 ? null : (m.done / m.total).clamp(0.0, 1.0));

    final (icon, title, subtitle) = switch (m.kind) {
      ModalKind.ncbiSearch => (FA.magnifyingGlass, l10n.progressSearchTitle, l10n.progressSearchingFor(m.subject)),
      ModalKind.ncbiFetch => (FA.download, l10n.progressFetchTitle, l10n.progressCounter('${m.done}', '${m.total}')),
      ModalKind.taxonomy => (FA.sitemap, l10n.progressTaxonomyTitle(m.subject), l10n.progressSpecies('${m.done}', '${m.total}')),
      ModalKind.speciesId => (FA.dna, l10n.speciesIdRunning, l10n.speciesIdGenesUsed(m.subject)),
      ModalKind.alignment => (FA.alignLeft, l10n.alignRunning, l10n.markersProgress('${m.done}', '${m.total}', _modeName(l10n, m.subject))),
      ModalKind.trimming => (FA.scissors, l10n.trimRunning, l10n.markersProgress('${m.done}', '${m.total}', _modeName(l10n, m.subject))),
      ModalKind.saturation => (FA.chartLine, l10n.saturationRunning, l10n.saturationProgress('${m.done}', '${m.total}')),
      ModalKind.matrix => (FA.tableCells, l10n.matrixRunning, m.current ?? l10n.matrixRunningBody),
      ModalKind.iqtree => (FA.diagramProject, l10n.iqtreeRunning, l10n.treeRunSubtitle(m.subject, _elapsed(m))),
      ModalKind.raxml => (FA.diagramProject, l10n.raxmlRunning, l10n.treeRunSubtitle(m.subject, _elapsed(m))),
      ModalKind.bothTrees => (FA.diagramProject, l10n.bothRunning, l10n.bothRunningSub(m.subject, _elapsed(m))),
      ModalKind.delimitation => (FA.users, l10n.delimRunning, m.subject),
    };
    final tall = tree != null || m.items.isNotEmpty;
    final finishedIcon = m.failed ? FA.circleXmark : FA.circleCheck;
    final finishedColor = m.failed ? SplaceTokens.danger : t.accent;

    return Stack(children: [
      const ModalBarrier(color: Color(0x660F172A), dismissible: false),
      Center(
        child: Material(
          color: Colors.transparent,
          child: Container(
            width: SplaceDialog.widthFor(MediaQuery.sizeOf(context)),
            // Runs with one card per marker (MAFFT, trimming, saturation,
            // trees) take the size of the other modals.
            height: tall ? SplaceDialog.heightFor(MediaQuery.sizeOf(context)) : null,
            decoration: BoxDecoration(color: t.surface, borderRadius: BorderRadius.circular(18), boxShadow: shadowLg, border: Border.all(color: t.border)),
            clipBehavior: Clip.antiAlias,
            child: Column(mainAxisSize: MainAxisSize.min, crossAxisAlignment: CrossAxisAlignment.stretch, children: [
              DialogHeader(
                leading: m.finished
                    ? DialogHeader.iconTile(context, finishedIcon, color: finishedColor)
                    : Container(
                        width: 40,
                        height: 40,
                        padding: const EdgeInsets.all(11),
                        decoration: BoxDecoration(
                          color: t.accentSoft,
                          borderRadius: BorderRadius.circular(10),
                          border: Border.all(color: t.accentSoftBorder),
                        ),
                        child: CircularProgressIndicator(strokeWidth: 2.4, color: t.accent),
                      ),
                title: title,
                subtitle: subtitle,
              ),
              _grow(
                tall,
                Padding(
                padding: const EdgeInsets.all(24),
                child: Column(mainAxisSize: tall ? MainAxisSize.max : MainAxisSize.min, crossAxisAlignment: CrossAxisAlignment.stretch, children: [
                  // GenBank downloads name the species (italic) and its
                  // voucher, as in the HTML app, instead of the request id.
                  if (m.kind == ModalKind.ncbiFetch && m.species != null)
                    Text.rich(
                      TextSpan(style: TextStyle(fontSize: 14, color: t.text), children: [
                        TextSpan(text: m.species, style: const TextStyle(fontStyle: FontStyle.italic, fontWeight: FontWeight.w600)),
                        if (m.voucher != null) TextSpan(text: ' (${m.voucher})', style: TextStyle(color: t.textMuted)),
                      ]),
                      maxLines: 2,
                      overflow: TextOverflow.ellipsis,
                    )
                  else if (m.trees.isNotEmpty)
                    // Both programs: one row each, as in their own modal.
                    Expanded(
                      child: SingleChildScrollView(
                        child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
                          for (final (i, e) in m.trees.entries.indexed) ...[
                            if (i > 0) const SizedBox(height: 14),
                            _ProgramRow(name: e.key, progress: e.value, tail: m.tails[e.key] ?? const [], item: m.items[e.key]),
                          ],
                        ]),
                      ),
                    )
                  else if (tree != null) ...[
                    _TreeStages(progress: tree, finished: m.finished && !m.failed),
                    const SizedBox(height: 12),
                    _treeDetail(context, tree),
                    if (m.tail.isNotEmpty) ...[
                      const SizedBox(height: 12),
                      _Tail(lines: m.tail),
                    ],
                    if (m.items.isNotEmpty) ...[
                      const SizedBox(height: 14),
                      Kicker(l10n.geneTreesProgress, size: 11),
                      const SizedBox(height: 8),
                      Expanded(child: _ItemList(items: m.items)),
                    ] else
                      const Spacer(),
                  ] else if (m.items.isNotEmpty)
                    Expanded(child: _ItemList(items: m.items, noteFor: m.kind == ModalKind.saturation ? (name) => _saturationNote(context, c, name) : null))
                  else
                    Text(
                      switch (m.kind) {
                        ModalKind.ncbiFetch => l10n.progressStatusIdle,
                        // Markers running now.
                        ModalKind.alignment || ModalKind.trimming =>
                          (m.current ?? '').isEmpty ? (m.finished ? '' : l10n.progressStatusIdle) : l10n.alignCurrent(m.current!),
                        _ => m.current ?? l10n.progressStatusIdle,
                      },
                      maxLines: 2,
                      overflow: TextOverflow.ellipsis,
                      style: TextStyle(fontSize: 14, color: t.text, fontStyle: m.kind == ModalKind.taxonomy ? FontStyle.italic : FontStyle.normal),
                    ),
                  if (ncbi) ...[
                    const SizedBox(height: 4),
                    Text(
                      '${c.ncbi.hasKey ? l10n.progressRateWithKey : l10n.progressRateDefault} · ${c.ncbi.hasKey ? l10n.progressApiEnabled : l10n.progressApiDisabled}',
                      style: TextStyle(fontSize: 13, color: t.textMuted),
                    ),
                  ],
                  const SizedBox(height: 14),
                  RainbowProgressBar(value: m.finished ? 1 : fraction),
                  const SizedBox(height: 6),
                  Row(children: [
                    Expanded(
                      child: Text(
                        // Search and identification report stages, not items.
                        m.result ??
                            (fraction == null || m.kind == ModalKind.speciesId || m.kind == ModalKind.ncbiSearch || tree != null || m.trees.isNotEmpty || m.kind == ModalKind.matrix
                                ? l10n.progressWorking
                                : l10n.progressCounter('${m.done}', '${m.total}')),
                        maxLines: 2,
                        overflow: TextOverflow.ellipsis,
                        style: TextStyle(fontSize: 13, color: m.finished ? finishedColor : t.textMuted),
                      ),
                    ),
                    Text(fraction == null ? '' : '${(fraction * 100).round()}%',
                        style: TextStyle(fontSize: 13, color: t.textMuted)),
                  ]),
                ]),
                ),
              ),
              Container(
                padding: const EdgeInsets.fromLTRB(20, 12, 14, 12),
                decoration: BoxDecoration(color: t.surface2, border: Border(top: BorderSide(color: t.border))),
                child: Row(children: [
                  Expanded(
                    child: Text(
                        switch (m.kind) {
                          ModalKind.taxonomy => l10n.progressTaxonomyFooter(m.subject),
                          ModalKind.speciesId => l10n.speciesIdFooter,
                          ModalKind.alignment => l10n.alignFooter,
                          ModalKind.trimming => l10n.trimFooter,
                          ModalKind.saturation || ModalKind.matrix || ModalKind.delimitation => l10n.localFooter,
                          ModalKind.iqtree || ModalKind.raxml || ModalKind.bothTrees => l10n.treeFooter,
                          _ => l10n.progressFooter,
                        },
                        style: TextStyle(fontSize: 12, color: t.textFaint)),
                  ),
                  const SizedBox(width: 10),
                  if (m.finished)
                    FilledButton(onPressed: c.closeModal, child: Text(l10n.close))
                  else
                    DialogCancelButton(onPressed: c.cancelTask),
                ]),
              ),
            ]),
          ),
        ),
      ),
    ]);
  }
}

/// Mean progress of the markers: finished ones count whole, running ones
/// by their stage, queued ones as zero.
double _itemsFraction(Map<String, ModalItem> items) {
  var sum = 0.0;
  for (final it in items.values) {
    sum += switch (it.state) {
      ItemState.queued => 0,
      ItemState.running => it.progress.clamp(0.0, 1.0),
      _ => 1,
    };
  }
  return sum / items.length;
}

/// [child] stretched over the free height of a tall modal.
Widget _grow(bool tall, Widget child) => tall ? Expanded(child: child) : child;

String _elapsed(TaskModal m) {
  final s = DateTime.now().difference(m.started).inSeconds;
  return s < 3600 ? '${s ~/ 60}:${(s % 60).toString().padLeft(2, '0')}' : '${s ~/ 3600}h${((s % 3600) ~/ 60).toString().padLeft(2, '0')}';
}

/// What the program is doing now, in words (models tested, iteration…).
Widget _treeDetail(BuildContext context, TreeProgress p) {
  final t = context.tokens;
  final l10n = context.l10n;
  String n(num? v) => v == null ? '?' : '$v';
  final lines = <String>[
    switch (p.stage) {
      TreeStage.prepare => l10n.treeDetailPrepare,
      TreeStage.modelFinder when p.moose => l10n.treeDetailMoose(n(p.modelsTotal)),
      TreeStage.modelFinder => p.modelsTotal == null ? l10n.treeDetailModelStart : l10n.treeDetailModels(n(p.modelsDone ?? 0), n(p.modelsTotal)),
      TreeStage.search => p.startsTotal != null
          ? l10n.treeDetailStarts(n(p.startsDone ?? 0), n(p.startsTotal))
          : p.iteration == null
              ? l10n.treeDetailSearchStart
              : l10n.treeDetailIteration(n(p.iteration)) +
                  (p.secondsLeft == null ? '' : ' · ${l10n.treeDetailLeft(_duration(p.secondsLeft!))}'),
      TreeStage.support => p.bootTotal != null && p.bootDone != null
          ? l10n.treeDetailBoot(n(p.bootDone), n(p.bootTotal))
          : l10n.treeDetailSupport,
      TreeStage.finish => l10n.treeDetailFinish,
    },
    if (p.bestModel != null) l10n.treeDetailBestModel(p.bestModel!),
    if ((p.bestLogL ?? p.logL) != null) 'lnL ${(p.bestLogL ?? p.logL)!.toStringAsFixed(3)}',
    if (p.bootTotal != null && p.stage == TreeStage.search && p.startsTotal == null) l10n.treeDetailUfboot(n(p.bootTotal)),
    if (p.bootConverged) l10n.treeDetailConverged(n(p.bootDone)),
  ];
  return Wrap(spacing: 8, runSpacing: 6, crossAxisAlignment: WrapCrossAlignment.center, children: [
    Text(lines.first, style: TextStyle(fontSize: 14, fontWeight: FontWeight.w600, color: t.text)),
    for (final l in lines.skip(1)) TagChip(l, color: t.accent),
  ]);
}

String _duration(int s) => s < 60 ? '${s}s' : (s < 3600 ? '${s ~/ 60} min ${s % 60}s' : '${s ~/ 3600} h ${(s % 3600) ~/ 60} min');

/// Saturation verdicts of a finished locus, in the user's language.
String? _saturationNote(BuildContext context, ProjectController c, String name) {
  final r = c.loci[name]?.saturation;
  if (r == null) return null;
  return r.subsets.entries.map((e) => '${e.key.label}: ${verdictLabel(context, e.value.verdict)}').join(' · ');
}

/// The stages of a tree inference, done / running / to do.
class _TreeStages extends StatelessWidget {
  const _TreeStages({required this.progress, required this.finished, this.failed = false});
  final TreeProgress progress;
  final bool finished;

  /// The program stopped at the current stage.
  final bool failed;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final l10n = context.l10n;
    final stages = progress.stages;
    final current = finished ? stages.length : stages.indexOf(progress.stage);
    String label(TreeStage s) => switch (s) {
          TreeStage.prepare => l10n.treeStagePrepare,
          TreeStage.modelFinder => l10n.treeStageModel,
          TreeStage.search => l10n.treeStageSearch,
          TreeStage.support => progress.raxml ? l10n.treeStageBootstrap : l10n.treeStageSupport,
          TreeStage.finish => l10n.treeStageFinish,
        };
    return Row(crossAxisAlignment: CrossAxisAlignment.start, children: [
      for (final (i, s) in stages.indexed) ...[
        Expanded(
          child: Column(children: [
            Row(children: [
              Expanded(child: Container(height: 2, color: i == 0 ? Colors.transparent : (i <= current ? t.accent : t.border))),
              Container(
                width: 26,
                height: 26,
                alignment: Alignment.center,
                decoration: BoxDecoration(
                  shape: BoxShape.circle,
                  color: i < current ? t.accent : (i == current ? t.accentSoft : t.surface2),
                  border: Border.all(color: i == current && failed ? SplaceTokens.danger : (i <= current ? t.accent : t.border), width: 1.5),
                ),
                child: i < current
                    ? const Icon(FA.check, size: 11, color: Colors.white)
                    : i == current && failed
                        ? const Icon(FA.xmark, size: 11, color: SplaceTokens.danger)
                        : i == current
                        ? SizedBox(width: 12, height: 12, child: CircularProgressIndicator(strokeWidth: 2, color: t.accent))
                        : Text('${i + 1}', style: TextStyle(fontSize: 11, fontWeight: FontWeight.w700, color: t.textFaint)),
              ),
              Expanded(child: Container(height: 2, color: i == stages.length - 1 ? Colors.transparent : (i < current ? t.accent : t.border))),
            ]),
            const SizedBox(height: 6),
            Text(label(s),
                textAlign: TextAlign.center,
                style: TextStyle(
                  fontSize: 12,
                  fontWeight: i == current ? FontWeight.w700 : FontWeight.w500,
                  color: i <= current ? t.text : t.textFaint,
                )),
          ]),
        ),
      ],
    ]);
  }
}

/// One of several programs running together: its name and state, its
/// stages, what it is doing and its last lines.
class _ProgramRow extends StatelessWidget {
  const _ProgramRow({required this.name, required this.progress, required this.tail, this.item});
  final String name;
  final TreeProgress progress;
  final List<String> tail;
  final ModalItem? item;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final state = item?.state ?? ItemState.running;
    final failed = state == ItemState.failed;
    final done = state == ItemState.done;
    final (icon, color) = switch (state) {
      ItemState.done || ItemState.cached => (FA.circleCheck, SplaceTokens.success),
      ItemState.failed => (FA.circleXmark, SplaceTokens.danger),
      ItemState.skipped => (FA.circleMinus, t.textFaint),
      _ => (FA.spinner, t.accent),
    };
    return Container(
      padding: const EdgeInsets.all(14),
      decoration: BoxDecoration(
        color: t.surface,
        borderRadius: BorderRadius.circular(12),
        border: Border.all(color: failed ? SplaceTokens.danger : t.border),
      ),
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        Row(children: [
          if (state == ItemState.running || state == ItemState.queued)
            SizedBox(width: 14, height: 14, child: CircularProgressIndicator(strokeWidth: 2, color: t.accent))
          else
            Icon(icon, size: 14, color: color),
          const SizedBox(width: 8),
          Text(name, style: TextStyle(fontSize: 14, fontWeight: FontWeight.w700, color: t.text)),
          const SizedBox(width: 10),
          Expanded(
            child: Text(
              item?.note ?? '',
              maxLines: 2,
              overflow: TextOverflow.ellipsis,
              style: TextStyle(fontSize: 12.5, color: failed ? SplaceTokens.danger : t.textMuted, fontWeight: failed ? FontWeight.w600 : FontWeight.w400),
            ),
          ),
        ]),
        const SizedBox(height: 10),
        // This program's own bar (the one below the list is the whole run).
        Row(children: [
          Expanded(child: RainbowProgressBar(value: done ? 1 : progress.fraction, height: 6)),
          const SizedBox(width: 8),
          SizedBox(
            width: 38,
            child: Text('${((done ? 1 : progress.fraction) * 100).round()}%',
                textAlign: TextAlign.right, style: TextStyle(fontSize: 12, fontWeight: FontWeight.w700, color: t.accent)),
          ),
        ]),
        const SizedBox(height: 12),
        _TreeStages(progress: progress, finished: done, failed: failed),
        if (!failed && !done) ...[const SizedBox(height: 10), _treeDetail(context, progress)],
        if (tail.isNotEmpty) ...[const SizedBox(height: 10), _Tail(lines: tail.length > 4 ? tail.sublist(tail.length - 4) : tail)],
      ]),
    );
  }
}

/// Last lines printed by the program (monospace).
class _Tail extends StatelessWidget {
  const _Tail({required this.lines});
  final List<String> lines;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    return Container(
      width: double.infinity,
      padding: const EdgeInsets.symmetric(horizontal: 12, vertical: 8),
      decoration: BoxDecoration(color: t.surface2, borderRadius: BorderRadius.circular(8), border: Border.all(color: t.border)),
      child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
        for (final l in lines)
          Text(l, maxLines: 1, overflow: TextOverflow.ellipsis, style: TextStyle(fontFamily: monoFont, fontSize: 11.5, color: t.textMuted, height: 1.45)),
      ]),
    );
  }
}

String _modeName(AppLocalizations l10n, String mode) => switch (mode) {
      'aminoAcid' => l10n.seqAaTitle,
      'codon' => l10n.seqCodonTitle,
      _ => l10n.seqNtTitle,
    };

/// One progress bar per marker (alignment and trimming), with its stage or
/// result.
class _ItemList extends StatelessWidget {
  const _ItemList({required this.items, this.noteFor});
  final Map<String, ModalItem> items;

  /// Localized note of a finished item (saturation verdicts).
  final String? Function(String name)? noteFor;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final l10n = context.l10n;
    final counts = <ItemState, int>{};
    for (final i in items.values) {
      counts.update(i.state, (v) => v + 1, ifAbsent: () => 1);
    }

    // Each marker keeps its place for the whole run (cards never move).
    final ordered = items.entries.toList();

    // One card per marker, three per row.
    Widget cell(String name, ModalItem it) {
      final (icon, color, status) = switch (it.state) {
        ItemState.queued => (FA.clock, t.textFaint, l10n.itemQueued),
        ItemState.running => (FA.spinner, t.accent, it.note ?? ''),
        ItemState.done => (FA.circleCheck, t.accent, noteFor?.call(name) ?? it.note ?? ''),
        ItemState.cached => (FA.clockRotateLeft, SplaceTokens.success, l10n.itemCached(it.note ?? '')),
        ItemState.skipped => (FA.circleMinus, SplaceTokens.warn, it.note ?? l10n.itemSkipped),
        ItemState.failed => (FA.circleXmark, SplaceTokens.danger, it.note ?? ''),
      };
      // Keyed by marker: when the cards reorder (running ones first) a
      // card never animates from another marker's bar.
      return Container(
        key: ValueKey(name),
        padding: const EdgeInsets.fromLTRB(10, 8, 10, 8),
        decoration: BoxDecoration(
          color: t.surface,
          borderRadius: BorderRadius.circular(10),
          border: Border.all(color: it.state == ItemState.running ? t.accentSoftBorder : t.border),
        ),
        child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
          Row(children: [
            Icon(icon, size: 12, color: color),
            const SizedBox(width: 7),
            Expanded(
              child: Text(name,
                  maxLines: 1, overflow: TextOverflow.ellipsis, style: TextStyle(fontSize: 13, fontWeight: FontWeight.w700, color: t.text)),
            ),
            if (it.state == ItemState.running)
              Text('${(it.progress * 100).round()}%', style: TextStyle(fontSize: 12, fontWeight: FontWeight.w700, color: t.accent)),
          ]),
          const SizedBox(height: 4),
          SplaceTooltip(
            message: status,
            child: Text(status,
                maxLines: 1,
                overflow: TextOverflow.ellipsis,
                style: TextStyle(fontSize: 12, color: it.state == ItemState.failed ? SplaceTokens.danger : t.textMuted)),
          ),
          const SizedBox(height: 6),
          // The marker's own bar, in the colours of every progress bar of
          // the app: its stages while it runs, full when done.
          RainbowProgressBar(
            height: 5,
            value: switch (it.state) {
              ItemState.queued => 0.0,
              ItemState.running => it.progress.clamp(0.0, 1.0),
              _ => 1.0,
            },
          ),
        ]),
      );
    }

    return Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
      Wrap(spacing: 8, runSpacing: 6, children: [
        if ((counts[ItemState.running] ?? 0) > 0) TagChip(l10n.itemsRunning('${counts[ItemState.running]}'), color: t.accent),
        if ((counts[ItemState.queued] ?? 0) > 0) TagChip(l10n.itemsQueued('${counts[ItemState.queued]}'), color: t.textMuted),
        if ((counts[ItemState.done] ?? 0) > 0) TagChip(l10n.itemsDone('${counts[ItemState.done]}'), color: t.accent),
        if ((counts[ItemState.cached] ?? 0) > 0) TagChip(l10n.itemsCached('${counts[ItemState.cached]}'), color: SplaceTokens.success),
        if ((counts[ItemState.skipped] ?? 0) > 0) TagChip(l10n.itemsSkipped('${counts[ItemState.skipped]}'), color: SplaceTokens.warn),
        if ((counts[ItemState.failed] ?? 0) > 0) TagChip(l10n.itemsFailed('${counts[ItemState.failed]}'), color: SplaceTokens.danger),
      ]),
      const SizedBox(height: 10),
      Expanded(
        child: Container(
        width: double.infinity,
        padding: const EdgeInsets.all(10),
        decoration: BoxDecoration(color: t.surface2, borderRadius: BorderRadius.circular(10), border: Border.all(color: t.border)),
        // Built lazily (only the visible cards), so hundreds of markers
        // stay smooth; the cards stay in the markers' order.
        child: GridView.builder(
          padding: const EdgeInsets.symmetric(vertical: 4),
          gridDelegate: const SliverGridDelegateWithFixedCrossAxisCount(
            crossAxisCount: 3,
            crossAxisSpacing: 10,
            mainAxisSpacing: 10,
            mainAxisExtent: 74,
          ),
          itemCount: ordered.length,
          itemBuilder: (context, i) => cell(ordered[i].key, ordered[i].value),
        ),
        ),
      ),
    ]);
  }
}
