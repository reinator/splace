import 'package:flutter/material.dart';

import '../../core/saturation.dart';
import '../../services/tool_locator.dart' show ExternalTool;
import '../../state/project_controller.dart';
import '../common.dart';
import '../fa.dart';
import '../theme.dart';
import '../widgets/examples.dart';
import '../widgets/citations.dart';
import '../widgets/kit.dart';
import '../widgets/saturation_plot.dart';
import '../widgets/splace_tooltip.dart';

/// Step 6 — saturation test by codon position (new in the Flutter app).
class SaturationSection extends StatefulWidget {
  const SaturationSection({super.key, required this.active});
  final bool active;

  @override
  State<SaturationSection> createState() => _SaturationSectionState();
}

/// "COI, CYTB" or "COI, CYTB +3".
String _names(Iterable<String> names) {
  final l = names.toList();
  return l.length <= 4 ? l.join(', ') : '${l.take(4).join(', ')} +${l.length - 4}';
}

class _SaturationSectionState extends State<SaturationSection> {
  String? _selected;
  CodonSubset _subset = CodonSubset.pos3;

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final reports = c.loci.values.where((l) => l.included).map((l) => l.saturation).whereType<SaturationReport>().toList();
    final pooled = c.pooledSaturation;
    final selected = reports.where((r) => r.name == _selected).firstOrNull ?? reports.firstOrNull;
    final coding = c.loci.values.any((l) => l.included && c.codonPositionsUsable(l) && l.alignment != null);
    // Loci without codons are tested on all their sites.
    final canRun = c.canTestSaturation;

    return StepSection(
      badge: '5',
      title: l10n.step6,
      description: l10n.saturationIntro,
      active: widget.active,
      child: Column(
        crossAxisAlignment: CrossAxisAlignment.stretch,
        children: [
          SubAnchor(id: '5a', child: SubHeading(id: '5a', title: l10n.subRunSaturation)),
          const _SaturationParams(),
          const SizedBox(height: 12),
          if (c.saturationStale) ...[
            TonePanel(tone: PanelTone.warn, icon: FA.rotateLeft, child: Text(l10n.satParamsStale, style: const TextStyle(fontSize: 13))),
            const SizedBox(height: 12),
          ],
          Align(
            alignment: Alignment.centerLeft,
            child: PrimaryButton(label: l10n.runSaturation, icon: FA.chartLine, onPressed: c.busy || !canRun ? null : c.runSaturation),
          ),
          if (reports.isEmpty && c.lociMode) ...[
            const SizedBox(height: 12),
            TonePanel(
              tone: PanelTone.info,
              icon: FA.circleInfo,
              child: Text(coding ? l10n.satOptionalLoci : l10n.satOptionalLociNoCoding, style: const TextStyle(fontSize: 13, height: 1.4)),
            ),
          ] else if (reports.isEmpty && c.partialLoci.isNotEmpty && !c.saturationRequired) ...[
            const SizedBox(height: 12),
            TonePanel(
              tone: PanelTone.info,
              icon: FA.circleInfo,
              child: Text(l10n.satOptionalPartial(_names(c.partialLoci.map((l) => l.name))), style: const TextStyle(fontSize: 13, height: 1.4)),
            ),
          ] else if (reports.isEmpty && coding) ...[
            const SizedBox(height: 12),
            TonePanel(tone: PanelTone.info, icon: FA.lock, child: Text(l10n.satUnlocksStep7, style: const TextStyle(fontSize: 13, height: 1.4))),
            if (c.partialLoci.isNotEmpty) ...[
              const SizedBox(height: 8),
              TonePanel(
                tone: PanelTone.info,
                icon: FA.circleInfo,
                child: Text(l10n.satPartialMixed(_names(c.partialLoci.map((l) => l.name))), style: const TextStyle(fontSize: 13, height: 1.4)),
              ),
            ],
          ],
          if (reports.isNotEmpty) ...[
            const SizedBox(height: 22),
            SubAnchor(
              id: '5b',
              child: SubHeading(id: '5b', title: l10n.saturationHeatmapTitle, body: l10n.saturationHeatmapHelp),
            ),
            const SizedBox(height: 8),
            _Heatmap(reports: reports, selected: selected?.name, onSelect: (n) => setState(() => _selected = n)),
          ],
          // Table, then the plots of the selected row, then the suggestion
          // (which leads into step 7: concatenation and partitions).
          if (selected != null) ...[const SizedBox(height: 14), _detail(context, selected)],
          if (pooled != null) ...[const SizedBox(height: 18), SubAnchor(id: '5c', child: _Recommendation(controller: c, pooled: pooled))],
        ],
      ),
    );
  }

  /// Which alignment of the locus was tested, and how trimming kept its
  /// codon positions.
  String? _source(BuildContext context, SaturationReport r) {
    final c = ProjectScope.read(context);
    final l10n = context.l10n;
    final l = c.loci[r.name];
    if (l == null || l.aligned == null) return null;
    if (!c.useTrimmed || l.trimmed == null) return l10n.satSourceMafft;
    if (l.proteinAligned != null && l.trimmedColumns == null) return l10n.satSourceClipkit;
    final cols = l.trimmedColumns;
    if (cols == null) return null;
    final total = (l.aligned!.first.length - l.frame) ~/ 3;
    final text = l10n.satSourceTrimmed('${cols.length ~/ 3}', '$total');
    return l.frame == 0 ? text : '$text ${l10n.satSourceFrame('${l.frame + 1}')}';
  }

  Widget _detail(BuildContext context, SaturationReport r) {
    final l10n = context.l10n;
    final t = context.tokens;
    // Non-coding loci have only "all sites".
    final subset = r.subsets.containsKey(_subset) ? _subset : (r.subsets.keys.firstOrNull ?? _subset);
    final s = r.subsets[subset];
    final codons = r.subsets.keys.any((k) => k != CodonSubset.all);
    final name = r.name;
    return Container(
      padding: const EdgeInsets.all(16),
      decoration: BoxDecoration(
        color: t.surface,
        borderRadius: BorderRadius.circular(16),
        border: Border.all(color: t.border),
        boxShadow: shadowSm,
      ),
      child: Column(
        crossAxisAlignment: CrossAxisAlignment.stretch,
        children: [
          Row(
            children: [
              Text(
                name,
                style: TextStyle(fontSize: 16, fontWeight: FontWeight.w600, color: t.text),
              ),
              const SizedBox(width: 10),
              Flexible(
                child: SplaceTooltip(
                  message: l10n.frameInfo,
                  child: Text.rich(
                    TextSpan(
                      children: [
                        TextSpan(text: '${r.taxa} taxa · ${r.length} bp${codons ? ' · ${l10n.frameLabel('${r.frame + 1}')} ' : ' · ${l10n.satAllSitesNonCoding} '}'),
                        WidgetSpan(
                          alignment: PlaceholderAlignment.middle,
                          child: Icon(FA.circleInfo, size: 11, color: t.textFaint),
                        ),
                      ],
                    ),
                    style: TextStyle(fontSize: 13, color: t.textMuted),
                  ),
                ),
              ),
            ],
          ),
          if (_source(context, r) case final source?) ...[
            const SizedBox(height: 8),
            Row(crossAxisAlignment: CrossAxisAlignment.start, children: [
              Padding(padding: const EdgeInsets.only(top: 3), child: Icon(FA.scissors, size: 11, color: t.textFaint)),
              const SizedBox(width: 8),
              Expanded(child: Text(source, style: TextStyle(fontSize: 12.5, color: t.textMuted, height: 1.45))),
            ]),
          ],
          const SizedBox(height: 12),
          Wrap(
            spacing: 6,
            runSpacing: 6,
            children: [
              for (final sub in r.subsets.keys)
                GestureDetector(
                  onTap: () => setState(() => _subset = sub),
                  child: MouseRegion(
                    cursor: SystemMouseCursors.click,
                    child: Container(
                      padding: const EdgeInsets.symmetric(horizontal: 13, vertical: 7),
                      decoration: BoxDecoration(
                        color: sub == subset ? t.accent : t.surface2,
                        borderRadius: BorderRadius.circular(8),
                        border: Border.all(color: sub == subset ? t.accent : t.border),
                      ),
                      child: Row(
                        mainAxisSize: MainAxisSize.min,
                        children: [
                          Container(
                            width: 8,
                            height: 8,
                            decoration: BoxDecoration(color: verdictColor(r.subsets[sub]!.verdict), shape: BoxShape.circle),
                          ),
                          const SizedBox(width: 6),
                          Text(
                            sub.label,
                            style: TextStyle(fontSize: 13, fontWeight: FontWeight.w600, color: sub == subset ? Colors.white : t.textMuted),
                          ),
                        ],
                      ),
                    ),
                  ),
                ),
            ],
          ),
          if (s != null) ...[
            const SizedBox(height: 10),
            Row(crossAxisAlignment: CrossAxisAlignment.start, children: [
              Container(
                margin: const EdgeInsets.only(top: 4),
                width: 9,
                height: 9,
                decoration: BoxDecoration(color: verdictColor(s.verdict), shape: BoxShape.circle),
              ),
              const SizedBox(width: 8),
              Expanded(
                child: Text.rich(
                  TextSpan(children: [
                    TextSpan(text: '${s.subset.label} · ${verdictLabel(context, s.verdict)}: ', style: const TextStyle(fontWeight: FontWeight.w700)),
                    TextSpan(text: verdictReason(context, s, ProjectScope.read(context).saturationThresholds)),
                  ]),
                  style: TextStyle(fontSize: 13, color: t.text, height: 1.4),
                ),
              ),
            ]),
          ],
          const SizedBox(height: 14),
          const SaturationReadingExample(),
          const SizedBox(height: 14),
          if (s != null)
            LayoutBuilder(builder: (context, constraints) {
              final plots = [
                ScatterPlot(
                  title: l10n.plotPhykit(fmt(s.slope)),
                  subtitle: l10n.plotPhykitBody,
                  xLabel: l10n.patristicDistance,
                  yLabel: 'p',
                  identityLine: true,
                  series: observedVsCorrected(s),
                ),
                if (r.positionPairs.isNotEmpty)
                  ScatterPlot(
                    title: l10n.plotPositionsPatristic,
                    subtitle: l10n.plotPositionsBody,
                    xLabel: '${l10n.patristicDistance} (${l10n.allPositions})',
                    yLabel: 'p',
                    identityLine: true,
                    series: positionsVsDistance(r, '1st', '2nd', '3rd'),
                  ),
              ];
              // All plots in one row (same size), or one per row on narrow
              // pages: never two plus a lonely third.
              if (constraints.maxWidth >= 720 && plots.length > 1) return EqualHeightRow(gap: 20, children: plots);
              return Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
                for (final (i, p) in plots.indexed) Padding(padding: EdgeInsets.only(top: i == 0 ? 0 : 18), child: p),
              ]);
            }),
          const SizedBox(height: 16),
          RecordsTable(
            columns: [
              '',
              l10n.metricSites,
              l10n.metricVariable,
              l10n.metricInformative,
              l10n.metricPairs,
              l10n.metricMeanP,
              l10n.metricMeanPatristic,
              l10n.metricPhykit,
              '|1 − s|',
              l10n.metricVerdict,
            ],
            flex: const [2, 2, 2, 2, 2, 2, 2, 2, 2, 3],
            rowCount: r.subsets.length,
            cellBuilder: (context, row, col) {
              final x = r.subsets.values.elementAt(row);
              return switch (col) {
                0 => Text(x.subset.label, style: const TextStyle(fontWeight: FontWeight.w700)),
                1 => Text('${x.sites}'),
                2 => Text('${x.variableSites}'),
                3 => Text('${x.informativeSites}'),
                4 => Text('${x.pairs}'),
                5 => Text(fmt(x.meanP, 3)),
                6 => Text(fmt(x.meanDistance, 3)),
                7 => Text(fmt(x.slope, 3), style: const TextStyle(fontWeight: FontWeight.w700)),
                8 => Text(fmt(x.slope == null ? null : (1 - x.slope!).abs(), 3)),
                _ => Align(
                  alignment: Alignment.centerLeft,
                  child: TagChip(verdictLabel(context, x.verdict), color: verdictColor(x.verdict)),
                ),
              };
            },
          ),
          if (r.frameStops != null)
            Padding(
              padding: const EdgeInsets.only(top: 8),
              child: Text('Internal stop codons per frame (+1/+2/+3): ${r.frameStops!.join(' / ')}', style: TextStyle(fontSize: 12, color: t.textFaint)),
            ),
        ],
      ),
    );
  }
}

/// What the test means for the matrix: the verdict of each codon position
/// over the coding genes (each gene tested on its own), the suggestion in
/// plain words, why, and the button that marks it in step 7b.
class _Recommendation extends StatelessWidget {
  const _Recommendation({required this.controller, required this.pooled});
  final ProjectController controller;
  final PooledSaturation pooled;

  @override
  Widget build(BuildContext context) {
    final l10n = context.l10n;
    final t = context.tokens;
    final c = controller;
    final rec = pooled.recommendation;
    final exclude = rec.exclude, caution = rec.caution;
    final single = pooled.reports.length == 1 ? pooled.reports.single : null;
    final (tone, icon, text) = switch (rec.kind) {
      RecommendationKind.keepAll => (PanelTone.success, FA.circleCheck, l10n.recKeepAll),
      RecommendationKind.excludePositions => (
          PanelTone.danger,
          FA.scissors,
          exclude.length == 1 ? l10n.recExcludeOne(codonPositions(context, exclude)) : l10n.recExcludeMany(codonPositions(context, exclude)),
        ),
      RecommendationKind.useAminoAcids => (PanelTone.warn, FA.triangleExclamation, l10n.recAminoAcids),
      RecommendationKind.undetermined => (PanelTone.info, FA.circleQuestion, l10n.recUndetermined),
    };
    final differing = [
      for (final r in pooled.reports)
        if (!_same(r.recommendation.exclude, exclude)) r.name,
    ];
    final applicable = rec.kind == RecommendationKind.excludePositions || rec.kind == RecommendationKind.keepAll;
    final appliedGlobal = c.concat.perLocusExclude.isEmpty && _same(c.concat.excludePositions, exclude);
    final appliedPerLocus = c.concat.perLocusExclude.isNotEmpty;
    // The labels of step 7b's "Codon positions" chips.
    String quoted(Iterable<int> ps) => (ps.toList()..sort()).map((p) => '“${l10n.excludePosition('$p')}”').join(' ${l10n.listAnd} ');

    Widget position(int n, CodonSubset subset) {
      final v = pooled.verdict(n);
      final color = verdictColor(v);
      final x = single?.subsets[subset];
      // One gene: why it got its verdict. Several: the genes of each verdict.
      final Widget body = single != null
          ? Text(x == null ? '—' : verdictReason(context, x, ProjectScope.read(context).saturationThresholds), style: TextStyle(fontSize: 12.5, color: t.textMuted, height: 1.4))
          : Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
              Text(l10n.recGenesCount('${pooled.count(n, SaturationVerdict.saturated)}', '${pooled.assessed(n)}'),
                  style: TextStyle(fontSize: 12.5, color: t.textMuted, height: 1.4)),
              for (final verdict in [SaturationVerdict.saturated, SaturationVerdict.moderate, SaturationVerdict.ok, SaturationVerdict.insufficient])
                if (pooled.genes[n]![verdict] case final genes? when genes.isNotEmpty)
                  Padding(
                    padding: const EdgeInsets.only(top: 3),
                    child: Text.rich(
                      TextSpan(children: [
                        TextSpan(text: '${verdictLabel(context, verdict)}: ', style: TextStyle(fontWeight: FontWeight.w700, color: verdictColor(verdict))),
                        TextSpan(text: _names(genes)),
                      ]),
                      style: TextStyle(fontSize: 12, color: t.textMuted, height: 1.35),
                    ),
                  ),
            ]);
      // A coloured strip on the left of a rounded card.
      return Container(
        clipBehavior: Clip.antiAlias,
        decoration: BoxDecoration(color: t.surface, borderRadius: BorderRadius.circular(10), border: Border.all(color: t.border)),
        child: IntrinsicHeight(
          child: Row(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
            Container(width: 4, color: color),
            Expanded(
              child: Padding(
                padding: const EdgeInsets.fromLTRB(12, 10, 12, 10),
                child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
                  Row(children: [
                    Expanded(
                      child: Text(l10n.codonPositionN(codonPositions(context, [n])),
                          style: TextStyle(fontSize: 13.5, fontWeight: FontWeight.w700, color: t.text)),
                    ),
                    TagChip(verdictLabel(context, v), color: color),
                  ]),
                  const SizedBox(height: 6),
                  body,
                ]),
              ),
            ),
          ]),
        ),
      );
    }

    return Container(
      padding: const EdgeInsets.all(16),
      decoration: BoxDecoration(color: t.surface2, borderRadius: BorderRadius.circular(16), border: Border.all(color: t.border)),
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        Row(children: [
          const StepBadge('5c', size: 26),
          const SizedBox(width: 10),
          Expanded(child: Text(l10n.recTitle, style: TextStyle(fontSize: 16, fontWeight: FontWeight.w700, color: t.text))),
          InfoDot(l10n.saturationThresholds),
        ]),
        const SizedBox(height: 4),
        Text(
          single != null ? l10n.recBasedOnGene(single.name, '${single.taxa}', '${single.length}') : l10n.recBasedOnGenes('${pooled.reports.length}'),
          style: TextStyle(fontSize: 12.5, color: t.textMuted),
        ),
        const SizedBox(height: 12),
        // Verdict of each position with its reason.
        LayoutBuilder(builder: (context, box) {
          final cells = [
            position(1, CodonSubset.pos1),
            position(2, CodonSubset.pos2),
            position(3, CodonSubset.pos3),
          ];
          if (box.maxWidth < 720) {
            return Column(children: [
              for (final (i, w) in cells.indexed) Padding(padding: EdgeInsets.only(top: i == 0 ? 0 : 8), child: w),
            ]);
          }
          return IntrinsicHeight(
            child: Row(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
              for (final (i, w) in cells.indexed) ...[if (i > 0) const SizedBox(width: 10), Expanded(child: w)],
            ]),
          );
        }),
        const SizedBox(height: 14),
        TonePanel(
          tone: tone,
          icon: icon,
          title: l10n.recWhatToDo,
          child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
            Text(text, style: const TextStyle(fontWeight: FontWeight.w600, fontSize: 14, height: 1.4)),
            if (rec.kind == RecommendationKind.excludePositions || rec.kind == RecommendationKind.useAminoAcids) ...[
              const SizedBox(height: 6),
              Text(l10n.recWhy, style: const TextStyle(fontSize: 13, height: 1.4)),
            ],
            if (caution.isNotEmpty) ...[
              const SizedBox(height: 6),
              Text(
                caution.length == 1 ? l10n.recCautionOne(codonPositions(context, caution)) : l10n.recCautionMany(codonPositions(context, caution)),
                style: const TextStyle(fontSize: 13, height: 1.4),
              ),
            ],
          ]),
        ),
        const SizedBox(height: 12),
        Wrap(spacing: 8, runSpacing: 8, crossAxisAlignment: WrapCrossAlignment.center, children: [
          if (applicable)
            ChipButton(
              label: exclude.isEmpty ? l10n.applyKeepAll : l10n.applyExclude(quoted(exclude)),
              icon: appliedGlobal ? FA.circleCheck : FA.check,
              style: ChipButtonStyle.accent,
              onPressed: appliedGlobal ? null : () => c.applySaturationRecommendation(perLocus: false),
            ),
          if (differing.isNotEmpty)
            SplaceTooltip(
              message: l10n.applyPerLocusInfo(differing.join(', ')),
              child: ChipButton(
                label: l10n.applyPerLocus,
                icon: appliedPerLocus ? FA.circleCheck : FA.listCheck,
                style: ChipButtonStyle.soft,
                onPressed: () => c.applySaturationRecommendation(perLocus: true),
              ),
            ),
          Text(_status(context), style: TextStyle(fontSize: 12.5, color: t.textMuted)),
        ]),
        if (differing.isNotEmpty) ...[
          const SizedBox(height: 8),
          Text(l10n.recLociDiffer('${differing.length}'), style: TextStyle(fontSize: 12.5, color: t.textMuted)),
        ],
      ]),
    );
  }

  static bool _same(Set<int> a, Set<int> b) => a.length == b.length && a.containsAll(b);

  /// What step 7 will do now.
  String _status(BuildContext context) {
    final l10n = context.l10n;
    final o = controller.concat;
    if (o.perLocusExclude.isNotEmpty) return l10n.recStatusPerLocus('${o.perLocusExclude.length}');
    if (o.excludePositions.isEmpty) return l10n.recStatusKeep;
    return l10n.recStatusExclude(codonPositions(context, o.excludePositions));
  }
}

/// Loci × codon position, coloured by verdict (slope and Iss in each
/// cell). Like the gene heatmaps: the bar on top describes the hovered
/// cell, its row and column are highlighted, and a click shows the locus.
class _Heatmap extends StatefulWidget {
  const _Heatmap({required this.reports, required this.selected, required this.onSelect});
  final List<SaturationReport> reports;
  final String? selected;
  final ValueChanged<String> onSelect;

  @override
  State<_Heatmap> createState() => _HeatmapState();
}

/// Heatmap columns: every codon subset, or only "all" when no locus has
/// codon positions (non-coding loci).
List<CodonSubset> _columns(List<SaturationReport> reports) =>
    reports.any((r) => r.subsets.keys.any((k) => k != CodonSubset.all)) ? CodonSubset.values : const [CodonSubset.all];

class _HeatmapState extends State<_Heatmap> {
  /// Hovered (row, column); either may be null (names or header).
  final _hover = ValueNotifier<(int?, int?)>((null, null));

  @override
  void dispose() {
    _hover.dispose();
    super.dispose();
  }

  String _name(BuildContext context, SaturationReport r) => r.name;

  String _describe(BuildContext context, (int?, int?) h) {
    final l10n = context.l10n;
    final subsets = _columns(widget.reports);
    final (row, col) = h;
    if (row == null && col == null) return l10n.selectLocusHint;
    if (row == null) {
      final s = subsets[col!];
      final counts = <SaturationVerdict, int>{};
      for (final r in widget.reports) {
        final v = r.subsets[s]?.verdict;
        if (v != null) counts.update(v, (n) => n + 1, ifAbsent: () => 1);
      }
      return '${s.label}: ${counts.entries.map((e) => '${verdictLabel(context, e.key)} ${e.value}').join(' · ')}';
    }
    final r = widget.reports[row];
    if (col == null) {
      return '${_name(context, r)}: ${r.taxa} taxa · ${r.length} bp · '
          '${r.subsets.entries.map((e) => '${e.key.label} ${verdictLabel(context, e.value.verdict)}').join(' · ')}';
    }
    final s = r.subsets[subsets[col]];
    if (s == null) return '${_name(context, r)} × ${subsets[col].label}: —';
    return '${_name(context, r)} × ${s.subset.label}: ${verdictLabel(context, s.verdict)}  ·  PhyKIT ${fmt(s.slope, 3)} · '
        '|1 − s| ${fmt(s.slope == null ? null : (1 - s.slope!).abs())} · ${s.pairs} ${l10n.metricPairs.toLowerCase()} · p ${fmt(s.meanP, 3)}';
  }

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final l10n = context.l10n;
    final subsets = _columns(widget.reports);
    const nameW = 190.0, rowH = 40.0, headH = 34.0;
    final reports = widget.reports;
    return LayoutBuilder(builder: (context, box) {
      // The cells share the width of the window (at least 96 px each).
      final cellW = ((box.maxWidth - nameW - 2) / subsets.length).clamp(96.0, double.infinity);
      final gridW = nameW + cellW * subsets.length;
      return Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        Container(
          decoration: BoxDecoration(borderRadius: BorderRadius.circular(12), border: Border.all(color: t.border), color: t.surface),
          clipBehavior: Clip.antiAlias,
          child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
            ValueListenableBuilder<(int?, int?)>(
              valueListenable: _hover,
              builder: (context, h, _) => Container(
                color: t.surface2,
                padding: const EdgeInsets.symmetric(horizontal: 12, vertical: 7),
                child: Text(
                  _describe(context, h),
                  maxLines: 1,
                  overflow: TextOverflow.ellipsis,
                  style: TextStyle(fontSize: 12.5, color: h == (null, null) ? t.textMuted : t.text),
                ),
              ),
            ),
            SingleChildScrollView(
              scrollDirection: Axis.horizontal,
              child: SizedBox(
                width: gridW,
                child: MouseRegion(
                  onExit: (_) => _hover.value = (null, null),
                  child: ValueListenableBuilder<(int?, int?)>(
                    valueListenable: _hover,
                    builder: (context, h, _) {
                      final (hr, hc) = h;
                      Widget header(int c) => MouseRegion(
                            onEnter: (_) => _hover.value = (null, c),
                            child: Container(
                              width: cellW,
                              height: headH,
                              alignment: Alignment.center,
                              color: hc == c ? t.accentSoft : null,
                              child: Text(
                                subsets[c].label,
                                style: TextStyle(fontSize: 12.5, fontWeight: hc == c ? FontWeight.w800 : FontWeight.w600, color: hc == c ? t.accentText : t.textMuted),
                              ),
                            ),
                          );
                      return Column(children: [
                        Row(children: [
                          SizedBox(
                            width: nameW,
                            height: headH,
                            child: Padding(
                              padding: const EdgeInsets.only(left: 12),
                              child: Align(alignment: Alignment.centerLeft, child: Kicker(l10n.colLocus)),
                            ),
                          ),
                          for (var c = 0; c < subsets.length; c++) header(c),
                        ]),
                        Divider(height: 1, color: t.border),
                        ConstrainedBox(
                          constraints: const BoxConstraints(maxHeight: 460),
                          child: ListView.builder(
                            shrinkWrap: true,
                            padding: const EdgeInsets.symmetric(vertical: 3),
                            itemCount: reports.length,
                            itemExtent: rowH,
                            itemBuilder: (context, i) {
                              final r = reports[i];
                              final selected = r.name == widget.selected;
                              final rowOn = hr == i;
                              return GestureDetector(
                                onTap: () => widget.onSelect(r.name),
                                child: MouseRegion(
                                  cursor: SystemMouseCursors.click,
                                  child: Container(
                                    decoration: BoxDecoration(
                                      color: selected ? t.accentSoft : (rowOn ? t.accentSoft.withValues(alpha: 0.5) : null),
                                      border: Border(left: BorderSide(color: selected ? t.accent : Colors.transparent, width: 3)),
                                    ),
                                    child: Row(children: [
                                      MouseRegion(
                                        onEnter: (_) => _hover.value = (i, null),
                                        child: SizedBox(
                                          width: nameW - 3,
                                          child: Padding(
                                            padding: const EdgeInsets.only(left: 9, right: 6),
                                            child: Text(
                                              _name(context, r),
                                              overflow: TextOverflow.ellipsis,
                                              style: TextStyle(
                                                fontSize: 13,
                                                fontWeight: rowOn ? FontWeight.w700 : FontWeight.w500,
                                                color: selected || rowOn ? t.accentText : t.text,
                                              ),
                                            ),
                                          ),
                                        ),
                                      ),
                                      for (var c = 0; c < subsets.length; c++)
                                        MouseRegion(
                                          onEnter: (_) => _hover.value = (i, c),
                                          child: Container(
                                            width: cellW,
                                            height: rowH,
                                            color: hc == c ? t.accentSoft.withValues(alpha: 0.6) : null,
                                            padding: const EdgeInsets.all(3),
                                            child: _cell(context, t, r.subsets[subsets[c]], outlined: rowOn && hc == c, outline: t.text),
                                          ),
                                        ),
                                    ]),
                                  ),
                                ),
                              );
                            },
                          ),
                        ),
                      ]);
                    },
                  ),
                ),
              ),
            ),
          ]),
        ),
        const SizedBox(height: 8),
        Wrap(spacing: 14, runSpacing: 6, crossAxisAlignment: WrapCrossAlignment.center, children: [
          Row(mainAxisSize: MainAxisSize.min, children: [
            Text(l10n.verdictOk, style: TextStyle(fontSize: 12, fontWeight: FontWeight.w600, color: t.textMuted)),
            const SizedBox(width: 6),
            Container(
              width: 220,
              height: 12,
              decoration: BoxDecoration(
                borderRadius: BorderRadius.circular(3),
                gradient: LinearGradient(colors: [for (final (_, c) in saturationScaleStops) c], stops: [for (final (x, _) in saturationScaleStops) x]),
              ),
            ),
            const SizedBox(width: 6),
            Text(l10n.verdictSaturated, style: TextStyle(fontSize: 12, fontWeight: FontWeight.w600, color: t.textMuted)),
          ]),
          Row(mainAxisSize: MainAxisSize.min, children: [
            Container(width: 12, height: 12, decoration: BoxDecoration(color: t.surface3, borderRadius: BorderRadius.circular(3), border: Border.all(color: t.border))),
            const SizedBox(width: 5),
            Text(l10n.verdictInsufficient, style: TextStyle(fontSize: 12, color: t.textMuted)),
          ]),
          Text(l10n.heatmapIntensity, style: TextStyle(fontSize: 12, color: t.textFaint)),
        ]),
      ]);
    });
  }

  Widget _cell(BuildContext context, SplaceTokens t, SubsetSaturation? s, {required bool outlined, required Color outline}) {
    if (s == null) return const SizedBox.shrink();
    final score = saturationScore(s);
    final color = score == null ? t.surface3 : saturationScale(score);
    return Container(
      decoration: BoxDecoration(
        color: color,
        borderRadius: BorderRadius.circular(6),
        border: outlined ? Border.all(color: outline, width: 2) : null,
      ),
      alignment: Alignment.center,
      child: Text(
        score == null ? verdictLabel(context, s.verdict) : fmt(s.slope),
        maxLines: 1,
        overflow: TextOverflow.ellipsis,
        style: TextStyle(color: score == null ? t.textMuted : onColor(color), fontSize: 11.5, fontWeight: FontWeight.w600),
      ),
    );
  }
}

/// Step 6a: how the test is made and its parameters, all visible; the
/// thresholds, the Iss replicates and the plot points can be changed
/// (they apply on the next run).
class _SaturationParams extends StatelessWidget {
  const _SaturationParams();

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    final th = c.saturationThresholds;
    SaturationThresholds with_({
      double? slopeModerate,
      double? slopeSaturated,
      int? minPairs,
      int? minSites,
      double? lowDivergence,
    }) =>
        SaturationThresholds(
          slopeModerate: slopeModerate ?? th.slopeModerate,
          slopeSaturated: slopeSaturated ?? th.slopeSaturated,
          minPairs: minPairs ?? th.minPairs,
          minSites: minSites ?? th.minSites,
          lowDivergence: lowDivergence ?? th.lowDivergence,
        );
    Widget number(String label, double value, double min, double max, ValueChanged<double> onChanged,
            {double step = 0.05, int decimals = 2, String? suffix, String? info}) =>
        // One line per title at any width (the whole text in the tooltip),
        // all fields the same size.
        Column(crossAxisAlignment: CrossAxisAlignment.stretch, mainAxisSize: MainAxisSize.min, children: [
          SizedBox(
            height: 20,
            child: Row(children: [
              Expanded(
                child: SplaceTooltip(
                  message: label,
                  child: Text(label,
                      maxLines: 1, softWrap: false, overflow: TextOverflow.ellipsis, style: TextStyle(fontSize: 13, fontWeight: FontWeight.w600, color: t.text)),
                ),
              ),
              if (info != null) ...[const SizedBox(width: 6), InfoDot(info)],
            ]),
          ),
          const SizedBox(height: 6),
          NumberField(value: value, min: min, max: max, step: step, decimals: decimals, suffix: suffix, width: double.infinity, onChanged: onChanged),
        ]);
    // Title (one line) and any field below it.
    Widget label(String text, String? info, Widget child) => Column(crossAxisAlignment: CrossAxisAlignment.stretch, mainAxisSize: MainAxisSize.min, children: [
          SizedBox(
            height: 20,
            child: Row(children: [
              Expanded(
                child: SplaceTooltip(
                  message: text,
                  child: Text(text,
                      maxLines: 1, softWrap: false, overflow: TextOverflow.ellipsis, style: TextStyle(fontSize: 13, fontWeight: FontWeight.w600, color: t.text)),
                ),
              ),
              if (info != null) ...[const SizedBox(width: 6), InfoDot(info)],
            ]),
          ),
          const SizedBox(height: 6),
          child,
        ]);
    Widget fact(String label, String value) => Padding(
          padding: const EdgeInsets.only(bottom: 4),
          child: Text.rich(TextSpan(style: TextStyle(fontSize: 13, color: t.text, height: 1.45), children: [
            TextSpan(text: '$label: ', style: TextStyle(color: t.textMuted)),
            TextSpan(style: const TextStyle(fontWeight: FontWeight.w600), children: citationSpans(value, fontSize: 12)),
          ])),
        );
    return Container(
      margin: const EdgeInsets.only(bottom: 4),
      padding: const EdgeInsets.all(16),
      decoration: BoxDecoration(color: t.surface, borderRadius: BorderRadius.circular(14), border: Border.all(color: t.border)),
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        Row(children: [
          Icon(FA.sliders, size: 14, color: t.accent),
          const SizedBox(width: 8),
          Text(l10n.satParamsTitle, style: TextStyle(fontSize: 15, fontWeight: FontWeight.w700, color: t.text)),
        ]),
        const SizedBox(height: 4),
        Text(l10n.satParamsBody, style: TextStyle(fontSize: 13, color: t.textMuted, height: 1.4)),
        const SizedBox(height: 12),
        Container(
          padding: const EdgeInsets.all(12),
          decoration: BoxDecoration(color: t.surface2, borderRadius: BorderRadius.circular(10), border: Border.all(color: t.border)),
          child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
            fact(l10n.satParamMethod, l10n.satParamMethodValue),
            fact(l10n.satParamDistance, l10n.satParamDistanceValue),
            fact(l10n.satParamSubsets, l10n.satParamSubsetsValue),
            fact(l10n.satParamVersions, [
              'PhyKIT ${c.tools[ExternalTool.phykit]?.version ?? '—'}',
              'IQ-TREE ${c.tools[ExternalTool.iqtree]?.version ?? '—'}',
            ].join(' · ')),
            fact(l10n.satParamFrame, l10n.satParamFrameValue),
          ]),
        ),
        const SizedBox(height: 14),
        Kicker(l10n.satParamsThresholds),
        const SizedBox(height: 8),
        AutoGrid(minItemWidth: 210, children: [
          number(l10n.satParamSlopeModerate, th.slopeModerate, 0, 1, (v) => c.setSaturationParams(thresholds: with_(slopeModerate: v)),
              info: l10n.satParamSlopeInfo),
          number(l10n.satParamSlopeSaturated, th.slopeSaturated, 0, 1, (v) => c.setSaturationParams(thresholds: with_(slopeSaturated: v))),
          number(l10n.satParamLowDivergence, th.lowDivergence, 0, 0.5, (v) => c.setSaturationParams(thresholds: with_(lowDivergence: v)),
              step: 0.01, decimals: 3, info: l10n.satParamLowDivergenceInfo),
          number(l10n.satParamMinSites, th.minSites.toDouble(), 1, 10000, (v) => c.setSaturationParams(thresholds: with_(minSites: v.round())),
              step: 5, decimals: 0),
          number(l10n.satParamMinPairs, th.minPairs.toDouble(), 1, 10000, (v) => c.setSaturationParams(thresholds: with_(minPairs: v.round())),
              step: 1, decimals: 0),
          number(l10n.satParamPlotPoints, c.satPlotPoints.toDouble(), 100, 50000, (v) => c.setSaturationParams(plotPoints: v.round()),
              step: 500, decimals: 0, info: l10n.satParamPlotPointsInfo),
          label(l10n.satParamModel, l10n.satParamModelInfo, SizedBox(
            height: 36,
            child: TextFormField(
              key: ValueKey('satModel-${c.satModel}'),
              initialValue: c.satModel,
              style: const TextStyle(fontFamily: monoFont, fontSize: 13),
              decoration: const InputDecoration(isDense: true, contentPadding: EdgeInsets.symmetric(horizontal: 10, vertical: 10)),
              onFieldSubmitted: (v) => c.setSaturationParams(model: v),
              onChanged: (v) => c.setSaturationParams(model: v),
            ),
          )),
          // A check box of the same height and frame as the other fields.
          label(
            l10n.satParamGaps,
            l10n.satParamGapsInfo,
            InkWell(
              borderRadius: BorderRadius.circular(8),
              onTap: () => c.setSaturationParams(excludeGaps: !c.satExcludeGaps),
              child: Container(
                height: 36,
                padding: const EdgeInsets.symmetric(horizontal: 10),
                decoration: BoxDecoration(
                  color: c.satExcludeGaps ? t.accentSoft : t.surface,
                  borderRadius: BorderRadius.circular(8),
                  border: Border.all(color: c.satExcludeGaps ? t.accentSoftBorder : t.border),
                ),
                child: Row(children: [
                  Container(
                    width: 16,
                    height: 16,
                    alignment: Alignment.center,
                    decoration: BoxDecoration(
                      color: c.satExcludeGaps ? t.accent : Colors.transparent,
                      borderRadius: BorderRadius.circular(4),
                      border: Border.all(color: c.satExcludeGaps ? t.accent : t.textFaint, width: 1.5),
                    ),
                    child: c.satExcludeGaps ? const Icon(FA.check, size: 9, color: Colors.white) : null,
                  ),
                  const SizedBox(width: 8),
                  Expanded(
                    child: Text(l10n.satParamGapsSwitch,
                        maxLines: 1,
                        overflow: TextOverflow.ellipsis,
                        style: TextStyle(fontSize: 13, color: c.satExcludeGaps ? t.accentText : t.text)),
                  ),
                ]),
              ),
            ),
          ),
        ]),
        if (c.toolsChecked && c.tools[ExternalTool.phykit] == null) ...[
          const SizedBox(height: 12),
          TonePanel(
            tone: PanelTone.warn,
            child: Row(children: [
              Expanded(child: Text(l10n.phykitRequired, style: const TextStyle(fontSize: 13))),
              const SizedBox(width: 10),
              ChipButton(
                label: l10n.installPhykit,
                icon: FA.download,
                style: ChipButtonStyle.accent,
                onPressed: c.busy || c.python == null ? null : c.installPhykit,
              ),
            ]),
          ),
        ],
        const SizedBox(height: 10),
        Align(
          alignment: Alignment.centerLeft,
          child: ChipButton(label: l10n.satParamsReset, icon: FA.rotateLeft, onPressed: c.busy ? null : c.resetSaturationParams),
        ),
      ]),
    );
  }
}
