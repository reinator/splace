import 'dart:math' as math;

import 'package:file_selector/file_selector.dart';
import 'package:flutter/material.dart';

import '../../core/genbank_markers.dart';
import '../../core/locus_organizer.dart';
import '../../core/marker_diagnostics.dart';
import '../../state/project_controller.dart';
import '../common.dart';
import '../fa.dart';
import '../theme.dart';
import '../widgets/gene_heatmap.dart';
import '../widgets/kit.dart';
import '../widgets/splace_tooltip.dart';

/// Step 4 — Select Markers (GenBank genes and FASTA loci), followed by the
/// Gene Presence and Marker Length heatmaps.
class MarkersSection extends StatelessWidget {
  const MarkersSection({super.key, required this.active});
  final bool active;

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    final markers = c.markers;
    final n = c.totalSamples;
    final selected = c.activeMarkers;
    final samples = c.samples;
    final org = c.organized;
    final cat = c.catalog;
    final hasGb = c.records.isNotEmpty && c.featureTypes.isNotEmpty;

    // Missing samples per selected marker.
    final missing = <String, List<SampleRef>>{};
    for (final name in selected) {
      final m = markers[name]!;
      final absent = [for (final s in samples) if (!m.lengths.containsKey(s.key)) s];
      if (absent.isNotEmpty) missing[name] = absent;
    }
    // GenBank duplicates (CDS/rRNA with more than one copy in a genome).
    final gbDup = <String, List<(int, List<GeneHit>)>>{};
    if (cat != null) {
      for (final type in const ['CDS', 'rRNA']) {
        if (!c.featureTypes.contains(type)) continue;
        for (final g in cat.byType[type]?.values ?? const <GeneSummary>[]) {
          if (!selected.contains(g.name)) continue;
          for (final e in g.byRecord.entries) {
            if (e.value.length > 1) (gbDup[g.name] ??= []).add((e.key, e.value));
          }
        }
      }
    }
    final faDup = <String, List<DuplicateGroup>>{};
    for (final d in org?.duplicates ?? const <DuplicateGroup>[]) {
      if (selected.contains(d.locus)) (faDup[d.locus] ??= []).add(d);
    }

    Widget legend(String pill, String text, ChipButtonStyle style) => Padding(
          padding: const EdgeInsets.symmetric(vertical: 3),
          child: Row(children: [
            IgnorePointer(child: ChipButton(label: pill, style: style, onPressed: () {})),
            const SizedBox(width: 8),
            Expanded(child: Text(text)),
          ]),
        );

    final data = HeatmapData(
      rows: [for (final s in samples) s.label],
      rowItalic: [for (final s in samples) s.italic],
      rowDetail: [for (final s in samples) s.detail],
      cols: markers.keys.toList(),
      values: [
        for (final s in samples)
          [
            for (final m in markers.values)
              () {
                final l = m.lengths[s.key];
                return l == null || l.isEmpty ? null : l.reduce(math.max);
              }(),
          ],
      ],
      colActive: [for (final m in markers.keys) selected.contains(m)],
    );
    final colIndex = {for (final (i, m) in data.cols.indexed) m: i};
    final selectedData = HeatmapData(
      rows: data.rows,
      rowItalic: data.rowItalic,
      rowDetail: data.rowDetail,
      cols: selected,
      values: [
        for (var r = 0; r < samples.length; r++) [for (final m in selected) data.values[r][colIndex[m]!]],
      ],
      colActive: [for (final _ in selected) true],
    );

    return Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
      StepSection(
        badge: '4',
        icon: FA.barcode,
        title: l10n.step4,
        count: selected.length,
        active: active,
        actions: [
          SplaceTooltip(message: l10n.tooltipAll, child: ChipButton(label: l10n.selectAll, onPressed: c.selectAllMarkers)),
          SplaceTooltip(message: l10n.tooltipNone, child: ChipButton(label: l10n.selectNone, onPressed: c.selectNoMarkers)),
          SplaceTooltip(message: l10n.tooltipShared, child: ChipButton(label: l10n.sharedOnly, style: ChipButtonStyle.soft, onPressed: c.selectSharedMarkers)),
          if (hasGb)
            SplaceTooltip(message: l10n.tooltipDefault, child: ChipButton(label: l10n.defaultPhylogeny, style: ChipButtonStyle.accent, onPressed: c.selectDefaultMarkers)),
        ],
        child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
          SubAnchor(id: '4a', child: SubHeading(id: '4a', title: l10n.subChooseMarkers)),
          TonePanel(
            tone: PanelTone.info,
            child: Row(crossAxisAlignment: CrossAxisAlignment.start, children: [
              Padding(padding: const EdgeInsets.only(top: 8), child: Icon(FA.circleInfo, size: 15, color: t.accent)),
              const SizedBox(width: 12),
              Expanded(
                child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
                  legend(l10n.selectAll, '— ${l10n.legendSelectAll}', ChipButtonStyle.neutral),
                  legend(l10n.selectNone, '— ${l10n.legendNone}', ChipButtonStyle.neutral),
                  legend(l10n.sharedOnly, '— ${l10n.legendShared}', ChipButtonStyle.soft),
                  if (hasGb) legend(l10n.defaultPhylogeny, '— ${l10n.legendDefault}', ChipButtonStyle.accent),
                ]),
              ),
            ]),
          ),
          const SizedBox(height: 14),
          const SizedBox(height: 6),
          if (markers.isEmpty)
            Text(l10n.noMarkers, style: TextStyle(fontStyle: FontStyle.italic, color: t.textFaint))
          else
            AutoGrid(minItemWidth: 230, spacing: 8, children: [
              for (final m in markers.values)
                SplaceTooltip(
                  message: _tooltip(context, m, n, org),
                  child: MarkerChip(
                    name: m.name,
                    active: selected.contains(m.name),
                    count: '${m.samples}/$n',
                    meta: m.meanLength > 0 ? '~${m.meanLength.round()} bp' : null,
                    onChanged: (v) => c.toggleMarker(m.name, v),
                  ),
                ),
            ]),
          if (missing.isNotEmpty) ...[
            const SizedBox(height: 16),
            TonePanel(
              tone: PanelTone.warn,
              title: l10n.missingTitle,
              icon: FA.circleExclamation,
              child: _MissingPanel(missing: missing, selected: selected.length, samples: samples),
            ),
          ],
          if (gbDup.isNotEmpty) ...[
            const SizedBox(height: 16),
            TonePanel(
              tone: PanelTone.duplicate,
              title: l10n.duplicatesTitle,
              icon: FA.clone,
              child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
                Text(
                  l10n.dupSummary('${gbDup.length}', '${gbDup.values.expand((l) => l).map((e) => e.$1).toSet().length}'),
                  style: const TextStyle(fontSize: 13.5, fontWeight: FontWeight.w600),
                ),
                const SizedBox(height: 4),
                Text(l10n.duplicateChoiceHelp, style: const TextStyle(fontSize: 13)),
                const SizedBox(height: 8),
                const _RelationLegend(),
                const SizedBox(height: 12),
                for (final e in gbDup.entries) _GbDuplicateBlock(gene: e.key, entries: e.value),
              ]),
            ),
          ],
          if (faDup.isNotEmpty) ...[
            const SizedBox(height: 16),
            TonePanel(
              tone: PanelTone.duplicate,
              title: '${l10n.duplicatesTitle} (FASTA)',
              icon: FA.clone,
              child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
                Text(l10n.dupFastaHelp, style: const TextStyle(fontSize: 13)),
                const SizedBox(height: 8),
                const _RelationLegend(),
                const SizedBox(height: 12),
                for (final e in faDup.entries.take(30)) _DuplicateBlock(locus: e.key, groups: e.value),
                if (faDup.length > 30) Text(l10n.andMore('${faDup.length - 30}')),
              ]),
            ),
          ],
          if (org != null && org.unassigned.isNotEmpty) ...[
            const SizedBox(height: 16),
            TonePanel(
              tone: PanelTone.danger,
              title: '${l10n.unassignedTitle} (${org.unassigned.length})',
              icon: FA.triangleExclamation,
              child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
                for (final u in org.unassigned.take(20))
                  Padding(
                    padding: const EdgeInsets.symmetric(vertical: 2),
                    child: Text.rich(TextSpan(children: [
                      TextSpan(text: u.header, style: const TextStyle(fontFamily: monoFont, fontSize: 12)),
                      TextSpan(text: '  · ${u.sourceFile} · ${[if (u.missingSample) l10n.missingSample, if (u.missingLocus) l10n.missingLocus].join(', ')}'),
                    ])),
                  ),
                if (org.unassigned.length > 20) Text(l10n.andMore('${org.unassigned.length - 20}')),
              ]),
            ),
          ],
          if (org != null) ...[
            const SizedBox(height: 18),
            Align(
              alignment: Alignment.centerLeft,
              child: OutlinedButton(
                onPressed: () async {
                  final dir = await getDirectoryPath();
                  if (dir != null) await c.exportLoci(targetDir: dir);
                },
                child: Text(l10n.exportLoci),
              ),
            ),
          ],
        ]),
      ),
      if (markers.isNotEmpty)
        _HeatmapPanel(
          badge: '4b',
          icon: FA.tableCells,
          title: l10n.heatmapTitle,
          summary: l10n.heatmapSummary('${markers.length}', '$n'),
          legend: [
            _Legend(color: GeneHeatmap.absent, label: l10n.heatmapAbsent),
            _Legend(color: GeneHeatmap.present, label: l10n.heatmapPresent),
          ],
          heatmap: () => GeneHeatmap(data: data, onToggleColumn: c.toggleMarker, hint: l10n.heatmapHoverHint),
        ),
      if (selected.isNotEmpty)
        _HeatmapPanel(
          badge: '4c',
          icon: FA.chartLine,
          title: l10n.lengthHeatmapTitle,
          summary: l10n.lengthHeatmapBody,
          legend: [
            Row(mainAxisSize: MainAxisSize.min, children: [
              Text('min', style: TextStyle(fontSize: 12, color: t.textMuted)),
              const SizedBox(width: 6),
              Container(
                width: 90,
                height: 10,
                decoration: BoxDecoration(
                  borderRadius: BorderRadius.circular(5),
                  gradient: const LinearGradient(colors: [GeneHeatmap.lengthLow, GeneHeatmap.lengthHigh]),
                ),
              ),
              const SizedBox(width: 6),
              Text('max', style: TextStyle(fontSize: 12, color: t.textMuted)),
            ]),
          ],
          heatmap: () => GeneHeatmap(data: selectedData, mode: HeatmapMode.length, hint: l10n.heatmapHoverHint),
        ),
    ]);
  }

  static String _tooltip(BuildContext context, MarkerSummary m, int n, OrganizerResult? org) {
    final l10n = context.l10n;
    final lengths = m.allLengths.toList();
    final aliases = org?.locusAliases[m.name];
    return [
      m.name,
      '${l10n.featureLabel}: ${m.sources.join(', ')}',
      '${l10n.recordsLabel}: ${m.samples}/$n',
      if (m.duplicated > 0) l10n.duplicatedIn('${m.duplicated}'),
      if (lengths.isNotEmpty) 'Avg ${m.meanLength.round()} bp · Min ${lengths.reduce(math.min)} bp · Max ${lengths.reduce(math.max)} bp',
      if (aliases != null && aliases.length > 1) aliases.join(', '),
    ].join('\n');
  }
}

class _Legend extends StatelessWidget {
  const _Legend({required this.color, required this.label});
  final Color color;
  final String label;

  @override
  Widget build(BuildContext context) => Row(mainAxisSize: MainAxisSize.min, children: [
        Text(label, style: TextStyle(fontSize: 13, fontWeight: FontWeight.w600, color: context.tokens.textMuted)),
        const SizedBox(width: 6),
        Container(width: 14, height: 14, decoration: BoxDecoration(color: color, borderRadius: BorderRadius.circular(3))),
      ]);
}

/// Short label and explanation of why a marker is missing.
(String, String) _reasonText(BuildContext context, MissingCause cause, String gene) {
  final l10n = context.l10n;
  final d = cause.detail ?? '';
  return switch (cause.reason) {
    MissingReason.typeNotSelected => (l10n.missTypeNotSelected(d), l10n.missTypeNotSelectedInfo(gene, d)),
    MissingReason.geneFeatureOnly => (l10n.missGeneOnly, l10n.missGeneOnlyInfo(gene, d)),
    MissingReason.noAnnotation => (l10n.missNoAnnotation(d), l10n.missNoAnnotationInfo(d)),
    MissingReason.unknownName => (l10n.missUnknownName, l10n.missUnknownNameInfo(d, gene)),
    MissingReason.partialGenome => (l10n.missPartial, l10n.missPartialInfo(d.split('|').first, d.split('|').last, gene)),
    MissingReason.notAnnotated => (l10n.missNotAnnotated, l10n.missNotAnnotatedInfo(gene)),
    MissingReason.notInFiles => (l10n.missNotInFiles, l10n.missNotInFilesInfo(gene)),
  };
}

/// Reasons the user can fix in SPLACE (warn colour) vs. data gaps (muted).
bool _fixable(MissingReason r) =>
    r == MissingReason.typeNotSelected || r == MissingReason.unknownName || r == MissingReason.geneFeatureOnly;

/// A small tag with the reason a marker is missing, explained on hover.
class _ReasonTag extends StatelessWidget {
  const _ReasonTag({required this.cause, required this.gene, this.prefix});
  final MissingCause cause;
  final String gene;
  final String? prefix;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final (label, info) = _reasonText(context, cause, gene);
    final color = _fixable(cause.reason) ? SplaceTokens.warn : t.textMuted;
    return SplaceTooltip(
      message: info,
      child: Container(
        padding: const EdgeInsets.symmetric(horizontal: 7, vertical: 2),
        decoration: BoxDecoration(
          color: color.withValues(alpha: 0.08),
          borderRadius: BorderRadius.circular(8),
          border: Border.all(color: color.withValues(alpha: 0.35)),
        ),
        child: Text(prefix == null ? label : '$prefix · $label', style: TextStyle(fontSize: 11.5, fontWeight: FontWeight.w600, color: color)),
      ),
    );
  }
}

/// Missing markers: a summary with what happens in the supermatrix, and
/// the gaps listed by marker or by sample, each with its likely reason.
class _MissingPanel extends StatefulWidget {
  const _MissingPanel({required this.missing, required this.selected, required this.samples});
  final Map<String, List<SampleRef>> missing;
  final int selected;
  final List<SampleRef> samples;

  @override
  State<_MissingPanel> createState() => _MissingPanelState();
}

class _MissingPanelState extends State<_MissingPanel> {
  /// The user's view; by default by sample when a few samples lack many
  /// markers (one partial genome would fill every marker's list).
  bool? _view;

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    final n = widget.samples.length;
    // sample key → missing markers.
    final bySample = <String, List<String>>{};
    for (final e in widget.missing.entries) {
      for (final s in e.value) {
        (bySample[s.key] ??= []).add(e.key);
      }
    }
    final causes = {
      for (final e in widget.missing.entries)
        for (final s in e.value) '${e.key}|${s.key}': c.missingCause(e.key, s),
    };
    final fixable = causes.values.where((x) => _fixable(x.reason)).length;
    final refs = {for (final s in widget.samples) s.key: s};
    final bySampleView = _view ?? bySample.length < widget.missing.length;

    Widget sampleName(SampleRef s) => Text.rich(
          TextSpan(children: [
            TextSpan(text: s.label.replaceAll('_', ' '), style: TextStyle(fontStyle: s.italic ? FontStyle.italic : FontStyle.normal, fontWeight: FontWeight.w600)),
            if (s.detail != null) TextSpan(text: '  ${s.detail}', style: TextStyle(fontSize: 12, color: t.textMuted)),
          ]),
          style: TextStyle(fontSize: 13, color: t.text),
          overflow: TextOverflow.ellipsis,
        );

    Widget coverage(int present, int total) => SizedBox(
          width: 90,
          child: RainbowProgressBar(value: total == 0 ? 0 : present / total, height: 7),
        );

    final byMarker = [
      for (final e in widget.missing.entries.take(40))
        Container(
          margin: const EdgeInsets.only(bottom: 8),
          padding: const EdgeInsets.fromLTRB(12, 10, 12, 10),
          decoration: BoxDecoration(color: t.surface, borderRadius: BorderRadius.circular(10), border: Border.all(color: t.border)),
          child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
            Row(children: [
              Text(e.key, style: TextStyle(fontWeight: FontWeight.w700, fontSize: 14, color: t.text)),
              const SizedBox(width: 10),
              coverage(n - e.value.length, n),
              const SizedBox(width: 8),
              Expanded(child: Text(l10n.missPresentIn('${n - e.value.length}', '$n', '${e.value.length}'), style: TextStyle(fontSize: 12.5, color: t.textMuted))),
            ]),
            const SizedBox(height: 8),
            for (final s in e.value.take(30))
              Padding(
                padding: const EdgeInsets.symmetric(vertical: 2),
                child: Row(children: [
                  Flexible(child: sampleName(s)),
                  const SizedBox(width: 8),
                  _ReasonTag(cause: causes['${e.key}|${s.key}']!, gene: e.key),
                ]),
              ),
            if (e.value.length > 30) Text(l10n.andMore('${e.value.length - 30}'), style: TextStyle(fontSize: 12.5, color: t.textMuted)),
          ]),
        ),
      if (widget.missing.length > 40) Text(l10n.andMore('${widget.missing.length - 40}')),
    ];

    final sampleRows = bySample.entries.toList()..sort((a, b) => b.value.length.compareTo(a.value.length));
    final sampleList = [
      for (final e in sampleRows.take(60))
        Container(
          margin: const EdgeInsets.only(bottom: 6),
          padding: const EdgeInsets.fromLTRB(12, 8, 12, 8),
          decoration: BoxDecoration(color: t.surface, borderRadius: BorderRadius.circular(10), border: Border.all(color: t.border)),
          child: Row(crossAxisAlignment: CrossAxisAlignment.start, children: [
            Expanded(
              flex: 4,
              child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
                sampleName(refs[e.key]!),
                const SizedBox(height: 4),
                Row(children: [
                  coverage(widget.selected - e.value.length, widget.selected),
                  const SizedBox(width: 8),
                  Flexible(
                    child: Text(l10n.missSampleCount('${e.value.length}', '${widget.selected}'),
                        style: TextStyle(fontSize: 12, color: e.value.length * 2 > widget.selected ? SplaceTokens.warn : t.textMuted)),
                  ),
                ]),
              ]),
            ),
            const SizedBox(width: 12),
            Expanded(
              flex: 7,
              child: Wrap(spacing: 6, runSpacing: 6, children: [
                for (final g in e.value) _ReasonTag(cause: causes['$g|${e.key}']!, gene: g, prefix: g),
              ]),
            ),
          ]),
        ),
      if (sampleRows.length > 60) Text(l10n.andMore('${sampleRows.length - 60}')),
    ];

    return Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
      Text(l10n.missSummary('${widget.missing.length}', '${widget.selected}', '${bySample.length}', '$n'),
          style: const TextStyle(fontSize: 13.5, fontWeight: FontWeight.w600)),
      const SizedBox(height: 4),
      Text(l10n.missMatrixNote, style: const TextStyle(fontSize: 13)),
      if (fixable > 0) ...[
        const SizedBox(height: 4),
        Text(l10n.missFixableNote('$fixable'), style: const TextStyle(fontSize: 13, fontWeight: FontWeight.w600, color: SplaceTokens.warn)),
      ],
      const SizedBox(height: 10),
      Row(children: [
        Seg<bool>(
          value: bySampleView,
          options: [(false, l10n.missByMarker), (true, l10n.missBySample)],
          onChanged: (v) => setState(() => _view = v),
        ),
        const SizedBox(width: 10),
        Expanded(child: Text(l10n.missHoverHint, style: TextStyle(fontSize: 12, color: t.textMuted))),
      ]),
      const SizedBox(height: 10),
      ...(bySampleView ? sampleList : byMarker),
    ]);
  }
}

/// Label and explanation of how the copies of a duplicated marker relate.
(String, String, Color) _relationText(BuildContext context, CopyRelation r) {
  final l10n = context.l10n;
  final t = context.tokens;
  return switch (r) {
    CopyRelation.identical => (l10n.dupIdentical, l10n.dupIdenticalInfo, t.accent),
    CopyRelation.contained => (l10n.dupContained, l10n.dupContainedInfo, SplaceTokens.warn),
    CopyRelation.different => (l10n.dupDifferent, l10n.dupDifferentInfo, SplacePalette.red500),
  };
}

class _RelationTag extends StatelessWidget {
  const _RelationTag(this.relation);
  final CopyRelation relation;

  @override
  Widget build(BuildContext context) {
    final (label, info, color) = _relationText(context, relation);
    return SplaceTooltip(message: info, child: TagChip(label, color: color));
  }
}

/// What "identical", "fragment" and "different" copies mean.
class _RelationLegend extends StatelessWidget {
  const _RelationLegend();

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    return Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
      for (final r in CopyRelation.values)
        Padding(
          padding: const EdgeInsets.symmetric(vertical: 2),
          child: Row(crossAxisAlignment: CrossAxisAlignment.start, children: [
            SizedBox(width: 210, child: Align(alignment: Alignment.centerLeft, child: _RelationTag(r))),
            Expanded(child: Text(_relationText(context, r).$2, style: TextStyle(fontSize: 12.5, color: t.textMuted))),
          ]),
        ),
    ]);
  }
}

/// GenBank genes with more than one copy in a genome: the user picks the
/// copy to use (the longest starts selected).
class _GbDuplicateBlock extends StatelessWidget {
  const _GbDuplicateBlock({required this.gene, required this.entries});
  final String gene;
  final List<(int, List<GeneHit>)> entries;

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    final dark = Theme.of(context).brightness == Brightness.dark;
    final head = TextStyle(fontSize: 12.5, fontWeight: FontWeight.w700, color: dark ? null : SplacePalette.red900);

    // Reference: the same gene in genomes with a single copy.
    final single = <int>[];
    for (final type in const ['CDS', 'rRNA']) {
      for (final hits in c.catalog?.byType[type]?[gene]?.byRecord.values ?? const <List<GeneHit>>[]) {
        if (hits.length == 1 && hits.first.length > 0) single.add(hits.first.length);
      }
    }
    final avg = single.isEmpty ? null : single.reduce((a, b) => a + b) / single.length;
    final lo = single.isEmpty ? 0 : single.reduce(math.min), hi = single.isEmpty ? 0 : single.reduce(math.max);

    return Padding(
      padding: const EdgeInsets.only(bottom: 14),
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        Row(children: [
          const Icon(FA.clone, size: 13, color: SplacePalette.red500),
          const SizedBox(width: 6),
          Text(gene, style: const TextStyle(fontWeight: FontWeight.w600, fontSize: 14)),
          const SizedBox(width: 6),
          Text('— ${l10n.duplicatedIn('${entries.length}')}', style: const TextStyle(fontSize: 13)),
        ]),
        const SizedBox(height: 4),
        Text(
          avg == null ? l10n.dupNoReference : l10n.dupReference(gene, '${single.length}', avg.round().toString(), '$lo', '$hi'),
          style: TextStyle(fontSize: 12.5, color: t.textMuted),
        ),
        const SizedBox(height: 6),
        Container(
          decoration: BoxDecoration(border: Border.all(color: dark ? const Color(0x33FDAE61) : SplacePalette.orange)),
          child: Column(children: [
            Container(
              color: dark ? const Color(0x1AFDAE61) : const Color(0x33FDAE61),
              padding: const EdgeInsets.symmetric(horizontal: 12, vertical: 8),
              child: Row(children: [
                Expanded(flex: 3, child: Text(l10n.colSpecies, style: head)),
                Expanded(flex: 2, child: Text(l10n.colAccession, style: head)),
                Expanded(flex: 7, child: Text(l10n.dupCopiesColumn, style: head)),
              ]),
            ),
            for (final (ri, hits) in entries)
              () {
                final r = c.records[ri];
                final chosen = chosenHit(hits, c.duplicateChoices[duplicateChoiceKey(r, gene)]);
                final longest = chosenHit(hits, null);
                // Copy whose length is closest to the other genomes.
                int? closest;
                if (avg != null) {
                  var best = double.infinity;
                  for (var i = 0; i < hits.length; i++) {
                    final d = (hits[i].length - avg).abs();
                    if (d < best) (best, closest) = (d, i);
                  }
                }
                return Padding(
                  padding: const EdgeInsets.symmetric(horizontal: 12, vertical: 6),
                  child: Row(crossAxisAlignment: CrossAxisAlignment.start, children: [
                    Expanded(
                      flex: 3,
                      child: Padding(
                        padding: const EdgeInsets.only(top: 8),
                        child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
                          Text(r.organismName, style: const TextStyle(fontStyle: FontStyle.italic, fontSize: 13), overflow: TextOverflow.ellipsis),
                          const SizedBox(height: 4),
                          _RelationTag(c.copyRelationOf(ri, hits)),
                        ]),
                      ),
                    ),
                    Expanded(flex: 2, child: Padding(padding: const EdgeInsets.only(top: 8), child: Text(r.accession, style: const TextStyle(fontSize: 13)))),
                    Expanded(
                      flex: 7,
                      child: Wrap(spacing: 8, runSpacing: 6, children: [
                        for (var i = 0; i < hits.length; i++)
                          _CopyOption(
                            length: hits[i].length,
                            location: r.features[hits[i].feature].location,
                            diff: avg == null ? null : (hits[i].length - avg).round(),
                            closest: i == closest,
                            longest: i == longest,
                            selected: i == chosen,
                            onTap: () => c.setDuplicateChoice(r, gene, i),
                          ),
                      ]),
                    ),
                  ]),
                );
              }(),
          ]),
        ),
      ]),
    );
  }
}

/// One copy of a duplicated gene: length, difference to the average of the
/// other genomes, location, and tags for the longest / closest copy.
class _CopyOption extends StatelessWidget {
  const _CopyOption({
    required this.length,
    required this.location,
    required this.diff,
    required this.closest,
    required this.longest,
    required this.selected,
    required this.onTap,
  });
  final int length;
  final String location;
  final int? diff;
  final bool closest, longest, selected;
  final VoidCallback onTap;

  @override
  Widget build(BuildContext context) {
    final l10n = context.l10n;
    final t = context.tokens;
    final d = diff;
    return InkWell(
      borderRadius: BorderRadius.circular(10),
      onTap: onTap,
      child: Container(
        padding: const EdgeInsets.symmetric(horizontal: 10, vertical: 6),
        decoration: BoxDecoration(
          color: selected ? t.accentSoft : t.surface,
          borderRadius: BorderRadius.circular(10),
          border: Border.all(color: selected ? t.accent : t.border, width: selected ? 2 : 1),
        ),
        child: Column(crossAxisAlignment: CrossAxisAlignment.start, mainAxisSize: MainAxisSize.min, children: [
          Row(mainAxisSize: MainAxisSize.min, children: [
            Text('$length bp', style: TextStyle(fontSize: 13, fontWeight: FontWeight.w700, color: selected ? t.accentText : t.text)),
            if (d != null) ...[
              const SizedBox(width: 6),
              Text(
                d == 0 ? l10n.dupSameAsAverage : l10n.dupVsAverage('${d > 0 ? '+' : '−'}${d.abs()}'),
                style: TextStyle(fontSize: 12, color: d.abs() <= 15 ? t.accentText : SplaceTokens.warn),
              ),
            ],
          ]),
          Text(location, maxLines: 1, overflow: TextOverflow.ellipsis, style: TextStyle(fontSize: 11.5, color: t.textMuted)),
          if (closest || longest)
            Padding(
              padding: const EdgeInsets.only(top: 3),
              child: Wrap(spacing: 4, children: [
                if (closest) TagChip(l10n.dupClosest),
                if (longest) TagChip(l10n.dupLongestTag, color: t.textMuted),
              ]),
            ),
        ]),
      ),
    );
  }
}

/// Duplicated sample/locus pairs (`.duplicate-genes` table).
class _DuplicateBlock extends StatelessWidget {
  const _DuplicateBlock({required this.locus, required this.groups});
  final String locus;
  final List<DuplicateGroup> groups;

  @override
  Widget build(BuildContext context) {
    final l10n = context.l10n;
    final dark = Theme.of(context).brightness == Brightness.dark;
    final head = TextStyle(fontSize: 12.5, fontWeight: FontWeight.w700, color: dark ? null : SplacePalette.red900);
    Widget row(List<Widget> cells, {bool header = false}) => Container(
          padding: const EdgeInsets.symmetric(horizontal: 12, vertical: 8),
          color: header ? (dark ? const Color(0x1AFDAE61) : const Color(0x33FDAE61)) : null,
          child: Row(children: [
            Expanded(flex: 3, child: cells[0]),
            Expanded(flex: 3, child: cells[1]),
            Expanded(flex: 2, child: cells[2]),
            Expanded(flex: 3, child: cells[3]),
          ]),
        );
    return Padding(
      padding: const EdgeInsets.only(bottom: 12),
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        Row(children: [
          const Icon(FA.clone, size: 13, color: SplacePalette.red500),
          const SizedBox(width: 6),
          Text(locus, style: const TextStyle(fontWeight: FontWeight.w600, fontSize: 14)),
          const SizedBox(width: 6),
          Text('— ${l10n.duplicatedIn('${groups.length}')}', style: const TextStyle(fontSize: 13)),
        ]),
        const SizedBox(height: 6),
        Container(
          decoration: BoxDecoration(border: Border.all(color: dark ? const Color(0x33FDAE61) : SplacePalette.orange)),
          child: Column(children: [
            row([
              Text(l10n.colSample, style: head),
              Text(l10n.colCopies, style: head),
              Text(l10n.colUsed, style: head),
              Text(l10n.colFile, style: head),
            ], header: true),
            for (final g in groups.take(8))
              row([
                Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
                  Text(g.sample, style: const TextStyle(fontStyle: FontStyle.italic, fontSize: 13), overflow: TextOverflow.ellipsis),
                  const SizedBox(height: 4),
                  _RelationTag(copyRelation([for (final x in g.candidates) x.record.sequence.replaceAll('-', '')])),
                ]),
                Text(g.candidates.map((x) => '${x.record.length} bp').join(', '), style: const TextStyle(fontSize: 13)),
                Text(g.kept.map((x) => '${x.record.length} bp').join(', '), style: const TextStyle(fontWeight: FontWeight.w700, fontSize: 13)),
                Text(g.candidates.first.sourceFile, style: const TextStyle(fontSize: 12), overflow: TextOverflow.ellipsis),
              ]),
            if (groups.length > 8)
              Padding(padding: const EdgeInsets.all(8), child: Text(l10n.andMore('${groups.length - 8}'))),
          ]),
        ),
      ]),
    );
  }
}

/// A heatmap section that stays collapsed until opened (large datasets
/// make long heatmaps), with its legend and a one-line summary.
class _HeatmapPanel extends StatefulWidget {
  const _HeatmapPanel({required this.badge, required this.icon, required this.title, required this.summary, required this.legend, required this.heatmap});
  final String badge;
  final IconData icon;
  final String title, summary;
  final List<Widget> legend;
  final Widget Function() heatmap;

  @override
  State<_HeatmapPanel> createState() => _HeatmapPanelState();
}

class _HeatmapPanelState extends State<_HeatmapPanel> {
  bool _open = false;

  @override
  Widget build(BuildContext context) {
    final l10n = context.l10n;
    final t = context.tokens;
    return SubAnchor(
      id: widget.badge,
      child: StepSection(
      badge: widget.badge,
      icon: widget.icon,
      title: widget.title,
      actions: [
        if (_open) ...widget.legend,
        ChipButton(
          label: _open ? l10n.heatmapHide : l10n.heatmapShow,
          style: _open ? ChipButtonStyle.neutral : ChipButtonStyle.accent,
          onPressed: () => setState(() => _open = !_open),
        ),
      ],
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        Text(widget.summary, style: TextStyle(fontSize: 13, color: t.textMuted)),
        if (_open) ...[
          const SizedBox(height: 10),
          widget.heatmap(),
        ],
      ]),
    ),
    );
  }
}
