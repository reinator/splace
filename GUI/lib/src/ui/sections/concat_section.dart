
import 'package:file_selector/file_selector.dart';
import 'package:flutter/material.dart';
import 'package:path/path.dart' as p;

import '../../core/supermatrix.dart';
import '../common.dart';
import '../fa.dart';
import '../theme.dart';
import '../widgets/examples.dart';
import '../widgets/kit.dart';
import '../widgets/matrix_table.dart';

/// Step 7 — Concatenate & Generate NEXUS (`#concatSection`).
class ConcatSection extends StatefulWidget {
  const ConcatSection({super.key, required this.active});
  final bool active;

  @override
  State<ConcatSection> createState() => _ConcatSectionState();
}

class _ConcatSectionState extends State<ConcatSection> {
  String _outgroupFilter = '';

  Future<void> _save(String key, String suggested) async {
    final c = ProjectScope.read(context);
    final loc = await getSaveLocation(suggestedName: suggested);
    if (loc != null) await c.saveMatrixFile(key, loc.path);
  }

  @override
  Widget build(BuildContext context) {
    final active = widget.active;
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final o = c.concat;
    final m = c.matrix;
    final anyTrimmed = c.loci.values.any((l) => l.trimmed != null);
    final includedLoci = c.loci.values.where((l) => l.included && l.aligned != null).length;
    // The markers of this project, for the partition pictures.
    final projectLoci = [
      for (final l in c.loci.values)
        if (l.included && l.aligned != null) (name: l.name, coding: l.coding && !l.isProtein && c.codonPositionsUsable(l)),
    ];

    ConcatOptions copy({PartitionScheme? scheme, Set<int>? exclude, Map<String, Set<int>>? perLocus, int? minLoci}) =>
        o.copyWith(scheme: scheme, excludePositions: exclude, perLocusExclude: perLocus, minLociPerTaxon: minLoci);
    final t = context.tokens;
    final taxa = <String>{
      for (final l in c.loci.values)
        if (l.included) ...(l.alignment ?? l.raw).map((r) => r.header),
    }.toList()
      ..sort();
    final q = _outgroupFilter.toLowerCase();
    final filtered = q.isEmpty ? taxa : taxa.where((x) => x.toLowerCase().contains(q)).toList();
    void toggleOutgroup(String x) => c.setConcat(o.copyWith(outgroups: o.outgroups.contains(x) ? (List.of(o.outgroups)..remove(x)) : [...o.outgroups, x]));

    final source = c.useTrimmed && anyTrimmed ? l10n.srcTrimmed : l10n.srcAligned;
    // Only the schemes that mean something for these alignments.
    final schemes = c.partitionSchemes;
    final scheme = c.effectiveScheme;
    final anyCodons = projectLoci.any((l) => l.coding);
    final schemeName = switch (scheme) {
      PartitionScheme.byLocusAndCodon => l10n.schemeByLocusCodon,
      PartitionScheme.byLocus => l10n.schemeByLocus,
      PartitionScheme.byCodonPosition => l10n.schemeByCodon,
      PartitionScheme.none => l10n.schemeNone,
    };
    final excluded = o.perLocusExclude.isNotEmpty
        ? l10n.perLocusExclusions(perLocusSummary(context, o.perLocusExclude))
        : (o.excludePositions.isEmpty ? l10n.concatNothingExcluded : codonPositions(context, o.excludePositions));

    return StepSection(
      badge: '6',
      title: l10n.step7,
      description: l10n.step7Intro,
      active: active,
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        _HowItWorks(loci: [for (final l in c.loci.values) if (l.included && l.aligned != null) l.name]),
        const SizedBox(height: 16),
        // 1. Which alignments.
        _Block(
          badge: '6a',
          title: l10n.concatBlockSource,
          body: l10n.concatBlockSourceBody,
          child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
            AutoGrid(minItemWidth: 300, children: [
              RadioCard<bool>(
                value: true,
                group: c.useTrimmed && anyTrimmed,
                onChanged: anyTrimmed ? c.setUseTrimmed : (_) {},
                title: '${l10n.srcTrimmed} ${l10n.recommendedTag}',
                subtitle: l10n.srcTrimmedBody,
              ),
              RadioCard<bool>(
                value: false,
                group: c.useTrimmed && anyTrimmed,
                onChanged: c.setUseTrimmed,
                title: l10n.srcAligned,
                subtitle: l10n.srcAlignedBody,
              ),
            ]),
            const SizedBox(height: 8),
            CheckCard(
              value: o.replaceGapsWithMissing,
              onChanged: (v) => c.setConcat(o.copyWith(replaceGapsWithMissing: v)),
              title: l10n.replaceGaps,
              subtitle: l10n.replaceGapsBody,
            ),
          ]),
        ),
        const SizedBox(height: 12),
        // 2. Partitions and codon positions.
        _Block(
          badge: '6b',
          title: l10n.concatBlockPartitions,
          body: l10n.concatBlockPartitionsBody,
          child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
            // Why some schemes are not offered (one marker, or no codon
            // positions known).
            if (schemes.length < 4) ...[
              TonePanel(
                tone: PanelTone.info,
                icon: FA.circleInfo,
                child: Text(
                  includedLoci <= 1
                      ? (anyCodons ? l10n.schemesOneCoding : l10n.schemesOneNoCodons)
                      : l10n.schemesNoCodons,
                  style: const TextStyle(fontSize: 13),
                ),
              ),
              const SizedBox(height: 10),
            ],
            for (final (s, title, body) in [
              (PartitionScheme.byLocusAndCodon, l10n.schemeByLocusCodon, l10n.schemeByLocusCodonBody),
              (PartitionScheme.byLocus, l10n.schemeByLocus, l10n.schemeByLocusBody),
              (PartitionScheme.byCodonPosition, l10n.schemeByCodon, l10n.schemeByCodonBody),
              (PartitionScheme.none, l10n.schemeNone, l10n.schemeNoneBody),
            ])
              if (schemes.contains(s))
                Padding(
                  padding: const EdgeInsets.only(bottom: 8),
                  child: _SchemeOption(
                    selected: scheme == s,
                    title: s == schemes.first ? '$title ${l10n.recommendedTag}' : title,
                    body: body,
                    onTap: () => c.setConcat(copy(scheme: s)),
                    picture: _SchemePicture(scheme: s, loci: projectLoci),
                  ),
                ),
            if (anyCodons) ...[
            const SizedBox(height: 14),
            Row(children: [
              Text(l10n.codonExclusion, style: TextStyle(fontSize: 14, fontWeight: FontWeight.w600, color: t.text)),
              const SizedBox(width: 6),
              InfoDot(l10n.codonExclusionInfo),
            ]),
            const SizedBox(height: 2),
            Text(l10n.codonExclusionHelp, style: TextStyle(fontSize: 12.5, color: t.textMuted)),
            const SizedBox(height: 8),
            Wrap(spacing: 8, runSpacing: 8, children: [
              for (final pos in [1, 2, 3])
                MarkerChipLike(
                  label: l10n.excludePosition('$pos'),
                  selected: o.excludePositions.contains(pos),
                  onTap: () => c.setConcat(copy(
                    exclude: o.excludePositions.contains(pos) ? ({...o.excludePositions}..remove(pos)) : {...o.excludePositions, pos},
                  )),
                ),
              if (o.perLocusExclude.isNotEmpty)
                ChipButton(
                  label: l10n.perLocusExclusions(perLocusSummary(context, o.perLocusExclude)),
                  icon: FA.xmark,
                  style: ChipButtonStyle.soft,
                  tooltip: o.perLocusExclude.length <= 3
                      ? l10n.clearExclusions
                      : '${l10n.clearExclusions}\n\n${o.perLocusExclude.entries.take(40).map((e) => '${e.key}: ${e.value.join('+')}').join('\n')}'
                          '${o.perLocusExclude.length > 40 ? '\n…' : ''}',
                  onPressed: () => c.setConcat(copy(perLocus: const {})),
                ),
            ]),
            const SizedBox(height: 10),
            CodonExclusionExample(excluded: o.excludePositions),
            ],
          ]),
        ),
        const SizedBox(height: 16),
        // 3. Outgroup, next to the partitions.
        _Block(
          badge: '6c',
          title: l10n.concatBlockOutgroup,
          body: l10n.concatBlockOutgroupBody,
          child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
            _OutgroupSuggestionPanel(onUse: (taxa) => c.setConcat(o.copyWith(outgroups: taxa))),
            const SizedBox(height: 10),
            Container(
              padding: const EdgeInsets.all(10),
              decoration: BoxDecoration(color: t.surface, borderRadius: BorderRadius.circular(10), border: Border.all(color: t.border)),
              child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
                TextField(
                  decoration: InputDecoration(hintText: l10n.outgroupHint, prefixIcon: Icon(FA.magnifyingGlass, size: 13, color: t.textFaint)),
                  onChanged: (v) => setState(() => _outgroupFilter = v),
                ),
                const SizedBox(height: 8),
                ConstrainedBox(
                  constraints: const BoxConstraints(maxHeight: 160),
                  child: SingleChildScrollView(
                    child: Wrap(spacing: 6, runSpacing: 6, children: [
                      // A few chips at a time: hundreds of taxa are found by typing.
                      for (final x in filtered.take(60))
                        FilterChip(
                          label: HighlightText(x, _outgroupFilter, style: const TextStyle(fontSize: 12, fontStyle: FontStyle.italic)),
                          selected: o.outgroups.contains(x),
                          onSelected: (_) => toggleOutgroup(x),
                          visualDensity: VisualDensity.compact,
                        ),
                    ]),
                  ),
                ),
                if (filtered.length > 60) ...[
                  const SizedBox(height: 6),
                  Text(l10n.resultsShowing('60', '${filtered.length}'), style: TextStyle(fontSize: 12, color: t.textMuted)),
                ],
                if (o.outgroups.isNotEmpty) ...[
                  const SizedBox(height: 8),
                  Text('${l10n.selectedLabel}: ${o.outgroups.join(', ')}', style: TextStyle(fontSize: 12, color: t.accentText, fontWeight: FontWeight.w600)),
                ],
              ]),
            ),
          ]),
        ),
        const SizedBox(height: 12),
        // What will be built, right above the button.
        SubAnchor(
          id: '6d',
          child: Container(
          padding: const EdgeInsets.all(14),
          decoration: BoxDecoration(color: t.accentSoft, borderRadius: BorderRadius.circular(12), border: Border.all(color: t.accentSoftBorder)),
          child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
            SubHeading(id: '6d', title: l10n.subBuildMatrix, body: l10n.concatSummaryTitle),
            Wrap(spacing: 18, runSpacing: 6, children: [
              for (final (k, v) in [
                (l10n.concatSummaryLoci, '$includedLoci'),
                (l10n.concatSummarySource, source),
                (l10n.outgroup, o.outgroups.isEmpty ? l10n.concatNoOutgroup : o.outgroups.join(', ')),
                (l10n.partitionScheme, schemeName),
                (l10n.excludedLabel, excluded),
              ])
                Text.rich(
                  TextSpan(children: [
                    TextSpan(text: '$k: ', style: TextStyle(color: t.textMuted)),
                    TextSpan(text: v, style: TextStyle(fontWeight: FontWeight.w600, color: t.text)),
                  ]),
                  style: const TextStyle(fontSize: 13),
                ),
            ]),
            const SizedBox(height: 12),
            Align(
              alignment: Alignment.centerLeft,
              child: PrimaryButton(label: l10n.buildMatrix, icon: FA.puzzlePiece, onPressed: c.busy ? null : c.buildMatrix),
            ),
          ]),
        ),
        ),
        if (m != null) ...[
          const SizedBox(height: 16),
          TonePanel(
            tone: PanelTone.success,
            title: l10n.matrixReady,
            icon: FA.circleCheck,
            child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
              Text(l10n.matrixStats('${m.taxa.length}', '${m.length}', '${m.partitions.length}', pct(m.missingFraction)),
                  style: const TextStyle(fontWeight: FontWeight.w600)),
              const SizedBox(height: 4),
              Text('${l10n.excludedLabel}: ${c.describeExclusions()}', style: const TextStyle(fontSize: 13)),
              const SizedBox(height: 10),
              Align(
                alignment: Alignment.centerLeft,
                child: ChipButton(
                  label: l10n.openFolder,
                  icon: FA.folderOpen,
                  onPressed: () => c.openFolder(p.dirname(c.matrixFiles.values.first)),
                ),
              ),
            ]),
          ),
          const SizedBox(height: 12),
          Wrap(spacing: 8, runSpacing: 8, children: [
            ChipButton(label: l10n.saveNexus, icon: FA.fileArrowDown, style: ChipButtonStyle.accent, onPressed: () => _save('nexus', 'concatenated.nex')),
            ChipButton(label: l10n.savePartition, icon: FA.fileArrowDown, onPressed: () => _save('partition', 'partitions.txt')),
            ChipButton(label: l10n.saveFasta, icon: FA.fileArrowDown, onPressed: () => _save('fasta', 'concatenated.fasta')),
            ChipButton(label: l10n.savePhylip, icon: FA.fileArrowDown, onPressed: () => _save('phylip', 'concatenated.phy')),
          ]),
          const SizedBox(height: 12),
          PreviewBox(label: l10n.outputFiles, text: c.matrixFiles.values.join('\n'), maxHeight: 120),
          if (c.overview != null) ...[
            const SizedBox(height: 16),
            Kicker(l10n.mtTitle),
            const SizedBox(height: 4),
            Text(l10n.mtIntro, style: TextStyle(fontSize: 12.5, color: t.textMuted, height: 1.4)),
            const SizedBox(height: 8),
            MatrixOverviewTable(overview: c.overview!),
          ],
          const SizedBox(height: 16),
          Kicker(l10n.mtPartitions),
          const SizedBox(height: 8),
          _PartitionTable(matrix: m),
        ],
      ]),
    );
  }
}

/// One lettered group of options ("7a"): title, what it decides, options.
class _Block extends StatelessWidget {
  const _Block({required this.badge, required this.title, required this.body, required this.child});
  final String badge;
  final String title, body;
  final Widget child;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    return SubAnchor(
      id: badge,
      child: Container(
      padding: const EdgeInsets.all(16),
      decoration: BoxDecoration(color: t.surface2, borderRadius: BorderRadius.circular(14), border: Border.all(color: t.border)),
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        Row(crossAxisAlignment: CrossAxisAlignment.start, children: [
          StepBadge(badge, size: 26),
          const SizedBox(width: 10),
          Expanded(
            child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
              Text(title, style: TextStyle(fontSize: 15.5, fontWeight: FontWeight.w700, color: t.text)),
              const SizedBox(height: 2),
              Text(body, style: TextStyle(fontSize: 13, color: t.textMuted, height: 1.4)),
            ]),
          ),
        ]),
        const SizedBox(height: 12),
        child,
      ]),
    ),
    );
  }
}

/// How the supermatrix is built: the loci side by side, one row per taxon,
/// "?" where a taxon lacks a locus, and partitions marking each locus.
/// "How the supermatrix is built": one line that opens the explainer, so
/// the options of step 7 follow the suggestion of step 6 closely.
class _HowItWorks extends StatefulWidget {
  const _HowItWorks({required this.loci});
  final List<String> loci;

  @override
  State<_HowItWorks> createState() => _HowItWorksState();
}

class _HowItWorksState extends State<_HowItWorks> {
  bool _open = false;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final l10n = context.l10n;
    if (_open) return _HowItWorksBody(loci: widget.loci, onClose: () => setState(() => _open = false));
    return Material(
      color: t.surface,
      shape: RoundedRectangleBorder(borderRadius: BorderRadius.circular(12), side: BorderSide(color: t.accentSoftBorder)),
      child: InkWell(
        borderRadius: BorderRadius.circular(12),
        onTap: () => setState(() => _open = true),
        child: Padding(
          padding: const EdgeInsets.symmetric(horizontal: 16, vertical: 11),
          child: Row(children: [
            Icon(FA.puzzlePiece, size: 13, color: t.accent),
            const SizedBox(width: 8),
            Text(l10n.concatHowTitle, style: TextStyle(fontSize: 14, fontWeight: FontWeight.w600, color: t.text)),
            const SizedBox(width: 10),
            Expanded(child: Text(l10n.concatHowOpen, overflow: TextOverflow.ellipsis, style: TextStyle(fontSize: 12.5, color: t.textMuted))),
            Icon(FA.chevronDown, size: 11, color: t.textMuted),
          ]),
        ),
      ),
    );
  }
}

class _HowItWorksBody extends StatelessWidget {
  const _HowItWorksBody({required this.loci, required this.onClose});
  final List<String> loci;
  final VoidCallback onClose;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final l10n = context.l10n;
    const colors = SplacePalette.categorical;
    final names = loci.isEmpty ? ['COI', 'CYTB', '16S'] : loci.take(5).toList();
    Widget block(int i, {bool missing = false}) => Expanded(
          flex: 3 + (names[i].length % 3),
          child: Container(
            height: 18,
            margin: const EdgeInsets.only(right: 3),
            alignment: Alignment.center,
            decoration: BoxDecoration(
              color: missing ? t.surface3 : colors[i % colors.length].withValues(alpha: 0.75),
              borderRadius: BorderRadius.circular(4),
            ),
            child: Text(missing ? '? ? ?' : '', style: TextStyle(fontSize: 10, color: t.textMuted, fontWeight: FontWeight.w700)),
          ),
        );
    Widget row(String taxon, {int? missing}) => Padding(
          padding: const EdgeInsets.only(bottom: 4),
          child: Row(children: [
            SizedBox(width: 74, child: Text(taxon, style: TextStyle(fontSize: 11.5, fontStyle: FontStyle.italic, color: t.textMuted))),
            for (var i = 0; i < names.length; i++) block(i, missing: i == missing),
          ]),
        );
    return Container(
      padding: const EdgeInsets.all(16),
      decoration: BoxDecoration(color: t.surface, borderRadius: BorderRadius.circular(14), border: Border.all(color: t.accentSoftBorder), boxShadow: shadowSm),
      child: LayoutBuilder(builder: (context, box) {
        final picture = Column(crossAxisAlignment: CrossAxisAlignment.stretch, mainAxisSize: MainAxisSize.min, children: [
          Row(children: [
            const SizedBox(width: 74),
            for (var i = 0; i < names.length; i++)
              Expanded(
                flex: 3 + (names[i].length % 3),
                child: Text(names[i], overflow: TextOverflow.ellipsis, style: TextStyle(fontSize: 11, fontWeight: FontWeight.w700, color: colors[i % colors.length])),
              ),
          ]),
          const SizedBox(height: 4),
          row('Taxon A'),
          row('Taxon B', missing: names.length > 1 ? 1 : null),
          row('Taxon C'),
          const SizedBox(height: 2),
          Row(children: [
            const SizedBox(width: 74),
            Expanded(child: Text(l10n.concatPictureCaption, style: TextStyle(fontSize: 11, color: t.textFaint))),
          ]),
        ]);
        final text = Column(crossAxisAlignment: CrossAxisAlignment.start, mainAxisSize: MainAxisSize.min, children: [
          Row(children: [
            Icon(FA.puzzlePiece, size: 14, color: t.accent),
            const SizedBox(width: 8),
            Expanded(child: Text(l10n.concatHowTitle, style: TextStyle(fontSize: 16, fontWeight: FontWeight.w700, color: t.text))),
            IconButton(
              tooltip: l10n.close,
              visualDensity: VisualDensity.compact,
              onPressed: onClose,
              icon: Icon(FA.chevronUp, size: 12, color: t.textMuted),
            ),
          ]),
          const SizedBox(height: 8),
          for (final (i, s) in [l10n.concatHow1, l10n.concatHow2, l10n.concatHow3].indexed)
            Padding(
              padding: const EdgeInsets.only(bottom: 5),
              child: Row(crossAxisAlignment: CrossAxisAlignment.start, children: [
                Text('${i + 1}. ', style: TextStyle(fontSize: 13, fontWeight: FontWeight.w700, color: t.accentText)),
                Expanded(child: Text(s, style: TextStyle(fontSize: 13, color: t.textMuted, height: 1.4))),
              ]),
            ),
        ]);
        if (box.maxWidth < 760) return Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [text, const SizedBox(height: 12), picture]);
        return Row(crossAxisAlignment: CrossAxisAlignment.center, children: [
          Expanded(flex: 5, child: text),
          const SizedBox(width: 24),
          Expanded(flex: 4, child: picture),
        ]);
      }),
    );
  }
}

/// Partition table (`#partitionTable`): name, range, length.
class _PartitionTable extends StatelessWidget {
  const _PartitionTable({required this.matrix});
  final Supermatrix matrix;

  @override
  Widget build(BuildContext context) {
    final l10n = context.l10n;
    final t = context.tokens;
    final parts = matrix.partitions;
    if (parts.isEmpty) return const SizedBox.shrink();
    int size(Partition p) => p.ranges.fold(0, (a, r) => a + ((r.end - r.start) ~/ r.stride) + 1);
    return RecordsTable(
      columns: [l10n.partitionName, l10n.partitionRange, l10n.partitionSites, l10n.partitionType],
      flex: const [3, 4, 1, 1],
      rowCount: parts.length,
      height: 300,
      cellBuilder: (context, r, col) {
        final p = parts[r];
        return switch (col) {
          0 => Text(p.name, style: const TextStyle(fontWeight: FontWeight.w600)),
          1 => Text(p.rangeString, style: const TextStyle(fontFamily: monoFont, fontSize: 12)),
          2 => Text('${size(p)}', style: TextStyle(color: t.textMuted)),
          _ => Text(matrix.isProtein ? 'AA' : 'DNA', style: TextStyle(color: t.textMuted)),
        };
      },
    );
  }
}

/// One partition scheme on its own line: radio, name, what it does, and a
/// picture of the partitions it makes.
class _SchemeOption extends StatelessWidget {
  const _SchemeOption({required this.selected, required this.title, required this.body, required this.onTap, required this.picture});
  final bool selected;
  final String title, body;
  final VoidCallback onTap;
  final Widget picture;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final text = Column(crossAxisAlignment: CrossAxisAlignment.start, mainAxisSize: MainAxisSize.min, children: [
      Text(title, style: TextStyle(fontSize: 14, fontWeight: FontWeight.w700, color: selected ? t.accentText : t.text)),
      const SizedBox(height: 3),
      Text(body, style: TextStyle(fontSize: 12.5, color: t.textMuted, height: 1.4)),
    ]);
    return MouseRegion(
      cursor: SystemMouseCursors.click,
      child: GestureDetector(
        onTap: onTap,
        child: AnimatedContainer(
          duration: const Duration(milliseconds: 150),
          padding: const EdgeInsets.all(14),
          decoration: BoxDecoration(
            color: selected ? t.accentSoft : t.surface,
            borderRadius: BorderRadius.circular(12),
            border: Border.all(color: selected ? SplaceBlue.c500 : t.border, width: selected ? 2 : 1),
          ),
          child: LayoutBuilder(builder: (context, box) {
            final radio = Container(
              width: 16,
              height: 16,
              margin: const EdgeInsets.only(top: 2, right: 12),
              decoration: BoxDecoration(shape: BoxShape.circle, border: Border.all(color: selected ? t.accent : t.textFaint, width: selected ? 5 : 1.5)),
            );
            if (box.maxWidth < 720) {
              return Row(crossAxisAlignment: CrossAxisAlignment.start, children: [
                radio,
                Expanded(child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [text, const SizedBox(height: 10), picture])),
              ]);
            }
            return Row(crossAxisAlignment: CrossAxisAlignment.center, children: [
              radio,
              Expanded(flex: 5, child: text),
              const SizedBox(width: 20),
              Expanded(flex: 4, child: picture),
            ]);
          }),
        ),
      ),
    );
  }
}

/// The partitions a scheme makes from the project's coding loci (or an
/// example): coloured blocks along the matrix and their names.
class _SchemePicture extends StatelessWidget {
  const _SchemePicture({required this.scheme, required this.loci});
  final PartitionScheme scheme;

  /// Markers of the project (name, codon positions known).
  final List<({String name, bool coding})> loci;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final l10n = context.l10n;
    const locusColors = SplacePalette.categorical;
    const posColors = codonPositionColors;
    // Up to three markers of this project, a coding one first when there
    // is one (so the codon split shows); examples before anything is aligned.
    final all = loci.isEmpty ? const [(name: 'ATP6', coding: true), (name: 'COI', coding: true)] : loci;
    final shown = <({String name, bool coding})>[
      ...all.where((l) => l.coding).take(1),
      ...all.where((l) => !l.coding).take(1),
    ];
    for (final l in all) {
      if (shown.length >= 3) break;
      if (!shown.contains(l)) shown.add(l);
    }
    final codingCount = all.where((l) => l.coding).length, otherCount = all.length - codingCount;
    // Each marker: its blocks (label, colour, flex).
    final blocks = <List<(String, Color, int)>>[];
    for (final (i, l) in shown.indexed) {
      final codons = l.coding && (scheme == PartitionScheme.byLocusAndCodon || scheme == PartitionScheme.byCodonPosition);
      blocks.add(switch (scheme) {
        PartitionScheme.byLocusAndCodon when codons => [
            for (var p = 1; p <= 3; p++) ('$p', Color.lerp(locusColors[i % 3], Colors.white, (p - 1) * 0.28)!, 1),
          ],
        PartitionScheme.byCodonPosition when codons => [for (var p = 1; p <= 3; p++) ('$p', posColors[p - 1], 1)],
        PartitionScheme.none => [(l.name, const Color(0xFF9CA3AF), 3)],
        _ => [(l.name, locusColors[i % 3], 3)],
      });
    }
    // Every partition of the matrix, in words (the first ones named).
    final names = <String>[
      for (final l in all)
        ...switch (scheme) {
          PartitionScheme.byLocusAndCodon when l.coding => ['${l.name}_pos1', '${l.name}_pos2', '${l.name}_pos3'],
          PartitionScheme.byCodonPosition when l.coding => const <String>[],
          PartitionScheme.none => const <String>[],
          _ => [l.name],
        },
    ];
    if (scheme == PartitionScheme.byCodonPosition && codingCount > 0) names.insertAll(0, ['codon_pos1', 'codon_pos2', 'codon_pos3']);
    if (scheme == PartitionScheme.none) names.add(l10n.schemeOneModel);
    final more = all.length > shown.length ? all.length - shown.length : 0;
    return Column(crossAxisAlignment: CrossAxisAlignment.stretch, mainAxisSize: MainAxisSize.min, children: [
      Row(children: [
        for (final l in shown)
          Expanded(
            child: Padding(
              padding: const EdgeInsets.only(right: 4),
              child: Text(l.name, overflow: TextOverflow.ellipsis, style: TextStyle(fontSize: 10.5, fontWeight: FontWeight.w700, color: t.textMuted)),
            ),
          ),
        if (more > 0) Text('+$more', style: TextStyle(fontSize: 10.5, fontWeight: FontWeight.w700, color: t.textFaint)),
      ]),
      const SizedBox(height: 3),
      Row(children: [
        for (final b in blocks)
          Expanded(
            child: Padding(
              padding: const EdgeInsets.only(right: 4),
              child: Row(children: [
                for (final (label, color, flex) in b)
                  Expanded(
                    flex: flex,
                    child: Container(
                      height: 20,
                      margin: const EdgeInsets.only(right: 1.5),
                      alignment: Alignment.center,
                      decoration: BoxDecoration(color: color, borderRadius: BorderRadius.circular(3)),
                      child: Text(label, maxLines: 1, overflow: TextOverflow.clip, style: const TextStyle(fontSize: 10, fontWeight: FontWeight.w700, color: Colors.white)),
                    ),
                  ),
              ]),
            ),
          ),
        if (more > 0) Text('…', style: TextStyle(fontSize: 12, fontWeight: FontWeight.w700, color: t.textFaint)),
      ]),
      const SizedBox(height: 5),
      Text(
        '${l10n.schemePartitionsCount('${names.length}')}: ${names.take(6).join(', ')}${names.length > 6 ? ', …' : ''}',
        maxLines: 2,
        overflow: TextOverflow.ellipsis,
        style: TextStyle(fontSize: 11, color: t.textMuted, fontFamily: monoFont),
      ),
      if (loci.isNotEmpty && codingCount > 0 && otherCount > 0 && scheme != PartitionScheme.none && scheme != PartitionScheme.byLocus)
        Text(l10n.schemeMixed('$codingCount', '$otherCount'), style: TextStyle(fontSize: 11, color: t.textFaint)),
    ]);
  }
}

/// Per-locus codon exclusions in a few words: the loci when there are up
/// to three, otherwise how many loci lose each set of positions (hundreds
/// of loci in one line ran out of the window).
String perLocusSummary(BuildContext context, Map<String, Set<int>> exclude) {
  if (exclude.length <= 3) return exclude.entries.map((e) => '${e.key}: ${e.value.join('+')}').join('; ');
  final groups = <String, int>{};
  for (final v in exclude.values) {
    final k = (v.toList()..sort()).join('+');
    groups[k] = (groups[k] ?? 0) + 1;
  }
  final sorted = groups.entries.toList()..sort((a, b) => b.value.compareTo(a.value));
  return sorted.map((e) => context.l10n.perLocusCount(e.key, '${e.value}')).join('; ');
}

/// Step 7c: the outgroup the alignment distances point to (computed in the
/// background once per set of alignments), with a button to use it.
class _OutgroupSuggestionPanel extends StatelessWidget {
  const _OutgroupSuggestionPanel({required this.onUse});
  final ValueChanged<List<String>> onUse;

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    if (c.loci.values.any((l) => l.included && l.aligned != null) && !c.outgroupSuggesting && !c.busy) {
      WidgetsBinding.instance.addPostFrameCallback((_) => c.suggestOutgroups());
    }
    final s = c.outgroupSuggestion;
    if (s == null) {
      return TonePanel(
        tone: PanelTone.info,
        child: Row(children: [
          if (c.outgroupSuggesting) ...[
            SizedBox(width: 14, height: 14, child: CircularProgressIndicator(strokeWidth: 2, color: t.accent)),
            const SizedBox(width: 10),
          ],
          Expanded(child: Text(c.outgroupSuggesting ? l10n.outgroupSuggesting : l10n.outgroupSuggestNone, style: const TextStyle(fontSize: 13))),
        ]),
      );
    }
    final applied = c.concat.outgroups.length == s.taxa.length && c.concat.outgroups.toSet().containsAll(s.taxa);
    String d(double v) => v.toStringAsFixed(3);
    final others = [for (final x in s.ranked) if (!s.taxa.contains(x)) x].take(5).toList();
    return TonePanel(
      tone: s.confident ? PanelTone.success : PanelTone.warn,
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        Row(children: [
          const Icon(FA.wandMagicSparkles, size: 13),
          const SizedBox(width: 8),
          Expanded(
            child: Text(s.confident ? l10n.outgroupSuggestTitle : l10n.outgroupSuggestWeak,
                style: const TextStyle(fontSize: 13.5, fontWeight: FontWeight.w700)),
          ),
          InfoDot(l10n.outgroupSuggestInfo),
        ]),
        const SizedBox(height: 8),
        Wrap(spacing: 6, runSpacing: 6, crossAxisAlignment: WrapCrossAlignment.center, children: [
          for (final x in s.taxa)
            TagChip('$x · ${d(s.meanDistance[x]!)}', color: t.accentText),
          const SizedBox(width: 4),
          ChipButton(
            label: applied ? l10n.outgroupSuggestApplied : l10n.outgroupSuggestUse,
            icon: applied ? FA.circleCheck : FA.check,
            keepIcon: true,
            style: ChipButtonStyle.accent,
            onPressed: applied ? null : () => onUse(s.taxa),
          ),
        ]),
        const SizedBox(height: 8),
        Text(
          l10n.outgroupSuggestStats(d(s.between), d(s.within), '${s.loci}', '${s.taxonCount}'),
          style: TextStyle(fontSize: 12.5, color: t.textMuted),
        ),
        if (others.isNotEmpty) ...[
          const SizedBox(height: 4),
          Text(
            '${l10n.outgroupSuggestNext}: ${others.map((x) => '$x (${d(s.meanDistance[x]!)})').join(', ')}',
            style: TextStyle(fontSize: 12, color: t.textMuted, fontStyle: FontStyle.italic),
          ),
        ],
      ]),
    );
  }
}
