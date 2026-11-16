import 'package:flutter/material.dart';

import '../../core/fasta.dart';
import '../../core/locus_organizer.dart';
import '../common.dart';
import '../fa.dart';
import '../theme.dart';
import '../widgets/kit.dart';
import 'feature_types_panel.dart';
import 'organizer_dialog.dart';

/// Step 3 — Feature Types (GenBank: CDS / rRNA / tRNA) and, for FASTA
/// files, where the sample and the locus names come from.
class FieldsSection extends StatelessWidget {
  const FieldsSection({super.key, required this.active});
  final bool active;

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final gb = c.records.isNotEmpty && c.taxonomyResolved;
    final fa = c.inputs.isNotEmpty;
    return StepSection(
      badge: '3',
      icon: FA.layerGroup,
      title: gb ? l10n.step3 : l10n.step3Fasta,
      description: gb ? l10n.featureTypesIntro : l10n.step3Intro,
      active: active,
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        if (gb)
          SubAnchor(
            id: '3a',
            child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
              SubHeading(id: '3a', title: l10n.subFeatureTypes),
              const FeatureTypesPanel(),
            ]),
          ),
        if (gb && fa) const SizedBox(height: 24),
        if (fa)
          SubAnchor(
            id: fastaFieldsSub(c),
            child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
              SubHeading(id: fastaFieldsSub(c), title: l10n.step3Fasta, body: gb ? l10n.step3Intro : null),
              const _FastaFieldsPanel(),
            ]),
          ),
      ]),
    );
  }
}

class _FastaFieldsPanel extends StatelessWidget {
  const _FastaFieldsPanel();

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    final cfg = c.organizer;
    final r = c.organized;
    // Same size as "Open folder" beside it (and as step 2's button).
    final button = ChipButton(
      label: l10n.orgOpen,
      icon: FA.tableColumns,
      style: ChipButtonStyle.accent,
      onPressed: c.busy ? null : () => showOrganizerDialog(context),
    );
    if (r == null) {
      return TonePanel(
        tone: c.organizeProblem == null ? PanelTone.info : PanelTone.danger,
        icon: c.organizeProblem == null ? FA.circleInfo : FA.triangleExclamation,
        title: c.organizeProblem == null ? l10n.orgPending : l10n.orgProblem(c.organizeProblem!),
        child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
          Text(l10n.orgPendingBody),
          const SizedBox(height: 12),
          button,
        ]),
      );
    }
    Widget fact(String label, String value, {Color? color}) => Padding(
          padding: const EdgeInsets.only(bottom: 6),
          child: Row(crossAxisAlignment: CrossAxisAlignment.start, children: [
            SizedBox(width: 170, child: Text(label, style: TextStyle(fontSize: 13, color: t.textMuted))),
            Expanded(child: Text(value, style: TextStyle(fontSize: 13, fontWeight: FontWeight.w600, color: color ?? t.text))),
          ]),
        );
    return OptionCard(
      kicker: l10n.optionsKicker,
      title: l10n.orgTitle,
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        fact(l10n.sampleField, cfg.sample.describe(), color: sampleColor),
        fact(l10n.locusField, cfg.locus.describe(), color: markerColor),
        fact(l10n.orgMatrix, l10n.orgSummary('${r.samples.length}', '${r.loci.length}', '${(r.occupancy * 100).round()}')),
        if (r.unassigned.isNotEmpty) fact(l10n.orgUnassignedLabel, '${r.unassigned.length}', color: SplaceTokens.danger),
        if (c.lociExportDir != null) fact(l10n.orgWritten, c.lociExportDir!),
        const SizedBox(height: 8),
        Wrap(spacing: 8, runSpacing: 8, children: [
          button,
          if (c.lociExportDir != null)
            ChipButton(label: l10n.openFolder, icon: FA.folderOpen, onPressed: () => c.openFolder(c.lociExportDir!)),
        ]),
      ]),
    );
  }
}

/// Edits one [FieldRule] and previews its value on an example record.
class FieldRuleEditor extends StatefulWidget {
  const FieldRuleEditor({super.key, required this.rule, required this.example, required this.onChanged});

  final FieldRule rule;
  final (String, SeqRecord) example;
  final ValueChanged<FieldRule> onChanged;

  @override
  State<FieldRuleEditor> createState() => _FieldRuleEditorState();
}

class _FieldRuleEditorState extends State<FieldRuleEditor> {
  late final _delimiters = TextEditingController(text: widget.rule.delimiters);
  late final _start = TextEditingController(text: '${widget.rule.start}');
  late final _end = TextEditingController(text: widget.rule.end?.toString() ?? '');
  late final _tag = TextEditingController(text: widget.rule.tag);
  late final _pattern = TextEditingController(text: widget.rule.pattern);

  @override
  void didUpdateWidget(FieldRuleEditor old) {
    super.didUpdateWidget(old);
    // Sync when the rule changes from outside (e.g. a new layout guess)
    // without disturbing the field being typed in.
    void sync(TextEditingController t, String v) {
      if (t.text != v) t.text = v;
    }

    final r = widget.rule;
    sync(_delimiters, r.delimiters);
    if (int.tryParse(_start.text) != r.start) sync(_start, '${r.start}');
    if (int.tryParse(_end.text) != r.end) sync(_end, r.end?.toString() ?? '');
    sync(_tag, r.tag);
    sync(_pattern, r.pattern);
  }

  @override
  void dispose() {
    for (final t in [_delimiters, _start, _end, _tag, _pattern]) {
      t.dispose();
    }
    super.dispose();
  }

  @override
  Widget build(BuildContext context) {
    final l10n = context.l10n;
    final t = context.tokens;
    final r = widget.rule;
    final (fileName, record) = widget.example;
    final preview = r.extract(fileName, record);
    final source = r.source == FieldSource.fileName ? fileStem(fileName) : record.header;
    final tokens = tokenize(source, r.delimiters);

    return Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
      Align(
        alignment: Alignment.centerLeft,
        child: Seg<FieldSource>(
          value: r.source,
          options: [
            (FieldSource.fileName, l10n.sourceFileName),
            (FieldSource.header, l10n.sourceHeader),
            (FieldSource.tag, l10n.sourceTag),
            (FieldSource.regex, l10n.sourceRegex),
          ],
          onChanged: (s) => widget.onChanged(r.copyWith(source: s)),
        ),
      ),
      const SizedBox(height: 14),
      if (r.source == FieldSource.fileName || r.source == FieldSource.header) ...[
        AutoGrid(minItemWidth: 110, children: [
          Field(
            label: l10n.delimiters,
            child: TextField(
              controller: _delimiters,
              style: const TextStyle(fontFamily: monoFont),
              decoration: const InputDecoration(hintText: '_ | - .'),
              onChanged: (v) => widget.onChanged(r.copyWith(delimiters: v)),
            ),
          ),
          Field(
            label: l10n.tokenStart,
            child: TextField(
              controller: _start,
              keyboardType: TextInputType.number,
              onChanged: (v) {
                final n = int.tryParse(v);
                if (n != null) widget.onChanged(r.copyWith(start: n));
              },
            ),
          ),
          Field(
            label: l10n.tokenEnd,
            child: TextField(
              controller: _end,
              keyboardType: TextInputType.number,
              onChanged: (v) => widget.onChanged(r.copyWith(end: () => int.tryParse(v))),
            ),
          ),
        ]),
        const SizedBox(height: 10),
        // Click a token to select it alone.
        Wrap(spacing: 6, runSpacing: 6, children: [
          for (var i = 0; i < tokens.length && i < 12; i++)
            ChipButton(
              label: '$i: ${tokens[i]}',
              style: (r.start == i && r.end == i + 1) ? ChipButtonStyle.soft : ChipButtonStyle.neutral,
              onPressed: () => widget.onChanged(r.copyWith(start: i, end: () => i + 1)),
            ),
        ]),
      ],
      if (r.source == FieldSource.tag)
        Field(
          label: l10n.tagKey,
          child: TextField(
            controller: _tag,
            decoration: const InputDecoration(hintText: 'gene, locus, organism, sample…'),
            onChanged: (v) => widget.onChanged(r.copyWith(tag: v)),
          ),
        ),
      if (r.source == FieldSource.regex) ...[
        Field(
          label: l10n.regexPattern,
          child: TextField(
            controller: _pattern,
            decoration: const InputDecoration(hintText: r'^([^_]+_[^_]+)'),
            style: const TextStyle(fontFamily: monoFont),
            onChanged: (v) => widget.onChanged(r.copyWith(pattern: v)),
          ),
        ),
        const SizedBox(height: 10),
        SwitchBox(label: l10n.regexOnFileName, value: r.regexOnFileName, onChanged: (v) => widget.onChanged(r.copyWith(regexOnFileName: v))),
      ],
      const SizedBox(height: 14),
      // `.step1-query-preview`
      Container(
        padding: const EdgeInsets.symmetric(horizontal: 14, vertical: 12),
        decoration: BoxDecoration(color: t.surface2, borderRadius: BorderRadius.circular(14), border: Border.all(color: t.border)),
        child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
          Kicker(l10n.previewLabel),
          const SizedBox(height: 6),
          Text.rich(TextSpan(style: const TextStyle(fontFamily: monoFont, fontSize: 12), children: [
            TextSpan(text: '>${record.header.length > 48 ? '${record.header.substring(0, 48)}…' : record.header}', style: TextStyle(color: t.textMuted)),
            TextSpan(text: '   →   ', style: TextStyle(color: t.textFaint)),
            TextSpan(
              text: preview ?? '∅',
              style: TextStyle(color: preview == null ? SplaceTokens.danger : t.accent),
            ),
          ])),
        ]),
      ),
    ]);
  }
}
