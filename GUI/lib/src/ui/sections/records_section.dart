import 'package:flutter/material.dart';

import '../../core/fasta.dart';
import '../../core/locus_organizer.dart';
import '../common.dart';
import '../fa.dart';
import '../theme.dart';
import '../widgets/kit.dart';
import 'genbank_records_panel.dart';
import 'organizer_dialog.dart';

/// Step 2 — Loaded Records: the GenBank records table (taxonomy, marker
/// counts) and the FASTA header explorer (sample and locus read from each
/// header).
class RecordsSection extends StatelessWidget {
  const RecordsSection({super.key, required this.active});
  final bool active;

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final gb = c.records.isNotEmpty;
    final fa = c.inputs.isNotEmpty;
    return StepSection(
      badge: '2',
      icon: FA.listCheck,
      title: l10n.step2,
      count: c.records.length + c.inputs.fold<int>(0, (a, f) => a + f.records.length),
      active: active,
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        if (gb) ...[
          if (fa) SectionHeading(icon: FA.dna, title: l10n.gbRecordsHeading),
          const GenBankRecordsPanel(),
        ],
        if (gb && fa) const SizedBox(height: 24),
        if (fa)
          SubAnchor(
            id: fastaHeadersSub(c),
            child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
              SubHeading(id: fastaHeadersSub(c), title: l10n.fastaHeadersHeading),
              const _FastaHeadersPanel(),
            ]),
          ),
      ]),
    );
  }
}

class _FastaHeadersPanel extends StatefulWidget {
  const _FastaHeadersPanel();

  @override
  State<_FastaHeadersPanel> createState() => _FastaHeadersPanelState();
}

class _FastaHeadersPanelState extends State<_FastaHeadersPanel> {
  String _filter = '';
  List<(String, SeqRecord)> _rows = const [];
  Object? _source;
  int _sourceCount = -1;

  List<(String, SeqRecord)> _allRows(List<FastaInput> inputs) {
    if (!identical(inputs, _source) || inputs.length != _sourceCount) {
      _source = inputs;
      _sourceCount = inputs.length;
      _rows = [for (final f in inputs) for (final r in f.records) (f.fileName, r)];
    }
    return _rows;
  }

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    final cfg = c.organizer;
    final all = _allRows(c.inputs);
    final q = _filter.toLowerCase();
    final rows = q.isEmpty
        ? all
        : all.where((r) => r.$1.toLowerCase().contains(q) || r.$2.header.toLowerCase().contains(q)).toList();
    final columns = headerTokenColumns(
      all.take(2000).map((r) => r.$2.header).toList(),
      cfg.locus.source == FieldSource.header ? cfg.locus.delimiters : '_ ',
    );
    final s = c.suggestion;
    final layout = switch (s?.kind) {
      LayoutKind.samplePerFile => l10n.layoutSamplePerFile,
      LayoutKind.locusPerFile => l10n.layoutLocusPerFile,
      LayoutKind.mixedHeader => l10n.layoutMixedHeader,
      LayoutKind.taggedHeader => l10n.layoutTaggedHeader,
      _ => l10n.layoutUnknown,
    };

    return Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        Row(children: [
          const Spacer(),
          if (c.organizing) ...[
            SizedBox(width: 12, height: 12, child: CircularProgressIndicator(strokeWidth: 2, color: t.accent)),
            const SizedBox(width: 8),
            Text(l10n.organizing, style: TextStyle(fontSize: 13, color: t.textMuted)),
            const SizedBox(width: 10),
          ],
          ChipButton(label: l10n.orgOpen, icon: FA.tableColumns, style: ChipButtonStyle.accent, onPressed: c.busy ? null : () => showOrganizerDialog(context)),
          const SizedBox(width: 8),
          ChipButton(label: l10n.clearFasta, icon: FA.trashCan, style: ChipButtonStyle.danger, onPressed: c.busy ? null : c.clearInputs),
        ]),
        const SizedBox(height: 10),
        // `.step2-intro-card`
        Container(
          padding: const EdgeInsets.symmetric(horizontal: 18, vertical: 16),
          decoration: BoxDecoration(
            borderRadius: BorderRadius.circular(16),
            border: Border.all(color: t.accentSoftBorder),
            gradient: Theme.of(context).brightness == Brightness.dark
                ? null
                : const LinearGradient(begin: Alignment.topLeft, end: Alignment.bottomRight, colors: [Color(0xEBECF1FF), Color(0xFAF7F9FF)]),
            color: Theme.of(context).brightness == Brightness.dark ? t.accentSoft : null,
          ),
          child: DefaultTextStyle.merge(
            style: TextStyle(color: t.accentText, fontSize: 14, height: 1.55),
            child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
              Text(l10n.step2Intro),
              if (s != null) ...[
                const SizedBox(height: 10),
                Row(children: [
                  Icon(FA.wandMagicSparkles, size: 14, color: t.accent),
                  const SizedBox(width: 8),
                  Flexible(child: Text('${l10n.detectedLayout}: $layout', style: const TextStyle(fontWeight: FontWeight.w700))),
                  const SizedBox(width: 8),
                  TagChip(l10n.confidence('${(s.confidence * 100).round()}')),
                ]),
                for (final e in s.evidence) Bullet(Text(e, style: const TextStyle(fontSize: 13))),
              ],
            ]),
          ),
        ),
        const SizedBox(height: 16),
        Kicker(l10n.tokenColumns),
        const SizedBox(height: 8),
        SizedBox(
          height: 84,
          child: ListView(scrollDirection: Axis.horizontal, children: [
            for (final col in columns)
              Container(
                width: 210,
                margin: const EdgeInsets.only(right: 8),
                padding: const EdgeInsets.symmetric(horizontal: 10, vertical: 8),
                decoration: BoxDecoration(
                  color: col.geneLikeFraction > 0.6 ? t.accentSoft : t.surface2,
                  borderRadius: BorderRadius.circular(10),
                  border: Border.all(color: col.geneLikeFraction > 0.6 ? t.accentSoftBorder : t.border),
                ),
                child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
                  Text(col.index < 0 ? '#last' : '#${col.index}', style: TextStyle(fontWeight: FontWeight.w700, fontSize: 13, color: t.text)),
                  Text(l10n.tokenColumnInfo('${col.distinct}', '${(col.geneLikeFraction * 100).round()}'),
                      style: TextStyle(fontSize: 11, color: t.textMuted)),
                  Text(col.examples.join(', '),
                      maxLines: 1, overflow: TextOverflow.ellipsis, style: TextStyle(fontSize: 11, fontFamily: monoFont, color: t.textFaint)),
                ]),
              ),
          ]),
        ),
        const SizedBox(height: 14),
        TextField(
          decoration: InputDecoration(
            prefixIcon: Icon(FA.magnifyingGlass, size: 14, color: t.textFaint),
            hintText: l10n.searchHeaders,
          ),
          onChanged: (v) => setState(() => _filter = v),
        ),
        const SizedBox(height: 12),
        RecordsTable(
          columns: [l10n.colFile, l10n.colHeader, l10n.colSample, l10n.colLocus, l10n.colLength],
          flex: const [3, 4, 3, 2, 1],
          rowCount: rows.length,
          cellBuilder: (context, r, col) {
            final (file, record) = rows[r];
            switch (col) {
              case 0:
                return HighlightText(file, _filter, style: TextStyle(fontSize: 12, color: t.textMuted));
              case 1:
                return HighlightText(record.header, _filter, style: TextStyle(fontFamily: monoFont, fontSize: 12, color: t.text));
              case 2:
                final v = cfg.sample.extract(file, record);
                return Text(v ?? '∅',
                    style: TextStyle(fontStyle: FontStyle.italic, color: v == null ? SplaceTokens.danger : t.text));
              case 3:
                final v = cfg.locus.extract(file, record);
                return Text(v ?? '∅',
                    style: TextStyle(fontWeight: FontWeight.w600, color: v == null ? SplaceTokens.danger : t.accent));
              default:
                return Text('${record.length} bp', style: TextStyle(color: t.textMuted));
            }
          },
        ),
      ]);
  }
}
