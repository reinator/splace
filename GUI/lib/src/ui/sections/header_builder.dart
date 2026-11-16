import 'package:flutter/material.dart';

import '../../core/genbank_markers.dart';
import '../common.dart';
import '../fa.dart';
import '../theme.dart';
import '../widgets/dialogs.dart';
import '../widgets/kit.dart';

/// "FASTA Header Builder" modal. Left: presets, the header's fields (drag
/// to reorder, with an example value), the fields that can be added and
/// the separator. Right: the first record's header built piece by piece,
/// the headers of the first records and whether they are unique.
Future<void> showHeaderBuilder(BuildContext context) => showDialog(
      context: context,
      builder: (_) => ProjectScope(controller: ProjectScope.read(context), child: const _HeaderBuilder()),
    );

class _HeaderBuilder extends StatefulWidget {
  const _HeaderBuilder();

  @override
  State<_HeaderBuilder> createState() => _HeaderBuilderState();
}

/// Colours of the field pieces in the preview, by position.
const _pieceColors = [
  SplaceBlue.c600, SplacePalette.red900, SplacePalette.blue600, SplacePalette.red700, SplacePalette.blue900, SplaceBlue.c400,
  SplaceRed.c700, SplacePalette.blue600, SplaceBlue.c800,
];

const _presets = [
  ['species', 'accession'],
  ['accession', 'species'],
  ['family', 'species', 'accession'],
  ['order', 'family', 'species', 'accession'],
  ['accession', 'family', 'genus', 'species'],
];

class _HeaderBuilderState extends State<_HeaderBuilder> {
  late HeaderTemplate _h = ProjectScope.read(context).header;

  String _label(BuildContext context, String key) {
    final l10n = context.l10n;
    return switch (key) {
      'species' => l10n.colSpecies,
      'accession' => l10n.colAccession,
      'kingdom' => l10n.colKingdom,
      'phylum' => l10n.colPhylum,
      'class' => l10n.colClass,
      'order' => l10n.colOrder,
      'family' => l10n.colFamily,
      'genus' => l10n.colGenus,
      _ => l10n.colAuthorship,
    };
  }

  void _set(List<String> fields) => setState(() => _h = _h.copyWith(fields: fields));

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    final first = c.records.isEmpty ? null : c.records.first;
    String example(String k) => first == null ? '' : headerFieldValue(first, k);
    final available = [for (final k in headerFieldKeys) if (!_h.fields.contains(k)) k];
    // Headers repeated inside a marker: joined on purpose with one sample
    // per species, an error otherwise (the sequences could not be told apart).
    final dups = c.headerConflicts(_h);
    final blocked = dups.isNotEmpty && !c.onePerTaxonMode;
    final preview = c.records.take(12).map(_h.build).toList();

    Widget block(String label, Widget child, {String? hint}) => Padding(
          padding: const EdgeInsets.only(bottom: 16),
          child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
            Kicker(label, size: 12),
            if (hint != null) ...[
              const SizedBox(height: 2),
              Text(hint, style: TextStyle(fontSize: 12, color: t.textMuted)),
            ],
            const SizedBox(height: 8),
            child,
          ]),
        );

    final left = Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
      block(
        l10n.headerPresets,
        Wrap(spacing: 8, runSpacing: 8, children: [
          for (final p in _presets)
            ChipButton(
              label: p.map((k) => _label(context, k)).join(_h.separator),
              style: _h.fields.join() == p.join() ? ChipButtonStyle.accent : ChipButtonStyle.neutral,
              onPressed: () => _set(List.of(p)),
            ),
        ]),
      ),
      block(
        l10n.headerFields,
        hint: l10n.headerDragHint,
        _h.fields.isEmpty
            ? Container(
                padding: const EdgeInsets.all(14),
                decoration: BoxDecoration(color: t.surface2, borderRadius: BorderRadius.circular(10), border: Border.all(color: t.border)),
                child: Text(l10n.headerEmpty, style: TextStyle(color: t.textFaint, fontStyle: FontStyle.italic)),
              )
            : ReorderableListView(
                shrinkWrap: true,
                physics: const NeverScrollableScrollPhysics(),
                buildDefaultDragHandles: false,
                // The new index already accounts for the removed item.
                onReorderItem: (from, to) {
                  final f = List.of(_h.fields);
                  final x = f.removeAt(from);
                  f.insert(to, x);
                  _set(f);
                },
                children: [
                  for (final (i, k) in _h.fields.indexed)
                    Container(
                      key: ValueKey(k),
                      margin: const EdgeInsets.only(bottom: 6),
                      padding: const EdgeInsets.symmetric(horizontal: 8, vertical: 6),
                      decoration: BoxDecoration(
                        color: t.surface,
                        borderRadius: BorderRadius.circular(10),
                        border: Border.all(color: t.border),
                      ),
                      child: Row(children: [
                        ReorderableDragStartListener(
                          index: i,
                          child: MouseRegion(
                            cursor: SystemMouseCursors.grab,
                            child: Padding(padding: const EdgeInsets.all(6), child: Icon(FA.gripVertical, size: 13, color: t.textFaint)),
                          ),
                        ),
                        Container(
                          width: 24,
                          height: 24,
                          alignment: Alignment.center,
                          decoration: BoxDecoration(color: _pieceColors[i % _pieceColors.length], shape: BoxShape.circle),
                          child: Text('${i + 1}', style: const TextStyle(fontSize: 12, fontWeight: FontWeight.w700, color: Colors.white)),
                        ),
                        const SizedBox(width: 10),
                        SizedBox(width: 120, child: Text(_label(context, k), style: TextStyle(fontWeight: FontWeight.w700, color: t.text))),
                        Expanded(
                          child: Text(example(k), maxLines: 1, overflow: TextOverflow.ellipsis, style: TextStyle(fontSize: 12.5, color: t.textMuted)),
                        ),
                        IconBtn(icon: FA.xmark, tooltip: l10n.remove, size: ChipButton.height, onPressed: () => _set(List.of(_h.fields)..removeAt(i))),
                      ]),
                    ),
                ],
              ),
      ),
      if (available.isNotEmpty)
        block(
          l10n.headerAvailable,
          hint: l10n.headerAddHint,
          Wrap(spacing: 8, runSpacing: 8, children: [
            for (final k in available)
              ChipButton(
                label: example(k).isEmpty ? _label(context, k) : '${_label(context, k)} · ${example(k)}',
                style: ChipButtonStyle.soft,
                onPressed: () => _set([..._h.fields, k]),
              ),
          ]),
        ),
      block(
        l10n.headerSeparator,
        Align(
          alignment: Alignment.centerLeft,
          child: Seg<String>(
            value: _h.separator,
            options: const [('_', '_  underscore'), ('-', '-  hyphen'), ('|', '|  pipe'), ('.', '.  dot')],
            onChanged: (v) => setState(() => _h = _h.copyWith(separator: v)),
          ),
        ),
      ),
    ]);

    // The first record's header, one coloured piece per field.
    final pieces = <Widget>[
      Text('>', style: TextStyle(fontSize: 15, fontWeight: FontWeight.w700, color: t.textMuted)),
      for (final (i, k) in _h.fields.indexed) ...[
        if (i > 0) Text(_h.separator, style: TextStyle(fontSize: 15, fontWeight: FontWeight.w700, color: t.textMuted)),
        Tooltip(
          message: _label(context, k),
          child: Container(
            padding: const EdgeInsets.symmetric(horizontal: 6, vertical: 3),
            decoration: BoxDecoration(
              color: _pieceColors[i % _pieceColors.length].withValues(alpha: 0.12),
              borderRadius: BorderRadius.circular(6),
              border: Border.all(color: _pieceColors[i % _pieceColors.length].withValues(alpha: 0.5)),
            ),
            child: Text(example(k).replaceAll(' ', '_'),
                style: TextStyle(fontSize: 14, fontWeight: FontWeight.w600, color: _pieceColors[i % _pieceColors.length])),
          ),
        ),
      ],
    ];

    final right = Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
      block(
        l10n.headerFirstRecord,
        Container(
          padding: const EdgeInsets.all(12),
          decoration: BoxDecoration(color: t.surface2, borderRadius: BorderRadius.circular(10), border: Border.all(color: t.border)),
          child: Wrap(spacing: 2, runSpacing: 6, crossAxisAlignment: WrapCrossAlignment.center, children: pieces),
        ),
      ),
      // Right under the example, so a blocking problem is seen first.
      if (dups.isEmpty)
        TonePanel(tone: PanelTone.success, child: Text(l10n.headerUnique, style: const TextStyle(fontSize: 13)))
      else
        TonePanel(
          tone: blocked ? PanelTone.danger : PanelTone.info,
          icon: blocked ? FA.circleExclamation : FA.circleInfo,
          child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
            Text(blocked ? l10n.headerConflictBlock('${dups.length}') : l10n.onePerTaxonJoined('${dups.length}'),
                style: const TextStyle(fontSize: 13, fontWeight: FontWeight.w600)),
            if (blocked && !_h.fields.contains('accession')) ...[
              const SizedBox(height: 8),
              ChipButton(label: l10n.headerAddAccession, icon: FA.plus, style: ChipButtonStyle.accent, onPressed: () => _set([..._h.fields, 'accession'])),
            ],
            const SizedBox(height: 6),
            for (final e in dups.entries.take(4))
              Text('${e.key}: ${e.value.join(', ')}', maxLines: 1, overflow: TextOverflow.ellipsis, style: const TextStyle(fontSize: 12.5)),
            if (dups.length > 4) Text(l10n.andMore('${dups.length - 4}'), style: const TextStyle(fontSize: 12.5)),
          ]),
        ),
      const SizedBox(height: 16),
      block(l10n.headerFirstRecords, PreviewBox(label: l10n.previewLabel, text: preview.map((e) => '>$e').join('\n'), maxHeight: 260)),
    ]);

    return SplaceDialog(
      icon: FA.tags,
      title: l10n.headerBuilderTitle,
      subtitle: l10n.headerBuilderSubtitle,
      actions: [
        DialogCancelButton(onPressed: () => Navigator.pop(context)),
        FilledButton(
          onPressed: _h.fields.isEmpty || blocked
              ? null
              : () {
                  c.setHeader(_h);
                  c.confirmHeader();
                  Navigator.pop(context);
                },
          child: Text(l10n.headerConfirm),
        ),
      ],
      child: Row(crossAxisAlignment: CrossAxisAlignment.start, children: [
        Expanded(flex: 11, child: left),
        const SizedBox(width: 24),
        Expanded(flex: 10, child: right),
      ]),
    );
  }
}
