import 'package:flutter/material.dart';

import '../../core/marker_search.dart';
import '../../core/syngenes.dart';
import '../../services/ncbi.dart';
import '../common.dart';
import '../fa.dart';
import '../theme.dart';
import '../widgets/dialogs.dart';
import '../widgets/kit.dart';
import '../widgets/splace_tooltip.dart';

/// Colours of the HTML review modal (`.taxonomy-stacked-bar-segment`):
/// SPLACE blue for complete, yellow for partial, SPLACE red for UNVERIFIED.
const _completeColor = SplaceBlue.c600;
const _partialColor = SplacePalette.yellow;
const _partialText = Color(0xFF1F2937);
const _unverifiedColor = SplaceRed.c600;
const _otherColor = Color(0xFF9CA3AF);

/// "Review genomes before importing" (Option 1), laid out like the HTML
/// modal: almost the whole window, a stacked bar of complete / partial /
/// UNVERIFIED genomes (all and selected), the selection shortcuts, a search
/// that highlights its matches and one switch per genome. Returns the
/// genomes to import, or null when cancelled.
Future<List<GenomeCandidate>?> showGenomeReviewDialog(BuildContext context, GenomeSearch search, List<GenomeCandidate> candidates,
        {MarkerSearch? markers}) =>
    showDialog<List<GenomeCandidate>>(
      context: context,
      builder: (_) => _GenomeReviewDialog(search: search, candidates: candidates, markers: markers),
    );

/// Marker searches: for each marker, the longest record of each species
/// (verified records first), so a species is not downloaded many times.
void selectOnePerSpecies(List<GenomeCandidate> all) {
  final keep = <GenomeCandidate>{};
  final markers = {for (final c in all) ...c.markers};
  for (final m in markers) {
    final best = <String, GenomeCandidate>{};
    for (final c in all.where((c) => c.markers.contains(m))) {
      final key = c.species.isEmpty ? c.accession : c.species.toLowerCase();
      final b = best[key];
      bool better(GenomeCandidate x, GenomeCandidate y) =>
          x.unverified != y.unverified ? !x.unverified : (x.length ?? 0) > (y.length ?? 0);
      if (b == null || better(c, b)) best[key] = c;
    }
    keep.addAll(best.values);
  }
  for (final c in all) {
    c.selected = keep.contains(c);
  }
}

class _GenomeReviewDialog extends StatefulWidget {
  const _GenomeReviewDialog({required this.search, required this.candidates, this.markers});
  final GenomeSearch search;
  final List<GenomeCandidate> candidates;

  /// Set for a search of individual markers.
  final MarkerSearch? markers;

  @override
  State<_GenomeReviewDialog> createState() => _GenomeReviewDialogState();
}

enum _Kind { complete, partial, unverified, other }

_Kind _kindOf(GenomeCandidate c) => c.unverified
    ? _Kind.unverified
    : switch (c.classification) {
        GenomeClass.complete => _Kind.complete,
        GenomeClass.partial => _Kind.partial,
        _ => _Kind.other,
      };

class _GenomeReviewDialogState extends State<_GenomeReviewDialog> {
  final _filter = TextEditingController();

  @override
  void initState() {
    super.initState();
    if (widget.markers != null) selectOnePerSpecies(widget.candidates);
  }

  @override
  void dispose() {
    _filter.dispose();
    super.dispose();
  }

  void _select(bool Function(GenomeCandidate) on) => setState(() {
        for (final c in widget.candidates) {
          c.selected = on(c);
        }
      });

  @override
  Widget build(BuildContext context) {
    final l10n = context.l10n;
    final t = context.tokens;
    final search = widget.search;
    final all = widget.candidates;
    final selected = all.where((c) => c.selected).toList();
    final visible = all.where((c) => c.matches(_filter.text)).toList();
    final query = _filter.text;

    final markers = widget.markers;
    final scope = markers != null
        ? [
            markers.taxon,
            markers.organelle == OrganelleType.cp ? l10n.organelleChloro : l10n.organelleMito,
            markers.genes.join(', '),
            if (markers.refseqOnly) 'RefSeq',
          ].join(' · ')
        : [
            search.taxon,
            search.organelle == Organelle.chloroplast ? l10n.organelleChloro : l10n.organelleMito,
            search.completeOnly ? l10n.scopeComplete : '${l10n.scopeComplete} + ${l10n.scopePartial}',
            if (search.refseqOnly) 'RefSeq',
          ].join(' · ');

    (String, Color, Color) badge(GenomeCandidate c) => switch (_kindOf(c)) {
          _Kind.unverified => (l10n.reviewUnverified, _unverifiedColor, _unverifiedColor),
          _Kind.complete => (l10n.reviewComplete, _completeColor, t.accent),
          _Kind.partial => (l10n.reviewPartial, _partialColor, _partialText),
          _Kind.other => (c.classification == GenomeClass.mixed ? l10n.reviewMixed : l10n.reviewUnknown, _otherColor, t.textMuted),
        };

    final body = TextStyle(fontSize: 13, color: t.text);

    return SplaceDialog(
      icon: FA.listCheck,
      title: l10n.reviewTitle,
      subtitle: l10n.reviewFound('${all.length}', scope),
      size: DialogSize.full,
      scroll: false,
      footerLeading: Text(l10n.reviewSelectedOf('${selected.length}', '${all.length}'), style: TextStyle(fontSize: 13, color: t.textMuted)),
      actions: [
        DialogCancelButton(onPressed: () => Navigator.pop(context)),
        FilledButton(
          onPressed: selected.isEmpty ? null : () => Navigator.pop(context, selected),
          child: Text(l10n.reviewImportN('${selected.length}')),
        ),
      ],
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        // Summary: what the search returned and what will be imported.
        Container(
          padding: const EdgeInsets.fromLTRB(16, 14, 16, 12),
          decoration: BoxDecoration(color: t.surface2, borderRadius: BorderRadius.circular(12), border: Border.all(color: t.border)),
          child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
            _StackedBar(label: '${l10n.reviewTotal} (${all.length})', genomes: all),
            const SizedBox(height: 10),
            _StackedBar(label: '${l10n.reviewSelected} (${selected.length})', genomes: selected),
            if (markers != null) ...[
              const SizedBox(height: 10),
              // Records and species per marker, found and selected.
              Wrap(spacing: 8, runSpacing: 6, children: [
                for (final g in markers.genes)
                  TagChip(
                    l10n.reviewMarkerCount(
                      g,
                      '${selected.where((c) => c.markers.contains(g)).length}',
                      '${all.where((c) => c.markers.contains(g)).length}',
                      '${{for (final c in selected.where((c) => c.markers.contains(g))) c.species}.length}',
                    ),
                    color: t.accent,
                  ),
              ]),
            ],
            const SizedBox(height: 10),
            Divider(height: 1, color: t.border),
            const SizedBox(height: 10),
            Wrap(alignment: WrapAlignment.center, spacing: 24, runSpacing: 6, children: [
              for (final (kind, label, color) in [
                (_Kind.complete, l10n.reviewComplete, _completeColor),
                (_Kind.partial, l10n.reviewPartial, _partialColor),
                (_Kind.unverified, l10n.reviewUnverified, _unverifiedColor),
                (_Kind.other, l10n.reviewUnknown, _otherColor),
              ])
                if (all.any((c) => _kindOf(c) == kind))
                  Row(mainAxisSize: MainAxisSize.min, children: [
                    Container(width: 12, height: 12, decoration: BoxDecoration(color: color, borderRadius: BorderRadius.circular(2))),
                    const SizedBox(width: 8),
                    Text('$label  ', style: TextStyle(fontSize: 13, color: t.textMuted)),
                    Text('${all.where((c) => _kindOf(c) == kind).length}', style: TextStyle(fontSize: 13, fontWeight: FontWeight.w700, color: t.text)),
                  ]),
            ]),
          ]),
        ),
        const SizedBox(height: 14),
        // Search and selection shortcuts (`.taxonomy-import-toolbar`).
        Row(children: [
          Expanded(
            child: TextField(
              controller: _filter,
              decoration: InputDecoration(
                hintText: l10n.reviewSearch,
                prefixIcon: Icon(FA.magnifyingGlass, size: 13, color: t.textFaint),
                prefixIconConstraints: const BoxConstraints(minWidth: 38),
                suffixIcon: _filter.text.isEmpty
                    ? null
                    : IconButton(onPressed: () => setState(_filter.clear), icon: Icon(FA.xmark, size: 12, color: t.textFaint)),
              ),
              onChanged: (_) => setState(() {}),
            ),
          ),
          const SizedBox(width: 12),
          Wrap(spacing: 8, runSpacing: 8, children: [
            if (markers != null)
              SplaceTooltip(
                message: l10n.reviewOnePerSpeciesInfo,
                child: ChipButton(
                  label: l10n.reviewOnePerSpecies,
                  style: ChipButtonStyle.soft,
                  onPressed: () => setState(() => selectOnePerSpecies(all)),
                ),
              ),
            ChipButton(label: l10n.reviewSelectAll, onPressed: () => _select((_) => true)),
            ChipButton(label: l10n.reviewSelectNone, onPressed: () => _select((_) => false)),
            ChipButton(label: l10n.reviewDisableUnverified, onPressed: () => _select((c) => c.selected && !c.unverified)),
            ChipButton(label: l10n.reviewOnlyComplete, onPressed: () => _select((c) => c.classification == GenomeClass.complete && !c.unverified)),
            ChipButton(label: l10n.reviewOnlyPartial, onPressed: () => _select((c) => c.classification == GenomeClass.partial && !c.unverified)),
          ]),
        ]),
        if (query.trim().isNotEmpty) ...[
          const SizedBox(height: 6),
          Text(l10n.reviewShowing('${visible.length}', '${all.length}'), style: TextStyle(fontSize: 12, color: t.textMuted)),
        ],
        const SizedBox(height: 12),
        Expanded(
          child: visible.isEmpty
              ? Center(child: Text(l10n.reviewNoMatch, style: TextStyle(color: t.textMuted)))
              : RecordsTable(
                  fill: true,
                  columns: [
                    l10n.reviewColImport,
                    l10n.reviewColAccession,
                    l10n.reviewColDescription,
                    if (markers != null) l10n.reviewColMarkers,
                    l10n.reviewColScope,
                  ],
                  flex: markers != null ? const [2, 3, 10, 4, 3] : const [2, 3, 12, 3],
                  rowCount: visible.length,
                  rowHeight: 52,
                  rowColor: (r) => visible[r].unverified
                      ? _unverifiedColor.withValues(alpha: 0.05)
                      : (visible[r].selected ? null : t.surface2.withValues(alpha: 0.6)),
                  onRowTap: (r) => setState(() => visible[r].selected = !visible[r].selected),
                  cellBuilder: (context, r, column) {
                    final c = visible[r];
                    // The markers column sits before the scope.
                    if (markers != null && column == 3) {
                      return Align(
                        alignment: Alignment.centerLeft,
                        child: Column(mainAxisAlignment: MainAxisAlignment.center, crossAxisAlignment: CrossAxisAlignment.start, children: [
                          Text(
                            [c.markers.join(', '), if (c.length != null) '${c.length} bp'].join(' · '),
                            maxLines: 1,
                            overflow: TextOverflow.ellipsis,
                            style: TextStyle(fontSize: 12.5, fontWeight: FontWeight.w600, color: t.accentText),
                          ),
                          // Where it was found (NCBI, BOLD or both), with the
                          // BOLD country and BIN.
                          SplaceTooltip(
                            message: [
                              if (c.boldId != null) 'BOLD ${c.boldId}',
                              if (c.bin != null) 'BIN ${c.bin}',
                              if (c.country != null) c.country!,
                              if (c.genbank != null) 'GenBank ${c.genbank}',
                            ].join(' · '),
                            child: Text(
                              [c.sources.join(' + '), ?c.country].join(' · '),
                              maxLines: 1,
                              overflow: TextOverflow.ellipsis,
                              style: TextStyle(fontSize: 11.5, color: c.sources.contains('BOLD') ? SplaceRed.c700 : t.textMuted, fontWeight: FontWeight.w600),
                            ),
                          ),
                        ]),
                      );
                    }
                    final col = markers != null && column == 4 ? 3 : column;
                    final (organism, rest) = _split(c);
                    final (label, dot, fg) = badge(c);
                    return switch (col) {
                      0 => Align(
                          alignment: Alignment.centerLeft,
                          child: SizedBox(
                            height: 26,
                            child: FittedBox(child: Switch(value: c.selected, onChanged: (v) => setState(() => c.selected = v))),
                          ),
                        ),
                      1 => HighlightText(c.accession, query, style: body),
                      2 => SplaceTooltip(
                          message: c.title,
                          child: Column(mainAxisAlignment: MainAxisAlignment.center, crossAxisAlignment: CrossAxisAlignment.start, children: [
                            Row(children: [
                              Flexible(child: HighlightText(organism, query, style: body.copyWith(fontStyle: FontStyle.italic, fontWeight: FontWeight.w700))),
                              // An old name kept in the title (NCBI updates only the organism).
                              if (c.filedAs != null)
                                Flexible(
                                  child: Text.rich(
                                    TextSpan(style: body.copyWith(fontSize: 11.5, color: t.textMuted), children: [
                                      TextSpan(text: '  ${context.l10n.reviewFiledAs} '),
                                      TextSpan(text: c.filedAs, style: const TextStyle(fontStyle: FontStyle.italic)),
                                    ]),
                                    overflow: TextOverflow.ellipsis,
                                    maxLines: 1,
                                  ),
                                ),
                            ]),
                            HighlightText(
                              rest,
                              query,
                              style: body.copyWith(
                                color: c.unverified ? _unverifiedColor : t.textMuted,
                                fontWeight: c.unverified ? FontWeight.w600 : FontWeight.w400,
                              ),
                            ),
                          ]),
                        ),
                      _ => Align(
                          alignment: Alignment.centerLeft,
                          child: Container(
                            padding: const EdgeInsets.symmetric(horizontal: 10, vertical: 4),
                            decoration: BoxDecoration(
                              color: dot.withValues(alpha: _kindOf(c) == _Kind.partial ? 0.45 : 0.12),
                              borderRadius: BorderRadius.circular(8),
                            ),
                            child: Text(label, style: TextStyle(fontSize: 12, fontWeight: FontWeight.w700, color: fg)),
                          ),
                        ),
                    };
                  },
                ),
        ),
      ]),
    );
  }
}

/// `.taxonomy-stacked-bar`: share of complete, partial, UNVERIFIED and
/// other genomes in [genomes], with the percentages written inside.
class _StackedBar extends StatelessWidget {
  const _StackedBar({required this.label, required this.genomes});
  final String label;
  final List<GenomeCandidate> genomes;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final total = genomes.length;
    final parts = [
      for (final (kind, color, fg) in [
        (_Kind.complete, _completeColor, Colors.white),
        (_Kind.partial, _partialColor, _partialText),
        (_Kind.unverified, _unverifiedColor, Colors.white),
        (_Kind.other, _otherColor, Colors.white),
      ])
        (genomes.where((c) => _kindOf(c) == kind).length, color, fg),
    ];
    return Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
      Kicker(label, size: 11),
      const SizedBox(height: 5),
      ClipRRect(
        borderRadius: BorderRadius.circular(6),
        child: Container(
          height: 28,
          color: t.surface3,
          child: total == 0
              ? null
              : Row(children: [
                  for (final (n, color, fg) in parts)
                    if (n > 0)
                      Expanded(
                        flex: n,
                        child: Container(
                          color: color,
                          alignment: Alignment.center,
                          child: n / total >= 0.06
                              ? Text('${(100 * n / total).round()}%  ($n)',
                                  maxLines: 1, overflow: TextOverflow.clip, style: TextStyle(fontSize: 12, fontWeight: FontWeight.w700, color: fg))
                              : null,
                        ),
                      ),
                ]),
        ),
      ),
    ]);
  }
}

/// Organism and the rest of the NCBI title, without the "UNVERIFIED:"
/// prefix (the classification column shows it).
(String, String) _split(GenomeCandidate c) {
  final desc = c.description.replaceFirst(RegExp(r'^UNVERIFIED:\s*', caseSensitive: false), '');
  final organism = c.species.isNotEmpty ? c.species : desc.split(RegExp(r'\s+')).take(2).join(' ').replaceAll(RegExp(r'[,;:.]+$'), '');
  final rest = desc.startsWith(organism) ? desc.substring(organism.length).replaceFirst(RegExp(r'^[\s,;:]+'), '') : desc;
  return (organism, rest.isEmpty ? c.title : rest);
}
