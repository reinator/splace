import 'package:flutter/material.dart';

import '../../services/species_id.dart';
import '../common.dart';
import '../fa.dart';
import '../theme.dart';
import '../widgets/dialogs.dart';
import '../widgets/kit.dart';
import '../widgets/splace_tooltip.dart';
import 'id_alignment_dialog.dart';
import 'species_id_params.dart';

/// Runs the whole identification: choose the genes, send them to GenBank
/// and BOLD, then let the user pick a name per genome.
Future<void> runSpeciesIdentification(BuildContext context) async {
  final c = ProjectScope.read(context);
  final genes = await showDialog<List<String>>(
    context: context,
    builder: (_) => ProjectScope(controller: c, child: const _SetupDialog()),
  );
  if (genes == null || genes.isEmpty || !context.mounted) return;
  final ok = await c.identifySpecies(genes);
  // Let the progress modal show its result before the list opens.
  await Future<void>.delayed(const Duration(milliseconds: 1500));
  if (ok && context.mounted) await showSpeciesIdResults(context);
}

Future<void> showSpeciesIdResults(BuildContext context) => showDialog(
      context: context,
      barrierDismissible: false,
      builder: (_) => ProjectScope(controller: ProjectScope.read(context), child: const _ResultsDialog()),
    );

/// Gene choice: every gene goes to GenBank; BOLD only takes COI.
class _SetupDialog extends StatefulWidget {
  const _SetupDialog();

  @override
  State<_SetupDialog> createState() => _SetupDialogState();
}

class _SetupDialogState extends State<_SetupDialog> {
  // The markers BOLD identifies with its database (COI; rbcL and matK for
  // plants) come selected.
  late final Set<String> _genes = {
    for (final g in ProjectScope.read(context).identificationMarkers.keys)
      if (ProjectScope.read(context).speciesIdOptions.boldAccepts(g)) g,
  };

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    final markers = c.identificationMarkers;
    final o = c.speciesIdOptions;
    final bold = _genes.any(o.boldAccepts);
    final cover = c.speciesIdCacheCoverage(_genes.toList());

    return SplaceDialog(
      icon: FA.dna,
      title: l10n.speciesIdSetupTitle,
      subtitle: l10n.speciesIdSetupSubtitle,
      actions: [
        DialogCancelButton(onPressed: () => Navigator.pop(context)),
        FilledButton(onPressed: _genes.isEmpty ? null : () => Navigator.pop(context, _genes.toList()), child: Text(l10n.speciesIdRun)),
      ],
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        Row(crossAxisAlignment: CrossAxisAlignment.start, children: [
          Expanded(
            child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
              Kicker(l10n.speciesIdGenesLabel, size: 12),
              const SizedBox(height: 4),
              Text(l10n.speciesIdGenesHintFor(o.boldGenes.join(', ')), style: TextStyle(fontSize: 13, color: t.textMuted)),
            ]),
          ),
          // Searches already answered are reused, so adding a gene only
          // sends the new ones (on the same line: the dialog never scrolls).
          if (cover.cached > 0) ...[
            const SizedBox(width: 12),
            Container(
              padding: const EdgeInsets.only(left: 10),
              decoration: BoxDecoration(color: SplaceTokens.success.withValues(alpha: 0.08), borderRadius: BorderRadius.circular(8)),
              child: Row(mainAxisSize: MainAxisSize.min, children: [
                Icon(FA.circleCheck, size: 12, color: SplaceTokens.success),
                const SizedBox(width: 6),
                Text(
                  cover.cached == cover.total
                      ? l10n.siCacheAll('${cover.total}')
                      : l10n.siCacheInfo('${cover.cached}', '${cover.total}', '${cover.total - cover.cached}'),
                  style: TextStyle(fontSize: 12, color: t.text),
                ),
                TextButton(onPressed: c.clearSpeciesIdCache, child: Text(l10n.siCacheClear)),
              ]),
            ),
          ],
        ]),
        const SizedBox(height: 8),
        Wrap(spacing: 8, runSpacing: 8, children: [
          for (final e in markers.entries)
            ChipButton(
              label: '${e.key}  ·  ${l10n.speciesIdInGenomes('${e.value}')}',
              style: _genes.contains(e.key) ? ChipButtonStyle.accent : ChipButtonStyle.neutral,
              onPressed: () => setState(() => _genes.contains(e.key) ? _genes.remove(e.key) : _genes.add(e.key)),
            ),
        ]),
        const SizedBox(height: 10),
        // Each source with its parameters, always in view (both cards the
        // same height); a source that will not run is dimmed.
        SpeciesIdParams(
          blastNote: l10n.speciesIdGenbankBodyFor(o.blastProgram, o.blastDatabase, '${o.hitlistSize}'),
          blastOn: _genes.isNotEmpty,
          boldNote: bold ? l10n.speciesIdBoldBodyFor(o.boldDatabase, o.boldGenes.join(', ')) : l10n.speciesIdBoldOffFor(o.boldGenes.join(', ')),
          boldOn: bold,
        ),
      ]),
    );
  }
}

/// One card per genome: a verdict in plain words and a table of the
/// species GenBank and BOLD point to, with the identity of each source; the
/// user keeps the current name or picks another.
class _ResultsDialog extends StatefulWidget {
  const _ResultsDialog();

  @override
  State<_ResultsDialog> createState() => _ResultsDialogState();
}

String _id(double v) => v.toStringAsFixed(1);

/// ≥ 98 % same species, 95–98 % close relative, below another species.
Color _identityColor(BuildContext context, double v) =>
    v >= 98 ? context.tokens.accent : (v >= 95 ? SplaceTokens.warn : context.tokens.textMuted);

class _ResultsDialogState extends State<_ResultsDialog> {
  final _choice = <String, String>{};
  bool _onlyDisagree = false;

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    final ids = c.speciesIds.values.toList()
      ..sort((a, b) => a.disagrees == b.disagrees ? a.currentName.compareTo(b.currentName) : (a.disagrees ? -1 : 1));
    final disagree = ids.where((r) => r.disagrees).length;
    final shown = _onlyDisagree ? ids.where((r) => r.disagrees).toList() : ids;
    final changes = {
      for (final e in _choice.entries)
        if (!sameSpecies(c.speciesIds[e.key]?.currentName ?? '', e.value)) e.key: e.value,
    };
    final head = TextStyle(fontSize: 12.5, fontWeight: FontWeight.w700, color: t.textMuted);
    // Use (fixed) · Species · GenBank · BOLD · Sequences (flex).
    const flex = [0, 5, 5, 2, 2];

    Widget identity(double v) => Text('${_id(v)}%', style: TextStyle(fontSize: 13, fontWeight: FontWeight.w700, color: _identityColor(context, v)));

    Widget row(RecordIdentification r, String species, SpeciesCandidate? s) {
      final chosen = sameSpecies(_choice[r.accession] ?? r.currentName, species);
      final isCurrent = r.isCurrent(species);
      return InkWell(
        onTap: () => setState(() => _choice[r.accession] = species),
        child: Container(
          padding: const EdgeInsets.symmetric(horizontal: 10, vertical: 7),
          decoration: BoxDecoration(
            color: chosen ? t.accentSoft : null,
            border: Border(top: BorderSide(color: t.borderSoft)),
          ),
          child: Row(crossAxisAlignment: CrossAxisAlignment.center, children: [
            SizedBox(width: 44, child: Icon(chosen ? FA.circleCheck : FA.circle, size: 15, color: chosen ? t.accent : t.surface3)),
            Expanded(
              flex: flex[1],
              child: Wrap(spacing: 6, crossAxisAlignment: WrapCrossAlignment.center, children: [
                Text(species, style: TextStyle(fontSize: 13.5, fontStyle: FontStyle.italic, fontWeight: FontWeight.w700, color: chosen ? t.accentText : t.text)),
                if (isCurrent) TagChip(l10n.siCurrentTag, color: t.textMuted),
                // The alignments of the query with this species' hits.
                if (s != null && s.hits.isNotEmpty)
                  IconButton(
                    tooltip: l10n.idAlnOpen,
                    visualDensity: VisualDensity.compact,
                    icon: Icon(FA.alignLeft, size: 13, color: t.accent),
                    onPressed: () => showIdAlignments(context, r, species: species),
                  ),
                // Hits that came under synonyms of this species.
                if (s != null && s.foundAs.isNotEmpty)
                  Text(l10n.siFoundAs(s.foundAs.join(', ')), style: TextStyle(fontSize: 11.5, fontStyle: FontStyle.italic, color: t.textMuted)),
              ]),
            ),
            // One line per gene, so any number of genes stays readable.
            Expanded(
              flex: flex[2],
              child: s == null || s.genbank.isEmpty
                  ? Text(s == null ? l10n.siNotFound : '—', style: TextStyle(fontSize: 12.5, color: t.textFaint))
                  : Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
                      for (final e in (s.genbank.entries.toList()..sort((a, b) => a.key.compareTo(b.key))))
                        Padding(
                          padding: const EdgeInsets.symmetric(vertical: 1),
                          child: Row(children: [
                            SizedBox(width: 64, child: Text(e.key, style: TextStyle(fontSize: 12, color: t.textMuted))),
                            if (e.value.onlySelf) ...[
                              // Only the record itself: shown, not evidence.
                              Text('${_id(e.value.selfIdentity!)}%', style: TextStyle(fontSize: 13, fontWeight: FontWeight.w700, color: t.textMuted)),
                              const SizedBox(width: 6),
                              SplaceTooltip(
                                message: l10n.siSelfInfo(e.value.selfAccessions.join(', ')),
                                child: TagChip(l10n.siSelfTag, color: t.textMuted),
                              ),
                            ] else ...[
                              identity(e.value.identity),
                              const SizedBox(width: 6),
                              Text('(${e.value.hits})', style: TextStyle(fontSize: 11.5, color: t.textFaint)),
                            ],
                          ]),
                        ),
                    ]),
            ),
            Expanded(
              flex: flex[3],
              child: s?.bold != null
                  ? identity(s!.bold!.identity)
                  : Text(r.sources.contains(IdSource.bold) ? '—' : l10n.siNotFound, style: TextStyle(fontSize: 12.5, color: t.textFaint)),
            ),
            Expanded(
              flex: flex[4],
              child: Text(s == null ? '0' : '${s.totalHits}', style: TextStyle(fontSize: 13, color: t.text)),
            ),
          ]),
        ),
      );
    }

    Widget genome(RecordIdentification r) {
      final cands = r.candidates;
      final best = cands.firstOrNull;
      final cur = r.current;
      final selfNote = r.onlySelfInGenbank ? ' ${l10n.siVerdictSelfOnly(r.selfAccessions.join(', '), r.currentName)}' : '';
      final (verdict, badge, color) = best == null || best.bestIdentity == 0
          ? (r.onlySelfInGenbank ? (selfNote.trim(), l10n.siNoReference, t.textMuted) : (l10n.siVerdictNone, l10n.speciesIdNoHits, t.textMuted))
          : r.disagrees
              ? cur != null && cur.bestIdentity > 0
                  ? (l10n.siVerdictMismatch(best.species, _id(best.bestIdentity), r.currentName, _id(cur.bestIdentity)), l10n.speciesIdMismatch, SplaceTokens.warn)
                  : (l10n.siVerdictMissing(best.species, _id(best.bestIdentity), r.currentName), l10n.speciesIdMismatch, SplaceTokens.warn)
              : r.isCurrent(best.species)
                  ? (l10n.siVerdictConfirmed(_id(best.bestIdentity)) + (r.currentAccepted == null ? '' : ' ${l10n.siVerdictAccepted(r.currentAccepted!)}'),
                      l10n.speciesIdMatch, t.accent)
                  // The closest match is only a relative: nothing to
                  // compare the name with.
                  : (l10n.siVerdictNoReference(r.currentName, best.species, _id(best.bestIdentity)) + selfNote, l10n.siNoReference, t.textMuted);
      return Container(
        margin: const EdgeInsets.only(bottom: 10),
        decoration: BoxDecoration(
          color: t.surface,
          borderRadius: BorderRadius.circular(12),
          border: Border.all(color: r.disagrees ? SplaceTokens.warn.withValues(alpha: 0.5) : t.border, width: r.disagrees ? 1.5 : 1),
        ),
        clipBehavior: Clip.antiAlias,
        child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
          Container(
            padding: const EdgeInsets.fromLTRB(12, 10, 12, 10),
            color: r.disagrees ? SplaceTokens.warn.withValues(alpha: 0.06) : t.surface2,
            child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
              Row(children: [
                Text(r.currentName, style: TextStyle(fontSize: 15, fontStyle: FontStyle.italic, fontWeight: FontWeight.w700, color: t.text)),
                const SizedBox(width: 8),
                Text(r.accession, style: TextStyle(fontSize: 12.5, color: t.textMuted)),
                const Spacer(),
                TagChip(badge, color: color),
              ]),
              const SizedBox(height: 4),
              // The species named in the verdict, always in italics.
              Text.rich(italicNames(verdict, [
                r.currentName,
                r.currentAccepted,
                best?.species,
                cur?.species,
                for (final k in r.shortlist()) k.species,
              ], style: TextStyle(fontSize: 13, height: 1.4, color: r.disagrees ? SplaceTokens.warn : t.text))),
            ]),
          ),
          ...[
            Padding(
              padding: const EdgeInsets.symmetric(horizontal: 10, vertical: 6),
              child: Row(children: [
                SizedBox(width: 44, child: Text(l10n.siColUse, style: head)),
                Expanded(flex: flex[1], child: Text(l10n.siColSpecies, style: head)),
                Expanded(flex: flex[2], child: Text(l10n.siColGenbank, style: head)),
                Expanded(flex: flex[3], child: Text(l10n.siColBold, style: head)),
                Expanded(flex: flex[4], child: Text(l10n.siColSeqs, style: head)),
              ]),
            ),
            for (final s in r.shortlist()) row(r, s.species, s),
            if (cur == null) row(r, r.currentName, null),
          ],
        ]),
      );
    }

    Widget legendDot(Color color, String text) => Row(mainAxisSize: MainAxisSize.min, children: [
          Container(width: 10, height: 10, decoration: BoxDecoration(color: color, shape: BoxShape.circle)),
          const SizedBox(width: 6),
          Text(text, style: const TextStyle(fontSize: 13)),
        ]);

    return SplaceDialog(
      icon: FA.dna,
      title: l10n.speciesIdResultsTitle,
      subtitle: l10n.speciesIdResultsSubtitle(c.speciesIdGenes.join(', ')),
      footerLeading: Text(l10n.speciesIdNote, style: TextStyle(fontSize: 12, color: t.textMuted)),
      actions: [
        DialogCancelButton(onPressed: () => Navigator.pop(context)),
        FilledButton(
          onPressed: changes.isEmpty
              ? null
              : () {
                  c.applySpeciesChoices(changes);
                  Navigator.pop(context);
                },
          child: Text(l10n.speciesIdApply('${changes.length}')),
        ),
      ],
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        TonePanel(
          tone: PanelTone.info,
          title: l10n.siHowTitle,
          icon: FA.circleInfo,
          child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
            Text(l10n.siHowBody, style: const TextStyle(fontSize: 13)),
            const SizedBox(height: 8),
            Wrap(spacing: 20, runSpacing: 6, children: [
              legendDot(t.accent, l10n.siLegendSame),
              legendDot(SplaceTokens.warn, l10n.siLegendClose),
              legendDot(t.textMuted, l10n.siLegendOther),
            ]),
          ]),
        ),
        // A source that failed is said plainly, not left as empty columns.
        for (final e in c.speciesIdFailures.entries) ...[
          const SizedBox(height: 10),
          TonePanel(
            tone: PanelTone.danger,
            child: Text(l10n.siSourceFailed(e.key == IdSource.genbank ? 'GenBank (BLAST)' : 'BOLD', e.value), style: const TextStyle(fontSize: 13)),
          ),
        ],
        const SizedBox(height: 12),
        Row(children: [
          Expanded(
            child: Text(l10n.speciesIdSummary('$disagree', '${ids.length}'),
                style: TextStyle(fontSize: 14, fontWeight: FontWeight.w600, color: disagree > 0 ? SplaceTokens.warn : t.accentText)),
          ),
          SizedBox(
            width: 340,
            child: SwitchBox(label: l10n.speciesIdOnlyDisagree, value: _onlyDisagree, onChanged: (v) => setState(() => _onlyDisagree = v)),
          ),
        ]),
        const SizedBox(height: 12),
        for (final r in shown) genome(r),
      ]),
    );
  }
}
