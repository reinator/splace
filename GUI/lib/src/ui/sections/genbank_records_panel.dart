import 'package:flutter/material.dart';
import 'package:url_launcher/url_launcher.dart';

import '../../core/genbank.dart';
import '../../core/genbank_markers.dart';
import '../../core/taxonomy_sources.dart';
import '../../services/species_id.dart' show sameSpecies;
import 'feature_types_panel.dart' show showGeneReviewDialog;
import '../common.dart';
import '../fa.dart';
import '../theme.dart';
import '../chrome/species_id_dialogs.dart';
import '../widgets/citations.dart';
import '../widgets/dialogs.dart';
import '../widgets/kit.dart';
import '../widgets/salve_widgets.dart';
import '../widgets/splace_tooltip.dart';

const recordColumns = [
  'accession', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'authorship', 'status', //
  ...salveColumns, 'pcgs', 'rrnas', 'trnas', 'source',
];

/// Columns filled from the SALVE sheet (shown once a record has one).
const salveColumns = [
  'criteria', 'commonNames', 'salveGroup', 'trend', 'generation', 'endemic', 'migratory', 'biomes', 'states', 'habitats', 'threats', //
  'officialList', 'ucs', 'assessed', 'doi', 'sheet',
];

/// GenBank part of step 2 (`#recordsSection`): taxonomy actions, the
/// taxonomic tree, column toggles, search and the records table.
class GenBankRecordsPanel extends StatefulWidget {
  const GenBankRecordsPanel({super.key});

  @override
  State<GenBankRecordsPanel> createState() => _GenBankRecordsPanelState();
}

class _GenBankRecordsPanelState extends State<GenBankRecordsPanel> {
  String _filter = '';
  bool _prompted = false;

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;

    // New records with unknown gene names or /gene vs /product conflicts
    // open the gene review once, before anything uses the names.
    if (!c.unknownGenesReviewed && !_prompted && c.catalog != null && !c.busy && c.modal == null) {
      _prompted = true;
      WidgetsBinding.instance.addPostFrameCallback((_) {
        if (mounted) showGeneReviewDialog(context);
      });
    }
    if (c.unknownGenesReviewed) _prompted = false;

    final sorted = List.of(c.records)..sort((a, b) => a.organismName.compareTo(b.organismName));
    final q = _filter.toLowerCase().trim();
    final rows = q.isEmpty
        ? sorted
        : sorted.where((r) {
            final tx = r.taxonomy;
            return '${r.accession} ${r.organismName} ${tx.family} ${tx.order} ${tx.klass} ${tx.phylum} ${tx.kingdom} ${r.fileName ?? ''}'
                .toLowerCase()
                .contains(q);
          }).toList();
    final counts = [for (final r in c.records) r.featureCounts];
    final pcg = MarkerCountStats.of(counts.map((e) => e.pcgs));
    final rrna = MarkerCountStats.of(counts.map((e) => e.rrnas));
    final trna = MarkerCountStats.of(counts.map((e) => e.trnas));

    String colLabel(String k) => switch (k) {
      'accession' => l10n.colAccession,
      'kingdom' => l10n.colKingdom,
      'phylum' => l10n.colPhylum,
      'class' => l10n.colClass,
      'order' => l10n.colOrder,
      'family' => l10n.colFamily,
      'genus' => l10n.colGenus,
      'species' => l10n.colSpecies,
      'authorship' => l10n.colAuthorship,
      'pcgs' => 'PCGs',
      'rrnas' => 'rRNAs',
      'trnas' => 'tRNAs',
      'status' => l10n.colStatus,
      'criteria' => l10n.colCriteria,
      'commonNames' => l10n.colCommonNames,
      'salveGroup' => l10n.colSalveGroup,
      'endemic' => l10n.colEndemic,
      'biomes' => l10n.colBiomes,
      'states' => l10n.colStates,
      'assessed' => l10n.colAssessed,
      'doi' => 'DOI',
      'sheet' => l10n.colSheet,
      'trend' => l10n.colTrend,
      'generation' => l10n.colGeneration,
      'migratory' => l10n.colMigratory,
      'habitats' => l10n.colHabitats,
      'threats' => l10n.colThreats,
      'officialList' => l10n.colOfficialList,
      'ucs' => l10n.colUcs,
      _ => l10n.colSource,
    };
    final hasStatus = c.records.any((r) => r.conservation != null);
    final hasSalve = c.records.any((r) => r.salve != null);
    final columns = [
      for (final k in recordColumns)
        if ((k != 'status' || hasStatus) && (!salveColumns.contains(k) || hasSalve)) k,
    ];
    final visible = [
      for (final k in columns)
        if (!c.hiddenColumns.contains(k)) k,
    ];
    int flexOf(String k) => switch (k) {
      'species' || 'authorship' || 'commonNames' || 'biomes' || 'states' || 'status' => 3,
      'source' || 'sheet' || 'endemic' => 1,
      _ => 2,
    };

    double minWidthOf(String k) => switch (k) {
      'species' || 'commonNames' => 190,
      'status' || 'biomes' || 'states' || 'doi' || 'authorship' => 170,
      'accession' || 'family' || 'order' || 'class' || 'phylum' || 'kingdom' || 'genus' || 'salveGroup' => 110,
      'criteria' || 'assessed' || 'endemic' || 'source' || 'migratory' || 'officialList' || 'ucs' => 104,
      'trend' || 'generation' => 130,
      'habitats' || 'threats' => 200,
      'sheet' => 64,
      _ => 70,
    };

    Widget cell(GbRecord r, String k) {
      final tx = r.taxonomy;
      final fc = r.featureCounts;
      String dash(String v) => v.isEmpty ? '—' : v;
      // Text cells highlight what the search box matched.
      Widget hl(String v, TextStyle style) => HighlightText(v, _filter, style: style);
      final muted = TextStyle(fontSize: 13, color: t.textMuted);
      return switch (k) {
        'accession' => hl(r.accession, TextStyle(fontSize: 13, color: t.text)),
        'kingdom' => hl(dash(tx.kingdom), muted),
        'phylum' => hl(dash(tx.phylum), muted),
        'class' => hl(dash(tx.klass), muted),
        'order' => hl(dash(tx.order), muted),
        'family' => hl(dash(tx.family), muted),
        'genus' => hl(tx.genus, muted.copyWith(fontStyle: FontStyle.italic)),
        'species' => SplaceTooltip(
          message: '${r.organismName}\n${r.definition}',
          child: hl(r.organismName, TextStyle(fontSize: 13, fontStyle: FontStyle.italic, color: t.text)),
        ),
        'authorship' => hl(dash(tx.authorship), muted),
        'status' => r.salve != null
            ? Align(
                alignment: Alignment.centerLeft,
                child: IucnBadge(code: r.salve!.code, name: r.salve!.categoryName, possiblyExtinct: r.salve!.possiblyExtinct),
              )
            : r.conservation == null
                ? Text('—', style: TextStyle(color: t.textMuted))
                : Align(alignment: Alignment.centerLeft, child: TagChip(r.conservation!, color: t.textMuted)),
        'criteria' => hl(dash(r.salve?.criteria ?? ''), muted),
        'commonNames' => SplaceTooltip(message: r.salve?.commonNames ?? '', child: hl(dash(r.salve?.commonNames ?? ''), muted)),
        'salveGroup' => hl(dash(r.salve?.group ?? ''), muted),
        'endemic' => hl(dash(r.salve?.endemic ?? ''), muted),
        'biomes' => SplaceTooltip(message: r.salve?.biomes ?? '', child: hl(dash(r.salve?.biomes ?? ''), muted)),
        'states' => SplaceTooltip(message: r.salve?.states ?? '', child: hl(dash(r.salve?.states ?? ''), muted)),
        'assessed' => hl(dash(r.salve?.assessed ?? ''), muted),
        'trend' => hl(dash(r.salve?.populationTrend ?? ''), muted),
        'generation' => hl(dash(r.salve?.generationTime ?? ''), muted),
        'migratory' => hl(dash(r.salve?.migratory ?? ''), muted),
        'officialList' => SplaceTooltip(message: (r.salve?.officialLists ?? const []).join('\n'), child: hl(dash(r.salve?.officialList ?? ''), muted)),
        'habitats' => SplaceTooltip(message: (r.salve?.habitats ?? const []).join('\n'), child: hl(dash((r.salve?.habitats ?? const []).join(', ')), muted)),
        'threats' => SplaceTooltip(
            message: (r.salve?.threats ?? const []).join('\n'),
            child: hl((r.salve?.threats ?? const []).isEmpty ? '—' : '${r.salve!.threats.length} · ${r.salve!.threats.join(', ')}', muted)),
        'ucs' => SplaceTooltip(
            message: (r.salve?.protectedAreas ?? const []).take(40).join('\n'),
            child: hl(r.salve == null ? '—' : '${r.salve!.protectedAreas.length}', muted)),
        'doi' => (r.salve?.doi ?? '').isEmpty
            ? Text('—', style: muted)
            : MouseRegion(
                cursor: SystemMouseCursors.click,
                child: GestureDetector(
                  onTap: () => launchUrl(Uri.parse('https://doi.org/${r.salve!.doi}')),
                  child: Text(r.salve!.doi, style: TextStyle(fontSize: 13, color: t.accentText, decoration: TextDecoration.underline)),
                ),
              ),
        'sheet' => r.salve == null
            ? Text('—', style: muted)
            : Align(
                alignment: Alignment.centerLeft,
                child: IconBtn(
                  icon: FA.circleInfo,
                  tooltip: l10n.salveSheetSubtitle,
                  size: 26,
                  onPressed: () => showSalveSheetDialog(context, r.organismName, r.salve!),
                ),
              ),
        'pcgs' => SplaceTooltip(
          message: 'PCGs: ${fc.pcgs} (mode ${pcg.mode})',
          child: Text('${fc.pcgs}', style: TextStyle(fontWeight: FontWeight.w700, color: CountCell.colors(context, fc.pcgs, pcg.mode, pcg.min, pcg.max).$2)),
        ),
        'rrnas' => SplaceTooltip(
          message: 'rRNAs: ${fc.rrnas} (mode ${rrna.mode})',
          child: Text('${fc.rrnas}', style: TextStyle(fontWeight: FontWeight.w700, color: CountCell.colors(context, fc.rrnas, rrna.mode, rrna.min, rrna.max).$2)),
        ),
        'trnas' => SplaceTooltip(
          message: 'tRNAs: ${fc.trnas} (mode ${trna.mode})',
          child: Text('${fc.trnas}', style: TextStyle(fontWeight: FontWeight.w700, color: CountCell.colors(context, fc.trnas, trna.mode, trna.min, trna.max).$2)),
        ),
        _ => Align(
          alignment: Alignment.centerLeft,
          child: TagChip(r.source == 'ncbi' ? 'NCBI' : l10n.sourceFile, color: r.source == 'ncbi' ? t.accent : t.textMuted),
        ),
      };
    }

    return Column(
      crossAxisAlignment: CrossAxisAlignment.stretch,
      children: [
        Row(
          children: [
            const Spacer(),
            if (c.recordsWithoutMarkers.isNotEmpty) ...[
              ChipButton(
                label: l10n.noMarkersButton('${c.recordsWithoutMarkers.length}'),
                style: ChipButtonStyle.danger,
                onPressed: c.busy ? null : () => _showNoMarkersDialog(context),
              ),
              const SizedBox(width: 8),
            ],
            ChipButton(
              label: l10n.clearAll,
              icon: FA.trashCan,
              style: ChipButtonStyle.danger,
              onPressed: c.busy
                  ? null
                  : () async {
                      if (await confirmDialog(context, title: l10n.clearTitle, body: l10n.clearBody, confirm: l10n.clearAll)) c.clearRecords();
                    },
            ),
          ],
        ),
        const SizedBox(height: 10),
        const SubAnchor(id: '2a', child: _SpeciesIdCard()),
        const SizedBox(height: 14),
        const SubAnchor(id: '2b', child: _TaxonomyCard()),
        if (c.records.any((r) => r.hasFetchedTaxonomy)) ...[
          const SizedBox(height: 14),
          _TaxonomyTree(root: buildTaxonomyTree(c.records), records: c.records),
        ],
        const SizedBox(height: 18),
        SubAnchor(id: '2c', child: SubHeading(id: '2c', title: l10n.subRecordsTable, body: l10n.subRecordsTableBody)),
        Wrap(
          spacing: 6,
          runSpacing: 6,
          crossAxisAlignment: WrapCrossAlignment.center,
          children: [
            Row(
              mainAxisSize: MainAxisSize.min,
              children: [
                Icon(FA.tableColumns, size: 12, color: t.textMuted),
                const SizedBox(width: 6),
                Text(
                  l10n.columnsLabel,
                  style: TextStyle(fontSize: 13, color: t.textMuted, fontWeight: FontWeight.w600),
                ),
              ],
            ),
            for (final k in columns)
              ChipButton(
                label: colLabel(k),
                icon: c.hiddenColumns.contains(k) ? FA.eyeSlash : FA.eye,
                style: c.hiddenColumns.contains(k) ? ChipButtonStyle.neutral : ChipButtonStyle.soft,
                keepIcon: true,
                onPressed: () => c.setColumnVisible(k, c.hiddenColumns.contains(k)),
              ),
          ],
        ),
        const SizedBox(height: 12),
        TextField(
          decoration: InputDecoration(
            prefixIcon: Icon(FA.magnifyingGlass, size: 14, color: t.textFaint),
            hintText: l10n.recordsSearch,
          ),
          onChanged: (v) => setState(() => _filter = v),
        ),
        const SizedBox(height: 12),
        RecordsTable(
          columns: [...visible.map(colLabel), ''],
          flex: [...visible.map(flexOf), 2],
          // Starting widths that fit their content; the table scrolls
          // sideways when they do not all fit.
          minWidths: [...visible.map(minWidthOf), 72],
          rowCount: rows.length,
          rowHeight: 38,
          height: 460,
          // Drag a header edge to widen a column and read all of it.
          resizable: true,
          // Gene counts colour the whole cell: red below the mode, blue above.
          cellColor: (i, col) {
            if (col >= visible.length) return null;
            final fc = rows[i].featureCounts;
            return switch (visible[col]) {
              'pcgs' => CountCell.colors(context, fc.pcgs, pcg.mode, pcg.min, pcg.max).$1,
              'rrnas' => CountCell.colors(context, fc.rrnas, rrna.mode, rrna.min, rrna.max).$1,
              'trnas' => CountCell.colors(context, fc.trnas, trna.mode, trna.min, trna.max).$1,
              _ => null,
            };
          },
          cellBuilder: (context, i, col) {
            final r = rows[i];
            if (col < visible.length) return cell(r, visible[col]);
            return Wrap(
              alignment: WrapAlignment.end,
              spacing: 2,
              children: [
                IconBtn(
                  icon: FA.penToSquare,
                  tooltip: l10n.editRecordTooltip(r.organismName, r.accession),
                  size: 26,
                  onPressed: () => showEditRecordDialog(context, r),
                ),
                IconBtn(
                  icon: FA.trashCan,
                  tooltip: l10n.removeRecordTooltip(r.organismName, r.accession),
                  size: 26,
                  onPressed: () async {
                    if (await confirmDialog(
                      context,
                      title: l10n.removeTitle,
                      body: '${r.organismName} (${r.accession}) — ${l10n.removeBody}',
                      confirm: l10n.remove,
                    )) {
                      c.removeRecord(r);
                    }
                  },
                ),
              ],
            );
          },
        ),
      ],
    );
  }
}

String taxonomySourceName(TaxonomySource s) => switch (s) {
  TaxonomySource.gbif => 'GBIF',
  TaxonomySource.ncbi => 'NCBI Taxonomy',
  TaxonomySource.salve => 'ICMBio SALVE',
};

/// Species-identification card of step 2: BLAST/BOLD check of the names
/// before the taxonomy, which depends on them. Optional.
class _SpeciesIdCard extends StatelessWidget {
  const _SpeciesIdCard();

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    final ids = c.speciesIds.values.toList();
    final flagged = ids.where((r) => r.disagrees).length;
    final confirmed = ids.where((r) {
      final best = r.candidates.firstOrNull;
      return !r.disagrees && best != null && best.bestIdentity > 0 && sameSpecies(best.species, r.currentName);
    }).length;
    final unverified = ids.length - flagged - confirmed;
    final pending = (c.catalog?.unknown.length ?? 0) + (c.geneReview?.conflicts.length ?? 0);

    Widget count(IconData icon, String text, Color color) => Row(mainAxisSize: MainAxisSize.min, children: [
          Icon(icon, size: 12, color: color),
          const SizedBox(width: 6),
          Text(text, style: TextStyle(fontSize: 13, fontWeight: FontWeight.w600, color: color)),
        ]);

    return _StepCard(
      icon: FA.barcode,
      badge: '2a',
      title: l10n.siCardTitle,
      subtitle: Text(l10n.siCardBody, style: TextStyle(fontSize: 13, color: t.textMuted, height: 1.6)),
      trailing: TagChip(l10n.siCardOptional, color: t.textMuted),
      body: ids.isEmpty
          ? null
          : Wrap(spacing: 18, runSpacing: 6, children: [
              count(FA.listCheck, l10n.siSumChecked('${ids.length}'), t.text),
              count(FA.circleCheck, l10n.siSumConfirmed('$confirmed'), t.accent),
              if (flagged > 0) count(FA.triangleExclamation, l10n.siSumFlagged('$flagged'), SplaceTokens.warn),
              if (unverified > 0) count(FA.circleMinus, l10n.siSumNoReference('$unverified'), t.textMuted),
            ]),
      footer: [
        Icon(FA.globe, size: 13, color: t.textFaint),
        const SizedBox(width: 8),
        Expanded(child: Text(l10n.internetRequired, maxLines: 2, overflow: TextOverflow.ellipsis, style: TextStyle(fontSize: 12, color: t.textMuted))),
        const SizedBox(width: 12),
        if (pending > 0) ...[
          ChipButton(
            label: c.unknownGenesReviewed ? l10n.geneReviewAgain('$pending') : l10n.geneReviewButtonN('$pending'),
            icon: c.unknownGenesReviewed ? FA.circleCheck : null,
            keepIcon: c.unknownGenesReviewed,
            style: c.unknownGenesReviewed ? ChipButtonStyle.soft : ChipButtonStyle.danger,
            onPressed: c.busy ? null : () => showGeneReviewDialog(context),
          ),
          const SizedBox(width: 10),
        ],
        if (ids.isNotEmpty) ...[
          OutlinedButton(onPressed: c.busy ? null : () => showSpeciesIdResults(context), child: Text(l10n.speciesIdShowResults)),
          const SizedBox(width: 10),
        ],
        PrimaryButton(
          label: ids.isEmpty ? l10n.speciesIdButton : l10n.siCardRerun,
          icon: FA.barcode,
          onPressed: c.busy ? null : () => runSpeciesIdentification(context),
        ),
      ],
    );
  }
}

/// Card of step 2: icon, title and subtitle, optional body, and a footer
/// bar with the actions.
class _StepCard extends StatelessWidget {
  const _StepCard({required this.badge, required this.icon, required this.title, required this.subtitle, this.trailing, this.body, this.bodyInset = true, required this.footer});
  final String badge;
  final IconData icon;
  final String title;
  final Widget subtitle;
  final Widget? trailing;
  final Widget? body;

  /// The body lines up with the title (true) or spans the card.
  final bool bodyInset;
  final List<Widget> footer;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final dark = Theme.of(context).brightness == Brightness.dark;
    return Container(
      decoration: BoxDecoration(
        borderRadius: BorderRadius.circular(18),
        border: Border.all(color: t.accentSoftBorder),
        gradient: dark
            ? null
            : const LinearGradient(
                begin: Alignment.topLeft,
                end: Alignment.bottomRight,
                colors: [Color(0xFFEEF0FB), Colors.white, Color(0xFFFEF2F2)],
                stops: [0, .55, 1],
              ),
        color: dark ? t.surface2 : null,
        boxShadow: shadowSm,
      ),
      clipBehavior: Clip.antiAlias,
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        Padding(
          padding: const EdgeInsets.fromLTRB(20, 18, 20, 0),
          child: Row(crossAxisAlignment: CrossAxisAlignment.start, children: [
            Container(
              width: 42,
              height: 42,
              decoration: BoxDecoration(color: t.accent, borderRadius: BorderRadius.circular(12)),
              child: Icon(icon, size: 17, color: Colors.white),
            ),
            const SizedBox(width: 14),
            Expanded(
              child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
                Row(children: [
                  StepBadge(badge, size: 24),
                  const SizedBox(width: 8),
                  Flexible(child: Text(title, style: TextStyle(fontSize: 16, fontWeight: FontWeight.w600, color: t.accentText))),
                ]),
                const SizedBox(height: 2),
                subtitle,
              ]),
            ),
            if (trailing != null) ...[const SizedBox(width: 10), trailing!],
          ]),
        ),
        if (body != null) Padding(padding: EdgeInsets.fromLTRB(bodyInset ? 76 : 20, bodyInset ? 12 : 16, 20, 0), child: body),
        const SizedBox(height: 16),
        Container(
          padding: const EdgeInsets.fromLTRB(20, 12, 20, 12),
          decoration: BoxDecoration(
            color: dark ? t.surface : Colors.white.withValues(alpha: 0.7),
            border: Border(top: BorderSide(color: t.accentSoftBorder)),
          ),
          child: Row(children: footer),
        ),
      ]),
    );
  }
}

/// Taxonomy card of step 2 (`.step2-intro-card` of the HTML app): choose
/// the source (GBIF, NCBI Taxonomy or an ICMBio SALVE table), fetch or
/// skip, with the dataFishing reference.
class _TaxonomyCard extends StatelessWidget {
  const _TaxonomyCard();

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    final withTaxId = c.records.where((r) => r.taxId != null).length;
    final canFetch = !c.busy;

    Widget source(TaxonomySource s, IconData icon, String body, String badge) {
      final selected = c.taxonomySource == s;
      return _SourceTile(selected: selected, icon: icon, title: taxonomySourceName(s), body: body, badge: badge, onTap: () => c.setTaxonomySource(s));
    }

    return _StepCard(
      icon: FA.sitemap,
      badge: '2b',
      title: l10n.taxonomyCardTitle,
      // One sentence, with dataFishing as an inline pill.
      subtitle: Text.rich(
        TextSpan(style: TextStyle(fontSize: 13, color: t.textMuted, height: 1.6), children: citationSpans(l10n.recordsReviewMethod)),
      ),
      trailing: c.taxonomyResolved ? TagChip(l10n.taxonomyResolved, icon: FA.check, color: SplaceTokens.success) : null,
      body: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        Kicker(l10n.taxonomySourceLabel, color: t.accentText),
        Padding(
          padding: const EdgeInsets.only(top: 8),
          // The three sources share the row equally, with the same height.
          child: IntrinsicHeight(
            child: Row(
              crossAxisAlignment: CrossAxisAlignment.stretch,
              children: [
                Expanded(child: source(TaxonomySource.gbif, FA.globe, l10n.sourceGbifBody, l10n.badgeOnline)),
                const SizedBox(width: 12),
                Expanded(child: source(TaxonomySource.ncbi, FA.dna, l10n.sourceNcbiBody('$withTaxId', '${c.records.length}'), l10n.badgeOnline)),
                const SizedBox(width: 12),
                Expanded(child: source(TaxonomySource.salve, FA.leaf, l10n.sourceSalveBody, l10n.badgeBrazil)),
              ],
            ),
          ),
        ),
        if (c.taxonomySource == TaxonomySource.salve)
          Padding(
            padding: const EdgeInsets.only(top: 10),
            // SALVE is the national assessment: Brazilian fauna only.
            child: TonePanel(
              tone: PanelTone.warn,
              child: Row(
                children: [
                  const Icon(FA.flag, size: 14),
                  const SizedBox(width: 10),
                  Expanded(child: MarkupText(l10n.salveHint, style: const TextStyle(fontSize: 13, height: 1.45))),
                ],
              ),
            ),
          ),
      ]),
      bodyInset: false,
      footer: [
        Icon(FA.globe, size: 13, color: t.textFaint),
        const SizedBox(width: 8),
        Expanded(
          child: Text(
            l10n.internetRequired,
            maxLines: 2,
            overflow: TextOverflow.ellipsis,
            style: TextStyle(fontSize: 12, color: t.textMuted),
          ),
        ),
        const SizedBox(width: 12),
        SplaceTooltip(
          message: '${l10n.skipTaxBody1} ${l10n.skipTaxBody2}',
          child: OutlinedButton(
            onPressed: c.busy || c.taxonomyResolved ? null : c.skipTaxonomy,
            child: Text(l10n.skipTaxonomy),
          ),
        ),
        const SizedBox(width: 10),
        PrimaryButton(
          label: l10n.fetchTaxonomyFrom(taxonomySourceName(c.taxonomySource)),
          icon: FA.sitemap,
          onPressed: canFetch ? c.fetchTaxonomy : null,
        ),
      ],
    );
  }
}

/// Selectable source card (Material 3 outlined card, filled when chosen).
class _SourceTile extends StatelessWidget {
  const _SourceTile({required this.selected, required this.icon, required this.title, required this.body, required this.badge, required this.onTap});
  final bool selected;
  final IconData icon;
  final String title, body, badge;
  final VoidCallback onTap;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    return Material(
      color: selected ? t.accentSoft : t.surface,
      shape: RoundedRectangleBorder(
        borderRadius: BorderRadius.circular(14),
        side: BorderSide(color: selected ? t.accent : t.border, width: selected ? 2 : 1),
      ),
      child: InkWell(
        borderRadius: BorderRadius.circular(14),
        onTap: onTap,
        child: Padding(
          padding: const EdgeInsets.all(12),
          child: Row(
            crossAxisAlignment: CrossAxisAlignment.start,
            children: [
              Container(
                width: 32,
                height: 32,
                decoration: BoxDecoration(color: selected ? t.accent : t.surface2, borderRadius: BorderRadius.circular(9)),
                child: Icon(icon, size: 14, color: selected ? Colors.white : t.accent),
              ),
              const SizedBox(width: 10),
              Expanded(
                child: Column(
                  crossAxisAlignment: CrossAxisAlignment.start,
                  children: [
                    Row(
                      children: [
                        Flexible(
                          child: Text(
                            title,
                            style: TextStyle(fontWeight: FontWeight.w600, fontSize: 14, color: selected ? t.accentText : t.text),
                          ),
                        ),
                        const SizedBox(width: 6),
                        TagChip(badge, color: selected ? t.accent : t.textMuted),
                      ],
                    ),
                    const SizedBox(height: 3),
                    Text(body, style: TextStyle(fontSize: 12, color: t.textMuted, height: 1.4)),
                  ],
                ),
              ),
              const SizedBox(width: 6),
              Icon(selected ? FA.circleCheck : FA.circleInfo, size: 15, color: selected ? t.accent : Colors.transparent),
            ],
          ),
        ),
      ),
    );
  }
}


/// `#taxonomyTreePanel`: Kingdom → … → Species with record counts and
/// marker chips.
class _TaxonomyTree extends StatefulWidget {
  const _TaxonomyTree({required this.root, required this.records});
  final TaxonNode root;
  final List<GbRecord> records;

  @override
  State<_TaxonomyTree> createState() => _TaxonomyTreeState();
}

/// Collapsed by default; the header opens it.
class _TaxonomyTreeState extends State<_TaxonomyTree> {
  bool _open = false;

  @override
  Widget build(BuildContext context) {
    final root = widget.root;
    final records = widget.records;
    final t = context.tokens;
    final dark = Theme.of(context).brightness == Brightness.dark;
    return Container(
      padding: const EdgeInsets.all(16),
      decoration: BoxDecoration(
        borderRadius: BorderRadius.circular(16),
        border: Border.all(color: t.accentSoftBorder),
        gradient: dark ? null : const LinearGradient(colors: [Color(0xFFEEF0FB), Colors.white, Color(0xFFFEF2F2)], stops: [0, .52, 1]),
        color: dark ? t.surface2 : null,
      ),
      child: Column(
        crossAxisAlignment: CrossAxisAlignment.stretch,
        children: [
          InkWell(
            borderRadius: BorderRadius.circular(8),
            onTap: () => setState(() => _open = !_open),
            child: Row(
              children: [
                Icon(FA.sitemap, size: 14, color: t.accent),
                const SizedBox(width: 8),
                Expanded(
                  child: Text.rich(TextSpan(children: [
                    TextSpan(text: context.l10n.taxonomyTree, style: TextStyle(fontWeight: FontWeight.w600, color: t.text, fontSize: 14)),
                    TextSpan(text: '   ${context.l10n.taxonomyTreeBody}', style: TextStyle(fontSize: 13, color: t.textMuted)),
                  ])),
                ),
                AnimatedRotation(
                  turns: _open ? 0.5 : 0,
                  duration: const Duration(milliseconds: 200),
                  child: Icon(FA.chevronDown, size: 13, color: t.accent),
                ),
              ],
            ),
          ),
          if (_open) ...[
          const SizedBox(height: 10),
          Container(
            constraints: const BoxConstraints(maxHeight: 300),
            decoration: BoxDecoration(
              color: t.surface.withValues(alpha: 0.8),
              borderRadius: BorderRadius.circular(12),
              border: Border.all(color: t.accentSoftBorder),
            ),
            child: SingleChildScrollView(
              padding: const EdgeInsets.all(8),
              child: Column(
                crossAxisAlignment: CrossAxisAlignment.stretch,
                children: [for (final n in root.sortedChildren) _TaxonTile(node: n, records: records, depth: 0)],
              ),
            ),
          ),
          ],
        ],
      ),
    );
  }
}

class _TaxonTile extends StatefulWidget {
  const _TaxonTile({required this.node, required this.records, required this.depth});
  final TaxonNode node;
  final List<GbRecord> records;
  final int depth;

  @override
  State<_TaxonTile> createState() => _TaxonTileState();
}

class _TaxonTileState extends State<_TaxonTile> {
  bool _open = false;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final n = widget.node;
    final italic = n.rank == 'Genus' || n.rank == 'Species';
    final recs = [for (final i in n.records) widget.records[i]];
    String mean(int Function(GbRecord) f) => recs.isEmpty ? '0' : (recs.map(f).reduce((a, b) => a + b) / recs.length).round().toString();
    final leaf = n.children.isEmpty;
    return Column(
      crossAxisAlignment: CrossAxisAlignment.stretch,
      children: [
        InkWell(
          onTap: leaf ? null : () => setState(() => _open = !_open),
          borderRadius: BorderRadius.circular(8),
          child: Padding(
            padding: EdgeInsets.fromLTRB(8.0 + widget.depth * 18, 5, 8, 5),
            child: Row(
              children: [
                SizedBox(width: 14, child: leaf ? null : Icon(_open ? FA.chevronDown : FA.chevronRight, size: 10, color: t.textMuted)),
                const SizedBox(width: 6),
                SizedBox(
                  width: 64,
                  child: Text(
                    _rankLabel(context, n.rank).toUpperCase(),
                    maxLines: 1,
                    overflow: TextOverflow.ellipsis,
                    style: TextStyle(fontSize: 11, fontWeight: FontWeight.w700, color: t.textFaint, letterSpacing: 0.6),
                  ),
                ),
                Flexible(
                  child: Text(
                    n.name,
                    overflow: TextOverflow.ellipsis,
                    style: TextStyle(
                      fontWeight: FontWeight.w600,
                      fontStyle: italic ? FontStyle.italic : FontStyle.normal,
                      color: t.text,
                      fontSize: 13,
                    ),
                  ),
                ),
                const SizedBox(width: 6),
                Text('(${n.records.length})', style: TextStyle(fontSize: 12, color: t.textFaint)),
                const Spacer(),
                TagChip('PCGs ~${mean((r) => r.markerLengths.pcgs)} bp'),
                const SizedBox(width: 4),
                TagChip('rRNAs ~${mean((r) => r.markerLengths.rrnas)} bp'),
                const SizedBox(width: 4),
                TagChip('tRNAs ~${mean((r) => r.markerLengths.trnas)} bp'),
              ],
            ),
          ),
        ),
        if (_open)
          for (final c in n.sortedChildren) _TaxonTile(node: c, records: widget.records, depth: widget.depth + 1),
      ],
    );
  }
}

String _rankLabel(BuildContext context, String rank) {
  final l10n = context.l10n;
  return switch (rank) {
    'Kingdom' => l10n.colKingdom,
    'Phylum' => l10n.colPhylum,
    'Class' => l10n.colClass,
    'Order' => l10n.colOrder,
    'Family' => l10n.colFamily,
    'Genus' => l10n.colGenus,
    _ => l10n.colSpecies,
  };
}

/// "Edit Record" modal: species name, GBIF lookup and taxonomy ranks.
Future<void> showEditRecordDialog(BuildContext context, GbRecord r) => showDialog(
  context: context,
  builder: (_) => ProjectScope(
    controller: ProjectScope.read(context),
    child: _EditRecordDialog(record: r),
  ),
);

class _EditRecordDialog extends StatefulWidget {
  const _EditRecordDialog({required this.record});
  final GbRecord record;

  @override
  State<_EditRecordDialog> createState() => _EditRecordDialogState();
}

class _EditRecordDialogState extends State<_EditRecordDialog> {
  late final GbTaxonomy _tx = widget.record.taxonomy;
  late final _species = TextEditingController(text: widget.record.organismName);
  late final _kingdom = TextEditingController(text: _tx.kingdom);
  late final _phylum = TextEditingController(text: _tx.phylum);
  late final _class = TextEditingController(text: _tx.klass);
  late final _order = TextEditingController(text: _tx.order);
  late final _family = TextEditingController(text: _tx.family);
  late String _authorship = _tx.authorship;
  bool _fetching = false;
  String? _gbifStatus;

  @override
  void dispose() {
    for (final c in [_species, _kingdom, _phylum, _class, _order, _family]) {
      c.dispose();
    }
    super.dispose();
  }

  String? _speciesError(BuildContext context) {
    final v = _species.text.trim();
    final l10n = context.l10n;
    if (v.isEmpty) return l10n.speciesRequired;
    if (!RegExp(r'^[A-Za-z\s.\-]+$').hasMatch(v)) return l10n.speciesLetters;
    final parts = v.split(RegExp(r'\s+'));
    if (parts.length < 2) return l10n.speciesBinomial;
    if (parts.first[0] != parts.first[0].toUpperCase()) return l10n.speciesGenusUpper;
    return null;
  }

  String? _familyError(BuildContext context) {
    final v = _family.text.trim();
    if (v.isEmpty || RegExp(r'(idae|aceae)$', caseSensitive: false).hasMatch(v)) return null;
    return context.l10n.familySuffix;
  }

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    final r = widget.record;
    final fc = r.featureCounts;
    final speciesError = _speciesError(context);
    final familyError = _familyError(context);
    final genus = _species.text.trim().split(RegExp(r'\s+')).first;

    Widget field(String label, TextEditingController ctl, {String? error}) => Field(
      label: label,
      child: TextField(
        controller: ctl,
        decoration: InputDecoration(errorText: error),
        onChanged: (_) => setState(() {}),
      ),
    );

    return SplaceDialog(
      icon: FA.penToSquare,
      title: l10n.editTitle,
      subtitle: l10n.editHint,
      actions: [
        DialogCancelButton(onPressed: () => Navigator.pop(context)),
        FilledButton(
          onPressed: speciesError != null || familyError != null
              ? null
              : () {
                  c.editRecord(
                    r,
                    organism: _species.text.trim(),
                    ranks: TaxonomyRanks(
                      kingdom: _kingdom.text.trim(),
                      phylum: _phylum.text.trim(),
                      klass: _class.text.trim(),
                      order: _order.text.trim(),
                      family: _family.text.trim(),
                      authorship: _authorship,
                    ),
                  );
                  Navigator.pop(context);
                },
          child: Text(l10n.save),
        ),
      ],
      child: Column(
        crossAxisAlignment: CrossAxisAlignment.stretch,
        children: [
          Container(
            padding: const EdgeInsets.all(12),
            decoration: BoxDecoration(
              color: t.surface2,
              borderRadius: BorderRadius.circular(10),
              border: Border.all(color: t.border),
            ),
            child: Column(
              crossAxisAlignment: CrossAxisAlignment.start,
              children: [
                Row(
                  children: [
                    Expanded(child: InfoLine(l10n.colAccession, r.accession)),
                    if (r.source == 'ncbi')
                      TextButton(
                        onPressed: () => launchUrl(Uri.parse('https://www.ncbi.nlm.nih.gov/nuccore/${r.accession}')),
                        child: const Text('NCBI'),
                      ),
                  ],
                ),
                InfoLine(l10n.colSource, r.source == 'ncbi' ? 'NCBI' : l10n.sourceFile),
                if (r.fileName != null) InfoLine(l10n.colFile, r.fileName!),
                if (r.fileSize != null) InfoLine(l10n.fileSize, '${(r.fileSize! / 1024).toStringAsFixed(1)} KB'),
                InfoLine(l10n.sequenceLabel, '${r.sequence.length} bp'),
                InfoLine(l10n.featuresLabel, '${fc.pcgs} PCGs, ${fc.rrnas} rRNAs, ${fc.trnas} tRNAs'),
                if (r.lineage.isNotEmpty) InfoLine(l10n.lineageLabel, r.lineage),
              ],
            ),
          ),
          const SizedBox(height: 14),
          Row(
            crossAxisAlignment: CrossAxisAlignment.start,
            children: [
              Expanded(flex: 3, child: field(l10n.colSpecies, _species, error: speciesError)),
              const SizedBox(width: 12),
              Expanded(
                flex: 2,
                child: Field(
                  label: '${l10n.colGenus} ${l10n.autoLabel}',
                  child: InputDecorator(
                    decoration: const InputDecoration(),
                    child: Text(genus, style: const TextStyle(fontStyle: FontStyle.italic)),
                  ),
                ),
              ),
            ],
          ),
          const SizedBox(height: 10),
          Row(
            children: [
              ChipButton(
                label: l10n.fetchTaxonomyFrom(taxonomySourceName(c.taxonomySource)),
                icon: _fetching ? FA.spinner : FA.globe,
                style: ChipButtonStyle.soft,
                onPressed: speciesError != null || _fetching
                    ? null
                    : () async {
                        setState(() => _fetching = true);
                        final tax = await c.lookupTaxonomy(r, _species.text.trim());
                        if (!mounted) return;
                        setState(() {
                          _fetching = false;
                          if (tax == null) {
                            _gbifStatus = l10n.taxonomyNoMatch(taxonomySourceName(c.taxonomySource));
                          } else {
                            _kingdom.text = tax.kingdom;
                            _phylum.text = tax.phylum;
                            _class.text = tax.klass;
                            _order.text = tax.order;
                            _family.text = tax.family;
                            _authorship = tax.authorship;
                            _gbifStatus = '${l10n.taxonomyMatched(tax.source)} ${tax.authorship}';
                          }
                        });
                      },
              ),
              const SizedBox(width: 10),
              if (_gbifStatus != null)
                Expanded(
                  child: Text(_gbifStatus!, style: TextStyle(fontSize: 12, color: t.textMuted)),
                ),
            ],
          ),
          const SizedBox(height: 12),
          AutoGrid(
            minItemWidth: 200,
            children: [
              field(l10n.colKingdom, _kingdom),
              field(l10n.colPhylum, _phylum),
              field(l10n.colClass, _class),
              field(l10n.colOrder, _order),
              field(l10n.colFamily, _family, error: familyError),
            ],
          ),
        ],
      ),
    );
  }
}


/// Lists the genomes without any usable PCG or rRNA, with the reason for
/// each, before removing them.
Future<void> _showNoMarkersDialog(BuildContext context) => showDialog(
      context: context,
      builder: (_) => ProjectScope(controller: ProjectScope.read(context), child: const _NoMarkersDialog()),
    );

class _NoMarkersDialog extends StatelessWidget {
  const _NoMarkersDialog();

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    final gone = c.recordsWithoutMarkers;
    String reason(({GbRecord record, int pcgs, int rrnas, int trnas}) g) => g.pcgs == 0 && g.rrnas == 0
        ? (g.trnas == 0 ? l10n.noMarkersReasonNone : l10n.noMarkersReasonTrna('${g.trnas}'))
        : l10n.noMarkersReasonUnnamed('${g.pcgs}', '${g.rrnas}');
    return SplaceDialog(
      icon: FA.trashCan,
      title: l10n.noMarkersTitle,
      subtitle: l10n.noMarkersSubtitle,
      scroll: false,
      actions: [
        DialogCancelButton(onPressed: () => Navigator.pop(context)),
        FilledButton(
          style: DialogCancelButton.style(context),
          onPressed: gone.isEmpty
              ? null
              : () {
                  c.removeRecords(gone.map((g) => g.record));
                  Navigator.pop(context);
                },
          child: Text(l10n.noMarkersRemove('${gone.length}')),
        ),
      ],
      child: RecordsTable(
        fill: true,
        columns: [l10n.colSpecies, l10n.colAccession, 'PCGs', 'rRNAs', 'tRNAs', l10n.colReason],
        flex: const [4, 3, 1, 1, 1, 9],
        rowCount: gone.length,
        rowHeight: 48,
        cellBuilder: (context, i, col) {
          final g = gone[i];
          return switch (col) {
            0 => Text(g.record.organismName, style: TextStyle(fontStyle: FontStyle.italic, color: t.text)),
            1 => Text(g.record.accession),
            2 => Text('${g.pcgs}'),
            3 => Text('${g.rrnas}'),
            4 => Text('${g.trnas}'),
            _ => Text(reason(g), maxLines: 2, style: TextStyle(fontSize: 12.5, color: t.textMuted)),
          };
        },
      ),
    );
  }
}
