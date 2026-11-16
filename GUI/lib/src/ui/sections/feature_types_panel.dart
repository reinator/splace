import 'dart:math' as math;

import 'package:flutter/material.dart';

import '../../core/gene_suggestions.dart';
import '../../core/genbank_markers.dart';
import '../../core/syngenes.dart';
import '../../core/syngenes_data.dart';
import '../common.dart';
import '../fa.dart';
import '../theme.dart';
import '../widgets/dialogs.dart';
import '../widgets/kit.dart';
import '../widgets/splace_tooltip.dart';

const synGenesCitation =
    'Rabelo et al. (2024). SynGenes: a Python class for standardizing nomenclatures of mitochondrial and chloroplast genes and a web form for enhancing searches for evolutionary analyses. BMC Bioinformatics 25:160.';

/// Step 3 for GenBank data (`#featureTypesSection`): CDS / rRNA / tRNA
/// switches and the SynGenes note, plus the unknown-gene review.
class FeatureTypesPanel extends StatefulWidget {
  const FeatureTypesPanel({super.key});

  @override
  State<FeatureTypesPanel> createState() => _FeatureTypesPanelState();
}

class _FeatureTypesPanelState extends State<FeatureTypesPanel> {
  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    final cat = c.catalog;
    if (cat == null) return const SizedBox.shrink();

    // The gene review itself runs in step 2, before the taxonomy; here it
    // can be reopened.
    final pending = cat.unknown.length + (c.geneReview?.conflicts.length ?? 0);

    String info(String type) => switch (type) {
          'CDS' => l10n.featureCdsInfo,
          'rRNA' => l10n.featureRrnaInfo,
          _ => l10n.featureTrnaInfo,
        };

    return Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
      TonePanel(
        tone: PanelTone.info,
        child: Row(children: [
          const Icon(FA.wandMagicSparkles, size: 13),
          const SizedBox(width: 8),
          Expanded(
            child: Text.rich(TextSpan(children: [
              TextSpan(text: '${l10n.standardizedPrefix} '),
              WidgetSpan(
                alignment: PlaceholderAlignment.baseline,
                baseline: TextBaseline.alphabetic,
                child: SplaceTooltip(
                  message: synGenesCitation,
                  // A link: same weight as the sentence, underlined.
                  child: Text('SynGenes (Rabelo et al. 2024)', style: TextStyle(decoration: TextDecoration.underline, color: t.accentText)),
                ),
              ),
              TextSpan(text: ' · ${cat.organelle == OrganelleType.cp ? l10n.organelleChloro : l10n.organelleMito}'),
            ])),
          ),
          // Names still to review: the review button. Reviewed (or all
          // recognized): a check, and the review can be reopened.
          if (pending > 0 && !c.unknownGenesReviewed)
            ChipButton(
              label: l10n.geneReviewButtonN('$pending'),
              style: ChipButtonStyle.danger,
              onPressed: () => showGeneReviewDialog(context),
            )
          else ...[
            TagChip(pending > 0 ? l10n.geneReviewDone : l10n.geneNamesAllKnown, icon: FA.circleCheck, color: SplaceTokens.success),
            if (pending > 0) ...[
              const SizedBox(width: 8),
              ChipButton(
                label: l10n.geneReviewAgain('$pending'),
                icon: FA.penToSquare,
                style: ChipButtonStyle.soft,
                onPressed: () => showGeneReviewDialog(context),
              ),
            ],
          ],
        ]),
      ),
      const SizedBox(height: 14),
      Wrap(spacing: 12, runSpacing: 12, children: [
        for (final type in markerFeatureTypes)
          if ((cat.featureTypeCounts[type] ?? 0) > 0)
            SplaceTooltip(
              message: info(type),
              child: GestureDetector(
                onTap: () => c.toggleFeatureType(type, !c.featureTypes.contains(type)),
                child: MouseRegion(
                  cursor: SystemMouseCursors.click,
                  child: Container(
                    padding: const EdgeInsets.fromLTRB(6, 4, 14, 4),
                    decoration: BoxDecoration(
                      color: c.featureTypes.contains(type) ? t.accentSoft : t.surface,
                      borderRadius: BorderRadius.circular(10),
                      border: Border.all(color: c.featureTypes.contains(type) ? t.accentSoftBorder : t.border),
                    ),
                    child: Row(mainAxisSize: MainAxisSize.min, children: [
                      Switch(value: c.featureTypes.contains(type), onChanged: (v) => c.toggleFeatureType(type, v)),
                      const SizedBox(width: 4),
                      Text(featureTypeLabels[type]!,
                          style: TextStyle(fontSize: 14, fontWeight: FontWeight.w600, color: c.featureTypes.contains(type) ? t.text : t.textMuted)),
                      const SizedBox(width: 6),
                      Text('(${cat.featureTypeCounts[type]})', style: TextStyle(fontSize: 13, color: t.textFaint)),
                    ]),
                  ),
                ),
              ),
            ),
      ]),
      if (c.featureTypes.isEmpty) ...[
        const SizedBox(height: 10),
        Text(l10n.featureTypesHint, style: TextStyle(fontSize: 13, color: t.textFaint, fontStyle: FontStyle.italic)),
      ],
    ]);
  }
}

/// Gene review, shown in step 2 before the taxonomy: names SynGenes did not
/// recognize (with the evidence behind each suggestion) and features whose
/// /gene and /product, both normalized with SynGenes, name different
/// genes. Applies global renames for unknown names and per-feature fixes
/// for the conflicts.
Future<void> showGeneReviewDialog(BuildContext context) => showDialog(
      context: context,
      barrierDismissible: false,
      builder: (_) => ProjectScope(controller: ProjectScope.read(context), child: const _GeneReviewDialog()),
    );

String _pct(double v) => '${(v * 100).round()}';

/// Strong / moderate / weak similarity colours (see grHowSequence).
Color _simColor(BuildContext context, double v) => v >= 0.5 ? context.tokens.accent : (v >= 0.35 ? SplaceTokens.warn : context.tokens.textFaint);

/// One "gene ███ 72%" line.
class _SimilarityBar extends StatelessWidget {
  const _SimilarityBar(this.label, this.value);
  final String label;
  final double value;
  static const width = 70.0;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final color = _simColor(context, value);
    return Padding(
      padding: const EdgeInsets.symmetric(vertical: 2),
      child: Row(children: [
        SizedBox(width: width, child: Text(label, maxLines: 1, overflow: TextOverflow.ellipsis, style: TextStyle(fontSize: 12, fontWeight: FontWeight.w600, color: t.text))),
        Expanded(
          child: ClipRRect(
            borderRadius: BorderRadius.circular(4),
            child: Container(
              height: 8,
              color: t.surface3,
              alignment: Alignment.centerLeft,
              child: FractionallySizedBox(widthFactor: value.clamp(0.0, 1.0), heightFactor: 1, child: ColoredBox(color: color)),
            ),
          ),
        ),
        SizedBox(width: 44, child: Text('${_pct(value)}%', textAlign: TextAlign.right, style: TextStyle(fontSize: 12, fontWeight: FontWeight.w700, color: color))),
      ]),
    );
  }
}

class _GeneReviewDialog extends StatefulWidget {
  const _GeneReviewDialog();

  @override
  State<_GeneReviewDialog> createState() => _GeneReviewDialogState();
}

class _GeneReviewDialogState extends State<_GeneReviewDialog> {
  /// Unknown name → standard name ('' = skip).
  final _choice = <String, String>{};

  /// featureKey → chosen name for each conflict.
  final _conflictChoice = <String, String>{};

  /// Which tag the sequence supports: 1 /gene, -1 /product, 0 unclear.
  static int _favours(QualifierConflict k) {
    final g = k.geneScore ?? 0, p = k.productScore ?? 0;
    if (math.max(g, p) < 0.35 || (g - p).abs() < 0.1) return 0;
    return g > p ? 1 : -1;
  }

  @override
  void initState() {
    super.initState();
    final c = ProjectScope.read(context);
    final review = c.geneReview;
    // Start from the suggestions; the user can still change or skip each.
    for (final s in review?.suggestions.values ?? const <GeneSuggestion>[]) {
      if (s.suggested != null) _choice[s.name] = s.suggested!;
    }
    // Conflicts start on what the sequence supports, else the name in use.
    for (final k in review?.conflicts ?? const <QualifierConflict>[]) {
      final f = _favours(k);
      _conflictChoice[featureKey(c.records[k.record], k.feature)] = f == 1 ? k.geneName : (f == -1 ? k.productName : k.used);
    }
  }

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    final unknown = c.catalog?.unknown ?? const <UnknownGene>[];
    final review = c.geneReview;
    final conflicts = review?.conflicts ?? const <QualifierConflict>[];
    final primary = c.organelle == OrganelleType.cp ? chloroplastGenes : mitochondrialGenes;
    final secondary = c.organelle == OrganelleType.cp ? mitochondrialGenes : chloroplastGenes;
    final options = [
      '',
      ...primary.keys.toList()..sort(),
      ...secondary.keys.where((k) => !primary.containsKey(k)).toList()..sort(),
    ];

    Widget section(String title) => Padding(
          padding: const EdgeInsets.only(bottom: 8, top: 10),
          child: Kicker(title, size: 12, color: t.accentText),
        );

    BoxDecoration card() => BoxDecoration(color: t.surface2, borderRadius: BorderRadius.circular(12), border: Border.all(color: t.border));

    // Evidence rows share one label column so values line up.
    const labelW = 92.0;
    Widget evidenceLabel(String text) => SizedBox(
          width: labelW,
          child: Text(text, style: TextStyle(fontSize: 12.5, fontWeight: FontWeight.w700, color: t.textMuted)),
        );

    // "/gene   “orf220”   ✗ not in SynGenes"
    Widget tagRow(String tag, Set<String> values, String? standard) => Padding(
          padding: const EdgeInsets.symmetric(vertical: 4),
          child: Row(crossAxisAlignment: CrossAxisAlignment.center, children: [
            evidenceLabel(tag),
            Expanded(
              child: Text(values.isEmpty ? l10n.grAbsent : values.map((v) => '“$v”').join(', '),
                  maxLines: 2, overflow: TextOverflow.ellipsis, style: TextStyle(fontSize: 13, color: values.isEmpty ? t.textFaint : t.text)),
            ),
            const SizedBox(width: 10),
            if (values.isNotEmpty)
              standard != null ? TagChip(l10n.grRecognized(standard), icon: FA.check) : TagChip(l10n.grNotRecognized, icon: FA.xmark, color: SplaceTokens.warn),
          ]),
        );

    Widget sequenceRow(List<(String, double)> ranking) => Padding(
          padding: const EdgeInsets.only(top: 6),
          child: Row(crossAxisAlignment: CrossAxisAlignment.start, children: [
            Padding(padding: const EdgeInsets.only(top: 2), child: evidenceLabel(l10n.grSequence)),
            Expanded(
              child: ranking.isEmpty
                  ? Text(l10n.grNoSequenceRef, style: TextStyle(fontSize: 12.5, color: t.textFaint))
                  : Column(children: [for (final (g, v) in ranking) _SimilarityBar(g, v)]),
            ),
          ]),
        );

    // Header, evidence (left) and decision (right, same height).
    Widget reviewCard({required Widget header, required List<Widget> evidence, required Widget decision}) => Container(
          margin: const EdgeInsets.only(bottom: 12),
          decoration: card(),
          clipBehavior: Clip.antiAlias,
          child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
            Container(
              padding: const EdgeInsets.fromLTRB(16, 12, 16, 12),
              decoration: BoxDecoration(color: t.surface, border: Border(bottom: BorderSide(color: t.border))),
              child: header,
            ),
            IntrinsicHeight(
              child: Row(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
                Expanded(
                  flex: 6,
                  child: Padding(
                    padding: const EdgeInsets.fromLTRB(16, 12, 16, 14),
                    child: Column(mainAxisAlignment: MainAxisAlignment.center, crossAxisAlignment: CrossAxisAlignment.stretch, children: evidence),
                  ),
                ),
                Container(width: 1, color: t.border),
                Expanded(
                  flex: 4,
                  child: Container(
                    color: t.surface,
                    padding: const EdgeInsets.all(16),
                    child: Column(mainAxisAlignment: MainAxisAlignment.center, crossAxisAlignment: CrossAxisAlignment.stretch, children: [decision]),
                  ),
                ),
              ]),
            ),
          ]),
        );

    Widget unknownCard(UnknownGene u) {
      final s = review?.suggestions[u.name];
      final m = s?.match;
      final (why, whyColor) = switch (s?.source) {
        SuggestionSource.product => (l10n.grWhyProduct(s!.products.first, s.suggested!), t.accentText),
        SuggestionSource.gene => (l10n.grWhyGene(s!.genes.first, s.suggested!), t.accentText),
        SuggestionSource.sequence => (l10n.grWhySequence(_pct(m!.score), m.gene, '${m.references}'), t.accentText),
        _ => (l10n.grWhyNone, SplaceTokens.warn),
      };
      final seqNote = s == null || s.source == SuggestionSource.sequence || m == null
          ? null
          : (m.gene == s.suggested ? l10n.grSeqAgrees(_pct(m.score)) : l10n.grSeqDisagrees(m.gene, _pct(m.score)));
      return reviewCard(
        header: Row(children: [
          Text(u.name, style: TextStyle(fontSize: 15, fontWeight: FontWeight.w700, color: t.text)),
          const SizedBox(width: 8),
          TagChip(featureTypeLabels[u.type] ?? u.type),
          const SizedBox(width: 10),
          Expanded(
            child: Text(
              '${l10n.grInGenomes('${u.records.length}')} · ${[
                for (final i in u.records.take(3)) '${c.records[i].organismName} (${c.records[i].accession})',
                if (u.records.length > 3) '+${u.records.length - 3} ${l10n.more}',
              ].join(', ')}',
              maxLines: 1,
              overflow: TextOverflow.ellipsis,
              style: TextStyle(fontSize: 12.5, color: t.textMuted, fontStyle: FontStyle.italic),
            ),
          ),
        ]),
        evidence: [
          tagRow('/gene', s?.genes ?? const {}, s?.byGene),
          tagRow('/product', s?.products ?? const {}, s?.byProduct),
          sequenceRow(s?.ranking ?? const []),
        ],
        decision: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
          Kicker(l10n.grUseName, size: 11),
          const SizedBox(height: 6),
          DropdownButtonFormField<String>(
            icon: dropdownChevron,
            initialValue: _choice[u.name] ?? '',
            isExpanded: true,
            menuMaxHeight: 360,
            items: [for (final o in options) DropdownMenuItem(value: o, child: Text(o.isEmpty ? l10n.unknownSkip : o))],
            onChanged: (v) => setState(() => _choice[u.name] = v ?? ''),
          ),
          const SizedBox(height: 10),
          Text(why, style: TextStyle(fontSize: 12.5, color: whyColor, height: 1.45)),
          if (seqNote != null) ...[
            const SizedBox(height: 4),
            Text(seqNote, style: TextStyle(fontSize: 12.5, height: 1.45, color: m!.gene == s!.suggested ? t.accentText : SplaceTokens.warn)),
          ],
        ]),
      );
    }

    Widget conflictCard(QualifierConflict k) {
      final r = c.records[k.record];
      final key = featureKey(r, k.feature);
      final fav = _favours(k);
      final chosen = _conflictChoice[key] ?? k.used;
      final explanation = switch (fav) {
        1 => l10n.grConflictGeneRight(_pct(k.geneScore!), k.geneName, _pct(k.productScore ?? 0), k.productName),
        -1 => l10n.grConflictProductRight(_pct(k.productScore!), k.productName, _pct(k.geneScore ?? 0), k.geneName),
        _ => l10n.grConflictUnclear,
      };
      // Bars: the closest markers plus both names the tags point to.
      final ranking = <(String, double)>[
        for (final e in k.scores.entries.take(3)) (e.key, e.value),
        for (final n in {k.geneName, k.productName})
          if (!k.scores.keys.take(3).contains(n) && k.scores.containsKey(n)) (n, k.scores[n]!),
      ];

      Widget option(String name, String tag, bool supported) {
        final selected = chosen == name;
        return Padding(
          padding: const EdgeInsets.only(bottom: 8),
          child: InkWell(
            borderRadius: BorderRadius.circular(10),
            onTap: () => setState(() => _conflictChoice[key] = name),
            child: Container(
              padding: const EdgeInsets.symmetric(horizontal: 12, vertical: 10),
              decoration: BoxDecoration(
                color: selected ? t.accentSoft : t.surface2,
                borderRadius: BorderRadius.circular(10),
                border: Border.all(color: selected ? t.accent : t.border, width: selected ? 2 : 1),
              ),
              child: Row(children: [
                Icon(selected ? FA.circleCheck : FA.circle, size: 15, color: selected ? t.accent : t.surface3),
                const SizedBox(width: 10),
                Expanded(
                  child: Text(l10n.grSays(tag, name), style: TextStyle(fontSize: 13.5, fontWeight: FontWeight.w700, color: selected ? t.accentText : t.text)),
                ),
                if (fav != 0) TagChip(supported ? l10n.grSupported : l10n.grNotSupported, color: supported ? t.accent : SplaceTokens.warn),
              ]),
            ),
          ),
        );
      }

      return reviewCard(
        header: Row(children: [
          Text(r.organismName, style: TextStyle(fontSize: 15, fontStyle: FontStyle.italic, fontWeight: FontWeight.w700, color: t.text)),
          const SizedBox(width: 8),
          TagChip(featureTypeLabels[k.type] ?? k.type),
          const SizedBox(width: 10),
          Expanded(
            child: Text('${r.accession} · ${r.features[k.feature].location}',
                maxLines: 1, overflow: TextOverflow.ellipsis, style: TextStyle(fontSize: 12.5, color: t.textMuted)),
          ),
        ]),
        evidence: [
          tagRow('/gene', {k.gene}, k.geneName),
          tagRow('/product', {k.product}, k.productName),
          sequenceRow(ranking),
        ],
        decision: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
          Kicker(l10n.grUseName, size: 11),
          const SizedBox(height: 6),
          option(k.geneName, '/gene', fav == 1),
          option(k.productName, '/product', fav == -1),
          const SizedBox(height: 2),
          Text(explanation, style: TextStyle(fontSize: 12.5, height: 1.45, fontWeight: FontWeight.w600, color: fav == 0 ? SplaceTokens.warn : t.accentText)),
        ]),
      );
    }

    return SplaceDialog(
      icon: FA.dna,
      title: l10n.geneReviewTitle,
      subtitle: l10n.geneReviewSubtitle,
      actions: [
        DialogCancelButton(
          label: l10n.unknownSkipAll,
          onPressed: () {
            c.dismissUnknownGenes();
            Navigator.pop(context);
          },
        ),
        FilledButton(
          onPressed: () {
            final perFeature = <String, String>{};
            for (final k in conflicts) {
              final key = featureKey(c.records[k.record], k.feature);
              final chosen = _conflictChoice[key];
              if (chosen != null && chosen != k.used) perFeature[key] = chosen;
            }
            c.applyGeneOverrides(_choice, perFeature: perFeature);
            Navigator.pop(context);
          },
          child: Text(l10n.unknownApply),
        ),
      ],
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        if (unknown.isEmpty && conflicts.isEmpty)
          TonePanel(tone: PanelTone.info, child: Text(l10n.geneReviewNothing))
        else
          TonePanel(
            tone: PanelTone.info,
            title: l10n.grHowTitle,
            icon: FA.circleInfo,
            child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
              Text(l10n.grHowSynGenes, style: const TextStyle(fontSize: 13)),
              const SizedBox(height: 4),
              Text(l10n.grHowSequence, style: const TextStyle(fontSize: 13)),
            ]),
          ),
        if (conflicts.isNotEmpty) ...[
          section(l10n.geneReviewConflictSection('${conflicts.length}')),
          for (final k in conflicts) conflictCard(k),
        ],
        if (unknown.isNotEmpty) ...[
          section(l10n.geneReviewUnknownSection('${unknown.length}')),
          for (final u in unknown) unknownCard(u),
        ],
      ]),
    );
  }
}
