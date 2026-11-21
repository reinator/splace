import 'package:file_selector/file_selector.dart';
import 'package:flutter/material.dart';
import 'package:path/path.dart' as p;

import '../../services/species_id.dart';
import '../../state/project_controller.dart';
import '../common.dart';
import '../fa.dart';
import '../theme.dart';
import '../widgets/kit.dart';
import '../widgets/splace_tooltip.dart';
import '../widgets/toast.dart';
import '../chrome/id_alignment_dialog.dart';
import '../chrome/species_id_params.dart';

/// Step 2 (FASTA): the barcoding switch and, when on, the identification of
/// each specimen (GenBank BLAST; BOLD for COI), the name kept for each one
/// and the BOLD submission spreadsheet (BOLD5 single-page template) with
/// the sequences.
class BarcodingCard extends StatefulWidget {
  const BarcodingCard({super.key});

  @override
  State<BarcodingCard> createState() => _BarcodingCardState();
}

class _BarcodingCardState extends State<BarcodingCard> {
  final Map<String, TextEditingController> _names = {};

  @override
  void dispose() {
    for (final c in _names.values) {
      c.dispose();
    }
    super.dispose();
  }

  TextEditingController _field(ProjectController c, String sample) {
    final ctrl = _names.putIfAbsent(sample, () => TextEditingController(text: c.barcodeNames[sample] ?? ''));
    final v = c.barcodeNames[sample] ?? '';
    if (ctrl.text != v) ctrl.text = v;
    return ctrl;
  }

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    final on = c.barcodingMode;
    final seqs = on ? c.barcodeSequences : const <String, String>{};
    final locus = c.barcodeLocusShown;
    final coi = locus != null && c.speciesIdOptions.boldAccepts(ProjectController.barcodeGene(locus));
    final named = c.barcodeNames.values.where((n) => n.trim().isNotEmpty).length;

    return Container(
      padding: const EdgeInsets.all(16),
      decoration: BoxDecoration(color: t.surface, borderRadius: BorderRadius.circular(12), border: Border.all(color: t.border), boxShadow: shadowSm),
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        SwitchBox(label: l10n.barcodingSwitch, detail: l10n.barcodingSwitchInfo, value: on, onChanged: c.busy ? null : c.setBarcodingMode),
        if (on) ...[
          const SizedBox(height: 14),
          Wrap(spacing: 14, runSpacing: 10, crossAxisAlignment: WrapCrossAlignment.end, children: [
            SizedBox(
              width: 220,
              child: Field(
                label: l10n.barcodingLocus,
                info: l10n.barcodingLocusInfo,
                child: SplaceSelect<String>(
                  value: locus ?? '',
                  options: [for (final l in c.barcodeLoci) SelectOption(l, l)],
                  onChanged: c.busy ? null : c.setBarcodeLocus,
                ),
              ),
            ),
            TagChip(l10n.barcodingCount('${seqs.length}'), color: t.textMuted),
            TagChip(coi ? l10n.barcodingSourcesCoi : l10n.barcodingSourcesGenbank, icon: FA.circleInfo, color: t.textMuted),
            // BLAST and BOLD parameters.
            OutlinedButton.icon(
              onPressed: c.busy ? null : () => showSpeciesIdParams(context),
              icon: const Icon(FA.sliders, size: 13),
              label: Text(l10n.idpTitle),
            ),
          ]),
          const SizedBox(height: 14),
          Wrap(spacing: 10, runSpacing: 8, crossAxisAlignment: WrapCrossAlignment.center, children: [
            PrimaryButton(
              label: c.barcodeIds.isEmpty ? l10n.barcodingIdentify : l10n.siCardRerun,
              icon: FA.barcode,
              onPressed: c.busy || seqs.isEmpty ? null : c.identifyBarcodes,
            ),
            SplaceTooltip(
              message: l10n.barcodingTaxonomyInfo,
              child: OutlinedButton(
                onPressed: c.busy || named == 0 ? null : c.fetchBarcodeTaxonomy,
                child: Text(l10n.barcodingTaxonomy),
              ),
            ),
            SplaceTooltip(
              message: l10n.barcodingExportInfo,
              child: OutlinedButton(
                onPressed: c.busy || seqs.isEmpty
                    ? null
                    : () async {
                        final loc = await getSaveLocation(
                          suggestedName: 'BOLD_submission.xlsx',
                          acceptedTypeGroups: const [XTypeGroup(label: 'Excel', extensions: ['xlsx'])],
                        );
                        if (loc == null) return;
                        final n = await c.exportBoldSubmission(loc.path);
                        if (context.mounted) {
                          showSplaceToast(context,
                              title: l10n.barcodingExported('$n'),
                              detail: p.basename(loc.path),
                              path: p.dirname(loc.path),
                              action: (label: l10n.openFolder, onPressed: () => c.openFolder(p.dirname(loc.path))));
                        }
                      },
                child: Text(l10n.barcodingExport),
              ),
            ),
          ]),
          const SizedBox(height: 8),
          Text(l10n.internetRequired, style: TextStyle(fontSize: 12, color: t.textFaint)),
          if (seqs.isNotEmpty) ...[
            const SizedBox(height: 12),
            Row(children: [
              SizedBox(width: 200, child: Kicker(l10n.barcodingColSample, size: 10.5)),
              Expanded(child: Kicker(l10n.barcodingColMatches, size: 10.5)),
              SizedBox(width: 300, child: Kicker(l10n.barcodingColName, size: 10.5)),
            ]),
            const SizedBox(height: 6),
            // Lazy: barcoding sets can have hundreds of specimens.
            SizedBox(
              height: (seqs.length * 48.0).clamp(48.0, 440.0),
              child: ListView.builder(
                itemCount: seqs.length,
                itemExtent: 48,
                itemBuilder: (context, i) {
                  final sample = seqs.keys.elementAt(i);
                  final id = c.barcodeIds[sample];
                  final cands = id?.shortlist(perSource: 3) ?? const <SpeciesCandidate>[];
                  String pct(double v) => '${v.toStringAsFixed(1)}%';
                  final summary = cands.isEmpty
                      ? (id == null ? '—' : l10n.barcodingNoMatch)
                      : cands.take(3).map((k) {
                          final src = [if (k.bold != null) 'BOLD ${pct(k.bold!.identity)}', if (k.inGenbank) 'GB ${pct(k.genbankIdentity)}'];
                          return '${binomialOf(k.species)} (${src.join(', ')})';
                        }).join(' · ');
                  return Row(children: [
                    SizedBox(width: 200, child: Text(sample, overflow: TextOverflow.ellipsis, style: TextStyle(fontFamily: monoFont, fontSize: 12.5, color: t.text))),
                    Expanded(
                      child: Text.rich(italicNames(summary, [for (final k in cands) binomialOf(k.species)], style: TextStyle(fontSize: 12.5, color: t.textMuted)),
                          overflow: TextOverflow.ellipsis, maxLines: 2),
                    ),
                    SizedBox(
                      width: 300,
                      child: Row(children: [
                        Expanded(
                          child: SizedBox(
                            height: 36,
                            child: TextField(
                              controller: _field(c, sample),
                              style: const TextStyle(fontSize: 13, fontStyle: FontStyle.italic),
                              decoration: InputDecoration(isDense: true, hintText: l10n.barcodingNameHint),
                              onChanged: (v) => c.setBarcodeName(sample, v),
                            ),
                          ),
                        ),
                        if (id != null && cands.any((k) => k.hits.isNotEmpty))
                          IconButton(
                            tooltip: l10n.idAlnOpen,
                            visualDensity: VisualDensity.compact,
                            icon: Icon(FA.alignLeft, size: 13, color: t.accent),
                            onPressed: () => showIdAlignments(context, id, name: sample),
                          ),
                        if (cands.isNotEmpty)
                          PopupMenuButton<String>(
                            tooltip: l10n.barcodingPick,
                            icon: Icon(FA.chevronDown, size: 12, color: t.textMuted),
                            onSelected: (v) => c.setBarcodeName(sample, v),
                            itemBuilder: (_) => [
                              for (final k in cands)
                                PopupMenuItem(
                                  value: binomialOf(k.species),
                                  child: Text.rich(TextSpan(children: [
                                    TextSpan(text: binomialOf(k.species), style: const TextStyle(fontStyle: FontStyle.italic)),
                                    TextSpan(text: ' · ${pct(k.bestIdentity)}'),
                                  ])),
                                ),
                            ],
                          ),
                      ]),
                    ),
                  ]);
                },
              ),
            ),
          ],
        ],
      ]),
    );
  }
}
