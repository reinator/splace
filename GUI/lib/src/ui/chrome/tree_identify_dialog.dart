import 'package:flutter/material.dart';

import '../../services/species_id.dart';
import '../../state/project_controller.dart';
import '../common.dart';
import '../fa.dart';
import '../theme.dart';
import '../widgets/dialogs.dart';
import '../widgets/kit.dart';
import 'id_alignment_dialog.dart';

/// "Identify by BLAST" on the tree: the [tips] chosen (a tip or a whole
/// clade, e.g. a genus that falls inside another one) are sent with the
/// sequence of one marker the user picks to GenBank BLAST (and BOLD for its
/// markers), with the options of step 2. Each tip gets the species of its
/// best match and a verdict against its current name. Gives back the tips
/// to rename in the figure (tip → best species), or null. With [results]
/// (identifications already in hand) it opens on them.
Future<Map<String, String>?> showTreeIdentifyDialog(BuildContext context, List<String> tips,
        {String Function(String tip)? label, String? marker, Map<String, RecordIdentification>? results}) =>
    showDialog<Map<String, String>>(
      context: context,
      barrierDismissible: false,
      builder: (_) => _TreeIdentifyDialog(tips: tips, label: label ?? (t) => t.replaceAll('_', ' '), marker: marker, results: results),
    );

/// Verdict of one tip.
enum TipVerdict {
  /// The best match is the tip's species.
  confirmed,

  /// Another species of the same genus at a same-species identity.
  otherSpecies,

  /// A species of another genus at a same-species identity.
  otherGenus,

  /// The tip's species and another one match as well (shared haplotypes).
  ambiguous,

  /// Nothing at a same-species identity: only relatives.
  noMatch,

  /// No hit at all.
  none,
}

/// The verdict of [r] (see [TipVerdict]).
TipVerdict tipVerdict(RecordIdentification? r) {
  final best = r?.candidates.firstOrNull;
  if (r == null || best == null) return TipVerdict.none;
  if (r.isCurrent(best.species)) return TipVerdict.confirmed;
  if (best.bestIdentity < sameSpeciesIdentity) return TipVerdict.noMatch;
  if (!r.disagrees) return TipVerdict.ambiguous;
  String genus(String s) => binomialOf(s).split(' ').first.toLowerCase();
  return genus(best.species) == genus(r.currentName) ? TipVerdict.otherSpecies : TipVerdict.otherGenus;
}

class _TreeIdentifyDialog extends StatefulWidget {
  const _TreeIdentifyDialog({required this.tips, required this.label, this.marker, this.results});
  final List<String> tips;
  final String Function(String tip) label;
  final String? marker;
  final Map<String, RecordIdentification>? results;

  @override
  State<_TreeIdentifyDialog> createState() => _TreeIdentifyDialogState();
}

class _TreeIdentifyDialogState extends State<_TreeIdentifyDialog> {
  late String? _marker = widget.marker;
  bool _running = false;
  double _share = 0;
  String _status = '';
  Map<String, RecordIdentification>? _ids;

  /// Tips to rename in the figure with their best match.
  final Set<String> _rename = {};

  @override
  void initState() {
    super.initState();
    if (widget.results != null) _show(widget.results);
  }

  /// Shows [ids]; the tips named as another species are ticked to rename.
  void _show(Map<String, RecordIdentification>? ids) {
    _ids = ids;
    _rename
      ..clear()
      ..addAll([
        for (final e in (ids ?? const <String, RecordIdentification>{}).entries)
          if (const {TipVerdict.otherSpecies, TipVerdict.otherGenus}.contains(tipVerdict(e.value))) e.key,
      ]);
  }

  Future<void> _run(ProjectController c) async {
    final marker = _marker;
    if (marker == null) return;
    setState(() {
      _running = true;
      _share = 0;
      _status = '';
    });
    final ids = await c.identifyTips(widget.tips, marker, onProgress: (share, status) {
      if (mounted) {
        setState(() {
          _share = share;
          _status = status;
        });
      }
    });
    if (!mounted) return;
    setState(() {
      _running = false;
      _show(ids);
    });
  }

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    final markers = c.tipMarkers(widget.tips);
    _marker ??= markers.firstOrNull?.marker;
    final ids = _ids;
    final o = c.speciesIdOptions;

    Widget body;
    if (ids != null) {
      body = _Results(tips: widget.tips, ids: ids, label: widget.label, rename: _rename, onToggle: (tip) => setState(() => _rename.contains(tip) ? _rename.remove(tip) : _rename.add(tip)));
    } else if (_running) {
      body = Padding(
        padding: const EdgeInsets.symmetric(vertical: 40),
        child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
          Text(l10n.treeIdRunning(_marker ?? '', '${widget.tips.length}'), style: TextStyle(fontSize: 15, fontWeight: FontWeight.w600, color: t.text)),
          const SizedBox(height: 12),
          RainbowProgressBar(value: _share <= 0 ? null : _share),
          const SizedBox(height: 8),
          Text(_status.isEmpty ? ' ' : _status, style: TextStyle(fontSize: 12.5, color: t.textMuted)),
        ]),
      );
    } else {
      body = Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        Text(l10n.treeIdIntro, style: TextStyle(fontSize: 13.5, color: t.textMuted, height: 1.5)),
        const SizedBox(height: 16),
        Kicker(l10n.treeIdMarker, size: 10.5),
        const SizedBox(height: 8),
        if (markers.isEmpty)
          TonePanel(tone: PanelTone.warn, child: Text(l10n.treeIdNoMarker, style: const TextStyle(fontSize: 13)))
        else
          AutoGrid(minItemWidth: 260, children: [
            for (final m in markers)
              RadioCard<String>(
                value: m.marker,
                group: _marker ?? '',
                onChanged: (v) => setState(() => _marker = v),
                title: m.marker,
                subtitle: l10n.treeIdMarkerTips('${m.tips}', '${widget.tips.length}'),
                trailing: o.boldAccepts(ProjectController.barcodeGene(m.marker)) ? TagChip('GenBank + BOLD', color: t.accentText) : TagChip('GenBank', color: t.textMuted),
              ),
          ]),
        const SizedBox(height: 16),
        Text(l10n.treeIdOptions(o.blastProgram, o.blastDatabase, '${(o.minCoverage * 100).round()}', sameSpeciesIdentity.toStringAsFixed(0)),
            style: TextStyle(fontSize: 12.5, color: t.textFaint, height: 1.4)),
      ]);
    }

    final failures = c.speciesIdFailures;
    return SplaceDialog(
      icon: FA.magnifyingGlass,
      title: l10n.treeIdTitle,
      subtitle: l10n.treeIdSubtitle('${widget.tips.length}'),
      footerLeading: ids != null && failures.isNotEmpty
          ? Text(failures.entries.map((e) => '${e.key == IdSource.bold ? 'BOLD' : 'GenBank'}: ${e.value}').join(' · '),
              maxLines: 2, overflow: TextOverflow.ellipsis, style: const TextStyle(fontSize: 12, color: SplaceTokens.danger))
          : null,
      actions: [
        if (_running)
          DialogCancelButton(onPressed: c.cancelTask)
        else ...[
          TextButton(onPressed: () => Navigator.of(context).pop(), child: Text(ids == null ? l10n.cancel : l10n.close)),
          if (ids == null)
            PrimaryButton(label: l10n.treeIdRun, icon: FA.magnifyingGlass, onPressed: _marker == null || c.busy ? null : () => _run(c))
          else ...[
            ChipButton(label: l10n.treeIdOtherMarker, onPressed: () => setState(() => _ids = null)),
            PrimaryButton(
              label: l10n.treeIdRename('${_rename.length}'),
              icon: FA.penToSquare,
              onPressed: _rename.isEmpty
                  ? null
                  : () => Navigator.of(context).pop({
                        for (final tip in _rename)
                          if (ids[tip]?.candidates.firstOrNull case final best?) tip: binomialOf(best.species),
                      }),
            ),
          ],
        ],
      ],
      child: body,
    );
  }
}

class _Results extends StatelessWidget {
  const _Results({required this.tips, required this.ids, required this.label, required this.rename, required this.onToggle});
  final List<String> tips;
  final Map<String, RecordIdentification> ids;
  final String Function(String tip) label;
  final Set<String> rename;
  final ValueChanged<String> onToggle;

  @override
  Widget build(BuildContext context) {
    final l10n = context.l10n;
    final t = context.tokens;
    // Problems first, then the tips that agree, then those not sent.
    const order = [TipVerdict.otherGenus, TipVerdict.otherSpecies, TipVerdict.ambiguous, TipVerdict.noMatch, TipVerdict.confirmed, TipVerdict.none];
    final rows = [...tips]..sort((a, b) {
        final va = order.indexOf(tipVerdict(ids[a])), vb = order.indexOf(tipVerdict(ids[b]));
        return va != vb ? va - vb : a.compareTo(b);
      });
    final counts = <TipVerdict, int>{};
    for (final tip in tips) {
      counts.update(tipVerdict(ids[tip]), (v) => v + 1, ifAbsent: () => 1);
    }
    (String, Color, IconData) look(TipVerdict v) => switch (v) {
          TipVerdict.confirmed => (l10n.treeIdConfirmed, SplaceTokens.success, FA.circleCheck),
          TipVerdict.otherSpecies => (l10n.treeIdOtherSpecies, SplaceTokens.warn, FA.triangleExclamation),
          TipVerdict.otherGenus => (l10n.treeIdOtherGenus, SplaceTokens.danger, FA.triangleExclamation),
          TipVerdict.ambiguous => (l10n.treeIdAmbiguous, SplaceTokens.warn, FA.circleInfo),
          TipVerdict.noMatch => (l10n.treeIdNoMatch(sameSpeciesIdentity.toStringAsFixed(0)), t.textMuted, FA.circleInfo),
          TipVerdict.none => (l10n.treeIdNone, t.textFaint, FA.minus),
        };
    String pct(double v) => '${v.toStringAsFixed(1)} %';

    return Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
      Wrap(spacing: 8, runSpacing: 6, children: [
        for (final v in order)
          if ((counts[v] ?? 0) > 0)
            () {
              final (text, color, icon) = look(v);
              return TagChip('$text · ${counts[v]}', icon: icon, color: color);
            }(),
      ]),
      const SizedBox(height: 14),
      for (final tip in rows)
        () {
          final r = ids[tip];
          final v = tipVerdict(r);
          final (text, color, icon) = look(v);
          final best = r?.candidates.firstOrNull;
          final wrong = v == TipVerdict.otherSpecies || v == TipVerdict.otherGenus;
          return Container(
            margin: const EdgeInsets.only(bottom: 8),
            padding: const EdgeInsets.fromLTRB(12, 10, 10, 10),
            decoration: BoxDecoration(
              color: t.surface,
              borderRadius: BorderRadius.circular(10),
              border: Border.all(color: wrong ? color.withValues(alpha: 0.5) : t.border, width: wrong ? 1.5 : 1),
            ),
            child: Row(children: [
              Icon(icon, size: 14, color: color),
              const SizedBox(width: 10),
              Expanded(
                flex: 5,
                child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
                  Text(label(tip), maxLines: 1, overflow: TextOverflow.ellipsis, style: TextStyle(fontSize: 13, fontWeight: FontWeight.w600, color: t.text)),
                  const SizedBox(height: 2),
                  Text.rich(
                    TextSpan(style: TextStyle(fontSize: 12, color: t.textMuted), children: [
                      TextSpan(text: '${l10n.treeIdCurrent} '),
                      TextSpan(text: r?.currentName ?? '—', style: const TextStyle(fontStyle: FontStyle.italic)),
                    ]),
                    maxLines: 1,
                    overflow: TextOverflow.ellipsis,
                  ),
                ]),
              ),
              const SizedBox(width: 12),
              Expanded(
                flex: 5,
                child: best == null
                    ? Text(text, style: TextStyle(fontSize: 12.5, color: color))
                    : Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
                        Text(binomialOf(best.species),
                            maxLines: 1, overflow: TextOverflow.ellipsis, style: TextStyle(fontSize: 13, fontStyle: FontStyle.italic, fontWeight: FontWeight.w600, color: wrong ? color : t.text)),
                        const SizedBox(height: 2),
                        Text(
                          [
                            text,
                            if (best.inGenbank) 'GenBank ${pct(best.genbankIdentity)}',
                            if (best.bold != null) 'BOLD ${pct(best.bold!.identity)}',
                            if (v != TipVerdict.confirmed && r?.current != null) l10n.treeIdCurrentAt(pct(r!.current!.bestIdentity)),
                          ].join(' · '),
                          maxLines: 2,
                          overflow: TextOverflow.ellipsis,
                          style: TextStyle(fontSize: 12, color: v == TipVerdict.confirmed ? t.textMuted : color),
                        ),
                      ]),
              ),
              const SizedBox(width: 8),
              SizedBox(
                width: 130,
                child: r != null && best != null
                    ? Align(
                        alignment: Alignment.centerRight,
                        child: ChipButton(label: l10n.treeIdHits, onPressed: () => showIdAlignments(context, r, species: best.species, name: label(tip))),
                      )
                    : null,
              ),
              const SizedBox(width: 8),
              // The same slot on every row, so the columns line up.
              SizedBox(
                width: 40,
                child: wrong
                    ? Tooltip(
                        message: l10n.treeIdRenameOne(binomialOf(best!.species)),
                        child: Checkbox(value: rename.contains(tip), onChanged: (_) => onToggle(tip)),
                      )
                    : null,
              ),
            ]),
          );
        }(),
    ]);
  }
}
