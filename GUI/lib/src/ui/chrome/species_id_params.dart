import 'package:flutter/material.dart';

import '../../services/species_id.dart';
import '../common.dart';
import '../fa.dart';
import '../theme.dart';
import '../widgets/dialogs.dart';
import '../widgets/kit.dart';

/// The parameters of the identification, GenBank BLAST beside the BOLD ID
/// engine; every change is applied (and saved) at once.
class SpeciesIdParams extends StatefulWidget {
  const SpeciesIdParams({super.key, this.blastNote, this.boldNote, this.blastOn = true, this.boldOn = true});

  /// What each source will do with the genes chosen (setup dialog), shown
  /// under its title; a source that will not run is dimmed.
  final String? blastNote, boldNote;
  final bool blastOn, boldOn;

  @override
  State<SpeciesIdParams> createState() => _SpeciesIdParamsState();
}

class _SpeciesIdParamsState extends State<SpeciesIdParams> {
  late final SpeciesIdOptions _start = ProjectScope.read(context).speciesIdOptions;
  late final _hits = TextEditingController(text: '${_start.hitlistSize}');
  late final _expect = TextEditingController(text: '${_start.expect}');
  late final _word = TextEditingController(text: _start.wordSize?.toString() ?? '');
  late final _cover = TextEditingController(text: (_start.minCoverage * 100).round().toString());
  late final _entrez = TextEditingController(text: _start.entrezQuery);
  late final _mi = TextEditingController(text: (_start.boldMinIdentity * 100).round().toString());
  late final _mo = TextEditingController(text: '${_start.boldMinOverlap}');
  late final _maxh = TextEditingController(text: '${_start.boldMaxHits}');

  @override
  void dispose() {
    for (final c in [_hits, _expect, _word, _cover, _entrez, _mi, _mo, _maxh]) {
      c.dispose();
    }
    super.dispose();
  }

  void _set(SpeciesIdOptions Function(SpeciesIdOptions o) f) {
    final c = ProjectScope.read(context);
    c.setSpeciesIdOptions(f(c.speciesIdOptions));
  }

  /// A number in [min, max], or null (left as it was) when it is not.
  num? _num(String v, num min, num max) {
    final n = num.tryParse(v.trim().replaceAll(',', '.'));
    return n == null || n < min || n > max ? null : n;
  }

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    final o = c.speciesIdOptions;
    Widget field(String label, String info, TextEditingController ctl, void Function(String) onChanged, {String? hint}) => Field(
          label: label,
          info: info,
          child: TextField(
            controller: ctl,
            decoration: InputDecoration(hintText: hint),
            onChanged: onChanged,
          ),
        );
    // Two fields side by side.
    Widget pair(Widget a, Widget b) => Row(crossAxisAlignment: CrossAxisAlignment.start, children: [
          Expanded(child: a),
          const SizedBox(width: 10),
          Expanded(child: b),
        ]);
    Widget card(String title, String? note, bool on, List<Widget> children) => Container(
          padding: const EdgeInsets.all(12),
          decoration: BoxDecoration(
            color: on ? t.accentSoft : t.surface2,
            borderRadius: BorderRadius.circular(10),
            border: Border.all(color: on ? t.accentSoftBorder : t.border),
          ),
          child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, mainAxisSize: MainAxisSize.min, children: [
            Text(title, style: TextStyle(fontWeight: FontWeight.w700, color: on ? t.accentText : t.textMuted)),
            if (note != null) Text(note, style: TextStyle(fontSize: 12, color: t.textMuted)),
            const SizedBox(height: 6),
            for (final (i, w) in children.indexed) ...[if (i > 0) const SizedBox(height: 6), w],
          ]),
        );

    // Both cards the same height, side by side.
    return EqualHeightRow(gap: 12, children: [
      card('GenBank (NCBI BLAST)', widget.blastNote, widget.blastOn, [
        pair(
          Field(
            label: l10n.idpDatabase,
            info: l10n.idpBlastDbInfo,
            child: SplaceSelect<String>(
              value: o.blastDatabase,
              options: [for (final d in SpeciesIdOptions.blastDatabases) SelectOption(d, d)],
              onChanged: (v) => _set((o) => o.copyWith(blastDatabase: v)),
            ),
          ),
          Field(
            label: l10n.idpProgram,
            info: l10n.idpProgramInfo,
            child: SplaceSelect<String>(
              value: o.blastProgram,
              options: const [SelectOption('megablast', 'megablast'), SelectOption('blastn', 'blastn')],
              onChanged: (v) => _set((o) => o.copyWith(blastProgram: v)),
            ),
          ),
        ),
        pair(
          field(l10n.idpHits, l10n.idpHitsInfo, _hits, (v) {
            final n = _num(v, 1, 5000);
            if (n != null) _set((o) => o.copyWith(hitlistSize: n.toInt()));
          }),
          field(l10n.idpExpect, l10n.idpExpectInfo, _expect, (v) {
            final n = _num(v, 0, 1000);
            if (n != null) _set((o) => o.copyWith(expect: n.toDouble()));
          }),
        ),
        pair(
          field(l10n.idpWord, l10n.idpWordInfo, _word, (v) {
            if (v.trim().isEmpty) return _set((o) => o.copyWith(clearWordSize: true));
            final n = _num(v, 7, 256);
            if (n != null) _set((o) => o.copyWith(wordSize: n.toInt()));
          }, hint: o.blastProgram == 'megablast' ? '28' : '11'),
          field(l10n.idpCoverage, l10n.idpCoverageInfo, _cover, (v) {
            final n = _num(v, 10, 100);
            if (n != null) _set((o) => o.copyWith(minCoverage: n / 100));
          }),
        ),
        field(l10n.idpEntrez, l10n.idpEntrezInfo, _entrez, (v) => _set((o) => o.copyWith(entrezQuery: v)), hint: 'Bufonidae[Organism]'),
      ]),
      card('BOLD Systems (ID engine)', widget.boldNote, widget.boldOn, [
        Field(
          label: l10n.idpDatabase,
          info: l10n.idpBoldDbInfo,
          child: SplaceSelect<String>(
            value: o.boldDatabase,
            options: [
              for (final e in SpeciesIdOptions.boldDatabases.entries) SelectOption(e.key, e.key, detail: e.value.join(', ')),
            ],
            onChanged: (v) => _set((o) => o.copyWith(boldDatabase: v)),
          ),
        ),
        pair(
          field(l10n.idpMinIdentity, l10n.idpMinIdentityInfo, _mi, (v) {
            final n = _num(v, 70, 100);
            if (n != null) _set((o) => o.copyWith(boldMinIdentity: n / 100));
          }),
          field(l10n.idpMinOverlap, l10n.idpMinOverlapInfo, _mo, (v) {
            final n = _num(v, 1, 1000);
            if (n != null) _set((o) => o.copyWith(boldMinOverlap: n.toInt()));
          }),
        ),
        field(l10n.idpMaxHits, l10n.idpMaxHitsInfo, _maxh, (v) {
          final n = _num(v, 1, 1000);
          if (n != null) _set((o) => o.copyWith(boldMaxHits: n.toInt()));
        }),
        Text(l10n.idpBoldMarkers(o.boldGenes.join(', ')), style: TextStyle(fontSize: 12, color: t.textMuted)),
      ]),
    ]);
  }
}

/// The parameters on their own (barcoding card).
Future<void> showSpeciesIdParams(BuildContext context) => showDialog(
      context: context,
      builder: (_) => ProjectScope(
        controller: ProjectScope.read(context),
        child: Builder(
          builder: (context) => SplaceDialog(
            icon: FA.sliders,
            title: context.l10n.idpTitle,
            subtitle: context.l10n.idpSubtitle,
            actions: [FilledButton(onPressed: () => Navigator.pop(context), child: Text(context.l10n.close))],
            child: const SpeciesIdParams(),
          ),
        ),
      ),
    );
