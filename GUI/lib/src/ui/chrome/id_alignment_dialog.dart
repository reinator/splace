import 'package:flutter/material.dart';
import 'package:flutter/services.dart';

import '../../services/species_id.dart';
import '../common.dart';
import '../fa.dart';
import '../theme.dart';
import '../widgets/dialogs.dart';
import '../widgets/kit.dart';

/// The alignments behind an identification: the sequence that was sent
/// (the query) and each hit of GenBank (BLAST) and BOLD aligned to it, with
/// identity, coverage, mismatches, gaps, e-value and bit score. [species]
/// limits the hits to one candidate; otherwise the best ones are shown.
Future<void> showIdAlignments(BuildContext context, RecordIdentification r, {String? species, String? name}) => showDialog(
      context: context,
      builder: (_) => ProjectScope(
        controller: ProjectScope.read(context),
        child: _AlignmentsDialog(r: r, species: species, name: name ?? r.currentName),
      ),
    );

class _AlignmentsDialog extends StatefulWidget {
  const _AlignmentsDialog({required this.r, required this.species, required this.name});
  final RecordIdentification r;
  final String? species;
  final String name;

  @override
  State<_AlignmentsDialog> createState() => _AlignmentsDialogState();
}

class _AlignmentsDialogState extends State<_AlignmentsDialog> {
  late String _gene = widget.r.queries.keys.firstOrNull ?? '';
  IdSource? _source;

  /// Hits of [_gene] (one candidate, or every candidate's best), best first.
  List<IdHit> get _hits {
    final cands = widget.r.candidates.where((c) => widget.species == null || sameSpecies(c.species, widget.species!));
    final all = [
      for (final c in cands)
        for (final h in c.hits[_gene] ?? const <IdHit>[])
          if (_source == null || h.source == _source) h,
    ]..sort((a, b) => b.identity != a.identity ? b.identity.compareTo(a.identity) : (b.bitScore ?? 0).compareTo(a.bitScore ?? 0));
    return widget.species == null ? all.take(30).toList() : all;
  }

  @override
  Widget build(BuildContext context) {
    final l10n = context.l10n;
    final t = context.tokens;
    final query = widget.r.queries[_gene] ?? '';
    final hits = _hits;
    return SplaceDialog(
      icon: FA.alignLeft,
      title: l10n.idAlnTitle,
      subtitle: [widget.name, widget.r.accession, if (widget.species != null) widget.species!].join(' · '),
      size: DialogSize.full,
      // The list scrolls itself.
      scroll: false,
      actions: [FilledButton(onPressed: () => Navigator.pop(context), child: Text(l10n.close))],
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        Wrap(spacing: 8, runSpacing: 8, crossAxisAlignment: WrapCrossAlignment.center, children: [
          for (final g in widget.r.queries.keys)
            ChipButton(label: g, style: g == _gene ? ChipButtonStyle.accent : ChipButtonStyle.neutral, onPressed: () => setState(() => _gene = g)),
          const SizedBox(width: 12),
          for (final (s, label) in [(null, l10n.idAlnAll), (IdSource.genbank, 'GenBank'), (IdSource.bold, 'BOLD')])
            ChipButton(label: label, style: s == _source ? ChipButtonStyle.soft : ChipButtonStyle.neutral, onPressed: () => setState(() => _source = s)),
        ]),
        const SizedBox(height: 12),
        Expanded(
          child: ListView(children: [
            // The sequence that was sent.
            _Card(
              title: Text.rich(TextSpan(children: [
                TextSpan(text: '${l10n.idAlnQuery}  ', style: TextStyle(fontWeight: FontWeight.w700, color: t.text)),
                TextSpan(text: '$_gene · ${query.length} bp', style: TextStyle(color: t.textMuted)),
              ])),
              trailing: IconButton(
                tooltip: l10n.idAlnCopy,
                icon: Icon(FA.copy, size: 14, color: t.textMuted),
                onPressed: () => Clipboard.setData(ClipboardData(text: '>${widget.r.accession}_$_gene\n$query')),
              ),
              child: _Sequence(query),
            ),
            if (hits.isEmpty)
              Padding(padding: const EdgeInsets.all(24), child: Center(child: Text(l10n.idAlnNone, style: TextStyle(color: t.textMuted))))
            else
              for (final (i, h) in hits.indexed) _HitCard(key: ValueKey('${h.source}${h.accession}${h.qFrom}'), hit: h, query: query, open: i < 3),
          ]),
        ),
      ]),
    );
  }
}

class _Card extends StatelessWidget {
  const _Card({required this.title, required this.child, this.trailing});
  final Widget title;
  final Widget? trailing;
  final Widget child;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    return Container(
      margin: const EdgeInsets.only(bottom: 10),
      padding: const EdgeInsets.fromLTRB(14, 10, 10, 12),
      decoration: BoxDecoration(color: t.surface, borderRadius: BorderRadius.circular(12), border: Border.all(color: t.border)),
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        Row(children: [Expanded(child: title), ?trailing]),
        const SizedBox(height: 8),
        child,
      ]),
    );
  }
}

/// A sequence in lines of 60 with the position of each line.
class _Sequence extends StatelessWidget {
  const _Sequence(this.seq);
  final String seq;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final style = TextStyle(fontFamily: monoFont, fontSize: 12, height: 1.45, color: t.text);
    return SelectionArea(
      child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
        for (var i = 0; i < seq.length; i += 60)
          Text('${'${i + 1}'.padLeft(5)}  ${seq.substring(i, (i + 60).clamp(0, seq.length))}', style: style),
      ]),
    );
  }
}

/// One hit: its numbers and its alignment with the query (fetched and
/// aligned on demand for BOLD).
class _HitCard extends StatefulWidget {
  const _HitCard({super.key, required this.hit, required this.query, required this.open});
  final IdHit hit;
  final String query;
  final bool open;

  @override
  State<_HitCard> createState() => _HitCardState();
}

class _HitCardState extends State<_HitCard> {
  late bool _open = widget.open;
  Future<void>? _loading;
  String? _error;

  void _ensure() {
    final h = widget.hit;
    if (h.hasAlignment || _loading != null) return;
    final c = ProjectScope.read(context);
    _loading = c.loadBoldAlignment(h, widget.query).then((_) {
      if (mounted) setState(() {});
    }).catchError((Object e) {
      if (mounted) setState(() => _error = '$e');
    });
  }

  @override
  Widget build(BuildContext context) {
    final l10n = context.l10n;
    final t = context.tokens;
    final h = widget.hit;
    if (_open) _ensure();
    final cover = widget.query.isEmpty || h.qTo == 0 ? null : 100 * (h.qTo - h.qFrom + 1) / widget.query.length;
    String num(double v) => v == v.roundToDouble() ? v.toStringAsFixed(0) : v.toStringAsFixed(1);
    final facts = [
      '${l10n.idAlnIdentity} ${num(h.identity)} %',
      if (cover != null) '${l10n.idAlnCoverage} ${num(cover)} %',
      if (h.alignLength > 0) '${h.alignLength} bp',
      if (h.mismatches != null) '${l10n.idAlnMismatches} ${h.mismatches}',
      '${l10n.idAlnGaps} ${h.gaps}',
      if (h.evalue != null) 'e ${h.evalue == 0 ? '0' : h.evalue!.toStringAsExponential(1)}',
      if (h.bitScore != null) '${h.bitScore!.round()} bits',
      if (h.qTo > 0) '${l10n.idAlnQuery} ${h.qFrom}–${h.qTo} · ${l10n.idAlnHit} ${h.sFrom}–${h.sTo}${h.reversed ? ' (−)' : ''}',
    ];
    return _Card(
      title: InkWell(
        onTap: () => setState(() => _open = !_open),
        child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
          Wrap(spacing: 8, crossAxisAlignment: WrapCrossAlignment.center, children: [
            Icon(_open ? FA.chevronDown : FA.chevronRight, size: 10, color: t.textMuted),
            TagChip(h.source == IdSource.bold ? 'BOLD' : 'GenBank', color: h.source == IdSource.bold ? SplaceRed.c700 : t.accent),
            Text(h.species, style: TextStyle(fontStyle: FontStyle.italic, fontWeight: FontWeight.w700, color: t.text)),
            Text(h.accession, style: TextStyle(fontFamily: monoFont, fontSize: 12, color: t.textMuted)),
          ]),
          const SizedBox(height: 4),
          Text(facts.join('   ·   '), style: TextStyle(fontSize: 12, color: t.textMuted)),
          if (h.title.isNotEmpty && h.source == IdSource.genbank)
            Text(h.title, maxLines: 1, overflow: TextOverflow.ellipsis, style: TextStyle(fontSize: 11.5, color: t.textFaint)),
        ]),
      ),
      child: !_open
          ? const SizedBox.shrink()
          : h.hasAlignment
              ? _Alignment(hit: h)
              : _error != null
                  ? Text(l10n.idAlnLoadError(_error!), style: const TextStyle(fontSize: 12.5, color: SplaceTokens.danger))
                  : Row(children: [
                      SizedBox(width: 12, height: 12, child: CircularProgressIndicator(strokeWidth: 2, color: t.accent)),
                      const SizedBox(width: 8),
                      Text(l10n.idAlnLoading, style: TextStyle(fontSize: 12.5, color: t.textMuted)),
                    ]),
    );
  }
}

/// Query and hit in blocks of 60 columns, the match line between them;
/// differences and gaps are marked.
class _Alignment extends StatelessWidget {
  const _Alignment({required this.hit});
  final IdHit hit;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final q = hit.qseq!, s = hit.hseq!;
    final base = TextStyle(fontFamily: monoFont, fontSize: 12, height: 1.35, color: t.text);
    final diff = base.copyWith(color: Colors.white, backgroundColor: SplaceRed.c600);
    final gap = base.copyWith(color: t.textFaint, backgroundColor: t.surface3);
    final blocks = <Widget>[];
    var qPos = hit.qFrom > 0 ? hit.qFrom : 1;
    var sPos = hit.sFrom > 0 ? hit.sFrom : 1;
    final step = hit.reversed ? -1 : 1;
    // Label and position, the same width on every line (monospaced).
    String prefix(String label, int pos) => '${label.padRight(6)}${'$pos'.padLeft(5)}  ';
    for (var i = 0; i < q.length; i += 60) {
      final end = (i + 60).clamp(0, q.length);
      final qs = q.substring(i, end), ss = s.substring(i, end);
      final qBases = qs.replaceAll('-', '').length, sBases = ss.replaceAll('-', '').length;
      final qEnd = qPos + qBases - 1, sEnd = sPos + step * (sBases - 1);
      List<InlineSpan> row(String a, String b) => [
            for (var k = 0; k < a.length; k++)
              TextSpan(
                text: a[k],
                style: a[k] == '-' ? gap : (a[k].toUpperCase() != b[k].toUpperCase() && b[k] != '-' ? diff : null),
              ),
          ];
      final mid = StringBuffer();
      for (var k = 0; k < qs.length; k++) {
        mid.write(qs[k].toUpperCase() == ss[k].toUpperCase() && qs[k] != '-' ? '|' : ' ');
      }
      blocks.add(Padding(
        padding: const EdgeInsets.only(bottom: 8),
        child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
          Text.rich(TextSpan(style: base, children: [TextSpan(text: prefix('Query', qPos)), ...row(qs, ss), TextSpan(text: '  $qEnd')])),
          // The match line starts under the first base: as wide as the
          // "Query     1  " prefix.
          Text('${' ' * prefix('', 0).length}${mid.toString()}', style: base.copyWith(color: t.textMuted)),
          Text.rich(TextSpan(style: base, children: [TextSpan(text: prefix('Hit', sPos)), ...row(ss, qs), TextSpan(text: '  $sEnd')])),
        ]),
      ));
      qPos = qEnd + 1;
      sPos = sEnd + step;
    }
    return SelectionArea(
      child: SingleChildScrollView(
        scrollDirection: Axis.horizontal,
        child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: blocks),
      ),
    );
  }
}
