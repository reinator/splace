import 'dart:async';
import 'dart:isolate';
import 'dart:math' as math;

import 'package:flutter/material.dart';

import '../../core/fasta.dart';
import '../../core/locus_organizer.dart';
import '../../state/project_controller.dart';
import '../common.dart';
import '../fa.dart';
import '../theme.dart';
import '../widgets/dialogs.dart';
import '../widgets/kit.dart';
import 'fields_section.dart' show FieldRuleEditor;

/// Colours of the two parts of a header.
const sampleColor = SplaceBlue.c500;
const markerColor = SplacePalette.red500;

/// "Sample and marker": after FASTA files are read, the user tells which
/// part of the headers is the sample (ID) and which the marker, with a
/// live preview of the samples × markers matrix. Confirming writes one
/// FASTA per marker holding every sample.
Future<void> showOrganizerDialog(BuildContext context) async {
  final c = ProjectScope.read(context);
  if (c.inputs.isEmpty) return;
  final config = await showDialog<OrganizerConfig>(
    context: context,
    barrierDismissible: false,
    builder: (_) => ProjectScope(controller: c, child: _OrganizerDialog(controller: c)),
  );
  if (config == null) {
    c.dismissOrganizerReview();
    return;
  }
  await c.applyOrganizer(config);
}

// Top level, so the isolate closure captures only its arguments (never
// the dialog state, which cannot be sent).
Future<OrganizerResult> _previewInIsolate(List<FastaInput> inputs, OrganizerConfig cfg) =>
    Isolate.run(() => organizeLoci(inputs, cfg));

/// Delimiters offered as switches (any other can be typed).
const _commonDelimiters = ['-', '_', '|', '.', ' ', ':', ';', '/'];

enum _Role { none, sample, marker }

class _OrganizerDialog extends StatefulWidget {
  const _OrganizerDialog({required this.controller});
  final ProjectController controller;

  @override
  State<_OrganizerDialog> createState() => _OrganizerDialogState();
}

class _OrganizerDialogState extends State<_OrganizerDialog> {
  late OrganizerConfig _cfg = widget.controller.organizer;
  late final _delims = TextEditingController(text: _headerDelimiters);

  /// Headers only (no sequences): the preview is computed in an isolate
  /// on every change, so it must stay small even for large data sets.
  late final List<FastaInput> _light = [
    for (final f in widget.controller.inputs) FastaInput(f.fileName, [for (final r in f.records) SeqRecord(r.header, '')]),
  ];
  late final int _records = _light.fold(0, (a, f) => a + f.records.length);

  /// A few headers: one of each format first (`x-x`, `x_x`, `x|x|x`…,
  /// so different delimiters show up), then from different files.
  late final List<(String, SeqRecord)> _examples = () {
    final files = _light.where((f) => f.records.isNotEmpty).toList();
    final out = <(String, SeqRecord)>[];
    final shapes = <String>{};
    var scanned = 0;
    for (final f in files) {
      for (final r in f.records) {
        if (++scanned > 50000 || out.length >= 4) break;
        if (shapes.add(_shape(r.header))) out.add((f.fileName, r));
      }
    }
    _formats = shapes.length;
    final seen = {for (final (_, r) in out) r};
    final all = files.length >= 5
        ? [for (var i = 0; i < 5; i++) () {
            final f = files[(i * (files.length - 1) / 4).round()];
            return (f.fileName, f.records[math.min(f.records.length - 1, i * 7)]);
          }()]
        : [for (final f in files) for (final r in f.records.take(20)) (f.fileName, r)];
    for (final e in all) {
      if (out.length >= 5) break;
      if (seen.add(e.$2)) out.add(e);
    }
    return out;
  }();

  /// Number of header formats found.
  int _formats = 1;

  /// Format of a header: letters and digits collapsed (`ITV01-L12` -> `x-x`).
  static String _shape(String header) => headerShape(header);

  /// Every header format: how many sequences, in which files, and the
  /// first ones (to show them), most common format first.
  late final List<({String shape, int count, Set<String> files, List<(String, SeqRecord)> records})> _shapes = () {
    final count = <String, int>{};
    final files = <String, Set<String>>{};
    final records = <String, List<(String, SeqRecord)>>{};
    for (final f in _light) {
      for (final r in f.records) {
        final sh = headerShape(r.header);
        count[sh] = (count[sh] ?? 0) + 1;
        (files[sh] ??= {}).add(f.fileName);
        final list = records[sh] ??= [];
        if (list.length < 200) list.add((f.fileName, r));
      }
    }
    final keys = count.keys.toList()..sort((a, b) => count[b]!.compareTo(count[a]!));
    return [for (final k in keys) (shape: k, count: count[k]!, files: files[k]!, records: records[k]!)];
  }();

  /// Formats whose sequences are listed.
  final Set<String> _openFormats = {};

  OrganizerResult? _preview;
  bool _computing = false;
  int _generation = 0;
  Timer? _debounce;
  bool _advanced = false;
  String _sampleQuery = '';

  @override
  void initState() {
    super.initState();
    _examples.length; // also counts the header formats
    _update(_cfg, now: true);
  }

  @override
  void dispose() {
    _debounce?.cancel();
    _delims.dispose();
    super.dispose();
  }

  String get _headerDelimiters =>
      _cfg.locus.source == FieldSource.header ? _cfg.locus.delimiters : (_cfg.sample.source == FieldSource.header ? _cfg.sample.delimiters : '-_|');

  void _update(OrganizerConfig cfg, {bool now = false}) {
    setState(() => _cfg = cfg);
    _debounce?.cancel();
    _debounce = Timer(now ? Duration.zero : const Duration(milliseconds: 180), () async {
      final generation = ++_generation;
      setState(() => _computing = true);
      // Headers only: tens of thousands take a few tens of ms; larger
      // sets are counted in an isolate so typing stays smooth.
      final r = _records <= 30000 ? organizeLoci(_light, cfg) : await _previewInIsolate(_light, cfg);
      if (!mounted || generation != _generation) return;
      setState(() {
        _preview = r;
        _computing = false;
      });
    });
  }

  // ---------------------------------------------------------------------
  // Header parts
  // ---------------------------------------------------------------------

  void _setDelimiters(String d) {
    if (_delims.text != d) _delims.text = d;
    FieldRule apply(FieldRule r) => r.source == FieldSource.header ? r.copyWith(delimiters: d) : r;
    _update(_cfg.copyWith(sample: apply(_cfg.sample), locus: apply(_cfg.locus)));
  }

  void _toggleDelimiter(String ch) {
    final d = _headerDelimiters;
    _setDelimiters(d.contains(ch) ? d.replaceAll(ch, '') : d + ch);
  }

  /// Token positions of [rule] in a header split into [n] tokens.
  static Set<int> _positions(FieldRule rule, int n) {
    if (rule.source != FieldSource.header || n == 0) return const {};
    final s = rule.start < 0 ? n + rule.start : rule.start;
    final e = rule.end == null ? n : (rule.end! < 0 ? n + rule.end! : rule.end!);
    if (s < 0 || s >= n || e <= s) return const {};
    return {for (var i = s; i < math.min(e, n); i++) i};
  }

  /// A header rule covering [positions] (made contiguous) of a header of
  /// [n] tokens: counted from the end when it reaches the last token, so
  /// headers with more or fewer parts still work.
  FieldRule _ruleFor(Set<int> positions, int n, String delimiters, {required bool otherAtEnd}) {
    final s = positions.reduce(math.min), e = positions.reduce(math.max) + 1;
    final toEnd = e == n;
    final start = toEnd && s > 0 ? s - n : s;
    final int? end = toEnd ? null : (s == 0 && otherAtEnd && n - e <= 2 ? e - n : e);
    return FieldRule(source: FieldSource.header, start: start, end: end, delimiters: delimiters);
  }

  /// Click on token [i] of a header of [n] tokens: none → sample → marker.
  void _cycle(int i, int n) {
    final (sample, locus) = _cycled(_cfg.sample, _cfg.locus, i, n);
    _update(_cfg.copyWith(sample: sample, locus: locus));
  }

  /// The same click on the example of one header format: that format gets
  /// its own parts (starting from the main rule).
  void _cycleFormat(String shape, int i, int n) {
    final f = _cfg.formats[shape];
    final (sample, locus) = _cycled(f?.sample ?? _cfg.sample, f?.locus ?? _cfg.locus, i, n);
    _update(_cfg.copyWith(formats: {..._cfg.formats, shape: FormatRule(sample: sample, locus: locus)}));
  }

  void _setFormat(String shape, FormatRule? rule) {
    final formats = {..._cfg.formats};
    if (rule == null) {
      formats.remove(shape);
    } else {
      formats[shape] = rule;
    }
    _update(_cfg.copyWith(formats: formats));
  }

  (FieldRule, FieldRule) _cycled(FieldRule sampleRule, FieldRule locusRule, int i, int n) {
    final d = _headerDelimiters;
    final sample = {..._positions(sampleRule, n)};
    final marker = {..._positions(locusRule, n)};
    final role = sample.contains(i) ? _Role.sample : (marker.contains(i) ? _Role.marker : _Role.none);
    switch (role) {
      case _Role.none:
        sample.add(i);
        marker.remove(i);
      case _Role.sample:
        sample.remove(i);
        marker.add(i);
      case _Role.marker:
        marker.remove(i);
    }
    // Ranges are contiguous: a part cannot hold tokens of the other.
    if (sample.isNotEmpty && marker.isNotEmpty) {
      final sMin = sample.reduce(math.min), sMax = sample.reduce(math.max);
      marker.removeWhere((m) => m >= sMin && m <= sMax);
      if (marker.isNotEmpty) {
        final mMin = marker.reduce(math.min), mMax = marker.reduce(math.max);
        if (role == _Role.sample) sample.removeWhere((s) => s >= mMin && s <= mMax);
      }
    }
    var newSample = sampleRule, newLocus = locusRule;
    final markerAtEnd = marker.contains(n - 1), sampleAtEnd = sample.contains(n - 1);
    if (sample.isNotEmpty) {
      newSample = _ruleFor(sample, n, d, otherAtEnd: markerAtEnd);
    } else if (sampleRule.source == FieldSource.header) {
      // No sample part left in the header: take it from the file name.
      newSample = const FieldRule.fileName(delimiters: '');
    }
    newLocus = marker.isNotEmpty
        ? _ruleFor(marker, n, d, otherAtEnd: sampleAtEnd)
        : FieldRule(source: FieldSource.header, start: 0, end: 0, delimiters: d);
    return (newSample, newLocus);
  }

  void _renameSample(String from, String to) => _renameSamples({from: to});

  /// Several renames at once (sample shown -> new name).
  void _renameSamples(Map<String, String> renames) {
    final aliases = {..._cfg.sampleAliases};
    for (final e in renames.entries) {
      _alias(aliases, e.key, e.value.trim());
    }
    _update(_cfg.copyWith(sampleAliases: aliases));
  }

  static void _alias(Map<String, String> aliases, String from, String to) {
    for (final k in aliases.keys.toList()) {
      if (aliases[k] == from) aliases[k] = to.isEmpty ? k : to;
    }
    if (to.isEmpty || to == from) {
      aliases.remove(from);
    } else {
      aliases[from] = to;
    }
    aliases.removeWhere((k, v) => k == v);
  }

  @override
  Widget build(BuildContext context) {
    final l10n = context.l10n;
    final t = context.tokens;
    final r = _preview;
    final plausible = r != null && ProjectController.plausibleOrganization(r);
    final canApply = r != null && plausible && r.loci.isNotEmpty && !_computing;
    final headerSample = _cfg.sample.source == FieldSource.header;
    final simple = (_cfg.sample.source == FieldSource.header || _cfg.sample.source == FieldSource.fileName) &&
        _cfg.locus.source == FieldSource.header;

    Widget section(String title, {String? info, required Widget child}) => Padding(
          padding: const EdgeInsets.only(bottom: 20),
          child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
            Row(children: [
              Kicker(title),
              if (info != null) ...[const SizedBox(width: 6), InfoDot(info)],
            ]),
            const SizedBox(height: 8),
            child,
          ]),
        );

    // 1. Delimiters.
    final delims = _headerDelimiters;
    final delimiterRow = Wrap(spacing: 6, runSpacing: 6, crossAxisAlignment: WrapCrossAlignment.center, children: [
      for (final ch in _commonDelimiters)
        ChipButton(
          label: ch == ' ' ? l10n.orgSpace : ch,
          style: delims.contains(ch) ? ChipButtonStyle.accent : ChipButtonStyle.neutral,
          onPressed: () => _toggleDelimiter(ch),
        ),
      SizedBox(
        width: 120,
        child: TextField(
          controller: _delims,
          style: const TextStyle(fontFamily: monoFont, fontSize: 13),
          decoration: InputDecoration(isDense: true, hintText: l10n.orgOther),
          onChanged: _setDelimiters,
        ),
      ),
    ]);

    // 2. Example headers, split, with the role of each part.
    Widget exampleRow(String file, SeqRecord record) {
      final tokens = tokenize(record.header, delims);
      final n = tokens.length;
      final sPos = _positions(_cfg.sample, n), mPos = _positions(_cfg.locus, n);
      final sample = _cfg.sample.extract(file, record);
      final marker = _cfg.locus.extract(file, record);
      return Container(
        padding: const EdgeInsets.symmetric(horizontal: 12, vertical: 8),
        margin: const EdgeInsets.only(bottom: 6),
        decoration: BoxDecoration(color: t.surface2, borderRadius: BorderRadius.circular(10), border: Border.all(color: t.border)),
        child: Row(children: [
          SizedBox(
            width: 150,
            child: Text(file, maxLines: 1, overflow: TextOverflow.ellipsis, style: TextStyle(fontSize: 11.5, color: t.textFaint)),
          ),
          Expanded(
            child: Wrap(spacing: 4, runSpacing: 4, crossAxisAlignment: WrapCrossAlignment.center, children: [
              for (var i = 0; i < n && i < 16; i++)
                _TokenChip(
                  text: tokens[i],
                  role: sPos.contains(i) ? _Role.sample : (mPos.contains(i) ? _Role.marker : _Role.none),
                  onTap: simple ? () => _cycle(i, n) : null,
                ),
            ]),
          ),
          const SizedBox(width: 12),
          SizedBox(
            width: 250,
            child: Text.rich(
              TextSpan(style: const TextStyle(fontSize: 12.5), children: [
                TextSpan(text: sample ?? '∅', style: TextStyle(color: sample == null ? SplaceTokens.danger : sampleColor, fontWeight: FontWeight.w700)),
                TextSpan(text: '  ·  ', style: TextStyle(color: t.textFaint)),
                TextSpan(text: marker ?? '∅', style: TextStyle(color: marker == null ? SplaceTokens.danger : markerColor, fontWeight: FontWeight.w700)),
              ]),
              maxLines: 1,
              overflow: TextOverflow.ellipsis,
            ),
          ),
        ]),
      );
    }

    // File names that differ from the sample read in their headers.
    final mismatches = <(String, String)>[];
    if (headerSample && r != null) {
      for (final e in r.sampleFiles.entries) {
        if (e.value.length != 1) continue;
        final stem = fileStem(e.value.first);
        String norm(String s) => stripLeadingZeros(s.toUpperCase());
        final onlySample = r.sampleFiles.values.where((f) => f.contains(e.value.first)).length == 1;
        // Written differently (case aside) but close: ITV01054.fasta with
        // headers ITV1054.
        if (onlySample && stem.toUpperCase() != e.key.toUpperCase() && stem.contains(RegExp(r'\d')) && _similar(norm(stem), norm(e.key))) {
          mismatches.add((e.value.first, e.key));
        }
      }
    }

    // 3. Header formats, when the headers are written in more than one way.
    Widget formatCard(({String shape, int count, Set<String> files, List<(String, SeqRecord)> records}) f) {
      final own = _cfg.formats[f.shape];
      final skip = own?.skip ?? false;
      final sampleRule = own?.sample ?? _cfg.sample, locusRule = own?.locus ?? _cfg.locus;
      final (file, record) = f.records.first;
      final tokens = tokenize(record.header, delims);
      final n = tokens.length;
      final sPos = _positions(sampleRule, n), mPos = _positions(locusRule, n);
      final sample = sampleRule.extract(file, record), marker = locusRule.extract(file, record);
      final open = _openFormats.contains(f.shape);
      return Container(
        margin: const EdgeInsets.only(bottom: 8),
        padding: const EdgeInsets.all(12),
        decoration: BoxDecoration(
          color: skip ? t.surface2 : t.surface,
          borderRadius: BorderRadius.circular(10),
          border: Border.all(color: !skip && (sample == null || marker == null) ? SplaceTokens.danger : t.border),
        ),
        child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
          Row(children: [
            Container(
              padding: const EdgeInsets.symmetric(horizontal: 8, vertical: 3),
              decoration: BoxDecoration(color: t.surface2, borderRadius: BorderRadius.circular(8), border: Border.all(color: t.border)),
              child: Text(f.shape, style: TextStyle(fontFamily: monoFont, fontSize: 12.5, fontWeight: FontWeight.w700, color: t.text)),
            ),
            const SizedBox(width: 10),
            Expanded(
              child: Text(
                l10n.orgFormatCount('${f.count}', '${f.files.length}'),
                style: TextStyle(fontSize: 12.5, color: t.textMuted),
              ),
            ),
            if (skip)
              TagChip(l10n.orgFormatSkipped, color: t.textFaint)
            else if (own != null)
              TagChip(l10n.orgFormatOwnRule, color: t.accent),
          ]),
          const SizedBox(height: 10),
          Opacity(
            opacity: skip ? 0.45 : 1,
            child: Row(children: [
              Expanded(
                child: Wrap(spacing: 4, runSpacing: 4, crossAxisAlignment: WrapCrossAlignment.center, children: [
                  for (var i = 0; i < n && i < 16; i++)
                    _TokenChip(
                      text: tokens[i],
                      role: sPos.contains(i) ? _Role.sample : (mPos.contains(i) ? _Role.marker : _Role.none),
                      onTap: simple && !skip ? () => _cycleFormat(f.shape, i, n) : null,
                    ),
                ]),
              ),
              const SizedBox(width: 12),
              SizedBox(
                width: 250,
                child: Text.rich(
                  TextSpan(style: const TextStyle(fontSize: 12.5), children: [
                    TextSpan(text: sample ?? '∅', style: TextStyle(color: sample == null ? SplaceTokens.danger : sampleColor, fontWeight: FontWeight.w700)),
                    TextSpan(text: '  ·  ', style: TextStyle(color: t.textFaint)),
                    TextSpan(text: marker ?? '∅', style: TextStyle(color: marker == null ? SplaceTokens.danger : markerColor, fontWeight: FontWeight.w700)),
                  ]),
                  maxLines: 1,
                  overflow: TextOverflow.ellipsis,
                ),
              ),
            ]),
          ),
          const SizedBox(height: 10),
          Wrap(spacing: 8, runSpacing: 6, children: [
            ChipButton(
              label: open ? l10n.orgFormatHide : l10n.orgFormatShow('${f.count}'),
              icon: open ? FA.chevronUp : FA.chevronDown,
              onPressed: () => setState(() => open ? _openFormats.remove(f.shape) : _openFormats.add(f.shape)),
            ),
            if (own != null && !skip)
              ChipButton(label: l10n.orgFormatUseMain, icon: FA.rotateLeft, onPressed: () => _setFormat(f.shape, null)),
            ChipButton(
              label: skip ? l10n.orgFormatUse : l10n.orgFormatSkip,
              icon: skip ? FA.check : FA.xmark,
              onPressed: () => _setFormat(f.shape, skip ? null : FormatRule(sample: sampleRule, locus: locusRule, skip: true)),
            ),
          ]),
          if (open) ...[
            const SizedBox(height: 8),
            Container(
              constraints: const BoxConstraints(maxHeight: 180),
              decoration: BoxDecoration(color: t.surface2, borderRadius: BorderRadius.circular(8), border: Border.all(color: t.border)),
              child: ListView.builder(
                shrinkWrap: true,
                padding: const EdgeInsets.symmetric(horizontal: 10, vertical: 6),
                itemCount: f.records.length + (f.count > f.records.length ? 1 : 0),
                itemExtent: 20,
                itemBuilder: (context, i) => i == f.records.length
                    ? Text('+ ${f.count - f.records.length}', style: TextStyle(fontSize: 12, color: t.textMuted))
                    : Row(children: [
                        SizedBox(
                          width: 160,
                          child: Text(f.records[i].$1, maxLines: 1, overflow: TextOverflow.ellipsis, style: TextStyle(fontSize: 11.5, color: t.textFaint)),
                        ),
                        Expanded(
                          child: Text('>${f.records[i].$2.header}',
                              maxLines: 1, overflow: TextOverflow.ellipsis, style: TextStyle(fontFamily: monoFont, fontSize: 12, color: t.text)),
                        ),
                      ]),
              ),
            ),
          ],
        ]),
      );
    }

    // 4. Sample names written with different structures (ITV1054 and
    // ITV01054): the user picks one and the others are converted to it.
    final structures = <String, List<String>>{};
    for (final sm in r?.samples ?? const <String>[]) {
      (structures[nameStructure(sm)] ??= []).add(sm);
    }
    final convertible = <String, Map<String, String>>{};
    if (structures.length > 1 && structures.length <= 12) {
      for (final target in structures.keys) {
        final renames = <String, String>{
          for (final e in structures.entries)
            if (e.key != target)
              for (final sm in e.value)
                if (conformName(sm, target) case final to? when to != sm) sm: to,
        };
        if (renames.isNotEmpty) convertible[target] = renames;
      }
    }
    final structureGroups = structures.entries.toList()..sort((a, b) => b.value.length.compareTo(a.value.length));

    final samples = r == null
        ? const <String>[]
        : [for (final s in r.samples) if (_sampleQuery.isEmpty || s.toLowerCase().contains(_sampleQuery.toLowerCase())) s];

    return SplaceDialog(
      icon: FA.tableColumns,
      title: l10n.orgTitle,
      subtitle: l10n.orgSubtitle('${_light.length}', '$_records'),
      size: DialogSize.full,
      footerLeading: r == null
          ? null
          : Text(
              l10n.orgSummary('${r.samples.length}', '${r.loci.length}', '${(r.occupancy * 100).round()}'),
              style: TextStyle(fontSize: 13, fontWeight: FontWeight.w700, color: plausible ? t.text : SplaceTokens.danger),
            ),
      actions: [
        DialogCancelButton(onPressed: () => Navigator.pop(context)),
        FilledButton.icon(
          onPressed: canApply ? () => Navigator.pop(context, _cfg) : null,
          icon: const Icon(FA.fileExport, size: 14),
          label: Text(l10n.orgApply('${r?.loci.length ?? 0}')),
        ),
      ],
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        Text(l10n.orgIntro, style: TextStyle(fontSize: 13.5, color: t.textMuted, height: 1.5)),
        const SizedBox(height: 16),
        AutoGrid(minItemWidth: 380, spacing: 20, children: [
          section(l10n.orgSampleFrom,
              info: l10n.orgSampleFromInfo,
              child: Align(
                alignment: Alignment.centerLeft,
                child: Seg<bool>(
                  value: headerSample,
                  options: [(true, l10n.orgSampleHeader), (false, l10n.orgSampleFile)],
                  onChanged: (v) => _update(_cfg.copyWith(
                    sample: v
                        ? FieldRule(source: FieldSource.header, start: 0, end: 1, delimiters: delims)
                        : const FieldRule.fileName(delimiters: ''),
                  )),
                ),
              )),
          section(l10n.delimiters, info: l10n.orgDelimitersInfo, child: delimiterRow),
        ]),
        section(
          l10n.orgParts,
          info: l10n.orgPartsInfo,
          child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
            Wrap(spacing: 16, runSpacing: 4, children: [
              if (_formats > 1) Text(l10n.orgFormats('$_formats'), style: TextStyle(fontSize: 12, fontWeight: FontWeight.w700, color: t.accentText)),
              _Legend(color: sampleColor, text: l10n.orgLegendSample),
              _Legend(color: markerColor, text: l10n.orgLegendMarker),
              _Legend(color: t.textFaint, text: l10n.orgLegendIgnored),
            ]),
            const SizedBox(height: 8),
            for (final (f, rec) in _examples) exampleRow(f, rec),
            if (!simple)
              Text(l10n.orgAdvancedActive, style: TextStyle(fontSize: 12, color: t.textMuted, fontStyle: FontStyle.italic)),
            const SizedBox(height: 4),
            Align(
              alignment: Alignment.centerLeft,
              child: ChipButton(
                label: l10n.orgAdvanced,
                icon: _advanced ? FA.chevronUp : FA.chevronDown,
                onPressed: () => setState(() => _advanced = !_advanced),
              ),
            ),
            if (_advanced) ...[
              const SizedBox(height: 10),
              AutoGrid(minItemWidth: 420, spacing: 16, children: [
                OptionCard(
                  title: l10n.sampleField,
                  info: l10n.sampleFieldInfo,
                  child: FieldRuleEditor(
                    rule: _cfg.sample,
                    example: _examples.isEmpty ? ('', SeqRecord('', '')) : _examples.first,
                    onChanged: (v) => _update(_cfg.copyWith(sample: v)),
                  ),
                ),
                OptionCard(
                  title: l10n.locusField,
                  info: l10n.locusFieldInfo,
                  child: FieldRuleEditor(
                    rule: _cfg.locus,
                    example: _examples.isEmpty ? ('', SeqRecord('', '')) : _examples.first,
                    onChanged: (v) => _update(_cfg.copyWith(locus: v)),
                  ),
                ),
              ]),
            ],
          ]),
        ),
        if (_shapes.length > 1)
          section(
            l10n.orgHeaderFormats('${_shapes.length}'),
            info: l10n.orgHeaderFormatsInfo,
            child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
              for (final f in _shapes.take(12)) formatCard(f),
              if (_shapes.length > 12) Text('+ ${_shapes.length - 12}', style: TextStyle(color: t.textMuted)),
            ]),
          ),
        if (convertible.isNotEmpty)
          section(
            l10n.orgNameStructures('${structures.length}'),
            info: l10n.orgNameStructuresInfo,
            child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
              for (final e in structureGroups)
                Container(
                  margin: const EdgeInsets.only(bottom: 8),
                  padding: const EdgeInsets.all(12),
                  decoration: BoxDecoration(color: t.surface, borderRadius: BorderRadius.circular(10), border: Border.all(color: t.border)),
                  child: Row(children: [
                    // Fixed widths: the names of every row start together.
                    SizedBox(
                      width: 150,
                      child: Align(
                        alignment: Alignment.centerLeft,
                        child: Container(
                          padding: const EdgeInsets.symmetric(horizontal: 8, vertical: 3),
                          decoration: BoxDecoration(color: t.surface2, borderRadius: BorderRadius.circular(8), border: Border.all(color: t.border)),
                          child: Text(e.key,
                              maxLines: 1,
                              overflow: TextOverflow.ellipsis,
                              style: TextStyle(fontFamily: monoFont, fontSize: 12.5, fontWeight: FontWeight.w700, color: t.text)),
                        ),
                      ),
                    ),
                    const SizedBox(width: 10),
                    SizedBox(
                      width: 100,
                      child: Text(l10n.orgStructureCount(e.value.length), style: TextStyle(fontSize: 12.5, color: t.textMuted)),
                    ),
                    Expanded(
                      child: Text(
                        e.value.take(8).join(', ') + (e.value.length > 8 ? ', …' : ''),
                        maxLines: 2,
                        overflow: TextOverflow.ellipsis,
                        style: TextStyle(fontFamily: monoFont, fontSize: 12, color: t.text),
                      ),
                    ),
                    const SizedBox(width: 10),
                    if (convertible[e.key] case final renames?)
                      ChipButton(
                        label: l10n.orgStandardize(e.key),
                        icon: FA.wandMagicSparkles,
                        style: ChipButtonStyle.soft,
                        tooltip: renames.entries.take(12).map((x) => '${x.key} → ${x.value}').join('\n') + (renames.length > 12 ? '\n…' : ''),
                        onPressed: () => _renameSamples(renames),
                      ),
                  ]),
                ),
            ]),
          ),
        // Problems first: they explain an odd count.
        if (r != null && !plausible)
          Padding(
            padding: const EdgeInsets.only(bottom: 16),
            child: TonePanel(
              tone: PanelTone.danger,
              icon: FA.triangleExclamation,
              title: l10n.orgImplausible('${r.loci.length}', '${r.records}'),
              child: Text(l10n.orgImplausibleBody),
            ),
          ),
        if (r != null && r.unassigned.isNotEmpty)
          Padding(
            padding: const EdgeInsets.only(bottom: 16),
            child: TonePanel(
              tone: PanelTone.warn,
              icon: FA.triangleExclamation,
              title: l10n.orgUnassigned('${r.unassigned.length}'),
              child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
                for (final u in r.unassigned.take(4))
                  Text('${u.sourceFile}: >${u.header}', maxLines: 1, overflow: TextOverflow.ellipsis, style: const TextStyle(fontFamily: monoFont, fontSize: 12)),
              ]),
            ),
          ),
        section(
          l10n.orgVariations,
          info: l10n.orgVariationsInfo,
          child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
            AutoGrid(minItemWidth: 320, children: [
              CheckCard(
                value: _cfg.ignoreCase,
                onChanged: (v) => _update(_cfg.copyWith(ignoreCase: v)),
                title: l10n.orgIgnoreCase,
                subtitle: 'itv01054 = ITV01054',
              ),
              CheckCard(
                value: _cfg.ignoreLeadingZeros,
                onChanged: (v) => _update(_cfg.copyWith(ignoreLeadingZeros: v)),
                title: l10n.orgLeadingZeros,
                subtitle: 'ITV1054 = ITV01054 · L01 = L1',
              ),
            ]),
            if (mismatches.isNotEmpty) ...[
              const SizedBox(height: 12),
              TonePanel(
                tone: PanelTone.info,
                icon: FA.circleInfo,
                title: l10n.orgFileMismatch('${mismatches.length}'),
                child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
                  for (final (file, sample) in mismatches.take(6))
                    Padding(
                      padding: const EdgeInsets.only(bottom: 4),
                      child: Row(children: [
                        Expanded(
                          child: Text('$file  →  $sample', style: const TextStyle(fontFamily: monoFont, fontSize: 12)),
                        ),
                        ChipButton(
                          label: l10n.orgRenameTo(fileStem(file)),
                          onPressed: () => _renameSample(sample, fileStem(file)),
                        ),
                      ]),
                    ),
                  if (mismatches.length > 6) Text('+ ${mismatches.length - 6}', style: TextStyle(color: t.textMuted)),
                  const SizedBox(height: 6),
                  Wrap(spacing: 8, runSpacing: 6, children: [
                    ChipButton(
                      label: l10n.orgRenameAll,
                      style: ChipButtonStyle.soft,
                      onPressed: () {
                        final aliases = {..._cfg.sampleAliases};
                        for (final (file, sample) in mismatches) {
                          aliases[sample] = fileStem(file);
                        }
                        _update(_cfg.copyWith(sampleAliases: aliases));
                      },
                    ),
                    ChipButton(
                      label: l10n.orgUseFileNames,
                      onPressed: () => _update(_cfg.copyWith(sample: const FieldRule.fileName(delimiters: ''))),
                    ),
                  ]),
                ]),
              ),
            ],
            if (r != null && r.sampleSpellings.isNotEmpty) ...[
              const SizedBox(height: 12),
              Text(l10n.orgMerged('${r.sampleSpellings.length}'), style: TextStyle(fontSize: 12.5, fontWeight: FontWeight.w700, color: t.text)),
              const SizedBox(height: 4),
              for (final e in r.sampleSpellings.entries.take(6))
                Text('${e.key}  ←  ${e.value.where((v) => v != e.key).join(', ')}', style: TextStyle(fontFamily: monoFont, fontSize: 12, color: t.textMuted)),
            ],
          ]),
        ),
        AutoGrid(minItemWidth: 420, spacing: 20, children: [
          section(
            l10n.orgSamples('${r?.samples.length ?? 0}'),
            info: l10n.orgSamplesInfo,
            child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
              TextField(
                decoration: InputDecoration(isDense: true, prefixIcon: Icon(FA.magnifyingGlass, size: 13, color: t.textFaint), hintText: l10n.orgFindSample),
                onChanged: (v) => setState(() => _sampleQuery = v),
              ),
              const SizedBox(height: 8),
              Container(
                height: 260,
                decoration: BoxDecoration(border: Border.all(color: t.border), borderRadius: BorderRadius.circular(10)),
                child: r == null
                    ? const SizedBox()
                    : ListView.builder(
                        itemCount: samples.length,
                        itemExtent: 40,
                        itemBuilder: (context, i) {
                          final s = samples[i];
                          final markers = r.lociFor(s);
                          return Padding(
                            padding: const EdgeInsets.symmetric(horizontal: 10),
                            child: Row(children: [
                              Expanded(
                                child: _RenameField(
                                  key: ValueKey(s),
                                  name: s,
                                  onSubmitted: (v) => _renameSample(s, v),
                                ),
                              ),
                              const SizedBox(width: 10),
                              SizedBox(
                                width: 130,
                                child: Text(l10n.orgSampleMarkers('$markers', '${r.loci.length}'),
                                    textAlign: TextAlign.right, style: TextStyle(fontSize: 12, color: t.textMuted)),
                              ),
                            ]),
                          );
                        },
                      ),
              ),
            ]),
          ),
          section(
            l10n.orgMarkers('${r?.loci.length ?? 0}'),
            info: l10n.orgMarkersInfo,
            child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
              Row(children: [
                Expanded(child: Text(l10n.duplicates, style: TextStyle(fontSize: 13, color: t.text))),
                SizedBox(
                  width: 230,
                  child: DropdownButtonFormField<DuplicatePolicy>(
                    icon: dropdownChevron,
                    initialValue: _cfg.duplicatePolicy,
                    isExpanded: true,
                    isDense: true,
                    items: [
                      DropdownMenuItem(value: DuplicatePolicy.keepLongest, child: Text(l10n.dupLongest)),
                      DropdownMenuItem(value: DuplicatePolicy.keepMostComplete, child: Text(l10n.dupMostComplete)),
                      DropdownMenuItem(value: DuplicatePolicy.keepFirst, child: Text(l10n.dupFirst)),
                      DropdownMenuItem(value: DuplicatePolicy.keepAllRenamed, child: Text(l10n.dupAll)),
                    ],
                    onChanged: (v) => _update(_cfg.copyWith(duplicatePolicy: v)),
                  ),
                ),
              ]),
              if (r != null && r.duplicates.isNotEmpty)
                Padding(
                  padding: const EdgeInsets.only(top: 4),
                  child: Text(l10n.orgDuplicates('${r.duplicates.length}'), style: TextStyle(fontSize: 12, color: t.textMuted)),
                ),
              const SizedBox(height: 12),
              if (r != null) _MarkerCoverage(result: r),
            ]),
          ),
        ]),
      ]),
    );
  }

  /// Two ids that differ only a little (a zero, a dash): a variation of
  /// the same name rather than another sample.
  static bool _similar(String a, String b) {
    String core(String s) => s.replaceAll(RegExp(r'[^A-Z0-9]'), '');
    final x = core(a), y = core(b);
    return x == y || x.contains(y) || y.contains(x);
  }
}

class _TokenChip extends StatelessWidget {
  const _TokenChip({required this.text, required this.role, this.onTap});
  final String text;
  final _Role role;
  final VoidCallback? onTap;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final color = switch (role) {
      _Role.sample => sampleColor,
      _Role.marker => markerColor,
      _Role.none => null,
    };
    return MouseRegion(
      cursor: onTap == null ? MouseCursor.defer : SystemMouseCursors.click,
      child: GestureDetector(
        onTap: onTap,
        child: Container(
          padding: const EdgeInsets.symmetric(horizontal: 8, vertical: 4),
          decoration: BoxDecoration(
            color: color?.withValues(alpha: 0.14) ?? t.surface,
            borderRadius: BorderRadius.circular(8),
            border: Border.all(color: color ?? t.border, width: color == null ? 1 : 1.5),
          ),
          child: Text(text,
              style: TextStyle(fontFamily: monoFont, fontSize: 12.5, fontWeight: color == null ? FontWeight.w400 : FontWeight.w700, color: color ?? t.textMuted)),
        ),
      ),
    );
  }
}

class _Legend extends StatelessWidget {
  const _Legend({required this.color, required this.text});
  final Color color;
  final String text;

  @override
  Widget build(BuildContext context) => Row(mainAxisSize: MainAxisSize.min, children: [
        Container(width: 10, height: 10, decoration: BoxDecoration(color: color.withValues(alpha: 0.2), border: Border.all(color: color, width: 1.5), borderRadius: BorderRadius.circular(3))),
        const SizedBox(width: 6),
        Text(text, style: TextStyle(fontSize: 12, color: context.tokens.textMuted)),
      ]);
}

/// Sample name that can be edited (Enter or leaving the field renames).
class _RenameField extends StatefulWidget {
  const _RenameField({super.key, required this.name, required this.onSubmitted});
  final String name;
  final ValueChanged<String> onSubmitted;

  @override
  State<_RenameField> createState() => _RenameFieldState();
}

class _RenameFieldState extends State<_RenameField> {
  late final _ctrl = TextEditingController(text: widget.name);
  final _focus = FocusNode();

  @override
  void initState() {
    super.initState();
    _focus.addListener(() {
      if (!_focus.hasFocus && _ctrl.text.trim() != widget.name) widget.onSubmitted(_ctrl.text);
    });
  }

  @override
  void dispose() {
    _ctrl.dispose();
    _focus.dispose();
    super.dispose();
  }

  @override
  Widget build(BuildContext context) => TextField(
        controller: _ctrl,
        focusNode: _focus,
        style: TextStyle(fontSize: 13, color: context.tokens.text),
        decoration: const InputDecoration(
          isDense: true,
          filled: false,
          border: InputBorder.none,
          enabledBorder: InputBorder.none,
          focusedBorder: UnderlineInputBorder(),
          contentPadding: EdgeInsets.symmetric(vertical: 8),
        ),
        onSubmitted: widget.onSubmitted,
      );
}

/// How many markers are in how many samples: a bar per sample count.
class _MarkerCoverage extends StatelessWidget {
  const _MarkerCoverage({required this.result});
  final OrganizerResult result;

  @override
  Widget build(BuildContext context) {
    final l10n = context.l10n;
    final t = context.tokens;
    final n = result.samples.length;
    if (n == 0 || result.loci.isEmpty) return const SizedBox();
    final counts = List.filled(n + 1, 0);
    for (final l in result.loci) {
      counts[math.min(n, result.samplesFor(l))]++;
    }
    final peak = counts.reduce(math.max);
    final complete = counts[n];
    final half = [for (var k = (n / 2).ceil(); k <= n; k++) counts[k]].fold(0, (a, b) => a + b);
    return Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
      Text(l10n.orgCoverage, style: TextStyle(fontSize: 12.5, fontWeight: FontWeight.w700, color: t.text)),
      const SizedBox(height: 6),
      SizedBox(
        height: 70,
        child: Row(crossAxisAlignment: CrossAxisAlignment.end, children: [
          for (var k = 1; k <= n; k++)
            Expanded(
              child: Tooltip(
                message: l10n.orgCoverageBar('${counts[k]}', '$k'),
                child: Container(
                  margin: const EdgeInsets.symmetric(horizontal: 0.5),
                  height: peak == 0 ? 0 : math.max(counts[k] == 0 ? 0 : 2, 70 * counts[k] / peak).toDouble(),
                  color: k == n ? SplaceBlue.c500 : SplaceBlue.c300,
                ),
              ),
            ),
        ]),
      ),
      const SizedBox(height: 2),
      Row(children: [
        Text('1', style: TextStyle(fontSize: 10.5, color: t.textFaint)),
        const Spacer(),
        Text(l10n.orgSamplesAxis, style: TextStyle(fontSize: 10.5, color: t.textFaint)),
        const Spacer(),
        Text('$n', style: TextStyle(fontSize: 10.5, color: t.textFaint)),
      ]),
      const SizedBox(height: 6),
      Text(l10n.orgCoverageSummary('$complete', '$half'), style: TextStyle(fontSize: 12, color: t.textMuted)),
    ]);
  }
}
