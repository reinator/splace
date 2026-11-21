import 'package:flutter/material.dart';

import '../../core/delimitation.dart';
import '../../services/tool_locator.dart';
import '../../state/project_controller.dart';
import '../common.dart';
import '../fa.dart';
import '../theme.dart';
import '../widgets/citations.dart';
import '../widgets/kit.dart';

/// Step 8e — species delimitation, by published programs: mPTP (PTP and
/// mPTP models) on the ML trees of the runs, ASAP and ABGD on the
/// alignments. IQ-TREE and RAxML-NG do not delimit species; mPTP reads
/// their trees.
class DelimitationCard extends StatefulWidget {
  const DelimitationCard({super.key, this.onViewTree});

  /// Opens step 9, where each delimitation has its button on the tree.
  final VoidCallback? onViewTree;

  @override
  State<DelimitationCard> createState() => _DelimitationCardState();
}

class _DelimitationCardState extends State<DelimitationCard> {
  final Set<PtpModel> _models = {PtpModel.multi, PtpModel.single};

  /// Runs left out by the user (all runs are delimited by default).
  final Set<TreeRun> _skipped = {};
  bool _autoMinbr = true, _crop = true, _gmyc = true;

  /// ASAP and ABGD: the alignments chosen (none by default: single-locus
  /// methods, the user picks the barcode marker) and their settings.
  final Set<String> _alignments = {};
  final Set<DelimMethod> _distanceMethods = {DelimMethod.asap, DelimMethod.abgd};
  DelimDistance _distance = DelimDistance.jc69;
  double _pMin = 0.001, _pMax = 0.1, _steps = 10, _slope = 1.5;

  /// Significance levels of the likelihood-ratio tests (the programs'
  /// defaults).
  double _mptpP = mptpDefaultPvalue, _gmycP = gmycDefaultPvalue;

  Widget _toggle(String label, bool on, VoidCallback onTap) => ChipButton(
        label: label,
        icon: on ? FA.check : null,
        keepIcon: on,
        style: on ? ChipButtonStyle.accent : ChipButtonStyle.neutral,
        onPressed: onTap,
      );

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    final runs = c.delimitableRuns;
    final chosenRuns = [for (final r in runs) if (!_skipped.contains(r)) r];
    final alignments = c.delimitableAlignments;
    final chosenAlignments = [for (final a in alignments) if (_alignments.contains(a)) a];
    final hasOutgroup = c.concat.outgroups.isNotEmpty;
    bool missing(ExternalTool tool) => c.toolsChecked && c.tools[tool] == null;
    String when(TreeRun r) => r.started.toIso8601String().substring(11, 16);
    final distanceMethods = {
      for (final m in _distanceMethods)
        if (!missing(m == DelimMethod.asap ? ExternalTool.asap : ExternalTool.abgd)) m,
    };
    // ABGD computes only JC69 and K2P distances.
    if (!_distance.abgd && !_distanceMethods.contains(DelimMethod.asap)) _distance = DelimDistance.jc69;
    final gmycReady = _gmyc && !missing(ExternalTool.lsd2) && !missing(ExternalTool.gmyc);
    final mptpReady = _models.isNotEmpty && !missing(ExternalTool.mptp);
    final treeWork = chosenRuns.isNotEmpty && (mptpReady || gmycReady);
    final alignmentWork = chosenAlignments.isNotEmpty && distanceMethods.isNotEmpty;

    Widget card(String title, IconData icon, List<Widget> children) => Container(
          padding: const EdgeInsets.all(16),
          decoration: BoxDecoration(color: t.surface, borderRadius: BorderRadius.circular(12), border: Border.all(color: t.border), boxShadow: shadowSm),
          child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
            Row(children: [
              Icon(icon, size: 14, color: t.accent),
              const SizedBox(width: 8),
              Expanded(child: Text(title, style: TextStyle(fontSize: 14, fontWeight: FontWeight.w700, color: t.accentText))),
            ]),
            const SizedBox(height: 10),
            ...children,
          ]),
        );
    Widget missingPanel(ExternalTool tool) => missing(tool)
        ? Padding(
            padding: const EdgeInsets.only(bottom: 10),
            child: TonePanel(tone: PanelTone.warn, child: Text('${l10n.toolMissing(tool.label)} — ${tool.installHint}', style: const TextStyle(fontSize: 13))),
          )
        : const SizedBox.shrink();
    // Every field the same: one line of label (the details are in the
    // info dot) above a 44 px input, in a grid of equal columns.
    Widget field(String label, Widget child, {String? info}) => Column(crossAxisAlignment: CrossAxisAlignment.stretch, mainAxisSize: MainAxisSize.min, children: [
          SizedBox(
            height: 20,
            child: Row(children: [
              Flexible(
                child: Text(label, maxLines: 1, overflow: TextOverflow.ellipsis, style: TextStyle(fontSize: 13.5, fontWeight: FontWeight.w600, color: t.text)),
              ),
              if (info != null) ...[const SizedBox(width: 6), InfoDot(info)],
            ]),
          ),
          const SizedBox(height: 6),
          child,
        ]);
    Widget number(double value, ValueChanged<double> set, {double min = 0.0001, double max = 0.5, double step = 0.001, int decimals = 4}) =>
        NumberField(value: value, min: min, max: max, step: step, decimals: decimals, width: double.infinity, onChanged: (v) => setState(() => set(v)));
    Widget part(int n, String title) => Padding(
          padding: const EdgeInsets.only(bottom: 6),
          child: Row(children: [
            StepBadge('$n', size: 18),
            const SizedBox(width: 8),
            Expanded(child: Kicker(title, size: 10.5)),
          ]),
        );

    final treeCard = card(l10n.delimOnTrees, FA.diagramProject, [
      Text.rich(TextSpan(style: TextStyle(fontSize: 13, color: t.textMuted, height: 1.5), children: citationSpans(l10n.delimTreesIntro))),
      const SizedBox(height: 10),
      missingPanel(ExternalTool.mptp),
      if (_gmyc) missingPanel(ExternalTool.lsd2),
      if (_gmyc) missingPanel(ExternalTool.gmyc),
      // One sequence per species: the tests of these methods tend to fail.
      if (c.oneSequencePerSpecies)
        Padding(
          padding: const EdgeInsets.only(bottom: 12),
          child: TonePanel(tone: PanelTone.warn, child: Text(l10n.delimOnePerSpecies, style: const TextStyle(fontSize: 13, height: 1.45))),
        ),
      part(1, l10n.delimTrees),
      Wrap(spacing: 8, runSpacing: 6, children: [
        for (final r in runs)
          _toggle('${r.engine} · ${when(r)}', !_skipped.contains(r), () => setState(() => _skipped.contains(r) ? _skipped.remove(r) : _skipped.add(r))),
      ]),
      const SizedBox(height: 14),
      part(2, l10n.delimModels),
      Wrap(spacing: 8, runSpacing: 6, children: [
        for (final m in PtpModel.values)
          _toggle(m == PtpModel.multi ? l10n.delimModelMulti : l10n.delimModelSingle, _models.contains(m),
              () => setState(() => _models.contains(m) ? _models.remove(m) : _models.add(m))),
        _toggle(l10n.delimModelGmyc, _gmyc, () => setState(() => _gmyc = !_gmyc)),
      ]),
      const SizedBox(height: 14),
      part(3, l10n.delimSettings),
      AutoGrid(minItemWidth: 170, children: [
        field(l10n.delimLrtLevel('PTP / mPTP'), number(_mptpP, (v) => _mptpP = v), info: l10n.delimLrtInfoMptp),
        field(l10n.delimLrtLevel('GMYC'), number(_gmycP, (v) => _gmycP = v, step: 0.01), info: l10n.delimLrtInfoGmyc),
      ]),
      const SizedBox(height: 12),
      SwitchBox(label: l10n.delimMinbr, detail: l10n.delimMinbrInfo, value: _autoMinbr, onChanged: (v) => setState(() => _autoMinbr = v)),
      if (hasOutgroup) ...[
        const SizedBox(height: 8),
        SwitchBox(label: l10n.delimCrop, detail: l10n.delimCropInfo, value: _crop, onChanged: (v) => setState(() => _crop = v)),
      ],
    ]);

    final alignmentCard = card(l10n.delimOnAlignments, FA.alignLeft, [
      Text.rich(TextSpan(style: TextStyle(fontSize: 13, color: t.textMuted, height: 1.5), children: citationSpans(l10n.delimAlignmentsIntro))),
      const SizedBox(height: 10),
      missingPanel(ExternalTool.asap),
      missingPanel(ExternalTool.abgd),
      part(1, l10n.delimAlignments),
      if (alignments.isEmpty)
        Text(l10n.delimNoAlignments, style: TextStyle(fontSize: 13, color: t.textFaint))
      else
        Wrap(spacing: 8, runSpacing: 6, children: [
          for (final a in alignments)
            _toggle(a == ProjectController.supermatrixAlignment ? l10n.delimSupermatrix : a, _alignments.contains(a),
                () => setState(() => _alignments.contains(a) ? _alignments.remove(a) : _alignments.add(a))),
        ]),
      const SizedBox(height: 14),
      part(2, l10n.delimModels),
      Wrap(spacing: 8, runSpacing: 6, children: [
        for (final m in [DelimMethod.asap, DelimMethod.abgd])
          _toggle(m.label, _distanceMethods.contains(m), () => setState(() => _distanceMethods.contains(m) ? _distanceMethods.remove(m) : _distanceMethods.add(m))),
      ]),
      const SizedBox(height: 14),
      part(3, l10n.delimSettings),
      AutoGrid(minItemWidth: 170, children: [
        field(
          l10n.delimDistance,
          SplaceSelect<DelimDistance>(
            value: _distance,
            options: [
              for (final d in DelimDistance.values)
                if (d.abgd || _distanceMethods.contains(DelimMethod.asap)) SelectOption(d, d.label),
            ],
            onChanged: (v) => setState(() => _distance = v),
          ),
          info: l10n.delimDistanceInfo,
        ),
        if (_distanceMethods.contains(DelimMethod.abgd)) ...[
          field(l10n.delimPMin, number(_pMin, (v) => _pMin = v), info: l10n.delimPriorInfo),
          field(l10n.delimPMax, number(_pMax, (v) => _pMax = v, step: 0.01), info: l10n.delimPriorInfo),
          field(l10n.delimSteps, number(_steps, (v) => _steps = v, min: 2, max: 50, step: 1, decimals: 0)),
          field(l10n.delimSlope, number(_slope, (v) => _slope = v, min: 0.1, max: 10, step: 0.1, decimals: 2), info: l10n.delimSlopeInfo),
        ],
      ]),
    ]);

    return Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
      Text.rich(TextSpan(style: TextStyle(fontSize: 13, color: t.textMuted, height: 1.5), children: citationSpans(l10n.delimIntro))),
      const SizedBox(height: 12),
      LayoutBuilder(builder: (context, box) {
        if (box.maxWidth < 900) {
          return Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [treeCard, const SizedBox(height: 12), alignmentCard]);
        }
        return EqualHeightRow(children: [treeCard, alignmentCard]);
      }),
      const SizedBox(height: 14),
      Wrap(spacing: 10, runSpacing: 8, crossAxisAlignment: WrapCrossAlignment.center, children: [
        PrimaryButton(
          label: l10n.delimRun,
          icon: FA.play,
          onPressed: c.busy || (!treeWork && !alignmentWork)
              ? null
              : () => c.runDelimitation(
                    runs: treeWork ? chosenRuns : const [],
                    models: mptpReady ? _models : const {},
                    gmyc: gmycReady,
                    mptpPvalue: _mptpP,
                    gmycPvalue: _gmycP,
                    autoMinbr: _autoMinbr,
                    cropOutgroup: _crop && hasOutgroup,
                    alignments: alignmentWork ? chosenAlignments : const [],
                    distanceMethods: distanceMethods,
                    distance: _distance,
                    abgdPMin: _pMin,
                    abgdPMax: _pMax,
                    abgdSteps: _steps.round(),
                    abgdSlope: _slope,
                  ),
        ),
        for (final tool in [ExternalTool.mptp, ExternalTool.lsd2, ExternalTool.gmyc, ExternalTool.asap, ExternalTool.abgd])
          if (c.tools[tool] case final r?) TagChip('${tool.label}${r.version == null ? '' : ' v${r.version}'}', color: t.textMuted),
      ]),
      const SizedBox(height: 12),
      if (c.delimitations.isEmpty)
        Text(l10n.delimNone, style: TextStyle(color: t.textFaint))
      else
        for (final d in c.delimitations) _DelimitationRow(key: ValueKey(d.id), d: d, onViewTree: widget.onViewTree),
    ]);
  }
}

class _DelimitationRow extends StatefulWidget {
  const _DelimitationRow({super.key, required this.d, this.onViewTree});
  final DelimitationRun d;
  final VoidCallback? onViewTree;

  @override
  State<_DelimitationRow> createState() => _DelimitationRowState();
}

class _DelimitationRowState extends State<_DelimitationRow> {
  bool _open = false;

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.read(context);
    final l10n = context.l10n;
    final t = context.tokens;
    final d = widget.d;
    final r = d.result;
    String fmt(double v) => v.abs() >= 0.001 || v == 0 ? v.toStringAsFixed(4).replaceFirst(RegExp(r'\.?0+$'), '') : v.toStringAsExponential(2);

    /// "12 species · ASAP score 2.5" / "12 species · P = 0.0129 · initial".
    String partitionLabel(int i) {
      final sp = d.partitions[i];
      final parts = <String>[l10n.delimSpecies('${sp.species.length}')];
      if (d.method == DelimMethod.asap) {
        parts.add('#${i + 1}');
        if (sp.score != null && sp.score!.isNotEmpty) parts.add(l10n.delimAsapScore(sp.score!));
      } else if (i < d.partitionNotes.length) {
        final note = d.partitionNotes[i].split('|');
        final prior = double.tryParse(note.last);
        if (prior != null) parts.add('P = ${fmt(prior)}');
        parts.add(note.first == 'rec' ? l10n.delimRecursive : l10n.delimInitial);
      }
      return parts.join(' · ');
    }

    return Container(
      margin: const EdgeInsets.only(bottom: 8),
      decoration: BoxDecoration(color: t.surface, borderRadius: BorderRadius.circular(12), border: Border.all(color: t.border)),
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        InkWell(
          onTap: () => setState(() => _open = !_open),
          child: Padding(
            padding: const EdgeInsets.all(12),
            child: Row(children: [
              StatusDot(state: stepState(d.status)),
              const SizedBox(width: 10),
              Text('${d.label} · ${d.time}', style: TextStyle(fontWeight: FontWeight.w600, fontSize: 14, color: t.text)),
              const SizedBox(width: 10),
              Expanded(
                child: Wrap(spacing: 6, runSpacing: 4, children: [
                  if (r != null) ...[
                    TagChip(l10n.delimSpecies('${r.species.length}'), icon: FA.users, color: r.collapsed ? SplaceTokens.warn : t.accentText),
                    if (r.singletons > 0 && !r.collapsed) TagChip(l10n.delimSingletons('${r.singletons}'), color: t.textMuted),
                    if (d.method.onTree) ...[
                      if (r.collapsed)
                        TagChip(l10n.delimLrtFailed, icon: FA.triangleExclamation, color: SplaceTokens.warn)
                      else if (r.pvalue != null)
                        TagChip('LRT p = ${fmt(r.pvalue!)}${r.alpha == null ? '' : ' < ${fmt(r.alpha!)}'}', color: t.textMuted),
                      if (d.method != DelimMethod.gmyc) TagChip('--minbr ${d.minbrAuto == null ? '0.0001' : fmt(d.minbrAuto!)}', color: t.textMuted),
                    ],
                    if (d.distance != null) TagChip(d.distance!.label, color: t.textMuted),
                  ],
                  if (d.status == StepStatus.error) TagChip(d.error ?? l10n.statusError, icon: FA.triangleExclamation, color: SplaceTokens.danger),
                ]),
              ),
              if (r != null) ...[
                ChipButton(
                  label: l10n.delimViewOnTree,
                  style: ChipButtonStyle.accent,
                  onPressed: () {
                    if (d.run != null) c.viewTree(c.mlTreeOf(d.run!));
                    widget.onViewTree?.call();
                  },
                ),
                const SizedBox(width: 8),
              ],
              IconBtn(icon: FA.folderOpen, tooltip: l10n.openFolder, size: ChipButton.height, bordered: true, onPressed: () => c.openFolder(d.directory)),
              const SizedBox(width: 8),
              Icon(_open ? FA.chevronUp : FA.chevronDown, size: 12, color: t.textMuted),
            ]),
          ),
        ),
        // The test failed: why every tip is one species, and what to try.
        if (r != null && r.collapsed)
          Padding(
            padding: const EdgeInsets.fromLTRB(12, 0, 12, 12),
            child: TonePanel(
              tone: PanelTone.warn,
              child: Text(
                [
                  l10n.delimCollapsed(d.method.label, r.pvalue == null ? '?' : fmt(r.pvalue!), r.alpha == null ? '?' : fmt(r.alpha!), '${r.species.fold<int>(0, (a, s) => a + s.length)}'),
                  if (r.bestSpecies != null && r.bestSpecies! > 1) l10n.delimCollapsedBest('${r.bestSpecies}'),
                  c.oneSequencePerSpecies ? l10n.delimCollapsedOnePerSpecies : l10n.delimCollapsedTry,
                ].join(' '),
                style: const TextStyle(fontSize: 13, height: 1.45),
              ),
            ),
          ),
        // ASAP and ABGD propose several partitions: the one shown on the
        // tree and in the methods text is chosen here.
        if (d.partitions.length > 1)
          Padding(
            padding: const EdgeInsets.fromLTRB(12, 0, 12, 12),
            child: Row(children: [
              Text(l10n.delimPartition, style: TextStyle(fontSize: 13, fontWeight: FontWeight.w600, color: t.text)),
              const SizedBox(width: 6),
              InfoDot(d.method == DelimMethod.asap ? l10n.delimPartitionAsapInfo : l10n.delimPartitionAbgdInfo),
              const SizedBox(width: 10),
              SizedBox(
                width: 380,
                child: SplaceSelect<int>(
                  value: d.chosen,
                  options: [for (var i = 0; i < d.partitions.length; i++) SelectOption(i, partitionLabel(i))],
                  onChanged: (i) => c.choosePartition(d, i),
                ),
              ),
            ]),
          ),
        if (_open)
          Padding(
            padding: const EdgeInsets.fromLTRB(12, 0, 12, 12),
            child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
              if (d.minbrCommand != null) ...[
                PreviewBox(label: d.method == DelimMethod.gmyc ? l10n.delimLsd2Command : l10n.delimMinbrCommand, text: d.minbrCommand!, maxHeight: 90),
                const SizedBox(height: 8),
              ],
              if (d.command != null) ...[
                PreviewBox(label: l10n.commandPreview, text: d.command!, maxHeight: 90),
                const SizedBox(height: 8),
              ],
              if (r != null)
                PreviewBox(
                  label: l10n.delimSpeciesList,
                  text: [for (final (i, s) in r.species.indexed) '${i + 1}. ${s.join(', ')}'].join('\n'),
                  maxHeight: 220,
                ),
            ]),
          ),
      ]),
    );
  }
}
