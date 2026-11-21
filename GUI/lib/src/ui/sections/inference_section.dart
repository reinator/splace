import 'dart:io';
import 'dart:math' as math;

import 'package:flutter/material.dart';
import 'package:flutter/services.dart';
import 'package:path/path.dart' as p;

import '../../core/tool_commands.dart';
import '../../services/tool_locator.dart';
import '../../state/project_controller.dart';
import '../common.dart';
import '../fa.dart';
import '../theme.dart';
import '../widgets/kit.dart';
import 'delimitation_card.dart';

const _iqModels = ['TEST', 'MFP', 'MFP+MERGE', 'GTR+G', 'GTR+G+I', 'HKY+G', 'custom'];
const _bootstraps = ['-B 1000', '-B 5000', '-b 100', ''];

/// Step 8 — Phylogenetic Inference with IQ-TREE3 (as in the HTML app) and
/// RAxML-NG, per-gene trees and the results ZIP.
class InferenceSection extends StatefulWidget {
  const InferenceSection({super.key, required this.active, this.onViewTree});
  final bool active;

  /// Opens step 9 (tree viewer).
  final VoidCallback? onViewTree;

  @override
  State<InferenceSection> createState() => _InferenceSectionState();
}

class _InferenceSectionState extends State<InferenceSection> {
  String _iqModel = 'MFP+MERGE';
  final _custom = TextEditingController(text: 'GTR+R4');
  String _bootstrap = '-B 1000';
  final _extra = TextEditingController();
  int _alrt = 1000;
  bool _bnni = true;
  double _nstop = 100, _pers = 0.5, _runs = 1;
  final _iqSeed = TextEditingController();
  late final _threads = TextEditingController(text: '${ProjectScope.read(context).threads}');
  bool _perGene = false;
  RaxmlModels _rxModels = RaxmlModels.moose;
  final _rxFixed = TextEditingController(text: 'GTR+G');
  String _criterion = 'BIC';
  bool _rxFast = false;
  final _rxBootstraps = TextEditingController(text: 'autoMRE{1000}');
  Set<RaxmlSupport> _metrics = {RaxmlSupport.fbp, RaxmlSupport.tbe};
  final _seed = TextEditingController(text: '12345');

  @override
  void initState() {
    super.initState();
    // A project opened from a file: the form shows the settings of its
    // last runs.
    final runs = ProjectScope.read(context).treeRuns;
    final iq = runs.where((r) => r.iqtree != null && !r.iqtree!.raxmlModelsOnly).firstOrNull;
    if (iq != null) {
      final o = iq.iqtree!;
      if (_iqModels.contains(o.model)) {
        _iqModel = o.model;
      } else {
        _iqModel = 'custom';
        _custom.text = o.model;
      }
      final b = o.classicBootstrap > 0 ? '-b ${o.classicBootstrap}' : (o.ufboot > 0 ? '-B ${o.ufboot}' : '');
      if (_bootstraps.contains(b)) _bootstrap = b;
      _alrt = o.alrt;
      _bnni = o.bnni;
      _nstop = (o.nstop ?? 100).toDouble();
      _pers = o.perturbation ?? 0.5;
      _runs = o.runs.toDouble();
      _iqSeed.text = o.seed?.toString() ?? '';
      _extra.text = o.extra;
      _perGene = iq.geneTrees.isNotEmpty;
    }
    final rx = runs.where((r) => r.raxml != null).firstOrNull;
    if (rx != null) {
      final o = rx.raxml!;
      _rxBootstraps.text = o.bootstraps;
      _metrics = {for (final m in RaxmlSupport.values) if (o.supportMetrics.split(',').contains(m.flag)) m};
      if (_metrics.isEmpty) _metrics = {RaxmlSupport.fbp, RaxmlSupport.tbe};
      _rxFast = o.fast;
      _criterion = o.criterion;
      _seed.text = '${o.seed}';
      _rxModels = rx.models;
      if (rx.fixedModel != null) _rxFixed.text = rx.fixedModel!;
    }
  }

  @override
  void dispose() {
    for (final x in [_custom, _extra, _threads, _rxFixed, _rxBootstraps, _seed, _iqSeed]) {
      x.dispose();
    }
    super.dispose();
  }

  IqtreeOptions _iqOptions(ProjectController c) {
    final b = _bootstrap.split(' ');
    return IqtreeOptions(
      model: _iqModel == 'custom' ? (_custom.text.trim().isEmpty ? 'GTR+G' : _custom.text.trim()) : _iqModel,
      ufboot: b.first == '-B' ? int.parse(b.last) : 0,
      classicBootstrap: b.first == '-b' ? int.parse(b.last) : 0,
      alrt: _alrt,
      bnni: _bnni,
      nstop: _nstop == 100 ? null : _nstop.round(),
      perturbation: _pers == 0.5 ? null : _pers,
      runs: _runs.round(),
      seed: int.tryParse(_iqSeed.text.trim()),
      threads: _threads.text.trim().isEmpty ? 'AUTO' : _threads.text.trim(),
      maxThreads: c.threads,
      outgroup: c.concat.outgroups.where((o) => c.matrix?.taxa.contains(o) ?? false).toList(),
      extra: _extra.text,
    );
  }

  RaxmlOptions _rxOptions(ProjectController c, List<String> outgroups) => RaxmlOptions(
        bootstraps: _rxBootstraps.text.trim().isEmpty ? 'autoMRE{1000}' : _rxBootstraps.text.trim(),
        supportMetrics: [for (final m in RaxmlSupport.values) if (_metrics.contains(m)) m.flag].join(','),
        fast: _rxFast,
        criterion: _criterion,
        maxThreads: c.threads,
        seed: int.tryParse(_seed.text) ?? 12345,
        outgroup: outgroups,
      );

  /// FBP and RBS exclude each other (RAxML-NG refuses both).
  void _toggleMetric(RaxmlSupport m, bool on) => setState(() {
        if (on) {
          _metrics.add(m);
          if (m == RaxmlSupport.fbp) _metrics.remove(RaxmlSupport.rbs);
          if (m == RaxmlSupport.rbs) _metrics.remove(RaxmlSupport.fbp);
        } else if (_metrics.length > 1) {
          _metrics.remove(m);
        }
        _metrics = {..._metrics};
      });

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    final m = c.matrix!;
    final outgroups = c.concat.outgroups.where(m.taxa.contains).toList();
    final iq = _iqOptions(c);
    final preview = c.iqtreeCommandPreview(iq);
    final rx = _rxOptions(c, outgroups);
    final needsReplicates = rx.needsReplicates;
    final rxVersion = c.tools[ExternalTool.raxmlng]?.version;
    final rxOld = (int.tryParse(rxVersion?.split('.').first ?? '') ?? 2) < 2;
    final rxPreview = c.raxmlCommandPreview(rx, models: _rxModels, fixedModel: _rxFixed.text.trim());

    Widget missing(ExternalTool tool) => c.toolsChecked && c.tools[tool] == null
        ? Padding(
            padding: const EdgeInsets.only(bottom: 12),
            child: TonePanel(
              tone: PanelTone.warn,
              child: Text('${l10n.toolMissing(tool.label)} — ${tool.installHint}', style: const TextStyle(fontSize: 13)),
            ),
          )
        : const SizedBox.shrink();

    String modelLabel(String v) => switch (v) {
          'TEST' => l10n.modelTest,
          'MFP' => l10n.modelMfp,
          'MFP+MERGE' => l10n.modelMfpMerge,
          'custom' => l10n.modelCustom,
          _ => v,
        };
    String bootLabel(String v) => switch (v) {
          '-B 1000' => l10n.bootUf1000,
          '-B 5000' => l10n.bootUf5000,
          '-b 100' => l10n.bootClassic,
          _ => l10n.bootNone,
        };

    return StepSection(
      badge: '7',
      title: l10n.step8,
      description: l10n.step8Intro,
      active: widget.active,
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        TonePanel(
          tone: PanelTone.info,
          child: Row(children: [
            const Icon(FA.sitemap, size: 13),
            const SizedBox(width: 8),
            Expanded(
              child: Text(
                '${l10n.outgroup}: ${outgroups.isEmpty ? l10n.noneSelected : outgroups.join(', ')} · ${l10n.outgroupFromStep7}',
                style: const TextStyle(fontSize: 13),
              ),
            ),
            Text(l10n.matrixStats('${m.taxa.length}', '${m.length}', '${m.partitions.length}', pct(m.missingFraction)),
                style: const TextStyle(fontSize: 12, fontWeight: FontWeight.w600)),
          ]),
        ),
        const SizedBox(height: 16),
        LayoutBuilder(builder: (context, box) {
          // Side by side on wide pages, with the same height (run buttons
          // aligned at the bottom); one above the other otherwise.
          final wide = box.maxWidth >= 2 * 480 + 16;
          // The command previews sit in the footers, right above the run
          // buttons, with the same number of lines: side by side they start
          // and end at the same height.
          final lines = math.max('\n'.allMatches(preview).length, '\n'.allMatches(rxPreview).length) + 1;
          Widget bottom(String command, Widget button) => Column(mainAxisSize: MainAxisSize.min, crossAxisAlignment: CrossAxisAlignment.stretch, children: [
                const SizedBox(height: 12),
                _CommandBox(text: command, lines: wide ? lines : null),
                const SizedBox(height: 16),
                button,
              ]);
          final iqCard = OptionCard(
            badge: '7a',
            title: 'IQ-TREE3 ${c.tools[ExternalTool.iqtree]?.version ?? ''}',
            info: l10n.iqtreeInfo,
            footer: bottom(preview, PrimaryButton(
              label: l10n.runIqtree,
              icon: FA.play,
              expand: true,
              onPressed: c.busy ? null : () => c.runIqtree(_iqOptions(c), perGene: _perGene),
            )),
            child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
              missing(ExternalTool.iqtree),
              AutoGrid(minItemWidth: 300, children: [
                Field(
                  label: l10n.substitutionModel,
                  child: DropdownButtonFormField<String>(
                    icon: dropdownChevron,
                    initialValue: _iqModel,
                    isExpanded: true,
                    items: [for (final v in _iqModels) DropdownMenuItem(value: v, child: Text(modelLabel(v)))],
                    onChanged: (v) => setState(() => _iqModel = v!),
                  ),
                ),
                Field(
                  label: l10n.bootstrapSupport,
                  child: DropdownButtonFormField<String>(
                    icon: dropdownChevron,
                    initialValue: _bootstrap,
                    isExpanded: true,
                    items: [for (final v in _bootstraps) DropdownMenuItem(value: v, child: Text(bootLabel(v)))],
                    onChanged: (v) => setState(() => _bootstrap = v!),
                  ),
                ),
                if (_iqModel == 'custom')
                  Field(
                    label: l10n.modelCustom,
                    child: TextField(controller: _custom, decoration: const InputDecoration(hintText: 'e.g. GTR+R4'), onChanged: (_) => setState(() {})),
                  ),
                Field(
                  label: l10n.iqAlrt,
                  info: l10n.iqAlrtInfo,
                  child: DropdownButtonFormField<int>(
                    icon: dropdownChevron,
                    initialValue: _alrt,
                    isExpanded: true,
                    items: [
                      DropdownMenuItem(value: 1000, child: Text(l10n.iqAlrt1000)),
                      DropdownMenuItem(value: 0, child: Text(l10n.bootNone)),
                    ],
                    onChanged: (v) => setState(() => _alrt = v!),
                  ),
                ),
                Field(
                  label: l10n.threadsLabel,
                  child: TextField(controller: _threads, decoration: const InputDecoration(hintText: 'AUTO'), onChanged: (_) => setState(() {})),
                ),
              ]),
              const SizedBox(height: 14),
              Row(children: [
                Kicker(l10n.iqSearchHeading),
                const SizedBox(width: 6),
                InfoDot(l10n.iqSearchInfo),
              ]),
              const SizedBox(height: 8),
              AutoGrid(minItemWidth: 220, children: [
                Field(
                  label: l10n.iqRuns,
                  info: l10n.iqRunsInfo,
                  child: NumberField(value: _runs, min: 1, max: 20, width: double.infinity, onChanged: (v) => setState(() => _runs = v)),
                ),
                Field(
                  label: l10n.iqNstop,
                  info: l10n.iqNstopInfo,
                  child: NumberField(value: _nstop, min: 10, max: 2000, step: 50, width: double.infinity, onChanged: (v) => setState(() => _nstop = v)),
                ),
                Field(
                  label: l10n.iqPers,
                  info: l10n.iqPersInfo,
                  child: NumberField(value: _pers, min: 0.1, max: 1, step: 0.1, decimals: 1, width: double.infinity, onChanged: (v) => setState(() => _pers = v)),
                ),
                Field(
                  label: l10n.seed,
                  info: l10n.iqSeedInfo,
                  child: TextField(controller: _iqSeed, keyboardType: TextInputType.number, decoration: InputDecoration(hintText: l10n.iqSeedRandom), onChanged: (_) => setState(() {})),
                ),
              ]),
              if (_bootstrap.startsWith('-B')) ...[
                const SizedBox(height: 12),
                CheckCard(value: _bnni, onChanged: (v) => setState(() => _bnni = v), title: l10n.iqBnni, subtitle: l10n.iqBnniBody),
              ],
              const SizedBox(height: 12),
              Field(
                label: '${l10n.extraFlags} ${l10n.spaceSeparated}',
                child: TextField(
                  controller: _extra,
                  style: const TextStyle(fontFamily: monoFont, fontSize: 13),
                  decoration: const InputDecoration(hintText: 'e.g. -mrate G,R,I+G --mset GTR,HKY'),
                  onChanged: (_) => setState(() {}),
                ),
              ),
              const SizedBox(height: 12),
              CheckCard(value: _perGene, onChanged: (v) => setState(() => _perGene = v), title: l10n.perGeneTrees, subtitle: l10n.perGeneTreesBody),
            ]),
          );
          final rxCard = OptionCard(
            badge: '7b',
            title: 'RAxML-NG ${c.tools[ExternalTool.raxmlng]?.version ?? ''}',
            info: l10n.raxmlInfo,
            footer: bottom(rxPreview, PrimaryButton(
              label: l10n.runRaxml,
              icon: FA.play,
              expand: true,
              onPressed: c.busy
                  ? null
                  : () => c.runRaxml(
                        _rxOptions(c, outgroups),
                        models: _rxModels,
                        fixedModel: _rxModels == RaxmlModels.fixed && _rxFixed.text.trim().isNotEmpty ? _rxFixed.text.trim() : null,
                      ),
            )),
            child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
              missing(ExternalTool.raxmlng),
              if (rxOld)
                Padding(
                  padding: const EdgeInsets.only(bottom: 12),
                  child: TonePanel(tone: PanelTone.warn, child: Text(l10n.raxmlOld(rxVersion ?? ''), style: const TextStyle(fontSize: 13))),
                ),
              Kicker(l10n.raxmlModelsHeading),
              const SizedBox(height: 8),
              RadioCard<RaxmlModels>(
                value: RaxmlModels.moose,
                group: _rxModels,
                onChanged: (v) => setState(() => _rxModels = v),
                title: l10n.raxmlMoose,
                subtitle: l10n.raxmlMooseBody,
              ),
              const SizedBox(height: 8),
              RadioCard<RaxmlModels>(
                value: RaxmlModels.modelFinder,
                group: _rxModels,
                onChanged: (v) => setState(() => _rxModels = v),
                title: l10n.raxmlModelFinder,
                subtitle: l10n.raxmlModelFinderBody,
              ),
              const SizedBox(height: 8),
              RadioCard<RaxmlModels>(
                value: RaxmlModels.fixed,
                group: _rxModels,
                onChanged: (v) => setState(() => _rxModels = v),
                title: l10n.raxmlFixedModel,
                subtitle: l10n.raxmlFixedBody,
              ),
              const SizedBox(height: 12),
              AutoGrid(minItemWidth: 180, children: [
                if (_rxModels == RaxmlModels.moose)
                  Field(
                    label: l10n.raxmlCriterion,
                    child: DropdownButtonFormField<String>(
                      icon: dropdownChevron,
                      initialValue: _criterion,
                      isExpanded: true,
                      items: [for (final v in const ['BIC', 'AIC', 'AICc']) DropdownMenuItem(value: v, child: Text(v))],
                      onChanged: (v) => setState(() => _criterion = v!),
                    ),
                  ),
                if (_rxModels == RaxmlModels.fixed)
                  Field(
                    label: l10n.substitutionModel,
                    child: TextField(controller: _rxFixed, decoration: const InputDecoration(hintText: 'GTR+G'), onChanged: (_) => setState(() {})),
                  ),
                Field(
                  label: l10n.raxmlSearch,
                  child: DropdownButtonFormField<bool>(
                    icon: dropdownChevron,
                    initialValue: _rxFast,
                    isExpanded: true,
                    items: [
                      DropdownMenuItem(value: false, child: Text(l10n.raxmlSearchAdaptive)),
                      DropdownMenuItem(value: true, child: Text(l10n.raxmlSearchFast)),
                    ],
                    onChanged: (v) => setState(() => _rxFast = v!),
                  ),
                ),
                if (needsReplicates) Field(label: l10n.bootstraps, child: TextField(controller: _rxBootstraps, onChanged: (_) => setState(() {}))),
                Field(label: l10n.seed, child: TextField(controller: _seed, keyboardType: TextInputType.number, onChanged: (_) => setState(() {}))),
              ]),
              const SizedBox(height: 12),
              Row(children: [
                Kicker(l10n.supportMetrics),
                const SizedBox(width: 6),
                InfoDot(l10n.raxmlSupportInfo),
              ]),
              const SizedBox(height: 8),
              AutoGrid(minItemWidth: 130, spacing: 8, children: [
                for (final m in RaxmlSupport.values)
                  MarkerChip(
                    name: m.label,
                    active: _metrics.contains(m),
                    onChanged: (v) => _toggleMetric(m, v),
                  ),
              ]),
            ]),
          );
          if (!wide) {
            return Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
              SubAnchor(id: '7a', child: iqCard),
              const SizedBox(height: 16),
              SubAnchor(id: '7b', child: rxCard),
            ]);
          }
          // Side by side: both letters lead to the row.
          return SubAnchor(id: '7a', child: EqualHeightRow(children: [iqCard, rxCard]));
        }),
        const SizedBox(height: 16),
        SubAnchor(
          id: '7c',
          child: _BothCard(
          threads: c.threads,
          busy: c.busy || c.tools[ExternalTool.iqtree] == null || c.tools[ExternalTool.raxmlng] == null,
          onRun: () => c.runBoth(
            _iqOptions(c),
            _rxOptions(c, outgroups),
            models: _rxModels,
            fixedModel: _rxModels == RaxmlModels.fixed && _rxFixed.text.trim().isNotEmpty ? _rxFixed.text.trim() : null,
          ),
        ),
        ),
        const SizedBox(height: 22),
        SubAnchor(
          id: '7d',
          child: SubHeading(
            id: '7d',
            title: l10n.runs,
            trailing: ChipButton(label: l10n.openFolder, icon: FA.folderOpen, onPressed: () => c.openFolder(c.workDir)),
          ),
        ),
        if (c.treeRuns.isEmpty) Text(l10n.noRuns, style: TextStyle(color: t.textFaint)) else for (final r in c.treeRuns) _RunCard(run: r, onViewTree: widget.onViewTree),
        // Species delimitation on the ML trees, once there are trees.
        if (c.delimitableRuns.isNotEmpty || c.delimitableAlignments.isNotEmpty) ...[
          const SizedBox(height: 22),
          SubAnchor(id: '7e', child: SubHeading(id: '7e', title: l10n.delimTitle)),
          DelimitationCard(onViewTree: widget.onViewTree),
        ],
      ]),
    );
  }
}

/// Running both programs at once: the threads split between them, drawn as
/// one square per thread.
class _BothCard extends StatelessWidget {
  const _BothCard({required this.threads, required this.busy, required this.onRun});
  final int threads;
  final bool busy;
  final VoidCallback onRun;

  @override
  Widget build(BuildContext context) {
    final l10n = context.l10n;
    final t = context.tokens;
    final total = threads < 2 ? 2 : threads;
    final iq = total ~/ 2, rx = total - total ~/ 2;
    const iqColor = SplacePalette.blue600, rxColor = SplacePalette.red500;
    Widget legend(Color c, String text) => Row(mainAxisSize: MainAxisSize.min, children: [
          Container(width: 12, height: 12, decoration: BoxDecoration(color: c, borderRadius: BorderRadius.circular(3))),
          const SizedBox(width: 6),
          Text(text, style: TextStyle(fontSize: 12.5, fontWeight: FontWeight.w600, color: t.text)),
        ]);
    return OptionCard(
      badge: '7c',
      title: l10n.bothTitle,
      info: l10n.bothInfo,
      child: LayoutBuilder(builder: (context, box) {
        final text = Column(crossAxisAlignment: CrossAxisAlignment.start, mainAxisSize: MainAxisSize.min, children: [
          Text(l10n.bothBody('$iq', '$rx'), style: TextStyle(fontSize: 13, color: t.textMuted, height: 1.45)),
          const SizedBox(height: 12),
          PrimaryButton(label: l10n.runBoth, icon: FA.play, onPressed: busy ? null : onRun),
        ]);
        final picture = Column(crossAxisAlignment: CrossAxisAlignment.start, mainAxisSize: MainAxisSize.min, children: [
          Kicker(l10n.bothThreads('$total')),
          const SizedBox(height: 6),
          Wrap(spacing: 3, runSpacing: 3, children: [
            for (var i = 0; i < total; i++)
              Container(
                width: 18,
                height: 18,
                alignment: Alignment.center,
                decoration: BoxDecoration(color: i < iq ? iqColor : rxColor, borderRadius: BorderRadius.circular(4)),
                child: Text('${i + 1}', style: const TextStyle(fontSize: 8.5, fontWeight: FontWeight.w700, color: Colors.white)),
              ),
          ]),
          const SizedBox(height: 8),
          Wrap(spacing: 14, runSpacing: 6, children: [legend(iqColor, 'IQ-TREE · $iq'), legend(rxColor, 'RAxML-NG · $rx')]),
        ]);
        if (box.maxWidth < 720) return Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [text, const SizedBox(height: 14), picture]);
        return Row(crossAxisAlignment: CrossAxisAlignment.start, children: [
          Expanded(flex: 5, child: text),
          const SizedBox(width: 24),
          Expanded(flex: 4, child: picture),
        ]);
      }),
    );
  }
}

/// Dark box with a command line and a copy button.
class _CommandBox extends StatelessWidget {
  const _CommandBox({required this.text, this.lines});
  final String text;

  /// Height in lines (the longer of two previews shown side by side);
  /// long lines then scroll sideways instead of wrapping.
  final int? lines;

  static const _style = TextStyle(fontFamily: monoFont, fontSize: 12, height: 1.45, color: Color(0xFFE2E8F0));

  @override
  Widget build(BuildContext context) {
    final l10n = context.l10n;
    return Container(
      padding: const EdgeInsets.all(12),
      decoration: BoxDecoration(color: const Color(0xFF0F172A), borderRadius: BorderRadius.circular(10)),
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        Row(children: [
          const Icon(FA.terminal, size: 11, color: Color(0xFF94A3B8)),
          const SizedBox(width: 6),
          Expanded(child: Kicker(l10n.commandPreview, color: const Color(0xFF94A3B8), size: 11)),
          IconButton(
            tooltip: l10n.logCopy,
            visualDensity: VisualDensity.compact,
            onPressed: () => Clipboard.setData(ClipboardData(text: text)),
            icon: const Icon(FA.copy, size: 12, color: Color(0xFF94A3B8)),
          ),
        ]),
        if (lines == null)
          SelectableText(text, style: _style)
        else
          SizedBox(
            height: lines! * 12 * 1.45 + 14,
            child: Scrollbar(
              child: SingleChildScrollView(
                scrollDirection: Axis.horizontal,
                padding: const EdgeInsets.only(bottom: 14),
                child: SelectableText(text, style: _style),
              ),
            ),
          ),
      ]),
    );
  }
}

class _RunCard extends StatefulWidget {
  const _RunCard({required this.run, this.onViewTree});
  final VoidCallback? onViewTree;
  final TreeRun run;

  @override
  State<_RunCard> createState() => _RunCardState();
}

class _RunCardState extends State<_RunCard> {
  bool _open = true;

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.read(context);
    final t = context.tokens;
    final l10n = context.l10n;
    final run = widget.run;
    final newick = _open && run.treeFile != null && File(run.treeFile!).existsSync() ? File(run.treeFile!).readAsStringSync().trim() : null;
    final files = Directory(run.directory).existsSync()
        ? (Directory(run.directory).listSync(recursive: true).whereType<File>().map((f) => p.relative(f.path, from: run.directory)).toList()..sort())
        : const <String>[];
    return Container(
      margin: const EdgeInsets.only(bottom: 10),
      decoration: BoxDecoration(color: t.surface, borderRadius: BorderRadius.circular(12), border: Border.all(color: t.border), boxShadow: shadowSm),
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        InkWell(
          onTap: () => setState(() => _open = !_open),
          child: Padding(
            padding: const EdgeInsets.all(14),
            child: Row(children: [
              StatusDot(state: stepState(run.status)),
              const SizedBox(width: 10),
              Text(run.engine, style: TextStyle(fontWeight: FontWeight.w600, fontSize: 14, color: t.text)),
              const SizedBox(width: 10),
              Expanded(
                child: Text(
                  [run.started.toIso8601String().substring(0, 19).replaceFirst('T', ' '), ...run.notes.take(2)].join(' · '),
                  overflow: TextOverflow.ellipsis,
                  style: TextStyle(fontSize: 13, color: t.textMuted),
                ),
              ),
              if (run.status == StepStatus.done) TagChip(l10n.analysisComplete, icon: FA.check, color: SplaceTokens.success),
              if (run.status == StepStatus.done && run.treeFile != null) ...[
                const SizedBox(width: 8),
                ChipButton(
                  label: l10n.viewTree,
                  style: ChipButtonStyle.accent,
                  onPressed: () {
                    c.viewTree(run.outputs['consensus'] ?? run.treeFile);
                    widget.onViewTree?.call();
                  },
                ),
              ],
              const SizedBox(width: 8),
              IconBtn(icon: FA.folderOpen, tooltip: l10n.openFolder, size: ChipButton.height, bordered: true, onPressed: () => c.openFolder(run.directory)),
              const SizedBox(width: 8),
              Icon(_open ? FA.chevronUp : FA.chevronDown, size: 12, color: t.textMuted),
            ]),
          ),
        ),
        if (_open)
          Padding(
            padding: const EdgeInsets.fromLTRB(14, 0, 14, 14),
            child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
              if (run.status == StepStatus.error) ...[
                TonePanel(
                  tone: PanelTone.danger,
                  icon: FA.triangleExclamation,
                  child: Text(run.error ?? l10n.runFailed, style: const TextStyle(fontSize: 13, fontWeight: FontWeight.w600)),
                ),
                const SizedBox(height: 8),
                if (run.tail.isNotEmpty) ...[
                  PreviewBox(label: l10n.runLastLines, text: run.tail.join('\n'), maxHeight: 220),
                  const SizedBox(height: 8),
                ],
                if (run.outputs['log'] != null || _logOf(run) != null)
                  Align(
                    alignment: Alignment.centerLeft,
                    child: ChipButton(
                      label: l10n.runOpenLog,
                      icon: FA.fileLines,
                      onPressed: () => c.openFolder(run.outputs['log'] ?? _logOf(run)!),
                    ),
                  ),
                const SizedBox(height: 8),
              ],
              Text('${l10n.outputDirectory}: ${run.directory}', style: TextStyle(fontFamily: monoFont, fontSize: 12, color: t.textMuted)),
              for (final n in run.notes.skip(2)) Text('• $n', style: TextStyle(fontSize: 13, color: t.textMuted)),
              if (run.command != null) ...[
                const SizedBox(height: 8),
                PreviewBox(label: l10n.commandPreview, text: run.commandDisplay ?? run.command!, maxHeight: 200),
              ],
              const SizedBox(height: 8),
              PreviewBox(label: l10n.outputFiles, text: files.join('\n'), maxHeight: 160),
              if (run.geneTrees.isNotEmpty) ...[
                const SizedBox(height: 8),
                Kicker(l10n.geneTreesLabel),
                const SizedBox(height: 6),
                Wrap(spacing: 6, runSpacing: 6, children: [
                  for (final e in run.geneTrees.entries)
                    TagChip(e.key, icon: e.value ? FA.check : FA.xmark, color: e.value ? SplaceTokens.success : SplaceTokens.danger),
                ]),
              ],
              if (newick != null) ...[
                const SizedBox(height: 10),
                PreviewBox(label: 'Newick', text: newick, maxHeight: 140),
              ],
            ]),
          ),
      ]),
    );
  }
}

/// The program's own log file in the run folder, if it wrote one.
String? _logOf(TreeRun run) {
  for (final name in ['splace.raxml.log', 'splace.log', 'modelfinder.log']) {
    final f = p.join(run.directory, name);
    if (File(f).existsSync()) return f;
  }
  return null;
}
