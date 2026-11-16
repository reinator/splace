import 'dart:io';

import 'package:file_selector/file_selector.dart';
import 'package:flutter/material.dart';
import 'package:flutter/services.dart';
import 'package:url_launcher/url_launcher.dart';

import '../core/genbank_markers.dart';
import '../core/saturation.dart';
import '../services/resource_monitor.dart';
import '../state/project_controller.dart';
import 'chrome/api_key_dialog.dart';
import 'chrome/progress_modal.dart';
import 'chrome/splash.dart';
import 'chrome/tools_dialog.dart';
import 'common.dart';
import 'fa.dart';
import 'sections/organizer_dialog.dart';
import 'sections/alignment_section.dart';
import 'sections/concat_section.dart';
import 'sections/fields_section.dart';
import 'sections/inference_section.dart';
import 'sections/tree_section.dart';
import 'sections/load_section.dart';
import 'sections/markers_section.dart';
import 'sections/records_section.dart';
import 'sections/saturation_section.dart';
import 'theme.dart';
import 'widgets/citations.dart';
import 'widgets/kit.dart';
import 'widgets/splace_tooltip.dart';

class _Sub {
  _Sub(this.id, this.title, this.done);
  final String id;
  final String title;
  final bool done;
}

class _Step {
  _Step(this.number, this.title, this.done, this.locked, {this.meta, this.subs = const []});
  final int number;
  final String title;
  final bool done;
  final bool locked;

  /// Short status under the title ("12 records", "8/8 aligned"…).
  final String? meta;
  final List<_Sub> subs;
}

/// Same layout as the HTML app (`.app-shell`): header, step navigator on
/// the left, one scrolling page with every step, activity log on the right
/// and a fixed footer.
class SplaceShell extends StatefulWidget {
  const SplaceShell({super.key});

  @override
  State<SplaceShell> createState() => _SplaceShellState();
}

class _SplaceShellState extends State<SplaceShell> {
  bool _minimapOpen = true;

  /// Decided on the first layout: open when the panel can be docked.
  bool? _logsOpen;
  bool _settingsOpen = false;
  int _active = 1;
  String? _activeSub;
  final _scroll = ScrollController();
  final _keys = {for (var i = 1; i <= 9; i++) i: GlobalKey()};

  /// The lettered parts of the steps ("2a", "7b"…), filled as they build.
  final _subKeys = <String, GlobalKey>{};

  @override
  void initState() {
    super.initState();
    _scroll.addListener(_trackActiveStep);
    WidgetsBinding.instance.addPostFrameCallback((_) {
      final c = ProjectScope.read(context);
      if (c.logs.isEmpty) c.log(context.l10n.logWelcome);
      if (c.autoDetectTools) c.detectTools();
    });
  }

  @override
  void dispose() {
    _scroll.dispose();
    super.dispose();
  }

  double? _top(GlobalKey key) {
    final box = key.currentContext?.findRenderObject() as RenderBox?;
    if (box == null || !box.attached) return null;
    return box.localToGlobal(Offset.zero).dy;
  }

  /// The active step is the last section whose top passed the header.
  void _trackActiveStep() {
    var active = 1;
    for (final e in _keys.entries) {
      final top = _top(e.value);
      if (top != null && top <= 180) active = e.key;
    }
    String? sub;
    var subTop = double.negativeInfinity;
    for (final e in _subKeys.entries) {
      if (!e.key.startsWith('$active')) continue;
      final top = _top(e.value);
      if (top != null && top <= 220 && top > subTop) {
        sub = e.key;
        subTop = top;
      }
    }
    if (active != _active || sub != _activeSub) {
      setState(() {
        _active = active;
        _activeSub = sub;
      });
    }
  }

  static bool _sameSet(Set<int> a, Set<int> b) => a.length == b.length && a.containsAll(b);

  void _goTo(int step, [String? sub]) {
    final ctx = (sub == null ? _keys[step] : _subKeys[sub])?.currentContext;
    if (ctx == null) return;
    Scrollable.ensureVisible(ctx, duration: const Duration(milliseconds: 350), curve: Curves.easeInOut);
    setState(() {
      _active = step;
      _activeSub = sub;
    });
  }

  List<_Step> _steps(ProjectController c) {
    final l10n = context.l10n;
    final gb = c.records.isNotEmpty;
    final fa = c.inputs.isNotEmpty;
    final loaded = gb || fa;
    final recordsDone = (!gb || c.taxonomyResolved) && (!fa || c.organized != null);
    final typesVisible = (gb && c.taxonomyResolved) || fa;
    final markersVisible = c.markers.isNotEmpty;
    final selected = c.activeMarkers.length;
    final included = c.loci.values.where((l) => l.included).toList();
    final aligned = c.loci.values.where((l) => l.aligned != null).length;
    final trimmed = c.loci.values.where((l) => l.trimmed != null).length;
    final allAligned = included.isNotEmpty && included.every((l) => l.aligned != null || l.alignStatus == StepStatus.skipped);
    final anyAligned = aligned > 0;
    final m = c.matrix;
    final sat = c.pooledSaturation;
    // Step 7 opens once the saturation test ran (its suggestion feeds the
    // matrix), or right away when no locus can be tested by codon position.
    final saturationRan = c.loci.values.any((l) => l.saturation != null);
    // Loci from FASTA files and partial markers: the test is optional.
    final saturationApplies = c.saturationRequired;
    final trees = c.treeRuns.where((r) => r.status == StepStatus.done).length;
    final used = c.loadOptionsUsed;
    bool ranWith(String engine) => c.treeRuns.any((r) => r.status == StepStatus.done && r.engine == engine);
    return [
      _Step(1, l10n.step1, loaded, false,
          meta: loaded ? [if (gb) l10n.metaRecords('${c.records.length}'), if (fa) l10n.metaFastaFiles('${c.inputs.length}')].join(' · ') : null,
          subs: [
            _Sub('1a', l10n.taxTitle, used.contains(1)),
            _Sub('1b', l10n.accTitle, used.contains(2)),
            _Sub('1c', l10n.gbFilesTitle, used.contains(3)),
            _Sub('1d', l10n.filesTitle, used.contains(4)),
          ]),
      _Step(2, l10n.step2, loaded && recordsDone, !loaded, meta: loaded ? l10n.metaSamples('${c.totalSamples}') : null, subs: [
        if (gb) ...[
          _Sub('2a', l10n.siCardTitle, c.speciesIds.isNotEmpty),
          _Sub('2b', l10n.taxonomyCardTitle, c.taxonomyResolved),
          _Sub('2c', l10n.subRecordsTable, c.taxonomyResolved),
        ],
        if (fa) _Sub(fastaHeadersSub(c), l10n.fastaHeadersHeading, c.organized != null),
      ]),
      _Step(3, gb ? l10n.step3 : l10n.step3Fasta, c.featureTypes.isNotEmpty || (!gb && c.organized != null), !typesVisible,
          meta: gb && c.featureTypes.isNotEmpty ? c.featureTypes.map((t) => featureTypeLabels[t] ?? t).join(' · ') : null,
          subs: [
            if (gb && c.taxonomyResolved) _Sub('3a', l10n.subFeatureTypes, c.featureTypes.isNotEmpty),
            if (fa) _Sub(fastaFieldsSub(c), l10n.step3Fasta, c.organized != null),
          ]),
      _Step(4, l10n.step4, selected > 0, !markersVisible,
          meta: markersVisible ? l10n.metaSelected('$selected', '${c.markers.length}') : null,
          subs: [
            _Sub('4a', l10n.subChooseMarkers, selected > 0),
            if (c.markers.isNotEmpty) _Sub('4b', l10n.heatmapTitle, selected > 0),
            if (selected > 0) _Sub('4c', l10n.lengthHeatmapTitle, selected > 0),
          ]),
      _Step(5, l10n.step5, allAligned, selected == 0,
          meta: c.loci.isEmpty ? null : l10n.metaAligned('$aligned', '${c.loci.length}') + (trimmed > 0 ? ' · ${l10n.metaTrimmed('$trimmed')}' : ''),
          subs: [
            if (gb) _Sub('5a', l10n.panelHeader, (c.headerConfirmed || anyAligned) && !c.headerBlocked),
            _Sub('5b', l10n.panelSeqType, anyAligned),
            _Sub('5c', l10n.panelMafft, anyAligned),
            _Sub('5d', l10n.panelRun, anyAligned),
            if (c.alignPanelVisible('5e')) _Sub('5e', l10n.panelTrimal, trimmed > 0 || c.trimSkipped),
          ]),
      _Step(6, l10n.step6, saturationRan, !anyAligned,
          subs: [
            _Sub('6a', l10n.subRunSaturation, saturationRan),
            if (saturationRan) _Sub('6b', l10n.saturationHeatmapTitle, true),
            if (sat != null) _Sub('6c', l10n.recTitle, c.concat.perLocusExclude.isNotEmpty || _sameSet(c.concat.excludePositions, sat.recommendation.exclude)),
          ],
          meta: sat == null
              ? null
              : switch (sat.recommendation.kind) {
                  RecommendationKind.keepAll => l10n.metaSatKeep,
                  RecommendationKind.excludePositions => l10n.metaSatExclude(sat.recommendation.exclude.join('+')),
                  RecommendationKind.useAminoAcids => l10n.metaSatAmino,
                  RecommendationKind.undetermined => l10n.verdictInsufficient,
                }),
      _Step(7, l10n.step7, m != null, m == null && (!anyAligned || (saturationApplies && !saturationRan)),
          meta: m == null ? null : l10n.metaMatrix('${m.taxa.length}', '${m.length}'),
          subs: [
            _Sub('7a', l10n.concatBlockSource, m != null),
            _Sub('7b', l10n.concatBlockPartitions, m != null),
            _Sub('7c', l10n.concatBlockOutgroup, c.concat.outgroups.isNotEmpty),
            _Sub('7d', l10n.subBuildMatrix, m != null),
          ]),
      _Step(8, l10n.step8, trees > 0, m == null, meta: c.treeRuns.isEmpty ? null : l10n.metaTrees('$trees'), subs: [
        _Sub('8a', 'IQ-TREE', ranWith('IQ-TREE')),
        _Sub('8b', 'RAxML-NG', ranWith('RAxML-NG')),
        _Sub('8c', l10n.bothTitle, ranWith('IQ-TREE') && ranWith('RAxML-NG')),
        _Sub('8d', l10n.runs, trees > 0),
      ]),
      _Step(9, l10n.step9, trees > 0, trees == 0, meta: trees == 0 ? null : l10n.metaViewable('${c.viewableTrees.length}'), subs: [
        _Sub('9a', l10n.subTree, trees > 0),
        if (c.treeRuns.isNotEmpty) _Sub('9b', l10n.zipTitle, c.zipSaved),
        _Sub('9c', l10n.howToCite, trees > 0),
      ]),
    ];
  }

  /// The "Sample and marker" dialog is open.
  bool _organizerOpen = false;

  /// New FASTA files (from any button or a drop): ask which part of the
  /// headers is the sample and which the marker.
  void _maybeOpenOrganizer(ProjectController c) {
    if (_organizerOpen || !c.organizerReview || c.busy || c.inputs.isEmpty) return;
    _organizerOpen = true;
    WidgetsBinding.instance.addPostFrameCallback((_) async {
      if (!mounted) return;
      await showOrganizerDialog(context);
      _organizerOpen = false;
    });
  }

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final t = context.tokens;
    _maybeOpenOrganizer(c);
    final steps = _steps(c);
    final progress = steps.where((s) => s.done).length / steps.length;

    return Scaffold(
      body: LayoutBuilder(builder: (context, constraints) {
        // The window is at least 1280×720. The log docks beside the page on
        // wide windows and floats over it otherwise, so the steps keep
        // enough width at the minimum size.
        final dockLogs = constraints.maxWidth >= 1500;
        // The log stays closed until the user opens it.
        final logsOpen = _logsOpen ??= false;
        final minimap = _Minimap(steps: steps, active: _active, activeSub: _activeSub, progress: progress, onTap: _goTo);
        final logs = LogsPanel(onClose: () => setState(() => _logsOpen = false));
        return Stack(children: [
          Column(children: [
            _Header(
              minimapOpen: _minimapOpen,
              logsOpen: logsOpen,
              settingsOpen: _settingsOpen,
              onMenu: () => setState(() => _minimapOpen = !_minimapOpen),
              onLogs: () => setState(() => _logsOpen = !logsOpen),
              onSettings: () => setState(() => _settingsOpen = !_settingsOpen),
            ),
            if (c.busy && c.modal == null) _TaskBar(controller: c),
            Expanded(
              child: Row(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
                if (_minimapOpen)
                  Container(
                    width: 264,
                    decoration: BoxDecoration(color: t.surface, border: Border(right: BorderSide(color: t.border))),
                    child: minimap,
                  ),
                Expanded(child: _main(c, steps)),
                if (logsOpen && dockLogs)
                  Container(
                    width: 340,
                    decoration: BoxDecoration(color: t.surface, border: Border(left: BorderSide(color: t.border))),
                    child: logs,
                  ),
              ]),
            ),
            const _Footer(),
          ]),
          if (logsOpen && !dockLogs)
            Positioned(
              top: 64,
              right: 0,
              bottom: 36,
              width: 400,
              child: Material(elevation: 12, color: t.surface, child: logs),
            ),
          if (_settingsOpen)
            Positioned.fill(
              child: GestureDetector(behavior: HitTestBehavior.opaque, onTap: () => setState(() => _settingsOpen = false)),
            ),
          if (_settingsOpen) const Positioned(top: 60, right: 120, child: _SettingsPopover()),
          if (c.modal != null) Positioned.fill(child: ProgressModal(controller: c)),
          if (c.splashVisible) const Positioned.fill(child: SplashOverlay()),
        ]);
      }),
    );
  }

  Widget _main(ProjectController c, List<_Step> steps) {
    final visible = {for (final s in steps) s.number: !s.locked};
    Widget section(int n, Widget child) => visible[n]! ? KeyedSubtree(key: _keys[n], child: child) : const SizedBox.shrink();
    return SubAnchors(
      keys: _subKeys,
      child: Scrollbar(
      controller: _scroll,
      child: SingleChildScrollView(
        controller: _scroll,
        // The steps use the whole width between the navigator and the log.
        child: Padding(
          padding: const EdgeInsets.fromLTRB(20, 20, 20, 32),
          child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
            section(1, LoadSection(active: _active == 1)),
            section(2, RecordsSection(active: _active == 2)),
            section(3, FieldsSection(active: _active == 3)),
            section(4, MarkersSection(active: _active == 4)),
            section(5, AlignmentSection(active: _active == 5)),
            section(6, SaturationSection(active: _active == 6)),
            section(7, ConcatSection(active: _active == 7)),
            section(8, InferenceSection(active: _active == 8, onViewTree: () => _goTo(9))),
            section(9, TreeSection(active: _active == 9)),
            // Saving everything and how to cite come after the trees (step 8
            // runs, step 9 views).
            if (visible[8]! && c.treeRuns.isNotEmpty) const SubAnchor(id: '9b', child: ResultsZipCard()),
            if (visible[8]!) const SubAnchor(id: '9c', child: CitationsCard()),
          ]),
        ),
      ),
    ),
    );
  }
}

/// Footer: credits on the left; on the right, what the app is using
/// (machine CPU, memory, free disk on the project drive, network).
class _Footer extends StatefulWidget {
  const _Footer();

  @override
  State<_Footer> createState() => _FooterState();
}

class _FooterState extends State<_Footer> {
  ResourceMonitor? _monitor;

  @override
  void initState() {
    super.initState();
    // Not under `flutter test` (no PowerShell loop in tests) unless asked.
    final env = Platform.environment;
    if (!env.containsKey('FLUTTER_TEST') || env['SPLACE_MONITOR'] == '1') {
      final c = ProjectScope.read(context);
      _monitor = ResourceMonitor(folder: () => c.workDir)..start();
    }
  }

  @override
  void dispose() {
    _monitor?.dispose();
    super.dispose();
  }

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    return Container(
      height: 36,
      padding: const EdgeInsets.symmetric(horizontal: 16),
      decoration: BoxDecoration(color: t.surface, border: Border(top: BorderSide(color: t.border))),
      child: Row(children: [
        Expanded(
          child: Text(
            'SPLACE v$appVersion · GPL-3.0 · Luan Rabelo & Renato Moreira-Oliveira',
            maxLines: 1,
            overflow: TextOverflow.ellipsis,
            style: TextStyle(fontSize: 12, color: t.textFaint),
          ),
        ),
        const SizedBox(width: 16),
        if (_monitor != null)
          ListenableBuilder(listenable: _monitor!, builder: (context, _) => _Meters(_monitor!.value))
        else
          const _Meters(ResourceSnapshot()),
      ]),
    );
  }
}

class _Meters extends StatelessWidget {
  const _Meters(this.v);
  final ResourceSnapshot v;

  @override
  Widget build(BuildContext context) {
    final l10n = context.l10n;
    final memUsed = v.memTotal == null || v.memFree == null ? null : v.memTotal! - v.memFree!;
    final diskUsed = v.diskTotal == null || v.diskFree == null ? null : v.diskTotal! - v.diskFree!;
    return Row(mainAxisSize: MainAxisSize.min, children: [
      _Meter(
        icon: FA.gaugeHigh,
        label: l10n.resCpu,
        value: v.cpu == null ? '—' : '${v.cpu!.round()}%',
        fraction: v.cpu == null ? null : v.cpu! / 100,
        tooltip: l10n.resCpuInfo,
      ),
      _Meter(
        icon: FA.microchip,
        label: l10n.resMemory,
        value: v.appMemory == 0 ? '—' : formatBytes(v.appMemory),
        fraction: memUsed == null ? null : memUsed / v.memTotal!,
        tooltip: memUsed == null
            ? l10n.resMemoryInfoApp(formatBytes(v.appMemory))
            : l10n.resMemoryInfo(formatBytes(v.appMemory), formatBytes(memUsed), formatBytes(v.memTotal!)),
      ),
      _Meter(
        icon: FA.hardDrive,
        label: l10n.resDisk,
        value: v.diskFree == null ? '—' : l10n.resFree(formatBytes(v.diskFree!)),
        fraction: diskUsed == null ? null : diskUsed / v.diskTotal!,
        tooltip: v.diskFree == null ? l10n.resDiskInfoNone : l10n.resDiskInfo(formatBytes(v.diskFree!), formatBytes(v.diskTotal!)),
      ),
      _Meter(
        icon: FA.wifi,
        label: l10n.resNetwork,
        value: '${formatBytes(v.netRate)}/s',
        tooltip: l10n.resNetworkInfo(formatBytes(v.netTotal)),
      ),
    ]);
  }
}

/// One resource: icon, name, value and a small bar coloured by load
/// (blue, orange from 70%, red from 90%).
class _Meter extends StatelessWidget {
  const _Meter({required this.icon, required this.label, required this.value, required this.tooltip, this.fraction});
  final IconData icon;
  final String label, value, tooltip;
  final double? fraction;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final f = fraction?.clamp(0.0, 1.0);
    final color = f == null ? t.accent : (f >= 0.9 ? SplaceRed.c600 : (f >= 0.7 ? const Color(0xFFF46D43) : t.accent));
    return SplaceTooltip(
      message: tooltip,
      child: Padding(
        padding: const EdgeInsets.only(left: 16),
        child: Row(mainAxisSize: MainAxisSize.min, children: [
          Icon(icon, size: 11, color: t.textFaint),
          const SizedBox(width: 5),
          Text(label, style: TextStyle(fontSize: 11.5, fontWeight: FontWeight.w600, color: t.textMuted)),
          const SizedBox(width: 5),
          if (f != null) ...[
            Container(
              width: 38,
              height: 5,
              decoration: BoxDecoration(color: t.surface3, borderRadius: BorderRadius.circular(3)),
              alignment: Alignment.centerLeft,
              child: FractionallySizedBox(
                widthFactor: f,
                child: Container(decoration: BoxDecoration(color: color, borderRadius: BorderRadius.circular(3))),
              ),
            ),
            const SizedBox(width: 5),
          ],
          Text(value, style: TextStyle(fontSize: 11.5, color: t.textMuted, fontFeatures: const [FontFeature.tabularFigures()])),
        ]),
      ),
    );
  }
}

// ---------------------------------------------------------------------------
// Header
// ---------------------------------------------------------------------------

class _Header extends StatelessWidget {
  const _Header({
    required this.minimapOpen,
    required this.logsOpen,
    required this.settingsOpen,
    required this.onMenu,
    required this.onLogs,
    required this.onSettings,
  });

  final bool minimapOpen, logsOpen, settingsOpen;
  final VoidCallback onMenu, onLogs, onSettings;

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final t = context.tokens;
    final l10n = context.l10n;
    final lang = c.locale?.languageCode ?? Localizations.localeOf(context).languageCode;
    final bold = TextStyle(fontWeight: FontWeight.w700, color: t.accentText);
    return Container(
      height: 64,
      padding: const EdgeInsets.symmetric(horizontal: 20),
      decoration: BoxDecoration(color: t.surface, border: Border(bottom: BorderSide(color: t.border))),
      child: Row(children: [
        IconBtn(icon: FA.bars, tooltip: l10n.tooltipMenu, active: minimapOpen, onPressed: onMenu),
        const SizedBox(width: 16),
        const SplaceLogo(height: 40, radius: 8, shadow: false),
        const SizedBox(width: 12),
        Column(mainAxisAlignment: MainAxisAlignment.center, crossAxisAlignment: CrossAxisAlignment.start, children: [
          Text('SPLACE', style: TextStyle(fontSize: 20, fontWeight: FontWeight.w700, color: t.text, height: 1.1, letterSpacing: -0.2)),
          const SizedBox(height: 2),
          Text.rich(
            TextSpan(style: TextStyle(fontSize: 12, color: t.textMuted, height: 1.1), children: [
              TextSpan(text: 'SPL', style: bold),
              const TextSpan(text: 'it, '),
              TextSpan(text: 'A', style: bold),
              const TextSpan(text: 'lign and '),
              TextSpan(text: 'C', style: bold),
              const TextSpan(text: 'oncatenat'),
              TextSpan(text: 'E', style: bold),
            ]),
          ),
        ]),
        const Spacer(),
        LangToggle(current: lang, onChanged: (code) => c.setLocale(Locale(code))),
        const SizedBox(width: 8),
        IconBtn(
          icon: FA.key,
          tooltip: l10n.apiKeyTitle,
          active: c.ncbiApiKey != null,
          onPressed: () => showApiKeyDialog(context),
        ),
        const SizedBox(width: 8),
        IconBtn(icon: FA.screwdriverWrench, tooltip: l10n.tooltipTools, onPressed: () => showToolsDialog(context)),
        const SizedBox(width: 8),
        IconBtn(icon: FA.sliders, tooltip: l10n.tooltipSettings, active: settingsOpen, onPressed: onSettings),
        const SizedBox(width: 8),
        IconBtn(icon: FA.rectangleList, tooltip: l10n.tooltipLogs, active: logsOpen, onPressed: onLogs),
        const SizedBox(width: 8),
        IconBtn(icon: FA.circleQuestion, tooltip: l10n.tooltipHelp, onPressed: () => c.setSplashVisible(true)),
        const SizedBox(width: 8),
        IconBtn(icon: FA.github, tooltip: 'GitHub', onPressed: () => launchUrl(Uri.parse(githubUrl))),
      ]),
    );
  }
}

class _TaskBar extends StatelessWidget {
  const _TaskBar({required this.controller});
  final ProjectController controller;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final p = controller.progress;
    return Container(
      decoration: BoxDecoration(color: t.accentSoft, border: Border(bottom: BorderSide(color: t.accentSoftBorder))),
      child: Column(children: [
        // Real progress only; tasks without a count show just the spinner.
        if (p != null) RainbowProgressBar(value: p, height: 3),
        Padding(
          padding: const EdgeInsets.symmetric(horizontal: 20, vertical: 6),
          child: Row(children: [
            SizedBox(width: 14, height: 14, child: CircularProgressIndicator(strokeWidth: 2, color: t.accent)),
            const SizedBox(width: 10),
            Expanded(
              child: Text(
                [
                  p == null ? controller.busyTask! : '${controller.busyTask!} — ${(p * 100).round()}%',
                  ?controller.progressDetail,
                ].join(' · '),
                maxLines: 1,
                overflow: TextOverflow.ellipsis,
                style: TextStyle(fontSize: 13, fontWeight: FontWeight.w600, color: t.accentText),
              ),
            ),
            ChipButton(label: context.l10n.cancel, icon: FA.stop, style: ChipButtonStyle.danger, onPressed: controller.cancelTask),
          ]),
        ),
      ]),
    );
  }
}

// ---------------------------------------------------------------------------
// Step navigator (`.app-minimap`)
// ---------------------------------------------------------------------------

/// The steps with their lettered parts: the active step lists them, the
/// others show them as small letters (green when done), so the progress of
/// every step is in view without a long list. Follows the active step.
class _Minimap extends StatefulWidget {
  const _Minimap({required this.steps, required this.active, required this.activeSub, required this.progress, required this.onTap});
  final List<_Step> steps;
  final int active;
  final String? activeSub;
  final double progress;
  final void Function(int step, [String? sub]) onTap;

  @override
  State<_Minimap> createState() => _MinimapState();
}

class _MinimapState extends State<_Minimap> {
  final _rows = {for (var i = 1; i <= 9; i++) i: GlobalKey()};

  @override
  void didUpdateWidget(_Minimap old) {
    super.didUpdateWidget(old);
    if (old.active == widget.active) return;
    WidgetsBinding.instance.addPostFrameCallback((_) {
      final ctx = _rows[widget.active]?.currentContext;
      if (ctx == null || !mounted) return;
      // Only when it is out of view (both calls do nothing otherwise).
      Scrollable.ensureVisible(ctx, duration: const Duration(milliseconds: 200), alignmentPolicy: ScrollPositionAlignmentPolicy.keepVisibleAtEnd);
      Scrollable.ensureVisible(ctx, duration: const Duration(milliseconds: 200), alignmentPolicy: ScrollPositionAlignmentPolicy.keepVisibleAtStart);
    });
  }

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final l10n = context.l10n;
    final steps = widget.steps, active = widget.active, activeSub = widget.activeSub, progress = widget.progress;
    final onTap = widget.onTap;
    return ListView(padding: EdgeInsets.zero, children: [
      Container(
        padding: const EdgeInsets.fromLTRB(20, 14, 20, 14),
        decoration: BoxDecoration(border: Border(bottom: BorderSide(color: t.borderSoft))),
        child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
          Row(children: [
            Expanded(child: Kicker(l10n.navProgress, color: t.textFaint, size: 11)),
            Text('${(progress * 100).round()}%', style: TextStyle(fontSize: 12, fontWeight: FontWeight.w700, color: t.accentText)),
          ]),
          const SizedBox(height: 8),
          RainbowProgressBar(value: progress, height: 4),
        ]),
      ),
      Container(
        padding: const EdgeInsets.fromLTRB(12, 12, 12, 12),
        child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
          Padding(
            padding: const EdgeInsets.only(left: 8, bottom: 6),
            child: Kicker(l10n.navSteps, color: t.textFaint, size: 11),
          ),
          for (final s in steps) ...[
            _MinimapStep(
              key: _rows[s.number],
              step: s,
              active: s.number == active && !s.locked,
              onTap: () => onTap(s.number),
              // Other steps: their parts as small letters.
              subs: s.locked || s.number == active || s.subs.isEmpty ? null : _SubDots(subs: s.subs, onTap: (id) => onTap(s.number, id)),
            ),
            if (!s.locked && s.number == active)
              for (final sub in s.subs)
                _MinimapSub(sub: sub, active: activeSub == sub.id, onTap: () => onTap(s.number, sub.id)),
          ],
        ]),
      ),
    ]);
  }
}

class _MinimapStep extends StatefulWidget {
  const _MinimapStep({super.key, required this.step, required this.active, required this.onTap, this.subs});
  final _Step step;
  final bool active;
  final VoidCallback onTap;

  /// The parts of a step that is not the active one.
  final Widget? subs;

  @override
  State<_MinimapStep> createState() => _MinimapStepState();
}

class _MinimapStepState extends State<_MinimapStep> {
  bool _hover = false;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final s = widget.step;
    final active = widget.active;
    final hover = _hover && !s.locked;
    final numFilled = active || s.done;
    final row = AnimatedContainer(
      duration: const Duration(milliseconds: 150),
      margin: const EdgeInsets.only(bottom: 2),
      padding: const EdgeInsets.symmetric(horizontal: 12, vertical: 6),
      decoration: BoxDecoration(
        color: active ? t.accentSoft : (hover ? t.surface2 : Colors.transparent),
        borderRadius: BorderRadius.circular(10),
        border: Border(left: BorderSide(color: active ? t.accent : Colors.transparent, width: 3)),
      ),
      child: Row(children: [
        Container(
          width: 26,
          height: 26,
          alignment: Alignment.center,
          decoration: BoxDecoration(
            color: s.done ? SplaceTokens.success : (numFilled ? t.accent : t.surface3),
            shape: BoxShape.circle,
          ),
          child: s.locked
              ? Icon(FA.lock, size: 10, color: t.textFaint)
              : s.done && !active
                  ? const Icon(FA.check, size: 11, color: Colors.white)
                  : Text('${s.number}', style: TextStyle(fontSize: 12, fontWeight: FontWeight.w700, color: numFilled ? Colors.white : t.textMuted)),
        ),
        const SizedBox(width: 12),
        Expanded(
          child: Column(crossAxisAlignment: CrossAxisAlignment.start, mainAxisSize: MainAxisSize.min, children: [
            Text(s.title,
                maxLines: 2,
                overflow: TextOverflow.ellipsis,
                style: TextStyle(
                  fontSize: 14,
                  height: 1.25,
                  fontWeight: active ? FontWeight.w600 : FontWeight.w500,
                  color: active ? t.accentText : (hover ? t.text : t.textMuted),
                )),
            if (s.meta != null && s.meta!.isNotEmpty)
              Padding(
                padding: const EdgeInsets.only(top: 2),
                child: Text(s.meta!, maxLines: 1, overflow: TextOverflow.ellipsis, style: TextStyle(fontSize: 11, color: t.textFaint)),
              ),
            if (widget.subs != null) Padding(padding: const EdgeInsets.only(top: 4), child: widget.subs),
          ]),
        ),
      ]),
    );
    return SplaceTooltip(
      message: s.locked ? context.l10n.minimapLocked : '',
      child: Opacity(
        opacity: s.locked ? 0.5 : 1,
        child: MouseRegion(
          cursor: s.locked ? SystemMouseCursors.forbidden : SystemMouseCursors.click,
          onEnter: (_) => setState(() => _hover = true),
          onExit: (_) => setState(() => _hover = false),
          child: GestureDetector(onTap: s.locked ? null : widget.onTap, child: row),
        ),
      ),
    );
  }
}

/// The parts of a step in one line: a small circle per letter, green with
/// a check when done; hover names it, a click goes there.
class _SubDots extends StatelessWidget {
  const _SubDots({required this.subs, required this.onTap});
  final List<_Sub> subs;
  final ValueChanged<String> onTap;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    return Wrap(spacing: 4, runSpacing: 4, children: [
      for (final sub in subs)
        SplaceTooltip(
          message: '${sub.id} · ${sub.title}',
          child: MouseRegion(
            cursor: SystemMouseCursors.click,
            child: GestureDetector(
              onTap: () => onTap(sub.id),
              child: Container(
                width: 17,
                height: 17,
                alignment: Alignment.center,
                decoration: BoxDecoration(
                  color: sub.done ? SplaceTokens.success : Colors.transparent,
                  shape: BoxShape.circle,
                  border: Border.all(color: sub.done ? SplaceTokens.success : t.border, width: 1.2),
                ),
                child: Text(sub.id.substring(1),
                    style: TextStyle(fontSize: 9.5, height: 1, fontWeight: FontWeight.w700, color: sub.done ? Colors.white : t.textFaint)),
              ),
            ),
          ),
        ),
    ]);
  }
}

/// Lettered parts of a step ("2a", "5c"…), indented under it.
class _MinimapSub extends StatelessWidget {
  const _MinimapSub({required this.sub, required this.active, required this.onTap});
  final _Sub sub;
  final bool active;
  final VoidCallback onTap;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    return MouseRegion(
      cursor: SystemMouseCursors.click,
      child: GestureDetector(
        onTap: onTap,
        child: Container(
          margin: const EdgeInsets.only(left: 26, bottom: 1),
          padding: const EdgeInsets.symmetric(horizontal: 10, vertical: 3),
          decoration: BoxDecoration(
            color: active ? t.accentSoft : Colors.transparent,
            borderRadius: BorderRadius.circular(8),
            border: Border(left: BorderSide(color: active ? t.accent : t.border, width: 2)),
          ),
          child: Row(children: [
            SizedBox(
              width: 22,
              child: Text(sub.id, style: TextStyle(fontSize: 11, fontWeight: FontWeight.w700, color: sub.done ? SplaceTokens.success : t.textFaint)),
            ),
            Expanded(
              child: Text(sub.title,
                  maxLines: 1,
                  overflow: TextOverflow.ellipsis,
                  style: TextStyle(fontSize: 12, color: active ? t.accentText : t.textMuted, fontWeight: active ? FontWeight.w600 : FontWeight.w400)),
            ),
            if (sub.done) const Icon(FA.check, size: 9, color: SplaceTokens.success),
          ]),
        ),
      ),
    );
  }
}

// ---------------------------------------------------------------------------
// Settings popover (`.settings-popover`)
// ---------------------------------------------------------------------------

class _SettingsPopover extends StatelessWidget {
  const _SettingsPopover();

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final t = context.tokens;
    final l10n = context.l10n;
    Widget row(String label, Widget child) => Padding(
          padding: const EdgeInsets.symmetric(vertical: 8),
          child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
            Kicker(label, size: 11),
            const SizedBox(height: 6),
            child,
          ]),
        );
    return Material(
      color: Colors.transparent,
      child: Container(
        width: 300,
        padding: const EdgeInsets.symmetric(horizontal: 16, vertical: 8),
        decoration: BoxDecoration(
          color: t.surface,
          borderRadius: BorderRadius.circular(14),
          border: Border.all(color: t.border),
          boxShadow: shadowLg,
        ),
        child: Column(mainAxisSize: MainAxisSize.min, crossAxisAlignment: CrossAxisAlignment.stretch, children: [
          row(
            l10n.settingsTheme,
            Seg<ThemeMode>(
              value: c.themeMode,
              options: [(ThemeMode.light, l10n.themeLight), (ThemeMode.dark, l10n.themeDark), (ThemeMode.system, l10n.themeAuto)],
              onChanged: c.setThemeMode,
            ),
          ),
          row(
            l10n.settingsDensity,
            Seg<bool>(
              value: c.compact,
              options: [(false, l10n.densityComfy), (true, l10n.densityCompact)],
              onChanged: c.setCompact,
            ),
          ),
          row(
            l10n.settingsFontSize,
            Row(children: [
              Text('A', style: TextStyle(fontSize: 11, color: t.textMuted)),
              Expanded(
                child: Slider(value: c.fontSize, min: 13, max: 20, divisions: 7, onChanged: c.setFontSize),
              ),
              Text('A', style: TextStyle(fontSize: 17, color: t.textMuted)),
              const SizedBox(width: 8),
              Text('${c.fontSize.round()}px', style: TextStyle(fontSize: 12, fontWeight: FontWeight.w700, color: t.text)),
              const SizedBox(width: 4),
              IconBtn(icon: FA.rotateLeft, tooltip: l10n.fontSizeReset, size: ChipButton.height, bordered: true, onPressed: () => c.setFontSize(16)),
            ]),
          ),
        ]),
      ),
    );
  }
}

// ---------------------------------------------------------------------------
// Activity log (`.app-logs`)
// ---------------------------------------------------------------------------

class LogsPanel extends StatelessWidget {
  const LogsPanel({super.key, required this.onClose});
  final VoidCallback onClose;

  String _text(ProjectController c) =>
      c.logs.map((e) => '${e.time.toIso8601String().substring(11, 19)}\t${e.level.name}\t[${e.source}] ${e.message}').join('\n');

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.read(context);
    final t = context.tokens;
    final l10n = context.l10n;
    return Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
      Container(
        padding: const EdgeInsets.symmetric(horizontal: 16, vertical: 12),
        decoration: BoxDecoration(border: Border(bottom: BorderSide(color: t.border))),
        child: Row(children: [
          Icon(FA.rectangleList, size: 16, color: t.accent),
          const SizedBox(width: 8),
          Expanded(child: Text(l10n.logTitle, style: TextStyle(fontSize: 14, fontWeight: FontWeight.w600, color: t.text))),
          IconBtn(icon: FA.copy, tooltip: l10n.logCopy, size: ChipButton.height, bordered: true, onPressed: () => Clipboard.setData(ClipboardData(text: _text(c)))),
          const SizedBox(width: 6),
          IconBtn(
            icon: FA.download,
            tooltip: l10n.logDownload,
            size: 28,
            bordered: true,
            onPressed: () async {
              final loc = await getSaveLocation(suggestedName: 'splace_log.txt');
              if (loc != null) await File(loc.path).writeAsString(_text(c));
            },
          ),
          const SizedBox(width: 6),
          IconBtn(icon: FA.broom, tooltip: l10n.logClear, size: ChipButton.height, bordered: true, onPressed: c.clearLogs),
          const SizedBox(width: 6),
          IconBtn(icon: FA.xmark, tooltip: l10n.logCollapse, size: ChipButton.height, bordered: true, onPressed: onClose),
        ]),
      ),
      Expanded(child: _LogBody(controller: c)),
    ]);
  }
}

/// Oldest line first, like the HTML log; follows new lines while the user
/// is at the bottom.
class _LogBody extends StatefulWidget {
  const _LogBody({required this.controller});
  final ProjectController controller;

  @override
  State<_LogBody> createState() => _LogBodyState();
}

class _LogBodyState extends State<_LogBody> {
  final _scroll = ScrollController();
  bool _follow = true;

  @override
  void initState() {
    super.initState();
    widget.controller.logTick.addListener(_onLog);
    _scroll.addListener(() {
      if (_scroll.hasClients) _follow = _scroll.position.extentAfter < 40;
    });
  }

  void _onLog() {
    setState(() {});
    if (_follow) {
      WidgetsBinding.instance.addPostFrameCallback((_) {
        if (_scroll.hasClients) _scroll.jumpTo(_scroll.position.maxScrollExtent);
      });
    }
  }

  @override
  void dispose() {
    widget.controller.logTick.removeListener(_onLog);
    _scroll.dispose();
    super.dispose();
  }

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final logs = widget.controller.logs;
    if (logs.isEmpty) {
      return Padding(
        padding: const EdgeInsets.all(24),
        child: Column(children: [
          Icon(FA.rectangleList, size: 24, color: t.textFaint.withValues(alpha: 0.4)),
          const SizedBox(height: 8),
          Text(context.l10n.logEmpty, textAlign: TextAlign.center, style: TextStyle(fontSize: 13, color: t.textFaint)),
        ]),
      );
    }
    return Scrollbar(
      controller: _scroll,
      child: ListView.builder(
        controller: _scroll,
        padding: const EdgeInsets.all(12),
        itemCount: logs.length,
        itemBuilder: (context, i) => _LogLine(entry: logs[i]),
      ),
    );
  }
}

class _LogLine extends StatelessWidget {
  const _LogLine({required this.entry});
  final LogEntry entry;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final (icon, color, bg) = switch (entry.level) {
      LogLevel.info => (FA.circleInfo, t.accent, Colors.transparent),
      LogLevel.success => (FA.circleCheck, SplaceTokens.success, Colors.transparent),
      LogLevel.warn => (FA.triangleExclamation, SplaceTokens.warn, const Color(0x0FF46D43)),
      LogLevel.error => (FA.circleXmark, SplaceTokens.danger, const Color(0x0FD73027)),
    };
    final mono = TextStyle(fontFamily: monoFont, fontSize: 12, color: t.text, height: 1.55);
    return Container(
      padding: const EdgeInsets.symmetric(horizontal: 4, vertical: 2),
      decoration: BoxDecoration(color: bg, borderRadius: BorderRadius.circular(4)),
      child: Row(crossAxisAlignment: CrossAxisAlignment.start, children: [
        SizedBox(
          width: 56,
          child: Padding(
            padding: const EdgeInsets.only(top: 2),
            child: Text(entry.time.toIso8601String().substring(11, 19), style: mono.copyWith(fontSize: 11, color: t.textFaint)),
          ),
        ),
        SizedBox(width: 18, child: Padding(padding: const EdgeInsets.only(top: 4), child: Icon(icon, size: 11, color: color))),
        const SizedBox(width: 8),
        Expanded(
          child: SelectableText.rich(TextSpan(style: mono, children: [
            if (entry.source != 'SPLACE') TextSpan(text: '[${entry.source}] ', style: TextStyle(color: t.textMuted)),
            TextSpan(text: entry.message),
          ])),
        ),
      ]),
    );
  }
}
