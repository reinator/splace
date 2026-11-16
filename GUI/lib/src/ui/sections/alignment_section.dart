import 'dart:io';
import 'dart:math' as math;

import 'package:flutter/material.dart';

import '../../core/fasta.dart';
import '../../core/genetic_code.dart';
import '../../core/tool_commands.dart';
import '../../services/tool_locator.dart';
import '../../state/project_controller.dart';
import '../common.dart';
import '../fa.dart';
import '../theme.dart';
import '../widgets/alignment_viewer.dart';
import '../widgets/examples.dart';
import '../widgets/kit.dart';
import 'header_builder.dart';
import '../widgets/splace_tooltip.dart';

const _mafftMethods = ['--auto', '--localpair', '--globalpair', '--genafpair', '--retree 2 --reorder', '--retree 1'];
const _trimalModes = ['-automated1', '-strict', '-strictplus', '-gappyout', '-nogaps', '-noallgaps', 'manual'];
const _clipkitModes = ['gappy', 'kpic', 'kpic-gappy', 'smart-gap'];

/// Step 5 — Sequence Alignment (`#nextStepSection`, panels 5a–5e and the
/// alignment results).
class AlignmentSection extends StatelessWidget {
  const AlignmentSection({super.key, required this.active});
  final bool active;

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    final gb = c.records.isNotEmpty;
    final results = c.loci.values.where((l) => l.aligned != null || l.alignStatus != StepStatus.idle).toList();

    Widget keyed(String id, Widget child) => SubAnchor(id: id, child: child);

    return Padding(
      padding: EdgeInsets.only(bottom: t.gapSection),
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        Padding(
          padding: const EdgeInsets.only(bottom: 16),
          child: Row(children: [
            StepBadge('5', active: active),
            const SizedBox(width: 8),
            Icon(FA.alignLeft, size: 16, color: t.accent),
            const SizedBox(width: 8),
            Text(l10n.step5, style: TextStyle(fontSize: 16, fontWeight: FontWeight.w600, color: t.text)),
            const SizedBox(width: 8),
            CountPill('${c.activeMarkers.length}'),
          ]),
        ),
        // 5a–5d show together with their defaults: one click aligns.
        if (gb) keyed('5a', const _HeaderPanel()),
        keyed('5b', const _SeqTypePanel()),
        keyed('5c', const _MafftPanel()),
        keyed('5d', const _RunPanel()),
        if (c.alignPanelVisible('5e')) keyed('5e', const _TrimPanel()),
        if (results.isNotEmpty)
          StepSection(
            icon: FA.circleCheck,
            title: l10n.resultsTitle,
            count: results.length,
            actions: [
              if (c.loci.values.any((l) => l.aligned != null))
                ChipButton(label: l10n.openFolder, icon: FA.folderOpen, onPressed: () => c.openFolder(c.dir('02_aligned'))),
            ],
            child: _ResultsList(results: results),
          ),
      ]),
    );
  }
}

// ---------------------------------------------------------------------------
// 5a Configure FASTA Header
// ---------------------------------------------------------------------------

class _HeaderPanel extends StatelessWidget {
  const _HeaderPanel();

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    final examples = c.records.take(3).map(c.header.build).toList();
    final dups = c.headerConflicts();
    return SubPanel(
      badge: '5a',
      title: l10n.panelHeader,
      body: l10n.panelHeaderBody,
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        Row(children: [
          PrimaryButton(label: l10n.openHeaderBuilder, icon: FA.tags, onPressed: () => showHeaderBuilder(context)),
          const SizedBox(width: 12),
          if (c.headerBlocked)
            TagChip(l10n.alignBlockedHeader, icon: FA.lock, color: SplaceTokens.danger)
          else if (c.headerConfirmed)
            TagChip(l10n.headerConfirmed, icon: FA.check, color: SplaceTokens.success)
          else
            TagChip(l10n.headerDefault, icon: FA.circleInfo, color: t.textMuted),
        ]),
        const SizedBox(height: 12),
        PreviewBox(label: l10n.previewLabel, text: examples.toSet().map((e) => '>$e').join('\n'), maxHeight: 90),
        if (dups.isNotEmpty) ...[
          const SizedBox(height: 10),
          TonePanel(
            tone: c.onePerTaxonMode ? PanelTone.info : PanelTone.danger,
            icon: c.onePerTaxonMode ? FA.circleInfo : FA.circleExclamation,
            child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
              Text(
                c.onePerTaxonMode ? l10n.onePerTaxonJoined('${dups.length}') : l10n.headerConflictBlock('${dups.length}'),
                style: TextStyle(fontSize: 13, color: t.text),
              ),
              if (!c.onePerTaxonMode) ...[
                if (!c.header.fields.contains('accession')) ...[
                  const SizedBox(height: 8),
                  ChipButton(
                    label: l10n.headerAddAccession,
                    icon: FA.plus,
                    style: ChipButtonStyle.accent,
                    onPressed: c.busy ? null : () => c.setHeader(c.header.copyWith(fields: [...c.header.fields, 'accession'])),
                  ),
                ],
                const SizedBox(height: 6),
                for (final e in dups.entries.take(4))
                  Text('${e.key}: ${e.value.join(', ')}', maxLines: 1, overflow: TextOverflow.ellipsis, style: TextStyle(fontSize: 12.5, color: t.textMuted)),
                if (dups.length > 4) Text(l10n.andMore('${dups.length - 4}'), style: TextStyle(fontSize: 12.5, color: t.textMuted)),
              ],
            ]),
          ),
        ],
      ]),
    );
  }
}

// ---------------------------------------------------------------------------
// 5b Sequence Type
// ---------------------------------------------------------------------------

class _SeqTypePanel extends StatelessWidget {
  const _SeqTypePanel();

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    final codon = c.seqMode == SeqMode.codon;
    final translating = c.seqMode != SeqMode.nucleotide;

    // One card per mode: choice, what it does and when to use it. The
    // three cards share one height and their "when to use" dividers sit on
    // one line: the backgrounds are drawn behind two rows of content (the
    // descriptions, then the recommendations), each row as tall as its
    // tallest cell.
    final modes = [
      (SeqMode.nucleotide, l10n.seqNtTitle, l10n.seqNtBody, FA.dna, l10n.recNtBody),
      (SeqMode.aminoAcid, l10n.seqAaTitle, l10n.seqAaBody, FA.puzzlePiece, l10n.recAaBody),
      (SeqMode.codon, l10n.seqCodonTitle, l10n.seqCodonBody, FA.shieldHalved, l10n.recCodonBody),
    ];
    Widget cells(Widget Function(SeqMode m, String title, String body, IconData icon, String rec, bool selected) cell) => IntrinsicHeight(
          child: Row(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
            for (final (i, (m, title, body, icon, rec)) in modes.indexed) ...[
              if (i > 0) const SizedBox(width: 12),
              Expanded(
                child: MouseRegion(
                  cursor: SystemMouseCursors.click,
                  child: GestureDetector(
                    behavior: HitTestBehavior.opaque,
                    onTap: () => c.setSeqMode(m),
                    child: cell(m, title, body, icon, rec, c.seqMode == m),
                  ),
                ),
              ),
            ],
          ]),
        );
    final cards = Stack(children: [
      Positioned.fill(
        child: Row(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
          for (final (i, (m, _, _, _, _)) in modes.indexed) ...[
            if (i > 0) const SizedBox(width: 12),
            Expanded(
              child: AnimatedContainer(
                duration: const Duration(milliseconds: 150),
                decoration: BoxDecoration(
                  color: c.seqMode == m ? t.accentSoft : t.surface,
                  borderRadius: BorderRadius.circular(12),
                  border: Border.all(color: c.seqMode == m ? SplaceBlue.c500 : t.border, width: 2),
                ),
              ),
            ),
          ],
        ]),
      ),
      Column(mainAxisSize: MainAxisSize.min, crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        cells((m, title, body, icon, rec, selected) => Padding(
              padding: const EdgeInsets.fromLTRB(16, 16, 16, 12),
              child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
                Row(children: [
                  Container(
                    width: 16,
                    height: 16,
                    decoration: BoxDecoration(shape: BoxShape.circle, border: Border.all(color: selected ? t.accent : t.textFaint, width: selected ? 5 : 1.5)),
                  ),
                  const SizedBox(width: 10),
                  Expanded(child: Text(title, style: TextStyle(fontSize: 15, fontWeight: FontWeight.w700, color: selected ? t.accentText : t.text))),
                ]),
                const SizedBox(height: 6),
                Text(body, style: TextStyle(fontSize: 13, color: t.textMuted, height: 1.4)),
              ]),
            )),
        cells((m, title, body, icon, rec, selected) => Padding(
              padding: const EdgeInsets.fromLTRB(16, 0, 16, 12),
              child: SeqModeExample(mode: m),
            )),
        cells((m, title, body, icon, rec, selected) => Padding(
              padding: const EdgeInsets.fromLTRB(16, 0, 16, 16),
              child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
                Divider(height: 1, color: selected ? t.accentSoftBorder : t.border),
                const SizedBox(height: 10),
                Row(children: [
                  Icon(icon, size: 12, color: t.accent),
                  const SizedBox(width: 6),
                  Expanded(child: Text(l10n.whenToUse, style: TextStyle(fontWeight: FontWeight.w700, fontSize: 12.5, color: t.text))),
                ]),
                const SizedBox(height: 4),
                Text(rec, style: TextStyle(fontSize: 12.5, color: t.textMuted, height: 1.4)),
              ]),
            )),
      ]),
    ]);

    return SubPanel(
      badge: '5b',
      title: l10n.panelSeqType,
      body: l10n.seqTypeBody,
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        cards,
        // The genetic code translates the CDS in the amino-acid and codon
        // modes. In the nucleotide mode the viewer picks its own code for
        // the amino-acid view.
        if (translating) ...[
          const SizedBox(height: 16),
          AutoGrid(minItemWidth: 340, children: [
            Field(
              label: l10n.geneticCode,
              info: l10n.geneticCodeInfo,
              child: DropdownButtonFormField<int>(
                icon: dropdownChevron,
                initialValue: c.geneticCode,
                isExpanded: true,
                items: [for (final e in geneticCodeNames.entries) DropdownMenuItem(value: e.key, child: Text(e.value))],
                onChanged: (v) => c.setGeneticCode(v!),
              ),
            ),
            Padding(
              padding: const EdgeInsets.only(top: 22),
              child: TonePanel(tone: PanelTone.info, child: Text(l10n.translationInfo, style: const TextStyle(fontSize: 13))),
            ),
          ]),
        ],
        if (codon) ...[
          const SizedBox(height: 12),
          TonePanel(tone: PanelTone.warn, child: Text(l10n.codonWarning, style: const TextStyle(fontSize: 13))),
          const SizedBox(height: 12),
          Kicker(l10n.codonWorkflow),
          const SizedBox(height: 8),
          // Five steps in one row of equal-height cards.
          IntrinsicHeight(
            child: Row(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
              for (final (i, (title, body)) in [
                (l10n.codonStep1, l10n.codonStep1Body),
                (l10n.codonStep2, l10n.codonStep2Body),
                (l10n.codonStep3, l10n.codonStep3Body),
                (l10n.codonStep4, l10n.codonStep4Body),
                (l10n.codonStep5, l10n.codonStep5Body),
              ].indexed) ...[
                if (i > 0) const SizedBox(width: 8),
                Expanded(
                  child: Container(
                    padding: const EdgeInsets.all(10),
                    decoration: BoxDecoration(color: t.accentSoft, borderRadius: BorderRadius.circular(10), border: Border.all(color: t.accentSoftBorder)),
                    child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
                      Text('${i + 1}. $title', style: TextStyle(fontWeight: FontWeight.w700, fontSize: 13, color: t.accentText)),
                      const SizedBox(height: 2),
                      Text(body, style: TextStyle(fontSize: 12, color: t.textMuted)),
                    ]),
                  ),
                ),
              ],
            ]),
          ),
          const SizedBox(height: 12),
          CheckCard(value: c.allowPseudo, onChanged: (v) => c.setCodonPolicy(pseudo: v), title: l10n.allowPseudo, subtitle: l10n.allowPseudoBody),
          const SizedBox(height: 8),
          CheckCard(
              value: c.allowInternalStops, onChanged: (v) => c.setCodonPolicy(internalStops: v), title: l10n.allowStops, subtitle: l10n.allowStopsBody),
          const SizedBox(height: 8),
          CheckCard(
            value: c.allowBacktranslationWarnings,
            onChanged: (v) => c.setCodonPolicy(backtranslationWarnings: v),
            title: l10n.allowBtWarnings,
            subtitle: l10n.allowBtWarningsBody,
          ),
        ],
      ]),
    );
  }
}

/// "Python tools for codon trimming" (ClipKIT status and installer).
class _PythonCard extends StatelessWidget {
  const _PythonCard();

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    final clip = c.tools[ExternalTool.clipkit];
    final status = !c.toolsChecked
        ? l10n.pythonNotChecked
        : c.python == null
            ? l10n.pythonMissing
            : '${l10n.pythonFound(c.python!.display, c.python!.version)} · ${clip == null ? l10n.clipkitMissing : l10n.clipkitFound(clip.version ?? '')}';
    return Container(
      padding: const EdgeInsets.all(14),
      decoration: BoxDecoration(color: t.surface2, borderRadius: BorderRadius.circular(12), border: Border.all(color: t.border)),
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        Row(children: [
          Icon(FA.python, size: 15, color: t.accent),
          const SizedBox(width: 8),
          Text(l10n.pythonTitle, style: TextStyle(fontWeight: FontWeight.w700, color: t.text)),
        ]),
        const SizedBox(height: 4),
        Text(l10n.pythonBody, style: TextStyle(fontSize: 13, color: t.textMuted)),
        const SizedBox(height: 8),
        Text(status, style: TextStyle(fontFamily: monoFont, fontSize: 12, color: clip != null ? SplaceTokens.success : t.textMuted)),
        const SizedBox(height: 10),
        Wrap(spacing: 8, children: [
          ChipButton(label: l10n.checkPython, icon: FA.arrowsRotate, onPressed: c.busy ? null : c.detectTools),
          ChipButton(
            label: l10n.installClipkit,
            icon: FA.download,
            style: ChipButtonStyle.soft,
            onPressed: c.busy || c.python == null ? null : c.installClipkit,
          ),
        ]),
      ]),
    );
  }
}

// ---------------------------------------------------------------------------
// 5c MAFFT Settings
// ---------------------------------------------------------------------------

class _MafftPanel extends StatefulWidget {
  const _MafftPanel();

  @override
  State<_MafftPanel> createState() => _MafftPanelState();
}

class _MafftPanelState extends State<_MafftPanel> {
  late final c = ProjectScope.read(context);
  late final _maxIt = TextEditingController(text: '${c.mafft.maxIterate}');
  late final _threads = TextEditingController(text: '${c.mafft.threads}');
  late final _op = TextEditingController(text: '${c.mafft.gapOpen}');
  late final _ep = TextEditingController(text: '${c.mafft.gapExtend}');
  late final _extra = TextEditingController(text: c.mafft.extra);

  @override
  void dispose() {
    for (final x in [_maxIt, _threads, _op, _ep, _extra]) {
      x.dispose();
    }
    super.dispose();
  }

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    final m = c.mafft;
    final cpus = Platform.numberOfProcessors;
    final jobs = math.max(1, (cpus - 2) ~/ math.max(1, m.threads));
    String methodLabel(String v) => switch (v) {
          '--auto' => l10n.mafftAuto,
          '--localpair' => l10n.mafftLocalpair,
          '--globalpair' => l10n.mafftGlobalpair,
          '--genafpair' => l10n.mafftGenafpair,
          '--retree 2 --reorder' => l10n.mafftRetree2,
          _ => l10n.mafftRetree1,
        };
    Widget num(TextEditingController ctl, void Function(String) onChanged, {String? help}) => TextField(
          controller: ctl,
          keyboardType: const TextInputType.numberWithOptions(decimal: true),
          decoration: InputDecoration(helperText: help),
          onChanged: onChanged,
        );
    return SubPanel(
      badge: '5c',
      title: l10n.panelMafft,
      body: l10n.panelMafftBody,
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        AutoGrid(minItemWidth: 260, children: [
          Field(
            label: l10n.mafftStrategy,
            child: DropdownButtonFormField<String>(
              icon: dropdownChevron,
              initialValue: _mafftMethods.contains(m.strategy) ? m.strategy : '--auto',
              isExpanded: true,
              items: [for (final v in _mafftMethods) DropdownMenuItem(value: v, child: Text(methodLabel(v)))],
              onChanged: (v) => c.setMafft(m.copyWith(strategy: v)),
            ),
          ),
          Field(
            label: '${l10n.mafftMaxIterate} (--maxiterate)',
            child: num(_maxIt, (v) => c.setMafft(m.copyWith(maxIterate: int.tryParse(v) ?? 0)), help: l10n.mafftMaxIterateHelp),
          ),
          Field(
            label: '${l10n.mafftThreads} ${l10n.mafftThreadsHelp}',
            child: num(_threads, (v) => c.setMafft(m.copyWith(threads: (int.tryParse(v) ?? 1).clamp(1, cpus))),
                help: l10n.mafftConcurrency('$cpus', '$jobs', '${m.threads}')),
          ),
          Field(label: '${l10n.mafftGapOpening} (--op)', child: num(_op, (v) => c.setMafft(m.copyWith(gapOpen: double.tryParse(v) ?? 1.53)))),
          Field(label: '${l10n.mafftGapExtension} (--ep)', child: num(_ep, (v) => c.setMafft(m.copyWith(gapExtend: double.tryParse(v) ?? 0.0)))),
          Field(
            label: '${l10n.mafftExtra} ${l10n.optional}',
            child: TextField(
              controller: _extra,
              style: const TextStyle(fontFamily: monoFont, fontSize: 13),
              decoration: const InputDecoration(hintText: 'e.g. --nomemsave'),
              onChanged: (v) => c.setMafft(m.copyWith(extra: v)),
            ),
          ),
        ]),
        const SizedBox(height: 12),
        AutoGrid(minItemWidth: 260, children: [
          SwitchBox(
            label: '--adjustdirection',
            detail: l10n.adjustDirectionHelp,
            mono: true,
            value: m.adjustDirection,
            onChanged: (v) => c.setMafft(m.copyWith(adjustDirection: v)),
          ),
          SwitchBox(label: '--reorder', detail: l10n.reorderHelp, mono: true, value: m.reorder, onChanged: (v) => c.setMafft(m.copyWith(reorder: v))),
          SwitchBox(
            label: '--preservecase',
            detail: l10n.preserveCaseHelp,
            mono: true,
            value: m.preserveCase,
            onChanged: (v) => c.setMafft(m.copyWith(preserveCase: v)),
          ),
        ]),
        if (c.seqMode != SeqMode.nucleotide) ...[
          const SizedBox(height: 10),
          Text(l10n.mafftAminoNote, style: TextStyle(fontSize: 12, color: t.textMuted, fontStyle: FontStyle.italic)),
        ],
      ]),
    );
  }
}

// ---------------------------------------------------------------------------
// 5d Run Alignment
// ---------------------------------------------------------------------------

class _RunPanel extends StatelessWidget {
  const _RunPanel();

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    final mafft = c.tools[ExternalTool.mafft];
    final aligned = c.loci.values.any((l) => l.aligned != null);
    return SubPanel(
      badge: '5d',
      title: l10n.panelRun,
      body: l10n.panelRunBody(c.workDir),
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        if (c.toolsChecked && mafft == null)
          Padding(
            padding: const EdgeInsets.only(bottom: 10),
            child: TonePanel(tone: PanelTone.danger, child: Text('${l10n.toolMissing('MAFFT')} — ${ExternalTool.mafft.installHint}')),
          ),
        Wrap(spacing: 10, runSpacing: 10, crossAxisAlignment: WrapCrossAlignment.center, children: [
          PrimaryButton(
            label: l10n.runMafftMode(switch (c.seqMode) {
              SeqMode.nucleotide => l10n.seqNtTitle,
              SeqMode.aminoAcid => l10n.seqAaTitle,
              SeqMode.codon => l10n.seqCodonTitle,
            }),
            icon: FA.play,
            onPressed: c.busy || c.activeMarkers.isEmpty || c.headerBlocked ? null : c.alignAll,
          ),
          if (c.headerBlocked) TagChip(l10n.alignBlockedHeader, icon: FA.lock, color: SplaceTokens.danger),
          if (c.inputLooksAligned && c.seqMode == SeqMode.nucleotide)
            OutlinedButton(onPressed: c.busy ? null : c.useInputAsAlignment, child: Text(l10n.useAsAlignment)),
          Text(l10n.markersToAlign('${c.activeMarkers.length}'), style: TextStyle(fontSize: 13, color: t.textMuted)),
        ]),
        if (c.trimPrompt) ...[
          const SizedBox(height: 16),
          Container(
            padding: const EdgeInsets.all(16),
            decoration: BoxDecoration(
              color: t.accentSoft,
              borderRadius: BorderRadius.circular(14),
              border: Border.all(color: t.accentSoftBorder, width: 1.5),
            ),
            child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
              Row(children: [
                Icon(FA.scissors, size: 22, color: t.accent),
                const SizedBox(width: 14),
                Expanded(
                  child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
                    Text(l10n.trimPromptTitle, style: TextStyle(fontWeight: FontWeight.w600, fontSize: 14, color: t.accentText)),
                    const SizedBox(height: 4),
                    Text(c.seqMode == SeqMode.codon ? l10n.trimPromptCodon : l10n.trimPromptBody, style: TextStyle(fontSize: 13, color: t.textMuted)),
                  ]),
                ),
                const SizedBox(width: 12),
                OutlinedButton(
                  onPressed: () => c.dismissTrimPrompt(skip: true),
                  child: Text(l10n.trimPromptNo),
                ),
                const SizedBox(width: 8),
                FilledButton(onPressed: c.openTrim, child: Text(l10n.trimPromptYes)),
              ]),
              // What trimming does, before choosing.
              const SizedBox(height: 12),
              _trimPreview(c),
            ]),
          ),
        ] else if (aligned && !c.alignPanelVisible('5e')) ...[
          // Trimming skipped (or not chosen yet): it can still be done.
          const SizedBox(height: 12),
          Wrap(spacing: 10, runSpacing: 8, crossAxisAlignment: WrapCrossAlignment.center, children: [
            if (c.trimSkipped) TagChip(l10n.trimSkippedTag, icon: FA.arrowRight, color: t.textMuted),
            ChipButton(label: l10n.configureTrimming, icon: FA.scissors, style: ChipButtonStyle.soft, onPressed: c.openTrim),
          ]),
          const SizedBox(height: 10),
          _trimPreview(c),
        ],
      ]),
    );
  }
}

/// The trimming animation with the method 5e starts with (ClipKIT in the
/// codon-aware mode), shown before trimming is configured.
Widget _trimPreview(ProjectController c) => c.seqMode == SeqMode.codon
    ? TrimExample(method: c.clipkitMode)
    : c.seqMode == SeqMode.nucleotide && _hasCoding(c)
        ? TrimFrameExample(
            method: c.trimal.method,
            gapThreshold: c.trimal.gapThreshold,
            codonMethod: c.trimTool == TrimTool.clipkitCodon ? c.clipkitMode : c.trimal.method,
            codonTool: c.trimTool == TrimTool.clipkitCodon ? 'ClipKIT --codon' : 'trimAl -backtrans',
            keepFrame: c.codonAwareTrimming)
        : TrimExample(method: c.trimal.method, gapThreshold: c.trimal.gapThreshold);

/// Coding markers are present (or the markers are not known yet): trimming
/// can break their reading frame.
bool _hasCoding(ProjectController c) => c.loci.isEmpty || c.loci.values.any((l) => l.included && l.coding);

// ---------------------------------------------------------------------------
// 5e Trimming
// ---------------------------------------------------------------------------

class _TrimPanel extends StatefulWidget {
  const _TrimPanel();

  @override
  State<_TrimPanel> createState() => _TrimPanelState();
}

class _TrimPanelState extends State<_TrimPanel> {
  // Created eagerly: lazy fields first read in dispose() would look up the
  // controller on a deactivated element.
  late final ProjectController c;
  late final TextEditingController _gt, _st, _cons, _w, _extra;

  @override
  void initState() {
    super.initState();
    c = ProjectScope.read(context);
    _gt = TextEditingController(text: '${c.trimal.gapThreshold}');
    _st = TextEditingController(text: c.trimal.similarity?.toString() ?? '');
    _cons = TextEditingController(text: c.trimal.consistency?.toString() ?? '');
    _w = TextEditingController(text: c.trimal.window?.toString() ?? '');
    _extra = TextEditingController(text: c.trimal.extra);
  }

  @override
  void dispose() {
    for (final x in [_gt, _st, _cons, _w, _extra]) {
      x.dispose();
    }
    super.dispose();
  }

  void _update() => c.setTrimal(TrimalOptions(
        method: c.trimal.method,
        gapThreshold: double.tryParse(_gt.text) ?? 0.9,
        similarity: double.tryParse(_st.text),
        consistency: int.tryParse(_cons.text),
        window: int.tryParse(_w.text),
        extra: _extra.text,
      ));

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    final tr = c.trimal;
    final alignedLoci = c.loci.values.where((l) => l.aligned != null).toList();
    final aligned = alignedLoci.isNotEmpty;
    final trimmed = alignedLoci.where((l) => l.trimmed != null).length;
    // The type of the existing alignments decides the tool; before
    // aligning, the mode chosen in 5b.
    final active = !aligned
        ? c.seqMode
        : alignedLoci.any((l) => l.proteinAligned != null && !l.isProtein)
            ? SeqMode.codon
            : alignedLoci.any((l) => l.isProtein)
                ? SeqMode.aminoAcid
                : SeqMode.nucleotide;
    final codon = active == SeqMode.codon;
    final clipkit = c.tools[ExternalTool.clipkit];
    // One tool per run: ClipKIT for codon-aware alignments, or when chosen
    // for a nucleotide alignment (coding markers with --codon); trimAl
    // otherwise.
    final dna = active != SeqMode.aminoAcid;
    // Codon-aware alignments always keep codons: trimAl or ClipKIT.
    final toolNow = codon && c.trimTool == TrimTool.trimalColumns ? TrimTool.trimalCodon : c.trimTool;
    final useClip = dna && toolNow == TrimTool.clipkitCodon;
    final backtrans = dna && toolNow == TrimTool.trimalCodon;
    String modeLabel(String v) => switch (v) {
          '-automated1' => l10n.trimAutomated1,
          '-strict' => l10n.trimStrict,
          '-strictplus' => l10n.trimStrictplus,
          '-gappyout' => l10n.trimGappyout,
          '-nogaps' => l10n.trimNogaps,
          '-noallgaps' => l10n.trimNoallgaps,
          _ => l10n.trimManual,
        };
    String modeDesc(String v) => switch (v) {
          '-automated1' => l10n.trimDescAutomated1,
          '-strict' => l10n.trimDescStrict,
          '-strictplus' => l10n.trimDescStrictplus,
          '-gappyout' => l10n.trimDescGappyout,
          '-nogaps' => l10n.trimDescNogaps,
          '-noallgaps' => l10n.trimDescNoallgaps,
          _ => l10n.trimDescManual,
        };
    String clipDesc(String v) => switch (v) {
          'kpic' => l10n.clipkitDescKpic,
          'kpic-gappy' => l10n.clipkitDescKpicGappy,
          'smart-gap' => l10n.clipkitDescSmartGap,
          _ => l10n.clipkitDescGappy,
        };
    Widget num(TextEditingController ctl, String hint) => TextField(
          controller: ctl,
          keyboardType: const TextInputType.numberWithOptions(decimal: true),
          decoration: InputDecoration(hintText: hint),
          onChanged: (_) => _update(),
        );
    Widget help(String text) => Padding(
          padding: const EdgeInsets.only(top: 6),
          child: Text(text, style: TextStyle(fontSize: 12.5, color: t.textMuted, height: 1.4)),
        );

    final tool = useClip ? 'ClipKIT --codon' : (backtrans ? 'trimAl -backtrans' : 'trimAl');
    Widget clipMode() => Field(
          label: l10n.clipkitMode,
          child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
            DropdownButtonFormField<String>(
              icon: dropdownChevron,
              initialValue: c.clipkitMode,
              isExpanded: true,
              items: [for (final m in _clipkitModes) DropdownMenuItem(value: m, child: Text(m))],
              onChanged: (v) => c.setClipkitMode(v!),
            ),
            help(clipDesc(c.clipkitMode)),
          ]),
        );
    return SubPanel(
      badge: '5e',
      title: l10n.panelTrimParams(tool),
      body: codon ? l10n.trimToolCodon : (active == SeqMode.nucleotide ? l10n.trimToolNt : l10n.trimToolAa),
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        // DNA alignments: which tool trims, one at a time.
        if (dna) ...[
          Field(
            label: l10n.trimToolLabel,
            info: l10n.trimCodonHandlingInfo,
            child: Align(
              alignment: Alignment.centerLeft,
              child: Seg<TrimTool>(
                value: toolNow,
                options: [
                  (TrimTool.trimalCodon, l10n.trimToolTrimalCodon),
                  (TrimTool.clipkitCodon, l10n.trimToolClipkit),
                  if (!codon) (TrimTool.trimalColumns, l10n.trimToolTrimal),
                ],
                onChanged: c.setTrimTool,
              ),
            ),
          ),
          help(switch (toolNow) {
            TrimTool.trimalCodon => l10n.trimToolTrimalCodonBody,
            TrimTool.clipkitCodon => l10n.trimToolClipkitBody,
            TrimTool.trimalColumns => l10n.trimToolTrimalBody,
          }),
          const SizedBox(height: 12),
        ],
        if (useClip) ...[
          AutoGrid(minItemWidth: 300, children: [clipMode()]),
          const SizedBox(height: 12),
          if (_hasCoding(c))
            TrimFrameExample(method: tr.method, gapThreshold: tr.gapThreshold, codonMethod: c.clipkitMode, codonTool: 'ClipKIT --codon', keepFrame: true)
          else
            TrimExample(method: c.clipkitMode),
          if (c.toolsChecked && clipkit == null) ...[
            const SizedBox(height: 12),
            TonePanel(tone: PanelTone.warn, child: Text(l10n.clipkitRequired, style: const TextStyle(fontSize: 13))),
            const SizedBox(height: 12),
            const _PythonCard(),
          ],
        ] else ...[
          AutoGrid(minItemWidth: 300, children: [
            Field(
              label: l10n.trimalMethod,
              child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
                DropdownButtonFormField<String>(
                  icon: dropdownChevron,
                  initialValue: _trimalModes.contains(tr.method) ? tr.method : '-automated1',
                  isExpanded: true,
                  items: [for (final v in _trimalModes) DropdownMenuItem(value: v, child: Text(modeLabel(v)))],
                  onChanged: (v) => c.setTrimal(TrimalOptions(
                    method: v!,
                    gapThreshold: tr.gapThreshold,
                    similarity: tr.similarity,
                    consistency: tr.consistency,
                    window: tr.window,
                    extra: tr.extra,
                  )),
                ),
                help(modeDesc(tr.method)),
              ]),
            ),
            Field(
              label: '${l10n.trimExtra} ${l10n.optional}',
              child: TextField(
                controller: _extra,
                style: const TextStyle(fontFamily: monoFont, fontSize: 13),
                decoration: const InputDecoration(hintText: 'e.g. -terminalonly'),
                onChanged: (_) => _update(),
              ),
            ),
          ]),
          // Only coding genes lose their codons; loci without coding
          // markers need no warning.
          if (active == SeqMode.nucleotide && !backtrans && c.loci.values.any((l) => l.included && l.coding)) ...[
            const SizedBox(height: 12),
            TonePanel(tone: PanelTone.warn, child: MarkupText(l10n.trimFrameWarning, style: const TextStyle(fontSize: 13))),
          ],
          if (tr.manual) ...[
            const SizedBox(height: 12),
            AutoGrid(minItemWidth: 200, children: [
              Field(label: '${l10n.trimGapThreshold} (-gt, 0–1)', child: num(_gt, '0.9')),
              Field(label: '${l10n.trimSimilarity} (-st, 0–1)', child: num(_st, l10n.optionalLabel)),
              Field(label: '${l10n.trimConsistency} (-cons, 0–100)', child: num(_cons, l10n.optionalLabel)),
              Field(label: '${l10n.trimWindow} (-w)', child: num(_w, l10n.optionalLabel)),
            ]),
          ],
          const SizedBox(height: 12),
          if (dna && _hasCoding(c))
            TrimFrameExample(
                method: tr.method, gapThreshold: tr.gapThreshold, codonMethod: tr.method, codonTool: 'trimAl -backtrans', keepFrame: backtrans)
          else
            TrimExample(method: tr.method, gapThreshold: tr.gapThreshold),
        ],
        const SizedBox(height: 14),
        Row(children: [
          PrimaryButton(
            label: l10n.runTrimWith(tool),
            icon: FA.scissors,
            onPressed: c.busy || !aligned || (useClip && c.toolsChecked && clipkit == null) ? null : c.trimAll,
          ),
          const SizedBox(width: 12),
          if (trimmed > 0) TagChip(l10n.trimmedCount('$trimmed', '${alignedLoci.length}'), icon: FA.check, color: SplaceTokens.success),
        ]),
      ]),
    );
  }
}

// ---------------------------------------------------------------------------
// Alignment results
// ---------------------------------------------------------------------------

/// The result cards. With many markers only a page is built (each card
/// holds stats, pills and possibly an alignment viewer), with a search.
class _ResultsList extends StatefulWidget {
  const _ResultsList({required this.results});
  final List<LocusState> results;

  @override
  State<_ResultsList> createState() => _ResultsListState();
}

class _ResultsListState extends State<_ResultsList> {
  static const _page = 25;
  int _shown = _page;
  String _query = '';

  @override
  Widget build(BuildContext context) {
    final l10n = context.l10n;
    final t = context.tokens;
    final all = widget.results;
    if (all.length <= _page) {
      return Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [for (final l in all) _LocusCard(key: ValueKey(l.name), locus: l)]);
    }
    final q = _query.toLowerCase();
    final visible = [for (final l in all) if (q.isEmpty || l.name.toLowerCase().contains(q)) l];
    final page = visible.take(_shown).toList();
    return Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
      Row(children: [
        Expanded(
          child: TextField(
            decoration: InputDecoration(isDense: true, prefixIcon: Icon(FA.magnifyingGlass, size: 13, color: t.textFaint), hintText: l10n.findMarker),
            onChanged: (v) => setState(() {
              _query = v;
              _shown = _page;
            }),
          ),
        ),
      ]),
      const SizedBox(height: 12),
      for (final l in page) _LocusCard(key: ValueKey(l.name), locus: l),
      if (visible.length > page.length)
        Padding(
          padding: const EdgeInsets.only(top: 4),
          child: Row(children: [
            Text(l10n.resultsShowing('${page.length}', '${visible.length}'), style: TextStyle(fontSize: 13, color: t.textMuted)),
            const Spacer(),
            ChipButton(label: l10n.showMoreN('$_page'), onPressed: () => setState(() => _shown += _page)),
            const SizedBox(width: 8),
            ChipButton(label: l10n.showAll, onPressed: () => setState(() => _shown = visible.length)),
          ]),
        ),
      if (visible.isEmpty) Text(l10n.noMatches, style: TextStyle(fontStyle: FontStyle.italic, color: t.textFaint)),
    ]);
  }
}

class _LocusCard extends StatefulWidget {
  const _LocusCard({super.key, required this.locus});
  final LocusState locus;

  @override
  State<_LocusCard> createState() => _LocusCardState();
}

class _LocusCardState extends State<_LocusCard> {
  bool _open = false;
  String _tab = 'aligned';

  /// Trimmed tab: the whole MAFFT alignment with the removed columns in
  /// grey (or only the kept columns).
  bool _showRemoved = true;

  /// Where the last view was scrolled to: the next tab opens there.
  AlignmentAnchor? _anchor;

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    final l = widget.locus;
    final lengths = l.raw.map((r) => r.length).toList();
    final mean = lengths.isEmpty ? 0 : lengths.fold<int>(0, (a, b) => a + b) ~/ lengths.length;
    final done = l.alignStatus == StepStatus.done;
    final alnLen = l.aligned?.first.length;
    final trimLen = l.trimmed?.first.length;
    final unit = l.isProtein ? 'aa' : 'bp';

    Widget stat(String label, String main, String? sub) => Expanded(
          child: Column(children: [
            Text(label, style: TextStyle(fontSize: 13, color: t.textMuted)),
            const SizedBox(height: 2),
            Text(main, style: TextStyle(fontSize: 13, fontWeight: FontWeight.w700, color: t.text)),
            if (sub != null) Text(sub, style: TextStyle(fontSize: 12, color: t.textFaint)),
          ]),
        );

    String codons(int len) => l.coding && !l.isProtein ? ' · ${len ~/ 3} codons' : '';
    final tabs = <String, List<dynamic>?>{'raw': l.raw, 'aligned': l.aligned, 'trimmed': l.trimmed, if (l.proteinAligned != null && !l.isProtein) 'protein': l.proteinAligned};
    if (tabs[_tab] == null) _tab = tabs.entries.lastWhere((e) => e.value != null, orElse: () => tabs.entries.first).key;
    // Columns of the MAFFT alignment removed by trimming (grey in the
    // trimmed tab), when the kept columns are known.
    final cols = l.trimmedColumns;
    List<bool>? mask;
    if (l.trimmed != null && cols != null && l.aligned != null && l.aligned!.length == l.trimmed!.length) {
      mask = List<bool>.filled(l.aligned!.first.length, true);
      for (final i in cols) {
        if (i < mask.length) mask[i] = false;
      }
    }
    final shown = switch (_tab) {
      'raw' => l.raw,
      'trimmed' => mask != null && _showRemoved ? l.aligned : l.trimmed,
      'protein' => l.proteinAligned,
      _ => l.aligned,
    };

    return Opacity(
      opacity: l.included ? 1 : 0.55,
      child: Container(
        margin: const EdgeInsets.only(bottom: 12),
        padding: const EdgeInsets.all(14),
        decoration: BoxDecoration(color: t.surface, borderRadius: BorderRadius.circular(12), border: Border.all(color: t.border), boxShadow: shadowSm),
        child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
          Row(children: [
            Container(
              width: 26,
              height: 26,
              decoration: BoxDecoration(color: done ? t.accent : t.surface3, shape: BoxShape.circle),
              child: l.alignStatus == StepStatus.running || l.trimStatus == StepStatus.running
                  ? const Padding(padding: EdgeInsets.all(6), child: CircularProgressIndicator(strokeWidth: 2, color: Colors.white))
                  : Icon(l.alignStatus == StepStatus.error ? FA.xmark : FA.check, size: 12, color: done ? Colors.white : t.textFaint),
            ),
            const SizedBox(width: 10),
            Text(l.name, style: TextStyle(fontSize: 14, fontWeight: FontWeight.w600, color: t.text)),
            const SizedBox(width: 10),
            Expanded(
              child: Text(
                alnLen == null ? '${l.raw.length} seq · ~$mean bp' : '${l.aligned!.length} seq · $alnLen $unit${codons(alnLen)}',
                overflow: TextOverflow.ellipsis,
                style: TextStyle(fontSize: 13, color: t.textMuted),
              ),
            ),
            if (!l.isProtein) ...[
              Text(l10n.colCoding, style: TextStyle(fontSize: 12, color: t.textMuted)),
              Transform.scale(scale: 0.75, child: Switch(value: l.coding, onChanged: (v) => c.setLocusCoding(l.name, coding: v))),
              // The frame can be changed only for DNA alignments not yet
              // trimmed; codon-aware ones always start at base 1.
              if (l.coding && l.trimmed == null && l.proteinAligned == null) ...[
                const SizedBox(width: 6),
                SplaceTooltip(
                  message: l10n.frameInfo,
                  child: Row(mainAxisSize: MainAxisSize.min, children: [
                    Text(l10n.colFrame, style: TextStyle(fontSize: 12, color: t.textMuted)),
                    const SizedBox(width: 4),
                    Icon(FA.circleInfo, size: 11, color: t.textFaint),
                  ]),
                ),
                const SizedBox(width: 6),
                DropdownButton<int>(
                  icon: dropdownChevron,
                  value: l.frame,
                  isDense: true,
                  underline: const SizedBox.shrink(),
                  items: [for (final f in [0, 1, 2]) DropdownMenuItem(value: f, child: Text(l10n.frameOption('${f + 1}')))],
                  onChanged: (v) => c.setLocusCoding(l.name, frame: v),
                ),
              ],
            ],
            const SizedBox(width: 10),
            SplaceTooltip(
              message: l10n.colInclude,
              child: SizedBox(width: 22, height: 22, child: Checkbox(value: l.included, onChanged: (v) => c.setLocusIncluded(l.name, v!))),
            ),
            const SizedBox(width: 6),
            IconBtn(
              icon: _open ? FA.chevronUp : FA.chevronDown,
              tooltip: l10n.viewAlignment,
              size: 28,
              bordered: true,
              onPressed: () => setState(() => _open = !_open),
            ),
          ]),
          const SizedBox(height: 10),
          // The steps this locus went through, in order, in the app's blue.
          _Pipeline(steps: [
            for (final s in l.sources) (s, s == 'FASTA' ? FA.fileLines : FA.dna),
            for (final p in l.pipeline) (p, null),
            if (l.trimMethod != null) (l.trimMethod!, FA.scissors),
          ]),
          const SizedBox(height: 6),
          Wrap(spacing: 6, runSpacing: 6, children: [
            if (l.coding && !l.isProtein && !(l.trimmed != null && l.trimmedFrameLost))
              SplaceTooltip(message: l10n.frameInfo, child: TagChip(l10n.frameTag('${l.alignmentFrame + 1}'), icon: FA.circleInfo, color: t.textMuted)),
            if (l.trimmed != null && l.trimmedFrameLost)
              SplaceTooltip(
                message: l10n.frameLostInfo,
                child: TagChip(l10n.frameLostTag, icon: FA.triangleExclamation, color: SplaceTokens.danger),
              ),
            if (l.frameShifted.isNotEmpty)
              SplaceTooltip(
                message: '${l10n.frameShiftInfo}\n\n${l.frameShifted.take(12).join('\n')}${l.frameShifted.length > 12 ? '\n…' : ''}',
                child: TagChip(l10n.frameShiftTag('${l.frameShifted.length}'), icon: FA.triangleExclamation, color: SplaceTokens.warn),
              ),
            if (l.alignStatus == StepStatus.skipped && l.message != null) TagChip(l.message!, icon: FA.circleInfo, color: SplaceTokens.warn),
            if (l.alignStatus == StepStatus.error) TagChip(l.message ?? l10n.statusError, icon: FA.triangleExclamation, color: SplaceTokens.danger),
            if (l.trimStatus == StepStatus.error) TagChip('${l10n.trimming}: ${l.message ?? l10n.statusError}', icon: FA.triangleExclamation, color: SplaceTokens.danger),
          ]),
          const SizedBox(height: 12),
          Row(children: [
            stat(l10n.statRaw, '${l.raw.length} seq · $mean bp', null),
            stat(l10n.statAligned, alnLen == null ? '—' : '${l.aligned!.length} seq · $alnLen $unit${codons(alnLen)}',
                l.stats != null && l.trimmed == null ? '${pct(l.stats!.gapFraction)}% gaps · ${l.stats!.informativeSites} PI' : null),
            stat(l10n.statTrimmed, trimLen == null ? l10n.noTrimming : '${l.trimmed!.length} seq · $trimLen $unit${codons(trimLen)}',
                l.stats != null && l.trimmed != null ? '${pct(l.stats!.gapFraction)}% gaps · ${l.stats!.informativeSites} PI' : null),
          ]),
          if (_open) ...[
            const SizedBox(height: 12),
            Wrap(spacing: 12, runSpacing: 8, crossAxisAlignment: WrapCrossAlignment.center, children: [
              Seg<String>(
                value: _tab,
                options: [
                  for (final e in tabs.entries)
                    if (e.value != null)
                      (
                        e.key,
                        switch (e.key) {
                          'raw' => l10n.statRaw,
                          'trimmed' => l10n.statTrimmed,
                          'protein' => l10n.tabProtein,
                          _ => l10n.statAligned,
                        }
                      ),
                ],
                onChanged: (v) => setState(() => _tab = v),
              ),
              if (_tab == 'trimmed' && mask != null)
                SwitchBox(label: l10n.showRemovedColumns, value: _showRemoved, onChanged: (v) => setState(() => _showRemoved = v)),
              // The trimmer kept every column: nothing to grey.
              if ((_tab == 'trimmed' || _tab == 'aligned') && mask != null && !mask.contains(true))
                TagChip(l10n.trimNothingRemoved(l.trimMethod ?? '', '${mask.length}'), icon: FA.circleInfo, color: t.textMuted),
            ]),
            const SizedBox(height: 10),
            if (shown != null) _viewer(c, l, shown, mask),
          ],
        ]),
      ),
    );
  }


  Widget _viewer(ProjectController c, LocusState l, List<SeqRecord> shown, List<bool>? mask) {
    final full = _tab == 'trimmed' && mask != null && _showRemoved;
    final compact = _tab == 'trimmed' && !full;
    final dna = !l.isProtein && _tab != 'protein';
    // Codon positions of the compact trimmed view are known from the column
    // map, or because whole codons were kept.
    final positions = l.coding && dna && _tab != 'raw' && !(compact && l.trimmedFrameLost && l.trimmedColumns == null);
    return AlignmentViewer(
      key: ValueKey('${l.name}-$_tab-$full'),
      records: shown,
      protein: _tab == 'protein' || (l.isProtein && _tab != 'raw'),
      // Only protein-coding genes are translated.
      translatable: l.coding && dna && _tab != 'raw' && !(compact && l.trimmedFrameLost),
      frame: compact ? l.alignmentFrame : l.frame,
      geneticCode: c.geneticCode,
      coding: positions,
      columnMap: compact ? l.trimmedColumns : null,
      mapFrame: l.frame,
      // The removed columns are grey in the MAFFT alignment too.
      removed: full || _tab == 'aligned' ? mask : null,
      space: switch (_tab) {
        'raw' => ViewSpace.raw,
        'protein' => ViewSpace.protein,
        _ => ViewSpace.mafft,
      },
      anchor: _anchor,
      onAnchor: (a) => _anchor = a,
    );
  }
}

/// Steps of a locus as one connected row: numbered tags in the style of
/// the step 4 buttons, joined by arrows.
class _Pipeline extends StatelessWidget {
  const _Pipeline({required this.steps});
  final List<(String, IconData?)> steps;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    return Wrap(spacing: 4, runSpacing: 6, crossAxisAlignment: WrapCrossAlignment.center, children: [
      for (final (i, (label, icon)) in steps.indexed) ...[
        if (i > 0) Icon(FA.arrowRight, size: 9, color: t.accent.withValues(alpha: 0.55)),
        TagChip(
          label,
          icon: icon,
          leading: Text('${i + 1}', style: TextStyle(fontSize: 12, fontWeight: FontWeight.w800, color: t.accent)),
        ),
      ],
    ]);
  }
}
