import 'dart:io';

import 'package:file_selector/file_selector.dart';

import 'package:flutter/material.dart';

import '../../state/project_controller.dart';
import '../common.dart';
import '../fa.dart';
import '../theme.dart';
import '../widgets/kit.dart';
import '../widgets/tree_viewer.dart';

/// Step 9 — the inferred trees: pick one (ML tree, consensus, RAxML-NG
/// support variants, gene trees) and look at it in the tree viewer.
class TreeSection extends StatelessWidget {
  const TreeSection({super.key, required this.active});
  final bool active;

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    final trees = c.viewableTrees;
    if (trees.isEmpty) return const SizedBox.shrink();
    // The consensus tree by default (of the newest run that has one, runs
    // are newest first); the
    // ML trees and the others when the user picks them.
    final chosen = trees.where((x) => x.path == c.treeToView).firstOrNull ??
        trees.where((x) => x.kind == 'consensus').firstOrNull ??
        trees.first;
    String label(({TreeRun run, String kind, String? gene, String path}) x) {
      final when = x.run.started.toIso8601String().substring(11, 16);
      final what = switch (x.kind) {
        'consensus' => l10n.treeKindConsensus,
        'treeFBP' => l10n.treeKindFbp,
        'treeTBE' => l10n.treeKindTbe,
        'treeRBS' => l10n.treeKindRbs,
        'treeSH' => l10n.treeKindSh,
        'treeEBG' => l10n.treeKindEbg,
        'bestTree' => l10n.treeKindBest,
        'gene' => l10n.treeKindGene(x.gene ?? ''),
        _ => l10n.treeKindMl,
      };
      return '${x.run.engine} · $when · $what';
    }

    final file = File(chosen.path);
    final newick = file.existsSync() ? file.readAsStringSync().trim() : '';
    final bootstrap = chosen.run.engine != 'IQ-TREE' || (chosen.run.command?.contains(' -b ') ?? false);
    return StepSection(
      badge: '9',
      icon: FA.diagramProject,
      title: l10n.step9,
      description: l10n.treeViewIntro,
      active: active,
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        SubAnchor(id: '9a', child: SubHeading(id: '9a', title: l10n.subTree)),
        Row(children: [
          Expanded(
            child: Field(
              label: l10n.treeChoose,
              child: DropdownButtonFormField<String>(
                key: ValueKey(chosen.path),
                icon: dropdownChevron,
                initialValue: chosen.path,
                isExpanded: true,
                items: [for (final x in trees) DropdownMenuItem(value: x.path, child: Text(label(x), overflow: TextOverflow.ellipsis))],
                onChanged: c.viewTree,
              ),
            ),
          ),
          const SizedBox(width: 12),
          Padding(
            padding: const EdgeInsets.only(top: 22),
            child: ChipButton(label: l10n.openFolder, onPressed: () => c.openFolder(chosen.run.directory)),
          ),
        ]),
        const SizedBox(height: 6),
        Text(chosen.path, style: TextStyle(fontFamily: monoFont, fontSize: 11.5, color: t.textFaint)),
        const SizedBox(height: 14),
        if (newick.isEmpty)
          Text(l10n.treeMissing, style: TextStyle(color: t.textFaint))
        else
          TreeViewer(
            key: ValueKey(chosen.path),
            newick: newick,
            outgroup: chosen.kind == 'gene' ? const [] : c.concat.outgroups,
            fileName: chosen.gene == null ? 'splace_${chosen.run.engine.toLowerCase().replaceAll('-', '')}_${chosen.kind}' : 'splace_${chosen.gene}_tree',
            bootstrapScale: bootstrap,
            tipInfo: c.tipInfo(),
            title: label(chosen),
            compareWith: [
              for (final x in trees)
                if (x.path != chosen.path)
                  TreeOption(
                    label: label(x),
                    load: () => File(x.path).readAsStringSync().trim(),
                    outgroup: x.kind == 'gene' ? const [] : c.concat.outgroups,
                  ),
            ],
          ),
      ]),
    );
  }
}

/// "Save everything": the whole project folder (sequences, alignments,
/// matrices, trees and logs) as one ZIP, after the trees.
class ResultsZipCard extends StatelessWidget {
  const ResultsZipCard({super.key});

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    return Container(
      margin: const EdgeInsets.only(bottom: 24),
      padding: const EdgeInsets.all(18),
      decoration: BoxDecoration(
        color: t.surface,
        borderRadius: BorderRadius.circular(16),
        border: Border.all(color: t.accentSoftBorder, width: 1.5),
        boxShadow: shadowSm,
      ),
      child: Wrap(spacing: 16, runSpacing: 12, crossAxisAlignment: WrapCrossAlignment.center, alignment: WrapAlignment.spaceBetween, children: [
        ConstrainedBox(
          constraints: const BoxConstraints(maxWidth: 640),
          child: Row(crossAxisAlignment: CrossAxisAlignment.start, children: [
            const Padding(padding: EdgeInsets.only(top: 2), child: StepBadge('9b', size: 26)),
            const SizedBox(width: 12),
            Expanded(
              child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
                Text(l10n.zipTitle, style: TextStyle(fontSize: 16, fontWeight: FontWeight.w700, color: t.text)),
                const SizedBox(height: 3),
                Text(l10n.zipBody, style: TextStyle(fontSize: 13, color: t.textMuted, height: 1.4)),
              ]),
            ),
          ]),
        ),
        Row(mainAxisSize: MainAxisSize.min, children: [
          ChipButton(label: l10n.openFolder, icon: FA.folderOpen, onPressed: () => c.openFolder(c.workDir)),
          const SizedBox(width: 8),
          PrimaryButton(
            label: l10n.saveZip,
            icon: FA.fileZipper,
            onPressed: c.busy
                ? null
                : () async {
                    final loc = await getSaveLocation(suggestedName: 'splace_results.zip', acceptedTypeGroups: const [XTypeGroup(label: 'ZIP', extensions: ['zip'])]);
                    if (loc != null) await c.saveZip(loc.path.endsWith('.zip') ? loc.path : '${loc.path}.zip');
                  },
          ),
        ]),
      ]),
    );
  }
}
