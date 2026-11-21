import 'dart:io';

import 'package:file_selector/file_selector.dart';
import 'package:flutter/material.dart';

import '../../services/tool_locator.dart';
import '../common.dart';
import '../fa.dart';
import '../theme.dart';
import '../widgets/dialogs.dart';
import '../widgets/kit.dart';


/// External tools, threads and project folder, in a modal (the genetic
/// code is chosen in step 5, where it is used),
/// styled like the HTML "FASTA Header Builder" dialog.
Future<void> showToolsDialog(BuildContext context) => showDialog(
      context: context,
      builder: (_) => ProjectScope(controller: ProjectScope.read(context), child: const _ToolsDialog()),
    );

class _ToolsDialog extends StatelessWidget {
  const _ToolsDialog();

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final t = context.tokens;
    final l10n = context.l10n;
    return SplaceDialog(
      icon: FA.screwdriverWrench,
      title: l10n.toolsTitle,
      subtitle: l10n.toolsBody,
      actions: [FilledButton(onPressed: () => Navigator.pop(context), child: Text(l10n.done))],
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        Kicker(l10n.tools),
        const SizedBox(height: 8),
        for (final tool in ExternalTool.values) _ToolRow(tool: tool),
        const SizedBox(height: 8),
        Align(
          alignment: Alignment.centerLeft,
          child: ChipButton(label: l10n.redetect, icon: FA.arrowsRotate, onPressed: c.detectTools),
        ),
        const SizedBox(height: 20),
        Kicker(l10n.analysisSettings),
        const SizedBox(height: 10),
        AutoGrid(minItemWidth: 320, children: [
          Field(
            label: l10n.threads('${c.threads}'),
            child: Slider(
              value: c.threads.toDouble(),
              min: 1,
              max: Platform.numberOfProcessors.toDouble(),
              divisions: Platform.numberOfProcessors > 1 ? Platform.numberOfProcessors - 1 : 1,
              label: '${c.threads}',
              onChanged: (v) => c.setThreads(v.round()),
            ),
          ),
        ]),
        const SizedBox(height: 16),
        Field(
          label: l10n.projectFolder,
          child: Row(children: [
            Expanded(
              child: Container(
                padding: const EdgeInsets.symmetric(horizontal: 12, vertical: 10),
                decoration: BoxDecoration(color: t.surface2, borderRadius: BorderRadius.circular(8), border: Border.all(color: t.border)),
                child: SelectableText(c.workDir, style: TextStyle(fontFamily: monoFont, fontSize: 12, color: t.text)),
              ),
            ),
            const SizedBox(width: 8),
            ChipButton(
              label: l10n.change,
              icon: FA.folder,
              onPressed: () async {
                final dir = await getDirectoryPath();
                if (dir != null) c.workDir = dir;
              },
            ),
            const SizedBox(width: 8),
            ChipButton(label: l10n.openFolder, icon: FA.folderOpen, onPressed: () => c.openFolder(c.dir(''))),
          ]),
        ),
      ]),
    );
  }
}

class _ToolRow extends StatefulWidget {
  const _ToolRow({required this.tool});
  final ExternalTool tool;

  @override
  State<_ToolRow> createState() => _ToolRowState();
}

class _ToolRowState extends State<_ToolRow> {
  final _path = TextEditingController();

  @override
  void dispose() {
    _path.dispose();
    super.dispose();
  }

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final t = context.tokens;
    final l10n = context.l10n;
    final r = c.tools[widget.tool];
    final checked = c.toolsChecked;
    return Container(
      margin: const EdgeInsets.only(bottom: 10),
      padding: const EdgeInsets.all(12),
      decoration: BoxDecoration(
        color: r != null ? t.accentSoft : t.surface,
        borderRadius: BorderRadius.circular(12),
        border: Border.all(color: r != null ? t.accentSoftBorder : t.border),
      ),
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        Row(children: [
          Icon(r != null ? FA.circleCheck : (checked ? FA.circleXmark : FA.spinner),
              size: 15, color: r != null ? SplaceTokens.success : (checked ? SplaceTokens.danger : t.textFaint)),
          const SizedBox(width: 10),
          Text(widget.tool.label, style: TextStyle(fontSize: 14, fontWeight: FontWeight.w600, color: t.text)),
          if (r?.version != null) ...[const SizedBox(width: 8), TagChip('v${r!.version}')],
          const SizedBox(width: 10),
          Expanded(
            child: Text(
              r?.description ?? (checked ? '${l10n.notFound} — ${widget.tool.installHint}' : '…'),
              overflow: TextOverflow.ellipsis,
              style: TextStyle(fontFamily: monoFont, fontSize: 12, color: t.textMuted),
            ),
          ),
        ]),
        const SizedBox(height: 8),
        Row(children: [
          Expanded(
            child: TextField(
              controller: _path,
              style: const TextStyle(fontSize: 13),
              decoration: InputDecoration(hintText: l10n.overridePath),
              onSubmitted: (v) => c.setToolOverride(widget.tool, v.trim()),
            ),
          ),
          const SizedBox(width: 8),
          ChipButton(
            label: l10n.browse,
            icon: FA.folderOpen,
            onPressed: () async {
              final f = await openFile();
              if (f != null) {
                _path.text = f.path;
                c.setToolOverride(widget.tool, f.path);
              }
            },
          ),
        ]),
      ]),
    );
  }
}
