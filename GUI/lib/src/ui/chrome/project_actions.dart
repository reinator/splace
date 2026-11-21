import 'package:file_selector/file_selector.dart';
import 'package:flutter/material.dart';
import 'package:path/path.dart' as p;

import '../../state/project_controller.dart';
import '../common.dart';
import '../fa.dart';
import '../widgets/dialogs.dart';
import '../widgets/toast.dart';

const _projectType = XTypeGroup(label: 'SPLACE project', extensions: [projectExtension]);

/// "Save project": to the file it was last saved to or opened from, or
/// (first time, or [saveAs]) to a file the user picks.
Future<void> saveProjectAction(BuildContext context, {bool saveAs = false}) async {
  final c = ProjectScope.read(context);
  final l10n = context.l10n;
  if (c.busy) return;
  var target = saveAs ? null : c.projectFile;
  if (target == null) {
    final stamp = DateTime.now().toIso8601String().substring(0, 10);
    final loc = await getSaveLocation(
      suggestedName: c.projectFile == null ? 'splace_project_$stamp.$projectExtension' : p.basename(c.projectFile!),
      acceptedTypeGroups: const [_projectType],
    );
    if (loc == null) return;
    target = loc.path;
  }
  final ok = await c.saveProject(target, appVersion: appVersion);
  if (!context.mounted) return;
  final file = c.projectFile ?? target;
  showSplaceToast(
    context,
    title: ok ? l10n.projectSavedTitle : l10n.projectSaveFailed,
    detail: p.basename(file),
    path: p.dirname(file),
    error: !ok,
    action: ok ? (label: l10n.openFolder, onPressed: () => c.openFolder(p.dirname(file))) : null,
  );
}

/// "Open project": asks first when this project has data (it is closed;
/// its files stay in its folder).
Future<void> openProjectAction(BuildContext context) async {
  final c = ProjectScope.read(context);
  final l10n = context.l10n;
  if (c.busy) return;
  final f = await openFile(acceptedTypeGroups: const [_projectType]);
  if (f == null || !context.mounted) return;
  final hasData = c.records.isNotEmpty || c.inputs.isNotEmpty;
  if (hasData) {
    final ok = await confirmDialog(
      context,
      title: l10n.openProjectTitle,
      body: l10n.openProjectBody(p.basename(f.path)),
      body2: l10n.newProjectKept(c.workDir),
      confirm: l10n.openProject,
      danger: false,
      icon: FA.folderOpen,
    );
    if (!ok) return;
  }
  await c.requestOpenProject(f.path);
}
