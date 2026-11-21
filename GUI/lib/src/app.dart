import 'package:flutter/material.dart';

import '../l10n/app_localizations.dart';
import 'services/window_title.dart';
import 'state/project_controller.dart';
import 'ui/common.dart';
import 'ui/shell.dart';
import 'ui/theme.dart';

class SplaceApp extends StatefulWidget {
  const SplaceApp({super.key, this.controller, this.openFile});

  /// Injected in tests; otherwise the app owns its controller.
  final ProjectController? controller;

  /// Project file (`.splc`) to open at start.
  final String? openFile;

  @override
  State<SplaceApp> createState() => _SplaceAppState();
}

class _SplaceAppState extends State<SplaceApp> {
  late ProjectController controller = (widget.controller ?? ProjectController(settingsFile: ProjectController.defaultSettingsFile))
    ..newProjectHandler = _newProject
    ..openProjectHandler = _openProject;

  @override
  void initState() {
    super.initState();
    // A project file given on the command line (double-click on a .splc).
    final file = widget.openFile;
    if (file != null) WidgetsBinding.instance.addPostFrameCallback((_) => _openProject(file));
  }

  /// "New project": a fresh controller (empty analysis, new project folder)
  /// with the tools already found, so nothing has to be detected again; the
  /// shell is rebuilt for it (its key is the controller).
  void _newProject() => _replace(_fresh());

  /// A fresh controller with the same settings and the tools already found,
  /// so nothing has to be detected again.
  ProjectController _fresh() {
    final old = controller;
    return ProjectController(settingsFile: old.settingsFile, locator: old.locator, autoDetectTools: false)
      ..tools = Map.of(old.tools)
      ..toolsChecked = old.toolsChecked
      ..python = old.python
      ..splashVisible = false
      ..newProjectHandler = _newProject
      ..openProjectHandler = _openProject;
  }

  /// The shell is rebuilt for [next] (its key is the controller).
  void _replace(ProjectController next) {
    final old = controller;
    setState(() => controller = next);
    WidgetsBinding.instance.addPostFrameCallback((_) {
      if (widget.controller == null || !identical(old, widget.controller)) old.dispose();
    });
  }

  /// "Open project": the file is read into a fresh controller, which
  /// replaces this one once it is ready; if it cannot be read, the current
  /// project stays and the error goes to its log.
  Future<void> _openProject(String file) async {
    final old = controller;
    final next = _fresh();
    // Tools of a project opened at start are detected as usual.
    if (!old.toolsChecked) next.detectTools();
    old.showBusy('Opening project');
    final ok = await next.openProject(file);
    old.showBusy(null);
    if (ok) {
      _replace(next);
    } else {
      for (final e in next.logs.where((e) => e.level == LogLevel.error)) {
        old.log(e.message, error: true);
      }
      next.dispose();
    }
  }

  @override
  void dispose() {
    if (!identical(controller, widget.controller)) controller.dispose();
    super.dispose();
  }


  @override
  Widget build(BuildContext context) {
    return ProjectScope(
      controller: controller,
      child: ListenableBuilder(
        listenable: controller,
        builder: (context, _) => MaterialApp(
          title: 'SPLACE',
          debugShowCheckedModeBanner: false,
          theme: buildSplaceTheme(Brightness.light, compact: controller.compact),
          darkTheme: buildSplaceTheme(Brightness.dark, compact: controller.compact),
          themeMode: controller.themeMode,
          // Text size slider of the settings popover (16 px = 100%).
          builder: (context, child) => MediaQuery(
            data: MediaQuery.of(context).copyWith(textScaler: TextScaler.linear(controller.fontSize / 16)),
            child: _WindowTitle(controller: controller, child: child!),
          ),
          locale: controller.locale,
          localizationsDelegates: AppLocalizations.localizationsDelegates,
          supportedLocales: AppLocalizations.supportedLocales,
          home: SplaceShell(key: ObjectKey(controller)),
        ),
      ),
    );
  }
}

/// The title bar of the window: the project file (full path) and whether
/// it is saved, e.g. `D:\\data\\fish.splc · saved — SPLACE` or
/// `● D:\\data\\fish.splc · unsaved changes — SPLACE`.
class _WindowTitle extends StatelessWidget {
  const _WindowTitle({required this.controller, required this.child});
  final ProjectController controller;
  final Widget child;

  @override
  Widget build(BuildContext context) {
    final l10n = context.l10n;
    return ListenableBuilder(
      listenable: Listenable.merge([controller, controller.unsaved]),
      builder: (context, child) {
        final file = controller.projectFile;
        final edited = controller.unsaved.value;
        final title = file == null
            ? (edited ? '● ${l10n.windowTitleNew} — SPLACE' : 'SPLACE')
            : edited
                ? '● $file · ${l10n.windowTitleUnsaved} — SPLACE'
                : '$file · ${l10n.windowTitleSaved} — SPLACE';
        setWindowTitle(title, file: file, edited: edited);
        return child!;
      },
      child: child,
    );
  }
}
