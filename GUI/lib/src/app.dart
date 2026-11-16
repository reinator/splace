import 'package:flutter/material.dart';

import '../l10n/app_localizations.dart';
import 'state/project_controller.dart';
import 'ui/common.dart';
import 'ui/shell.dart';
import 'ui/theme.dart';

class SplaceApp extends StatefulWidget {
  const SplaceApp({super.key, this.controller});

  /// Injected in tests; otherwise the app owns its controller.
  final ProjectController? controller;

  @override
  State<SplaceApp> createState() => _SplaceAppState();
}

class _SplaceAppState extends State<SplaceApp> {
  late final controller = widget.controller ?? ProjectController(settingsFile: ProjectController.defaultSettingsFile);

  @override
  void dispose() {
    if (widget.controller == null) controller.dispose();
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
            child: child!,
          ),
          locale: controller.locale,
          localizationsDelegates: AppLocalizations.localizationsDelegates,
          supportedLocales: AppLocalizations.supportedLocales,
          home: const SplaceShell(),
        ),
      ),
    );
  }
}
