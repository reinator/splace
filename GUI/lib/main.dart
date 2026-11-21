import 'package:flutter/foundation.dart';
import 'package:flutter/material.dart';

import 'src/app.dart';
import 'src/services/file_association.dart';

/// `splace project.splc` (or a double-click on a project file) opens it.
void main(List<String> args) {
  final file = args.where((a) => a.toLowerCase().endsWith('.splc')).firstOrNull;
  // Windows: double-clicking a .splc opens this SPLACE (the released
  // program, not a development build).
  if (kReleaseMode) registerProjectFileType();
  runApp(SplaceApp(openFile: file));
}
