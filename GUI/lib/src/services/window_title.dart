import 'package:flutter/services.dart';

const _channel = MethodChannel('splace/window');
String? _last;

/// Sets the title bar of the window (Windows, macOS and Linux runners):
/// [title], the project [file] (macOS shows its icon) and whether it has
/// unsaved changes ([edited], the dot of the macOS close button). Sent only
/// when something changed; nothing happens where there is no runner
/// (tests).
Future<void> setWindowTitle(String title, {String? file, bool edited = false}) async {
  final key = '$title\u0000$file\u0000$edited';
  if (key == _last) return;
  _last = key;
  try {
    await _channel.invokeMethod<void>('setTitle', {'title': title, 'file': file, 'edited': edited});
  } on MissingPluginException {
    // No runner (widget tests).
  } on PlatformException {
    // An older runner without the channel.
  }
}
