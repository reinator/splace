import 'package:flutter/material.dart';

/// A small part of HTML as text spans, for texts people edit by hand (the
/// tooltips of the example buttons): `<b>`/`<strong>` bold, `<i>`/`<em>`
/// italic, `<u>` underline, `<br>` a new line, `&amp; &lt; &gt; &nbsp;
/// &quot;`. Other tags are left out, their text kept.
InlineSpan htmlSpan(String html, {TextStyle? style}) {
  final spans = <InlineSpan>[];
  var bold = 0, italic = 0, underline = 0;
  String decode(String s) => s
      .replaceAll('&nbsp;', ' ')
      .replaceAll('&lt;', '<')
      .replaceAll('&gt;', '>')
      .replaceAll('&quot;', '"')
      .replaceAll('&#39;', "'")
      .replaceAll('&amp;', '&');
  void text(String s) {
    if (s.isEmpty) return;
    spans.add(TextSpan(
      text: decode(s.replaceAll(RegExp(r'\s*\n\s*'), ' ')),
      style: TextStyle(
        fontWeight: bold > 0 ? FontWeight.w700 : null,
        fontStyle: italic > 0 ? FontStyle.italic : null,
        decoration: underline > 0 ? TextDecoration.underline : null,
      ),
    ));
  }

  var at = 0;
  for (final m in RegExp(r'<\s*(/?)\s*([a-zA-Z]+)[^>]*?(/?)\s*>').allMatches(html)) {
    text(html.substring(at, m.start));
    at = m.end;
    final closing = m.group(1) == '/';
    final d = closing ? -1 : 1;
    switch (m.group(2)!.toLowerCase()) {
      case 'b' || 'strong':
        bold = (bold + d).clamp(0, 99);
      case 'i' || 'em':
        italic = (italic + d).clamp(0, 99);
      case 'u':
        underline = (underline + d).clamp(0, 99);
      case 'br':
        spans.add(const TextSpan(text: '\n'));
      case 'p' when closing:
        spans.add(const TextSpan(text: '\n\n'));
    }
  }
  text(html.substring(at));
  return TextSpan(style: style, children: spans);
}
