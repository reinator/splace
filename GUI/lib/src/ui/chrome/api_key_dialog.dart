import 'dart:math' as math;

import 'package:flutter/material.dart';
import 'package:flutter/services.dart';
import 'package:intl/intl.dart';
import 'package:url_launcher/url_launcher.dart';

import '../../services/ncbi.dart';
import '../common.dart';
import '../fa.dart';
import '../theme.dart';
import '../widgets/dialogs.dart';
import '../widgets/kit.dart';

/// "NCBI API key" modal (key button in the top bar): status, the request
/// rate with and without a key, how to get one and the key field.
Future<void> showApiKeyDialog(BuildContext context) => showDialog(
      context: context,
      builder: (_) => ProjectScope(controller: ProjectScope.read(context), child: const _ApiKeyDialog()),
    );

class _ApiKeyDialog extends StatefulWidget {
  const _ApiKeyDialog();

  @override
  State<_ApiKeyDialog> createState() => _ApiKeyDialogState();
}

class _ApiKeyDialogState extends State<_ApiKeyDialog> {
  late final _key = TextEditingController(text: ProjectScope.read(context).ncbiApiKey ?? '');
  bool _obscure = true;
  String? _error;

  static final _account = Uri.parse('https://account.ncbi.nlm.nih.gov/settings/');

  @override
  void dispose() {
    _key.dispose();
    super.dispose();
  }

  void _save() {
    if (!looksLikeNcbiApiKey(_key.text)) {
      setState(() => _error = context.l10n.apiKeyInvalid);
      return;
    }
    ProjectScope.read(context).setNcbiApiKey(_key.text.trim());
    Navigator.pop(context);
  }

  Future<void> _paste() async {
    final data = await Clipboard.getData(Clipboard.kTextPlain);
    final text = data?.text?.trim();
    if (text == null || text.isEmpty) return;
    setState(() {
      _key.text = text;
      _error = null;
    });
  }

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    final has = c.ncbiApiKey != null;
    final valid = looksLikeNcbiApiKey(_key.text);

    Widget step(int n, String text) => Padding(
          padding: const EdgeInsets.only(bottom: 10),
          child: Row(crossAxisAlignment: CrossAxisAlignment.start, children: [
            Container(
              width: 22,
              height: 22,
              alignment: Alignment.center,
              decoration: BoxDecoration(color: t.accentSoft, shape: BoxShape.circle, border: Border.all(color: t.accentSoftBorder)),
              child: Text('$n', style: TextStyle(fontSize: 12, fontWeight: FontWeight.w600, color: t.accentText)),
            ),
            const SizedBox(width: 12),
            Expanded(
              child: Padding(
                padding: const EdgeInsets.only(top: 2),
                child: Text(text, style: TextStyle(fontSize: 13, color: t.text, height: 1.45)),
              ),
            ),
          ]),
        );

    return SplaceDialog(
      icon: FA.key,
      title: l10n.apiKeyTitle,
      subtitle: l10n.apiKeyUsedFor,
      footerLeading: Text(l10n.apiKeyStored, style: TextStyle(fontSize: 12, color: t.textMuted)),
      actions: [
        if (has)
          FilledButton(
            style: DialogCancelButton.style(context),
            onPressed: () async {
              final ok = await confirmDialog(context, title: l10n.apiClearTitle, body: l10n.apiClearBody, confirm: l10n.apiKeyClear);
              if (ok && context.mounted) {
                c.setNcbiApiKey(null);
                Navigator.pop(context);
              }
            },
            child: Text(l10n.apiKeyClear),
          ),
        FilledButton(onPressed: () => Navigator.pop(context), child: Text(l10n.close)),
      ],
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        // Starts with the time comparison (the "Current" tag tells whether a
        // key is in use), so the steps and the NCBI link fit below.
        _SpeedComparison(hasKey: has),
        const SizedBox(height: 20),
        // Key field with the save button beside it, as in the HTML modal.
        Text(l10n.apiKeyFieldLabel, style: TextStyle(fontSize: 14, fontWeight: FontWeight.w600, color: t.text)),
        const SizedBox(height: 8),
        Row(crossAxisAlignment: CrossAxisAlignment.start, children: [
          Expanded(
            child: TextField(
              controller: _key,
              obscureText: _obscure,
              autofocus: !has,
              style: const TextStyle(fontFamily: monoFont, fontSize: 13, letterSpacing: 0.4),
              decoration: InputDecoration(
                hintText: '0123456789abcdef0123456789abcdef0123',
                errorText: _error,
                prefixIcon: Icon(FA.key, size: 13, color: valid ? t.accent : t.textFaint),
                prefixIconConstraints: const BoxConstraints(minWidth: 38),
                suffixIcon: Row(mainAxisSize: MainAxisSize.min, children: [
                  IconBtn(icon: FA.paste, tooltip: l10n.apiKeyPaste, size: ChipButton.height, onPressed: _paste),
                  IconBtn(icon: _obscure ? FA.eye : FA.eyeSlash, tooltip: '', size: ChipButton.height, onPressed: () => setState(() => _obscure = !_obscure)),
                  const SizedBox(width: 6),
                ]),
              ),
              onChanged: (_) => setState(() => _error = null),
              onSubmitted: (_) => _save(),
            ),
          ),
          const SizedBox(width: 10),
          FilledButton(
            onPressed: valid ? _save : null,
            child: Text(l10n.apiKeySave),
          ),
        ]),
        if (_key.text.trim().isNotEmpty && _error == null) ...[
          const SizedBox(height: 8),
          Row(children: [
            Icon(valid ? FA.circleCheck : FA.triangleExclamation, size: 12, color: valid ? SplaceTokens.success : SplaceTokens.warn),
            const SizedBox(width: 6),
            Text(valid ? l10n.apiKeyLooksValid : l10n.apiKeyInvalid, style: TextStyle(fontSize: 12, color: valid ? SplaceTokens.success : SplaceTokens.warn)),
          ]),
        ],
        const SizedBox(height: 16),
        Divider(color: t.border),
        const SizedBox(height: 10),
        // How to get one; the NCBI link sits on the title row so it is
        // always in view next to the steps.
        Row(children: [
          Expanded(child: Text(l10n.apiKeyHowTitle, style: TextStyle(fontSize: 14, fontWeight: FontWeight.w600, color: t.text))),
          OutlinedButton(
            onPressed: () => launchUrl(_account),
            child: Text(l10n.apiKeyOpenNcbi),
          ),
        ]),
        const SizedBox(height: 10),
        step(1, l10n.apiKeyHow1),
        step(2, l10n.apiKeyHow2),
        step(3, l10n.apiKeyHow3),
      ]),
    );
  }
}

/// Download time without and with a key for 50, 200 and 1000 genomes,
/// from [NcbiClient.estimateFetch]: one row per size with both times and a
/// bar whose blue part is the time with a key (the whole bar is the time
/// without one).
class _SpeedComparison extends StatelessWidget {
  const _SpeedComparison({required this.hasKey});
  final bool hasKey;

  static const _sizes = [50, 200, 1000];

  static String _time(Duration d) {
    final s = math.max(1, (d.inMilliseconds / 1000).round());
    if (s < 60) return '$s s';
    final m = s ~/ 60, r = s % 60;
    return r == 0 ? '$m min' : '$m min $r s';
  }

  @override
  Widget build(BuildContext context) {
    final l10n = context.l10n;
    final t = context.tokens;
    final dark = Theme.of(context).brightness == Brightness.dark;
    final lang = Localizations.localeOf(context).languageCode;
    final extraColor = dark ? t.accentSoftBorder : SplaceBlue.c100;
    final head = TextStyle(fontSize: 12, fontWeight: FontWeight.w700, color: t.textMuted);
    const numW = 112.0, gainW = 56.0, nameW = 110.0;

    Widget cell(double w, Widget child) => SizedBox(width: w, child: Align(alignment: Alignment.centerRight, child: child));

    Widget row(int n) {
      final without = NcbiClient.estimateFetch(n, withKey: false);
      final keyed = NcbiClient.estimateFetch(n, withKey: true);
      final share = keyed.inMilliseconds / without.inMilliseconds;
      final factor = NumberFormat('0.0', lang).format(without.inMilliseconds / keyed.inMilliseconds);
      return Padding(
        padding: const EdgeInsets.symmetric(vertical: 5),
        child: Row(children: [
          SizedBox(
            width: nameW,
            child: Text(l10n.apiKeyEstimateGenomes(NumberFormat.decimalPattern(lang).format(n)),
                style: TextStyle(fontSize: 13, fontWeight: FontWeight.w600, color: t.text)),
          ),
          Expanded(
            child: ClipRRect(
              borderRadius: BorderRadius.circular(6),
              child: SizedBox(
                height: 14,
                child: LayoutBuilder(
                  builder: (context, box) => Row(children: [
                    Container(width: box.maxWidth * share, color: t.accent),
                    Expanded(child: Container(color: extraColor)),
                  ]),
                ),
              ),
            ),
          ),
          cell(numW, Text('≈ ${_time(without)}', style: TextStyle(fontSize: 13, color: t.text))),
          cell(numW, Text('≈ ${_time(keyed)}', style: TextStyle(fontSize: 13, fontWeight: FontWeight.w700, color: t.accent))),
          cell(gainW, Text('$factor×', style: TextStyle(fontSize: 13, fontWeight: FontWeight.w600, color: t.accentText))),
        ]),
      );
    }

    Widget legend({required bool keyed, required Color color}) => Row(crossAxisAlignment: CrossAxisAlignment.start, children: [
          Container(
            margin: const EdgeInsets.only(top: 4),
            width: 10,
            height: 10,
            decoration: BoxDecoration(color: color, borderRadius: BorderRadius.circular(3), border: Border.all(color: t.accentSoftBorder)),
          ),
          const SizedBox(width: 8),
          Flexible(
            child: Text.rich(TextSpan(style: TextStyle(fontSize: 12, color: t.textMuted), children: [
              TextSpan(text: keyed ? l10n.apiKeyWith : l10n.apiKeyWithout, style: TextStyle(fontWeight: FontWeight.w700, color: t.text)),
              TextSpan(text: ' · ${l10n.apiKeyEstimateDetail('${NcbiClient.fetchBatch(withKey: keyed)}', '${NcbiClient.rateLimit(withKey: keyed)}')}'),
            ])),
          ),
          if (keyed == hasKey) ...[const SizedBox(width: 6), TagChip(l10n.apiKeyCurrent)],
        ]);

    final big = _sizes.last;
    final saved = NcbiClient.estimateFetch(big, withKey: false) - NcbiClient.estimateFetch(big, withKey: true);

    return Container(
      padding: const EdgeInsets.all(16),
      decoration: BoxDecoration(color: t.surface2, borderRadius: BorderRadius.circular(12), border: Border.all(color: t.border)),
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        Row(children: [
          Icon(FA.hourglassHalf, size: 12, color: t.textMuted),
          const SizedBox(width: 8),
          Flexible(child: Kicker(l10n.apiKeyEstimateTitle, size: 12)),
          const SizedBox(width: 6),
          InfoDot(l10n.apiKeyEstimateNote),
        ]),
        const SizedBox(height: 10),
        Row(children: [
          SizedBox(width: nameW, child: Text(l10n.apiKeyEstimateColSize, style: head)),
          const Spacer(),
          cell(numW, Text(l10n.apiKeyWithout, style: head)),
          cell(numW, Text(l10n.apiKeyWith, style: head)),
          cell(gainW, Text(l10n.apiKeyEstimateColGain, style: head)),
        ]),
        for (final n in _sizes) row(n),
        const SizedBox(height: 10),
        legend(keyed: true, color: t.accent),
        const SizedBox(height: 4),
        legend(keyed: false, color: extraColor),
        const SizedBox(height: 10),
        Row(children: [
          Icon(FA.bolt, size: 12, color: t.accent),
          const SizedBox(width: 8),
          Expanded(
            child: Text(l10n.apiKeyEstimateSaved(NumberFormat.decimalPattern(lang).format(big), _time(saved)),
                style: TextStyle(fontSize: 13, fontWeight: FontWeight.w600, color: t.accentText)),
          ),
        ]),
      ]),
    );
  }
}
