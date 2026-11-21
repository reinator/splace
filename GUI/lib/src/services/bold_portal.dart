import 'dart:async';
import 'dart:convert';
import 'dart:io';

import '../core/bold.dart';
import 'process_runner.dart' show CancelToken;
import 'resource_monitor.dart';

/// The public records of BOLD Systems through the v5 portal API
/// (portal.boldsystems.org/api): taxon names resolved by the query
/// preprocessor (any rank: species, genus, family…), the marker counts of
/// a query and its records (TSV), kept for the markers asked for.
class BoldPortal {
  BoldPortal({HttpClient? http}) : _http = http ?? (HttpClient()..connectionTimeout = const Duration(seconds: 30));

  static const base = 'https://portal.boldsystems.org/api';

  /// Above this many records the query is too broad to download.
  static const maxRecords = 60000;

  final HttpClient _http;

  /// The API takes queries of at most 250 characters.
  static const maxQuery = 240;

  /// Busy answers (429, 5xx) are retried after 2, 5, 10 and 20 s: the
  /// portal turns away requests that come in quick succession.
  Future<String> _get(String path, Map<String, String> query, {CancelToken? cancel}) async {
    final uri = Uri.parse('$base$path').replace(queryParameters: query);
    for (var attempt = 0;; attempt++) {
      if (cancel?.isCancelled ?? false) throw const BoldCancelled();
      final req = await _http.getUrl(uri).timeout(const Duration(seconds: 60));
      final res = await req.close().timeout(const Duration(minutes: 5));
      final body = await res.transform(utf8.decoder).join();
      NetMeter.add(received: body.length, sent: uri.toString().length);
      if (res.statusCode == 200) return body;
      final busy = res.statusCode == 429 || res.statusCode >= 500;
      if (busy && attempt < retryDelays.length) {
        await Future<void>.delayed(retryDelays[attempt]);
        continue;
      }
      throw HttpException('BOLD returned ${res.statusCode}: ${body.length > 200 ? body.substring(0, 200) : body}', uri: uri);
    }
  }

  /// Waits between the attempts of a busy request.
  List<Duration> retryDelays = const [Duration(seconds: 2), Duration(seconds: 5), Duration(seconds: 10), Duration(seconds: 20)];


  /// The BOLD taxon terms (`tax:species:Rhinella marina`,
  /// `tax:family:Bufonidae`) of [names]; names BOLD does not hold as a taxon
  /// (some synonyms) are left out.
  Future<List<String>> taxonTerms(List<String> names, {CancelToken? cancel}) async {
    final clean = [for (final n in names) n.replaceAll(';', ' ').trim()].where((n) => n.isNotEmpty && n.length <= maxQuery).toSet().toList();
    final terms = <String>{};
    for (final batch in queryBatches(clean)) {
      final body = await _get('/query/preprocessor', {'query': batch}, cancel: cancel);
      for (final t in ((jsonDecode(body) as Map<String, dynamic>)['successful_terms'] as List? ?? const [])) {
        for (final m in '${(t as Map)['matched'] ?? ''}'.split(';')) {
          if (m.startsWith('tax:')) terms.add(m);
        }
      }
    }
    return terms.toList();
  }

  /// Records per BOLD marker code of the query.
  Future<Map<String, int>> markerCounts(List<String> terms, {CancelToken? cancel}) async {
    final out = <String, int>{};
    for (final batch in queryBatches(terms)) {
      final body = await _get('/summary', {'query': batch, 'fields': 'marker_code', 'reduce_operation': 'count'}, cancel: cancel);
      final m = (jsonDecode(body) as Map<String, dynamic>)['marker_code'];
      if (m is Map) {
        for (final e in m.entries) {
          out.update('${e.key}', (v) => v + (e.value as num).toInt(), ifAbsent: () => (e.value as num).toInt());
        }
      }
    }
    return out;
  }

  /// The public records of [terms] for the marker codes [codes] (with a
  /// sequence).
  Future<List<BoldRecord>> records(List<String> terms, Set<String> codes, {CancelToken? cancel}) async {
    final out = <String, BoldRecord>{};
    for (final batch in queryBatches(terms)) {
      final q = jsonDecode(await _get('/query', {'query': batch, 'extent': 'full'}, cancel: cancel)) as Map<String, dynamic>;
      final id = q['query_id'] as String?;
      if (id == null) throw HttpException('BOLD did not accept the query $batch');
      final tsv = await _get('/documents/$id/download', {'format': 'tsv'}, cancel: cancel);
      // A specimen found by two terms of different batches comes once.
      for (final r in parseBoldTsv(tsv, codes)) {
        out['${r.processId}|${r.markerCode}'] = r;
      }
    }
    return out.values.toList();
  }
}

/// [items] joined by `;` in queries of at most [BoldPortal.maxQuery]
/// characters.
List<String> queryBatches(List<String> items) {
  final out = <String>[];
  var cur = '';
  for (final it in items) {
    if (cur.isNotEmpty && cur.length + 1 + it.length > BoldPortal.maxQuery) {
      out.add(cur);
      cur = '';
    }
    cur = cur.isEmpty ? it : '$cur;$it';
  }
  if (cur.isNotEmpty) out.add(cur);
  return out;
}

/// The rows of a portal TSV for the marker codes [codes] (all when empty),
/// with a sequence.
List<BoldRecord> parseBoldTsv(String tsv, Set<String> codes) {
  final lines = const LineSplitter().convert(tsv);
  if (lines.isEmpty) return const [];
  final head = lines.first.split('\t');
  final want = {for (final c in codes) c.toLowerCase()};
  final out = <BoldRecord>[];
  for (final line in lines.skip(1)) {
    if (line.trim().isEmpty) continue;
    final cells = line.split('\t');
    final row = {for (var i = 0; i < head.length && i < cells.length; i++) head[i]: cells[i]};
    final r = BoldRecord.fromRow(row);
    if (r.sequence.isEmpty || r.processId.isEmpty) continue;
    if (want.isNotEmpty && !want.contains(r.markerCode.toLowerCase())) continue;
    out.add(r);
  }
  return out;
}

class BoldCancelled implements Exception {
  const BoldCancelled();
  @override
  String toString() => 'Cancelled';
}
