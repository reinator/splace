import 'dart:async';
import 'dart:convert';
import 'dart:io';
import 'dart:isolate';
import 'dart:math' as math;

import 'package:path/path.dart' as p;

import '../core/genbank.dart';
import '../core/taxonomy_sources.dart';
import 'process_runner.dart' show CancelToken;
import 'resource_monitor.dart';

const ncbiExampleAccessions =
    'NC_024026.1, KJ184305.1, NC_024184.1, KJ556976.1, OP056886.1, PV818392.1, PX528615.1, NC_022707.1, KF356397.1, KJ642220.1, NC_057648.1, NC_082556.1, NC_023954.1';

/// NCBI accession numbers: 1-6 letters (optional underscore), 4+ digits,
/// optional version.
final accessionPattern = RegExp(r'\b[A-Za-z]{1,6}_?\d{4,}(?:\.\d+)?\b');

bool looksLikeNcbiApiKey(String? v) => v != null && RegExp(r'^[A-Za-z0-9_-]{20,}$').hasMatch(v.trim());

enum Organelle { mitochondrion, chloroplast }

class GenomeSearch {
  const GenomeSearch({required this.taxon, this.organelle = Organelle.mitochondrion, this.completeOnly = false, this.refseqOnly = false});
  final String taxon;
  final Organelle organelle;
  final bool completeOnly;
  final bool refseqOnly;

  /// Entrez query, identical to the HTML app's query preview.
  String get term {
    final scopes = completeOnly ? ['complete genome'] : ['complete genome', 'partial genome'];
    return [
      '"${taxon.trim()}"[Organism]',
      scopes.length == 1 ? '"${scopes.first}"[Title]' : '(${scopes.map((s) => '"$s"[Title]').join(' OR ')})',
      organelle == Organelle.chloroplast
          ? '(chloroplast[Title] OR plastid[Title] OR chloroplast[Filter])'
          : '(mitochondrion[Title] OR mitochondrial[Title] OR mitochondrion[Filter])',
      'biomol_genomic[PROP]',
      if (refseqOnly) 'srcdb_refseq[PROP]',
    ].join(' AND ');
  }
}

enum GenomeClass { complete, partial, mixed, unknown }

/// One NCBI hit shown in the "Review genomes before importing" dialog.
class GenomeCandidate {
  GenomeCandidate({required this.id, required this.accession, required this.title, required this.species, required this.description})
      : classification = _classify(title),
        unverified = RegExp(r'\bunverified\b', caseSensitive: false).hasMatch(title);

  final String id;
  final String accession;
  final String title;
  final String species;
  final String description;
  final GenomeClass classification;
  final bool unverified;
  bool selected = true;

  /// Sequence length (esummary `slen`).
  int? length;

  /// Marker searches: the markers whose query found this record.
  final Set<String> markers = {};

  static GenomeClass _classify(String title) {
    final t = title.toLowerCase();
    final complete = t.contains('complete');
    final partial = t.contains('partial') || t.contains('incomplete');
    if (complete && partial) return GenomeClass.mixed;
    if (complete) return GenomeClass.complete;
    if (partial) return GenomeClass.partial;
    return GenomeClass.unknown;
  }

  factory GenomeCandidate.fromSummary(String id, Map<String, dynamic> s) {
    final title = (s['title'] as String?) ?? '';
    final accession = (s['accessionversion'] as String?) ?? (s['caption'] as String?) ?? id;
    final species = _speciesFromTitle(title);
    return GenomeCandidate(id: id, accession: accession, title: title, species: species, description: _cleanDescription(title, accession, species))
      ..length = s['slen'] is int ? s['slen'] as int : int.tryParse('${s['slen'] ?? ''}');
  }

  bool matches(String query) {
    final q = _fold(query);
    return q.isEmpty || _fold('$accession $id $species $description $title').contains(q);
  }
}

String _fold(String s) => foldForSearch(s).trim();

/// Lower case without accents, one character per character, so match
/// positions in the folded text are positions in the original.
String foldForSearch(String s) => s
    .toLowerCase()
    .replaceAll(RegExp('[áàâãä]'), 'a')
    .replaceAll(RegExp('[éèêë]'), 'e')
    .replaceAll(RegExp('[íìîï]'), 'i')
    .replaceAll(RegExp('[óòôõö]'), 'o')
    .replaceAll(RegExp('[úùûü]'), 'u')
    .replaceAll('ç', 'c');

String _speciesFromTitle(String title) {
  final clean = title
      .replaceFirst(RegExp(r'\b(?:complete|partial|incomplete|unverified)\b.*$', caseSensitive: false), '')
      .replaceFirst(RegExp(r'\b(?:mitochondrion|mitochondrial|chloroplast|plastid)\b.*$', caseSensitive: false), '')
      .trim();
  final m = RegExp(r'^([A-Z][a-zA-Z.-]+\s+(?:cf\.\s+|aff\.\s+|sp\.\s+)?[a-z][a-zA-Z.-]+)').firstMatch(clean);
  return m?.group(1)?.replaceAll(RegExp(r'\s+'), ' ').trim() ?? '';
}

String _cleanDescription(String title, String accession, String species) {
  var d = title;
  for (final v in {accession, accession.replaceFirst(RegExp(r'\.\d+$'), '')}) {
    if (v.isNotEmpty) d = d.replaceAll(RegExp('\\b${RegExp.escape(v)}\\b'), '');
  }
  d = d
      .replaceAll(RegExp(r'\b(voucher|isolate|strain|specimen)\s+[^,;:]+', caseSensitive: false), '')
      .replaceAll(RegExp(r'\s{2,}'), ' ')
      .replaceAllMapped(RegExp(r'\s+([,;:])'), (m) => m.group(1)!)
      .replaceAll(RegExp(r'^[,;:\s]+|[,;:\s]+$'), '')
      .trim();
  if (species.isNotEmpty && !d.toLowerCase().startsWith(species.toLowerCase())) d = '$species $d'.trim();
  return d.isEmpty ? (species.isEmpty ? '-' : species) : d;
}

/// NCBI E-utilities (esearch, esummary, efetch) and the GBIF species API,
/// with an on-disk GenBank cache (7 days) and the NCBI request rate
/// (3 requests/s, 10 with an API key).
class NcbiClient {
  NcbiClient({required this.cacheDir, this.apiKey});

  static const _eutils = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils';
  static const cacheTtl = Duration(days: 7);

  final String cacheDir;
  String? apiKey;
  final _http = HttpClient()..userAgent = 'SPLACE/4 (https://github.com/luanrabelo/splace)';
  DateTime _last = DateTime.fromMillisecondsSinceEpoch(0);

  bool get hasKey => looksLikeNcbiApiKey(apiKey);
  Duration get _gap => requestGap(withKey: hasKey);

  /// NCBI's limit: 3 requests/s without a key, 10 with one.
  static int rateLimit({required bool withKey}) => withKey ? 10 : 3;

  /// Minimum time between two request starts, a little under [rateLimit]
  /// so network jitter does not push NCBI's count over it (at 350 ms NCBI
  /// still answered some requests with 429 in Nov 2026 tests).
  static Duration requestGap({required bool withKey}) => Duration(milliseconds: withKey ? 120 : 400);

  /// GenBank records per efetch request.
  static int fetchBatch({required bool withKey}) => withKey ? 10 : 5;

  /// efetch requests running at the same time. Starts are still spaced by
  /// [requestGap], so the rate never passes [rateLimit]; running several
  /// at once is what lets SPLACE reach it, since one response takes ~1 s.
  static int maxInFlight({required bool withKey}) => withKey ? 12 : 4;

  /// Average time between two efetch requests in a long download,
  /// occasional 429 retries included, and the time of the first response.
  /// Measured on Nov 8 2026 with Testudines mitogenomes (about 75 kB of
  /// GenBank text each, empty cache):
  ///
  /// | records | with a key (10/request) | without (5/request) |
  /// |---------|-------------------------|---------------------|
  /// | 50      | 1.6 s (5 requests)      | 5.5 s (10 + 2 retries) |
  /// | 200     | 4.4 s (20)              | 16.9 s (40)         |
  /// | 700     | 10.1 s (70 + 2 retries) | 64.6 s (140 + 20 retries) |
  static int _perRequestMs({required bool withKey}) => withKey ? 140 : 450;
  static const _firstResponseMs = 1200;

  /// Rough time to download [records] GenBank records that are not in the
  /// cache: requests are paced by the rate limit, so the time grows with
  /// the number of requests.
  static Duration estimateFetch(int records, {required bool withKey}) {
    final requests = (records / fetchBatch(withKey: withKey)).ceil();
    if (requests == 0) return Duration.zero;
    return Duration(milliseconds: requests * _perRequestMs(withKey: withKey) + _firstResponseMs);
  }

  /// Requests sent so far and HTTP codes that were retried (to check the
  /// rate).
  int requestCount = 0;
  final retries = <int>[];
  Future<void> _gate = Future.value();

  /// Waits for this request's turn: starts are queued and spaced by the
  /// gap, also when several requests run concurrently.
  Future<void> _throttle() {
    final turn = _gate.then((_) async {
      final wait = _last.add(_gap).difference(DateTime.now());
      if (!wait.isNegative) await Future<void>.delayed(wait);
      _last = DateTime.now();
      requestCount++;
    });
    _gate = turn;
    return turn;
  }

  Uri _uri(String tool, Map<String, String> params) =>
      Uri.parse('$_eutils/$tool.fcgi').replace(queryParameters: {...params, if (hasKey) 'api_key': apiKey!.trim(), 'tool': 'splace'});

  /// [post] sends the parameters in the body (long queries, such as the
  /// marker searches with every SynGenes synonym, are too long for a URL).
  Future<String> _get(Uri uri, {CancelToken? cancel, bool post = false}) async {
    for (var attempt = 0;; attempt++) {
      if (cancel?.isCancelled ?? false) throw const CancelledException();
      await _throttle();
      final HttpClientRequest req;
      if (post) {
        req = await _http.postUrl(uri.replace(query: ''));
        req.headers.contentType = ContentType('application', 'x-www-form-urlencoded', charset: 'utf-8');
        req.write(uri.query);
      } else {
        req = await _http.getUrl(uri);
      }
      final res = await req.close();
      final body = await res.transform(utf8.decoder).join();
      NetMeter.add(received: body.length, sent: uri.toString().length);
      if (res.statusCode == 200) return body;
      // 429: too many requests; retry with backoff.
      if ((res.statusCode == 429 || res.statusCode >= 500) && attempt < 4) {
        retries.add(res.statusCode);
        // Too many requests: hold the whole queue for a second, not only
        // this request, so the others do not hit the limit too.
        if (res.statusCode == 429) _last = DateTime.now().add(const Duration(seconds: 1));
        await Future<void>.delayed(Duration(milliseconds: 600 * (attempt + 1)));
        continue;
      }
      throw HttpException('HTTP ${res.statusCode}: ${body.length > 200 ? body.substring(0, 200) : body}', uri: uri);
    }
  }

  /// E-utilities JSON sometimes carries raw tabs/newlines inside strings
  /// (error messages), which strict JSON parsers reject.
  static Map<String, dynamic> _json(String body) =>
      jsonDecode(body.replaceAll(RegExp(r'[\x00-\x1F]'), ' ')) as Map<String, dynamic>;

  /// Ids of the genomes matching [s], at most [retmax] (esearch returns up
  /// to 10,000 per request; larger searches are paged).
  Future<List<String>> searchGenomes(GenomeSearch s, {int retmax = 100000}) => searchTerm(s.term, retmax: retmax);

  /// Ids of the nucleotide records matching an Entrez [term].
  Future<List<String>> searchTerm(String term, {int retmax = 100000, CancelToken? cancel}) async {
    final ids = <String>[];
    while (ids.length < retmax) {
      final page = await _searchPage(term, retstart: ids.length, retmax: math.min(10000, retmax - ids.length), cancel: cancel);
      ids.addAll(page.ids);
      if (page.ids.isEmpty || ids.length >= page.count) break;
    }
    return ids;
  }

  Future<({List<String> ids, int count})> _searchPage(String term, {required int retstart, required int retmax, CancelToken? cancel}) async {
    for (var attempt = 0;; attempt++) {
      final body = await _get(_uri('esearch', {'db': 'nucleotide', 'retmode': 'json', 'retstart': '$retstart', 'retmax': '$retmax', 'term': term}),
          cancel: cancel, post: true);
      final result = _json(body)['esearchresult'] as Map<String, dynamic>?;
      final error = result?['ERROR'];
      if (error != null) {
        // "Search Backend failed" is transient on NCBI's side.
        if (attempt < 3) {
          await Future<void>.delayed(Duration(seconds: 1 + attempt));
          continue;
        }
        throw HttpException('NCBI search failed: $error');
      }
      final ids = ((result?['idlist'] as List?) ?? const []).cast<String>();
      return (ids: ids, count: int.tryParse('${result?['count'] ?? ''}') ?? ids.length);
    }
  }

  Future<List<GenomeCandidate>> summaries(List<String> ids, {void Function(int done, int total)? onProgress}) async {
    final out = <GenomeCandidate>[];
    final chunks = (ids.length / 200).ceil();
    for (var i = 0; i < ids.length; i += 200) {
      onProgress?.call(i ~/ 200, chunks);
      final chunk = ids.sublist(i, (i + 200).clamp(0, ids.length));
      final body = await _get(_uri('esummary', {'db': 'nucleotide', 'retmode': 'json', 'id': chunk.join(',')}));
      final result = (_json(body)['result'] as Map<String, dynamic>?) ?? const {};
      for (final id in chunk) {
        final s = result[id];
        if (s is Map<String, dynamic>) out.add(GenomeCandidate.fromSummary(id, s));
      }
    }
    return out;
  }

  File _cacheFile(String id) => File(p.join(cacheDir, '${id.replaceAll(RegExp(r'[^A-Za-z0-9_.-]'), '_')}.gb'));

  String? _cached(String id) {
    final f = _cacheFile(id);
    if (!f.existsSync()) return null;
    if (DateTime.now().difference(f.lastModifiedSync()) > cacheTtl) return null;
    return f.readAsStringSync();
  }

  /// Downloads GenBank records (accessions or UIDs) in batches, reporting
  /// progress per record. Failed ids are returned in `failed`.
  Future<({List<GbRecord> records, List<String> failed})> fetchGenBank(
    List<String> ids, {
    CancelToken? cancel,
    void Function(int done, int total, String current)? onProgress,
    void Function(GbRecord record)? onRecord,
    void Function(String message)? log,
  }) async {
    final records = <GbRecord>[];
    final failed = <String>[];
    var done = 0;
    final pending = <String>[];
    // Cached records are parsed together in one isolate.
    final cachedTexts = <String>[];
    for (final id in ids) {
      final cached = _cached(id);
      if (cached != null) {
        cachedTexts.add(cached);
      } else {
        pending.add(id);
      }
    }
    if (cachedTexts.isNotEmpty) {
      final parsed = await _parseInIsolate(cachedTexts.join('\n'));
      records.addAll(parsed);
      parsed.forEach(onRecord ?? (_) {});
      done += cachedTexts.length;
      onProgress?.call(done, ids.length, parsed.isEmpty ? '' : parsed.last.accession);
      log?.call('${cachedTexts.length} record(s) loaded from the local cache');
    }
    final batch = fetchBatch(withKey: hasKey);
    Directory(cacheDir).createSync(recursive: true);
    final chunks = [for (var i = 0; i < pending.length; i += batch) pending.sublist(i, (i + batch).clamp(0, pending.length))];

    Future<void> fetchChunk(List<String> chunk) async {
      try {
        final text = await _get(_uri('efetch', {'db': 'nucleotide', 'id': chunk.join(','), 'rettype': 'gb', 'retmode': 'text'}), cancel: cancel);
        final parsed = await _parseInIsolate(text);
        // Cache per record, under both the requested id and the accession.
        final raws = text.split(RegExp(r'^//\s*$', multiLine: true));
        for (var k = 0; k < parsed.length; k++) {
          final raw = k < raws.length ? '${raws[k].trim()}\n//\n' : null;
          if (raw == null) continue;
          _cacheFile(parsed[k].accession).writeAsStringSync(raw);
          if (k < chunk.length && parsed.length == chunk.length) _cacheFile(chunk[k]).writeAsStringSync(raw);
        }
        records.addAll(parsed);
        if (parsed.length < chunk.length) {
          final got = parsed.map((r) => r.accession.split('.').first).toSet();
          failed.addAll(chunk.where((c) => !got.contains(c.split('.').first) && !RegExp(r'^\d+$').hasMatch(c)));
        }
        for (final r in parsed) {
          onRecord?.call(r);
          log?.call('Fetched ${r.organism} (${r.accession})');
        }
      } on CancelledException {
        return;
      } catch (e) {
        log?.call('Failed to fetch ${chunk.join(', ')}: $e');
        failed.addAll(chunk);
      }
      done += chunk.length;
      onProgress?.call(done.clamp(0, ids.length), ids.length, chunk.last);
    }

    // A small pool of workers takes the batches in order; the throttle
    // spaces their requests.
    var next = 0;
    Future<void> worker() async {
      while (next < chunks.length && !(cancel?.isCancelled ?? false)) {
        await fetchChunk(chunks[next++]);
      }
    }

    await Future.wait([for (var w = 0; w < math.min(maxInFlight(withKey: hasKey), chunks.length); w++) worker()]);
    return (records: records, failed: failed);
  }

  /// NCBI Taxonomy id of a scientific name (records without
  /// `/db_xref="taxon:…"`).
  Future<String?> ncbiTaxId(String species) async {
    final body = await _get(_uri('esearch', {'db': 'taxonomy', 'retmode': 'json', 'term': '"${speciesKey(species)}"[Scientific Name]'}));
    final ids = (_json(body)['esearchresult']?['idlist'] as List?) ?? const [];
    return ids.isEmpty ? null : ids.first as String;
  }

  /// Ranks of NCBI Taxonomy ids (batches of 100).
  Future<Map<String, TaxonInfo>> ncbiTaxonomy(List<String> taxIds) async {
    final out = <String, TaxonInfo>{};
    for (var i = 0; i < taxIds.length; i += 100) {
      final chunk = taxIds.sublist(i, (i + 100).clamp(0, taxIds.length));
      final xml = await _get(_uri('efetch', {'db': 'taxonomy', 'id': chunk.join(','), 'retmode': 'xml'}));
      out.addAll(parseNcbiTaxonomyXml(xml));
    }
    return out;
  }

  Future<TaxonInfo?> gbif(String species) async {
    final uri = Uri.parse('https://api.gbif.org/v1/species').replace(queryParameters: {'name': species});
    final req = await _http.getUrl(uri);
    req.headers.set(HttpHeaders.acceptHeader, 'application/json');
    final res = await req.close();
    final body = await res.transform(utf8.decoder).join();
    NetMeter.add(received: body.length, sent: uri.toString().length);
    if (res.statusCode != 200) throw HttpException('GBIF returned ${res.statusCode}', uri: uri);
    final results = (jsonDecode(body)['results'] as List?)?.cast<Map<String, dynamic>>() ?? const [];
    if (results.isEmpty) return null;
    final a = results.firstWhere((r) => r['taxonomicStatus'] == 'ACCEPTED', orElse: () => results.first);
    String s(String k) => (a[k] as String?) ?? '';
    return TaxonInfo(
        kingdom: s('kingdom'), phylum: s('phylum'), klass: s('class'), order: s('order'), family: s('family'), authorship: s('authorship'), source: 'GBIF');
  }

  static const _salveApi = 'https://salve.icmbio.gov.br/salve-api/public';

  Future<dynamic> _salveJson(String endpoint, Map<String, String> params, {CancelToken? cancel}) async {
    final uri = Uri.parse('$_salveApi/$endpoint').replace(queryParameters: params);
    for (var attempt = 0;; attempt++) {
      if (cancel?.isCancelled ?? false) throw const CancelledException();
      final req = await _http.getUrl(uri);
      req.headers.set(HttpHeaders.acceptHeader, 'application/json');
      final res = await req.close();
      final body = await res.transform(utf8.decoder).join();
      NetMeter.add(received: body.length, sent: uri.toString().length);
      if (res.statusCode == 200) return _json(body)['data'];
      if ((res.statusCode == 429 || res.statusCode >= 500) && attempt < 3) {
        await Future<void>.delayed(Duration(milliseconds: 800 * (attempt + 1)));
        continue;
      }
      throw HttpException('SALVE returned ${res.statusCode}', uri: uri);
    }
  }

  /// ICMBio SALVE public API (salve.icmbio.gov.br/salve-api/public): the
  /// species sheet from `search`, then every section of the sheet from
  /// `fichaHtml` (see [salveSections]), together. Null when the species has no sheet: SALVE only
  /// assesses the Brazilian fauna.
  Future<TaxonInfo?> salve(String species, {CancelToken? cancel}) async {
    final key = speciesKey(species);
    if (key.length < 3) return null;
    final hits = await _salveJson('search', {'q': key, 'paginationPageSize': '50'}, cancel: cancel);
    if (hits is! List || hits.isEmpty) return null;
    final row = salvePickSheet(hits.whereType<Map<String, dynamic>>().toList(), key);
    if (row == null) return null;
    final id = (row['id_ficha_publicada'] ?? row['id_ficha']) as String?;
    Future<Map<String, dynamic>?> section(String name) async {
      if (id == null) return null;
      try {
        final d = await _salveJson('fichaHtml', {'idFicha': id, 'section': name}, cancel: cancel);
        return d is Map<String, dynamic> ? d : null;
      } on CancelledException {
        rethrow;
      } catch (_) {
        // The sheet is still useful without one of its sections.
        return null;
      }
    }

    // Every section of the sheet, together.
    final parts = await Future.wait(salveSections.map(section));
    return salveTaxonInfo(species, row, sections: {for (final (i, k) in salveSections.indexed) k: parts[i]});
  }

  void close() => _http.close(force: true);
}



/// Top level so the isolate closure never captures the client.
Future<List<GbRecord>> _parseInIsolate(String text) => Isolate.run(() => parseGenBankFile(text, source: 'ncbi'));

class CancelledException implements Exception {
  const CancelledException();
  @override
  String toString() => 'cancelled';
}
