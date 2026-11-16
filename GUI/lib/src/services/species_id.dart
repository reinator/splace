/// Species identification of loaded genomes against GenBank (NCBI BLAST)
/// and BOLD (ID engine), so wrong names can be fixed before the taxonomy.
///
/// - GenBank: the NCBI BLAST URL API (`Blast.cgi`, CMD=Put / Get), blastn
///   megablast on core_nt, following NCBI's usage rules: one submission at
///   a time and at most one status check per minute for each RID.
/// - BOLD: the v5 Identification Engine at id.boldsystems.org (the v4
///   Ids_xml service was retired in August 2025). It is the API the BOLD
///   web form uses (`/submission`, `/submission/status/{id}`,
///   `/submission/results/{id}`); BOLD does not document it, so it may
///   change. Only COI is sent to BOLD.
library;

import 'dart:async';
import 'dart:convert';
import 'dart:io';
import 'dart:math' as math;

import '../core/fasta.dart' show SeqRecord, fingerprint;
import '../core/genbank.dart';
import 'process_runner.dart' show CancelToken;
import 'ncbi.dart' show CancelledException;
import 'resource_monitor.dart';

enum IdSource { genbank, bold }

/// One query sent to BLAST/BOLD: a gene of a record.
class IdQuery {
  const IdQuery(this.record, this.gene, this.sequence);
  final String record; // accession
  final String gene;
  final String sequence;

  String get id => '$record|$gene';
}

/// Best evidence for one species from one source and gene.
class IdEvidence {
  IdEvidence(this.identity, this.hits);

  /// Best percent identity (0–100) among other GenBank/BOLD records.
  double identity;

  /// Hits of this species at that gene (other records).
  int hits;

  /// Identity of the record itself in the results (BLAST finds the query's
  /// own GenBank record, or the INSDC record a RefSeq copy came from). It
  /// shows the data is in GenBank but does not confirm the name.
  double? selfIdentity;

  /// Accessions of those self hits.
  final Set<String> selfAccessions = {};

  /// Only the record itself was found for this species.
  bool get onlySelf => hits == 0 && selfIdentity != null;
}

/// A species proposed for a record, with the evidence of each source.
class SpeciesCandidate {
  SpeciesCandidate(this.species);
  final String species;

  /// gene → evidence from GenBank BLAST.
  final Map<String, IdEvidence> genbank = {};

  /// COI evidence from BOLD.
  IdEvidence? bold;

  /// Genes with evidence from other GenBank records (not only the record
  /// itself).
  Iterable<IdEvidence> get _independent => genbank.values.where((e) => e.hits > 0);

  bool get inGenbank => _independent.isNotEmpty;

  double get genbankIdentity => _independent.fold(0.0, (a, e) => math.max(a, e.identity));

  double get bestIdentity => math.max(genbankIdentity, bold?.identity ?? 0);

  /// Sources (GenBank, BOLD) that found this species at a same-species
  /// identity: a 93 % GenBank hit does not back a 99.8 % BOLD one.
  int get sourceCount =>
      (genbankIdentity >= sameSpeciesIdentity ? 1 : 0) + ((bold?.identity ?? 0) >= sameSpeciesIdentity ? 1 : 0);

  int get totalHits => genbank.values.fold(0, (a, e) => a + e.hits) + (bold?.hits ?? 0);
}

class RecordIdentification {
  RecordIdentification(this.accession, this.currentName);
  final String accession;
  final String currentName;
  final Map<String, SpeciesCandidate> _candidates = {};

  /// Sources that answered for this record.
  final Set<IdSource> sources = {};

  SpeciesCandidate _candidate(String species) => _candidates.putIfAbsent(species, () => SpeciesCandidate(species));

  /// Counts one hit of [species] at [identity] % from [source] ([gene] for
  /// GenBank; BOLD is always COI). A hit on the record itself ([self], its
  /// accession [acc]) is kept apart: shown, but not counted as evidence.
  void addEvidence(IdSource source, String gene, String species, double identity, {bool self = false, String acc = ''}) {
    sources.add(source);
    final c = _candidate(species);
    final ev = source == IdSource.genbank ? c.genbank.putIfAbsent(gene, () => IdEvidence(0, 0)) : (c.bold ??= IdEvidence(0, 0));
    if (self) {
      ev.selfIdentity = math.max(ev.selfIdentity ?? 0, identity);
      if (acc.isNotEmpty) ev.selfAccessions.add(acc);
      return;
    }
    ev
      ..identity = math.max(ev.identity, identity)
      ..hits += 1;
  }

  /// GenBank holds this record's sequence (BLAST found it) but no other
  /// record of the same species.
  bool get onlySelfInGenbank {
    final c = current;
    return c != null && !c.inGenbank && c.genbank.values.any((e) => e.onlySelf);
  }

  /// Candidates, best first: identity; within 0.5 % the species backed by
  /// both GenBank and BOLD comes first, then the one with more hits.
  List<SpeciesCandidate> get candidates => _candidates.values.toList()
    ..sort((a, b) {
      if ((a.bestIdentity - b.bestIdentity).abs() >= 0.5) return b.bestIdentity.compareTo(a.bestIdentity);
      final bySources = b.sourceCount.compareTo(a.sourceCount);
      if (bySources != 0) return bySources;
      final byId = b.bestIdentity.compareTo(a.bestIdentity);
      return byId != 0 ? byId : b.totalHits.compareTo(a.totalHits);
    });

  /// What the results table shows: the best species of each source (BOLD
  /// alone can return many near-identical names that would push the
  /// GenBank matches out of view) and the current name when found.
  List<SpeciesCandidate> shortlist({int perSource = 4}) {
    final all = candidates;
    final gb = all.where((c) => c.inGenbank).toList()..sort((a, b) => b.genbankIdentity.compareTo(a.genbankIdentity));
    final bold = all.where((c) => c.bold != null).toList()..sort((a, b) => b.bold!.identity.compareTo(a.bold!.identity));
    final pick = <SpeciesCandidate>{...gb.take(perSource), ...bold.take(perSource), ?current};
    return all.where(pick.contains).toList();
  }

  /// The candidate with the record's current name, if found.
  SpeciesCandidate? get current => candidates.where((c) => sameSpecies(c.species, currentName)).firstOrNull;

  /// The best candidate names another species than the record, at an
  /// identity that means the same species. Below that the closest match is
  /// only a relative (expected when no other record of the species exists),
  /// which does not contest the name.
  bool get disagrees {
    final c = candidates;
    if (c.isEmpty || _norm(c.first.species) == _norm(currentName) || c.first.bestIdentity < sameSpeciesIdentity) return false;
    // The current name matching as well (identical haplotypes of close
    // species) is ambiguous, not a misidentification.
    final cur = current;
    return cur == null || cur.bestIdentity < c.first.bestIdentity;
  }

  /// Accessions of the record found by BLAST (itself or its INSDC source).
  List<String> get selfAccessions => {
        for (final c in _candidates.values)
          for (final e in c.genbank.values) ...e.selfAccessions,
      }.toList()
        ..sort();
}

/// Identity (%) from which a match is taken as the same species.
const sameSpeciesIdentity = 98.0;

/// Same species name, ignoring case, punctuation and spacing.
bool sameSpecies(String a, String b) => _norm(a) == _norm(b);

String _norm(String s) => s.toLowerCase().replaceAll(RegExp(r'[^a-z ]'), '').replaceAll(RegExp(r'\s+'), ' ').trim();

/// Genus + species of a name ("Caranx hippos voucher X" → "Caranx hippos").
String binomialOf(String name) {
  final parts = name.trim().split(RegExp(r'\s+'));
  if (parts.length < 2) return name.trim();
  if (parts[1] == 'sp.' || parts[1] == 'cf.' || parts[1] == 'aff.') return parts.take(math.min(3, parts.length)).join(' ');
  return '${parts[0]} ${parts[1]}';
}

/// Minimum share of the query that a hit must cover to count.
const minCoverage = 0.8;

class SpeciesIdentifier {
  SpeciesIdentifier({this.onStatus, this.onProgress});

  /// Progress messages for the modal.
  final void Function(String message)? onStatus;

  /// Share of each source's work done (0–1). BLAST advances by stage
  /// (submitted, waiting against NCBI's own time estimate, results); BOLD
  /// by sequences completed.
  final void Function(IdSource source, double fraction)? onProgress;
  final _http = HttpClient()..userAgent = 'SPLACE/4 (https://github.com/luanrabelo/splace)';

  static const _blast = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi';
  static const _bold = 'https://id.boldsystems.org';

  /// Queries per BLAST submission and per BOLD submission.
  static const blastChunk = 40;
  static const boldChunk = 100;

  void close() => _http.close(force: true);

  /// One request; network errors and 5xx/429 answers are retried (3
  /// tries, growing pauses).
  Future<String> _send(String method, Uri uri, {List<int>? body, String? contentType}) async {
    for (var attempt = 1;; attempt++) {
      try {
        final req = await _http.openUrl(method, uri);
        if (contentType != null) req.headers.set(HttpHeaders.contentTypeHeader, contentType);
        if (body != null) req.add(body);
        final res = await req.close();
        final text = await res.transform(const Utf8Decoder(allowMalformed: true)).join();
        NetMeter.add(received: text.length, sent: uri.toString().length + (body?.length ?? 0));
        if (res.statusCode == 200) return text;
        final error = HttpException('HTTP ${res.statusCode}: ${text.length > 200 ? text.substring(0, 200) : text}', uri: uri);
        if ((res.statusCode < 500 && res.statusCode != 429) || attempt >= 3) throw error;
      } on SocketException {
        if (attempt >= 3) rethrow;
      } on HttpException catch (e) {
        if (attempt >= 3 || !(e.message.contains('Connection') || e.message.startsWith('HTTP 5') || e.message.startsWith('HTTP 429'))) rethrow;
      }
      await Future<void>.delayed(Duration(seconds: 3 * attempt));
    }
  }

  /// Why each source failed in the last [identify], if it did.
  final Map<IdSource, String> failures = {};

  static Future<void> _wait(Duration d, CancelToken? cancel, [void Function()? tick]) async {
    final end = DateTime.now().add(d);
    var last = DateTime.now();
    while (DateTime.now().isBefore(end)) {
      if (cancel?.isCancelled ?? false) throw const CancelledException();
      await Future<void>.delayed(const Duration(milliseconds: 250));
      if (tick != null && DateTime.now().difference(last).inMilliseconds >= 1000) {
        last = DateTime.now();
        tick();
      }
    }
  }

  static String _fasta(List<IdQuery> qs) => qs.map((q) => '>${q.id}\n${q.sequence}').join('\n');

  // -------------------------------------------------------------------------
  // GenBank (BLAST)
  // -------------------------------------------------------------------------

  /// BLAST hits per query id: (species, identity %, accession).
  Future<Map<String, List<(String, double, String)>>> blast(List<IdQuery> queries, {CancelToken? cancel}) async {
    final out = <String, List<(String, double, String)>>{};
    for (var i = 0; i < queries.length; i += blastChunk) {
      final chunk = queries.sublist(i, math.min(i + blastChunk, queries.length));
      final part = i ~/ blastChunk + 1, parts = (queries.length / blastChunk).ceil();
      void report(double within) => onProgress?.call(IdSource.genbank, (part - 1 + within) / parts);
      report(0.02);
      final body = Uri(queryParameters: {
        'CMD': 'Put',
        'PROGRAM': 'blastn',
        'MEGABLAST': 'on',
        'DATABASE': 'core_nt',
        'HITLIST_SIZE': '25',
        'TOOL': 'splace',
        'QUERY': _fasta(chunk),
      }).query;
      final put = await _send('POST', Uri.parse(_blast), body: utf8.encode(body), contentType: 'application/x-www-form-urlencoded');
      final rid = RegExp(r'RID = (\w+)').firstMatch(put)?.group(1);
      final rtoe = int.tryParse(RegExp(r'RTOE = (\d+)').firstMatch(put)?.group(1) ?? '') ?? 30;
      if (rid == null) throw const HttpException('BLAST did not return a request id (RID)');
      onStatus?.call('GenBank BLAST $part/$parts · RID $rid · ~${rtoe}s');
      // While NCBI works the bar follows its own estimate (RTOE), up to 90 %.
      final started = DateTime.now();
      final expected = math.max(rtoe, 20) * 1.5;
      void tick() => report(0.05 + 0.85 * math.min(1, DateTime.now().difference(started).inSeconds / expected));
      report(0.05);
      // First check after the estimate, then at most once a minute.
      await _wait(Duration(seconds: math.max(rtoe, 20)), cancel, tick);
      for (var polls = 0;; polls++) {
        final info = await _send('GET', Uri.parse('$_blast?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=$rid'));
        final status = RegExp(r'Status=(\w+)').firstMatch(info)?.group(1) ?? 'UNKNOWN';
        if (status == 'READY') break;
        if (status != 'WAITING' || polls > 30) throw HttpException('BLAST $rid: $status');
        onStatus?.call('GenBank BLAST $part/$parts · RID $rid · waiting for NCBI…');
        await _wait(const Duration(seconds: 60), cancel, tick);
      }
      final json = jsonDecode(await _send('GET', Uri.parse('$_blast?CMD=Get&RID=$rid&FORMAT_TYPE=JSON2_S'))) as Map<String, dynamic>;
      for (final o in (json['BlastOutput2'] as List? ?? const [])) {
        final search = (o as Map)['report']?['results']?['search'] as Map?;
        if (search == null) continue;
        final id = search['query_title'] as String? ?? '';
        final qlen = (search['query_len'] as num?)?.toInt() ?? 0;
        final hits = out.putIfAbsent(id, () => []);
        for (final h in (search['hits'] as List? ?? const [])) {
          // core_nt merges identical sequences into one hit with several
          // descriptions (e.g. a RefSeq record and the GenBank record it
          // was copied from, or two specimens with the same haplotype).
          final descs = ((h as Map)['description'] as List?)?.cast<Map>() ?? const <Map>[];
          final hsp = (h['hsps'] as List?)?.cast<Map>().firstOrNull;
          if (descs.isEmpty || hsp == null) continue;
          final len = (hsp['align_len'] as num).toInt();
          if (qlen > 0 && len / qlen < minCoverage) continue;
          final identity = 100 * (hsp['identity'] as num) / len;
          for (final desc in descs) {
            final name = (desc['sciname'] as String?) ?? '';
            if (name.isEmpty) continue;
            hits.add((binomialOf(name), identity, (desc['accession'] as String?) ?? ''));
          }
        }
      }
      report(1);
      // NCBI asks for no more than one submission every 10 seconds.
      if (i + blastChunk < queries.length) await _wait(const Duration(seconds: 10), cancel);
    }
    return out;
  }

  // -------------------------------------------------------------------------
  // BOLD (ID engine v5)
  // -------------------------------------------------------------------------

  /// BOLD hits per query id: (species, identity %, process id).
  Future<Map<String, List<(String, double, String)>>> bold(List<IdQuery> queries, {CancelToken? cancel}) async {
    final out = <String, List<(String, double, String)>>{};
    for (var i = 0; i < queries.length; i += boldChunk) {
      final chunk = queries.sublist(i, math.min(i + boldChunk, queries.length));
      final part = i ~/ boldChunk, parts = (queries.length / boldChunk).ceil();
      void report(double within) => onProgress?.call(IdSource.bold, (part + within) / parts);
      report(0.02);
      const boundary = 'splace-bold-boundary';
      final body = utf8.encode('--$boundary\r\n'
          'Content-Disposition: form-data; name="fasta_file"; filename="submitted.fas"\r\n'
          'Content-Type: text/plain\r\n\r\n'
          '${_fasta(chunk)}\r\n'
          '--$boundary--\r\n');
      final sub = jsonDecode(await _send(
        'POST',
        Uri.parse('$_bold/submission?db=public.tax-derep&mi=0.9&mo=100&maxh=25&order=3'),
        body: body,
        contentType: 'multipart/form-data; boundary=$boundary',
      )) as Map<String, dynamic>;
      final id = sub['sub_id'] as String?;
      if (id == null) throw const HttpException('BOLD did not return a submission id');
      final path = Uri.encodeComponent(id);
      for (var polls = 0;; polls++) {
        await _wait(Duration(seconds: polls == 0 ? 5 : 10), cancel);
        final st = jsonDecode(await _send('GET', Uri.parse('$_bold/submission/status/$path'))) as Map<String, dynamic>;
        if (st['status'] == 'complete') break;
        final done = st['completed'] ?? 0, total = (st['queued'] ?? 0) + (st['processing'] ?? 0) + done;
        onStatus?.call('BOLD · $done/$total ${total == 0 ? '' : 'sequences'}');
        if (total > 0) report(0.05 + 0.9 * done / total);
        if (polls > 90) throw const HttpException('BOLD did not finish in 15 minutes');
      }
      final jsonl = await _send('GET', Uri.parse('$_bold/submission/results/$path'));
      report(1);
      for (final line in const LineSplitter().convert(jsonl)) {
        if (line.trim().isEmpty) continue;
        final r = jsonDecode(line) as Map<String, dynamic>;
        final hits = out.putIfAbsent(r['seqid'] as String? ?? '', () => []);
        final results = (r['results'] as Map?) ?? const {};
        for (final e in results.entries) {
          final v = e.value as Map;
          final species = ((v['taxonomy'] as Map?)?['species'] as String?) ?? '';
          if (species.isEmpty) continue;
          hits.add((binomialOf(species), (v['pident'] as num).toDouble(), (e.key as String).split('|').first));
        }
      }
    }
    return out;
  }

  /// Queries answered from the cache in the last [identify].
  int reused = 0;

  /// Queries sent to BLAST/BOLD in the last [identify].
  int sent = 0;

  /// Cache key of one query at one source: record, gene and sequence.
  /// GenBank keys carry a version: answers saved before every description
  /// of merged BLAST hits was read (`genbank|…`) are not reused.
  static String cacheKey(IdSource source, String record, String gene, String sequence) =>
      '${source == IdSource.genbank ? 'genbank2' : source.name}|$record|$gene|${fingerprint([SeqRecord('', sequence)])}';

  /// Cache entries no longer produced by [cacheKey].
  static bool staleCacheKey(String key) => key.startsWith('genbank|');

  /// Identifies [records] with the [genes] given (name → sequence per
  /// record): every gene goes to GenBank, COI also to BOLD. Queries found in
  /// [cache] (see [cacheKey]) are not sent again; new answers are added.
  Future<Map<String, RecordIdentification>> identify(
    List<GbRecord> records,
    Map<String, Map<String, String>> sequences, {
    CancelToken? cancel,
    void Function(String message)? log,
    Map<String, List<(String, double, String)>>? cache,
  }) async {
    failures.clear();
    final store = cache ?? {};
    final ids = {for (final r in records) r.accession: RecordIdentification(r.accession, r.organismName)};
    String base(String acc) => acc.split('.').first;
    // The record and the INSDC records a RefSeq copy came from.
    final self = {
      for (final r in records) r.accession: {base(r.accession), for (final a in r.sourceAccessions) base(a)},
    };
    final queries = [
      for (final r in records)
        for (final e in (sequences[r.accession] ?? const <String, String>{}).entries)
          if (e.value.length >= 100) IdQuery(r.accession, e.key, e.value),
    ];

    void add(IdSource source, List<IdQuery> qs, Map<String, List<(String, double, String)>> hits) {
      for (final q in qs) {
        final list = hits[q.id];
        if (list == null) continue;
        final rid = ids[q.record]!..sources.add(source);
        for (final (species, identity, acc) in list) {
          // The record itself is shown but is not evidence.
          final isSelf = source == IdSource.genbank && (self[q.record]?.contains(base(acc)) ?? false);
          rid.addEvidence(source, q.gene, species, identity, self: isSelf, acc: acc);
        }
      }
    }

    String key(IdSource s, IdQuery q) => cacheKey(s, q.record, q.gene, q.sequence);

    /// Answers [qs] from the cache and sends the rest with [run].
    Future<void> source(IdSource s, List<IdQuery> qs, Future<Map<String, List<(String, double, String)>>> Function(List<IdQuery>) run) async {
      final known = [for (final q in qs) if (store.containsKey(key(s, q))) q];
      final todo = [for (final q in qs) if (!store.containsKey(key(s, q))) q];
      reused += known.length;
      sent += todo.length;
      add(s, known, {for (final q in known) q.id: store[key(s, q)]!});
      if (todo.isEmpty) {
        onProgress?.call(s, 1);
        return;
      }
      final hits = await run(todo);
      for (final q in todo) {
        store[key(s, q)] = hits[q.id] ?? const [];
      }
      add(s, todo, hits);
    }

    reused = 0;
    sent = 0;
    // Both sources run at the same time; a failure of one keeps the other.
    final coi = queries.where((q) => q.gene == 'COI').toList();
    await Future.wait([
      source(IdSource.genbank, queries, (qs) => blast(qs, cancel: cancel)).catchError((Object e) {
        if (e is CancelledException) throw e;
        failures[IdSource.genbank] = '$e';
        log?.call('GenBank BLAST failed: $e');
      }),
      if (coi.isNotEmpty)
        source(IdSource.bold, coi, (qs) => bold(qs, cancel: cancel)).catchError((Object e) {
          if (e is CancelledException) throw e;
          failures[IdSource.bold] = '$e';
          log?.call('BOLD failed: $e');
        }),
    ]);
    return ids;
  }
}
