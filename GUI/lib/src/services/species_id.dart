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

/// Parameters of the identification: NCBI BLAST (URL API) and the BOLD v5
/// ID engine (id.boldsystems.org).
class SpeciesIdOptions {
  const SpeciesIdOptions({
    this.blastDatabase = 'core_nt',
    this.blastProgram = 'megablast',
    this.expect = 0.05,
    this.wordSize,
    this.hitlistSize = 25,
    this.entrezQuery = '',
    this.minCoverage = 0.8,
    this.boldDatabase = 'public.tax-derep',
    this.boldMinIdentity = 0.9,
    this.boldMinOverlap = 100,
    this.boldMaxHits = 25,
  });

  /// BLAST databases offered (core_nt: the default nucleotide database of
  /// web BLAST since 2024; nt: the complete one).
  static const blastDatabases = ['core_nt', 'nt', 'refseq_select_rna', 'refseq_rna'];
  static const blastPrograms = ['megablast', 'blastn'];

  /// BOLD ID engine databases and the markers each one identifies.
  static const boldDatabases = {
    'public.tax-derep': ['COI'],
    'species': ['COI'],
    'all.tax-derep': ['COI'],
    'all.animal-alt': ['COI'],
    'public.plants': ['rbcL', 'matK'],
    'public.fungi': ['ITS'],
  };

  final String blastDatabase, blastProgram, entrezQuery, boldDatabase;
  final double expect, minCoverage, boldMinIdentity;
  final int? wordSize;
  final int hitlistSize, boldMinOverlap, boldMaxHits;

  /// Markers the chosen BOLD database identifies.
  List<String> get boldGenes => boldDatabases[boldDatabase] ?? const ['COI'];
  bool boldAccepts(String gene) => boldGenes.any((g) => g.toLowerCase() == gene.toLowerCase());

  /// What changes the answers: part of the cache keys.
  String get blastSignature => [blastDatabase, blastProgram, expect, wordSize ?? '', hitlistSize, entrezQuery.trim()].join('~');
  String get boldSignature => [boldDatabase, boldMinIdentity, boldMinOverlap, boldMaxHits].join('~');

  SpeciesIdOptions copyWith({
    String? blastDatabase,
    String? blastProgram,
    double? expect,
    int? wordSize,
    bool clearWordSize = false,
    int? hitlistSize,
    String? entrezQuery,
    double? minCoverage,
    String? boldDatabase,
    double? boldMinIdentity,
    int? boldMinOverlap,
    int? boldMaxHits,
  }) =>
      SpeciesIdOptions(
        blastDatabase: blastDatabase ?? this.blastDatabase,
        blastProgram: blastProgram ?? this.blastProgram,
        expect: expect ?? this.expect,
        wordSize: clearWordSize ? null : (wordSize ?? this.wordSize),
        hitlistSize: hitlistSize ?? this.hitlistSize,
        entrezQuery: entrezQuery ?? this.entrezQuery,
        minCoverage: minCoverage ?? this.minCoverage,
        boldDatabase: boldDatabase ?? this.boldDatabase,
        boldMinIdentity: boldMinIdentity ?? this.boldMinIdentity,
        boldMinOverlap: boldMinOverlap ?? this.boldMinOverlap,
        boldMaxHits: boldMaxHits ?? this.boldMaxHits,
      );

  Map<String, dynamic> toJson() => {
        'blastDatabase': blastDatabase,
        'blastProgram': blastProgram,
        'expect': expect,
        'wordSize': wordSize,
        'hitlistSize': hitlistSize,
        'entrezQuery': entrezQuery,
        'minCoverage': minCoverage,
        'boldDatabase': boldDatabase,
        'boldMinIdentity': boldMinIdentity,
        'boldMinOverlap': boldMinOverlap,
        'boldMaxHits': boldMaxHits,
      };

  factory SpeciesIdOptions.fromJson(Map<String, dynamic> j) {
    const d = SpeciesIdOptions();
    T v<T>(String k, T def) => j[k] is T ? j[k] as T : def;
    return SpeciesIdOptions(
      blastDatabase: blastDatabases.contains(j['blastDatabase']) ? j['blastDatabase'] as String : d.blastDatabase,
      blastProgram: blastPrograms.contains(j['blastProgram']) ? j['blastProgram'] as String : d.blastProgram,
      expect: (j['expect'] as num?)?.toDouble() ?? d.expect,
      wordSize: (j['wordSize'] as num?)?.toInt(),
      hitlistSize: (j['hitlistSize'] as num?)?.toInt() ?? d.hitlistSize,
      entrezQuery: v('entrezQuery', ''),
      minCoverage: (j['minCoverage'] as num?)?.toDouble() ?? d.minCoverage,
      boldDatabase: boldDatabases.containsKey(j['boldDatabase']) ? j['boldDatabase'] as String : d.boldDatabase,
      boldMinIdentity: (j['boldMinIdentity'] as num?)?.toDouble() ?? d.boldMinIdentity,
      boldMinOverlap: (j['boldMinOverlap'] as num?)?.toInt() ?? d.boldMinOverlap,
      boldMaxHits: (j['boldMaxHits'] as num?)?.toInt() ?? d.boldMaxHits,
    );
  }
}

/// One hit of a query: who it is, how similar, where it aligns and (when
/// the source gives it) the aligned sequences.
class IdHit {
  IdHit({
    required this.source,
    required this.species,
    required this.accession,
    required this.identity,
    this.title = '',
    this.alignLength = 0,
    this.qFrom = 0,
    this.qTo = 0,
    this.sFrom = 0,
    this.sTo = 0,
    this.mismatches,
    this.gaps = 0,
    this.evalue,
    this.bitScore,
    this.qseq,
    this.hseq,
  });

  final IdSource source;
  final String species, accession, title;
  final double identity;
  final int alignLength, qFrom, qTo, sFrom, sTo, gaps;
  final int? mismatches;
  final double? evalue, bitScore;

  /// Aligned query and hit (with gaps): BLAST gives them; for BOLD the hit
  /// sequence is fetched and aligned when the alignment is opened.
  String? qseq, hseq;

  bool get hasAlignment => qseq != null && hseq != null && qseq!.length == hseq!.length;

  /// The hit runs on the reverse strand of the query.
  bool get reversed => sTo < sFrom;

  Map<String, dynamic> toJson({bool withAlignment = true}) => {
        's': source.name,
        'sp': species,
        'acc': accession,
        'id': identity,
        if (title.isNotEmpty) 't': title,
        'len': alignLength,
        'q': [qFrom, qTo],
        'h': [sFrom, sTo],
        if (mismatches != null) 'mm': mismatches,
        'gaps': gaps,
        if (evalue != null) 'e': evalue,
        if (bitScore != null) 'bits': bitScore,
        if (withAlignment && qseq != null) 'qseq': qseq,
        if (withAlignment && hseq != null) 'hseq': hseq,
      };

  factory IdHit.fromJson(Map<String, dynamic> j) {
    List<int> pair(String k) => ((j[k] as List?) ?? const [0, 0]).map((e) => (e as num).toInt()).toList();
    final q = pair('q'), h = pair('h');
    return IdHit(
      source: j['s'] == 'bold' ? IdSource.bold : IdSource.genbank,
      species: j['sp'] as String? ?? '',
      accession: j['acc'] as String? ?? '',
      identity: (j['id'] as num?)?.toDouble() ?? 0,
      title: j['t'] as String? ?? '',
      alignLength: (j['len'] as num?)?.toInt() ?? 0,
      qFrom: q[0],
      qTo: q[1],
      sFrom: h[0],
      sTo: h[1],
      mismatches: (j['mm'] as num?)?.toInt(),
      gaps: (j['gaps'] as num?)?.toInt() ?? 0,
      evalue: (j['e'] as num?)?.toDouble(),
      bitScore: (j['bits'] as num?)?.toDouble(),
      qseq: j['qseq'] as String?,
      hseq: j['hseq'] as String?,
    );
  }
}

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

  Map<String, dynamic> toJson() => {
        'id': identity,
        'n': hits,
        if (selfIdentity != null) 'self': selfIdentity,
        if (selfAccessions.isNotEmpty) 'selfAcc': selfAccessions.toList(),
      };

  factory IdEvidence.fromJson(Map<String, dynamic> j) => IdEvidence((j['id'] as num?)?.toDouble() ?? 0, (j['n'] as num?)?.toInt() ?? 0)
    ..selfIdentity = (j['self'] as num?)?.toDouble()
    ..selfAccessions.addAll([...?(j['selfAcc'] as List?)?.cast<String>()]);
}

/// A species proposed for a record, with the evidence of each source.
class SpeciesCandidate {
  SpeciesCandidate(this.species);
  final String species;

  /// Other names (synonyms in GBIF) the hits of this species came under,
  /// e.g. Bufo marinus and Chaunus marinus for Rhinella marina.
  final Set<String> foundAs = {};

  /// The hits behind this species (gene → hits, best first), for the
  /// alignment viewer.
  final Map<String, List<IdHit>> hits = {};

  void _keep(String gene, IdHit h) {
    final list = hits.putIfAbsent(gene, () => [])..add(h);
    list.sort((a, b) => b.identity.compareTo(a.identity));
  }

  /// Joins the evidence of [other] (a synonym) into this one.
  void absorb(SpeciesCandidate other) {
    foundAs
      ..add(other.species)
      ..addAll(other.foundAs);
    for (final e in other.genbank.entries) {
      genbank.update(e.key, (mine) => _merge(mine, e.value), ifAbsent: () => e.value);
    }
    if (other.bold != null) bold = bold == null ? other.bold : _merge(bold!, other.bold!);
    for (final e in other.hits.entries) {
      for (final h in e.value) {
        _keep(e.key, h);
      }
    }
  }

  static IdEvidence _merge(IdEvidence a, IdEvidence b) {
    final out = IdEvidence(math.max(a.identity, b.identity), a.hits + b.hits)
      ..selfIdentity = (a.selfIdentity == null && b.selfIdentity == null) ? null : math.max(a.selfIdentity ?? 0, b.selfIdentity ?? 0);
    out.selfAccessions
      ..addAll(a.selfAccessions)
      ..addAll(b.selfAccessions);
    return out;
  }

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

  Map<String, dynamic> toJson() => {
        'species': species,
        if (foundAs.isNotEmpty) 'foundAs': foundAs.toList(),
        'genbank': {for (final e in genbank.entries) e.key: e.value.toJson()},
        if (bold != null) 'bold': bold!.toJson(),
        'hits': {
          for (final e in hits.entries) e.key: [for (final (i, h) in e.value.indexed) h.toJson(withAlignment: i < 5)],
        },
      };

  factory SpeciesCandidate.fromJson(Map<String, dynamic> j) {
    final c = SpeciesCandidate(j['species'] as String)..foundAs.addAll([...?(j['foundAs'] as List?)?.cast<String>()]);
    for (final e in ((j['genbank'] as Map?) ?? const {}).entries) {
      c.genbank[e.key as String] = IdEvidence.fromJson((e.value as Map).cast<String, dynamic>());
    }
    if (j['bold'] is Map) c.bold = IdEvidence.fromJson((j['bold'] as Map).cast<String, dynamic>());
    for (final e in ((j['hits'] as Map?) ?? const {}).entries) {
      c.hits[e.key as String] = [for (final h in e.value as List) IdHit.fromJson((h as Map).cast<String, dynamic>())];
    }
    return c;
  }
}

class RecordIdentification {
  RecordIdentification(this.accession, this.currentName);
  final String accession;
  final String currentName;

  /// The sequence sent for each gene (the query of the alignments).
  final Map<String, String> queries = {};
  final Map<String, SpeciesCandidate> _candidates = {};

  /// The accepted name (GBIF) of [currentName] when it is a synonym.
  String? currentAccepted;

  /// [name] is the record's species: its name or, when that is a synonym,
  /// the accepted one.
  bool isCurrent(String name) => sameSpecies(name, currentName) || (currentAccepted != null && sameSpecies(name, currentAccepted!));

  /// Groups the candidates by accepted name: [acceptedOf] maps a name
  /// (any case) to its accepted name in GBIF; a synonym's evidence joins
  /// the accepted species, which remembers the names it came under.
  void mergeSynonyms(Map<String, String> acceptedOf) {
    final byNorm = {for (final e in acceptedOf.entries) _norm(e.key): e.value};
    String? accepted(String name) => byNorm[_norm(binomialOf(name))];
    final cur = accepted(currentName);
    if (cur != null && !sameSpecies(cur, binomialOf(currentName))) currentAccepted = cur;
    for (final c in _candidates.values.toList()) {
      final a = accepted(c.species);
      if (a == null || sameSpecies(a, binomialOf(c.species))) continue;
      _candidates.remove(c.species);
      final target = _candidates.values.where((x) => sameSpecies(x.species, a)).firstOrNull ?? _candidate(a);
      target.absorb(c);
    }
  }

  /// Sources that answered for this record.
  final Set<IdSource> sources = {};

  SpeciesCandidate _candidate(String species) => _candidates.putIfAbsent(species, () => SpeciesCandidate(species));

  /// Counts one hit of [species] at [identity] % from [source] ([gene] for
  /// GenBank; BOLD is always COI). A hit on the record itself ([self], its
  /// accession [acc]) is kept apart: shown, but not counted as evidence.
  void addEvidence(IdSource source, String gene, String species, double identity, {bool self = false, String acc = '', IdHit? hit}) {
    sources.add(source);
    final c = _candidate(species);
    if (hit != null) c._keep(gene, hit);
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
  SpeciesCandidate? get current => candidates.where((c) => isCurrent(c.species)).firstOrNull;

  /// The best candidate names another species than the record, at an
  /// identity that means the same species. Below that the closest match is
  /// only a relative (expected when no other record of the species exists),
  /// which does not contest the name.
  bool get disagrees {
    final c = candidates;
    if (c.isEmpty || isCurrent(c.first.species) || c.first.bestIdentity < sameSpeciesIdentity) return false;
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

  Map<String, dynamic> toJson() => {
        'accession': accession,
        'currentName': currentName,
        if (currentAccepted != null) 'currentAccepted': currentAccepted,
        'queries': queries,
        'sources': [for (final x in sources) x.name],
        'candidates': [for (final c in _candidates.values) c.toJson()],
      };

  factory RecordIdentification.fromJson(Map<String, dynamic> j) {
    final r = RecordIdentification(j['accession'] as String, j['currentName'] as String? ?? '')
      ..currentAccepted = j['currentAccepted'] as String?;
    r.queries.addAll(((j['queries'] as Map?) ?? const {}).cast<String, String>());
    for (final x in (j['sources'] as List?) ?? const []) {
      if (IdSource.values.where((v) => v.name == x).firstOrNull case final v?) r.sources.add(v);
    }
    for (final c in (j['candidates'] as List?) ?? const []) {
      final k = SpeciesCandidate.fromJson((c as Map).cast<String, dynamic>());
      r._candidates[k.species] = k;
    }
    return r;
  }
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


class SpeciesIdentifier {
  SpeciesIdentifier({this.onStatus, this.onProgress, this.options = const SpeciesIdOptions()});

  final SpeciesIdOptions options;

  /// Progress messages for the modal.
  final void Function(String message)? onStatus;

  /// Share of each source's work done (0–1). BLAST advances by stage
  /// (submitted, waiting against NCBI's own time estimate, results); BOLD
  /// by sequences completed.
  final void Function(IdSource source, double fraction)? onProgress;
  final _http = HttpClient()..userAgent = 'SPLACE/4 (https://github.com/itvgenomics/splace)';

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

  /// BLAST hits per query id, with their alignments.
  Future<Map<String, List<IdHit>>> blast(List<IdQuery> queries, {CancelToken? cancel}) async {
    final out = <String, List<IdHit>>{};
    final o = options;
    for (var i = 0; i < queries.length; i += blastChunk) {
      final chunk = queries.sublist(i, math.min(i + blastChunk, queries.length));
      final part = i ~/ blastChunk + 1, parts = (queries.length / blastChunk).ceil();
      void report(double within) => onProgress?.call(IdSource.genbank, (part - 1 + within) / parts);
      report(0.02);
      final body = Uri(queryParameters: {
        'CMD': 'Put',
        'PROGRAM': 'blastn',
        if (o.blastProgram == 'megablast') 'MEGABLAST': 'on',
        'DATABASE': o.blastDatabase,
        'HITLIST_SIZE': '${o.hitlistSize}',
        'EXPECT': '${o.expect}',
        if (o.wordSize != null) 'WORD_SIZE': '${o.wordSize}',
        if (o.entrezQuery.trim().isNotEmpty) 'ENTREZ_QUERY': o.entrezQuery.trim(),
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
          if (qlen > 0 && len / qlen < 0.3) continue;
          final ident = (hsp['identity'] as num).toInt();
          final identity = 100 * ident / len;
          final gaps = (hsp['gaps'] as num?)?.toInt() ?? 0;
          for (final desc in descs) {
            final name = (desc['sciname'] as String?) ?? '';
            if (name.isEmpty) continue;
            hits.add(IdHit(
              source: IdSource.genbank,
              species: binomialOf(name),
              accession: (desc['accession'] as String?) ?? '',
              title: (desc['title'] as String?) ?? '',
              identity: identity,
              alignLength: len,
              qFrom: (hsp['query_from'] as num?)?.toInt() ?? 0,
              qTo: (hsp['query_to'] as num?)?.toInt() ?? 0,
              sFrom: (hsp['hit_from'] as num?)?.toInt() ?? 0,
              sTo: (hsp['hit_to'] as num?)?.toInt() ?? 0,
              mismatches: len - ident - gaps,
              gaps: gaps,
              evalue: (hsp['evalue'] as num?)?.toDouble(),
              bitScore: (hsp['bit_score'] as num?)?.toDouble(),
              qseq: hsp['qseq'] as String?,
              hseq: hsp['hseq'] as String?,
            ));
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

  /// BOLD hits per query id. BOLD gives the coordinates, not the aligned
  /// sequences: the viewer fetches the hit and aligns it.
  Future<Map<String, List<IdHit>>> bold(List<IdQuery> queries, {CancelToken? cancel}) async {
    final out = <String, List<IdHit>>{};
    final o = options;
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
        Uri.parse('$_bold/submission').replace(queryParameters: {
          'db': o.boldDatabase,
          'mi': '${o.boldMinIdentity}',
          'mo': '${o.boldMinOverlap}',
          'maxh': '${o.boldMaxHits}',
          'order': '3',
        }),
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
          // "ABTVC001-05|COI-5P|BOLD:AAB1186|…"
          final key = (e.key as String).split('|');
          int n(String k) => (v[k] as num?)?.toInt() ?? 0;
          hits.add(IdHit(
            source: IdSource.bold,
            species: binomialOf(species),
            accession: key.first,
            title: [key.first, if (key.length > 1) key[1], if (key.length > 2 && key[2].startsWith('BOLD:')) key[2]].join(' · '),
            identity: (v['pident'] as num).toDouble(),
            alignLength: n('length'),
            qFrom: n('qstart'),
            qTo: n('qend'),
            sFrom: n('sstart'),
            sTo: n('send'),
            mismatches: n('mismatch'),
            gaps: n('gaps'),
            evalue: (v['evalue'] as num?)?.toDouble(),
            bitScore: (v['bitscore'] as num?)?.toDouble(),
          ));
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
  /// The parameters of the source are part of the key: changing them
  /// sends the queries again.
  static String cacheKey(IdSource source, String record, String gene, String sequence, [SpeciesIdOptions options = const SpeciesIdOptions()]) =>
      '${source == IdSource.genbank ? 'genbank3~${options.blastSignature}' : 'bold3~${options.boldSignature}'}|$record|$gene|${fingerprint([SeqRecord('', sequence)])}';

  /// Cache entries no longer produced by [cacheKey] (no hit details).
  static bool staleCacheKey(String key) => !(key.startsWith('genbank3~') || key.startsWith('bold3~'));

  /// Identifies [records] with the [genes] given (name → sequence per
  /// record): every gene goes to GenBank, COI also to BOLD. Queries found in
  /// [cache] (see [cacheKey]) are not sent again; new answers are added.
  Future<Map<String, RecordIdentification>> identify(
    List<GbRecord> records,
    Map<String, Map<String, String>> sequences, {
    CancelToken? cancel,
    void Function(String message)? log,
    Map<String, List<IdHit>>? cache,
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

    for (final q in queries) {
      ids[q.record]!.queries[q.gene] = q.sequence;
    }

    void add(IdSource source, List<IdQuery> qs, Map<String, List<IdHit>> hits) {
      for (final q in qs) {
        final list = hits[q.id];
        if (list == null) continue;
        final rid = ids[q.record]!..sources.add(source);
        for (final h in list) {
          // Short hits (a conserved stretch) do not name the species: the
          // minimum coverage is applied here, so changing it needs no new
          // search.
          if (source == IdSource.genbank && q.sequence.isNotEmpty && h.alignLength / q.sequence.length < options.minCoverage) continue;
          // The record itself is shown but is not evidence.
          final isSelf = source == IdSource.genbank && (self[q.record]?.contains(base(h.accession)) ?? false);
          rid.addEvidence(source, q.gene, h.species, h.identity, self: isSelf, acc: h.accession, hit: h);
        }
      }
    }

    String key(IdSource s, IdQuery q) => cacheKey(s, q.record, q.gene, q.sequence, options);

    /// Answers [qs] from the cache and sends the rest with [run].
    Future<void> source(IdSource s, List<IdQuery> qs, Future<Map<String, List<IdHit>>> Function(List<IdQuery>) run) async {
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
    // BOLD: the markers its chosen database identifies (COI; rbcL and matK
    // in the plant database; ITS in the fungal one).
    final coi = queries.where((q) => options.boldAccepts(q.gene)).toList();
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
