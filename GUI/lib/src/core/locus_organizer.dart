/// Reorganizes FASTA records into one file per locus.
///
/// Input data often comes "sample-centric" (one file per specimen holding
/// several loci, e.g. `NC_030596_Oryza_brachyantha.fasta` with `>psbA`,
/// `>matK`, ...) or with sample and locus mixed in the header
/// (`>Sample12_COI`). Phylogenetic pipelines need "locus-centric" files
/// whose headers are identical across loci so they can be aligned and
/// concatenated. This module extracts a sample name and a locus name from
/// every record according to user-chosen [FieldRule]s, normalizes locus
/// synonyms, resolves duplicates and builds the locus × sample matrix.
library;

import 'dart:math' as math;

import 'fasta.dart';

class FastaInput {
  FastaInput(this.fileName, this.records);

  /// File name (not the full path); used as a field source.
  final String fileName;
  final List<SeqRecord> records;

  String get stem => fileStem(fileName);
}

const _fastaExtensions = [
  '.fasta', '.fas', '.fa', '.fna', '.ffn', '.faa', '.fsa', '.txt', '.aln',
];

/// File name without directory and without FASTA extension.
String fileStem(String fileName) {
  var name = fileName.split(RegExp(r'[\\/]')).last;
  final lower = name.toLowerCase();
  for (final ext in _fastaExtensions) {
    if (lower.endsWith(ext)) {
      name = name.substring(0, name.length - ext.length);
      break;
    }
  }
  return name;
}

/// Splits [text] on any character contained in [delimiters]. Empty tokens
/// are dropped so `a__b` yields `[a, b]`.
List<String> tokenize(String text, String delimiters) {
  if (delimiters.isEmpty) return [text];
  final tokens = <String>[];
  final current = StringBuffer();
  for (final rune in text.runes) {
    final ch = String.fromCharCode(rune);
    if (delimiters.contains(ch)) {
      if (current.isNotEmpty) tokens.add(current.toString());
      current.clear();
    } else {
      current.write(ch);
    }
  }
  if (current.isNotEmpty) tokens.add(current.toString());
  return tokens;
}

final _bracketTag = RegExp(r'\[([A-Za-z_][\w.-]*)=([^\]]*)\]');
final _plainTag = RegExp(r'(?:^|[\s|;])([A-Za-z_][\w.-]*)=([^\s|;\[\]]+)');

/// Extracts `key=value` annotations such as NCBI's `[gene=COX1]` or
/// `locus=COI|sample=X`. Keys are lower-cased; bracketed values may contain
/// spaces, plain ones may not.
Map<String, String> parseHeaderTags(String header) {
  final tags = <String, String>{};
  for (final m in _bracketTag.allMatches(header)) {
    tags.putIfAbsent(m.group(1)!.toLowerCase(), () => m.group(2)!.trim());
  }
  for (final m in _plainTag.allMatches(header.replaceAll(_bracketTag, ' '))) {
    tags.putIfAbsent(m.group(1)!.toLowerCase(), () => m.group(2)!.trim());
  }
  return tags;
}

enum FieldSource { fileName, header, tag, regex }

/// Describes how to extract one field (sample or locus) from a record.
///
/// For [FieldSource.fileName] and [FieldSource.header] the text is split on
/// [delimiters] and the tokens `[start, end)` are joined with `_`. Negative
/// indices count from the end, `end == null` means "to the last token"; so
/// `start: -1` is the last token and `start: 0, end: -1` is everything but
/// the last token.
class FieldRule {
  const FieldRule({
    required this.source,
    this.start = 0,
    this.end,
    this.delimiters = '_',
    this.tag = 'gene',
    this.pattern = '',
    this.regexOnFileName = false,
  });

  const FieldRule.fileName({int start = 0, int? end, String delimiters = '_'})
      : this(source: FieldSource.fileName, start: start, end: end, delimiters: delimiters);

  const FieldRule.header({int start = 0, int? end, String delimiters = ' '})
      : this(source: FieldSource.header, start: start, end: end, delimiters: delimiters);

  const FieldRule.tag(String tag) : this(source: FieldSource.tag, tag: tag);

  const FieldRule.regex(String pattern, {bool onFileName = false})
      : this(source: FieldSource.regex, pattern: pattern, regexOnFileName: onFileName);

  final FieldSource source;
  final int start;
  final int? end;
  final String delimiters;
  final String tag;

  /// Regular expression whose first capture group (or whole match when there
  /// is no group) becomes the field value.
  final String pattern;
  final bool regexOnFileName;

  FieldRule copyWith({
    FieldSource? source,
    int? start,
    int? Function()? end,
    String? delimiters,
    String? tag,
    String? pattern,
    bool? regexOnFileName,
  }) =>
      FieldRule(
        source: source ?? this.source,
        start: start ?? this.start,
        end: end != null ? end() : this.end,
        delimiters: delimiters ?? this.delimiters,
        tag: tag ?? this.tag,
        pattern: pattern ?? this.pattern,
        regexOnFileName: regexOnFileName ?? this.regexOnFileName,
      );

  /// Returns null when the field cannot be extracted.
  String? extract(String fileName, SeqRecord record) {
    switch (source) {
      case FieldSource.fileName:
        return _slice(tokenize(fileStem(fileName), delimiters));
      case FieldSource.header:
        return _slice(tokenize(record.header, delimiters));
      case FieldSource.tag:
        final v = parseHeaderTags(record.header)[tag.toLowerCase()];
        return (v == null || v.isEmpty) ? null : v;
      case FieldSource.regex:
        if (pattern.isEmpty) return null;
        final RegExp re;
        try {
          re = RegExp(pattern);
        } on FormatException {
          return null;
        }
        final m = re.firstMatch(regexOnFileName ? fileStem(fileName) : record.header);
        if (m == null) return null;
        final v = m.groupCount >= 1 ? m.group(1) : m.group(0);
        return (v == null || v.isEmpty) ? null : v;
    }
  }

  String? _slice(List<String> tokens) {
    final n = tokens.length;
    if (n == 0) return null;
    final s = start < 0 ? n + start : start;
    final e = end == null ? n : (end! < 0 ? n + end! : end!);
    if (s < 0 || s >= n || e <= s) return null;
    return tokens.sublist(s, math.min(e, n)).join('_');
  }

  String describe() {
    switch (source) {
      case FieldSource.fileName:
      case FieldSource.header:
        final where = source == FieldSource.fileName ? 'file name' : 'header';
        final range = end == null ? '[$start:]' : '[$start:$end]';
        return '$where tokens $range split on "${delimiters.replaceAll(' ', '␣')}"';
      case FieldSource.tag:
        return 'header tag "$tag="';
      case FieldSource.regex:
        return 'regex /$pattern/ on ${regexOnFileName ? 'file name' : 'header'}';
    }
  }
}

/// Canonical names for common animal mitochondrial loci. Keys are
/// normalization keys (see [locusKey]).
const Map<String, String> mitochondrialSynonyms = {
  'COI': 'COI', 'CO1': 'COI', 'COX1': 'COI', 'COXI': 'COI',
  'COII': 'COII', 'CO2': 'COII', 'COX2': 'COII', 'COXII': 'COII',
  'COIII': 'COIII', 'CO3': 'COIII', 'COX3': 'COIII', 'COXIII': 'COIII',
  'CYTB': 'CYTB', 'COB': 'CYTB', 'CYB': 'CYTB',
  'ND1': 'ND1', 'NAD1': 'ND1', 'NADH1': 'ND1',
  'ND2': 'ND2', 'NAD2': 'ND2', 'NADH2': 'ND2',
  'ND3': 'ND3', 'NAD3': 'ND3', 'NADH3': 'ND3',
  'ND4': 'ND4', 'NAD4': 'ND4', 'NADH4': 'ND4',
  'ND4L': 'ND4L', 'NAD4L': 'ND4L', 'NADH4L': 'ND4L',
  'ND5': 'ND5', 'NAD5': 'ND5', 'NADH5': 'ND5',
  'ND6': 'ND6', 'NAD6': 'ND6', 'NADH6': 'ND6',
  'ATP6': 'ATP6', 'ATPASE6': 'ATP6',
  'ATP8': 'ATP8', 'ATPASE8': 'ATP8',
  '12S': '12S', '12SRRNA': '12S', 'RRNS': '12S', 'SRRNA': '12S',
  '16S': '16S', '16SRRNA': '16S', 'RRNL': '16S', 'LRRNA': '16S',
};

/// Case/punctuation-insensitive key used to group locus names.
String locusKey(String name) =>
    name.toUpperCase().replaceAll(RegExp(r'[\s_\-.]'), '');

final _geneLike = RegExp(
  r'^(COI{1,3}|CO[123]|COX[123I]+|CYTB|COB|ND[1-6]L?|NAD[1-9]L?|ATP\d+|ATPASE\d|'
  r'\d{1,2}S(RRNA)?|RRN\d*[SL]?|ITS\d?|5\.8S|18S|28S|H3|EF1A?|RAG[12]|'
  r'(PSA|PSB|PET|ATP|NDH|RPS|RPL|RPO|YCF|MAT|RBC|ACC|CCS|CEM|CLP|INF|LHB|RRN|TRN|CHL|PAF|PBF)[A-Z0-9]*)$',
  caseSensitive: false,
);

/// Heuristic: does [token] look like a gene/locus name?
bool looksLikeLocus(String token) =>
    _geneLike.hasMatch(token.replaceAll(RegExp(r'[\s_\-]'), ''));

class LocusNormalizer {
  LocusNormalizer({Map<String, String>? synonyms, this.useSynonyms = true})
      : synonyms = synonyms ?? mitochondrialSynonyms;

  final Map<String, String> synonyms;
  final bool useSynonyms;

  /// Grouping key for a raw locus name.
  String keyOf(String raw) {
    final key = locusKey(raw);
    if (useSynonyms) {
      final canonical = synonyms[key];
      if (canonical != null) return locusKey(canonical);
    }
    return key;
  }

  String? canonicalOf(String raw) =>
      useSynonyms ? synonyms[locusKey(raw)] : null;
}

enum DuplicatePolicy { keepLongest, keepMostComplete, keepFirst, keepAllRenamed }

class OrganizedSeq {
  OrganizedSeq({
    required this.sample,
    required this.locus,
    required this.record,
    required this.sourceFile,
  });

  final String sample;
  final String locus;
  final SeqRecord record;
  final String sourceFile;

  /// Count of unambiguous nucleotides; used to rank duplicates.
  int get informativeLength =>
      RegExp('[ACGTUacgtu]').allMatches(record.sequence).length;
}

class DuplicateGroup {
  DuplicateGroup(this.sample, this.locus, this.candidates, this.kept);
  final String sample;
  final String locus;
  final List<OrganizedSeq> candidates;
  final List<OrganizedSeq> kept;
}

class UnassignedRecord {
  UnassignedRecord(this.sourceFile, this.header, {required this.missingSample, required this.missingLocus});
  final String sourceFile;
  final String header;
  final bool missingSample;
  final bool missingLocus;
}

final _word = RegExp(r'[A-Za-z0-9]+');

/// Structure of a header: its words collapsed, the delimiters kept
/// (`ITV01-L12` -> `x-x`, `L12|ITV01|2020` -> `x|x|x`).
String headerShape(String header) => header.replaceAll(_word, 'x');

/// Structure of a name: letters kept (upper case), each digit as `#`
/// (`ITV01054` -> `ITV#####`, `ITV1054` -> `ITV####`).
String nameStructure(String name) => name.toUpperCase().replaceAll(RegExp(r'\d'), '#');

/// [name] written with the structure [target] when they differ only in how
/// many digits a number has: zeros are added or removed on the left
/// (`ITV1054` -> `ITV01054` for `ITV#####`). Null when it cannot be done.
String? conformName(String name, String target) {
  final parts = RegExp(r'\d+|\D+').allMatches(name).map((m) => m[0]!).toList();
  final want = RegExp(r'#+|[^#]+').allMatches(target).map((m) => m[0]!).toList();
  if (parts.length != want.length) return null;
  final out = StringBuffer();
  for (var i = 0; i < parts.length; i++) {
    final p = parts[i], w = want[i];
    final digits = RegExp(r'^\d+$').hasMatch(p);
    if (digits != w.startsWith('#')) return null;
    if (!digits) {
      if (p.toUpperCase() != w) return null;
      out.write(p);
      continue;
    }
    final core = p.replaceFirst(RegExp(r'^0+(?=\d)'), '');
    if (core.length > w.length) return null;
    out.write(core.padLeft(w.length, '0'));
  }
  return out.toString();
}

/// How the headers of one structure ([headerShape]) are read, when it is
/// not the main rule: its own sample and marker parts, or left out.
class FormatRule {
  const FormatRule({required this.sample, required this.locus, this.skip = false});
  final FieldRule sample;
  final FieldRule locus;

  /// The sequences with this structure are not used.
  final bool skip;
}

class OrganizerConfig {
  const OrganizerConfig({
    required this.sample,
    required this.locus,
    this.duplicatePolicy = DuplicatePolicy.keepLongest,
    this.useSynonyms = true,
    this.sanitizeNames = true,
    this.ignoreCase = true,
    this.ignoreLeadingZeros = false,
    this.sampleAliases = const {},
    this.formats = const {},
  });

  /// Header structure ([headerShape]) -> how those headers are read.
  final Map<String, FormatRule> formats;

  final FieldRule sample;
  final FieldRule locus;
  final DuplicatePolicy duplicatePolicy;
  final bool useSynonyms;
  final bool sanitizeNames;

  /// `itv01054` and `ITV01054` are the same sample (marker names are always
  /// compared without case).
  final bool ignoreCase;

  /// Numbers are compared without leading zeros, in samples and markers:
  /// `ITV1054` = `ITV01054`, `L01` = `L1`.
  final bool ignoreLeadingZeros;

  /// Sample names typed by the user: name read from the data -> name to use.
  /// Two names given the same new name are merged.
  final Map<String, String> sampleAliases;

  OrganizerConfig copyWith({
    FieldRule? sample,
    FieldRule? locus,
    DuplicatePolicy? duplicatePolicy,
    bool? useSynonyms,
    bool? sanitizeNames,
    bool? ignoreCase,
    bool? ignoreLeadingZeros,
    Map<String, String>? sampleAliases,
    Map<String, FormatRule>? formats,
  }) =>
      OrganizerConfig(
        sample: sample ?? this.sample,
        locus: locus ?? this.locus,
        duplicatePolicy: duplicatePolicy ?? this.duplicatePolicy,
        useSynonyms: useSynonyms ?? this.useSynonyms,
        sanitizeNames: sanitizeNames ?? this.sanitizeNames,
        ignoreCase: ignoreCase ?? this.ignoreCase,
        ignoreLeadingZeros: ignoreLeadingZeros ?? this.ignoreLeadingZeros,
        sampleAliases: sampleAliases ?? this.sampleAliases,
        formats: formats ?? this.formats,
      );
}

final _chunks = RegExp(r'\d+|\D+');

/// Case-insensitive order that reads numbers as numbers: L1, L2, L10, L100.
int naturalCompare(String a, String b) {
  final x = _chunks.allMatches(a.toLowerCase()).map((m) => m[0]!).toList();
  final y = _chunks.allMatches(b.toLowerCase()).map((m) => m[0]!).toList();
  for (var i = 0; i < x.length && i < y.length; i++) {
    final p = int.tryParse(x[i]), q = int.tryParse(y[i]);
    final c = p != null && q != null ? (p != q ? p.compareTo(q) : x[i].length.compareTo(y[i].length)) : x[i].compareTo(y[i]);
    if (c != 0) return c;
  }
  final c = x.length.compareTo(y.length);
  return c != 0 ? c : a.compareTo(b);
}

final _leadingZeros = RegExp(r'(?<!\d)0+(?=\d)');

/// [text] with the leading zeros of its numbers removed (`L007` -> `L7`).
String stripLeadingZeros(String text) => text.replaceAll(_leadingZeros, '');

class OrganizerResult {
  OrganizerResult({
    required this.samples,
    required this.loci,
    required this.byLocus,
    required this.duplicates,
    required this.unassigned,
    required this.locusAliases,
    this.sampleSpellings = const {},
    this.sampleFiles = const {},
    this.records = 0,
    this.skipped = 0,
  });

  /// Records left out with their header structure ([FormatRule.skip]).
  final int skipped;

  /// Sample -> the spellings merged into it (only when there is more than
  /// one, e.g. ITV01054 <- ITV1054).
  final Map<String, Set<String>> sampleSpellings;

  /// Sample -> files its sequences came from.
  final Map<String, Set<String>> sampleFiles;


  /// Records read (assigned or not).
  final int records;

  /// Sorted sample names.
  final List<String> samples;

  /// Loci sorted by decreasing sample coverage, then name.
  final List<String> loci;

  /// locus -> sample -> kept sequence(s). Usually one per cell; more only
  /// with [DuplicatePolicy.keepAllRenamed] (then the sample name is
  /// suffixed and the extra entries become their own rows).
  final Map<String, Map<String, OrganizedSeq>> byLocus;
  final List<DuplicateGroup> duplicates;
  final List<UnassignedRecord> unassigned;

  /// canonical locus -> raw spellings merged into it (e.g. COI <- cox1, CO1).
  final Map<String, Set<String>> locusAliases;

  int samplesFor(String locus) => byLocus[locus]?.length ?? 0;

  int lociFor(String sample) =>
      byLocus.values.where((m) => m.containsKey(sample)).length;

  bool has(String locus, String sample) => byLocus[locus]?.containsKey(sample) ?? false;

  /// Fraction of filled cells in the locus × sample matrix.
  double get occupancy {
    if (samples.isEmpty || loci.isEmpty) return 0;
    final filled = byLocus.values.fold<int>(0, (a, m) => a + m.length);
    return filled / (samples.length * loci.length);
  }

  /// FASTA text per locus with headers = sample name, restricted to
  /// [loci] (default: all) and to loci present in at least [minSamples].
  Map<String, String> locusFastas({Iterable<String>? loci, int minSamples = 1}) {
    final out = <String, String>{};
    for (final locus in loci ?? this.loci) {
      final cells = byLocus[locus];
      if (cells == null || cells.length < minSamples) continue;
      final names = cells.keys.toList()..sort();
      out[locus] = writeFasta(names.map((s) => SeqRecord(s, cells[s]!.record.sequence)));
    }
    return out;
  }

  /// Mapping table original record -> new location, for provenance.
  String mappingTsv() {
    final b = StringBuffer('source_file\toriginal_header\tsample\tlocus\tlength\tstatus\n');
    final keptSet = <OrganizedSeq>{};
    for (final m in byLocus.values) {
      keptSet.addAll(m.values);
    }
    for (final m in byLocus.values) {
      for (final s in m.values) {
        b.writeln('${s.sourceFile}\t${s.record.header}\t${s.sample}\t${s.locus}\t${s.record.length}\tkept');
      }
    }
    for (final d in duplicates) {
      for (final c in d.candidates.where((c) => !keptSet.contains(c))) {
        b.writeln('${c.sourceFile}\t${c.record.header}\t${c.sample}\t${c.locus}\t${c.record.length}\tdropped_duplicate');
      }
    }
    for (final u in unassigned) {
      b.writeln('${u.sourceFile}\t${u.header}\t\t\t\tunassigned');
    }
    return b.toString();
  }
}

OrganizerResult organizeLoci(List<FastaInput> inputs, OrganizerConfig config) {
  final normalizer = LocusNormalizer(useSynonyms: config.useSynonyms);
  final unassigned = <UnassignedRecord>[];
  // key -> raw spellings with counts, to choose the display name.
  final spellings = <String, Map<String, int>>{};
  final sampleSpell = <String, Map<String, int>>{};
  final sampleFiles = <String, Set<String>>{};
  final cells = <String, Map<String, List<OrganizedSeq>>>{};
  var total = 0, skipped = 0;

  String sampleKey(String s) {
    var k = config.ignoreCase ? s.toUpperCase() : s;
    if (config.ignoreLeadingZeros) k = stripLeadingZeros(k);
    return k;
  }

  String locusGroup(String raw) => normalizer.keyOf(config.ignoreLeadingZeros ? stripLeadingZeros(raw) : raw);

  for (final input in inputs) {
    for (final record in input.records) {
      total++;
      var sampleRule = config.sample, locusRule = config.locus;
      if (config.formats.isNotEmpty) {
        final f = config.formats[headerShape(record.header)];
        if (f != null) {
          if (f.skip) {
            skipped++;
            continue;
          }
          sampleRule = f.sample;
          locusRule = f.locus;
        }
      }
      var sample = sampleRule.extract(input.fileName, record);
      final rawLocus = locusRule.extract(input.fileName, record);
      if (sample == null || rawLocus == null) {
        unassigned.add(UnassignedRecord(input.fileName, record.header,
            missingSample: sample == null, missingLocus: rawLocus == null));
        continue;
      }
      if (config.sanitizeNames) sample = sanitizeTaxonName(sample);
      sample = config.sampleAliases[sample] ?? sample;
      final sk = sampleKey(sample);
      (sampleSpell[sk] ??= {}).update(sample, (c) => c + 1, ifAbsent: () => 1);
      (sampleFiles[sk] ??= {}).add(input.fileName);
      final key = locusGroup(rawLocus);
      (spellings[key] ??= {}).update(rawLocus, (c) => c + 1, ifAbsent: () => 1);
      ((cells[key] ??= {})[sk] ??= []).add(OrganizedSeq(
        sample: sk,
        locus: key,
        record: record,
        sourceFile: input.fileName,
      ));
    }
  }

  String mostFrequent(Map<String, int> raw) {
    final ranked = raw.entries.toList()
      ..sort((a, b) => b.value != a.value ? b.value.compareTo(a.value) : a.key.compareTo(b.key));
    return ranked.first.key;
  }

  String displayName(String key) {
    final raw = spellings[key]!;
    final canonical = normalizer.canonicalOf(raw.keys.first);
    if (canonical != null) return canonical;
    // Most frequent spelling wins; ties broken alphabetically.
    final name = mostFrequent(raw);
    return config.sanitizeNames ? sanitizeTaxonName(name) : name;
  }

  // Sample keys -> the name shown (and written in the FASTA headers).
  final sampleName = {for (final e in sampleSpell.entries) e.key: mostFrequent(e.value)};

  final byLocus = <String, Map<String, OrganizedSeq>>{};
  final duplicates = <DuplicateGroup>[];
  final aliases = <String, Set<String>>{};

  for (final entry in cells.entries) {
    final locus = displayName(entry.key);
    (aliases[locus] ??= {}).addAll(spellings[entry.key]!.keys);
    final row = byLocus[locus] ??= {};
    for (final cell in entry.value.entries) {
      final name = sampleName[cell.key]!;
      final candidates = cell.value
          .map((s) => OrganizedSeq(sample: name, locus: locus, record: s.record, sourceFile: s.sourceFile))
          .toList();
      final kept = _resolveDuplicates(candidates, config.duplicatePolicy);
      if (candidates.length > 1) {
        duplicates.add(DuplicateGroup(name, locus, candidates, kept));
      }
      for (var i = 0; i < kept.length; i++) {
        row[i == 0 ? name : '${name}_dup${i + 1}'] = kept[i];
      }
    }
  }

  final samples = <String>{for (final row in byLocus.values) ...row.keys};

  final loci = byLocus.keys.toList()
    ..sort((a, b) {
      final c = byLocus[b]!.length.compareTo(byLocus[a]!.length);
      return c != 0 ? c : naturalCompare(a, b);
    });

  return OrganizerResult(
    samples: samples.toList()..sort(naturalCompare),
    loci: loci,
    byLocus: byLocus,
    duplicates: duplicates,
    unassigned: unassigned,
    locusAliases: aliases,
    sampleSpellings: {
      for (final e in sampleSpell.entries)
        if (e.value.length > 1) sampleName[e.key]!: e.value.keys.toSet(),
    },
    sampleFiles: {for (final e in sampleFiles.entries) sampleName[e.key]!: e.value},
    records: total,
    skipped: skipped,
  );
}

List<OrganizedSeq> _resolveDuplicates(List<OrganizedSeq> c, DuplicatePolicy policy) {
  if (c.length <= 1) return c;
  switch (policy) {
    case DuplicatePolicy.keepFirst:
      return [c.first];
    case DuplicatePolicy.keepAllRenamed:
      return c;
    case DuplicatePolicy.keepLongest:
      return [c.reduce((a, b) => b.record.length > a.record.length ? b : a)];
    case DuplicatePolicy.keepMostComplete:
      return [c.reduce((a, b) => b.informativeLength > a.informativeLength ? b : a)];
  }
}

// ---------------------------------------------------------------------------
// Layout detection and header exploration
// ---------------------------------------------------------------------------

enum LayoutKind {
  /// One file per sample, header = locus (needs reorganizing).
  samplePerFile,

  /// One file per locus, header = sample (already in the right shape).
  locusPerFile,

  /// Sample and locus are both encoded in each header.
  mixedHeader,

  /// Locus comes from a `[gene=...]`-style tag.
  taggedHeader,
  unknown,
}

class LayoutSuggestion {
  LayoutSuggestion(this.kind, this.config, this.confidence, this.evidence);
  final LayoutKind kind;
  final OrganizerConfig config;

  /// 0..1, how sure the heuristic is.
  final double confidence;

  /// Human-readable justification (English; shown in the UI log).
  final List<String> evidence;
}

double _fraction(Iterable<String> values, bool Function(String) test) {
  var n = 0, hits = 0;
  for (final v in values) {
    n++;
    if (test(v)) hits++;
  }
  return n == 0 ? 0 : hits / n;
}

/// Guesses which part of file names / headers holds the sample and which
/// holds the locus. The UI always lets the user override the result.
LayoutSuggestion suggestLayout(List<FastaInput> inputs) {
  final evidence = <String>[];
  final records = [for (final f in inputs) ...f.records];
  if (records.isEmpty) {
    return LayoutSuggestion(
      LayoutKind.unknown,
      const OrganizerConfig(sample: FieldRule.fileName(), locus: FieldRule.header()),
      0,
      ['No records'],
    );
  }

  // 1. Explicit tags win.
  final tagged = _fraction(records.take(500).map((r) => r.header),
      (h) => parseHeaderTags(h).containsKey('gene') || parseHeaderTags(h).containsKey('locus'));
  if (tagged > 0.8) {
    final key = records.take(500).any((r) => parseHeaderTags(r.header).containsKey('gene')) ? 'gene' : 'locus';
    final hasOrganism = _fraction(records.take(500).map((r) => r.header),
            (h) => parseHeaderTags(h).containsKey('organism')) > 0.8;
    evidence.add('${(tagged * 100).round()}% of headers carry a [$key=...] tag');
    return LayoutSuggestion(
      LayoutKind.taggedHeader,
      OrganizerConfig(
        sample: hasOrganism ? const FieldRule.tag('organism') : const FieldRule.fileName(),
        locus: FieldRule.tag(key),
      ),
      0.9,
      evidence,
    );
  }

  final ids = records.map((r) => r.id).toList();
  final distinctIds = ids.toSet();
  final idFiles = <String, int>{};
  for (final f in inputs) {
    for (final id in f.records.map((r) => r.id).toSet()) {
      idFiles.update(id, (v) => v + 1, ifAbsent: () => 1);
    }
  }
  final recurring = distinctIds.where((id) => (idFiles[id] ?? 0) > 1).length / distinctIds.length;

  // 2. Several files whose header ids recur across files.
  if (inputs.length > 1 && recurring > 0.5) {
    final idGeneLike = _fraction(distinctIds, looksLikeLocus);
    final stemGeneLike = _fraction(inputs.map((f) => f.stem), looksLikeLocus);
    evidence.add('${(recurring * 100).round()}% of header IDs appear in more than one file');
    evidence.add('${(idGeneLike * 100).round()}% of header IDs look like locus names; '
        '${(stemGeneLike * 100).round()}% of file names do');
    final samplePerFile = idGeneLike > stemGeneLike ||
        (idGeneLike == stemGeneLike && distinctIds.length <= inputs.length);
    if (samplePerFile) {
      return LayoutSuggestion(
        LayoutKind.samplePerFile,
        const OrganizerConfig(sample: FieldRule.fileName(delimiters: ''), locus: FieldRule.header()),
        0.5 + 0.5 * (idGeneLike - stemGeneLike).abs().clamp(0, 1),
        evidence,
      );
    }
    return LayoutSuggestion(
      LayoutKind.locusPerFile,
      const OrganizerConfig(sample: FieldRule.header(), locus: FieldRule.fileName(delimiters: '')),
      0.5 + 0.5 * (stemGeneLike - idGeneLike).abs().clamp(0, 1),
      evidence,
    );
  }

  // 3. Sample and marker in each header, told apart by how they behave:
  // the sample part is the same within a file and the marker part changes
  // within a file and repeats across files (`ERR2040879-L10`).
  final structural = _structuralGuess(inputs);
  if (structural != null) return structural;

  // 4. Sample and locus mixed in the header: find a gene-like token column.
  for (final delimiters in ['_', '|', '_|', ' _|-']) {
    final columns = headerTokenColumns(records.take(2000).map((r) => r.id).toList(), delimiters);
    final geneColumn = columns.where((c) => c.geneLikeFraction > 0.6).toList()
      ..sort((a, b) => b.geneLikeFraction.compareTo(a.geneLikeFraction));
    if (geneColumn.isEmpty) continue;
    final col = geneColumn.first;
    evidence.add('Token ${col.index} (split on "$delimiters") looks like a locus in '
        '${(col.geneLikeFraction * 100).round()}% of headers');
    final locus = FieldRule(source: FieldSource.header, start: col.index, end: col.index == -1 ? null : col.index + 1, delimiters: delimiters);
    final sample = col.index == 0
        ? FieldRule(source: FieldSource.header, start: 1, delimiters: delimiters)
        : col.index == -1
            ? FieldRule(source: FieldSource.header, start: 0, end: -1, delimiters: delimiters)
            : FieldRule(source: FieldSource.header, start: 0, end: col.index, delimiters: delimiters);
    return LayoutSuggestion(
      LayoutKind.mixedHeader,
      OrganizerConfig(sample: sample, locus: locus),
      col.geneLikeFraction * 0.9,
      evidence,
    );
  }

  evidence.add('Could not recognize locus names automatically; please choose the fields manually');
  return LayoutSuggestion(
    LayoutKind.unknown,
    const OrganizerConfig(sample: FieldRule.fileName(delimiters: ''), locus: FieldRule.header()),
    0.2,
    evidence,
  );
}

/// Delimiter sets tried when guessing the parts of a header.
const headerDelimiterCandidates = ['-', '_', '|', ' ', '.', '-_', '_|', '-_|', '-_| .'];

LayoutSuggestion? _structuralGuess(List<FastaInput> inputs) {
  final files = [for (final f in inputs) if (f.records.length > 1) f];
  if (files.length < 2) return null;
  // A sample of the records keeps this fast on large data sets.
  final perFile = math.max(20, 4000 ~/ files.length);
  for (final delimiters in headerDelimiterCandidates) {
    final tokens = [
      for (final f in files.take(200)) [for (final r in f.records.take(perFile)) tokenize(r.header, delimiters)],
    ];
    final lengths = {for (final f in tokens) for (final t in f) t.length};
    if (lengths.any((n) => n < 2)) continue;
    final fixed = lengths.length == 1;
    final width = lengths.reduce(math.min);
    // Positions from the start when every header has as many parts,
    // otherwise the first part and the last one.
    final positions = fixed ? [for (var i = 0; i < math.min(width, 8); i++) i] : [0, -1];
    String at(List<String> t, int i) => t[i < 0 ? t.length + i : i];
    int? sample, locus;
    var sampleScore = 0.0, locusScore = 0.0;
    for (final i in positions) {
      var constant = 0;
      var variety = 0.0;
      final filesOf = <String, int>{};
      final perFileValues = <String>{};
      for (final f in tokens) {
        final values = {for (final t in f) at(t, i)};
        if (values.length == 1) {
          constant++;
          perFileValues.add(values.first);
        }
        variety += values.length / f.length;
        for (final v in values) {
          filesOf.update(v, (n) => n + 1, ifAbsent: () => 1);
        }
      }
      // Constant within files and different between them: the sample.
      final c = constant / tokens.length * (perFileValues.length / math.max(1, constant));
      if (c > sampleScore) {
        sampleScore = c;
        sample = i;
      }
      // Changing within files and shared between them: the marker.
      final recurring = filesOf.values.where((n) => n > 1).length / filesOf.length;
      final l = (variety / tokens.length) * recurring;
      if (l > locusScore) {
        locusScore = l;
        locus = i;
      }
    }
    if (sample == null || locus == null || sample == locus || sampleScore < 0.8 || locusScore < 0.3) continue;
    final last = !fixed || locus == width - 1;
    final locusRule = FieldRule(source: FieldSource.header, start: last ? -1 : locus, end: last ? null : locus + 1, delimiters: delimiters);
    final sampleRule = last && sample == 0
        // Everything before the marker, so sample names that contain the
        // delimiter stay whole.
        ? FieldRule(source: FieldSource.header, start: 0, end: -1, delimiters: delimiters)
        : FieldRule(source: FieldSource.header, start: sample, end: sample == -1 ? null : sample + 1, delimiters: delimiters);
    return LayoutSuggestion(
      LayoutKind.mixedHeader,
      OrganizerConfig(sample: sampleRule, locus: locusRule),
      (0.5 + 0.5 * math.min(sampleScore, locusScore)).clamp(0.0, 0.95),
      [
        'Split on "${delimiters.replaceAll(' ', '␣')}": one part is the same in all headers of a file '
            '(sample) and another changes within files and repeats across them (marker)',
      ],
    );
  }
  return null;
}

class TokenColumn {
  TokenColumn(this.index, this.distinct, this.coverage, this.geneLikeFraction, this.examples);

  /// Token index; negative values count from the end (-1 = last token).
  final int index;
  final int distinct;

  /// Fraction of headers that have this token.
  final double coverage;
  final double geneLikeFraction;
  final List<String> examples;
}

/// Per-position statistics of header tokens; helps the user choose which
/// token is the sample and which is the locus. Returns positions from the
/// start (0, 1, ...) plus the last token (-1).
List<TokenColumn> headerTokenColumns(List<String> headers, String delimiters, {int maxColumns = 8}) {
  final tokenized = headers.map((h) => tokenize(h, delimiters)).toList();
  final maxLen = tokenized.fold<int>(0, (a, t) => math.max(a, t.length));
  final columns = <TokenColumn>[];
  for (final index in [for (var i = 0; i < math.min(maxLen, maxColumns); i++) i, if (maxLen > 1) -1]) {
    final values = <String>[];
    for (final t in tokenized) {
      final i = index < 0 ? t.length + index : index;
      if (i >= 0 && i < t.length) values.add(t[i]);
    }
    if (values.isEmpty) continue;
    final distinct = values.toSet();
    columns.add(TokenColumn(
      index,
      distinct.length,
      values.length / headers.length,
      _fraction(values, looksLikeLocus),
      distinct.take(5).toList(),
    ));
  }
  return columns;
}
