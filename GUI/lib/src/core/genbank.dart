/// GenBank flat-file parsing, feature locations and taxonomy, ported from
/// the HTML app (docs/splace.js) with a recursive location parser that also
/// handles `join(complement(..),..)` and `order(..)`.
///
/// Pure Dart so parsing can run in isolates.
library;

class GbFeature {
  GbFeature(this.type, this.location, this.qualifiers);
  final String type;
  final String location;
  final Map<String, String> qualifiers;

  String? q(String key) {
    final v = qualifiers[key]?.replaceAll(RegExp(r'\s+'), ' ').trim();
    return v == null || v.isEmpty ? null : v;
  }

  bool get isPseudo => qualifiers.containsKey('pseudo') || qualifiers.containsKey('pseudogene');
}

class TaxonomyRanks {
  const TaxonomyRanks({
    this.source = '',
    this.kingdom = '',
    this.phylum = '',
    this.klass = '',
    this.order = '',
    this.family = '',
    this.authorship = '',
  });
  final String kingdom, phylum, klass, order, family, authorship;

  /// Where the ranks came from: GBIF, NCBI, SALVE or manual.
  final String source;

  bool get isEmpty =>
      kingdom.isEmpty && phylum.isEmpty && klass.isEmpty && order.isEmpty && family.isEmpty && authorship.isEmpty;
}

class GbRecord {
  GbRecord({
    required this.accession,
    required this.organism,
    required this.lineage,
    required this.definition,
    required this.features,
    required this.sequence,
    this.source = 'file',
    this.fileName,
    this.fileSize,
    this.fileDate,
    this.sourceAccessions = const [],
  });

  final String accession;

  /// INSDC records a RefSeq record was copied from (its COMMENT: "The
  /// reference sequence is identical to KM201334" / "was derived from …").
  /// Same sequence data, so not independent evidence in species checks.
  final List<String> sourceAccessions;
  final String organism;

  /// Semicolon-separated lineage from the ORGANISM block.
  final String lineage;
  final String definition;
  final List<GbFeature> features;

  /// Lower-case nucleotide sequence from ORIGIN.
  final String sequence;

  /// `ncbi` or `file`.
  String source;
  final String? fileName;
  final int? fileSize;
  final DateTime? fileDate;

  // Editable in the "Edit record" dialog / filled from GBIF.
  String? editedOrganism;
  String? editedFamily;
  TaxonomyRanks ranks = const TaxonomyRanks();

  String get organismName => (editedOrganism ?? organism).trim();

  GbTaxonomy get taxonomy => GbTaxonomy.of(this);

  ({int pcgs, int rrnas, int trnas}) get featureCounts {
    var p = 0, r = 0, t = 0;
    for (final f in features) {
      switch (f.type) {
        case 'CDS':
          p++;
        case 'rRNA':
          r++;
        case 'tRNA':
          t++;
      }
    }
    return (pcgs: p, rrnas: r, trnas: t);
  }

  /// Total annotated length (bp) per marker class.
  ({int pcgs, int rrnas, int trnas}) get markerLengths {
    var p = 0, r = 0, t = 0;
    for (final f in features) {
      final len = locationLength(f.location);
      switch (f.type) {
        case 'CDS':
          p += len;
        case 'rRNA':
          r += len;
        case 'tRNA':
          t += len;
      }
    }
    return (pcgs: p, rrnas: r, trnas: t);
  }

  bool get hasFetchedTaxonomy => !ranks.isEmpty;

  /// Voucher of this genome as written in its title (DEFINITION), e.g.
  /// "Caranx hippos voucher MZUSP 1234 mitochondrion, complete genome" →
  /// "voucher MZUSP 1234"; null when the title names none.
  String? get voucher => titleVoucher(definition);

  /// NCBI Taxonomy id from the source feature (`/db_xref="taxon:8384"`).
  String? get taxId {
    for (final f in features) {
      if (f.type != 'source') continue;
      final m = RegExp(r'taxon:(\d+)').firstMatch(f.qualifiers['db_xref'] ?? '');
      if (m != null) return m.group(1);
    }
    return null;
  }

  /// Conservation category (from ICMBio SALVE).
  String? conservation;

  /// The species' SALVE sheet (Brazilian fauna only), when fetched.
  SalveSheet? salve;
}

/// What the ICMBio SALVE public API says about a species of the Brazilian
/// fauna: the national extinction-risk assessment and every section of its
/// sheet (taxonomy, distribution, natural history, population, threats,
/// uses, conservation, research, references).
class SalveSheet {
  const SalveSheet({
    required this.category,
    this.code = '',
    this.criteria = '',
    this.commonNames = '',
    this.group = '',
    this.biomes = '',
    this.states = '',
    this.watersheds = '',
    this.endemic = '',
    this.possiblyExtinct = false,
    this.assessed = '',
    this.published = '',
    this.situation = '',
    this.doi = '',
    this.synonyms = '',
    this.assessors = '',
    this.citation = '',
    this.justification = '',
    this.url = '',
    this.populationTrend = '',
    this.generationTime = '',
    this.sexRatio = '',
    this.migratory = '',
    this.habitats = const [],
    this.foodHabits = const [],
    this.reproduction = '',
    this.interactions = const [],
    this.threats = const [],
    this.uses = const [],
    this.officialList = '',
    this.officialLists = const [],
    this.history = const [],
    this.actions = const [],
    this.protectedAreas = const [],
    this.researchNeeds = const [],
    this.references = 0,
    this.texts = const {},
  });

  /// "Vulnerável (VU)".
  final String category;

  /// IUCN code: EX, EW, RE, CR, EN, VU, NT, LC, DD, NA, NE.
  final String code;
  final String criteria, commonNames, group, biomes, states, watersheds;

  /// Endemic to Brazil ("Sim" / "Não"), as SALVE writes it.
  final String endemic;
  final bool possiblyExtinct;

  /// Assessment period ("01/2023") and publication date of the sheet.
  final String assessed, published;

  /// Sheet status ("PUBLICADA"…).
  final String situation;
  final String doi, synonyms, assessors, citation, justification;

  /// The sheet on the SALVE site.
  final String url;

  /// Population: trend ("Declinando"…), generation time, sex ratio.
  final String populationTrend, generationTime, sexRatio;

  /// Natural history: migratory ("Sim"/"Não"), habitats, diet, the
  /// reproduction summary and the species it interacts with.
  final String migratory;
  final List<String> habitats, foodHabits;
  final String reproduction;
  final List<String> interactions;

  /// Threats (IUCN classification, most specific level) and uses.
  final List<String> threats, uses;

  /// In an official list ("Sim"/"Não"), which lists, the assessment
  /// history ("2023 · Nacional Brasil · VU A2cd"), conservation actions
  /// (PANs…), protected areas (UCs) where it occurs, research needs.
  final String officialList;
  final List<String> officialLists, history, actions, protectedAreas, researchNeeds;

  /// Bibliographic references of the sheet.
  final int references;

  /// Long texts of the sheet by section key (naturalHistory, habitat, food,
  /// reproduction, interactions, population, genetics, threats, uses,
  /// conservation, research, globalDistribution, nationalDistribution).
  final Map<String, String> texts;

  /// "Vulnerável (VU)" → "Vulnerável".
  String get categoryName => category.replaceAll(RegExp(r'\s*\([A-Z]{2}\)\s*$'), '').trim();
}

/// Ranks shown in the records table (`extractTaxonomy` in the HTML app).
class GbTaxonomy {
  const GbTaxonomy(this.kingdom, this.phylum, this.klass, this.order, this.family, this.genus, this.species, this.authorship);

  factory GbTaxonomy.of(GbRecord r) {
    final parts = r.organismName.split(RegExp(r'\s+'));
    final genus = parts.isEmpty ? '' : parts.first;
    final species = parts.length > 1 ? parts.skip(1).join(' ') : '';
    var family = r.ranks.family.isNotEmpty ? r.ranks.family : (r.editedFamily ?? '');
    if (family.isEmpty && r.lineage.isNotEmpty) {
      for (final part in r.lineage.split(RegExp(r';\s*'))) {
        final s = part.trim();
        if (RegExp(r'(idae|aceae)$', caseSensitive: false).hasMatch(s)) {
          family = s;
        } else if (family.isEmpty && RegExp(r'ales$', caseSensitive: false).hasMatch(s)) {
          family = s;
        }
      }
    }
    return GbTaxonomy(r.ranks.kingdom, r.ranks.phylum, r.ranks.klass, r.ranks.order, family, genus, species, r.ranks.authorship);
  }

  final String kingdom, phylum, klass, order, family, genus, species, authorship;

  String get binomial => species.isEmpty ? genus : '$genus $species';
}

/// Parses every record (`//`-terminated) of a GenBank file.
List<GbRecord> parseGenBankFile(String text, {String source = 'file', String? fileName, int? fileSize, DateTime? fileDate}) {
  final normalized = text.replaceAll('\r\n', '\n').replaceAll('\r', '\n');
  final out = <GbRecord>[];
  final chunk = StringBuffer();
  for (final line in normalized.split('\n')) {
    chunk.writeln(line);
    if (line.startsWith('//')) {
      final r = parseGenBank(chunk.toString(), source: source, fileName: fileName, fileSize: fileSize, fileDate: fileDate);
      if (r != null) out.add(r);
      chunk.clear();
    }
  }
  if (chunk.toString().trim().isNotEmpty) {
    final r = parseGenBank(chunk.toString(), source: source, fileName: fileName, fileSize: fileSize, fileDate: fileDate);
    if (r != null) out.add(r);
  }
  return out;
}

/// Parses one GenBank record; null when there is no LOCUS line.
GbRecord? parseGenBank(String text, {String source = 'file', String? fileName, int? fileSize, DateTime? fileDate}) {
  final lines = text.replaceAll('\r\n', '\n').replaceAll('\r', '\n').split('\n');
  String? accession;
  var organism = '';
  var definition = '';
  final lineage = <String>[];
  final features = <GbFeature>[];
  final seq = StringBuffer();
  final sources = <String>[];
  var hasLocus = false;

  var i = 0;
  while (i < lines.length) {
    final line = lines[i];
    if (line.startsWith('LOCUS')) {
      hasLocus = true;
      final m = RegExp(r'^LOCUS\s+(\S+)').firstMatch(line);
      accession ??= m?.group(1);
      i++;
    } else if (line.startsWith('DEFINITION')) {
      final b = StringBuffer(line.substring(10).trim());
      i++;
      while (i < lines.length && lines[i].startsWith('            ')) {
        b.write(' ${lines[i].trim()}');
        i++;
      }
      definition = b.toString();
    } else if (line.startsWith('VERSION')) {
      final m = RegExp(r'^VERSION\s+(\S+)').firstMatch(line);
      if (m != null) accession = m.group(1);
      i++;
    } else if (RegExp(r'^\s+ORGANISM\s+').hasMatch(line)) {
      organism = line.replaceFirst(RegExp(r'^\s+ORGANISM\s+'), '').trim();
      i++;
      while (i < lines.length && RegExp(r'^\s{10,}').hasMatch(lines[i]) && !lines[i].trimLeft().startsWith('/')) {
        final l = lines[i].trim();
        if (!(l.contains(';') || l.endsWith('.'))) break;
        lineage.add(l.replaceFirst(RegExp(r'\.$'), ''));
        i++;
      }
    } else if (line.startsWith('COMMENT')) {
      final b = StringBuffer(line.substring(7));
      i++;
      while (i < lines.length && lines[i].startsWith('            ')) {
        b.write(' ${lines[i].trim()}');
        i++;
      }
      final m = RegExp(r'(?:identical to|derived from)\s+([A-Z0-9_.,\s]+?)(?:\.\s|\.$|;|$)').firstMatch(b.toString());
      if (m != null) sources.addAll(RegExp(r'[A-Z]{1,2}_?\d{5,}(?:\.\d+)?').allMatches(m.group(1)!).map((x) => x.group(0)!));
    } else if (line.startsWith('FEATURES')) {
      i++;
      final block = <String>[];
      while (i < lines.length && (lines[i].startsWith(' ') || lines[i].isEmpty)) {
        block.add(lines[i]);
        i++;
      }
      features.addAll(_parseFeatures(block));
    } else if (line.startsWith('ORIGIN')) {
      i++;
      while (i < lines.length && !lines[i].startsWith('//')) {
        seq.write(lines[i].replaceAll(RegExp(r'[\s\d]'), ''));
        i++;
      }
    } else {
      i++;
    }
  }
  if (!hasLocus || accession == null) return null;
  return GbRecord(
    accession: accession,
    organism: organism,
    lineage: lineage.join(' '),
    definition: definition,
    features: features,
    sequence: seq.toString().toLowerCase(),
    source: source,
    fileName: fileName,
    fileSize: fileSize,
    fileDate: fileDate,
    sourceAccessions: sources,
  );
}

final _labelSuffix = RegExp(
  r'(\s*-\s*(truncated|partial|pseudo(gene)?)\s*)$|(\s+(cds|gene|trna|rrna|mrna|exon|intron|orf|ncrna|misc_rna|misc_feature)\s*)$',
  caseSensitive: false,
);
final _orfLabel = RegExp(r'^orf\s+\d+\s*\(frame', caseSensitive: false);

/// Gene name in a /label qualifier, written by editors such as Geneious when
/// a feature is curated by hand: `ycf1 CDS`, `trnC-GCA tRNA`,
/// `ycf2 gene - Truncated` and `rps12` all name the gene alone. Automatic
/// ORF-finder labels (`ORF 4 (frame 2)`) are not genes.
String? geneFromLabel(String? label) {
  var v = (label ?? '').replaceAll(RegExp(r'\s+'), ' ').trim();
  if (v.isEmpty || _orfLabel.hasMatch(v)) return null;
  for (var i = 0; i < 3; i++) {
    final next = v.replaceFirst(_labelSuffix, '').trim();
    if (next == v) break;
    v = next;
  }
  return v.isEmpty ? null : v;
}

final _featureStart = RegExp(r'^ {5}(\S+)\s+(\S.*)$');
final _qualifier = RegExp(r'^ {21}/([^=\s]+)(?:=(.*))?$');

List<GbFeature> _parseFeatures(List<String> lines) {
  final out = <GbFeature>[];
  String? type;
  final loc = StringBuffer();
  var quals = <String, String>{};
  String? lastKey;
  var inLocation = false;

  void flush() {
    final t = type;
    if (t != null) {
      // Curated records (Geneious and similar editors) often name a feature
      // only with /label="ycf1 CDS": treat it as /gene when /gene is absent.
      if ((quals['gene'] ?? '').trim().isEmpty) {
        final g = geneFromLabel(quals['label']);
        if (g != null) quals['gene'] = g;
      }
      out.add(GbFeature(t, loc.toString().replaceAll(' ', ''), quals));
    }
    type = null;
    loc.clear();
    quals = {};
    lastKey = null;
  }

  for (final line in lines) {
    if (line.trim().isEmpty) continue;
    final f = _featureStart.firstMatch(line);
    if (f != null) {
      flush();
      type = f.group(1);
      loc.write(f.group(2)!.trim());
      inLocation = true;
      continue;
    }
    if (type == null) continue;
    final q = _qualifier.firstMatch(line);
    if (q != null) {
      inLocation = false;
      lastKey = q.group(1);
      var v = q.group(2) ?? '';
      if (v.startsWith('"')) v = v.substring(1);
      if (v.endsWith('"')) v = v.substring(0, v.length - 1);
      // /db_xref repeats (taxon:, BOLD:…): keep every value.
      final previous = quals[lastKey!];
      quals[lastKey!] = lastKey == 'db_xref' && previous != null ? '$previous|$v' : v;
    } else if (inLocation) {
      loc.write(line.trim());
    } else if (lastKey != null) {
      var v = line.trim();
      if (v.endsWith('"')) v = v.substring(0, v.length - 1);
      // Translations and sequences continue without spaces.
      final sep = lastKey == 'translation' ? '' : ' ';
      quals[lastKey!] = '${quals[lastKey]}$sep$v';
    }
  }
  flush();
  return out;
}

// ---------------------------------------------------------------------------
// Locations
// ---------------------------------------------------------------------------

/// Splits `a,b(c,d),e` on top-level commas.
List<String> _splitTop(String s) {
  final parts = <String>[];
  var depth = 0, start = 0;
  for (var i = 0; i < s.length; i++) {
    final ch = s[i];
    if (ch == '(') depth++;
    if (ch == ')') depth--;
    if (ch == ',' && depth == 0) {
      parts.add(s.substring(start, i));
      start = i + 1;
    }
  }
  parts.add(s.substring(start));
  return parts;
}

String _inner(String loc, String fn) => loc.substring(fn.length + 1, loc.lastIndexOf(')'));

/// Extracts the feature sequence (upper case) following the location.
String extractLocation(String sequence, String location) => _extract(sequence, location.trim()).toUpperCase();

String _extract(String seq, String loc) {
  if (loc.startsWith('complement(')) return reverseComplement(_extract(seq, _inner(loc, 'complement')));
  if (loc.startsWith('join(')) return _splitTop(_inner(loc, 'join')).map((p) => _extract(seq, p.trim())).join();
  if (loc.startsWith('order(')) return _splitTop(_inner(loc, 'order')).map((p) => _extract(seq, p.trim())).join();
  if (loc.contains(':')) return ''; // remote entry (other accession)
  final range = RegExp(r'^<?(\d+)\.\.>?(\d+)$').firstMatch(loc);
  if (range != null) {
    final start = int.parse(range.group(1)!) - 1;
    final end = int.parse(range.group(2)!);
    if (start < 0 || end > seq.length || start >= end) {
      // Circular genomes: a range crossing the origin is written start..end with end < start.
      if (start >= 0 && end < start && start < seq.length && end <= seq.length) return seq.substring(start) + seq.substring(0, end);
      return '';
    }
    return seq.substring(start, end);
  }
  final single = RegExp(r'^[<>]?(\d+)$').firstMatch(loc);
  if (single != null) {
    final pos = int.parse(single.group(1)!) - 1;
    return pos >= 0 && pos < seq.length ? seq[pos] : '';
  }
  return ''; // between-bases (12^13) or unsupported
}

/// Number of bases covered by a location.
int locationLength(String location) => _length(location.trim());

int _length(String loc) {
  if (loc.startsWith('complement(')) return _length(_inner(loc, 'complement'));
  if (loc.startsWith('join(')) return _splitTop(_inner(loc, 'join')).fold(0, (a, p) => a + _length(p.trim()));
  if (loc.startsWith('order(')) return _splitTop(_inner(loc, 'order')).fold(0, (a, p) => a + _length(p.trim()));
  if (loc.contains(':')) return 0;
  final range = RegExp(r'^<?(\d+)\.\.>?(\d+)$').firstMatch(loc);
  if (range != null) {
    final len = int.parse(range.group(2)!) - int.parse(range.group(1)!) + 1;
    return len > 0 ? len : 0;
  }
  return RegExp(r'^[<>]?\d+$').hasMatch(loc) ? 1 : 0;
}

const _complement = {
  'A': 'T', 'T': 'A', 'U': 'A', 'G': 'C', 'C': 'G', 'N': 'N', //
  'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K',
  'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D',
  'a': 't', 't': 'a', 'u': 'a', 'g': 'c', 'c': 'g', 'n': 'n',
  'r': 'y', 'y': 'r', 's': 's', 'w': 'w', 'k': 'm', 'm': 'k',
  'b': 'v', 'v': 'b', 'd': 'h', 'h': 'd',
};

String reverseComplement(String seq) {
  final out = StringBuffer();
  for (var i = seq.length - 1; i >= 0; i--) {
    out.write(_complement[seq[i]] ?? seq[i]);
  }
  return out.toString();
}

/// "voucher X", "isolate X", "strain X" or "specimen X" from a GenBank
/// title, stopping before the organelle / completeness words.
String? titleVoucher(String title) {
  final m = RegExp(
    r'\b((?:specimen\s+)?voucher|isolate|strain|specimen)\s+(.+?)(?=\s+(?:mitochondri\w*|chloroplast|plastid|complete|partial|genomic|DNA|genome)\b|[,;]|$)',
    caseSensitive: false,
  ).firstMatch(title);
  if (m == null) return null;
  final v = m.group(2)!.trim();
  return v.isEmpty ? null : '${m.group(1)!.toLowerCase()} $v';
}
