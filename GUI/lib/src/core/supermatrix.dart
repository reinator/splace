/// Concatenation of per-locus alignments into a supermatrix, with optional
/// codon-position exclusion and partition files for RAxML-NG and IQ-TREE.
library;

import 'dart:typed_data';

import 'fasta.dart';
import 'saturation.dart';

class LocusAlignment {
  LocusAlignment(this.name, this.records, {this.coding = true, this.frame = 0});

  final String name;
  final List<SeqRecord> records;

  /// Protein-coding loci get codon partitions and position exclusion.
  final bool coding;

  /// Column (0..2) where the first complete codon starts.
  final int frame;

  int get length => records.isEmpty ? 0 : records.first.length;
}

final _nonCodingName = RegExp(r'^(rrn|trn|rRNA|tRNA|\d+S|ITS|5\.8S|18S|28S|12S|16S|D-?loop|CR)', caseSensitive: false);

/// Guesses whether an alignment is protein-coding: rRNA/tRNA-like names are
/// not; otherwise a locus is coding when some frame has (almost) no internal
/// stop codons. Returns the coding flag and best frame.
({bool coding, int frame}) guessCoding(String name, List<SeqRecord> alignment, int geneticCode) {
  if (_nonCodingName.hasMatch(name)) return (coding: false, frame: 0);
  final stops = [for (var f = 0; f < 3; f++) countInternalStops(alignment, f, geneticCode)];
  var best = 0;
  for (var f = 1; f < 3; f++) {
    if (stops[f] < stops[best]) best = f;
  }
  // Tolerate a few stops (sequencing errors, pseudogenes) relative to size.
  final tolerated = (alignment.length * 0.05).ceil();
  return (coding: stops[best] <= tolerated, frame: best);
}

/// Names of protein-coding genes of organelle genomes (SynGenes style):
/// they are coding even when stop codons point the other way (wrong code,
/// pseudogene copies, sequencing errors).
final _codingName = RegExp(
    r'^(CO[I123]+|COX\d|CYTB|COB|ND\d|ND4L|NAD\d|NAD4L|ATP\d|ATP[689]|rbcL|matK|ndh[A-K]|psa[A-J]|psb[A-Z]|pet[A-N]|rp[ls]\d+|rpo[ABC]\d?|ycf\d+|ccsA|cemA|clpP|accD|inf[A1])\b',
    caseSensitive: false);

/// A known protein-coding gene name (cox1, rbcL, ndhF…).
bool knownCodingName(String name) => _codingName.hasMatch(name);

/// Coding detection of a DNA alignment with the project's genetic code
/// only (no guessing of another code, which could mislead): rRNA/tRNA-like
/// names are not coding; known protein-coding names always are; other
/// loci are coding when some frame has (almost) no internal stops. The
/// frame is the one with the fewest stops.
({bool coding, int frame, int stops}) detectCoding(String name, List<SeqRecord> alignment, int code) {
  if (_nonCodingName.hasMatch(name)) return (coding: false, frame: 0, stops: 0);
  final tolerated = (alignment.length * 0.05).ceil();
  final stops = [for (var f = 0; f < 3; f++) countInternalStops(alignment, f, code)];
  var best = 0;
  for (var f = 1; f < 3; f++) {
    if (stops[f] < stops[best]) best = f;
  }
  final coding = stops[best] <= tolerated || _codingName.hasMatch(name);
  return (coding: coding, frame: best, stops: stops[best]);
}

enum MissingDataPolicy {
  /// Keep all loci; taxa missing a locus get `?` for that block.
  fill,

  /// Keep only loci present in every taxon.
  completeLociOnly,
}

enum PartitionScheme { none, byLocus, byLocusAndCodon, byCodonPosition }

enum SeqDataType { dna, protein }

class ConcatOptions {
  const ConcatOptions({
    this.missingPolicy = MissingDataPolicy.fill,
    this.excludePositions = const {},
    this.perLocusExclude = const {},
    this.scheme = PartitionScheme.byLocusAndCodon,
    this.minLociPerTaxon = 1,
    this.replaceGapsWithMissing = false,
    this.outgroups = const [],
    this.dataType = SeqDataType.dna,
  });

  final MissingDataPolicy missingPolicy;

  /// Codon positions (1..3) removed from every coding locus.
  final Set<int> excludePositions;

  /// Locus-specific exclusions (e.g. from per-locus saturation results);
  /// merged with [excludePositions].
  final Map<String, Set<int>> perLocusExclude;
  final PartitionScheme scheme;

  /// Taxa with fewer loci than this are dropped.
  final int minLociPerTaxon;

  /// Writes alignment gaps as `?` (missing) in the exported matrices.
  final bool replaceGapsWithMissing;

  /// Written to the NEXUS `assumptions` block.
  final List<String> outgroups;
  final SeqDataType dataType;

  ConcatOptions copyWith({
    MissingDataPolicy? missingPolicy,
    Set<int>? excludePositions,
    Map<String, Set<int>>? perLocusExclude,
    PartitionScheme? scheme,
    int? minLociPerTaxon,
    bool? replaceGapsWithMissing,
    List<String>? outgroups,
    SeqDataType? dataType,
  }) =>
      ConcatOptions(
        missingPolicy: missingPolicy ?? this.missingPolicy,
        excludePositions: excludePositions ?? this.excludePositions,
        perLocusExclude: perLocusExclude ?? this.perLocusExclude,
        scheme: scheme ?? this.scheme,
        minLociPerTaxon: minLociPerTaxon ?? this.minLociPerTaxon,
        replaceGapsWithMissing: replaceGapsWithMissing ?? this.replaceGapsWithMissing,
        outgroups: outgroups ?? this.outgroups,
        dataType: dataType ?? this.dataType,
      );
}

/// 1-based inclusive range with optional stride (`start-end\stride`).
class CharRange {
  const CharRange(this.start, this.end, [this.stride = 1]);
  final int start;
  final int end;
  final int stride;

  String format() => stride == 1 ? '$start-$end' : '$start-$end\\$stride';
}

class Partition {
  Partition(this.name, this.ranges, {this.coding = false});
  final String name;
  final List<CharRange> ranges;
  final bool coding;

  String get rangeString => ranges.map((r) => r.format()).join(', ');
}

class Supermatrix {
  Supermatrix({
    required this.taxa,
    required this.sequences,
    required this.partitions,
    required this.locusRanges,
    required this.droppedLoci,
    required this.droppedTaxa,
    this.dataType = SeqDataType.dna,
    this.outgroups = const [],
  });

  final List<String> taxa;
  final Map<String, String> sequences;
  final List<Partition> partitions;

  /// locus -> 1-based inclusive range in the matrix (after exclusions).
  final Map<String, CharRange> locusRanges;
  final List<String> droppedLoci;
  final List<String> droppedTaxa;
  final SeqDataType dataType;
  final List<String> outgroups;

  bool get isProtein => dataType == SeqDataType.protein;

  int get length => sequences.isEmpty ? 0 : sequences.values.first.length;

  /// Gaps and unknown bases over all cells. Computed once: the screens
  /// read it on every rebuild and the matrix can hold millions of cells.
  late final double missingFraction = () {
    if (sequences.isEmpty || length == 0) return 0.0;
    var missing = 0;
    for (final s in sequences.values) {
      for (var i = 0; i < s.length; i++) {
        switch (s.codeUnitAt(i)) {
          case 0x2D || 0x3F || 0x4E || 0x6E || 0x58 || 0x78: // - ? N n X x
            missing++;
        }
      }
    }
    return missing / (sequences.length * length);
  }();

  String toFasta() => writeFasta(taxa.map((t) => SeqRecord(t, sequences[t]!)));

  /// Relaxed PHYLIP (names may exceed 10 characters), sequential.
  String toPhylip() {
    final b = StringBuffer('${taxa.length} $length\n');
    final width = taxa.fold<int>(0, (a, t) => t.length > a ? t.length : a) + 2;
    for (final t in taxa) {
      b.writeln('${t.padRight(width)}${sequences[t]}');
    }
    return b.toString();
  }

  String toNexus({bool includeSets = true}) {
    final b = StringBuffer('#NEXUS\n\nBEGIN DATA;\n')
      ..writeln('  DIMENSIONS NTAX=${taxa.length} NCHAR=$length;')
      ..writeln('  FORMAT DATATYPE=${isProtein ? 'PROTEIN' : 'DNA'} GAP=- MISSING=?;')
      ..writeln('  MATRIX');
    final width = taxa.fold<int>(0, (a, t) => t.length > a ? t.length : a) + 4;
    for (final t in taxa) {
      b.writeln('    ${"'$t'".padRight(width)}${sequences[t]}');
    }
    b.writeln('  ;\nEND;');
    if (includeSets && partitions.isNotEmpty) {
      b.writeln('\nBEGIN SETS;');
      for (final p in partitions) {
        b.writeln('  CHARSET ${p.name} = ${p.ranges.map((r) => r.format()).join(' ')};');
      }
      b.writeln('END;');
    }
    final og = outgroups.where(taxa.contains).toList();
    if (og.isNotEmpty) {
      b.writeln('\nBEGIN ASSUMPTIONS;\n  OUTGROUP ${og.map((t) => "'$t'").join(' ')};\nEND;');
    }
    return b.toString();
  }

  /// Partition table in the RAxML/IQ-TREE `DNA, name = range` style (the
  /// "Save partition table" file of the HTML app).
  String toPartitionTable() =>
      partitions.map((p) => '${isProtein ? 'AA' : 'DNA'}, ${p.name} = ${p.rangeString}').join('\n');

  /// RAxML-NG partition file. [models] maps partition name to a RAxML-NG
  /// model string; missing entries use [defaultModel].
  String toRaxmlPartitions({Map<String, String> models = const {}, String? defaultModel}) {
    final b = StringBuffer();
    final fallback = defaultModel ?? (isProtein ? 'LG+G' : 'GTR+G');
    for (final p in partitions) {
      b.writeln('${models[p.name] ?? fallback}, ${p.name} = ${p.rangeString}');
    }
    return b.toString();
  }

  /// NEXUS partition file for IQ-TREE (`-p`/`-q`/`-Q`). With [models] a
  /// charpartition fixes the model of each partition; otherwise IQ-TREE's
  /// ModelFinder chooses them.
  String toIqtreePartitions({Map<String, String> models = const {}}) {
    final b = StringBuffer('#nexus\nbegin sets;\n');
    for (final p in partitions) {
      b.writeln('  charset ${p.name} = ${p.ranges.map((r) => r.format()).join(' ')};');
    }
    if (models.isNotEmpty) {
      final parts = partitions.map((p) => '${models[p.name] ?? (isProtein ? 'LG+G4' : 'GTR+F+G4')}: ${p.name}').join(', ');
      b.writeln('  charpartition splace = $parts;');
    }
    b.writeln('end;');
    return b.toString();
  }
}

/// Builds the supermatrix. Taxon names must already match across loci
/// (see the locus organizer).
Supermatrix buildSupermatrix(List<LocusAlignment> loci, [ConcatOptions options = const ConcatOptions()]) {
  final allTaxa = <String>{for (final l in loci) ...l.records.map((r) => r.header)};
  final dropped = <String>[];
  final kept = <LocusAlignment>[];
  for (final l in loci) {
    if (l.records.isEmpty || !isAligned(l.records)) {
      dropped.add(l.name);
      continue;
    }
    if (options.missingPolicy == MissingDataPolicy.completeLociOnly &&
        l.records.map((r) => r.header).toSet().length < allTaxa.length) {
      dropped.add(l.name);
      continue;
    }
    kept.add(l);
  }

  // Drop poorly sampled taxa.
  final lociPerTaxon = <String, int>{};
  for (final l in kept) {
    for (final r in l.records) {
      lociPerTaxon.update(r.header, (v) => v + 1, ifAbsent: () => 1);
    }
  }
  final droppedTaxa = [
    for (final t in allTaxa)
      if ((lociPerTaxon[t] ?? 0) < options.minLociPerTaxon) t,
  ]..sort();
  final taxa = (lociPerTaxon.keys.where((t) => !droppedTaxa.contains(t)).toList())..sort();

  final buffers = {for (final t in taxa) t: StringBuffer()};
  final partitions = <Partition>[];
  final locusRanges = <String, CharRange>{};
  final globalCodon = <int, List<CharRange>>{1: [], 2: [], 3: []};
  var offset = 0;

  final protein = options.dataType == SeqDataType.protein;
  for (final l in kept) {
    final coding = l.coding && !protein;
    final excluded = coding
        ? {...options.excludePositions, ...?options.perLocusExclude[l.name]}
        : <int>{};
    // Output columns and their codon position (0 = non-coding).
    final columns = <int>[];
    final positions = <int>[];
    for (var c = 0; c < l.length; c++) {
      final pos = coding ? codonPositionOf(c, l.frame) : 0;
      if (excluded.contains(pos)) continue;
      columns.add(c);
      positions.add(pos);
    }
    if (columns.isEmpty) {
      dropped.add(l.name);
      continue;
    }
    final byTaxon = {for (final r in l.records) r.header: r.sequence};
    for (final t in taxa) {
      final seq = byTaxon[t];
      final buf = buffers[t]!;
      if (seq == null) {
        buf.write('?' * columns.length);
      } else if (columns.length == seq.length) {
        buf.write(seq);
      } else {
        for (final c in columns) {
          buf.writeCharCode(seq.codeUnitAt(c));
        }
      }
    }

    final name = sanitizeTaxonName(l.name);
    locusRanges[l.name] = CharRange(offset + 1, offset + columns.length);
    final usesCodons = coding && (options.scheme == PartitionScheme.byLocusAndCodon || options.scheme == PartitionScheme.byCodonPosition);
    if (usesCodons) {
      for (final pos in {1, 2, 3}.difference(excluded)) {
        final ranges = _rangesFor(positions, pos, offset);
        if (ranges.isEmpty) continue;
        if (options.scheme == PartitionScheme.byLocusAndCodon) {
          partitions.add(Partition('${name}_pos$pos', ranges, coding: true));
        } else {
          globalCodon[pos]!.addAll(ranges);
        }
      }
    } else if (options.scheme != PartitionScheme.none) {
      partitions.add(Partition(name, [locusRanges[l.name]!]));
    }
    offset += columns.length;
  }

  if (options.scheme == PartitionScheme.byCodonPosition) {
    partitions.insertAll(0, [
      for (final pos in [1, 2, 3])
        if (globalCodon[pos]!.isNotEmpty) Partition('codon_pos$pos', globalCodon[pos]!, coding: true),
    ]);
  }

  String finish(String s) => options.replaceGapsWithMissing ? s.replaceAll('-', '?') : s;
  return Supermatrix(
    taxa: taxa,
    sequences: {for (final t in taxa) t: finish(buffers[t].toString())},
    dataType: options.dataType,
    outgroups: options.outgroups,
    partitions: partitions,
    locusRanges: locusRanges,
    droppedLoci: dropped,
    droppedTaxa: droppedTaxa,
  );
}

/// Ranges (1-based, matrix coordinates) of the columns with codon position
/// [pos]. Columns of one position are periodic, so a single strided range
/// is enough; non-periodic input falls back to explicit runs.
List<CharRange> _rangesFor(List<int> positions, int pos, int offset) {
  final idx = [for (var i = 0; i < positions.length; i++) if (positions[i] == pos) i];
  if (idx.isEmpty) return const [];
  if (idx.length == 1) return [CharRange(offset + idx.first + 1, offset + idx.first + 1)];
  final stride = idx[1] - idx[0];
  var periodic = true;
  for (var i = 2; i < idx.length; i++) {
    if (idx[i] - idx[i - 1] != stride) {
      periodic = false;
      break;
    }
  }
  if (periodic) return [CharRange(offset + idx.first + 1, offset + idx.last + 1, stride)];
  final ranges = <CharRange>[];
  var start = idx.first, prev = idx.first;
  for (final i in idx.skip(1)) {
    if (i != prev + 1) {
      ranges.add(CharRange(offset + start + 1, offset + prev + 1));
      start = i;
    }
    prev = i;
  }
  ranges.add(CharRange(offset + start + 1, offset + prev + 1));
  return ranges;
}

/// What each taxon has in each locus of a supermatrix (like the main table
/// of SequenceMatrix): bases and gaps per cell, a missing locus has no base.
class MatrixOverview {
  MatrixOverview._(this.taxa, this.loci, this.lengths, this.bases, this.gaps);

  final List<String> taxa;

  /// Loci in matrix order, with their length in the matrix.
  final List<String> loci;
  final List<int> lengths;

  /// [taxon][locus]: characters that are not `-` or `?`, and gaps (`-`).
  final List<Int32List> bases;
  final List<Int32List> gaps;

  bool present(int taxon, int locus) => bases[taxon][locus] > 0;

  int taxaWith(int locus) => [for (final b in bases) if (b[locus] > 0) 1].length;

  int lociOf(int taxon) => bases[taxon].where((b) => b > 0).length;

  int basesOf(int taxon) => bases[taxon].fold(0, (a, b) => a + b);

  int get length => lengths.fold(0, (a, b) => a + b);

  /// Gaps among the taxa that have the locus.
  double gapFraction(int locus) {
    var g = 0, n = 0;
    for (var t = 0; t < taxa.length; t++) {
      if (bases[t][locus] == 0) continue;
      g += gaps[t][locus];
      n += lengths[locus];
    }
    return n == 0 ? 0 : g / n;
  }

  /// One row per taxon × locus, for spreadsheets.
  String toTsv() {
    final b = StringBuffer('taxon\tlocus\tlocus_length\tbases\tgaps\tstatus\n');
    for (var t = 0; t < taxa.length; t++) {
      for (var l = 0; l < loci.length; l++) {
        b.writeln('${taxa[t]}\t${loci[l]}\t${lengths[l]}\t${bases[t][l]}\t${gaps[t][l]}\t${present(t, l) ? 'present' : 'missing'}');
      }
    }
    return b.toString();
  }
}

MatrixOverview matrixOverview(Supermatrix m) {
  final ranges = m.locusRanges.entries.toList()..sort((a, b) => a.value.start.compareTo(b.value.start));
  final bases = <Int32List>[], gaps = <Int32List>[];
  for (final taxon in m.taxa) {
    final s = m.sequences[taxon] ?? '';
    final b = Int32List(ranges.length), g = Int32List(ranges.length);
    for (final (i, e) in ranges.indexed) {
      for (var c = e.value.start - 1; c < e.value.end && c < s.length; c++) {
        final u = s.codeUnitAt(c);
        if (u == 0x2D) {
          g[i]++;
        } else if (u != 0x3F) {
          b[i]++;
        }
      }
    }
    bases.add(b);
    gaps.add(g);
  }
  return MatrixOverview._(
    m.taxa,
    [for (final e in ranges) e.key],
    [for (final e in ranges) e.value.end - e.value.start + 1],
    bases,
    gaps,
  );
}
