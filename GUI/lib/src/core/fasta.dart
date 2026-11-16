/// FASTA parsing and writing.
///
/// Pure Dart (no `dart:io`, no Flutter) so it can run inside isolates and be
/// unit-tested in isolation.
library;

class SeqRecord {
  SeqRecord(this.header, this.sequence);

  /// Full header line without the leading `>`.
  final String header;

  /// Sequence with whitespace removed. Case is preserved.
  final String sequence;

  /// First whitespace-delimited token of the header.
  String get id {
    final space = header.indexOf(RegExp(r'\s'));
    return space < 0 ? header : header.substring(0, space);
  }

  int get length => sequence.length;

  SeqRecord copyWith({String? header, String? sequence}) =>
      SeqRecord(header ?? this.header, sequence ?? this.sequence);
}

/// Parses FASTA text. Blank lines and `;` comment lines are ignored.
/// Text before the first header is ignored as well.
List<SeqRecord> parseFasta(String text) {
  final records = <SeqRecord>[];
  String? header;
  final buffer = StringBuffer();

  void flush() {
    if (header != null) records.add(SeqRecord(header, buffer.toString()));
    buffer.clear();
  }

  for (final rawLine in text.split('\n')) {
    final line = rawLine.endsWith('\r')
        ? rawLine.substring(0, rawLine.length - 1)
        : rawLine;
    if (line.isEmpty || line.startsWith(';')) continue;
    if (line.startsWith('>')) {
      flush();
      header = line.substring(1).trim();
    } else if (header != null) {
      buffer.write(line.replaceAll(RegExp(r'\s'), ''));
    }
  }
  flush();
  return records;
}

String writeFasta(Iterable<SeqRecord> records, {int lineWidth = 0}) {
  final out = StringBuffer();
  for (final r in records) {
    out.writeln('>${r.header}');
    if (lineWidth <= 0) {
      out.writeln(r.sequence);
    } else {
      for (var i = 0; i < r.sequence.length; i += lineWidth) {
        final end = i + lineWidth > r.sequence.length
            ? r.sequence.length
            : i + lineWidth;
        out.writeln(r.sequence.substring(i, end));
      }
    }
  }
  return out.toString();
}

/// Makes a taxon name safe for PHYLIP/NEXUS/Newick and for RAxML-NG, which
/// rejects names containing spaces, colons, commas, parentheses, brackets
/// or semicolons.
String sanitizeTaxonName(String name) {
  final cleaned = name
      .trim()
      .replaceAll(RegExp(r"[\s:,;()\[\]'" '"' r'=]+'), '_')
      .replaceAll(RegExp(r'_+'), '_')
      .replaceAll(RegExp(r'^_|_$'), '');
  return cleaned.isEmpty ? 'unnamed' : cleaned;
}

/// True when all sequences have the same length (an alignment).
bool isAligned(Iterable<SeqRecord> records) {
  int? len;
  for (final r in records) {
    len ??= r.length;
    if (r.length != len) return false;
  }
  return len != null;
}

/// Short stable fingerprint of [records] (headers and sequences, FNV-1a),
/// used to tell whether a marker changed since its last alignment.
String fingerprint(List<SeqRecord> records) {
  var h = 0x811c9dc5;
  void add(String s) {
    for (var i = 0; i < s.length; i++) {
      h = ((h ^ s.codeUnitAt(i)) * 0x01000193) & 0xFFFFFFFF;
    }
    h = ((h ^ 0x0a) * 0x01000193) & 0xFFFFFFFF;
  }

  for (final r in records) {
    add(r.header);
    add(r.sequence);
  }
  return '${records.length}:${h.toRadixString(16)}';
}
