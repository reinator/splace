/// Writes a GenBank record back as a flat file, corrected: the organism and
/// lineage from the species identification and the taxonomy (GBIF, NCBI,
/// SALVE), and the gene names standardized by SynGenes. The original header
/// (references, comments) is kept; every feature and qualifier is written
/// again, with the original name of a renamed gene in a /note.
library;

import 'genbank.dart';

/// Qualifiers written without quotes (INSDC feature table).
const _unquoted = {
  'codon_start', 'transl_table', 'number', 'anticodon', 'transl_except', 'rpt_type', 'rpt_unit_range',
  'direction', 'mod_base', 'estimated_length', 'citation', 'compare', 'tag_peptide', 'rpt_unit_seq',
};

/// Lineage line of a record: the ranks found (kingdom to family) and the
/// genus, as GenBank writes it ("Animalia; Chordata; …; Caranx."); the
/// original lineage when no rank is known.
String correctedLineage(GbRecord r) {
  final t = r.taxonomy;
  final ranks = [t.kingdom, t.phylum, t.klass, t.order, t.family, t.genus].where((x) => x.trim().isNotEmpty).toList();
  if (r.ranks.isEmpty || ranks.length < 2) return r.lineage.replaceFirst(RegExp(r'\.$'), '');
  return ranks.join('; ');
}

/// [r] as a GenBank flat file. [geneNames] gives the standard name of
/// features by index (SynGenes, as in the gene catalog); [note] is added to
/// the COMMENT (what was corrected and by which program).
String writeGenBank(GbRecord r, {Map<int, String> geneNames = const {}, String? note}) {
  final organism = r.organismName;
  final lineage = correctedLineage(r);
  final b = StringBuffer();
  b.write(_header(r, organism, lineage, note));

  // Original → standard name of each renamed gene, so the `gene` features
  // (and any other feature naming it) follow their CDS, rRNA or tRNA.
  final renamed = <String, String>{};
  for (final e in geneNames.entries) {
    final old = r.features[e.key].q('gene');
    if (old != null && old != e.value) renamed[old] = e.value;
  }

  b.writeln('FEATURES             Location/Qualifiers');
  for (var fi = 0; fi < r.features.length; fi++) {
    final f = r.features[fi];
    final quals = Map<String, String>.of(f.qualifiers);
    if (f.type == 'source') {
      quals['organism'] = organism;
    }
    final old = quals['gene']?.replaceAll(RegExp(r'\s+'), ' ').trim();
    final std = geneNames[fi] ?? (old == null ? null : renamed[old]);
    if (std != null && std != old) {
      quals['gene'] = std;
      final was = old == null || old.isEmpty ? '' : '; original annotation: $old';
      final n = 'gene name standardized by SynGenes$was';
      quals['note'] = (quals['note'] ?? '').isEmpty ? n : '${quals['note']}; $n';
    }
    b.write(_feature(f.type, f.location));
    for (final e in quals.entries) {
      // /db_xref values were joined with | when read.
      final values = e.key == 'db_xref' ? e.value.split('|') : [e.value];
      for (final v in values) {
        b.write(_qualifier(e.key, v));
      }
    }
  }
  b.writeln('ORIGIN      ');
  final seq = r.sequence.toLowerCase();
  for (var i = 0; i < seq.length; i += 60) {
    final line = seq.substring(i, (i + 60).clamp(0, seq.length));
    final groups = [for (var k = 0; k < line.length; k += 10) line.substring(k, (k + 10).clamp(0, line.length))];
    b.writeln('${'${i + 1}'.padLeft(9)} ${groups.join(' ')}');
  }
  b.writeln('//');
  return b.toString();
}

/// The original header with the new SOURCE/ORGANISM block and the old
/// species name replaced in the DEFINITION; a minimal header when the
/// record has none.
String _header(GbRecord r, String organism, String lineage, String? note) {
  final orgBlock = [
    'SOURCE      $organism',
    '  ORGANISM  $organism',
    ..._wrap('$lineage.', 68).map((l) => '            $l'),
  ];
  final commentLines = note == null ? <String>[] : _wrap(note, 68);
  if (r.header.trim().isEmpty) {
    return [
      'LOCUS       ${r.accession.padRight(16)} ${r.sequence.length} bp    DNA     linear',
      'DEFINITION  ${r.definition.isEmpty ? organism : r.definition.replaceAll(r.organism, organism)}',
      'ACCESSION   ${r.accession.split('.').first}',
      'VERSION     ${r.accession}',
      ...orgBlock,
      if (commentLines.isNotEmpty) ...['COMMENT     ${commentLines.first}', ...commentLines.skip(1).map((l) => '            $l')],
      '',
    ].join('\n');
  }
  final lines = r.header.replaceAll('\r\n', '\n').split('\n');
  if (lines.isNotEmpty && lines.last.isEmpty) lines.removeLast();
  final out = <String>[];
  var replacedSource = false, inDefinition = false, addedComment = false;
  for (var i = 0; i < lines.length; i++) {
    final l = lines[i];
    if (l.startsWith('SOURCE')) {
      // Skip the old SOURCE block: its continuation and the ORGANISM lines.
      var j = i + 1;
      while (j < lines.length && !RegExp(r'^[A-Z]').hasMatch(lines[j])) {
        j++;
      }
      out.addAll(orgBlock);
      replacedSource = true;
      i = j - 1;
      continue;
    }
    if (l.startsWith('DEFINITION')) inDefinition = true;
    if (RegExp(r'^[A-Z]').hasMatch(l) && !l.startsWith('DEFINITION')) inDefinition = false;
    if (l.startsWith('COMMENT') && commentLines.isNotEmpty && !addedComment) {
      // Ours first, then the original comment on the following lines.
      out.add('COMMENT     ${commentLines.first}');
      out.addAll(commentLines.skip(1).map((c) => '            $c'));
      out.add('            ${l.substring(7).trim()}');
      addedComment = true;
      continue;
    }
    out.add(inDefinition && r.organism.isNotEmpty && organism != r.organism ? l.replaceAll(r.organism, organism) : l);
  }
  if (!replacedSource) {
    final at = out.indexWhere((l) => l.startsWith('REFERENCE') || l.startsWith('COMMENT'));
    out.insertAll(at < 0 ? out.length : at, orgBlock);
  }
  if (commentLines.isNotEmpty && !addedComment) {
    out.add('COMMENT     ${commentLines.first}');
    out.addAll(commentLines.skip(1).map((c) => '            $c'));
  }
  return '${out.join('\n')}\n';
}

String _feature(String type, String location) {
  final lines = _wrapAt(location, 58, ',');
  return [
    '     ${type.padRight(16)}${lines.first}',
    ...lines.skip(1).map((l) => '${' ' * 21}$l'),
  ].map((l) => '$l\n').join();
}

String _qualifier(String key, String value) {
  final v = value.replaceAll(RegExp(r'\s+'), key == 'translation' ? '' : ' ').trim();
  final text = value.isEmpty && key != 'note'
      ? '/$key'
      : _unquoted.contains(key) && !v.contains(' ')
          ? '/$key=$v'
          : '/$key="${v.replaceAll('"', "''")}"';
  final lines = key == 'translation' ? _chunks(text, 58) : _wrap(text, 58);
  return lines.map((l) => '${' ' * 21}$l\n').join();
}

/// Splits [text] into lines of at most [width] characters at spaces (a
/// word longer than [width] is cut).
List<String> _wrap(String text, int width) {
  final out = <String>[];
  var line = '';
  for (final word in text.split(' ')) {
    var w = word;
    while (w.length > width) {
      if (line.isNotEmpty) {
        out.add(line);
        line = '';
      }
      out.add(w.substring(0, width));
      w = w.substring(width);
    }
    if (line.isEmpty) {
      line = w;
    } else if (line.length + 1 + w.length <= width) {
      line = '$line $w';
    } else {
      out.add(line);
      line = w;
    }
  }
  if (line.isNotEmpty || out.isEmpty) out.add(line);
  return out;
}

List<String> _chunks(String text, int width) => [for (var i = 0; i < text.length; i += width) text.substring(i, (i + width).clamp(0, text.length))];

/// Splits a location after [sep] so each line has at most [width]
/// characters.
List<String> _wrapAt(String text, int width, String sep) {
  if (text.length <= width) return [text];
  final out = <String>[];
  var line = '';
  for (final part in text.split(sep)) {
    final piece = line.isEmpty ? part : '$line$sep$part';
    if (piece.length + 1 > width && line.isNotEmpty) {
      out.add('$line$sep');
      line = part;
    } else {
      line = piece;
    }
  }
  if (line.isNotEmpty) out.add(line);
  return out;
}
