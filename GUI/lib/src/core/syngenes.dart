/// Gene-name standardization with the SynGenes tables (Rabelo et al. 2024),
/// ported from the HTML app.
library;

import 'genbank.dart';
import 'syngenes_data.dart';

enum OrganelleType { mt, cp }

const mtDefaultGenes = ['COI', 'COII', 'COIII', 'CYTB', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6', 'ATP6', 'ATP8'];
const cpDefaultGenes = [
  'rbcL', 'matK', 'ndhF', 'atpB', 'psaA', 'psbA', 'psbB', 'psbC', 'psbD', 'psbE', 'psbF', 'psbH', //
  'psbI', 'psbJ', 'psbK', 'psbL', 'psbM', 'psbN', 'psbT',
];

class _Lookup {
  _Lookup(Map<String, List<String>> table) {
    for (final e in table.entries) {
      keysUpper.putIfAbsent(e.key.toUpperCase(), () => e.key);
      for (final s in e.value) {
        synonyms.putIfAbsent(s.toLowerCase().trim(), () => e.key);
      }
    }
  }
  final keysUpper = <String, String>{};
  final synonyms = <String, String>{};
}

final _mt = _Lookup(mitochondrialGenes);
final _cp = _Lookup(chloroplastGenes);

bool isKnownGeneName(String? name) =>
    name != null && (mitochondrialGenes.containsKey(name) || chloroplastGenes.containsKey(name));

/// Standard name for a CDS/rRNA from its /gene and /product qualifiers:
/// product synonyms first, then the gene against the keys and synonyms of
/// the organelle's table, then the other organelle. Falls back to the raw
/// name. [overrides] holds names the user mapped by hand.
String? standardizeGeneName(String? rawName, String? product, OrganelleType type, [Map<String, String> overrides = const {}]) {
  final overrideKey = rawName ?? product;
  if (overrideKey != null && overrides.containsKey(overrideKey)) return overrides[overrideKey];
  final primary = type == OrganelleType.cp ? _cp : _mt;
  final alt = type == OrganelleType.cp ? _mt : _cp;
  if (product != null) {
    final p = product.toLowerCase().trim();
    final hit = primary.synonyms[p] ?? alt.synonyms[p];
    if (hit != null) return hit;
  }
  if (rawName != null) {
    final upper = rawName.toUpperCase();
    final key = primary.keysUpper[upper] ?? alt.keysUpper[upper];
    if (key != null) return key;
    final raw = rawName.toLowerCase().trim();
    final hit = primary.synonyms[raw] ?? alt.synonyms[raw];
    if (hit != null) return hit;
  }
  return rawName ?? product;
}

const _aaNames = {
  'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys', 'E': 'Glu', 'Q': 'Gln', //
  'G': 'Gly', 'H': 'His', 'I': 'Ile', 'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
  'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',
};

String _cap(String s) => s.isEmpty ? s : '${s[0].toUpperCase()}${s.substring(1).toLowerCase()}';

/// tRNA-Xxx names, splitting Leu1/Leu2 and Ser1/Ser2 by anticodon.
String standardizeTrnaName(GbFeature f) {
  final product = f.q('product') ?? '';
  final gene = f.q('gene') ?? '';
  final note = f.q('note') ?? '';
  final anticodon = f.q('anticodon') ?? '';
  final raw = product.isNotEmpty ? product : gene;
  if (raw.isEmpty) {
    final m = RegExp(r'tRNA-(\w+)', caseSensitive: false).firstMatch(note) ??
        RegExp(r'transfer\s+RNA[- ]+(\w+)', caseSensitive: false).firstMatch(note);
    return m != null ? 'tRNA-${_cap(m.group(1)!)}' : 'tRNA-Unknown';
  }
  final rawLower = raw.toLowerCase();
  final geneLower = gene.toLowerCase();
  final combined = '$anticodon $product $note $gene'.toLowerCase();

  if (rawLower.contains('leu') || geneLower.startsWith('trnl')) {
    if (geneLower == 'trnl1' || geneLower.contains('trnl(uaa)') || geneLower.contains('trnl(taa)')) return 'tRNA-Leu1';
    if (geneLower == 'trnl2' || geneLower.contains('trnl(uag)') || geneLower.contains('trnl(tag)')) return 'tRNA-Leu2';
    if (combined.contains('taa') || combined.contains('uaa')) return 'tRNA-Leu1';
    if (combined.contains('tag') || combined.contains('uag') || combined.contains('cun')) return 'tRNA-Leu2';
    return 'tRNA-Leu1';
  }
  if (rawLower.contains('ser') || geneLower.startsWith('trns')) {
    if (geneLower == 'trns1' || geneLower.contains('trns(uga)') || geneLower.contains('trns(tga)')) return 'tRNA-Ser1';
    if (geneLower == 'trns2' || geneLower.contains('trns(gcu)') || geneLower.contains('trns(gct)')) return 'tRNA-Ser2';
    if (combined.contains('tga') || combined.contains('uga') || combined.contains('ucn')) return 'tRNA-Ser1';
    if (['gct', 'gcu', 'agn', 'agy'].any(combined.contains)) return 'tRNA-Ser2';
    return 'tRNA-Ser1';
  }
  final aa = RegExp(r'tRNA-(\w+)', caseSensitive: false).firstMatch(raw);
  if (aa != null) return 'tRNA-${_cap(aa.group(1)!)}';
  final trn = RegExp(r'^trn([A-Za-z])$', caseSensitive: false).firstMatch(raw);
  if (trn != null) {
    final letter = trn.group(1)!.toUpperCase();
    return 'tRNA-${_aaNames[letter] ?? letter}';
  }
  final transfer = RegExp(r'transfer\s+RNA[- ]+(\w+)', caseSensitive: false).firstMatch(raw);
  if (transfer != null) return 'tRNA-${_cap(transfer.group(1)!)}';
  return raw;
}

/// Mitochondrion or plastid, by counting genes that match each table.
OrganelleType detectOrganelleType(Iterable<GbRecord> records) {
  var mt = 0, cp = 0;
  final mtKeys = {for (final k in mitochondrialGenes.keys) k.toLowerCase()};
  final cpKeys = {for (final k in chloroplastGenes.keys) k.toLowerCase()};
  for (final r in records) {
    for (final f in r.features) {
      final g = (f.q('gene') ?? '').toLowerCase();
      final p = (f.q('product') ?? '').toLowerCase();
      if (mtKeys.contains(g) || mtKeys.contains(p)) mt++;
      if (cpKeys.contains(g) || cpKeys.contains(p)) cp++;
    }
  }
  return cp > mt ? OrganelleType.cp : OrganelleType.mt;
}

/// Gene name of a CDS used in codon mode, also trying /note, /locus_tag
/// and /gene_synonym when /gene and /product are missing.
String? resolveCodonGeneName(GbFeature f, OrganelleType type, [Map<String, String> overrides = const {}]) {
  final gene = f.q('gene');
  final product = f.q('product');
  final direct = standardizeGeneName(gene, product, type, overrides);
  if (direct != null) return direct;
  for (final fallback in [f.q('note'), f.q('locus_tag'), f.q('gene_synonym')].whereType<String>()) {
    final hit = standardizeGeneName(gene ?? fallback, fallback, type, overrides);
    if (hit != null) return hit;
  }
  return null;
}
