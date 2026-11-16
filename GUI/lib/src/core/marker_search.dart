/// Searching NCBI for individual markers (COI, 12S, CYTB, rbcL…) instead of
/// whole genomes: one Entrez query per marker, built from the names the
/// SynGenes tables (Rabelo et al. 2024) know for it.
library;

import 'syngenes.dart';
import 'syngenes_data.dart';

/// Markers SynGenes knows for an organelle: protein-coding genes and rRNAs.
List<String> syngenesMarkers(OrganelleType type) =>
    (type == OrganelleType.cp ? chloroplastGenes : mitochondrialGenes).keys.where((k) => k != 'CR').toList();

/// Gene symbols (searched in the [Gene] field) and product names (searched
/// in the title) of [gene] in SynGenes. tRNA names that the tables list
/// among rRNA synonyms (`trnS` for 12S) and symbols of one or two letters
/// or a roman numeral (`I` for COI) are left out: they would bring other
/// genes.
({List<String> symbols, List<String> phrases}) markerNames(String gene, OrganelleType type) {
  final table = type == OrganelleType.cp ? chloroplastGenes : mitochondrialGenes;
  final names = [gene, ...?table[gene]];
  final symbols = <String>[], phrases = <String>[];
  final seen = <String>{};
  for (final raw in names) {
    final n = raw.trim();
    if (n.isEmpty || n.length > 70 || !seen.add(n.toLowerCase())) continue;
    if (RegExp(r'^trn', caseSensitive: false).hasMatch(n)) continue;
    // Entrez phrases cannot hold quotes or brackets.
    if (RegExp(r'["\[\]()]').hasMatch(n)) continue;
    if (n.contains(' ')) {
      phrases.add(n);
    } else if (n.length >= 3 && !RegExp(r'^([IVX]+|\d+)$', caseSensitive: false).hasMatch(n)) {
      // Symbols such as "I" (listed for COI) would match far too much.
      symbols.add(n);
    }
  }
  return (symbols: symbols, phrases: phrases);
}

class MarkerSearch {
  const MarkerSearch({
    required this.taxon,
    required this.genes,
    this.organelle = OrganelleType.mt,
    this.includeGenomes = false,
    this.refseqOnly = false,
    this.maxLength = 10000,
  });

  final String taxon;

  /// SynGenes names of the markers (`COI`, `12S`, `CYTB`).
  final List<String> genes;
  final OrganelleType organelle;

  /// Also records longer than [maxLength] (complete genomes and long
  /// fragments, which hold the markers too).
  final bool includeGenomes;
  final bool refseqOnly;
  final int maxLength;

  /// Entrez query of one marker.
  String termFor(String gene) {
    final n = markerNames(gene, organelle);
    // Symbols in the gene field (/gene), product names in the title.
    final names = [
      for (final s in n.symbols) '$s[Gene]',
      for (final p in n.phrases) '"$p"[Title]',
    ];
    return [
      '"${taxon.trim()}"[Organism]',
      '(${names.join(' OR ')})',
      organelle == OrganelleType.cp
          ? '(chloroplast[Filter] OR chloroplast[Title] OR plastid[Title])'
          : '(mitochondrion[Filter] OR mitochondrial[Title] OR mitochondrion[Title])',
      'biomol_genomic[PROP]',
      if (!includeGenomes) '1:$maxLength[SLEN]',
      if (refseqOnly) 'srcdb_refseq[PROP]',
    ].join(' AND ');
  }

  /// The queries shown before searching, one line per marker.
  String get preview => [for (final g in genes) '$g: ${termFor(g)}'].join('\n\n');
}
