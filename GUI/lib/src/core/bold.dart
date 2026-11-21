/// BOLD Systems (Ratnasingham & Hebert 2007) as a second source of
/// individual markers: the public records of the BOLD v5 portal API
/// (portal.boldsystems.org/api) for the markers BOLD holds — COI and the
/// plant barcodes of the chloroplast (rbcL, matK…) — turned into GenBank
/// records, so they go through the same pipeline as the NCBI ones.
library;

import 'genbank.dart';
import 'genetic_code.dart';
import 'syngenes.dart' show OrganelleType;

/// BOLD marker codes of each SynGenes marker that BOLD holds.
const boldMarkerCodes = <String, List<String>>{
  'COI': ['COI-5P', 'COI-3P'],
  'rbcL': ['rbcL', 'rbcLa'],
  'matK': ['matK'],
  'rpoB': ['rpoB'],
  'rpoC1': ['rpoC1'],
  'ycf1': ['ycf1'],
};

/// Product names written on the BOLD records (and read back by SynGenes).
const _products = <String, String>{
  'COI': 'cytochrome c oxidase subunit I',
  'rbcL': 'ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit',
  'matK': 'maturase K',
  'rpoB': 'RNA polymerase beta subunit',
  'rpoC1': "RNA polymerase beta' subunit",
  'ycf1': 'hypothetical chloroplast RF1',
};

bool isBoldMarker(String gene) => boldMarkerCodes.containsKey(gene);

/// The SynGenes marker of a BOLD marker code (`COI-5P` → `COI`).
String? geneOfBoldCode(String code) {
  for (final e in boldMarkerCodes.entries) {
    if (e.value.any((c) => c.toLowerCase() == code.toLowerCase())) return e.key;
  }
  return null;
}

/// One public BOLD record (one specimen and one marker).
class BoldRecord {
  const BoldRecord({
    required this.processId,
    required this.markerCode,
    required this.sequence,
    this.sampleId = '',
    this.species = '',
    this.identification = '',
    this.kingdom = '',
    this.phylum = '',
    this.klass = '',
    this.order = '',
    this.family = '',
    this.genus = '',
    this.insdc = '',
    this.country = '',
    this.bin = '',
    this.institution = '',
    this.boldName = '',
    this.province = '',
    this.region = '',
    this.site = '',
    this.lat,
    this.lon,
    this.elevation = '',
  });

  /// The same record under [accepted] (the species searched), when BOLD
  /// files it under a synonym; [boldName] keeps the BOLD name.
  BoldRecord renamed(String accepted) => BoldRecord(
        processId: processId,
        markerCode: markerCode,
        sequence: sequence,
        sampleId: sampleId,
        species: accepted,
        identification: identification,
        kingdom: kingdom,
        phylum: phylum,
        klass: klass,
        order: order,
        family: family,
        genus: accepted.split(' ').first,
        insdc: insdc,
        country: country,
        bin: bin,
        institution: institution,
        boldName: name,
        province: province,
        region: region,
        site: site,
        lat: lat,
        lon: lon,
        elevation: elevation,
      );

  /// A row of the portal's TSV, by column name.
  factory BoldRecord.fromRow(Map<String, String> r) {
    String v(String k) {
      final s = (r[k] ?? '').trim();
      // "Unrecoverable" and the like mean no value.
      return s == 'None' || s == 'null' || RegExp(r'^(unrecoverable|unknown|not (specified|available|collected)|n/?a)$', caseSensitive: false).hasMatch(s) ? '' : s;
    }

    // "[4.773, -73.037]"
    final coord = RegExp(r'(-?\d+(?:\.\d+)?)\s*,\s*(-?\d+(?:\.\d+)?)').firstMatch(v('coord'));

    return BoldRecord(
      processId: v('processid'),
      markerCode: v('marker_code'),
      sequence: cleanBoldSequence(v('nuc')),
      sampleId: v('sampleid'),
      species: v('species'),
      identification: v('identification'),
      kingdom: v('kingdom'),
      phylum: v('phylum'),
      klass: v('class'),
      order: v('order'),
      family: v('family'),
      genus: v('genus'),
      insdc: v('insdc_acs'),
      country: v('country/ocean'),
      province: v('province/state'),
      region: v('region'),
      site: v('site'),
      lat: coord == null ? null : double.tryParse(coord.group(1)!),
      lon: coord == null ? null : double.tryParse(coord.group(2)!),
      elevation: v('elev').replaceFirst(RegExp(r'\.0$'), ''),
      bin: v('bin_uri'),
      institution: v('inst'),
    );
  }

  final String processId, markerCode, sequence, sampleId, species, identification;
  final String kingdom, phylum, klass, order, family, genus;

  /// GenBank accession of the same sequence, when BOLD mirrored it.
  final String insdc;
  final String country, bin, institution;

  /// The name in BOLD when the record was renamed to the accepted one.
  final String boldName;

  /// Where the specimen was collected (any of it may be missing).
  final String province, region, site, elevation;
  final double? lat, lon;

  /// INSDC /geo_loc_name: "Colombia: Meta, San Juan de Arama" (the country
  /// alone, or empty, when that is all BOLD has).
  String get geoLocName {
    if (country.isEmpty) return '';
    final rest = [province, region].where((x) => x.isNotEmpty).join(', ');
    return rest.isEmpty ? country : '$country: $rest';
  }

  /// INSDC /lat_lon: "4.7730 N 73.0370 W".
  String get latLon => lat == null || lon == null
      ? ''
      : '${lat!.abs().toStringAsFixed(4)} ${lat! < 0 ? 'S' : 'N'} ${lon!.abs().toStringAsFixed(4)} ${lon! < 0 ? 'W' : 'E'}';

  /// The name of the specimen: the species, else the identification
  /// (e.g. "Rhinella sp.", a genus).
  String get name => species.isNotEmpty ? species : (identification.isNotEmpty ? identification : genus);

  String? get gene => geneOfBoldCode(markerCode);
}

/// Sequence without alignment gaps, blanks and the Ns at the ends.
String cleanBoldSequence(String nuc) =>
    nuc.replaceAll(RegExp(r'[\s\-.]'), '').toLowerCase().replaceFirst(RegExp(r'^n+'), '').replaceFirst(RegExp(r'n+$'), '');

/// Mitochondrial genetic code of a BOLD lineage (NCBI tables).
int boldMitoCode(String phylum, String klass) {
  final ph = phylum.toLowerCase(), cl = klass.toLowerCase();
  if (ph == 'chordata') return cl == 'ascidiacea' ? 13 : 2;
  if (const {'echinodermata', 'hemichordata', 'platyhelminthes'}.contains(ph)) return 9;
  if (const {'cnidaria', 'porifera', 'ctenophora', 'placozoa'}.contains(ph)) return 4;
  if (const {'arthropoda', 'mollusca', 'annelida', 'nematoda', 'bryozoa', 'brachiopoda', 'rotifera', 'tardigrada', 'onychophora', 'nemertea'}.contains(ph)) {
    return 5;
  }
  // Fungi and protists.
  return 4;
}

/// Reading frame (0, 1 or 2) with the fewest internal stops under [code].
int readingFrame(String dna, int code) {
  final stops = stopCodonsFor(code);
  var best = 0, bestStops = 1 << 30;
  for (var f = 0; f < 3; f++) {
    var n = 0;
    final s = dna.toUpperCase();
    // The last codon may be a real stop: not counted.
    for (var i = f; i + 6 <= s.length; i += 3) {
      if (stops.contains(s.substring(i, i + 3))) n++;
    }
    if (n < bestStops) {
      best = f;
      bestStops = n;
    }
  }
  return best;
}

/// A BOLD record as a GenBank record: a source feature (voucher, country,
/// BIN) and one CDS of the marker, partial at both ends, in the reading
/// frame found under the genetic code of its lineage (11 in plastids).
GbRecord boldToGbRecord(BoldRecord b, {String? accession}) {
  final gene = b.gene ?? b.markerCode;
  final seq = b.sequence;
  final plastid = !const {'COI'}.contains(gene);
  final code = plastid ? 11 : boldMitoCode(b.phylum, b.klass);
  final frame = readingFrame(seq, code);
  final end = frame + 3 * ((seq.length - frame) ~/ 3);
  final acc = accession ?? b.processId;
  final lineage = [b.kingdom, b.phylum, b.klass, b.order, b.family, b.genus].where((x) => x.isNotEmpty).join('; ');
  final organelle = plastid ? 'chloroplast' : 'mitochondrion';
  return GbRecord(
    accession: acc,
    organism: b.name,
    lineage: lineage,
    definition: [
      b.name,
      if (b.sampleId.isNotEmpty) 'voucher ${b.sampleId}',
      '${_products[gene] ?? gene} (${b.markerCode}) gene, partial cds; $organelle; BOLD ${b.processId}',
    ].join(' '),
    features: [
      GbFeature('source', '1..${seq.length}', {
        'organism': b.name,
        'organelle': organelle,
        'mol_type': 'genomic DNA',
        if (b.sampleId.isNotEmpty) 'specimen_voucher': b.sampleId,
        if (b.geoLocName.isNotEmpty) 'geo_loc_name': b.geoLocName,
        if (b.latLon.isNotEmpty) 'lat_lon': b.latLon,
        'db_xref': 'BOLD:${b.processId}',
        'note': [
          'BOLD Systems record ${b.processId} (${b.markerCode})',
          if (b.boldName.isNotEmpty) 'identified in BOLD as ${b.boldName}',
          if (b.site.isNotEmpty) 'collection site: ${b.site}',
          if (b.elevation.isNotEmpty) 'elevation: ${b.elevation} m',
          if (b.bin.isNotEmpty) 'BIN ${b.bin}',
          if (b.insdc.isNotEmpty) 'GenBank ${b.insdc}',
          if (b.institution.isNotEmpty) b.institution,
        ].join('; '),
      }),
      GbFeature('CDS', '<1..>$end', {
        'gene': gene,
        'product': _products[gene] ?? gene,
        'codon_start': '${frame + 1}',
        'transl_table': '$code',
      }),
    ],
    sequence: seq,
    source: 'bold',
  );
}

/// The organelle of a marker search that BOLD can serve.
bool boldServes(OrganelleType organelle, Iterable<String> genes) => genes.any((g) => isBoldMarker(g) && (g == 'COI') == (organelle == OrganelleType.mt));
