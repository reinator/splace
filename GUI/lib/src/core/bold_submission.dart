/// BOLD Systems specimen submission: the columns of the official BOLD5
/// "BCDM Specimen Submission – Single Page" template (sheet "Submissions",
/// boldsystems.org/resource-hub/templates, 2025), filled for barcode
/// sequences identified in SPLACE. Fields SPLACE cannot know (voucher,
/// collection data, collectors…) stay empty for the user to complete.
library;

/// Columns of the template, in its order.
const boldColumns = [
  'SAMPLEID', 'FIELDID', 'MUSEUMID', 'COLLECTION_CODE', 'INST', 'FUNDING_SRC',
  'PHYLUM', 'CLASS', 'ORDER', 'FAMILY', 'SUBFAMILY', 'TRIBE', 'GENUS', 'SPECIES', 'SUBSPECIES',
  'IDENTIFIED_BY', 'IDENTIFICATION_METHOD', 'TAXONOMY_NOTES',
  'SEX', 'REPRODUCTION', 'LIFE_STAGE', 'SHORT_NOTE', 'NOTES', 'VOUCHER_TYPE', 'TISSUE_TYPE', 'SPECIMEN_LINKOUT',
  'ASSOCIATED_TAXA', 'ASSOCIATED_SPECIMENS',
  'COLLECTORS', 'COLLECTION_DATE_START', 'COLLECTION_DATE_END', 'COUNTRY/OCEAN', 'PROVINCE/STATE', 'REGION', 'SECTOR', 'SITE',
  'COORD:LAT', 'COORD:LON', 'ELEV', 'DEPTH', 'ELEV_ACCURACY', 'DEPTH_ACCURACY', 'COORD_SOURCE', 'COORD_ACCURACY',
  'COLLECTION_TIME', 'HABITAT', 'SAMPLING_PROTOCOL', 'COLLECTION_NOTES', 'SITE_CODE', 'COLLECTION_EVENT_ID',
];

/// Name of the template's only sheet.
const boldSheet = 'Submissions';

/// One specimen: its sample id (also the sequence header) and what SPLACE
/// knows of it.
class BoldSpecimen {
  const BoldSpecimen({
    required this.sampleId,
    this.phylum = '',
    this.klass = '',
    this.order = '',
    this.family = '',
    this.genus = '',
    this.species = '',
    this.identificationMethod = '',
    this.taxonomyNotes = '',
  });

  final String sampleId, phylum, klass, order, family, genus, species, identificationMethod, taxonomyNotes;

  /// The row in [boldColumns] order.
  List<String> get row {
    final values = {
      'SAMPLEID': sampleId,
      'PHYLUM': phylum,
      'CLASS': klass,
      'ORDER': order,
      'FAMILY': family,
      'GENUS': genus,
      'SPECIES': species,
      'IDENTIFICATION_METHOD': identificationMethod,
      'TAXONOMY_NOTES': taxonomyNotes,
    };
    return [for (final c in boldColumns) values[c] ?? ''];
  }
}

/// Header and rows of the submission sheet.
List<List<String>> boldSheetRows(List<BoldSpecimen> specimens) => [boldColumns, for (final s in specimens) s.row];

/// The sequences for BOLD's sequence upload: one record per specimen, the
/// header being its sample id (the SAMPLEID column), gaps removed.
String boldFasta(Map<String, String> sequences) => [
      for (final e in sequences.entries) '>${e.key}\n${e.value.replaceAll(RegExp(r'[-.\s]'), '').toUpperCase()}',
    ].join('\n');
