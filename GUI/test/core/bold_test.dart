import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/core/bold.dart';
import 'package:splace/src/core/genbank_markers.dart';
import 'package:splace/src/core/genetic_code.dart';
import 'package:splace/src/core/marker_search.dart';
import 'package:splace/src/core/syngenes.dart';
import 'package:splace/src/services/bold_portal.dart';

/// A BOLD portal TSV (real record ABTVC001-05, Rhinella marina, COI-5P).
const _tsv = 'processid	sampleid	marker_code	nuc	insdc_acs	species	identification	kingdom	phylum	class	order	family	genus	country/ocean	bin_uri	inst\nABTVC001-05	05-PARA-0001	COI-5P	ACCCTGTACTTAATCTTTGGTGCCTGAGCAGGAATAGTGGGAACTGCCCTCAGTCTTCTCATCCGAGCCGAACTGAGCCAACCAGGCTTCCTCCTTGGTGACGATCAAATCTATAATGTTATCGTTACCGCACATGCTTTTGTTATAATTTTCTTTATGGTTATGCCAATCCTAATCGGCGGCTTTGGTAACTGACTTGTGCCTTTAATAATTGGAGCCCCCGATATGGCTTTCCCTCGTATAAACAACATAAGCTTTTGATTACTGCCCCCATCTTTCCTCCTCCTTCTGGCCTCCGCTGGAGTAGAAGCCGGAGCGGGGACCGGCTGAACTGTTTATCCCCCATTAGCCGGAAATCTCGCCCATGCAGGGCCGTCTGTAGACTTAACTATTTTTTCCCTTCATTTAGCAGGTGTTTCATCCATTCTAGGAGCAATTAATTTTATTACCACAACCCTGAATATGAAACCTCCATCAATAACTCAATATCAAACACCTTTATTTGTCTGGTCTGTGTTGATTACCGCTGTTCTACTACTTCTCTCCCTACCAGTTCTCGCAGCTGGAATTACAATACTTTTAACTGATCGCAACTTAAATACAACATTCTTTGATCCCGCTGGAGGAGGAGACCCTATCTTATATCAACATCTTTTC		Rhinella marina	Rhinella marina	Animalia	Chordata	Amphibia	Anura	Bufonidae	Rhinella	Costa Rica	BOLD:AAB1186	Centre for Biodiversity Genomics\nX2		16S	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT		Rhinella marina										\nX3		COI-5P	-CCCTGTACTTAATCTTTGGTGCCTGAGCAGGAATAGTGGGAACTGCCCTCAGTCTTCTCATCCGAGCCGAACTGAGCCAACCAGGCTTCCTCCTTGGTGACGATCAAATCTATAATGTTATCGTTACCGCACATGCTTTTGTTATAATTTTCTTTATGGTTATGCCAATCCTAATCGGCGGCTTTGGTAACTGACTTGTGCCTTTAATAATTGGAGCCCCCGATATGGCTTTCCCTCGTATAAACAACATAAGCTTTTGATTACTGCCCCCATCTTTCCTCCTCCTTCTGGCCTCCGCTGGAGTAGAAGCCGGAGCGGGGACCGGCTGAACTGTTTATCCCCCATTAGCCGGAAATCTCGCCCATGCAGGGCCGTCTGTAGACTTAACTATTTTTTCCCTTCATTTAGCAGGTGTTTCATCCATTCTAGGAGCAATTAATTTTATTACCACAACCCTGAATATGAAACCTCCATCAATAACTCAATATCAAACACCTTTATTTGTCTGGTCTGTGTTGATTACCGCTGTTCTACTACTTCTCTCCCTACCAGTTCTCGCAGCTGGAATTACAATACTTTTAACTGATCGCAACTTAAATACAACATTCTTTGATCCCGCTGGAGGAGGAGACCCTATCTTATATCAACATCTTTTC--	MN123456	Rhinella marina			Chordata	Amphibia						\n';

void main() {
  test('BOLD TSV: only the markers asked for, with a sequence; gaps removed', () {
    final rows = parseBoldTsv(_tsv, {'COI-5P', 'COI-3P'});
    expect(rows.map((r) => r.processId), ['ABTVC001-05', 'X3']);
    expect(rows.first.gene, 'COI');
    expect(rows.first.name, 'Rhinella marina');
    expect(rows.first.country, 'Costa Rica');
    expect(rows.first.bin, 'BOLD:AAB1186');
    expect(rows.last.sequence.contains('-'), isFalse);
    expect(rows.last.insdc, 'MN123456');
  });

  test('a BOLD record becomes a GenBank record SynGenes reads as COI, in frame', () {
    final b = parseBoldTsv(_tsv, {'COI-5P'}).first;
    final r = boldToGbRecord(b);
    expect(r.source, 'bold');
    expect(r.accession, 'ABTVC001-05');
    expect(r.lineage, contains('Bufonidae'));
    final cds = r.features.firstWhere((f) => f.type == 'CDS');
    expect(cds.q('transl_table'), '2');
    final frame = int.parse(cds.q('codon_start')!) - 1;
    // No internal stop in the vertebrate mitochondrial code.
    final protein = translate(r.sequence.substring(frame, frame + 3 * ((r.sequence.length - frame) ~/ 3)), 2);
    expect(protein.substring(0, protein.length - 1).contains('*'), isFalse);
    final catalog = buildGeneCatalog([r]);
    expect(catalog.byType['CDS']!.keys, contains('COI'));
  });

  test('BOLD codes, organelles and the search', () {
    expect(geneOfBoldCode('rbcLa'), 'rbcL');
    expect(geneOfBoldCode('COI-3P'), 'COI');
    expect(boldServes(OrganelleType.mt, ['COI', 'CYTB']), isTrue);
    expect(boldServes(OrganelleType.mt, ['CYTB']), isFalse);
    expect(boldServes(OrganelleType.cp, ['matK']), isTrue);
    expect(boldMitoCode('Arthropoda', 'Insecta'), 5);
    expect(boldMitoCode('Echinodermata', ''), 9);
    final s = MarkerSearch(taxon: 'Rhinella marina', genes: const ['COI', 'CYTB'], synonyms: const ['Rhinella marina', 'Bufo marinus'], bold: true);
    expect(s.boldGenes, ['COI']);
    expect(s.boldCodes, {'COI-5P', 'COI-3P'});
    expect(s.boldNames, ['Rhinella marina', 'Bufo marinus']);
    expect(s.preview, contains('BOLD Systems'));
    expect(MarkerSearch(taxon: 'Rhinella', genes: const ['COI']).boldGenes, isEmpty);
  });

  test('a record BOLD files under a synonym takes the accepted name; placeholder countries are dropped', () {
    final b = BoldRecord.fromRow({'processid': 'X9', 'marker_code': 'COI-5P', 'nuc': 'ACGT' * 60, 'species': 'Chaunus marinus', 'country/ocean': 'Unrecoverable'});
    expect(b.country, '');
    final r = b.renamed('Rhinella marina');
    expect(r.name, 'Rhinella marina');
    expect(r.genus, 'Rhinella');
    expect(r.boldName, 'Chaunus marinus');
    expect(boldToGbRecord(r).features.first.q('note'), contains('identified in BOLD as Chaunus marinus'));
  });

  test('the locality of a BOLD record goes to /geo_loc_name and /lat_lon', () {
    final b = BoldRecord.fromRow({
      'processid': 'X10', 'marker_code': 'COI-5P', 'nuc': 'ACGT' * 60, 'species': 'Rhinella marina',
      'country/ocean': 'Colombia', 'province/state': 'Meta', 'region': 'San Juan de Arama', 'site': 'Finca Los Almendros',
      'coord': '[3.228, -73.852]', 'elev': '325.0',
    });
    expect(b.geoLocName, 'Colombia: Meta, San Juan de Arama');
    expect(b.latLon, '3.2280 N 73.8520 W');
    final r = boldToGbRecord(b);
    expect(r.locality, 'Colombia: Meta, San Juan de Arama');
    expect(r.latLon, '3.2280 N 73.8520 W');
    expect(r.features.first.q('note'), allOf(contains('Finca Los Almendros'), contains('325 m')));
    // Nothing known: nothing written (the locality is optional).
    final bare = boldToGbRecord(BoldRecord.fromRow({'processid': 'X11', 'marker_code': 'COI-5P', 'nuc': 'ACGT' * 60, 'species': 'Rhinella marina'}));
    expect(bare.locality, isNull);
    expect(bare.latLon, isNull);
  });
}
