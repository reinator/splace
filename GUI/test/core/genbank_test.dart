import 'dart:io';

import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/core/genbank.dart';
import 'package:splace/src/core/genbank_markers.dart';
import 'package:splace/src/core/syngenes.dart';

/// 60 bp toy mitogenome: a forward CDS with /codon_start=2, a CDS split
/// over a join on the reverse strand, two tRNA-Leu and an rRNA.
String _record(String accession, String organism, {bool duplicateCox = false}) => '''
LOCUS       ${accession.split('.').first}               60 bp    DNA     circular VRT 01-JAN-2024
DEFINITION  $organism mitochondrion, complete genome.
ACCESSION   ${accession.split('.').first}
VERSION     $accession
SOURCE      mitochondrion $organism
  ORGANISM  $organism
            Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Amphibia;
            Anura; Bufonidae; Rhinella.
FEATURES             Location/Qualifiers
     source          1..60
                     /organism="$organism"
     CDS             1..13
                     /gene="COX1"
                     /codon_start=2
                     /transl_table=2
                     /product="cytochrome c oxidase subunit I"
${duplicateCox ? '''     CDS             1..10
                     /gene="COX1"
                     /product="cytochrome c oxidase subunit I"
''' : ''}     CDS             complement(join(20..25,30..32))
                     /gene="nad6"
                     /product="NADH dehydrogenase subunit 6"
     tRNA            40..45
                     /product="tRNA-Leu"
                     /anticodon=(pos:41..43,aa:Leu,seq:taa)
     tRNA            46..50
                     /product="tRNA-Leu"
                     /note="codons recognized: CUN"
     rRNA            51..60
                     /product="16S ribosomal RNA"
ORIGIN
        1 aatggccatt aaaacccggg tttaaacccc ggggttttaa aaccccgggg ttttaaaacc
//
''';

void main() {
  test('parses records, locations and taxonomy', () {
    final records = parseGenBankFile(_record('NC_000001.1', 'Rhinella marina') + _record('NC_000002.1', 'Rhinella major'));
    expect(records, hasLength(2));
    final r = records.first;
    expect(r.accession, 'NC_000001.1');
    expect(r.organism, 'Rhinella marina');
    expect(r.sequence, hasLength(60));
    expect(r.featureCounts, (pcgs: 2, rrnas: 1, trnas: 2));
    expect(r.taxonomy.family, 'Bufonidae');
    expect(r.taxonomy.genus, 'Rhinella');
    expect(r.taxonomy.species, 'marina');

    final cds = r.features.firstWhere((f) => f.q('gene') == 'nad6');
    // complement(join(20..25,30..32)) = revcomp(seq[19..25] + seq[29..32])
    final forward = r.sequence.substring(19, 25) + r.sequence.substring(29, 32);
    expect(extractLocation(r.sequence, cds.location), reverseComplement(forward).toUpperCase());
    expect(locationLength(cds.location), 9);
    expect(locationLength('join(complement(1..3),complement(10..12))'), 6);
  });

  test('standardizes gene names with SynGenes', () {
    expect(standardizeGeneName('COX1', 'cytochrome c oxidase subunit I', OrganelleType.mt), 'COI');
    expect(standardizeGeneName('nad6', 'NADH dehydrogenase subunit 6', OrganelleType.mt), 'ND6');
    expect(standardizeGeneName(null, '16S ribosomal RNA', OrganelleType.mt), '16S');
    expect(standardizeGeneName('weird1', null, OrganelleType.mt), 'weird1');
    expect(standardizeGeneName('weird1', null, OrganelleType.mt, {'weird1': 'COII'}), 'COII');

    final r = parseGenBankFile(_record('NC_000001.1', 'Rhinella marina')).first;
    final trnas = r.features.where((f) => f.type == 'tRNA').map(standardizeTrnaName).toList();
    expect(trnas, ['tRNA-Leu1', 'tRNA-Leu2']);
    expect(detectOrganelleType([r]), OrganelleType.mt);
  });

  test('catalog, headers, duplicates and codon candidates', () {
    final records = parseGenBankFile(
      _record('NC_000001.1', 'Rhinella marina', duplicateCox: true) + _record('NC_000002.1', 'Rhinella major'),
    );
    final cat = buildGeneCatalog(records);
    expect(cat.featureTypeCounts, {'CDS': 5, 'tRNA': 4, 'rRNA': 2});
    final visible = cat.visible({'CDS'});
    expect(visible.keys, ['COI', 'ND6']);
    expect(visible['COI']!.byRecord[0], hasLength(2), reason: 'duplicated COX1 in the first genome');
    expect(cat.unknown, isEmpty);

    const header = HeaderTemplate(fields: ['species', 'accession']);
    expect(header.build(records.first), 'Rhinella_marina_NC_000001.1');
    expect(header.duplicates(records), isEmpty);
    expect(const HeaderTemplate(fields: ['genus']).duplicates(records), {'Rhinella': ['NC_000001.1', 'NC_000002.1']});

    // Longest copy by default, the user's choice otherwise.
    var fastas = extractGeneFastas(records, cat, {'CDS'}, {'COI'}, header, const {});
    expect(fastas['COI']!.first.sequence, hasLength(13));
    fastas = extractGeneFastas(records, cat, {'CDS'}, {'COI'}, header, {duplicateChoiceKey(records.first, 'COI'): 1});
    expect(fastas['COI']!.first.sequence, hasLength(10));

    final codon = extractCodonCandidates(records, cat, {'COI'}, header, const {}, fallbackCode: 1);
    final c = codon['COI']!.first;
    expect(c.dna, records.first.sequence.substring(1, 13).toUpperCase(), reason: '/codon_start=2 drops one base');
    expect(c.code, 2, reason: '/transl_table=2');
    expect(c.codeFromFeature, isTrue);

    final tree = buildTaxonomyTree(records..first.ranks = const TaxonomyRanks(kingdom: 'Animalia', family: 'Bufonidae'));
    expect(tree.children.keys, ['Kingdom:Animalia']);
  });

  test('partial features (barcode fragments) are counted per gene', () {
    expect(isPartialLocation('<1..>654'), isTrue);
    expect(isPartialLocation('complement(<20..32)'), isTrue);
    expect(isPartialLocation('1..13'), isFalse);
    final partial = _record('MK000001.1', 'Caranx hippos').replaceFirst('     CDS             1..13', '     CDS             <1..>13');
    final records = parseGenBankFile(partial + _record('NC_000002.1', 'Rhinella major'));
    final cat = buildGeneCatalog(records);
    expect(partialCounts(records, cat, {'CDS'}, {'COI', 'ND6'}, const {}), {'COI': 1});
  });

  test('curated /label names a feature without /gene', () {
    expect(geneFromLabel('ycf1 CDS'), 'ycf1');
    expect(geneFromLabel(' trnC-GCA tRNA'), 'trnC-GCA');
    expect(geneFromLabel('ycf2 gene - Truncated'), 'ycf2');
    expect(geneFromLabel('rps14 cds'), 'rps14');
    expect(geneFromLabel('rps12'), 'rps12');
    expect(geneFromLabel('rrn4.5 rRNA'), 'rrn4.5');
    expect(geneFromLabel('ORF 4 (frame 2)'), isNull);
    expect(geneFromLabel(''), isNull);

    final record = parseGenBankFile(_record('NC_000001.1', 'Rhinella marina').replaceFirst(
      '''     tRNA            40..45
                     /product="tRNA-Leu"''',
      '''     CDS             40..45
                     /created_by="C0566"
                     /label="ycf1 CDS"
     CDS             40..45
                     /gene="clpP1"
                     /label="clpP CDS"''',
    )).first;
    final cds = record.features.where((f) => f.type == 'CDS').toList();
    expect(cds[2].q('gene'), 'ycf1', reason: 'only /label');
    expect(cds[3].q('gene'), 'clpP1', reason: '/gene wins over /label');
    expect(cds[2].q('label'), 'ycf1 CDS', reason: 'the label stays');

    // Real curated plastomes (GeSeq + Geneious), when the folder is present.
    final file = File('teste/Plastomas_Uleiorchis/Uleiorchis_ITV01433_25_Aguas.gb');
    if (file.existsSync()) {
      final real = parseGenBankFile(file.readAsStringSync());
      final cat = buildGeneCatalog(real, organelle: OrganelleType.cp);
      final genes = cat.visible({'CDS', 'tRNA'}).keys.toSet();
      expect(genes, containsAll(['ycf1', 'ycf2', 'rps12', 'trnC-GCA', 'trnQ-UUG']));
    }
  });

  test('marker count statistics use the mode', () {
    final s = MarkerCountStats.of([13, 13, 12, 14, 13]);
    expect((s.min, s.max, s.mode), (12, 14, 13));
  });
}
