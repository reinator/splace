import 'dart:convert';

import 'package:archive/archive.dart';
import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/core/bold_submission.dart';
import 'package:splace/src/core/genbank.dart';
import 'package:splace/src/core/genbank_writer.dart';
import 'package:splace/src/core/xlsx.dart';

void main() {
  const gb = '''LOCUS       NEW001     60 bp    DNA     circular VRT 01-JAN-2026
DEFINITION  Caranx sp. mitochondrion, complete genome.
ACCESSION   NEW001
VERSION     NEW001.1
SOURCE      mitochondrion Caranx sp.
  ORGANISM  Caranx sp.
            Eukaryota; Metazoa; Chordata.
REFERENCE   1  (bases 1 to 60)
  AUTHORS   Someone A.
  TITLE     Direct Submission
COMMENT     Assembled with MitoFinder.
FEATURES             Location/Qualifiers
     source          1..60
                     /organism="Caranx sp."
                     /mol_type="genomic DNA"
                     /db_xref="taxon:1"
     gene            1..30
                     /gene="cox1"
     CDS             1..30
                     /gene="cox1"
                     /codon_start=1
                     /transl_table=2
                     /product="cytochrome c oxidase subunit I"
     rRNA            31..60
                     /product="16S ribosomal RNA"
ORIGIN      
        1 atggcatgga tggcatggat ggcatggatg gcatggatgg catggatggc atggatggca
//
''';

  test('corrected GenBank: organism, lineage, standard gene names, header kept', () {
    final r = parseGenBank(gb)!
      ..editedOrganism = 'Caranx hippos'
      ..ranks = const TaxonomyRanks(kingdom: 'Animalia', phylum: 'Chordata', klass: 'Actinopterygii', order: 'Carangiformes', family: 'Carangidae', source: 'GBIF');
    expect(r.header, contains('REFERENCE   1'));
    final text = writeGenBank(r, geneNames: {2: 'COI', 3: '16S'}, note: 'Corrected with SPLACE.');
    final back = parseGenBank(text)!;
    expect(back.organism, 'Caranx hippos');
    expect(back.lineage, 'Animalia; Chordata; Actinopterygii; Carangiformes; Carangidae; Caranx');
    expect(back.definition, 'Caranx hippos mitochondrion, complete genome.');
    expect(back.sequence, r.sequence);
    expect(back.features.length, r.features.length);
    expect(back.features[0].q('organism'), 'Caranx hippos');
    // The gene feature follows its CDS; the original name in a note.
    expect(back.features[1].q('gene'), 'COI');
    expect(back.features[2].q('gene'), 'COI');
    expect(back.features[2].q('note'), contains('original annotation: cox1'));
    expect(back.features[2].q('codon_start'), '1');
    expect(back.features[3].q('gene'), '16S');
    expect(text, contains('/codon_start=1\n'));
    expect(text, contains('REFERENCE   1'));
    expect(text, contains('COMMENT     Corrected with SPLACE.'));
    expect(text, contains('Assembled with MitoFinder.'));
    expect(text, contains('db_xref="taxon:1"'));
    expect(text.trim().endsWith('//'), isTrue);
  });

  test('BOLD5 single-page template: the 50 columns in order', () {
    expect(boldColumns, hasLength(50));
    expect(boldColumns.first, 'SAMPLEID');
    expect(boldColumns[13], 'SPECIES');
    expect(boldColumns.last, 'COLLECTION_EVENT_ID');
    const s = BoldSpecimen(sampleId: 'ITV001', phylum: 'Chordata', family: 'Carangidae', genus: 'Caranx', species: 'Caranx hippos');
    final row = s.row;
    expect(row, hasLength(50));
    expect(row[boldColumns.indexOf('GENUS')], 'Caranx');
    expect(boldFasta({'ITV001': 'acg-t'}), '>ITV001\nACGT');
  });

  test('xlsx: one sheet, inline strings, header in bold', () {
    final bytes = buildXlsx(boldSheet, boldSheetRows([const BoldSpecimen(sampleId: '0012 & <x>')]));
    final zip = ZipDecoder().decodeBytes(bytes);
    String read(String n) => utf8.decode(zip.findFile(n)!.content as List<int>);
    expect(read('xl/workbook.xml'), contains('name="Submissions"'));
    final sheet = read('xl/worksheets/sheet1.xml');
    expect(sheet, contains('<c r="A1" t="inlineStr" s="1"><is><t xml:space="preserve">SAMPLEID</t>'));
    expect(sheet, contains('<c r="AX1"'));
    expect(sheet, contains('0012 &amp; &lt;x&gt;'));
    expect(read('[Content_Types].xml'), contains('/xl/worksheets/sheet1.xml'));
    expect(xlsxColumn(0), 'A');
    expect(xlsxColumn(25), 'Z');
    expect(xlsxColumn(26), 'AA');
    expect(xlsxColumn(49), 'AX');
  });
}
