import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/core/fasta.dart';
import 'package:splace/src/core/genbank_markers.dart';
import 'package:splace/src/core/marker_search.dart';
import 'package:splace/src/core/syngenes.dart';

void main() {
  test('SynGenes names become an Entrez query per marker', () {
    expect(syngenesMarkers(OrganelleType.mt), containsAll(['COI', '12S', 'CYTB', '16S', 'ND2']));
    expect(syngenesMarkers(OrganelleType.mt), isNot(contains('CR')));
    expect(syngenesMarkers(OrganelleType.cp), containsAll(['rbcL', 'matK']));

    final rrna = markerNames('12S', OrganelleType.mt);
    expect(rrna.symbols, isNot(contains('trnS')), reason: 'a tRNA listed among the 12S names');
    expect(rrna.phrases, contains('small subunit ribosomal RNA'));

    const s = MarkerSearch(taxon: 'Astyanax', genes: ['COI', 'CYTB']);
    final t = s.termFor('COI');
    expect(t, startsWith('"Astyanax"[Organism] AND ('));
    expect(t, contains('COX1[Gene]'));
    expect(t.toLowerCase(), contains('"cytochrome c oxidase subunit i"[title]'));
    expect(RegExp(r'[(\s]I\[Gene\]').hasMatch(t), isFalse, reason: 'one-letter symbols match other genes');
    expect(t, isNot(contains('COI[Title]')));
    expect(t, contains('1:10000[SLEN]'));
    expect(t, contains('mitochondrion[Filter]'));
    expect(const MarkerSearch(taxon: 'x', genes: ['rbcL'], organelle: OrganelleType.cp, includeGenomes: true).termFor('rbcL'),
        allOf(contains('chloroplast[Filter]'), isNot(contains('[SLEN]'))));
    expect(s.preview.split('\n\n'), hasLength(2));
  });

  test('records of one species are joined, the longest kept', () {
    final r = onePerTaxon({
      'COI': [SeqRecord('Bufonidae_Bufo_bufo', 'ACGTAC'), SeqRecord('Bufonidae_Bufo_bufo', 'ACGTACGTAA'), SeqRecord('Ranidae_Rana_x', 'ACG')],
      'CYTB': [SeqRecord('Ranidae_Rana_x', 'TTTT')],
    });
    expect(r.fastas['COI']!.map((s) => s.sequence), ['ACGTACGTAA', 'ACG']);
    expect(r.dropped, {'COI': 1});
    expect(taxonHeader.fields, ['family', 'species']);
  });
}
