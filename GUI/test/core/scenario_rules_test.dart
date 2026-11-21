import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/core/fasta.dart';
import 'package:splace/src/core/genbank.dart';
import 'package:splace/src/core/genbank_markers.dart';
import 'package:splace/src/core/supermatrix.dart';
import 'package:splace/src/state/project_controller.dart';

LocusState _locus(String name, {required bool coding}) => LocusState(name, [SeqRecord('A', 'ATG'), SeqRecord('B', 'ATG')])
  ..aligned = [SeqRecord('A', 'ATGAAA'), SeqRecord('B', 'ATGAAG')]
  ..coding = coding;

void main() {
  test('partition schemes follow the scenario', () {
    final c = ProjectController(autoDetectTools: false);
    addTearDown(c.dispose);
    c.setConcat(c.concat.copyWith(scheme: PartitionScheme.byLocusAndCodon));

    // One coding gene: by codon position or none ("by locus" = one partition).
    c.loci['COI'] = _locus('COI', coding: true);
    expect(c.partitionSchemes, [PartitionScheme.byCodonPosition, PartitionScheme.none]);
    expect(c.effectiveScheme, PartitionScheme.byCodonPosition);

    // One gene without known codon positions: one partition only.
    c.loci['COI'] = _locus('COI', coding: false);
    expect(c.partitionSchemes, [PartitionScheme.none]);
    expect(c.effectiveScheme, PartitionScheme.none);

    // Several genes, none coding: by locus or none.
    c.loci['16S'] = _locus('16S', coding: false);
    expect(c.partitionSchemes, [PartitionScheme.byLocus, PartitionScheme.none]);
    expect(c.effectiveScheme, PartitionScheme.byLocus);

    // Several genes, one coding: every scheme, the choice kept.
    c.loci['CYTB'] = _locus('CYTB', coding: true);
    expect(c.partitionSchemes, hasLength(4));
    expect(c.effectiveScheme, PartitionScheme.byLocusAndCodon);
  });

  test('a fragment partial at the 3\' end keeps only whole codons in codon mode', () {
    final r = GbRecord(
      accession: 'X1',
      organism: 'Caranx hippos',
      lineage: 'Eukaryota; Chordata; Actinopteri; Carangiformes; Carangidae; Caranx',
      definition: 'Caranx hippos COI gene, partial cds; mitochondrial',
      features: [
        GbFeature('source', '1..20', {'organism': 'Caranx hippos', 'organelle': 'mitochondrion'}),
        GbFeature('CDS', '<1..>20', {'gene': 'COX1', 'codon_start': '2', 'transl_table': '2'}),
      ],
      sequence: 'gatgacccacctactaaaaa',
    );
    final cat = buildGeneCatalog([r]);
    final out = extractCodonCandidates([r], cat, {'COI'}, const HeaderTemplate(fields: ['species']), const {}, fallbackCode: 2);
    final dna = out['COI']!.single.dna;
    // 20 bp, codon_start 2: 19 bp, the last (incomplete) base left out.
    expect(dna.length, 18);
    expect(dna.length % 3, 0);
  });
}
