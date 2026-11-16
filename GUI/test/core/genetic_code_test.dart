import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/core/fasta.dart';
import 'package:splace/src/core/genetic_code.dart';
import 'package:splace/src/core/supermatrix.dart';

/// Same cases as `runCodonAwareSelfTest` in the HTML app.
void main() {
  test('CDS validation: terminal and incomplete stops', () {
    final terminal = validateCds('T', 'terminal-stop', 'ATGTAG', 2)!;
    expect((terminal.validDna, terminal.protein, terminal.terminalStopRemoved), ('ATG', 'M', true));

    final incomplete = validateCds('T', 'incomplete-stop', 'ATGTA', 2)!;
    expect((incomplete.validDna, incomplete.protein, incomplete.normalized), ('ATG', 'M', true));

    expect(validateCds('T', 'not-multiple', 'ATGCC', 2), isNull);
  });

  test('NCBI tables 2 and 11', () {
    expect(validateCds('T', 'mt', 'ATATGATAA', 2)!.protein, 'MW');
    expect(validateCds('T', 'cp', 'ATATTTTAA', 11)!.protein, 'IF');
    expect(stopCodonsFor(2), {'TAA', 'TAG', 'AGA', 'AGG'});
    expect(stopCodonsFor(1), {'TAA', 'TAG', 'TGA'});
    expect(translate('ATG---TGA', 2), 'M-W');
  });

  test('internal stop policy', () {
    const dna = 'ATGTAAATGTAA'; // M * M *
    expect(validateCds('T', 'x', dna, 1), isNull);
    final kept = validateCds('T', 'x', dna, 1, stopPolicy: InternalStopPolicy.allowWithWarning)!;
    expect((kept.protein, kept.internalStops), ('M*M', 1));
  });

  test('back-translation keeps codons and gaps in triplets', () {
    final bt = backTranslate('M-F', 'ATGTTT', 1);
    expect(bt.dna, 'ATG---TTT');
    expect(bt.mismatches, 0);
    expect(bt.frameDisrupted, isFalse);

    final mismatch = backTranslate('MW', 'ATGTTT', 1);
    expect(mismatch.mismatches, 1);

    final integrity = codonIntegrity([SeqRecord('a', bt.dna), SeqRecord('b', 'ATA---TTT')]);
    expect(integrity.ok, isTrue);
    expect(codonIntegrity([SeqRecord('a', 'AT--GTTTT')]).gapsInTriplets, isFalse);

    final aln = backTranslateAlignment(
      [SeqRecord('a', 'M-F'), SeqRecord('b', 'MLF')],
      {'a': 'ATGTTT', 'b': 'ATGCTGTTT'},
      const {},
      1,
    );
    expect(aln.dna.map((r) => r.sequence), ['ATG---TTT', 'ATGCTGTTT']);
  });

  test('codon-position removal and partitions', () {
    expect(removeCodonPositions([SeqRecord('a', 'ATGCCCTTT')], {3}).single.sequence, 'ATCCTT');
    final m = buildSupermatrix(
      [LocusAlignment('geneA', [SeqRecord('tax1', 'ATGCCCTTT'), SeqRecord('tax2', 'ATA---TTT')])],
      const ConcatOptions(excludePositions: {3}, replaceGapsWithMissing: true, outgroups: ['tax2']),
    );
    expect(m.sequences['tax2'], 'AT??TT');
    expect(m.partitions.map((p) => '${p.name}=${p.rangeString}'), ['geneA_pos1=1-5\\2', 'geneA_pos2=2-6\\2']);
    expect(m.toNexus(), contains("OUTGROUP 'tax2';"));
    expect(m.toPartitionTable(), 'DNA, geneA_pos1 = 1-5\\2\nDNA, geneA_pos2 = 2-6\\2');

    final protein = buildSupermatrix(
      [LocusAlignment('COI', [SeqRecord('a', 'MLF'), SeqRecord('b', 'M-F')])],
      const ConcatOptions(dataType: SeqDataType.protein),
    );
    expect(protein.partitions.single.name, 'COI', reason: 'no codon partitions for proteins');
    expect(protein.toNexus(), contains('DATATYPE=PROTEIN'));
    expect(protein.toRaxmlPartitions(), 'LG+G, COI = 1-3\n');
  });
}
