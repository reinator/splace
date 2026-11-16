import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/core/alignment_tools.dart';
import 'package:splace/src/core/fasta.dart';
import 'package:splace/src/core/supermatrix.dart';

void main() {
  // A mitochondrial CDS: TGA is Trp in code 2 and a stop in code 1.
  const mito = 'ATAGGTTGACTTCTGAAAATTAGTTCTGAGGTTCGA';

  test('only the project genetic code is used (no code detection)', () {
    final aln = [for (var i = 0; i < 6; i++) SeqRecord('s$i', mito)];
    // Code 2: TGA is Trp, no internal stop.
    final r2 = detectCoding('COI', aln, 2);
    expect(r2.coding, isTrue);
    expect(r2.frame, 0);
    expect(r2.stops, 0);
    // Code 1: TGA stays a stop (another code is not tried); COI is still
    // coding by name.
    final r1 = detectCoding('COI', aln, 1);
    expect(r1.coding, isTrue);
    expect(detectCoding('unknownLocus', aln, 1).stops, greaterThan(0));
    // rRNA names are never coding.
    expect(detectCoding('16S', aln, 1).coding, isFalse);
  });

  test('known protein-coding names stay coding even with stops', () {
    final aln = [for (var i = 0; i < 4; i++) SeqRecord('s$i', 'AACTCCGTATGAGATGAAATAGTATAGGTGAATAGC')];
    expect(detectCoding('ND4', aln, 2).coding, isTrue);
    expect(detectCoding('mystery', aln, 2).coding, isFalse);
  });

  test('gaps that are not multiples of 3 are reported as frame shifts', () {
    final aln = [
      SeqRecord('ok', 'ATGAAA---CCCGGGTTT'),
      SeqRecord('partial', '--GAAACCCCCCGGGTTT'),
      SeqRecord('shift', 'ATGAA-ACCCCCGGGTTT'),
    ];
    expect(frameShiftedSequences(aln), ['shift']);
  });

  test('a locus without known codon positions is one partition', () {
    final m = buildSupermatrix([
      LocusAlignment('COI', [SeqRecord('a', 'ATGAAACCC'), SeqRecord('b', 'ATGAAACCC')]),
      LocusAlignment('ND1', [SeqRecord('a', 'ATGAAACC'), SeqRecord('b', 'ATGAAACC')], coding: false),
    ]);
    expect(m.toPartitionTable(), 'DNA, COI_pos1 = 1-7\\3\nDNA, COI_pos2 = 2-8\\3\nDNA, COI_pos3 = 3-9\\3\nDNA, ND1 = 10-17');
  });
}
