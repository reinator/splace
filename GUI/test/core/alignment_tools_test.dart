import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/core/alignment_tools.dart';
import 'package:splace/src/core/fasta.dart';

void main() {
  test('trimAl column map', () {
    // Verbatim trimAl 1.5 output with -colnumbering.
    expect(parseTrimalColumnMap('#ColumnsMap\t0, 1, 2, 6, 7, 8, 9, 10, 11\n'), [0, 1, 2, 6, 7, 8, 9, 10, 11]);
    expect(parseTrimalColumnMap('nothing'), isNull);
  });

  test('ClipKIT log: kept columns (--codon trims whole codons)', () {
    // Verbatim ClipKIT 2.14 `-m gappy -g 0.5 --codon -l` output.
    const log = '1 keep constant 0.0\n2 keep constant 0.0\n3 keep constant 0.0\n4 trim parsimony-informative 0.0\n'
        '5 trim other 0.8\n6 trim other 0.8\n7 keep constant 0.0\n8 keep singleton 0.0\n9 keep constant 0.4\n';
    expect(parseClipkitLog(log), [0, 1, 2, 6, 7, 8]);
    expect(parseClipkitLog('nothing here'), isNull);
  });

  test('MAFFT reversed prefix', () {
    final r = stripMafftReversedPrefix([SeqRecord('_R_s1', 'AC'), SeqRecord('s2', 'AC')]);
    expect(r.records.map((s) => s.header), ['s1', 's2']);
    expect(r.reversed, ['s1']);
  });

  test('alignment statistics', () {
    final s = alignmentStats([
      SeqRecord('a', 'AAC-'),
      SeqRecord('b', 'AGC-'),
      SeqRecord('c', 'GGCA'),
      SeqRecord('d', 'GGTA'),
    ]);
    expect(s.variableSites, 3);
    expect(s.informativeSites, 1, reason: 'only column 0 has two states each in ≥ 2 taxa');
    expect(s.gapFraction, closeTo(2 / 16, 1e-9));
  });

  test('keptColumns recovers the columns ClipKIT kept', () {
    final full = [SeqRecord('a', 'MK-LV'), SeqRecord('b', 'MKALV'), SeqRecord('c', 'MR-LI')];
    final trimmed = [SeqRecord('a', 'MKLV'), SeqRecord('b', 'MKLV'), SeqRecord('c', 'MRLI')];
    expect(keptColumns(full, trimmed), [0, 1, 3, 4]);
    expect(keptColumns(full, [SeqRecord('a', 'X'), SeqRecord('b', 'X'), SeqRecord('c', 'X')]), isNull);
  });
}
