import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/core/newick.dart';
import 'package:splace/src/core/tool_commands.dart';

void main() {
  group('multilineCommand', () {
    test('one option per line with its values, ended by a backslash', () {
      final s = multilineCommand('iqtree3', ['-s', 'concatenated.phy', '-p', 'partitions.nex', '-m', 'MFP+MERGE', '-B', '1000', '-redo']);
      expect(s.split('\n'), [
        r'iqtree3 \',
        r'  -s concatenated.phy \',
        r'  -p partitions.nex \',
        r'  -m MFP+MERGE \',
        r'  -B 1000 \',
        '  -redo',
      ]);
    });

    test('quotes paths with spaces and keeps negative numbers as values', () {
      final s = multilineCommand('raxml-ng', ['--msa', 'My Data/a.phy', '--brlen', '-1']);
      expect(s, "raxml-ng \\\n  --msa 'My Data/a.phy' \\\n  --brlen -1");
    });
  });

  group('TipInfo.guess', () {
    test('finds species, genus, family and accession in a SPLACE header', () {
      final i = TipInfo.guess('NC_012345_Ranidae_Rana_temporaria');
      expect(i.species, 'Rana temporaria');
      expect(i.genus, 'Rana');
      expect(i.family, 'Ranidae');
      expect(i.accession, 'NC_012345');
    });

    test('leaves out what the name does not have', () {
      final i = TipInfo.guess('sample_17');
      expect(i.species, '');
      expect(i.accession, '');
    });
  });

  group('pruneNode', () {
    test('removes a tip and joins the lengths of its parent and sibling', () {
      final t = parseNewick('((A:1,B:2):3,(C:1,D:1):1);');
      final b = t.tips.firstWhere((n) => n.name == 'B');
      final r = pruneNode(t, b);
      expect(toNewick(r), '(A:4,(C:1,D:1):1);');
    });

    test('a root left with one child gives its place to that child', () {
      final t = parseNewick('((A:1,B:1):1,C:2);');
      final r = pruneNode(t, t.tips.firstWhere((n) => n.name == 'C'));
      expect(r.tips.map((n) => n.name), ['A', 'B']);
      expect(r.parent, isNull);
    });

    test('keeps at least two tips', () {
      final t = parseNewick('(A:1,B:1);');
      expect(identical(pruneNode(t, t.children.first), t), isTrue);
      expect(t.tipCount, 2);
    });
  });
}
