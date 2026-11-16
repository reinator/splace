import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/core/newick.dart';
import 'package:splace/src/core/tree_compare.dart';

List<String> order(TreeNode t) => [for (final n in t.tips) tipKey(n.name)];

void main() {
  group('compareTrees', () {
    test('same topology, other rooting and rotation: distance 0', () {
      final a = parseNewick('((A,B),(C,(D,E)));');
      final b = parseNewick('(((E,D),C),(B,A));');
      final r = compareTrees(a, b);
      expect(r.rf, 0);
      expect(r.normalizedRf, 0);
      expect(r.conflictA, isEmpty);
      expect(r.shared, r.splitsA.length);
    });

    test('one moved tip gives the conflicting clades on both sides', () {
      final a = parseNewick('(((A,B),C),(D,E),F);');
      final b = parseNewick('(((A,C),B),(D,E),F);');
      final r = compareTrees(a, b);
      expect(r.rf, 2);
      expect(r.conflictA.single.tips.map((t) => t.name).toSet(), {'A', 'B'});
      expect(r.conflictB.single.tips.map((t) => t.name).toSet(), {'A', 'C'});
      expect(r.maxRf, 6);
    });

    test('tips in one tree only are left out of the clades', () {
      final a = parseNewick('((A,B),(C,(D,X)));');
      final b = parseNewick('((A_,B),(C,D));'.replaceFirst('A_', 'a'));
      final r = compareTrees(a, b);
      expect(r.common, {'a', 'b', 'c', 'd'});
      expect(r.onlyA, {'x'});
      expect(r.rf, 0);
    });
  });

  group('untangle', () {
    test('rotates the second tree to follow the first one', () {
      final a = parseNewick('((A,B),(C,(D,E)));');
      final b = parseNewick('(((E,D),C),(B,A));');
      expect(crossings(order(a), order(b)), greaterThan(0));
      untangle(b, {for (final (i, k) in order(a).indexed) k: i});
      expect(order(b), order(a));
      expect(crossings(order(a), order(b)), 0);
    });
  });
}
