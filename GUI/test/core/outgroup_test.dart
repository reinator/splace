import 'dart:math' as math;

import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/core/fasta.dart';
import 'package:splace/src/core/outgroup.dart';
import 'package:splace/src/core/supermatrix.dart';

void main() {
  // An ingroup of 8 close taxa and two outgroup taxa far from them (and
  // close to each other), in two loci; one taxon is missing a locus.
  final rnd = math.Random(3);
  String mutate(String s, double rate) => String.fromCharCodes(s.codeUnits.map((c) {
        if (rnd.nextDouble() >= rate) return c;
        final others = 'ACGT'.codeUnits.where((b) => b != c).toList();
        return others[rnd.nextInt(3)];
      }));
  List<LocusAlignment> data() {
    final loci = <LocusAlignment>[];
    for (final name in ['COI', 'CYTB']) {
      final root = List.generate(300, (_) => 'ACGT'[rnd.nextInt(4)]).join();
      final ingroup = mutate(root, 0.02), out = mutate(root, 0.25);
      loci.add(LocusAlignment(name, [
        for (var i = 0; i < 8; i++)
          if (!(name == 'CYTB' && i == 3)) SeqRecord('In_$i', mutate(ingroup, 0.02)),
        SeqRecord('Out_a', mutate(out, 0.02)),
        SeqRecord('Out_b', mutate(out, 0.02)),
      ]));
    }
    return loci;
  }

  test('the group farthest from the rest is suggested as outgroup', () {
    final s = suggestOutgroup(data())!;
    expect(s.confident, isTrue);
    expect(s.taxa.toSet(), {'Out_a', 'Out_b'});
    expect(s.between, greaterThan(3 * s.within));
    expect(s.ranked.take(2).toSet(), {'Out_a', 'Out_b'});
    expect(s.taxonCount, 10);
    expect(s.loci, 2);
  });

  test('gaps are not compared and too few taxa give no suggestion', () {
    expect(suggestOutgroup([LocusAlignment('x', [SeqRecord('a', 'ACGT'), SeqRecord('b', 'ACGA')])]), isNull);
    final s = suggestOutgroup([
      LocusAlignment('x', [
        SeqRecord('a', 'A' * 40),
        SeqRecord('b', 'A' * 40),
        SeqRecord('c', 'A' * 39 + 'C'),
        SeqRecord('d', '-' * 20 + 'C' * 20),
      ]),
    ])!;
    // d differs from the others only where it has bases.
    expect(s.meanDistance['d'], closeTo(1.0, 0.06));
    expect(s.taxa, ['d']);
  });
}
