import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/core/taxonomy_sources.dart';
import 'package:splace/src/services/species_id.dart';

void main() {
  test('hits under a synonym count for the accepted species', () {
    final r = RecordIdentification('X1', 'Rhinella marina')
      ..addEvidence(IdSource.genbank, 'COI', 'Bufo marinus', 99.6)
      ..addEvidence(IdSource.bold, 'COI', 'Chaunus marinus', 99.8)
      ..addEvidence(IdSource.genbank, 'COI', 'Rhinella marina', 99.1)
      ..addEvidence(IdSource.genbank, 'COI', 'Rhinella horribilis', 97.0);
    // Before: the best name is a synonym and the record looks misidentified.
    expect(r.candidates.first.species, 'Chaunus marinus');
    expect(r.disagrees, isTrue);
    r.mergeSynonyms({'Bufo marinus': 'Rhinella marina', 'Chaunus marinus': 'Rhinella marina', 'Rhinella marina': 'Rhinella marina'});
    final best = r.candidates.first;
    expect(best.species, 'Rhinella marina');
    expect(best.foundAs, {'Bufo marinus', 'Chaunus marinus'});
    expect(best.bestIdentity, 99.8);
    expect(best.genbank['COI']!.hits, 2);
    expect(best.sourceCount, 2, reason: 'GenBank and BOLD now back the same species');
    expect(r.disagrees, isFalse);
    expect(r.current, same(best));
  });

  test('a record named with a synonym is confirmed by its accepted name', () {
    final r = RecordIdentification('X2', 'Bufo marinus')..addEvidence(IdSource.genbank, 'COI', 'Rhinella marina', 99.9);
    expect(r.disagrees, isTrue);
    r.mergeSynonyms({'Bufo marinus': 'Rhinella marina', 'Rhinella marina': 'Rhinella marina'});
    expect(r.currentAccepted, 'Rhinella marina');
    expect(r.isCurrent('Rhinella marina'), isTrue);
    expect(r.disagrees, isFalse);
  });

  test('names identified to species', () {
    expect(isDeterminateSpecies('Rhinella marina'), isTrue);
    expect(isDeterminateSpecies('Rhinella sp.'), isFalse);
    expect(isDeterminateSpecies('Rhinella cf. marina'), isFalse);
    expect(isDeterminateSpecies('Rhinella'), isFalse);
  });
}
