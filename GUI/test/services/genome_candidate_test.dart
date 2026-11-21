import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/services/ncbi.dart';

void main() {
  test('the current NCBI organism names the record; an old name in the title is kept as "filed as"', () {
    final c = GenomeCandidate.fromSummary('1', {
      'accessionversion': 'AB123456.1',
      'title': 'Bufo marinus cytochrome oxidase subunit I (COI) gene, partial cds; mitochondrial',
      'organism': 'Rhinella marina',
      'slen': 658,
    });
    expect(c.species, 'Rhinella marina');
    expect(c.filedAs, 'Bufo marinus');
    expect(c.description, startsWith('Rhinella marina cytochrome oxidase'));
    expect(c.matches('Bufo marinus'), isTrue);

    final same = GenomeCandidate.fromSummary('2', {
      'accessionversion': 'PV403566.1',
      'title': 'Rhinella marina voucher FPI0001257 cytochrome oxidase subunit 1 (COI) gene, partial cds; mitochondrial',
      'organism': 'Rhinella marina',
    });
    expect(same.species, 'Rhinella marina');
    expect(same.filedAs, isNull);

    // No organism in the summary: the name of the title.
    final old = GenomeCandidate.fromSummary('3', {'caption': 'X1', 'title': 'Rhinella marina mitochondrion, complete genome'});
    expect(old.species, 'Rhinella marina');
    expect(old.filedAs, isNull);
  });
}
