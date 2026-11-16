import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/core/fasta.dart';
import 'package:splace/src/core/mafft_progress.dart';
import 'package:splace/src/core/genbank.dart';
import 'package:splace/src/services/species_id.dart';

void main() {
  test('MAFFT progress follows the stages and never goes back', () {
    final p = MafftProgress(46);
    final values = <double>[];
    final lines = [
      'All-to-all alignment.',
      for (var i = 0; i < 46; i += 10) '   $i / 46 (by thread   ${i % 2}) ',
      'Constructing a UPGMA tree ... ',
      for (var i = 0; i < 46; i += 10) '   $i / 46',
      'Progressive alignment ... ',
      for (final i in [1, 2, 4, 3, 20, 45]) 'STEP    $i /45 (thread    0) ',
      'Segment   1/  1    1-1552',
      for (var j = 0; j < 90; j += 10) '001-${j.toString().padLeft(4, '0')}-1 (thread    2) identical     ',
      '002-0010-0 (thread    1) identical     ',
    ];
    for (final l in lines) {
      p.add(l);
      values.add(p.value);
    }
    for (var i = 1; i < values.length; i++) {
      expect(values[i], greaterThanOrEqualTo(values[i - 1]));
    }
    // Before the progressive alignment: under 35 %.
    expect(values[lines.indexOf('Progressive alignment ... ')], lessThan(0.35));
    // End of the progressive alignment: 70 %.
    expect(values[lines.indexWhere((l) => l.startsWith('STEP    45'))], closeTo(0.70, 0.001));
    expect(p.value, greaterThan(0.84));
    expect(p.value, lessThanOrEqualTo(0.98));
    expect(MafftProgress.isProgressLine('STEP    3 /45 (thread    0)'), isTrue);
    expect(MafftProgress.isProgressLine('Strategy:'), isFalse);
  });

  test('FFT-NS-2 second progressive alignment goes past 70 %', () {
    final p = MafftProgress(10);
    for (final l in ['Progressive alignment 1/2... ', 'STEP     9 / 9 ', '    5 / 10', 'Progressive alignment 2/2... ', 'STEP     5 / 9 ']) {
      p.add(l);
    }
    expect(p.value, greaterThan(0.70));
  });

  test('fingerprint changes with the sequences', () {
    final a = [SeqRecord('a', 'ACGT'), SeqRecord('b', 'ACGA')];
    expect(fingerprint(a), fingerprint([SeqRecord('a', 'ACGT'), SeqRecord('b', 'ACGA')]));
    expect(fingerprint(a), isNot(fingerprint([SeqRecord('a', 'ACGT'), SeqRecord('b', 'ACGG')])));
  });

  test('species identification reuses cached answers without sending anything', () async {
    final seq = 'ACGT' * 40;
    final records = [
      GbRecord(accession: 'NC_1.1', organism: 'Caranx ruber', lineage: '', definition: '', features: const [], sequence: seq),
    ];
    final cache = {
      SpeciesIdentifier.cacheKey(IdSource.genbank, 'NC_1.1', 'COI', seq): [('Caranx bartholomaei', 99.5, 'X1.1')],
      SpeciesIdentifier.cacheKey(IdSource.bold, 'NC_1.1', 'COI', seq): [('Caranx bartholomaei', 100.0, 'P1')],
    };
    final id = SpeciesIdentifier();
    final out = await id.identify(records, {
      'NC_1.1': {'COI': seq},
    }, cache: cache);
    id.close();
    expect(id.reused, 2);
    expect(id.sent, 0);
    final r = out['NC_1.1']!;
    expect(r.candidates.first.species, 'Caranx bartholomaei');
    expect(r.candidates.first.bold?.identity, 100.0);
  });
}
