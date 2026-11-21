import 'dart:io';

import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/core/bold.dart';
import 'package:splace/src/services/bold_portal.dart';

/// The live BOLD v5 portal API (SPLACE_NETWORK=1).
void main() {
  test('queries are split under the 250-character limit of the API', () {
    final names = List.generate(25, (i) => 'Bufo marinus subspecies$i');
    final batches = queryBatches(names);
    expect(batches.length, greaterThan(1));
    expect(batches.every((b) => b.length <= BoldPortal.maxQuery), isTrue);
    expect(batches.expand((b) => b.split(';')).toList(), names);
  });

  final network = Platform.environment['SPLACE_NETWORK'] == '1';

  test('species, genus and family are resolved; COI records come back', () async {
    final bold = BoldPortal();
    // More names than one query holds (the GBIF synonyms of R. marina).
    final terms = await bold.taxonTerms([
      'Rhinella marina', 'Bombinator horridus', 'Bombinator maculatus', 'Bufo albicans', 'Bufo brasiliensis', 'Bufo horridus', 'Bufo humeralis',
      'Bufo lazarus', 'Bufo maculiventris', 'Bufo marinis', 'Bufo marinus', 'Bufo marinus fluminensis', 'Bufo marinus horribilis', 'Bufo marinus marinus',
      'Bufo marinus napensis', 'Bufo pithecodactylus', 'Chaunus marinus', 'Docidophryne lazarus', 'Rana brasiliana', 'Rana marina', 'Bufonidae',
    ]);
    expect(terms, containsAll(['tax:species:Rhinella marina', 'tax:family:Bufonidae']));
    expect(terms.any((t) => t.contains('Bufo marinus')), isFalse);
    final counts = await bold.markerCounts(['tax:species:Rhinella marina']);
    expect(counts['COI-5P'], greaterThan(0));
    final rows = await bold.records(['tax:species:Rhinella marina'], {'COI-5P'});
    expect(rows, isNotEmpty);
    expect(rows.every((r) => r.markerCode == 'COI-5P' && r.sequence.length > 400), isTrue);
  }, skip: network ? false : 'set SPLACE_NETWORK=1', timeout: const Timeout(Duration(minutes: 3)));

  test('chloroplast barcodes: matK and rbcL of a genus become in-frame CDS', () async {
    final bold = BoldPortal();
    final terms = await bold.taxonTerms(['Quercus robur']);
    final rows = await bold.records(terms, {...boldMarkerCodes['matK']!, ...boldMarkerCodes['rbcL']!});
    expect(rows.map((r) => r.gene).toSet(), containsAll(['matK', 'rbcL']));
    final rec = boldToGbRecord(rows.firstWhere((r) => r.gene == 'rbcL'));
    final cds = rec.features.last;
    expect(cds.q('transl_table'), '11');
    expect(rec.features.first.q('organelle'), 'chloroplast');
  }, skip: network ? false : 'set SPLACE_NETWORK=1', timeout: const Timeout(Duration(minutes: 3)));
}
