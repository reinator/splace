@Tags(['network'])
library;

import 'dart:io';

import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/core/genbank_markers.dart';
import 'package:splace/src/core/marker_search.dart';
import 'package:splace/src/core/syngenes.dart';
import 'package:splace/src/services/ncbi.dart';
import 'package:splace/src/ui/chrome/genome_review_dialog.dart' show selectOnePerSpecies;

/// Live NCBI/GBIF requests. Opt-in: SPLACE_NETWORK=1 flutter test test/integration/ncbi_test.dart
void main() {
  final enabled = Platform.environment['SPLACE_NETWORK'] == '1';
  late Directory cache;
  late NcbiClient ncbi;

  setUp(() {
    cache = Directory.systemTemp.createTempSync('splace_ncbi');
    ncbi = NcbiClient(cacheDir: cache.path);
  });
  tearDown(() {
    ncbi.close();
    cache.deleteSync(recursive: true);
  });

  test('taxonomy search, review summaries and import', () async {
    const search = GenomeSearch(taxon: 'Bufonidae', completeOnly: true, refseqOnly: true);
    final ids = await ncbi.searchGenomes(search, retmax: 4);
    expect(ids, isNotEmpty);
    final candidates = await ncbi.summaries(ids);
    expect(candidates.first.accession, startsWith('NC_'));
    expect(candidates.first.classification, GenomeClass.complete);
    expect(candidates.first.species, isNotEmpty);

    final r = await ncbi.fetchGenBank(ids.take(2).toList());
    expect(r.records, hasLength(2));
    final cat = buildGeneCatalog(r.records);
    expect(cat.organelle, OrganelleType.mt);
    expect(cat.visible({'CDS'}).keys, containsAll(mtDefaultGenes));
  }, skip: !enabled);

  test('individual markers: one query per marker, joined by species', () async {
    const search = MarkerSearch(taxon: 'Astyanax lacustris', genes: ['COI', 'CYTB']);
    final byId = <String, Set<String>>{};
    for (final g in search.genes) {
      final ids = await ncbi.searchTerm(search.termFor(g), retmax: 400);
      expect(ids, isNotEmpty, reason: g);
      for (final id in ids) {
        (byId[id] ??= {}).add(g);
      }
    }
    final candidates = await ncbi.summaries(byId.keys.toList());
    for (final c in candidates) {
      c.markers.addAll(byId[c.id]!);
    }
    expect(candidates.every((c) => (c.length ?? 0) <= 10000), isTrue, reason: 'no genomes');
    selectOnePerSpecies(candidates);
    final selected = candidates.where((c) => c.selected).toList();
    expect(selected.length, lessThanOrEqualTo(4), reason: 'one per species and marker');
    final r = await ncbi.fetchGenBank([for (final c in selected) c.id]);
    final cat = buildGeneCatalog(r.records);
    final genes = extractGeneFastas(r.records, cat, {'CDS'}, {'COI', 'CYTB'}, taxonHeader, const {});
    final joined = onePerTaxon(genes).fastas;
    expect(joined.keys, containsAll(['COI', 'CYTB']));
    expect(joined['COI']!.map((s) => s.header).toSet().length, joined['COI']!.length);
    expect(joined['COI']!.first.header, contains('Astyanax_lacustris'));
  }, skip: !enabled);

  test('accessions are cached on disk', () async {
    final firstLogs = <String>[];
    final first = await ncbi.fetchGenBank(['NC_024026.1'], log: firstLogs.add);
    expect(first.records.map((r) => r.accession), ['NC_024026.1'], reason: firstLogs.join('\n'));
    final logs = <String>[];
    final again = await ncbi.fetchGenBank(['NC_024026.1'], log: logs.add);
    expect(again.records.single.accession, 'NC_024026.1');
    expect(logs, [contains('local cache')]);
  }, skip: !enabled);

  test('NCBI Taxonomy by id and by name', () async {
    final id = await ncbi.ncbiTaxId('Rhinella marina');
    expect(id, isNotNull);
    final t = (await ncbi.ncbiTaxonomy([id!]))[id]!;
    expect((t.klass, t.order, t.family), ('Amphibia', 'Anura', 'Bufonidae'));
  }, skip: !enabled);

  test('GBIF taxonomy', () async {
    final t = await ncbi.gbif('Rhinella marina');
    expect(t!.family, 'Bufonidae');
    expect(t.order, 'Anura');
  }, skip: !enabled);
}
