import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/core/genbank.dart';
import 'package:splace/src/services/species_id.dart';

GbRecord _rec(String acc, String organism, {List<String> sources = const []}) => GbRecord(
      accession: acc,
      organism: organism,
      lineage: '',
      definition: '',
      features: const [],
      sequence: '',
      sourceAccessions: sources,
    );

void main() {
  final seq = 'ACGT' * 200;

  /// Runs [SpeciesIdentifier.identify] from the cache only (no network).
  Future<RecordIdentification> identify(GbRecord r, List<(String, double, String)> genbank, {List<(String, double, String)> bold = const []}) async {
    List<IdHit> hits(IdSource src, List<(String, double, String)> l) =>
        [for (final (sp, id, acc) in l) IdHit(source: src, species: sp, accession: acc, identity: id, alignLength: seq.length)];
    final cache = {
      SpeciesIdentifier.cacheKey(IdSource.genbank, r.accession, 'COI', seq): hits(IdSource.genbank, genbank),
      SpeciesIdentifier.cacheKey(IdSource.bold, r.accession, 'COI', seq): hits(IdSource.bold, bold),
    };
    final id = SpeciesIdentifier();
    final out = await id.identify([r], {
      r.accession: {'COI': seq},
    }, cache: cache);
    id.close();
    expect(id.sent, 0);
    return out[r.accession]!;
  }

  test('the record and its RefSeq source are shown, not counted', () async {
    // BLAST merges identical sequences: NC_025644 and KM201334 (the record
    // it was copied from) come as one hit.
    final r = await identify(_rec('NC_025644.1', 'Carangoides equula', sources: ['KM201334']), [
      ('Carangoides equula', 100, 'NC_025644'),
      ('Carangoides equula', 100, 'KM201334'),
      ('Pseudocaranx georgianus', 92.4, 'OZ556862'),
    ]);
    final cur = r.current!;
    expect(cur.inGenbank, isFalse);
    expect(cur.genbank['COI']!.onlySelf, isTrue);
    expect(cur.genbank['COI']!.selfIdentity, 100);
    expect(r.selfAccessions, ['KM201334', 'NC_025644']);
    expect(r.onlySelfInGenbank, isTrue);
    // The closest other species is only a relative: the name stands.
    expect(r.disagrees, isFalse);
    expect(r.shortlist().map((c) => c.species), contains('Carangoides equula'));
  });

  test('other records of the species confirm the name', () async {
    final r = await identify(_rec('NC_088027.1', 'Selene vomer'), [
      ('Selene vomer', 100, 'NC_088027'),
      ('Selene vomer', 100, 'PX069997'),
      ('Selene brevoortii', 94.9, 'PZ310832'),
    ]);
    expect(r.current!.inGenbank, isTrue);
    expect(r.current!.genbank['COI']!.hits, 1);
    expect(r.onlySelfInGenbank, isFalse);
    expect(r.disagrees, isFalse);
  });

  test('a same-species match to another name is a possible misidentification', () async {
    final r = await identify(_rec('NC_000001.1', 'Caranx latus'), [
      ('Caranx hippos', 99.9, 'PZ428442'),
      ('Caranx latus', 92.8, 'PP032957'),
    ]);
    expect(r.disagrees, isTrue);
    expect(r.candidates.first.species, 'Caranx hippos');
  });

  test('old GenBank cache keys are stale', () {
    expect(SpeciesIdentifier.staleCacheKey('genbank|NC_1.1|COI|1:abc'), isTrue);
    expect(SpeciesIdentifier.staleCacheKey(SpeciesIdentifier.cacheKey(IdSource.genbank, 'NC_1.1', 'COI', 'ACGT')), isFalse);
    expect(SpeciesIdentifier.staleCacheKey(SpeciesIdentifier.cacheKey(IdSource.bold, 'NC_1.1', 'COI', 'ACGT')), isFalse);
  });

  test('RefSeq COMMENT gives the source accessions', () {
    final r = parseGenBank('''
LOCUS       NC_025644              16580 bp    DNA     circular VRT 05-SEP-2014
DEFINITION  Carangoides equula mitochondrion, complete genome.
ACCESSION   NC_025644
VERSION     NC_025644.1
COMMENT     REVIEWED REFSEQ: This record has been curated by NCBI staff. The
            reference sequence is identical to KM201334.
FEATURES             Location/Qualifiers
ORIGIN
        1 acgtacgt
//
''')!;
    expect(r.accession, 'NC_025644.1');
    expect(r.sourceAccessions, ['KM201334']);
  });
}
