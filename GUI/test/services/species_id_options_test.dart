import 'package:flutter/material.dart';
import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/core/pairwise.dart';
import 'package:splace/src/services/species_id.dart';
import 'package:splace/src/ui/widgets/html_text.dart';

void main() {
  test('pairwise alignment puts the gap where the base is missing', () {
    final (a, b) = alignPair('ACGTTGCAACGT', 'ACGTGCAACGT');
    expect(a.length, b.length);
    expect(a.replaceAll('-', ''), 'ACGTTGCAACGT');
    expect(b.replaceAll('-', ''), 'ACGTGCAACGT');
    expect('-'.allMatches(b).length, 1);
  });

  test('hits keep their alignment through the cache; parameters change the key', () {
    final h = IdHit(
      source: IdSource.genbank,
      species: 'Rhinella marina',
      accession: 'MN1.1',
      identity: 99.2,
      alignLength: 650,
      qFrom: 1,
      qTo: 650,
      sFrom: 20,
      sTo: 669,
      mismatches: 5,
      evalue: 0,
      bitScore: 1180,
      qseq: 'ACGT',
      hseq: 'ACGA',
    );
    final back = IdHit.fromJson(h.toJson());
    expect(back.hasAlignment, isTrue);
    expect(back.sFrom, 20);
    expect(IdHit.fromJson(h.toJson(withAlignment: false)).hasAlignment, isFalse);

    const a = SpeciesIdOptions();
    final b = a.copyWith(hitlistSize: 100);
    final k1 = SpeciesIdentifier.cacheKey(IdSource.genbank, 'X', 'COI', 'ACGT', a);
    expect(k1, isNot(SpeciesIdentifier.cacheKey(IdSource.genbank, 'X', 'COI', 'ACGT', b)));
    // Coverage is applied after the cache: the same key.
    expect(k1, SpeciesIdentifier.cacheKey(IdSource.genbank, 'X', 'COI', 'ACGT', a.copyWith(minCoverage: 0.5)));
    expect(SpeciesIdentifier.staleCacheKey('genbank2|X|COI|f'), isTrue);
    expect(SpeciesIdentifier.staleCacheKey(k1), isFalse);
    expect(SpeciesIdOptions.fromJson(b.toJson()).hitlistSize, 100);
  });

  test('BOLD identifies the markers of its database', () {
    const o = SpeciesIdOptions();
    expect(o.boldAccepts('COI'), isTrue);
    expect(o.boldAccepts('rbcL'), isFalse);
    final plants = o.copyWith(boldDatabase: 'public.plants');
    expect(plants.boldAccepts('rbcL'), isTrue);
    expect(plants.boldAccepts('matK'), isTrue);
    expect(plants.boldAccepts('COI'), isFalse);
  });

  test('tooltip HTML: bold, italic and line breaks', () {
    final span = htmlSpan('A <b>bold</b> and <i>italic</i> name &amp; more<br>next') as TextSpan;
    final parts = span.children!.cast<TextSpan>();
    expect(parts.map((p) => p.text).join(), 'A bold and italic name & more\nnext');
    expect(parts.firstWhere((p) => p.text == 'bold').style!.fontWeight, FontWeight.w700);
    expect(parts.firstWhere((p) => p.text == 'italic').style!.fontStyle, FontStyle.italic);
  });
}
