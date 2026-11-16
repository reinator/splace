import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/core/genbank.dart';
import 'package:splace/src/core/taxonomy_sources.dart';

void main() {
  test('NCBI Taxonomy XML: lineage ranks and authority', () {
    const xml = '''<?xml version="1.0" ?>
<TaxaSet><Taxon>
    <TaxId>8384</TaxId>
    <ScientificName>Rhinella marina</ScientificName>
    <OtherNames><Name><ClassCDE>authority</ClassCDE><DispName>Rhinella marina (Linnaeus, 1758)</DispName></Name></OtherNames>
    <Rank>species</Rank>
    <LineageEx>
        <Taxon><TaxId>2759</TaxId><ScientificName>Eukaryota</ScientificName><Rank>superkingdom</Rank></Taxon>
        <Taxon><TaxId>33208</TaxId><ScientificName>Metazoa</ScientificName><Rank>kingdom</Rank></Taxon>
        <Taxon><TaxId>7711</TaxId><ScientificName>Chordata</ScientificName><Rank>phylum</Rank></Taxon>
        <Taxon><TaxId>8292</TaxId><ScientificName>Amphibia</ScientificName><Rank>class</Rank></Taxon>
        <Taxon><TaxId>8342</TaxId><ScientificName>Anura</ScientificName><Rank>order</Rank></Taxon>
        <Taxon><TaxId>8382</TaxId><ScientificName>Bufonidae</ScientificName><Rank>family</Rank></Taxon>
    </LineageEx>
</Taxon></TaxaSet>''';
    final t = parseNcbiTaxonomyXml(xml)['8384']!;
    expect([t.kingdom, t.phylum, t.klass, t.order, t.family, t.authorship], ['Metazoa', 'Chordata', 'Amphibia', 'Anura', 'Bufonidae', '(Linnaeus, 1758)']);
    expect(t.source, 'NCBI');
  });

  test('SALVE API: sheet, ranks and the whole assessment', () {
    final rows = <Map<String, dynamic>>[
      {'st_excluida': true, 'nm_cientifico': '<i>Panthera onca</i>', 'cd_categoria_final': 'CR'},
      {
        'st_excluida': false,
        'cd_situacao_ficha': 'PUBLICADA',
        'nm_cientifico_atual': '<i>Panthera onca</i>&nbsp;<span>(Linnaeus, 1758)</span>',
        'de_categoria_final_completa': 'Vulnerável (VU)',
        'cd_categoria_final': 'VU',
        'ds_criterio_aval_iucn_final': 'A2cd',
        'no_comum': 'Onça-pintada',
        'ds_grupo_salve': 'Mamíferos',
        'no_bioma': 'Amazônia, Cerrado',
        'sg_estado': 'AC, AM',
        'de_periodo_avaliacao': '01/2023',
        'ds_doi': '10.37002/salve.ficha.14015.2',
        'st_possivelmente_extinta': '',
      },
      {'st_excluida': false, 'nm_cientifico': '<i>Panthera leo</i>'},
    ];
    final row = salvePickSheet(rows, speciesKey('Panthera_onca'))!;
    final t = salveTaxonInfo('Panthera onca', row, sections: {'taxonomicClassification': {
      'tree': [
        {'name': 'Reino', 'value': 'Animalia'},
        {'name': 'Ordem', 'value': 'Carnivora'},
        {'name': 'Família', 'value': 'Felidae'},
      ],
      'oldNames': 'Felis onca',
    }, 'distribution': {
      'endemicOfBrazil': 'Não'
    }, 'header': {
      'justificative': '<p>Declínio <em>projetado</em>.</p>',
      'authorship': 'Morato, R.G.',
    }, 'population': {
      'populationTrend': 'Declinando',
      'generationalTime': '6.84 Ano(s)',
    }, 'threats': {
      'table': [
        {'type': {'description': '4 - Transportes|4.1 - Estradas e ferrovias'}},
        {'type': {'description': '2 - Agropecuária|2.1 - Culturas anuais'}},
      ],
    }, 'conservation': {
      'officialListPresence': 'Sim',
      'historyTable': [
        {'year': 2023, 'type': 'Nacional Brasil', 'category': 'Vulnerável (VU)', 'criteria': 'A2cd'},
      ],
      'UCPresenceTable': [
        {'name': 'APA Cavernas do Peruaçu'},
        {'name': 'PARNA Iguaçu'},
      ],
    }, 'bibliographicReferences': {
      'table': {'refsTaxon': [1, 2, 3], 'refsOccurrences': []},
    }});
    expect([t.kingdom, t.order, t.family, t.authorship], ['Animalia', 'Carnivora', 'Felidae', '(Linnaeus, 1758)']);
    final s = t.salve!;
    expect([s.code, s.categoryName, s.criteria, s.endemic, s.synonyms], ['VU', 'Vulnerável', 'A2cd', 'Não', 'Felis onca']);
    expect(s.justification, 'Declínio projetado.');
    expect(s.possiblyExtinct, isFalse);
    expect([s.populationTrend, s.generationTime, s.officialList], ['Declinando', '6.84 Ano(s)', 'Sim']);
    expect(s.threats, ['Estradas e ferrovias', 'Culturas anuais']);
    expect(s.history, ['2023 · Nacional Brasil · Vulnerável (VU) · A2cd']);
    expect(s.protectedAreas.length, 2);
    expect(s.references, 3);
    expect(salvePickSheet(rows, 'puma concolor'), isNull);
  });

  test('taxon id from the source feature', () {
    final r = parseGenBank('''LOCUS       X1  10 bp DNA
VERSION     X1.1
FEATURES             Location/Qualifiers
     source          1..10
                     /organism="Rhinella marina"
                     /db_xref="BOLD:ABC123"
                     /db_xref="taxon:8384"
ORIGIN
        1 acgtacgtac
//
''')!;
    expect(r.taxId, '8384');
  });
}
