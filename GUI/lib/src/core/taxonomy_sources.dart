/// Taxonomy from GBIF, NCBI Taxonomy and the ICMBio SALVE API, reduced to
/// the ranks shown in the records table.
library;

import 'genbank.dart';

enum TaxonomySource { gbif, ncbi, salve }

class TaxonInfo {
  const TaxonInfo({
    this.kingdom = '',
    this.phylum = '',
    this.klass = '',
    this.order = '',
    this.family = '',
    this.authorship = '',
    this.status,
    this.salve,
    required this.source,
  });

  final String kingdom, phylum, klass, order, family, authorship;

  /// Conservation category (SALVE).
  final String? status;

  /// The whole SALVE sheet.
  final SalveSheet? salve;
  final String source;

  TaxonomyRanks toRanks() =>
      TaxonomyRanks(kingdom: kingdom, phylum: phylum, klass: klass, order: order, family: family, authorship: authorship, source: source);
}

/// "Rhinella marina (Linnaeus, 1758)" → "rhinella marina".
String speciesKey(String name) {
  final words = name.trim().replaceAll('_', ' ').split(RegExp(r'\s+')).where((w) => w.isNotEmpty).toList();
  return words.take(2).join(' ').toLowerCase();
}

// ---------------------------------------------------------------------------
// NCBI Taxonomy (efetch db=taxonomy, XML)
// ---------------------------------------------------------------------------

String _unescape(String s) =>
    s.replaceAll('&amp;', '&').replaceAll('&lt;', '<').replaceAll('&gt;', '>').replaceAll('&quot;', '"').replaceAll('&apos;', "'");

/// Top-level `<Taxon>` blocks of a TaxaSet (nested lineage taxa skipped).
List<String> _topLevelTaxa(String xml) {
  final out = <String>[];
  final tag = RegExp(r'<(/?)Taxon>');
  var depth = 0, start = -1;
  for (final m in tag.allMatches(xml)) {
    if (m.group(1)!.isEmpty) {
      if (depth == 0) start = m.end;
      depth++;
    } else {
      depth--;
      if (depth == 0 && start >= 0) out.add(xml.substring(start, m.start));
    }
  }
  return out;
}

/// TaxId → ranks from an NCBI Taxonomy efetch response.
Map<String, TaxonInfo> parseNcbiTaxonomyXml(String xml) {
  final out = <String, TaxonInfo>{};
  for (final block in _topLevelTaxa(xml)) {
    final id = RegExp(r'<TaxId>(\d+)</TaxId>').firstMatch(block)?.group(1);
    if (id == null) continue;
    final name = _unescape(RegExp(r'<ScientificName>([^<]*)</ScientificName>').firstMatch(block)?.group(1) ?? '');
    final ranks = <String, String>{};
    final lineage = RegExp(r'<LineageEx>([\s\S]*?)</LineageEx>').firstMatch(block)?.group(1) ?? '';
    for (final m in RegExp(r'<ScientificName>([^<]*)</ScientificName>\s*<Rank>([^<]*)</Rank>').allMatches(lineage)) {
      ranks.putIfAbsent(m.group(2)!, () => _unescape(m.group(1)!));
    }
    var authority = _unescape(
        RegExp(r'<ClassCDE>authority</ClassCDE>\s*<DispName>([^<]*)</DispName>').firstMatch(block)?.group(1) ?? '');
    if (authority.toLowerCase().startsWith(name.toLowerCase())) authority = authority.substring(name.length).trim();
    out[id] = TaxonInfo(
      // NCBI calls animals "Metazoa" (kingdom) and plants "Viridiplantae".
      kingdom: ranks['kingdom'] ?? ranks['superkingdom'] ?? ranks['domain'] ?? '',
      phylum: ranks['phylum'] ?? '',
      klass: ranks['class'] ?? '',
      order: ranks['order'] ?? '',
      family: ranks['family'] ?? '',
      authorship: authority,
      source: 'NCBI',
    );
  }
  return out;
}

// ---------------------------------------------------------------------------
// ICMBio SALVE public API (salve.icmbio.gov.br/salve-api/public)
// ---------------------------------------------------------------------------

/// SALVE returns names and texts as HTML ("<i>Puma concolor</i>&nbsp;…").
String salveText(Object? html) => '${html ?? ''}'
    .replaceAll(RegExp(r'<br\s*/?>|</p>\s*<p>', caseSensitive: false), '\n')
    // Inline tags (italics of names…) vanish; the others separate words.
    .replaceAll(RegExp(r'</?(i|em|b|strong|u|sup|sub|a)(\s[^>]*)?>', caseSensitive: false), '')
    .replaceAll(RegExp(r'<[^>]+>'), ' ')
    .replaceAll('&nbsp;', ' ')
    .replaceAll('&amp;', '&')
    .replaceAll('&lt;', '<')
    .replaceAll('&gt;', '>')
    .replaceAll('&quot;', '"')
    .replaceAll(RegExp(r'[ \t]+'), ' ')
    .replaceAll(RegExp(r' *\n *'), '\n')
    .trim();

/// The sheet of [key] (see [speciesKey]) among the rows of a `search`
/// response: the exact name (current, as assessed or former), the
/// published sheet first; excluded sheets are left out.
Map<String, dynamic>? salvePickSheet(List<Map<String, dynamic>> rows, String key) {
  Map<String, dynamic>? best;
  for (final h in rows) {
    if (h['st_excluida'] == true) continue;
    final names = [h['nm_cientifico_atual'], h['nm_cientifico'], h['nm_cientifico_anterior']].whereType<String>().map(salveText);
    if (!names.any((n) => speciesKey(n) == key)) continue;
    if (best == null || (h['cd_situacao_ficha'] == 'PUBLICADA' && best['cd_situacao_ficha'] != 'PUBLICADA')) best = h;
  }
  return best;
}

/// Sections of a SALVE sheet (`fichaHtml?section=`).
const salveSections = [
  'header', 'taxonomicClassification', 'distribution', 'naturalHistory', 'population', 'threats', 'uses', 'conservation', 'research',
  'bibliographicReferences',
];

/// Ranks and sheet from the `search` row and the `fichaHtml` [sections]
/// (by name, see [salveSections]; missing ones are skipped).
TaxonInfo salveTaxonInfo(String species, Map<String, dynamic> row, {Map<String, Map<String, dynamic>?> sections = const {}}) {
  Map<String, dynamic>? sec(String k) => sections[k];
  final taxonomy = sec('taxonomicClassification'), distribution = sec('distribution'), header = sec('header');
  final natural = sec('naturalHistory'), population = sec('population'), threats = sec('threats'), uses = sec('uses');
  final conservation = sec('conservation'), research = sec('research'), refs = sec('bibliographicReferences');
  String s(Map<String, dynamic>? m, String k) => m == null ? '' : salveText(m[k]);
  Map<String, dynamic>? map(Object? v) => v is Map<String, dynamic> ? v : null;
  List<Map<String, dynamic>> table(Object? v) {
    final t = v is Map ? v['table'] : v;
    return t is List ? t.whereType<Map<String, dynamic>>().toList() : const [];
  }

  List<String> unique(Iterable<String> xs) => xs.map((x) => x.trim()).where((x) => x.isNotEmpty).toSet().toList();
  final ranks = <String, String>{};
  final tree = taxonomy?['tree'];
  if (tree is List) {
    for (final n in tree.whereType<Map<String, dynamic>>()) {
      ranks[salveText(n['name']).toLowerCase()] = salveText(n['value']);
    }
  }
  final rawName = '${row['nm_cientifico_atual'] ?? row['nm_cientifico'] ?? ''}';
  final author = RegExp(r'<span>(.*?)</span>').firstMatch(rawName)?.group(1);
  final category = s(row, 'de_categoria_final_completa').isNotEmpty ? s(row, 'de_categoria_final_completa') : s(header, 'category');
  var code = s(row, 'cd_categoria_final').toUpperCase();
  if (code.isEmpty) code = RegExp(r'\(([A-Z]{2})\)').firstMatch(category)?.group(1) ?? '';
  final doi = s(row, 'ds_doi').isNotEmpty ? s(row, 'ds_doi') : s(header, 'DOI');

  // Natural history.
  final food = map(natural?['foodHabit']);
  final repro = map(natural?['reproduction']);
  final reproParts = [
    if (s(repro, 'gestationTime').isNotEmpty) 'gestação ${s(repro, 'gestationTime')}',
    if (s(repro, 'offspringSize').isNotEmpty) 'prole ${s(repro, 'offspringSize')}',
    if (s(repro, 'birthInterval').isNotEmpty) 'intervalo ${s(repro, 'birthInterval')}',
    for (final r in table(repro))
      if (s(r, 'male').isNotEmpty || s(r, 'female').isNotEmpty)
        '${s(r, 'name')}: ${[if (s(r, 'male').isNotEmpty) '♂ ${s(r, 'male')}', if (s(r, 'female').isNotEmpty) '♀ ${s(r, 'female')}'].join(', ')}',
  ];
  final interactions = map(natural?['interactionsWithOtherSpecies']);
  // Threats: "2 - Agropecuária|2.1 - Culturas…" → the most specific level.
  String lastLevel(String v) => v.split('|').last.replaceFirst(RegExp(r'^[\d.]+\s*-\s*'), '').trim();
  final threatList = unique(table(threats).map((t) => lastLevel(salveText(map(t['type'])?['description'] ?? t['type']))));
  final useList = unique(table(uses).map((u) => salveText(u['type']).replaceFirst(RegExp(r'^[\d.]+\s*-\s*'), '')));
  final history = [
    for (final h in table(conservation?['historyTable']))
      [s(h, 'year'), s(h, 'type'), s(h, 'category'), s(h, 'criteria')].where((x) => x.isNotEmpty).join(' · '),
  ];
  final actions = [
    for (final a in table(conservation?['conservationActionsTable']))
      [s(a, 'actionPanName').isNotEmpty ? s(a, 'actionPanName') : s(a, 'name'), s(a, 'situation')].where((x) => x.isNotEmpty).join(' — '),
  ];
  final officialLists = unique([
    for (final l in (conservation?['officialLists'] is List ? conservation!['officialLists'] as List : const []).whereType<Map<String, dynamic>>())
      s(l, 'de_titulo_formatado'),
    for (final l in table(conservation?['conventionPresenceTable'])) s(l, 'name'),
  ]);
  final ucs = unique(table(conservation?['UCPresenceTable']).map((u) => s(u, 'name')));
  final researchNeeds = unique(table(research).map((r) => [for (final v in r.values) if (v is String) salveText(v)].where((x) => x.isNotEmpty).join(' · ')));
  final refsTable = map(refs?['table']);
  final references = ((refsTable?['refsTaxon'] as List?)?.length ?? 0) + ((refsTable?['refsOccurrences'] as List?)?.length ?? 0);
  final texts = <String, String>{
    'naturalHistory': s(natural, 'description'),
    'habitat': s(natural, 'habitatObservations'),
    'food': s(food, 'observations'),
    'reproduction': s(repro, 'observations'),
    'interactions': s(interactions, 'observations'),
    'population': s(population, 'populationObservations'),
    'genetics': s(population, 'geneticCharacteristics'),
    'threats': s(threats, 'description'),
    'uses': s(uses, 'description'),
    'conservation': s(conservation, 'description'),
    'research': s(research, 'description'),
    'globalDistribution': s(distribution, 'globalDistribution'),
    'nationalDistribution': s(distribution, 'nationalDistribution'),
  }..removeWhere((_, v) => v.isEmpty);

  final sheet = SalveSheet(
    category: category.isEmpty ? code : category,
    code: code,
    criteria: s(row, 'ds_criterio_aval_iucn_final'),
    commonNames: s(row, 'no_comum').isNotEmpty ? s(row, 'no_comum') : s(taxonomy, 'commonNames'),
    group: s(row, 'ds_grupo_salve'),
    biomes: s(row, 'no_bioma').isNotEmpty ? s(row, 'no_bioma') : s(distribution, 'biomes'),
    states: s(row, 'sg_estado').isNotEmpty ? s(row, 'sg_estado') : s(distribution, 'states'),
    watersheds: s(distribution, 'watersheds'),
    endemic: s(distribution, 'endemicOfBrazil'),
    possiblyExtinct: const {'S', 'SIM', 'TRUE'}.contains(s(row, 'st_possivelmente_extinta').toUpperCase()),
    assessed: s(row, 'de_periodo_avaliacao').isNotEmpty ? s(row, 'de_periodo_avaliacao') : s(header, 'lastAvaliation'),
    published: s(header, 'publicationDate'),
    situation: s(row, 'cd_situacao_ficha'),
    doi: doi,
    synonyms: s(taxonomy, 'oldNames'),
    assessors: s(header, 'authorship'),
    citation: s(header, 'citation'),
    justification: s(header, 'justificative'),
    url: 'https://salve.icmbio.gov.br/#/?q=${Uri.encodeComponent(speciesKey(species))}',
    populationTrend: s(population, 'populationTrend'),
    generationTime: s(population, 'generationalTime'),
    sexRatio: s(population, 'sexRatio'),
    migratory: s(natural, 'migratorySpecie'),
    habitats: unique(table(natural?['habitats']).map((h) => s(h, 'habitat').replaceFirst(RegExp(r'^[A-Z0-9.]+\s*-\s*'), ''))),
    foodHabits: unique(table(food).map((f) => s(f, 'type'))),
    reproduction: reproParts.join('; '),
    interactions: unique(table(interactions).map((i) => [s(i, 'taxon'), s(i, 'type'), s(i, 'category')].where((x) => x.isNotEmpty).join(' · '))),
    threats: threatList,
    uses: useList,
    officialList: s(conservation, 'officialListPresence'),
    officialLists: officialLists,
    history: history,
    actions: actions,
    protectedAreas: ucs,
    researchNeeds: researchNeeds,
    references: references,
    texts: texts,
  );
  return TaxonInfo(
    kingdom: ranks['reino'] ?? '',
    phylum: ranks['filo'] ?? '',
    klass: ranks['classe'] ?? '',
    order: ranks['ordem'] ?? '',
    family: ranks['família'] ?? ranks['familia'] ?? '',
    authorship: author == null ? '' : salveText(author),
    status: sheet.category.isEmpty ? null : sheet.category,
    salve: sheet,
    source: 'SALVE',
  );
}
