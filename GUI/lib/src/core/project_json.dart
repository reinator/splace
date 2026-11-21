/// JSON of the analysis settings and results kept in a project file
/// (`.splc`): every option of every step, so a project can be reopened where
/// it stopped and its parameters read by whoever receives it.
library;

import 'delimitation.dart';
import 'fasta.dart';
import 'genbank.dart';
import 'genbank_markers.dart';
import 'locus_organizer.dart';
import 'marker_search.dart';
import 'saturation.dart';
import 'supermatrix.dart';
import 'syngenes.dart' show OrganelleType;
import 'tool_commands.dart';

typedef Json = Map<String, dynamic>;

/// [name] as a value of [values], or [fallback].
T enumOf<T extends Enum>(List<T> values, Object? name, T fallback) => values.where((v) => v.name == name).firstOrNull ?? fallback;

Json asJson(Object? o) => o is Map ? o.cast<String, dynamic>() : <String, dynamic>{};
List<String> strings(Object? o) => [...?(o as List?)?.map((e) => '$e')];
double? dbl(Object? o) => (o as num?)?.toDouble();
int? integer(Object? o) => (o as num?)?.toInt();

// --- sequences ------------------------------------------------------------

List<List<String>> seqsToJson(List<SeqRecord> s) => [for (final r in s) [r.header, r.sequence]];
List<SeqRecord> seqsFromJson(Object? o) => [for (final r in (o as List?) ?? const []) SeqRecord('${(r as List)[0]}', '${r[1]}')];

// --- step 1 ---------------------------------------------------------------

Json markerSearchToJson(MarkerSearch m) => {
      'taxon': m.taxon,
      'genes': m.genes,
      'organelle': m.organelle.name,
      'includeGenomes': m.includeGenomes,
      'refseqOnly': m.refseqOnly,
      'maxLength': m.maxLength,
      'synonyms': m.synonyms,
      'bold': m.bold,
    };

MarkerSearch markerSearchFromJson(Json j) => MarkerSearch(
      taxon: j['taxon'] as String? ?? '',
      genes: strings(j['genes']),
      organelle: enumOf(OrganelleType.values, j['organelle'], OrganelleType.mt),
      includeGenomes: j['includeGenomes'] == true,
      refseqOnly: j['refseqOnly'] == true,
      maxLength: integer(j['maxLength']) ?? 10000,
      synonyms: strings(j['synonyms']),
      bold: j['bold'] == true,
    );

Json fieldRuleToJson(FieldRule r) => {
      'source': r.source.name,
      'start': r.start,
      if (r.end != null) 'end': r.end,
      'delimiters': r.delimiters,
      'tag': r.tag,
      'pattern': r.pattern,
      'regexOnFileName': r.regexOnFileName,
    };

FieldRule fieldRuleFromJson(Json j) => FieldRule(
      source: enumOf(FieldSource.values, j['source'], FieldSource.header),
      start: integer(j['start']) ?? 0,
      end: integer(j['end']),
      delimiters: j['delimiters'] as String? ?? '_',
      tag: j['tag'] as String? ?? 'gene',
      pattern: j['pattern'] as String? ?? '',
      regexOnFileName: j['regexOnFileName'] == true,
    );

Json organizerToJson(OrganizerConfig c) => {
      'sample': fieldRuleToJson(c.sample),
      'locus': fieldRuleToJson(c.locus),
      'duplicatePolicy': c.duplicatePolicy.name,
      'useSynonyms': c.useSynonyms,
      'sanitizeNames': c.sanitizeNames,
      'ignoreCase': c.ignoreCase,
      'ignoreLeadingZeros': c.ignoreLeadingZeros,
      'sampleAliases': c.sampleAliases,
      'formats': {
        for (final e in c.formats.entries)
          e.key: {'sample': fieldRuleToJson(e.value.sample), 'locus': fieldRuleToJson(e.value.locus), 'skip': e.value.skip},
      },
    };

OrganizerConfig organizerFromJson(Json j) => OrganizerConfig(
      sample: fieldRuleFromJson(asJson(j['sample'])),
      locus: fieldRuleFromJson(asJson(j['locus'])),
      duplicatePolicy: enumOf(DuplicatePolicy.values, j['duplicatePolicy'], DuplicatePolicy.keepLongest),
      useSynonyms: j['useSynonyms'] != false,
      sanitizeNames: j['sanitizeNames'] != false,
      ignoreCase: j['ignoreCase'] != false,
      ignoreLeadingZeros: j['ignoreLeadingZeros'] == true,
      sampleAliases: asJson(j['sampleAliases']).cast<String, String>(),
      formats: {
        for (final e in asJson(j['formats']).entries)
          e.key: FormatRule(
            sample: fieldRuleFromJson(asJson(asJson(e.value)['sample'])),
            locus: fieldRuleFromJson(asJson(asJson(e.value)['locus'])),
            skip: asJson(e.value)['skip'] == true,
          ),
      },
    );

Json fastaInputToJson(FastaInput f) => {'fileName': f.fileName, 'records': seqsToJson(f.records)};
FastaInput fastaInputFromJson(Json j) => FastaInput(j['fileName'] as String? ?? '', seqsFromJson(j['records']));

// --- step 4 ---------------------------------------------------------------

Json headerToJson(HeaderTemplate h) => {'fields': h.fields, 'separator': h.separator};
HeaderTemplate headerFromJson(Json j) => HeaderTemplate(
      fields: j['fields'] == null ? const ['accession', 'family', 'genus', 'species'] : strings(j['fields']),
      separator: j['separator'] as String? ?? '_',
    );

Json mafftToJson(MafftOptions m) => {
      'strategy': m.strategy,
      'threads': m.threads,
      'adjustDirection': m.adjustDirection,
      'maxIterate': m.maxIterate,
      'gapOpen': m.gapOpen,
      'gapExtend': m.gapExtend,
      'reorder': m.reorder,
      'preserveCase': m.preserveCase,
      'extra': m.extra,
      'amino': m.amino,
    };

MafftOptions mafftFromJson(Json j) => MafftOptions(
      strategy: j['strategy'] as String? ?? '--auto',
      threads: integer(j['threads']) ?? 1,
      adjustDirection: j['adjustDirection'] == true,
      maxIterate: integer(j['maxIterate']) ?? 0,
      gapOpen: dbl(j['gapOpen']) ?? 1.53,
      gapExtend: dbl(j['gapExtend']) ?? 0.0,
      reorder: j['reorder'] == true,
      preserveCase: j['preserveCase'] == true,
      extra: j['extra'] as String? ?? '',
      amino: j['amino'] == true,
    );

Json trimalToJson(TrimalOptions t) => {
      'method': t.method,
      'gapThreshold': t.gapThreshold,
      if (t.similarity != null) 'similarity': t.similarity,
      if (t.consistency != null) 'consistency': t.consistency,
      if (t.window != null) 'window': t.window,
      'extra': t.extra,
    };

TrimalOptions trimalFromJson(Json j) => TrimalOptions(
      method: j['method'] as String? ?? '-automated1',
      gapThreshold: dbl(j['gapThreshold']) ?? 0.9,
      similarity: dbl(j['similarity']),
      consistency: integer(j['consistency']),
      window: integer(j['window']),
      extra: j['extra'] as String? ?? '',
    );

// --- step 5 ---------------------------------------------------------------

Json thresholdsToJson(SaturationThresholds t) => {
      'slopeModerate': t.slopeModerate,
      'slopeSaturated': t.slopeSaturated,
      'undefinedModerate': t.undefinedModerate,
      'undefinedSaturated': t.undefinedSaturated,
      'issSaturated': t.issSaturated,
      'minPairs': t.minPairs,
      'minSites': t.minSites,
      'lowDivergence': t.lowDivergence,
    };

SaturationThresholds thresholdsFromJson(Json j, SaturationThresholds d) => SaturationThresholds(
      slopeModerate: dbl(j['slopeModerate']) ?? d.slopeModerate,
      slopeSaturated: dbl(j['slopeSaturated']) ?? d.slopeSaturated,
      undefinedModerate: dbl(j['undefinedModerate']) ?? d.undefinedModerate,
      undefinedSaturated: dbl(j['undefinedSaturated']) ?? d.undefinedSaturated,
      issSaturated: dbl(j['issSaturated']) ?? d.issSaturated,
      minPairs: integer(j['minPairs']) ?? d.minPairs,
      minSites: integer(j['minSites']) ?? d.minSites,
      lowDivergence: dbl(j['lowDivergence']) ?? d.lowDivergence,
    );

Json saturationToJson(SaturationReport r) => {
      'name': r.name,
      'taxa': r.taxa,
      'length': r.length,
      'frame': r.frame,
      if (r.frameStops != null) 'frameStops': r.frameStops,
      'pairs': [for (final x in r.positionPairs) [x.distance, x.p1, x.p2, x.p3]],
      'subsets': {
        for (final e in r.subsets.entries)
          e.key.name: {
            'sites': e.value.sites,
            'pairs': e.value.pairs,
            'undefinedPairs': e.value.undefinedPairs,
            'meanP': e.value.meanP,
            'meanDistance': e.value.meanDistance,
            'slope': e.value.slope,
            'r2': e.value.r2,
            'tsTv': e.value.tsTvRatio,
            'iss': e.value.iss,
            'base': e.value.baseFrequencies,
            'model': e.value.model,
            'verdict': e.value.verdict.name,
            'variable': e.value.variableSites,
            'informative': e.value.informativeSites,
            if (e.value.composition case final c?) 'composition': {'chi2': c.chi2, 'df': c.df, 'p': c.p, 'failing': c.failing, 'taxa': c.taxa},
            'points': [for (final x in e.value.points) [x.distance, x.p, x.ts, x.tv]],
          },
      },
    };

SaturationReport saturationFromJson(Json j) {
  List<double> nums(Object? o) => [for (final x in (o as List?) ?? const []) (x as num).toDouble()];
  final subsets = <CodonSubset, SubsetSaturation>{};
  for (final e in asJson(j['subsets']).entries) {
    final s = asJson(e.value);
    final c = s['composition'] == null ? null : asJson(s['composition']);
    final subset = enumOf(CodonSubset.values, e.key, CodonSubset.all);
    subsets[subset] = SubsetSaturation(
      subset: subset,
      sites: integer(s['sites']) ?? 0,
      pairs: integer(s['pairs']) ?? 0,
      undefinedPairs: integer(s['undefinedPairs']) ?? 0,
      meanP: dbl(s['meanP']) ?? 0,
      meanDistance: dbl(s['meanDistance']) ?? 0,
      slope: dbl(s['slope']),
      r2: dbl(s['r2']),
      tsTvRatio: dbl(s['tsTv']),
      iss: dbl(s['iss']),
      baseFrequencies: nums(s['base']),
      model: s['model'] as String? ?? '',
      verdict: enumOf(SaturationVerdict.values, s['verdict'], SaturationVerdict.insufficient),
      variableSites: integer(s['variable']) ?? 0,
      informativeSites: integer(s['informative']) ?? 0,
      composition: c == null
          ? null
          : CompositionTest(chi2: dbl(c['chi2']) ?? 0, df: integer(c['df']) ?? 0, p: dbl(c['p']) ?? 1, failing: integer(c['failing']) ?? 0, taxa: integer(c['taxa']) ?? 0),
      points: [
        for (final x in (s['points'] as List?) ?? const [])
          SaturationPoint(dbl((x as List)[0]) ?? 0, dbl(x[1]) ?? 0, dbl(x[2]) ?? 0, dbl(x[3]) ?? 0),
      ],
    );
  }
  return SaturationReport(
    j['name'] as String? ?? '',
    integer(j['taxa']) ?? 0,
    integer(j['length']) ?? 0,
    integer(j['frame']) ?? 0,
    subsets,
    j['frameStops'] == null ? null : [for (final x in j['frameStops'] as List) (x as num).toInt()],
    [
      for (final x in (j['pairs'] as List?) ?? const [])
        PositionPair(dbl((x as List)[0]) ?? 0, dbl(x[1]) ?? 0, dbl(x[2]) ?? 0, dbl(x[3]) ?? 0),
    ],
  );
}

// --- step 6 ---------------------------------------------------------------

Json concatToJson(ConcatOptions c) => {
      'missingPolicy': c.missingPolicy.name,
      'excludePositions': c.excludePositions.toList(),
      'perLocusExclude': {for (final e in c.perLocusExclude.entries) e.key: e.value.toList()},
      'scheme': c.scheme.name,
      'minLociPerTaxon': c.minLociPerTaxon,
      'replaceGapsWithMissing': c.replaceGapsWithMissing,
      'outgroups': c.outgroups,
      'dataType': c.dataType.name,
    };

ConcatOptions concatFromJson(Json j) {
  Set<int> ints(Object? o) => {for (final x in (o as List?) ?? const []) (x as num).toInt()};
  return ConcatOptions(
    missingPolicy: enumOf(MissingDataPolicy.values, j['missingPolicy'], MissingDataPolicy.fill),
    excludePositions: ints(j['excludePositions']),
    perLocusExclude: {for (final e in asJson(j['perLocusExclude']).entries) e.key: ints(e.value)},
    scheme: enumOf(PartitionScheme.values, j['scheme'], PartitionScheme.byLocusAndCodon),
    minLociPerTaxon: integer(j['minLociPerTaxon']) ?? 1,
    replaceGapsWithMissing: j['replaceGapsWithMissing'] == true,
    outgroups: strings(j['outgroups']),
    dataType: enumOf(SeqDataType.values, j['dataType'], SeqDataType.dna),
  );
}

// --- step 7 ---------------------------------------------------------------

Json iqtreeToJson(IqtreeOptions o) => {
      'model': o.model,
      'ufboot': o.ufboot,
      'alrt': o.alrt,
      'classicBootstrap': o.classicBootstrap,
      'threads': o.threads,
      if (o.maxThreads != null) 'maxThreads': o.maxThreads,
      if (o.seed != null) 'seed': o.seed,
      'outgroup': o.outgroup,
      'raxmlModelsOnly': o.raxmlModelsOnly,
      'bnni': o.bnni,
      if (o.nstop != null) 'nstop': o.nstop,
      if (o.perturbation != null) 'perturbation': o.perturbation,
      'runs': o.runs,
      'extra': o.extra,
    };

IqtreeOptions iqtreeFromJson(Json j) => IqtreeOptions(
      model: j['model'] as String? ?? 'MFP+MERGE',
      ufboot: integer(j['ufboot']) ?? 1000,
      alrt: integer(j['alrt']) ?? 1000,
      classicBootstrap: integer(j['classicBootstrap']) ?? 0,
      threads: j['threads'] as String? ?? 'AUTO',
      maxThreads: integer(j['maxThreads']),
      seed: integer(j['seed']),
      outgroup: strings(j['outgroup']),
      raxmlModelsOnly: j['raxmlModelsOnly'] == true,
      bnni: j['bnni'] == true,
      nstop: integer(j['nstop']),
      perturbation: dbl(j['perturbation']),
      runs: integer(j['runs']) ?? 1,
      extra: j['extra'] as String? ?? '',
    );

Json raxmlToJson(RaxmlOptions o) => {
      'bootstraps': o.bootstraps,
      'startTrees': o.startTrees,
      'maxThreads': o.maxThreads,
      'seed': o.seed,
      'outgroup': o.outgroup,
      'supportMetrics': o.supportMetrics,
      'fast': o.fast,
      'criterion': o.criterion,
      'extra': o.extra,
    };

RaxmlOptions raxmlFromJson(Json j) => RaxmlOptions(
      bootstraps: j['bootstraps'] as String? ?? 'autoMRE{1000}',
      startTrees: j['startTrees'] as String? ?? '',
      maxThreads: integer(j['maxThreads']) ?? 4,
      seed: integer(j['seed']) ?? 12345,
      outgroup: strings(j['outgroup']),
      supportMetrics: j['supportMetrics'] as String? ?? 'fbp,tbe',
      fast: j['fast'] == true,
      criterion: j['criterion'] as String? ?? 'BIC',
      extra: j['extra'] as String? ?? '',
    );

Json delimitationToJson(Delimitation d) => {
      if (d.model != null) 'model': d.model!.name,
      'species': d.species,
      if (d.nullLogL != null) 'nullLogL': d.nullLogL,
      if (d.logL != null) 'logL': d.logL,
      if (d.pvalue != null) 'pvalue': d.pvalue,
      if (d.lrtPassed != null) 'lrtPassed': d.lrtPassed,
      if (d.minbr != null) 'minbr': d.minbr,
      if (d.bestSpecies != null) 'bestSpecies': d.bestSpecies,
      if (d.alpha != null) 'alpha': d.alpha,
    };

List<List<String>> _groups(Object? o) => [for (final s in (o as List?) ?? const []) strings(s)];

Delimitation delimitationFromJson(Json j) => Delimitation(
      model: j['model'] == null ? null : enumOf(PtpModel.values, j['model'], PtpModel.multi),
      species: _groups(j['species']),
      nullLogL: dbl(j['nullLogL']),
      logL: dbl(j['logL']),
      pvalue: dbl(j['pvalue']),
      lrtPassed: j['lrtPassed'] as bool?,
      minbr: dbl(j['minbr']),
      bestSpecies: integer(j['bestSpecies']),
      alpha: dbl(j['alpha']),
    );

Json spartitionToJson(Spartition s) => {'name': s.name, 'species': s.species, if (s.score != null) 'score': s.score};
Spartition spartitionFromJson(Json j) => Spartition(j['name'] as String? ?? '', _groups(j['species']), j['score'] as String?);

/// The GenBank records, in JSON (with the user's edits).
List<Json> recordsToJson(List<GbRecord> records) => [for (final r in records) r.toJson()];
List<GbRecord> recordsFromJson(Object? o) => [for (final r in (o as List?) ?? const []) GbRecord.fromJson(asJson(r))];
