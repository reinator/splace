part of 'project_controller.dart';

/// Extension of SPLACE project files: a ZIP with
///
/// - `project.json`: every parameter of every step, readable by whoever
///   receives the file (format, version, what the data are);
/// - `state.json`: the records with their edits, the FASTA inputs, the
///   loci with their alignments, the identifications, saturation reports,
///   tree and delimitation runs and the log;
/// - `work/`: the project folder (loci, alignments, matrices, trees…).
///
/// Opening it unpacks `work/` into a new project folder and puts the
/// analysis back as it was saved, so no step has to be run again.
const projectExtension = 'splc';
const _projectFormat = 'splace-project';
const _projectVersion = 1;

/// What [_readProjectFile] gives back: the two JSON documents (decoded off
/// the UI thread) and the folder the work files went to.
class _ProjectData {
  _ProjectData(this.params, this.state, this.workDir);
  final Map<String, dynamic> params;
  final Map<String, dynamic> state;
  final String workDir;
}

/// Writes [params], [state] and the files of [workDir] to [target] (a
/// temporary file first, renamed at the end, so a failed save never
/// breaks the previous one).
Future<void> _writeProjectFile(String target, Map<String, dynamic> params, Map<String, dynamic> state, String workDir) => Isolate.run(() {
      final tmp = '$target.part';
      final encoder = ZipFileEncoder()..create(tmp);
      encoder.addArchiveFile(ArchiveFile.string('project.json', const JsonEncoder.withIndent('  ').convert(params)));
      encoder.addArchiveFile(ArchiveFile.string('state.json', jsonEncode(state)));
      final root = Directory(workDir);
      if (root.existsSync()) {
        for (final f in root.listSync(recursive: true).whereType<File>()) {
          encoder.addFileSync(f, 'work/${p.relative(f.path, from: workDir).replaceAll('\\', '/')}');
        }
      }
      encoder.closeSync();
      final out = File(target);
      if (out.existsSync()) out.deleteSync();
      File(tmp).renameSync(target);
    });

/// Reads a project file: the JSON documents, and the work files unpacked
/// into [workDir] (never outside it, whatever names the ZIP holds).
Future<_ProjectData> _readProjectFile(String file, String workDir) => Isolate.run(() {
      final input = InputFileStream(file);
      try {
        final archive = ZipDecoder().decodeStream(input);
        String? params, state;
        final root = p.normalize(p.absolute(workDir));
        for (final f in archive.files) {
          if (!f.isFile) continue;
          final name = f.name.replaceAll('\\', '/');
          if (name == 'project.json') {
            params = utf8.decode(f.content);
          } else if (name == 'state.json') {
            state = utf8.decode(f.content);
          } else if (name.startsWith('work/')) {
            final out = p.normalize(p.join(root, p.joinAll(name.substring(5).split('/'))));
            if (!p.isWithin(root, out)) continue;
            Directory(p.dirname(out)).createSync(recursive: true);
            final os = OutputFileStream(out);
            f.writeContent(os);
            os.closeSync();
          }
        }
        if (params == null || state == null) throw const FormatException('not a SPLACE project (project.json or state.json missing)');
        final j = jsonDecode(params) as Map<String, dynamic>;
        if (j['format'] != _projectFormat) throw const FormatException('not a SPLACE project');
        if ((j['version'] as num? ?? 0) > _projectVersion) {
          throw FormatException('made by a newer SPLACE (project version ${j['version']}); update SPLACE to open it');
        }
        return _ProjectData(j, jsonDecode(state) as Map<String, dynamic>, workDir);
      } finally {
        input.closeSync();
      }
    });

int _jsonHash(Object? o) => switch (o) {
      Map m => Object.hashAll([for (final e in m.entries) Object.hash(e.key, _jsonHash(e.value))]),
      Iterable l => Object.hashAll(l.map(_jsonHash)),
      _ => o.hashCode,
    };

Map<String, dynamic> _locusToJson(LocusState l, String Function(StepStatus) status) => {
      'name': l.name,
      'sources': l.sources.toList(),
      'raw': seqsToJson(l.raw),
      if (l.aligned != null) 'aligned': seqsToJson(l.aligned!),
      if (l.trimmed != null) 'trimmed': seqsToJson(l.trimmed!),
      if (l.proteinAligned != null) 'protein': seqsToJson(l.proteinAligned!),
      if (l.cds.isNotEmpty) 'cds': l.cds,
      if (l.codes.isNotEmpty) 'codes': l.codes,
      'coding': l.coding,
      'isProtein': l.isProtein,
      'partial': l.partial,
      'frame': l.frame,
      if (l.frameShifted.isNotEmpty) 'frameShifted': l.frameShifted,
      'trimmedFrameLost': l.trimmedFrameLost,
      if (l.trimmedColumns != null) 'trimmedColumns': l.trimmedColumns,
      'alignStatus': status(l.alignStatus),
      'trimStatus': status(l.trimStatus),
      if (l.message != null) 'message': l.message,
      if (l.trimMethod != null) 'trimMethod': l.trimMethod,
      'pipeline': l.pipeline,
      if (l.saturation != null) 'saturation': saturationToJson(l.saturation!),
      'included': l.included,
    };

LocusState _locusFromJson(Map<String, dynamic> j) {
  final l = LocusState(j['name'] as String, seqsFromJson(j['raw']), sources: strings(j['sources']).toSet())
    ..aligned = j['aligned'] == null ? null : seqsFromJson(j['aligned'])
    ..trimmed = j['trimmed'] == null ? null : seqsFromJson(j['trimmed'])
    ..proteinAligned = j['protein'] == null ? null : seqsFromJson(j['protein'])
    ..cds = asJson(j['cds']).cast<String, String>()
    ..codes = {for (final e in asJson(j['codes']).entries) e.key: (e.value as num).toInt()}
    ..coding = j['coding'] != false
    ..isProtein = j['isProtein'] == true
    ..partial = integer(j['partial']) ?? 0
    ..frame = integer(j['frame']) ?? 0
    ..frameShifted = strings(j['frameShifted'])
    ..trimmedFrameLost = j['trimmedFrameLost'] == true
    ..trimmedColumns = j['trimmedColumns'] == null ? null : [for (final x in j['trimmedColumns'] as List) (x as num).toInt()]
    ..alignStatus = enumOf(StepStatus.values, j['alignStatus'], StepStatus.idle)
    ..trimStatus = enumOf(StepStatus.values, j['trimStatus'], StepStatus.idle)
    ..message = j['message'] as String?
    ..trimMethod = j['trimMethod'] as String?
    ..saturation = j['saturation'] == null ? null : saturationFromJson(asJson(j['saturation']))
    ..included = j['included'] != false;
  l.pipeline.addAll(strings(j['pipeline']));
  final aln = l.alignment;
  if (aln != null && aln.isNotEmpty) l.stats = alignmentStats(aln);
  return l;
}

extension ProjectFileX on ProjectController {
  /// A path inside the project folder, relative to it (the folder changes
  /// when the project is opened elsewhere); other paths as they are.
  String? _rel(String? path) =>
      path == null ? null : (p.isWithin(workDir, path) || p.equals(workDir, path) ? p.relative(path, from: workDir).replaceAll('\\', '/') : path);

  String? _abs(String? path) => path == null ? null : (p.isAbsolute(path) ? path : p.normalize(p.join(workDir, p.joinAll(path.split('/')))));

  Map<String, String> _relMap(Map<String, String> m) => {for (final e in m.entries) e.key: _rel(e.value)!};
  Map<String, String> _absMap(Object? m) => {for (final e in asJson(m).entries) e.key: _abs(e.value as String)!};

  /// A step that was running when the project was saved did not finish.
  String _status(StepStatus s) => (s == StepStatus.running ? StepStatus.idle : s).name;

  // -------------------------------------------------------------------
  // Parameters (project.json)
  // -------------------------------------------------------------------

  /// Every parameter of the analysis, step by step (also what a reader of
  /// the file sees first).
  Map<String, dynamic> projectParameters(String appVersion) => {
        'format': _projectFormat,
        'version': _projectVersion,
        'app': 'SPLACE $appVersion',
        'saved': DateTime.now().toIso8601String(),
        'data': {
          'genbankRecords': records.length,
          'fastaFiles': inputs.length,
          'fastaSequences': inputs.fold<int>(0, (a, f) => a + f.records.length),
          'markers': loci.length,
          'markersIncluded': [for (final l in loci.values) if (l.included) l.name],
          'treeRuns': treeRuns.length,
          'delimitations': delimitations.length,
        },
        '1_load': {
          'sourcesUsed': loadOptionsUsed.toList()..sort(),
          'markerSearches': [for (final m in markerSearches) {...markerSearchToJson(m.search), 'records': m.records}],
          if (inputs.isNotEmpty) 'fastaOrganizer': organizerToJson(organizer),
          'organizerConfirmed': organizerConfirmed,
          'barcoding': barcodingMode,
          if (barcodeLocus != null) 'barcodeLocus': barcodeLocus,
        },
        '2_records': {
          'taxonomySource': taxonomySource.name,
          'taxonomyResolved': taxonomyResolved,
          'speciesIdentification': {'options': speciesIdOptions.toJson(), 'genes': speciesIdGenes},
        },
        '3_features_and_markers': {
          'featureTypes': featureTypes.toList(),
          'selectedMarkers': selectedMarkers.toList(),
          'duplicateChoices': duplicateChoices,
          'geneOverrides': geneOverrides,
          'featureOverrides': featureOverrides,
          'unknownGenesReviewed': unknownGenesReviewed,
        },
        '4_alignment': {
          'header': headerToJson(header),
          'headerConfirmed': headerConfirmed,
          'onePerTaxon': onePerTaxonMode,
          'sequenceType': seqMode.name,
          'geneticCode': geneticCode,
          'geneticCodeSetByUser': _geneticCodeTouched,
          'allowPseudo': allowPseudo,
          'allowInternalStops': allowInternalStops,
          'allowBacktranslationWarnings': allowBacktranslationWarnings,
          'mafft': mafftToJson(mafft),
          'trimTool': trimTool.name,
          'trimal': trimalToJson(trimal),
          'clipkitMode': clipkitMode,
          'trimSkipped': trimSkipped,
        },
        '5_saturation': {
          'thresholds': thresholdsToJson(saturationThresholds),
          'issReplicates': issReplicates,
          'plotPoints': satPlotPoints,
          'model': satModel,
          'excludeGaps': satExcludeGaps,
          if (_satParamsUsed != null) 'resultsMadeWith': _satParamsUsed,
        },
        '6_concatenation': {'options': concatToJson(concat), 'useTrimmed': useTrimmed},
        '7_inference': [
          for (final r in treeRuns)
            {
              'engine': r.engine,
              'started': r.started.toIso8601String(),
              'status': _status(r.status),
              if (r.command != null) 'command': r.command,
              if (r.iqtree != null) 'iqtree': iqtreeToJson(r.iqtree!),
              if (r.raxml != null) 'raxml': raxmlToJson(r.raxml!),
              if (r.bestModel != null) 'bestModel': r.bestModel,
              if (r.mooseModels.isNotEmpty) 'mooseModels': r.mooseModels,
            },
        ],
        '7_delimitation': [
          for (final d in delimitations)
            {
              'method': d.method.name,
              'source': d.source,
              'status': _status(d.status),
              if (d.command != null) 'command': d.command,
              if (d.distance != null) 'distance': d.distance!.name,
              if (d.result != null) 'species': d.result!.species.length,
            },
        ],
      };

  // -------------------------------------------------------------------
  // State (state.json)
  // -------------------------------------------------------------------

  Map<String, dynamic> _treeRunToJson(TreeRun r) => {
        'engine': r.engine,
        'directory': _rel(r.directory),
        'started': r.started.toIso8601String(),
        'status': _status(r.status),
        if (r.treeFile != null) 'treeFile': _rel(r.treeFile),
        if (r.command != null) 'command': r.command,
        if (r.commandDisplay != null) 'commandDisplay': r.commandDisplay,
        'outputs': _relMap(r.outputs),
        'notes': r.notes,
        'geneTrees': r.geneTrees,
        if (r.iqtree != null) 'iqtree': iqtreeToJson(r.iqtree!),
        if (r.raxml != null) 'raxml': raxmlToJson(r.raxml!),
        'models': r.models.name,
        if (r.fixedModel != null) 'fixedModel': r.fixedModel,
        'mooseModels': r.mooseModels,
        if (r.bestModel != null) 'bestModel': r.bestModel,
        if (r.error != null) 'error': r.error,
        'tail': r.tail,
      };

  TreeRun _treeRunFromJson(Map<String, dynamic> j) {
    final r = TreeRun(j['engine'] as String? ?? '', _abs(j['directory'] as String?) ?? workDir, started: DateTime.tryParse(j['started'] as String? ?? ''))
      ..status = enumOf(StepStatus.values, j['status'], StepStatus.idle)
      ..treeFile = _abs(j['treeFile'] as String?)
      ..command = j['command'] as String?
      ..commandDisplay = j['commandDisplay'] as String?
      ..iqtree = j['iqtree'] == null ? null : iqtreeFromJson(asJson(j['iqtree']))
      ..raxml = j['raxml'] == null ? null : raxmlFromJson(asJson(j['raxml']))
      ..models = enumOf(RaxmlModels.values, j['models'], RaxmlModels.moose)
      ..fixedModel = j['fixedModel'] as String?
      ..bestModel = j['bestModel'] as String?
      ..error = j['error'] as String?;
    r.outputs.addAll(_absMap(j['outputs']));
    r.notes.addAll(strings(j['notes']));
    r.geneTrees.addAll({for (final e in asJson(j['geneTrees']).entries) e.key: e.value == true});
    r.mooseModels.addAll(strings(j['mooseModels']));
    r.tail.addAll(strings(j['tail']));
    return r;
  }

  Map<String, dynamic> _delimitationToJson(DelimitationRun d) => {
        'method': d.method.name,
        'source': d.source,
        'directory': _rel(d.directory),
        if (d.run != null) 'run': treeRuns.indexOf(d.run!),
        if (d.tree != null) 'tree': _rel(d.tree),
        'started': d.started.toIso8601String(),
        'status': _status(d.status),
        if (d.result != null) 'result': delimitationToJson(d.result!),
        if (d.command != null) 'command': d.command,
        if (d.minbrCommand != null) 'minbrCommand': d.minbrCommand,
        if (d.error != null) 'error': d.error,
        if (d.minbrAuto != null) 'minbrAuto': d.minbrAuto,
        'outgroup': d.outgroup,
        'cropped': d.cropped,
        'partitions': [for (final s in d.partitions) spartitionToJson(s)],
        'chosen': d.chosen,
        if (d.distance != null) 'distance': d.distance!.name,
        'partitionNotes': d.partitionNotes,
        'outputs': _relMap(d.outputs),
      };

  DelimitationRun _delimitationFromJson(Map<String, dynamic> j) {
    final run = integer(j['run']);
    final d = DelimitationRun(
      method: enumOf(DelimMethod.values, j['method'], DelimMethod.mptp),
      source: j['source'] as String? ?? '',
      directory: _abs(j['directory'] as String?) ?? workDir,
      run: run != null && run >= 0 && run < treeRuns.length ? treeRuns[run] : null,
      tree: _abs(j['tree'] as String?),
      started: DateTime.tryParse(j['started'] as String? ?? ''),
    )
      ..status = enumOf(StepStatus.values, j['status'], StepStatus.idle)
      ..result = j['result'] == null ? null : delimitationFromJson(asJson(j['result']))
      ..command = j['command'] as String?
      ..minbrCommand = j['minbrCommand'] as String?
      ..error = j['error'] as String?
      ..minbrAuto = dbl(j['minbrAuto'])
      ..outgroup = strings(j['outgroup'])
      ..cropped = j['cropped'] == true
      ..chosen = integer(j['chosen']) ?? 0
      ..distance = j['distance'] == null ? null : enumOf(DelimDistance.values, j['distance'], DelimDistance.jc69);
    d.partitions.addAll([for (final s in (j['partitions'] as List?) ?? const []) spartitionFromJson(asJson(s))]);
    d.partitionNotes.addAll(strings(j['partitionNotes']));
    d.outputs.addAll(_absMap(j['outputs']));
    return d;
  }

  Map<String, dynamic> _projectState() => {
        'records': recordsToJson(records),
        'inputs': [for (final f in inputs) fastaInputToJson(f)],
        'seenFastaLoci': _seenFastaLoci.toList(),
        if (lociExportDir != null) 'lociExportDir': _rel(lociExportDir),
        if (ncbiStatus != null) 'ncbiStatus': ncbiStatus,
        'loci': [for (final l in loci.values) _locusToJson(l, _status)],
        'speciesIds': [for (final r in speciesIds.values) r.toJson()],
        'speciesIdFailures': {for (final e in speciesIdFailures.entries) e.key.name: e.value},
        'acceptedNames': _acceptedNames,
        'barcodeIds': {for (final e in barcodeIds.entries) e.key: e.value.toJson()},
        'barcodeNames': barcodeNames,
        'barcodeTaxa': {for (final e in barcodeTaxa.entries) e.key: e.value.toJson()},
        'trimPrompt': trimPrompt,
        'trimOpen': trimOpen,
        'matrix': matrix != null,
        'matrixFiles': _relMap(matrixFiles),
        'treeRuns': [for (final r in treeRuns) _treeRunToJson(r)],
        'delimitations': [for (final d in delimitations) _delimitationToJson(d)],
        if (treeToView != null) 'treeToView': _rel(treeToView),
        'treeFigures': {for (final e in treeFigures.entries) _rel(e.key)!: e.value},
        'zipSaved': zipSaved,
        'logs': [
          for (final e in logs) [e.time.toIso8601String(), e.source, e.level.name, e.message],
        ],
      };

  /// A hash of everything a save writes but the log and the date: equal
  /// when nothing changed. Long strings (sequences) cost nothing after the
  /// first time (the VM keeps a string's hash).
  int _projectFingerprint() {
    final params = projectParameters('')..remove('saved');
    final state = _projectState()
      ..remove('logs')
      ..remove('zipSaved');
    return _jsonHash([params, state]);
  }

  // -------------------------------------------------------------------
  // Save and open
  // -------------------------------------------------------------------

  /// Saves the whole project to [target] (`.splc`): parameters, results and
  /// the project folder. Later saves without a path go to the same file.
  Future<bool> saveProject(String target, {required String appVersion}) async {
    final path = p.extension(target).toLowerCase() == '.$projectExtension' ? target : '$target.$projectExtension';
    final ok = await _task('Saving project', (_) async {
      // Every file in the project folder goes in, so a step's files must
      // be written before it is saved.
      Directory(workDir).createSync(recursive: true);
      await _writeProjectFile(path, projectParameters(appVersion), _projectState(), workDir);
      projectFile = path;
      _markSaved();
      log('Project saved to $path', level: LogLevel.success);
      return true;
    });
    return ok ?? false;
  }

  /// Opens a project saved by [saveProject]: its files are unpacked into a
  /// new project folder ([into], or one named after the file in the
  /// projects folder) and every step is put back as it was.
  Future<bool> openProject(String file, {String? into}) async {
    final ok = await _task('Opening project', (_) async {
      final stamp = DateTime.now().toIso8601String().replaceAll(RegExp(r'[:.]'), '-').substring(0, 19);
      final stem = p.basenameWithoutExtension(file).replaceAll(RegExp(r'[^\w.-]+'), '_');
      final target = into ?? p.join(documentsDir, 'projects', '${stem}_$stamp');
      final data = await _readProjectFile(file, target);
      _workDir = data.workDir;
      await _applyProject(data.params, data.state);
      projectFile = file;
      _markSaved();
      log('Opened project $file (${data.params['app'] ?? 'SPLACE'}, saved ${data.params['saved'] ?? '?'}); files in $workDir',
          level: LogLevel.success);
      return true;
    });
    return ok ?? false;
  }

  Future<void> _applyProject(Map<String, dynamic> params, Map<String, dynamic> st) async {
    final load = asJson(params['1_load']);
    final rec = asJson(params['2_records']);
    final feat = asJson(params['3_features_and_markers']);
    final aln = asJson(params['4_alignment']);
    final sat = asJson(params['5_saturation']);
    final cat = asJson(params['6_concatenation']);

    // The log of the project first: what follows is added to it.
    logs
      ..clear()
      ..addAll([
        for (final e in (st['logs'] as List?) ?? const [])
          LogEntry('${(e as List)[1]}', '${e[3]}', level: enumOf(LogLevel.values, e[2], LogLevel.info), time: DateTime.tryParse('${e[0]}')),
      ]);
    logTick.value++;

    // 1. Data.
    records
      ..clear()
      ..addAll(recordsFromJson(st['records']));
    loadOptionsUsed
      ..clear()
      ..addAll([for (final x in (load['sourcesUsed'] as List?) ?? const []) (x as num).toInt()]);
    markerSearches
      ..clear()
      ..addAll([
        for (final m in (load['markerSearches'] as List?) ?? const [])
          (search: markerSearchFromJson(asJson(m)), records: integer(asJson(m)['records']) ?? 0),
      ]);
    ncbiStatus = st['ncbiStatus'] as String?;
    barcodingMode = load['barcoding'] == true;
    barcodeLocus = load['barcodeLocus'] as String?;

    // 3. Features and markers (the catalog reads the gene overrides).
    featureTypes
      ..clear()
      ..addAll(strings(feat['featureTypes']));
    selectedMarkers
      ..clear()
      ..addAll(strings(feat['selectedMarkers']));
    duplicateChoices
      ..clear()
      ..addAll({for (final e in asJson(feat['duplicateChoices']).entries) e.key: (e.value as num).toInt()});
    geneOverrides
      ..clear()
      ..addAll(asJson(feat['geneOverrides']).cast<String, String>());
    featureOverrides
      ..clear()
      ..addAll(asJson(feat['featureOverrides']).cast<String, String>());
    _rebuildCatalog();
    unknownGenesReviewed = feat['unknownGenesReviewed'] != false;

    // 2. Records: taxonomy and identification.
    taxonomySource = enumOf(TaxonomySource.values, rec['taxonomySource'], taxonomySource);
    taxonomyResolved = rec['taxonomyResolved'] == true;
    final id = asJson(rec['speciesIdentification']);
    if (id['options'] is Map) speciesIdOptions = SpeciesIdOptions.fromJson(asJson(id['options']));
    speciesIdGenes = strings(id['genes']);
    speciesIds = {
      for (final r in (st['speciesIds'] as List?) ?? const [])
        if (RecordIdentification.fromJson(asJson(r)) case final x) x.accession: x,
    };
    speciesIdFailures = {
      for (final e in asJson(st['speciesIdFailures']).entries) enumOf(IdSource.values, e.key, IdSource.genbank): '${e.value}',
    };
    _acceptedNames
      ..clear()
      ..addAll(asJson(st['acceptedNames']).cast<String, String?>());
    barcodeIds = {for (final e in asJson(st['barcodeIds']).entries) e.key: RecordIdentification.fromJson(asJson(e.value))};
    barcodeNames
      ..clear()
      ..addAll(asJson(st['barcodeNames']).cast<String, String>());
    barcodeTaxa
      ..clear()
      ..addAll({for (final e in asJson(st['barcodeTaxa']).entries) e.key: TaxonInfo.fromJson(asJson(e.value))});

    // 4. Alignment settings.
    header = headerFromJson(asJson(aln['header']));
    headerConfirmed = aln['headerConfirmed'] == true;
    onePerTaxonMode = aln['onePerTaxon'] == true;
    seqMode = enumOf(SeqMode.values, aln['sequenceType'], SeqMode.nucleotide);
    geneticCode = integer(aln['geneticCode']) ?? geneticCode;
    _geneticCodeTouched = aln['geneticCodeSetByUser'] == true;
    allowPseudo = aln['allowPseudo'] == true;
    allowInternalStops = aln['allowInternalStops'] == true;
    allowBacktranslationWarnings = aln['allowBacktranslationWarnings'] == true;
    if (aln['mafft'] is Map) {
      final m = mafftFromJson(asJson(aln['mafft']));
      mafft = m.copyWith(threads: m.threads.clamp(1, Platform.numberOfProcessors));
    }
    trimTool = enumOf(TrimTool.values, aln['trimTool'], TrimTool.trimalCodon);
    if (aln['trimal'] is Map) trimal = trimalFromJson(asJson(aln['trimal']));
    clipkitMode = aln['clipkitMode'] as String? ?? clipkitMode;
    trimSkipped = aln['trimSkipped'] == true;
    trimPrompt = st['trimPrompt'] == true;
    trimOpen = st['trimOpen'] == true;

    // 5-6.
    saturationThresholds = thresholdsFromJson(asJson(sat['thresholds']), phykitThresholds);
    issReplicates = integer(sat['issReplicates']) ?? issReplicates;
    satPlotPoints = integer(sat['plotPoints']) ?? satPlotPoints;
    satModel = sat['model'] as String? ?? satModel;
    satExcludeGaps = sat['excludeGaps'] == true;
    _satParamsUsed = sat['resultsMadeWith'] as String?;
    if (cat['options'] is Map) concat = concatFromJson(asJson(cat['options']));
    useTrimmed = cat['useTrimmed'] != false;

    // FASTA input: the split is computed again from the files and the
    // confirmed rule (markers already seen keep the user's selection).
    inputs
      ..clear()
      ..addAll([for (final f in (st['inputs'] as List?) ?? const []) fastaInputFromJson(asJson(f))]);
    if (load['fastaOrganizer'] is Map) organizer = organizerFromJson(asJson(load['fastaOrganizer']));
    organizerConfirmed = load['organizerConfirmed'] == true;
    organizerReview = false;
    lociExportDir = _abs(st['lociExportDir'] as String?);
    _seenFastaLoci
      ..clear()
      ..addAll(strings(st['seenFastaLoci']));
    if (inputs.isNotEmpty) await reorganize(immediate: true, confirmed: true);
    _markersCache = null;

    // The loci as saved, with their alignments: the caches make running a
    // step again reuse what did not change.
    loci
      ..clear()
      ..addAll({
        for (final l in (st['loci'] as List?) ?? const [])
          if (_locusFromJson(asJson(l)) case final x) x.name: x,
      });
    _lociBuiltFrom = loci.isEmpty ? null : _lociSignature;
    _alignCache.clear();
    _trimCache.clear();
    for (final l in loci.values) {
      if (l.aligned == null || l.alignStatus != StepStatus.done) continue;
      _alignCache[l.name] = _AlignSnapshot(l, _alignKey(l));
      final key = _trimKey(l);
      if (l.trimmed != null && key != null) {
        _trimCache[l.name] = (key: key, trimmed: l.trimmed!, method: l.trimMethod, frameLost: l.trimmedFrameLost, columns: l.trimmedColumns);
      }
    }

    // 6. The supermatrix is computed again (the same input gives the same
    // matrix); its files came with the project.
    matrixFiles
      ..clear()
      ..addAll(_absMap(st['matrixFiles']));
    matrix = null;
    overview = null;
    if (st['matrix'] == true) {
      try {
        final included = loci.values.where((l) => l.included && _aln(l) != null).toList();
        final options = concat.copyWith(dataType: proteinData ? SeqDataType.protein : SeqDataType.dna, scheme: effectiveScheme);
        matrix = await _matrixInIsolate([
          for (final l in included) LocusAlignment(l.name, _aln(l)!, coding: codonPositionsUsable(l), frame: _frame(l)),
        ], options);
        overview = await _overviewInIsolate(matrix!);
      } catch (e) {
        log('The supermatrix could not be rebuilt: $e', level: LogLevel.warn);
      }
    }

    // 7-8.
    treeRuns
      ..clear()
      ..addAll([for (final r in (st['treeRuns'] as List?) ?? const []) _treeRunFromJson(asJson(r))]);
    delimitations
      ..clear()
      ..addAll([for (final d in (st['delimitations'] as List?) ?? const []) _delimitationFromJson(asJson(d))]);
    treeToView = _abs(st['treeToView'] as String?);
    treeFigures
      ..clear()
      ..addAll({
        for (final e in ((st['treeFigures'] as Map?) ?? const {}).entries)
          if (e.value is Map) _abs(e.key as String)!: (e.value as Map).cast<String, dynamic>(),
      });
    zipSaved = st['zipSaved'] == true;
    _tipInfo = null;
    outgroupSuggestion = null;
    _outgroupKey = null;
    splashVisible = false;
    _changed();
  }
}
