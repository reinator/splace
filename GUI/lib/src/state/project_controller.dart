import 'dart:async';
import 'dart:convert';
import 'dart:io';
import 'dart:isolate';
import 'dart:math' as math;
import 'dart:ui' show Locale;

import 'package:archive/archive_io.dart';
import 'package:flutter/foundation.dart';
import 'package:flutter/material.dart' show ThemeMode;
import 'package:path/path.dart' as p;

import '../core/alignment_tools.dart';
import '../core/bold.dart';
import '../core/bold_submission.dart';
import '../core/delimitation.dart';
import '../core/fasta.dart';
import '../core/gene_suggestions.dart';
import '../core/genbank.dart';
import '../core/genbank_markers.dart';
import '../core/genbank_writer.dart';
import '../core/marker_search.dart';
import '../core/genetic_code.dart';
import '../core/locus_organizer.dart';
import '../core/marker_diagnostics.dart';
import '../core/mafft_progress.dart';
import '../core/tree_progress.dart';
import '../core/saturation.dart';
import '../core/supermatrix.dart';
import '../core/syngenes.dart';
import '../core/taxonomy_sources.dart';
import '../core/newick.dart';
import '../core/outgroup.dart';
import '../core/pairwise.dart' show alignPair;
import '../core/phykit.dart';
import '../core/tool_commands.dart';
import '../core/xlsx.dart';
import '../services/bold_portal.dart';
import '../services/ncbi.dart';
import '../services/process_runner.dart';
import '../services/species_id.dart';
import '../services/tool_locator.dart';
import '../core/project_json.dart';

part 'project_file.dart';

enum StepStatus { idle, running, done, error, skipped }

enum LogLevel { info, warn, error, success }

/// Step 5b "Sequence Type" of the HTML app.
enum SeqMode { nucleotide, aminoAcid, codon }

/// Trimming of DNA alignments: trimAl keeping whole codons (`-backtrans`),
/// ClipKIT keeping whole codons (`--codon`), trimAl column by column.
enum TrimTool { trimalCodon, clipkitCodon, trimalColumns }

class LogEntry {
  LogEntry(this.source, this.message, {this.level = LogLevel.info, DateTime? time}) : time = time ?? DateTime.now();
  final DateTime time;
  final String source;
  final String message;
  final LogLevel level;
  bool get error => level == LogLevel.error;
}

class LocusState {
  LocusState(this.name, this.raw, {this.sources = const {}});

  final String name;
  final List<SeqRecord> raw;

  /// CDS/rRNA/tRNA (GenBank feature types) and/or FASTA.
  final Set<String> sources;
  List<SeqRecord>? aligned;
  List<SeqRecord>? trimmed;

  /// Codon-aware mode: protein alignment and the validated CDS used for
  /// back-translation (header → CDS / genetic code).
  List<SeqRecord>? proteinAligned;
  Map<String, String> cds = const {};
  Map<String, int> codes = const {};
  bool coding = true;
  bool isProtein = false;

  /// GenBank records of this marker whose feature is partial (`<`/`>` in
  /// the location: fragments such as barcodes). Their saturation test is
  /// optional.
  int partial = 0;

  /// Most sequences are fragments: the saturation test is optional. A few
  /// partial genes in complete genomes do not count.
  bool get mostlyPartial => partial * 2 > raw.length;

  /// Frame of [aligned]; trimmed coding alignments are always frame 0.
  int frame = 0;

  /// Sequences of a coding DNA alignment with gaps that are not multiples
  /// of three (their codon positions are shifted from that gap on).
  List<String> frameShifted = const [];

  /// trimAl removed single columns of this coding DNA alignment, so the
  /// codon positions of [trimmed] are no longer known.
  bool trimmedFrameLost = false;

  /// Column of [aligned] that each column of [trimmed] comes from (trimAl;
  /// null for ClipKIT): the viewer and the saturation test refer to the
  /// MAFFT columns with it.
  List<int>? trimmedColumns;
  StepStatus alignStatus = StepStatus.idle;
  StepStatus trimStatus = StepStatus.idle;
  String? message;
  String? trimMethod;
  final List<String> pipeline = [];
  AlignmentStats? stats;
  SaturationReport? saturation;
  bool included = true;

  /// Best available alignment (trimmed > aligned) or null.
  List<SeqRecord>? get alignment => trimmed ?? aligned;

  /// Codon positions of [alignment] are known (coding DNA whose frame was
  /// not broken by trimming).
  bool get codonPositionsKnown => coding && !isProtein && !(trimmed != null && trimmedFrameLost);
  int get alignmentFrame => trimmed != null && coding ? 0 : frame;
}

/// A marker's alignment kept for re-runs, with the key of its input and
/// settings.
class _AlignSnapshot {
  _AlignSnapshot(LocusState l, this.key)
      : aligned = l.aligned!,
        proteinAligned = l.proteinAligned,
        cds = l.cds,
        codes = l.codes,
        isProtein = l.isProtein,
        coding = l.coding,
        frame = l.frame,
        frameShifted = l.frameShifted,
        pipeline = List.of(l.pipeline);
  final String key;
  final List<SeqRecord> aligned;
  final List<SeqRecord>? proteinAligned;
  final Map<String, String> cds;
  final Map<String, int> codes;
  final bool isProtein, coding;
  final int frame;
  final List<String> frameShifted;
  final List<String> pipeline;

  void restore(LocusState l) {
    l
      ..aligned = aligned
      ..proteinAligned = proteinAligned
      ..cds = cds
      ..codes = codes
      ..isProtein = isProtein
      ..coding = coding
      ..frame = frame
      ..frameShifted = frameShifted
      ..trimmedFrameLost = false
      ..trimmedColumns = null
      ..stats = alignmentStats(aligned)
      ..trimmed = null
      ..trimMethod = null
      ..trimStatus = StepStatus.idle
      ..saturation = null
      ..message = null
      ..alignStatus = StepStatus.done
      ..pipeline.clear()
      ..pipeline.addAll([...pipeline, 'cache']);
  }
}

class TreeRun {
  TreeRun(this.engine, this.directory, {DateTime? started}) : started = started ?? DateTime.now();
  final String engine;
  final String directory;
  final DateTime started;
  StepStatus status = StepStatus.running;
  String? treeFile;
  String? command;

  /// [command] with one option per line (for display).
  String? commandDisplay;
  final Map<String, String> outputs = {};
  final List<String> notes = [];

  /// Per-gene trees: gene → ok.
  final Map<String, bool> geneTrees = {};

  /// Settings of the run, for the methods text of "How to cite".
  IqtreeOptions? iqtree;
  RaxmlOptions? raxml;

  /// RAxML-NG: where the models came from ([fixedModel] when fixed).
  RaxmlModels models = RaxmlModels.moose;
  String? fixedModel;

  /// RAxML-NG: models chosen by MOOSE, one per partition.
  final List<String> mooseModels = [];

  /// Best model reported (IQ-TREE ModelFinder).
  String? bestModel;

  /// Why the run failed, and the last lines the program printed.
  String? error;
  final List<String> tail = [];

  void addTail(String line) {
    if (line.trim().isEmpty) return;
    tail.add(line.trimRight());
    if (tail.length > 40) tail.removeAt(0);
  }
}

/// One species delimitation: PTP or mPTP (mPTP program) on the ML tree of a
/// tree run, or ASAP or ABGD on an alignment.
class DelimitationRun {
  DelimitationRun({required this.method, required this.source, required this.directory, this.run, this.tree, DateTime? started})
      : started = started ?? DateTime.now();
  final DelimMethod method;

  /// Tree run whose ML tree was read (PTP, mPTP).
  final TreeRun? run;

  /// The program of the tree (IQ-TREE, RAxML-NG) or the alignment (a
  /// marker or the supermatrix).
  final String source;

  /// Tree file given to mPTP (the run's ML tree, as mPTP reads it).
  final String? tree;
  final String directory;
  final DateTime started;
  StepStatus status = StepStatus.running;
  Delimitation? result;
  String? command, minbrCommand, error;

  /// Minimum branch length detected from the p-distances (`--minbr_auto`),
  /// or null when mPTP's default (0.0001) was used.
  double? minbrAuto;
  List<String> outgroup = const [];
  bool cropped = false;

  /// ASAP and ABGD: every partition the program proposed (ASAP: best score
  /// first; ABGD: initial and recursive at each prior), the one shown,
  /// and the settings.
  final List<Spartition> partitions = [];
  int chosen = 0;
  DelimDistance? distance;

  /// ABGD: the prior of each partition of [partitions].
  final List<String> partitionNotes = [];

  /// Output files: delimitation text, SVG, SPART and tip → species table.
  final Map<String, String> outputs = {};

  /// Time of the tree run (tree methods) or of this run.
  String get time => (run?.started ?? started).toIso8601String().substring(11, 16);

  /// "mPTP · IQ-TREE", "ASAP · COI".
  String get label => '${method.label} · $source';

  /// Stable id (tree viewer buttons).
  String get id => '${run?.directory ?? directory}|${method.name}|$source';
}

enum ModalKind { ncbiSearch, ncbiFetch, taxonomy, speciesId, alignment, trimming, saturation, matrix, iqtree, raxml, bothTrees, delimitation }

/// State of the progress modal shown during NCBI/GBIF requests.
class TaskModal {
  TaskModal(this.kind, this.subject);
  final ModalKind kind;

  /// Taxon searched, number of accessions or taxonomy source.
  final String subject;
  int done = 0;
  int total = 0;
  String? current;

  /// Last GenBank record downloaded: shown as the species in italics with
  /// its voucher (or accession) in parentheses.
  String? species;
  String? voucher;

  /// One line per item (alignment/trimming: per marker), in order.
  final Map<String, ModalItem> items = {};

  /// Tree inference: stage and numbers read from IQ-TREE / RAxML-NG.
  TreeProgress? tree;

  /// Last lines printed by the running program.
  final List<String> tail = [];
  final DateTime started = DateTime.now();

  /// Several programs at once (IQ-TREE and RAxML-NG): the stages and the
  /// last lines of each, by name.
  final Map<String, TreeProgress> trees = {};
  final Map<String, List<String>> tails = {};

  void addTail(String line, [String? program]) {
    final l = line.trimRight();
    // Blank lines and the boxes IQ-TREE draws with - and |.
    if (l.trim().isEmpty || RegExp(r'^[\s|=*-]*$').hasMatch(l) || RegExp(r'^\s*\|.*\|\s*$').hasMatch(l)) return;
    final list = program != null && trees.containsKey(program) ? tails.putIfAbsent(program, () => []) : tail;
    list.add(l);
    if (list.length > 6) list.removeAt(0);
  }

  /// Final message; the modal shows it briefly before closing.
  String? result;
  bool failed = false;
  bool get finished => result != null;
}

enum ItemState { queued, running, done, cached, skipped, failed }

/// Progress of one item of a [TaskModal].
class ModalItem {
  double progress = 0;
  ItemState state = ItemState.queued;

  /// Stage or result, e.g. "Progressive alignment" or "1 552 bp".
  String? note;
}

/// One marker of step 4, merged across GenBank genes and FASTA loci.
class MarkerSummary {
  MarkerSummary(this.name);
  final String name;

  /// CDS/rRNA/tRNA and/or FASTA.
  final Set<String> sources = {};

  /// Sample key (`gb:<record>` or `fa:<sample>`) → sequence lengths.
  final Map<String, List<int>> lengths = {};
  int get samples => lengths.length;
  Iterable<int> get allLengths => lengths.values.expand((l) => l).where((l) => l > 0);
  int get duplicated => lengths.values.where((l) => l.length > 1).length;
  double get meanLength {
    final l = allLengths.toList();
    return l.isEmpty ? 0 : l.reduce((a, b) => a + b) / l.length;
  }
}

/// A sample row of the heatmaps: a GenBank record or a FASTA sample.
class SampleRef {
  const SampleRef(this.key, this.label, {this.detail, this.italic = false});
  final String key;
  final String label;
  final String? detail;
  final bool italic;
}

// Isolate entry points are top-level so closures capture only their
// arguments (never the controller, which is not sendable).

Future<({List<FastaInput> inputs, LayoutSuggestion suggestion})> _loadInIsolate(
        List<String> paths, List<FastaInput> existing) =>
    Isolate.run(() {
      final loaded = [
        for (final path in paths)
          FastaInput(p.basename(path), parseFasta(utf8.decode(File(path).readAsBytesSync(), allowMalformed: true))),
      ];
      final all = [...existing, ...loaded];
      return (inputs: loaded, suggestion: suggestLayout(all));
    });

const fastaExtensions = {'.fasta', '.fas', '.fa', '.fna', '.ffn', '.fsa', '.aln'};
const genbankExtensions = {'.gb', '.gbk', '.genbank', '.gbff'};

Future<List<String>> _listFiles(String dir, Set<String> exts) => Isolate.run(() => Directory(dir)
    .listSync(recursive: true)
    .whereType<File>()
    .map((f) => f.path)
    .where((path) => exts.contains(p.extension(path).toLowerCase()))
    .toList()
  ..sort());

Future<List<GbRecord>> _loadGenBankInIsolate(List<String> paths) => Isolate.run(() => [
      for (final path in paths)
        ...parseGenBankFile(
          utf8.decode(File(path).readAsBytesSync(), allowMalformed: true),
          fileName: p.basename(path),
          fileSize: File(path).lengthSync(),
          fileDate: File(path).lastModifiedSync(),
        ),
    ]);

Future<int> _writeLociInIsolate(OrganizerResult result, List<String> loci, String out) => Isolate.run(() {
      Directory(out).createSync(recursive: true);
      final fastas = result.locusFastas(loci: loci);
      for (final e in fastas.entries) {
        File(p.join(out, '${sanitizeTaxonName(e.key)}.fasta')).writeAsStringSync(e.value);
      }
      File(p.join(out, 'splace_locus_mapping.tsv')).writeAsStringSync(result.mappingTsv());
      return fastas.length;
    });

Future<OrganizerResult> _organizeInIsolate(List<FastaInput> inputs, OrganizerConfig config) =>
    Isolate.run(() => organizeLoci(inputs, config));

Future<({List<SeqRecord> records, List<String> reversed, AlignmentStats stats, bool coding, int frame, List<String> shifted})>
    _postAlignInIsolate(String fasta, String name, int code) => Isolate.run(() {
          final stripped = stripMafftReversedPrefix(parseFasta(fasta));
          final guess = detectCoding(name, stripped.records, code);
          return (
            records: stripped.records,
            reversed: stripped.reversed,
            stats: alignmentStats(stripped.records),
            coding: guess.coding,
            frame: guess.frame,
            shifted: guess.coding ? frameShiftedSequences(stripped.records) : const <String>[],
          );
        });

Future<MatrixOverview> _overviewInIsolate(Supermatrix m) => Isolate.run(() => matrixOverview(m));

Future<Supermatrix> _matrixInIsolate(List<LocusAlignment> loci, ConcatOptions options) =>
    Isolate.run(() => buildSupermatrix(loci, options));


Future<Map<String, List<SeqRecord>>> _geneFastasInIsolate(List<GbRecord> records, GeneCatalog cat, Set<String> types, Set<String> genes,
        HeaderTemplate header, Map<String, int> choices) =>
    Isolate.run(() => extractGeneFastas(records, cat, types, genes, header, choices));

Future<Map<String, List<CodonCandidate>>> _codonCandidatesInIsolate(List<GbRecord> records, GeneCatalog cat, Set<String> genes,
        HeaderTemplate header, Map<String, int> choices, int code, bool pseudo) =>
    Isolate.run(() => extractCodonCandidates(records, cat, genes, header, choices, fallbackCode: code, allowPseudo: pseudo));

/// Project subfolders → folders inside the ZIP.
Future<void> _zipInIsolate(String root, Map<String, String> folders, String target) => Isolate.run(() {
      final encoder = ZipFileEncoder()..create(target);
      for (final e in folders.entries) {
        final d = Directory(p.join(root, e.key));
        if (!d.existsSync()) continue;
        for (final f in d.listSync(recursive: true).whereType<File>()) {
          encoder.addFileSync(f, p.join(e.value, p.relative(f.path, from: d.path)).replaceAll('\\', '/'));
        }
      }
      encoder.closeSync();
    });

/// Holds the whole analysis state and runs every step without blocking the
/// UI: parsing and statistics run in isolates, external tools run as
/// asynchronous processes that stream their logs and can be cancelled.
class ProjectController extends ChangeNotifier {
  ProjectController({ToolLocator? locator, this.autoDetectTools = true, this.settingsFile, bool? watchUnsaved})
      : locator = locator ?? ToolLocator(),
        _watchUnsaved = watchUnsaved ?? settingsFile != null {
    _loadSettings();
    if (_watchUnsaved) addListener(_scheduleUnsavedCheck);
  }

  /// [unsaved] follows every change (the app: a controller with settings);
  /// tests call [checkUnsaved] themselves, so no timer is left behind.
  final bool _watchUnsaved;

  /// JSON file with preferences (theme, language, NCBI key, tool paths).
  /// Null in tests: nothing is persisted.
  final String? settingsFile;

  static String get defaultSettingsFile {
    final home = Platform.environment['USERPROFILE'] ?? Platform.environment['HOME'] ?? Directory.current.path;
    return p.join(home, 'Documents', 'SPLACE', 'settings.json');
  }

  // Display preferences (settings popover, same options as the HTML app).
  ThemeMode themeMode = ThemeMode.light;
  bool compact = false;
  double fontSize = 16;
  bool splashVisible = true;

  void setThemeMode(ThemeMode m) {
    themeMode = m;
    _changed(save: true);
  }

  void setCompact(bool v) {
    compact = v;
    _changed(save: true);
  }

  void setFontSize(double v) {
    fontSize = v.clamp(13, 20);
    _changed(save: true);
  }

  void setSplashVisible(bool v) {
    splashVisible = v;
    notifyListeners();
  }

  /// The tutorial (one page per step) opens after "Start App".
  bool tutorialAtStart = true;

  void setTutorialAtStart(bool v) {
    tutorialAtStart = v;
    _changed(save: true);
  }

  /// Set by the app: replaces this controller with a fresh one (a new
  /// project in a new folder); see [requestNewProject].
  void Function()? newProjectHandler;

  /// Set by the app: opens a project file in a fresh controller (see
  /// [requestOpenProject]).
  Future<void> Function(String file)? openProjectHandler;

  /// Project file (`.splc`) this project was last saved to or opened from;
  /// "Save project" writes to it again.
  String? projectFile;

  /// The project has changes that are not in [projectFile] (or, never
  /// saved, it has data): the window title says so. Checked a moment after
  /// each change, comparing a fingerprint of what would be saved.
  final ValueNotifier<bool> unsaved = ValueNotifier(false);

  /// Fingerprint of the project when it was last saved or opened.
  int? _savedPrint;
  Timer? _unsavedTimer;

  void _scheduleUnsavedCheck() {
    _unsavedTimer?.cancel();
    _unsavedTimer = Timer(const Duration(milliseconds: 700), checkUnsaved);
  }

  /// Updates [unsaved] now (a running step is checked when it ends).
  void checkUnsaved() {
    _unsavedTimer?.cancel();
    if (busy) return;
    final hasData = records.isNotEmpty || inputs.isNotEmpty;
    unsaved.value = projectFile == null ? hasData : _projectFingerprint() != _savedPrint;
  }

  /// What was just saved or opened is the reference.
  void _markSaved() {
    _savedPrint = _projectFingerprint();
    unsaved.value = false;
  }

  /// Opens [file] as a new project (this one stays in [workDir]).
  Future<void> requestOpenProject(String file) async {
    if (busy) return;
    await openProjectHandler?.call(file);
  }

  /// Marks this project busy while something outside it works (another
  /// project being opened); null clears it.
  void showBusy(String? label) {
    busyTask = label;
    notifyListeners();
  }

  /// Starts a new project: the app builds a new controller (same settings
  /// and tools, empty analysis, new project folder) and disposes this one.
  /// The files of this project stay in [workDir].
  void requestNewProject() {
    if (busy) return;
    log('New project; the previous one stays in $workDir');
    newProjectHandler?.call();
  }

  final ToolLocator locator;

  /// Probe external tools when the UI starts.
  final bool autoDetectTools;
  Map<ExternalTool, ResolvedTool?> tools = {};
  bool toolsChecked = false;
  PythonCommand? python;

  int threads = math.max(1, Platform.numberOfProcessors - 1);
  int geneticCode = 1;
  bool _geneticCodeTouched = false;
  String? _workDir;
  Locale? locale;

  // Logs are high-frequency (IQ-TREE prints thousands of lines), so they
  // notify through their own listenable instead of rebuilding every page.
  final List<LogEntry> logs = [];
  final ValueNotifier<int> logTick = ValueNotifier(0);
  static const _maxLogs = 20000;

  String? busyTask;
  double? progress;
  String? progressDetail;
  CancelToken? _cancel;
  bool get busy => busyTask != null;

  // ---------------------------------------------------------------------
  // Settings
  // ---------------------------------------------------------------------

  String? ncbiApiKey;
  bool apiHintDismissed = false;

  void _loadSettings() {
    final f = settingsFile;
    if (f == null || !File(f).existsSync()) return;
    try {
      final j = jsonDecode(File(f).readAsStringSync()) as Map<String, dynamic>;
      themeMode = ThemeMode.values.firstWhere((m) => m.name == j['theme'], orElse: () => ThemeMode.light);
      compact = j['compact'] == true;
      fontSize = ((j['fontSize'] as num?) ?? 16).toDouble().clamp(13, 20);
      final lang = j['locale'] as String?;
      locale = lang == null ? null : Locale(lang);
      ncbiApiKey = j['ncbiApiKey'] as String?;
      apiHintDismissed = j['apiHintDismissed'] == true;
      tutorialAtStart = j['tutorialAtStart'] != false;
      if (j['speciesId'] is Map) speciesIdOptions = SpeciesIdOptions.fromJson((j['speciesId'] as Map).cast<String, dynamic>());
      taxonomySource = TaxonomySource.values.firstWhere((v) => v.name == j['taxonomySource'], orElse: () => TaxonomySource.gbif);
      final hidden = j['hiddenColumns'];
      if (hidden is List) {
        hiddenColumns
          ..clear()
          ..addAll(hidden.cast<String>());
      }
      final overrides = j['toolOverrides'];
      if (overrides is Map) {
        for (final e in overrides.entries) {
          final tool = ExternalTool.values.where((t) => t.name == e.key).firstOrNull;
          if (tool != null) locator.overrides[tool] = '${e.value}';
        }
      }
    } catch (_) {
      // A corrupt settings file must never block the app.
    }
  }

  void _saveSettings() {
    final f = settingsFile;
    if (f == null) return;
    try {
      File(f)
        ..parent.createSync(recursive: true)
        ..writeAsStringSync(const JsonEncoder.withIndent('  ').convert({
          'theme': themeMode.name,
          'compact': compact,
          'fontSize': fontSize,
          'locale': locale?.languageCode,
          'ncbiApiKey': ncbiApiKey,
          'apiHintDismissed': apiHintDismissed,
          'tutorialAtStart': tutorialAtStart,
          'speciesId': speciesIdOptions.toJson(),
          'hiddenColumns': hiddenColumns.toList(),
          'taxonomySource': taxonomySource.name,
          'toolOverrides': {for (final e in locator.overrides.entries) e.key.name: e.value},
        }));
    } catch (e) {
      log('Could not save settings: $e', level: LogLevel.warn);
    }
  }

  void _changed({bool save = false}) {
    if (save) _saveSettings();
    notifyListeners();
  }

  void setNcbiApiKey(String? key) {
    ncbiApiKey = key?.trim().isEmpty ?? true ? null : key!.trim();
    _ncbi?.apiKey = ncbiApiKey;
    log(ncbiApiKey == null ? 'NCBI API key removed' : 'NCBI API key saved (10 requests/s)', level: LogLevel.success);
    _changed(save: true);
  }

  void dismissApiHint() {
    apiHintDismissed = true;
    _changed(save: true);
  }

  // ---------------------------------------------------------------------
  // Infrastructure
  // ---------------------------------------------------------------------

  void log(String message, {String source = 'SPLACE', bool error = false, LogLevel? level}) {
    logs.add(LogEntry(source, message, level: level ?? (error ? LogLevel.error : LogLevel.info)));
    if (logs.length > _maxLogs) logs.removeRange(0, logs.length - _maxLogs);
    logTick.value++;
  }

  void clearLogs() {
    logs.clear();
    logTick.value++;
  }

  String get documentsDir {
    final home = Platform.environment['USERPROFILE'] ?? Platform.environment['HOME'] ?? Directory.current.path;
    return p.join(home, 'Documents', 'SPLACE');
  }

  String get workDir {
    if (_workDir != null) return _workDir!;
    final stamp = DateTime.now().toIso8601String().replaceAll(RegExp(r'[:.]'), '-').substring(0, 19);
    _workDir = p.join(documentsDir, 'projects', stamp);
    return _workDir!;
  }

  set workDir(String value) {
    _workDir = value;
    notifyListeners();
  }

  String dir(String sub) {
    final d = p.join(workDir, sub);
    Directory(d).createSync(recursive: true);
    return d;
  }

  /// BOLD Systems public records (individual markers).
  BoldPortal boldPortal = BoldPortal();

  NcbiClient? _ncbi;
  NcbiClient get ncbi =>
      _ncbi ??= NcbiClient(cacheDir: p.join(documentsDir, 'cache', 'genbank'), apiKey: ncbiApiKey);

  void setThreads(int value) {
    threads = value.clamp(1, Platform.numberOfProcessors);
    notifyListeners();
  }

  void setGeneticCode(int code) {
    geneticCode = code;
    _geneticCodeTouched = true;
    for (final l in loci.values) {
      l.saturation = null;
    }
    notifyListeners();
  }

  void setLocale(Locale? value) {
    locale = value;
    _changed(save: true);
  }

  Future<void> detectTools() async {
    locator.clearCache();
    tools = await locator.resolveAll();
    // ClipKIT and PhyKIT ship inside the app (standalone, with their own
    // Python); a Python of the system is only looked for when one is
    // missing (development runs), to install it with pip.
    final needsPython = ExternalTool.values.any((t) => t.pythonModule != null && tools[t] == null);
    python = needsPython ? await findPython() : null;
    toolsChecked = true;
    for (final e in tools.entries) {
      final t = e.value;
      if (t == null && e.key.optional) continue;
      log(t == null
          ? '${e.key.label}: not found. ${e.key.installHint}'
          : '${e.key.label} ${t.version ?? ''} → ${t.description}', level: t == null ? LogLevel.warn : LogLevel.info);
    }
    notifyListeners();
  }

  void setToolOverride(ExternalTool tool, String? path) {
    if (path == null || path.isEmpty) {
      locator.overrides.remove(tool);
    } else {
      locator.overrides[tool] = path;
    }
    _saveSettings();
    detectTools();
  }

  /// Installs ClipKIT with pip (`--user --upgrade`), as the Electron app.
  Future<void> installClipkit() => _installPython(ExternalTool.clipkit);

  /// Installs PhyKIT (saturation test, step 6) with pip.
  Future<void> installPhykit() => _installPython(ExternalTool.phykit);

  Future<void> _installPython(ExternalTool target) => _task('Installing ${target.label}', (cancel) async {
        final py = python ?? await findPython();
        if (py == null) throw StateError('Python 3 was not found in PATH.');
        final tool = ResolvedTool(target, py.executable, ToolRuntime.native, prefixArgs: py.args);
        final out = await runTool(tool, ['-m', 'pip', 'install', '--user', '--upgrade', target.pythonModule!],
            workingDirectory: Directory.systemTemp.path, cancel: cancel, onLog: (m) => log(m, source: 'pip'));
        if (!out.ok) throw StateError('pip exited with code ${out.exitCode}');
        locator.clearCache();
        tools[target] = await locator.resolve(target);
      });

  /// Progress modal of NCBI/GBIF requests (`#progressModal` of the HTML app).
  TaskModal? modal;
  Timer? _modalTimer;

  void closeModal() {
    _modalTimer?.cancel();
    modal = null;
    notifyListeners();
  }

  Future<T?> _task<T>(String label, Future<T> Function(CancelToken cancel) body, {TaskModal? modal}) async {
    if (busy) {
      log('Another task is running: $busyTask', level: LogLevel.warn);
      return null;
    }
    busyTask = label;
    progress = null;
    progressDetail = null;
    _cancel = CancelToken();
    _modalTimer?.cancel();
    this.modal = modal;
    notifyListeners();
    final watch = Stopwatch()..start();
    try {
      final result = await body(_cancel!);
      if (!_cancel!.isCancelled) {
        log('$label finished in ${(watch.elapsedMilliseconds / 1000).toStringAsFixed(1)} s', level: LogLevel.success);
      }
      modal?.result ??= _cancel!.isCancelled ? 'cancelled' : 'done';
      return result;
    } catch (e, st) {
      log('$label failed: $e', error: true);
      debugPrintStack(stackTrace: st);
      modal
        ?..result = '$e'
        ..failed = true;
      return null;
    } finally {
      // Like the HTML modal: show the outcome briefly, then close.
      // Per-marker lists stay open when a marker failed or was skipped, so
      // the user can read why.
      final stay = modal != null && modal.items.values.any((i) => i.state == ItemState.failed || i.state == ItemState.skipped);
      if (modal != null && identical(this.modal, modal) && !stay) {
        _modalTimer = Timer(Duration(milliseconds: modal.failed ? 4000 : modal.items.isEmpty ? 1400 : 2500), () {
          if (identical(this.modal, modal)) closeModal();
        });
      }
      busyTask = null;
      progress = null;
      progressDetail = null;
      _cancel = null;
      notifyListeners();
    }
  }

  void cancelTask() {
    if (_cancel == null) return;
    log('Cancelling: $busyTask', level: LogLevel.warn);
    _cancel!.cancel();
  }

  void _progress(int done, int total, [String? detail]) {
    progress = total == 0 ? null : done / total;
    if (detail != null) progressDetail = detail;
    modal
      ?..done = done
      ..total = total
      ..current = detail ?? modal!.current;
    // The last step is shown at once; the others are gathered.
    if (done >= total) {
      notifyListeners();
    } else {
      _notifySoon();
    }
  }

  Future<ResolvedTool> _require(ExternalTool tool) async {
    if (!toolsChecked) await detectTools();
    final t = tools[tool];
    if (t == null) throw StateError('${tool.label} not found. ${tool.installHint}');
    return t;
  }

  /// Runs [task] over [items] with at most [concurrency] in flight. With
  /// [label], the progress modal lists the items running.
  /// With [reportProgress] false the task reports the overall progress.
  Future<void> _pool<T>(List<T> items, int concurrency, CancelToken cancel, Future<void> Function(T) task,
      {String Function(T)? label, bool reportProgress = true}) async {
    var next = 0, done = 0;
    final running = <String>[];
    if (reportProgress) _progress(0, items.length);
    Future<void> worker() async {
      while (next < items.length && !cancel.isCancelled) {
        final item = items[next++];
        final name = label?.call(item);
        if (name != null) {
          running.add(name);
          modal?.current = running.join(', ');
          notifyListeners();
        }
        await task(item);
        if (name != null) {
          running.remove(name);
          modal?.current = running.join(', ');
        }
        done++;
        if (reportProgress) _progress(done, items.length);
      }
    }

    await Future.wait(List.generate(math.max(1, math.min(concurrency, items.length)), (_) => worker()));
  }

  // ---------------------------------------------------------------------
  // 1. Load data: NCBI, GenBank files, FASTA files
  // ---------------------------------------------------------------------

  final List<GbRecord> records = [];
  GeneCatalog? catalog;

  GeneReview? _geneReview;
  GeneCatalog? _geneReviewOf;

  /// Name suggestions for unknown genes and /gene vs /product conflicts,
  /// computed once per catalog.
  GeneReview? get geneReview {
    final cat = catalog;
    if (cat == null) return null;
    if (!identical(_geneReviewOf, cat)) {
      _geneReview = reviewGeneNames(records, cat, geneticCode: geneticCode, overrides: geneOverrides, featureOverrides: featureOverrides);
      _geneReviewOf = cat;
    }
    return _geneReview;
  }

  /// Step 2 is resolved once taxonomy was fetched (GBIF) or skipped.
  bool taxonomyResolved = false;
  final Set<String> hiddenColumns = {'source', 'accession', 'kingdom', 'phylum'};
  String? ncbiStatus;

  OrganelleType get organelle => catalog?.organelle ?? OrganelleType.mt;

  int _mergeRecords(List<GbRecord> incoming) {
    final seen = {for (final r in records) r.accession.trim()};
    var added = 0;
    for (final r in incoming) {
      if (r.accession.trim().isNotEmpty && seen.contains(r.accession.trim())) continue;
      records.add(r);
      seen.add(r.accession.trim());
      added++;
    }
    if (added > 0) _rebuildCatalog(newRecords: true);
    return added;
  }

  void _rebuildCatalog({bool newRecords = false}) {
    if (records.isEmpty) {
      catalog = null;
      _markersCache = null;
      return;
    }
    catalog = buildGeneCatalog(records, overrides: geneOverrides, featureOverrides: featureOverrides);
    if (!_geneticCodeTouched) geneticCode = catalog!.organelle == OrganelleType.cp ? 11 : 2;
    // New records need the gene review (unknown names and /gene vs
    // /product conflicts) before the taxonomy.
    if (newRecords) unknownGenesReviewed = catalog!.unknown.isEmpty && (geneReview?.conflicts.isEmpty ?? true);
    _markersCache = null;
  }

  Future<void> addGenBankFiles(List<String> paths) => _task('Reading ${paths.length} GenBank file(s)', (_) async {
        if (paths.isEmpty) return;
        final loaded = await _loadGenBankInIsolate(paths);
        final added = _mergeRecords(loaded);
        if (added > 0) loadOptionsUsed.add(3);
        log('Loaded $added GenBank record(s) from ${paths.length} file(s)'
            '${loaded.length > added ? ' (${loaded.length - added} duplicate accession(s) skipped)' : ''}');
        if (loaded.isEmpty) log('No GenBank records found in the selected files', level: LogLevel.warn);
      });

  /// Files and/or folders (drag and drop): GenBank and FASTA are told apart
  /// by extension.
  Future<void> addPaths(List<String> paths) async {
    final gb = <String>[], fa = <String>[];
    for (final path in paths) {
      if (Directory(path).existsSync()) {
        gb.addAll(await _listFiles(path, genbankExtensions));
        fa.addAll(await _listFiles(path, fastaExtensions));
      } else if (genbankExtensions.contains(p.extension(path).toLowerCase())) {
        gb.add(path);
      } else {
        fa.add(path);
      }
    }
    if (gb.isNotEmpty) await addGenBankFiles(gb);
    if (fa.isNotEmpty) await addFiles(fa);
  }

  Future<void> addGenBankDirectory(String path) async {
    final files = await _listFiles(path, genbankExtensions);
    if (files.isEmpty) {
      log('No GenBank files (.gb, .gbk, .genbank) found in $path', error: true);
      return;
    }
    await addGenBankFiles(files);
  }

  /// Option 2: fetch accession numbers from NCBI.
  Future<void> fetchAccessions(String text) {
    final ids = text.split(RegExp(r'[\s,;]+')).map((s) => s.trim()).where((s) => s.isNotEmpty).toSet().toList();
    if (ids.isEmpty) return Future.value();
    return _task('Fetching GenBank records from NCBI', (cancel) async {
      await _download(ids, cancel, label: ids.join(', '));
      if (records.isNotEmpty) loadOptionsUsed.add(2);
    },
        modal: TaskModal(ModalKind.ncbiFetch, '${ids.length}')..total = ids.length);
  }

  Future<void> _download(List<String> ids, CancelToken cancel, {required String label}) async {
    log('Fetching ${ids.length} record(s) from NCBI${ncbi.hasKey ? ' with API key' : ''}…', source: 'NCBI');
    final r = await ncbi.fetchGenBank(
      ids,
      cancel: cancel,
      onProgress: (d, t, current) => _progress(d, t, current),
      onRecord: (r) => modal
        ?..species = r.organism
        ..voucher = [?r.voucher, r.accession].join(' · '),
      log: (m) => log(m, source: 'NCBI'),
    );
    final added = _mergeRecords(r.records);
    ncbiStatus = 'Fetched ${r.records.length}/${ids.length} records';
    modal?.result = ncbiStatus;
    modal?.failed = r.records.isEmpty;
    log('Loaded $added new record(s); ${r.failed.isEmpty ? 'no failures' : 'failed: ${r.failed.join(', ')}'}',
        source: 'NCBI', level: r.failed.isEmpty ? LogLevel.success : LogLevel.warn);
  }

  /// Option 1: search NCBI by taxon; returns the genomes to review before
  /// importing (null when the search failed).
  Future<List<GenomeCandidate>?> searchGenomes(GenomeSearch search) =>
      _task(modal: TaskModal(ModalKind.ncbiSearch, search.taxon)..total = 2, 'Searching NCBI for ${search.taxon}', (cancel) async {
        log('Query: ${search.term}', source: 'NCBI');
        // Steps: the search, then one summary request per 200 ids.
        final ids = await ncbi.searchGenomes(search);
        _progress(1, 1 + math.max(1, (ids.length / 200).ceil()).toInt());
        if (ids.isEmpty) {
          ncbiStatus = 'No records found for ${search.taxon}';
          log(ncbiStatus!, source: 'NCBI', level: LogLevel.warn);
          modal
            ?..result = ncbiStatus
            ..failed = true;
          return <GenomeCandidate>[];
        }
        modal?.current = '${ids.length} ids';
        final steps = 1 + math.max(1, (ids.length / 200).ceil()).toInt();
        final candidates = await ncbi.summaries(ids, onProgress: (d, _) => _progress(1 + d, steps));
        _progress(steps, steps);
        final complete = candidates.where((c) => c.classification == GenomeClass.complete).length;
        ncbiStatus = 'Found ${candidates.length} genome(s) ($complete complete); review before importing';
        log(ncbiStatus!, source: 'NCBI');
        // The review dialog replaces the modal right away.
        _modalTimer?.cancel();
        modal = null;
        return candidates;
      });

  /// Option 1 for individual markers: one NCBI query per marker (SynGenes
  /// names) and, for COI and the chloroplast barcodes, the public records
  /// of BOLD Systems; the records found, each tagged with its markers and
  /// its source, to review.
  Future<List<GenomeCandidate>?> searchMarkers(MarkerSearch search) => _task(
          modal: TaskModal(ModalKind.ncbiSearch, '${search.taxon} · ${search.genes.join(', ')}')
            ..total = search.genes.length + 1 + (search.boldGenes.isEmpty ? 0 : 1),
          'Searching NCBI${search.boldGenes.isEmpty ? '' : ' and BOLD'} for ${search.genes.join(', ')} of ${search.taxon}', (cancel) async {
        final steps = search.genes.length + 1 + (search.boldGenes.isEmpty ? 0 : 1);
        final byId = <String, Set<String>>{};
        final counts = <String, int>{};
        for (final (i, gene) in search.genes.indexed) {
          if (cancel.isCancelled) break;
          modal?.current = gene;
          final ids = await ncbi.searchTerm(search.termFor(gene), retmax: 20000, cancel: cancel);
          counts[gene] = ids.length;
          for (final id in ids) {
            (byId[id] ??= {}).add(gene);
          }
          log('$gene: ${ids.length} record(s)', source: 'NCBI');
          _progress(i + 1, steps, gene);
        }
        final candidates = <GenomeCandidate>[];
        if (byId.isNotEmpty) {
          final ids = byId.keys.toList();
          candidates.addAll(await ncbi.summaries(ids, onProgress: (d, n) => modal?.current = 'summaries ${d * 200}/${ids.length}'));
          for (final c in candidates) {
            c.markers.addAll(byId[c.id] ?? const {});
          }
        }
        _progress(search.genes.length + 1, steps);
        var boldNote = '';
        if (search.boldGenes.isNotEmpty && !cancel.isCancelled) {
          modal?.current = 'BOLD Systems';
          try {
            final bold = await _searchBold(search, candidates, cancel);
            boldNote = '; BOLD: ${bold.added} only in BOLD, ${bold.shared} also in GenBank';
          } catch (e) {
            if (e is BoldCancelled) rethrow;
            boldNote = '; BOLD could not be searched';
            log('BOLD Systems could not be searched: $e', source: 'BOLD', level: LogLevel.warn);
          }
          _progress(steps, steps);
        }
        if (candidates.isEmpty) {
          ncbiStatus = 'No records of ${search.genes.join(', ')} found for ${search.taxon}';
          log(ncbiStatus!, source: 'NCBI', level: LogLevel.warn);
          modal
            ?..result = ncbiStatus
            ..failed = true;
          return <GenomeCandidate>[];
        }
        ncbiStatus = 'Found ${candidates.length} record(s): ${counts.entries.map((e) => '${e.key} ${e.value}').join(', ')}$boldNote; review before importing';
        log(ncbiStatus!, source: 'NCBI');
        _modalTimer?.cancel();
        modal = null;
        return candidates;
      });

  /// Adds the BOLD records of [search] to [candidates]: a record BOLD
  /// mirrors from GenBank marks the NCBI one (NCBI + BOLD); the others
  /// become candidates of their own, imported from BOLD.
  Future<({int added, int shared})> _searchBold(MarkerSearch search, List<GenomeCandidate> candidates, CancelToken cancel) async {
    final terms = await boldPortal.taxonTerms(search.boldNames, cancel: cancel);
    final skipped = search.boldNames.where((n) => !terms.any((t) => t.toLowerCase().endsWith(':${n.toLowerCase()}'))).toList();
    if (skipped.isNotEmpty) log('BOLD does not hold as a taxon: ${skipped.join(', ')}', source: 'BOLD');
    if (terms.isEmpty) {
      log('BOLD Systems: no taxon matches ${search.taxon}', source: 'BOLD', level: LogLevel.warn);
      return (added: 0, shared: 0);
    }
    log('BOLD query: ${terms.join(';')} · markers ${search.boldCodes.join(', ')}', source: 'BOLD');
    final counts = await boldPortal.markerCounts(terms, cancel: cancel);
    final total = counts.values.fold(0, (a, b) => a + b);
    final wanted = search.boldCodes.fold(0, (a, c) => a + (counts[c] ?? 0));
    if (wanted == 0) {
      log('BOLD Systems: no public ${search.boldCodes.join('/')} records for ${search.taxon}', source: 'BOLD');
      return (added: 0, shared: 0);
    }
    if (total > BoldPortal.maxRecords) {
      log('BOLD Systems: $total records for ${search.taxon} is too many to download; search a lower rank', source: 'BOLD', level: LogLevel.warn);
      return (added: 0, shared: 0);
    }
    final rows = await boldPortal.records(terms, search.boldCodes, cancel: cancel);
    final byAccession = {for (final c in candidates) c.accession.replaceFirst(RegExp(r'\.\d+$'), ''): c};
    // A specimen with several markers (rbcL and matK): one record each.
    final perSpecimen = <String, int>{};
    for (final r in rows) {
      perSpecimen.update(r.processId, (v) => v + 1, ifAbsent: () => 1);
    }
    var added = 0, shared = 0;
    // A species searched with its synonyms: the records BOLD files under a
    // synonym take the accepted name (as NCBI does), or each name would
    // become a species of its own.
    final accepted = isSpeciesName(search.taxon) ? search.taxon.trim() : null;
    final synonyms = {for (final n in search.synonyms) n.toLowerCase()};
    for (final found in rows) {
      final r = accepted != null && found.name.toLowerCase() != accepted.toLowerCase() && synonyms.contains(found.name.toLowerCase())
          ? found.renamed(accepted)
          : found;
      final gene = r.gene;
      if (gene == null || !search.boldGenes.contains(gene) || r.sequence.length < 100) continue;
      final mirror = byAccession[r.insdc.replaceFirst(RegExp(r'\.\d+$'), '')];
      if (r.insdc.isNotEmpty && mirror != null) {
        mirror
          ..sources.add('BOLD')
          ..boldId ??= r.processId
          ..bin ??= r.bin.isEmpty ? null : r.bin
          ..country ??= r.country.isEmpty ? null : r.country;
        shared++;
        continue;
      }
      final acc = (perSpecimen[r.processId] ?? 1) > 1 ? '${r.processId}_${r.markerCode}' : r.processId;
      final rec = boldToGbRecord(r, accession: acc);
      candidates.add(GenomeCandidate(
        id: 'bold:$acc',
        accession: acc,
        title: rec.definition,
        species: r.name,
        description: rec.definition,
        filedAs: r.boldName.isEmpty ? null : r.boldName,
      )
        ..length = r.sequence.length
        ..markers.add(gene)
        ..sources.remove('NCBI')
        ..sources.add('BOLD')
        ..genbank = r.insdc.isEmpty ? null : r.insdc
        ..boldRecord = rec
        ..boldId = r.processId
        ..bin = r.bin.isEmpty ? null : r.bin
        ..country = r.country.isEmpty ? null : r.country);
      added++;
    }
    log('BOLD Systems: ${rows.length} record(s) of ${search.boldCodes.join('/')}; $added only in BOLD, $shared already in GenBank', source: 'BOLD');
    return (added: added, shared: shared);
  }

  /// Options of step 1 that loaded data (1 taxon search, 2 accessions,
  /// 3 GenBank files, 4 FASTA files), for the step navigator.
  final Set<int> loadOptionsUsed = {};

  /// Marker searches whose records were imported (methods text).
  final List<({MarkerSearch search, int records})> markerSearches = [];

  /// Imports the records of a marker search and joins them by taxon.
  Future<void> importMarkers(List<GenomeCandidate> selected, MarkerSearch search) async {
    final fromBold = [for (final c in selected) ?c.boldRecord];
    final fromNcbi = [for (final c in selected) if (c.boldRecord == null) c.id];
    await _task('Importing ${selected.length} record(s) of ${search.taxon}', (cancel) async {
      if (fromNcbi.isNotEmpty) await _download(fromNcbi, cancel, label: search.taxon);
      if (fromBold.isNotEmpty) {
        // Saved as GenBank files, so the project keeps what was imported.
        final dir = Directory(p.join(this.dir('00_ncbi'), 'bold'))..createSync(recursive: true);
        for (final r in fromBold) {
          File(p.join(dir.path, '${r.accession.replaceAll(RegExp(r'[^A-Za-z0-9_.-]'), '_')}.gb')).writeAsStringSync(writeGenBank(r));
        }
        final added = _mergeRecords(fromBold);
        log('Loaded $added record(s) from BOLD Systems → ${dir.path}', source: 'BOLD', level: LogLevel.success);
      }
    }, modal: TaskModal(ModalKind.ncbiFetch, search.taxon)..total = selected.length);
    if (records.isEmpty) return;
    loadOptionsUsed.add(1);
    // The queries come from SynGenes: kept for the methods text and saved
    // with the project so they can be reported or repeated.
    markerSearches.add((search: search, records: selected.length));
    try {
      final file = File(p.join(dir('00_ncbi'), 'marker_queries.txt'));
      file.writeAsStringSync(
        '# ${DateTime.now().toIso8601String()} · ${search.taxon} · ${selected.length} record(s) imported\n'
        '# Entrez (nuccore) queries built from the SynGenes names of each marker\n${search.preview}\n\n',
        mode: FileMode.append,
      );
      log('Marker queries saved → ${file.path}');
    } catch (e) {
      log('Could not save the marker queries: $e', level: LogLevel.warn);
    }
    // Steps 3 and 4 start with what was searched: CDS and/or rRNA, and
    // the searched markers selected.
    final rna = RegExp(r'^(12S|16S|rrn)', caseSensitive: false);
    for (final g in search.genes) {
      featureTypes.add(rna.hasMatch(g) ? 'rRNA' : 'CDS');
    }
    _markersCache = null;
    selectedMarkers.addAll(search.genes.where((g) => catalog?.byType.values.any((t) => t.containsKey(g)) ?? false));
    notifyListeners();
    if (!onePerTaxonMode) {
      setOnePerTaxon(true);
      log('Individual markers: one sample per species (FASTA header family_species; the longest sequence of each marker is kept)',
          level: LogLevel.success);
    }
  }

  Future<void> importGenomes(List<GenomeCandidate> selected, String taxon) =>
      _task('Importing ${selected.length} genome(s) of $taxon', (cancel) async {
        await _download([for (final c in selected) c.id], cancel, label: taxon);
        if (records.isNotEmpty) loadOptionsUsed.add(1);
      }, modal: TaskModal(ModalKind.ncbiFetch, taxon)..total = selected.length);

  void clearRecords() {
    records.clear();
    markerSearches.clear();
    loadOptionsUsed.removeAll(const [1, 2, 3]);
    selectedMarkers.removeWhere((m) => !(organized?.loci.contains(m) ?? false));
    duplicateChoices.clear();
    geneOverrides.clear();
    featureOverrides.clear();
    featureTypes.clear();
    taxonomyResolved = false;
    if (onePerTaxonMode) setOnePerTaxon(false);
    _rebuildCatalog();
    log('All GenBank records cleared');
    notifyListeners();
  }

  void removeRecord(GbRecord r) {
    records.remove(r);
    _rebuildCatalog();
    log('Removed ${r.organismName} (${r.accession})');
    notifyListeners();
  }

  /// Records without any protein-coding or rRNA gene SPLACE can name, with
  /// their feature counts (to explain why).
  List<({GbRecord record, int pcgs, int rrnas, int trnas})> get recordsWithoutMarkers {
    final cat = catalog;
    if (cat == null) return const [];
    return [
      for (var ri = 0; ri < records.length; ri++)
        if (!cat.recordGenes[ri].keys.any((fi) => const {'CDS', 'rRNA'}.contains(records[ri].features[fi].type)))
          (record: records[ri], pcgs: records[ri].featureCounts.pcgs, rrnas: records[ri].featureCounts.rrnas, trnas: records[ri].featureCounts.trnas),
    ];
  }

  void removeRecords(Iterable<GbRecord> gone) {
    final set = gone.toSet();
    records.removeWhere(set.contains);
    _rebuildCatalog();
    log('Removed ${set.length} record(s) without markers: ${set.map((r) => r.accession).join(', ')}');
    notifyListeners();
  }

  void editRecord(GbRecord r, {required String organism, TaxonomyRanks? ranks}) {
    r.editedOrganism = organism.trim() == r.organism ? null : organism.trim();
    if (ranks != null) {
      r.ranks = ranks;
      if (ranks.family.isNotEmpty) r.editedFamily = ranks.family;
    }
    _markersCache = null;
    log('Edited ${r.accession}: ${r.organismName}');
    notifyListeners();
  }

  /// Where "Fetch Taxonomy" reads the ranks from.
  TaxonomySource taxonomySource = TaxonomySource.gbif;

  void setTaxonomySource(TaxonomySource v) {
    taxonomySource = v;
    _changed(save: true);
  }

  String get _sourceLabel => switch (taxonomySource) {
        TaxonomySource.gbif => 'GBIF',
        TaxonomySource.ncbi => 'NCBI',
        TaxonomySource.salve => 'SALVE',
      };

  void _applyTaxon(GbRecord r, TaxonInfo t) {
    r.ranks = t.toRanks();
    if (t.family.isNotEmpty) r.editedFamily ??= t.family;
    if (t.status != null) r.conservation = t.status;
    if (t.salve != null) r.salve = t.salve;
  }

  /// "Fetch Taxonomy": ranks for every distinct species name from GBIF,
  /// NCBI Taxonomy (the records' taxon ids) or a SALVE table.
  Future<void> fetchTaxonomy() => _task('Fetching taxonomy from $_sourceLabel', (cancel) async {
        final source = _sourceLabel;
        final bySpecies = <String, List<GbRecord>>{};
        for (final r in records) {
          if (r.organismName.isNotEmpty) bySpecies.putIfAbsent(r.organismName, () => []).add(r);
        }
        final names = bySpecies.keys.toList();
        var ok = 0, failed = 0, done = 0;
        void miss(String name, [Object? e]) {
          failed++;
          log(
              e == null
                  ? (taxonomySource == TaxonomySource.salve ? 'No SALVE sheet for $name (not assessed as Brazilian fauna)' : 'No $source match for $name')
                  : '$source failed for $name: $e',
              source: source,
              level: LogLevel.warn);
        }

        switch (taxonomySource) {
          case TaxonomySource.gbif:
            for (var i = 0; i < names.length; i += 5) {
              if (cancel.isCancelled) break;
              await Future.wait(names.sublist(i, math.min(i + 5, names.length)).map((name) async {
                try {
                  final t = await ncbi.gbif(name);
                  if (t == null) {
                    miss(name);
                  } else {
                    ok++;
                    for (final r in bySpecies[name]!) {
                      _applyTaxon(r, t);
                    }
                  }
                } catch (e) {
                  miss(name, e);
                }
                _progress(++done, names.length, name);
              }));
            }
          case TaxonomySource.ncbi:
            // Taxon ids from /db_xref, a name search otherwise.
            final idOf = <String, String>{};
            for (final name in names) {
              if (cancel.isCancelled) break;
              var id = bySpecies[name]!.map((r) => r.taxId).whereType<String>().firstOrNull;
              try {
                id ??= await ncbi.ncbiTaxId(name);
              } catch (e) {
                log('NCBI Taxonomy search failed for $name: $e', source: source, level: LogLevel.warn);
              }
              if (id != null) idOf[name] = id;
              _progress(++done, names.length * 2, name);
            }
            final info = await ncbi.ncbiTaxonomy(idOf.values.toSet().toList());
            for (final name in names) {
              final t = info[idOf[name]];
              if (t == null) {
                miss(name);
              } else {
                ok++;
                for (final r in bySpecies[name]!) {
                  _applyTaxon(r, t);
                }
              }
              _progress(++done, names.length * 2, name);
            }
          case TaxonomySource.salve:
            // Brazilian fauna only: species without a sheet are reported.
            // Three species at a time (each sheet is ten requests).
            for (var i = 0; i < names.length; i += 3) {
              if (cancel.isCancelled) break;
              await Future.wait(names.sublist(i, math.min(i + 3, names.length)).map((name) async {
                try {
                  final t = await ncbi.salve(name, cancel: cancel);
                  if (t == null) {
                    miss(name);
                  } else {
                    ok++;
                    for (final r in bySpecies[name]!) {
                      _applyTaxon(r, t);
                    }
                  }
                } on CancelledException {
                  return;
                } catch (e) {
                  miss(name, e);
                }
                _progress(++done, names.length, name);
              }));
            }
        }
        taxonomyResolved = true;
        _markersCache = null;
        final summary = 'Taxonomy fetched: $ok/${names.length} species${failed > 0 ? ', $failed not found' : ''}'
            '${failed > 0 && taxonomySource == TaxonomySource.salve ? ' (SALVE only assesses species of the Brazilian fauna)' : ''}'
            '${names.length < records.length ? ' (${names.length} unique names for ${records.length} records)' : ''}';
        modal?.result = summary;
        log(summary, source: source, level: LogLevel.success);
      }, modal: TaskModal(ModalKind.taxonomy, _sourceLabel));

  /// Taxonomy of one name with the selected source (edit dialog).
  Future<TaxonInfo?> lookupTaxonomy(GbRecord r, String species) async {
    try {
      switch (taxonomySource) {
        case TaxonomySource.gbif:
          return await ncbi.gbif(species);
        case TaxonomySource.ncbi:
          final id = species == r.organism ? (r.taxId ?? await ncbi.ncbiTaxId(species)) : await ncbi.ncbiTaxId(species);
          return id == null ? null : (await ncbi.ncbiTaxonomy([id]))[id];
        case TaxonomySource.salve:
          return await ncbi.salve(species);
      }
    } catch (e) {
      log('$_sourceLabel failed for $species: $e', source: _sourceLabel, level: LogLevel.warn);
      return null;
    }
  }

  void skipTaxonomy() {
    taxonomyResolved = true;
    log('Taxonomy skipped; using the names in the GenBank records');
    notifyListeners();
  }

  void setColumnVisible(String column, bool visible) {
    visible ? hiddenColumns.remove(column) : hiddenColumns.add(column);
    _changed(save: true);
  }

  // ---------------------------------------------------------------------
  // 2-3. FASTA import and organization
  // ---------------------------------------------------------------------

  final List<FastaInput> inputs = [];
  LayoutSuggestion? suggestion;
  OrganizerConfig organizer =
      const OrganizerConfig(sample: FieldRule.fileName(delimiters: ''), locus: FieldRule.header());
  OrganizerResult? organized;
  bool organizing = false;
  int _organizeGeneration = 0;
  Timer? _organizeDebounce;
  final Set<String> _seenFastaLoci = {};

  /// New FASTA files were read: the user should confirm (in the "Sample
  /// and marker" dialog) which part of the headers is the sample and which
  /// the marker. Cleared by [applyOrganizer] or [dismissOrganizerReview].
  bool organizerReview = false;

  /// Why the last organization was not used (too many markers: the sample
  /// / marker split is most likely wrong).
  String? organizeProblem;

  /// The user confirmed the split at least once.
  bool organizerConfirmed = false;

  /// Folder of the FASTA files written per marker after the split was
  /// confirmed.
  String? lociExportDir;

  Future<void> addFiles(List<String> paths) => _task('Reading ${paths.length} FASTA file(s)', (_) async {
        if (paths.isEmpty) return;
        final r = await _loadInIsolate(paths, List.of(inputs));
        inputs.addAll(r.inputs);
        if (r.inputs.isNotEmpty) loadOptionsUsed.add(4);
        final n = r.inputs.fold<int>(0, (a, f) => a + f.records.length);
        log('Loaded ${r.inputs.length} FASTA file(s), $n sequence(s)');
        // Keep the user's split when files are added to a confirmed one.
        if (!organizerConfirmed) _applySuggestion(r.suggestion);
        organizerReview = r.inputs.isNotEmpty;
      }).then((_) => reorganize(immediate: true));

  /// Uses [config] (confirmed in the dialog), organizes the sequences and
  /// writes one FASTA per marker, every sample in it, to `01_loci`.
  /// Returns that folder (null when nothing was organized).
  Future<String?> applyOrganizer(OrganizerConfig config) async {
    organizer = config;
    organizerConfirmed = true;
    organizerReview = false;
    await reorganize(immediate: true, confirmed: true);
    if (organized == null || organized!.loci.isEmpty) return null;
    // Loci (target capture, UCEs…) are mostly non-coding: trimming does not
    // keep codons unless the user turns it on.
    if (lociMode) codonAwareTrimming = false;
    lociExportDir = await exportLoci(selectedOnly: false);
    notifyListeners();
    return lociExportDir;
  }

  void dismissOrganizerReview() {
    organizerReview = false;
    notifyListeners();
  }

  /// A split that gives about one marker per sequence is wrong (e.g. the
  /// whole header taken as the marker): showing thousands of one-sequence
  /// markers would freeze every step, so such a result is not used until
  /// the user fixes the split.
  static bool plausibleOrganization(OrganizerResult r) =>
      r.loci.length <= 20000 && !(r.records >= 200 && r.loci.length > r.records * 0.5);

  Future<void> addDirectory(String path) async {
    final files = await _listFiles(path, fastaExtensions);
    if (files.isEmpty) {
      log('No FASTA files found in $path', error: true);
      return;
    }
    await addFiles(files);
  }

  void _applySuggestion(LayoutSuggestion s) {
    suggestion = s;
    organizer = s.config;
    log('Layout guess: ${s.kind.name} (confidence ${(s.confidence * 100).round()}%)');
    for (final e in s.evidence) {
      log('  • $e');
    }
  }

  void clearInputs() {
    inputs.clear();
    loadOptionsUsed.remove(4);
    organized = null;
    suggestion = null;
    organizerReview = false;
    organizerConfirmed = false;
    organizeProblem = null;
    lociExportDir = null;
    selectedMarkers.removeWhere((m) => !(catalog?.byType.values.any((t) => t.containsKey(m)) ?? false));
    _seenFastaLoci.clear();
    _markersCache = null;
    notifyListeners();
  }

  void updateOrganizer(OrganizerConfig config) {
    organizer = config;
    notifyListeners();
    reorganize();
  }

  /// Recomputes the locus × sample matrix in an isolate. Rapid edits are
  /// debounced and stale results discarded, so typing a regex never
  /// freezes the UI.
  Future<void> reorganize({bool immediate = false, bool confirmed = false}) async {
    _organizeDebounce?.cancel();
    if (!immediate) {
      _organizeDebounce = Timer(const Duration(milliseconds: 250), () => reorganize(immediate: true));
      return;
    }
    if (inputs.isEmpty) return;
    final generation = ++_organizeGeneration;
    organizing = true;
    notifyListeners();
    try {
      final result = await _organizeInIsolate(List.of(inputs), organizer);
      if (generation != _organizeGeneration) return;
      if (!plausibleOrganization(result)) {
        organized = null;
        _markersCache = null;
        organizeProblem = '${result.loci.length} markers for ${result.records} sequences';
        log('The sample / marker split gives ${result.loci.length} markers for ${result.records} sequences; '
            'choose which part of the headers is the sample and which the marker', level: LogLevel.warn);
        if (!confirmed) organizerReview = true;
        return;
      }
      organizeProblem = null;
      organized = result;
      _markersCache = null;
      // New FASTA loci start selected; the user's later choices are kept.
      for (final l in result.loci) {
        if (_seenFastaLoci.add(l)) selectedMarkers.add(l);
      }
    } catch (e) {
      log('Organizing failed: $e', error: true);
    } finally {
      if (generation == _organizeGeneration) {
        organizing = false;
        notifyListeners();
      }
    }
  }

  /// Writes one FASTA per FASTA-derived locus plus a provenance table.
  Future<String?> exportLoci({String? targetDir, bool selectedOnly = true}) async {
    final result = organized;
    if (result == null) return null;
    final out = targetDir ?? dir('01_loci');
    final loci = [for (final l in result.loci) if (!selectedOnly || selectedMarkers.contains(l)) l];
    // Hundreds of markers: built and written in one go off the UI thread.
    final n = await _writeLociInIsolate(result, loci, out);
    log('Wrote $n locus FASTA file(s) and splace_locus_mapping.tsv to $out');
    return out;
  }

  // ---------------------------------------------------------------------
  // 3-4. Feature types and markers
  // ---------------------------------------------------------------------

  final Set<String> featureTypes = {};
  final Set<String> selectedMarkers = {};
  final Map<String, int> duplicateChoices = {};
  final Map<String, String> geneOverrides = {};
  bool unknownGenesReviewed = true;

  void toggleFeatureType(String type, bool on) {
    on ? featureTypes.add(type) : featureTypes.remove(type);
    _markersCache = null;
    notifyListeners();
  }

  /// Gene names mapped by hand in the "unknown genes" dialog.
  // -------------------------------------------------------------------------
  // Species identification (GenBank BLAST + BOLD), before the taxonomy
  // -------------------------------------------------------------------------

  /// Last identification, by accession.
  Map<String, RecordIdentification> speciesIds = {};

  /// Genes used in the last identification.
  List<String> speciesIdGenes = [];

  /// Sources that failed in the last identification, with the error.
  Map<IdSource, String> speciesIdFailures = {};

  /// Standard CDS/rRNA markers found in the records, with how many records
  /// carry each, for the identification setup.
  Map<String, int> get identificationMarkers {
    final cat = catalog;
    if (cat == null) return const {};
    final out = <String, int>{};
    for (final type in const ['CDS', 'rRNA']) {
      for (final g in cat.byType[type]?.values ?? const <GeneSummary>[]) {
        if (isKnownGeneName(g.name)) out[g.name] = g.byRecord.length;
      }
    }
    return Map.fromEntries(out.entries.toList()..sort((a, b) => b.value != a.value ? b.value - a.value : a.key.compareTo(b.key)));
  }

  /// Nucleotide sequence of [gene] in each record (first copy).
  Map<String, Map<String, String>> _geneSequences(List<String> genes) {
    final cat = catalog!;
    final out = <String, Map<String, String>>{};
    for (var ri = 0; ri < records.length; ri++) {
      final r = records[ri];
      for (final e in cat.recordGenes[ri].entries) {
        if (!genes.contains(e.value)) continue;
        final byGene = out.putIfAbsent(r.accession, () => {});
        if (byGene.containsKey(e.value)) continue;
        final seq = extractLocation(r.sequence, r.features[e.key].location);
        if (seq.isNotEmpty) byGene[e.value] = seq;
      }
    }
    return out;
  }

  /// BLAST/BOLD answers already received (see [SpeciesIdentifier.cacheKey]),
  /// kept on disk: running again with another gene only sends what is new.
  Map<String, List<IdHit>>? _idCache;

  File get _idCacheFile => File(p.join(documentsDir, 'cache', 'species_id.json'));

  Map<String, List<IdHit>> get _speciesIdCache {
    if (_idCache != null) return _idCache!;
    final cache = <String, List<IdHit>>{};
    try {
      if (_idCacheFile.existsSync()) {
        final json = jsonDecode(_idCacheFile.readAsStringSync()) as Map<String, dynamic>;
        for (final e in json.entries) {
          if (SpeciesIdentifier.staleCacheKey(e.key)) continue;
          cache[e.key] = [for (final h in e.value as List) IdHit.fromJson((h as Map).cast<String, dynamic>())];
        }
      }
    } catch (e) {
      debugPrint('Species identification cache could not be read: $e');
    }
    return _idCache = cache;
  }

  /// Hits kept with their aligned sequences in the cache, per query (the
  /// others keep the numbers only, so the file stays small).
  static const _cachedAlignments = 15;

  Future<void> _saveSpeciesIdCache() async {
    final cache = _idCache;
    if (cache == null) return;
    await _idCacheFile.parent.create(recursive: true);
    await _idCacheFile.writeAsString(jsonEncode({
      for (final e in cache.entries)
        e.key: [
          for (final (i, h) in (List.of(e.value)..sort((a, b) => b.identity.compareTo(a.identity))).indexed) h.toJson(withAlignment: i < _cachedAlignments),
        ],
    }));
  }

  /// The alignment of a BOLD hit with [query]: BOLD returns where they
  /// align, not the aligned sequences, so the hit's sequence is fetched
  /// from the BOLD portal (by process id) and that stretch is aligned to
  /// the query here. Kept in the hit (and the cache) once done.
  Future<void> loadBoldAlignment(IdHit hit, String query) async {
    if (hit.hasAlignment || hit.source != IdSource.bold) return;
    final rows = await boldPortal.records(['ids:processid:${hit.accession}'], const {});
    if (rows.isEmpty) throw StateError('BOLD has no public sequence for ${hit.accession}');
    final code = hit.title.split(' · ').elementAtOrNull(1) ?? '';
    final rec = rows.where((r) => r.markerCode == code).firstOrNull ?? (rows.toList()..sort((a, b) => b.sequence.length - a.sequence.length)).first;
    final subject = rec.sequence.toUpperCase();
    String cut(String s, int from, int to) {
      final a = (from - 1).clamp(0, s.length), b = to.clamp(a, s.length);
      return from > 0 && to >= from ? s.substring(a, b) : s;
    }

    final q = cut(query.toUpperCase(), hit.qFrom, hit.qTo);
    var h = hit.reversed ? reverseComplement(cut(subject, hit.sTo, hit.sFrom)) : cut(subject, hit.sFrom, hit.sTo);
    if (h.isEmpty) h = subject;
    final (qa, ha) = await Isolate.run(() => alignPair(q, h));
    hit
      ..qseq = qa
      ..hseq = ha;
    await _saveSpeciesIdCache();
  }

  /// Parameters of BLAST and the BOLD ID engine (saved with the settings).
  SpeciesIdOptions speciesIdOptions = const SpeciesIdOptions();

  void setSpeciesIdOptions(SpeciesIdOptions o) {
    speciesIdOptions = o;
    _saveSettings();
    notifyListeners();
  }

  /// The marker name a barcode locus is sent under: COI (any of its names),
  /// the plastid barcodes (rbcL, matK) by their SynGenes name, else as is
  /// (ITS…).
  static String barcodeGene(String locus) {
    if (isCoiLocus(locus)) return 'COI';
    final cp = standardizeGeneName(locus, null, OrganelleType.cp);
    return cp == 'rbcL' || cp == 'matK' ? cp! : (RegExp(r'^its', caseSensitive: false).hasMatch(locus.trim()) ? 'ITS' : locus);
  }

  /// How many queries for [genes] are already in the cache, per source
  /// (GenBank, BOLD) and how many would be sent.
  ({int cached, int total}) speciesIdCacheCoverage(List<String> genes) {
    final cache = _speciesIdCache;
    var cached = 0, total = 0;
    final seqs = catalog == null ? const <String, Map<String, String>>{} : _geneSequences(genes);
    for (final e in seqs.entries) {
      for (final g in e.value.entries) {
        if (g.value.length < 100) continue;
        for (final s in [IdSource.genbank, if (speciesIdOptions.boldAccepts(g.key)) IdSource.bold]) {
          total++;
          if (cache.containsKey(SpeciesIdentifier.cacheKey(s, e.key, g.key, g.value, speciesIdOptions))) cached++;
        }
      }
    }
    return (cached: cached, total: total);
  }

  int get speciesIdCacheSize => _speciesIdCache.length;

  Future<void> clearSpeciesIdCache() async {
    _idCache = {};
    try {
      if (_idCacheFile.existsSync()) await _idCacheFile.delete();
    } catch (_) {}
    log('Species identification cache cleared');
    notifyListeners();
  }

  /// Sends [genes] of every record to GenBank (BLAST) and COI to BOLD;
  /// answers already in the cache are reused.
  Future<bool> identifySpecies(List<String> genes) async {
    if (catalog == null || genes.isEmpty) return false;
    final o = speciesIdOptions;
    final toBold = genes.any(o.boldAccepts);
    final ok = await _task('Identifying species with GenBank${toBold ? ' and BOLD' : ''} (${genes.join(', ')})', (cancel) async {
      // The bar is the mean of the sources in use (BOLD only with its markers).
      final share = {IdSource.genbank: 0.0, if (toBold) IdSource.bold: 0.0};
      final id = SpeciesIdentifier(
        options: o,
        onStatus: (m) {
          modal?.current = m;
          notifyListeners();
        },
        onProgress: (source, f) {
          share[source] = f;
          final mean = share.values.fold(0.0, (a, v) => a + v) / share.length;
          _progress((mean * 1000).round(), 1000);
        },
      );
      try {
        log('Species identification: ${records.length} record(s), genes ${genes.join(', ')}; '
            'GenBank BLAST (${o.blastProgram}, ${o.blastDatabase}, e ≤ ${o.expect}, ${o.hitlistSize} hits${o.entrezQuery.trim().isEmpty ? '' : ', ${o.entrezQuery.trim()}'})'
            '${toBold ? ' and BOLD ID engine (${o.boldDatabase}, ≥ ${(o.boldMinIdentity * 100).round()} %, ${o.boldMaxHits} hits)' : ''}');
        final cache = _speciesIdCache;
        try {
          speciesIds = await id.identify(records, _geneSequences(genes),
              cancel: cancel, log: (m) => log(m, level: LogLevel.warn), cache: cache);
        } finally {
          // Answers received before a cancel or failure are kept too.
          await _saveSpeciesIdCache();
        }
        modal?.current = 'GBIF synonyms';
        await _mergeIdSynonyms(speciesIds.values);
        if (id.reused > 0) log('Species identification: ${id.reused} quer${id.reused == 1 ? 'y' : 'ies'} reused from the cache, ${id.sent} sent');
        speciesIdGenes = genes;
        speciesIdFailures = Map.of(id.failures);
        final disagree = speciesIds.values.where((r) => r.disagrees).length;
        modal?.result = '$disagree of ${speciesIds.length} record(s) may be misidentified'
            '${id.reused == 0 ? '' : ' · ${id.reused} result(s) reused from the cache'}';
        log(modal!.result!, level: disagree == 0 ? LogLevel.success : LogLevel.warn);
        return true;
      } finally {
        id.close();
      }
    }, modal: TaskModal(ModalKind.speciesId, genes.join(', ')));
    return ok ?? false;
  }

  /// Writes every record as a corrected GenBank file (one per record,
  /// `<accession>_<Genus_species>.gb`) in [folder] and in the project's
  /// `00_genbank_corrected`: the organism as identified or edited, the
  /// lineage from the taxonomy, and the gene names standardized by SynGenes
  /// (the original name kept in a /note). [version] is SPLACE's, for the
  /// COMMENT. Returns the number of files written.
  Future<int> exportCorrectedGenBank(String folder, {required String version}) async {
    final cat = catalog;
    final what = <String>[
      if (records.any((r) => r.editedOrganism != null)) 'organism name from the species identification or edited by the user',
      if (taxonomyResolved) 'lineage from $_sourceLabel (dataFishing method)',
      if (cat != null) 'gene names standardized with SynGenes',
    ];
    final note = 'Corrected with SPLACE v$version${what.isEmpty ? '' : ': ${what.join('; ')}'}.';
    final project = dir('00_genbank_corrected');
    var n = 0;
    final used = <String>{};
    for (var ri = 0; ri < records.length; ri++) {
      final r = records[ri];
      final text = writeGenBank(r, geneNames: cat == null || ri >= cat.recordGenes.length ? const {} : cat.recordGenes[ri], note: note);
      String safe(String v) => v.trim().replaceAll(RegExp(r'[^A-Za-z0-9._-]+'), '_').replaceAll(RegExp(r'^_+|_+$'), '');
      var name = [safe(r.accession), safe(r.organismName)].where((x) => x.isNotEmpty).join('_');
      if (name.isEmpty) name = 'record_${ri + 1}';
      while (!used.add(name)) {
        name = '${name}_${ri + 1}';
      }
      for (final d in {folder, project}) {
        await File(p.join(d, '$name.gb')).writeAsString(text);
      }
      n++;
    }
    log('Corrected GenBank files: $n written to $folder (copy in $project)', level: LogLevel.success);
    return n;
  }

  // -------------------------------------------------------------------------
  // DNA barcoding (FASTA): identification and BOLD submission
  // -------------------------------------------------------------------------

  /// The FASTA input is DNA barcodes (one specimen per sequence): step 2
  /// offers the identification and the BOLD submission spreadsheet.
  bool barcodingMode = false;

  void setBarcodingMode(bool v) {
    barcodingMode = v;
    notifyListeners();
  }

  /// Marker identified (a locus of the organized FASTA).
  String? barcodeLocus;

  List<String> get barcodeLoci => organized?.loci ?? const [];

  /// [barcodeLocus], or the COI-like locus, or the first one.
  String? get barcodeLocusShown {
    final loci = barcodeLoci;
    if (loci.isEmpty) return null;
    if (barcodeLocus != null && loci.contains(barcodeLocus)) return barcodeLocus;
    return loci.where(isCoiLocus).firstOrNull ?? loci.first;
  }

  void setBarcodeLocus(String? v) {
    barcodeLocus = v;
    notifyListeners();
  }

  /// COI under any of its names (SynGenes): BOLD identifies only COI.
  static bool isCoiLocus(String name) {
    if (RegExp(r'^(coi|cox1|co1|cox-?i|mt-?co1|coxi)$', caseSensitive: false).hasMatch(name.trim())) return true;
    return standardizeGeneName(name, null, OrganelleType.mt) == 'COI';
  }

  /// Sample → sequence (gaps removed) of the barcode locus.
  Map<String, String> get barcodeSequences {
    final locus = barcodeLocusShown;
    final cells = locus == null ? null : organized?.byLocus[locus];
    if (cells == null) return const {};
    return {for (final e in cells.entries) e.key: e.value.record.sequence.replaceAll(RegExp(r'[-.\s]'), '')};
  }

  /// Last identification of the barcodes, by sample.
  Map<String, RecordIdentification> barcodeIds = {};

  /// Accepted names (GBIF) already looked up: name → accepted (null: GBIF
  /// does not know it).
  final Map<String, String?> _acceptedNames = {};

  /// Synonyms in the identification: the record's name and the best
  /// species of each record are looked up in GBIF, and the hits that came
  /// under a synonym count for the accepted species (Bufo marinus and
  /// Chaunus marinus hits for Rhinella marina). Without GBIF the names stay
  /// as they are.
  Future<void> _mergeIdSynonyms(Iterable<RecordIdentification> ids) async {
    final names = <String>{
      for (final r in ids) ...[
        if (isSpeciesName(binomialOf(r.currentName))) binomialOf(r.currentName),
        for (final c in r.shortlist(perSource: 6)) if (isSpeciesName(binomialOf(c.species))) binomialOf(c.species),
      ],
    }.where((n) => !_acceptedNames.containsKey(n)).toList();
    var failed = 0;
    // Four at a time: GBIF answers quickly and has no key.
    for (var i = 0; i < names.length; i += 4) {
      await Future.wait([
        for (final n in names.skip(i).take(4))
          ncbi.gbifAccepted(n).then((a) => _acceptedNames[n] = a).catchError((Object _) {
            failed++;
            return null;
          }),
      ]);
    }
    if (failed > 0) log('GBIF could not be reached for $failed name(s): their synonyms are not grouped', level: LogLevel.warn);
    final map = <String, String>{
      for (final e in _acceptedNames.entries)
        if (e.value != null) e.key.toLowerCase(): e.value!,
    };
    var merged = 0;
    for (final r in ids) {
      final before = r.candidates.length;
      r.mergeSynonyms(map);
      merged += before - r.candidates.length;
    }
    final renamed = ids.where((r) => r.currentAccepted != null).length;
    log('Species identification: synonyms checked in GBIF (${names.length} name(s)); $merged synonym match(es) joined to the accepted species'
        '${renamed == 0 ? '' : ', $renamed record(s) named with a synonym'}');
  }

  /// Name chosen for each sample (empty: not identified).
  final Map<String, String> barcodeNames = {};

  /// Taxonomy (GBIF) of the chosen names, by name.
  final Map<String, TaxonInfo> barcodeTaxa = {};

  /// Sends every barcode to GenBank BLAST and, for COI, to BOLD; each
  /// sample gets the best match at a same-species identity (98 %) as its
  /// name, which the user can change.
  Future<bool> identifyBarcodes() async {
    final seqs = barcodeSequences;
    final locus = barcodeLocusShown;
    if (seqs.isEmpty || locus == null) return false;
    final gene = barcodeGene(locus);
    final o = speciesIdOptions;
    final toBold = o.boldAccepts(gene);
    final ok = await _task('Identifying ${seqs.length} barcode(s) with GenBank${toBold ? ' and BOLD' : ''}', (cancel) async {
      final share = {IdSource.genbank: 0.0, if (toBold) IdSource.bold: 0.0};
      final id = SpeciesIdentifier(
        options: o,
        onStatus: (m) {
          modal?.current = m;
          notifyListeners();
        },
        onProgress: (source, f) {
          share[source] = f;
          _progress((share.values.fold(0.0, (a, v) => a + v) / share.length * 1000).round(), 1000);
        },
      );
      try {
        // Each sample as a minimal record: identify() needs only its id
        // and name.
        final records = [
          for (final sample in seqs.keys)
            GbRecord(accession: sample, organism: barcodeNames[sample] ?? '', lineage: '', definition: '', features: const [], sequence: seqs[sample]!),
        ];
        final cache = _speciesIdCache;
        try {
          barcodeIds = await id.identify(records, {for (final e in seqs.entries) e.key: {gene: e.value}},
              cancel: cancel, log: (m) => log(m, level: LogLevel.warn), cache: cache);
        } finally {
          await _saveSpeciesIdCache();
        }
        modal?.current = 'GBIF synonyms';
        await _mergeIdSynonyms(barcodeIds.values);
        speciesIdFailures = Map.of(id.failures);
        var named = 0;
        for (final e in barcodeIds.entries) {
          final best = e.value.candidates.firstOrNull;
          if ((barcodeNames[e.key] ?? '').isEmpty && best != null && best.bestIdentity >= sameSpeciesIdentity) {
            barcodeNames[e.key] = binomialOf(best.species);
            named++;
          }
        }
        modal?.result = '$named of ${seqs.length} barcode(s) matched a species at ≥ ${sameSpeciesIdentity.toStringAsFixed(0)} %';
        log(modal!.result!, level: LogLevel.success);
        return true;
      } finally {
        id.close();
      }
    }, modal: TaskModal(ModalKind.speciesId, gene));
    return ok ?? false;
  }

  // -------------------------------------------------------------------------
  // Identification of tips from the tree viewer
  // -------------------------------------------------------------------------

  /// Identifications started from the tree, by marker and tip.
  final Map<String, Map<String, RecordIdentification>> tipIds = {};

  /// Sequence of [tip] in [marker] (its alignment without gaps, else its
  /// raw sequence), or null when the tip has none. Tip names are the
  /// alignment headers (also with the characters a tree program changes).
  String? tipSequence(String marker, String tip) {
    final l = loci[marker];
    if (l == null || l.isProtein) return null;
    String norm(String s) => s.replaceAll(RegExp(r'[^A-Za-z0-9.]+'), '_');
    for (final list in [l.aligned, l.raw]) {
      if (list == null) continue;
      final r = list.where((x) => x.header == tip).firstOrNull ?? list.where((x) => norm(x.header) == norm(tip)).firstOrNull;
      if (r == null) continue;
      final seq = r.sequence.replaceAll(RegExp(r'[-.?\s]'), '');
      if (seq.isNotEmpty) return seq;
    }
    return null;
  }

  /// DNA markers with a sequence of at least 100 bases for some of [tips],
  /// with how many tips have it (most first).
  List<({String marker, int tips})> tipMarkers(List<String> tips) {
    final out = [
      for (final l in loci.values)
        if (!l.isProtein)
          (marker: l.name, tips: tips.where((t) => (tipSequence(l.name, t)?.length ?? 0) >= 100).length),
    ].where((m) => m.tips > 0).toList()
      ..sort((a, b) => b.tips != a.tips ? b.tips - a.tips : a.marker.compareTo(b.marker));
    return out;
  }

  /// Identifies [tips] with their [marker] sequence: GenBank BLAST and, for
  /// the markers BOLD identifies, BOLD (the options of step 2; answers in
  /// the cache are reused). The record a tip came from is not evidence.
  /// [onProgress] gets the share done and what is being done (the tree
  /// editor covers the progress window). Null when it could not run.
  Future<Map<String, RecordIdentification>?> identifyTips(List<String> tips, String marker, {void Function(double share, String status)? onProgress}) async {
    final seqs = {for (final t in tips) if (tipSequence(marker, t) case final s? when s.length >= 100) t: s};
    if (seqs.isEmpty) return null;
    final gene = barcodeGene(marker);
    final o = speciesIdOptions;
    final toBold = o.boldAccepts(gene);
    final info = tipInfo();
    final byAccession = {for (final r in records) r.accession: r};
    var status = '';
    return _task('Identifying ${seqs.length} tip(s) of the tree with $marker (GenBank${toBold ? ' and BOLD' : ''})', (cancel) async {
      final share = {IdSource.genbank: 0.0, if (toBold) IdSource.bold: 0.0};
      double mean() => share.values.fold(0.0, (a, v) => a + v) / share.length;
      final id = SpeciesIdentifier(
        options: o,
        onStatus: (m) {
          status = m;
          onProgress?.call(mean(), m);
        },
        onProgress: (source, f) {
          share[source] = f;
          _progress((mean() * 1000).round(), 1000);
          onProgress?.call(mean(), status);
        },
      );
      try {
        // Each tip as a minimal record named after it, with the record it
        // came from (and that one's INSDC sources) as "itself".
        final queryRecords = [
          for (final tip in seqs.keys)
            () {
              final i = info[tip];
              final r = i == null ? null : byAccession[i.accession];
              final name = (i?.species ?? '').isNotEmpty ? i!.species : _plainTip(tip);
              return GbRecord(
                accession: tip,
                organism: name,
                lineage: '',
                definition: '',
                features: const [],
                sequence: seqs[tip]!,
                sourceAccessions: [if (r != null) r.accession, ...?r?.sourceAccessions],
              );
            }(),
        ];
        log('Tree identification: ${seqs.length} tip(s), $marker; GenBank BLAST (${o.blastProgram}, ${o.blastDatabase})${toBold ? ' and BOLD (${o.boldDatabase})' : ''}');
        final cache = _speciesIdCache;
        Map<String, RecordIdentification> ids;
        try {
          ids = await id.identify(queryRecords, {for (final e in seqs.entries) e.key: {gene: e.value}},
              cancel: cancel, log: (m) => log(m, level: LogLevel.warn), cache: cache);
        } finally {
          await _saveSpeciesIdCache();
        }
        onProgress?.call(1, 'GBIF synonyms');
        await _mergeIdSynonyms(ids.values);
        speciesIdFailures = Map.of(id.failures);
        tipIds.putIfAbsent(marker, () => {}).addAll(ids);
        final disagree = ids.values.where((r) => r.disagrees).length;
        log('Tree identification ($marker): $disagree of ${ids.length} tip(s) may be another species'
            '${id.reused == 0 ? '' : ' · ${id.reused} reused from the cache'}', level: disagree == 0 ? LogLevel.success : LogLevel.warn);
        return ids;
      } finally {
        id.close();
      }
    });
  }

  /// A tip name as words: `NC_0123.1_Family_Genus_species` → the part
  /// after the accession, with spaces.
  static String _plainTip(String tip) =>
      tip.replaceFirst(RegExp(r'^[A-Z]{1,6}_?\d+(\.\d+)?_'), '').replaceAll('_', ' ').trim();

  void setBarcodeName(String sample, String name) {
    barcodeNames[sample] = name.trim();
    notifyListeners();
  }

  /// Ranks (phylum to family) of the chosen names from GBIF, for the BOLD
  /// sheet (dataFishing method).
  Future<void> fetchBarcodeTaxonomy() => _task('Fetching taxonomy of the barcodes from GBIF', (cancel) async {
        final names = {for (final n in barcodeNames.values) if (n.trim().isNotEmpty) n.trim()}.toList();
        var done = 0;
        for (final name in names) {
          if (cancel.isCancelled) break;
          try {
            final t = await ncbi.gbif(name);
            if (t != null) barcodeTaxa[name] = t;
          } catch (e) {
            log('GBIF failed for $name: $e', level: LogLevel.warn, source: 'GBIF');
          }
          _progress(++done, names.length, name);
        }
        log('Taxonomy of ${barcodeTaxa.length} of ${names.length} name(s) from GBIF', level: LogLevel.success);
      });

  /// Writes the BOLD5 single-page specimen sheet ([path], .xlsx) and, next
  /// to it, the sequences with the sample ids as headers (`.fas`) for
  /// BOLD's sequence upload. Returns the number of specimens.
  Future<int> exportBoldSubmission(String path) async {
    final seqs = barcodeSequences;
    final locus = barcodeLocusShown ?? '';
    final specimens = <BoldSpecimen>[];
    for (final sample in seqs.keys) {
      final name = (barcodeNames[sample] ?? '').trim();
      final words = name.split(RegExp(r'\s+'));
      final t = barcodeTaxa[name];
      final id = barcodeIds[sample];
      final best = id?.candidates.where((c) => sameSpecies(c.species, name)).firstOrNull;
      String pct(double v) => '${v.toStringAsFixed(2)} %';
      final evidence = [
        if (best != null && best.bold != null) 'BOLD ID Engine ${pct(best.bold!.identity)}',
        if (best != null && best.inGenbank) 'GenBank BLAST ${pct(best.genbankIdentity)}',
      ];
      specimens.add(BoldSpecimen(
        sampleId: sample,
        phylum: t?.phylum ?? '',
        klass: t?.klass ?? '',
        order: t?.order ?? '',
        family: t?.family ?? '',
        genus: name.isEmpty ? '' : words.first,
        species: words.length >= 2 ? name : '',
        identificationMethod: evidence.isEmpty ? '' : 'DNA barcode ($locus) similarity: ${evidence.join('; ')} (SPLACE)',
        taxonomyNotes: evidence.isEmpty ? '' : 'Best match ${binomialOf(best!.species)}${t == null ? '' : '; ranks from GBIF'}',
      ));
    }
    final target = path.toLowerCase().endsWith('.xlsx') ? path : '$path.xlsx';
    await File(target).writeAsBytes(buildXlsx(boldSheet, boldSheetRows(specimens)));
    final fas = '${target.substring(0, target.length - 5)}_sequences.fas';
    await File(fas).writeAsString('${boldFasta(seqs)}\n');
    log('BOLD submission: ${specimens.length} specimen(s) in $target; sequences in $fas', level: LogLevel.success);
    return specimens.length;
  }

  /// Renames the records the user chose to (accession → species).
  void applySpeciesChoices(Map<String, String> choices) {
    for (final e in choices.entries) {
      final r = records.where((r) => r.accession == e.key).firstOrNull;
      if (r == null || r.organismName == e.value) continue;
      log('Species of ${r.accession}: ${r.organismName} → ${e.value} (GenBank/BOLD identification)');
      r.editedOrganism = e.value == r.organism ? null : e.value;
    }
    _markersCache = null;
    notifyListeners();
  }

  /// Names fixed for single features (see [featureKey]).
  final Map<String, String> featureOverrides = {};

  /// Applies the gene review: [mapping] renames raw names everywhere,
  /// [perFeature] fixes single features.
  void applyGeneOverrides(Map<String, String> mapping, {Map<String, String> perFeature = const {}}) {
    geneOverrides.addAll({
      for (final e in mapping.entries)
        if (e.value.isNotEmpty) e.key: e.value,
    });
    featureOverrides.addAll(perFeature);
    for (final e in perFeature.entries) {
      log('Gene of ${e.key.replaceFirst('#', ' feature ')} → ${e.value}');
    }
    unknownGenesReviewed = true;
    _rebuildCatalog();
    for (final e in mapping.entries) {
      if (e.value.isNotEmpty) log('Gene name ${e.key} → ${e.value}');
    }
    notifyListeners();
  }

  void dismissUnknownGenes() {
    unknownGenesReviewed = true;
    notifyListeners();
  }

  Map<String, MarkerSummary>? _markersCache;

  /// Markers of step 4: GenBank genes of the selected feature types plus
  /// FASTA loci, merged by name.
  Map<String, MarkerSummary> get markers {
    final cached = _markersCache;
    if (cached != null) return cached;
    final out = <String, MarkerSummary>{};
    final cat = catalog;
    if (cat != null) {
      for (final g in cat.visible(featureTypes).values) {
        final m = out.putIfAbsent(g.name, () => MarkerSummary(g.name))..sources.addAll(g.types);
        for (final e in g.byRecord.entries) {
          m.lengths['gb:${e.key}'] = [for (final h in e.value) h.length];
        }
      }
    }
    final org = organized;
    if (org != null) {
      for (final l in org.loci) {
        final m = out.putIfAbsent(l, () => MarkerSummary(l))..sources.add('FASTA');
        for (final e in org.byLocus[l]!.entries) {
          m.lengths['fa:${e.key}'] = [e.value.record.sequence.replaceAll(RegExp('[-?]'), '').length];
        }
      }
    }
    return _markersCache = Map.fromEntries(out.entries.toList()..sort((a, b) => naturalCompare(a.key, b.key)));
  }

  /// Rows of the heatmaps (GenBank records sorted by organism, then FASTA
  /// samples).
  List<SampleRef> get samples {
    final order = List.generate(records.length, (i) => i)
      ..sort((a, b) => records[a].organismName.compareTo(records[b].organismName));
    return [
      for (final i in order) SampleRef('gb:$i', records[i].organismName, detail: records[i].accession, italic: true),
      for (final s in organized?.samples ?? const <String>[]) SampleRef('fa:$s', s),
    ];
  }

  int get totalSamples => records.length + (organized?.samples.length ?? 0);

  /// Why [gene] is missing from sample [s] (step 4's missing panel).
  MissingCause missingCause(String gene, SampleRef s) {
    final cat = catalog;
    if (!s.key.startsWith('gb:') || cat == null) return const MissingCause(MissingReason.notInFiles);
    return diagnoseMissing(
      records: records,
      index: int.parse(s.key.substring(3)),
      gene: gene,
      catalog: cat,
      types: featureTypes,
      overrides: geneOverrides,
      typical: typicalRecordLength(records),
    );
  }

  /// How the copies of a duplicated gene in record [ri] relate.
  CopyRelation copyRelationOf(int ri, List<GeneHit> hits) {
    final r = records[ri];
    return copyRelation([for (final h in hits) extractLocation(r.sequence, r.features[h.feature].location)]);
  }

  void toggleMarker(String name, [bool? on]) {
    final value = on ?? !selectedMarkers.contains(name);
    value ? selectedMarkers.add(name) : selectedMarkers.remove(name);
    notifyListeners();
  }

  void selectAllMarkers() {
    selectedMarkers.addAll(markers.keys);
    notifyListeners();
  }

  void selectNoMarkers() {
    selectedMarkers.clear();
    notifyListeners();
  }

  /// "Shared Only": markers present in every sample.
  void selectSharedMarkers() {
    selectedMarkers
      ..clear()
      ..addAll([for (final m in markers.values) if (m.samples >= totalSamples) m.name]);
    notifyListeners();
  }

  /// "Default to Phylogeny": the 13 mitochondrial PCGs or the plastid set.
  void selectDefaultMarkers() {
    featureTypes.add('CDS');
    _markersCache = null;
    final defaults = organelle == OrganelleType.cp ? cpDefaultGenes : mtDefaultGenes;
    selectedMarkers
      ..clear()
      ..addAll(defaults.where(markers.containsKey));
    notifyListeners();
  }

  /// Selected markers that exist (a marker may vanish when its feature
  /// type is switched off).
  List<String> get activeMarkers => [for (final m in markers.keys) if (selectedMarkers.contains(m)) m];

  void setDuplicateChoice(GbRecord r, String gene, int index) {
    duplicateChoices[duplicateChoiceKey(r, gene)] = index;
    notifyListeners();
  }

  // ---------------------------------------------------------------------
  // 5. Alignment
  // ---------------------------------------------------------------------

  final Map<String, LocusState> loci = {};
  HeaderTemplate header = const HeaderTemplate();

  /// Individual-marker data (separate records per marker): one sample per
  /// taxon, the longest sequence of each marker kept. Turned on by a
  /// marker search; the FASTA header then names the taxon only.
  bool onePerTaxonMode = false;

  void setOnePerTaxon(bool v) {
    onePerTaxonMode = v;
    if (v) {
      header = taxonHeader;
    } else if (header.fields.join(',') == taxonHeader.fields.join(',')) {
      header = const HeaderTemplate();
    }
    headerConfirmed = false;
    notifyListeners();
  }
  bool headerConfirmed = false;
  SeqMode seqMode = SeqMode.nucleotide;
  bool allowPseudo = false;
  bool allowInternalStops = false;
  bool allowBacktranslationWarnings = false;
  MafftOptions mafft = const MafftOptions(threads: 4);
  TrimalOptions trimal = const TrimalOptions();
  String clipkitMode = 'gappy';

  /// How DNA alignments are trimmed: trimAl keeping whole codons
  /// (`-backtrans`), ClipKIT keeping whole codons (`--codon`), or trimAl
  /// column by column.
  TrimTool trimTool = TrimTool.trimalCodon;

  /// Coding genes keep whole codons when trimmed.
  bool get codonAwareTrimming => trimTool != TrimTool.trimalColumns;
  set codonAwareTrimming(bool v) => trimTool = v ? (trimTool == TrimTool.trimalColumns ? TrimTool.trimalCodon : trimTool) : TrimTool.trimalColumns;

  /// "Alignment complete — trim poorly-aligned regions?" card.
  bool trimPrompt = false;
  bool trimSkipped = false;

  /// The user chose to trim: the trimming parameters (5d) are shown.
  bool trimOpen = false;
  /// The coding genes' saturation results, summed up per codon position
  /// (the suggestion of step 6).
  PooledSaturation? get pooledSaturation => PooledSaturation.of(loci.values.where((l) => l.included).map((l) => l.saturation).whereType<SaturationReport>());
  ConcatOptions concat = const ConcatOptions();
  Supermatrix? matrix;

  /// Bases and gaps of each taxon in each locus of [matrix] (step 7 table).
  MatrixOverview? overview;
  final Map<String, String> matrixFiles = {};
  final List<TreeRun> treeRuns = [];

  void setHeader(HeaderTemplate h) {
    header = h;
    headerConfirmed = false;
    notifyListeners();
  }

  void confirmHeader() {
    final dups = headerConflicts();
    if (dups.isNotEmpty) {
      headerConfirmed = false;
      log('FASTA header ${header.fields.join(header.separator)}: ${dups.length} header(s) repeat inside the same marker '
          '(${dups.keys.take(5).join(', ')}${dups.length > 5 ? '…' : ''}) — add Accession or Voucher to tell them apart',
          level: LogLevel.error);
      notifyListeners();
      return;
    }
    headerConfirmed = true;
    log('FASTA header: ${header.fields.join(header.separator)}'
        '${dups.isEmpty ? '' : ' — ${dups.length} taxa with several records, joined into one sample each'}',
        level: LogLevel.success);
    notifyListeners();
  }

  /// FASTA headers that two records would share inside the same marker
  /// (the selected markers, or all while none is selected): MAFFT, trimAl
  /// and the tree programs could not tell those sequences apart. The same
  /// header in different markers is fine: it is how a sample's markers are
  /// joined in the matrix. header → accessions.
  ///
  /// With one sample per taxon, records of the same species sharing a
  /// header is the point (their markers are joined); records not named to
  /// species (Rhinella sp., cf., aff., a genus) are different specimens and
  /// still conflict. A FASTA sample named like a GenBank header conflicts
  /// too.
  Map<String, List<String>> headerConflicts([HeaderTemplate? template]) {
    final h = template ?? header;
    final cat = catalog;
    final visible = cat?.visible(featureTypes) ?? const <String, GeneSummary>{};
    if (visible.isEmpty) return onePerTaxonMode ? const {} : h.duplicates(records);
    final genes = selectedMarkers.isEmpty ? visible.keys : selectedMarkers.where(visible.containsKey);
    final names = <int, String>{};
    final out = <String, Set<String>>{};
    for (final g in genes) {
      final seen = <String, List<int>>{};
      for (final ri in visible[g]!.byRecord.keys) {
        seen.putIfAbsent(names[ri] ??= h.build(records[ri]), () => []).add(ri);
      }
      for (final e in seen.entries) {
        final joinable = onePerTaxonMode && e.value.every((ri) => isDeterminateSpecies(records[ri].organismName));
        if (e.value.length > 1 && !joinable) (out[e.key] ??= {}).addAll(e.value.map((ri) => records[ri].accession));
      }
      // FASTA samples of the same marker with a GenBank header's name.
      final cells = organized?.byLocus[g];
      if (cells != null) {
        for (final sample in cells.keys) {
          final clash = seen[sample];
          if (clash != null) (out[sample] ??= {}).addAll([...clash.map((ri) => records[ri].accession), 'FASTA $sample']);
        }
      }
    }
    return {for (final e in out.entries) e.key: e.value.toList()};
  }

  /// Records in a name conflict that are not named to species (Caranx sp.,
  /// cf., aff.): with one sample per species they cannot be joined.
  List<GbRecord> get undeterminedConflicts {
    final acc = {for (final v in headerConflicts().values) ...v};
    return [for (final r in records) if (acc.contains(r.accession) && !isDeterminateSpecies(r.organismName)) r];
  }

  /// Leaves out [undeterminedConflicts] (one click instead of removing them
  /// one by one in the records table).
  void removeUndeterminedConflicts() {
    final gone = undeterminedConflicts;
    if (gone.isEmpty) return;
    records.removeWhere(gone.toSet().contains);
    _rebuildCatalog();
    loci.clear();
    log('Left out ${gone.length} record(s) not identified to species: '
        '${gone.take(8).map((r) => '${r.organismName} (${r.accession})').join(', ')}${gone.length > 8 ? '…' : ''}', level: LogLevel.warn);
    notifyListeners();
  }

  /// Steps that use the sequence names (alignment onwards) stop while two
  /// sequences of one marker would share a name.
  bool _namesOk(String what) {
    if (!headerBlocked) return true;
    final dups = headerConflicts();
    log('$what blocked: ${dups.length} FASTA name(s) repeat inside the same marker (${dups.keys.take(5).join(', ')}${dups.length > 5 ? '…' : ''}); '
        'fix the names in 4a (add Accession or Voucher) or in step 2 before going on',
        level: LogLevel.error);
    notifyListeners();
    return false;
  }

  /// Repeated headers stop the alignment, unless one sample per taxon joins
  /// them on purpose.
  bool get headerBlocked => headerRequired && headerConflicts().isNotEmpty;

  /// GenBank data uses a FASTA header (5a); the default one
  /// (accession_family_genus_species) is used until the user changes it.
  bool get headerRequired => records.isNotEmpty;

  bool get _anyAligned => loci.values.any((l) => l.aligned != null);

  /// Panels 5a–5c show together with their defaults, so aligning takes
  /// one click (the run button ends 5c); 5d (trimming parameters) appears
  /// when the user chooses to trim after aligning, or once something was
  /// trimmed.
  bool alignPanelVisible(String id) => switch (id) {
        '4a' => headerRequired,
        // The rest waits while 5a has repeated names.
        _ when headerBlocked => false,
        '4b' || '4c' => true,
        _ => _anyAligned && (trimOpen || loci.values.any((l) => l.trimmed != null)),
      };

  /// Shows the trimming parameters (5d).
  void openTrim() {
    trimOpen = true;
    trimPrompt = false;
    trimSkipped = false;
    notifyListeners();
  }

  void setSeqMode(SeqMode m) {
    seqMode = m;
    notifyListeners();
  }

  void setCodonPolicy({bool? pseudo, bool? internalStops, bool? backtranslationWarnings}) {
    allowPseudo = pseudo ?? allowPseudo;
    allowInternalStops = internalStops ?? allowInternalStops;
    allowBacktranslationWarnings = backtranslationWarnings ?? allowBacktranslationWarnings;
    notifyListeners();
  }

  void setMafft(MafftOptions o) {
    mafft = o;
    notifyListeners();
  }

  void setTrimal(TrimalOptions o) {
    trimal = o;
    notifyListeners();
  }

  void setClipkitMode(String mode) {
    clipkitMode = mode;
    notifyListeners();
  }

  void setCodonAwareTrimming(bool v) {
    codonAwareTrimming = v;
    notifyListeners();
  }

  void setTrimTool(TrimTool t) {
    trimTool = t;
    notifyListeners();
  }

  void setLocusIncluded(String locus, bool included) {
    loci[locus]!.included = included;
    notifyListeners();
  }

  void setLocusCoding(String locus, {bool? coding, int? frame}) {
    final l = loci[locus]!;
    if (coding != null) l.coding = coding;
    if (frame != null) l.frame = frame;
    l.saturation = null;
    notifyListeners();
  }

  void dismissTrimPrompt({required bool skip}) {
    trimPrompt = false;
    trimSkipped = skip;
    if (skip) log('Trimming skipped; the MAFFT alignments will be concatenated');
    notifyListeners();
  }

  /// Builds the per-marker sequence sets from GenBank records (with the
  /// confirmed header) and FASTA loci. Called before aligning.
  Future<void> prepareLoci() async {
    _lociBuiltFrom = _lociSignature;
    final names = activeMarkers;
    final cat = catalog;
    // Copies: isolate closures must not capture the controller.
    final recs = List.of(records), types = Set.of(featureTypes), h = header, choices = Map.of(duplicateChoices);
    final wanted = names.toSet();
    var gb = cat == null
        ? const <String, List<SeqRecord>>{}
        : await _geneFastasInIsolate(recs, cat, types, wanted, h, choices);
    if (onePerTaxonMode) {
      final joined = onePerTaxon(gb);
      gb = joined.fastas;
      if (joined.dropped.isNotEmpty) {
        log('One sample per taxon: kept the longest sequence of each marker '
            '(${joined.dropped.entries.map((e) => '${e.key}: ${e.value} shorter left out').join('; ')})');
      }
    }
    final org = organized;
    loci.clear();
    matrix = null;
    overview = null;
    trimPrompt = false;
    trimSkipped = false;
    for (final name in names) {
      final seqs = <SeqRecord>[...?gb[name]];
      final cells = org?.byLocus[name];
      if (cells != null) {
        for (final s in cells.keys.toList()..sort()) {
          seqs.add(SeqRecord(s, cells[s]!.record.sequence.replaceAll('-', '')));
        }
      }
      if (seqs.isEmpty) continue;
      loci[name] = LocusState(name, _uniqueHeaders(seqs, name), sources: markers[name]?.sources ?? const {});
    }
    if (cat != null) {
      final partial = partialCounts(recs, cat, types, wanted, choices);
      for (final e in partial.entries) {
        final l = loci[e.key];
        if (l == null) continue;
        l.partial = math.min(e.value, l.raw.length);
        log('${l.name}: ${e.value} partial record(s) (incomplete feature, < or > in the location)'
            '${l.mostlyPartial ? ' — mostly fragments, the saturation test is optional for it' : ''}');
      }
    }
    final out = dir('01_loci');
    for (final l in loci.values) {
      await File(p.join(out, '${sanitizeTaxonName(l.name)}.fasta')).writeAsString(writeFasta(l.raw));
    }
    if (org != null) await File(p.join(out, 'splace_locus_mapping.tsv')).writeAsString(org.mappingTsv());
    log('${loci.length} marker(s) ready for alignment → $out');
  }

  /// Duplicated headers get `_2`, `_3`… (MAFFT and IQ-TREE need unique
  /// names). The alignment is blocked before it gets here ([headerBlocked]);
  /// a repeat that still reaches it is reported, never silent.
  List<SeqRecord> _uniqueHeaders(List<SeqRecord> seqs, [String locus = '']) {
    final count = <String, int>{};
    final out = [
      for (final s in seqs)
        () {
          final n = count.update(s.header, (v) => v + 1, ifAbsent: () => 1);
          return n == 1 ? s : SeqRecord('${s.header}_$n', s.sequence);
        }(),
    ];
    final repeated = [for (final e in count.entries) if (e.value > 1) e.key];
    if (repeated.isNotEmpty) {
      log('$locus: repeated FASTA name(s) renamed with _2, _3…: ${repeated.take(5).join(', ')}${repeated.length > 5 ? '…' : ''}', level: LogLevel.warn);
    }
    return out;
  }

  /// FASTA-only projects: use the organized loci (all selected).
  Future<void> useOrganizedLoci() async {
    await prepareLoci();
    notifyListeners();
  }

  /// FASTA-only projects whose selected loci are already aligned (every
  /// sequence of a locus has the same length): MAFFT can be skipped.
  bool get inputLooksAligned {
    final org = organized;
    if (org == null || records.isNotEmpty) return false;
    final names = activeMarkers.where(org.byLocus.containsKey).toList();
    return names.isNotEmpty &&
        names.every((m) {
          final seqs = [for (final c in org.byLocus[m]!.values) c.record];
          return seqs.length > 1 && isAligned(seqs);
        });
  }

  /// Uses the imported sequences, gaps kept, as alignments (skips MAFFT).
  Future<void> useInputAsAlignment() async {
    if (loci.isEmpty || _lociBuiltFrom != _lociSignature) await prepareLoci();
    await _task('Reading alignments', (_) async {
      final org = organized;
      for (final l in loci.values) {
        final cells = org?.byLocus[l.name];
        final input = cells == null ? l.raw : _uniqueHeaders([for (final s in cells.keys.toList()..sort()) SeqRecord(s, cells[s]!.record.sequence)]);
        final r = await _postAlignInIsolate(writeFasta(input), l.name, geneticCode);
        final coding = l.sources.contains('CDS') || (r.coding && _mayGuessCoding(l));
        l
          ..aligned = r.records
          ..stats = r.stats
          ..coding = coding
          ..frame = r.frame
          ..frameShifted = coding ? r.shifted : const []
          ..alignStatus = StepStatus.skipped;
      }
      notifyListeners();
    });
  }

  /// Loci read only from FASTA files have no annotation: short non-coding
  /// loci (UCEs, target capture) often have a frame without stop codons by
  /// chance, so only known protein-coding gene names are taken as coding;
  /// the user can mark others as coding in step 5.
  bool _mayGuessCoding(LocusState l) => l.sources.any((t) => t != 'FASTA') || knownCodingName(l.name);

  /// FASTA-only project (loci from per-sample or per-marker FASTA files).
  bool get lociMode => organized != null && records.isEmpty;

  /// Last alignment of each marker with the key of its input and settings:
  /// running again (e.g. after adding a marker) reuses what did not change.
  final Map<String, _AlignSnapshot> _alignCache = {};

  /// Last trimming of each marker (alignment key + trimming settings).
  final Map<String, ({String key, List<SeqRecord> trimmed, String? method, bool frameLost, List<int>? columns})> _trimCache = {};

  String _alignKey(LocusState l) {
    final args = mafft.args('').where((a) => a.isNotEmpty).toList();
    final t = args.indexOf('--thread');
    if (t >= 0) args.removeRange(t, t + 2);
    return [seqMode.name, args.join(' '), geneticCode, allowPseudo, allowInternalStops, allowBacktranslationWarnings, fingerprint(l.raw)].join('|');
  }

  String? _trimKey(LocusState l) {
    final a = _alignCache[l.name];
    if (a == null) return null;
    return [
      a.key,
      trimmedByClipkit(l)
          ? 'clipkit $clipkitMode ${codonTrimmed(l) ? '--codon ${l.frame}' : ''}'
          : 'trimal ${trimal.args('', '').join(' ')}${trimmedByBacktrans(l) ? ' -backtrans ${l.frame} $geneticCode' : ''}',
    ].join('|');
  }

  void _item(String name, {double? progress, ItemState? state, String? note}) {
    final it = modal?.items[name];
    if (it == null) return;
    if (state != null) it.state = state;
    // One pass per marker: while it runs its bar only moves forward
    // through all of its stages.
    if (progress != null) it.progress = it.state == ItemState.running ? math.max(it.progress, progress) : progress;
    if (note != null) it.note = note;
    _notifySoon();
  }

  Timer? _notifyTimer;

  /// Progress of many markers at once (parallel jobs, each reporting its
  /// stages) would rebuild the whole page dozens of times a second: with
  /// more than 24 markers the updates are gathered and shown about ten
  /// times a second.
  void _notifySoon() {
    if ((modal?.items.length ?? 0) <= 24) {
      notifyListeners();
      return;
    }
    if (_notifyTimer?.isActive ?? false) return;
    _notifyTimer = Timer(const Duration(milliseconds: 100), notifyListeners);
  }

  /// Runs MAFFT for every included marker in the selected sequence mode.
  /// Markers whose sequences and settings did not change since their last
  /// alignment are taken from the cache.
  Future<void> alignAll() async {
    if (headerRequired && !headerConfirmed) confirmHeader();
    if (!_namesOk('Alignment')) return;
    log('Sequence type: ${seqMode.name}${seqMode == SeqMode.nucleotide ? '' : ', genetic code $geneticCode'}');
    if (loci.isEmpty || _lociBuiltFrom != _lociSignature) await prepareLoci();
    final included = loci.values.where((l) => l.included).toList();
    final modal = TaskModal(ModalKind.alignment, seqMode.name)..total = included.length;
    for (final l in included) {
      modal.items[l.name] = ModalItem();
    }
    await _task(modal: modal, 'MAFFT alignment (${seqMode.name})', (cancel) async {
      final tool = await _require(ExternalTool.mafft);
      final keys = {for (final l in included) l.name: _alignKey(l)};
      // Unchanged markers come from the cache.
      final cached = [
        for (final l in included)
          if (_alignCache[l.name]?.key == keys[l.name]) l,
      ];
      for (final l in cached) {
        _alignCache[l.name]!.restore(l);
        _restoreTrim(l);
        await File(p.join(dir('02_aligned'), '${sanitizeTaxonName(l.name)}_aligned.fasta')).writeAsString(writeFasta(l.aligned!));
        _item(l.name, progress: 1, state: ItemState.cached, note: _alignedNote(l));
      }
      final items = included.where((l) => !cached.contains(l)).toList();
      final cores = math.max(1, Platform.numberOfProcessors - 2);
      // Many markers: one thread each and as many at once as there are
      // cores (MAFFT gains little from threads on small alignments, and
      // each run has a fixed start-up cost).
      final perJob = items.length >= cores ? 1 : math.max(1, mafft.threads);
      _mafftThreads = perJob;
      final jobs = math.max(1, math.min(math.max(1, items.length), cores ~/ perJob));
      if (cached.isNotEmpty) {
        log('${cached.length} marker(s) unchanged since the last alignment, reused'
            '${cached.length <= 20 ? ': ${cached.map((l) => l.name).join(', ')}' : ''}', source: 'MAFFT');
      }
      log('Aligning ${items.length} marker(s): $jobs parallel job(s) × $perJob thread(s)', source: 'MAFFT');
      final codonLog = <String>[];
      var finished = cached.length;
      _progress(finished, included.length);
      await _pool(items, jobs, cancel, (l) async {
        l
          ..alignStatus = StepStatus.running
          ..trimmed = null
          ..trimmedColumns = null
          ..trimStatus = StepStatus.idle
          ..trimMethod = null
          ..saturation = null
          ..message = null
          ..pipeline.clear();
        _item(l.name, state: ItemState.running, progress: 0.01);
        void stage(double f, String note) => _item(l.name, progress: f, note: note);
        try {
          switch (seqMode) {
            case SeqMode.nucleotide:
              await _alignNucleotide(tool, l, cancel, stage);
            case SeqMode.aminoAcid:
              await _alignProtein(tool, l, cancel, backTranslate: false, logLines: codonLog, stage: stage);
            case SeqMode.codon:
              await _alignProtein(tool, l, cancel, backTranslate: true, logLines: codonLog, stage: stage);
          }
          if (l.alignStatus == StepStatus.done) {
            _alignCache[l.name] = _AlignSnapshot(l, keys[l.name]!);
            _item(l.name, progress: 1, state: ItemState.done, note: _alignedNote(l));
          } else {
            _item(l.name, progress: 1, state: ItemState.skipped, note: l.message);
          }
        } catch (e) {
          l
            ..alignStatus = StepStatus.error
            ..message = '$e';
          log('${l.name}: $e', source: 'MAFFT', error: true);
          _item(l.name, progress: 1, state: cancel.isCancelled ? ItemState.skipped : ItemState.failed, note: '$e');
        }
        _progress(++finished, included.length);
      }, label: (l) => l.name, reportProgress: false);
      if (codonLog.isNotEmpty) {
        await File(p.join(dir('02_aligned'), 'codon_validation.log')).writeAsString(codonLog.join('\n'));
      }
      final ok = included.where((l) => l.alignStatus == StepStatus.done).length;
      modal
        ..result = 'Aligned $ok of ${included.length} marker(s)${cached.isEmpty ? '' : ' (${cached.length} from the cache)'}'
        ..failed = ok == 0;
      trimPrompt = ok > 0 && included.any((l) => l.trimmed == null && l.alignStatus == StepStatus.done);
      log('Aligned $ok/${included.length} marker(s) → ${dir('02_aligned')}', source: 'MAFFT', level: ok == included.length ? LogLevel.success : LogLevel.warn);
    });
  }

  String _alignedNote(LocusState l) => '${l.aligned!.length} seq · ${l.aligned!.first.length} ${l.isProtein ? 'aa' : 'bp'}';

  /// Puts back the cached trimming of [l] when it matches its alignment.
  void _restoreTrim(LocusState l) {
    final t = _trimCache[l.name];
    if (t == null || t.key != _trimKey(l)) return;
    l
      ..trimmed = t.trimmed
      ..trimMethod = t.method
      ..trimmedFrameLost = t.frameLost
      ..trimmedColumns = t.columns
      ..stats = alignmentStats(t.trimmed)
      ..trimStatus = StepStatus.done;
  }

  /// What the loci were built from: the selected markers, the FASTA
  /// header of each record, the chosen copies of duplicated genes and the
  /// FASTA loci. Aligning rebuilds the loci when it changed.
  String get _lociSignature => [
        activeMarkers.join(','),
        [for (final r in records) header.build(r)].join(','),
        featureTypes.join(','),
        (duplicateChoices.entries.map((e) => '${e.key}=${e.value}').toList()..sort()).join(','),
        identityHashCode(catalog),
        identityHashCode(organized),
      ].join('|');

  String? _lociBuiltFrom;

  /// Runs MAFFT on [input] ([sequences] records); [onProgress] gets the
  /// fraction read from MAFFT's counters. Only notable lines are logged
  /// (the strategy, warnings and errors), or the last lines on failure.
  /// Threads of each MAFFT run of the current alignment (see [alignAll]).
  int? _mafftThreads;

  Future<String> _runMafft(ResolvedTool tool, LocusState l, String input, CancelToken cancel,
      {bool amino = false, int sequences = 2, void Function(double)? onProgress}) async {
    final progress = MafftProgress(sequences);
    final tail = <String>[];
    var strategy = false;
    final out = await runTool(tool, mafft.copyWith(amino: amino, threads: _mafftThreads ?? mafft.threads).args(input),
        workingDirectory: dir('02_aligned'), captureStdout: true, cancel: cancel, onLog: (m) {
      if (MafftProgress.isProgressLine(m)) {
        if (progress.add(m)) onProgress?.call(progress.value);
        return;
      }
      tail.add(m);
      if (tail.length > 25) tail.removeAt(0);
      final notable = m.startsWith(r'$ ') || strategy || RegExp(r'strategy|error|warning|cannot|not found', caseSensitive: false).hasMatch(m);
      strategy = m.trim() == 'Strategy:';
      if (notable && m.trim().isNotEmpty) log(m.trim(), source: 'MAFFT ${l.name}');
    });
    if (!out.ok || out.stdout.trim().isEmpty) {
      for (final m in tail) {
        log(m, source: 'MAFFT ${l.name}', level: LogLevel.warn);
      }
      throw StateError(out.cancelled ? 'cancelled' : 'MAFFT exit code ${out.exitCode}');
    }
    return out.stdout;
  }

  Future<void> _alignNucleotide(ResolvedTool tool, LocusState l, CancelToken cancel, void Function(double, String) stage) async {
    final name = sanitizeTaxonName(l.name);
    final input = p.join(dir('01_loci'), '$name.fasta');
    await File(input).writeAsString(writeFasta(l.raw));
    stage(0.03, 'MAFFT');
    final text = await _runMafft(tool, l, input, cancel, sequences: l.raw.length, onProgress: (f) => stage(0.03 + 0.9 * f, 'MAFFT'));
    stage(0.95, 'checking frame and orientation');
    await File(p.join(dir('02_aligned'), '${name}_aligned.fasta')).writeAsString(text);
    final r = await _postAlignInIsolate(text, l.name, geneticCode);
    if (r.reversed.isNotEmpty) {
      log('${l.name}: reverse-complemented ${r.reversed.length} sequence(s): ${r.reversed.take(5).join(', ')}', source: 'MAFFT');
    }
    final types = l.sources;
    final rnaOnly = types.isNotEmpty && types.every((t) => t == 'rRNA' || t == 'tRNA');
    // Fragments (partial CDS, e.g. barcodes) start at different codon
    // positions: their codon positions are only used when the alignment
    // reads in one frame without stops. Otherwise the marker is trimmed
    // column by column and stays one partition (no codon positions).
    final partialCds = types.contains('CDS') && l.partial > 0;
    final coding = !rnaOnly && (partialCds ? r.coding : (types.contains('CDS') || (r.coding && _mayGuessCoding(l))));
    if (partialCds && !r.coding) {
      log('${l.name}: the partial sequences do not read in one frame without stop codons in the nucleotide alignment, so its codon '
          'positions are unknown: it is trimmed column by column and kept as one partition. The Codon-aware mode (4b) aligns the '
          'translated fragments and keeps the codon positions.', source: 'MAFFT', level: LogLevel.warn);
    }
    l
      ..aligned = r.records
      ..isProtein = false
      ..stats = r.stats
      ..coding = coding
      ..frame = r.frame
      ..frameShifted = coding ? r.shifted : const []
      ..alignStatus = StepStatus.done
      ..pipeline.addAll(['MAFFT ${mafft.strategy.split(' ').first}', 'DNA']);
    if (l.frameShifted.isNotEmpty) {
      log('${l.name}: ${l.frameShifted.length} sequence(s) have gaps that are not multiples of 3 (frame shifts): '
          '${l.frameShifted.take(5).join(', ')}${l.frameShifted.length > 5 ? '…' : ''}. Their codon positions are shifted in the '
          'nucleotide alignment; the Codon-aware mode avoids this.', source: 'MAFFT', level: LogLevel.warn);
    }
  }

  /// Amino-acid and codon-aware modes: validate and translate the CDS,
  /// align the proteins with `mafft --amino`, and (codon mode)
  /// back-translate to a codon-preserving DNA alignment.
  Future<void> _alignProtein(ResolvedTool tool, LocusState l, CancelToken cancel,
      {required bool backTranslate, required List<String> logLines, required void Function(double, String) stage}) async {
    final name = sanitizeTaxonName(l.name);
    void clog(String m) {
      logLines.add(m);
      log(m, source: 'Codon ${l.name}');
    }

    stage(0.02, 'CDS validation');
    l.frameShifted = const [];
    final cdsList = await _cdsFor(l, clog);
    if (cdsList == null) {
      l
        ..alignStatus = StepStatus.skipped
        ..message = 'not protein-coding (skipped in ${backTranslate ? 'codon-aware' : 'amino-acid'} mode)';
      log('${l.name}: ${l.message}', source: 'MAFFT', level: LogLevel.warn);
      return;
    }
    final policy = allowInternalStops ? InternalStopPolicy.allowWithWarning : InternalStopPolicy.exclude;
    final proteins = <SeqRecord>[];
    final cds = <String, String>{};
    final codes = <String, int>{};
    var skipped = 0, normalized = 0, stops = 0;
    for (final c in cdsList) {
      final v = validateCds(l.name, c.header, c.dna, c.code, stopPolicy: policy, log: clog);
      if (v == null) {
        skipped++;
        continue;
      }
      if (v.normalized) normalized++;
      if (v.terminalStopRemoved) stops++;
      proteins.add(SeqRecord(c.header, v.protein));
      cds[c.header] = v.validDna;
      codes[c.header] = v.code;
    }
    clog('[CDS validation] [${l.name}] ${cdsList.length} raw | ${proteins.length} valid | $skipped skipped | '
        '$normalized incomplete stop(s) completed | $stops terminal stop(s) removed');
    if (proteins.length < 2) throw StateError('fewer than 2 valid CDS after validation');
    final aligned = dir('02_aligned');
    await File(p.join(aligned, '${name}_cds_validated.fasta')).writeAsString(writeFasta([for (final e in cds.entries) SeqRecord(e.key, e.value)]));
    final protIn = p.join(aligned, '${name}_protein.fasta');
    await File(protIn).writeAsString(writeFasta(proteins));
    stage(0.1, 'MAFFT --amino');
    final span = backTranslate ? 0.8 : 0.88;
    final protAln = parseFasta(await _runMafft(tool, l, protIn, cancel,
        amino: true, sequences: proteins.length, onProgress: (f) => stage(0.1 + span * f, 'MAFFT --amino')));
    await File(p.join(aligned, '${name}_protein_aligned.fasta')).writeAsString(writeFasta(protAln));
    l
      ..proteinAligned = protAln
      ..cds = cds
      ..codes = codes
      ..pipeline.addAll(['CDS ${proteins.length}/${cdsList.length}', 'translate (code ${codes.values.toSet().join('/')})', 'MAFFT --amino']);
    if (!backTranslate) {
      l
        ..aligned = protAln
        ..isProtein = true
        ..coding = false
        ..frame = 0
        ..stats = alignmentStats(protAln)
        ..alignStatus = StepStatus.done;
      return;
    }
    stage(0.92, 'back-translation');
    final bt = backTranslateAlignment(protAln, cds, codes, geneticCode);
    final integrity = codonIntegrity(bt.dna);
    final problems = [...integrity.issues(frameDisrupted: bt.frameDisrupted), if (bt.mismatches > 0) '${bt.mismatches} amino-acid/codon mismatch(es)'];
    for (final i in bt.issues.take(10)) {
      clog('[Back-translation] [${l.name}] $i');
    }
    if (problems.isNotEmpty) {
      final msg = problems.join('; ');
      if (!allowBacktranslationWarnings || !integrity.ok) throw StateError('codon integrity check failed: $msg');
      log('${l.name}: back-translation warnings kept by policy: $msg', source: 'Codon', level: LogLevel.warn);
    }
    await File(p.join(aligned, '${name}_aligned.fasta')).writeAsString(writeFasta(bt.dna));
    l
      ..aligned = bt.dna
      ..isProtein = false
      ..coding = true
      ..frame = 0
      ..stats = alignmentStats(bt.dna)
      ..alignStatus = StepStatus.done
      ..pipeline.add('back-translate');
  }

  /// CDS for a marker: GenBank CDS features (codon_start, transl_table)
  /// and FASTA sequences of coding loci, shifted to their best frame.
  Future<List<CodonCandidate>?> _cdsFor(LocusState l, void Function(String) clog) async {
    final out = <CodonCandidate>[];
    final cat = catalog;
    if (cat != null && l.sources.contains('CDS')) {
      final names = {l.name};
      final recs = List.of(records), h = header, choices = Map.of(duplicateChoices);
      final code = geneticCode, pseudo = allowPseudo;
      final c = await _codonCandidatesInIsolate(recs, cat, names, h, choices, code, pseudo);
      out.addAll(onePerTaxonMode ? oneCdsPerTaxon(c[l.name] ?? const []) : (c[l.name] ?? const []));
    }
    final org = organized;
    if (org != null && l.sources.contains('FASTA')) {
      final cells = org.byLocus[l.name] ?? const {};
      final seqs = [for (final e in cells.entries) SeqRecord(e.key, e.value.record.sequence.replaceAll('-', ''))];
      final guess = guessCoding(l.name, seqs, geneticCode);
      if (guess.coding) {
        for (final s in seqs) {
          out.add(CodonCandidate(s.header, s.sequence.substring(math.min(guess.frame, s.sequence.length)), geneticCode, false));
        }
        if (guess.frame > 0) clog('[CDS extraction] [${l.name}] FASTA sequences read from frame +${guess.frame + 1}');
      }
    }
    if (out.isEmpty) return null;
    final seen = <String, int>{};
    return [
      for (final c in out)
        () {
          final n = seen.update(c.header, (v) => v + 1, ifAbsent: () => 1);
          return n == 1 ? c : CodonCandidate('${c.header}_$n', c.dna, c.code, c.codeFromFeature);
        }(),
    ];
  }

  // ---------------------------------------------------------------------
  // 5d. Trimming (trimAl / ClipKIT)
  // ---------------------------------------------------------------------

  /// Trims every aligned marker: trimAl for nucleotide and amino-acid
  /// alignments, ClipKIT (on the protein, then back-translated) for
  /// codon-aware ones. Unchanged markers come from the cache.
  Future<void> trimAll() async {
    if (!_namesOk('Trimming')) return;
    final items = loci.values.where((l) => l.included && l.aligned != null).toList();
    final modal = TaskModal(ModalKind.trimming, seqMode.name)..total = items.length;
    for (final l in items) {
      modal.items[l.name] = ModalItem();
    }
    await _task(modal: modal, 'Trimming', (cancel) async {
      trimPrompt = false;
      trimSkipped = false;
      final trimalTool = tools[ExternalTool.trimal] ?? await locator.resolve(ExternalTool.trimal);
      var finished = 0, cached = 0;
      _progress(0, items.length);
      final todo = <LocusState>[];
      for (final l in items) {
        final key = _trimKey(l);
        final t = _trimCache[l.name];
        if (key != null && t != null && t.key == key && l.trimmed != null) {
          cached++;
          _item(l.name, progress: 1, state: ItemState.cached, note: '${t.method} · ${l.trimmed!.first.length} ${l.isProtein ? 'aa' : 'bp'}');
          _progress(++finished, items.length);
        } else {
          todo.add(l);
        }
      }
      await _pool(todo, math.max(1, Platform.numberOfProcessors - 2), cancel, (l) async {
        l.trimStatus = StepStatus.running;
        final clip = trimmedByClipkit(l), codon = codonTrimmed(l);
        _item(l.name,
            state: ItemState.running,
            progress: 0.05,
            note: clip ? 'ClipKIT $clipkitMode${codon ? ' --codon' : ''}' : 'trimAl ${trimal.label}${trimmedByBacktrans(l) ? ' -backtrans' : ''}');
        try {
          if (clip) {
            await _trimClipkit(l, cancel, (f, [note]) => _item(l.name, progress: f, note: note), codon: codon);
          } else if (trimmedByBacktrans(l)) {
            if (trimalTool == null) throw StateError('trimAl not found. ${ExternalTool.trimal.installHint}');
            await _trimBacktrans(trimalTool, l, cancel, (f, [note]) => _item(l.name, progress: f, note: note));
          } else {
            if (trimalTool == null) throw StateError('trimAl not found. ${ExternalTool.trimal.installHint}');
            await _trimPlain(trimalTool, l, cancel, (f, [note]) => _item(l.name, progress: f, note: note));
          }
          l
            ..stats = alignmentStats(l.trimmed!)
            ..saturation = null
            ..trimStatus = StepStatus.done;
          final key = _trimKey(l);
          if (key != null) _trimCache[l.name] = (key: key, trimmed: l.trimmed!, method: l.trimMethod, frameLost: l.trimmedFrameLost, columns: l.trimmedColumns);
          final before = l.aligned!.first.length, after = l.trimmed!.first.length;
          _item(l.name, progress: 1, state: ItemState.done, note: '$before → $after ${l.isProtein ? 'aa' : 'bp'}');
        } catch (e) {
          l
            ..trimStatus = StepStatus.error
            ..message = '$e';
          log('${l.name}: $e', source: 'trimming', error: true);
          _item(l.name, progress: 1, state: cancel.isCancelled ? ItemState.skipped : ItemState.failed, note: '$e');
        }
        _progress(++finished, items.length);
      }, label: (l) => l.name, reportProgress: false);
      final ok = items.where((l) => l.trimStatus == StepStatus.done).length;
      final removed = items
          .where((l) => l.trimStatus == StepStatus.done && l.trimmed != null && l.aligned != null)
          .fold<int>(0, (a, l) => a + l.aligned!.first.length - l.trimmed!.first.length);
      if (ok > 0 && removed == 0) {
        log('No column was removed: the alignments have no gappy or poorly aligned regions for this method '
            '(the kpic modes of ClipKIT also remove sites without phylogenetic information).', source: 'trimming', level: LogLevel.warn);
      }
      modal
        ..result = 'Trimmed $ok of ${items.length} marker(s)${cached == 0 ? '' : ' ($cached from the cache)'}'
            '${ok > 0 && removed == 0 ? ' — no column removed (nothing to trim for this method)' : ' — $removed column(s) removed'}'
        ..failed = ok == 0 && items.isNotEmpty;
    });
  }

  /// Coding DNA trimmed in whole codons by ClipKIT `--codon`: the codon-aware
  /// alignments (5b), and the coding markers of a nucleotide alignment when
  /// "Keep the reading frame" is on.
  bool codonTrimmed(LocusState l) =>
      !l.isProtein && codonAwareTrimming && ((l.proteinAligned != null && l.cds.isNotEmpty) || l.coding);

  /// The marker is trimmed by ClipKIT (with `--codon` when [codonTrimmed]):
  /// codon-aware alignments, and every DNA marker when ClipKIT is the
  /// trimming tool of a nucleotide alignment ([codonAwareTrimming]).
  /// trimAl trims the others (and amino-acid alignments).
  bool trimmedByClipkit(LocusState l) => !l.isProtein && trimTool == TrimTool.clipkitCodon;

  /// Coding DNA trimmed by trimAl in whole codons (`-backtrans`).
  bool trimmedByBacktrans(LocusState l) => trimTool == TrimTool.trimalCodon && codonTrimmed(l);

  /// Programs that print no progress: the bar eases from [from] towards
  /// [to] at a pace set by the size of the alignment, until [run] ends.
  Future<T> _eased<T>(LocusState l, double from, double to, void Function(double f, [String? note]) stage, Future<T> Function() run) async {
    final cells = l.aligned!.length * (l.aligned!.firstOrNull?.length ?? 0);
    final tau = (cells / 4e5).clamp(0.4, 30.0);
    final started = DateTime.now();
    final ticker = Timer.periodic(const Duration(milliseconds: 150), (_) {
      final s = DateTime.now().difference(started).inMilliseconds / 1000;
      stage(from + (to - from) * (1 - math.exp(-s / tau)));
    });
    try {
      return await run();
    } finally {
      ticker.cancel();
    }
  }

  /// trimAl on one alignment, column by column. In a coding DNA alignment
  /// this can remove single bases of a codon, so the frame is marked lost.
  Future<void> _trimPlain(ResolvedTool tool, LocusState l, CancelToken cancel, void Function(double f, [String? note]) stage) async {
    final name = sanitizeTaxonName(l.name);
    final input = p.join(dir('02_aligned'), '${name}_aligned.fasta');
    final output = p.join(dir('03_trimmed'), '${name}_trimmed.fasta');
    await File(input).writeAsString(writeFasta(l.aligned!));
    l.trimmedFrameLost = false;
    stage(0.1, 'trimAl ${trimal.label}');
    // -colnumbering: the kept columns, to map the trimmed alignment back to
    // the MAFFT one.
    final out = await _eased(
        l,
        0.1,
        0.85,
        stage,
        () => runTool(tool, [...trimal.args(input, output), '-colnumbering'],
            workingDirectory: dir('03_trimmed'), captureStdout: true, cancel: cancel, onLog: (m) => log(m, source: 'trimAl ${l.name}')));
    if (!out.ok) throw StateError('trimAl exit code ${out.exitCode}');
    stage(0.9, 'reading the trimmed alignment');
    final kept = parseTrimalColumnMap(out.stdout);
    l
      ..trimmed = parseFasta(await File(output).readAsString())
      ..trimmedColumns = kept != null && l.trimmed != null && kept.length == l.trimmed!.firstOrNull?.length ? kept : null
      ..trimMethod = 'trimAl ${trimal.label}';
    // trimAl removes single columns: in a coding DNA alignment this shifts
    // the reading frame, so codon positions are unknown from here.
    if (l.coding && !l.isProtein) {
      l.trimmedFrameLost = true;
      log('${l.name}: trimAl removed single columns of a protein-coding DNA alignment; the reading frame is lost '
          '(${l.trimmed!.first.length} bp${l.trimmed!.first.length % 3 == 0 ? '' : ', not a multiple of 3'}). It will be one '
          'partition without codon positions and is left out of the saturation test. Turn on "Keep the reading frame" to avoid this.',
          source: 'trimAl', level: LogLevel.warn);
    }
  }

  /// Whole-codon trimming by trimAl (`-backtrans`, Capella-Gutiérrez et al.
  /// 2009): the DNA alignment is read codon by codon from its frame and
  /// translated with the project's genetic code; trimAl trims that protein
  /// alignment with the chosen method and writes back the codons of the
  /// kept amino acids, so whole codons go or stay. Codons broken by a gap
  /// (frameshifts) go in as X / NNN (logged). trimAl checks stops with the
  /// standard code, hence `-ignorestopcodon`. The kept columns come from a
  /// second run of the same method with `-colnumbering`.
  Future<void> _trimBacktrans(ResolvedTool tool, LocusState l, CancelToken cancel, void Function(double f, [String? note]) stage) async {
    final name = sanitizeTaxonName(l.name);
    final work = dir('03_trimmed');
    final aln = l.aligned!;
    final frame = l.proteinAligned != null ? 0 : l.frame;
    final length = aln.first.length;
    final codons = (length - frame) ~/ 3;
    if (codons <= 0) throw StateError('no complete codon');
    final table = codonTable(geneticCode);
    final ids = [for (var i = 0; i < aln.length; i++) 'S${i + 1}'];
    final prot = <SeqRecord>[], cds = <SeqRecord>[];
    var masked = 0;
    for (final (i, r) in aln.indexed) {
      final pb = StringBuffer(), cb = StringBuffer();
      for (var k = 0; k < codons; k++) {
        final at = frame + 3 * k;
        final codon = (at + 3 <= r.sequence.length ? r.sequence.substring(at, at + 3) : '---').toUpperCase().replaceAll('U', 'T').replaceAll('.', '-');
        if (codon == '---') {
          pb.write('-');
        } else if (codon.contains('-')) {
          pb.write('X');
          cb.write('NNN');
          masked++;
        } else {
          final aa = table[codon] ?? 'X';
          pb.write(aa);
          cb.write(codon);
        }
      }
      prot.add(SeqRecord(ids[i], pb.toString()));
      cds.add(SeqRecord(ids[i], cb.toString()));
    }
    final protIn = p.join(work, '${name}_backtrans_protein.fasta');
    final cdsIn = p.join(work, '${name}_backtrans_cds.fasta');
    final protOut = p.join(work, '${name}_backtrans_protein_trimmed.fasta');
    final output = p.join(work, '${name}_trimmed.fasta');
    await File(protIn).writeAsString(writeFasta(prot));
    await File(cdsIn).writeAsString(writeFasta(cds));
    if (masked > 0 || frame > 0 || frame + 3 * codons < length) {
      log('${l.name}: ${masked > 0 ? '$masked codon(s) broken by gaps (frameshifts) given to trimAl as X/NNN' : ''}'
          '${masked > 0 && (frame > 0 || frame + 3 * codons < length) ? '; ' : ''}'
          '${frame > 0 || frame + 3 * codons < length ? '${length - 3 * codons} column(s) of incomplete codons at the ends left out' : ''}',
          source: 'trimAl', level: masked > 0 ? LogLevel.warn : LogLevel.info);
    }
    stage(0.1, 'trimAl ${trimal.label} -backtrans');
    // Kept amino-acid columns (same method, no back-translation).
    final map = await _eased(
        l,
        0.1,
        0.45,
        stage,
        () => runTool(tool, [...trimal.args(protIn, protOut), '-colnumbering'],
            workingDirectory: work, captureStdout: true, cancel: cancel, onLog: (m) => log(m, source: 'trimAl ${l.name}')));
    if (!map.ok) throw StateError(map.cancelled ? 'cancelled' : 'trimAl exit code ${map.exitCode}');
    final kept = parseTrimalColumnMap(map.stdout);
    // The codon alignment.
    final out = await _eased(
        l,
        0.5,
        0.88,
        stage,
        () => runTool(tool, [...trimal.args(protIn, output, backtranslateFrom: cdsIn), '-ignorestopcodon'],
            workingDirectory: work, cancel: cancel, onLog: (m) => log(m, source: 'trimAl ${l.name}')));
    if (!out.ok || !File(output).existsSync()) throw StateError(out.cancelled ? 'cancelled' : 'trimAl -backtrans exit code ${out.exitCode}');
    stage(0.92, 'reading the trimmed alignment');
    final byId = {for (final r in parseFasta(await File(output).readAsString())) r.header.split(RegExp(r'\s')).first: r.sequence};
    final trimmed = [for (final (i, r) in aln.indexed) SeqRecord(r.header, byId[ids[i]] ?? '')];
    if (trimmed.any((r) => r.sequence.isEmpty) || trimmed.first.length % 3 != 0) {
      throw StateError('trimAl -backtrans returned ${byId.length} of ${aln.length} sequence(s), ${trimmed.first.length} bp');
    }
    await File(output).writeAsString(writeFasta(trimmed));
    final columns = kept == null ? null : [for (final k in kept) ...[frame + 3 * k, frame + 3 * k + 1, frame + 3 * k + 2]];
    l
      ..trimmed = trimmed
      ..trimmedFrameLost = false
      ..trimmedColumns = columns != null && columns.length == trimmed.first.length ? columns : null
      ..trimMethod = 'trimAl ${trimal.label} -backtrans';
    log('${l.name}: trimAl -backtrans kept ${trimmed.first.length ~/ 3}/$codons codons', source: 'trimAl');
  }

  /// Whole-codon trimming by ClipKIT `--codon` (Steenwyk et al. 2020): when
  /// one position of a codon is to be trimmed, ClipKIT removes the whole
  /// codon. The alignment goes in from its first complete codon; the kept
  /// columns come from ClipKIT's own log (`-l`).
  Future<void> _trimClipkit(LocusState l, CancelToken cancel, void Function(double f, [String? note]) stage, {required bool codon}) async {
    final name = sanitizeTaxonName(l.name);
    final trimmedDir = dir('03_trimmed');
    final clip = tools[ExternalTool.clipkit] ?? await locator.resolve(ExternalTool.clipkit);
    if (clip == null) throw StateError('ClipKIT is not installed (step 4d: Install ClipKIT).');
    final aln = l.aligned!;
    // Without codons (non-coding markers) the whole alignment goes in.
    final frame = !codon || l.proteinAligned != null ? 0 : l.frame;
    final length = aln.first.length;
    final end = codon ? frame + 3 * ((length - frame) ~/ 3) : length;
    if (end <= frame) throw StateError('no complete codon');
    final input = p.join(trimmedDir, '${name}_codon_input.fasta');
    final output = p.join(trimmedDir, '${name}_codon_trimmed.fasta');
    await File(input).writeAsString(writeFasta([for (final r in aln) SeqRecord(r.header, r.sequence.substring(frame, end))]));
    if (codon && (frame > 0 || end < length)) {
      log('${l.name}: ${frame + length - end} column(s) of incomplete codons at the ends left out before ClipKIT --codon', source: 'ClipKIT');
    }
    stage(0.1, 'ClipKIT $clipkitMode${codon ? ' --codon' : ''}');
    final out = await _eased(
        l,
        0.1,
        0.85,
        stage,
        () => runTool(clip, [input, '-o', output, '-m', clipkitMode, if (codon) '--codon', '-l'],
            workingDirectory: trimmedDir, cancel: cancel, onLog: (m) => log(m, source: 'ClipKIT ${l.name}')));
    if (!out.ok || !File(output).existsSync()) throw StateError(out.cancelled ? 'cancelled' : 'ClipKIT exit code ${out.exitCode}');
    stage(0.9, 'reading the trimmed alignment');
    final trimmed = parseFasta(await File(output).readAsString());
    final logFile = File('$output.log');
    final kept = logFile.existsSync() ? parseClipkitLog(await logFile.readAsString()) : null;
    final columns = kept == null ? null : [for (final c in kept) c + frame];
    if (trimmed.isEmpty || trimmed.length != aln.length || (codon && trimmed.first.length % 3 != 0)) {
      throw StateError('ClipKIT returned ${trimmed.length} sequence(s) of ${trimmed.firstOrNull?.length ?? 0} bp');
    }
    l
      ..trimmed = trimmed
      ..trimmedFrameLost = false
      ..trimmedColumns = columns != null && columns.length == trimmed.first.length ? columns : null
      ..trimMethod = 'ClipKIT $clipkitMode${codon ? ' --codon' : ''}';
    log(codon
        ? '${l.name}: ClipKIT --codon kept ${trimmed.first.length ~/ 3}/${(end - frame) ~/ 3} codons'
        : '${l.name}: ClipKIT kept ${trimmed.first.length}/$length columns', source: 'ClipKIT');
  }

  // ---------------------------------------------------------------------
  // 6. Saturation
  // ---------------------------------------------------------------------

  /// DNA loci without codon positions (non-coding, or trimmed without
  /// keeping the frame): their saturation is tested on all sites.
  List<LocusState> get _wholeSiteLoci =>
      loci.values.where((l) => l.included && !l.isProtein && !codonPositionsUsable(l) && _aln(l) != null).toList();

  /// Step 7 waits for the saturation test: GenBank data with a complete
  /// protein-coding marker. Loci from FASTA files and markers with partial
  /// sequences (individual-marker searches bring fragments) can skip it.
  bool get saturationRequired =>
      !lociMode && loci.values.any((l) => l.included && !l.mostlyPartial && codonPositionsUsable(l) && l.alignment != null);

  /// Coding markers made mostly of partial sequences (their test is optional).
  List<LocusState> get partialLoci =>
      [for (final l in loci.values) if (l.included && l.mostlyPartial && codonPositionsUsable(l) && l.alignment != null) l];

  /// Something can be tested in step 6.
  bool get canTestSaturation => loci.values.any((l) => l.included && !l.isProtein && _aln(l) != null);

  /// Parameters of the saturation test (step 6a): the verdict thresholds,
  /// the Iss replicates and the points kept for the plots.
  SaturationThresholds saturationThresholds = phykitThresholds;
  int issReplicates = 400, satPlotPoints = 3000;

  /// IQ-TREE model of the gene trees and of the branch lengths of each
  /// codon subset, and PhyKIT's `--exclude_gaps`.
  String satModel = 'GTR+G';
  bool satExcludeGaps = false;

  /// The results on screen were made with other parameters.
  bool get saturationStale => loci.values.any((l) => l.saturation != null) && _satParamsUsed != _satParams;
  String? _satParamsUsed;
  String get _satParams {
    final th = saturationThresholds;
    return [th.slopeModerate, th.slopeSaturated, th.minPairs, th.minSites, th.lowDivergence, satPlotPoints, satModel, satExcludeGaps].join('|');
  }

  void setSaturationParams({SaturationThresholds? thresholds, int? replicates, int? plotPoints, String? model, bool? excludeGaps}) {
    if (thresholds != null) saturationThresholds = thresholds;
    if (replicates != null) issReplicates = replicates;
    if (plotPoints != null) satPlotPoints = plotPoints;
    if (model != null && model.trim().isNotEmpty) satModel = model.trim();
    if (excludeGaps != null) satExcludeGaps = excludeGaps;
    notifyListeners();
  }

  void resetSaturationParams() {
    saturationThresholds = phykitThresholds;
    issReplicates = 400;
    satPlotPoints = 3000;
    satModel = 'GTR+G';
    satExcludeGaps = false;
    notifyListeners();
  }

  Future<void> runSaturation() {
    if (!_namesOk('Saturation test')) return Future.value();
    final coding = loci.values.where((l) => l.included && codonPositionsUsable(l) && _aln(l) != null).toList();
    final whole = _wholeSiteLoci;
    final steps = coding.length + whole.length;
    final modal = TaskModal(ModalKind.saturation, '${coding.length + whole.length}')..total = steps;
    for (final l in [...coding, ...whole]) {
      modal.items[l.name] = ModalItem();
    }
    return _task(modal: modal, 'Saturation test', (cancel) async {
        if (coding.isEmpty && whole.isEmpty) throw StateError('No aligned DNA loci to test. Align first.');
        final skipped = loci.values.where((l) => l.included && l.coding && !l.isProtein && _aln(l) != null && !codonPositionsUsable(l));
        if (skipped.isNotEmpty) {
          log('Saturation: ${skipped.map((l) => l.name).join(', ')} left out — trimming removed single columns, so codon positions are unknown',
              level: LogLevel.warn);
        }
        final outDir = dir('04_saturation');
        _satParamsUsed = _satParams;
        final th = saturationThresholds;
        final iqtree = await _require(ExternalTool.iqtree);
        final phykit = tools[ExternalTool.phykit] ?? await locator.resolve(ExternalTool.phykit);
        if (phykit == null) throw StateError('PhyKIT is not installed (step 5a: Install PhyKIT). ${ExternalTool.phykit.installHint}');
        tools[ExternalTool.phykit] = phykit;
        log('Saturation by PhyKIT ${phykit.version ?? ''} (Philippe et al. 2011): gene trees and branch lengths of each codon subset by IQ-TREE '
            '($satModel), thresholds: moderate below ${th.slopeModerate}, saturated below ${th.slopeSaturated}; low divergence p < ${th.lowDivergence}; '
            'min. ${th.minSites} sites and ${th.minPairs} pairs${satExcludeGaps ? '; sites with gaps excluded' : ''}', source: 'PhyKIT');
        var done = 0;
        String note(SaturationReport r) => r.subsets.entries.map((e) => '${e.key.label}: ${e.value.verdict.name}').join(' · ');
        final wholeSet = whole.toSet();
        await _pool([...coding, ...whole], math.max(1, Platform.numberOfProcessors - 2), cancel, (l) async {
          final aln = _aln(l)!;
          _item(l.name, state: ItemState.running, progress: 0.02, note: '${aln.length} seq · ${aln.first.length} bp');
          try {
            l.saturation = await _phykitSaturation(l, aln, wholeSet.contains(l) ? null : _frame(l), iqtree, phykit, outDir, cancel,
                (f, [note]) => _item(l.name, progress: f, note: note));
            _item(l.name, state: ItemState.done, progress: 1, note: note(l.saturation!));
          } catch (e) {
            log('${l.name}: $e', source: 'PhyKIT', error: true);
            _item(l.name, state: cancel.isCancelled ? ItemState.skipped : ItemState.failed, progress: 1, note: '$e');
          }
          _progress(++done, steps, l.name);
        }, label: (l) => l.name, reportProgress: false);
        if (cancel.isCancelled) return;
        if (coding.isEmpty) {
          await File(p.join(outDir, 'saturation_summary.tsv')).writeAsString(saturationTsv());
          final saturated = whole.where((l) => l.saturation?.subsets[CodonSubset.all]?.verdict == SaturationVerdict.saturated).length;
          log('Saturation (all sites): $saturated of ${whole.length} loci saturated');
          modal.result = 'All sites of ${whole.length} loci: $saturated saturated';
          notifyListeners();
          return;
        }
        await File(p.join(outDir, 'saturation_summary.tsv')).writeAsString(saturationTsv());
        final pooled = pooledSaturation;
        if (pooled == null) return;
        final rec = pooled.recommendation;
        log('Saturation: ${rec.kind.name}${rec.exclude.isEmpty ? '' : ' → exclude codon position(s) ${rec.exclude.join(', ')}'}');
        modal.result = 'Recommendation: ${switch (rec.kind) {
          RecommendationKind.keepAll => 'keep all codon positions',
          RecommendationKind.excludePositions => 'exclude codon position(s) ${rec.exclude.join('+')}',
          RecommendationKind.useAminoAcids => 'use amino acids',
          RecommendationKind.undetermined => 'undetermined',
        }}';
        notifyListeners();
      });
  }

  /// Saturation of one marker by PhyKIT: an IQ-TREE tree of all its sites
  /// ([satModel], fast search), the branch lengths of each codon subset on
  /// that topology (`-te`), then `phykit saturation -v --json` on each
  /// subset with its tree. [frame] null: non-coding, all sites only. Taxa
  /// are renamed S1, S2… in the files so every program reads the same
  /// names.
  Future<SaturationReport> _phykitSaturation(LocusState l, List<SeqRecord> aln, int? frame, ResolvedTool iqtree, ResolvedTool phykit,
      String outDir, CancelToken cancel, void Function(double f, [String? note]) stage) async {
    final work = Directory(p.join(outDir, sanitizeTaxonName(l.name)))..createSync(recursive: true);
    final length = aln.first.length;
    final ids = [for (var i = 0; i < aln.length; i++) 'S${i + 1}'];
    await File(p.join(work.path, 'taxa.tsv')).writeAsString([for (final (i, r) in aln.indexed) '${ids[i]}\t${r.header}'].join('\n'));
    final subsets = frame == null ? const [CodonSubset.all] : CodonSubset.values;
    Future<String> write(CodonSubset sub) async {
      final cols = subsetColumns(sub, length, frame ?? 0);
      final path = p.join(work.path, '${sub.name}.fasta');
      await File(path).writeAsString(writeFasta([
        for (final (i, r) in aln.indexed) SeqRecord(ids[i], String.fromCharCodes(cols.map((c) => c < r.sequence.length ? r.sequence.codeUnitAt(c) : 0x2d))),
      ]));
      return path;
    }

    Future<void> run(ResolvedTool tool, List<String> args, String what) async {
      final out = await runTool(tool, args, workingDirectory: work.path, cancel: cancel, captureStdout: tool == phykit, onLog: (_) {});
      if (!out.ok) throw StateError(out.cancelled ? 'cancelled' : '$what exit code ${out.exitCode}');
    }

    // PhyKIT's saturation value (`--json`) and its pairs (`-v --json`),
    // both saved next to the subset.
    Future<PhykitSaturation> saturation(CodonSubset sub, String fasta, String tree) async {
      Future<String> call(bool verbose) async {
        final out = await runTool(phykit, ['saturation', '-a', fasta, '-t', tree, if (verbose) '-v', '--json', if (satExcludeGaps) '-e'],
            workingDirectory: work.path, cancel: cancel, captureStdout: true, onLog: (_) {});
        if (!out.ok) throw StateError(out.cancelled ? 'cancelled' : 'PhyKIT exit code ${out.exitCode}');
        await File(p.join(work.path, '${sub.name}_phykit${verbose ? '_pairs' : ''}.json')).writeAsString(out.stdout);
        return out.stdout;
      }

      return parsePhykitSaturation(await call(false), await call(true));
    }

    // 1. Tree of all sites.
    stage(0.05, 'IQ-TREE $satModel');
    final allFasta = await write(CodonSubset.all);
    await run(iqtree, ['-s', allFasta, '-m', satModel, '-fast', '-nt', '1', '--prefix', 'all', '-redo', '-quiet'], 'IQ-TREE');
    final allTree = p.join(work.path, 'all.treefile');
    final runs = <CodonSubset, PhykitSaturation>{};
    final cells = <CodonSubset, int>{};
    final per = 0.9 / subsets.length;
    var at = 0.3;
    stage(at, 'PhyKIT');
    for (final sub in subsets) {
      if (cancel.isCancelled) throw StateError('cancelled');
      final fasta = sub == CodonSubset.all ? allFasta : await write(sub);
      var tree = allTree;
      if (sub != CodonSubset.all) {
        // 2. Branch lengths of this subset on the same topology.
        stage(at, 'IQ-TREE ${sub.label}');
        await run(iqtree, ['-s', fasta, '-te', allTree, '-m', satModel, '-nt', '1', '--prefix', sub.name, '-redo', '-quiet'], 'IQ-TREE');
        tree = p.join(work.path, '${sub.name}.treefile');
      }
      // 3. PhyKIT.
      stage(at + per * 0.5, 'PhyKIT ${sub.label}');
      runs[sub] = await saturation(sub, fasta, tree);
      cells[sub] = subsetColumns(sub, length, frame ?? 0).length;
      at = math.min(0.97, at + per * 0.7);
      stage(at);
    }
    final seqs = [for (final r in aln) r.sequence];
    final out = <CodonSubset, SubsetSaturation>{
      for (final sub in subsets)
        sub: () {
          final (variable, informative) = variableInformative(seqs, subsetColumns(sub, length, frame ?? 0));
          return phykitSubset(sub, cells[sub]!, runs[sub]!, saturationThresholds,
              model: satModel, maxPoints: satPlotPoints, variableSites: variable, informativeSites: informative);
        }(),
    };
    final positions = frame == null
        ? const <PositionPair>[]
        : phykitPositionPairs(runs[CodonSubset.all]!, runs[CodonSubset.pos1]!, runs[CodonSubset.pos2]!, runs[CodonSubset.pos3]!);
    return SaturationReport(l.name, aln.length, length, frame ?? 0, out, null, positions);
  }

  String saturationTsv() {
    final b = StringBuffer('locus\tsubset\tsites\tpairs\tmean_p\tmean_patristic\tphykit_saturation\tone_minus_saturation\tmodel\tverdict\t'
        'variable_sites\tinformative_sites\n');
    final reports = [
      ...loci.values.map((l) => l.saturation).whereType<SaturationReport>(),
    ];
    String f(double? v) => v == null ? 'NA' : v.toStringAsFixed(4);
    for (final r in reports) {
      for (final s in r.subsets.values) {
        b.writeln('${r.name}\t${s.subset.label}\t${s.sites}\t${s.pairs}\t${f(s.meanP)}\t${f(s.meanDistance)}\t'
            '${f(s.slope)}\t${f(s.slope == null ? null : (1 - s.slope!).abs())}\t${s.model}\t${s.verdict.name}\t'
            '${s.variableSites}\t${s.informativeSites}');
      }
    }
    return b.toString();
  }

  /// Applies the saturation recommendation to the concatenation settings,
  /// globally or per locus.
  void applySaturationRecommendation({required bool perLocus}) {
    concat = perLocus
        ? concat.copyWith(excludePositions: const {}, perLocusExclude: {
            for (final l in loci.values)
              if (l.saturation != null && l.saturation!.recommendation.exclude.isNotEmpty) l.name: l.saturation!.recommendation.exclude,
          })
        : concat.copyWith(perLocusExclude: const {}, excludePositions: pooledSaturation?.recommendation.exclude ?? const {});
    log('Concatenation will exclude: ${describeExclusions()}');
    notifyListeners();
  }

  /// Taxonomy, accession and voucher of every tip name the GenBank records
  /// produce (the FASTA header of 5a, with the `_2`, `_3`… suffixes of
  /// duplicated headers), for the tree viewer's labels and row colours.
  Map<String, TipInfo> tipInfo() {
    // Same records, header and taxonomy: the same map (the viewer keeps its
    // caches while the identity holds).
    final key = Object.hash(header.separator, Object.hashAll(header.fields),
        Object.hashAll(records.map((r) => '${r.accession}|${r.organismName}|${r.ranks.family}|${r.ranks.order}|${r.ranks.klass}|${r.editedFamily}|${r.salve?.code}')));
    if (key == _tipInfoKey && _tipInfo != null) return _tipInfo!;
    final out = <String, TipInfo>{};
    final count = <String, int>{};
    for (final r in records) {
      final h = header.build(r);
      final n = count.update(h, (v) => v + 1, ifAbsent: () => 1);
      final t = r.taxonomy;
      out[n == 1 ? h : '${h}_$n'] = TipInfo(
        species: t.binomial,
        genus: t.genus,
        family: t.family,
        order: t.order,
        klass: t.klass,
        phylum: t.phylum,
        kingdom: t.kingdom,
        accession: r.accession,
        voucher: (r.voucher ?? '').replaceFirst(RegExp(r'^(specimen\s+)?(voucher|isolate|strain|specimen)\s+', caseSensitive: false), ''),
        locality: r.locality ?? '',
        latLon: r.latLon ?? '',
        status: r.salve?.code ?? '',
        statusName: r.salve?.categoryName ?? '',
        commonName: (r.salve?.commonNames ?? '').split(',').first.trim(),
        endemic: r.salve?.endemic ?? '',
        biomes: r.salve?.biomes ?? '',
        trend: r.salve?.populationTrend ?? '',
        officialList: r.salve?.officialList ?? '',
        threats: r.salve?.threats.length ?? 0,
      );
    }
    _tipInfoKey = key;
    return _tipInfo = out;
  }

  Map<String, TipInfo>? _tipInfo;
  int? _tipInfoKey;

  String describeExclusions() {
    if (concat.excludePositions.isEmpty && concat.perLocusExclude.isEmpty) return 'nothing';
    final parts = <String>[
      if (concat.excludePositions.isNotEmpty) 'position(s) ${concat.excludePositions.join('+')} in all coding loci',
      for (final e in concat.perLocusExclude.entries) '${e.key}: ${e.value.join('+')}',
    ];
    return parts.join('; ');
  }

  // ---------------------------------------------------------------------
  // 7. Supermatrix
  // ---------------------------------------------------------------------

  /// Concatenate the trimmed alignments (recommended) or the MAFFT output.
  bool useTrimmed = true;

  void setUseTrimmed(bool v) {
    useTrimmed = v;
    matrix = null;
    overview = null;
    notifyListeners();
  }

  List<SeqRecord>? _aln(LocusState l) => useTrimmed ? l.alignment : l.aligned;
  int _frame(LocusState l) => useTrimmed ? l.alignmentFrame : l.frame;

  /// Codon positions of the alignment used downstream are known.
  bool codonPositionsUsable(LocusState l) => useTrimmed ? l.codonPositionsKnown : l.coding && !l.isProtein;

  bool get proteinData => loci.values.any((l) => l.included && l.isProtein && l.aligned != null);

  /// Outgroup suggested from the alignment distances (step 7c), computed
  /// once per set of alignments.
  OutgroupSuggestion? outgroupSuggestion;
  bool outgroupSuggesting = false;
  int? _outgroupKey;

  int get _alignmentsKey => Object.hashAll([
        useTrimmed,
        for (final l in loci.values)
          if (l.included && _aln(l) != null) identityHashCode(_aln(l)),
      ]);

  /// Computes [outgroupSuggestion] when the alignments changed.
  Future<void> suggestOutgroups() async {
    final key = _alignmentsKey;
    if (outgroupSuggesting || key == _outgroupKey) return;
    final input = [
      for (final l in loci.values)
        if (l.included && _aln(l) != null) LocusAlignment(l.name, _aln(l)!),
    ];
    if (input.isEmpty) return;
    outgroupSuggesting = true;
    notifyListeners();
    try {
      outgroupSuggestion = await Isolate.run(() => suggestOutgroup(input));
      _outgroupKey = key;
      final s = outgroupSuggestion;
      if (s != null) {
        log('Outgroup suggested from the p-distances of ${s.loci} alignment(s): ${s.taxa.join(', ')} '
            '(mean distance to the ingroup ${s.between.toStringAsFixed(3)}, within the ingroup ${s.within.toStringAsFixed(3)})'
            '${s.confident ? '' : ' — no clearly separate group; the most distant taxon is offered'}');
      }
    } catch (e) {
      log('Outgroup suggestion failed: $e', level: LogLevel.warn);
    } finally {
      outgroupSuggesting = false;
      notifyListeners();
    }
  }

  void setConcat(ConcatOptions o) {
    concat = o;
    notifyListeners();
  }

  /// The partition schemes that mean something for these alignments, the
  /// recommended first:
  /// - one marker: by codon position when its codon positions are known
  ///   (by locus would be a single partition), else one partition;
  /// - several markers without known codon positions (non-coding, protein,
  ///   fragments without a confirmed frame, frame lost in trimming): by
  ///   locus, or none;
  /// - several markers, some coding: every scheme.
  List<PartitionScheme> get partitionSchemes {
    final inc = loci.values.where((l) => l.included && (l.aligned != null)).toList();
    final codons = inc.any((l) => codonPositionsUsable(l));
    if (inc.length <= 1) return codons ? const [PartitionScheme.byCodonPosition, PartitionScheme.none] : const [PartitionScheme.none];
    return codons
        ? const [PartitionScheme.byLocusAndCodon, PartitionScheme.byLocus, PartitionScheme.byCodonPosition, PartitionScheme.none]
        : const [PartitionScheme.byLocus, PartitionScheme.none];
  }

  /// The scheme used: the one chosen when it applies, else the recommended
  /// one for these alignments (one gene chosen "by locus and codon" is
  /// partitioned by codon position).
  PartitionScheme get effectiveScheme {
    final allowed = partitionSchemes;
    if (allowed.contains(concat.scheme)) return concat.scheme;
    // The nearest meaning: codons kept when asked for and possible.
    if (concat.scheme == PartitionScheme.byLocusAndCodon && allowed.contains(PartitionScheme.byCodonPosition)) return PartitionScheme.byCodonPosition;
    if (concat.scheme == PartitionScheme.byCodonPosition && allowed.contains(PartitionScheme.byLocus)) return PartitionScheme.byLocus;
    return allowed.first;
  }

  Future<void> buildMatrix() async {
    if (!_namesOk('Supermatrix')) return;
    await _buildMatrix();
  }

  Future<void> _buildMatrix() => _task(modal: TaskModal(ModalKind.matrix, '')..total = 2, 'Building supermatrix', (_) async {
        final included = loci.values.where((l) => l.included && _aln(l) != null).toList();
        if (included.isEmpty) throw StateError('No alignments available.');
        final protein = proteinData;
        if (protein && included.any((l) => !l.isProtein)) {
          throw StateError('Cannot concatenate protein and DNA alignments together. Re-align all markers in one mode.');
        }
        final scheme = effectiveScheme;
        if (scheme != concat.scheme) log('Partitions: ${concat.scheme.name} does not apply to these alignments; ${scheme.name} is used');
        final options = concat.copyWith(dataType: protein ? SeqDataType.protein : SeqDataType.dna, scheme: scheme);
        final m = await _matrixInIsolate([
          // Loci whose frame was lost in trimming are one partition.
          for (final l in included) LocusAlignment(l.name, _aln(l)!, coding: codonPositionsUsable(l), frame: _frame(l)),
        ], options);
        matrix = m;
        _progress(1, 2, 'writing PHYLIP, FASTA, NEXUS and partition files');
        final out = dir('05_matrix');
        final files = {
          'phylip': p.join(out, 'concatenated.phy'),
          'fasta': p.join(out, 'concatenated.fasta'),
          'nexus': p.join(out, 'concatenated.nex'),
          'partition': p.join(out, 'partitions.txt'),
          'raxml': p.join(out, 'partitions_raxml.txt'),
          'iqtree': p.join(out, 'partitions_iqtree.nex'),
        };
        await File(files['phylip']!).writeAsString(m.toPhylip());
        await File(files['fasta']!).writeAsString(m.toFasta());
        await File(files['nexus']!).writeAsString(m.toNexus());
        await File(files['partition']!).writeAsString(m.toPartitionTable());
        await File(files['raxml']!).writeAsString(m.toRaxmlPartitions());
        await File(files['iqtree']!).writeAsString(m.toIqtreePartitions());
        // Bases and gaps of every taxon in every locus (step 7 table).
        files['overview'] = p.join(out, 'matrix_overview.tsv');
        overview = await _overviewInIsolate(m);
        await File(files['overview']!).writeAsString(overview!.toTsv());
        matrixFiles
          ..clear()
          ..addAll(files);
        _progress(2, 2);
        modal?.result = '${m.taxa.length} taxa × ${m.length} sites · ${m.partitions.length} partition(s)';
        log('Supermatrix: ${m.taxa.length} taxa × ${m.length} sites, ${m.partitions.length} partition(s), '
            '${(m.missingFraction * 100).toStringAsFixed(1)}% missing. Excluded: ${describeExclusions()}', level: LogLevel.success);
        if (m.droppedLoci.isNotEmpty) log('Dropped loci: ${m.droppedLoci.join(', ')}', level: LogLevel.warn);
        if (m.droppedTaxa.isNotEmpty) log('Dropped taxa: ${m.droppedTaxa.join(', ')}', level: LogLevel.warn);
        notifyListeners();
      });

  Future<void> saveMatrixFile(String key, String target) async {
    final src = matrixFiles[key];
    if (src == null) return;
    await File(src).copy(target);
    log('Saved $target', level: LogLevel.success);
  }

  // ---------------------------------------------------------------------
  // 8. Trees
  // ---------------------------------------------------------------------

  String _runDir(String engine) {
    final stamp = DateTime.now().toIso8601String().replaceAll(RegExp(r'[:.]'), '-').substring(0, 19);
    return dir(p.join('06_trees', '${engine}_$stamp'));
  }

  int get _partitionCount => _partitioned ? matrix!.partitions.length : 0;

  bool get _partitioned => matrix != null && matrix!.partitions.isNotEmpty && effectiveScheme != PartitionScheme.none;

  String iqtreeCommandPreview(IqtreeOptions options) {
    final exe = tools[ExternalTool.iqtree]?.executable ?? 'iqtree3';
    return multilineCommand(
      p.basename(exe),
      options.args(
        alignment: matrixFiles['phylip'] ?? 'concatenated.phy',
        partitions: _partitioned ? (matrixFiles['iqtree'] ?? 'partitions_iqtree.nex') : null,
        prefix: '06_trees/iqtree/splace',
        partitionCount: _partitionCount,
      ),
    );
  }

  /// The RAxML-NG tree search as it will run (the `--parse` check and,
  /// with ModelFinder, IQ-TREE's model selection come before it).
  String raxmlCommandPreview(RaxmlOptions options, {required RaxmlModels models, String? fixedModel}) {
    final exe = tools[ExternalTool.raxmlng]?.executable ?? 'raxml-ng';
    const prefix = '06_trees/raxmlng/splace';
    final args = switch (models) {
      RaxmlModels.moose => options.allArgs(
          msa: matrixFiles['phylip'] ?? 'concatenated.phy',
          model: _partitioned ? 'auto=${matrixFiles['raxml'] ?? 'partitions_raxml.txt'}' : (matrix?.isProtein ?? false ? 'AA' : 'DNA'),
          prefix: prefix,
          moose: true,
          dataType: matrix?.isProtein ?? false ? 'AA' : 'DNA'),
      _ => options.allArgs(msa: '$prefix.raxml.rba', prefix: prefix),
    };
    return multilineCommand(p.basename(exe), args);
  }

  /// Feeds one line of IQ-TREE/RAxML-NG output to the progress modal.
  void _treeLine(String line, [TreeProgress? progress, String? item]) {
    final m = modal;
    if (m == null) return;
    m.addTail(line, item);
    final p = progress ?? m.tree;
    if (p == null) return;
    if (p.add(line) || item == null) {
      if (item != null) _item(item, progress: p.fraction);
      notifyListeners();
    }
  }

  static bool _usesModelFinder(String model) => RegExp(r'(^|\+)(MFP|MF|TEST\w*)(\+|$)').hasMatch(model) || model.isEmpty;

  List<SeqRecord> _geneTreeInput(LocusState l) => _aln(l)!;

  TreeProgress _iqProgress(IqtreeOptions o) => TreeProgress(stages: [
        TreeStage.prepare,
        if (_usesModelFinder(o.model)) TreeStage.modelFinder,
        TreeStage.search,
        if (o.ufboot > 0 || o.alrt > 0 || o.classicBootstrap > 0) TreeStage.support,
        TreeStage.finish,
      ]);

  Future<void> runIqtree(IqtreeOptions options, {bool perGene = false}) {
    if (!_namesOk('IQ-TREE')) return Future.value();
    final modal = TaskModal(ModalKind.iqtree, options.model)..tree = _iqProgress(options);
    if (perGene) {
      for (final l in loci.values.where((l) => l.included && _aln(l) != null && _geneTreeInput(l).length >= 4)) {
        modal.items[l.name] = ModalItem();
      }
    }
    return _task(modal: modal, 'IQ-TREE', (cancel) async {
        final run = await _iqtreeRun(options, cancel, modal.tree!);
        final t = modal.tree!;
        modal.result = 'ML tree done${t.bestModel == null ? '' : ' · ${t.bestModel}'}${t.bestLogL == null ? '' : ' · lnL ${t.bestLogL}'}';
        notifyListeners();
        if (perGene && !cancel.isCancelled) await _geneTrees(await _require(ExternalTool.iqtree), run, options, cancel);
      });
  }

  /// One IQ-TREE analysis of the supermatrix; [item] names the line of a
  /// shared modal that shows its progress (both programs at once).
  Future<TreeRun> _iqtreeRun(IqtreeOptions options, CancelToken cancel, TreeProgress progress, {String? item}) async {
    final tool = await _require(ExternalTool.iqtree);
    if (matrix == null) throw StateError('Build the supermatrix first.');
    final run = TreeRun('IQ-TREE', _runDir('iqtree'))..iqtree = options;
    treeRuns.insert(0, run);
    notifyListeners();
    return _failRun(run, () => _iqtreeBody(tool, run, options, cancel, progress, item: item));
  }

  /// Marks [run] as failed, with the reason, when [body] throws.
  Future<TreeRun> _failRun(TreeRun run, Future<TreeRun> Function() body) async {
    try {
      return await body();
    } catch (e) {
      run
        ..status = StepStatus.error
        ..error ??= _errorText(e);
      notifyListeners();
      rethrow;
    }
  }

  static String _errorText(Object e) => e is StateError ? e.message : '$e';

  Future<TreeRun> _iqtreeBody(ResolvedTool tool, TreeRun run, IqtreeOptions options, CancelToken cancel, TreeProgress progress, {String? item}) async {
    final args = options.args(
      alignment: matrixFiles['phylip']!,
      partitions: _partitioned ? matrixFiles['iqtree'] : null,
      prefix: p.join(run.directory, 'splace'),
      partitionCount: _partitionCount,
    );
    run.command = '${tool.description} ${args.join(' ')}';
    run.commandDisplay = multilineCommand(tool.description, args);
    final out = await runTool(tool, args, workingDirectory: run.directory, cancel: cancel, onLog: (m) {
      log(m, source: 'IQ-TREE');
      run.addTail(m);
      _treeLine(m, progress, item);
    });
    progress.finish();
    _collectOutputs(run, {
      'tree': 'splace.treefile',
      'consensus': 'splace.contree',
      'report': 'splace.iqtree',
      'log': 'splace.log',
    });
    run.status = out.ok && run.treeFile != null ? StepStatus.done : StepStatus.error;
    notifyListeners();
    if (run.status == StepStatus.error) {
      throw StateError(out.cancelled ? 'cancelled' : (out.ok ? 'IQ-TREE wrote no tree' : 'IQ-TREE ${out.describe}'));
    }
    run.bestModel = progress.bestModel;
    return run;
  }

  /// IQ-TREE and RAxML-NG at the same time on the supermatrix, each with
  /// half of the threads, so the computer is not overloaded.
  Future<void> runBoth(IqtreeOptions iq, RaxmlOptions rx, {RaxmlModels models = RaxmlModels.moose, String? fixedModel}) {
    if (!_namesOk('Tree inference')) return Future.value();
    final total = math.max(2, threads);
    final forIq = total ~/ 2, forRx = total - total ~/ 2;
    final iqo = iq.copyWith(threads: '$forIq', maxThreads: forIq);
    final rxo = rx.copyWith(maxThreads: forRx);
    final modal = TaskModal(ModalKind.bothTrees, '$forIq + $forRx')
      ..items['IQ-TREE'] = ModalItem()
      ..items['RAxML-NG'] = ModalItem();
    return _task(modal: modal, 'IQ-TREE + RAxML-NG', (cancel) async {
        log('Running IQ-TREE ($forIq threads) and RAxML-NG ($forRx threads) together', source: 'Trees');
        final pi = _iqProgress(iqo), pr = _rxProgress(rxo, models);
        modal.trees
          ..['IQ-TREE'] = pi
          ..['RAxML-NG'] = pr;
        notifyListeners();
        var failed = 0;
        Future<void> one(String name, TreeProgress p, Future<TreeRun> Function() body) async {
          _item(name, state: ItemState.running, progress: 0.01);
          try {
            await body();
            _item(name, state: ItemState.done, progress: 1, note: [?p.bestModel, if (p.bestLogL != null) 'lnL ${p.bestLogL}'].join(' · '));
          } catch (e) {
            failed++;
            _item(name, state: cancel.isCancelled ? ItemState.skipped : ItemState.failed, note: _errorText(e));
            notifyListeners();
          }
        }

        await Future.wait([
          one('IQ-TREE', pi, () => _iqtreeRun(iqo, cancel, pi, item: 'IQ-TREE')),
          one('RAxML-NG', pr, () => _raxmlRun(rxo, cancel, pr, models: models, fixedModel: fixedModel, item: 'RAxML-NG', threads: forRx)),
        ]);
        if (failed == 2) throw StateError(cancel.isCancelled ? 'cancelled' : 'both analyses failed');
        modal.result = failed == 0 ? 'Both trees done' : 'One analysis failed';
        notifyListeners();
      });
  }

  /// "Also infer per-gene trees": one IQ-TREE run per marker alignment.
  Future<void> _geneTrees(ResolvedTool tool, TreeRun run, IqtreeOptions options, CancelToken cancel) async {
    final genesDir = p.join(run.directory, 'gene_trees');
    Directory(genesDir).createSync(recursive: true);
    final genes = loci.values.where((l) => l.included && _aln(l) != null && _aln(l)!.length >= 4).toList();
    var done = 0;
    for (final l in genes) {
      if (cancel.isCancelled) break;
      final name = sanitizeTaxonName(l.name);
      final fasta = p.join(genesDir, '$name.fasta');
      await File(fasta).writeAsString(writeFasta(_aln(l)!));
      final geneOptions = IqtreeOptions(
        model: options.model.contains('MERGE') ? options.model.replaceAll('+MERGE', '') : options.model,
        ufboot: options.ufboot,
        alrt: options.alrt,
        classicBootstrap: options.classicBootstrap,
        threads: options.threads,
        maxThreads: options.maxThreads,
        seed: options.seed,
        outgroup: options.outgroup.where((o) => _aln(l)!.any((r) => r.header == o)).toList(),
        bnni: options.bnni,
        nstop: options.nstop,
        perturbation: options.perturbation,
        extra: options.extra,
      );
      final progress = TreeProgress(stages: [
        TreeStage.prepare,
        if (_usesModelFinder(geneOptions.model)) TreeStage.modelFinder,
        TreeStage.search,
        if (geneOptions.ufboot > 0 || geneOptions.alrt > 0 || geneOptions.classicBootstrap > 0) TreeStage.support,
        TreeStage.finish,
      ]);
      _item(l.name, state: ItemState.running, progress: 0.01, note: 'IQ-TREE');
      final out = await runTool(tool, geneOptions.args(alignment: fasta, prefix: p.join(genesDir, '${name}_tree')),
          workingDirectory: genesDir, cancel: cancel, onLog: (m) {
        log(m, source: 'IQ-TREE ${l.name}');
        _treeLine(m, progress, l.name);
      });
      run.geneTrees[l.name] = out.ok;
      _item(l.name,
          progress: 1,
          state: out.ok ? ItemState.done : (cancel.isCancelled ? ItemState.skipped : ItemState.failed),
          note: out.ok ? (progress.bestModel ?? 'ML tree') : 'exit code ${out.exitCode}');
      _progress(++done, genes.length, l.name);
      notifyListeners();
    }
    log('Per-gene trees: ${run.geneTrees.values.where((v) => v).length}/${genes.length} → $genesDir',
        source: 'IQ-TREE', level: LogLevel.success);
  }

  /// RAxML-NG with the models chosen by its own MOOSE (RAxML-NG 2), by
  /// IQ-TREE's ModelFinder (restricted to models RAxML-NG supports) or one
  /// fixed model for every partition.
  Future<void> runRaxml(RaxmlOptions options, {RaxmlModels models = RaxmlModels.moose, String? fixedModel}) {
    if (!_namesOk('RAxML-NG')) return Future.value();
    final modal = TaskModal(ModalKind.raxml, switch (models) {
      RaxmlModels.moose => 'MOOSE',
      RaxmlModels.modelFinder => 'ModelFinder',
      RaxmlModels.fixed => fixedModel ?? 'GTR+G',
    })
      ..tree = _rxProgress(options, models);
    return _task(modal: modal, 'RAxML-NG', (cancel) async {
        await _raxmlRun(options, cancel, modal.tree!, models: models, fixedModel: fixedModel);
        final t = modal.tree!;
        modal.result = 'ML tree done${t.bestLogL == null ? '' : ' · lnL ${t.bestLogL}'}${t.bootDone == null ? '' : ' · ${t.bootDone} bootstrap replicates'}';
        notifyListeners();
      });
  }

  TreeProgress _rxProgress(RaxmlOptions options, RaxmlModels models) {
    final boots = options.bootstraps.trim();
    final support = options.needsReplicates && boots.isNotEmpty && boots != '0';
    return TreeProgress(raxml: true, stages: [
      TreeStage.prepare,
      if (models != RaxmlModels.fixed) TreeStage.modelFinder,
      TreeStage.search,
      if (support || options.metrics.isNotEmpty) TreeStage.support,
      TreeStage.finish,
    ]);
  }

  /// One RAxML-NG analysis; [item] names its line in a shared modal and
  /// [threads] caps IQ-TREE's ModelFinder when both programs run together.
  Future<TreeRun> _raxmlRun(RaxmlOptions options, CancelToken cancel, TreeProgress progress,
      {required RaxmlModels models, String? fixedModel, String? item, int? threads}) async {
    final raxml = await _require(ExternalTool.raxmlng);
    if (matrix == null) throw StateError('Build the supermatrix first.');
    final major = int.tryParse(raxml.version?.split('.').first ?? '');
    if (major != null && major < 2 && (models == RaxmlModels.moose || options.fast || options.metrics.any((m) => m.index > 1))) {
      throw StateError('RAxML-NG ${raxml.version} found (${raxml.description}); MOOSE, --fast and the SH-aLRT/EBG/RBS supports need RAxML-NG 2.0 or newer.');
    }
    final defaultModel = fixedModel ?? (matrix!.isProtein ? 'LG+G' : 'GTR+G');
    final run = TreeRun('RAxML-NG', _runDir('raxmlng'))
      ..raxml = options
      ..models = models
      ..fixedModel = models == RaxmlModels.fixed ? defaultModel : null;
    treeRuns.insert(0, run);
    notifyListeners();
    return _failRun(run, () => _raxmlBody(raxml, run, options, cancel, progress,
        models: models, defaultModel: defaultModel, item: item, threads: threads));
  }

  Future<TreeRun> _raxmlBody(ResolvedTool raxml, TreeRun run, RaxmlOptions options, CancelToken cancel, TreeProgress progress,
      {required RaxmlModels models, required String defaultModel, String? item, int? threads}) async {
    void rlog(String m) {
      log(m, source: 'RAxML-NG');
      run.addTail(m);
    }

    var model = _partitioned ? matrixFiles['raxml']! : defaultModel;
    if (_partitioned && models == RaxmlModels.fixed) {
      model = p.join(run.directory, 'partitions_fixed.txt');
      await File(model).writeAsString(matrix!.toRaxmlPartitions(defaultModel: defaultModel));
    }

    if (models == RaxmlModels.modelFinder) {
      final iq = await _require(ExternalTool.iqtree);
      final mfPrefix = p.join(run.directory, 'modelfinder');
      rlog('Selecting partitioning scheme and models with IQ-TREE ModelFinder (--mset raxml)…');
      final mf = await runTool(
        iq,
        IqtreeOptions(model: _partitioned ? 'MF+MERGE' : 'MF', raxmlModelsOnly: true, threads: '${threads ?? this.threads}')
            .args(alignment: matrixFiles['phylip']!, partitions: _partitioned ? matrixFiles['iqtree'] : null, prefix: mfPrefix, partitionCount: _partitionCount),
        workingDirectory: run.directory,
        cancel: cancel,
        onLog: (m) {
          log(m, source: 'ModelFinder');
          run.addTail(m);
          _treeLine(m, progress, item);
        },
      );
      if (!mf.ok) throw StateError(mf.cancelled ? 'cancelled' : 'ModelFinder failed (${mf.describe}).');
      final schemeFile = File('$mfPrefix.best_scheme.nex');
      if (schemeFile.existsSync()) {
        final converted = schemeToRaxml(parseIqtreeNexusScheme(await schemeFile.readAsString()));
        model = p.join(run.directory, 'partitions_modelfinder.txt');
        await File(model).writeAsString(converted.partitions);
        run.notes.addAll(converted.replaced);
        converted.replaced.forEach(rlog);
        rlog('Scheme for RAxML-NG:\n${converted.partitions}');
      } else {
        // Unpartitioned: read the single best model from the report.
        final report = await File('$mfPrefix.iqtree').readAsString();
        final best = RegExp(r'Best-fit model according to \w+:\s*(\S+)').firstMatch(report)?.group(1);
        model = (best == null ? null : iqtreeModelToRaxml(best)) ?? defaultModel;
        rlog('Best model: $best → RAxML-NG: $model');
      }
    }

    // MOOSE: `DNA`/`AA` for one partition, `auto=<file>` for several
    // (the models in the file are replaced by the selected ones).
    final moose = models == RaxmlModels.moose;
    if (moose) model = _partitioned ? 'auto=$model' : (matrix!.isProtein ? 'AA' : 'DNA');

    final prefix = p.join(run.directory, 'splace');
    final dataType = matrix!.isProtein ? 'AA' : 'DNA';
    final parse = await runTool(raxml, options.parseArgs(msa: matrixFiles['phylip']!, model: model, prefix: prefix, dataType: moose ? dataType : null),
        workingDirectory: run.directory, cancel: cancel, onLog: rlog);
    if (!parse.ok) throw StateError(parse.cancelled ? 'cancelled' : 'RAxML-NG could not read the alignment/partitions (${parse.describe}).');
    final logFile = File('$prefix.raxml.log');
    if (logFile.existsSync()) {
      final s = parseRaxmlLog(await logFile.readAsString());
      rlog('RAxML-NG recommends ${s.recommendedThreads ?? '?'} thread(s); memory ${s.memory ?? '?'}');
    }

    // The .rba keeps the models it was parsed with, so MOOSE reads the
    // PHYLIP file again.
    final allArgs = moose
        ? options.allArgs(msa: matrixFiles['phylip']!, model: model, prefix: prefix, moose: true, dataType: dataType)
        : options.allArgs(msa: '$prefix.raxml.rba', prefix: prefix);
    run.command = '${raxml.description} ${allArgs.join(' ')}';
    run.commandDisplay = multilineCommand(raxml.description, allArgs);
    final output = StringBuffer();
    void onLog(String m) {
      rlog(m);
      output.writeln(m);
      _treeLine(m, progress, item);
    }

    var all = await runTool(raxml, allArgs, workingDirectory: run.directory, cancel: cancel, onLog: onLog);
    // A crash (seen once after 40 min of bootstrapping, 0xC0000005 on
    // Windows) does not lose the work: RAxML-NG resumes from its
    // checkpoint when run again without --redo.
    for (var attempt = 1; !all.ok && !all.cancelled && attempt <= 2 && File('$prefix.raxml.ckp').existsSync(); attempt++) {
      final why = 'RAxML-NG stopped: ${all.describe}. Resuming from its checkpoint (attempt $attempt of 2)…';
      log(why, source: 'RAxML-NG', level: LogLevel.warn);
      run
        ..addTail(why)
        ..notes.add('Resumed from checkpoint after ${all.describe}');
      all = await runTool(raxml, [for (final a in allArgs) if (a != '--redo') a],
          workingDirectory: run.directory, cancel: cancel, onLog: onLog);
    }
    progress.finish();
    _collectOutputs(run, {
      'tree': 'splace.raxml.support',
      for (final m in RaxmlSupport.values) m.outputKey: 'splace.raxml.${m.suffix}',
      'bestTree': 'splace.raxml.bestTree',
      'bootstraps': 'splace.raxml.bootstraps',
      'model': 'splace.raxml.bestModel',
      'log': 'splace.raxml.log',
    });
    run.treeFile = run.outputs['tree'] ??
        [for (final m in options.metrics) run.outputs[m.outputKey]].nonNulls.firstOrNull ??
        run.outputs['bestTree'];
    // Majority-rule consensus of the bootstrap replicates (the tree step 9
    // shows first, as IQ-TREE's .contree).
    final reps = run.outputs['bootstraps'];
    if (all.ok && reps != null && !cancel.isCancelled) {
      final cPrefix = p.join(run.directory, 'splace_consensus');
      final cons = await runTool(raxml, ['--consense', 'MR', '--tree', reps, '--prefix', cPrefix, '--redo'],
          workingDirectory: run.directory, cancel: cancel, onLog: rlog);
      final file = File('$cPrefix.raxml.consensusTreeMR');
      if (cons.ok && file.existsSync()) {
        run.outputs['consensus'] = file.path;
        run.notes.add('Majority-rule consensus of the bootstrap replicates (--consense MR)');
      } else {
        rlog('RAxML-NG could not build the consensus tree (${cons.describe})');
      }
    }
    if (moose) {
      final chosen = parseMooseModels(output.toString());
      run.mooseModels.addAll(chosen);
      if (chosen.length == 1) {
        run.bestModel = chosen.single;
        run.notes.add('MOOSE: ${chosen.single}');
      } else if (chosen.isNotEmpty) {
        run.notes.add('MOOSE: ${chosen.length} partition models');
      }
    }
    if (run.outputs['log'] != null) {
      final s = parseRaxmlLog(await File(run.outputs['log']!).readAsString());
      if (s.finalLogLikelihood != null) run.notes.add('Final lnL ${s.finalLogLikelihood}');
      if (s.bootstrapTrees != null) run.notes.add('Bootstrap converged after ${s.bootstrapTrees} replicates');
    }
    run.status = all.ok && run.treeFile != null ? StepStatus.done : StepStatus.error;
    if (run.status == StepStatus.error) {
      throw StateError(all.cancelled ? 'cancelled' : (all.ok ? 'RAxML-NG wrote no tree' : 'RAxML-NG ${all.describe}'));
    }
    notifyListeners();
    return run;
  }

  /// Trees that can be shown in step 9: (run, kind, gene, path). Kind is
  /// `tree`, `consensus`, `treeFBP`/`treeTBE`/`treeRBS`/`treeSH`/`treeEBG`,
  /// `bestTree` or `gene`.
  List<({TreeRun run, String kind, String? gene, String path})> get viewableTrees => [
        for (final r in treeRuns.where((r) => r.status == StepStatus.done)) ...[
          for (final k in ['consensus', 'tree', for (final m in RaxmlSupport.values) m.outputKey, 'bestTree'])
            if (r.outputs[k] != null) (run: r, kind: k, gene: null, path: r.outputs[k]!),
          for (final g in r.geneTrees.entries)
            if (g.value && File(p.join(r.directory, 'gene_trees', '${sanitizeTaxonName(g.key)}_tree.treefile')).existsSync())
              (run: r, kind: 'gene', gene: g.key, path: p.join(r.directory, 'gene_trees', '${sanitizeTaxonName(g.key)}_tree.treefile')),
        ],
      ];

  // ---------------------------------------------------------------------
  // 8e. Species delimitation (mPTP on trees; ASAP and ABGD on alignments)
  // ---------------------------------------------------------------------

  /// Delimitations, newest first (one per tree or alignment and method).
  final List<DelimitationRun> delimitations = [];

  /// Finished delimitations, for the buttons of the tree viewer.
  List<DelimitationRun> get delimitationsDone => [for (final d in delimitations) if (d.status == StepStatus.done && d.result != null) d];

  /// ML tree of a finished run (branch lengths in substitutions per site):
  /// IQ-TREE's `.treefile`, RAxML-NG's best tree.
  String? mlTreeOf(TreeRun r) => r.engine == 'RAxML-NG' ? (r.outputs['bestTree'] ?? r.treeFile) : (r.outputs['tree'] ?? r.treeFile);

  /// Every tip of the supermatrix is a different species (one sequence
  /// per species): PTP, mPTP and GMYC look for the change from branching
  /// between species to branching within them, so without several
  /// individuals of a species their test tends to fail.
  bool get oneSequencePerSpecies {
    final taxa = matrix?.taxa;
    if (taxa == null || taxa.length < 4) return false;
    final info = tipInfo();
    final seen = <String>{};
    for (final t in taxa) {
      final s = info[t]?.species ?? '';
      if (s.isEmpty || !seen.add(s.toLowerCase())) return false;
    }
    return true;
  }

  /// Runs whose ML tree can be delimited.
  List<TreeRun> get delimitableRuns => [
        for (final r in treeRuns)
          if (r.status == StepStatus.done && mlTreeOf(r) != null && File(mlTreeOf(r)!).existsSync()) r,
      ];

  /// Name of the supermatrix among the alignments ASAP and ABGD can read.
  static const supermatrixAlignment = 'supermatrix';

  /// DNA alignments ASAP and ABGD can read: the supermatrix and each
  /// included marker (as used downstream: trimmed or MAFFT).
  List<String> get delimitableAlignments => [
        if (matrix != null && !matrix!.isProtein && matrixFiles['fasta'] != null) supermatrixAlignment,
        for (final l in loci.values)
          if (l.included && !l.isProtein && _aln(l) != null && _aln(l)!.length >= 3) l.name,
      ];

  Future<List<(String, String)>> _alignmentForDelimitation(String name) async {
    if (name == supermatrixAlignment) {
      return [for (final r in parseFasta(await File(matrixFiles['fasta']!).readAsString())) (r.header.split(RegExp(r'\s')).first, r.sequence)];
    }
    return [for (final r in _aln(loci[name]!)!) (r.header.split(RegExp(r'\s')).first, r.sequence)];
  }

  /// Delimits species on the ML tree of each of [runs] with mPTP, with
  /// each of [models]. With [autoMinbr], the minimum branch length comes
  /// from the p-distances of the supermatrix (`--minbr_auto`); with
  /// [cropOutgroup], the outgroups root the tree and are left out.
  /// With [alignments] and [distanceMethods] (ASAP, ABGD), each alignment
  /// is also delimited from its pairwise [distance]s.
  Future<void> runDelimitation({
    List<TreeRun> runs = const [],
    Set<PtpModel> models = const {},
    bool gmyc = false,
    double gmycPvalue = gmycDefaultPvalue,
    double mptpPvalue = mptpDefaultPvalue,
    bool autoMinbr = true,
    bool cropOutgroup = true,
    List<String> alignments = const [],
    Set<DelimMethod> distanceMethods = const {},
    DelimDistance distance = DelimDistance.jc69,
    double abgdPMin = 0.001,
    double abgdPMax = 0.1,
    int abgdSteps = 10,
    double abgdSlope = 1.5,
  }) {
    if (!_namesOk('Species delimitation')) return Future.value();
    final methods = [for (final m in PtpModel.values) if (models.contains(m)) m.label, if (gmyc) 'GMYC', for (final m in distanceMethods) m.label];
    final modal = TaskModal(ModalKind.delimitation, methods.join(' + '));
    for (final r in runs) {
      for (final m in PtpModel.values.where(models.contains)) {
        modal.items['${m.label} · ${r.engine} ${_runTime(r)}'] = ModalItem();
      }
      if (gmyc) modal.items['GMYC · ${r.engine} ${_runTime(r)}'] = ModalItem();
    }
    for (final a in alignments) {
      for (final m in distanceMethods) {
        modal.items['${m.label} · $a'] = ModalItem();
      }
    }
    return _task(modal: modal, 'Species delimitation', (cancel) async {
      var failed = 0;
      if (runs.isNotEmpty && models.isNotEmpty) {
        final tool = await _require(ExternalTool.mptp);
        for (final r in runs) {
          if (cancel.isCancelled) break;
          try {
            await _delimitRun(tool, r, models, autoMinbr: autoMinbr, cropOutgroup: cropOutgroup, pvalue: mptpPvalue, cancel: cancel);
          } catch (e) {
            failed++;
            log('mPTP on the ${r.engine} tree failed: ${_errorText(e)}', error: true);
          }
        }
      }
      if (runs.isNotEmpty && gmyc) {
        for (final r in runs) {
          if (cancel.isCancelled) break;
          try {
            await _gmycRun(r, gmycPvalue, cropOutgroup: cropOutgroup, cancel: cancel);
          } catch (e) {
            failed++;
            _item('GMYC · ${r.engine} ${_runTime(r)}', state: ItemState.failed, note: _errorText(e));
            log('GMYC on the ${r.engine} tree failed: ${_errorText(e)}', error: true);
          }
        }
      }
      for (final a in alignments) {
        for (final m in [DelimMethod.asap, DelimMethod.abgd].where(distanceMethods.contains)) {
          if (cancel.isCancelled) break;
          try {
            await _delimitAlignment(a, m, distance, cancel,
                pMin: abgdPMin, pMax: abgdPMax, steps: abgdSteps, slope: abgdSlope);
          } catch (e) {
            failed++;
            _item('${m.label} · $a', state: ItemState.failed, note: _errorText(e));
            log('${m.label} on $a failed: ${_errorText(e)}', error: true);
          }
        }
      }
      final done = delimitations.where((d) => d.status == StepStatus.done).length;
      modal.result = '$done delimitation(s)${failed > 0 ? ' · $failed failed' : ''}';
      notifyListeners();
    });
  }

  static String _runTime(TreeRun r) => r.started.toIso8601String().substring(11, 16);

  Future<void> _delimitRun(ResolvedTool tool, TreeRun r, Set<PtpModel> models,
      {required bool autoMinbr, required bool cropOutgroup, required double pvalue, required CancelToken cancel}) async {
    final out = p.join(r.directory, 'delimitation');
    Directory(out).createSync(recursive: true);
    final tree = parseNewick(await File(mlTreeOf(r)!).readAsString());
    final tips = {for (final t in tree.tips) t.name};
    final input = p.join(out, 'ml_tree_for_mptp.nwk');
    await File(input).writeAsString('${mptpNewick(tree)}\n');
    final outgroup = [for (final o in concat.outgroups) if (tips.contains(o)) o];
    // Cropping needs tips left on both sides.
    final crop = cropOutgroup && outgroup.isNotEmpty && outgroup.length < tips.length - 2;
    void mlog(String m) => log(m, source: 'mPTP');

    // Minimum branch length from the p-distances of the supermatrix.
    double? minbr;
    String? minbrCommand;
    final fasta = matrixFiles['fasta'];
    if (autoMinbr && fasta != null && File(fasta).existsSync() && !(matrix?.isProtein ?? false)) {
      final prefix = p.join(out, 'minbr');
      final args = mptpMinbrArgs(tree: input, fasta: fasta, prefix: prefix, outgroup: outgroup);
      minbrCommand = '${tool.description} ${args.join(' ')}';
      final res = await runTool(tool, args, workingDirectory: out, cancel: cancel, onLog: mlog);
      final file = File('$prefix.txt');
      if (res.ok && file.existsSync()) minbr = parseMinbr(await file.readAsString());
      if (minbr == null) log('mPTP could not detect the minimum branch length; its default (0.0001) is used', level: LogLevel.warn, source: 'mPTP');
    }

    for (final model in PtpModel.values.where(models.contains)) {
      if (cancel.isCancelled) return;
      final method = DelimMethod.ofModel(model);
      final item = '${model.label} · ${r.engine} ${_runTime(r)}';
      _item(item, state: ItemState.running, progress: 0.2);
      delimitations.removeWhere((d) => identical(d.run, r) && d.method == method);
      final d = DelimitationRun(method: method, source: r.engine, run: r, tree: input, directory: out)
        ..minbrAuto = minbr
        ..minbrCommand = minbrCommand
        ..outgroup = outgroup
        ..cropped = crop;
      delimitations.insert(0, d);
      notifyListeners();
      final prefix = p.join(out, model == PtpModel.single ? 'ptp' : 'mptp');
      final args = mptpArgs(tree: input, prefix: prefix, model: model, outgroup: outgroup, cropOutgroup: crop, minbr: minbr, pvalue: pvalue);
      d.command = '${tool.description} ${args.join(' ')}';
      final res = await runTool(tool, args, workingDirectory: out, cancel: cancel, onLog: mlog);
      final txt = File('$prefix.txt');
      if (!res.ok || !txt.existsSync()) {
        d
          ..status = StepStatus.error
          ..error = res.cancelled ? 'cancelled' : 'mPTP ${res.describe}';
        _item(item, state: res.cancelled ? ItemState.skipped : ItemState.failed, note: d.error);
        notifyListeners();
        continue;
      }
      final result = parseMptp(await txt.readAsString(), model, alpha: pvalue)..minbr = minbr;
      final tsv = '$prefix.species.tsv';
      await File(tsv).writeAsString('${delimitationTsv(result)}\n');
      d
        ..result = result
        ..status = StepStatus.done
        ..outputs['txt'] = txt.path
        ..outputs['tsv'] = tsv;
      if (File('$prefix.svg').existsSync()) d.outputs['svg'] = '$prefix.svg';
      _item(item, state: ItemState.done, progress: 1, note: '${result.species.length} species');
      log('${d.label}: ${result.species.length} species (${result.singletons} with one tip)'
          '${result.lrtPassed == false ? ' — LRT failed, the null model (one species) is preferred' : ''}', level: LogLevel.success, source: 'mPTP');
      notifyListeners();
    }
  }

  /// GMYC on the ML tree of [r]: LSD2 dates the tree (relative dates, so it
  /// becomes ultrametric), rooted on the outgroups and without them when
  /// there are outgroups; then the pGMYC command line delimits the species.
  Future<void> _gmycRun(TreeRun r, double pvalue, {required bool cropOutgroup, required CancelToken cancel}) async {
    final lsd2 = await _require(ExternalTool.lsd2);
    final gmyc = await _require(ExternalTool.gmyc);
    final item = 'GMYC · ${r.engine} ${_runTime(r)}';
    _item(item, state: ItemState.running, progress: 0.1);
    final out = p.join(r.directory, 'delimitation', 'gmyc');
    Directory(out).createSync(recursive: true);
    final tree = parseNewick(await File(mlTreeOf(r)!).readAsString());
    final tips = {for (final t in tree.tips) t.name};
    final input = p.join(out, 'ml_tree.nwk');
    await File(input).writeAsString('${mptpNewick(tree)}\n');
    final outgroup = [for (final o in concat.outgroups) if (tips.contains(o)) o];
    final useOutgroup = cropOutgroup && outgroup.isNotEmpty && outgroup.length < tips.length - 2;
    String? ogFile;
    if (useOutgroup) {
      ogFile = p.join(out, 'outgroup.txt');
      await File(ogFile).writeAsString(lsd2OutgroupFile(outgroup));
    }
    delimitations.removeWhere((d) => identical(d.run, r) && d.method == DelimMethod.gmyc);
    final d = DelimitationRun(method: DelimMethod.gmyc, source: r.engine, run: r, tree: input, directory: out)
      ..outgroup = outgroup
      ..cropped = useOutgroup;
    delimitations.insert(0, d);
    notifyListeners();
    void fail(String why) {
      d
        ..status = StepStatus.error
        ..error = why;
      _item(item, state: cancel.isCancelled ? ItemState.skipped : ItemState.failed, note: why);
      notifyListeners();
    }

    // 1. LSD2: ultrametric tree with relative dates.
    final prefix = p.join(out, 'lsd2');
    final lsdArgs = lsd2Args(tree: input, prefix: prefix, seqLength: matrix?.length ?? 1000, outgroupFile: ogFile);
    d.minbrCommand = '${lsd2.description} ${lsdArgs.join(' ')}';
    final lsdRes = await runTool(lsd2, lsdArgs, workingDirectory: out, cancel: cancel, onLog: (m) => log(m, source: 'LSD2'));
    final nexus = File('$prefix.date.nexus');
    final um = nexus.existsSync() ? ultrametricFromNexus(await nexus.readAsString()) : null;
    if (!lsdRes.ok || um == null) return fail(lsdRes.cancelled ? 'cancelled' : 'LSD2 ${lsdRes.ok ? 'wrote no dated tree' : lsdRes.describe}');
    final umFile = p.join(out, 'ultrametric_tree.nwk');
    await File(umFile).writeAsString('$um\n');
    d.outputs['ultrametric'] = umFile;
    _item(item, progress: 0.4);

    // 2. GMYC (single threshold).
    final json = p.join(out, 'gmyc.json');
    final args = gmycArgs(tree: umFile, out: json, pvalue: pvalue);
    d.command = '${gmyc.description} ${args.join(' ')}';
    final res = await runTool(gmyc, args, workingDirectory: out, cancel: cancel, captureStdout: true, onLog: (m) => log(m, source: 'GMYC'));
    if (!res.ok || !File(json).existsSync()) return fail(res.cancelled ? 'cancelled' : 'GMYC ${res.describe}');
    final result = parseGmyc(await File(json).readAsString(), alpha: pvalue);
    final tsv = p.join(out, 'gmyc.species.tsv');
    await File(tsv).writeAsString('${delimitationTsv(result)}\n');
    d
      ..result = result
      ..status = StepStatus.done
      ..outputs['json'] = json
      ..outputs['tsv'] = tsv;
    _item(item, state: ItemState.done, progress: 1, note: '${result.species.length} species');
    log('${d.label}: ${result.species.length} species${result.lrtPassed == false ? ' — LRT failed, one species preferred' : ''}', level: LogLevel.success, source: 'GMYC');
    notifyListeners();
  }

  /// ASAP or ABGD on one alignment: the sequences are written with short
  /// ids, the program runs in its own folder of `07_delimitation`, and its
  /// SPART partitions are read back with the sample names.
  Future<void> _delimitAlignment(String name, DelimMethod method, DelimDistance distance, CancelToken cancel,
      {required double pMin, required double pMax, required int steps, required double slope}) async {
    final tool = await _require(method == DelimMethod.asap ? ExternalTool.asap : ExternalTool.abgd);
    final item = '${method.label} · $name';
    _item(item, state: ItemState.running, progress: 0.2);
    final stamp = DateTime.now().toIso8601String().replaceAll(RegExp(r'[:.]'), '-').substring(0, 19);
    final out = dir(p.join('07_delimitation', '${method.name}_${sanitizeTaxonName(name)}_$stamp'));
    final (fasta, names) = delimitationFasta(await _alignmentForDelimitation(name));
    await File(p.join(out, 'input.fasta')).writeAsString(fasta);
    await File(p.join(out, 'input_names.tsv')).writeAsString('${['id\tname', for (final e in names.entries) '${e.key}\t${e.value}'].join('\n')}\n');
    final abgdDistance = distance.abgd ? distance : DelimDistance.jc69;
    final args = method == DelimMethod.asap
        ? asapArgs(fasta: 'input.fasta', outDir: 'asap', distance: distance)
        : abgdArgs(fasta: 'input.fasta', outDir: 'abgd', distance: abgdDistance, pMin: pMin, pMax: pMax, steps: steps, slope: slope);
    delimitations.removeWhere((d) => d.run == null && d.method == method && d.source == name);
    final d = DelimitationRun(method: method, source: name, directory: out)
      ..distance = method == DelimMethod.abgd ? abgdDistance : distance
      ..command = '${tool.description} ${args.join(' ')}';
    delimitations.insert(0, d);
    notifyListeners();
    final output = StringBuffer();
    final res = await runTool(tool, args, workingDirectory: out, cancel: cancel, onLog: (m) {
      output.writeln(m);
      log(m, source: method.label);
    });
    final sub = Directory(p.join(out, method == DelimMethod.asap ? 'asap' : 'abgd'));
    final sparts = sub.existsSync()
        ? (sub.listSync().whereType<File>().where((f) => f.path.endsWith('.spart')).toList()..sort((a, b) => a.path.compareTo(b.path)))
        : <File>[];
    if (res.cancelled) {
      d
        ..status = StepStatus.error
        ..error = 'cancelled';
      _item(item, state: ItemState.skipped);
      notifyListeners();
      return;
    }
    if (method == DelimMethod.asap) {
      if (sparts.isNotEmpty) {
        d.partitions.addAll(parseSpart(await sparts.first.readAsString(), names));
        d.outputs['spart'] = sparts.first.path;
      }
    } else {
      final priors = abgdPriors(pMin, pMax, steps);
      // ABGD writes "<name>.spart" (initial partitions) and
      // "<name>.rec.spart" (recursive partitions).
      for (final f in sparts) {
        final recursive = f.path.endsWith('.rec.spart');
        final parts = parseSpart(await f.readAsString(), names);
        for (final (i, sp) in parts.indexed) {
          d.partitions.add(sp);
          d.partitionNotes.add('${recursive ? 'rec' : 'init'}|${i < priors.length ? priors[i] : ''}');
        }
        d.outputs[recursive ? 'spartRec' : 'spart'] = f.path;
      }
      // Only one group at every prior: ABGD writes no SPART.
      if (d.partitions.isEmpty && res.ok && output.toString().contains('Only one partition')) {
        d.partitions.add(Spartition('one group', [names.values.toList()]));
        d.partitionNotes.add('init|');
      }
      // Shown first: the initial partition at the prior nearest 0.01.
      var best = 0;
      var gap = double.infinity;
      for (final (i, n) in d.partitionNotes.indexed) {
        final parts = n.split('|');
        final prior = double.tryParse(parts.last);
        if (parts.first != 'init' || prior == null) continue;
        final g = (math.log(prior) - math.log(0.01)).abs();
        if (g < gap) {
          gap = g;
          best = i;
        }
      }
      d.chosen = best;
    }
    if (!res.ok || d.partitions.isEmpty) {
      d
        ..status = StepStatus.error
        ..error = res.ok ? '${method.label} wrote no partition' : '${method.label} ${res.describe}';
      _item(item, state: ItemState.failed, note: d.error);
      notifyListeners();
      return;
    }
    _applyPartition(d);
    d.status = StepStatus.done;
    _item(item, state: ItemState.done, progress: 1, note: '${d.result!.species.length} species');
    log('${d.label}: ${d.partitions.length} partition(s); shown: ${d.result!.species.length} species', level: LogLevel.success, source: method.label);
    notifyListeners();
  }

  void _applyPartition(DelimitationRun d) {
    final sp = d.partitions[d.chosen.clamp(0, d.partitions.length - 1)];
    d.result = Delimitation(species: sp.species);
    final tsv = p.join(d.directory, '${d.method.name}.species.tsv');
    File(tsv).writeAsStringSync('${delimitationTsv(d.result!)}\n');
    d.outputs['tsv'] = tsv;
  }

  /// Shows partition [index] of an ASAP or ABGD run (on the tree, in the
  /// methods text and in its tip → species table).
  void choosePartition(DelimitationRun d, int index) {
    if (index < 0 || index >= d.partitions.length) return;
    d.chosen = index;
    _applyPartition(d);
    notifyListeners();
  }

  /// Tree shown in step 9 (path), or null for the newest.
  String? treeToView;

  void viewTree(String? path) {
    treeToView = path;
    notifyListeners();
  }

  /// Figure of each tree in the viewer (by tree file): its options and
  /// edits, as the viewer reports them. Kept while the viewer is closed
  /// and saved in the project file.
  final Map<String, Map<String, dynamic>> treeFigures = {};

  /// Keeps the figure of [tree] (no rebuild: the viewer already shows it).
  void setTreeFigure(String tree, Map<String, dynamic> figure) {
    treeFigures[tree] = figure;
    if (_watchUnsaved) _scheduleUnsavedCheck();
  }

  void _collectOutputs(TreeRun run, Map<String, String> candidates) {
    for (final e in candidates.entries) {
      final f = p.join(run.directory, e.value);
      if (File(f).existsSync()) run.outputs[e.key] = f;
    }
    run.treeFile ??= run.outputs['tree'];
  }

  /// "Save Results as ZIP": raw/, aligned/, trimmed/, saturation/,
  /// concatenated/ and trees/ from the project folder.
  /// The results were saved as a ZIP (step navigator).
  bool zipSaved = false;

  Future<void> saveZip(String target) => _task('Saving ZIP', (_) async {
        final root = workDir;
        const map = {
          '01_loci': 'raw',
          '02_aligned': 'aligned',
          '03_trimmed': 'trimmed',
          '04_saturation': 'saturation',
          '05_matrix': 'concatenated',
          '06_trees': 'trees',
          '07_delimitation': 'delimitation',
        };
        await _zipInIsolate(root, map, target);
        zipSaved = true;
        log('Saved $target', level: LogLevel.success);
      });

  Future<void> openFolder(String path) async {
    final cmd = Platform.isWindows ? 'explorer' : (Platform.isMacOS ? 'open' : 'xdg-open');
    await Process.run(cmd, [path]);
  }

  @override
  void dispose() {
    _modalTimer?.cancel();
    _notifyTimer?.cancel();
    _organizeDebounce?.cancel();
    _unsavedTimer?.cancel();
    _ncbi?.close();
    logTick.dispose();
    unsaved.dispose();
    super.dispose();
  }
}
