import 'dart:convert';
import 'dart:io';

import 'package:archive/archive_io.dart';
import 'package:flutter_test/flutter_test.dart';
import 'package:path/path.dart' as p;
import 'package:splace/src/core/delimitation.dart';
import 'package:splace/src/core/fasta.dart';
import 'package:splace/src/core/genbank.dart';
import 'package:splace/src/core/genbank_markers.dart';
import 'package:splace/src/core/saturation.dart';
import 'package:splace/src/core/tool_commands.dart';
import 'package:splace/src/services/species_id.dart';
import 'package:splace/src/state/project_controller.dart';

/// project.json and state.json of a project file (without what changes on
/// every save: the time and the log).
({Map<String, dynamic> params, Map<String, dynamic> state}) _read(String file) {
  final archive = ZipDecoder().decodeBytes(File(file).readAsBytesSync());
  Map<String, dynamic> json(String name) => jsonDecode(utf8.decode(archive.findFile(name)!.content)) as Map<String, dynamic>;
  return (params: json('project.json')..remove('saved'), state: json('state.json')..remove('logs'));
}

void main() {
  late Directory tmp;

  setUp(() => tmp = Directory.systemTemp.createTempSync('splace_spc_'));
  tearDown(() => tmp.deleteSync(recursive: true));

  test('a saved project opens where it stopped, with the same parameters and results', () async {
    final c = ProjectController(autoDetectTools: false)..workDir = p.join(tmp.path, 'a');
    addTearDown(c.dispose);
    await c.addGenBankFiles([for (final n in ['01', '02', '03', '04']) 'assets/examples/genbank/$n.gb']);
    expect(c.records, hasLength(4));
    c.skipTaxonomy();
    c.toggleFeatureType('CDS', true);
    // Edits made in step 2.
    c.records.first
      ..editedOrganism = 'Edited species'
      ..ranks = const TaxonomyRanks(source: 'manual', family: 'Bufonidae', order: 'Anura');
    final markers = c.markers.keys.where((m) => const ['CYTB', 'ND2'].contains(m)).toList();
    expect(markers, hasLength(2));
    for (final m in markers) {
      c.toggleMarker(m, true);
    }
    expect(c.activeMarkers, markers);
    c.setHeader(const HeaderTemplate(fields: ['accession', 'species']));
    c.confirmHeader();
    await c.prepareLoci();
    // Alignments as MAFFT would leave them (padded), and a saturation test.
    for (final l in c.loci.values) {
      final width = l.raw.map((r) => r.sequence.length).reduce((a, b) => a > b ? a : b);
      l
        ..aligned = [for (final r in l.raw) SeqRecord(r.header, r.sequence.padRight(width, '-'))]
        ..alignStatus = StepStatus.done
        ..pipeline.addAll(['MAFFT --auto', 'DNA']);
    }
    final first = c.loci.values.first;
    first.saturation = testSaturation(first.name, first.aligned!, const SaturationOptions(frame: 0));
    c.mafft = c.mafft.copyWith(strategy: '--localpair', maxIterate: 1000);
    c.satModel = 'HKY+G';
    c.concat = c.concat.copyWith(outgroups: [first.aligned!.first.header]);
    await c.buildMatrix();
    expect(c.matrix, isNotNull);
    // A tree run with its files, and a delimitation on its tree.
    final runDir = Directory(p.join(c.workDir, '06_trees', 'iqtree_x'))..createSync(recursive: true);
    final treeFile = File(p.join(runDir.path, 'splace.treefile'))..writeAsStringSync('(a:1,b:1,c:1);\n');
    final run = TreeRun('IQ-TREE', runDir.path)
      ..status = StepStatus.done
      ..command = 'iqtree3 -s concatenated.phy -m MFP+MERGE'
      ..iqtree = const IqtreeOptions(model: 'MFP+MERGE', seed: 7)
      ..bestModel = 'GTR+F+I+G4'
      ..treeFile = treeFile.path;
    run.outputs['tree'] = treeFile.path;
    c.treeRuns.add(run);
    // The figure of that tree, as the viewer reports it.
    c.setTreeFigure(treeFile.path, {
      'version': 1,
      'title': 'Figure 1',
      'font': 10.5,
      'collapsed': [2],
      'renamed': {'a': {'text': 'A. alpha', 'italic': true, 'bold': false}},
    });
    c.delimitations.add(DelimitationRun(method: DelimMethod.mptp, source: 'IQ-TREE', directory: runDir.path, run: run, tree: treeFile.path)
      ..status = StepStatus.done
      ..result = Delimitation(model: PtpModel.multi, species: [
        ['a', 'b'],
        ['c'],
      ]));
    c.speciesIds = {
      c.records.first.accession: RecordIdentification(c.records.first.accession, 'Edited species')
        ..addEvidence(IdSource.genbank, markers.first, 'Other species', 99.4)
        ..addEvidence(IdSource.bold, markers.first, 'Edited species', 98.2),
    };

    final file = p.join(tmp.path, 'study.splc');
    // Never saved, with data: unsaved (the window title says so).
    c.checkUnsaved();
    expect(c.unsaved.value, isTrue);
    expect(await c.saveProject(file, appVersion: 'test'), isTrue);
    expect(c.projectFile, file);
    c.checkUnsaved();
    expect(c.unsaved.value, isFalse);
    // A change makes it unsaved; undoing it makes it saved again.
    c.satModel = 'GTR+G';
    c.checkUnsaved();
    expect(c.unsaved.value, isTrue);
    c.satModel = 'HKY+G';
    c.checkUnsaved();
    expect(c.unsaved.value, isFalse);
    final saved = _read(file);
    expect(saved.params['format'], 'splace-project');
    expect(saved.params['4_alignment']['mafft']['strategy'], '--localpair');
    expect(saved.params['7_inference'].single['command'], contains('iqtree3'));

    final d = ProjectController(autoDetectTools: false);
    addTearDown(d.dispose);
    expect(await d.openProject(file, into: p.join(tmp.path, 'b')), isTrue);
    expect(d.workDir, p.join(tmp.path, 'b'));
    expect(d.projectFile, file);
    d.checkUnsaved();
    expect(d.unsaved.value, isFalse, reason: 'a project just opened is saved');
    expect(d.records, hasLength(4));
    expect(d.records.first.organismName, 'Edited species');
    expect(d.records.first.ranks.family, 'Bufonidae');
    expect(d.records.first.features.length, c.records.first.features.length);
    expect(d.activeMarkers, markers);
    expect(d.headerConfirmed, isTrue);
    expect(d.loci.keys, c.loci.keys);
    for (final l in c.loci.values) {
      expect(d.loci[l.name]!.aligned!.map((r) => r.sequence), l.aligned!.map((r) => r.sequence));
    }
    expect(d.loci[first.name]!.saturation!.subsets[CodonSubset.pos3]!.slope, first.saturation!.subsets[CodonSubset.pos3]!.slope);
    expect(d.mafft.strategy, '--localpair');
    expect(d.satModel, 'HKY+G');
    expect(d.concat.outgroups, c.concat.outgroups);
    expect(d.matrix!.taxa, c.matrix!.taxa);
    expect(File(d.matrixFiles['phylip']!).existsSync(), isTrue);
    expect(p.isWithin(d.workDir, d.matrixFiles['phylip']!), isTrue);
    // The tree and its files moved with the project.
    final r = d.treeRuns.single;
    expect(r.iqtree!.seed, 7);
    expect(r.bestModel, 'GTR+F+I+G4');
    expect(r.started, run.started);
    expect(p.isWithin(d.workDir, r.treeFile!), isTrue);
    expect(File(r.treeFile!).readAsStringSync(), '(a:1,b:1,c:1);\n');
    expect(d.treeFigures[r.treeFile!], c.treeFigures[treeFile.path]);
    expect(identical(d.delimitations.single.run, r), isTrue);
    expect(d.delimitations.single.result!.species, hasLength(2));
    final id = d.speciesIds[c.records.first.accession]!;
    expect(id.candidates.map((x) => x.species), c.speciesIds.values.single.candidates.map((x) => x.species));
    expect(id.disagrees, c.speciesIds.values.single.disagrees);
    // The log of the project came too.
    expect(d.logs.any((e) => e.message.contains('Supermatrix')), isTrue);

    // Saved again, the project is the same (paths are kept relative).
    final again = p.join(tmp.path, 'again.splc');
    expect(await d.saveProject(again, appVersion: 'test'), isTrue);
    final resaved = _read(again);
    expect(resaved.params, saved.params);
    expect(resaved.state, saved.state);
  });

  test('a project file cannot write outside the project folder', () async {
    final file = p.join(tmp.path, 'evil.splc');
    final archive = Archive()
      ..addFile(ArchiveFile.string('project.json', jsonEncode({'format': 'splace-project', 'version': 1})))
      ..addFile(ArchiveFile.string('state.json', '{}'))
      ..addFile(ArchiveFile.string('work/../../escaped.txt', 'x'))
      ..addFile(ArchiveFile.string('work/01_loci/ok.fasta', '>a\nACGT\n'));
    File(file).writeAsBytesSync(ZipEncoder().encode(archive));
    final d = ProjectController(autoDetectTools: false);
    addTearDown(d.dispose);
    final into = p.join(tmp.path, 'x', 'project');
    expect(await d.openProject(file, into: into), isTrue);
    expect(File(p.join(into, '01_loci', 'ok.fasta')).existsSync(), isTrue);
    expect(File(p.join(tmp.path, 'escaped.txt')).existsSync(), isFalse);
    expect(File(p.join(tmp.path, 'x', 'escaped.txt')).existsSync(), isFalse);
  });

  test('a file that is not a SPLACE project is refused', () async {
    final file = p.join(tmp.path, 'other.splc');
    File(file).writeAsBytesSync(ZipEncoder().encode(Archive()..addFile(ArchiveFile.string('readme.txt', 'hi'))));
    final d = ProjectController(autoDetectTools: false);
    addTearDown(d.dispose);
    expect(await d.openProject(file, into: p.join(tmp.path, 'y')), isFalse);
    expect(d.logs.any((e) => e.error && e.message.contains('not a SPLACE project')), isTrue);
  });
}
