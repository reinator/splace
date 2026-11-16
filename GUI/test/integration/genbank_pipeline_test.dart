@Tags(['integration', 'network'])
@Timeout(Duration(minutes: 15))
library;

import 'dart:io';
import 'dart:math' as math;

import 'package:flutter_test/flutter_test.dart';
import 'package:path/path.dart' as p;
import 'package:splace/src/core/genbank_markers.dart';
import 'package:splace/src/core/genetic_code.dart';
import 'package:splace/src/core/tool_commands.dart';
import 'package:splace/src/services/tool_locator.dart';
import 'package:splace/src/state/project_controller.dart';

/// The HTML app's GenBank workflow end to end, with live NCBI/GBIF and the
/// real tools: accessions → taxonomy → CDS → default markers → header →
/// codon-aware MAFFT → trimAl -backtrans → saturation → NEXUS → IQ-TREE
/// (+ per-gene trees) → ZIP.
///
///   SPLACE_NETWORK=1 flutter test test/integration/genbank_pipeline_test.dart
void main() {
  final enabled = Platform.environment['SPLACE_NETWORK'] == '1';

  test('GenBank accessions → codon-aware alignment → tree', () async {
    final tools = await ToolLocator().resolveAll();
    // Codon-aware trimming uses ClipKIT.
    for (final t in [ExternalTool.mafft, ExternalTool.clipkit, ExternalTool.iqtree]) {
      if (tools[t] == null) {
        markTestSkipped('${t.label} not available');
        return;
      }
    }
    final tmp = Directory.systemTemp.createTempSync('splace_gb_');
    addTearDown(() => tmp.deleteSync(recursive: true));
    final c = ProjectController(autoDetectTools: false)..workDir = p.join(tmp.path, 'project');
    addTearDown(c.dispose);
    c.setThreads(math.min(4, Platform.numberOfProcessors));
    await c.detectTools();
    String why() => c.logs.reversed.take(25).map((e) => '[${e.source}] ${e.message}').toList().reversed.join('\n');

    // Five RefSeq Bufonidae mitogenomes found by Option 1 (NCBI UIDs).
    await c.fetchAccessions('2216740843 2301278816 1915405715 1834962159 1834962145');
    expect(c.records, hasLength(5), reason: why());
    expect(c.organelle.name, 'mt');
    expect(c.geneticCode, 2, reason: 'vertebrate mitochondrial code for mitogenomes');

    await c.fetchTaxonomy();
    expect(c.taxonomyResolved, isTrue);
    expect(c.records.where((r) => r.taxonomy.family == 'Bufonidae'), isNotEmpty, reason: why());

    c.toggleFeatureType('CDS', true);
    c.selectDefaultMarkers();
    expect(c.activeMarkers, hasLength(13));
    // Keep the test quick: three genes.
    for (final m in c.activeMarkers.where((m) => !{'COI', 'CYTB', 'ND2'}.contains(m)).toList()) {
      c.toggleMarker(m, false);
    }
    c.setHeader(const HeaderTemplate(fields: ['species', 'accession']));
    c.confirmHeader();

    c.setSeqMode(SeqMode.codon);
    await c.alignAll();
    for (final l in c.loci.values) {
      expect(l.alignStatus, StepStatus.done, reason: '${l.name}: ${l.message}\n${why()}');
      expect(codonIntegrity(l.aligned!).ok, isTrue);
      expect(l.proteinAligned, isNotNull);
      expect(l.aligned!.first.header, matches(RegExp(r'^[A-Z][a-z]+_[a-z]+.*_NC_\d+\.\d$')));
    }
    expect(c.trimPrompt, isTrue);

    await c.trimAll();
    for (final l in c.loci.values) {
      expect(l.trimStatus, StepStatus.done, reason: '${l.name}: ${l.message}\n${why()}');
      expect(codonIntegrity(l.trimmed!).ok, isTrue, reason: l.trimMethod);
      expect(l.trimMethod, contains('ClipKIT'), reason: why());
    }

    await c.runSaturation();
    expect(c.pooledSaturation, isNotNull, reason: why());

    c.setConcat(c.concat.copyWith(outgroups: [c.loci.values.first.trimmed!.first.header]));
    await c.buildMatrix();
    final m = c.matrix!;
    expect(m.taxa, hasLength(5));
    expect(File(c.matrixFiles['nexus']!).readAsStringSync(), contains('OUTGROUP'));
    expect(File(c.matrixFiles['partition']!).readAsStringSync(), contains('DNA, COI_pos1'));

    await c.runIqtree(IqtreeOptions(model: 'GTR+F+G4', ufboot: 1000, alrt: 0, threads: '2', seed: 1, outgroup: c.concat.outgroups), perGene: true);
    final run = c.treeRuns.first;
    expect(run.status, StepStatus.done, reason: why());
    expect(run.geneTrees, hasLength(3));
    expect(run.geneTrees.values.every((ok) => ok), isTrue, reason: why());

    final zip = p.join(tmp.path, 'results.zip');
    await c.saveZip(zip);
    expect(File(zip).lengthSync(), greaterThan(1000));
  }, skip: !enabled);
}
