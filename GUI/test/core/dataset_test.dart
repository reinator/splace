@Tags(['dataset'])
library;

import 'dart:io';

import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/core/fasta.dart';
import 'package:splace/src/core/locus_organizer.dart';

/// Runs the organizer on the repository dataset (datasets/inputs.tar.gz,
/// one chloroplast genome per file, headers = gene names). Extract it and
/// point SPLACE_DATASET to the `inputs` folder; skipped otherwise.
void main() {
  final dir = Platform.environment['SPLACE_DATASET'];
  test('repository dataset is reorganized into locus files', () {
    final files = Directory(dir!)
        .listSync()
        .whereType<File>()
        .where((f) => f.path.endsWith('.fasta'))
        .toList();
    final inputs = [
      for (final f in files) FastaInput(f.uri.pathSegments.last, parseFasta(f.readAsStringSync())),
    ];
    final s = suggestLayout(inputs);
    expect(s.kind, LayoutKind.samplePerFile);

    final watch = Stopwatch()..start();
    final result = organizeLoci(inputs, s.config);
    // ignore: avoid_print
    print('${inputs.length} files → ${result.loci.length} loci × ${result.samples.length} samples, '
        'occupancy ${(result.occupancy * 100).toStringAsFixed(1)}%, '
        '${result.duplicates.length} duplicated cells, ${result.unassigned.length} unassigned, '
        '${watch.elapsedMilliseconds} ms');
    // ignore: avoid_print
    print('Top loci: ${result.loci.take(8).map((l) => '$l(${result.samplesFor(l)})').join(', ')}');
    expect(result.samples.length, inputs.length, reason: 'one sample per file');
    expect(result.unassigned, isEmpty);
    expect(result.loci, contains('rbcL'));
    final rbcL = parseFasta(result.locusFastas()['rbcL']!);
    expect(rbcL.map((r) => r.header).toSet().length, rbcL.length, reason: 'unique sample names per locus');
  }, skip: dir == null ? 'set SPLACE_DATASET to run' : false);
}
