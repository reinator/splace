import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/core/fasta.dart';
import 'package:splace/src/core/supermatrix.dart';
import 'package:splace/src/core/tool_commands.dart';

void main() {
  final coi = LocusAlignment('COI', [
    SeqRecord('A', 'ATGAAACCC'),
    SeqRecord('B', 'ATGAAGCCT'),
  ]);
  final rrna = LocusAlignment('16S', [
    SeqRecord('A', 'GGGG'),
    SeqRecord('C', 'GGGA'),
  ], coding: false);

  test('missing taxa are filled with ?', () {
    final m = buildSupermatrix([coi, rrna]);
    expect(m.taxa, ['A', 'B', 'C']);
    expect(m.sequences['B'], 'ATGAAGCCT????');
    expect(m.sequences['C'], '?????????GGGA');
    expect(m.length, 13);
  });

  test('codon partitions use strided ranges', () {
    final m = buildSupermatrix([coi, rrna]);
    expect(m.partitions.map((p) => '${p.name}=${p.rangeString}'),
        ['COI_pos1=1-7\\3', 'COI_pos2=2-8\\3', 'COI_pos3=3-9\\3', '16S=10-13']);
    expect(m.toRaxmlPartitions(), contains('GTR+G, COI_pos1 = 1-7\\3'));
    expect(m.toIqtreePartitions(), contains('charset COI_pos3 = 3-9\\3;'));
  });

  test('excluding 3rd positions keeps 1st/2nd with stride 2', () {
    final m = buildSupermatrix([coi, rrna], const ConcatOptions(excludePositions: {3}));
    expect(m.sequences['A'], 'ATAACC' 'GGGG');
    expect(m.partitions.map((p) => '${p.name}=${p.rangeString}'),
        ['COI_pos1=1-5\\2', 'COI_pos2=2-6\\2', '16S=7-10']);
  });

  test('frame offset is respected when excluding positions', () {
    final shifted = LocusAlignment('ND1', [SeqRecord('A', 'CATGAAA'), SeqRecord('B', 'CATGAAG')], frame: 1);
    final m = buildSupermatrix([shifted], const ConcatOptions(excludePositions: {3}));
    // Column 0 is a 3rd position (tail of a previous codon) and goes away.
    expect(m.sequences['A'], 'ATAA');
  });

  test('complete-loci-only policy drops incomplete loci', () {
    final m = buildSupermatrix([coi, rrna], const ConcatOptions(missingPolicy: MissingDataPolicy.completeLociOnly));
    expect(m.droppedLoci, ['COI', '16S']);
  });

  test('global codon-position scheme merges loci', () {
    final nd1 = LocusAlignment('ND1', [SeqRecord('A', 'ATGCCC'), SeqRecord('B', 'ATGCCA')]);
    final m = buildSupermatrix([coi, nd1], const ConcatOptions(scheme: PartitionScheme.byCodonPosition));
    expect(m.partitions.first.name, 'codon_pos1');
    expect(m.partitions.first.rangeString, '1-7\\3, 10-13\\3');
  });

  test('guessCoding', () {
    expect(guessCoding('16S', rrna.records, 1).coding, isFalse);
    expect(guessCoding('COI', coi.records, 1), (coding: true, frame: 0));
  });

  group('RAxML-NG interoperability', () {
    test('IQ-TREE model names are translated', () {
      expect(iqtreeModelToRaxml('TIM2+F+I+G4'), 'TIM2uf+F+I+G4');
      expect(iqtreeModelToRaxml('TIM2e+G4'), 'TIM2+G4');
      expect(iqtreeModelToRaxml('K2P+R3'), 'K80+R3');
      expect(iqtreeModelToRaxml('TNe'), 'TN93ef');
      expect(iqtreeModelToRaxml('HKY+F+G'), 'HKY+F+G4');
      expect(iqtreeModelToRaxml('LG+G4'), isNull);
    });

    test('best_scheme conversion', () {
      const scheme = 'TIM2+F+G4, COI_pos1 = 1-7\\3\nHKY+F+I, COI_pos2_COI_pos3 = 2-8\\3, 3-9\\3\nWEIRD, x = 10-13\n';
      final parsed = parseRaxmlStyleScheme(scheme);
      expect(parsed, hasLength(3));
      final out = schemeToRaxml(parsed);
      expect(out.partitions, contains('TIM2uf+F+G4, COI_pos1 = 1-7\\3'));
      expect(out.partitions, contains('HKY+F+I, COI_pos2_COI_pos3 = 2-8\\3, 3-9\\3'));
      expect(out.replaced.single, startsWith('x:'));
    });

    test('IQ-TREE 3 best_scheme.nex / best_model.nex parsing', () {
      // Verbatim output of IQ-TREE 3.1.4 `-m MF+MERGE --mset raxml`.
      const schemeNex = '#nexus\nbegin sets;\n'
          '  charset COI_pos1_COI_pos2_ND1_pos1_ND1_pos2 = 1-298\\3  2-299\\3  301-598\\3  302-599\\3;\n'
          '  charset COI_pos3_ND1_pos3 = 3-300\\3  303-600\\3;\n'
          '  charpartition mymodels =\n'
          '    GTR: COI_pos1_COI_pos2_ND1_pos1_ND1_pos2,\n'
          '    GTR+ASC: COI_pos3_ND1_pos3;\nend;\n';
      final entries = parseIqtreeNexusScheme(schemeNex);
      expect(entries.map((e) => e.model), ['GTR', 'GTR+ASC']);
      expect(entries.last.ranges, '3-300\\3, 303-600\\3');
      final out = schemeToRaxml(entries);
      expect(out.partitions, 'GTR, COI_pos1_COI_pos2_ND1_pos1_ND1_pos2 = 1-298\\3, 2-299\\3, 301-598\\3, 302-599\\3\n'
          'GTR, COI_pos3_ND1_pos3 = 3-300\\3, 303-600\\3\n');
      expect(out.replaced.single, contains('ascertainment correction dropped'));

      const modelNex = '#nexus\nbegin sets;\n  charset a = 1-298\\3;\n  charpartition mymodels =\n'
          '    GTR{0.87,0.78,1.14,1.17,0.84}+F{0.23,0.27,0.26,0.22}+G4{0.5}: a{0.571004};\nend;\n';
      expect(parseIqtreeNexusScheme(modelNex).single.model, 'GTR+F+G4');
    });

    test('raxml log parsing', () {
      const log = '* Recommended number of threads / MPI processes: 6\n'
          'Bootstrapping converged after 350 replicates.\n'
          'Final LogLikelihood: -12345.6789\n';
      final s = parseRaxmlLog(log);
      expect(s.recommendedThreads, 6);
      expect(s.bootstrapTrees, 350);
      expect(s.finalLogLikelihood, closeTo(-12345.6789, 1e-6));
    });

    test('command lines', () {
      const r = RaxmlOptions(maxThreads: 8, outgroup: ['A']);
      expect(r.allArgs(msa: 'x.rba', prefix: 'out/t'),
          containsAllInOrder(['--all', '--msa', 'x.rba', '--threads', 'auto{8}', '--outgroup', 'A']));
      const fast = RaxmlOptions(fast: true, supportMetrics: 'sh,ebg');
      final moose = fast.allArgs(msa: 'm.phy', model: 'auto=p.txt', prefix: 't', moose: true);
      expect(moose, containsAllInOrder(['--all', '--fast', '--msa', 'm.phy', '--model', 'auto=p.txt', '--moose-options', 'criterion=bic']));
      expect(moose, isNot(contains('--bs-trees')), reason: 'SH-aLRT and EBG need no replicates');
      expect(moose, isNot(contains('--tree')), reason: 'adaptive starting trees by default');
      expect(RaxmlSupport.parse('fbp, SH').map((m) => m.outputKey), ['treeFBP', 'treeSH']);
      expect(parseMooseModels('Best model(s):\nPartition #0: HKY+FO+I+G4m (LogLH = -1782.1  BIC = 3672.7)\nPartition #1: TIM2+FO+G4m (LogLH = -2.0  BIC = 4.0)'),
          ['HKY+FO+I+G4m', 'TIM2+FO+G4m']);
      const search = IqtreeOptions(model: 'GTR+G', ufboot: 1000, alrt: 1000, bnni: true, nstop: 500, perturbation: 0.2, runs: 3, seed: 7);
      expect(search.args(alignment: 'm.phy', prefix: 'p'),
          containsAllInOrder(['-B', '1000', '--alrt', '1000', '--bnni', '-nstop', '500', '-pers', '0.2', '--runs', '3', '--seed', '7']));
      expect(const IqtreeOptions(ufboot: 0, bnni: true).args(alignment: 'm', prefix: 'p'), isNot(contains('--bnni')), reason: '--bnni needs UFBoot');
      // Hundreds of partitions: MERGE uses the relaxed clustering.
      const merge = IqtreeOptions();
      expect(merge.args(alignment: 'm', partitions: 'p.nex', prefix: 'p', partitionCount: 1347), containsAllInOrder(['--rcluster', '10']));
      expect(merge.args(alignment: 'm', partitions: 'p.nex', prefix: 'p', partitionCount: 12), isNot(contains('--rcluster')));
      expect(const IqtreeOptions(model: 'GTR+G').args(alignment: 'm', partitions: 'p.nex', prefix: 'p', partitionCount: 1347), isNot(contains('--rcluster')));
      expect(search.copyWith(threads: '3', maxThreads: 3).args(alignment: 'm', prefix: 'p'), containsAllInOrder(['--runs', '3', '-T', '3', '--threads-max', '3']));
      const iq = IqtreeOptions(model: 'MF+MERGE', raxmlModelsOnly: true);
      final args = iq.args(alignment: 'm.phy', partitions: 'p.nex', prefix: 'mf');
      expect(args, containsAllInOrder(['-m', 'MF+MERGE', '--mset', 'raxml']));
      expect(args, isNot(contains('-B')), reason: 'no bootstrap without tree inference');
    });
  });

  test('overview: bases, gaps and missing loci per taxon', () {
    final a = LocusAlignment('A', [SeqRecord('t1', 'ACG-T'), SeqRecord('t2', 'AC--T')]);
    final b = LocusAlignment('B', [SeqRecord('t1', 'GGG')]);
    final o = matrixOverview(buildSupermatrix([a, b], const ConcatOptions(scheme: PartitionScheme.byLocus)));
    expect(o.loci, ['A', 'B']);
    expect(o.lengths, [5, 3]);
    final t2 = o.taxa.indexOf('t2');
    expect(o.bases[t2], [3, 0], reason: 't2 has no B: missing');
    expect(o.gaps[t2][0], 2);
    expect(o.taxaWith(1), 1);
    expect(o.gapFraction(0), closeTo(3 / 10, 1e-9));
    expect(o.toTsv(), contains('t2\tB\t3\t0\t0\tmissing'));
  });

  test('sanitizeTaxonName', () {
    expect(sanitizeTaxonName("Bufo bufo (voucher: 12)"), 'Bufo_bufo_voucher_12');
  });
}
