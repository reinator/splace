import 'package:flutter_test/flutter_test.dart';
import 'package:splace/src/core/fasta.dart';
import 'package:splace/src/core/locus_organizer.dart';

FastaInput _file(String name, Map<String, String> records) =>
    FastaInput(name, [for (final e in records.entries) SeqRecord(e.key, e.value)]);

void main() {
  test('tokenize drops empty tokens and supports several delimiters', () {
    expect(tokenize('a__b|c d', '_| '), ['a', 'b', 'c', 'd']);
    expect(tokenize('abc', ''), ['abc']);
  });

  test('header tags, bracketed and plain', () {
    final tags = parseHeaderTags('lcl|NC_1 [gene=COX1] [protein=cytochrome c oxidase subunit I] sample=S12;locus=COI');
    expect(tags['gene'], 'COX1');
    expect(tags['protein'], 'cytochrome c oxidase subunit I');
    expect(tags['sample'], 'S12');
    expect(tags['locus'], 'COI');
  });

  test('field rule slicing with negative indices', () {
    final r = SeqRecord('Genus_species_voucher12_COI', 'ACGT');
    expect(const FieldRule.header(start: -1, delimiters: '_').extract('x.fasta', r), 'COI');
    expect(const FieldRule.header(start: 0, end: -1, delimiters: '_').extract('x.fasta', r), 'Genus_species_voucher12');
    expect(const FieldRule.fileName(start: 2).extract('NC_030596_Oryza_brachyantha.fasta', r), 'Oryza_brachyantha');
    expect(const FieldRule.header(start: 9, delimiters: '_').extract('x.fasta', r), isNull);
    expect(const FieldRule.regex(r'voucher(\d+)').extract('x.fasta', r), '12');
  });

  group('sample-per-file data (one specimen per file, header = locus)', () {
    final inputs = [
      _file('NC_030596_Oryza_brachyantha.fasta', {'psbA': 'AAAA', 'matK': 'CCCC', 'rbcL': 'GGGG'}),
      _file('NC_015401_Olea_europaea.fasta', {'psbA': 'AAAT', 'matK': 'CCCT', 'rbcL': 'GGGTT', 'rbcL ': 'GG'}),
      _file('NC_007898_Solanum_lycopersicum.fasta', {'psbA': 'AATT', 'MATK': 'CCTT'}),
    ];

    test('layout is detected', () {
      final s = suggestLayout(inputs);
      expect(s.kind, LayoutKind.samplePerFile);
      expect(s.config.sample.source, FieldSource.fileName);
      expect(s.config.locus.source, FieldSource.header);
    });

    test('reorganized into locus files with sample headers', () {
      final result = organizeLoci(inputs, suggestLayout(inputs).config);
      expect(result.loci, ['matK', 'psbA', 'rbcL'],
          reason: 'MATK merged with matK; sorted by coverage, then name');
      expect(result.samples, hasLength(3));
      expect(result.samplesFor('rbcL'), 2);
      expect(result.duplicates, hasLength(1));
      // keepLongest keeps the 5 bp rbcL copy.
      expect(result.byLocus['rbcL']!['NC_015401_Olea_europaea']!.record.sequence, 'GGGTT');

      final fastas = result.locusFastas();
      final psbA = parseFasta(fastas['psbA']!);
      expect(psbA.map((r) => r.header), containsAll(['NC_030596_Oryza_brachyantha', 'NC_007898_Solanum_lycopersicum']));
      expect(result.locusFastas(minSamples: 3).keys, ['matK', 'psbA']);
    });

    test('mapping table records dropped duplicates', () {
      final tsv = organizeLoci(inputs, suggestLayout(inputs).config).mappingTsv();
      expect(tsv, contains('dropped_duplicate'));
    });

    test('occupancy', () {
      final result = organizeLoci(inputs, suggestLayout(inputs).config);
      expect(result.occupancy, closeTo(8 / 9, 1e-9));
    });
  });

  test('locus-per-file data is recognized', () {
    final inputs = [
      _file('COI.fasta', {'Bufo_bufo': 'AC', 'Rhinella_marina': 'AG'}),
      _file('CYTB.fasta', {'Bufo_bufo': 'AC', 'Rhinella_marina': 'AG'}),
    ];
    expect(suggestLayout(inputs).kind, LayoutKind.locusPerFile);
  });

  test('mixed headers in a single file, with synonym merging', () {
    final inputs = [
      _file('all.fasta', {
        'Bufo_bufo_COX1': 'AAAA',
        'Bufo_bufo_CYTB': 'CCCC',
        'Rhinella_marina_COI': 'AAAT',
        'Rhinella_marina_cob': 'CCCT',
      }),
    ];
    final s = suggestLayout(inputs);
    expect(s.kind, LayoutKind.mixedHeader);
    final result = organizeLoci(inputs, s.config);
    expect(result.loci.toSet(), {'COI', 'CYTB'});
    expect(result.samples, ['Bufo_bufo', 'Rhinella_marina']);
    expect(result.locusAliases['COI'], containsAll(['COX1', 'COI']));
  });

  test('NCBI tagged headers', () {
    final inputs = [
      _file('seqs.fasta', {
        'lcl|MN1_cds_1 [gene=COX1] [organism=Bufo bufo]': 'AAAA',
        'lcl|MN2_cds_1 [gene=CYTB] [organism=Bufo bufo]': 'CCCC',
      }),
    ];
    final s = suggestLayout(inputs);
    expect(s.kind, LayoutKind.taggedHeader);
    final result = organizeLoci(inputs, s.config);
    expect(result.samples, ['Bufo_bufo']);
  });

  test('unassigned records are reported, not silently dropped', () {
    final inputs = [_file('a.fasta', {'onlyonetoken': 'AC'})];
    final result = organizeLoci(
      inputs,
      const OrganizerConfig(sample: FieldRule.header(start: 0, delimiters: '_'), locus: FieldRule.header(start: 1, delimiters: '_')),
    );
    expect(result.unassigned, hasLength(1));
    expect(result.unassigned.single.missingLocus, isTrue);
  });

  group('many markers in per-sample files (target capture)', () {
    // `>ITV01908-L12`: sample and marker in each header, one file per
    // sample; one file writes its sample without the zero (ITV1054).
    List<FastaInput> capture() => [
          for (final (file, id) in [('ITV01908.fasta', 'ITV01908'), ('ITV01054.fasta', 'ITV1054'), ('ERR2040879.fasta', 'ERR2040879')])
            FastaInput(file, [
              for (var l = 1; l <= 30; l++)
                if (!(id == 'ERR2040879' && l > 25)) SeqRecord('$id-L$l', 'ACGT' * (5 + l)),
            ]),
        ];

    test('the layout guess finds sample and marker by how they behave', () {
      final s = suggestLayout(capture());
      expect(s.kind, LayoutKind.mixedHeader);
      final r = organizeLoci(capture(), s.config);
      expect(r.samples, ['ERR2040879', 'ITV1054', 'ITV01908'], reason: 'numbers in natural order');
      expect(r.loci, hasLength(30));
      expect(r.samplesFor('L1'), 3);
      expect(r.samplesFor('L30'), 2);
      expect(r.sampleFiles['ITV1054'], {'ITV01054.fasta'});
    });

    test('name variations and renaming', () {
      final cfg = suggestLayout(capture()).config;
      // Leading zeros: L01 = L1 and ITV1054 = ITV01054.
      final inputs = capture()
        ..add(FastaInput('extra.fasta', [SeqRecord('ITV01054-L01', 'ACGTACGTAC'), SeqRecord('itv01908-L002', 'ACGT')]));
      final zeros = organizeLoci(inputs, cfg.copyWith(ignoreLeadingZeros: true));
      expect(zeros.loci, hasLength(30), reason: 'L01 joins L1, L002 joins L2');
      expect(zeros.samples, hasLength(3), reason: 'ITV01054 = ITV1054, itv01908 = ITV01908');
      expect(zeros.sampleSpellings.values.expand((v) => v), containsAll(['ITV1054', 'ITV01054']));
      expect(zeros.duplicates, isNotEmpty);

      final renamed = organizeLoci(capture(), cfg.copyWith(sampleAliases: {'ITV1054': 'ITV01054'}));
      expect(renamed.samples, contains('ITV01054'));
      expect(renamed.samples, isNot(contains('ITV1054')));
    });

    test('a split with one marker per sequence is not plausible', () {
      final bad = organizeLoci(capture(), const OrganizerConfig(sample: FieldRule.fileName(delimiters: ''), locus: FieldRule.header()));
      expect(bad.loci.length, bad.records);
      expect(bad.records, 85);
      expect(stripLeadingZeros('ITV01054-L007 0.5'), 'ITV1054-L7 0.5');
    });

    test('name structures are standardized by adding or removing zeros', () {
      expect(nameStructure('ITV01054'), 'ITV#####');
      expect(nameStructure('itv1054'), 'ITV####');
      expect(conformName('ITV1054', 'ITV#####'), 'ITV01054');
      expect(conformName('ITV001054', 'ITV#####'), 'ITV01054');
      expect(conformName('ITV123456', 'ITV#####'), isNull, reason: 'too many digits');
      expect(conformName('MPEG1054', 'ITV#####'), isNull, reason: 'other letters');
      expect(conformName('ITV1054A', 'ITV#####'), isNull);
      expect(headerShape('ITV01054-L12'), 'x-x');
      expect(headerShape('L12|ITV01054|2020'), 'x|x|x');
    });

    test('each header format can have its own parts or be left out', () {
      final inputs = [
        FastaInput('a.fasta', [SeqRecord('S1-L1', 'ACGT'), SeqRecord('S1-L2', 'ACGT')]),
        FastaInput('b.fasta', [SeqRecord('L1|S2', 'ACGT'), SeqRecord('L2|S2', 'ACGT')]),
        FastaInput('c.fasta', [SeqRecord('junk header', 'ACGT')]),
      ];
      const main = OrganizerConfig(sample: FieldRule.header(start: 0, end: 1, delimiters: '-|'), locus: FieldRule.header(start: 1, end: 2, delimiters: '-|'));
      final plain = organizeLoci(inputs, main);
      expect(plain.samples, containsAll(['S1', 'L1', 'L2']), reason: 'the second format is read backwards');

      final fixed = organizeLoci(
        inputs,
        main.copyWith(formats: {
          'x|x': const FormatRule(sample: FieldRule.header(start: 1, end: 2, delimiters: '-|'), locus: FieldRule.header(start: 0, end: 1, delimiters: '-|')),
          'x x': const FormatRule(sample: FieldRule.header(), locus: FieldRule.header(), skip: true),
        }),
      );
      expect(fixed.samples, ['S1', 'S2']);
      expect(fixed.loci, ['L1', 'L2']);
      expect(fixed.samplesFor('L1'), 2);
      expect(fixed.skipped, 1);
      expect(fixed.unassigned, isEmpty);
    });
  });
}
