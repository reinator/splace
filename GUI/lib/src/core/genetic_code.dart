/// NCBI genetic codes, CDS validation, translation and protein-guided
/// back-translation (codon-aware alignment), ported from the HTML app.
library;

import 'fasta.dart';

const _standard = {
  'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', //
  'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
  'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
  'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
  'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
  'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
  'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
  'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
};

/// Differences from the standard code (NCBI translation tables).
const _changes = <int, Map<String, String>>{
  1: {},
  2: {'ATA': 'M', 'TGA': 'W', 'AGA': '*', 'AGG': '*'},
  3: {'ATA': 'M', 'TGA': 'W', 'CTT': 'T', 'CTC': 'T', 'CTA': 'T', 'CTG': 'T'},
  4: {'TGA': 'W'},
  5: {'ATA': 'M', 'TGA': 'W', 'AGA': 'S', 'AGG': 'S'},
  6: {'TAA': 'Q', 'TAG': 'Q'},
  9: {'AAA': 'N', 'AGA': 'S', 'AGG': 'S', 'TGA': 'W'},
  11: {},
  13: {'AGA': 'G', 'AGG': 'G', 'ATA': 'M', 'TGA': 'W'},
  14: {'AAA': 'N', 'AGA': 'S', 'AGG': 'S', 'TAA': 'Y', 'TGA': 'W'},
  16: {'TAG': 'L'},
  21: {'ATA': 'M', 'TGA': 'W', 'AGA': 'S', 'AGG': 'S', 'AAA': 'N'},
  24: {'AGA': 'S', 'AGG': 'K', 'TGA': 'W'},
};

/// Names shown in the genetic-code selectors (same list as the HTML app).
const geneticCodeNames = {
  1: 'Standard (NCBI 1)',
  2: 'Vertebrate mitochondrial (NCBI 2)',
  3: 'Yeast mitochondrial (NCBI 3)',
  4: 'Mold mitochondrial (NCBI 4)',
  5: 'Invertebrate mitochondrial (NCBI 5)',
  6: 'Ciliate nuclear (NCBI 6)',
  9: 'Echinoderm mitochondrial (NCBI 9)',
  11: 'Bacterial, archaeal and plant plastid (NCBI 11)',
  13: 'Ascidian mitochondrial (NCBI 13)',
  14: 'Alternative flatworm mitochondrial (NCBI 14)',
  16: 'Chlorophycean mitochondrial (NCBI 16)',
  21: 'Trematode mitochondrial (NCBI 21)',
  24: 'Pterobranchia mitochondrial (NCBI 24)',
};

final _tables = <int, Map<String, String>>{};

/// Codon → amino acid table for an NCBI code (standard code if unknown).
Map<String, String> codonTable(int code) =>
    _tables.putIfAbsent(code, () => {..._standard, ...?_changes[code] ?? _changes[1]});

bool isKnownGeneticCode(int code) => _changes.containsKey(code);

String geneticCodeName(int code) => geneticCodeNames[code] ?? 'NCBI $code';

Set<String> stopCodonsFor(int code) => {
      for (final e in codonTable(code).entries)
        if (e.value == '*') e.key,
    };

String translateCodon(String codon, int code) => codonTable(code)[codon.toUpperCase().replaceAll('U', 'T')] ?? 'X';

/// Translates a DNA sequence from position 0. Gapped codons become `-`.
String translate(String dna, int code) {
  final seq = dna.toUpperCase().replaceAll('U', 'T');
  final out = StringBuffer();
  for (var i = 0; i + 3 <= seq.length; i += 3) {
    final codon = seq.substring(i, i + 3);
    out.write(codon == '---' ? '-' : translateCodon(codon, code));
  }
  return out.toString();
}

enum InternalStopPolicy { exclude, allowWithWarning }

class CdsValidation {
  CdsValidation({
    required this.validDna,
    required this.protein,
    required this.normalized,
    required this.terminalStopRemoved,
    required this.internalStops,
    required this.code,
  });

  /// DNA without the terminal stop (protein.length * 3 == validDna.length).
  final String validDna;
  final String protein;

  /// An incomplete terminal stop (T or TA, completed by polyadenylation in
  /// mitochondria) was removed.
  final bool normalized;
  final bool terminalStopRemoved;
  final int internalStops;
  final int code;
}

/// Validates one CDS (already shifted by /codon_start) and translates it.
/// Returns null (and logs why) when it cannot be used in codon mode.
CdsValidation? validateCds(
  String gene,
  String header,
  String dna,
  int code, {
  InternalStopPolicy stopPolicy = InternalStopPolicy.exclude,
  void Function(String)? log,
}) {
  void say(String m) => log?.call(m);
  if (dna.isEmpty) {
    say('[CDS validation] [$gene] SKIP "$header": empty sequence');
    return null;
  }
  if (!isKnownGeneticCode(code)) {
    say('[Translation] [$gene] WARN "$header": unknown transl_table $code, using the standard code');
    code = 1;
  }
  var work = dna.toUpperCase().replaceAll('U', 'T');
  var normalized = false;
  final remainder = work.length % 3;
  if (remainder != 0) {
    final terminal = work.substring(work.length - remainder);
    if ((remainder == 1 && terminal == 'T') || (remainder == 2 && terminal == 'TA')) {
      say('[CDS validation] [$gene] "$header": removed incomplete terminal stop fragment $terminal '
          '(${dna.length} bp → ${dna.length - remainder} bp)');
      work = work.substring(0, work.length - remainder);
      normalized = true;
    } else {
      say('[CDS validation] [$gene] SKIP "$header": length ${dna.length} bp is not a multiple of 3 '
          '(terminal "$terminal" is not an incomplete stop fragment)');
      return null;
    }
  }
  if (work.isEmpty) {
    say('[CDS validation] [$gene] SKIP "$header": empty after normalization');
    return null;
  }
  final table = codonTable(code);
  final protein = StringBuffer();
  var internalStops = 0;
  for (var i = 0; i < work.length - 3; i += 3) {
    final aa = table[work.substring(i, i + 3)] ?? 'X';
    if (aa == '*') internalStops++;
    protein.write(aa);
  }
  final last = work.substring(work.length - 3);
  final lastAa = table[last] ?? 'X';
  final terminalStop = lastAa == '*';
  final validDna = terminalStop ? work.substring(0, work.length - 3) : work;
  if (!terminalStop) protein.write(lastAa);
  say(terminalStop
      ? '[CDS validation] [$gene] "$header": removed terminal stop codon $last (${work.length ~/ 3} codons, code $code)'
      : '[CDS validation] [$gene] "$header": ${work.length ~/ 3} codons, no terminal stop (code $code)');
  if (validDna.isEmpty) {
    say('[CDS validation] [$gene] SKIP "$header": nothing left after removing the terminal stop');
    return null;
  }
  if (internalStops > 0) {
    if (stopPolicy == InternalStopPolicy.exclude) {
      say('[CDS validation] [$gene] SKIP "$header": $internalStops internal stop codon(s) (policy: exclude)');
      return null;
    }
    say('[Translation] [$gene] WARN "$header": $internalStops internal stop codon(s) (policy: allow with warning)');
  }
  return CdsValidation(
    validDna: validDna,
    protein: protein.toString(),
    normalized: normalized,
    terminalStopRemoved: terminalStop,
    internalStops: internalStops,
    code: code,
  );
}

class BackTranslation {
  BackTranslation(this.dna, this.mismatches, this.missingCodons, this.leftoverCodons, this.issues);
  final String dna;
  final int mismatches;
  final int missingCodons;
  final int leftoverCodons;
  final List<String> issues;
  bool get frameDisrupted => missingCodons > 0 || leftoverCodons > 0;
}

/// Back-translates one aligned protein using the original (unaligned) CDS:
/// every residue takes the next codon, every gap becomes `---`.
BackTranslation backTranslate(String alignedProtein, String cds, int code) {
  final dna = cds.toUpperCase().replaceAll('U', 'T').replaceAll(RegExp('[^ACGTNRYSWKMBDHV]'), '');
  final out = StringBuffer();
  final issues = <String>[];
  var idx = 0, mismatches = 0, missing = 0;
  final codons = dna.length ~/ 3;
  for (var pos = 0; pos < alignedProtein.length; pos++) {
    final aa = alignedProtein[pos].toUpperCase();
    if (aa == '-' || aa == '.') {
      out.write('---');
      continue;
    }
    if (idx >= codons) {
      out.write('NNN');
      missing++;
      issues.add('AA position ${pos + 1}: no CDS codon left for $aa');
      continue;
    }
    final codon = dna.substring(idx * 3, idx * 3 + 3);
    idx++;
    final translated = translateCodon(codon, code);
    if (aa != 'X' && aa != '?' && aa != translated) {
      mismatches++;
      if (issues.length < 20) issues.add('AA position ${pos + 1}: protein $aa but $codon is $translated (code $code)');
    }
    out.write(codon);
  }
  final leftover = codons - idx;
  if (leftover > 0) issues.add('$leftover CDS codon(s) left unused after back-translation');
  return BackTranslation(out.toString(), mismatches, missing, leftover, issues);
}

class CodonIntegrity {
  CodonIntegrity({required this.equalLengths, required this.lengthMultipleOf3, required this.gapsInTriplets, required this.length});
  final bool equalLengths;
  final bool lengthMultipleOf3;
  final bool gapsInTriplets;
  final int length;
  bool get ok => equalLengths && lengthMultipleOf3 && gapsInTriplets;

  List<String> issues({bool frameDisrupted = false}) => [
        if (!equalLengths) 'sequences have unequal lengths',
        if (!lengthMultipleOf3) 'alignment length $length is not divisible by 3',
        if (!gapsInTriplets) 'gap runs are not multiples of 3',
        if (frameDisrupted) 'codon frame disruption detected during back-translation',
      ];
}

/// Checks that a back-translated alignment still respects the codon frame.
CodonIntegrity codonIntegrity(List<SeqRecord> aln) {
  final first = aln.isEmpty ? 0 : aln.first.length;
  final gapRun = RegExp('-+');
  return CodonIntegrity(
    equalLengths: aln.every((r) => r.length == first),
    lengthMultipleOf3: aln.every((r) => r.length % 3 == 0),
    gapsInTriplets: aln.every((r) => gapRun.allMatches(r.sequence).every((m) => (m.end - m.start) % 3 == 0)),
    length: first,
  );
}

/// Back-translates a whole protein alignment. Sequences are matched by
/// header; proteins without a CDS are dropped and reported.
({List<SeqRecord> dna, List<String> issues, int mismatches, bool frameDisrupted}) backTranslateAlignment(
  List<SeqRecord> proteinAln,
  Map<String, String> cdsByHeader,
  Map<String, int> codeByHeader,
  int fallbackCode,
) {
  final out = <SeqRecord>[];
  final issues = <String>[];
  var mismatches = 0;
  var disrupted = false;
  for (final p in proteinAln) {
    final cds = cdsByHeader[p.header];
    if (cds == null) {
      issues.add('${p.header}: no CDS to back-translate');
      continue;
    }
    final bt = backTranslate(p.sequence, cds, codeByHeader[p.header] ?? fallbackCode);
    mismatches += bt.mismatches;
    disrupted |= bt.frameDisrupted;
    for (final i in bt.issues.take(3)) {
      issues.add('${p.header}: $i');
    }
    out.add(SeqRecord(p.header, bt.dna));
  }
  return (dna: out, issues: issues, mismatches: mismatches, frameDisrupted: disrupted);
}

/// Removes the given codon positions (1-3) from an in-frame alignment.
List<SeqRecord> removeCodonPositions(List<SeqRecord> aln, Set<int> positions, {int frame = 0}) => [
      for (final r in aln)
        SeqRecord(r.header, [
          for (var i = 0; i < r.length; i++)
            if (!positions.contains((i - frame) % 3 + 1)) r.sequence[i],
        ].join()),
    ];
