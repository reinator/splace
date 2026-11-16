import 'package:flutter/material.dart';
import 'package:flutter/services.dart';
import 'package:url_launcher/url_launcher.dart';

import '../common.dart';
import '../fa.dart';
import '../../services/tool_locator.dart' show ExternalTool, ExternalToolInfo;
import '../methods_text.dart';
import '../theme.dart';
import 'kit.dart';
import 'splace_tooltip.dart';

class Citation {
  const Citation(this.name, this.authors, this.year, this.title, this.journal, this.volume, this.doi, {this.url});
  final String name, authors, year, title, journal, volume, doi;

  /// Source code (GitHub) of the tool, when it has one.
  final String? url;

  String get text => '$authors ($year) $title. $journal $volume. https://doi.org/$doi${url == null ? '' : ' · $url'}';
}

/// Tools used by SPLACE (docs/citations.js, plus RAxML-NG, ClipKIT, the
/// IQ-TREE methods, the distances and Xia's saturation index).
const citations = [
  Citation('SPLACE', 'Oliveira RRM, Vasconcelos S, Oliveira G', '2022',
      'SPLACE: A tool to automatically SPLit, Align, and ConcatenatE genes for phylogenomic inference of several organisms', 'Frontiers in Bioinformatics',
      '2:1074802', '10.3389/fbinf.2022.1074802',
      url: 'https://github.com/luanrabelo/splace'),
  Citation('SynGenes', 'Rabelo LP, Sodré D, de Sousa RPC, Watanabe L, Gomes G, Sampaio I, Vallinoto M', '2024',
      'SynGenes: a Python class for standardizing nomenclatures of mitochondrial and chloroplast genes and a web form for enhancing searches for evolutionary analyses',
      'BMC Bioinformatics', '25:160', '10.1186/s12859-024-05781-y',
      url: 'https://github.com/luanrabelo/SynGenes'),
  Citation('dataFishing', 'Rabelo L, Sodré D, Balcázar ODA, do Rosário MF, Guimarães-Costa AJ, Gomes G, Sampaio I, Vallinoto M', '2025',
      'dataFishing: An efficient Python tool and user-friendly web-form for mining mitochondrial and chloroplast sequences, taxonomic, and biodiversity data',
      'Ecological Informatics', '85:102970', '10.1016/j.ecoinf.2024.102970',
      url: 'https://github.com/luanrabelo/dataFishing'),
  Citation('MAFFT', 'Katoh K, Standley DM', '2013', 'MAFFT Multiple Sequence Alignment Software Version 7: Improvements in Performance and Usability',
      'Molecular Biology and Evolution', '30(4):772–780', '10.1093/molbev/mst010'),
  Citation('trimAl', 'Capella-Gutiérrez S, Silla-Martínez JM, Gabaldón T', '2009',
      'trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses', 'Bioinformatics', '25(15):1972–1973', '10.1093/bioinformatics/btp348'),
  Citation('ClipKIT', 'Steenwyk JL, Buida TJ, Li Y, Shen XX, Rokas A', '2020',
      'ClipKIT: A multiple sequence alignment trimming software for accurate phylogenomic inference', 'PLoS Biology', '18(12):e3001007',
      '10.1371/journal.pbio.3001007'),
  Citation('IQ-TREE 3',
      'Wong TKF, Ly-Trong N, Ren H, Demotte P, Baños H, Roger AJ, Susko E, Bielow C, De Maio N, Goldman N, Hahn MW, dos Reis M, Vinh LS, Huttley G, Lanfear R, Minh BQ',
      '2026', 'IQ-TREE 3: phylogenomic inference software using complex evolutionary models', 'Molecular Biology and Evolution', '43(5):msag117',
      '10.1093/molbev/msag117'),
  Citation('IQ-TREE 2', 'Minh BQ, Schmidt HA, Chernomor O, Schrempf D, Woodhams MD, von Haeseler A, Lanfear R', '2020',
      'IQ-TREE 2: New Models and Methods for Phylogenetic Inference', 'Molecular Biology and Evolution', '37(5):1530–1534', '10.1093/molbev/msaa015'),
  Citation('RAxML-NG', 'Kozlov AM, Darriba D, Flouri T, Morel B, Stamatakis A', '2019',
      'RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference', 'Bioinformatics', '35(21):4453–4455',
      '10.1093/bioinformatics/btz305'),
  Citation('Adaptive RAxML-NG', 'Togkousidis A, Kozlov OM, Haag J, Höhler D, Stamatakis A', '2023',
      'Adaptive RAxML-NG: Accelerating Phylogenetic Inference under Maximum Likelihood using Dataset Difficulty', 'Molecular Biology and Evolution', '40(10):msad227',
      '10.1093/molbev/msad227'),
  Citation('RAxML-NG fast', 'Togkousidis A, Stamatakis A, Gascuel O', '2025',
      'Accelerating Maximum Likelihood Phylogenetic Inference via Early Stopping to Evade (Over-)optimization', 'Systematic Biology', '74(6):1020–1028',
      '10.1093/sysbio/syaf043'),
  Citation('RBS', 'Stamatakis A, Hoover P, Rougemont J', '2008', 'A Rapid Bootstrap Algorithm for the RAxML Web Servers', 'Systematic Biology', '57(5):758–771',
      '10.1080/10635150802429642'),
  Citation('EBG', 'Wiegert J, Höhler D, Haag J, Stamatakis A', '2024', 'Predicting Phylogenetic Bootstrap Values via Machine Learning',
      'Molecular Biology and Evolution', '41(10):msae215', '10.1093/molbev/msae215'),
  Citation('ModelFinder', 'Kalyaanamoorthy S, Minh BQ, Wong TKF, von Haeseler A, Jermiin LS', '2017',
      'ModelFinder: fast model selection for accurate phylogenetic estimates', 'Nature Methods', '14:587–589', '10.1038/nmeth.4285'),
  Citation('rcluster', 'Lanfear R, Calcott B, Kainer D, Mayer C, Stamatakis A', '2014', 'Selecting optimal partitioning schemes for phylogenomic datasets',
      'BMC Evolutionary Biology', '14:82', '10.1186/1471-2148-14-82'),
  Citation('Partition merging', 'Chernomor O, von Haeseler A, Minh BQ', '2016', 'Terrace aware data structure for phylogenomic inference from supermatrices',
      'Systematic Biology', '65(6):997–1008', '10.1093/sysbio/syw037'),
  Citation('UFBoot2', 'Hoang DT, Chernomor O, von Haeseler A, Minh BQ, Vinh LS', '2018', 'UFBoot2: Improving the Ultrafast Bootstrap Approximation',
      'Molecular Biology and Evolution', '35(2):518–522', '10.1093/molbev/msx281'),
  Citation('SH-aLRT', 'Guindon S, Dufayard JF, Lefort V, Anisimova M, Hordijk W, Gascuel O', '2010',
      'New algorithms and methods to estimate maximum-likelihood phylogenies: assessing the performance of PhyML 3.0', 'Systematic Biology', '59(3):307–321',
      '10.1093/sysbio/syq010'),
  Citation('TBE', 'Lemoine F, Domelevo Entfellner JB, Wilkinson E, Correia D, Dávila Felipe M, De Oliveira T, Gascuel O', '2018',
      "Renewing Felsenstein's phylogenetic bootstrap in the era of big data", 'Nature', '556:452–456', '10.1038/s41586-018-0043-0'),
  Citation('TN93', 'Tamura K, Nei M', '1993',
      'Estimation of the number of nucleotide substitutions in the control region of mitochondrial DNA in humans and chimpanzees', 'Molecular Biology and Evolution',
      '10(3):512–526', '10.1093/oxfordjournals.molbev.a040023'),
  Citation('K2P', 'Kimura M', '1980',
      'A simple method for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences', 'Journal of Molecular Evolution',
      '16:111–120', '10.1007/BF01731581'),
  Citation('PhyKIT', 'Steenwyk JL, Buida TJ, Labella AL, Li Y, Shen XX, Rokas A', '2021',
      'PhyKIT: a broadly applicable UNIX shell toolkit for processing and analyzing phylogenomic data', 'Bioinformatics', '37(16):2325–2331',
      '10.1093/bioinformatics/btab096',
      url: 'https://github.com/JLSteenwyk/PhyKIT'),
  Citation('Philippe 2011', 'Philippe H, Brinkmann H, Lavrov DV, Littlewood DTJ, Manuel M, Wörheide G, Baurain D', '2011',
      'Resolving difficult phylogenetic questions: why more sequences are not enough', 'PLoS Biology', '9(3):e1000602', '10.1371/journal.pbio.1000602'),
  Citation('Iss (saturation)', 'Xia X, Xie Z, Salemi M, Chen L, Wang Y', '2003', 'An index of substitution saturation and its application',
      'Molecular Phylogenetics and Evolution', '26(1):1–7', '10.1016/S1055-7903(02)00326-3'),
];

/// Works always listed in "How to cite" (the tools SPLACE is built on).
const _coreCitations = ['SPLACE', 'SynGenes', 'dataFishing', 'MAFFT', 'trimAl', 'ClipKIT', 'PhyKIT', 'IQ-TREE 3', 'RAxML-NG'];

/// "How to Cite" card shown after the tree: a methods text written from
/// what was done in the project (academic language, tools cited), ready to
/// copy, and the references — those cited in the text, or the tools of
/// SPLACE while nothing has been run.
class CitationsCard extends StatelessWidget {
  const CitationsCard({super.key});

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final l10n = context.l10n;
    final c = ProjectScope.of(context);
    final methods = buildMethodsText(c, Localizations.localeOf(context).languageCode);
    final refs = methods.cited.isEmpty ? _coreCitations : methods.cited;
    return Container(
      padding: const EdgeInsets.all(16),
      decoration: BoxDecoration(color: t.accentSoft, borderRadius: BorderRadius.circular(14), border: Border.all(color: t.accentSoftBorder)),
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        Row(children: [
          const StepBadge('9c', size: 26),
          const SizedBox(width: 10),
          Text(l10n.howToCite, style: TextStyle(fontWeight: FontWeight.w600, fontSize: 14, color: t.accentText)),
          const SizedBox(width: 12),
          Expanded(child: Text(l10n.howToCiteHint, textAlign: TextAlign.right, style: TextStyle(fontSize: 12, color: t.textMuted))),
        ]),
        if (!methods.isEmpty) ...[
          const SizedBox(height: 12),
          Container(
            padding: const EdgeInsets.fromLTRB(16, 12, 12, 14),
            decoration: BoxDecoration(color: t.surface, borderRadius: BorderRadius.circular(10), border: Border.all(color: t.accentSoftBorder)),
            child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
              Row(children: [
                Icon(FA.fileLines, size: 13, color: t.accent),
                const SizedBox(width: 8),
                Expanded(child: Text(l10n.methodsTitle, style: TextStyle(fontWeight: FontWeight.w700, fontSize: 14, color: t.accentText))),
                ChipButton(
                  label: l10n.methodsCopy,
                  style: ChipButtonStyle.accent,
                  onPressed: () => Clipboard.setData(ClipboardData(text: methods.plain)),
                ),
              ]),
              const SizedBox(height: 4),
              Text(l10n.methodsHelp, style: TextStyle(fontSize: 12, color: t.textMuted, height: 1.4)),
              const SizedBox(height: 10),
              SelectionArea(
                child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
                  for (final p in methods.paragraphs)
                    Padding(
                      padding: const EdgeInsets.only(bottom: 10),
                      child: Text(p, textAlign: TextAlign.justify, style: TextStyle(fontSize: 13.5, color: t.text, height: 1.6)),
                    ),
                ]),
              ),
            ]),
          ),
        ],
        const SizedBox(height: 14),
        // The versions that ran here, as they should be reported.
        Kicker(l10n.citeVersions, size: 10.5),
        const SizedBox(height: 8),
        Wrap(spacing: 8, runSpacing: 6, children: [
          TagChip('SPLACE v$appVersion', color: t.accentText),
          for (final tool in ExternalTool.values)
            if (c.tools[tool] case final r?) TagChip('${tool.label} ${r.version == null ? '' : 'v${r.version}'}'.trim(), color: t.accentText),
        ]),
        const SizedBox(height: 14),
        Kicker(methods.cited.isEmpty ? l10n.citeTools : l10n.citeReferences, size: 10.5),
        const SizedBox(height: 8),
        for (final name in refs)
          if (citations.where((x) => x.name == name).firstOrNull case final c?)
            Container(
              margin: const EdgeInsets.only(bottom: 8),
              padding: const EdgeInsets.all(12),
              decoration: BoxDecoration(color: t.surface, borderRadius: BorderRadius.circular(10), border: Border.all(color: t.accentSoftBorder)),
              child: Row(children: [
                Expanded(
                  child: Column(crossAxisAlignment: CrossAxisAlignment.start, children: [
                    Text(c.name, style: TextStyle(fontWeight: FontWeight.w700, color: t.accentText)),
                    const SizedBox(height: 2),
                    Text.rich(TextSpan(style: TextStyle(fontSize: 13, color: t.text, height: 1.45), children: [
                      TextSpan(text: '${c.authors} (${c.year}) '),
                      TextSpan(text: c.title, style: const TextStyle(fontStyle: FontStyle.italic)),
                      TextSpan(text: '. ${c.journal} ${c.volume}. '),
                      WidgetSpan(
                        child: InkWell(
                          onTap: () => launchUrl(Uri.parse('https://doi.org/${c.doi}')),
                          child: Text('DOI: ${c.doi}', style: TextStyle(fontSize: 13, color: t.accent, decoration: TextDecoration.underline)),
                        ),
                      ),
                      if (c.url != null) ...[
                        const TextSpan(text: '  ·  '),
                        WidgetSpan(
                          child: InkWell(
                            onTap: () => launchUrl(Uri.parse(c.url!)),
                            child: Row(mainAxisSize: MainAxisSize.min, children: [
                              Icon(FA.github, size: 12, color: t.accent),
                              const SizedBox(width: 4),
                              Text(c.url!.replaceFirst('https://', ''),
                                  style: TextStyle(fontSize: 13, color: t.accent, decoration: TextDecoration.underline)),
                            ]),
                          ),
                        ),
                      ],
                    ])),
                  ]),
                ),
                const SizedBox(width: 10),
                ChipButton(label: l10n.copy, onPressed: () => Clipboard.setData(ClipboardData(text: c.text))),
              ]),
            ),
      ]),
    );
  }
}

/// "Katoh & Standley", "Rabelo et al.": first author's surname, as in the
/// HTML references.
String shortAuthors(String authors) {
  if (authors.contains('et al')) return authors.replaceAll(RegExp(r'\.$'), '.');
  final names = authors.split(',').map((a) => a.trim().split(RegExp(r'\s+')).first).toList();
  return switch (names.length) { 1 => names.first, 2 => '${names[0]} & ${names[1]}', _ => '${names.first} et al.' };
}

/// Inline reference to a tool of [citations] ("**SynGenes** (Rabelo et al.
/// 2024)"): full citation in the tooltip (after [note]), DOI on click.
/// [label] replaces the tool name; with [pill], drawn as a small rounded
/// badge.
class CitationRef extends StatelessWidget {
  const CitationRef(this.name, {super.key, this.fontSize = 13, this.note, this.label, this.pill = false});
  final String name;
  final double fontSize;
  final String? note;
  final String? label;
  final bool pill;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final c = citations.firstWhere((c) => c.name == name);
    Widget text = Text.rich(TextSpan(style: TextStyle(fontSize: fontSize, color: t.accentText, height: 1.2), children: [
      TextSpan(text: label ?? c.name, style: const TextStyle(fontWeight: FontWeight.w700)),
      TextSpan(text: ' (${shortAuthors(c.authors)} ${c.year})'),
    ]));
    if (pill) {
      text = Container(
        padding: const EdgeInsets.symmetric(horizontal: 8, vertical: 2),
        decoration: BoxDecoration(color: t.accentSoftBorder, borderRadius: BorderRadius.circular(8), border: Border.all(color: t.accentSoftBorder)),
        child: text,
      );
    }
    return SplaceTooltip(
      message: note == null ? c.text : '$note ${c.name}.\n\n${c.text}',
      child: MouseRegion(
        cursor: SystemMouseCursors.click,
        child: GestureDetector(onTap: () => launchUrl(Uri.parse('https://doi.org/${c.doi}')), child: text),
      ),
    );
  }
}

/// [text] with `[[Tool]]` or `[[Tool|label]]` markers turned into
/// [CitationRef] pills (names from [citations]).
List<InlineSpan> citationSpans(String text, {double fontSize = 12}) {
  final out = <InlineSpan>[];
  var start = 0;
  for (final m in RegExp(r'\[\[([^\]|]+)(?:\|([^\]]+))?\]\]').allMatches(text)) {
    if (m.start > start) out.add(TextSpan(text: text.substring(start, m.start)));
    out.add(WidgetSpan(
      alignment: PlaceholderAlignment.middle,
      child: CitationRef(m.group(1)!, label: m.group(2), pill: true, fontSize: fontSize),
    ));
    start = m.end;
  }
  if (start < text.length) out.add(TextSpan(text: text.substring(start)));
  return out;
}
