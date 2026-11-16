import 'package:flutter/material.dart';
import 'package:url_launcher/url_launcher.dart';

import '../../core/genbank.dart';
import '../common.dart';
import '../fa.dart';
import '../theme.dart';
import 'dialogs.dart';
import 'kit.dart';
import 'splace_tooltip.dart';

/// IUCN Red List colours of each category (background, text), as in the
/// IUCN and ICMBio publications. RE (regionally extinct) and NA (not
/// applicable) are regional categories used by SALVE.
(Color, Color) iucnColors(String code) => switch (code.toUpperCase()) {
      'EX' => (const Color(0xFF000000), Colors.white),
      'EW' => (const Color(0xFF542344), Colors.white),
      'RE' => (const Color(0xFF9B4F96), Colors.white),
      'CR' => (const Color(0xFFD81E05), Colors.white),
      'EN' => (const Color(0xFFFC7F3F), Colors.black),
      'VU' => (const Color(0xFFF9E814), Colors.black),
      'NT' => (const Color(0xFFCCE226), Colors.black),
      'LC' => (const Color(0xFF60C659), Colors.black),
      'DD' => (const Color(0xFFD1D1C6), Colors.black),
      'NA' => (const Color(0xFFC1B5A5), Colors.black),
      _ => (Colors.white, Colors.black), // NE and unknown
    };

/// The category code on its IUCN colour, with the name after it.
class IucnBadge extends StatelessWidget {
  const IucnBadge({super.key, required this.code, this.name = '', this.possiblyExtinct = false});
  final String code, name;
  final bool possiblyExtinct;

  @override
  Widget build(BuildContext context) {
    final t = context.tokens;
    final (bg, fg) = iucnColors(code);
    final label = code.isEmpty ? '?' : code.toUpperCase();
    return SplaceTooltip(
      message: [if (name.isNotEmpty) '$name ($label)' else label, if (possiblyExtinct) context.l10n.salvePossiblyExtinct].join('\n'),
      child: Row(mainAxisSize: MainAxisSize.min, children: [
        Container(
          padding: const EdgeInsets.symmetric(horizontal: 7, vertical: 2),
          decoration: BoxDecoration(
            color: bg,
            borderRadius: BorderRadius.circular(6),
            border: Border.all(color: bg == Colors.white ? t.border : bg),
          ),
          child: Text(possiblyExtinct ? '$label (PEX)' : label,
              style: TextStyle(fontSize: 11.5, fontWeight: FontWeight.w800, color: fg, letterSpacing: 0.3)),
        ),
        if (name.isNotEmpty) ...[
          const SizedBox(width: 6),
          Flexible(child: Text(name, maxLines: 1, overflow: TextOverflow.ellipsis, style: TextStyle(fontSize: 12.5, color: t.textMuted))),
        ],
      ]),
    );
  }
}

/// Everything the SALVE sheet of [species] says.
Future<void> showSalveSheetDialog(BuildContext context, String species, SalveSheet s) {
  final l10n = context.l10n;
  return showDialog(
    context: context,
    builder: (context) {
      final t = context.tokens;
      Widget line(String label, String value) => value.isEmpty ? const SizedBox.shrink() : InfoLine(label, value);
      Widget block(String label, String text) =>
          text.isEmpty ? const SizedBox.shrink() : DialogBlock(label: label, child: SelectableText(text, style: TextStyle(fontSize: 13, height: 1.5, color: t.text)));
      Widget items(String label, List<String> xs) => xs.isEmpty
          ? const SizedBox.shrink()
          : DialogBlock(
              label: '$label (${xs.length})',
              child: SelectableText(xs.map((x) => '• $x').join('\n'), style: TextStyle(fontSize: 13, height: 1.5, color: t.text)),
            );
      String text(String k) => s.texts[k] ?? '';
      return SplaceDialog(
        icon: FA.leaf,
        title: species,
        subtitle: l10n.salveSheetSubtitle,
        actions: [
          if (s.doi.isNotEmpty) OutlinedButton(onPressed: () => launchUrl(Uri.parse('https://doi.org/${s.doi}')), child: const Text('DOI')),
          FilledButton.icon(
            onPressed: () => launchUrl(Uri.parse(s.url)),
            icon: const Icon(FA.arrowUpRightFromSquare, size: 12),
            label: Text(l10n.salveOpenSheet),
          ),
        ],
        child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
          Row(children: [
            IucnBadge(code: s.code, name: s.categoryName, possiblyExtinct: s.possiblyExtinct),
            if (s.criteria.isNotEmpty) ...[const SizedBox(width: 10), TagChip('${l10n.colCriteria}: ${s.criteria}', color: t.textMuted)],
          ]),
          const SizedBox(height: 14),
          // Assessment and taxonomy.
          line(l10n.colCommonNames, s.commonNames),
          line(l10n.colSalveGroup, s.group),
          line(l10n.salveSynonyms, s.synonyms),
          line(l10n.colAssessed, s.assessed),
          line(l10n.salvePublished, s.published),
          line(l10n.salveSituation, s.situation),
          line(l10n.colOfficialList, s.officialList),
          line('DOI', s.doi),
          line(l10n.salveReferences, s.references == 0 ? '' : '${s.references}'),
          const SizedBox(height: 10),
          block(l10n.salveJustification, s.justification),
          items(l10n.salveHistory, s.history),
          items(l10n.salveLists, s.officialLists),
          // Distribution.
          Kicker(l10n.salveSecDistribution),
          const SizedBox(height: 8),
          line(l10n.colEndemic, s.endemic),
          line(l10n.colBiomes, s.biomes),
          line(l10n.colStates, s.states),
          line(l10n.salveWatersheds, s.watersheds),
          block(l10n.salveNational, text('nationalDistribution')),
          block(l10n.salveGlobal, text('globalDistribution')),
          // Natural history and population.
          Kicker(l10n.salveSecNatural),
          const SizedBox(height: 8),
          line(l10n.colMigratory, s.migratory),
          line(l10n.colHabitats, s.habitats.join(', ')),
          line(l10n.salveDiet, s.foodHabits.join(', ')),
          line(l10n.salveReproduction, s.reproduction),
          line(l10n.colTrend, s.populationTrend),
          line(l10n.colGeneration, s.generationTime),
          line(l10n.salveSexRatio, s.sexRatio),
          block(l10n.salveNaturalHistory, text('naturalHistory')),
          block(l10n.colHabitats, text('habitat')),
          block(l10n.salveDiet, text('food')),
          block(l10n.salveReproduction, text('reproduction')),
          items(l10n.salveInteractions, s.interactions),
          block(l10n.salveInteractions, text('interactions')),
          block(l10n.salvePopulation, text('population')),
          block(l10n.salveGenetics, text('genetics')),
          // Threats, uses, conservation, research.
          Kicker(l10n.salveSecThreats),
          const SizedBox(height: 8),
          items(l10n.colThreats, s.threats),
          block(l10n.colThreats, text('threats')),
          items(l10n.salveUses, s.uses),
          block(l10n.salveUses, text('uses')),
          items(l10n.salveActions, s.actions),
          block(l10n.salveConservation, text('conservation')),
          items(l10n.colUcs, s.protectedAreas),
          items(l10n.salveResearch, s.researchNeeds),
          block(l10n.salveResearch, text('research')),
          // Authorship.
          block(l10n.salveAssessors, s.assessors),
          block(l10n.salveCitation, s.citation),
        ]),
      );
    },
  );
}
