import 'dart:async';

import 'package:desktop_drop/desktop_drop.dart';
import 'package:file_selector/file_selector.dart';
import 'package:flutter/material.dart';

import '../../core/taxonomy_sources.dart' show isSpeciesName;
import '../../core/bold.dart' show boldMarkerCodes, boldServes;
import '../../core/marker_search.dart';
import '../../core/syngenes.dart';
import '../../services/ncbi.dart';
import '../chrome/api_key_dialog.dart';
import '../chrome/genome_review_dialog.dart';
import '../common.dart';
import '../fa.dart';
import '../theme.dart';
import '../widgets/example_data.dart';
import '../widgets/html_text.dart';
import '../widgets/kit.dart';
import '../widgets/splace_tooltip.dart';

const fastaTypes = XTypeGroup(label: 'FASTA', extensions: ['fasta', 'fas', 'fa', 'fna', 'ffn', 'fsa', 'aln', 'txt']);
const genbankTypes = XTypeGroup(label: 'GenBank', extensions: ['gb', 'gbk', 'genbank', 'gbff']);

/// Step 1 — Load Data (`#step1Section`): NCBI taxonomy search, accession
/// numbers, local GenBank files and local FASTA files.
class LoadSection extends StatefulWidget {
  const LoadSection({super.key, required this.active});
  final bool active;

  @override
  State<LoadSection> createState() => _LoadSectionState();
}

class _LoadSectionState extends State<LoadSection> {
  final _taxon = TextEditingController();
  final _accessions = TextEditingController();
  Organelle _organelle = Organelle.mitochondrion;
  bool _completeOnly = false;
  bool _refseqOnly = false;

  /// Option 1 searches individual markers (COI, 12S, CYTB…) instead of
  /// genomes.
  bool _markerMode = false;
  final Set<String> _genes = {'COI', '12S', 'CYTB'};
  bool _includeGenomes = false;

  /// Individual markers: BOLD Systems too, for COI and the plant barcodes.
  bool _bold = true;

  static const _defaultGenes = {Organelle.mitochondrion: {'COI', '12S', 'CYTB'}, Organelle.chloroplast: {'rbcL', 'matK'}};

  OrganelleType get _type => _organelle == Organelle.chloroplast ? OrganelleType.cp : OrganelleType.mt;

  MarkerSearch get _markerSearch => MarkerSearch(
        taxon: _taxon.text.trim(),
        genes: [for (final g in syngenesMarkers(_type)) if (_genes.contains(g)) g],
        organelle: _type,
        includeGenomes: _includeGenomes,
        refseqOnly: _refseqOnly,
        synonyms: _synonymsInUse,
        bold: _bold,
      );
  bool _dragGb = false;
  bool _dragFa = false;

  /// A species name in 1a: searched with its synonyms from GBIF when
  /// [_withSynonyms] is on.
  bool _withSynonyms = false;
  List<String>? _synonyms;
  String? _synonymsFor;
  bool _synonymsLoading = false;
  String? _synonymsError;
  Timer? _synonymsDebounce;

  bool get _isSpecies => isSpeciesName(_taxon.text);

  /// Example files bundled in assets/examples, by kind.
  Map<ExampleKind, List<String>> _examples = const {};

  @override
  void initState() {
    super.initState();
    listExampleAssets().then((e) {
      if (mounted) setState(() => _examples = e);
    });
    loadExampleTooltips().then((t) {
      if (mounted) setState(() => _tips = t);
    });
  }

  /// Tooltips of the example buttons (assets/examples/tooltips.json, HTML).
  Map<String, Map<String, String>> _tips = const {};

  /// Synonyms in use for the current name (empty when off or not found).
  List<String> get _synonymsInUse => _withSynonyms && _isSpecies && _synonymsFor == _taxon.text.trim() ? (_synonyms ?? const []) : const [];

  /// Fetches the synonyms of the typed species from GBIF.
  Future<void> _loadSynonyms() async {
    final name = _taxon.text.trim();
    if (!_withSynonyms || !isSpeciesName(name) || _synonymsFor == name) return;
    setState(() {
      _synonymsLoading = true;
      _synonymsError = null;
    });
    final c = ProjectScope.read(context);
    try {
      final names = await c.ncbi.gbifSynonyms(name);
      if (!mounted || _taxon.text.trim() != name) return;
      setState(() {
        _synonyms = names;
        _synonymsFor = name;
      });
      c.log(names.isEmpty ? 'GBIF does not know $name: searching the name alone' : 'GBIF synonyms of $name: ${names.join(', ')}', source: 'GBIF');
    } catch (e) {
      if (mounted) setState(() => _synonymsError = '$e');
    } finally {
      if (mounted) setState(() => _synonymsLoading = false);
    }
  }

  void _onTaxonChanged() {
    setState(() {});
    _synonymsDebounce?.cancel();
    if (_withSynonyms && _isSpecies) _synonymsDebounce = Timer(const Duration(milliseconds: 700), _loadSynonyms);
  }

  @override
  void dispose() {
    _synonymsDebounce?.cancel();
    _taxon.dispose();
    _accessions.dispose();
    super.dispose();
  }

  GenomeSearch get _search =>
      GenomeSearch(taxon: _taxon.text.trim(), organelle: _organelle, completeOnly: _completeOnly, refseqOnly: _refseqOnly, synonyms: _synonymsInUse);

  Future<void> _runMarkerSearch() async {
    final c = ProjectScope.read(context);
    final search = _markerSearch;
    final candidates = await c.searchMarkers(search);
    if (!mounted || candidates == null || candidates.isEmpty) return;
    final cancelled = context.l10n.reviewCancelled;
    final selected = await showGenomeReviewDialog(context, _search, candidates, markers: search);
    if (selected == null || selected.isEmpty) {
      c.log(cancelled, source: 'NCBI');
      return;
    }
    await c.importMarkers(selected, search);
  }

  Future<void> _runSearch() async {
    if (_markerMode) return _runMarkerSearch();
    final c = ProjectScope.read(context);
    final search = _search;
    final candidates = await c.searchGenomes(search);
    if (!mounted || candidates == null || candidates.isEmpty) return;
    final cancelled = context.l10n.reviewCancelled;
    final selected = await showGenomeReviewDialog(context, search, candidates);
    if (selected == null || selected.isEmpty) {
      c.log(cancelled, source: 'NCBI');
      return;
    }
    await c.importGenomes(selected, search.taxon);
  }

  @override
  Widget build(BuildContext context) {
    final c = ProjectScope.of(context);
    final l10n = context.l10n;
    final t = context.tokens;
    final sequences = c.inputs.fold<int>(0, (a, f) => a + f.records.length);
    final taxonOk = _taxon.text.trim().length >= 3 && (!_markerMode || _genes.isNotEmpty);
    final accessionOk = accessionPattern.hasMatch(_accessions.text);

    Future<void> browseGb() async {
      final files = await openFiles(acceptedTypeGroups: [genbankTypes]);
      if (files.isNotEmpty) await c.addGenBankFiles(files.map((f) => f.path).toList());
    }

    Future<void> browseFa() async {
      final files = await openFiles(acceptedTypeGroups: [fastaTypes]);
      if (files.isNotEmpty) await c.addFiles(files.map((f) => f.path).toList());
    }

    Future<void> folderGb() async {
      final dir = await getDirectoryPath();
      if (dir != null) await c.addGenBankDirectory(dir);
    }

    Future<void> folderFa() async {
      final dir = await getDirectoryPath();
      if (dir != null) await c.addDirectory(dir);
    }

    Widget example(VoidCallback onTap, String tooltip) => SplaceTooltip(
          message: tooltip,
          child: ChipButton(label: l10n.example, icon: FA.wandMagicSparkles, style: ChipButtonStyle.soft, onPressed: onTap),
        );
    final lang = Localizations.localeOf(context).languageCode;

    /// An example button whose tooltip comes from tooltips.json ([key]),
    /// in HTML; [fallback] when the file has none.
    Widget editableExample(String key, String label, VoidCallback? onTap, String fallback) {
      final html = exampleTooltip(_tips, key, lang);
      return SplaceTooltip(
        message: html == null ? fallback : '',
        richMessage: html == null ? null : htmlSpan(html),
        child: ChipButton(label: label, icon: FA.wandMagicSparkles, style: ChipButtonStyle.soft, onPressed: onTap),
      );
    }

    /// Example data bundled in assets/examples/<kind>: the button loads it
    /// (disabled, saying where to put the files, while the folder is empty).
    Widget dataExample(ExampleKind k, String key, String label) {
      final files = _examples[k] ?? const <String>[];
      return editableExample(
        key,
        label,
        c.busy || files.isEmpty ? null : () => loadExample(c, k, files),
        files.isEmpty ? l10n.examplesEmpty(k.source) : l10n.examplesFiles('${files.length}', k.source),
      );
    }


    final option1 = OptionCard(
      badge: '1a',
      title: l10n.taxTitle,
      info: l10n.taxInfo,
      text: l10n.taxBody,
      action: example(() {
        setState(() {
          _taxon.text = 'Carangidae';
          _organelle = Organelle.mitochondrion;
          _refseqOnly = true;
          _completeOnly = true;
        });
      }, l10n.taxExampleTooltip),
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        if (c.ncbiApiKey == null && !c.apiHintDismissed) ...[
          TonePanel(
            tone: PanelTone.info,
            child: Row(children: [
              const Icon(FA.key, size: 12),
              const SizedBox(width: 8),
              Expanded(child: Text(l10n.taxApiHint, style: const TextStyle(fontSize: 13))),
              TextButton(onPressed: () => showApiKeyDialog(context), child: Text(l10n.apiKeyTitle)),
              IconBtn(icon: FA.xmark, tooltip: l10n.close, size: ChipButton.height, onPressed: c.dismissApiHint),
            ]),
          ),
          const SizedBox(height: 12),
        ],
        // What to search: whole genomes or individual markers.
        Align(
          alignment: Alignment.centerLeft,
          child: Seg<bool>(
            value: _markerMode,
            options: [(false, l10n.searchGenomes), (true, l10n.searchMarkers)],
            onChanged: (v) => setState(() => _markerMode = v),
          ),
        ),
        const SizedBox(height: 12),
        // Name, genome type, scope and database side by side; labels share
        // the same top and every control is 44 px high.
        Row(crossAxisAlignment: CrossAxisAlignment.start, children: [
          Expanded(
            flex: 3,
            child: Field(
              label: l10n.taxName,
              child: TextField(
                controller: _taxon,
                // A species (two names) is written in italics.
                style: _isSpecies ? const TextStyle(fontStyle: FontStyle.italic) : null,
                decoration: InputDecoration(
                  hintText: 'Carangidae',
                  prefixIcon: Icon(FA.magnifyingGlass, size: 13, color: t.textFaint),
                  prefixIconConstraints: const BoxConstraints(minWidth: 38),
                ),
                onChanged: (_) => _onTaxonChanged(),
                onSubmitted: (_) => taxonOk && !c.busy ? _runSearch() : null,
              ),
            ),
          ),
          const SizedBox(width: 12),
          Expanded(
            flex: 3,
            child: Field(
              label: l10n.taxOrganelle,
              child: SplaceSelect<Organelle>(
                value: _organelle,
                options: [
                  SelectOption(Organelle.mitochondrion, l10n.organelleMito, detail: 'mtDNA'),
                  SelectOption(Organelle.chloroplast, l10n.organelleChloro, detail: 'cpDNA'),
                ],
                onChanged: (v) => setState(() {
                  _organelle = v;
                  _genes
                    ..clear()
                    ..addAll(_defaultGenes[v]!);
                }),
              ),
            ),
          ),
          const SizedBox(width: 12),
          Expanded(
            flex: 3,
            child: _markerMode
                ? Field(
                    label: l10n.taxScope,
                    info: l10n.markerGenomesInfo,
                    child: SwitchBox(label: l10n.markerIncludeGenomes, value: _includeGenomes, onChanged: (v) => setState(() => _includeGenomes = v)),
                  )
                : Field(
                    label: l10n.taxScope,
                    info: l10n.taxScopeInfo,
                    child: SwitchBox(label: l10n.taxCompleteOnly, value: _completeOnly, onChanged: (v) => setState(() => _completeOnly = v)),
                  ),
          ),
          const SizedBox(width: 12),
          Expanded(
            flex: 3,
            child: Field(
              label: l10n.taxDatabase,
              child: SwitchBox(label: l10n.taxRefseqOnly, value: _refseqOnly, onChanged: (v) => setState(() => _refseqOnly = v)),
            ),
          ),
        ]),
        // A species: it can also be searched under its other names.
        if (_isSpecies) ...[
          const SizedBox(height: 12),
          Container(
            padding: const EdgeInsets.all(12),
            decoration: BoxDecoration(color: t.surface2, borderRadius: BorderRadius.circular(10), border: Border.all(color: t.border)),
            child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
              Text.rich(TextSpan(style: TextStyle(fontSize: 13, color: t.text), children: [
                TextSpan(text: '${l10n.taxSpeciesRecognized} ', style: TextStyle(color: t.textMuted)),
                TextSpan(text: _taxon.text.trim(), style: const TextStyle(fontStyle: FontStyle.italic, fontWeight: FontWeight.w600)),
              ])),
              const SizedBox(height: 10),
              SwitchBox(
                label: l10n.taxSynonyms,
                detail: l10n.taxSynonymsInfo,
                value: _withSynonyms,
                onChanged: (v) {
                  setState(() => _withSynonyms = v);
                  if (v) _loadSynonyms();
                },
              ),
              if (_withSynonyms) ...[
                const SizedBox(height: 10),
                if (_synonymsLoading)
                  Row(children: [
                    SizedBox(width: 12, height: 12, child: CircularProgressIndicator(strokeWidth: 2, color: t.accent)),
                    const SizedBox(width: 8),
                    Text(l10n.taxSynonymsLoading, style: TextStyle(fontSize: 12.5, color: t.textMuted)),
                  ])
                else if (_synonymsError != null)
                  Text(l10n.taxSynonymsError(_synonymsError!), style: const TextStyle(fontSize: 12.5, color: SplaceTokens.danger))
                else if (_synonymsFor == _taxon.text.trim())
                  _synonymsInUse.isEmpty
                      ? Text(l10n.taxSynonymsNone, style: TextStyle(fontSize: 12.5, color: t.textMuted))
                      : Wrap(spacing: 6, runSpacing: 6, crossAxisAlignment: WrapCrossAlignment.center, children: [
                          Text(l10n.taxSynonymsFound('${_synonymsInUse.length}'), style: TextStyle(fontSize: 12.5, color: t.textMuted)),
                          for (final n in _synonymsInUse)
                            Container(
                              padding: const EdgeInsets.symmetric(horizontal: 8, vertical: 3),
                              decoration: BoxDecoration(color: t.accentSoft, borderRadius: BorderRadius.circular(8), border: Border.all(color: t.accentSoftBorder)),
                              child: Text(n, style: TextStyle(fontSize: 12, fontStyle: FontStyle.italic, color: t.accentText)),
                            ),
                        ]),
              ],
            ]),
          ),
        ],
        if (_markerMode) ...[
          const SizedBox(height: 12),
          Field(
            label: l10n.markerPick('${_genes.length}'),
            info: l10n.markerPickInfo,
            child: Wrap(spacing: 6, runSpacing: 6, children: [
              for (final g in syngenesMarkers(_type))
                ChipButton(
                  label: g,
                  style: _genes.contains(g) ? ChipButtonStyle.accent : ChipButtonStyle.neutral,
                  onPressed: () => setState(() => _genes.contains(g) ? _genes.remove(g) : _genes.add(g)),
                ),
            ]),
          ),
          const SizedBox(height: 12),
          // BOLD Systems: COI and the plant barcodes of the chloroplast.
          SwitchBox(
            label: l10n.markerBold,
            detail: boldServes(_type, _genes)
                ? l10n.markerBoldInfo([for (final g in syngenesMarkers(_type)) if (_genes.contains(g) && boldServes(_type, [g])) '$g (${boldMarkerCodes[g]!.join(', ')})'].join(', '))
                : l10n.markerBoldNone(_type == OrganelleType.cp ? 'rbcL, matK, rpoB, rpoC1, ycf1' : 'COI'),
            value: _bold && boldServes(_type, _genes),
            onChanged: boldServes(_type, _genes) ? (v) => setState(() => _bold = v) : null,
          ),
        ],
        const SizedBox(height: 12),
        PreviewBox(
          label: l10n.taxQueryPreview,
          text: _taxon.text.trim().isEmpty ? '' : (_markerMode ? _markerSearch.preview : _search.term),
          maxHeight: _markerMode ? 120 : 80,
        ),
        const SizedBox(height: 12),
        Row(children: [
          PrimaryButton(label: l10n.search, icon: FA.magnifyingGlass, onPressed: taxonOk && !c.busy ? _runSearch : null),
          const SizedBox(width: 12),
          if (c.ncbiStatus != null)
            Expanded(child: Text(c.ncbiStatus!, overflow: TextOverflow.ellipsis, style: TextStyle(fontSize: 13, color: t.textMuted))),
        ]),
      ]),
    );

    final option2 = OptionCard(
      badge: '1b',
      title: l10n.accTitle,
      info: l10n.accInfo,
      text: l10n.accBody,
      action: editableExample('1b', l10n.example, () => setState(() => _accessions.text = ncbiExampleAccessions), l10n.accExampleTooltip),
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        TextField(
          controller: _accessions,
          minLines: 4,
          maxLines: 8,
          style: const TextStyle(fontFamily: monoFont, fontSize: 13),
          decoration: const InputDecoration(hintText: 'NC_024026.1, KJ184305.1, NC_024184.1, KJ556976.1…'),
          onChanged: (_) => setState(() {}),
        ),
        const SizedBox(height: 12),
        Align(
          alignment: Alignment.centerLeft,
          child: SplaceTooltip(
            message: l10n.accFetchTooltip,
            child: PrimaryButton(
              label: l10n.fetch,
              icon: FA.download,
              onPressed: accessionOk && !c.busy ? () => c.fetchAccessions(_accessions.text) : null,
            ),
          ),
        ),
      ]),
    );

    final option3 = OptionCard(
      badge: '1c',
      title: l10n.gbFilesTitle,
      info: l10n.gbFilesInfo,
      text: l10n.gbFilesBody,
      action: dataExample(ExampleKind.genbank, '1c', l10n.example),
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        DropTarget(
          onDragEntered: (_) => setState(() => _dragGb = true),
          onDragExited: (_) => setState(() => _dragGb = false),
          onDragDone: (d) {
            setState(() => _dragGb = false);
            c.addPaths(d.files.map((f) => f.path).toList());
          },
          child: DropZone(title: l10n.gbFilesDrop, subtitle: l10n.gbFilesBrowse, highlight: _dragGb, onTap: browseGb),
        ),
        const SizedBox(height: 10),
        Align(
          alignment: Alignment.centerLeft,
          child: ChipButton(label: l10n.addFolder, icon: FA.folderOpen, onPressed: c.busy ? null : folderGb),
        ),
      ]),
    );

    final option4 = OptionCard(
      badge: '1d',
      title: l10n.filesTitle,
      info: l10n.filesInfo,
      text: l10n.filesBody,
      // Two kinds of FASTA example: loci and DNA barcodes.
      action: Wrap(spacing: 8, runSpacing: 6, children: [
        dataExample(ExampleKind.loci, '1d_loci', l10n.exampleLociShort),
        dataExample(ExampleKind.barcoding, '1d_barcoding', l10n.exampleBarcodingShort),
      ]),
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        DropTarget(
          onDragEntered: (_) => setState(() => _dragFa = true),
          onDragExited: (_) => setState(() => _dragFa = false),
          onDragDone: (d) {
            setState(() => _dragFa = false);
            c.addPaths(d.files.map((f) => f.path).toList());
          },
          child: DropZone(title: l10n.filesDrop, subtitle: l10n.filesBrowse, highlight: _dragFa, onTap: browseFa),
        ),
        const SizedBox(height: 10),
        Align(
          alignment: Alignment.centerLeft,
          child: ChipButton(label: l10n.addFolder, icon: FA.folderOpen, onPressed: c.busy ? null : folderFa),
        ),
      ]),
    );

    return StepSection(
      badge: '1',
      title: l10n.step1,
      description: l10n.step1Intro,
      active: widget.active,
      child: Column(crossAxisAlignment: CrossAxisAlignment.stretch, children: [
        // `.step1-tips-accordion`: intro, the "Example" button note and the
        // tips list with the red warning markers; key terms in bold.
        TipsAccordion(kicker: l10n.tipsKicker, title: l10n.tipsTitle, children: [
          MarkupText(l10n.tipsIntro),
          Padding(
            padding: const EdgeInsets.symmetric(vertical: 14),
            child: Text.rich(TextSpan(children: [
              TextSpan(text: '${l10n.tipsExamplePrefix} '),
              WidgetSpan(
                alignment: PlaceholderAlignment.middle,
                child: IgnorePointer(
                  child: ChipButton(label: l10n.example, icon: FA.wandMagicSparkles, style: ChipButtonStyle.soft, onPressed: () {}),
                ),
              ),
              TextSpan(text: ' ${l10n.tipsExampleSuffix}'),
            ])),
          ),
          for (final tip in [
            l10n.tipComplete,
            l10n.tipMultipleSearches,
            l10n.tipCombine,
            '**${l10n.tipSampleTitle}:** ${l10n.tipSampleBody}',
            '**${l10n.tipLocusTitle}:** ${l10n.tipLocusBody}',
            '**${l10n.tipMixedTitle}:** ${l10n.tipMixedBody}',
            '**${l10n.tipAlignedTitle}** ${l10n.tipAlignedBody}',
          ])
            Bullet(MarkupText(tip), icon: FA.triangleExclamation, gap: 10),
        ]),
        // One option per row, as in the HTML app.
        for (final (i, card) in [option1, option2, option3, option4].indexed) ...[
          if (i > 0) const SizedBox(height: 16),
          SubAnchor(id: '1${'abcd'[i]}', child: card),
        ],
        if (c.records.isNotEmpty || c.inputs.isNotEmpty) ...[
          const SizedBox(height: 16),
          TonePanel(
            tone: PanelTone.success,
            child: Row(children: [
              const Icon(FA.circleCheck, size: 15, color: SplaceTokens.success),
              const SizedBox(width: 10),
              Expanded(
                child: Text(
                  [
                    if (c.records.isNotEmpty) l10n.gbSummary('${c.records.length}'),
                    if (c.inputs.isNotEmpty) l10n.inputSummary('${c.inputs.length}', '$sequences'),
                  ].join(' · '),
                  style: const TextStyle(fontWeight: FontWeight.w600),
                ),
              ),
              if (c.inputs.isNotEmpty)
                ChipButton(label: l10n.clearFasta, icon: FA.trashCan, style: ChipButtonStyle.danger, onPressed: c.busy ? null : c.clearInputs),
            ]),
          ),
        ],
      ]),
    );
  }
}
