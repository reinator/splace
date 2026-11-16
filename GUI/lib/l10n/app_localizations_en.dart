// ignore: unused_import
import 'package:intl/intl.dart' as intl;

import 'app_localizations.dart';

// ignore_for_file: type=lint

/// The translations for English (`en`).
class AppLocalizationsEn extends AppLocalizations {
  AppLocalizationsEn([String locale = 'en']) : super(locale);

  @override
  String get navProgress => 'Progress';

  @override
  String get navSteps => 'Pipeline steps';

  @override
  String get logTitle => 'Activity log';

  @override
  String get logEmpty =>
      'No activity yet. Logs from fetching, alignment and analyses will appear here.';

  @override
  String get logCopy => 'Copy all';

  @override
  String get logDownload => 'Download as .txt';

  @override
  String get logClear => 'Clear';

  @override
  String get logCollapse => 'Hide log panel';

  @override
  String get logWelcome => 'SPLACE ready. Start by loading data in step 1.';

  @override
  String get settingsTheme => 'Theme';

  @override
  String get themeLight => 'Light';

  @override
  String get themeDark => 'Dark';

  @override
  String get themeAuto => 'Auto';

  @override
  String get settingsDensity => 'Density';

  @override
  String get densityComfy => 'Comfortable';

  @override
  String get densityCompact => 'Compact';

  @override
  String get settingsFontSize => 'Text size';

  @override
  String get fontSizeReset => 'Reset';

  @override
  String get minimapLocked => 'Complete the previous step to unlock';

  @override
  String get tooltipMenu => 'Toggle pipeline navigator';

  @override
  String get tooltipLogs => 'Activity log';

  @override
  String get tooltipSettings => 'Display & language';

  @override
  String get tooltipHelp => 'Show welcome screen';

  @override
  String get tooltipTools => 'External tools and analysis settings';

  @override
  String get cancel => 'Cancel';

  @override
  String get done => 'Done';

  @override
  String andMore(String n) {
    return '…and $n more';
  }

  @override
  String get splashSummary =>
      'A complete toolkit for phylogenomics: fetch GenBank records or organize FASTA loci, align, test saturation, concatenate and build trees.';

  @override
  String get splashStart => 'Start App';

  @override
  String get splashStep1 =>
      'Load genomes or individual markers (COI, 12S, CYTB…) from NCBI, or local GenBank and FASTA files. [[SynGenes]] builds the marker searches from every synonym of each gene and standardizes the marker names for the next steps.';

  @override
  String get splashStep2 =>
      'Inspect imported records, refine species names and recover taxonomy from GBIF, NCBI or SALVE with [[dataFishing]] before marker selection.';

  @override
  String get splashStep3 =>
      'Choose which feature types (CDS, rRNA, tRNA) enter the analysis, or show SPLACE where FASTA headers keep the sample and locus names.';

  @override
  String get splashStep4 =>
      'Pick the markers for the matrix and review missing and duplicated sequences in the gene presence and length heatmaps.';

  @override
  String get splashStep5 =>
      'Align each marker with [[MAFFT]] and trim it with [[trimAl]] or [[ClipKIT]], keeping whole codons in protein-coding loci.';

  @override
  String get splashStep6 =>
      'Test substitution saturation by codon position with the index of [[Iss (saturation)|Xia]] and decide which positions to keep.';

  @override
  String get splashStep7 =>
      'Concatenate the aligned loci into a supermatrix and export NEXUS, PHYLIP and partition files, optionally split by codon position.';

  @override
  String get splashStep8 =>
      'Estimate the maximum-likelihood phylogeny with [[IQ-TREE 3|IQ-TREE]] or [[RAxML-NG]], with branch support and optional per-gene trees.';

  @override
  String get step1 => 'Load Data';

  @override
  String get step2 => 'Loaded Records';

  @override
  String get step3 => 'Feature Types';

  @override
  String get step4 => 'Select Markers';

  @override
  String get step5 => 'Sequence Alignment';

  @override
  String get step6 => 'Saturation Test';

  @override
  String get step7 => 'Concatenate and NEXUS';

  @override
  String get step8 => 'Phylogenetic Inference';

  @override
  String get step1Intro =>
      'Choose how to load data into SPLACE: GenBank records from NCBI or local files, or your own FASTA files.';

  @override
  String get tipsKicker => 'Tips';

  @override
  String get tipsTitle => 'Before you start';

  @override
  String get tipsIntro =>
      'Build the **taxonomy table** gradually and combine **complementary sources** when needed.';

  @override
  String get tipsExamplePrefix => 'Each option below includes an';

  @override
  String get tipsExampleSuffix => 'button to help you get started quickly.';

  @override
  String get tipSampleTitle => 'One file per sample';

  @override
  String get tipSampleBody =>
      'each header is a **locus** (e.g. NC_030596_Oryza.fasta with >psbA, >matK). SPLACE regroups them into **one file per locus**.';

  @override
  String get tipLocusTitle => 'One file per locus';

  @override
  String get tipLocusBody =>
      'each header is a **sample** (e.g. COI.fasta). Already in the right shape.';

  @override
  String get tipMixedTitle => 'Mixed headers';

  @override
  String get tipMixedBody =>
      'sample and locus in the **same header** (>Bufo_bufo_COI) or **NCBI tags** ([gene=COX1]).';

  @override
  String get tipAlignedTitle => 'Already aligned?';

  @override
  String get tipAlignedBody =>
      'You can **skip MAFFT** in step 5 and go straight to trimming, saturation and trees.';

  @override
  String optionN(String n) {
    return 'Option $n';
  }

  @override
  String get filesTitle => 'Use local FASTA files';

  @override
  String get filesInfo =>
      'Import .fasta, .fas, .fa or .fna files that are already on your computer.';

  @override
  String get filesBody =>
      'Import your own .fasta, .fas, .fa or .fna files from disk.';

  @override
  String get filesDrop => 'Drag and drop FASTA files or folders here';

  @override
  String get filesBrowse => 'or click to browse (.fasta, .fas, .fa, .fna)';

  @override
  String get folderTitle => 'Import a whole folder';

  @override
  String get folderInfo =>
      'Reads every FASTA file inside the folder and its subfolders.';

  @override
  String get folderBody => 'Useful for hundreds of files, one per sample.';

  @override
  String get addFolder => 'Choose folder';

  @override
  String inputSummary(String files, String sequences) {
    return '$files files · $sequences sequences loaded';
  }

  @override
  String get clearRecords => 'Clear records';

  @override
  String get step2Intro =>
      'Every FASTA header is listed below with the sample and locus SPLACE reads from it. Adjust the rules in step 3 until every record has both.';

  @override
  String get detectedLayout => 'Detected layout';

  @override
  String get layoutSamplePerFile => 'one file per sample — headers are loci';

  @override
  String get layoutLocusPerFile => 'one file per locus — headers are samples';

  @override
  String get layoutMixedHeader => 'sample and locus are both in the header';

  @override
  String get layoutTaggedHeader => 'locus taken from header tags ([gene=…])';

  @override
  String get layoutUnknown => 'not recognized — choose the fields manually';

  @override
  String confidence(String percent) {
    return 'confidence $percent%';
  }

  @override
  String get organizing => 'Updating…';

  @override
  String get tokenColumns => 'Header parts (split on the locus delimiters)';

  @override
  String tokenColumnInfo(String distinct, String geneLike) {
    return '$distinct distinct · $geneLike% look like loci';
  }

  @override
  String get searchHeaders => 'Filter headers or files';

  @override
  String get colFile => 'File';

  @override
  String get colHeader => 'Header';

  @override
  String get colSample => 'Sample';

  @override
  String get colLocus => 'Locus';

  @override
  String get colLength => 'Length';

  @override
  String get step3Intro =>
      'Choose where the sample name and the locus name come from. The preview and the table in step 2 update as you type.';

  @override
  String fieldN(String n) {
    return 'Field $n';
  }

  @override
  String get sampleField => 'Sample name';

  @override
  String get sampleFieldInfo =>
      'Becomes the sequence name in every locus file, so it must be identical across loci for concatenation.';

  @override
  String get locusField => 'Locus name';

  @override
  String get locusFieldInfo =>
      'Defines which file each sequence goes to. Synonyms such as COX1/COI/CO1 can be merged.';

  @override
  String get sourceFileName => 'File name';

  @override
  String get sourceHeader => 'Header parts';

  @override
  String get sourceTag => 'Header tag';

  @override
  String get sourceRegex => 'Regex';

  @override
  String get delimiters => 'Delimiters';

  @override
  String get tokenStart => 'First part';

  @override
  String get tokenEnd => 'End (blank = last)';

  @override
  String get tagKey => 'Tag key';

  @override
  String get regexPattern => 'Pattern (first group is used)';

  @override
  String get regexOnFileName => 'Apply to file name';

  @override
  String get previewLabel => 'Preview';

  @override
  String get optionsKicker => 'Options';

  @override
  String get duplicates => 'Duplicates and names';

  @override
  String get duplicatesInfo =>
      'What to do when the same sample has the same locus more than once (e.g. inverted-repeat copies in plastomes).';

  @override
  String get dupLongest => 'Keep longest';

  @override
  String get dupMostComplete => 'Keep fewest N/gaps';

  @override
  String get dupFirst => 'Keep first';

  @override
  String get dupAll => 'Keep all (rename)';

  @override
  String get mergeSynonyms => 'Merge locus synonyms (COX1 = COI = CO1…)';

  @override
  String get sanitizeNames => 'Names safe for tree software';

  @override
  String get selectAll => 'Select All';

  @override
  String get selectNone => 'None';

  @override
  String get sharedOnly => 'Shared Only';

  @override
  String get legendSelectAll => 'selects every visible marker';

  @override
  String get legendNone => 'clears all selections';

  @override
  String get legendShared => 'markers present in all samples';

  @override
  String get missingTitle => 'Missing markers';

  @override
  String missingIn(String n) {
    return 'missing in $n sample(s)';
  }

  @override
  String get duplicatesTitle => 'Duplicate markers';

  @override
  String duplicatedIn(String n) {
    return 'duplicated in $n sample(s)';
  }

  @override
  String get colCopies => 'Copies (sizes)';

  @override
  String get colUsed => 'Used';

  @override
  String get unassignedTitle => 'Records without sample or locus';

  @override
  String get missingSample => 'no sample';

  @override
  String get missingLocus => 'no locus';

  @override
  String useLoci(String n) {
    return 'Use $n markers';
  }

  @override
  String get exportLoci => 'Export FASTA per locus…';

  @override
  String get heatmapTitle => 'Gene Presence Heatmap';

  @override
  String get heatmapAbsent => 'Absent';

  @override
  String get heatmapPresent => 'Present';

  @override
  String matrixSummary(String loci, String samples, String occupancy) {
    return '$loci markers × $samples samples · $occupancy% filled. Click a marker name to include or exclude it.';
  }

  @override
  String get panelMafft => 'MAFFT Settings';

  @override
  String get panelMafftBody =>
      'Each marker is aligned separately; several markers run in parallel.';

  @override
  String get mafftStrategy => 'Alignment Method';

  @override
  String get directionLabel => 'Sequence direction';

  @override
  String get adjustDirection => 'Fix reversed sequences';

  @override
  String get panelTrimal => 'Trimming (trimAl / ClipKIT)';

  @override
  String get trimalMethod => 'Trimming Mode';

  @override
  String get sequenceType => 'Sequence Type';

  @override
  String get codonAwareTrimming => 'Keep the reading frame';

  @override
  String get codonAwareHint =>
      'Removing single columns from protein-coding alignments shifts the reading frame and invalidates codon partitions and the saturation test.';

  @override
  String get panelRun => 'Run Alignment';

  @override
  String panelRunBody(String dir) {
    return 'Align all selected markers with MAFFT. Output is saved to $dir.';
  }

  @override
  String get runMafft => 'Align (MAFFT)';

  @override
  String get runTrimal => 'Trim (trimAl)';

  @override
  String get useAsAlignment => 'Input is already aligned';

  @override
  String get colCoding => 'Coding';

  @override
  String get colFrame => 'Reading frame';

  @override
  String get colInclude => 'Use this marker';

  @override
  String get statRaw => 'Raw';

  @override
  String get statAligned => 'Aligned';

  @override
  String get statTrimmed => 'Trimmed';

  @override
  String get statusError => 'Error';

  @override
  String get saturationIntro =>
      'Each marker is tested with PhyKIT (Philippe et al. 2011): IQ-TREE infers the tree of the gene, then, for the 1st, 2nd and 3rd codon positions (and 1st+2nd, all sites), PhyKIT regresses the uncorrected p-distance of every pair of taxa on their distance in the tree. A value near 1 means no saturation; the closer to 0, the more the observed differences hide repeated substitutions.';

  @override
  String get saturationThresholds =>
      'Heuristic thresholds, set in 6a (test parameters): moderate below the moderate limit of PhyKIT\'s saturation value, saturated below the saturated limit; too little divergence (mean p) counts as OK.';

  @override
  String get runSaturation => 'Run saturation test';

  @override
  String get recommendation => 'Recommendation';

  @override
  String get recKeepAll =>
      'Keep all three codon positions. In step 7b, partition the matrix by codon position so each one gets its own substitution model.';

  @override
  String get recAminoAcids =>
      'All three positions are saturated: at this divergence the nucleotides say little about the relationships. Align as amino acids (step 5b) or use slower markers.';

  @override
  String get recUndetermined =>
      'There are too few sites or sequence pairs to judge saturation. Keep all positions.';

  @override
  String get applyPerLocus => 'Use each gene\'s own result';

  @override
  String get concatenated => 'Concatenated';

  @override
  String get saturationHeatmapTitle => 'Saturation by locus and codon position';

  @override
  String get selectLocusHint =>
      'Hover a cell to see its values; click a row to see its plots.';

  @override
  String get verdictOk => 'OK';

  @override
  String get verdictModerate => 'Moderate';

  @override
  String get verdictSaturated => 'Saturated';

  @override
  String get verdictInsufficient => 'Too little data';

  @override
  String plotObservedVsCorrected(String model) {
    return 'Observed p vs $model distance';
  }

  @override
  String plotTsTv(String model) {
    return 'Transitions and transversions vs $model distance';
  }

  @override
  String get transitions => 'Transitions';

  @override
  String get transversions => 'Transversions';

  @override
  String get metricSites => 'Sites';

  @override
  String get metricPairs => 'Pairs';

  @override
  String get metricMeanP => 'Mean p';

  @override
  String get metricSlope => 'Slope';

  @override
  String get metricTsTv => 'Ts/Tv';

  @override
  String get metricIss => 'Iss';

  @override
  String get metricUndefined => 'Uncorrectable';

  @override
  String get metricVerdict => 'Verdict';

  @override
  String frameLabel(String frame) {
    return 'codons start at base $frame';
  }

  @override
  String get step7Intro =>
      'Build a supermatrix from aligned genes and export for phylogenetic analysis.';

  @override
  String get srcTrimmed => 'Trimmed sequences';

  @override
  String get recommendedTag => '(recommended)';

  @override
  String get srcTrimmedBody =>
      'Use trimAl output — removes poorly aligned positions.';

  @override
  String get srcAligned => 'Aligned sequences';

  @override
  String get srcAlignedBody => 'Use MAFFT output without trimming.';

  @override
  String get partitionScheme => 'Partitions';

  @override
  String get schemeNone => 'None';

  @override
  String get schemeByLocus => 'By locus';

  @override
  String get schemeByLocusCodon => 'By locus and codon position';

  @override
  String get schemeByCodon => 'By codon position (all loci)';

  @override
  String get codonExclusion => 'Codon positions';

  @override
  String get codonExclusionInfo =>
      'Excluded positions are removed from coding loci only. The saturation test (step 6) can fill this in.';

  @override
  String excludePosition(String n) {
    return 'Exclude position $n';
  }

  @override
  String perLocusExclusions(String list) {
    return 'Per locus: $list';
  }

  @override
  String get clearExclusions => 'Clear per-locus exclusions';

  @override
  String get buildMatrix => 'Build supermatrix';

  @override
  String get matrixReady => 'Supermatrix ready';

  @override
  String matrixStats(
    String taxa,
    String sites,
    String partitions,
    String missing,
  ) {
    return '$taxa taxa × $sites sites · $partitions partitions · $missing% missing';
  }

  @override
  String get excludedLabel => 'Excluded';

  @override
  String get outputFiles => 'Output files';

  @override
  String get openFolder => 'Open folder';

  @override
  String get partitionPreview => 'RAxML-NG partitions';

  @override
  String get step8Intro =>
      'Maximum-likelihood tree estimation from the concatenated supermatrix with IQ-TREE3 or RAxML-NG.';

  @override
  String get outgroup => 'Outgroup';

  @override
  String get outgroupInfo =>
      'Selected species are written to the NEXUS file and used to root IQ-TREE and RAxML-NG trees.';

  @override
  String get outgroupHint => 'Search species…';

  @override
  String engineN(String n) {
    return 'Engine $n';
  }

  @override
  String get iqtreeInfo =>
      'ModelFinder chooses the models (and can merge partitions); ultrafast bootstrap and SH-aLRT give support values.';

  @override
  String get raxmlInfo =>
      'RAxML-NG 2 chooses the models itself (MOOSE), adapts the search to how hard the alignment is and offers fast support measures. Like IQ-TREE it is a maximum-likelihood program (not Bayesian).';

  @override
  String get model => 'Substitution Model';

  @override
  String get ufboot => 'Ultrafast bootstrap (-B)';

  @override
  String get alrt => 'SH-aLRT (--alrt)';

  @override
  String get runIqtree => 'Run IQ-TREE3';

  @override
  String get raxmlModelFinder =>
      'ModelFinder (IQ-TREE, RAxML-compatible models)';

  @override
  String get bootstraps => 'Bootstrap trees';

  @override
  String get supportMetrics => 'Support';

  @override
  String get seed => 'Seed';

  @override
  String get runRaxml => 'Run RAxML-NG';

  @override
  String get runs => 'Runs';

  @override
  String get noRuns => 'No trees yet.';

  @override
  String toolMissing(String tool) {
    return '$tool not found';
  }

  @override
  String get toolsTitle => 'External tools and analysis settings';

  @override
  String get toolsBody =>
      'SPLACE looks for bundled binaries, then PATH, then WSL. Set a custom path if needed.';

  @override
  String get tools => 'Tools';

  @override
  String get redetect => 'Detect again';

  @override
  String get notFound => 'Not found';

  @override
  String get overridePath => 'Custom path (or wsl:/path/inside/wsl)';

  @override
  String get browse => 'Browse…';

  @override
  String get analysisSettings => 'Analysis';

  @override
  String threads(String n) {
    return 'Threads: $n';
  }

  @override
  String get geneticCode => 'Genetic code';

  @override
  String get projectFolder => 'Project folder';

  @override
  String get change => 'Change…';

  @override
  String get tipComplete =>
      'Prefer **complete genomes** whenever possible because they usually improve downstream analyses.';

  @override
  String get tipMultipleSearches =>
      'Run **multiple taxonomic searches** to add **ingroups** and **outgroups**. New results are appended to the same taxonomy table.';

  @override
  String get tipCombine =>
      'You can **combine** taxonomic searches, accession imports, and local **.gb, .gbk or .genbank** files in the same dataset.';

  @override
  String get example => 'Example';

  @override
  String get taxTitle => 'Search by taxonomic name';

  @override
  String get taxInfo =>
      'Search NCBI by taxonomic name, genome type, database and completeness before importing the chosen genomes.';

  @override
  String get taxBody =>
      'Search NCBI and download matching GenBank records automatically into the tool.';

  @override
  String get taxExampleTooltip =>
      'Fill the taxonomic search with Carangidae, mitochondrion, complete genomes and RefSeq only.';

  @override
  String get taxApiHint =>
      'Add an NCBI API key to speed up imports (key button in the top bar).';

  @override
  String get taxName => 'Taxonomic name';

  @override
  String get taxOrganelle => 'Genome type';

  @override
  String get organelleMito => 'Mitochondrion';

  @override
  String get organelleChloro => 'Chloroplast';

  @override
  String get taxScope => 'Genome scope';

  @override
  String get taxCompleteOnly => 'Complete only';

  @override
  String get taxScopeInfo =>
      'Complete only searches NCBI for records with “complete genome” in the title. When off, records with “partial genome” in the title are included too.';

  @override
  String get scopeComplete => 'Complete genomes';

  @override
  String get scopePartial => 'Partial genomes';

  @override
  String get taxDatabase => 'Database';

  @override
  String get taxRefseqOnly => 'RefSeq only';

  @override
  String get taxQueryPreview => 'Query preview';

  @override
  String get search => 'Search';

  @override
  String get accTitle => 'Paste accession numbers';

  @override
  String get accInfo =>
      'Paste one or more known accessions and fetch those exact GenBank records from NCBI.';

  @override
  String get accBody =>
      'Enter accession numbers separated by commas, spaces or line breaks.';

  @override
  String get accExampleTooltip =>
      'Load example accessions from Balcázar et al. (2026)';

  @override
  String get accFetchTooltip =>
      'Fetch GenBank records from NCBI for the entered accession numbers';

  @override
  String get fetch => 'Fetch';

  @override
  String get gbFilesTitle => 'Use local GenBank files';

  @override
  String get gbFilesInfo =>
      'Import local .gb, .gbk or .genbank files when the records are already on your computer.';

  @override
  String get gbFilesBody =>
      'Import your own .gb, .gbk or .genbank files from disk.';

  @override
  String get gbFilesDrop => 'Drag and drop GenBank files or folders here';

  @override
  String get gbFilesBrowse => 'or click to browse (.gb, .gbk, .genbank)';

  @override
  String gbSummary(String n) {
    return '$n GenBank records loaded';
  }

  @override
  String get clearFasta => 'Clear FASTA files';

  @override
  String get close => 'Close';

  @override
  String get apiKeyTitle => 'NCBI API key';

  @override
  String get apiKeyHelper =>
      'Used by Option 1 and Option 2 for NCBI requests. Without an API key, SPLACE runs at about 3 requests per second. With a saved valid key, SPLACE uses up to 10 requests per second.';

  @override
  String get apiKeySave => 'Save key';

  @override
  String get apiKeyClear => 'Clear key';

  @override
  String get apiKeyGetLink => 'Get an API key from your NCBI account settings';

  @override
  String get apiKeySaved =>
      'API key saved. Option 1 and Option 2 now use up to 10 requests per second.';

  @override
  String get apiKeyInvalid =>
      'The API key looks invalid. Check it before saving.';

  @override
  String get apiKeyNone =>
      'No API key saved. Option 1 and Option 2 are using up to 3 requests per second.';

  @override
  String get apiClearTitle => 'Remove saved API key?';

  @override
  String get apiClearBody =>
      'Option 1 and Option 2 will go back to the default NCBI request rate.';

  @override
  String get reviewTitle => 'Review genomes before importing';

  @override
  String get reviewSubtitle =>
      'Choose which genomes should be imported. All switches start enabled.';

  @override
  String get reviewDetailTaxon => 'Taxonomic query';

  @override
  String get reviewDetailGenomeType => 'Genome type';

  @override
  String get reviewDetailSearchTypes => 'Search types';

  @override
  String get reviewDetailRefseq => 'RefSeq';

  @override
  String get reviewRefseqOn => 'Used';

  @override
  String get reviewRefseqOff => 'Not used';

  @override
  String get reviewSelectAll => 'Enable all';

  @override
  String get reviewSelectNone => 'Disable all';

  @override
  String get reviewDisableUnverified => 'Disable UNVERIFIED';

  @override
  String get reviewOnlyComplete => 'Only complete';

  @override
  String get reviewOnlyPartial => 'Only partial';

  @override
  String get reviewColImport => 'Import';

  @override
  String get reviewColAccession => 'Accession';

  @override
  String get reviewColDescription => 'Genome description';

  @override
  String get reviewColScope => 'Classification';

  @override
  String get reviewImport => 'Import selected genomes';

  @override
  String get reviewTotal => 'Total';

  @override
  String get reviewSelected => 'Selected';

  @override
  String get reviewComplete => 'Complete';

  @override
  String get reviewPartial => 'Partial';

  @override
  String get reviewMixed => 'Mixed title';

  @override
  String get reviewUnknown => 'Unclassified';

  @override
  String get reviewUnverified => 'UNVERIFIED';

  @override
  String reviewUnverifiedCount(String count) {
    return 'UNVERIFIED genomes found: $count';
  }

  @override
  String get reviewUnverifiedAdvice =>
      'We recommend removing UNVERIFIED genomes before importing.';

  @override
  String get reviewCancelled => 'Genome selection cancelled.';

  @override
  String get reviewSearch => 'Filter by accession, species or description';

  @override
  String reviewShowing(String shown, String total) {
    return 'Showing $shown/$total genomes matching the search.';
  }

  @override
  String get gbRecordsHeading => 'GenBank records';

  @override
  String get fastaHeadersHeading => 'FASTA headers';

  @override
  String get clearAll => 'Clear all';

  @override
  String get clearTitle => 'Clear all records?';

  @override
  String get clearBody =>
      'This action cannot be undone. All loaded records and selections will be removed.';

  @override
  String get recordsReview =>
      'Review and correct the species names below using the edit button.';

  @override
  String get recordsReviewMethod =>
      'Review and correct the species names below using the edit button, then fetch the taxonomy with the method of [[dataFishing]].';

  @override
  String get recordsFetchPrefix => 'Then click';

  @override
  String get recordsFetchSuffix =>
      'to retrieve classification ranks from GBIF via';

  @override
  String get internetRequired => 'Internet required';

  @override
  String get fetchTaxonomy => 'Fetch Taxonomy';

  @override
  String get skipTaxonomy => 'Skip taxonomy and proceed';

  @override
  String get skipTaxBody1 =>
      'If species names contain special characters or there are duplicate entries, errors may occur in future analyses.';

  @override
  String get skipTaxBody2 =>
      'It is recommended to correct species names and fetch taxonomy for best results.';

  @override
  String get taxonomyResolved => 'Taxonomy ready';

  @override
  String get taxonomyTree => 'Taxonomic tree';

  @override
  String get taxonomyTreeBody =>
      'Indented taxonomy summary using the same records shown below.';

  @override
  String get columnsLabel => 'Columns:';

  @override
  String get colAccession => 'Accession';

  @override
  String get colKingdom => 'Kingdom';

  @override
  String get colPhylum => 'Phylum';

  @override
  String get colClass => 'Class';

  @override
  String get colOrder => 'Order';

  @override
  String get colFamily => 'Family';

  @override
  String get colGenus => 'Genus';

  @override
  String get colSpecies => 'Species';

  @override
  String get colAuthorship => 'Authorship';

  @override
  String get colSource => 'Source';

  @override
  String get sourceFile => 'File';

  @override
  String get recordsSearch => 'Search accession, species, family...';

  @override
  String editRecordTooltip(String species, String accession) {
    return 'Edit $species ($accession)';
  }

  @override
  String removeRecordTooltip(String species, String accession) {
    return 'Remove $species ($accession)';
  }

  @override
  String get removeTitle => 'Remove record?';

  @override
  String get removeBody => 'This action cannot be undone.';

  @override
  String get remove => 'Remove';

  @override
  String get editTitle => 'Edit Record';

  @override
  String get editHint =>
      'Review the current record and update species and taxonomy fields below.';

  @override
  String get editGbif => 'Fetch Taxonomy from GBIF';

  @override
  String get save => 'Save';

  @override
  String get autoLabel => '(auto)';

  @override
  String get fileSize => 'Size';

  @override
  String get sequenceLabel => 'Sequence';

  @override
  String get featuresLabel => 'Features';

  @override
  String get lineageLabel => 'Lineage';

  @override
  String get speciesRequired => 'Species name is required';

  @override
  String get speciesLetters => 'Only letters, spaces and hyphens are allowed';

  @override
  String get speciesBinomial => 'Use genus and epithet (e.g. Rhinella marina)';

  @override
  String get speciesGenusUpper =>
      'The genus must start with an uppercase letter';

  @override
  String get familySuffix => 'Family names end in -idae or -aceae';

  @override
  String get gbifNoMatch => 'No GBIF match for this name.';

  @override
  String get gbifMatched => 'GBIF match found.';

  @override
  String get step3Fasta => 'Sample and Locus Fields';

  @override
  String get featureTypesIntro =>
      'Select the gene categories you want to include in your analysis. Toggle each category to activate or deactivate it.';

  @override
  String get standardizedPrefix => 'Gene names have been standardized using';

  @override
  String get featureCdsInfo =>
      'Protein-coding sequences (genes translated into proteins)';

  @override
  String get featureRrnaInfo => 'Ribosomal RNA genes (12S, 16S, etc.)';

  @override
  String get featureTrnaInfo => 'Transfer RNA genes';

  @override
  String get featureTypesHint =>
      'Turn on at least one category to list the markers in step 4.';

  @override
  String unknownGenesButton(String n) {
    return '$n unknown names';
  }

  @override
  String get unknownTitle =>
      'Gene Names Not Found in SynGenes (Rabelo et al. 2024)';

  @override
  String get unknownSubtitle =>
      'Select the correct standardized name for each gene found in your files.';

  @override
  String unknownCount(String count) {
    return '$count gene(s) could not be automatically standardized.';
  }

  @override
  String get unknownHowto =>
      'Select the standardized name from the dropdown in the Correct Name column. Choose (skip) to exclude the gene from analysis.';

  @override
  String get unknownSkip => '(skip)';

  @override
  String geneSuggested(String name) {
    return 'Suggested: $name';
  }

  @override
  String geneFromProduct(String value) {
    return 'from /product “$value”';
  }

  @override
  String geneFromGene(String value) {
    return 'from /gene “$value”';
  }

  @override
  String geneFromSequence(String percent, String name, String n) {
    return 'sequence $percent% similar to $name in $n other genomes';
  }

  @override
  String geneSequenceAgrees(String percent) {
    return 'the sequence agrees ($percent%)';
  }

  @override
  String geneSequenceDisagrees(String name, String percent) {
    return 'but the sequence looks like $name ($percent%)';
  }

  @override
  String get geneNoSuggestion =>
      'No match in SynGenes or in the other genomes; check the annotation.';

  @override
  String geneQualifierMismatch(String gene, String product) {
    return '/gene “$gene” and /product “$product” differ';
  }

  @override
  String conflictsBanner(String n) {
    return '$n feature(s) have /gene and /product naming different genes.';
  }

  @override
  String get conflictsReview => 'Review';

  @override
  String get conflictsTitle => '/gene and /product disagree';

  @override
  String get conflictsSubtitle =>
      'SPLACE uses the /product name. Each sequence is compared with the same marker in the other genomes to show which qualifier is right.';

  @override
  String conflictUsed(String name) {
    return 'Used: $name';
  }

  @override
  String conflictSeqGene(String name, String percent) {
    return 'Sequence matches $name ($percent%): /gene looks right';
  }

  @override
  String conflictSeqUsed(String name, String percent) {
    return 'Sequence matches $name ($percent%): the name used looks right';
  }

  @override
  String conflictSeqOther(String name, String percent) {
    return 'Sequence matches $name ($percent%), neither qualifier';
  }

  @override
  String get conflictSeqNone => 'Sequence inconclusive';

  @override
  String conflictUseGene(String name) {
    return 'Use $name';
  }

  @override
  String get speciesIdRunning => 'Identifying species';

  @override
  String speciesIdGenesUsed(String genes) {
    return 'GenBank BLAST and BOLD · $genes';
  }

  @override
  String get speciesIdFooter =>
      'BLAST at NCBI usually takes 1–3 minutes; BOLD answers in seconds.';

  @override
  String get speciesIdButton => 'Identify species (GenBank + BOLD)';

  @override
  String get speciesIdShowResults => 'Identification results';

  @override
  String get speciesIdSetupTitle => 'Identify species with GenBank and BOLD';

  @override
  String get speciesIdSetupSubtitle =>
      'Compare the genes of each genome with GenBank (BLAST) and BOLD to check the species names before fetching the taxonomy.';

  @override
  String get speciesIdGenesLabel => 'Genes to compare';

  @override
  String get speciesIdGenesHint =>
      'Choose one or more genes. Every gene is searched in GenBank with BLAST; BOLD accepts COI only.';

  @override
  String speciesIdInGenomes(String n) {
    return 'in $n genomes';
  }

  @override
  String get speciesIdSources => 'Databases';

  @override
  String get speciesIdGenbankBody =>
      'BLAST (megablast) against core_nt, for every selected gene.';

  @override
  String get speciesIdBoldBody =>
      'BOLD ID engine, public species library, COI only.';

  @override
  String get speciesIdBoldOff => 'Select COI to also search BOLD.';

  @override
  String get speciesIdRun => 'Identify';

  @override
  String get speciesIdResultsTitle => 'Species identification';

  @override
  String speciesIdResultsSubtitle(String genes) {
    return 'GenBank BLAST and BOLD · $genes. Choose the name to keep for each genome.';
  }

  @override
  String get speciesIdOnlyDisagree => 'Only possible misidentifications';

  @override
  String speciesIdSummary(String disagree, String total) {
    return '$disagree of $total genomes may be misidentified';
  }

  @override
  String get speciesIdColGenome => 'Genome';

  @override
  String get speciesIdColCandidates => 'Species found (choose one)';

  @override
  String get speciesIdMismatch => 'Possible misidentification';

  @override
  String get speciesIdMatch => 'Confirmed';

  @override
  String get speciesIdNoHits => 'No match in GenBank or BOLD';

  @override
  String speciesIdKeepCurrent(String name) {
    return 'Keep $name';
  }

  @override
  String speciesIdHits(String n) {
    return '$n hits';
  }

  @override
  String speciesIdApply(String n) {
    return 'Apply $n change(s)';
  }

  @override
  String get speciesIdNote =>
      'Hits on the genome itself (or on the GenBank record its RefSeq NC_ copy came from) appear as “this record” and do not count as evidence.';

  @override
  String get geneReviewButton => 'Review gene names';

  @override
  String geneReviewButtonN(String n) {
    return 'Review gene names ($n)';
  }

  @override
  String get geneReviewDone => 'Gene names reviewed';

  @override
  String get geneReviewTitle => 'Review gene names';

  @override
  String get geneReviewSubtitle =>
      'Names SynGenes did not recognize and genes whose /gene and /product disagree. Review them before the taxonomy.';

  @override
  String geneReviewUnknownSection(String n) {
    return 'Names not found in SynGenes ($n)';
  }

  @override
  String geneReviewConflictSection(String n) {
    return '/gene and /product name different genes ($n)';
  }

  @override
  String get geneReviewConflictHelp =>
      'Both qualifiers were normalized with SynGenes. Choose which name this feature gets; only this genome changes.';

  @override
  String get geneReviewNothing =>
      'All gene names were standardized and every /gene agrees with its /product.';

  @override
  String conflictPickGene(String name) {
    return '/gene → $name';
  }

  @override
  String conflictPickProduct(String name) {
    return '/product → $name';
  }

  @override
  String get grHowTitle => 'How to read this review';

  @override
  String get grHowSynGenes =>
      'SynGenes recognizes a tag when it is a standard gene name or one of its known synonyms.';

  @override
  String get grHowSequence =>
      'Similarity is the share of the sequence (protein for CDS, DNA for rRNA) also found in the same gene of the other genomes in this dataset: above 50% is strong, 35–50% moderate, below 35% weak.';

  @override
  String get grSequence => 'Sequence';

  @override
  String grRecognized(String name) {
    return 'SynGenes: $name';
  }

  @override
  String get grNotRecognized => 'not in SynGenes';

  @override
  String get grAbsent => 'absent';

  @override
  String grInGenomes(String n) {
    return 'in $n genome(s)';
  }

  @override
  String get grUseName => 'Name to use';

  @override
  String get grNoSequenceRef => 'no other genome to compare';

  @override
  String grWhyProduct(String value, String name) {
    return 'The /product “$value” is a SynGenes name of $name.';
  }

  @override
  String grWhyGene(String value, String name) {
    return 'The /gene “$value” is a SynGenes name of $name.';
  }

  @override
  String grWhySequence(String percent, String name, String n) {
    return 'Neither tag is in SynGenes, but the sequence is $percent% similar to $name in $n other genomes.';
  }

  @override
  String get grWhyNone =>
      'Neither tag is in SynGenes and the sequence does not clearly match any marker. Check the annotation or skip this name.';

  @override
  String grSeqAgrees(String percent) {
    return 'The sequence agrees ($percent% similar).';
  }

  @override
  String grSeqDisagrees(String name, String percent) {
    return 'Careful: the sequence looks more like $name ($percent%).';
  }

  @override
  String grSays(String tag, String name) {
    return '$tag says $name';
  }

  @override
  String grSimilarityTo(String name) {
    return 'Similarity to $name in the other genomes';
  }

  @override
  String get grSupported => 'Supported by the sequence';

  @override
  String get grNotSupported => 'Not supported by the sequence';

  @override
  String grConflictGeneRight(
    String gp,
    String gene,
    String pp,
    String product,
  ) {
    return 'The sequence is $gp% similar to $gene and $pp% to $product in the other genomes, so the /gene tag ($gene) is probably right and the /product tag is wrong.';
  }

  @override
  String grConflictProductRight(
    String pp,
    String product,
    String gp,
    String gene,
  ) {
    return 'The sequence is $pp% similar to $product and $gp% to $gene in the other genomes, so the /product tag ($product) is probably right and the /gene tag is wrong.';
  }

  @override
  String get grConflictUnclear =>
      'The sequence does not clearly favour either gene. Check the annotation of this genome before choosing.';

  @override
  String get grConflictIntro =>
      'The /gene and /product tags of this feature name different genes. Click the name to use; only this genome changes.';

  @override
  String get siHowTitle => 'How to read the identification';

  @override
  String get siHowBody =>
      'The genes of each genome were compared with GenBank (BLAST) and BOLD. Identity is the share of identical bases with the closest sequence found.';

  @override
  String get siLegendSame => '≥ 98%: same species';

  @override
  String get siLegendClose => '95–98%: close relative';

  @override
  String get siLegendOther => '< 95%: another species';

  @override
  String get siColUse => 'Use';

  @override
  String get siColSpecies => 'Species';

  @override
  String get siColGenbank => 'GenBank (BLAST)';

  @override
  String get siColBold => 'BOLD (COI)';

  @override
  String get siColSeqs => 'Matching sequences';

  @override
  String get siCurrentTag => 'current name';

  @override
  String get siNotFound => 'not found';

  @override
  String siVerdictConfirmed(String identity) {
    return 'The matches confirm the current name (best identity $identity%).';
  }

  @override
  String siVerdictMismatch(
    String best,
    String identity,
    String current,
    String currentIdentity,
  ) {
    return 'The closest match is $best ($identity%), while the current name $current reaches $currentIdentity%. The genome may be misidentified.';
  }

  @override
  String siVerdictMissing(String best, String identity, String current) {
    return 'The closest match is $best ($identity%) and $current is not among the matches. The genome may be misidentified, or the name may be a synonym.';
  }

  @override
  String get siVerdictNone =>
      'No sequence in GenBank or BOLD matched the genes of this genome.';

  @override
  String noMarkersButton(String n) {
    return 'Remove genomes without markers ($n)';
  }

  @override
  String get noMarkersTitle => 'Genomes without markers';

  @override
  String get noMarkersSubtitle =>
      'These genomes have no protein-coding or rRNA gene that SPLACE can use, so they add nothing to the matrix. Check the list before removing them.';

  @override
  String get noMarkersReasonNone =>
      'No annotated genes (0 PCGs, 0 rRNAs, 0 tRNAs): the GenBank record has no gene annotation.';

  @override
  String noMarkersReasonTrna(String n) {
    return 'Only $n tRNA(s): no PCGs or rRNAs.';
  }

  @override
  String noMarkersReasonUnnamed(String pcgs, String rrnas) {
    return '$pcgs PCG(s) and $rrnas rRNA(s) annotated, but none has a /gene or /product name to identify it.';
  }

  @override
  String noMarkersRemove(String n) {
    return 'Remove $n genome(s)';
  }

  @override
  String get heatmapShow => 'Show heatmap';

  @override
  String get heatmapHide => 'Hide heatmap';

  @override
  String get heatmapHoverHint =>
      'Point at a cell, a marker or a sample to see its details.';

  @override
  String get headerPresets => 'Presets';

  @override
  String get headerDragHint =>
      'Drag to reorder. The example comes from the first record.';

  @override
  String get headerAddHint =>
      'Click a field to add it to the end of the header.';

  @override
  String get headerFirstRecord => 'First record';

  @override
  String get headerFirstRecords => 'Headers of the first records';

  @override
  String siSourceFailed(String source, String error) {
    return '$source did not answer, so its column has no results: $error';
  }

  @override
  String get whenToUse => 'When to use';

  @override
  String get trimCodonHandling => 'Protein-coding genes';

  @override
  String get trimCodonHandlingInfo =>
      'Codon positions are used later for the partitions and the saturation test (step 6), so coding genes should keep whole codons. trimAl does it with -backtrans (it trims the translated alignment and writes back the codons; -ignorestopcodon because it checks stops with the standard code) and ClipKIT with --codon. Excluding a codon position from the supermatrix is decided in step 7b, not here.';

  @override
  String get alignRunning => 'Aligning with MAFFT';

  @override
  String get trimRunning => 'Trimming alignments';

  @override
  String markersProgress(String done, String total, String mode) {
    return '$done of $total markers · $mode';
  }

  @override
  String get alignFooter =>
      'MAFFT runs on this computer, several markers at a time.';

  @override
  String get trimFooter =>
      'trimAl/ClipKIT run on this computer, several markers at a time.';

  @override
  String alignCurrent(String markers) {
    return 'Running: $markers';
  }

  @override
  String get colReason => 'Why';

  @override
  String get unknownSkipAll => 'Skip All';

  @override
  String get unknownApply => 'Apply Corrections';

  @override
  String get more => 'more';

  @override
  String get defaultPhylogeny => 'Default to Phylogeny';

  @override
  String get legendDefault => 'recommended set for phylogenetic analysis';

  @override
  String get tooltipAll => 'Select every visible marker';

  @override
  String get tooltipNone => 'Clear every selected marker';

  @override
  String get tooltipShared =>
      'Select only genes present in all genomes (no missing data)';

  @override
  String get tooltipDefault =>
      'Select the default gene set for phylogenetic analysis';

  @override
  String get duplicateChoiceHelp =>
      'The supermatrix uses one copy per genome: click the copy to use. The longest starts selected; “closest” marks the copy whose length is nearest the genomes with a single copy, usually the right one when the copies differ.';

  @override
  String get defaultLabel => 'default';

  @override
  String dupReference(
    String gene,
    String n,
    String avg,
    String min,
    String max,
  ) {
    return '$gene in $n genomes with a single copy: average $avg bp (from $min to $max bp). Pick the copy closest to it.';
  }

  @override
  String get dupNoReference =>
      'No other genome has a single copy of this gene to compare lengths.';

  @override
  String get dupCopiesColumn => 'Copies (click the one to use)';

  @override
  String dupVsAverage(String diff) {
    return '$diff bp vs average';
  }

  @override
  String get dupSameAsAverage => 'same as average';

  @override
  String get dupClosest => 'closest to average';

  @override
  String get dupLongestTag => 'longest';

  @override
  String get noMarkers => 'No markers found for the selected feature types.';

  @override
  String get featureLabel => 'Feature';

  @override
  String get recordsLabel => 'Records';

  @override
  String heatmapSummary(String markers, String samples) {
    return '$markers markers × $samples samples. Click a marker name to select it.';
  }

  @override
  String get lengthHeatmapTitle => 'Marker Length Heatmap';

  @override
  String get lengthHeatmapBody =>
      'Each selected marker is scaled independently by its own minimum and maximum length.';

  @override
  String get panelHeader => 'Configure FASTA Header';

  @override
  String get panelHeaderBody =>
      'Sequence identifier of the alignment files. The default (accession_family_genus_species) is used unless you change it in the Header Builder.';

  @override
  String get openHeaderBuilder => 'Open Header Builder';

  @override
  String get headerBuilderTitle => 'FASTA Header Builder';

  @override
  String get headerBuilderSubtitle =>
      'Choose, order and separate the fields of every sequence name.';

  @override
  String get headerFields => 'Header fields (in order)';

  @override
  String get headerAvailable => 'Available fields';

  @override
  String get headerSeparator => 'Separator';

  @override
  String get headerEmpty => 'Add at least one field.';

  @override
  String get headerConfirm => 'Confirm header';

  @override
  String get headerConfirmed => 'Header confirmed';

  @override
  String get headerUnique => 'Every header is unique.';

  @override
  String headerDuplicates(String n) {
    return '$n header(s) are shared by more than one record; numeric suffixes (_2, _3…) will be added. Add Accession to make them unique.';
  }

  @override
  String get moveLeft => 'Move earlier';

  @override
  String get moveRight => 'Move later';

  @override
  String get panelSeqType => 'Sequence Type';

  @override
  String get seqTypeBody =>
      'Choose how to align sequences. Codon-aware mode aligns protein translations via MAFFT then back-translates to DNA, preserving reading frame.';

  @override
  String get seqNtTitle => 'Nucleotide';

  @override
  String get seqNtBody => 'Align DNA/RNA sequences directly';

  @override
  String get seqAaTitle => 'Amino Acid';

  @override
  String get seqAaBody =>
      'Translate (codon_start, transl_table), then align proteins';

  @override
  String get seqCodonTitle => 'Codon-aware';

  @override
  String get seqCodonBody => 'Protein-guided codon alignment (CDS only)';

  @override
  String get recNtTitle => 'Nucleotide, direct DNA alignment';

  @override
  String get recNtBody =>
      'Best for rRNAs, tRNAs, non-coding regions and exploratory datasets. The FASTA goes directly to MAFFT.';

  @override
  String get recAaTitle => 'Amino Acid, protein alignment';

  @override
  String get recAaBody =>
      'Best when coding genes are divergent. DNA is translated first, MAFFT aligns proteins and the result is analysed as protein.';

  @override
  String get recCodonTitle => 'Codon-aware, CDS phylogenetics';

  @override
  String get recCodonBody =>
      'Recommended for CDS matrices with DNA output and codons preserved: validates CDS length, removes terminal stops, aligns proteins and back-translates.';

  @override
  String get geneticCodeInfo =>
      'Used for translation, frame detection and the saturation test. GenBank /transl_table qualifiers take precedence.';

  @override
  String get translationInfo =>
      'Sequences are translated before alignment. MAFFT runs in --amino mode.';

  @override
  String get codonWarning =>
      'Codon-aware mode only processes CDS features. Genes with no annotated CDS (rRNA, tRNA, etc.) will be skipped. The resulting alignment preserves reading frame and supports codon-position partitioning.';

  @override
  String get codonWorkflow => 'Codon-aware workflow';

  @override
  String get codonStep1 => 'Extract CDS';

  @override
  String get codonStep1Body =>
      'Only annotated CDS are used; rRNAs and tRNAs are skipped.';

  @override
  String get codonStep2 => 'Validate frame';

  @override
  String get codonStep2Body =>
      'Length divisible by 3 after codon_start and terminal stop handling.';

  @override
  String get codonStep3 => 'Translate';

  @override
  String get codonStep3Body => 'With the NCBI genetic code of each CDS.';

  @override
  String get codonStep4 => 'MAFFT protein';

  @override
  String get codonStep4Body => 'Proteins aligned with mafft --amino.';

  @override
  String get codonStep5 => 'Back-translate';

  @override
  String get codonStep5Body =>
      'The protein alignment guides a codon-preserving DNA alignment.';

  @override
  String get allowPseudo => 'Allow /pseudo and /pseudogene CDS';

  @override
  String get allowPseudoBody =>
      'Disabled by default. When off, pseudo CDS are excluded before translation.';

  @override
  String get allowStops => 'Allow internal stop codons';

  @override
  String get allowStopsBody =>
      'Disabled by default. When on, sequences with internal stops stay in the run and are reported with warnings.';

  @override
  String get allowBtWarnings =>
      'Continue when amino acid and codon checks report warnings';

  @override
  String get allowBtWarningsBody =>
      'Disabled by default. When off, any amino acid/codon mismatch after MAFFT fails the marker; when on, the warnings are logged.';

  @override
  String get pythonTitle => 'ClipKIT (codon-aware trimming)';

  @override
  String get pythonBody =>
      'ClipKIT runs with Python 3. SPLACE uses it for codon-aware alignments: it trims the protein alignment and SPLACE back-translates the result to codons.';

  @override
  String get pythonNotChecked => 'Python status was not checked yet.';

  @override
  String get pythonMissing => 'Python 3 was not found in PATH.';

  @override
  String pythonFound(String exe, String version) {
    return '$exe (Python $version)';
  }

  @override
  String clipkitFound(String version) {
    return 'ClipKIT $version available';
  }

  @override
  String get clipkitMissing => 'ClipKIT not installed';

  @override
  String get checkPython => 'Check Python';

  @override
  String get installClipkit => 'Install ClipKIT';

  @override
  String get mafftAuto => '--auto (automatic, recommended)';

  @override
  String get mafftLocalpair => '--localpair (L-INS-i, most accurate, slow)';

  @override
  String get mafftGlobalpair => '--globalpair (G-INS-i, global, slow)';

  @override
  String get mafftGenafpair => '--genafpair (E-INS-i, long gaps, slow)';

  @override
  String get mafftRetree2 => '--retree 2 (FFT-NS-2, fast)';

  @override
  String get mafftRetree1 => '--retree 1 (FFT-NS-1, fastest)';

  @override
  String get mafftMaxIterate => 'Max Iterations';

  @override
  String get mafftMaxIterateHelp => '0 = no refinement; 1000 = max';

  @override
  String get mafftThreads => 'Threads';

  @override
  String get mafftThreadsHelp => '(per job)';

  @override
  String mafftConcurrency(String cpus, String jobs, String threads) {
    return '$cpus logical CPUs: $jobs job(s) × $threads thread(s)';
  }

  @override
  String get mafftGapOpening => 'Gap Opening';

  @override
  String get mafftGapExtension => 'Gap Extension';

  @override
  String get mafftExtra => 'Additional Parameters';

  @override
  String get adjustDirectionHelp => 'Flips reverse-strand sequences';

  @override
  String get reorderHelp => 'Orders the output by similarity';

  @override
  String get optional => '(optional)';

  @override
  String get optionalLabel => 'Optional';

  @override
  String get mafftAminoNote =>
      'Amino-acid and codon-aware modes add --amino automatically (--adjustdirection does not apply to proteins).';

  @override
  String runMafftMode(String mode) {
    return 'Run MAFFT Alignment ($mode)';
  }

  @override
  String markersToAlign(String n) {
    return '$n marker(s) selected';
  }

  @override
  String get trimPromptTitle =>
      'Alignment complete - trim poorly aligned regions?';

  @override
  String get trimPromptBody =>
      'trimAl removes columns with excessive gaps or low conservation, improving phylogenetic inference quality.';

  @override
  String get trimPromptCodon =>
      'In Codon-aware mode, ClipKIT trims the protein alignment and SPLACE back-translates it to codon-preserving DNA.';

  @override
  String get trimPromptNo => 'No, skip trimming';

  @override
  String get trimPromptYes => 'Yes, set up trimming';

  @override
  String get trimAutomated1 => '-automated1 (heuristic, recommended)';

  @override
  String get trimStrict => '-strict (strict automated)';

  @override
  String get trimStrictplus => '-strictplus (strict+)';

  @override
  String get trimGappyout => '-gappyout (remove gappy columns)';

  @override
  String get trimNogaps => '-nogaps (remove all gap columns)';

  @override
  String get trimNoallgaps => '-noallgaps (only all-gap columns)';

  @override
  String get trimManual => 'Manual thresholds';

  @override
  String get trimGapThreshold => 'Gap Threshold';

  @override
  String get trimSimilarity => 'Similarity';

  @override
  String get trimConsistency => 'Consistency';

  @override
  String get trimWindow => 'Window Size';

  @override
  String get trimExtra => 'Additional Parameters';

  @override
  String get trimCodonNote =>
      'For codon-aware runs MAFFT aligned proteins: the trimmer trims the protein alignment, SPLACE back-translates it to codon-preserving DNA and validates the frame.';

  @override
  String get codonTrimBackend => 'Codon trimming backend';

  @override
  String get codonTrimAuto => 'Auto: trimAl, then ClipKIT if needed';

  @override
  String get codonTrimTrimal =>
      'trimAl, protein alignment with back-translation';

  @override
  String get codonTrimClipkit => 'ClipKIT (Python required)';

  @override
  String get clipkitMode => 'ClipKIT mode';

  @override
  String get codonTrimmingTitle => 'Codon trimming';

  @override
  String get codonTrimmingBody =>
      'Removes poorly aligned codon blocks before concatenation and before any codon-position filtering.';

  @override
  String get trimSkippedTag => 'Trimming skipped';

  @override
  String get trimming => 'Trimming';

  @override
  String get resultsTitle => 'Alignment Results';

  @override
  String get noTrimming => 'No trimming';

  @override
  String get viewAlignment => 'Show alignment';

  @override
  String get tabProtein => 'Protein';

  @override
  String get viewDna => 'DNA view';

  @override
  String get viewTranslate => 'Translate to amino acids';

  @override
  String possibleStops(String n) {
    return '$n possible internal stop(s)';
  }

  @override
  String get conservation => 'CONSERVATION';

  @override
  String get stopLegend => 'Stop';

  @override
  String get gapLegend => 'Gap';

  @override
  String get replaceGaps => 'Replace alignment gaps with ?';

  @override
  String get replaceGapsBody =>
      'Use when the downstream program should treat alignment gaps as missing data. Partition coordinates are preserved.';

  @override
  String get selectedLabel => 'Selected';

  @override
  String get saveNexus => 'Save NEXUS';

  @override
  String get savePartition => 'Save partition table';

  @override
  String get saveFasta => 'Save FASTA';

  @override
  String get savePhylip => 'Save PHYLIP';

  @override
  String get partitionName => 'Partition';

  @override
  String get partitionRange => 'Range';

  @override
  String get partitionSites => 'Sites';

  @override
  String get partitionType => 'Type';

  @override
  String get noneSelected => 'none selected';

  @override
  String get outgroupFromStep7 => 'chosen in step 7';

  @override
  String get substitutionModel => 'Substitution Model';

  @override
  String get modelTest => 'ModelFinder: auto best-fit (TEST)';

  @override
  String get modelMfp => 'MFP: ModelFinder Plus';

  @override
  String get modelMfpMerge => 'MFP+MERGE: partition model finder (recommended)';

  @override
  String get modelCustom => 'Custom…';

  @override
  String get bootstrapSupport => 'Bootstrap Support';

  @override
  String get bootUf1000 => 'UFBoot2: 1000 replicates (recommended)';

  @override
  String get bootUf5000 => 'UFBoot2: 5000 replicates';

  @override
  String get bootClassic => 'Classical bootstrap: 100 replicates (slow)';

  @override
  String get bootNone => 'No bootstrap';

  @override
  String get threadsLabel => 'Threads';

  @override
  String get extraFlags => 'Extra flags';

  @override
  String get spaceSeparated => '(space-separated)';

  @override
  String get perGeneTrees => 'Also infer per-gene trees';

  @override
  String get perGeneTreesBody =>
      'Runs a separate IQ-TREE analysis for each marker (saved in gene_trees/).';

  @override
  String get commandPreview => 'Command preview';

  @override
  String get saveZip => 'Save Results as ZIP';

  @override
  String get analysisComplete => 'Analysis complete';

  @override
  String get outputDirectory => 'Output directory';

  @override
  String get geneTreesLabel => 'Per-gene trees';

  @override
  String get howToCite => 'How to Cite';

  @override
  String get howToCiteHint =>
      'If you use these results, cite SPLACE and the tools it ran';

  @override
  String get copy => 'Copy';

  @override
  String metaGenomes(String n) {
    return '$n genomes';
  }

  @override
  String metaFastaFiles(String n) {
    return '$n FASTA files';
  }

  @override
  String metaSamples(String n) {
    return '$n samples';
  }

  @override
  String metaSelected(String selected, String total) {
    return '$selected of $total selected';
  }

  @override
  String metaAligned(String done, String total) {
    return '$done/$total aligned';
  }

  @override
  String metaTrimmed(String n) {
    return '$n trimmed';
  }

  @override
  String metaMatrix(String taxa, String sites) {
    return '$taxa taxa × $sites sites';
  }

  @override
  String reviewSelectedCount(String selected, String total) {
    return '$selected of $total genomes enabled';
  }

  @override
  String get progressSearchTitle => 'Searching NCBI by taxonomic name';

  @override
  String progressSearchingFor(String name) {
    return 'Searching for $name';
  }

  @override
  String get progressFetchTitle => 'Importing GenBank records';

  @override
  String progressCounter(String current, String total) {
    return '$current of $total';
  }

  @override
  String progressSpecies(String current, String total) {
    return '$current of $total species';
  }

  @override
  String get progressRateLabel => 'Request rate';

  @override
  String get progressRateWithKey => 'Up to 10 requests/s';

  @override
  String get progressRateDefault => 'Up to 3 requests/s';

  @override
  String get progressApiEnabled => 'Saved and active';

  @override
  String get progressApiDisabled => 'Not provided';

  @override
  String get progressCurrentLabel => 'Current item';

  @override
  String get progressCurrentSpecies => 'Current species';

  @override
  String get progressStatusIdle => 'Waiting to start';

  @override
  String get progressFooter =>
      'SPLACE is contacting NCBI and updating the progress continuously.';

  @override
  String progressTaxonomyTitle(String source) {
    return 'Fetching taxonomy from $source';
  }

  @override
  String progressTaxonomyFooter(String source) {
    return 'SPLACE is querying $source and updating the records as matches are found.';
  }

  @override
  String get progressWorking => 'Working…';

  @override
  String get colStatus => 'Status (SALVE)';

  @override
  String fetchTaxonomyFrom(String source) {
    return 'Fetch taxonomy ($source)';
  }

  @override
  String taxonomyMatched(String source) {
    return '$source match found.';
  }

  @override
  String taxonomyNoMatch(String source) {
    return 'No $source match for this name.';
  }

  @override
  String get taxonomyCardTitle => 'Taxonomy of the records';

  @override
  String get taxonomySourceLabel => 'Taxonomy source';

  @override
  String get sourceGbifBody =>
      'Accepted names and authorship from the Global Biodiversity Information Facility backbone.';

  @override
  String sourceNcbiBody(String withId, String total) {
    return 'Lineage of the taxon id in each GenBank record ($withId of $total records have one; the others are searched by name).';
  }

  @override
  String get sourceSalveBody =>
      'Brazilian fauna only. ICMBio\'s national assessment: taxonomy, extinction-risk category and the species sheet, from the SALVE public API.';

  @override
  String get badgeOnline => 'online';

  @override
  String get salveHint =>
      '**Only species of the Brazilian fauna.** SALVE holds ICMBio\'s national extinction-risk assessments; species from elsewhere have no sheet and are listed as not found (use GBIF or NCBI for them). Every section of each sheet is read online: ranks, category in the IUCN colours and criteria, common names, distribution (endemism, biomes, states, watersheds), natural history (habitats, diet, reproduction, migration), population (trend, generation time), threats, uses, official lists, assessment history, conservation actions, protected areas, research, assessors, citation and DOI.';

  @override
  String get poweredByDataFishing => 'Method of';

  @override
  String get metaSatKeep => 'keep all positions';

  @override
  String metaSatExclude(String positions) {
    return 'remove position $positions';
  }

  @override
  String get metaSatAmino => 'use amino acids';

  @override
  String get apiKeyActive => 'Active';

  @override
  String get apiKeyInactive => 'Not set';

  @override
  String get apiKeyWithout => 'Without a key';

  @override
  String get apiKeyWith => 'With a key';

  @override
  String get apiKeyUsedFor =>
      'Used by the NCBI search, accession downloads and NCBI Taxonomy.';

  @override
  String get apiKeyHowTitle => 'How to get a key';

  @override
  String get apiKeyHow1 => 'Sign in to your NCBI account (free).';

  @override
  String get apiKeyHow2 =>
      'Open Account settings → API Key Management and click “Create an API Key”.';

  @override
  String get apiKeyHow3 => 'Copy the key and paste it in the field below.';

  @override
  String get apiKeyOpenNcbi => 'Open NCBI account';

  @override
  String get apiKeyPaste => 'Paste';

  @override
  String get apiKeyFieldLabel => 'Your API key';

  @override
  String get apiKeyLooksValid => 'Format looks valid';

  @override
  String get apiKeyStored =>
      'Stored only on this computer (Documents/SPLACE/settings.json).';

  @override
  String get apiKeyEstimateTitle => 'Estimated download time';

  @override
  String apiKeyEstimateGenomes(String n) {
    return '$n genomes';
  }

  @override
  String apiKeyEstimateDetail(String batch, String rate) {
    return '$batch genomes per request · up to $rate requests/s';
  }

  @override
  String apiKeyEstimateGain(String factor, String saved) {
    return '$factor× faster with a key · $saved less';
  }

  @override
  String get apiKeyEstimateNote =>
      'Download of mitochondrial genomes not yet in the local cache, parsing included, from times measured with this app on Nov 8 2026 (Testudines, about 75 kB per record): 50, 200 and 700 genomes took 1.6, 4.4 and 10.1 s with a key and 5.5, 16.9 and 64.6 s without one. The search and the list of results add a few seconds; records already in the cache load in a moment. The actual time depends on your connection and on NCBI\'s load.';

  @override
  String get apiKeyCurrent => 'Current';

  @override
  String reviewFound(String n, String query) {
    return '$n genomes found · $query';
  }

  @override
  String get reviewColOrganism => 'Organism';

  @override
  String get reviewColStatus => 'Status';

  @override
  String reviewUnverifiedNote(String n) {
    return '$n genomes are flagged UNVERIFIED by NCBI. We recommend leaving them out.';
  }

  @override
  String get reviewDeselectUnverified => 'Deselect UNVERIFIED';

  @override
  String get reviewSelectMenu => 'Select';

  @override
  String reviewImportN(String n) {
    return 'Import $n';
  }

  @override
  String get reviewNoMatch => 'No genome matches the filter.';

  @override
  String reviewSelectedOf(String selected, String total) {
    return '$selected of $total selected';
  }

  @override
  String metaTrees(String n) {
    return '$n tree(s)';
  }

  @override
  String get itemQueued => 'Waiting';

  @override
  String itemCached(String note) {
    return 'From the cache · $note';
  }

  @override
  String get itemSkipped => 'Skipped';

  @override
  String itemsRunning(String n) {
    return '$n running';
  }

  @override
  String itemsQueued(String n) {
    return '$n waiting';
  }

  @override
  String itemsDone(String n) {
    return '$n done';
  }

  @override
  String itemsCached(String n) {
    return '$n from the cache (unchanged)';
  }

  @override
  String itemsSkipped(String n) {
    return '$n skipped';
  }

  @override
  String itemsFailed(String n) {
    return '$n failed';
  }

  @override
  String get trimToolNt =>
      'Choose one tool: trimAl or ClipKIT keeping whole codons in coding genes (recommended), or trimAl column by column (splits codons).';

  @override
  String get trimToolAa => 'trimAl trims the protein alignment.';

  @override
  String get trimToolCodon =>
      'The codon alignment is trimmed keeping whole codons: by trimAl (-backtrans) or ClipKIT (--codon).';

  @override
  String runTrimWith(String tool) {
    return 'Trim with $tool';
  }

  @override
  String get trimDescAutomated1 =>
      'trimAl looks at the alignment and picks the best of -gappyout and -strict by itself. Recommended for maximum-likelihood trees.';

  @override
  String get trimDescStrict =>
      'Removes gap-rich columns and also columns with low similarity between sequences.';

  @override
  String get trimDescStrictplus =>
      'Like -strict, but tuned for neighbour-joining trees; removes more.';

  @override
  String get trimDescGappyout =>
      'Removes only the gap-rich columns, with a threshold computed from the gap distribution of the alignment.';

  @override
  String get trimDescNogaps =>
      'Removes every column that has at least one gap. Can remove a lot.';

  @override
  String get trimDescNoallgaps =>
      'Removes only columns made entirely of gaps. The lightest option.';

  @override
  String get trimDescManual =>
      'You set the gap, similarity and consistency thresholds yourself.';

  @override
  String get clipkitDescGappy =>
      'Removes columns with more gaps than the threshold (90 %). Simple and conservative.';

  @override
  String get clipkitDescKpic =>
      'Keeps only parsimony-informative and constant columns; removes the rest.';

  @override
  String get clipkitDescKpicGappy => 'kpic, and also removes gap-rich columns.';

  @override
  String get clipkitDescSmartGap =>
      'Chooses the gap threshold from the alignment itself (ClipKIT default).';

  @override
  String get clipkitRequired =>
      'Codon-aware trimming uses ClipKIT, which is not installed yet. Install it below (needs Python 3).';

  @override
  String get frameInfo =>
      'Where the first complete codon starts in this alignment. Base 1 (+1): the alignment starts at the 1st position of a codon, as in most CDS. +2 or +3: it starts 1 or 2 bases into a codon. SPLACE picks the frame with the fewest stop codons; it is used to translate, to trim whole codons and to split the 1st, 2nd and 3rd codon positions (saturation test and partitions).';

  @override
  String frameOption(String n) {
    return 'Base $n (+$n)';
  }

  @override
  String frameTag(String n) {
    return 'Codons start at base $n';
  }

  @override
  String get viewerHint =>
      'Hover a base to see its position · drag with the right button to scroll · drag the edge of the names to resize, double-click to fit.';

  @override
  String viewerColumn(String n) {
    return 'Column $n';
  }

  @override
  String viewerSeqPos(String n) {
    return 'position $n in the sequence (gaps not counted)';
  }

  @override
  String get viewerGap => 'gap';

  @override
  String viewerIdentical(String pct) {
    return '$pct% of sequences share this base';
  }

  @override
  String viewerCodonPos(String n) {
    return 'codon position $n';
  }

  @override
  String get viewerNames => 'Sequences';

  @override
  String siCacheInfo(String cached, String total, String todo) {
    return '$cached of $total searches were already done before and will be reused; only $todo will be sent.';
  }

  @override
  String siCacheAll(String total) {
    return 'All $total searches were done before: the results appear right away, nothing is sent.';
  }

  @override
  String get siCacheClear => 'Clear saved results';

  @override
  String get trimFrameWarning =>
      '**trimAl removes single columns, not codons.** In a protein-coding gene aligned as DNA, removing one or two bases of a codon shifts the reading frame of everything after it: 1st, 2nd and 3rd codon positions get mixed. Without \"Keep the reading frame\" these genes become one partition without codon positions and are left out of the saturation test.';

  @override
  String get frameLostTag => 'Reading frame lost in trimming';

  @override
  String get frameLostInfo =>
      'trimAl removed single columns of this coding gene, so its length may not be a multiple of 3 and the codon positions are unknown. It is concatenated as one partition and left out of the saturation test. Turn on \"Keep the reading frame\" in 5e and trim again to keep codon partitions.';

  @override
  String frameShiftTag(String n) {
    return '$n sequence(s) with frame shifts';
  }

  @override
  String get frameShiftInfo =>
      'MAFFT placed gaps that are not multiples of 3 inside these sequences: from that gap on, their bases sit in the wrong codon position of the column. This is common when coding genes are aligned as nucleotides. Codon partitions and the saturation test are approximate for them; the Codon-aware mode (5b) avoids it.';

  @override
  String get saturationRunning => 'Saturation test';

  @override
  String saturationProgress(String done, String total) {
    return '$done of $total alignments';
  }

  @override
  String get matrixRunning => 'Building the supermatrix';

  @override
  String get matrixRunningBody => 'Concatenating the alignments';

  @override
  String get iqtreeRunning => 'IQ-TREE: maximum-likelihood tree';

  @override
  String get raxmlRunning => 'RAxML-NG: maximum-likelihood tree';

  @override
  String treeRunSubtitle(String model, String elapsed) {
    return 'Model $model · elapsed $elapsed';
  }

  @override
  String get geneTreesProgress => 'Gene trees';

  @override
  String get localFooter => 'Runs on this computer.';

  @override
  String get treeFooter =>
      'Runs on this computer. Long analyses can take minutes or hours; you can follow each stage here.';

  @override
  String get treeDetailPrepare => 'Reading the alignment and partitions';

  @override
  String get treeDetailModelStart => 'Testing substitution models';

  @override
  String treeDetailModels(String done, String total) {
    return 'Model $done of up to $total tested';
  }

  @override
  String treeDetailStarts(String done, String total) {
    return 'ML search from starting tree $done of $total';
  }

  @override
  String get treeDetailSearchStart => 'Building starting trees';

  @override
  String treeDetailIteration(String n) {
    return 'Tree search: iteration $n';
  }

  @override
  String treeDetailLeft(String time) {
    return '~$time left';
  }

  @override
  String treeDetailBoot(String done, String total) {
    return 'Bootstrap replicate $done of $total';
  }

  @override
  String get treeDetailSupport => 'Computing branch support';

  @override
  String get treeDetailFinish => 'Writing the results';

  @override
  String treeDetailBestModel(String model) {
    return 'Best model: $model';
  }

  @override
  String treeDetailUfboot(String n) {
    return '$n UFBoot replicates';
  }

  @override
  String treeDetailConverged(String n) {
    return 'Bootstrap converged after $n';
  }

  @override
  String get treeStagePrepare => 'Alignment';

  @override
  String get treeStageModel => 'Model selection';

  @override
  String get treeStageSearch => 'ML tree search';

  @override
  String get treeStageSupport => 'Branch support';

  @override
  String get treeStageBootstrap => 'Bootstrap';

  @override
  String get treeStageFinish => 'Results';

  @override
  String get step9 => 'Tree Viewer';

  @override
  String metaViewable(String n) {
    return '$n tree(s)';
  }

  @override
  String get viewTree => 'View tree';

  @override
  String get treeViewIntro =>
      'Look at the inferred trees. Root them on the outgroup, at the midpoint or on any branch, check the support of each clade, find species and export the figure.';

  @override
  String get treeChoose => 'Tree';

  @override
  String get treeMissing => 'The tree file was not found.';

  @override
  String get treeKindMl => 'ML tree with support';

  @override
  String get treeKindConsensus => 'bootstrap consensus';

  @override
  String get treeKindFbp => 'ML tree, bootstrap support (FBP)';

  @override
  String get treeKindTbe => 'ML tree, transfer bootstrap (TBE)';

  @override
  String get treeKindBest => 'best ML tree (no support)';

  @override
  String treeKindGene(String gene) {
    return 'gene tree $gene';
  }

  @override
  String treeParseError(String error) {
    return 'This tree could not be read: $error';
  }

  @override
  String get treeRootAsIs => 'As inferred';

  @override
  String get treeRootOutgroup => 'Outgroup';

  @override
  String get treeRootMidpoint => 'Midpoint';

  @override
  String get treeRootManual => 'Rooted on the chosen branch';

  @override
  String get treePhylogram => 'Phylogram';

  @override
  String get treeCladogram => 'Cladogram';

  @override
  String get treeShowSupport => 'Support';

  @override
  String get treeItalic => 'Italic';

  @override
  String get treeSearch => 'Find a species';

  @override
  String get treeFont => 'Font size';

  @override
  String treeStats(String tips, String length) {
    return '$tips tips · tree length $length';
  }

  @override
  String get treeHint =>
      'Click a branch or a name to select it: root the tree there, rotate, collapse or shorten it. Ctrl+Z undoes.';

  @override
  String treeClade(String n) {
    return 'Clade with $n tips';
  }

  @override
  String treeSupportValue(String value) {
    return 'support $value';
  }

  @override
  String get treeRootHere => 'Root here';

  @override
  String get treeCollapse => 'Collapse';

  @override
  String get treeExpand => 'Expand';

  @override
  String get treeClearSelection => 'Clear';

  @override
  String get treeSupportLegend => 'Branch support:';

  @override
  String get treeSupportInfoUfboot =>
      'IQ-TREE labels are SH-aLRT/UFBoot (or UFBoot alone); colours use the last value. UFBoot ≥ 95 is strong support (Hoang et al. 2018); SH-aLRT ≥ 80 together with it makes the clade well supported.';

  @override
  String get treeSupportInfoBootstrap =>
      'Standard (non-parametric) bootstrap: ≥ 70 is usually taken as good support (Hillis & Bull 1993).';

  @override
  String get siSelfTag => 'this record';

  @override
  String siSelfInfo(String accessions) {
    return 'BLAST found this genome\'s own sequence in GenBank ($accessions). The data is there, but a record cannot confirm its own name, so it does not count as evidence. Other records of the species would.';
  }

  @override
  String siVerdictSelfOnly(String accessions, String current) {
    return 'GenBank has this genome\'s sequence ($accessions) but no other record of $current, so GenBank cannot confirm the name.';
  }

  @override
  String siVerdictNoReference(String current, String best, String identity) {
    return 'No other sequence of $current was found. The closest is $best ($identity%), below 98% — the usual distance between species — so the current name is not contested.';
  }

  @override
  String get siNoReference => 'No reference';

  @override
  String get apiKeyEstimateColSize => 'Download';

  @override
  String get apiKeyEstimateColGain => 'Gain';

  @override
  String apiKeyEstimateSaved(String n, String saved) {
    return 'For $n genomes the key saves about $saved.';
  }

  @override
  String missTypeNotSelected(String type) {
    return 'annotated as $type';
  }

  @override
  String missTypeNotSelectedInfo(String gene, String type) {
    return '$gene is annotated in this record as $type, a feature type not selected in step 3. Select $type there to include it.';
  }

  @override
  String get missGeneOnly => 'only a \'gene\' feature';

  @override
  String missGeneOnlyInfo(String gene, String location) {
    return 'The record has a \'gene\' feature for $gene ($location) but no CDS/rRNA feature, which is what SPLACE reads. The gene is in the sequence; the annotation is incomplete. Load a better-annotated record, or annotate the genome (e.g. with MitoS or MitoFinder).';
  }

  @override
  String missNoAnnotation(String types) {
    return 'no $types annotated';
  }

  @override
  String missNoAnnotationInfo(String types) {
    return 'This record has no $types annotated at all, so no gene can be taken from it. It is probably a raw submission; annotate the genome (e.g. with MitoS or MitoFinder) or remove it.';
  }

  @override
  String get missUnknownName => 'unrecognized name';

  @override
  String missUnknownNameInfo(String names, String gene) {
    return 'The record has genes with names SPLACE did not recognize ($names). One of them may be $gene: fix the name in step 3 (unknown names).';
  }

  @override
  String get missPartial => 'partial genome';

  @override
  String missPartialInfo(String length, String typical, String gene) {
    return 'The record has $length bp, well below the typical genome of this dataset ($typical bp): it is a partial genome and does not cover $gene.';
  }

  @override
  String get missNotAnnotated => 'not annotated';

  @override
  String missNotAnnotatedInfo(String gene) {
    return '$gene is not among the genes annotated in this record: it is absent from the sequence or the authors did not annotate it.';
  }

  @override
  String get missNotInFiles => 'not in the files';

  @override
  String missNotInFilesInfo(String gene) {
    return 'No FASTA file has a $gene sequence for this sample.';
  }

  @override
  String missPresentIn(String present, String total, String missing) {
    return 'in $present of $total samples · missing in $missing';
  }

  @override
  String missSampleCount(String missing, String total) {
    return 'missing $missing of $total markers';
  }

  @override
  String missSummary(
    String markers,
    String selected,
    String samples,
    String total,
  ) {
    return '$markers of the $selected selected markers are missing in at least one sample; $samples of $total samples lack at least one marker.';
  }

  @override
  String get missMatrixNote =>
      'In the supermatrix (step 7) a missing marker becomes \'?\' for that sample when “Allow missing sequences” is on; otherwise markers absent from any taxon are left out. A few gaps are normal; a sample missing most markers adds little and may be worth removing.';

  @override
  String missFixableNote(String n) {
    return '$n gap(s) can be fixed in SPLACE (orange tags). Hover a tag to see how.';
  }

  @override
  String get missByMarker => 'By marker';

  @override
  String get missBySample => 'By sample';

  @override
  String get missHoverHint =>
      'Each tag says why the marker is missing; hover it for details.';

  @override
  String dupSummary(String genes, String genomes) {
    return '$genes gene(s) have more than one copy annotated in $genomes genome(s).';
  }

  @override
  String get dupFastaHelp =>
      'These samples have more than one sequence of the same locus in the FASTA files. The “Used” column shows what goes into the alignment, following the “Duplicates and names” setting.';

  @override
  String get dupIdentical => 'identical copies';

  @override
  String get dupIdenticalInfo =>
      'Same sequence (e.g. inverted-repeat copies in plastomes): any copy gives the same result.';

  @override
  String get dupContained => 'fragment of another copy';

  @override
  String get dupContainedInfo =>
      'The shorter copy is part of the longer one: usually a duplicated or truncated annotation. Use the complete copy.';

  @override
  String get dupDifferent => 'different copies';

  @override
  String get dupDifferentInfo =>
      'The copies differ in sequence (paralogues, NUMT-like copies or annotation errors). The choice changes the alignment: prefer the copy closest to the other genomes.';

  @override
  String get siCardTitle => 'Species identification';

  @override
  String get siCardBody =>
      'Compares a barcode gene (COI, 16S…) of each genome with GenBank (BLAST) and BOLD to catch misidentified records before fetching the taxonomy, which depends on the names.';

  @override
  String get siCardOptional => 'Optional';

  @override
  String get siCardRerun => 'Identify again';

  @override
  String siSumChecked(String n) {
    return '$n genomes checked';
  }

  @override
  String siSumConfirmed(String n) {
    return '$n names confirmed';
  }

  @override
  String siSumFlagged(String n) {
    return '$n possibly misidentified';
  }

  @override
  String siSumNoReference(String n) {
    return '$n without a reference';
  }

  @override
  String get headerDefault => 'Default header · used when aligning';

  @override
  String viewerMafftColumn(String n) {
    return 'column $n of the MAFFT alignment';
  }

  @override
  String get satSourceMafft => 'Tested on the MAFFT alignment (not trimmed).';

  @override
  String get satSourceClipkit =>
      'Tested on the trimmed alignment (ClipKIT on the proteins, back-translated): only whole codons were removed, so every base keeps its codon position.';

  @override
  String satSourceTrimmed(String kept, String total) {
    return 'Tested on the trimmed alignment: trimAl kept $kept of the $total MAFFT codons, and SPLACE removed only whole codons, so every column keeps the codon position it had in the MAFFT alignment (hover a trimmed column in step 5 to see which MAFFT column it is).';
  }

  @override
  String satSourceFrame(String frame) {
    return 'Codons now start at base 1 because the incomplete codon before base $frame of the MAFFT alignment was dropped; the positions are the same.';
  }

  @override
  String get preserveCaseHelp => 'Keeps the input letter case';

  @override
  String get configureTrimming => 'Set up trimming';

  @override
  String panelTrimParams(String tool) {
    return '$tool parameters';
  }

  @override
  String trimmedCount(String n, String total) {
    return '$n of $total markers trimmed';
  }

  @override
  String get showRemovedColumns => 'Show removed regions';

  @override
  String get nextStop => 'Next stop';

  @override
  String get translationTable => 'Translation table';

  @override
  String translationTableTitle(String code, String name) {
    return 'Genetic code $code · $name';
  }

  @override
  String get translationTableBody =>
      'Codons of the genetic code chosen above for the amino-acid view (it is never guessed from the data). Stop codons are red; the badge counts the internal stops of this alignment that use each one.';

  @override
  String codonSecondBase(String base) {
    return '2nd $base';
  }

  @override
  String codonDiffers(String codon, String aa, String standard) {
    return '$codon → $aa in this code ($standard in the standard code)';
  }

  @override
  String get codonDiffersLegend => 'Differs from the standard code';

  @override
  String stopCodonsFound(String list) {
    return 'Internal stops in this alignment: $list';
  }

  @override
  String removedColumns(String n, String total) {
    return '$n of $total columns removed';
  }

  @override
  String get removedColumnsInfo =>
      'The whole MAFFT alignment is shown; grey columns were removed by trimming. Turn off \"Show removed regions\" to see only the kept columns.';

  @override
  String get removedLegend => 'Removed';

  @override
  String get viewerRemoved => 'removed by trimming';

  @override
  String viewerCodon(String codon) {
    return 'codon $codon';
  }

  @override
  String get plotObservedVsCorrectedBody =>
      'Each point is a pair of sequences. Points on the dashed line: every substitution is visible; points levelling off below it: multiple hits hide substitutions (saturation).';

  @override
  String get plotTsTvBody =>
      'Transitions (red) usually saturate first: when they stop growing while transversions (blue) keep rising, the subset is saturated.';

  @override
  String plotPositions(String model) {
    return 'p of each codon position × $model distance';
  }

  @override
  String get plotPositionsBody =>
      'Same pairs, one colour per codon position, against the distance of the whole gene. A position whose points flatten while the others keep rising is saturated (usually the 3rd).';

  @override
  String get allPositions => 'all positions';

  @override
  String get compositionTitle => 'Composition and variability';

  @override
  String get compositionInfo =>
      'Chi-square test of base-frequency homogeneity across taxa (contingency test, as in PAUP*), and the number of taxa whose own composition differs from the subset (df = 3, p < 0.05, as IQ-TREE reports). Saturated positions often drift in composition (e.g. GC3), which can also bias the tree.';

  @override
  String get metricVariable => 'Variable';

  @override
  String get metricInformative => 'Parsimony-informative';

  @override
  String get metricTaxaFailing => 'Taxa failing';

  @override
  String get uncorrectablePairs => 'pairs beyond correction';

  @override
  String get heatmapIntensity =>
      'The shade inside each colour follows the slope and the uncorrectable pairs.';

  @override
  String get treeRootLabel => 'Root';

  @override
  String get treeShapeLabel => 'Shape';

  @override
  String get treeSupportValues => 'Values';

  @override
  String get treeSupportDots => 'Dots';

  @override
  String get treeSupportHidden => 'Hide';

  @override
  String get treeDisplay => 'Display';

  @override
  String get treeAlignLabels => 'Aligned on the right';

  @override
  String get treeExport => 'Export';

  @override
  String treeMatches(String n) {
    return '$n found';
  }

  @override
  String treeHoverTip(String name, String length) {
    return '$name · branch length $length';
  }

  @override
  String treeHoverClade(String n, String support, String length) {
    return 'Clade with $n tips · support $support · branch length $length — click to select it';
  }

  @override
  String get treeOutgroupLegend => 'Outgroup';

  @override
  String get ord1 => '1st';

  @override
  String get ord2 => '2nd';

  @override
  String get ord3 => '3rd';

  @override
  String get listAnd => 'and';

  @override
  String codonPositionN(String n) {
    return '$n position';
  }

  @override
  String reasonInsufficient(String sites, String pairs) {
    return 'Too little data to judge ($sites sites, $pairs sequence pairs).';
  }

  @override
  String reasonUndefined(String pct) {
    return '$pct% of the sequence pairs differ so much that the distance cannot be corrected.';
  }

  @override
  String reasonLowDivergence(String p) {
    return 'The sequences differ very little (mean p = $p), so saturation is not a concern.';
  }

  @override
  String reasonSlopeOk(String slope) {
    return 'PhyKIT saturation $slope: the observed differences still grow with the tree distance.';
  }

  @override
  String reasonSlopeModerate(String slope, String limit) {
    return 'PhyKIT saturation $slope (below $limit): the observed differences grow more slowly than the tree distance.';
  }

  @override
  String reasonSlopeSaturated(String slope, String limit) {
    return 'PhyKIT saturation $slope (below $limit): the observed differences have almost stopped growing with the tree distance; many changes were overwritten by later ones.';
  }

  @override
  String reasonIss(String iss) {
    return 'Iss $iss (0.8 or more): the sites look close to random.';
  }

  @override
  String get recTitle => 'What to do with the codon positions';

  @override
  String recBasedOn(String taxa, String length) {
    return 'Based on the concatenated matrix ($taxa taxa, $length bp). Each position is tested on its own.';
  }

  @override
  String get recWhatToDo => 'Suggestion';

  @override
  String recExcludeOne(String positions) {
    return 'Remove the $positions codon position from the matrix (step 7b) and keep the other two, partitioned by codon position.';
  }

  @override
  String recExcludeMany(String positions) {
    return 'Remove the $positions codon positions from the matrix (step 7b) and keep the remaining one.';
  }

  @override
  String get recWhy =>
      'Why: a saturated position has changed so many times that its differences no longer show how long ago the lineages split. It mostly adds noise (homoplasy) and can pull long branches together.';

  @override
  String recCautionOne(String positions) {
    return 'The $positions position shows moderate saturation: keep it in its own partition and check whether the tree changes without it.';
  }

  @override
  String recCautionMany(String positions) {
    return 'The $positions positions show moderate saturation: keep them in their own partitions and check whether the tree changes without them.';
  }

  @override
  String applyExclude(String positions) {
    return 'Mark $positions in step 7b';
  }

  @override
  String get applyKeepAll =>
      'Leave every \"Exclude position\" unmarked in step 7b';

  @override
  String applyPerLocusInfo(String loci) {
    return 'Each gene excludes only the positions saturated in it. Genes whose result differs from the overall suggestion: $loci.';
  }

  @override
  String recLociDiffer(String n) {
    return '$n gene(s) have a result different from the overall suggestion; see their rows in the heatmap.';
  }

  @override
  String get recStatusKeep => 'Step 7b now keeps all positions.';

  @override
  String recStatusExclude(String positions) {
    return 'Step 7b now excludes the $positions position(s).';
  }

  @override
  String recStatusPerLocus(String n) {
    return 'Step 7b now excludes positions per gene ($n genes).';
  }

  @override
  String get saturationHeatmapHelp =>
      'Each cell shows PhyKIT\'s saturation value (slope of p on the tree distance) — 1 means no saturation, near 0 saturated. Blue is OK, red is saturated.';

  @override
  String get rankSpecies => 'Species';

  @override
  String get rankGenus => 'Genus';

  @override
  String get rankFamily => 'Family';

  @override
  String get rankOrder => 'Order';

  @override
  String get rankClass => 'Class';

  @override
  String get rankPhylum => 'Phylum';

  @override
  String get rankKingdom => 'Kingdom';

  @override
  String get tipFieldOriginal => 'Name as in the tree';

  @override
  String get tipFieldAccession => 'Accession';

  @override
  String get tipFieldVoucher => 'Voucher';

  @override
  String get treeBold => 'Bold';

  @override
  String get treeParens => 'In parentheses';

  @override
  String get treeRemovePart => 'Remove this part';

  @override
  String get treeAddPart => 'Add field';

  @override
  String get treeMovePart => 'Move left';

  @override
  String get treeLabelFormat => 'Label format';

  @override
  String get treeLabelParts => 'Parts of the label';

  @override
  String get treeLabelPartsInfo =>
      'Each label is built from these fields, in this order. I = italic, B = bold, ( ) = in parentheses. The values come from the GenBank records and the taxonomy fetched in step 2; a field a tip does not have is left out.';

  @override
  String get treeTabTree => 'Tree';

  @override
  String get treeTabBranches => 'Branches and support';

  @override
  String get treeTabLabels => 'Labels';

  @override
  String get treeTabColours => 'Colours';

  @override
  String get treeOrderLabel => 'Order and rotation';

  @override
  String get treeLadderUp => 'Ladderize ↑';

  @override
  String get treeLadderDown => 'Ladderize ↓';

  @override
  String get treeFlip => 'Flip tree';

  @override
  String get treeRotateClade => 'Rotate clade';

  @override
  String get treeRotateInfo =>
      'Rotating a clade swaps the order of its branches; the relationships do not change. Select a clade by clicking its branch, then rotate it. Flip tree mirrors the whole tree.';

  @override
  String get treeLineWidth => 'Line width';

  @override
  String get treePointSize => 'Point size';

  @override
  String get treePoints => 'Points';

  @override
  String get treeTipDots => 'Points on the tips';

  @override
  String get treeSupportMin => 'Show support from';

  @override
  String get treeSupportMinInfo =>
      'Only nodes with support at or above this value get a number or a point (for example, points only on nodes with UFBoot ≥ 95).';

  @override
  String get treeSupportAll => 'All';

  @override
  String get treeLongBranches => 'Long branches';

  @override
  String get treeLongBranchesInfo =>
      'A shortened branch is drawn shorter than it is and marked with //, so the rest of the tree has room. Branches over 5× the median length and over 30% of the tree height count as long. You can also select any branch and shorten it. The real length stays in the info bar and in the exported Newick.';

  @override
  String treeShortenLong(String n) {
    return 'Shorten long branches ($n)';
  }

  @override
  String get treeRestoreLengths => 'Restore all lengths';

  @override
  String get treeShortenBranch => 'Shorten branch //';

  @override
  String get treeRestoreBranch => 'Restore length';

  @override
  String get treeShortenedLegend => 'shortened branch (not to scale)';

  @override
  String get treeColorRows => 'Colour rows by';

  @override
  String get treeColorRowsInfo =>
      'Each taxon of the chosen rank gets a background colour across the whole row, with its name on the right. The ranks come from the taxonomy fetched in step 2 (or the species and genus in the tip names).';

  @override
  String get treeColorNone => 'None';

  @override
  String get treeNoTaxonomy =>
      'No taxonomy for these tips. Fetch it in step 2.';

  @override
  String get treeExportInfo =>
      'SVG and PNG are the figure as you see it (labels, colours, legend, shortened branches), without the selection or search marks. Newick has the topology and the real branch lengths.';

  @override
  String get treeUndo => 'Undo';

  @override
  String get treeRedo => 'Redo';

  @override
  String get treeReset => 'Start over';

  @override
  String get zipTitle => 'Save everything';

  @override
  String get zipBody =>
      'One ZIP with the sequences, alignments, saturation tables, matrices, partitions, trees and logs of this project.';

  @override
  String get concatHowTitle => 'How the supermatrix is built';

  @override
  String get concatHow1 =>
      'The alignment of each locus is placed side by side, one row per taxon (matched by the FASTA header).';

  @override
  String get concatHow2 =>
      'A taxon that lacks a locus gets ? in that stretch (missing data), so every row has the same length.';

  @override
  String get concatHow3 =>
      'The partition file tells IQ-TREE / RAxML-NG where each locus (and each codon position) starts and ends, so each part can get its own substitution model.';

  @override
  String get concatPictureCaption => '? = taxon without that locus';

  @override
  String get concatBlockSource => 'Which alignments to join';

  @override
  String get concatBlockSourceBody =>
      'The trimmed alignments of step 5 (without the badly aligned columns) or the MAFFT output as it is.';

  @override
  String get concatBlockPartitions => 'How to split the matrix into partitions';

  @override
  String get concatBlockPartitionsBody =>
      'Each partition gets its own substitution model in the tree search. Codon positions can also be removed here.';

  @override
  String get schemeByLocusCodonBody =>
      'Each locus is a partition and, in protein-coding loci, each codon position (1st, 2nd, 3rd) gets its own. Codon positions evolve at very different rates; with MFP+MERGE, IQ-TREE joins the partitions that behave alike.';

  @override
  String get schemeByLocusBody =>
      'One partition per locus: each gene gets its own model, and its three codon positions share it. Good for rRNAs, non-coding loci or loci with few sites.';

  @override
  String get schemeByCodonBody =>
      'Three partitions shared by all coding loci (all 1st positions together, all 2nd, all 3rd); each non-coding locus stays a partition of its own. Fewer models than by locus and position.';

  @override
  String get schemeNoneBody =>
      'One model for the whole matrix (no partition file). Simplest, but usually fits worse when genes or positions evolve differently.';

  @override
  String get codonExclusionHelp =>
      'Removed positions leave the coding loci only. Step 6 (saturation) can fill this in.';

  @override
  String get concatSummaryTitle => 'What will be built';

  @override
  String get concatSummaryLoci => 'Loci';

  @override
  String get concatSummarySource => 'Alignments';

  @override
  String get concatNoOutgroup => 'none (midpoint)';

  @override
  String get concatNothingExcluded => 'nothing';

  @override
  String get treeXScale => 'Horizontal scale';

  @override
  String get treeXScaleInfo =>
      'Stretches (above 1) or shrinks (below 1) the branches on the x axis. Above 1 the figure grows wider than the set width.';

  @override
  String get treeYScale => 'Vertical scale';

  @override
  String get treeYScaleInfo =>
      'Space between the rows on the y axis (1 = 1.7 times the label size).';

  @override
  String get treeBranchLengths => 'Branch lengths';

  @override
  String get treeShowLengths => 'Show as labels';

  @override
  String get treeSupportFont => 'Support text size';

  @override
  String get treeProportional => 'Size by support value';

  @override
  String get treeShowLabels => 'Show tip labels';

  @override
  String get treeTabFigure => 'Figure and export';

  @override
  String get treeScaleBar => 'Scale bar';

  @override
  String get treeShowScaleBar => 'Show';

  @override
  String get treeScaleBarLength => 'Scale bar length';

  @override
  String get treeScaleBarLengthInfo =>
      'In substitutions per site. 0 picks a round value about a fifth of the tree.';

  @override
  String get treeFigureLegend => 'Legend';

  @override
  String get treeLegendInFigure => 'Draw in the figure';

  @override
  String get treeFigureLegendInfo =>
      'Support colours, outgroup, shortened branches and row colours are drawn under the tree, so they go with the PNG and SVG.';

  @override
  String get treeAtTips => 'At the tips';

  @override
  String get treeLegendMarks => 'Marks';

  @override
  String get treeBranchColor => 'Branch colour';

  @override
  String get treeSupportPosition => 'Support position';

  @override
  String get treeAboveBranch => 'Above the branch';

  @override
  String get treeBelowBranch => 'Below the branch';

  @override
  String get treeSupportColor => 'Colour';

  @override
  String get treeFontFamily => 'Font';

  @override
  String get treeLabelColor => 'Label colour';

  @override
  String get treeLabelPosition => 'Label position';

  @override
  String get treeLabelGap => 'Distance from the tip';

  @override
  String get treeLabelGapInfo =>
      'Space between the end of the branch (or the aligned column) and the label.';

  @override
  String get treeColorApply => 'Apply the colours to';

  @override
  String get treeColorApplyRows => 'Row background';

  @override
  String get treeColorApplyBranches => 'Branches';

  @override
  String get treeColorApplyLabels => 'Labels';

  @override
  String get treeCladeColors => 'Clade colours';

  @override
  String get treeCladeColorsInfo =>
      'Select a clade by clicking its branch and give it a colour: its branches take it. Undo works here too.';

  @override
  String get treeColorClade => 'Colour clade';

  @override
  String get treeRemoveColor => 'Remove colour';

  @override
  String treeClearCladeColors(String n) {
    return 'Remove all ($n)';
  }

  @override
  String get treeBackground => 'Figure background';

  @override
  String get treeLegendPosition => 'Position';

  @override
  String get treeLegendBelow => 'Below the tree';

  @override
  String get treeLegendRight => 'Right of the tree';

  @override
  String get treeLegendStyle => 'Style';

  @override
  String get treeLegendBox => 'In a box';

  @override
  String get treeLegendTitle => 'Legend title';

  @override
  String get treeLegendTitleHint => 'e.g. Legend';

  @override
  String get treeLegendFont => 'Legend text size';

  @override
  String get treeFigureTitle => 'Figure title';

  @override
  String get treeFigureTitleHint => 'Optional, drawn above the tree';

  @override
  String get treeDpi => 'PNG resolution';

  @override
  String get treeDpiInfo =>
      'Dots per inch of the PNG (300 for most journals, 600 for line art). The DPI is written in the file. Very large figures are capped at 16 000 px on the longest side.';

  @override
  String get treePngBackground => 'PNG background';

  @override
  String get treeTransparent => 'Transparent';

  @override
  String get treePickColor => 'Choose a colour';

  @override
  String get concatBlockOutgroup => 'Which taxa root the tree (outgroup)';

  @override
  String get concatBlockOutgroupBody =>
      'Pick one or more taxa outside the group you study; the trees of steps 8 and 9 are rooted on them. Without an outgroup the tree viewer roots at the midpoint.';

  @override
  String get schemeOneModel => 'one model for everything';

  @override
  String schemePartitionsCount(String n) {
    return '$n partition(s)';
  }

  @override
  String get resCpu => 'CPU';

  @override
  String get resCpuInfo =>
      'CPU load of the machine, including MAFFT, trimAl, IQ-TREE and the other tools SPLACE runs.';

  @override
  String get resMemory => 'RAM';

  @override
  String resMemoryInfo(String app, String used, String total) {
    return 'SPLACE uses $app. Machine: $used of $total in use (the bar).';
  }

  @override
  String resMemoryInfoApp(String app) {
    return 'SPLACE uses $app.';
  }

  @override
  String get resDisk => 'Disk';

  @override
  String resFree(String v) {
    return '$v free';
  }

  @override
  String resDiskInfo(String free, String total) {
    return '$free free of $total on the drive of the project folder (the bar is the space in use).';
  }

  @override
  String get resDiskInfoNone => 'Space on the drive of the project folder.';

  @override
  String get resNetwork => 'Network';

  @override
  String resNetworkInfo(String total) {
    return 'Data SPLACE exchanges with NCBI, GBIF and BLAST, per second. Since the app opened: $total.';
  }

  @override
  String get methodsTitle => 'Methods text';

  @override
  String get methodsCopy => 'Copy text and references';

  @override
  String get methodsHelp =>
      'Written from what was done in this project (data, alignment, trimming, saturation, supermatrix and trees). Review it and adapt it to your manuscript.';

  @override
  String get citeReferences => 'References cited in the text';

  @override
  String get citeTools => 'Tools used by SPLACE';

  @override
  String get treeBranchStyle => 'Branch style';

  @override
  String get treeBranchSquare => 'Right angles';

  @override
  String get treeBranchRounded => 'Rounded corners';

  @override
  String get treeBranchSlanted => 'Slanted';

  @override
  String get treeZoom => 'Zoom (screen only)';

  @override
  String get treeZoomFit => 'Fit the window';

  @override
  String get treeZoomInfo =>
      'Enlarges or reduces the view of the figure on this screen. The exported figure keeps the width set in Figure and export, whatever the monitor.';

  @override
  String get treeTipDotsLabel => 'Tips';

  @override
  String get treeLabelCustom => 'Custom';

  @override
  String get treeMarkOutgroup => 'Grey box behind';

  @override
  String get treeRenamed => 'Edited labels';

  @override
  String get treeRenamedInfo =>
      'Click a tip and choose Edit label to write its own text (with italic and bold). Clear the text to go back to the built label.';

  @override
  String get treeRenameTip => 'Edit label';

  @override
  String treeRestoreNames(String n) {
    return 'Restore edited labels ($n)';
  }

  @override
  String get treeBarsByRank => 'Bars with taxon names';

  @override
  String get treeBarsByRankInfo =>
      'Draws a bar beside each group of tips that share the taxon of this rank (e.g. one bar per family), with its name. Taxa split in the tree get one bar per part.';

  @override
  String get treeHighlightOpacity => 'Box opacity';

  @override
  String get treeHighlightClade => 'Highlight clade';

  @override
  String get treeRemoveHighlight => 'Remove highlight';

  @override
  String get treeCladeLabel => 'Name clade';

  @override
  String treeClearHighlights(String n) {
    return 'Remove highlights ($n)';
  }

  @override
  String treeClearCladeLabels(String n) {
    return 'Remove clade names ($n)';
  }

  @override
  String get treeRemoveFromFigure => 'Remove from tree';

  @override
  String get treeFigureWidth => 'Figure width';

  @override
  String get treeFigureWidthInfo =>
      'Real width of the exported figure. Font sizes are in points, so a 9 pt label prints at 9 pt at this width. The height follows from the number of tips and the vertical scale.';

  @override
  String get treeWidthPreset => 'Journal width';

  @override
  String get treeWidthOneCol => '1 column (8.5 cm)';

  @override
  String get treeWidthOneHalf => '1.5 columns (12 cm)';

  @override
  String get treeWidthTwoCol => 'Full page (17.5 cm)';

  @override
  String get treeWidthCustom => 'Custom';

  @override
  String get treeMargin => 'Margin';

  @override
  String get treeTitleFont => 'Title size';

  @override
  String treeFigureSize(String w, String h, String pw, String ph, String dpi) {
    return 'Figure $w × $h cm · PNG $pw × $ph px at $dpi dpi';
  }

  @override
  String get treeTextHint => 'Text';

  @override
  String get treeTextEmptyHint => 'Leave it empty to remove it.';

  @override
  String get concatHowOpen => 'Click to see how the loci are put side by side';

  @override
  String get treeOpenEditor => 'Open the tree editor';

  @override
  String get treeEditorTitle => 'Tree editor';

  @override
  String get treeEditorHint =>
      'Preview of the figure at its real size. Click it, or the button above, to edit it in the full window.';

  @override
  String get treeCompare => 'Compare trees';

  @override
  String get treeModeEdit => 'Edit figure';

  @override
  String get treeDone => 'Done';

  @override
  String get treeExportMenu => 'Export';

  @override
  String get treeZoomIn => 'Zoom in (Ctrl + wheel)';

  @override
  String get treeZoomOut => 'Zoom out (Ctrl + wheel)';

  @override
  String get treeFitWidth => 'Fit width';

  @override
  String get treeFitPage => 'Whole figure';

  @override
  String get treeCanvasHint =>
      'Drag to move · wheel scrolls, Ctrl + wheel zooms · click selects, Shift + click selects the clade joining two tips · right click for options · Del removes · Esc clears';

  @override
  String get treeTabCompare => 'Compare';

  @override
  String get treeCompareWith => 'Compare with';

  @override
  String get treeUntangle => 'Untangle';

  @override
  String get treeUntangleText => 'Rotate the other tree';

  @override
  String get treeUntangleInfo =>
      'Rotates the clades of the other tree (its topology does not change) so each tip faces its pair in the tree being edited, with as few crossing lines as possible.';

  @override
  String get treeConflicts => 'Different clades';

  @override
  String get treeConflictsText => 'Mark them';

  @override
  String get treeConflictsInfo =>
      'Branches leading to a clade that the other tree does not have are drawn thicker, in this colour.';

  @override
  String get treeConflictColor => 'Colour of different clades';

  @override
  String get treeConnectors => 'Lines between the trees';

  @override
  String get treeConnectorsGroup => 'Row colours';

  @override
  String get treeConnectorsInfo =>
      'With rows coloured by a rank (Colours tab), each line takes the colour of its taxon; otherwise it is grey.';

  @override
  String get treeCompareGap => 'Space between the trees';

  @override
  String get treeCompareResult => 'Result';

  @override
  String treeCompareCommon(String n) {
    return '$n tips in both trees';
  }

  @override
  String treeCompareOnly(String a, String b) {
    return '$a only on the left · $b only on the right';
  }

  @override
  String treeCompareRf(String rf, String nrf) {
    return 'Robinson–Foulds distance $rf (normalised $nrf)';
  }

  @override
  String treeCompareShared(String n, String a, String b) {
    return '$n clades in common · $a only on the left · $b only on the right';
  }

  @override
  String get treeCompareRfInfo =>
      'Counts the clades found in one tree and not in the other (unrooted, over the tips both trees have). 0 means the same topology; the normalised value divides by the largest possible distance.';

  @override
  String get treeCompareLegend => 'Clade not found in the other tree';

  @override
  String get treeCompareTooFew =>
      'The trees share fewer than 4 tips: their clades cannot be compared.';

  @override
  String treeCompareError(String e) {
    return 'The other tree could not be read: $e';
  }

  @override
  String get treeCompareIntro =>
      'The tree being edited and the one chosen here are drawn face to face, one mirrored, with a line joining each tip to its pair. Styles (fonts, branches, support, colours) come from the other tabs.';

  @override
  String get treeSelectParent => 'Parent clade';

  @override
  String get treeLeft => 'Left';

  @override
  String treeDetailMoose(String n) {
    return 'MOOSE is choosing the models ($n partition(s))';
  }

  @override
  String get treeKindRbs => 'ML tree, rapid bootstrap (RBS)';

  @override
  String get treeKindSh => 'ML tree, SH-aLRT';

  @override
  String get treeKindEbg => 'ML tree, predicted bootstrap (EBG)';

  @override
  String raxmlOld(String version) {
    return 'RAxML-NG $version found. MOOSE, the fast search and the SH-aLRT, EBG and RBS supports need version 2.0 or newer (SPLACE ships 2.0.4).';
  }

  @override
  String get raxmlModelsHeading => 'Substitution models';

  @override
  String get raxmlMoose => 'MOOSE (RAxML-NG\'s own model selection)';

  @override
  String get raxmlMooseBody =>
      'One model per partition, chosen inside the same run. Recommended.';

  @override
  String get raxmlModelFinderBody =>
      'IQ-TREE picks the models and can merge partitions; RAxML-NG then builds the tree.';

  @override
  String get raxmlFixedModel => 'One fixed model';

  @override
  String get raxmlFixedBody =>
      'The same model for every partition (e.g. GTR+G).';

  @override
  String get raxmlCriterion => 'Criterion';

  @override
  String get raxmlSearch => 'Tree search';

  @override
  String get raxmlSearchAdaptive => 'Adaptive (default)';

  @override
  String get raxmlSearchFast => 'Fast (--fast)';

  @override
  String get raxmlSupportInfo =>
      'FBP, TBE and RBS come from bootstrap replicate trees. SH-aLRT and EBG are computed on the ML tree only and take seconds. FBP and RBS cannot be used together.';

  @override
  String get treeLayoutLabel => 'Layout';

  @override
  String get treeLayoutRect => 'Rectangular';

  @override
  String get treeLayoutCircle => 'Circular';

  @override
  String get treeLayoutInfo =>
      'Circular: the root in the centre and the tips around it, good for trees with many tips. Names follow the radius; taxon colours become a ring outside the names and clade bars become arcs.';

  @override
  String get treeArc => 'Opening';

  @override
  String get treeArcInfo =>
      'How much of the circle the tips take: 360° closes it, 180° makes a half circle (fan).';

  @override
  String get treeRotation => 'Rotation';

  @override
  String get treeSwapSides => 'Swap sides';

  @override
  String get treeSwapSidesText => 'This tree on the right';

  @override
  String get treeSwapSidesInfo =>
      'Puts the tree being edited on the right and the other one on the left.';

  @override
  String get iqAlrt => 'SH-aLRT test';

  @override
  String get iqAlrtInfo =>
      'A second, fast support value for each branch (shown as SH-aLRT/UFBoot). Branches with SH-aLRT ≥ 80 and UFBoot ≥ 95 are well supported.';

  @override
  String get iqAlrt1000 => '1000 replicates (recommended)';

  @override
  String get iqSearchHeading => 'Tree search';

  @override
  String get iqSearchInfo =>
      'The defaults suit most data. More runs or more iterations take longer but explore more trees, useful when results change between runs.';

  @override
  String get iqRuns => 'Independent runs';

  @override
  String get iqRunsInfo =>
      '--runs: repeats the whole search from different starting trees and keeps the best tree.';

  @override
  String get iqNstop => 'Stop after (iterations)';

  @override
  String get iqNstopInfo =>
      '-nstop: the search stops after this many iterations without a better tree (default 100). Use more for large or difficult data.';

  @override
  String get iqPers => 'Perturbation strength';

  @override
  String get iqPersInfo =>
      '-pers: how much the current tree is changed at each iteration (default 0.5). Smaller values help with short alignments.';

  @override
  String get iqSeedInfo =>
      'The same seed gives the same result again. Empty: a random seed (written in the log).';

  @override
  String get iqSeedRandom => 'random';

  @override
  String get iqBnni => 'Refine the UFBoot trees (--bnni)';

  @override
  String get iqBnniBody =>
      'Optimises each bootstrap tree with NNI moves. Recommended: it reduces overestimated support when the model does not fit the data perfectly.';

  @override
  String get bothTitle => 'Run IQ-TREE and RAxML-NG together';

  @override
  String get bothInfo =>
      'Both are maximum-likelihood methods. Comparing their trees (step 9, Compare trees) shows which clades do not depend on the program.';

  @override
  String bothBody(String iq, String rx) {
    return 'Both analyses start at the same time with the options chosen in the two cards above. The threads are split so the computer stays usable: IQ-TREE gets $iq and RAxML-NG $rx. Change the total in the settings (top bar).';
  }

  @override
  String bothThreads(String n) {
    return '$n threads';
  }

  @override
  String get runBoth => 'Run both';

  @override
  String get bothRunning => 'IQ-TREE and RAxML-NG';

  @override
  String bothRunningSub(String threads, String time) {
    return 'Side by side · threads $threads · $time';
  }

  @override
  String get satUnlocksStep7 =>
      'Run the test to continue: its suggestion (which codon positions to keep) is used to build the matrix in step 7, which opens afterwards.';

  @override
  String get satAllGenes => 'All genes together';

  @override
  String get satAllGenesInfo =>
      '\"All genes together\" joins the coding genes only for this test (the matrix itself is built in step 7). It is the row that decides the suggestion of which codon positions to remove from the whole matrix; each gene\'s row allows a choice per gene.';

  @override
  String get exampleLabel => 'Example';

  @override
  String get exampleNt =>
      'The same alignment in the three modes. Here the gap of species B falls inside codons 3 and 4 and shifts the reading frame.';

  @override
  String get exampleAa =>
      'The same genes translated: one letter per codon, so fewer and more conserved sites.';

  @override
  String get exampleCodon =>
      'The gap of species B is the whole codon 3: the reading frame is kept.';

  @override
  String exampleTrim(String n, String total) {
    return 'Trimming removes unreliable columns (mostly gaps, or sites without phylogenetic signal) before the tree is built. Here $n of $total columns are removed (×).';
  }

  @override
  String get exampleCodonKeepAll => 'All three codon positions are kept.';

  @override
  String exampleCodonKept(String positions, String pct) {
    return 'Kept: $positions positions ($pct% of the coding sites).';
  }

  @override
  String get exampleSatLabel => 'How to read the plots';

  @override
  String get exampleSatOk => 'Not saturated';

  @override
  String get exampleSatOkBody =>
      'Observed differences (p) keep rising with distance (d), close to the dashed y = x line.';

  @override
  String get exampleSatBad => 'Saturated';

  @override
  String get exampleSatBadBody =>
      'p levels off while d keeps growing: new changes overwrite old ones and the signal is lost.';

  @override
  String get plotTrend => 'trend';

  @override
  String exampleSpecies(String letter) {
    return 'Species $letter';
  }

  @override
  String get examplePlayAgain => 'Play again';

  @override
  String get trimStepOriginal => 'Alignment';

  @override
  String get trimStepMarked => 'Columns removed';

  @override
  String get trimStepResult => 'Trimmed alignment';

  @override
  String get treeSupportArrowLabel => 'Arrows';

  @override
  String get treeSupportArrow => 'Value with an arrow to the branch';

  @override
  String get treeSupportArrowInfo =>
      'Moves each support value a little away from the branches (up or down, as in Support position) and joins it to its node with a small arrow. Useful when values fall on short branches or on top of other lines.';

  @override
  String get orgOpen => 'Sample and marker';

  @override
  String get orgTitle => 'Sample and marker in the headers';

  @override
  String orgSubtitle(String files, String seqs) {
    return '$files FASTA file(s) · $seqs sequences';
  }

  @override
  String orgSummary(String samples, String markers, String pct) {
    return '$samples samples × $markers markers · $pct% filled';
  }

  @override
  String orgApply(String n) {
    return 'Create $n marker FASTA files';
  }

  @override
  String get orgIntro =>
      'Each header must say which sample (ID) and which marker the sequence belongs to. Choose the characters that separate the parts of the header, then click the parts below: once for the sample, twice for the marker, a third time to ignore it. SPLACE then writes one FASTA per marker holding every sample, ready to align.';

  @override
  String get orgSampleFrom => 'Sample (ID) from';

  @override
  String get orgSampleFromInfo =>
      'Part of the header, or the file name when each file holds one sample and the headers only name the marker.';

  @override
  String get orgSampleHeader => 'Header';

  @override
  String get orgSampleFile => 'File name';

  @override
  String get orgDelimitersInfo =>
      'Every selected character splits the header. Several can be on at once when files use different separators (ITV01-L1 and ITV02_L1).';

  @override
  String get orgSpace => 'space';

  @override
  String get orgOther => 'others';

  @override
  String get orgParts => 'Parts of the header (examples from different files)';

  @override
  String get orgPartsInfo =>
      'Positions are kept for every header. A part that reaches the end of the header is counted from the end, so headers with more or fewer parts still work.';

  @override
  String get orgLegendSample => 'Sample (ID)';

  @override
  String get orgLegendMarker => 'Marker';

  @override
  String get orgLegendIgnored => 'Ignored';

  @override
  String get orgAdvanced => 'Advanced: tags and regular expressions';

  @override
  String get orgAdvancedActive =>
      'A tag or regular expression is in use (see Advanced); clicking the parts is off.';

  @override
  String orgImplausible(String markers, String seqs) {
    return '$markers markers for $seqs sequences';
  }

  @override
  String get orgImplausibleBody =>
      'Almost every sequence became its own marker, so the marker part is probably wrong (for example the whole header). Pick the delimiter and click the part that names the marker.';

  @override
  String orgUnassigned(String n) {
    return '$n sequence(s) without a sample or a marker (left out)';
  }

  @override
  String get orgUnassignedLabel => 'Without sample or marker';

  @override
  String get orgVariations => 'Name variations';

  @override
  String get orgVariationsInfo =>
      'The same sample or marker can be written in slightly different ways. These options join them; you can also rename any sample in the list below.';

  @override
  String get orgIgnoreCase => 'Ignore upper and lower case';

  @override
  String get orgLeadingZeros => 'Ignore zeros at the start of numbers';

  @override
  String orgFileMismatch(String n) {
    return '$n file(s) whose name differs from the sample in its headers';
  }

  @override
  String orgRenameTo(String name) {
    return 'Call it $name';
  }

  @override
  String get orgRenameAll => 'Use the file names for these';

  @override
  String get orgUseFileNames => 'Take every sample from the file name';

  @override
  String orgMerged(String n) {
    return '$n sample(s) written in more than one way, joined:';
  }

  @override
  String orgSamples(String n) {
    return 'Samples ($n)';
  }

  @override
  String get orgSamplesInfo =>
      'Click a name to rename it (Enter to apply). Giving two samples the same name joins them.';

  @override
  String get orgFindSample => 'Find a sample';

  @override
  String orgSampleMarkers(String n, String total) {
    return '$n of $total markers';
  }

  @override
  String orgMarkers(String n) {
    return 'Markers ($n)';
  }

  @override
  String get orgMarkersInfo =>
      'How many samples each marker has, and what to do when a sample has more than one sequence of the same marker. Every marker goes to step 4, where you choose which to use.';

  @override
  String orgDuplicates(String n) {
    return '$n sample(s) with more than one sequence of a marker';
  }

  @override
  String get orgCoverage => 'Markers by number of samples';

  @override
  String orgCoverageBar(String markers, String samples) {
    return '$markers markers in $samples samples';
  }

  @override
  String get orgSamplesAxis => 'samples';

  @override
  String orgCoverageSummary(String complete, String half) {
    return '$complete in every sample · $half in at least half of them';
  }

  @override
  String get orgPending =>
      'Tell SPLACE which part of the headers is the sample and which the marker';

  @override
  String orgProblem(String problem) {
    return 'The split is not used: $problem';
  }

  @override
  String get orgPendingBody =>
      'Then one FASTA per marker is written with every sample, and the markers appear in step 4.';

  @override
  String get orgMatrix => 'Matrix';

  @override
  String get orgWritten => 'FASTA per marker';

  @override
  String get findMarker => 'Find a marker';

  @override
  String resultsShowing(String n, String total) {
    return 'Showing $n of $total';
  }

  @override
  String showMoreN(String n) {
    return 'Show $n more';
  }

  @override
  String get showAll => 'Show all';

  @override
  String get noMatches => 'Nothing matches.';

  @override
  String orgFormats(String n) {
    return '$n header formats';
  }

  @override
  String perLocusCount(String positions, String n) {
    return '$positions in $n markers';
  }

  @override
  String get satOptionalLoci =>
      'Optional for loci: step 7 is already open. The markers marked as protein-coding are tested by codon position; the others on all their sites.';

  @override
  String get satOptionalLociNoCoding =>
      'Optional for loci: step 7 is already open. No marker is protein-coding, so each locus is tested on all its sites (there are no codon positions): the test shows which loci are saturated.';

  @override
  String exampleTrimNone(String total) {
    return 'None of the $total columns is removed: a MAFFT alignment has no column made only of gaps. This option only matters for alignments edited by hand or subsets of sequences.';
  }

  @override
  String get exampleGapsRow => 'gaps';

  @override
  String get runFailed => 'The analysis failed.';

  @override
  String get runLastLines => 'Last lines of the program';

  @override
  String get runOpenLog => 'Open the log';

  @override
  String orgHeaderFormats(String n) {
    return 'Header formats ($n)';
  }

  @override
  String get orgHeaderFormatsInfo =>
      'The headers are written in more than one way (other delimiters, more or fewer parts). Each format shows its sequences and can have its own sample and marker parts: click the parts of its example. Every sequence then ends up named the same way: the sample in the header, one FASTA file per marker.';

  @override
  String orgFormatCount(String n, String files) {
    return '$n sequences in $files files';
  }

  @override
  String get orgFormatOwnRule => 'Own parts';

  @override
  String get orgFormatSkipped => 'Left out';

  @override
  String orgFormatShow(String n) {
    return 'Show the $n sequences';
  }

  @override
  String get orgFormatHide => 'Hide the sequences';

  @override
  String get orgFormatUseMain => 'Same parts as the main rule';

  @override
  String get orgFormatSkip => 'Leave these sequences out';

  @override
  String get orgFormatUse => 'Use these sequences';

  @override
  String orgNameStructures(String n) {
    return 'Sample name structures ($n)';
  }

  @override
  String get orgNameStructuresInfo =>
      'Sample names grouped by how they are written: letters as they are, each digit as #. Names of one collection written with more or fewer digits (ITV1054 and ITV01054) can be standardized to the structure you choose: zeros are added or removed on the left. Hover the button to see the renames.';

  @override
  String orgStructureCount(int n) {
    String _temp0 = intl.Intl.pluralLogic(
      n,
      locale: localeName,
      other: '$n samples',
      one: '1 sample',
    );
    return '$_temp0';
  }

  @override
  String orgStandardize(String structure) {
    return 'Standardize to $structure';
  }

  @override
  String get satAllSitesNonCoding => 'all sites (non-coding)';

  @override
  String schemeMixed(String coding, String other) {
    return '$coding coding markers split by codon; $other without codons stay whole';
  }

  @override
  String get exampleCodonGene => 'Gene';

  @override
  String get exampleCodonMatrix => 'In the matrix';

  @override
  String get mtTitle => 'Taxa × markers';

  @override
  String get mtIntro =>
      'As in SequenceMatrix: each cell is the number of bases the taxon has in that marker (darker green = more complete), \"—\" when it has no sequence. The header gives each alignment\'s length, how many taxa have it and its gaps. Hover a cell or a header for details; the same table is saved as matrix_overview.tsv.';

  @override
  String get mtTaxon => 'Taxon';

  @override
  String get mtLoci => 'Markers';

  @override
  String get mtBases => 'Bases';

  @override
  String get mtMissing => 'Missing';

  @override
  String get mtPartitions => 'Partitions';

  @override
  String get mtLegend =>
      'Cell: bases of the taxon in the marker · — no sequence · header: length · taxa with the marker · % gaps';

  @override
  String mtLocusInfo(
    String locus,
    String length,
    String n,
    String total,
    String gaps,
  ) {
    return '$locus: $length bp, in $n of $total taxa, $gaps% gaps';
  }

  @override
  String get mtMissingCell => 'no sequence (filled with ?)';

  @override
  String mtCellInfo(String bases, String length, String gaps) {
    return '$bases of $length bp, $gaps gaps';
  }

  @override
  String get searchGenomes => 'Genomes';

  @override
  String get searchMarkers => 'Individual markers';

  @override
  String get markerIncludeGenomes => 'Also genomes';

  @override
  String get markerGenomesInfo =>
      'Off: only records up to 10 kb (gene sequences such as barcodes). On: complete genomes and long fragments too, which also hold the markers (and are much larger to download).';

  @override
  String markerPick(String n) {
    return 'Markers (SynGenes) · $n selected';
  }

  @override
  String get markerPickInfo =>
      'Protein-coding genes and rRNAs of the SynGenes tables. Each marker is searched with every name SynGenes knows for it (gene symbols and product names, e.g. COI, COX1, cytochrome c oxidase subunit I), so records annotated in different ways are found. Records of the same species are joined into one sample, and species missing a marker get it as missing data.';

  @override
  String reviewMarkerCount(
    String marker,
    String selected,
    String found,
    String species,
  ) {
    return '$marker: $selected/$found records · $species species';
  }

  @override
  String get reviewOnePerSpecies => 'One per species and marker';

  @override
  String get reviewOnePerSpeciesInfo =>
      'Keeps, for each marker, the longest record of each species (verified records first). It is the selection the review starts with.';

  @override
  String get reviewColMarkers => 'Markers';

  @override
  String get onePerTaxon => 'One sample per species (individual markers)';

  @override
  String get onePerTaxonBody =>
      'Joins records of the same species (a COI from one voucher, a CYTB from another) into one row of the matrix: the header is family_species and, for each marker, the longest sequence of the species is kept. Species without a marker get it as missing data (?).';

  @override
  String onePerTaxonJoined(String n) {
    return '$n species have several records: they are joined into one sample each.';
  }

  @override
  String satOptionalPartial(String markers) {
    return 'Optional here: step 7 is already open. $markers came from individual-marker records with partial sequences (fragments such as barcodes, marked < or > in GenBank), so the test is a guide, not a requirement.';
  }

  @override
  String satPartialMixed(String markers) {
    return '$markers have partial sequences: the test is optional for them. Step 7 waits for the complete markers.';
  }

  @override
  String metaRecords(String n) {
    return '$n records';
  }

  @override
  String headerConflictBlock(String n) {
    return '$n header(s) repeat inside the same marker, so those sequences could not be told apart in the alignment or the tree. Add a field that tells them apart (Accession or Voucher), or turn on \"One sample per species\" to join them. The alignment stays locked until every header is unique.';
  }

  @override
  String get headerAddAccession => 'Add Accession';

  @override
  String get alignBlockedHeader => 'Locked: repeated FASTA headers (5a)';

  @override
  String recGenesCount(String n, String total) {
    return 'Saturated in $n of $total genes.';
  }

  @override
  String recBasedOnGene(String gene, String taxa, String length) {
    return 'Result of $gene ($taxa taxa, $length bp). Each position is tested separately.';
  }

  @override
  String recBasedOnGenes(String n) {
    return '$n coding genes, each tested on its own (nothing is concatenated for the test). A position counts as saturated when it is saturated in more than half of the genes.';
  }

  @override
  String rec7bPreview(String section) {
    return 'How \"$section\" in step 7b will look:';
  }

  @override
  String get frameExampleLabel => 'Example · keep the reading frame';

  @override
  String get frameExampleCaption =>
      'trimAl removes single columns. With this option, when trimAl removes a column (×), SPLACE removes the whole codon it belongs to (+): the gene keeps whole codons and the 1st, 2nd and 3rd positions stay in order (the coloured ruler).';

  @override
  String get frameStepOriginal => 'Codons 1-2-3';

  @override
  String get frameStepColumn => 'trimAl removes a column';

  @override
  String get frameStepCodon => 'ClipKIT --codon takes the whole codon';

  @override
  String get frameStepResult => 'Frame kept';

  @override
  String get trimCasesLabel => 'Example · trimming a coding gene, two ways';

  @override
  String get trimCasesStepMark => 'Columns marked';

  @override
  String get trimCaseColumns => 'Column by column (trimAl)';

  @override
  String get trimCaseCodons => 'Keep the reading frame (ClipKIT --codon)';

  @override
  String get trimCaseSelected => 'selected';

  @override
  String trimCaseFrameLost(String n) {
    return 'Frame lost: $n columns out of codon position';
  }

  @override
  String trimCaseFrameKept(String k, String total) {
    return 'Frame kept: $k of $total codons, positions 1-2-3 in order';
  }

  @override
  String get trimCaseFrameKeptCols => 'Frame kept: no column out of position';

  @override
  String get trimCasesCaption =>
      'Same alignment, two tools. trimAl removes the columns it marks (×) one by one, so the codon positions (coloured ruler) fall out of order (!). ClipKIT with --codon removes the whole codon of every column it marks (+): the 1st, 2nd and 3rd positions stay in order for the partitions and the saturation test.';

  @override
  String get subRecordsTable => 'Records table';

  @override
  String get subRecordsTableBody =>
      'Every record with its taxonomy and gene counts; edit or remove records here.';

  @override
  String get subFeatureTypes => 'Annotations to use (CDS, rRNA, tRNA)';

  @override
  String get subChooseMarkers => 'Choose the markers';

  @override
  String get subRunSaturation => 'Run the test';

  @override
  String get subBuildMatrix => 'Build the matrix';

  @override
  String get subTree => 'Tree and figure';

  @override
  String get viewCodons => 'Codons';

  @override
  String get viewDots => 'Dots';

  @override
  String get viewDotsInfo =>
      'Shows a dot in place of the bases in columns where every sequence has the same base (the first sequence keeps its letters), so the variable columns stand out, as in MEGA and other viewers. In the codon view, codons that are the same in every sequence become dots.';

  @override
  String codonsVarying(String n, String total, String nonsyn) {
    return '$n of $total codons vary · $nonsyn change the amino acid';
  }

  @override
  String get codonVariation => 'Codon variation';

  @override
  String viewerCodonN(String n) {
    return 'codon $n';
  }

  @override
  String codonColumnInfo(String codons, String aa) {
    return '$codons different codon(s) here, $aa amino acid(s)';
  }

  @override
  String get codonLegendAa => 'colour = amino acid';

  @override
  String get codonLegendSyn => 'synonymous change';

  @override
  String get codonLegendNonSyn => 'amino acid changes';

  @override
  String get mtResizeNames =>
      'Drag to widen the names; double-click to fit the longest one';

  @override
  String treeLegendUnmarked(String min) {
    return 'no mark: support < $min';
  }

  @override
  String get treeStrongFrom => 'Strong support from';

  @override
  String get treeModerateFrom => 'Moderate support from';

  @override
  String get treeSupportClassesInfo =>
      'Each class colours the points (and values) from its threshold up to the next one. Pick a preset — automatic for the kind of support, 2 to 5 classes, or cut at the tertiles of this tree\'s own values — or set the thresholds and colours yourself, adding or removing classes.';

  @override
  String get treeLegendSupport => 'Branch support';

  @override
  String get treeLegendSupportShow => 'In the legend';

  @override
  String get treeLegendSupportTitle => 'Support heading';

  @override
  String get badgeBrazil => 'Brazil only';

  @override
  String get salvePossiblyExtinct => 'Possibly extinct';

  @override
  String get salveSheetSubtitle =>
      'SALVE sheet (ICMBio) · national extinction-risk assessment';

  @override
  String get salveOpenSheet => 'Open the sheet in SALVE';

  @override
  String get salveWatersheds => 'Watersheds';

  @override
  String get salvePublished => 'Published';

  @override
  String get salveSituation => 'Sheet status';

  @override
  String get salveSynonyms => 'Former names';

  @override
  String get salveJustification => 'Justification';

  @override
  String get salveAssessors => 'Assessors';

  @override
  String get salveCitation => 'How to cite';

  @override
  String get colCriteria => 'Criteria';

  @override
  String get colCommonNames => 'Common names';

  @override
  String get colSalveGroup => 'SALVE group';

  @override
  String get colEndemic => 'Endemic to Brazil';

  @override
  String get colBiomes => 'Biomes';

  @override
  String get colStates => 'States';

  @override
  String get colAssessed => 'Assessed';

  @override
  String get colSheet => 'Sheet';

  @override
  String get geneNamesAllKnown => 'All gene names recognized';

  @override
  String geneReviewAgain(String n) {
    return 'Review again ($n)';
  }

  @override
  String get clipkitCodonMode => 'Coding markers: ClipKIT mode (--codon)';

  @override
  String get clipkitCodonInfo =>
      'To keep the reading frame, the coding markers are trimmed by ClipKIT with its --codon option: when one position of a codon should go, ClipKIT removes the whole codon (Steenwyk et al. 2020). Non-coding markers are trimmed by trimAl with the method above.';

  @override
  String get outgroupSuggesting =>
      'Measuring the distances between the alignments\' taxa to suggest an outgroup…';

  @override
  String get outgroupSuggestNone =>
      'Align the markers to get an outgroup suggested from the distances.';

  @override
  String get outgroupSuggestTitle => 'Suggested outgroup (alignment distances)';

  @override
  String get outgroupSuggestWeak =>
      'No clearly separate group: the most distant taxon is suggested';

  @override
  String get outgroupSuggestInfo =>
      'p-distances between every pair of taxa over the sites both have in each alignment, summed across markers; an average-linkage (UPGMA) tree of them; the smaller of the two groups joined last is the one farthest from the rest. The value beside each taxon is its mean distance to all the others. Check it against what you know of the group: a long branch (fast-evolving or contaminated sequence) can also look distant.';

  @override
  String get outgroupSuggestApplied => 'In use';

  @override
  String get outgroupSuggestUse => 'Use as outgroup';

  @override
  String outgroupSuggestStats(
    String between,
    String within,
    String loci,
    String taxa,
  ) {
    return 'Mean p-distance to the ingroup $between; within the ingroup $within · $loci alignments, $taxa taxa.';
  }

  @override
  String get outgroupSuggestNext => 'Next most distant';

  @override
  String get treeArrowShape => 'Arrow shape';

  @override
  String get treeArrowStraight => 'Straight';

  @override
  String get treeArrowCurved => 'Curved';

  @override
  String get treeArrowShapeInfo =>
      'Values never overlap: when a value\'s place is taken it moves away and an arrow points to its node. The arrows can be straight or curved.';

  @override
  String get treeSupportClasses => 'Support classes';

  @override
  String treeClassesAuto(String kind) {
    return 'Automatic ($kind)';
  }

  @override
  String treeClassesN(String n) {
    return '$n classes';
  }

  @override
  String get treeClassesTertiles => 'Tertiles of this tree';

  @override
  String get treeClassHidden => 'below \"show support from\"';

  @override
  String get treeClassRemove => 'Remove this class';

  @override
  String get treeClassAdd => 'Add a class';

  @override
  String get tipFieldConservation => 'Category (SALVE)';

  @override
  String get treeHeatStatus => 'Category';

  @override
  String get treeHeatStatusLegend => 'Extinction risk (SALVE / ICMBio)';

  @override
  String get treeHeatYes => 'yes';

  @override
  String get treeHeatNo => 'no';

  @override
  String get treeHeatmap => 'SALVE heatmap';

  @override
  String get treeHeatmapInfo =>
      'Coloured columns beside the labels from the species\' SALVE sheets: category in the IUCN colours, endemic to Brazil, and one column per biome. Rectangular layout only. The same information can go into the labels as text (label parts above).';

  @override
  String get satParamsTitle => 'Test parameters';

  @override
  String get satParamsBody =>
      'How the saturation of each marker is measured and when it is called moderate or saturated. Changes apply on the next run.';

  @override
  String get satParamMethod => 'Method';

  @override
  String get satParamMethodValue =>
      'PhyKIT saturation (Steenwyk et al. 2021; Philippe et al. 2011): slope of the uncorrected p-distance of every pair of taxa on their patristic distance in the tree';

  @override
  String get satParamDistance => 'Trees';

  @override
  String get satParamDistanceValue =>
      'IQ-TREE, one per gene (fast search) with the model below; the branch lengths of each codon subset are estimated on that topology (-te)';

  @override
  String get satParamSubsets => 'Sites tested';

  @override
  String get satParamSubsetsValue =>
      '1st, 2nd and 3rd codon positions, 1st+2nd and all sites (non-coding markers: all sites)';

  @override
  String get satParamCode => 'Genetic code';

  @override
  String get satParamFrame => 'Reading frame';

  @override
  String get satParamFrameValue =>
      'from the alignment (CDS start, or the frame with fewest internal stops)';

  @override
  String get satParamsThresholds => 'Thresholds and settings';

  @override
  String get satParamSlopeModerate => 'Saturation · moderate <';

  @override
  String get satParamSlopeSaturated => 'Saturation · saturated <';

  @override
  String get satParamSlopeInfo =>
      'PhyKIT\'s saturation value: slope of the p-distance on the patristic distance. 1 = no saturation; the closer to 0, the more saturated. The defaults (0.5 and 0.3) are heuristic: the 2nd codon positions usually sit near 1 and saturated 3rd positions well below 0.5.';

  @override
  String get satParamIss => 'Iss · saturated from';

  @override
  String get satParamIssInfo =>
      'Xia\'s index: mean site entropy over the entropy expected under full saturation. Near 1, the sites look random.';

  @override
  String get satParamUndefinedModerate => 'No correction · moderate from';

  @override
  String get satParamUndefinedSaturated => 'No correction · saturated from';

  @override
  String get satParamUndefinedInfo =>
      'Pairs whose corrected distance cannot be computed (beyond what the model corrects): a strong sign of saturation.';

  @override
  String get satParamLowDivergence => 'Low divergence · p below';

  @override
  String get satParamLowDivergenceInfo =>
      'With a mean p-distance this low, saturation does not matter and the subset is called OK.';

  @override
  String get satParamMinSites => 'Minimum sites';

  @override
  String get satParamMinPairs => 'Minimum pairs';

  @override
  String get satParamReplicates => 'Iss · replicates';

  @override
  String get satParamReplicatesInfo =>
      'Random sites drawn to estimate the entropy expected under full saturation, for each number of sequences. More replicates, steadier Iss, longer run.';

  @override
  String get satParamPlotPoints => 'Plot points';

  @override
  String get satParamPlotPointsInfo =>
      'Pairs of sequences drawn in the plots (a sample when there are more); the statistics use every pair.';

  @override
  String get satParamsReset => 'Restore defaults';

  @override
  String get satParamsStale =>
      'The parameters changed since the last run: run the test again to update the results.';

  @override
  String get colTrend => 'Population trend';

  @override
  String get colGeneration => 'Generation time';

  @override
  String get colMigratory => 'Migratory';

  @override
  String get colHabitats => 'Habitats';

  @override
  String get colThreats => 'Threats';

  @override
  String get colOfficialList => 'Official list';

  @override
  String get colUcs => 'Protected areas';

  @override
  String get salveReferences => 'References in the sheet';

  @override
  String get salveHistory => 'Assessment history';

  @override
  String get salveLists => 'Official lists and conventions';

  @override
  String get salveSecDistribution => 'Distribution';

  @override
  String get salveNational => 'In Brazil';

  @override
  String get salveGlobal => 'Worldwide';

  @override
  String get salveSecNatural => 'Natural history and population';

  @override
  String get salveDiet => 'Diet';

  @override
  String get salveReproduction => 'Reproduction';

  @override
  String get salveSexRatio => 'Sex ratio';

  @override
  String get salveNaturalHistory => 'Natural history';

  @override
  String get salveInteractions => 'Interactions with other species';

  @override
  String get salvePopulation => 'Population';

  @override
  String get salveGenetics => 'Genetic characteristics';

  @override
  String get salveSecThreats => 'Threats, uses, conservation and research';

  @override
  String get salveUses => 'Uses';

  @override
  String get salveActions => 'Conservation actions';

  @override
  String get salveConservation => 'Conservation';

  @override
  String get salveResearch => 'Research';

  @override
  String get trendDeclining => 'declining';

  @override
  String get trendStable => 'stable';

  @override
  String get trendIncreasing => 'increasing';

  @override
  String get trendUnknown => 'unknown';

  @override
  String get treeHeatFew => 'few';

  @override
  String get treeHeatMany => 'many';

  @override
  String trimNothingRemoved(String method, String n) {
    return '$method kept all $n columns: this alignment has no gappy region to remove, so nothing is grey (the kpic modes also remove uninformative sites)';
  }

  @override
  String get treeBranchWidth => 'Branch width';

  @override
  String get treeBranchWidthInfo =>
      'Width of the tree itself (branches), in cm. 0 = automatic: what is left beside labels, heatmap, bars and legend, never under 45 % of the figure — with many extras the figure grows wider instead of squeezing the tree. Set a width to fix it; the figure grows to fit the rest.';

  @override
  String get treeHeatCell => 'Heatmap cell';

  @override
  String get treeHeatCellInfo =>
      'Side of each heatmap cell, in pt. 0 = automatic (from the row height).';

  @override
  String get citeVersions => 'Versions used';

  @override
  String get satParamVersions => 'Programs';

  @override
  String get satParamModel => 'IQ-TREE model';

  @override
  String get satParamModelInfo =>
      'Model of the gene trees and of the branch lengths of each codon subset (IQ-TREE -m), e.g. GTR+G, HKY+G, GTR+F+I+G4.';

  @override
  String get satParamGaps => 'Gaps';

  @override
  String get satParamGapsSwitch => 'Ignore sites with gaps (-e)';

  @override
  String get satParamGapsInfo =>
      'PhyKIT --exclude_gaps: sites with a gap in a pair are left out of its p-distance.';

  @override
  String get phykitRequired =>
      'The saturation test uses PhyKIT, which is not installed yet. Install it here (needs Python 3) or with: pip install phykit';

  @override
  String get installPhykit => 'Install PhyKIT';

  @override
  String plotPhykit(String s) {
    return 'p × tree distance · PhyKIT $s';
  }

  @override
  String get plotPhykitBody =>
      'Each point is a pair of taxa: uncorrected p-distance against their patristic distance in the IQ-TREE tree of these sites. Points far below the diagonal hide repeated substitutions.';

  @override
  String get plotPositionsPatristic => 'Codon positions × tree distance';

  @override
  String get patristicDistance => 'Patristic distance';

  @override
  String get metricMeanPatristic => 'Mean patristic';

  @override
  String get metricPhykit => 'PhyKIT saturation';

  @override
  String get trimToolLabel => 'Trimming tool';

  @override
  String get trimToolTrimal => 'trimAl · column by column';

  @override
  String get trimToolClipkit => 'ClipKIT · whole codons';

  @override
  String get trimToolTrimalBody =>
      'trimAl trims every marker column by column. In coding genes this can cut codons apart: their reading frame is lost, they become one partition without codon positions and stay out of the saturation test.';

  @override
  String get trimToolClipkitBody =>
      'ClipKIT trims every marker. In the coding genes it runs with --codon: when one position of a codon should go, the whole codon goes, so the reading frame and the codon positions are kept. Non-coding markers are trimmed by ClipKIT without --codon. trimAl is not used.';

  @override
  String get mtShowBases => 'Bases';

  @override
  String get mtShowPercent => '% of the marker';

  @override
  String get mtOrderMatrix => 'Matrix order';

  @override
  String get mtOrderComplete => 'Most complete first';

  @override
  String get mtScaleComplete => 'Completeness';

  @override
  String get mtScaleMissing => 'Missing / gaps';

  @override
  String get treeMarkShape => 'Point shape';

  @override
  String get treeShapeCircle => 'Circle';

  @override
  String get treeShapeSquare => 'Square';

  @override
  String get treeShapeDiamond => 'Diamond';

  @override
  String get treeShapeTriangle => 'Triangle';

  @override
  String get treeBarWidth => 'Clade bar width';

  @override
  String get treeBarWidthInfo =>
      'Width of the bars beside clades, in pt. 0 = automatic (from the branch width).';

  @override
  String get treeBarFont => 'Clade label size';

  @override
  String get treeBarFontInfo =>
      'Size of the names beside the clade bars, in pt. 0 = the tip label size.';

  @override
  String get treeHeatGap => 'Heatmap column gap';

  @override
  String get trimToolTrimalCodon => 'trimAl · whole codons';

  @override
  String get trimToolTrimalCodonBody =>
      'trimAl trims every marker. Coding genes are read codon by codon and translated with the project\'s genetic code; trimAl trims that protein alignment and, with -backtrans, writes back the codons of the kept amino acids, so whole codons go or stay and the codon positions are kept. Non-coding markers are trimmed column by column.';

  @override
  String frameStepCodonTool(String tool) {
    return '$tool takes the whole codon';
  }

  @override
  String trimCaseCodonsTool(String tool) {
    return 'Keep the reading frame ($tool)';
  }

  @override
  String trimCasesCaptionTool(String tool) {
    return 'Same alignment, two ways. Column by column, trimAl removes the columns it marks (×) one by one, so the codon positions (coloured ruler) fall out of order (!). With $tool the whole codon of every marked column goes (+): the 1st, 2nd and 3rd positions stay in order for the partitions and the saturation test.';
  }
}
