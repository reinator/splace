import 'dart:async';

import 'package:flutter/foundation.dart';
import 'package:flutter/widgets.dart';
import 'package:flutter_localizations/flutter_localizations.dart';
import 'package:intl/intl.dart' as intl;

import 'app_localizations_en.dart';
import 'app_localizations_es.dart';
import 'app_localizations_pt.dart';

// ignore_for_file: type=lint

/// Callers can lookup localized strings with an instance of AppLocalizations
/// returned by `AppLocalizations.of(context)`.
///
/// Applications need to include `AppLocalizations.delegate()` in their app's
/// `localizationDelegates` list, and the locales they support in the app's
/// `supportedLocales` list. For example:
///
/// ```dart
/// import 'l10n/app_localizations.dart';
///
/// return MaterialApp(
///   localizationsDelegates: AppLocalizations.localizationsDelegates,
///   supportedLocales: AppLocalizations.supportedLocales,
///   home: MyApplicationHome(),
/// );
/// ```
///
/// ## Update pubspec.yaml
///
/// Please make sure to update your pubspec.yaml to include the following
/// packages:
///
/// ```yaml
/// dependencies:
///   # Internationalization support.
///   flutter_localizations:
///     sdk: flutter
///   intl: any # Use the pinned version from flutter_localizations
///
///   # Rest of dependencies
/// ```
///
/// ## iOS Applications
///
/// iOS applications define key application metadata, including supported
/// locales, in an Info.plist file that is built into the application bundle.
/// To configure the locales supported by your app, you’ll need to edit this
/// file.
///
/// First, open your project’s ios/Runner.xcworkspace Xcode workspace file.
/// Then, in the Project Navigator, open the Info.plist file under the Runner
/// project’s Runner folder.
///
/// Next, select the Information Property List item, select Add Item from the
/// Editor menu, then select Localizations from the pop-up menu.
///
/// Select and expand the newly-created Localizations item then, for each
/// locale your application supports, add a new item and select the locale
/// you wish to add from the pop-up menu in the Value field. This list should
/// be consistent with the languages listed in the AppLocalizations.supportedLocales
/// property.
abstract class AppLocalizations {
  AppLocalizations(String locale)
    : localeName = intl.Intl.canonicalizedLocale(locale.toString());

  final String localeName;

  static AppLocalizations of(BuildContext context) {
    return Localizations.of<AppLocalizations>(context, AppLocalizations)!;
  }

  static const LocalizationsDelegate<AppLocalizations> delegate =
      _AppLocalizationsDelegate();

  /// A list of this localizations delegate along with the default localizations
  /// delegates.
  ///
  /// Returns a list of localizations delegates containing this delegate along with
  /// GlobalMaterialLocalizations.delegate, GlobalCupertinoLocalizations.delegate,
  /// and GlobalWidgetsLocalizations.delegate.
  ///
  /// Additional delegates can be added by appending to this list in
  /// MaterialApp. This list does not have to be used at all if a custom list
  /// of delegates is preferred or required.
  static const List<LocalizationsDelegate<dynamic>> localizationsDelegates =
      <LocalizationsDelegate<dynamic>>[
        delegate,
        GlobalMaterialLocalizations.delegate,
        GlobalCupertinoLocalizations.delegate,
        GlobalWidgetsLocalizations.delegate,
      ];

  /// A list of this localizations delegate's supported locales.
  static const List<Locale> supportedLocales = <Locale>[
    Locale('en'),
    Locale('es'),
    Locale('pt'),
  ];

  /// No description provided for @navProgress.
  ///
  /// In en, this message translates to:
  /// **'Progress'**
  String get navProgress;

  /// No description provided for @navSteps.
  ///
  /// In en, this message translates to:
  /// **'Pipeline steps'**
  String get navSteps;

  /// No description provided for @logTitle.
  ///
  /// In en, this message translates to:
  /// **'Activity log'**
  String get logTitle;

  /// No description provided for @logEmpty.
  ///
  /// In en, this message translates to:
  /// **'No activity yet. Logs from fetching, alignment and analyses will appear here.'**
  String get logEmpty;

  /// No description provided for @logCopy.
  ///
  /// In en, this message translates to:
  /// **'Copy all'**
  String get logCopy;

  /// No description provided for @logDownload.
  ///
  /// In en, this message translates to:
  /// **'Download as .txt'**
  String get logDownload;

  /// No description provided for @logClear.
  ///
  /// In en, this message translates to:
  /// **'Clear'**
  String get logClear;

  /// No description provided for @logCollapse.
  ///
  /// In en, this message translates to:
  /// **'Hide log panel'**
  String get logCollapse;

  /// No description provided for @logWelcome.
  ///
  /// In en, this message translates to:
  /// **'SPLACE ready. Start by loading data in step 1.'**
  String get logWelcome;

  /// No description provided for @settingsTheme.
  ///
  /// In en, this message translates to:
  /// **'Theme'**
  String get settingsTheme;

  /// No description provided for @themeLight.
  ///
  /// In en, this message translates to:
  /// **'Light'**
  String get themeLight;

  /// No description provided for @themeDark.
  ///
  /// In en, this message translates to:
  /// **'Dark'**
  String get themeDark;

  /// No description provided for @themeAuto.
  ///
  /// In en, this message translates to:
  /// **'Auto'**
  String get themeAuto;

  /// No description provided for @settingsDensity.
  ///
  /// In en, this message translates to:
  /// **'Density'**
  String get settingsDensity;

  /// No description provided for @densityComfy.
  ///
  /// In en, this message translates to:
  /// **'Comfortable'**
  String get densityComfy;

  /// No description provided for @densityCompact.
  ///
  /// In en, this message translates to:
  /// **'Compact'**
  String get densityCompact;

  /// No description provided for @settingsFontSize.
  ///
  /// In en, this message translates to:
  /// **'Text size'**
  String get settingsFontSize;

  /// No description provided for @fontSizeReset.
  ///
  /// In en, this message translates to:
  /// **'Reset'**
  String get fontSizeReset;

  /// No description provided for @minimapLocked.
  ///
  /// In en, this message translates to:
  /// **'Complete the previous step to unlock'**
  String get minimapLocked;

  /// No description provided for @tooltipMenu.
  ///
  /// In en, this message translates to:
  /// **'Toggle pipeline navigator'**
  String get tooltipMenu;

  /// No description provided for @tooltipLogs.
  ///
  /// In en, this message translates to:
  /// **'Activity log'**
  String get tooltipLogs;

  /// No description provided for @tooltipSettings.
  ///
  /// In en, this message translates to:
  /// **'Display & language'**
  String get tooltipSettings;

  /// No description provided for @tooltipHelp.
  ///
  /// In en, this message translates to:
  /// **'Show welcome screen'**
  String get tooltipHelp;

  /// No description provided for @tooltipTools.
  ///
  /// In en, this message translates to:
  /// **'External tools and analysis settings'**
  String get tooltipTools;

  /// No description provided for @cancel.
  ///
  /// In en, this message translates to:
  /// **'Cancel'**
  String get cancel;

  /// No description provided for @done.
  ///
  /// In en, this message translates to:
  /// **'Done'**
  String get done;

  /// No description provided for @andMore.
  ///
  /// In en, this message translates to:
  /// **'…and {n} more'**
  String andMore(String n);

  /// No description provided for @splashSummary.
  ///
  /// In en, this message translates to:
  /// **'A complete toolkit for phylogenomics: fetch GenBank records or organize FASTA loci, align, test saturation, concatenate and build trees.'**
  String get splashSummary;

  /// No description provided for @splashStart.
  ///
  /// In en, this message translates to:
  /// **'Start App'**
  String get splashStart;

  /// No description provided for @splashStep1.
  ///
  /// In en, this message translates to:
  /// **'Load genomes or individual markers (COI, 12S, CYTB…) from NCBI, or local GenBank and FASTA files. [[SynGenes]] builds the marker searches from every synonym of each gene and standardizes the marker names for the next steps.'**
  String get splashStep1;

  /// No description provided for @splashStep2.
  ///
  /// In en, this message translates to:
  /// **'Inspect imported records, refine species names and recover taxonomy from GBIF, NCBI or SALVE with [[dataFishing]] before marker selection.'**
  String get splashStep2;

  /// No description provided for @splashStep3.
  ///
  /// In en, this message translates to:
  /// **'Choose which feature types (CDS, rRNA, tRNA) enter the analysis, or show SPLACE where FASTA headers keep the sample and locus names.'**
  String get splashStep3;

  /// No description provided for @splashStep4.
  ///
  /// In en, this message translates to:
  /// **'Pick the markers for the matrix and review missing and duplicated sequences in the gene presence and length heatmaps.'**
  String get splashStep4;

  /// No description provided for @splashStep5.
  ///
  /// In en, this message translates to:
  /// **'Align each marker with [[MAFFT]] and trim it with [[trimAl]] or [[ClipKIT]], keeping whole codons in protein-coding loci.'**
  String get splashStep5;

  /// No description provided for @splashStep6.
  ///
  /// In en, this message translates to:
  /// **'Test substitution saturation of each coding gene at each codon position with [[PhyKIT]] and decide which positions to keep.'**
  String get splashStep6;

  /// No description provided for @splashStep7.
  ///
  /// In en, this message translates to:
  /// **'Concatenate the aligned loci into a supermatrix and export NEXUS, PHYLIP and partition files, optionally split by codon position.'**
  String get splashStep7;

  /// No description provided for @splashStep8.
  ///
  /// In en, this message translates to:
  /// **'Estimate the maximum-likelihood phylogeny with [[IQ-TREE 3|IQ-TREE]] or [[RAxML-NG]], with branch support and optional per-gene trees, and delimit species with [[mPTP]], [[GMYC]], [[ASAP]] and [[ABGD]].'**
  String get splashStep8;

  /// No description provided for @step1.
  ///
  /// In en, this message translates to:
  /// **'Load Data'**
  String get step1;

  /// No description provided for @step2.
  ///
  /// In en, this message translates to:
  /// **'Loaded Records'**
  String get step2;

  /// No description provided for @step3.
  ///
  /// In en, this message translates to:
  /// **'Feature Types'**
  String get step3;

  /// No description provided for @step4.
  ///
  /// In en, this message translates to:
  /// **'Select Markers'**
  String get step4;

  /// No description provided for @step5.
  ///
  /// In en, this message translates to:
  /// **'Sequence Alignment'**
  String get step5;

  /// No description provided for @step6.
  ///
  /// In en, this message translates to:
  /// **'Saturation Test'**
  String get step6;

  /// No description provided for @step7.
  ///
  /// In en, this message translates to:
  /// **'Concatenate and NEXUS'**
  String get step7;

  /// No description provided for @step8.
  ///
  /// In en, this message translates to:
  /// **'Phylogenetic Inference'**
  String get step8;

  /// No description provided for @step1Intro.
  ///
  /// In en, this message translates to:
  /// **'Choose how to load data into SPLACE: GenBank records from NCBI or local files, or your own FASTA files.'**
  String get step1Intro;

  /// No description provided for @tipsKicker.
  ///
  /// In en, this message translates to:
  /// **'Tips'**
  String get tipsKicker;

  /// No description provided for @tipsTitle.
  ///
  /// In en, this message translates to:
  /// **'Before you start'**
  String get tipsTitle;

  /// No description provided for @tipsIntro.
  ///
  /// In en, this message translates to:
  /// **'Build the **taxonomy table** gradually and combine **complementary sources** when needed.'**
  String get tipsIntro;

  /// No description provided for @tipsExamplePrefix.
  ///
  /// In en, this message translates to:
  /// **'Each option below includes an'**
  String get tipsExamplePrefix;

  /// No description provided for @tipsExampleSuffix.
  ///
  /// In en, this message translates to:
  /// **'button to help you get started quickly.'**
  String get tipsExampleSuffix;

  /// No description provided for @tipSampleTitle.
  ///
  /// In en, this message translates to:
  /// **'One file per sample'**
  String get tipSampleTitle;

  /// No description provided for @tipSampleBody.
  ///
  /// In en, this message translates to:
  /// **'each header is a **locus** (e.g. NC_030596_Oryza.fasta with >psbA, >matK). SPLACE regroups them into **one file per locus**.'**
  String get tipSampleBody;

  /// No description provided for @tipLocusTitle.
  ///
  /// In en, this message translates to:
  /// **'One file per locus'**
  String get tipLocusTitle;

  /// No description provided for @tipLocusBody.
  ///
  /// In en, this message translates to:
  /// **'each header is a **sample** (e.g. COI.fasta). Already in the right shape.'**
  String get tipLocusBody;

  /// No description provided for @tipMixedTitle.
  ///
  /// In en, this message translates to:
  /// **'Mixed headers'**
  String get tipMixedTitle;

  /// No description provided for @tipMixedBody.
  ///
  /// In en, this message translates to:
  /// **'sample and locus in the **same header** (>Bufo_bufo_COI) or **NCBI tags** ([gene=COX1]).'**
  String get tipMixedBody;

  /// No description provided for @tipAlignedTitle.
  ///
  /// In en, this message translates to:
  /// **'Already aligned?'**
  String get tipAlignedTitle;

  /// No description provided for @tipAlignedBody.
  ///
  /// In en, this message translates to:
  /// **'You can **skip MAFFT** in step 4 and go straight to trimming, saturation and trees.'**
  String get tipAlignedBody;

  /// No description provided for @optionN.
  ///
  /// In en, this message translates to:
  /// **'Option {n}'**
  String optionN(String n);

  /// No description provided for @filesTitle.
  ///
  /// In en, this message translates to:
  /// **'Use local FASTA files'**
  String get filesTitle;

  /// No description provided for @filesInfo.
  ///
  /// In en, this message translates to:
  /// **'Import .fasta, .fas, .fa or .fna files that are already on your computer.'**
  String get filesInfo;

  /// No description provided for @filesBody.
  ///
  /// In en, this message translates to:
  /// **'Import your own .fasta, .fas, .fa or .fna files from disk.'**
  String get filesBody;

  /// No description provided for @filesDrop.
  ///
  /// In en, this message translates to:
  /// **'Drag and drop FASTA files or folders here'**
  String get filesDrop;

  /// No description provided for @filesBrowse.
  ///
  /// In en, this message translates to:
  /// **'or click to browse (.fasta, .fas, .fa, .fna)'**
  String get filesBrowse;

  /// No description provided for @folderTitle.
  ///
  /// In en, this message translates to:
  /// **'Import a whole folder'**
  String get folderTitle;

  /// No description provided for @folderInfo.
  ///
  /// In en, this message translates to:
  /// **'Reads every FASTA file inside the folder and its subfolders.'**
  String get folderInfo;

  /// No description provided for @folderBody.
  ///
  /// In en, this message translates to:
  /// **'Useful for hundreds of files, one per sample.'**
  String get folderBody;

  /// No description provided for @addFolder.
  ///
  /// In en, this message translates to:
  /// **'Choose folder'**
  String get addFolder;

  /// No description provided for @inputSummary.
  ///
  /// In en, this message translates to:
  /// **'{files} files · {sequences} sequences loaded'**
  String inputSummary(String files, String sequences);

  /// No description provided for @clearRecords.
  ///
  /// In en, this message translates to:
  /// **'Clear records'**
  String get clearRecords;

  /// No description provided for @step2Intro.
  ///
  /// In en, this message translates to:
  /// **'Every FASTA header is listed below with the sample and locus SPLACE reads from it. Adjust the rules in step 3 until every record has both.'**
  String get step2Intro;

  /// No description provided for @detectedLayout.
  ///
  /// In en, this message translates to:
  /// **'Detected layout'**
  String get detectedLayout;

  /// No description provided for @layoutSamplePerFile.
  ///
  /// In en, this message translates to:
  /// **'one file per sample — headers are loci'**
  String get layoutSamplePerFile;

  /// No description provided for @layoutLocusPerFile.
  ///
  /// In en, this message translates to:
  /// **'one file per locus — headers are samples'**
  String get layoutLocusPerFile;

  /// No description provided for @layoutMixedHeader.
  ///
  /// In en, this message translates to:
  /// **'sample and locus are both in the header'**
  String get layoutMixedHeader;

  /// No description provided for @layoutTaggedHeader.
  ///
  /// In en, this message translates to:
  /// **'locus taken from header tags ([gene=…])'**
  String get layoutTaggedHeader;

  /// No description provided for @layoutUnknown.
  ///
  /// In en, this message translates to:
  /// **'not recognized — choose the fields manually'**
  String get layoutUnknown;

  /// No description provided for @confidence.
  ///
  /// In en, this message translates to:
  /// **'confidence {percent}%'**
  String confidence(String percent);

  /// No description provided for @organizing.
  ///
  /// In en, this message translates to:
  /// **'Updating…'**
  String get organizing;

  /// No description provided for @tokenColumns.
  ///
  /// In en, this message translates to:
  /// **'Header parts (split on the locus delimiters)'**
  String get tokenColumns;

  /// No description provided for @tokenColumnInfo.
  ///
  /// In en, this message translates to:
  /// **'{distinct} distinct · {geneLike}% look like loci'**
  String tokenColumnInfo(String distinct, String geneLike);

  /// No description provided for @searchHeaders.
  ///
  /// In en, this message translates to:
  /// **'Filter headers or files'**
  String get searchHeaders;

  /// No description provided for @colFile.
  ///
  /// In en, this message translates to:
  /// **'File'**
  String get colFile;

  /// No description provided for @colHeader.
  ///
  /// In en, this message translates to:
  /// **'Header'**
  String get colHeader;

  /// No description provided for @colSample.
  ///
  /// In en, this message translates to:
  /// **'Sample'**
  String get colSample;

  /// No description provided for @colLocus.
  ///
  /// In en, this message translates to:
  /// **'Locus'**
  String get colLocus;

  /// No description provided for @colLength.
  ///
  /// In en, this message translates to:
  /// **'Length'**
  String get colLength;

  /// No description provided for @step3Intro.
  ///
  /// In en, this message translates to:
  /// **'Choose where the sample name and the locus name come from. The preview and the table in step 2 update as you type.'**
  String get step3Intro;

  /// No description provided for @fieldN.
  ///
  /// In en, this message translates to:
  /// **'Field {n}'**
  String fieldN(String n);

  /// No description provided for @sampleField.
  ///
  /// In en, this message translates to:
  /// **'Sample name'**
  String get sampleField;

  /// No description provided for @sampleFieldInfo.
  ///
  /// In en, this message translates to:
  /// **'Becomes the sequence name in every locus file, so it must be identical across loci for concatenation.'**
  String get sampleFieldInfo;

  /// No description provided for @locusField.
  ///
  /// In en, this message translates to:
  /// **'Locus name'**
  String get locusField;

  /// No description provided for @locusFieldInfo.
  ///
  /// In en, this message translates to:
  /// **'Defines which file each sequence goes to. Synonyms such as COX1/COI/CO1 can be merged.'**
  String get locusFieldInfo;

  /// No description provided for @sourceFileName.
  ///
  /// In en, this message translates to:
  /// **'File name'**
  String get sourceFileName;

  /// No description provided for @sourceHeader.
  ///
  /// In en, this message translates to:
  /// **'Header parts'**
  String get sourceHeader;

  /// No description provided for @sourceTag.
  ///
  /// In en, this message translates to:
  /// **'Header tag'**
  String get sourceTag;

  /// No description provided for @sourceRegex.
  ///
  /// In en, this message translates to:
  /// **'Regex'**
  String get sourceRegex;

  /// No description provided for @delimiters.
  ///
  /// In en, this message translates to:
  /// **'Delimiters'**
  String get delimiters;

  /// No description provided for @tokenStart.
  ///
  /// In en, this message translates to:
  /// **'First part'**
  String get tokenStart;

  /// No description provided for @tokenEnd.
  ///
  /// In en, this message translates to:
  /// **'End (blank = last)'**
  String get tokenEnd;

  /// No description provided for @tagKey.
  ///
  /// In en, this message translates to:
  /// **'Tag key'**
  String get tagKey;

  /// No description provided for @regexPattern.
  ///
  /// In en, this message translates to:
  /// **'Pattern (first group is used)'**
  String get regexPattern;

  /// No description provided for @regexOnFileName.
  ///
  /// In en, this message translates to:
  /// **'Apply to file name'**
  String get regexOnFileName;

  /// No description provided for @previewLabel.
  ///
  /// In en, this message translates to:
  /// **'Preview'**
  String get previewLabel;

  /// No description provided for @optionsKicker.
  ///
  /// In en, this message translates to:
  /// **'Options'**
  String get optionsKicker;

  /// No description provided for @duplicates.
  ///
  /// In en, this message translates to:
  /// **'Duplicates and names'**
  String get duplicates;

  /// No description provided for @duplicatesInfo.
  ///
  /// In en, this message translates to:
  /// **'What to do when the same sample has the same locus more than once (e.g. inverted-repeat copies in plastomes).'**
  String get duplicatesInfo;

  /// No description provided for @dupLongest.
  ///
  /// In en, this message translates to:
  /// **'Keep longest'**
  String get dupLongest;

  /// No description provided for @dupMostComplete.
  ///
  /// In en, this message translates to:
  /// **'Keep fewest N/gaps'**
  String get dupMostComplete;

  /// No description provided for @dupFirst.
  ///
  /// In en, this message translates to:
  /// **'Keep first'**
  String get dupFirst;

  /// No description provided for @dupAll.
  ///
  /// In en, this message translates to:
  /// **'Keep all (rename)'**
  String get dupAll;

  /// No description provided for @mergeSynonyms.
  ///
  /// In en, this message translates to:
  /// **'Merge locus synonyms (COX1 = COI = CO1…)'**
  String get mergeSynonyms;

  /// No description provided for @sanitizeNames.
  ///
  /// In en, this message translates to:
  /// **'Names safe for tree software'**
  String get sanitizeNames;

  /// No description provided for @selectAll.
  ///
  /// In en, this message translates to:
  /// **'Select All'**
  String get selectAll;

  /// No description provided for @selectNone.
  ///
  /// In en, this message translates to:
  /// **'None'**
  String get selectNone;

  /// No description provided for @sharedOnly.
  ///
  /// In en, this message translates to:
  /// **'Shared Only'**
  String get sharedOnly;

  /// No description provided for @missingTitle.
  ///
  /// In en, this message translates to:
  /// **'Missing markers'**
  String get missingTitle;

  /// No description provided for @missingIn.
  ///
  /// In en, this message translates to:
  /// **'missing in {n} sample(s)'**
  String missingIn(String n);

  /// No description provided for @duplicatesTitle.
  ///
  /// In en, this message translates to:
  /// **'Duplicate markers'**
  String get duplicatesTitle;

  /// No description provided for @duplicatedIn.
  ///
  /// In en, this message translates to:
  /// **'duplicated in {n} sample(s)'**
  String duplicatedIn(String n);

  /// No description provided for @colCopies.
  ///
  /// In en, this message translates to:
  /// **'Copies (sizes)'**
  String get colCopies;

  /// No description provided for @colUsed.
  ///
  /// In en, this message translates to:
  /// **'Used'**
  String get colUsed;

  /// No description provided for @unassignedTitle.
  ///
  /// In en, this message translates to:
  /// **'Records without sample or locus'**
  String get unassignedTitle;

  /// No description provided for @missingSample.
  ///
  /// In en, this message translates to:
  /// **'no sample'**
  String get missingSample;

  /// No description provided for @missingLocus.
  ///
  /// In en, this message translates to:
  /// **'no locus'**
  String get missingLocus;

  /// No description provided for @useLoci.
  ///
  /// In en, this message translates to:
  /// **'Use {n} markers'**
  String useLoci(String n);

  /// No description provided for @exportLoci.
  ///
  /// In en, this message translates to:
  /// **'Export FASTA per locus…'**
  String get exportLoci;

  /// No description provided for @heatmapTitle.
  ///
  /// In en, this message translates to:
  /// **'Gene Presence Heatmap'**
  String get heatmapTitle;

  /// No description provided for @heatmapAbsent.
  ///
  /// In en, this message translates to:
  /// **'Absent'**
  String get heatmapAbsent;

  /// No description provided for @heatmapPresent.
  ///
  /// In en, this message translates to:
  /// **'Present'**
  String get heatmapPresent;

  /// No description provided for @matrixSummary.
  ///
  /// In en, this message translates to:
  /// **'{loci} markers × {samples} samples · {occupancy}% filled. Click a marker name to include or exclude it.'**
  String matrixSummary(String loci, String samples, String occupancy);

  /// No description provided for @panelMafft.
  ///
  /// In en, this message translates to:
  /// **'MAFFT Settings'**
  String get panelMafft;

  /// No description provided for @panelMafftBody.
  ///
  /// In en, this message translates to:
  /// **'Each marker is aligned separately; several markers run in parallel.'**
  String get panelMafftBody;

  /// No description provided for @mafftStrategy.
  ///
  /// In en, this message translates to:
  /// **'Alignment Method'**
  String get mafftStrategy;

  /// No description provided for @directionLabel.
  ///
  /// In en, this message translates to:
  /// **'Sequence direction'**
  String get directionLabel;

  /// No description provided for @adjustDirection.
  ///
  /// In en, this message translates to:
  /// **'Fix reversed sequences'**
  String get adjustDirection;

  /// No description provided for @panelTrimal.
  ///
  /// In en, this message translates to:
  /// **'Trimming (trimAl / ClipKIT)'**
  String get panelTrimal;

  /// No description provided for @trimalMethod.
  ///
  /// In en, this message translates to:
  /// **'Trimming Mode'**
  String get trimalMethod;

  /// No description provided for @sequenceType.
  ///
  /// In en, this message translates to:
  /// **'Sequence Type'**
  String get sequenceType;

  /// No description provided for @codonAwareTrimming.
  ///
  /// In en, this message translates to:
  /// **'Keep the reading frame'**
  String get codonAwareTrimming;

  /// No description provided for @codonAwareHint.
  ///
  /// In en, this message translates to:
  /// **'Removing single columns from protein-coding alignments shifts the reading frame and invalidates codon partitions and the saturation test.'**
  String get codonAwareHint;

  /// No description provided for @panelRun.
  ///
  /// In en, this message translates to:
  /// **'Run Alignment'**
  String get panelRun;

  /// No description provided for @panelRunBody.
  ///
  /// In en, this message translates to:
  /// **'Align all selected markers with MAFFT. Output is saved to {dir}.'**
  String panelRunBody(String dir);

  /// No description provided for @runMafft.
  ///
  /// In en, this message translates to:
  /// **'Align (MAFFT)'**
  String get runMafft;

  /// No description provided for @runTrimal.
  ///
  /// In en, this message translates to:
  /// **'Trim (trimAl)'**
  String get runTrimal;

  /// No description provided for @useAsAlignment.
  ///
  /// In en, this message translates to:
  /// **'Input is already aligned'**
  String get useAsAlignment;

  /// No description provided for @colCoding.
  ///
  /// In en, this message translates to:
  /// **'Coding'**
  String get colCoding;

  /// No description provided for @colFrame.
  ///
  /// In en, this message translates to:
  /// **'Reading frame'**
  String get colFrame;

  /// No description provided for @colInclude.
  ///
  /// In en, this message translates to:
  /// **'Use this marker'**
  String get colInclude;

  /// No description provided for @statRaw.
  ///
  /// In en, this message translates to:
  /// **'Raw'**
  String get statRaw;

  /// No description provided for @statAligned.
  ///
  /// In en, this message translates to:
  /// **'Aligned'**
  String get statAligned;

  /// No description provided for @statTrimmed.
  ///
  /// In en, this message translates to:
  /// **'Trimmed'**
  String get statTrimmed;

  /// No description provided for @statusError.
  ///
  /// In en, this message translates to:
  /// **'Error'**
  String get statusError;

  /// No description provided for @saturationIntro.
  ///
  /// In en, this message translates to:
  /// **'Each marker is tested with [[PhyKIT]], method of [[Philippe 2011|]]: IQ-TREE infers the tree of the gene, then, for the 1st, 2nd and 3rd codon positions (and 1st+2nd, all sites), PhyKIT regresses the uncorrected p-distance of every pair of taxa on their distance in the tree. A value near 1 means no saturation; the closer to 0, the more the observed differences hide repeated substitutions.'**
  String get saturationIntro;

  /// No description provided for @saturationThresholds.
  ///
  /// In en, this message translates to:
  /// **'Heuristic thresholds, set in 5a (test parameters): moderate below the moderate limit of PhyKIT\'s saturation value, saturated below the saturated limit; too little divergence (mean p) counts as OK.'**
  String get saturationThresholds;

  /// No description provided for @runSaturation.
  ///
  /// In en, this message translates to:
  /// **'Run saturation test'**
  String get runSaturation;

  /// No description provided for @recommendation.
  ///
  /// In en, this message translates to:
  /// **'Recommendation'**
  String get recommendation;

  /// No description provided for @recKeepAll.
  ///
  /// In en, this message translates to:
  /// **'Keep all three codon positions. In step 6b, partition the matrix by codon position so each one gets its own substitution model.'**
  String get recKeepAll;

  /// No description provided for @recAminoAcids.
  ///
  /// In en, this message translates to:
  /// **'All three positions are saturated: at this divergence the nucleotides say little about the relationships. Align as amino acids (step 4b) or use slower markers.'**
  String get recAminoAcids;

  /// No description provided for @recUndetermined.
  ///
  /// In en, this message translates to:
  /// **'There are too few sites or sequence pairs to judge saturation. Keep all positions.'**
  String get recUndetermined;

  /// No description provided for @applyPerLocus.
  ///
  /// In en, this message translates to:
  /// **'Use each gene\'s own result'**
  String get applyPerLocus;

  /// No description provided for @concatenated.
  ///
  /// In en, this message translates to:
  /// **'Concatenated'**
  String get concatenated;

  /// No description provided for @saturationHeatmapTitle.
  ///
  /// In en, this message translates to:
  /// **'Saturation by locus and codon position'**
  String get saturationHeatmapTitle;

  /// No description provided for @selectLocusHint.
  ///
  /// In en, this message translates to:
  /// **'Hover a cell to see its values; click a row to see its plots.'**
  String get selectLocusHint;

  /// No description provided for @verdictOk.
  ///
  /// In en, this message translates to:
  /// **'OK'**
  String get verdictOk;

  /// No description provided for @verdictModerate.
  ///
  /// In en, this message translates to:
  /// **'Moderate'**
  String get verdictModerate;

  /// No description provided for @verdictSaturated.
  ///
  /// In en, this message translates to:
  /// **'Saturated'**
  String get verdictSaturated;

  /// No description provided for @verdictInsufficient.
  ///
  /// In en, this message translates to:
  /// **'Too little data'**
  String get verdictInsufficient;

  /// No description provided for @plotObservedVsCorrected.
  ///
  /// In en, this message translates to:
  /// **'Observed p vs {model} distance'**
  String plotObservedVsCorrected(String model);

  /// No description provided for @plotTsTv.
  ///
  /// In en, this message translates to:
  /// **'Transitions and transversions vs {model} distance'**
  String plotTsTv(String model);

  /// No description provided for @transitions.
  ///
  /// In en, this message translates to:
  /// **'Transitions'**
  String get transitions;

  /// No description provided for @transversions.
  ///
  /// In en, this message translates to:
  /// **'Transversions'**
  String get transversions;

  /// No description provided for @metricSites.
  ///
  /// In en, this message translates to:
  /// **'Sites'**
  String get metricSites;

  /// No description provided for @metricPairs.
  ///
  /// In en, this message translates to:
  /// **'Pairs'**
  String get metricPairs;

  /// No description provided for @metricMeanP.
  ///
  /// In en, this message translates to:
  /// **'Mean p'**
  String get metricMeanP;

  /// No description provided for @metricSlope.
  ///
  /// In en, this message translates to:
  /// **'Slope'**
  String get metricSlope;

  /// No description provided for @metricTsTv.
  ///
  /// In en, this message translates to:
  /// **'Ts/Tv'**
  String get metricTsTv;

  /// No description provided for @metricIss.
  ///
  /// In en, this message translates to:
  /// **'Iss'**
  String get metricIss;

  /// No description provided for @metricUndefined.
  ///
  /// In en, this message translates to:
  /// **'Uncorrectable'**
  String get metricUndefined;

  /// No description provided for @metricVerdict.
  ///
  /// In en, this message translates to:
  /// **'Verdict'**
  String get metricVerdict;

  /// No description provided for @frameLabel.
  ///
  /// In en, this message translates to:
  /// **'codons start at base {frame}'**
  String frameLabel(String frame);

  /// No description provided for @step7Intro.
  ///
  /// In en, this message translates to:
  /// **'Build a supermatrix from aligned genes and export for phylogenetic analysis.'**
  String get step7Intro;

  /// No description provided for @srcTrimmed.
  ///
  /// In en, this message translates to:
  /// **'Trimmed sequences'**
  String get srcTrimmed;

  /// No description provided for @recommendedTag.
  ///
  /// In en, this message translates to:
  /// **'(recommended)'**
  String get recommendedTag;

  /// No description provided for @srcTrimmedBody.
  ///
  /// In en, this message translates to:
  /// **'Use trimAl output — removes poorly aligned positions.'**
  String get srcTrimmedBody;

  /// No description provided for @srcAligned.
  ///
  /// In en, this message translates to:
  /// **'Aligned sequences'**
  String get srcAligned;

  /// No description provided for @srcAlignedBody.
  ///
  /// In en, this message translates to:
  /// **'Use MAFFT output without trimming.'**
  String get srcAlignedBody;

  /// No description provided for @partitionScheme.
  ///
  /// In en, this message translates to:
  /// **'Partitions'**
  String get partitionScheme;

  /// No description provided for @schemeNone.
  ///
  /// In en, this message translates to:
  /// **'None'**
  String get schemeNone;

  /// No description provided for @schemeByLocus.
  ///
  /// In en, this message translates to:
  /// **'By locus'**
  String get schemeByLocus;

  /// No description provided for @schemeByLocusCodon.
  ///
  /// In en, this message translates to:
  /// **'By locus and codon position'**
  String get schemeByLocusCodon;

  /// No description provided for @schemeByCodon.
  ///
  /// In en, this message translates to:
  /// **'By codon position (all loci)'**
  String get schemeByCodon;

  /// No description provided for @codonExclusion.
  ///
  /// In en, this message translates to:
  /// **'Codon positions'**
  String get codonExclusion;

  /// No description provided for @codonExclusionInfo.
  ///
  /// In en, this message translates to:
  /// **'Excluded positions are removed from coding loci only. The saturation test (step 5) can fill this in.'**
  String get codonExclusionInfo;

  /// No description provided for @excludePosition.
  ///
  /// In en, this message translates to:
  /// **'Exclude position {n}'**
  String excludePosition(String n);

  /// No description provided for @perLocusExclusions.
  ///
  /// In en, this message translates to:
  /// **'Per locus: {list}'**
  String perLocusExclusions(String list);

  /// No description provided for @clearExclusions.
  ///
  /// In en, this message translates to:
  /// **'Clear per-locus exclusions'**
  String get clearExclusions;

  /// No description provided for @buildMatrix.
  ///
  /// In en, this message translates to:
  /// **'Build supermatrix'**
  String get buildMatrix;

  /// No description provided for @matrixReady.
  ///
  /// In en, this message translates to:
  /// **'Supermatrix ready'**
  String get matrixReady;

  /// No description provided for @matrixStats.
  ///
  /// In en, this message translates to:
  /// **'{taxa} taxa × {sites} sites · {partitions} partitions · {missing}% missing'**
  String matrixStats(
    String taxa,
    String sites,
    String partitions,
    String missing,
  );

  /// No description provided for @excludedLabel.
  ///
  /// In en, this message translates to:
  /// **'Excluded'**
  String get excludedLabel;

  /// No description provided for @outputFiles.
  ///
  /// In en, this message translates to:
  /// **'Output files'**
  String get outputFiles;

  /// No description provided for @openFolder.
  ///
  /// In en, this message translates to:
  /// **'Open folder'**
  String get openFolder;

  /// No description provided for @partitionPreview.
  ///
  /// In en, this message translates to:
  /// **'RAxML-NG partitions'**
  String get partitionPreview;

  /// No description provided for @step8Intro.
  ///
  /// In en, this message translates to:
  /// **'Maximum-likelihood tree estimation from the concatenated supermatrix with IQ-TREE3 or RAxML-NG.'**
  String get step8Intro;

  /// No description provided for @outgroup.
  ///
  /// In en, this message translates to:
  /// **'Outgroup'**
  String get outgroup;

  /// No description provided for @outgroupInfo.
  ///
  /// In en, this message translates to:
  /// **'Selected species are written to the NEXUS file and used to root IQ-TREE and RAxML-NG trees.'**
  String get outgroupInfo;

  /// No description provided for @outgroupHint.
  ///
  /// In en, this message translates to:
  /// **'Search species…'**
  String get outgroupHint;

  /// No description provided for @engineN.
  ///
  /// In en, this message translates to:
  /// **'Engine {n}'**
  String engineN(String n);

  /// No description provided for @iqtreeInfo.
  ///
  /// In en, this message translates to:
  /// **'ModelFinder chooses the models (and can merge partitions); ultrafast bootstrap and SH-aLRT give support values.'**
  String get iqtreeInfo;

  /// No description provided for @raxmlInfo.
  ///
  /// In en, this message translates to:
  /// **'RAxML-NG 2 chooses the models itself (MOOSE), adapts the search to how hard the alignment is and offers fast support measures. Like IQ-TREE it is a maximum-likelihood program (not Bayesian).'**
  String get raxmlInfo;

  /// No description provided for @model.
  ///
  /// In en, this message translates to:
  /// **'Substitution Model'**
  String get model;

  /// No description provided for @ufboot.
  ///
  /// In en, this message translates to:
  /// **'Ultrafast bootstrap (-B)'**
  String get ufboot;

  /// No description provided for @alrt.
  ///
  /// In en, this message translates to:
  /// **'SH-aLRT (--alrt)'**
  String get alrt;

  /// No description provided for @runIqtree.
  ///
  /// In en, this message translates to:
  /// **'Run IQ-TREE3'**
  String get runIqtree;

  /// No description provided for @raxmlModelFinder.
  ///
  /// In en, this message translates to:
  /// **'ModelFinder (IQ-TREE, RAxML-compatible models)'**
  String get raxmlModelFinder;

  /// No description provided for @bootstraps.
  ///
  /// In en, this message translates to:
  /// **'Bootstrap trees'**
  String get bootstraps;

  /// No description provided for @supportMetrics.
  ///
  /// In en, this message translates to:
  /// **'Support'**
  String get supportMetrics;

  /// No description provided for @seed.
  ///
  /// In en, this message translates to:
  /// **'Seed'**
  String get seed;

  /// No description provided for @runRaxml.
  ///
  /// In en, this message translates to:
  /// **'Run RAxML-NG'**
  String get runRaxml;

  /// No description provided for @runs.
  ///
  /// In en, this message translates to:
  /// **'Runs'**
  String get runs;

  /// No description provided for @noRuns.
  ///
  /// In en, this message translates to:
  /// **'No trees yet.'**
  String get noRuns;

  /// No description provided for @toolMissing.
  ///
  /// In en, this message translates to:
  /// **'{tool} not found'**
  String toolMissing(String tool);

  /// No description provided for @toolsTitle.
  ///
  /// In en, this message translates to:
  /// **'External tools and analysis settings'**
  String get toolsTitle;

  /// No description provided for @toolsBody.
  ///
  /// In en, this message translates to:
  /// **'SPLACE looks for bundled binaries, then PATH, then WSL. Set a custom path if needed.'**
  String get toolsBody;

  /// No description provided for @tools.
  ///
  /// In en, this message translates to:
  /// **'Tools'**
  String get tools;

  /// No description provided for @redetect.
  ///
  /// In en, this message translates to:
  /// **'Detect again'**
  String get redetect;

  /// No description provided for @notFound.
  ///
  /// In en, this message translates to:
  /// **'Not found'**
  String get notFound;

  /// No description provided for @overridePath.
  ///
  /// In en, this message translates to:
  /// **'Custom path (or wsl:/path/inside/wsl)'**
  String get overridePath;

  /// No description provided for @browse.
  ///
  /// In en, this message translates to:
  /// **'Browse…'**
  String get browse;

  /// No description provided for @analysisSettings.
  ///
  /// In en, this message translates to:
  /// **'Analysis'**
  String get analysisSettings;

  /// No description provided for @threads.
  ///
  /// In en, this message translates to:
  /// **'Threads: {n}'**
  String threads(String n);

  /// No description provided for @geneticCode.
  ///
  /// In en, this message translates to:
  /// **'Genetic code'**
  String get geneticCode;

  /// No description provided for @projectFolder.
  ///
  /// In en, this message translates to:
  /// **'Project folder'**
  String get projectFolder;

  /// No description provided for @change.
  ///
  /// In en, this message translates to:
  /// **'Change…'**
  String get change;

  /// No description provided for @tipComplete.
  ///
  /// In en, this message translates to:
  /// **'Prefer **complete genomes** whenever possible because they usually improve downstream analyses.'**
  String get tipComplete;

  /// No description provided for @tipMultipleSearches.
  ///
  /// In en, this message translates to:
  /// **'Run **multiple taxonomic searches** to add **ingroups** and **outgroups**. New results are appended to the same taxonomy table.'**
  String get tipMultipleSearches;

  /// No description provided for @tipCombine.
  ///
  /// In en, this message translates to:
  /// **'You can **combine** taxonomic searches, accession imports, and local **.gb, .gbk or .genbank** files in the same dataset.'**
  String get tipCombine;

  /// No description provided for @example.
  ///
  /// In en, this message translates to:
  /// **'Example'**
  String get example;

  /// No description provided for @taxTitle.
  ///
  /// In en, this message translates to:
  /// **'Search by taxonomic name'**
  String get taxTitle;

  /// No description provided for @taxInfo.
  ///
  /// In en, this message translates to:
  /// **'Search NCBI by taxonomic name, genome type, database and completeness before importing the chosen genomes.'**
  String get taxInfo;

  /// No description provided for @taxBody.
  ///
  /// In en, this message translates to:
  /// **'Search NCBI and download matching GenBank records automatically into the tool.'**
  String get taxBody;

  /// No description provided for @taxExampleTooltip.
  ///
  /// In en, this message translates to:
  /// **'Fill the taxonomic search with Carangidae, mitochondrion, complete genomes and RefSeq only.'**
  String get taxExampleTooltip;

  /// No description provided for @taxApiHint.
  ///
  /// In en, this message translates to:
  /// **'Add an NCBI API key to speed up imports (key button in the top bar).'**
  String get taxApiHint;

  /// No description provided for @taxName.
  ///
  /// In en, this message translates to:
  /// **'Taxonomic name'**
  String get taxName;

  /// No description provided for @taxOrganelle.
  ///
  /// In en, this message translates to:
  /// **'Genome type'**
  String get taxOrganelle;

  /// No description provided for @organelleMito.
  ///
  /// In en, this message translates to:
  /// **'Mitochondrion'**
  String get organelleMito;

  /// No description provided for @organelleChloro.
  ///
  /// In en, this message translates to:
  /// **'Chloroplast'**
  String get organelleChloro;

  /// No description provided for @taxScope.
  ///
  /// In en, this message translates to:
  /// **'Genome scope'**
  String get taxScope;

  /// No description provided for @taxCompleteOnly.
  ///
  /// In en, this message translates to:
  /// **'Complete only'**
  String get taxCompleteOnly;

  /// No description provided for @taxScopeInfo.
  ///
  /// In en, this message translates to:
  /// **'Complete only searches NCBI for records with “complete genome” in the title. When off, records with “partial genome” in the title are included too.'**
  String get taxScopeInfo;

  /// No description provided for @scopeComplete.
  ///
  /// In en, this message translates to:
  /// **'Complete genomes'**
  String get scopeComplete;

  /// No description provided for @scopePartial.
  ///
  /// In en, this message translates to:
  /// **'Partial genomes'**
  String get scopePartial;

  /// No description provided for @taxDatabase.
  ///
  /// In en, this message translates to:
  /// **'Database'**
  String get taxDatabase;

  /// No description provided for @taxRefseqOnly.
  ///
  /// In en, this message translates to:
  /// **'RefSeq only'**
  String get taxRefseqOnly;

  /// No description provided for @taxQueryPreview.
  ///
  /// In en, this message translates to:
  /// **'Query preview'**
  String get taxQueryPreview;

  /// No description provided for @search.
  ///
  /// In en, this message translates to:
  /// **'Search'**
  String get search;

  /// No description provided for @accTitle.
  ///
  /// In en, this message translates to:
  /// **'Paste accession numbers'**
  String get accTitle;

  /// No description provided for @accInfo.
  ///
  /// In en, this message translates to:
  /// **'Paste one or more known accessions and fetch those exact GenBank records from NCBI.'**
  String get accInfo;

  /// No description provided for @accBody.
  ///
  /// In en, this message translates to:
  /// **'Enter accession numbers separated by commas, spaces or line breaks.'**
  String get accBody;

  /// No description provided for @accExampleTooltip.
  ///
  /// In en, this message translates to:
  /// **'Load example accessions from Balcázar et al. (2026)'**
  String get accExampleTooltip;

  /// No description provided for @accFetchTooltip.
  ///
  /// In en, this message translates to:
  /// **'Fetch GenBank records from NCBI for the entered accession numbers'**
  String get accFetchTooltip;

  /// No description provided for @fetch.
  ///
  /// In en, this message translates to:
  /// **'Fetch'**
  String get fetch;

  /// No description provided for @gbFilesTitle.
  ///
  /// In en, this message translates to:
  /// **'Use local GenBank files'**
  String get gbFilesTitle;

  /// No description provided for @gbFilesInfo.
  ///
  /// In en, this message translates to:
  /// **'Import local .gb, .gbk or .genbank files when the records are already on your computer.'**
  String get gbFilesInfo;

  /// No description provided for @gbFilesBody.
  ///
  /// In en, this message translates to:
  /// **'Import your own .gb, .gbk or .genbank files from disk.'**
  String get gbFilesBody;

  /// No description provided for @gbFilesDrop.
  ///
  /// In en, this message translates to:
  /// **'Drag and drop GenBank files or folders here'**
  String get gbFilesDrop;

  /// No description provided for @gbFilesBrowse.
  ///
  /// In en, this message translates to:
  /// **'or click to browse (.gb, .gbk, .genbank)'**
  String get gbFilesBrowse;

  /// No description provided for @gbSummary.
  ///
  /// In en, this message translates to:
  /// **'{n} GenBank records loaded'**
  String gbSummary(String n);

  /// No description provided for @clearFasta.
  ///
  /// In en, this message translates to:
  /// **'Clear FASTA files'**
  String get clearFasta;

  /// No description provided for @close.
  ///
  /// In en, this message translates to:
  /// **'Close'**
  String get close;

  /// No description provided for @apiKeyTitle.
  ///
  /// In en, this message translates to:
  /// **'NCBI API key'**
  String get apiKeyTitle;

  /// No description provided for @apiKeyHelper.
  ///
  /// In en, this message translates to:
  /// **'Used by Option 1 and Option 2 for NCBI requests. Without an API key, SPLACE runs at about 3 requests per second. With a saved valid key, SPLACE uses up to 10 requests per second.'**
  String get apiKeyHelper;

  /// No description provided for @apiKeySave.
  ///
  /// In en, this message translates to:
  /// **'Save key'**
  String get apiKeySave;

  /// No description provided for @apiKeyClear.
  ///
  /// In en, this message translates to:
  /// **'Clear key'**
  String get apiKeyClear;

  /// No description provided for @apiKeyGetLink.
  ///
  /// In en, this message translates to:
  /// **'Get an API key from your NCBI account settings'**
  String get apiKeyGetLink;

  /// No description provided for @apiKeySaved.
  ///
  /// In en, this message translates to:
  /// **'API key saved. Option 1 and Option 2 now use up to 10 requests per second.'**
  String get apiKeySaved;

  /// No description provided for @apiKeyInvalid.
  ///
  /// In en, this message translates to:
  /// **'The API key looks invalid. Check it before saving.'**
  String get apiKeyInvalid;

  /// No description provided for @apiKeyNone.
  ///
  /// In en, this message translates to:
  /// **'No API key saved. Option 1 and Option 2 are using up to 3 requests per second.'**
  String get apiKeyNone;

  /// No description provided for @apiClearTitle.
  ///
  /// In en, this message translates to:
  /// **'Remove saved API key?'**
  String get apiClearTitle;

  /// No description provided for @apiClearBody.
  ///
  /// In en, this message translates to:
  /// **'Option 1 and Option 2 will go back to the default NCBI request rate.'**
  String get apiClearBody;

  /// No description provided for @reviewTitle.
  ///
  /// In en, this message translates to:
  /// **'Review genomes before importing'**
  String get reviewTitle;

  /// No description provided for @reviewSubtitle.
  ///
  /// In en, this message translates to:
  /// **'Choose which genomes should be imported. All switches start enabled.'**
  String get reviewSubtitle;

  /// No description provided for @reviewDetailTaxon.
  ///
  /// In en, this message translates to:
  /// **'Taxonomic query'**
  String get reviewDetailTaxon;

  /// No description provided for @reviewDetailGenomeType.
  ///
  /// In en, this message translates to:
  /// **'Genome type'**
  String get reviewDetailGenomeType;

  /// No description provided for @reviewDetailSearchTypes.
  ///
  /// In en, this message translates to:
  /// **'Search types'**
  String get reviewDetailSearchTypes;

  /// No description provided for @reviewDetailRefseq.
  ///
  /// In en, this message translates to:
  /// **'RefSeq'**
  String get reviewDetailRefseq;

  /// No description provided for @reviewRefseqOn.
  ///
  /// In en, this message translates to:
  /// **'Used'**
  String get reviewRefseqOn;

  /// No description provided for @reviewRefseqOff.
  ///
  /// In en, this message translates to:
  /// **'Not used'**
  String get reviewRefseqOff;

  /// No description provided for @reviewSelectAll.
  ///
  /// In en, this message translates to:
  /// **'Enable all'**
  String get reviewSelectAll;

  /// No description provided for @reviewSelectNone.
  ///
  /// In en, this message translates to:
  /// **'Disable all'**
  String get reviewSelectNone;

  /// No description provided for @reviewDisableUnverified.
  ///
  /// In en, this message translates to:
  /// **'Disable UNVERIFIED'**
  String get reviewDisableUnverified;

  /// No description provided for @reviewOnlyComplete.
  ///
  /// In en, this message translates to:
  /// **'Only complete'**
  String get reviewOnlyComplete;

  /// No description provided for @reviewOnlyPartial.
  ///
  /// In en, this message translates to:
  /// **'Only partial'**
  String get reviewOnlyPartial;

  /// No description provided for @reviewColImport.
  ///
  /// In en, this message translates to:
  /// **'Import'**
  String get reviewColImport;

  /// No description provided for @reviewColAccession.
  ///
  /// In en, this message translates to:
  /// **'Accession'**
  String get reviewColAccession;

  /// No description provided for @reviewColDescription.
  ///
  /// In en, this message translates to:
  /// **'Genome description'**
  String get reviewColDescription;

  /// No description provided for @reviewColScope.
  ///
  /// In en, this message translates to:
  /// **'Classification'**
  String get reviewColScope;

  /// No description provided for @reviewImport.
  ///
  /// In en, this message translates to:
  /// **'Import selected genomes'**
  String get reviewImport;

  /// No description provided for @reviewTotal.
  ///
  /// In en, this message translates to:
  /// **'Total'**
  String get reviewTotal;

  /// No description provided for @reviewSelected.
  ///
  /// In en, this message translates to:
  /// **'Selected'**
  String get reviewSelected;

  /// No description provided for @reviewComplete.
  ///
  /// In en, this message translates to:
  /// **'Complete'**
  String get reviewComplete;

  /// No description provided for @reviewPartial.
  ///
  /// In en, this message translates to:
  /// **'Partial'**
  String get reviewPartial;

  /// No description provided for @reviewMixed.
  ///
  /// In en, this message translates to:
  /// **'Mixed title'**
  String get reviewMixed;

  /// No description provided for @reviewUnknown.
  ///
  /// In en, this message translates to:
  /// **'Unclassified'**
  String get reviewUnknown;

  /// No description provided for @reviewUnverified.
  ///
  /// In en, this message translates to:
  /// **'UNVERIFIED'**
  String get reviewUnverified;

  /// No description provided for @reviewUnverifiedCount.
  ///
  /// In en, this message translates to:
  /// **'UNVERIFIED genomes found: {count}'**
  String reviewUnverifiedCount(String count);

  /// No description provided for @reviewUnverifiedAdvice.
  ///
  /// In en, this message translates to:
  /// **'We recommend removing UNVERIFIED genomes before importing.'**
  String get reviewUnverifiedAdvice;

  /// No description provided for @reviewCancelled.
  ///
  /// In en, this message translates to:
  /// **'Genome selection cancelled.'**
  String get reviewCancelled;

  /// No description provided for @reviewSearch.
  ///
  /// In en, this message translates to:
  /// **'Filter by accession, species or description'**
  String get reviewSearch;

  /// No description provided for @reviewShowing.
  ///
  /// In en, this message translates to:
  /// **'Showing {shown}/{total} genomes matching the search.'**
  String reviewShowing(String shown, String total);

  /// No description provided for @gbRecordsHeading.
  ///
  /// In en, this message translates to:
  /// **'GenBank records'**
  String get gbRecordsHeading;

  /// No description provided for @fastaHeadersHeading.
  ///
  /// In en, this message translates to:
  /// **'FASTA headers'**
  String get fastaHeadersHeading;

  /// No description provided for @clearAll.
  ///
  /// In en, this message translates to:
  /// **'Clear all'**
  String get clearAll;

  /// No description provided for @clearTitle.
  ///
  /// In en, this message translates to:
  /// **'Clear all records?'**
  String get clearTitle;

  /// No description provided for @clearBody.
  ///
  /// In en, this message translates to:
  /// **'This action cannot be undone. All loaded records and selections will be removed.'**
  String get clearBody;

  /// No description provided for @recordsReview.
  ///
  /// In en, this message translates to:
  /// **'Review and correct the species names below using the edit button.'**
  String get recordsReview;

  /// No description provided for @recordsReviewMethod.
  ///
  /// In en, this message translates to:
  /// **'Review and correct the species names below using the edit button, then fetch the taxonomy with the method of [[dataFishing]].'**
  String get recordsReviewMethod;

  /// No description provided for @recordsFetchPrefix.
  ///
  /// In en, this message translates to:
  /// **'Then click'**
  String get recordsFetchPrefix;

  /// No description provided for @recordsFetchSuffix.
  ///
  /// In en, this message translates to:
  /// **'to retrieve classification ranks from GBIF via'**
  String get recordsFetchSuffix;

  /// No description provided for @internetRequired.
  ///
  /// In en, this message translates to:
  /// **'Internet required'**
  String get internetRequired;

  /// No description provided for @fetchTaxonomy.
  ///
  /// In en, this message translates to:
  /// **'Fetch Taxonomy'**
  String get fetchTaxonomy;

  /// No description provided for @skipTaxonomy.
  ///
  /// In en, this message translates to:
  /// **'Skip taxonomy and proceed'**
  String get skipTaxonomy;

  /// No description provided for @skipTaxBody1.
  ///
  /// In en, this message translates to:
  /// **'If species names contain special characters or there are duplicate entries, errors may occur in future analyses.'**
  String get skipTaxBody1;

  /// No description provided for @skipTaxBody2.
  ///
  /// In en, this message translates to:
  /// **'It is recommended to correct species names and fetch taxonomy for best results.'**
  String get skipTaxBody2;

  /// No description provided for @taxonomyResolved.
  ///
  /// In en, this message translates to:
  /// **'Taxonomy ready'**
  String get taxonomyResolved;

  /// No description provided for @taxonomyTree.
  ///
  /// In en, this message translates to:
  /// **'Taxonomic tree'**
  String get taxonomyTree;

  /// No description provided for @taxonomyTreeBody.
  ///
  /// In en, this message translates to:
  /// **'Indented taxonomy summary using the same records shown below.'**
  String get taxonomyTreeBody;

  /// No description provided for @columnsLabel.
  ///
  /// In en, this message translates to:
  /// **'Columns:'**
  String get columnsLabel;

  /// No description provided for @colAccession.
  ///
  /// In en, this message translates to:
  /// **'Accession'**
  String get colAccession;

  /// No description provided for @colKingdom.
  ///
  /// In en, this message translates to:
  /// **'Kingdom'**
  String get colKingdom;

  /// No description provided for @colPhylum.
  ///
  /// In en, this message translates to:
  /// **'Phylum'**
  String get colPhylum;

  /// No description provided for @colClass.
  ///
  /// In en, this message translates to:
  /// **'Class'**
  String get colClass;

  /// No description provided for @colOrder.
  ///
  /// In en, this message translates to:
  /// **'Order'**
  String get colOrder;

  /// No description provided for @colFamily.
  ///
  /// In en, this message translates to:
  /// **'Family'**
  String get colFamily;

  /// No description provided for @colGenus.
  ///
  /// In en, this message translates to:
  /// **'Genus'**
  String get colGenus;

  /// No description provided for @colSpecies.
  ///
  /// In en, this message translates to:
  /// **'Species'**
  String get colSpecies;

  /// No description provided for @colAuthorship.
  ///
  /// In en, this message translates to:
  /// **'Authorship'**
  String get colAuthorship;

  /// No description provided for @colSource.
  ///
  /// In en, this message translates to:
  /// **'Source'**
  String get colSource;

  /// No description provided for @sourceFile.
  ///
  /// In en, this message translates to:
  /// **'File'**
  String get sourceFile;

  /// No description provided for @recordsSearch.
  ///
  /// In en, this message translates to:
  /// **'Search accession, species, family...'**
  String get recordsSearch;

  /// No description provided for @editRecordTooltip.
  ///
  /// In en, this message translates to:
  /// **'Edit {species} ({accession})'**
  String editRecordTooltip(String species, String accession);

  /// No description provided for @removeRecordTooltip.
  ///
  /// In en, this message translates to:
  /// **'Remove {species} ({accession})'**
  String removeRecordTooltip(String species, String accession);

  /// No description provided for @removeTitle.
  ///
  /// In en, this message translates to:
  /// **'Remove record?'**
  String get removeTitle;

  /// No description provided for @removeBody.
  ///
  /// In en, this message translates to:
  /// **'This action cannot be undone.'**
  String get removeBody;

  /// No description provided for @remove.
  ///
  /// In en, this message translates to:
  /// **'Remove'**
  String get remove;

  /// No description provided for @editTitle.
  ///
  /// In en, this message translates to:
  /// **'Edit Record'**
  String get editTitle;

  /// No description provided for @editHint.
  ///
  /// In en, this message translates to:
  /// **'Review the current record and update species and taxonomy fields below.'**
  String get editHint;

  /// No description provided for @editGbif.
  ///
  /// In en, this message translates to:
  /// **'Fetch Taxonomy from GBIF'**
  String get editGbif;

  /// No description provided for @save.
  ///
  /// In en, this message translates to:
  /// **'Save'**
  String get save;

  /// No description provided for @autoLabel.
  ///
  /// In en, this message translates to:
  /// **'(auto)'**
  String get autoLabel;

  /// No description provided for @fileSize.
  ///
  /// In en, this message translates to:
  /// **'Size'**
  String get fileSize;

  /// No description provided for @sequenceLabel.
  ///
  /// In en, this message translates to:
  /// **'Sequence'**
  String get sequenceLabel;

  /// No description provided for @featuresLabel.
  ///
  /// In en, this message translates to:
  /// **'Features'**
  String get featuresLabel;

  /// No description provided for @lineageLabel.
  ///
  /// In en, this message translates to:
  /// **'Lineage'**
  String get lineageLabel;

  /// No description provided for @speciesRequired.
  ///
  /// In en, this message translates to:
  /// **'Species name is required'**
  String get speciesRequired;

  /// No description provided for @speciesLetters.
  ///
  /// In en, this message translates to:
  /// **'Only letters, spaces and hyphens are allowed'**
  String get speciesLetters;

  /// No description provided for @speciesBinomial.
  ///
  /// In en, this message translates to:
  /// **'Use genus and epithet (e.g. Rhinella marina)'**
  String get speciesBinomial;

  /// No description provided for @speciesGenusUpper.
  ///
  /// In en, this message translates to:
  /// **'The genus must start with an uppercase letter'**
  String get speciesGenusUpper;

  /// No description provided for @familySuffix.
  ///
  /// In en, this message translates to:
  /// **'Family names end in -idae or -aceae'**
  String get familySuffix;

  /// No description provided for @gbifNoMatch.
  ///
  /// In en, this message translates to:
  /// **'No GBIF match for this name.'**
  String get gbifNoMatch;

  /// No description provided for @gbifMatched.
  ///
  /// In en, this message translates to:
  /// **'GBIF match found.'**
  String get gbifMatched;

  /// No description provided for @step3Fasta.
  ///
  /// In en, this message translates to:
  /// **'Sample and Locus Fields'**
  String get step3Fasta;

  /// No description provided for @featureTypesIntro.
  ///
  /// In en, this message translates to:
  /// **'Select the gene categories you want to include in your analysis. Toggle each category to activate or deactivate it.'**
  String get featureTypesIntro;

  /// No description provided for @standardizedPrefix.
  ///
  /// In en, this message translates to:
  /// **'Gene names have been standardized using'**
  String get standardizedPrefix;

  /// No description provided for @featureCdsInfo.
  ///
  /// In en, this message translates to:
  /// **'Protein-coding sequences (genes translated into proteins)'**
  String get featureCdsInfo;

  /// No description provided for @featureRrnaInfo.
  ///
  /// In en, this message translates to:
  /// **'Ribosomal RNA genes (12S, 16S, etc.)'**
  String get featureRrnaInfo;

  /// No description provided for @featureTrnaInfo.
  ///
  /// In en, this message translates to:
  /// **'Transfer RNA genes'**
  String get featureTrnaInfo;

  /// No description provided for @featureTypesHint.
  ///
  /// In en, this message translates to:
  /// **'Turn on at least one category to list the markers below.'**
  String get featureTypesHint;

  /// No description provided for @unknownGenesButton.
  ///
  /// In en, this message translates to:
  /// **'{n} unknown names'**
  String unknownGenesButton(String n);

  /// No description provided for @unknownTitle.
  ///
  /// In en, this message translates to:
  /// **'Gene Names Not Found in SynGenes (Rabelo et al. 2024)'**
  String get unknownTitle;

  /// No description provided for @unknownSubtitle.
  ///
  /// In en, this message translates to:
  /// **'Select the correct standardized name for each gene found in your files.'**
  String get unknownSubtitle;

  /// No description provided for @unknownCount.
  ///
  /// In en, this message translates to:
  /// **'{count} gene(s) could not be automatically standardized.'**
  String unknownCount(String count);

  /// No description provided for @unknownHowto.
  ///
  /// In en, this message translates to:
  /// **'Select the standardized name from the dropdown in the Correct Name column. Choose (skip) to exclude the gene from analysis.'**
  String get unknownHowto;

  /// No description provided for @unknownSkip.
  ///
  /// In en, this message translates to:
  /// **'(skip)'**
  String get unknownSkip;

  /// No description provided for @geneSuggested.
  ///
  /// In en, this message translates to:
  /// **'Suggested: {name}'**
  String geneSuggested(String name);

  /// No description provided for @geneFromProduct.
  ///
  /// In en, this message translates to:
  /// **'from /product “{value}”'**
  String geneFromProduct(String value);

  /// No description provided for @geneFromGene.
  ///
  /// In en, this message translates to:
  /// **'from /gene “{value}”'**
  String geneFromGene(String value);

  /// No description provided for @geneFromSequence.
  ///
  /// In en, this message translates to:
  /// **'sequence {percent}% similar to {name} in {n} other genomes'**
  String geneFromSequence(String percent, String name, String n);

  /// No description provided for @geneSequenceAgrees.
  ///
  /// In en, this message translates to:
  /// **'the sequence agrees ({percent}%)'**
  String geneSequenceAgrees(String percent);

  /// No description provided for @geneSequenceDisagrees.
  ///
  /// In en, this message translates to:
  /// **'but the sequence looks like {name} ({percent}%)'**
  String geneSequenceDisagrees(String name, String percent);

  /// No description provided for @geneNoSuggestion.
  ///
  /// In en, this message translates to:
  /// **'No match in SynGenes or in the other genomes; check the annotation.'**
  String get geneNoSuggestion;

  /// No description provided for @geneQualifierMismatch.
  ///
  /// In en, this message translates to:
  /// **'/gene “{gene}” and /product “{product}” differ'**
  String geneQualifierMismatch(String gene, String product);

  /// No description provided for @conflictsBanner.
  ///
  /// In en, this message translates to:
  /// **'{n} feature(s) have /gene and /product naming different genes.'**
  String conflictsBanner(String n);

  /// No description provided for @conflictsReview.
  ///
  /// In en, this message translates to:
  /// **'Review'**
  String get conflictsReview;

  /// No description provided for @conflictsTitle.
  ///
  /// In en, this message translates to:
  /// **'/gene and /product disagree'**
  String get conflictsTitle;

  /// No description provided for @conflictsSubtitle.
  ///
  /// In en, this message translates to:
  /// **'SPLACE uses the /product name. Each sequence is compared with the same marker in the other genomes to show which qualifier is right.'**
  String get conflictsSubtitle;

  /// No description provided for @conflictUsed.
  ///
  /// In en, this message translates to:
  /// **'Used: {name}'**
  String conflictUsed(String name);

  /// No description provided for @conflictSeqGene.
  ///
  /// In en, this message translates to:
  /// **'Sequence matches {name} ({percent}%): /gene looks right'**
  String conflictSeqGene(String name, String percent);

  /// No description provided for @conflictSeqUsed.
  ///
  /// In en, this message translates to:
  /// **'Sequence matches {name} ({percent}%): the name used looks right'**
  String conflictSeqUsed(String name, String percent);

  /// No description provided for @conflictSeqOther.
  ///
  /// In en, this message translates to:
  /// **'Sequence matches {name} ({percent}%), neither qualifier'**
  String conflictSeqOther(String name, String percent);

  /// No description provided for @conflictSeqNone.
  ///
  /// In en, this message translates to:
  /// **'Sequence inconclusive'**
  String get conflictSeqNone;

  /// No description provided for @conflictUseGene.
  ///
  /// In en, this message translates to:
  /// **'Use {name}'**
  String conflictUseGene(String name);

  /// No description provided for @speciesIdRunning.
  ///
  /// In en, this message translates to:
  /// **'Identifying species'**
  String get speciesIdRunning;

  /// No description provided for @speciesIdGenesUsed.
  ///
  /// In en, this message translates to:
  /// **'GenBank BLAST and BOLD · {genes}'**
  String speciesIdGenesUsed(String genes);

  /// No description provided for @speciesIdFooter.
  ///
  /// In en, this message translates to:
  /// **'BLAST at NCBI usually takes 1–3 minutes; BOLD answers in seconds.'**
  String get speciesIdFooter;

  /// No description provided for @speciesIdButton.
  ///
  /// In en, this message translates to:
  /// **'Identify species (GenBank + BOLD)'**
  String get speciesIdButton;

  /// No description provided for @speciesIdShowResults.
  ///
  /// In en, this message translates to:
  /// **'Identification results'**
  String get speciesIdShowResults;

  /// No description provided for @speciesIdSetupTitle.
  ///
  /// In en, this message translates to:
  /// **'Identify species with GenBank and BOLD'**
  String get speciesIdSetupTitle;

  /// No description provided for @speciesIdSetupSubtitle.
  ///
  /// In en, this message translates to:
  /// **'Compare the genes of each genome with GenBank (BLAST) and BOLD to check the species names before fetching the taxonomy.'**
  String get speciesIdSetupSubtitle;

  /// No description provided for @speciesIdGenesLabel.
  ///
  /// In en, this message translates to:
  /// **'Genes to compare'**
  String get speciesIdGenesLabel;

  /// No description provided for @speciesIdGenesHint.
  ///
  /// In en, this message translates to:
  /// **'Choose one or more genes. Every gene is searched in GenBank with BLAST; BOLD accepts COI only.'**
  String get speciesIdGenesHint;

  /// No description provided for @speciesIdInGenomes.
  ///
  /// In en, this message translates to:
  /// **'in {n} genomes'**
  String speciesIdInGenomes(String n);

  /// No description provided for @speciesIdSources.
  ///
  /// In en, this message translates to:
  /// **'Databases'**
  String get speciesIdSources;

  /// No description provided for @speciesIdGenbankBody.
  ///
  /// In en, this message translates to:
  /// **'BLAST (megablast) against core_nt, for every selected gene.'**
  String get speciesIdGenbankBody;

  /// No description provided for @speciesIdBoldBody.
  ///
  /// In en, this message translates to:
  /// **'BOLD ID engine, public species library, COI only.'**
  String get speciesIdBoldBody;

  /// No description provided for @speciesIdBoldOff.
  ///
  /// In en, this message translates to:
  /// **'Select COI to also search BOLD.'**
  String get speciesIdBoldOff;

  /// No description provided for @speciesIdRun.
  ///
  /// In en, this message translates to:
  /// **'Identify'**
  String get speciesIdRun;

  /// No description provided for @speciesIdResultsTitle.
  ///
  /// In en, this message translates to:
  /// **'Species identification'**
  String get speciesIdResultsTitle;

  /// No description provided for @speciesIdResultsSubtitle.
  ///
  /// In en, this message translates to:
  /// **'GenBank BLAST and BOLD · {genes}. Choose the name to keep for each genome.'**
  String speciesIdResultsSubtitle(String genes);

  /// No description provided for @speciesIdOnlyDisagree.
  ///
  /// In en, this message translates to:
  /// **'Only possible misidentifications'**
  String get speciesIdOnlyDisagree;

  /// No description provided for @speciesIdSummary.
  ///
  /// In en, this message translates to:
  /// **'{disagree} of {total} genomes may be misidentified'**
  String speciesIdSummary(String disagree, String total);

  /// No description provided for @speciesIdColGenome.
  ///
  /// In en, this message translates to:
  /// **'Genome'**
  String get speciesIdColGenome;

  /// No description provided for @speciesIdColCandidates.
  ///
  /// In en, this message translates to:
  /// **'Species found (choose one)'**
  String get speciesIdColCandidates;

  /// No description provided for @speciesIdMismatch.
  ///
  /// In en, this message translates to:
  /// **'Possible misidentification'**
  String get speciesIdMismatch;

  /// No description provided for @speciesIdMatch.
  ///
  /// In en, this message translates to:
  /// **'Confirmed'**
  String get speciesIdMatch;

  /// No description provided for @speciesIdNoHits.
  ///
  /// In en, this message translates to:
  /// **'No match in GenBank or BOLD'**
  String get speciesIdNoHits;

  /// No description provided for @speciesIdKeepCurrent.
  ///
  /// In en, this message translates to:
  /// **'Keep {name}'**
  String speciesIdKeepCurrent(String name);

  /// No description provided for @speciesIdHits.
  ///
  /// In en, this message translates to:
  /// **'{n} hits'**
  String speciesIdHits(String n);

  /// No description provided for @speciesIdApply.
  ///
  /// In en, this message translates to:
  /// **'Apply {n} change(s)'**
  String speciesIdApply(String n);

  /// No description provided for @speciesIdNote.
  ///
  /// In en, this message translates to:
  /// **'Hits on the genome itself (or on the GenBank record its RefSeq NC_ copy came from) appear as “this record” and do not count as evidence.'**
  String get speciesIdNote;

  /// No description provided for @geneReviewButton.
  ///
  /// In en, this message translates to:
  /// **'Review gene names'**
  String get geneReviewButton;

  /// No description provided for @geneReviewButtonN.
  ///
  /// In en, this message translates to:
  /// **'Review gene names ({n})'**
  String geneReviewButtonN(String n);

  /// No description provided for @geneReviewDone.
  ///
  /// In en, this message translates to:
  /// **'Gene names reviewed'**
  String get geneReviewDone;

  /// No description provided for @geneReviewTitle.
  ///
  /// In en, this message translates to:
  /// **'Review gene names'**
  String get geneReviewTitle;

  /// No description provided for @geneReviewSubtitle.
  ///
  /// In en, this message translates to:
  /// **'Names SynGenes did not recognize and genes whose /gene and /product disagree. Review them before the taxonomy.'**
  String get geneReviewSubtitle;

  /// No description provided for @geneReviewUnknownSection.
  ///
  /// In en, this message translates to:
  /// **'Names not found in SynGenes ({n})'**
  String geneReviewUnknownSection(String n);

  /// No description provided for @geneReviewConflictSection.
  ///
  /// In en, this message translates to:
  /// **'/gene and /product name different genes ({n})'**
  String geneReviewConflictSection(String n);

  /// No description provided for @geneReviewConflictHelp.
  ///
  /// In en, this message translates to:
  /// **'Both qualifiers were normalized with SynGenes. Choose which name this feature gets; only this genome changes.'**
  String get geneReviewConflictHelp;

  /// No description provided for @geneReviewNothing.
  ///
  /// In en, this message translates to:
  /// **'All gene names were standardized and every /gene agrees with its /product.'**
  String get geneReviewNothing;

  /// No description provided for @conflictPickGene.
  ///
  /// In en, this message translates to:
  /// **'/gene → {name}'**
  String conflictPickGene(String name);

  /// No description provided for @conflictPickProduct.
  ///
  /// In en, this message translates to:
  /// **'/product → {name}'**
  String conflictPickProduct(String name);

  /// No description provided for @grHowTitle.
  ///
  /// In en, this message translates to:
  /// **'How to read this review'**
  String get grHowTitle;

  /// No description provided for @grHowSynGenes.
  ///
  /// In en, this message translates to:
  /// **'SynGenes recognizes a tag when it is a standard gene name or one of its known synonyms.'**
  String get grHowSynGenes;

  /// No description provided for @grHowSequence.
  ///
  /// In en, this message translates to:
  /// **'Similarity is the share of the sequence (protein for CDS, DNA for rRNA) also found in the same gene of the other genomes in this dataset: above 50% is strong, 35–50% moderate, below 35% weak.'**
  String get grHowSequence;

  /// No description provided for @grSequence.
  ///
  /// In en, this message translates to:
  /// **'Sequence'**
  String get grSequence;

  /// No description provided for @grRecognized.
  ///
  /// In en, this message translates to:
  /// **'SynGenes: {name}'**
  String grRecognized(String name);

  /// No description provided for @grNotRecognized.
  ///
  /// In en, this message translates to:
  /// **'not in SynGenes'**
  String get grNotRecognized;

  /// No description provided for @grAbsent.
  ///
  /// In en, this message translates to:
  /// **'absent'**
  String get grAbsent;

  /// No description provided for @grInGenomes.
  ///
  /// In en, this message translates to:
  /// **'in {n} genome(s)'**
  String grInGenomes(String n);

  /// No description provided for @grUseName.
  ///
  /// In en, this message translates to:
  /// **'Name to use'**
  String get grUseName;

  /// No description provided for @grNoSequenceRef.
  ///
  /// In en, this message translates to:
  /// **'no other genome to compare'**
  String get grNoSequenceRef;

  /// No description provided for @grWhyProduct.
  ///
  /// In en, this message translates to:
  /// **'The /product “{value}” is a SynGenes name of {name}.'**
  String grWhyProduct(String value, String name);

  /// No description provided for @grWhyGene.
  ///
  /// In en, this message translates to:
  /// **'The /gene “{value}” is a SynGenes name of {name}.'**
  String grWhyGene(String value, String name);

  /// No description provided for @grWhySequence.
  ///
  /// In en, this message translates to:
  /// **'Neither tag is in SynGenes, but the sequence is {percent}% similar to {name} in {n} other genomes.'**
  String grWhySequence(String percent, String name, String n);

  /// No description provided for @grWhyNone.
  ///
  /// In en, this message translates to:
  /// **'Neither tag is in SynGenes and the sequence does not clearly match any marker. Check the annotation or skip this name.'**
  String get grWhyNone;

  /// No description provided for @grSeqAgrees.
  ///
  /// In en, this message translates to:
  /// **'The sequence agrees ({percent}% similar).'**
  String grSeqAgrees(String percent);

  /// No description provided for @grSeqDisagrees.
  ///
  /// In en, this message translates to:
  /// **'Careful: the sequence looks more like {name} ({percent}%).'**
  String grSeqDisagrees(String name, String percent);

  /// No description provided for @grSays.
  ///
  /// In en, this message translates to:
  /// **'{tag} says {name}'**
  String grSays(String tag, String name);

  /// No description provided for @grSimilarityTo.
  ///
  /// In en, this message translates to:
  /// **'Similarity to {name} in the other genomes'**
  String grSimilarityTo(String name);

  /// No description provided for @grSupported.
  ///
  /// In en, this message translates to:
  /// **'Supported by the sequence'**
  String get grSupported;

  /// No description provided for @grNotSupported.
  ///
  /// In en, this message translates to:
  /// **'Not supported by the sequence'**
  String get grNotSupported;

  /// No description provided for @grConflictGeneRight.
  ///
  /// In en, this message translates to:
  /// **'The sequence is {gp}% similar to {gene} and {pp}% to {product} in the other genomes, so the /gene tag ({gene}) is probably right and the /product tag is wrong.'**
  String grConflictGeneRight(String gp, String gene, String pp, String product);

  /// No description provided for @grConflictProductRight.
  ///
  /// In en, this message translates to:
  /// **'The sequence is {pp}% similar to {product} and {gp}% to {gene} in the other genomes, so the /product tag ({product}) is probably right and the /gene tag is wrong.'**
  String grConflictProductRight(
    String pp,
    String product,
    String gp,
    String gene,
  );

  /// No description provided for @grConflictUnclear.
  ///
  /// In en, this message translates to:
  /// **'The sequence does not clearly favour either gene. Check the annotation of this genome before choosing.'**
  String get grConflictUnclear;

  /// No description provided for @grConflictIntro.
  ///
  /// In en, this message translates to:
  /// **'The /gene and /product tags of this feature name different genes. Click the name to use; only this genome changes.'**
  String get grConflictIntro;

  /// No description provided for @siHowTitle.
  ///
  /// In en, this message translates to:
  /// **'How to read the identification'**
  String get siHowTitle;

  /// No description provided for @siHowBody.
  ///
  /// In en, this message translates to:
  /// **'The genes of each genome were compared with GenBank (BLAST) and BOLD. Identity is the share of identical bases with the closest sequence found.'**
  String get siHowBody;

  /// No description provided for @siLegendSame.
  ///
  /// In en, this message translates to:
  /// **'≥ 98%: same species'**
  String get siLegendSame;

  /// No description provided for @siLegendClose.
  ///
  /// In en, this message translates to:
  /// **'95–98%: close relative'**
  String get siLegendClose;

  /// No description provided for @siLegendOther.
  ///
  /// In en, this message translates to:
  /// **'< 95%: another species'**
  String get siLegendOther;

  /// No description provided for @siColUse.
  ///
  /// In en, this message translates to:
  /// **'Use'**
  String get siColUse;

  /// No description provided for @siColSpecies.
  ///
  /// In en, this message translates to:
  /// **'Species'**
  String get siColSpecies;

  /// No description provided for @siColGenbank.
  ///
  /// In en, this message translates to:
  /// **'GenBank (BLAST)'**
  String get siColGenbank;

  /// No description provided for @siColBold.
  ///
  /// In en, this message translates to:
  /// **'BOLD (COI)'**
  String get siColBold;

  /// No description provided for @siColSeqs.
  ///
  /// In en, this message translates to:
  /// **'Matching sequences'**
  String get siColSeqs;

  /// No description provided for @siCurrentTag.
  ///
  /// In en, this message translates to:
  /// **'current name'**
  String get siCurrentTag;

  /// No description provided for @siNotFound.
  ///
  /// In en, this message translates to:
  /// **'not found'**
  String get siNotFound;

  /// No description provided for @siVerdictConfirmed.
  ///
  /// In en, this message translates to:
  /// **'The matches confirm the current name (best identity {identity}%).'**
  String siVerdictConfirmed(String identity);

  /// No description provided for @siVerdictMismatch.
  ///
  /// In en, this message translates to:
  /// **'The closest match is {best} ({identity}%), while the current name {current} reaches {currentIdentity}%. The genome may be misidentified.'**
  String siVerdictMismatch(
    String best,
    String identity,
    String current,
    String currentIdentity,
  );

  /// No description provided for @siVerdictMissing.
  ///
  /// In en, this message translates to:
  /// **'The closest match is {best} ({identity}%) and {current} is not among the matches. The genome may be misidentified, or the name may be a synonym.'**
  String siVerdictMissing(String best, String identity, String current);

  /// No description provided for @siVerdictNone.
  ///
  /// In en, this message translates to:
  /// **'No sequence in GenBank or BOLD matched the genes of this genome.'**
  String get siVerdictNone;

  /// No description provided for @noMarkersButton.
  ///
  /// In en, this message translates to:
  /// **'Remove genomes without markers ({n})'**
  String noMarkersButton(String n);

  /// No description provided for @noMarkersTitle.
  ///
  /// In en, this message translates to:
  /// **'Genomes without markers'**
  String get noMarkersTitle;

  /// No description provided for @noMarkersSubtitle.
  ///
  /// In en, this message translates to:
  /// **'These genomes have no protein-coding or rRNA gene that SPLACE can use, so they add nothing to the matrix. Check the list before removing them.'**
  String get noMarkersSubtitle;

  /// No description provided for @noMarkersReasonNone.
  ///
  /// In en, this message translates to:
  /// **'No annotated genes (0 PCGs, 0 rRNAs, 0 tRNAs): the GenBank record has no gene annotation.'**
  String get noMarkersReasonNone;

  /// No description provided for @noMarkersReasonTrna.
  ///
  /// In en, this message translates to:
  /// **'Only {n} tRNA(s): no PCGs or rRNAs.'**
  String noMarkersReasonTrna(String n);

  /// No description provided for @noMarkersReasonUnnamed.
  ///
  /// In en, this message translates to:
  /// **'{pcgs} PCG(s) and {rrnas} rRNA(s) annotated, but none has a /gene or /product name to identify it.'**
  String noMarkersReasonUnnamed(String pcgs, String rrnas);

  /// No description provided for @noMarkersRemove.
  ///
  /// In en, this message translates to:
  /// **'Remove {n} genome(s)'**
  String noMarkersRemove(String n);

  /// No description provided for @heatmapShow.
  ///
  /// In en, this message translates to:
  /// **'Show heatmap'**
  String get heatmapShow;

  /// No description provided for @heatmapHide.
  ///
  /// In en, this message translates to:
  /// **'Hide heatmap'**
  String get heatmapHide;

  /// No description provided for @heatmapHoverHint.
  ///
  /// In en, this message translates to:
  /// **'Point at a cell, a marker or a sample to see its details.'**
  String get heatmapHoverHint;

  /// No description provided for @headerPresets.
  ///
  /// In en, this message translates to:
  /// **'Presets'**
  String get headerPresets;

  /// No description provided for @headerDragHint.
  ///
  /// In en, this message translates to:
  /// **'Drag to reorder. The example comes from the first record.'**
  String get headerDragHint;

  /// No description provided for @headerAddHint.
  ///
  /// In en, this message translates to:
  /// **'Click a field to add it to the end of the header.'**
  String get headerAddHint;

  /// No description provided for @headerFirstRecord.
  ///
  /// In en, this message translates to:
  /// **'First record'**
  String get headerFirstRecord;

  /// No description provided for @headerFirstRecords.
  ///
  /// In en, this message translates to:
  /// **'Headers of the first records'**
  String get headerFirstRecords;

  /// No description provided for @siSourceFailed.
  ///
  /// In en, this message translates to:
  /// **'{source} did not answer, so its column has no results: {error}'**
  String siSourceFailed(String source, String error);

  /// No description provided for @whenToUse.
  ///
  /// In en, this message translates to:
  /// **'When to use'**
  String get whenToUse;

  /// No description provided for @trimCodonHandling.
  ///
  /// In en, this message translates to:
  /// **'Protein-coding genes'**
  String get trimCodonHandling;

  /// No description provided for @trimCodonHandlingInfo.
  ///
  /// In en, this message translates to:
  /// **'Codon positions are used later for the partitions and the saturation test (step 5), so coding genes should keep whole codons. trimAl does it with -backtrans (it trims the translated alignment and writes back the codons; -ignorestopcodon because it checks stops with the standard code) and ClipKIT with --codon. Excluding a codon position from the supermatrix is decided in step 6b, not here.'**
  String get trimCodonHandlingInfo;

  /// No description provided for @alignRunning.
  ///
  /// In en, this message translates to:
  /// **'Aligning with MAFFT'**
  String get alignRunning;

  /// No description provided for @trimRunning.
  ///
  /// In en, this message translates to:
  /// **'Trimming alignments'**
  String get trimRunning;

  /// No description provided for @markersProgress.
  ///
  /// In en, this message translates to:
  /// **'{done} of {total} markers · {mode}'**
  String markersProgress(String done, String total, String mode);

  /// No description provided for @alignFooter.
  ///
  /// In en, this message translates to:
  /// **'MAFFT runs on this computer, several markers at a time.'**
  String get alignFooter;

  /// No description provided for @trimFooter.
  ///
  /// In en, this message translates to:
  /// **'trimAl/ClipKIT run on this computer, several markers at a time.'**
  String get trimFooter;

  /// No description provided for @alignCurrent.
  ///
  /// In en, this message translates to:
  /// **'Running: {markers}'**
  String alignCurrent(String markers);

  /// No description provided for @colReason.
  ///
  /// In en, this message translates to:
  /// **'Why'**
  String get colReason;

  /// No description provided for @unknownSkipAll.
  ///
  /// In en, this message translates to:
  /// **'Skip All'**
  String get unknownSkipAll;

  /// No description provided for @unknownApply.
  ///
  /// In en, this message translates to:
  /// **'Apply Corrections'**
  String get unknownApply;

  /// No description provided for @more.
  ///
  /// In en, this message translates to:
  /// **'more'**
  String get more;

  /// No description provided for @defaultPhylogeny.
  ///
  /// In en, this message translates to:
  /// **'Default to Phylogeny'**
  String get defaultPhylogeny;

  /// No description provided for @tooltipAll.
  ///
  /// In en, this message translates to:
  /// **'Select every visible marker'**
  String get tooltipAll;

  /// No description provided for @tooltipNone.
  ///
  /// In en, this message translates to:
  /// **'Clear every selected marker'**
  String get tooltipNone;

  /// No description provided for @tooltipShared.
  ///
  /// In en, this message translates to:
  /// **'Select only genes present in all genomes (no missing data)'**
  String get tooltipShared;

  /// No description provided for @tooltipDefault.
  ///
  /// In en, this message translates to:
  /// **'Select the default gene set for phylogenetic analysis'**
  String get tooltipDefault;

  /// No description provided for @duplicateChoiceHelp.
  ///
  /// In en, this message translates to:
  /// **'The supermatrix uses one copy per genome: click the copy to use. The longest starts selected; “closest” marks the copy whose length is nearest the genomes with a single copy, usually the right one when the copies differ.'**
  String get duplicateChoiceHelp;

  /// No description provided for @defaultLabel.
  ///
  /// In en, this message translates to:
  /// **'default'**
  String get defaultLabel;

  /// No description provided for @dupReference.
  ///
  /// In en, this message translates to:
  /// **'{gene} in {n} genomes with a single copy: average {avg} bp (from {min} to {max} bp). Pick the copy closest to it.'**
  String dupReference(
    String gene,
    String n,
    String avg,
    String min,
    String max,
  );

  /// No description provided for @dupNoReference.
  ///
  /// In en, this message translates to:
  /// **'No other genome has a single copy of this gene to compare lengths.'**
  String get dupNoReference;

  /// No description provided for @dupCopiesColumn.
  ///
  /// In en, this message translates to:
  /// **'Copies (click the one to use)'**
  String get dupCopiesColumn;

  /// No description provided for @dupVsAverage.
  ///
  /// In en, this message translates to:
  /// **'{diff} bp vs average'**
  String dupVsAverage(String diff);

  /// No description provided for @dupSameAsAverage.
  ///
  /// In en, this message translates to:
  /// **'same as average'**
  String get dupSameAsAverage;

  /// No description provided for @dupClosest.
  ///
  /// In en, this message translates to:
  /// **'closest to average'**
  String get dupClosest;

  /// No description provided for @dupLongestTag.
  ///
  /// In en, this message translates to:
  /// **'longest'**
  String get dupLongestTag;

  /// No description provided for @noMarkers.
  ///
  /// In en, this message translates to:
  /// **'No markers found for the selected feature types.'**
  String get noMarkers;

  /// No description provided for @featureLabel.
  ///
  /// In en, this message translates to:
  /// **'Feature'**
  String get featureLabel;

  /// No description provided for @recordsLabel.
  ///
  /// In en, this message translates to:
  /// **'Records'**
  String get recordsLabel;

  /// No description provided for @heatmapSummary.
  ///
  /// In en, this message translates to:
  /// **'{markers} markers × {samples} samples. Click a marker name to select it.'**
  String heatmapSummary(String markers, String samples);

  /// No description provided for @lengthHeatmapTitle.
  ///
  /// In en, this message translates to:
  /// **'Marker Length Heatmap'**
  String get lengthHeatmapTitle;

  /// No description provided for @lengthHeatmapBody.
  ///
  /// In en, this message translates to:
  /// **'Each selected marker is scaled independently by its own minimum and maximum length.'**
  String get lengthHeatmapBody;

  /// No description provided for @panelHeader.
  ///
  /// In en, this message translates to:
  /// **'Configure FASTA Header'**
  String get panelHeader;

  /// No description provided for @panelHeaderBody.
  ///
  /// In en, this message translates to:
  /// **'Sequence identifier of the alignment files. The default (accession_family_genus_species) is used unless you change it in the Header Builder.'**
  String get panelHeaderBody;

  /// No description provided for @openHeaderBuilder.
  ///
  /// In en, this message translates to:
  /// **'Open Header Builder'**
  String get openHeaderBuilder;

  /// No description provided for @headerBuilderTitle.
  ///
  /// In en, this message translates to:
  /// **'FASTA Header Builder'**
  String get headerBuilderTitle;

  /// No description provided for @headerBuilderSubtitle.
  ///
  /// In en, this message translates to:
  /// **'Choose, order and separate the fields of every sequence name.'**
  String get headerBuilderSubtitle;

  /// No description provided for @headerFields.
  ///
  /// In en, this message translates to:
  /// **'Header fields (in order)'**
  String get headerFields;

  /// No description provided for @headerAvailable.
  ///
  /// In en, this message translates to:
  /// **'Available fields'**
  String get headerAvailable;

  /// No description provided for @headerSeparator.
  ///
  /// In en, this message translates to:
  /// **'Separator'**
  String get headerSeparator;

  /// No description provided for @headerEmpty.
  ///
  /// In en, this message translates to:
  /// **'Add at least one field.'**
  String get headerEmpty;

  /// No description provided for @headerConfirm.
  ///
  /// In en, this message translates to:
  /// **'Confirm header'**
  String get headerConfirm;

  /// No description provided for @headerConfirmed.
  ///
  /// In en, this message translates to:
  /// **'Header confirmed'**
  String get headerConfirmed;

  /// No description provided for @headerUnique.
  ///
  /// In en, this message translates to:
  /// **'Every header is unique.'**
  String get headerUnique;

  /// No description provided for @headerDuplicates.
  ///
  /// In en, this message translates to:
  /// **'{n} header(s) are shared by more than one record; numeric suffixes (_2, _3…) will be added. Add Accession to make them unique.'**
  String headerDuplicates(String n);

  /// No description provided for @moveLeft.
  ///
  /// In en, this message translates to:
  /// **'Move earlier'**
  String get moveLeft;

  /// No description provided for @moveRight.
  ///
  /// In en, this message translates to:
  /// **'Move later'**
  String get moveRight;

  /// No description provided for @panelSeqType.
  ///
  /// In en, this message translates to:
  /// **'Sequence Type'**
  String get panelSeqType;

  /// No description provided for @seqTypeBody.
  ///
  /// In en, this message translates to:
  /// **'Choose how to align sequences. Codon-aware mode aligns protein translations via MAFFT then back-translates to DNA, preserving reading frame.'**
  String get seqTypeBody;

  /// No description provided for @seqNtTitle.
  ///
  /// In en, this message translates to:
  /// **'Nucleotide'**
  String get seqNtTitle;

  /// No description provided for @seqNtBody.
  ///
  /// In en, this message translates to:
  /// **'Align DNA/RNA sequences directly'**
  String get seqNtBody;

  /// No description provided for @seqAaTitle.
  ///
  /// In en, this message translates to:
  /// **'Amino Acid'**
  String get seqAaTitle;

  /// No description provided for @seqAaBody.
  ///
  /// In en, this message translates to:
  /// **'Translate (codon_start, transl_table), then align proteins'**
  String get seqAaBody;

  /// No description provided for @seqCodonTitle.
  ///
  /// In en, this message translates to:
  /// **'Codon-aware'**
  String get seqCodonTitle;

  /// No description provided for @seqCodonBody.
  ///
  /// In en, this message translates to:
  /// **'Protein-guided codon alignment (CDS only)'**
  String get seqCodonBody;

  /// No description provided for @recNtTitle.
  ///
  /// In en, this message translates to:
  /// **'Nucleotide, direct DNA alignment'**
  String get recNtTitle;

  /// No description provided for @recNtBody.
  ///
  /// In en, this message translates to:
  /// **'Best for rRNAs, tRNAs, non-coding regions and exploratory datasets. The FASTA goes directly to MAFFT.'**
  String get recNtBody;

  /// No description provided for @recAaTitle.
  ///
  /// In en, this message translates to:
  /// **'Amino Acid, protein alignment'**
  String get recAaTitle;

  /// No description provided for @recAaBody.
  ///
  /// In en, this message translates to:
  /// **'Best when coding genes are divergent. DNA is translated first, MAFFT aligns proteins and the result is analysed as protein.'**
  String get recAaBody;

  /// No description provided for @recCodonTitle.
  ///
  /// In en, this message translates to:
  /// **'Codon-aware, CDS phylogenetics'**
  String get recCodonTitle;

  /// No description provided for @recCodonBody.
  ///
  /// In en, this message translates to:
  /// **'Recommended for CDS matrices with DNA output and codons preserved: validates CDS length, removes terminal stops, aligns proteins and back-translates.'**
  String get recCodonBody;

  /// No description provided for @geneticCodeInfo.
  ///
  /// In en, this message translates to:
  /// **'Used for translation, frame detection and the saturation test. GenBank /transl_table qualifiers take precedence.'**
  String get geneticCodeInfo;

  /// No description provided for @translationInfo.
  ///
  /// In en, this message translates to:
  /// **'Sequences are translated before alignment. MAFFT runs in --amino mode.'**
  String get translationInfo;

  /// No description provided for @codonWarning.
  ///
  /// In en, this message translates to:
  /// **'Codon-aware mode only processes CDS features. Genes with no annotated CDS (rRNA, tRNA, etc.) will be skipped. The resulting alignment preserves reading frame and supports codon-position partitioning.'**
  String get codonWarning;

  /// No description provided for @codonWorkflow.
  ///
  /// In en, this message translates to:
  /// **'Codon-aware workflow'**
  String get codonWorkflow;

  /// No description provided for @codonStep1.
  ///
  /// In en, this message translates to:
  /// **'Extract CDS'**
  String get codonStep1;

  /// No description provided for @codonStep1Body.
  ///
  /// In en, this message translates to:
  /// **'Only annotated CDS are used; rRNAs and tRNAs are skipped.'**
  String get codonStep1Body;

  /// No description provided for @codonStep2.
  ///
  /// In en, this message translates to:
  /// **'Validate frame'**
  String get codonStep2;

  /// No description provided for @codonStep2Body.
  ///
  /// In en, this message translates to:
  /// **'Length divisible by 3 after codon_start and terminal stop handling.'**
  String get codonStep2Body;

  /// No description provided for @codonStep3.
  ///
  /// In en, this message translates to:
  /// **'Translate'**
  String get codonStep3;

  /// No description provided for @codonStep3Body.
  ///
  /// In en, this message translates to:
  /// **'With the NCBI genetic code of each CDS.'**
  String get codonStep3Body;

  /// No description provided for @codonStep4.
  ///
  /// In en, this message translates to:
  /// **'MAFFT protein'**
  String get codonStep4;

  /// No description provided for @codonStep4Body.
  ///
  /// In en, this message translates to:
  /// **'Proteins aligned with mafft --amino.'**
  String get codonStep4Body;

  /// No description provided for @codonStep5.
  ///
  /// In en, this message translates to:
  /// **'Back-translate'**
  String get codonStep5;

  /// No description provided for @codonStep5Body.
  ///
  /// In en, this message translates to:
  /// **'The protein alignment guides a codon-preserving DNA alignment.'**
  String get codonStep5Body;

  /// No description provided for @allowPseudo.
  ///
  /// In en, this message translates to:
  /// **'Allow /pseudo and /pseudogene CDS'**
  String get allowPseudo;

  /// No description provided for @allowPseudoBody.
  ///
  /// In en, this message translates to:
  /// **'Disabled by default. When off, pseudo CDS are excluded before translation.'**
  String get allowPseudoBody;

  /// No description provided for @allowStops.
  ///
  /// In en, this message translates to:
  /// **'Allow internal stop codons'**
  String get allowStops;

  /// No description provided for @allowStopsBody.
  ///
  /// In en, this message translates to:
  /// **'Disabled by default. When on, sequences with internal stops stay in the run and are reported with warnings.'**
  String get allowStopsBody;

  /// No description provided for @allowBtWarnings.
  ///
  /// In en, this message translates to:
  /// **'Continue when amino acid and codon checks report warnings'**
  String get allowBtWarnings;

  /// No description provided for @allowBtWarningsBody.
  ///
  /// In en, this message translates to:
  /// **'Disabled by default. When off, any amino acid/codon mismatch after MAFFT fails the marker; when on, the warnings are logged.'**
  String get allowBtWarningsBody;

  /// No description provided for @pythonTitle.
  ///
  /// In en, this message translates to:
  /// **'ClipKIT (codon-aware trimming)'**
  String get pythonTitle;

  /// No description provided for @pythonBody.
  ///
  /// In en, this message translates to:
  /// **'ClipKIT runs with Python 3. SPLACE uses it for codon-aware alignments: it trims the protein alignment and SPLACE back-translates the result to codons.'**
  String get pythonBody;

  /// No description provided for @pythonNotChecked.
  ///
  /// In en, this message translates to:
  /// **'Python status was not checked yet.'**
  String get pythonNotChecked;

  /// No description provided for @pythonMissing.
  ///
  /// In en, this message translates to:
  /// **'Python 3 was not found in PATH.'**
  String get pythonMissing;

  /// No description provided for @pythonFound.
  ///
  /// In en, this message translates to:
  /// **'{exe} (Python {version})'**
  String pythonFound(String exe, String version);

  /// No description provided for @clipkitFound.
  ///
  /// In en, this message translates to:
  /// **'ClipKIT {version} available'**
  String clipkitFound(String version);

  /// No description provided for @clipkitMissing.
  ///
  /// In en, this message translates to:
  /// **'ClipKIT not installed'**
  String get clipkitMissing;

  /// No description provided for @checkPython.
  ///
  /// In en, this message translates to:
  /// **'Check Python'**
  String get checkPython;

  /// No description provided for @installClipkit.
  ///
  /// In en, this message translates to:
  /// **'Install ClipKIT'**
  String get installClipkit;

  /// No description provided for @mafftAuto.
  ///
  /// In en, this message translates to:
  /// **'--auto (automatic, recommended)'**
  String get mafftAuto;

  /// No description provided for @mafftLocalpair.
  ///
  /// In en, this message translates to:
  /// **'--localpair (L-INS-i, most accurate, slow)'**
  String get mafftLocalpair;

  /// No description provided for @mafftGlobalpair.
  ///
  /// In en, this message translates to:
  /// **'--globalpair (G-INS-i, global, slow)'**
  String get mafftGlobalpair;

  /// No description provided for @mafftGenafpair.
  ///
  /// In en, this message translates to:
  /// **'--genafpair (E-INS-i, long gaps, slow)'**
  String get mafftGenafpair;

  /// No description provided for @mafftRetree2.
  ///
  /// In en, this message translates to:
  /// **'--retree 2 (FFT-NS-2, fast)'**
  String get mafftRetree2;

  /// No description provided for @mafftRetree1.
  ///
  /// In en, this message translates to:
  /// **'--retree 1 (FFT-NS-1, fastest)'**
  String get mafftRetree1;

  /// No description provided for @mafftMaxIterate.
  ///
  /// In en, this message translates to:
  /// **'Max Iterations'**
  String get mafftMaxIterate;

  /// No description provided for @mafftMaxIterateHelp.
  ///
  /// In en, this message translates to:
  /// **'0 = no refinement; 1000 = max'**
  String get mafftMaxIterateHelp;

  /// No description provided for @mafftThreads.
  ///
  /// In en, this message translates to:
  /// **'Threads'**
  String get mafftThreads;

  /// No description provided for @mafftThreadsHelp.
  ///
  /// In en, this message translates to:
  /// **'(per job)'**
  String get mafftThreadsHelp;

  /// No description provided for @mafftConcurrency.
  ///
  /// In en, this message translates to:
  /// **'{cpus} logical CPUs: {jobs} job(s) × {threads} thread(s)'**
  String mafftConcurrency(String cpus, String jobs, String threads);

  /// No description provided for @mafftGapOpening.
  ///
  /// In en, this message translates to:
  /// **'Gap Opening'**
  String get mafftGapOpening;

  /// No description provided for @mafftGapExtension.
  ///
  /// In en, this message translates to:
  /// **'Gap Extension'**
  String get mafftGapExtension;

  /// No description provided for @mafftExtra.
  ///
  /// In en, this message translates to:
  /// **'Additional Parameters'**
  String get mafftExtra;

  /// No description provided for @adjustDirectionHelp.
  ///
  /// In en, this message translates to:
  /// **'Flips reverse-strand sequences'**
  String get adjustDirectionHelp;

  /// No description provided for @reorderHelp.
  ///
  /// In en, this message translates to:
  /// **'Orders the output by similarity'**
  String get reorderHelp;

  /// No description provided for @optional.
  ///
  /// In en, this message translates to:
  /// **'(optional)'**
  String get optional;

  /// No description provided for @optionalLabel.
  ///
  /// In en, this message translates to:
  /// **'Optional'**
  String get optionalLabel;

  /// No description provided for @mafftAminoNote.
  ///
  /// In en, this message translates to:
  /// **'Amino-acid and codon-aware modes add --amino automatically (--adjustdirection does not apply to proteins).'**
  String get mafftAminoNote;

  /// No description provided for @runMafftMode.
  ///
  /// In en, this message translates to:
  /// **'Run MAFFT Alignment ({mode})'**
  String runMafftMode(String mode);

  /// No description provided for @markersToAlign.
  ///
  /// In en, this message translates to:
  /// **'{n} marker(s) selected'**
  String markersToAlign(String n);

  /// No description provided for @trimPromptTitle.
  ///
  /// In en, this message translates to:
  /// **'Alignment complete - trim poorly aligned regions?'**
  String get trimPromptTitle;

  /// No description provided for @trimPromptBody.
  ///
  /// In en, this message translates to:
  /// **'trimAl removes columns with excessive gaps or low conservation, improving phylogenetic inference quality.'**
  String get trimPromptBody;

  /// No description provided for @trimPromptCodon.
  ///
  /// In en, this message translates to:
  /// **'In Codon-aware mode, ClipKIT trims the protein alignment and SPLACE back-translates it to codon-preserving DNA.'**
  String get trimPromptCodon;

  /// No description provided for @trimPromptNo.
  ///
  /// In en, this message translates to:
  /// **'No, skip trimming'**
  String get trimPromptNo;

  /// No description provided for @trimPromptYes.
  ///
  /// In en, this message translates to:
  /// **'Yes, set up trimming'**
  String get trimPromptYes;

  /// No description provided for @trimAutomated1.
  ///
  /// In en, this message translates to:
  /// **'-automated1 (heuristic, recommended)'**
  String get trimAutomated1;

  /// No description provided for @trimStrict.
  ///
  /// In en, this message translates to:
  /// **'-strict (strict automated)'**
  String get trimStrict;

  /// No description provided for @trimStrictplus.
  ///
  /// In en, this message translates to:
  /// **'-strictplus (strict+)'**
  String get trimStrictplus;

  /// No description provided for @trimGappyout.
  ///
  /// In en, this message translates to:
  /// **'-gappyout (remove gappy columns)'**
  String get trimGappyout;

  /// No description provided for @trimNogaps.
  ///
  /// In en, this message translates to:
  /// **'-nogaps (remove all gap columns)'**
  String get trimNogaps;

  /// No description provided for @trimNoallgaps.
  ///
  /// In en, this message translates to:
  /// **'-noallgaps (only all-gap columns)'**
  String get trimNoallgaps;

  /// No description provided for @trimManual.
  ///
  /// In en, this message translates to:
  /// **'Manual thresholds'**
  String get trimManual;

  /// No description provided for @trimGapThreshold.
  ///
  /// In en, this message translates to:
  /// **'Gap Threshold'**
  String get trimGapThreshold;

  /// No description provided for @trimSimilarity.
  ///
  /// In en, this message translates to:
  /// **'Similarity'**
  String get trimSimilarity;

  /// No description provided for @trimConsistency.
  ///
  /// In en, this message translates to:
  /// **'Consistency'**
  String get trimConsistency;

  /// No description provided for @trimWindow.
  ///
  /// In en, this message translates to:
  /// **'Window Size'**
  String get trimWindow;

  /// No description provided for @trimExtra.
  ///
  /// In en, this message translates to:
  /// **'Additional Parameters'**
  String get trimExtra;

  /// No description provided for @trimCodonNote.
  ///
  /// In en, this message translates to:
  /// **'For codon-aware runs MAFFT aligned proteins: the trimmer trims the protein alignment, SPLACE back-translates it to codon-preserving DNA and validates the frame.'**
  String get trimCodonNote;

  /// No description provided for @codonTrimBackend.
  ///
  /// In en, this message translates to:
  /// **'Codon trimming backend'**
  String get codonTrimBackend;

  /// No description provided for @codonTrimAuto.
  ///
  /// In en, this message translates to:
  /// **'Auto: trimAl, then ClipKIT if needed'**
  String get codonTrimAuto;

  /// No description provided for @codonTrimTrimal.
  ///
  /// In en, this message translates to:
  /// **'trimAl, protein alignment with back-translation'**
  String get codonTrimTrimal;

  /// No description provided for @codonTrimClipkit.
  ///
  /// In en, this message translates to:
  /// **'ClipKIT (Python required)'**
  String get codonTrimClipkit;

  /// No description provided for @clipkitMode.
  ///
  /// In en, this message translates to:
  /// **'ClipKIT mode'**
  String get clipkitMode;

  /// No description provided for @codonTrimmingTitle.
  ///
  /// In en, this message translates to:
  /// **'Codon trimming'**
  String get codonTrimmingTitle;

  /// No description provided for @codonTrimmingBody.
  ///
  /// In en, this message translates to:
  /// **'Removes poorly aligned codon blocks before concatenation and before any codon-position filtering.'**
  String get codonTrimmingBody;

  /// No description provided for @trimSkippedTag.
  ///
  /// In en, this message translates to:
  /// **'Trimming skipped'**
  String get trimSkippedTag;

  /// No description provided for @trimming.
  ///
  /// In en, this message translates to:
  /// **'Trimming'**
  String get trimming;

  /// No description provided for @resultsTitle.
  ///
  /// In en, this message translates to:
  /// **'Alignment Results'**
  String get resultsTitle;

  /// No description provided for @noTrimming.
  ///
  /// In en, this message translates to:
  /// **'No trimming'**
  String get noTrimming;

  /// No description provided for @viewAlignment.
  ///
  /// In en, this message translates to:
  /// **'Show alignment'**
  String get viewAlignment;

  /// No description provided for @tabProtein.
  ///
  /// In en, this message translates to:
  /// **'Protein'**
  String get tabProtein;

  /// No description provided for @viewDna.
  ///
  /// In en, this message translates to:
  /// **'DNA view'**
  String get viewDna;

  /// No description provided for @viewTranslate.
  ///
  /// In en, this message translates to:
  /// **'Translate to amino acids'**
  String get viewTranslate;

  /// No description provided for @possibleStops.
  ///
  /// In en, this message translates to:
  /// **'{n} possible internal stop(s)'**
  String possibleStops(String n);

  /// No description provided for @conservation.
  ///
  /// In en, this message translates to:
  /// **'CONSERVATION'**
  String get conservation;

  /// No description provided for @stopLegend.
  ///
  /// In en, this message translates to:
  /// **'Stop'**
  String get stopLegend;

  /// No description provided for @gapLegend.
  ///
  /// In en, this message translates to:
  /// **'Gap'**
  String get gapLegend;

  /// No description provided for @replaceGaps.
  ///
  /// In en, this message translates to:
  /// **'Replace alignment gaps with ?'**
  String get replaceGaps;

  /// No description provided for @replaceGapsBody.
  ///
  /// In en, this message translates to:
  /// **'Use when the downstream program should treat alignment gaps as missing data. Partition coordinates are preserved.'**
  String get replaceGapsBody;

  /// No description provided for @selectedLabel.
  ///
  /// In en, this message translates to:
  /// **'Selected'**
  String get selectedLabel;

  /// No description provided for @saveNexus.
  ///
  /// In en, this message translates to:
  /// **'Save NEXUS'**
  String get saveNexus;

  /// No description provided for @savePartition.
  ///
  /// In en, this message translates to:
  /// **'Save partition table'**
  String get savePartition;

  /// No description provided for @saveFasta.
  ///
  /// In en, this message translates to:
  /// **'Save FASTA'**
  String get saveFasta;

  /// No description provided for @savePhylip.
  ///
  /// In en, this message translates to:
  /// **'Save PHYLIP'**
  String get savePhylip;

  /// No description provided for @partitionName.
  ///
  /// In en, this message translates to:
  /// **'Partition'**
  String get partitionName;

  /// No description provided for @partitionRange.
  ///
  /// In en, this message translates to:
  /// **'Range'**
  String get partitionRange;

  /// No description provided for @partitionSites.
  ///
  /// In en, this message translates to:
  /// **'Sites'**
  String get partitionSites;

  /// No description provided for @partitionType.
  ///
  /// In en, this message translates to:
  /// **'Type'**
  String get partitionType;

  /// No description provided for @noneSelected.
  ///
  /// In en, this message translates to:
  /// **'none selected'**
  String get noneSelected;

  /// No description provided for @outgroupFromStep7.
  ///
  /// In en, this message translates to:
  /// **'chosen in step 6'**
  String get outgroupFromStep7;

  /// No description provided for @substitutionModel.
  ///
  /// In en, this message translates to:
  /// **'Substitution Model'**
  String get substitutionModel;

  /// No description provided for @modelTest.
  ///
  /// In en, this message translates to:
  /// **'ModelFinder: auto best-fit (TEST)'**
  String get modelTest;

  /// No description provided for @modelMfp.
  ///
  /// In en, this message translates to:
  /// **'MFP: ModelFinder Plus'**
  String get modelMfp;

  /// No description provided for @modelMfpMerge.
  ///
  /// In en, this message translates to:
  /// **'MFP+MERGE: partition model finder (recommended)'**
  String get modelMfpMerge;

  /// No description provided for @modelCustom.
  ///
  /// In en, this message translates to:
  /// **'Custom…'**
  String get modelCustom;

  /// No description provided for @bootstrapSupport.
  ///
  /// In en, this message translates to:
  /// **'Bootstrap Support'**
  String get bootstrapSupport;

  /// No description provided for @bootUf1000.
  ///
  /// In en, this message translates to:
  /// **'UFBoot2: 1000 replicates (recommended)'**
  String get bootUf1000;

  /// No description provided for @bootUf5000.
  ///
  /// In en, this message translates to:
  /// **'UFBoot2: 5000 replicates'**
  String get bootUf5000;

  /// No description provided for @bootClassic.
  ///
  /// In en, this message translates to:
  /// **'Classical bootstrap: 100 replicates (slow)'**
  String get bootClassic;

  /// No description provided for @bootNone.
  ///
  /// In en, this message translates to:
  /// **'No bootstrap'**
  String get bootNone;

  /// No description provided for @threadsLabel.
  ///
  /// In en, this message translates to:
  /// **'Threads'**
  String get threadsLabel;

  /// No description provided for @extraFlags.
  ///
  /// In en, this message translates to:
  /// **'Extra flags'**
  String get extraFlags;

  /// No description provided for @spaceSeparated.
  ///
  /// In en, this message translates to:
  /// **'(space-separated)'**
  String get spaceSeparated;

  /// No description provided for @perGeneTrees.
  ///
  /// In en, this message translates to:
  /// **'Also infer per-gene trees'**
  String get perGeneTrees;

  /// No description provided for @perGeneTreesBody.
  ///
  /// In en, this message translates to:
  /// **'Runs a separate IQ-TREE analysis for each marker (saved in gene_trees/).'**
  String get perGeneTreesBody;

  /// No description provided for @commandPreview.
  ///
  /// In en, this message translates to:
  /// **'Command preview'**
  String get commandPreview;

  /// No description provided for @saveZip.
  ///
  /// In en, this message translates to:
  /// **'Save Results as ZIP'**
  String get saveZip;

  /// No description provided for @analysisComplete.
  ///
  /// In en, this message translates to:
  /// **'Analysis complete'**
  String get analysisComplete;

  /// No description provided for @outputDirectory.
  ///
  /// In en, this message translates to:
  /// **'Output directory'**
  String get outputDirectory;

  /// No description provided for @geneTreesLabel.
  ///
  /// In en, this message translates to:
  /// **'Per-gene trees'**
  String get geneTreesLabel;

  /// No description provided for @howToCite.
  ///
  /// In en, this message translates to:
  /// **'How to Cite'**
  String get howToCite;

  /// No description provided for @howToCiteHint.
  ///
  /// In en, this message translates to:
  /// **'If you use these results, cite SPLACE and the tools it ran'**
  String get howToCiteHint;

  /// No description provided for @copy.
  ///
  /// In en, this message translates to:
  /// **'Copy'**
  String get copy;

  /// No description provided for @metaGenomes.
  ///
  /// In en, this message translates to:
  /// **'{n} genomes'**
  String metaGenomes(String n);

  /// No description provided for @metaFastaFiles.
  ///
  /// In en, this message translates to:
  /// **'{n} FASTA files'**
  String metaFastaFiles(String n);

  /// No description provided for @metaSamples.
  ///
  /// In en, this message translates to:
  /// **'{n} samples'**
  String metaSamples(String n);

  /// No description provided for @metaSelected.
  ///
  /// In en, this message translates to:
  /// **'{selected} of {total} selected'**
  String metaSelected(String selected, String total);

  /// No description provided for @metaAligned.
  ///
  /// In en, this message translates to:
  /// **'{done}/{total} aligned'**
  String metaAligned(String done, String total);

  /// No description provided for @metaTrimmed.
  ///
  /// In en, this message translates to:
  /// **'{n} trimmed'**
  String metaTrimmed(String n);

  /// No description provided for @metaMatrix.
  ///
  /// In en, this message translates to:
  /// **'{taxa} taxa × {sites} sites'**
  String metaMatrix(String taxa, String sites);

  /// No description provided for @reviewSelectedCount.
  ///
  /// In en, this message translates to:
  /// **'{selected} of {total} genomes enabled'**
  String reviewSelectedCount(String selected, String total);

  /// No description provided for @progressSearchTitle.
  ///
  /// In en, this message translates to:
  /// **'Searching NCBI by taxonomic name'**
  String get progressSearchTitle;

  /// No description provided for @progressSearchingFor.
  ///
  /// In en, this message translates to:
  /// **'Searching for {name}'**
  String progressSearchingFor(String name);

  /// No description provided for @progressFetchTitle.
  ///
  /// In en, this message translates to:
  /// **'Importing GenBank records'**
  String get progressFetchTitle;

  /// No description provided for @progressCounter.
  ///
  /// In en, this message translates to:
  /// **'{current} of {total}'**
  String progressCounter(String current, String total);

  /// No description provided for @progressSpecies.
  ///
  /// In en, this message translates to:
  /// **'{current} of {total} species'**
  String progressSpecies(String current, String total);

  /// No description provided for @progressRateLabel.
  ///
  /// In en, this message translates to:
  /// **'Request rate'**
  String get progressRateLabel;

  /// No description provided for @progressRateWithKey.
  ///
  /// In en, this message translates to:
  /// **'Up to 10 requests/s'**
  String get progressRateWithKey;

  /// No description provided for @progressRateDefault.
  ///
  /// In en, this message translates to:
  /// **'Up to 3 requests/s'**
  String get progressRateDefault;

  /// No description provided for @progressApiEnabled.
  ///
  /// In en, this message translates to:
  /// **'Saved and active'**
  String get progressApiEnabled;

  /// No description provided for @progressApiDisabled.
  ///
  /// In en, this message translates to:
  /// **'Not provided'**
  String get progressApiDisabled;

  /// No description provided for @progressCurrentLabel.
  ///
  /// In en, this message translates to:
  /// **'Current item'**
  String get progressCurrentLabel;

  /// No description provided for @progressCurrentSpecies.
  ///
  /// In en, this message translates to:
  /// **'Current species'**
  String get progressCurrentSpecies;

  /// No description provided for @progressStatusIdle.
  ///
  /// In en, this message translates to:
  /// **'Waiting to start'**
  String get progressStatusIdle;

  /// No description provided for @progressFooter.
  ///
  /// In en, this message translates to:
  /// **'SPLACE is contacting NCBI and updating the progress continuously.'**
  String get progressFooter;

  /// No description provided for @progressTaxonomyTitle.
  ///
  /// In en, this message translates to:
  /// **'Fetching taxonomy from {source}'**
  String progressTaxonomyTitle(String source);

  /// No description provided for @progressTaxonomyFooter.
  ///
  /// In en, this message translates to:
  /// **'SPLACE is querying {source} and updating the records as matches are found.'**
  String progressTaxonomyFooter(String source);

  /// No description provided for @progressWorking.
  ///
  /// In en, this message translates to:
  /// **'Working…'**
  String get progressWorking;

  /// No description provided for @colStatus.
  ///
  /// In en, this message translates to:
  /// **'Status (SALVE)'**
  String get colStatus;

  /// No description provided for @fetchTaxonomyFrom.
  ///
  /// In en, this message translates to:
  /// **'Fetch taxonomy ({source})'**
  String fetchTaxonomyFrom(String source);

  /// No description provided for @taxonomyMatched.
  ///
  /// In en, this message translates to:
  /// **'{source} match found.'**
  String taxonomyMatched(String source);

  /// No description provided for @taxonomyNoMatch.
  ///
  /// In en, this message translates to:
  /// **'No {source} match for this name.'**
  String taxonomyNoMatch(String source);

  /// No description provided for @taxonomyCardTitle.
  ///
  /// In en, this message translates to:
  /// **'Taxonomy of the records'**
  String get taxonomyCardTitle;

  /// No description provided for @taxonomySourceLabel.
  ///
  /// In en, this message translates to:
  /// **'Taxonomy source'**
  String get taxonomySourceLabel;

  /// No description provided for @sourceGbifBody.
  ///
  /// In en, this message translates to:
  /// **'Accepted names and authorship from the Global Biodiversity Information Facility backbone.'**
  String get sourceGbifBody;

  /// No description provided for @sourceNcbiBody.
  ///
  /// In en, this message translates to:
  /// **'Lineage of the taxon id in each GenBank record ({withId} of {total} records have one; the others are searched by name).'**
  String sourceNcbiBody(String withId, String total);

  /// No description provided for @sourceSalveBody.
  ///
  /// In en, this message translates to:
  /// **'Brazilian fauna only. ICMBio\'s national assessment: taxonomy, extinction-risk category and the species sheet, from the SALVE public API.'**
  String get sourceSalveBody;

  /// No description provided for @badgeOnline.
  ///
  /// In en, this message translates to:
  /// **'online'**
  String get badgeOnline;

  /// No description provided for @salveHint.
  ///
  /// In en, this message translates to:
  /// **'**Only species of the Brazilian fauna.** SALVE holds ICMBio\'s national extinction-risk assessments; species from elsewhere have no sheet and are listed as not found (use GBIF or NCBI for them). Every section of each sheet is read online: ranks, category in the IUCN colours and criteria, common names, distribution (endemism, biomes, states, watersheds), natural history (habitats, diet, reproduction, migration), population (trend, generation time), threats, uses, official lists, assessment history, conservation actions, protected areas, research, assessors, citation and DOI.'**
  String get salveHint;

  /// No description provided for @poweredByDataFishing.
  ///
  /// In en, this message translates to:
  /// **'Method of'**
  String get poweredByDataFishing;

  /// No description provided for @metaSatKeep.
  ///
  /// In en, this message translates to:
  /// **'keep all positions'**
  String get metaSatKeep;

  /// No description provided for @metaSatExclude.
  ///
  /// In en, this message translates to:
  /// **'remove position {positions}'**
  String metaSatExclude(String positions);

  /// No description provided for @metaSatAmino.
  ///
  /// In en, this message translates to:
  /// **'use amino acids'**
  String get metaSatAmino;

  /// No description provided for @apiKeyActive.
  ///
  /// In en, this message translates to:
  /// **'Active'**
  String get apiKeyActive;

  /// No description provided for @apiKeyInactive.
  ///
  /// In en, this message translates to:
  /// **'Not set'**
  String get apiKeyInactive;

  /// No description provided for @apiKeyWithout.
  ///
  /// In en, this message translates to:
  /// **'Without a key'**
  String get apiKeyWithout;

  /// No description provided for @apiKeyWith.
  ///
  /// In en, this message translates to:
  /// **'With a key'**
  String get apiKeyWith;

  /// No description provided for @apiKeyUsedFor.
  ///
  /// In en, this message translates to:
  /// **'Used by the NCBI search, accession downloads and NCBI Taxonomy.'**
  String get apiKeyUsedFor;

  /// No description provided for @apiKeyHowTitle.
  ///
  /// In en, this message translates to:
  /// **'How to get a key'**
  String get apiKeyHowTitle;

  /// No description provided for @apiKeyHow1.
  ///
  /// In en, this message translates to:
  /// **'Sign in to your NCBI account (free).'**
  String get apiKeyHow1;

  /// No description provided for @apiKeyHow2.
  ///
  /// In en, this message translates to:
  /// **'Open Account settings → API Key Management and click “Create an API Key”.'**
  String get apiKeyHow2;

  /// No description provided for @apiKeyHow3.
  ///
  /// In en, this message translates to:
  /// **'Copy the key and paste it in the field below.'**
  String get apiKeyHow3;

  /// No description provided for @apiKeyOpenNcbi.
  ///
  /// In en, this message translates to:
  /// **'Open NCBI account'**
  String get apiKeyOpenNcbi;

  /// No description provided for @apiKeyPaste.
  ///
  /// In en, this message translates to:
  /// **'Paste'**
  String get apiKeyPaste;

  /// No description provided for @apiKeyFieldLabel.
  ///
  /// In en, this message translates to:
  /// **'Your API key'**
  String get apiKeyFieldLabel;

  /// No description provided for @apiKeyLooksValid.
  ///
  /// In en, this message translates to:
  /// **'Format looks valid'**
  String get apiKeyLooksValid;

  /// No description provided for @apiKeyStored.
  ///
  /// In en, this message translates to:
  /// **'Stored only on this computer (Documents/SPLACE/settings.json).'**
  String get apiKeyStored;

  /// No description provided for @apiKeyEstimateTitle.
  ///
  /// In en, this message translates to:
  /// **'Estimated download time'**
  String get apiKeyEstimateTitle;

  /// No description provided for @apiKeyEstimateGenomes.
  ///
  /// In en, this message translates to:
  /// **'{n} genomes'**
  String apiKeyEstimateGenomes(String n);

  /// No description provided for @apiKeyEstimateDetail.
  ///
  /// In en, this message translates to:
  /// **'{batch} genomes per request · up to {rate} requests/s'**
  String apiKeyEstimateDetail(String batch, String rate);

  /// No description provided for @apiKeyEstimateGain.
  ///
  /// In en, this message translates to:
  /// **'{factor}× faster with a key · {saved} less'**
  String apiKeyEstimateGain(String factor, String saved);

  /// No description provided for @apiKeyEstimateNote.
  ///
  /// In en, this message translates to:
  /// **'Download of mitochondrial genomes not yet in the local cache, parsing included, from times measured with this app on Nov 8 2026 (Testudines, about 75 kB per record): 50, 200 and 700 genomes took 1.6, 4.4 and 10.1 s with a key and 5.5, 16.9 and 64.6 s without one. The search and the list of results add a few seconds; records already in the cache load in a moment. The actual time depends on your connection and on NCBI\'s load.'**
  String get apiKeyEstimateNote;

  /// No description provided for @apiKeyCurrent.
  ///
  /// In en, this message translates to:
  /// **'Current'**
  String get apiKeyCurrent;

  /// No description provided for @reviewFound.
  ///
  /// In en, this message translates to:
  /// **'{n} genomes found · {query}'**
  String reviewFound(String n, String query);

  /// No description provided for @reviewColOrganism.
  ///
  /// In en, this message translates to:
  /// **'Organism'**
  String get reviewColOrganism;

  /// No description provided for @reviewColStatus.
  ///
  /// In en, this message translates to:
  /// **'Status'**
  String get reviewColStatus;

  /// No description provided for @reviewUnverifiedNote.
  ///
  /// In en, this message translates to:
  /// **'{n} genomes are flagged UNVERIFIED by NCBI. We recommend leaving them out.'**
  String reviewUnverifiedNote(String n);

  /// No description provided for @reviewDeselectUnverified.
  ///
  /// In en, this message translates to:
  /// **'Deselect UNVERIFIED'**
  String get reviewDeselectUnverified;

  /// No description provided for @reviewSelectMenu.
  ///
  /// In en, this message translates to:
  /// **'Select'**
  String get reviewSelectMenu;

  /// No description provided for @reviewImportN.
  ///
  /// In en, this message translates to:
  /// **'Import {n}'**
  String reviewImportN(String n);

  /// No description provided for @reviewNoMatch.
  ///
  /// In en, this message translates to:
  /// **'No genome matches the filter.'**
  String get reviewNoMatch;

  /// No description provided for @reviewSelectedOf.
  ///
  /// In en, this message translates to:
  /// **'{selected} of {total} selected'**
  String reviewSelectedOf(String selected, String total);

  /// No description provided for @metaTrees.
  ///
  /// In en, this message translates to:
  /// **'{n} tree(s)'**
  String metaTrees(String n);

  /// No description provided for @itemQueued.
  ///
  /// In en, this message translates to:
  /// **'Waiting'**
  String get itemQueued;

  /// No description provided for @itemCached.
  ///
  /// In en, this message translates to:
  /// **'From the cache · {note}'**
  String itemCached(String note);

  /// No description provided for @itemSkipped.
  ///
  /// In en, this message translates to:
  /// **'Skipped'**
  String get itemSkipped;

  /// No description provided for @itemsRunning.
  ///
  /// In en, this message translates to:
  /// **'{n} running'**
  String itemsRunning(String n);

  /// No description provided for @itemsQueued.
  ///
  /// In en, this message translates to:
  /// **'{n} waiting'**
  String itemsQueued(String n);

  /// No description provided for @itemsDone.
  ///
  /// In en, this message translates to:
  /// **'{n} done'**
  String itemsDone(String n);

  /// No description provided for @itemsCached.
  ///
  /// In en, this message translates to:
  /// **'{n} from the cache (unchanged)'**
  String itemsCached(String n);

  /// No description provided for @itemsSkipped.
  ///
  /// In en, this message translates to:
  /// **'{n} skipped'**
  String itemsSkipped(String n);

  /// No description provided for @itemsFailed.
  ///
  /// In en, this message translates to:
  /// **'{n} failed'**
  String itemsFailed(String n);

  /// No description provided for @trimToolNt.
  ///
  /// In en, this message translates to:
  /// **'Choose one tool: trimAl or ClipKIT keeping whole codons in coding genes (recommended), or trimAl column by column (splits codons).'**
  String get trimToolNt;

  /// No description provided for @trimToolAa.
  ///
  /// In en, this message translates to:
  /// **'trimAl trims the protein alignment.'**
  String get trimToolAa;

  /// No description provided for @trimToolCodon.
  ///
  /// In en, this message translates to:
  /// **'The codon alignment is trimmed keeping whole codons: by trimAl (-backtrans) or ClipKIT (--codon).'**
  String get trimToolCodon;

  /// No description provided for @runTrimWith.
  ///
  /// In en, this message translates to:
  /// **'Trim with {tool}'**
  String runTrimWith(String tool);

  /// No description provided for @trimDescAutomated1.
  ///
  /// In en, this message translates to:
  /// **'trimAl looks at the alignment and picks the best of -gappyout and -strict by itself. Recommended for maximum-likelihood trees.'**
  String get trimDescAutomated1;

  /// No description provided for @trimDescStrict.
  ///
  /// In en, this message translates to:
  /// **'Removes gap-rich columns and also columns with low similarity between sequences.'**
  String get trimDescStrict;

  /// No description provided for @trimDescStrictplus.
  ///
  /// In en, this message translates to:
  /// **'Like -strict, but tuned for neighbour-joining trees; removes more.'**
  String get trimDescStrictplus;

  /// No description provided for @trimDescGappyout.
  ///
  /// In en, this message translates to:
  /// **'Removes only the gap-rich columns, with a threshold computed from the gap distribution of the alignment.'**
  String get trimDescGappyout;

  /// No description provided for @trimDescNogaps.
  ///
  /// In en, this message translates to:
  /// **'Removes every column that has at least one gap. Can remove a lot.'**
  String get trimDescNogaps;

  /// No description provided for @trimDescNoallgaps.
  ///
  /// In en, this message translates to:
  /// **'Removes only columns made entirely of gaps. The lightest option.'**
  String get trimDescNoallgaps;

  /// No description provided for @trimDescManual.
  ///
  /// In en, this message translates to:
  /// **'You set the gap, similarity and consistency thresholds yourself.'**
  String get trimDescManual;

  /// No description provided for @clipkitDescGappy.
  ///
  /// In en, this message translates to:
  /// **'Removes columns with more gaps than the threshold (90 %). Simple and conservative.'**
  String get clipkitDescGappy;

  /// No description provided for @clipkitDescKpic.
  ///
  /// In en, this message translates to:
  /// **'Keeps only parsimony-informative and constant columns; removes the rest.'**
  String get clipkitDescKpic;

  /// No description provided for @clipkitDescKpicGappy.
  ///
  /// In en, this message translates to:
  /// **'kpic, and also removes gap-rich columns.'**
  String get clipkitDescKpicGappy;

  /// No description provided for @clipkitDescSmartGap.
  ///
  /// In en, this message translates to:
  /// **'Chooses the gap threshold from the alignment itself (ClipKIT default).'**
  String get clipkitDescSmartGap;

  /// No description provided for @clipkitRequired.
  ///
  /// In en, this message translates to:
  /// **'Codon-aware trimming uses ClipKIT, which is not installed yet. Install it below (needs Python 3).'**
  String get clipkitRequired;

  /// No description provided for @frameInfo.
  ///
  /// In en, this message translates to:
  /// **'Where the first complete codon starts in this alignment. Base 1 (+1): the alignment starts at the 1st position of a codon, as in most CDS. +2 or +3: it starts 1 or 2 bases into a codon. SPLACE picks the frame with the fewest stop codons; it is used to translate, to trim whole codons and to split the 1st, 2nd and 3rd codon positions (saturation test and partitions).'**
  String get frameInfo;

  /// No description provided for @frameOption.
  ///
  /// In en, this message translates to:
  /// **'Base {n} (+{n})'**
  String frameOption(String n);

  /// No description provided for @frameTag.
  ///
  /// In en, this message translates to:
  /// **'Codons start at base {n}'**
  String frameTag(String n);

  /// No description provided for @viewerHint.
  ///
  /// In en, this message translates to:
  /// **'Hover a base to see its position · drag with the right button to scroll · drag the edge of the names to resize, double-click to fit.'**
  String get viewerHint;

  /// No description provided for @viewerColumn.
  ///
  /// In en, this message translates to:
  /// **'Column {n}'**
  String viewerColumn(String n);

  /// No description provided for @viewerSeqPos.
  ///
  /// In en, this message translates to:
  /// **'position {n} in the sequence (gaps not counted)'**
  String viewerSeqPos(String n);

  /// No description provided for @viewerGap.
  ///
  /// In en, this message translates to:
  /// **'gap'**
  String get viewerGap;

  /// No description provided for @viewerIdentical.
  ///
  /// In en, this message translates to:
  /// **'{pct}% of sequences share this base'**
  String viewerIdentical(String pct);

  /// No description provided for @viewerCodonPos.
  ///
  /// In en, this message translates to:
  /// **'codon position {n}'**
  String viewerCodonPos(String n);

  /// No description provided for @viewerNames.
  ///
  /// In en, this message translates to:
  /// **'Sequences'**
  String get viewerNames;

  /// No description provided for @siCacheInfo.
  ///
  /// In en, this message translates to:
  /// **'{cached} of {total} searches were already done before and will be reused; only {todo} will be sent.'**
  String siCacheInfo(String cached, String total, String todo);

  /// No description provided for @siCacheAll.
  ///
  /// In en, this message translates to:
  /// **'All {total} searches were done before: the results appear right away, nothing is sent.'**
  String siCacheAll(String total);

  /// No description provided for @siCacheClear.
  ///
  /// In en, this message translates to:
  /// **'Clear saved results'**
  String get siCacheClear;

  /// No description provided for @trimFrameWarning.
  ///
  /// In en, this message translates to:
  /// **'**trimAl removes single columns, not codons.** In a protein-coding gene aligned as DNA, removing one or two bases of a codon shifts the reading frame of everything after it: 1st, 2nd and 3rd codon positions get mixed. Without \"Keep the reading frame\" these genes become one partition without codon positions and are left out of the saturation test.'**
  String get trimFrameWarning;

  /// No description provided for @frameLostTag.
  ///
  /// In en, this message translates to:
  /// **'Reading frame lost in trimming'**
  String get frameLostTag;

  /// No description provided for @frameLostInfo.
  ///
  /// In en, this message translates to:
  /// **'trimAl removed single columns of this coding gene, so its length may not be a multiple of 3 and the codon positions are unknown. It is concatenated as one partition and left out of the saturation test. Turn on \"Keep the reading frame\" in 4d and trim again to keep codon partitions.'**
  String get frameLostInfo;

  /// No description provided for @frameShiftTag.
  ///
  /// In en, this message translates to:
  /// **'{n} sequence(s) with frame shifts'**
  String frameShiftTag(String n);

  /// No description provided for @frameShiftInfo.
  ///
  /// In en, this message translates to:
  /// **'MAFFT placed gaps that are not multiples of 3 inside these sequences: from that gap on, their bases sit in the wrong codon position of the column. This is common when coding genes are aligned as nucleotides. Codon partitions and the saturation test are approximate for them; the Codon-aware mode (4b) avoids it.'**
  String get frameShiftInfo;

  /// No description provided for @saturationRunning.
  ///
  /// In en, this message translates to:
  /// **'Saturation test'**
  String get saturationRunning;

  /// No description provided for @saturationProgress.
  ///
  /// In en, this message translates to:
  /// **'{done} of {total} alignments'**
  String saturationProgress(String done, String total);

  /// No description provided for @matrixRunning.
  ///
  /// In en, this message translates to:
  /// **'Building the supermatrix'**
  String get matrixRunning;

  /// No description provided for @matrixRunningBody.
  ///
  /// In en, this message translates to:
  /// **'Concatenating the alignments'**
  String get matrixRunningBody;

  /// No description provided for @iqtreeRunning.
  ///
  /// In en, this message translates to:
  /// **'IQ-TREE: maximum-likelihood tree'**
  String get iqtreeRunning;

  /// No description provided for @raxmlRunning.
  ///
  /// In en, this message translates to:
  /// **'RAxML-NG: maximum-likelihood tree'**
  String get raxmlRunning;

  /// No description provided for @treeRunSubtitle.
  ///
  /// In en, this message translates to:
  /// **'Model {model} · elapsed {elapsed}'**
  String treeRunSubtitle(String model, String elapsed);

  /// No description provided for @geneTreesProgress.
  ///
  /// In en, this message translates to:
  /// **'Gene trees'**
  String get geneTreesProgress;

  /// No description provided for @localFooter.
  ///
  /// In en, this message translates to:
  /// **'Runs on this computer.'**
  String get localFooter;

  /// No description provided for @treeFooter.
  ///
  /// In en, this message translates to:
  /// **'Runs on this computer. Long analyses can take minutes or hours; you can follow each stage here.'**
  String get treeFooter;

  /// No description provided for @treeDetailPrepare.
  ///
  /// In en, this message translates to:
  /// **'Reading the alignment and partitions'**
  String get treeDetailPrepare;

  /// No description provided for @treeDetailModelStart.
  ///
  /// In en, this message translates to:
  /// **'Testing substitution models'**
  String get treeDetailModelStart;

  /// No description provided for @treeDetailModels.
  ///
  /// In en, this message translates to:
  /// **'Model {done} of up to {total} tested'**
  String treeDetailModels(String done, String total);

  /// No description provided for @treeDetailStarts.
  ///
  /// In en, this message translates to:
  /// **'ML search from starting tree {done} of {total}'**
  String treeDetailStarts(String done, String total);

  /// No description provided for @treeDetailSearchStart.
  ///
  /// In en, this message translates to:
  /// **'Building starting trees'**
  String get treeDetailSearchStart;

  /// No description provided for @treeDetailIteration.
  ///
  /// In en, this message translates to:
  /// **'Tree search: iteration {n}'**
  String treeDetailIteration(String n);

  /// No description provided for @treeDetailLeft.
  ///
  /// In en, this message translates to:
  /// **'~{time} left'**
  String treeDetailLeft(String time);

  /// No description provided for @treeDetailBoot.
  ///
  /// In en, this message translates to:
  /// **'Bootstrap replicate {done} of {total}'**
  String treeDetailBoot(String done, String total);

  /// No description provided for @treeDetailSupport.
  ///
  /// In en, this message translates to:
  /// **'Computing branch support'**
  String get treeDetailSupport;

  /// No description provided for @treeDetailFinish.
  ///
  /// In en, this message translates to:
  /// **'Writing the results'**
  String get treeDetailFinish;

  /// No description provided for @treeDetailBestModel.
  ///
  /// In en, this message translates to:
  /// **'Best model: {model}'**
  String treeDetailBestModel(String model);

  /// No description provided for @treeDetailUfboot.
  ///
  /// In en, this message translates to:
  /// **'{n} UFBoot replicates'**
  String treeDetailUfboot(String n);

  /// No description provided for @treeDetailConverged.
  ///
  /// In en, this message translates to:
  /// **'Bootstrap converged after {n}'**
  String treeDetailConverged(String n);

  /// No description provided for @treeStagePrepare.
  ///
  /// In en, this message translates to:
  /// **'Alignment'**
  String get treeStagePrepare;

  /// No description provided for @treeStageModel.
  ///
  /// In en, this message translates to:
  /// **'Model selection'**
  String get treeStageModel;

  /// No description provided for @treeStageSearch.
  ///
  /// In en, this message translates to:
  /// **'ML tree search'**
  String get treeStageSearch;

  /// No description provided for @treeStageSupport.
  ///
  /// In en, this message translates to:
  /// **'Branch support'**
  String get treeStageSupport;

  /// No description provided for @treeStageBootstrap.
  ///
  /// In en, this message translates to:
  /// **'Bootstrap'**
  String get treeStageBootstrap;

  /// No description provided for @treeStageFinish.
  ///
  /// In en, this message translates to:
  /// **'Results'**
  String get treeStageFinish;

  /// No description provided for @step9.
  ///
  /// In en, this message translates to:
  /// **'Tree Viewer'**
  String get step9;

  /// No description provided for @metaViewable.
  ///
  /// In en, this message translates to:
  /// **'{n} tree(s)'**
  String metaViewable(String n);

  /// No description provided for @viewTree.
  ///
  /// In en, this message translates to:
  /// **'View tree'**
  String get viewTree;

  /// No description provided for @treeViewIntro.
  ///
  /// In en, this message translates to:
  /// **'Look at the inferred trees. Root them on the outgroup, at the midpoint or on any branch, check the support of each clade, find species and export the figure.'**
  String get treeViewIntro;

  /// No description provided for @treeChoose.
  ///
  /// In en, this message translates to:
  /// **'Tree'**
  String get treeChoose;

  /// No description provided for @treeMissing.
  ///
  /// In en, this message translates to:
  /// **'The tree file was not found.'**
  String get treeMissing;

  /// No description provided for @treeKindMl.
  ///
  /// In en, this message translates to:
  /// **'ML tree with support'**
  String get treeKindMl;

  /// No description provided for @treeKindConsensus.
  ///
  /// In en, this message translates to:
  /// **'bootstrap consensus'**
  String get treeKindConsensus;

  /// No description provided for @treeKindFbp.
  ///
  /// In en, this message translates to:
  /// **'ML tree, bootstrap support (FBP)'**
  String get treeKindFbp;

  /// No description provided for @treeKindTbe.
  ///
  /// In en, this message translates to:
  /// **'ML tree, transfer bootstrap (TBE)'**
  String get treeKindTbe;

  /// No description provided for @treeKindBest.
  ///
  /// In en, this message translates to:
  /// **'best ML tree (no support)'**
  String get treeKindBest;

  /// No description provided for @treeKindGene.
  ///
  /// In en, this message translates to:
  /// **'gene tree {gene}'**
  String treeKindGene(String gene);

  /// No description provided for @treeParseError.
  ///
  /// In en, this message translates to:
  /// **'This tree could not be read: {error}'**
  String treeParseError(String error);

  /// No description provided for @treeRootAsIs.
  ///
  /// In en, this message translates to:
  /// **'As inferred'**
  String get treeRootAsIs;

  /// No description provided for @treeRootOutgroup.
  ///
  /// In en, this message translates to:
  /// **'Outgroup'**
  String get treeRootOutgroup;

  /// No description provided for @treeRootMidpoint.
  ///
  /// In en, this message translates to:
  /// **'Midpoint'**
  String get treeRootMidpoint;

  /// No description provided for @treeRootManual.
  ///
  /// In en, this message translates to:
  /// **'Rooted on the chosen branch'**
  String get treeRootManual;

  /// No description provided for @treePhylogram.
  ///
  /// In en, this message translates to:
  /// **'Phylogram'**
  String get treePhylogram;

  /// No description provided for @treeCladogram.
  ///
  /// In en, this message translates to:
  /// **'Cladogram'**
  String get treeCladogram;

  /// No description provided for @treeShowSupport.
  ///
  /// In en, this message translates to:
  /// **'Support'**
  String get treeShowSupport;

  /// No description provided for @treeItalic.
  ///
  /// In en, this message translates to:
  /// **'Italic'**
  String get treeItalic;

  /// No description provided for @treeSearch.
  ///
  /// In en, this message translates to:
  /// **'Find a species'**
  String get treeSearch;

  /// No description provided for @treeFont.
  ///
  /// In en, this message translates to:
  /// **'Font size'**
  String get treeFont;

  /// No description provided for @treeStats.
  ///
  /// In en, this message translates to:
  /// **'{tips} tips · tree length {length}'**
  String treeStats(String tips, String length);

  /// No description provided for @treeHint.
  ///
  /// In en, this message translates to:
  /// **'Click a branch or a name to select it: root the tree there, rotate, collapse or shorten it. Ctrl+Z undoes.'**
  String get treeHint;

  /// No description provided for @treeClade.
  ///
  /// In en, this message translates to:
  /// **'Clade with {n} tips'**
  String treeClade(String n);

  /// No description provided for @treeSupportValue.
  ///
  /// In en, this message translates to:
  /// **'support {value}'**
  String treeSupportValue(String value);

  /// No description provided for @treeRootHere.
  ///
  /// In en, this message translates to:
  /// **'Root here'**
  String get treeRootHere;

  /// No description provided for @treeCollapse.
  ///
  /// In en, this message translates to:
  /// **'Collapse'**
  String get treeCollapse;

  /// No description provided for @treeExpand.
  ///
  /// In en, this message translates to:
  /// **'Expand'**
  String get treeExpand;

  /// No description provided for @treeClearSelection.
  ///
  /// In en, this message translates to:
  /// **'Clear'**
  String get treeClearSelection;

  /// No description provided for @treeSupportLegend.
  ///
  /// In en, this message translates to:
  /// **'Branch support:'**
  String get treeSupportLegend;

  /// No description provided for @treeSupportInfoUfboot.
  ///
  /// In en, this message translates to:
  /// **'IQ-TREE labels are SH-aLRT/UFBoot (or UFBoot alone); colours use the last value. UFBoot ≥ 95 is strong support (Hoang et al. 2018); SH-aLRT ≥ 80 together with it makes the clade well supported.'**
  String get treeSupportInfoUfboot;

  /// No description provided for @treeSupportInfoBootstrap.
  ///
  /// In en, this message translates to:
  /// **'Standard (non-parametric) bootstrap: ≥ 70 is usually taken as good support (Hillis & Bull 1993).'**
  String get treeSupportInfoBootstrap;

  /// No description provided for @siSelfTag.
  ///
  /// In en, this message translates to:
  /// **'this record'**
  String get siSelfTag;

  /// No description provided for @siSelfInfo.
  ///
  /// In en, this message translates to:
  /// **'BLAST found this genome\'s own sequence in GenBank ({accessions}). The data is there, but a record cannot confirm its own name, so it does not count as evidence. Other records of the species would.'**
  String siSelfInfo(String accessions);

  /// No description provided for @siVerdictSelfOnly.
  ///
  /// In en, this message translates to:
  /// **'GenBank has this genome\'s sequence ({accessions}) but no other record of {current}, so GenBank cannot confirm the name.'**
  String siVerdictSelfOnly(String accessions, String current);

  /// No description provided for @siVerdictNoReference.
  ///
  /// In en, this message translates to:
  /// **'No other sequence of {current} was found. The closest is {best} ({identity}%), below 98% — the usual distance between species — so the current name is not contested.'**
  String siVerdictNoReference(String current, String best, String identity);

  /// No description provided for @siNoReference.
  ///
  /// In en, this message translates to:
  /// **'No reference'**
  String get siNoReference;

  /// No description provided for @apiKeyEstimateColSize.
  ///
  /// In en, this message translates to:
  /// **'Download'**
  String get apiKeyEstimateColSize;

  /// No description provided for @apiKeyEstimateColGain.
  ///
  /// In en, this message translates to:
  /// **'Gain'**
  String get apiKeyEstimateColGain;

  /// No description provided for @apiKeyEstimateSaved.
  ///
  /// In en, this message translates to:
  /// **'For {n} genomes the key saves about {saved}.'**
  String apiKeyEstimateSaved(String n, String saved);

  /// No description provided for @missTypeNotSelected.
  ///
  /// In en, this message translates to:
  /// **'annotated as {type}'**
  String missTypeNotSelected(String type);

  /// No description provided for @missTypeNotSelectedInfo.
  ///
  /// In en, this message translates to:
  /// **'{gene} is annotated in this record as {type}, a feature type not selected in step 3. Select {type} there to include it.'**
  String missTypeNotSelectedInfo(String gene, String type);

  /// No description provided for @missGeneOnly.
  ///
  /// In en, this message translates to:
  /// **'only a \'gene\' feature'**
  String get missGeneOnly;

  /// No description provided for @missGeneOnlyInfo.
  ///
  /// In en, this message translates to:
  /// **'The record has a \'gene\' feature for {gene} ({location}) but no CDS/rRNA feature, which is what SPLACE reads. The gene is in the sequence; the annotation is incomplete. Load a better-annotated record, or annotate the genome (e.g. with MitoS or MitoFinder).'**
  String missGeneOnlyInfo(String gene, String location);

  /// No description provided for @missNoAnnotation.
  ///
  /// In en, this message translates to:
  /// **'no {types} annotated'**
  String missNoAnnotation(String types);

  /// No description provided for @missNoAnnotationInfo.
  ///
  /// In en, this message translates to:
  /// **'This record has no {types} annotated at all, so no gene can be taken from it. It is probably a raw submission; annotate the genome (e.g. with MitoS or MitoFinder) or remove it.'**
  String missNoAnnotationInfo(String types);

  /// No description provided for @missUnknownName.
  ///
  /// In en, this message translates to:
  /// **'unrecognized name'**
  String get missUnknownName;

  /// No description provided for @missUnknownNameInfo.
  ///
  /// In en, this message translates to:
  /// **'The record has genes with names SPLACE did not recognize ({names}). One of them may be {gene}: fix the name in step 3 (unknown names).'**
  String missUnknownNameInfo(String names, String gene);

  /// No description provided for @missPartial.
  ///
  /// In en, this message translates to:
  /// **'partial genome'**
  String get missPartial;

  /// No description provided for @missPartialInfo.
  ///
  /// In en, this message translates to:
  /// **'The record has {length} bp, well below the typical genome of this dataset ({typical} bp): it is a partial genome and does not cover {gene}.'**
  String missPartialInfo(String length, String typical, String gene);

  /// No description provided for @missNotAnnotated.
  ///
  /// In en, this message translates to:
  /// **'not annotated'**
  String get missNotAnnotated;

  /// No description provided for @missNotAnnotatedInfo.
  ///
  /// In en, this message translates to:
  /// **'{gene} is not among the genes annotated in this record: it is absent from the sequence or the authors did not annotate it.'**
  String missNotAnnotatedInfo(String gene);

  /// No description provided for @missNotInFiles.
  ///
  /// In en, this message translates to:
  /// **'not in the files'**
  String get missNotInFiles;

  /// No description provided for @missNotInFilesInfo.
  ///
  /// In en, this message translates to:
  /// **'No FASTA file has a {gene} sequence for this sample.'**
  String missNotInFilesInfo(String gene);

  /// No description provided for @missPresentIn.
  ///
  /// In en, this message translates to:
  /// **'in {present} of {total} samples · missing in {missing}'**
  String missPresentIn(String present, String total, String missing);

  /// No description provided for @missSampleCount.
  ///
  /// In en, this message translates to:
  /// **'missing {missing} of {total} markers'**
  String missSampleCount(String missing, String total);

  /// No description provided for @missSummary.
  ///
  /// In en, this message translates to:
  /// **'{markers} of the {selected} selected markers are missing in at least one sample; {samples} of {total} samples lack at least one marker.'**
  String missSummary(
    String markers,
    String selected,
    String samples,
    String total,
  );

  /// No description provided for @missMatrixNote.
  ///
  /// In en, this message translates to:
  /// **'In the supermatrix (step 6) a missing marker becomes \'?\' for that sample when “Allow missing sequences” is on; otherwise markers absent from any taxon are left out. A few gaps are normal; a sample missing most markers adds little and may be worth removing.'**
  String get missMatrixNote;

  /// No description provided for @missFixableNote.
  ///
  /// In en, this message translates to:
  /// **'{n} gap(s) can be fixed in SPLACE (orange tags). Hover a tag to see how.'**
  String missFixableNote(String n);

  /// No description provided for @missByMarker.
  ///
  /// In en, this message translates to:
  /// **'By marker'**
  String get missByMarker;

  /// No description provided for @missBySample.
  ///
  /// In en, this message translates to:
  /// **'By sample'**
  String get missBySample;

  /// No description provided for @missHoverHint.
  ///
  /// In en, this message translates to:
  /// **'Each tag says why the marker is missing; hover it for details.'**
  String get missHoverHint;

  /// No description provided for @dupSummary.
  ///
  /// In en, this message translates to:
  /// **'{genes} gene(s) have more than one copy annotated in {genomes} genome(s).'**
  String dupSummary(String genes, String genomes);

  /// No description provided for @dupFastaHelp.
  ///
  /// In en, this message translates to:
  /// **'These samples have more than one sequence of the same locus in the FASTA files. The “Used” column shows what goes into the alignment, following the “Duplicates and names” setting.'**
  String get dupFastaHelp;

  /// No description provided for @dupIdentical.
  ///
  /// In en, this message translates to:
  /// **'identical copies'**
  String get dupIdentical;

  /// No description provided for @dupIdenticalInfo.
  ///
  /// In en, this message translates to:
  /// **'Same sequence (e.g. inverted-repeat copies in plastomes): any copy gives the same result.'**
  String get dupIdenticalInfo;

  /// No description provided for @dupContained.
  ///
  /// In en, this message translates to:
  /// **'fragment of another copy'**
  String get dupContained;

  /// No description provided for @dupContainedInfo.
  ///
  /// In en, this message translates to:
  /// **'The shorter copy is part of the longer one: usually a duplicated or truncated annotation. Use the complete copy.'**
  String get dupContainedInfo;

  /// No description provided for @dupDifferent.
  ///
  /// In en, this message translates to:
  /// **'different copies'**
  String get dupDifferent;

  /// No description provided for @dupDifferentInfo.
  ///
  /// In en, this message translates to:
  /// **'The copies differ in sequence (paralogues, NUMT-like copies or annotation errors). The choice changes the alignment: prefer the copy closest to the other genomes.'**
  String get dupDifferentInfo;

  /// No description provided for @siCardTitle.
  ///
  /// In en, this message translates to:
  /// **'Species identification'**
  String get siCardTitle;

  /// No description provided for @siCardBody.
  ///
  /// In en, this message translates to:
  /// **'Compares a barcode gene (COI, 16S…) of each genome with GenBank (BLAST) and BOLD to catch misidentified records before fetching the taxonomy, which depends on the names.'**
  String get siCardBody;

  /// No description provided for @siCardOptional.
  ///
  /// In en, this message translates to:
  /// **'Optional'**
  String get siCardOptional;

  /// No description provided for @siCardRerun.
  ///
  /// In en, this message translates to:
  /// **'Identify again'**
  String get siCardRerun;

  /// No description provided for @siSumChecked.
  ///
  /// In en, this message translates to:
  /// **'{n} genomes checked'**
  String siSumChecked(String n);

  /// No description provided for @siSumConfirmed.
  ///
  /// In en, this message translates to:
  /// **'{n} names confirmed'**
  String siSumConfirmed(String n);

  /// No description provided for @siSumFlagged.
  ///
  /// In en, this message translates to:
  /// **'{n} possibly misidentified'**
  String siSumFlagged(String n);

  /// No description provided for @siSumNoReference.
  ///
  /// In en, this message translates to:
  /// **'{n} without a reference'**
  String siSumNoReference(String n);

  /// No description provided for @headerDefault.
  ///
  /// In en, this message translates to:
  /// **'Default header · used when aligning'**
  String get headerDefault;

  /// No description provided for @viewerMafftColumn.
  ///
  /// In en, this message translates to:
  /// **'column {n} of the MAFFT alignment'**
  String viewerMafftColumn(String n);

  /// No description provided for @satSourceMafft.
  ///
  /// In en, this message translates to:
  /// **'Tested on the MAFFT alignment (not trimmed).'**
  String get satSourceMafft;

  /// No description provided for @satSourceClipkit.
  ///
  /// In en, this message translates to:
  /// **'Tested on the trimmed alignment (ClipKIT on the proteins, back-translated): only whole codons were removed, so every base keeps its codon position.'**
  String get satSourceClipkit;

  /// No description provided for @satSourceTrimmed.
  ///
  /// In en, this message translates to:
  /// **'Tested on the trimmed alignment: trimAl kept {kept} of the {total} MAFFT codons, and SPLACE removed only whole codons, so every column keeps the codon position it had in the MAFFT alignment (hover a trimmed column in step 4 to see which MAFFT column it is).'**
  String satSourceTrimmed(String kept, String total);

  /// No description provided for @satSourceFrame.
  ///
  /// In en, this message translates to:
  /// **'Codons now start at base 1 because the incomplete codon before base {frame} of the MAFFT alignment was dropped; the positions are the same.'**
  String satSourceFrame(String frame);

  /// No description provided for @preserveCaseHelp.
  ///
  /// In en, this message translates to:
  /// **'Keeps the input letter case'**
  String get preserveCaseHelp;

  /// No description provided for @configureTrimming.
  ///
  /// In en, this message translates to:
  /// **'Set up trimming'**
  String get configureTrimming;

  /// No description provided for @panelTrimParams.
  ///
  /// In en, this message translates to:
  /// **'{tool} parameters'**
  String panelTrimParams(String tool);

  /// No description provided for @trimmedCount.
  ///
  /// In en, this message translates to:
  /// **'{n} of {total} markers trimmed'**
  String trimmedCount(String n, String total);

  /// No description provided for @showRemovedColumns.
  ///
  /// In en, this message translates to:
  /// **'Show removed regions'**
  String get showRemovedColumns;

  /// No description provided for @nextStop.
  ///
  /// In en, this message translates to:
  /// **'Next stop'**
  String get nextStop;

  /// No description provided for @translationTable.
  ///
  /// In en, this message translates to:
  /// **'Translation table'**
  String get translationTable;

  /// No description provided for @translationTableTitle.
  ///
  /// In en, this message translates to:
  /// **'Genetic code {code} · {name}'**
  String translationTableTitle(String code, String name);

  /// No description provided for @translationTableBody.
  ///
  /// In en, this message translates to:
  /// **'Codons of the genetic code chosen above for the amino-acid view (it is never guessed from the data). Stop codons are red; the badge counts the internal stops of this alignment that use each one.'**
  String get translationTableBody;

  /// No description provided for @codonSecondBase.
  ///
  /// In en, this message translates to:
  /// **'2nd {base}'**
  String codonSecondBase(String base);

  /// No description provided for @codonDiffers.
  ///
  /// In en, this message translates to:
  /// **'{codon} → {aa} in this code ({standard} in the standard code)'**
  String codonDiffers(String codon, String aa, String standard);

  /// No description provided for @codonDiffersLegend.
  ///
  /// In en, this message translates to:
  /// **'Differs from the standard code'**
  String get codonDiffersLegend;

  /// No description provided for @stopCodonsFound.
  ///
  /// In en, this message translates to:
  /// **'Internal stops in this alignment: {list}'**
  String stopCodonsFound(String list);

  /// No description provided for @removedColumns.
  ///
  /// In en, this message translates to:
  /// **'{n} of {total} columns removed'**
  String removedColumns(String n, String total);

  /// No description provided for @removedColumnsInfo.
  ///
  /// In en, this message translates to:
  /// **'The whole MAFFT alignment is shown; grey columns were removed by trimming. Turn off \"Show removed regions\" to see only the kept columns.'**
  String get removedColumnsInfo;

  /// No description provided for @removedLegend.
  ///
  /// In en, this message translates to:
  /// **'Removed'**
  String get removedLegend;

  /// No description provided for @viewerRemoved.
  ///
  /// In en, this message translates to:
  /// **'removed by trimming'**
  String get viewerRemoved;

  /// No description provided for @viewerCodon.
  ///
  /// In en, this message translates to:
  /// **'codon {codon}'**
  String viewerCodon(String codon);

  /// No description provided for @plotObservedVsCorrectedBody.
  ///
  /// In en, this message translates to:
  /// **'Each point is a pair of sequences. Points on the dashed line: every substitution is visible; points levelling off below it: multiple hits hide substitutions (saturation).'**
  String get plotObservedVsCorrectedBody;

  /// No description provided for @plotTsTvBody.
  ///
  /// In en, this message translates to:
  /// **'Transitions (red) usually saturate first: when they stop growing while transversions (blue) keep rising, the subset is saturated.'**
  String get plotTsTvBody;

  /// No description provided for @plotPositions.
  ///
  /// In en, this message translates to:
  /// **'p of each codon position × {model} distance'**
  String plotPositions(String model);

  /// No description provided for @plotPositionsBody.
  ///
  /// In en, this message translates to:
  /// **'Same pairs, one colour per codon position, against the distance of the whole gene. A position whose points flatten while the others keep rising is saturated (usually the 3rd).'**
  String get plotPositionsBody;

  /// No description provided for @allPositions.
  ///
  /// In en, this message translates to:
  /// **'all positions'**
  String get allPositions;

  /// No description provided for @compositionTitle.
  ///
  /// In en, this message translates to:
  /// **'Composition and variability'**
  String get compositionTitle;

  /// No description provided for @compositionInfo.
  ///
  /// In en, this message translates to:
  /// **'Chi-square test of base-frequency homogeneity across taxa (contingency test, as in PAUP*), and the number of taxa whose own composition differs from the subset (df = 3, p < 0.05, as IQ-TREE reports). Saturated positions often drift in composition (e.g. GC3), which can also bias the tree.'**
  String get compositionInfo;

  /// No description provided for @metricVariable.
  ///
  /// In en, this message translates to:
  /// **'Variable'**
  String get metricVariable;

  /// No description provided for @metricInformative.
  ///
  /// In en, this message translates to:
  /// **'Parsimony-informative'**
  String get metricInformative;

  /// No description provided for @metricTaxaFailing.
  ///
  /// In en, this message translates to:
  /// **'Taxa failing'**
  String get metricTaxaFailing;

  /// No description provided for @uncorrectablePairs.
  ///
  /// In en, this message translates to:
  /// **'pairs beyond correction'**
  String get uncorrectablePairs;

  /// No description provided for @heatmapIntensity.
  ///
  /// In en, this message translates to:
  /// **'The shade inside each colour follows the slope and the uncorrectable pairs.'**
  String get heatmapIntensity;

  /// No description provided for @treeRootLabel.
  ///
  /// In en, this message translates to:
  /// **'Root'**
  String get treeRootLabel;

  /// No description provided for @treeShapeLabel.
  ///
  /// In en, this message translates to:
  /// **'Shape'**
  String get treeShapeLabel;

  /// No description provided for @treeSupportValues.
  ///
  /// In en, this message translates to:
  /// **'Values'**
  String get treeSupportValues;

  /// No description provided for @treeSupportDots.
  ///
  /// In en, this message translates to:
  /// **'Dots'**
  String get treeSupportDots;

  /// No description provided for @treeSupportHidden.
  ///
  /// In en, this message translates to:
  /// **'Hide'**
  String get treeSupportHidden;

  /// No description provided for @treeDisplay.
  ///
  /// In en, this message translates to:
  /// **'Display'**
  String get treeDisplay;

  /// No description provided for @treeAlignLabels.
  ///
  /// In en, this message translates to:
  /// **'Aligned on the right'**
  String get treeAlignLabels;

  /// No description provided for @treeExport.
  ///
  /// In en, this message translates to:
  /// **'Export'**
  String get treeExport;

  /// No description provided for @treeMatches.
  ///
  /// In en, this message translates to:
  /// **'{n} found'**
  String treeMatches(String n);

  /// No description provided for @treeHoverTip.
  ///
  /// In en, this message translates to:
  /// **'{name} · branch length {length}'**
  String treeHoverTip(String name, String length);

  /// No description provided for @treeHoverClade.
  ///
  /// In en, this message translates to:
  /// **'Clade with {n} tips · support {support} · branch length {length} — click to select it'**
  String treeHoverClade(String n, String support, String length);

  /// No description provided for @treeOutgroupLegend.
  ///
  /// In en, this message translates to:
  /// **'Outgroup'**
  String get treeOutgroupLegend;

  /// No description provided for @ord1.
  ///
  /// In en, this message translates to:
  /// **'1st'**
  String get ord1;

  /// No description provided for @ord2.
  ///
  /// In en, this message translates to:
  /// **'2nd'**
  String get ord2;

  /// No description provided for @ord3.
  ///
  /// In en, this message translates to:
  /// **'3rd'**
  String get ord3;

  /// No description provided for @listAnd.
  ///
  /// In en, this message translates to:
  /// **'and'**
  String get listAnd;

  /// No description provided for @codonPositionN.
  ///
  /// In en, this message translates to:
  /// **'{n} position'**
  String codonPositionN(String n);

  /// No description provided for @reasonInsufficient.
  ///
  /// In en, this message translates to:
  /// **'Too little data to judge ({sites} sites, {pairs} sequence pairs).'**
  String reasonInsufficient(String sites, String pairs);

  /// No description provided for @reasonUndefined.
  ///
  /// In en, this message translates to:
  /// **'{pct}% of the sequence pairs differ so much that the distance cannot be corrected.'**
  String reasonUndefined(String pct);

  /// No description provided for @reasonLowDivergence.
  ///
  /// In en, this message translates to:
  /// **'The sequences differ very little (mean p = {p}), so saturation is not a concern.'**
  String reasonLowDivergence(String p);

  /// No description provided for @reasonSlopeOk.
  ///
  /// In en, this message translates to:
  /// **'PhyKIT saturation {slope}: the observed differences still grow with the tree distance.'**
  String reasonSlopeOk(String slope);

  /// No description provided for @reasonSlopeModerate.
  ///
  /// In en, this message translates to:
  /// **'PhyKIT saturation {slope} (below {limit}): the observed differences grow more slowly than the tree distance.'**
  String reasonSlopeModerate(String slope, String limit);

  /// No description provided for @reasonSlopeSaturated.
  ///
  /// In en, this message translates to:
  /// **'PhyKIT saturation {slope} (below {limit}): the observed differences have almost stopped growing with the tree distance; many changes were overwritten by later ones.'**
  String reasonSlopeSaturated(String slope, String limit);

  /// No description provided for @reasonIss.
  ///
  /// In en, this message translates to:
  /// **'Iss {iss} (0.8 or more): the sites look close to random.'**
  String reasonIss(String iss);

  /// No description provided for @recTitle.
  ///
  /// In en, this message translates to:
  /// **'What to do with the codon positions'**
  String get recTitle;

  /// No description provided for @recBasedOn.
  ///
  /// In en, this message translates to:
  /// **'Based on the concatenated matrix ({taxa} taxa, {length} bp). Each position is tested on its own.'**
  String recBasedOn(String taxa, String length);

  /// No description provided for @recWhatToDo.
  ///
  /// In en, this message translates to:
  /// **'Suggestion'**
  String get recWhatToDo;

  /// No description provided for @recExcludeOne.
  ///
  /// In en, this message translates to:
  /// **'Remove the {positions} codon position from the matrix (step 6b) and keep the other two, partitioned by codon position.'**
  String recExcludeOne(String positions);

  /// No description provided for @recExcludeMany.
  ///
  /// In en, this message translates to:
  /// **'Remove the {positions} codon positions from the matrix (step 6b) and keep the remaining one.'**
  String recExcludeMany(String positions);

  /// No description provided for @recWhy.
  ///
  /// In en, this message translates to:
  /// **'Why: a saturated position has changed so many times that its differences no longer show how long ago the lineages split. It mostly adds noise (homoplasy) and can pull long branches together.'**
  String get recWhy;

  /// No description provided for @recCautionOne.
  ///
  /// In en, this message translates to:
  /// **'The {positions} position shows moderate saturation: keep it in its own partition and check whether the tree changes without it.'**
  String recCautionOne(String positions);

  /// No description provided for @recCautionMany.
  ///
  /// In en, this message translates to:
  /// **'The {positions} positions show moderate saturation: keep them in their own partitions and check whether the tree changes without them.'**
  String recCautionMany(String positions);

  /// No description provided for @applyExclude.
  ///
  /// In en, this message translates to:
  /// **'Mark {positions} in step 6b'**
  String applyExclude(String positions);

  /// No description provided for @applyKeepAll.
  ///
  /// In en, this message translates to:
  /// **'Leave every \"Exclude position\" unmarked in step 6b'**
  String get applyKeepAll;

  /// No description provided for @applyPerLocusInfo.
  ///
  /// In en, this message translates to:
  /// **'Each gene excludes only the positions saturated in it. Genes whose result differs from the overall suggestion: {loci}.'**
  String applyPerLocusInfo(String loci);

  /// No description provided for @recLociDiffer.
  ///
  /// In en, this message translates to:
  /// **'{n} gene(s) have a result different from the overall suggestion; see their rows in the heatmap.'**
  String recLociDiffer(String n);

  /// No description provided for @recStatusKeep.
  ///
  /// In en, this message translates to:
  /// **'Step 6b now keeps all positions.'**
  String get recStatusKeep;

  /// No description provided for @recStatusExclude.
  ///
  /// In en, this message translates to:
  /// **'Step 6b now excludes the {positions} position(s).'**
  String recStatusExclude(String positions);

  /// No description provided for @recStatusPerLocus.
  ///
  /// In en, this message translates to:
  /// **'Step 6b now excludes positions per gene ({n} genes).'**
  String recStatusPerLocus(String n);

  /// No description provided for @saturationHeatmapHelp.
  ///
  /// In en, this message translates to:
  /// **'Each cell shows PhyKIT\'s saturation value (slope of p on the tree distance) — 1 means no saturation, near 0 saturated. Blue is OK, red is saturated.'**
  String get saturationHeatmapHelp;

  /// No description provided for @rankSpecies.
  ///
  /// In en, this message translates to:
  /// **'Species'**
  String get rankSpecies;

  /// No description provided for @rankGenus.
  ///
  /// In en, this message translates to:
  /// **'Genus'**
  String get rankGenus;

  /// No description provided for @rankFamily.
  ///
  /// In en, this message translates to:
  /// **'Family'**
  String get rankFamily;

  /// No description provided for @rankOrder.
  ///
  /// In en, this message translates to:
  /// **'Order'**
  String get rankOrder;

  /// No description provided for @rankClass.
  ///
  /// In en, this message translates to:
  /// **'Class'**
  String get rankClass;

  /// No description provided for @rankPhylum.
  ///
  /// In en, this message translates to:
  /// **'Phylum'**
  String get rankPhylum;

  /// No description provided for @rankKingdom.
  ///
  /// In en, this message translates to:
  /// **'Kingdom'**
  String get rankKingdom;

  /// No description provided for @tipFieldOriginal.
  ///
  /// In en, this message translates to:
  /// **'Name as in the tree'**
  String get tipFieldOriginal;

  /// No description provided for @tipFieldAccession.
  ///
  /// In en, this message translates to:
  /// **'Accession'**
  String get tipFieldAccession;

  /// No description provided for @tipFieldVoucher.
  ///
  /// In en, this message translates to:
  /// **'Voucher'**
  String get tipFieldVoucher;

  /// No description provided for @treeBold.
  ///
  /// In en, this message translates to:
  /// **'Bold'**
  String get treeBold;

  /// No description provided for @treeParens.
  ///
  /// In en, this message translates to:
  /// **'In parentheses'**
  String get treeParens;

  /// No description provided for @treeRemovePart.
  ///
  /// In en, this message translates to:
  /// **'Remove this part'**
  String get treeRemovePart;

  /// No description provided for @treeAddPart.
  ///
  /// In en, this message translates to:
  /// **'Add field'**
  String get treeAddPart;

  /// No description provided for @treeMovePart.
  ///
  /// In en, this message translates to:
  /// **'Move left'**
  String get treeMovePart;

  /// No description provided for @treeLabelFormat.
  ///
  /// In en, this message translates to:
  /// **'Label format'**
  String get treeLabelFormat;

  /// No description provided for @treeLabelParts.
  ///
  /// In en, this message translates to:
  /// **'Parts of the label'**
  String get treeLabelParts;

  /// No description provided for @treeLabelPartsInfo.
  ///
  /// In en, this message translates to:
  /// **'Each label is built from these fields, in this order. I = italic, B = bold, ( ) = in parentheses. The values come from the GenBank records and the taxonomy fetched in step 2; a field a tip does not have is left out.'**
  String get treeLabelPartsInfo;

  /// No description provided for @treeTabTree.
  ///
  /// In en, this message translates to:
  /// **'Tree'**
  String get treeTabTree;

  /// No description provided for @treeTabBranches.
  ///
  /// In en, this message translates to:
  /// **'Branches and support'**
  String get treeTabBranches;

  /// No description provided for @treeTabLabels.
  ///
  /// In en, this message translates to:
  /// **'Labels'**
  String get treeTabLabels;

  /// No description provided for @treeTabColours.
  ///
  /// In en, this message translates to:
  /// **'Colours'**
  String get treeTabColours;

  /// No description provided for @treeOrderLabel.
  ///
  /// In en, this message translates to:
  /// **'Order and rotation'**
  String get treeOrderLabel;

  /// No description provided for @treeLadderUp.
  ///
  /// In en, this message translates to:
  /// **'Ladderize ↑'**
  String get treeLadderUp;

  /// No description provided for @treeLadderDown.
  ///
  /// In en, this message translates to:
  /// **'Ladderize ↓'**
  String get treeLadderDown;

  /// No description provided for @treeFlip.
  ///
  /// In en, this message translates to:
  /// **'Flip tree'**
  String get treeFlip;

  /// No description provided for @treeRotateClade.
  ///
  /// In en, this message translates to:
  /// **'Rotate clade'**
  String get treeRotateClade;

  /// No description provided for @treeRotateInfo.
  ///
  /// In en, this message translates to:
  /// **'Rotating a clade swaps the order of its branches; the relationships do not change. Select a clade by clicking its branch, then rotate it. Flip tree mirrors the whole tree.'**
  String get treeRotateInfo;

  /// No description provided for @treeLineWidth.
  ///
  /// In en, this message translates to:
  /// **'Line width'**
  String get treeLineWidth;

  /// No description provided for @treePointSize.
  ///
  /// In en, this message translates to:
  /// **'Point size'**
  String get treePointSize;

  /// No description provided for @treePoints.
  ///
  /// In en, this message translates to:
  /// **'Points'**
  String get treePoints;

  /// No description provided for @treeTipDots.
  ///
  /// In en, this message translates to:
  /// **'Points on the tips'**
  String get treeTipDots;

  /// No description provided for @treeSupportMin.
  ///
  /// In en, this message translates to:
  /// **'Show support from'**
  String get treeSupportMin;

  /// No description provided for @treeSupportMinInfo.
  ///
  /// In en, this message translates to:
  /// **'Only nodes with support at or above this value get a number or a point (for example, points only on nodes with UFBoot ≥ 95).'**
  String get treeSupportMinInfo;

  /// No description provided for @treeSupportAll.
  ///
  /// In en, this message translates to:
  /// **'All'**
  String get treeSupportAll;

  /// No description provided for @treeLongBranches.
  ///
  /// In en, this message translates to:
  /// **'Long branches'**
  String get treeLongBranches;

  /// No description provided for @treeLongBranchesInfo.
  ///
  /// In en, this message translates to:
  /// **'A shortened branch is drawn shorter than it is and marked with //, so the rest of the tree has room. Branches over 5× the median length and over 30% of the tree height count as long. You can also select any branch and shorten it. The real length stays in the info bar and in the exported Newick.'**
  String get treeLongBranchesInfo;

  /// No description provided for @treeShortenLong.
  ///
  /// In en, this message translates to:
  /// **'Shorten long branches ({n})'**
  String treeShortenLong(String n);

  /// No description provided for @treeRestoreLengths.
  ///
  /// In en, this message translates to:
  /// **'Restore all lengths'**
  String get treeRestoreLengths;

  /// No description provided for @treeShortenBranch.
  ///
  /// In en, this message translates to:
  /// **'Shorten branch //'**
  String get treeShortenBranch;

  /// No description provided for @treeRestoreBranch.
  ///
  /// In en, this message translates to:
  /// **'Restore length'**
  String get treeRestoreBranch;

  /// No description provided for @treeShortenedLegend.
  ///
  /// In en, this message translates to:
  /// **'shortened branch (not to scale)'**
  String get treeShortenedLegend;

  /// No description provided for @treeColorRows.
  ///
  /// In en, this message translates to:
  /// **'Colour rows by'**
  String get treeColorRows;

  /// No description provided for @treeColorRowsInfo.
  ///
  /// In en, this message translates to:
  /// **'Each taxon of the chosen rank gets a background colour across the whole row, with its name on the right. The ranks come from the taxonomy fetched in step 2 (or the species and genus in the tip names).'**
  String get treeColorRowsInfo;

  /// No description provided for @treeColorNone.
  ///
  /// In en, this message translates to:
  /// **'None'**
  String get treeColorNone;

  /// No description provided for @treeNoTaxonomy.
  ///
  /// In en, this message translates to:
  /// **'No taxonomy for these tips. Fetch it in step 2.'**
  String get treeNoTaxonomy;

  /// No description provided for @treeExportInfo.
  ///
  /// In en, this message translates to:
  /// **'SVG and PNG are the figure as you see it (labels, colours, legend, shortened branches), without the selection or search marks. Newick has the topology and the real branch lengths.'**
  String get treeExportInfo;

  /// No description provided for @treeUndo.
  ///
  /// In en, this message translates to:
  /// **'Undo'**
  String get treeUndo;

  /// No description provided for @treeRedo.
  ///
  /// In en, this message translates to:
  /// **'Redo'**
  String get treeRedo;

  /// No description provided for @treeReset.
  ///
  /// In en, this message translates to:
  /// **'Start over'**
  String get treeReset;

  /// No description provided for @zipTitle.
  ///
  /// In en, this message translates to:
  /// **'Save everything'**
  String get zipTitle;

  /// No description provided for @zipBody.
  ///
  /// In en, this message translates to:
  /// **'One ZIP with the sequences, alignments, saturation tables, matrices, partitions, trees and logs of this project.'**
  String get zipBody;

  /// No description provided for @concatHowTitle.
  ///
  /// In en, this message translates to:
  /// **'How the supermatrix is built'**
  String get concatHowTitle;

  /// No description provided for @concatHow1.
  ///
  /// In en, this message translates to:
  /// **'The alignment of each locus is placed side by side, one row per taxon (matched by the FASTA header).'**
  String get concatHow1;

  /// No description provided for @concatHow2.
  ///
  /// In en, this message translates to:
  /// **'A taxon that lacks a locus gets ? in that stretch (missing data), so every row has the same length.'**
  String get concatHow2;

  /// No description provided for @concatHow3.
  ///
  /// In en, this message translates to:
  /// **'The partition file tells IQ-TREE / RAxML-NG where each locus (and each codon position) starts and ends, so each part can get its own substitution model.'**
  String get concatHow3;

  /// No description provided for @concatPictureCaption.
  ///
  /// In en, this message translates to:
  /// **'? = taxon without that locus'**
  String get concatPictureCaption;

  /// No description provided for @concatBlockSource.
  ///
  /// In en, this message translates to:
  /// **'Which alignments to join'**
  String get concatBlockSource;

  /// No description provided for @concatBlockSourceBody.
  ///
  /// In en, this message translates to:
  /// **'The trimmed alignments of step 4 (without the badly aligned columns) or the MAFFT output as it is.'**
  String get concatBlockSourceBody;

  /// No description provided for @concatBlockPartitions.
  ///
  /// In en, this message translates to:
  /// **'How to split the matrix into partitions'**
  String get concatBlockPartitions;

  /// No description provided for @concatBlockPartitionsBody.
  ///
  /// In en, this message translates to:
  /// **'Each partition gets its own substitution model in the tree search. Codon positions can also be removed here.'**
  String get concatBlockPartitionsBody;

  /// No description provided for @schemeByLocusCodonBody.
  ///
  /// In en, this message translates to:
  /// **'Each locus is a partition and, in protein-coding loci, each codon position (1st, 2nd, 3rd) gets its own. Codon positions evolve at very different rates; with MFP+MERGE, IQ-TREE joins the partitions that behave alike.'**
  String get schemeByLocusCodonBody;

  /// No description provided for @schemeByLocusBody.
  ///
  /// In en, this message translates to:
  /// **'One partition per locus: each gene gets its own model, and its three codon positions share it. Good for rRNAs, non-coding loci or loci with few sites.'**
  String get schemeByLocusBody;

  /// No description provided for @schemeByCodonBody.
  ///
  /// In en, this message translates to:
  /// **'Three partitions shared by all coding loci (all 1st positions together, all 2nd, all 3rd); each non-coding locus stays a partition of its own. Fewer models than by locus and position.'**
  String get schemeByCodonBody;

  /// No description provided for @schemeNoneBody.
  ///
  /// In en, this message translates to:
  /// **'One model for the whole matrix (no partition file). Simplest, but usually fits worse when genes or positions evolve differently.'**
  String get schemeNoneBody;

  /// No description provided for @codonExclusionHelp.
  ///
  /// In en, this message translates to:
  /// **'Removed positions leave the coding loci only. Step 5 (saturation) can fill this in.'**
  String get codonExclusionHelp;

  /// No description provided for @concatSummaryTitle.
  ///
  /// In en, this message translates to:
  /// **'What will be built'**
  String get concatSummaryTitle;

  /// No description provided for @concatSummaryLoci.
  ///
  /// In en, this message translates to:
  /// **'Loci'**
  String get concatSummaryLoci;

  /// No description provided for @concatSummarySource.
  ///
  /// In en, this message translates to:
  /// **'Alignments'**
  String get concatSummarySource;

  /// No description provided for @concatNoOutgroup.
  ///
  /// In en, this message translates to:
  /// **'none (midpoint)'**
  String get concatNoOutgroup;

  /// No description provided for @concatNothingExcluded.
  ///
  /// In en, this message translates to:
  /// **'nothing'**
  String get concatNothingExcluded;

  /// No description provided for @treeXScale.
  ///
  /// In en, this message translates to:
  /// **'Horizontal scale'**
  String get treeXScale;

  /// No description provided for @treeXScaleInfo.
  ///
  /// In en, this message translates to:
  /// **'Stretches (above 1) or shrinks (below 1) the branches on the x axis. Above 1 the figure grows wider than the set width.'**
  String get treeXScaleInfo;

  /// No description provided for @treeYScale.
  ///
  /// In en, this message translates to:
  /// **'Vertical scale'**
  String get treeYScale;

  /// No description provided for @treeYScaleInfo.
  ///
  /// In en, this message translates to:
  /// **'Space between the rows on the y axis (1 = 1.7 times the label size).'**
  String get treeYScaleInfo;

  /// No description provided for @treeBranchLengths.
  ///
  /// In en, this message translates to:
  /// **'Branch lengths'**
  String get treeBranchLengths;

  /// No description provided for @treeShowLengths.
  ///
  /// In en, this message translates to:
  /// **'Show as labels'**
  String get treeShowLengths;

  /// No description provided for @treeSupportFont.
  ///
  /// In en, this message translates to:
  /// **'Support text size'**
  String get treeSupportFont;

  /// No description provided for @treeProportional.
  ///
  /// In en, this message translates to:
  /// **'Size by support value'**
  String get treeProportional;

  /// No description provided for @treeShowLabels.
  ///
  /// In en, this message translates to:
  /// **'Show tip labels'**
  String get treeShowLabels;

  /// No description provided for @treeTabFigure.
  ///
  /// In en, this message translates to:
  /// **'Figure and export'**
  String get treeTabFigure;

  /// No description provided for @treeScaleBar.
  ///
  /// In en, this message translates to:
  /// **'Scale bar'**
  String get treeScaleBar;

  /// No description provided for @treeShowScaleBar.
  ///
  /// In en, this message translates to:
  /// **'Show'**
  String get treeShowScaleBar;

  /// No description provided for @treeScaleBarLength.
  ///
  /// In en, this message translates to:
  /// **'Scale bar length'**
  String get treeScaleBarLength;

  /// No description provided for @treeScaleBarLengthInfo.
  ///
  /// In en, this message translates to:
  /// **'In substitutions per site. 0 picks a round value about a fifth of the tree.'**
  String get treeScaleBarLengthInfo;

  /// No description provided for @treeFigureLegend.
  ///
  /// In en, this message translates to:
  /// **'Legend'**
  String get treeFigureLegend;

  /// No description provided for @treeLegendInFigure.
  ///
  /// In en, this message translates to:
  /// **'Draw in the figure'**
  String get treeLegendInFigure;

  /// No description provided for @treeFigureLegendInfo.
  ///
  /// In en, this message translates to:
  /// **'Support colours, outgroup, shortened branches and row colours are drawn under the tree, so they go with the PNG and SVG.'**
  String get treeFigureLegendInfo;

  /// No description provided for @treeAtTips.
  ///
  /// In en, this message translates to:
  /// **'At the tips'**
  String get treeAtTips;

  /// No description provided for @treeLegendMarks.
  ///
  /// In en, this message translates to:
  /// **'Marks'**
  String get treeLegendMarks;

  /// No description provided for @treeBranchColor.
  ///
  /// In en, this message translates to:
  /// **'Branch colour'**
  String get treeBranchColor;

  /// No description provided for @treeSupportPosition.
  ///
  /// In en, this message translates to:
  /// **'Support position'**
  String get treeSupportPosition;

  /// No description provided for @treeAboveBranch.
  ///
  /// In en, this message translates to:
  /// **'Above the branch'**
  String get treeAboveBranch;

  /// No description provided for @treeBelowBranch.
  ///
  /// In en, this message translates to:
  /// **'Below the branch'**
  String get treeBelowBranch;

  /// No description provided for @treeSupportColor.
  ///
  /// In en, this message translates to:
  /// **'Colour'**
  String get treeSupportColor;

  /// No description provided for @treeFontFamily.
  ///
  /// In en, this message translates to:
  /// **'Font'**
  String get treeFontFamily;

  /// No description provided for @treeLabelColor.
  ///
  /// In en, this message translates to:
  /// **'Label colour'**
  String get treeLabelColor;

  /// No description provided for @treeLabelPosition.
  ///
  /// In en, this message translates to:
  /// **'Label position'**
  String get treeLabelPosition;

  /// No description provided for @treeLabelGap.
  ///
  /// In en, this message translates to:
  /// **'Distance from the tip'**
  String get treeLabelGap;

  /// No description provided for @treeLabelGapInfo.
  ///
  /// In en, this message translates to:
  /// **'Space between the end of the branch (or the aligned column) and the label.'**
  String get treeLabelGapInfo;

  /// No description provided for @treeColorApply.
  ///
  /// In en, this message translates to:
  /// **'Apply the colours to'**
  String get treeColorApply;

  /// No description provided for @treeColorApplyRows.
  ///
  /// In en, this message translates to:
  /// **'Row background'**
  String get treeColorApplyRows;

  /// No description provided for @treeColorApplyBranches.
  ///
  /// In en, this message translates to:
  /// **'Branches'**
  String get treeColorApplyBranches;

  /// No description provided for @treeColorApplyLabels.
  ///
  /// In en, this message translates to:
  /// **'Labels'**
  String get treeColorApplyLabels;

  /// No description provided for @treeCladeColors.
  ///
  /// In en, this message translates to:
  /// **'Clade colours'**
  String get treeCladeColors;

  /// No description provided for @treeCladeColorsInfo.
  ///
  /// In en, this message translates to:
  /// **'Select a clade by clicking its branch and give it a colour: its branches take it. Undo works here too.'**
  String get treeCladeColorsInfo;

  /// No description provided for @treeColorClade.
  ///
  /// In en, this message translates to:
  /// **'Colour clade'**
  String get treeColorClade;

  /// No description provided for @treeRemoveColor.
  ///
  /// In en, this message translates to:
  /// **'Remove colour'**
  String get treeRemoveColor;

  /// No description provided for @treeClearCladeColors.
  ///
  /// In en, this message translates to:
  /// **'Remove all ({n})'**
  String treeClearCladeColors(String n);

  /// No description provided for @treeBackground.
  ///
  /// In en, this message translates to:
  /// **'Figure background'**
  String get treeBackground;

  /// No description provided for @treeLegendPosition.
  ///
  /// In en, this message translates to:
  /// **'Position'**
  String get treeLegendPosition;

  /// No description provided for @treeLegendBelow.
  ///
  /// In en, this message translates to:
  /// **'Below the tree'**
  String get treeLegendBelow;

  /// No description provided for @treeLegendRight.
  ///
  /// In en, this message translates to:
  /// **'Right of the tree'**
  String get treeLegendRight;

  /// No description provided for @treeLegendStyle.
  ///
  /// In en, this message translates to:
  /// **'Style'**
  String get treeLegendStyle;

  /// No description provided for @treeLegendBox.
  ///
  /// In en, this message translates to:
  /// **'In a box'**
  String get treeLegendBox;

  /// No description provided for @treeLegendTitle.
  ///
  /// In en, this message translates to:
  /// **'Legend title'**
  String get treeLegendTitle;

  /// No description provided for @treeLegendTitleHint.
  ///
  /// In en, this message translates to:
  /// **'e.g. Legend'**
  String get treeLegendTitleHint;

  /// No description provided for @treeLegendFont.
  ///
  /// In en, this message translates to:
  /// **'Legend text size'**
  String get treeLegendFont;

  /// No description provided for @treeFigureTitle.
  ///
  /// In en, this message translates to:
  /// **'Figure title'**
  String get treeFigureTitle;

  /// No description provided for @treeFigureTitleHint.
  ///
  /// In en, this message translates to:
  /// **'Optional, drawn above the tree'**
  String get treeFigureTitleHint;

  /// No description provided for @treeDpi.
  ///
  /// In en, this message translates to:
  /// **'PNG resolution'**
  String get treeDpi;

  /// No description provided for @treeDpiInfo.
  ///
  /// In en, this message translates to:
  /// **'Dots per inch of the PNG (300 for most journals, 600 for line art). The DPI is written in the file. Very large figures are capped at 16 000 px on the longest side.'**
  String get treeDpiInfo;

  /// No description provided for @treePngBackground.
  ///
  /// In en, this message translates to:
  /// **'PNG background'**
  String get treePngBackground;

  /// No description provided for @treeTransparent.
  ///
  /// In en, this message translates to:
  /// **'Transparent'**
  String get treeTransparent;

  /// No description provided for @treePickColor.
  ///
  /// In en, this message translates to:
  /// **'Choose a colour'**
  String get treePickColor;

  /// No description provided for @concatBlockOutgroup.
  ///
  /// In en, this message translates to:
  /// **'Which taxa root the tree (outgroup)'**
  String get concatBlockOutgroup;

  /// No description provided for @concatBlockOutgroupBody.
  ///
  /// In en, this message translates to:
  /// **'Pick one or more taxa outside the group you study; the trees of steps 7 and 8 are rooted on them. Without an outgroup the tree viewer roots at the midpoint.'**
  String get concatBlockOutgroupBody;

  /// No description provided for @schemeOneModel.
  ///
  /// In en, this message translates to:
  /// **'one model for everything'**
  String get schemeOneModel;

  /// No description provided for @schemePartitionsCount.
  ///
  /// In en, this message translates to:
  /// **'{n} partition(s)'**
  String schemePartitionsCount(String n);

  /// No description provided for @resCpu.
  ///
  /// In en, this message translates to:
  /// **'CPU'**
  String get resCpu;

  /// No description provided for @resCpuInfo.
  ///
  /// In en, this message translates to:
  /// **'CPU load of the machine, including MAFFT, trimAl, IQ-TREE and the other tools SPLACE runs.'**
  String get resCpuInfo;

  /// No description provided for @resMemory.
  ///
  /// In en, this message translates to:
  /// **'RAM'**
  String get resMemory;

  /// No description provided for @resMemoryInfo.
  ///
  /// In en, this message translates to:
  /// **'SPLACE uses {app}. Machine: {used} of {total} in use (the bar).'**
  String resMemoryInfo(String app, String used, String total);

  /// No description provided for @resMemoryInfoApp.
  ///
  /// In en, this message translates to:
  /// **'SPLACE uses {app}.'**
  String resMemoryInfoApp(String app);

  /// No description provided for @resDisk.
  ///
  /// In en, this message translates to:
  /// **'Disk'**
  String get resDisk;

  /// No description provided for @resFree.
  ///
  /// In en, this message translates to:
  /// **'{v} free'**
  String resFree(String v);

  /// No description provided for @resDiskInfo.
  ///
  /// In en, this message translates to:
  /// **'{free} free of {total} on the drive of the project folder (the bar is the space in use).'**
  String resDiskInfo(String free, String total);

  /// No description provided for @resDiskInfoNone.
  ///
  /// In en, this message translates to:
  /// **'Space on the drive of the project folder.'**
  String get resDiskInfoNone;

  /// No description provided for @resNetwork.
  ///
  /// In en, this message translates to:
  /// **'Network'**
  String get resNetwork;

  /// No description provided for @resNetworkInfo.
  ///
  /// In en, this message translates to:
  /// **'Data SPLACE exchanges with NCBI, GBIF and BLAST, per second. Since the app opened: {total}.'**
  String resNetworkInfo(String total);

  /// No description provided for @methodsTitle.
  ///
  /// In en, this message translates to:
  /// **'Methods text'**
  String get methodsTitle;

  /// No description provided for @methodsCopy.
  ///
  /// In en, this message translates to:
  /// **'Copy text and references'**
  String get methodsCopy;

  /// No description provided for @methodsHelp.
  ///
  /// In en, this message translates to:
  /// **'Written from what was done in this project (data, alignment, trimming, saturation, supermatrix and trees). Review it and adapt it to your manuscript.'**
  String get methodsHelp;

  /// No description provided for @citeReferences.
  ///
  /// In en, this message translates to:
  /// **'References cited in the text'**
  String get citeReferences;

  /// No description provided for @citeTools.
  ///
  /// In en, this message translates to:
  /// **'Tools used by SPLACE'**
  String get citeTools;

  /// No description provided for @treeBranchStyle.
  ///
  /// In en, this message translates to:
  /// **'Branch style'**
  String get treeBranchStyle;

  /// No description provided for @treeBranchSquare.
  ///
  /// In en, this message translates to:
  /// **'Right angles'**
  String get treeBranchSquare;

  /// No description provided for @treeBranchRounded.
  ///
  /// In en, this message translates to:
  /// **'Rounded corners'**
  String get treeBranchRounded;

  /// No description provided for @treeBranchSlanted.
  ///
  /// In en, this message translates to:
  /// **'Slanted'**
  String get treeBranchSlanted;

  /// No description provided for @treeZoom.
  ///
  /// In en, this message translates to:
  /// **'Zoom (screen only)'**
  String get treeZoom;

  /// No description provided for @treeZoomFit.
  ///
  /// In en, this message translates to:
  /// **'Fit the window'**
  String get treeZoomFit;

  /// No description provided for @treeZoomInfo.
  ///
  /// In en, this message translates to:
  /// **'Enlarges or reduces the view of the figure on this screen. The exported figure keeps the width set in Figure and export, whatever the monitor.'**
  String get treeZoomInfo;

  /// No description provided for @treeTipDotsLabel.
  ///
  /// In en, this message translates to:
  /// **'Tips'**
  String get treeTipDotsLabel;

  /// No description provided for @treeLabelCustom.
  ///
  /// In en, this message translates to:
  /// **'Custom'**
  String get treeLabelCustom;

  /// No description provided for @treeMarkOutgroup.
  ///
  /// In en, this message translates to:
  /// **'Grey box behind'**
  String get treeMarkOutgroup;

  /// No description provided for @treeRenamed.
  ///
  /// In en, this message translates to:
  /// **'Edited labels'**
  String get treeRenamed;

  /// No description provided for @treeRenamedInfo.
  ///
  /// In en, this message translates to:
  /// **'Click a tip and choose Edit label to write its own text (with italic and bold). Clear the text to go back to the built label.'**
  String get treeRenamedInfo;

  /// No description provided for @treeRenameTip.
  ///
  /// In en, this message translates to:
  /// **'Edit label'**
  String get treeRenameTip;

  /// No description provided for @treeRestoreNames.
  ///
  /// In en, this message translates to:
  /// **'Restore edited labels ({n})'**
  String treeRestoreNames(String n);

  /// No description provided for @treeBarsByRank.
  ///
  /// In en, this message translates to:
  /// **'Bars with taxon names'**
  String get treeBarsByRank;

  /// No description provided for @treeBarsByRankInfo.
  ///
  /// In en, this message translates to:
  /// **'Draws a bar beside each group of tips that share the taxon of this rank (e.g. one bar per family), with its name. Taxa split in the tree get one bar per part.'**
  String get treeBarsByRankInfo;

  /// No description provided for @treeHighlightOpacity.
  ///
  /// In en, this message translates to:
  /// **'Box opacity'**
  String get treeHighlightOpacity;

  /// No description provided for @treeHighlightClade.
  ///
  /// In en, this message translates to:
  /// **'Highlight clade'**
  String get treeHighlightClade;

  /// No description provided for @treeRemoveHighlight.
  ///
  /// In en, this message translates to:
  /// **'Remove highlight'**
  String get treeRemoveHighlight;

  /// No description provided for @treeCladeLabel.
  ///
  /// In en, this message translates to:
  /// **'Name clade'**
  String get treeCladeLabel;

  /// No description provided for @treeClearHighlights.
  ///
  /// In en, this message translates to:
  /// **'Remove highlights ({n})'**
  String treeClearHighlights(String n);

  /// No description provided for @treeClearCladeLabels.
  ///
  /// In en, this message translates to:
  /// **'Remove clade names ({n})'**
  String treeClearCladeLabels(String n);

  /// No description provided for @treeRemoveFromFigure.
  ///
  /// In en, this message translates to:
  /// **'Remove from tree'**
  String get treeRemoveFromFigure;

  /// No description provided for @treeFigureWidth.
  ///
  /// In en, this message translates to:
  /// **'Figure width'**
  String get treeFigureWidth;

  /// No description provided for @treeFigureWidthInfo.
  ///
  /// In en, this message translates to:
  /// **'Real width of the exported figure. Font sizes are in points, so a 9 pt label prints at 9 pt at this width. The height follows from the number of tips and the vertical scale.'**
  String get treeFigureWidthInfo;

  /// No description provided for @treeWidthPreset.
  ///
  /// In en, this message translates to:
  /// **'Journal width'**
  String get treeWidthPreset;

  /// No description provided for @treeWidthOneCol.
  ///
  /// In en, this message translates to:
  /// **'1 column (8.5 cm)'**
  String get treeWidthOneCol;

  /// No description provided for @treeWidthOneHalf.
  ///
  /// In en, this message translates to:
  /// **'1.5 columns (12 cm)'**
  String get treeWidthOneHalf;

  /// No description provided for @treeWidthTwoCol.
  ///
  /// In en, this message translates to:
  /// **'Full page (17.5 cm)'**
  String get treeWidthTwoCol;

  /// No description provided for @treeWidthCustom.
  ///
  /// In en, this message translates to:
  /// **'Custom'**
  String get treeWidthCustom;

  /// No description provided for @treeMargin.
  ///
  /// In en, this message translates to:
  /// **'Margin'**
  String get treeMargin;

  /// No description provided for @treeTitleFont.
  ///
  /// In en, this message translates to:
  /// **'Title size'**
  String get treeTitleFont;

  /// No description provided for @treeFigureSize.
  ///
  /// In en, this message translates to:
  /// **'Figure {w} × {h} cm · PNG {pw} × {ph} px at {dpi} dpi'**
  String treeFigureSize(String w, String h, String pw, String ph, String dpi);

  /// No description provided for @treeTextHint.
  ///
  /// In en, this message translates to:
  /// **'Text'**
  String get treeTextHint;

  /// No description provided for @treeTextEmptyHint.
  ///
  /// In en, this message translates to:
  /// **'Leave it empty to remove it.'**
  String get treeTextEmptyHint;

  /// No description provided for @concatHowOpen.
  ///
  /// In en, this message translates to:
  /// **'Click to see how the loci are put side by side'**
  String get concatHowOpen;

  /// No description provided for @treeOpenEditor.
  ///
  /// In en, this message translates to:
  /// **'Open the tree editor'**
  String get treeOpenEditor;

  /// No description provided for @treeEditorTitle.
  ///
  /// In en, this message translates to:
  /// **'Tree editor'**
  String get treeEditorTitle;

  /// No description provided for @treeEditorHint.
  ///
  /// In en, this message translates to:
  /// **'Preview of the figure at its real size. Click it, or the button above, to edit it in the full window.'**
  String get treeEditorHint;

  /// No description provided for @treeCompare.
  ///
  /// In en, this message translates to:
  /// **'Compare trees'**
  String get treeCompare;

  /// No description provided for @treeModeEdit.
  ///
  /// In en, this message translates to:
  /// **'Edit figure'**
  String get treeModeEdit;

  /// No description provided for @treeDone.
  ///
  /// In en, this message translates to:
  /// **'Done'**
  String get treeDone;

  /// No description provided for @treeExportMenu.
  ///
  /// In en, this message translates to:
  /// **'Export'**
  String get treeExportMenu;

  /// No description provided for @treeZoomIn.
  ///
  /// In en, this message translates to:
  /// **'Zoom in (Ctrl + wheel)'**
  String get treeZoomIn;

  /// No description provided for @treeZoomOut.
  ///
  /// In en, this message translates to:
  /// **'Zoom out (Ctrl + wheel)'**
  String get treeZoomOut;

  /// No description provided for @treeFitWidth.
  ///
  /// In en, this message translates to:
  /// **'Fit width'**
  String get treeFitWidth;

  /// No description provided for @treeFitPage.
  ///
  /// In en, this message translates to:
  /// **'Whole figure'**
  String get treeFitPage;

  /// No description provided for @treeCanvasHint.
  ///
  /// In en, this message translates to:
  /// **'Drag to move · wheel scrolls, Ctrl + wheel zooms · click selects, Shift + click selects the clade joining two tips · right click for options · Del removes · Esc clears'**
  String get treeCanvasHint;

  /// No description provided for @treeTabCompare.
  ///
  /// In en, this message translates to:
  /// **'Compare'**
  String get treeTabCompare;

  /// No description provided for @treeCompareWith.
  ///
  /// In en, this message translates to:
  /// **'Compare with'**
  String get treeCompareWith;

  /// No description provided for @treeUntangle.
  ///
  /// In en, this message translates to:
  /// **'Untangle'**
  String get treeUntangle;

  /// No description provided for @treeUntangleText.
  ///
  /// In en, this message translates to:
  /// **'Rotate the other tree'**
  String get treeUntangleText;

  /// No description provided for @treeUntangleInfo.
  ///
  /// In en, this message translates to:
  /// **'Rotates the clades of the other tree (its topology does not change) so each tip faces its pair in the tree being edited, with as few crossing lines as possible.'**
  String get treeUntangleInfo;

  /// No description provided for @treeConflicts.
  ///
  /// In en, this message translates to:
  /// **'Different clades'**
  String get treeConflicts;

  /// No description provided for @treeConflictsText.
  ///
  /// In en, this message translates to:
  /// **'Mark them'**
  String get treeConflictsText;

  /// No description provided for @treeConflictsInfo.
  ///
  /// In en, this message translates to:
  /// **'Branches leading to a clade that the other tree does not have are drawn thicker, in this colour.'**
  String get treeConflictsInfo;

  /// No description provided for @treeConflictColor.
  ///
  /// In en, this message translates to:
  /// **'Colour of different clades'**
  String get treeConflictColor;

  /// No description provided for @treeConnectors.
  ///
  /// In en, this message translates to:
  /// **'Lines between the trees'**
  String get treeConnectors;

  /// No description provided for @treeConnectorsGroup.
  ///
  /// In en, this message translates to:
  /// **'Row colours'**
  String get treeConnectorsGroup;

  /// No description provided for @treeConnectorsInfo.
  ///
  /// In en, this message translates to:
  /// **'With rows coloured by a rank (Colours tab), each line takes the colour of its taxon; otherwise it is grey.'**
  String get treeConnectorsInfo;

  /// No description provided for @treeCompareGap.
  ///
  /// In en, this message translates to:
  /// **'Space between the trees'**
  String get treeCompareGap;

  /// No description provided for @treeCompareResult.
  ///
  /// In en, this message translates to:
  /// **'Result'**
  String get treeCompareResult;

  /// No description provided for @treeCompareCommon.
  ///
  /// In en, this message translates to:
  /// **'{n} tips in both trees'**
  String treeCompareCommon(String n);

  /// No description provided for @treeCompareOnly.
  ///
  /// In en, this message translates to:
  /// **'{a} only on the left · {b} only on the right'**
  String treeCompareOnly(String a, String b);

  /// No description provided for @treeCompareRf.
  ///
  /// In en, this message translates to:
  /// **'Robinson–Foulds distance {rf} (normalised {nrf})'**
  String treeCompareRf(String rf, String nrf);

  /// No description provided for @treeCompareShared.
  ///
  /// In en, this message translates to:
  /// **'{n} clades in common · {a} only on the left · {b} only on the right'**
  String treeCompareShared(String n, String a, String b);

  /// No description provided for @treeCompareRfInfo.
  ///
  /// In en, this message translates to:
  /// **'Counts the clades found in one tree and not in the other (unrooted, over the tips both trees have). 0 means the same topology; the normalised value divides by the largest possible distance.'**
  String get treeCompareRfInfo;

  /// No description provided for @treeCompareLegend.
  ///
  /// In en, this message translates to:
  /// **'Clade not found in the other tree'**
  String get treeCompareLegend;

  /// No description provided for @treeCompareTooFew.
  ///
  /// In en, this message translates to:
  /// **'The trees share fewer than 4 tips: their clades cannot be compared.'**
  String get treeCompareTooFew;

  /// No description provided for @treeCompareError.
  ///
  /// In en, this message translates to:
  /// **'The other tree could not be read: {e}'**
  String treeCompareError(String e);

  /// No description provided for @treeCompareIntro.
  ///
  /// In en, this message translates to:
  /// **'The tree being edited and the one chosen here are drawn face to face, one mirrored, with a line joining each tip to its pair. Styles (fonts, branches, support, colours) come from the other tabs.'**
  String get treeCompareIntro;

  /// No description provided for @treeSelectParent.
  ///
  /// In en, this message translates to:
  /// **'Parent clade'**
  String get treeSelectParent;

  /// No description provided for @treeLeft.
  ///
  /// In en, this message translates to:
  /// **'Left'**
  String get treeLeft;

  /// No description provided for @treeDetailMoose.
  ///
  /// In en, this message translates to:
  /// **'MOOSE is choosing the models ({n} partition(s))'**
  String treeDetailMoose(String n);

  /// No description provided for @treeKindRbs.
  ///
  /// In en, this message translates to:
  /// **'ML tree, rapid bootstrap (RBS)'**
  String get treeKindRbs;

  /// No description provided for @treeKindSh.
  ///
  /// In en, this message translates to:
  /// **'ML tree, SH-aLRT'**
  String get treeKindSh;

  /// No description provided for @treeKindEbg.
  ///
  /// In en, this message translates to:
  /// **'ML tree, predicted bootstrap (EBG)'**
  String get treeKindEbg;

  /// No description provided for @raxmlOld.
  ///
  /// In en, this message translates to:
  /// **'RAxML-NG {version} found. MOOSE, the fast search and the SH-aLRT, EBG and RBS supports need version 2.0 or newer (SPLACE ships 2.0.4).'**
  String raxmlOld(String version);

  /// No description provided for @raxmlModelsHeading.
  ///
  /// In en, this message translates to:
  /// **'Substitution models'**
  String get raxmlModelsHeading;

  /// No description provided for @raxmlMoose.
  ///
  /// In en, this message translates to:
  /// **'MOOSE (RAxML-NG\'s own model selection)'**
  String get raxmlMoose;

  /// No description provided for @raxmlMooseBody.
  ///
  /// In en, this message translates to:
  /// **'One model per partition, chosen inside the same run. Recommended.'**
  String get raxmlMooseBody;

  /// No description provided for @raxmlModelFinderBody.
  ///
  /// In en, this message translates to:
  /// **'IQ-TREE picks the models and can merge partitions; RAxML-NG then builds the tree.'**
  String get raxmlModelFinderBody;

  /// No description provided for @raxmlFixedModel.
  ///
  /// In en, this message translates to:
  /// **'One fixed model'**
  String get raxmlFixedModel;

  /// No description provided for @raxmlFixedBody.
  ///
  /// In en, this message translates to:
  /// **'The same model for every partition (e.g. GTR+G).'**
  String get raxmlFixedBody;

  /// No description provided for @raxmlCriterion.
  ///
  /// In en, this message translates to:
  /// **'Criterion'**
  String get raxmlCriterion;

  /// No description provided for @raxmlSearch.
  ///
  /// In en, this message translates to:
  /// **'Tree search'**
  String get raxmlSearch;

  /// No description provided for @raxmlSearchAdaptive.
  ///
  /// In en, this message translates to:
  /// **'Adaptive (default)'**
  String get raxmlSearchAdaptive;

  /// No description provided for @raxmlSearchFast.
  ///
  /// In en, this message translates to:
  /// **'Fast (--fast)'**
  String get raxmlSearchFast;

  /// No description provided for @raxmlSupportInfo.
  ///
  /// In en, this message translates to:
  /// **'FBP, TBE and RBS come from bootstrap replicate trees. SH-aLRT and EBG are computed on the ML tree only and take seconds. FBP and RBS cannot be used together.'**
  String get raxmlSupportInfo;

  /// No description provided for @treeLayoutLabel.
  ///
  /// In en, this message translates to:
  /// **'Layout'**
  String get treeLayoutLabel;

  /// No description provided for @treeLayoutRect.
  ///
  /// In en, this message translates to:
  /// **'Rectangular'**
  String get treeLayoutRect;

  /// No description provided for @treeLayoutCircle.
  ///
  /// In en, this message translates to:
  /// **'Circular'**
  String get treeLayoutCircle;

  /// No description provided for @treeLayoutInfo.
  ///
  /// In en, this message translates to:
  /// **'Circular: the root in the centre and the tips around it, good for trees with many tips. Names follow the radius; taxon colours become a ring outside the names and clade bars become arcs.'**
  String get treeLayoutInfo;

  /// No description provided for @treeArc.
  ///
  /// In en, this message translates to:
  /// **'Opening'**
  String get treeArc;

  /// No description provided for @treeArcInfo.
  ///
  /// In en, this message translates to:
  /// **'How much of the circle the tips take: 360° closes it, 180° makes a half circle (fan).'**
  String get treeArcInfo;

  /// No description provided for @treeRotation.
  ///
  /// In en, this message translates to:
  /// **'Rotation'**
  String get treeRotation;

  /// No description provided for @treeSwapSides.
  ///
  /// In en, this message translates to:
  /// **'Swap sides'**
  String get treeSwapSides;

  /// No description provided for @treeSwapSidesText.
  ///
  /// In en, this message translates to:
  /// **'This tree on the right'**
  String get treeSwapSidesText;

  /// No description provided for @treeSwapSidesInfo.
  ///
  /// In en, this message translates to:
  /// **'Puts the tree being edited on the right and the other one on the left.'**
  String get treeSwapSidesInfo;

  /// No description provided for @iqAlrt.
  ///
  /// In en, this message translates to:
  /// **'SH-aLRT test'**
  String get iqAlrt;

  /// No description provided for @iqAlrtInfo.
  ///
  /// In en, this message translates to:
  /// **'A second, fast support value for each branch (shown as SH-aLRT/UFBoot). Branches with SH-aLRT ≥ 80 and UFBoot ≥ 95 are well supported.'**
  String get iqAlrtInfo;

  /// No description provided for @iqAlrt1000.
  ///
  /// In en, this message translates to:
  /// **'1000 replicates (recommended)'**
  String get iqAlrt1000;

  /// No description provided for @iqSearchHeading.
  ///
  /// In en, this message translates to:
  /// **'Tree search'**
  String get iqSearchHeading;

  /// No description provided for @iqSearchInfo.
  ///
  /// In en, this message translates to:
  /// **'The defaults suit most data. More runs or more iterations take longer but explore more trees, useful when results change between runs.'**
  String get iqSearchInfo;

  /// No description provided for @iqRuns.
  ///
  /// In en, this message translates to:
  /// **'Independent runs'**
  String get iqRuns;

  /// No description provided for @iqRunsInfo.
  ///
  /// In en, this message translates to:
  /// **'--runs: repeats the whole search from different starting trees and keeps the best tree.'**
  String get iqRunsInfo;

  /// No description provided for @iqNstop.
  ///
  /// In en, this message translates to:
  /// **'Stop after (iterations)'**
  String get iqNstop;

  /// No description provided for @iqNstopInfo.
  ///
  /// In en, this message translates to:
  /// **'-nstop: the search stops after this many iterations without a better tree (default 100). Use more for large or difficult data.'**
  String get iqNstopInfo;

  /// No description provided for @iqPers.
  ///
  /// In en, this message translates to:
  /// **'Perturbation strength'**
  String get iqPers;

  /// No description provided for @iqPersInfo.
  ///
  /// In en, this message translates to:
  /// **'-pers: how much the current tree is changed at each iteration (default 0.5). Smaller values help with short alignments.'**
  String get iqPersInfo;

  /// No description provided for @iqSeedInfo.
  ///
  /// In en, this message translates to:
  /// **'The same seed gives the same result again. Empty: a random seed (written in the log).'**
  String get iqSeedInfo;

  /// No description provided for @iqSeedRandom.
  ///
  /// In en, this message translates to:
  /// **'random'**
  String get iqSeedRandom;

  /// No description provided for @iqBnni.
  ///
  /// In en, this message translates to:
  /// **'Refine the UFBoot trees (--bnni)'**
  String get iqBnni;

  /// No description provided for @iqBnniBody.
  ///
  /// In en, this message translates to:
  /// **'Optimises each bootstrap tree with NNI moves. Recommended: it reduces overestimated support when the model does not fit the data perfectly.'**
  String get iqBnniBody;

  /// No description provided for @bothTitle.
  ///
  /// In en, this message translates to:
  /// **'Run IQ-TREE and RAxML-NG together'**
  String get bothTitle;

  /// No description provided for @bothInfo.
  ///
  /// In en, this message translates to:
  /// **'Both are maximum-likelihood methods. Comparing their trees (step 8, Compare trees) shows which clades do not depend on the program.'**
  String get bothInfo;

  /// No description provided for @bothBody.
  ///
  /// In en, this message translates to:
  /// **'Both analyses start at the same time with the options chosen in the two cards above. The threads are split so the computer stays usable: IQ-TREE gets {iq} and RAxML-NG {rx}. Change the total in the settings (top bar).'**
  String bothBody(String iq, String rx);

  /// No description provided for @bothThreads.
  ///
  /// In en, this message translates to:
  /// **'{n} threads'**
  String bothThreads(String n);

  /// No description provided for @runBoth.
  ///
  /// In en, this message translates to:
  /// **'Run both'**
  String get runBoth;

  /// No description provided for @bothRunning.
  ///
  /// In en, this message translates to:
  /// **'IQ-TREE and RAxML-NG'**
  String get bothRunning;

  /// No description provided for @bothRunningSub.
  ///
  /// In en, this message translates to:
  /// **'Side by side · threads {threads} · {time}'**
  String bothRunningSub(String threads, String time);

  /// No description provided for @satUnlocksStep7.
  ///
  /// In en, this message translates to:
  /// **'Run the test to continue: its suggestion (which codon positions to keep) is used to build the matrix in step 6, which opens afterwards.'**
  String get satUnlocksStep7;

  /// No description provided for @satAllGenes.
  ///
  /// In en, this message translates to:
  /// **'All genes together'**
  String get satAllGenes;

  /// No description provided for @satAllGenesInfo.
  ///
  /// In en, this message translates to:
  /// **'\"All genes together\" joins the coding genes only for this test (the matrix itself is built in step 6). It is the row that decides the suggestion of which codon positions to remove from the whole matrix; each gene\'s row allows a choice per gene.'**
  String get satAllGenesInfo;

  /// No description provided for @exampleLabel.
  ///
  /// In en, this message translates to:
  /// **'Example'**
  String get exampleLabel;

  /// No description provided for @exampleNt.
  ///
  /// In en, this message translates to:
  /// **'The same alignment in the three modes. Here the gap of species B falls inside codons 3 and 4 and shifts the reading frame.'**
  String get exampleNt;

  /// No description provided for @exampleAa.
  ///
  /// In en, this message translates to:
  /// **'The same genes translated: one letter per codon, so fewer and more conserved sites.'**
  String get exampleAa;

  /// No description provided for @exampleCodon.
  ///
  /// In en, this message translates to:
  /// **'The gap of species B is the whole codon 3: the reading frame is kept.'**
  String get exampleCodon;

  /// No description provided for @exampleTrim.
  ///
  /// In en, this message translates to:
  /// **'Trimming removes unreliable columns (mostly gaps, or sites without phylogenetic signal) before the tree is built. Here {n} of {total} columns are removed (×).'**
  String exampleTrim(String n, String total);

  /// No description provided for @exampleCodonKeepAll.
  ///
  /// In en, this message translates to:
  /// **'All three codon positions are kept.'**
  String get exampleCodonKeepAll;

  /// No description provided for @exampleCodonKept.
  ///
  /// In en, this message translates to:
  /// **'Kept: {positions} positions ({pct}% of the coding sites).'**
  String exampleCodonKept(String positions, String pct);

  /// No description provided for @exampleSatLabel.
  ///
  /// In en, this message translates to:
  /// **'How to read the plots'**
  String get exampleSatLabel;

  /// No description provided for @exampleSatOk.
  ///
  /// In en, this message translates to:
  /// **'Not saturated'**
  String get exampleSatOk;

  /// No description provided for @exampleSatOkBody.
  ///
  /// In en, this message translates to:
  /// **'Observed differences (p) keep rising with distance (d), close to the dashed y = x line.'**
  String get exampleSatOkBody;

  /// No description provided for @exampleSatBad.
  ///
  /// In en, this message translates to:
  /// **'Saturated'**
  String get exampleSatBad;

  /// No description provided for @exampleSatBadBody.
  ///
  /// In en, this message translates to:
  /// **'p levels off while d keeps growing: new changes overwrite old ones and the signal is lost.'**
  String get exampleSatBadBody;

  /// No description provided for @plotTrend.
  ///
  /// In en, this message translates to:
  /// **'trend'**
  String get plotTrend;

  /// No description provided for @exampleSpecies.
  ///
  /// In en, this message translates to:
  /// **'Species {letter}'**
  String exampleSpecies(String letter);

  /// No description provided for @examplePlayAgain.
  ///
  /// In en, this message translates to:
  /// **'Play again'**
  String get examplePlayAgain;

  /// No description provided for @trimStepOriginal.
  ///
  /// In en, this message translates to:
  /// **'Alignment'**
  String get trimStepOriginal;

  /// No description provided for @trimStepMarked.
  ///
  /// In en, this message translates to:
  /// **'Columns removed'**
  String get trimStepMarked;

  /// No description provided for @trimStepResult.
  ///
  /// In en, this message translates to:
  /// **'Trimmed alignment'**
  String get trimStepResult;

  /// No description provided for @treeSupportArrowLabel.
  ///
  /// In en, this message translates to:
  /// **'Arrows'**
  String get treeSupportArrowLabel;

  /// No description provided for @treeSupportArrow.
  ///
  /// In en, this message translates to:
  /// **'Value with an arrow to the branch'**
  String get treeSupportArrow;

  /// No description provided for @treeSupportArrowInfo.
  ///
  /// In en, this message translates to:
  /// **'Moves each support value a little away from the branches (up or down, as in Support position) and joins it to its node with a small arrow. Useful when values fall on short branches or on top of other lines.'**
  String get treeSupportArrowInfo;

  /// No description provided for @orgOpen.
  ///
  /// In en, this message translates to:
  /// **'Sample and marker'**
  String get orgOpen;

  /// No description provided for @orgTitle.
  ///
  /// In en, this message translates to:
  /// **'Sample and marker in the headers'**
  String get orgTitle;

  /// No description provided for @orgSubtitle.
  ///
  /// In en, this message translates to:
  /// **'{files} FASTA file(s) · {seqs} sequences'**
  String orgSubtitle(String files, String seqs);

  /// No description provided for @orgSummary.
  ///
  /// In en, this message translates to:
  /// **'{samples} samples × {markers} markers · {pct}% filled'**
  String orgSummary(String samples, String markers, String pct);

  /// No description provided for @orgApply.
  ///
  /// In en, this message translates to:
  /// **'Create {n} marker FASTA files'**
  String orgApply(String n);

  /// No description provided for @orgIntro.
  ///
  /// In en, this message translates to:
  /// **'Each header must say which sample (ID) and which marker the sequence belongs to. Choose the characters that separate the parts of the header, then click the parts below: once for the sample, twice for the marker, a third time to ignore it. SPLACE then writes one FASTA per marker holding every sample, ready to align.'**
  String get orgIntro;

  /// No description provided for @orgSampleFrom.
  ///
  /// In en, this message translates to:
  /// **'Sample (ID) from'**
  String get orgSampleFrom;

  /// No description provided for @orgSampleFromInfo.
  ///
  /// In en, this message translates to:
  /// **'Part of the header, or the file name when each file holds one sample and the headers only name the marker.'**
  String get orgSampleFromInfo;

  /// No description provided for @orgSampleHeader.
  ///
  /// In en, this message translates to:
  /// **'Header'**
  String get orgSampleHeader;

  /// No description provided for @orgSampleFile.
  ///
  /// In en, this message translates to:
  /// **'File name'**
  String get orgSampleFile;

  /// No description provided for @orgDelimitersInfo.
  ///
  /// In en, this message translates to:
  /// **'Every selected character splits the header. Several can be on at once when files use different separators (ITV01-L1 and ITV02_L1).'**
  String get orgDelimitersInfo;

  /// No description provided for @orgSpace.
  ///
  /// In en, this message translates to:
  /// **'space'**
  String get orgSpace;

  /// No description provided for @orgOther.
  ///
  /// In en, this message translates to:
  /// **'others'**
  String get orgOther;

  /// No description provided for @orgParts.
  ///
  /// In en, this message translates to:
  /// **'Parts of the header (examples from different files)'**
  String get orgParts;

  /// No description provided for @orgPartsInfo.
  ///
  /// In en, this message translates to:
  /// **'Positions are kept for every header. A part that reaches the end of the header is counted from the end, so headers with more or fewer parts still work.'**
  String get orgPartsInfo;

  /// No description provided for @orgLegendSample.
  ///
  /// In en, this message translates to:
  /// **'Sample (ID)'**
  String get orgLegendSample;

  /// No description provided for @orgLegendMarker.
  ///
  /// In en, this message translates to:
  /// **'Marker'**
  String get orgLegendMarker;

  /// No description provided for @orgLegendIgnored.
  ///
  /// In en, this message translates to:
  /// **'Ignored'**
  String get orgLegendIgnored;

  /// No description provided for @orgAdvanced.
  ///
  /// In en, this message translates to:
  /// **'Advanced: tags and regular expressions'**
  String get orgAdvanced;

  /// No description provided for @orgAdvancedActive.
  ///
  /// In en, this message translates to:
  /// **'A tag or regular expression is in use (see Advanced); clicking the parts is off.'**
  String get orgAdvancedActive;

  /// No description provided for @orgImplausible.
  ///
  /// In en, this message translates to:
  /// **'{markers} markers for {seqs} sequences'**
  String orgImplausible(String markers, String seqs);

  /// No description provided for @orgImplausibleBody.
  ///
  /// In en, this message translates to:
  /// **'Almost every sequence became its own marker, so the marker part is probably wrong (for example the whole header). Pick the delimiter and click the part that names the marker.'**
  String get orgImplausibleBody;

  /// No description provided for @orgUnassigned.
  ///
  /// In en, this message translates to:
  /// **'{n} sequence(s) without a sample or a marker (left out)'**
  String orgUnassigned(String n);

  /// No description provided for @orgUnassignedLabel.
  ///
  /// In en, this message translates to:
  /// **'Without sample or marker'**
  String get orgUnassignedLabel;

  /// No description provided for @orgVariations.
  ///
  /// In en, this message translates to:
  /// **'Name variations'**
  String get orgVariations;

  /// No description provided for @orgVariationsInfo.
  ///
  /// In en, this message translates to:
  /// **'The same sample or marker can be written in slightly different ways. These options join them; you can also rename any sample in the list below.'**
  String get orgVariationsInfo;

  /// No description provided for @orgIgnoreCase.
  ///
  /// In en, this message translates to:
  /// **'Ignore upper and lower case'**
  String get orgIgnoreCase;

  /// No description provided for @orgLeadingZeros.
  ///
  /// In en, this message translates to:
  /// **'Ignore zeros at the start of numbers'**
  String get orgLeadingZeros;

  /// No description provided for @orgFileMismatch.
  ///
  /// In en, this message translates to:
  /// **'{n} file(s) whose name differs from the sample in its headers'**
  String orgFileMismatch(String n);

  /// No description provided for @orgRenameTo.
  ///
  /// In en, this message translates to:
  /// **'Call it {name}'**
  String orgRenameTo(String name);

  /// No description provided for @orgRenameAll.
  ///
  /// In en, this message translates to:
  /// **'Use the file names for these'**
  String get orgRenameAll;

  /// No description provided for @orgUseFileNames.
  ///
  /// In en, this message translates to:
  /// **'Take every sample from the file name'**
  String get orgUseFileNames;

  /// No description provided for @orgMerged.
  ///
  /// In en, this message translates to:
  /// **'{n} sample(s) written in more than one way, joined:'**
  String orgMerged(String n);

  /// No description provided for @orgSamples.
  ///
  /// In en, this message translates to:
  /// **'Samples ({n})'**
  String orgSamples(String n);

  /// No description provided for @orgSamplesInfo.
  ///
  /// In en, this message translates to:
  /// **'Click a name to rename it (Enter to apply). Giving two samples the same name joins them.'**
  String get orgSamplesInfo;

  /// No description provided for @orgFindSample.
  ///
  /// In en, this message translates to:
  /// **'Find a sample'**
  String get orgFindSample;

  /// No description provided for @orgSampleMarkers.
  ///
  /// In en, this message translates to:
  /// **'{n} of {total} markers'**
  String orgSampleMarkers(String n, String total);

  /// No description provided for @orgMarkers.
  ///
  /// In en, this message translates to:
  /// **'Markers ({n})'**
  String orgMarkers(String n);

  /// No description provided for @orgMarkersInfo.
  ///
  /// In en, this message translates to:
  /// **'How many samples each marker has, and what to do when a sample has more than one sequence of the same marker. Every marker goes to step 3, where you choose which to use.'**
  String get orgMarkersInfo;

  /// No description provided for @orgDuplicates.
  ///
  /// In en, this message translates to:
  /// **'{n} sample(s) with more than one sequence of a marker'**
  String orgDuplicates(String n);

  /// No description provided for @orgCoverage.
  ///
  /// In en, this message translates to:
  /// **'Markers by number of samples'**
  String get orgCoverage;

  /// No description provided for @orgCoverageBar.
  ///
  /// In en, this message translates to:
  /// **'{markers} markers in {samples} samples'**
  String orgCoverageBar(String markers, String samples);

  /// No description provided for @orgSamplesAxis.
  ///
  /// In en, this message translates to:
  /// **'samples'**
  String get orgSamplesAxis;

  /// No description provided for @orgCoverageSummary.
  ///
  /// In en, this message translates to:
  /// **'{complete} in every sample · {half} in at least half of them'**
  String orgCoverageSummary(String complete, String half);

  /// No description provided for @orgPending.
  ///
  /// In en, this message translates to:
  /// **'Tell SPLACE which part of the headers is the sample and which the marker'**
  String get orgPending;

  /// No description provided for @orgProblem.
  ///
  /// In en, this message translates to:
  /// **'The split is not used: {problem}'**
  String orgProblem(String problem);

  /// No description provided for @orgPendingBody.
  ///
  /// In en, this message translates to:
  /// **'Then one FASTA per marker is written with every sample, and the markers appear in step 3.'**
  String get orgPendingBody;

  /// No description provided for @orgMatrix.
  ///
  /// In en, this message translates to:
  /// **'Matrix'**
  String get orgMatrix;

  /// No description provided for @orgWritten.
  ///
  /// In en, this message translates to:
  /// **'FASTA per marker'**
  String get orgWritten;

  /// No description provided for @findMarker.
  ///
  /// In en, this message translates to:
  /// **'Find a marker'**
  String get findMarker;

  /// No description provided for @resultsShowing.
  ///
  /// In en, this message translates to:
  /// **'Showing {n} of {total}'**
  String resultsShowing(String n, String total);

  /// No description provided for @showMoreN.
  ///
  /// In en, this message translates to:
  /// **'Show {n} more'**
  String showMoreN(String n);

  /// No description provided for @showAll.
  ///
  /// In en, this message translates to:
  /// **'Show all'**
  String get showAll;

  /// No description provided for @noMatches.
  ///
  /// In en, this message translates to:
  /// **'Nothing matches.'**
  String get noMatches;

  /// No description provided for @orgFormats.
  ///
  /// In en, this message translates to:
  /// **'{n} header formats'**
  String orgFormats(String n);

  /// No description provided for @perLocusCount.
  ///
  /// In en, this message translates to:
  /// **'{positions} in {n} markers'**
  String perLocusCount(String positions, String n);

  /// No description provided for @satOptionalLoci.
  ///
  /// In en, this message translates to:
  /// **'Optional for loci: step 6 is already open. The markers marked as protein-coding are tested by codon position; the others on all their sites.'**
  String get satOptionalLoci;

  /// No description provided for @satOptionalLociNoCoding.
  ///
  /// In en, this message translates to:
  /// **'Optional for loci: step 6 is already open. No marker is protein-coding, so each locus is tested on all its sites (there are no codon positions): the test shows which loci are saturated.'**
  String get satOptionalLociNoCoding;

  /// No description provided for @exampleTrimNone.
  ///
  /// In en, this message translates to:
  /// **'None of the {total} columns is removed: a MAFFT alignment has no column made only of gaps. This option only matters for alignments edited by hand or subsets of sequences.'**
  String exampleTrimNone(String total);

  /// No description provided for @exampleGapsRow.
  ///
  /// In en, this message translates to:
  /// **'gaps'**
  String get exampleGapsRow;

  /// No description provided for @runFailed.
  ///
  /// In en, this message translates to:
  /// **'The analysis failed.'**
  String get runFailed;

  /// No description provided for @runLastLines.
  ///
  /// In en, this message translates to:
  /// **'Last lines of the program'**
  String get runLastLines;

  /// No description provided for @runOpenLog.
  ///
  /// In en, this message translates to:
  /// **'Open the log'**
  String get runOpenLog;

  /// No description provided for @orgHeaderFormats.
  ///
  /// In en, this message translates to:
  /// **'Header formats ({n})'**
  String orgHeaderFormats(String n);

  /// No description provided for @orgHeaderFormatsInfo.
  ///
  /// In en, this message translates to:
  /// **'The headers are written in more than one way (other delimiters, more or fewer parts). Each format shows its sequences and can have its own sample and marker parts: click the parts of its example. Every sequence then ends up named the same way: the sample in the header, one FASTA file per marker.'**
  String get orgHeaderFormatsInfo;

  /// No description provided for @orgFormatCount.
  ///
  /// In en, this message translates to:
  /// **'{n} sequences in {files} files'**
  String orgFormatCount(String n, String files);

  /// No description provided for @orgFormatOwnRule.
  ///
  /// In en, this message translates to:
  /// **'Own parts'**
  String get orgFormatOwnRule;

  /// No description provided for @orgFormatSkipped.
  ///
  /// In en, this message translates to:
  /// **'Left out'**
  String get orgFormatSkipped;

  /// No description provided for @orgFormatShow.
  ///
  /// In en, this message translates to:
  /// **'Show the {n} sequences'**
  String orgFormatShow(String n);

  /// No description provided for @orgFormatHide.
  ///
  /// In en, this message translates to:
  /// **'Hide the sequences'**
  String get orgFormatHide;

  /// No description provided for @orgFormatUseMain.
  ///
  /// In en, this message translates to:
  /// **'Same parts as the main rule'**
  String get orgFormatUseMain;

  /// No description provided for @orgFormatSkip.
  ///
  /// In en, this message translates to:
  /// **'Leave these sequences out'**
  String get orgFormatSkip;

  /// No description provided for @orgFormatUse.
  ///
  /// In en, this message translates to:
  /// **'Use these sequences'**
  String get orgFormatUse;

  /// No description provided for @orgNameStructures.
  ///
  /// In en, this message translates to:
  /// **'Sample name structures ({n})'**
  String orgNameStructures(String n);

  /// No description provided for @orgNameStructuresInfo.
  ///
  /// In en, this message translates to:
  /// **'Sample names grouped by how they are written: letters as they are, each digit as #. Names of one collection written with more or fewer digits (ITV1054 and ITV01054) can be standardized to the structure you choose: zeros are added or removed on the left. Hover the button to see the renames.'**
  String get orgNameStructuresInfo;

  /// No description provided for @orgStructureCount.
  ///
  /// In en, this message translates to:
  /// **'{n, plural, =1{1 sample} other{{n} samples}}'**
  String orgStructureCount(int n);

  /// No description provided for @orgStandardize.
  ///
  /// In en, this message translates to:
  /// **'Standardize to {structure}'**
  String orgStandardize(String structure);

  /// No description provided for @satAllSitesNonCoding.
  ///
  /// In en, this message translates to:
  /// **'all sites (non-coding)'**
  String get satAllSitesNonCoding;

  /// No description provided for @schemeMixed.
  ///
  /// In en, this message translates to:
  /// **'{coding} coding markers split by codon; {other} without codons stay whole'**
  String schemeMixed(String coding, String other);

  /// No description provided for @exampleCodonGene.
  ///
  /// In en, this message translates to:
  /// **'Gene'**
  String get exampleCodonGene;

  /// No description provided for @exampleCodonMatrix.
  ///
  /// In en, this message translates to:
  /// **'In the matrix'**
  String get exampleCodonMatrix;

  /// No description provided for @mtTitle.
  ///
  /// In en, this message translates to:
  /// **'Taxa × markers'**
  String get mtTitle;

  /// No description provided for @mtIntro.
  ///
  /// In en, this message translates to:
  /// **'As in SequenceMatrix: each cell is the number of bases the taxon has in that marker (darker green = more complete), \"—\" when it has no sequence. The header gives each alignment\'s length, how many taxa have it and its gaps. Hover a cell or a header for details; the same table is saved as matrix_overview.tsv.'**
  String get mtIntro;

  /// No description provided for @mtTaxon.
  ///
  /// In en, this message translates to:
  /// **'Taxon'**
  String get mtTaxon;

  /// No description provided for @mtLoci.
  ///
  /// In en, this message translates to:
  /// **'Markers'**
  String get mtLoci;

  /// No description provided for @mtBases.
  ///
  /// In en, this message translates to:
  /// **'Bases'**
  String get mtBases;

  /// No description provided for @mtMissing.
  ///
  /// In en, this message translates to:
  /// **'Missing'**
  String get mtMissing;

  /// No description provided for @mtPartitions.
  ///
  /// In en, this message translates to:
  /// **'Partitions'**
  String get mtPartitions;

  /// No description provided for @mtLegend.
  ///
  /// In en, this message translates to:
  /// **'Cell: bases of the taxon in the marker · — no sequence · header: length · taxa with the marker · % gaps'**
  String get mtLegend;

  /// No description provided for @mtLocusInfo.
  ///
  /// In en, this message translates to:
  /// **'{locus}: {length} bp, in {n} of {total} taxa, {gaps}% gaps'**
  String mtLocusInfo(
    String locus,
    String length,
    String n,
    String total,
    String gaps,
  );

  /// No description provided for @mtMissingCell.
  ///
  /// In en, this message translates to:
  /// **'no sequence (filled with ?)'**
  String get mtMissingCell;

  /// No description provided for @mtCellInfo.
  ///
  /// In en, this message translates to:
  /// **'{bases} of {length} bp, {gaps} gaps'**
  String mtCellInfo(String bases, String length, String gaps);

  /// No description provided for @searchGenomes.
  ///
  /// In en, this message translates to:
  /// **'Genomes'**
  String get searchGenomes;

  /// No description provided for @searchMarkers.
  ///
  /// In en, this message translates to:
  /// **'Individual markers'**
  String get searchMarkers;

  /// No description provided for @markerIncludeGenomes.
  ///
  /// In en, this message translates to:
  /// **'Also genomes'**
  String get markerIncludeGenomes;

  /// No description provided for @markerGenomesInfo.
  ///
  /// In en, this message translates to:
  /// **'Off: only records up to 10 kb (gene sequences such as barcodes). On: complete genomes and long fragments too, which also hold the markers (and are much larger to download).'**
  String get markerGenomesInfo;

  /// No description provided for @markerPick.
  ///
  /// In en, this message translates to:
  /// **'Markers (SynGenes) · {n} selected'**
  String markerPick(String n);

  /// No description provided for @markerPickInfo.
  ///
  /// In en, this message translates to:
  /// **'Protein-coding genes and rRNAs of the SynGenes tables. Each marker is searched with every name SynGenes knows for it (gene symbols and product names, e.g. COI, COX1, cytochrome c oxidase subunit I), so records annotated in different ways are found. Records of the same species are joined into one sample, and species missing a marker get it as missing data.'**
  String get markerPickInfo;

  /// No description provided for @reviewMarkerCount.
  ///
  /// In en, this message translates to:
  /// **'{marker}: {selected}/{found} records · {species} species'**
  String reviewMarkerCount(
    String marker,
    String selected,
    String found,
    String species,
  );

  /// No description provided for @reviewOnePerSpecies.
  ///
  /// In en, this message translates to:
  /// **'One per species and marker'**
  String get reviewOnePerSpecies;

  /// No description provided for @reviewOnePerSpeciesInfo.
  ///
  /// In en, this message translates to:
  /// **'Keeps, for each marker, the longest record of each species (verified records first). It is the selection the review starts with.'**
  String get reviewOnePerSpeciesInfo;

  /// No description provided for @reviewColMarkers.
  ///
  /// In en, this message translates to:
  /// **'Markers'**
  String get reviewColMarkers;

  /// No description provided for @onePerTaxon.
  ///
  /// In en, this message translates to:
  /// **'One sample per species (individual markers)'**
  String get onePerTaxon;

  /// No description provided for @onePerTaxonBody.
  ///
  /// In en, this message translates to:
  /// **'Joins records of the same species (a COI from one voucher, a CYTB from another) into one row of the matrix: the header is family_species and, for each marker, the longest sequence of the species is kept. Species without a marker get it as missing data (?).'**
  String get onePerTaxonBody;

  /// No description provided for @onePerTaxonJoined.
  ///
  /// In en, this message translates to:
  /// **'{n} species have several records: they are joined into one sample each.'**
  String onePerTaxonJoined(String n);

  /// No description provided for @satOptionalPartial.
  ///
  /// In en, this message translates to:
  /// **'Optional here: step 6 is already open. {markers} came from individual-marker records with partial sequences (fragments such as barcodes, marked < or > in GenBank), so the test is a guide, not a requirement.'**
  String satOptionalPartial(String markers);

  /// No description provided for @satPartialMixed.
  ///
  /// In en, this message translates to:
  /// **'{markers} have partial sequences: the test is optional for them. Step 6 waits for the complete markers.'**
  String satPartialMixed(String markers);

  /// No description provided for @metaRecords.
  ///
  /// In en, this message translates to:
  /// **'{n} records'**
  String metaRecords(String n);

  /// No description provided for @headerConflictBlock.
  ///
  /// In en, this message translates to:
  /// **'{n} header(s) repeat inside the same marker, so those sequences could not be told apart in the alignment or the tree. Add a field that tells them apart (Accession or Voucher), or turn on \"One sample per species\" to join them. The alignment and every step after it (saturation, matrix, trees, delimitation) stay locked until every header is unique.'**
  String headerConflictBlock(String n);

  /// No description provided for @headerAddAccession.
  ///
  /// In en, this message translates to:
  /// **'Add Accession'**
  String get headerAddAccession;

  /// No description provided for @alignBlockedHeader.
  ///
  /// In en, this message translates to:
  /// **'Locked: repeated FASTA headers (4a)'**
  String get alignBlockedHeader;

  /// No description provided for @recGenesCount.
  ///
  /// In en, this message translates to:
  /// **'Saturated in {n} of {total} genes.'**
  String recGenesCount(String n, String total);

  /// No description provided for @recBasedOnGene.
  ///
  /// In en, this message translates to:
  /// **'Result of {gene} ({taxa} taxa, {length} bp). Each position is tested separately.'**
  String recBasedOnGene(String gene, String taxa, String length);

  /// No description provided for @recBasedOnGenes.
  ///
  /// In en, this message translates to:
  /// **'{n} coding genes, each tested on its own (nothing is concatenated for the test). A position counts as saturated when it is saturated in more than half of the genes.'**
  String recBasedOnGenes(String n);

  /// No description provided for @rec7bPreview.
  ///
  /// In en, this message translates to:
  /// **'How \"{section}\" in step 6b will look:'**
  String rec7bPreview(String section);

  /// No description provided for @frameExampleLabel.
  ///
  /// In en, this message translates to:
  /// **'Example · keep the reading frame'**
  String get frameExampleLabel;

  /// No description provided for @frameExampleCaption.
  ///
  /// In en, this message translates to:
  /// **'trimAl removes single columns. With this option, when trimAl removes a column (×), SPLACE removes the whole codon it belongs to (+): the gene keeps whole codons and the 1st, 2nd and 3rd positions stay in order (the coloured ruler).'**
  String get frameExampleCaption;

  /// No description provided for @frameStepOriginal.
  ///
  /// In en, this message translates to:
  /// **'Codons 1-2-3'**
  String get frameStepOriginal;

  /// No description provided for @frameStepColumn.
  ///
  /// In en, this message translates to:
  /// **'trimAl removes a column'**
  String get frameStepColumn;

  /// No description provided for @frameStepCodon.
  ///
  /// In en, this message translates to:
  /// **'ClipKIT --codon takes the whole codon'**
  String get frameStepCodon;

  /// No description provided for @frameStepResult.
  ///
  /// In en, this message translates to:
  /// **'Frame kept'**
  String get frameStepResult;

  /// No description provided for @trimCasesLabel.
  ///
  /// In en, this message translates to:
  /// **'Example · trimming a coding gene, two ways'**
  String get trimCasesLabel;

  /// No description provided for @trimCasesStepMark.
  ///
  /// In en, this message translates to:
  /// **'Columns marked'**
  String get trimCasesStepMark;

  /// No description provided for @trimCaseColumns.
  ///
  /// In en, this message translates to:
  /// **'Column by column (trimAl)'**
  String get trimCaseColumns;

  /// No description provided for @trimCaseCodons.
  ///
  /// In en, this message translates to:
  /// **'Keep the reading frame (ClipKIT --codon)'**
  String get trimCaseCodons;

  /// No description provided for @trimCaseSelected.
  ///
  /// In en, this message translates to:
  /// **'selected'**
  String get trimCaseSelected;

  /// No description provided for @trimCaseFrameLost.
  ///
  /// In en, this message translates to:
  /// **'Frame lost: {n} columns out of codon position'**
  String trimCaseFrameLost(String n);

  /// No description provided for @trimCaseFrameKept.
  ///
  /// In en, this message translates to:
  /// **'Frame kept: {k} of {total} codons, positions 1-2-3 in order'**
  String trimCaseFrameKept(String k, String total);

  /// No description provided for @trimCaseFrameKeptCols.
  ///
  /// In en, this message translates to:
  /// **'Frame kept: no column out of position'**
  String get trimCaseFrameKeptCols;

  /// No description provided for @trimCasesCaption.
  ///
  /// In en, this message translates to:
  /// **'Same alignment, two tools. trimAl removes the columns it marks (×) one by one, so the codon positions (coloured ruler) fall out of order (!). ClipKIT with --codon removes the whole codon of every column it marks (+): the 1st, 2nd and 3rd positions stay in order for the partitions and the saturation test.'**
  String get trimCasesCaption;

  /// No description provided for @subRecordsTable.
  ///
  /// In en, this message translates to:
  /// **'Records table'**
  String get subRecordsTable;

  /// No description provided for @subRecordsTableBody.
  ///
  /// In en, this message translates to:
  /// **'Every record with its taxonomy and gene counts; edit or remove records here.'**
  String get subRecordsTableBody;

  /// No description provided for @subFeatureTypes.
  ///
  /// In en, this message translates to:
  /// **'Annotations to use (CDS, rRNA, tRNA)'**
  String get subFeatureTypes;

  /// No description provided for @subChooseMarkers.
  ///
  /// In en, this message translates to:
  /// **'Choose the markers'**
  String get subChooseMarkers;

  /// No description provided for @subRunSaturation.
  ///
  /// In en, this message translates to:
  /// **'Run the test'**
  String get subRunSaturation;

  /// No description provided for @subBuildMatrix.
  ///
  /// In en, this message translates to:
  /// **'Build the matrix'**
  String get subBuildMatrix;

  /// No description provided for @subTree.
  ///
  /// In en, this message translates to:
  /// **'Tree and figure'**
  String get subTree;

  /// No description provided for @viewCodons.
  ///
  /// In en, this message translates to:
  /// **'Codons'**
  String get viewCodons;

  /// No description provided for @viewDots.
  ///
  /// In en, this message translates to:
  /// **'Dots'**
  String get viewDots;

  /// No description provided for @viewDotsInfo.
  ///
  /// In en, this message translates to:
  /// **'Shows a dot in place of the bases in columns where every sequence has the same base (the first sequence keeps its letters), so the variable columns stand out, as in MEGA and other viewers. In the codon view, codons that are the same in every sequence become dots.'**
  String get viewDotsInfo;

  /// No description provided for @codonsVarying.
  ///
  /// In en, this message translates to:
  /// **'{n} of {total} codons vary · {nonsyn} change the amino acid'**
  String codonsVarying(String n, String total, String nonsyn);

  /// No description provided for @codonVariation.
  ///
  /// In en, this message translates to:
  /// **'Codon variation'**
  String get codonVariation;

  /// No description provided for @viewerCodonN.
  ///
  /// In en, this message translates to:
  /// **'codon {n}'**
  String viewerCodonN(String n);

  /// No description provided for @codonColumnInfo.
  ///
  /// In en, this message translates to:
  /// **'{codons} different codon(s) here, {aa} amino acid(s)'**
  String codonColumnInfo(String codons, String aa);

  /// No description provided for @codonLegendAa.
  ///
  /// In en, this message translates to:
  /// **'colour = amino acid'**
  String get codonLegendAa;

  /// No description provided for @codonLegendSyn.
  ///
  /// In en, this message translates to:
  /// **'synonymous change'**
  String get codonLegendSyn;

  /// No description provided for @codonLegendNonSyn.
  ///
  /// In en, this message translates to:
  /// **'amino acid changes'**
  String get codonLegendNonSyn;

  /// No description provided for @mtResizeNames.
  ///
  /// In en, this message translates to:
  /// **'Drag to widen the names; double-click to fit the longest one'**
  String get mtResizeNames;

  /// No description provided for @treeLegendUnmarked.
  ///
  /// In en, this message translates to:
  /// **'no mark: support < {min}'**
  String treeLegendUnmarked(String min);

  /// No description provided for @treeStrongFrom.
  ///
  /// In en, this message translates to:
  /// **'Strong support from'**
  String get treeStrongFrom;

  /// No description provided for @treeModerateFrom.
  ///
  /// In en, this message translates to:
  /// **'Moderate support from'**
  String get treeModerateFrom;

  /// No description provided for @treeSupportClassesInfo.
  ///
  /// In en, this message translates to:
  /// **'Each class colours the points (and values) from its threshold up to the next one. Pick a preset — automatic for the kind of support, 2 to 5 classes, or cut at the tertiles of this tree\'s own values — or set the thresholds and colours yourself, adding or removing classes.'**
  String get treeSupportClassesInfo;

  /// No description provided for @treeLegendSupport.
  ///
  /// In en, this message translates to:
  /// **'Branch support'**
  String get treeLegendSupport;

  /// No description provided for @treeLegendSupportShow.
  ///
  /// In en, this message translates to:
  /// **'In the legend'**
  String get treeLegendSupportShow;

  /// No description provided for @treeLegendSupportTitle.
  ///
  /// In en, this message translates to:
  /// **'Support heading'**
  String get treeLegendSupportTitle;

  /// No description provided for @badgeBrazil.
  ///
  /// In en, this message translates to:
  /// **'Brazil only'**
  String get badgeBrazil;

  /// No description provided for @salvePossiblyExtinct.
  ///
  /// In en, this message translates to:
  /// **'Possibly extinct'**
  String get salvePossiblyExtinct;

  /// No description provided for @salveSheetSubtitle.
  ///
  /// In en, this message translates to:
  /// **'SALVE sheet (ICMBio) · national extinction-risk assessment'**
  String get salveSheetSubtitle;

  /// No description provided for @salveOpenSheet.
  ///
  /// In en, this message translates to:
  /// **'Open the sheet in SALVE'**
  String get salveOpenSheet;

  /// No description provided for @salveWatersheds.
  ///
  /// In en, this message translates to:
  /// **'Watersheds'**
  String get salveWatersheds;

  /// No description provided for @salvePublished.
  ///
  /// In en, this message translates to:
  /// **'Published'**
  String get salvePublished;

  /// No description provided for @salveSituation.
  ///
  /// In en, this message translates to:
  /// **'Sheet status'**
  String get salveSituation;

  /// No description provided for @salveSynonyms.
  ///
  /// In en, this message translates to:
  /// **'Former names'**
  String get salveSynonyms;

  /// No description provided for @salveJustification.
  ///
  /// In en, this message translates to:
  /// **'Justification'**
  String get salveJustification;

  /// No description provided for @salveAssessors.
  ///
  /// In en, this message translates to:
  /// **'Assessors'**
  String get salveAssessors;

  /// No description provided for @salveCitation.
  ///
  /// In en, this message translates to:
  /// **'How to cite'**
  String get salveCitation;

  /// No description provided for @colCriteria.
  ///
  /// In en, this message translates to:
  /// **'Criteria'**
  String get colCriteria;

  /// No description provided for @colCommonNames.
  ///
  /// In en, this message translates to:
  /// **'Common names'**
  String get colCommonNames;

  /// No description provided for @colSalveGroup.
  ///
  /// In en, this message translates to:
  /// **'SALVE group'**
  String get colSalveGroup;

  /// No description provided for @colEndemic.
  ///
  /// In en, this message translates to:
  /// **'Endemic to Brazil'**
  String get colEndemic;

  /// No description provided for @colBiomes.
  ///
  /// In en, this message translates to:
  /// **'Biomes'**
  String get colBiomes;

  /// No description provided for @colStates.
  ///
  /// In en, this message translates to:
  /// **'States'**
  String get colStates;

  /// No description provided for @colAssessed.
  ///
  /// In en, this message translates to:
  /// **'Assessed'**
  String get colAssessed;

  /// No description provided for @colSheet.
  ///
  /// In en, this message translates to:
  /// **'Sheet'**
  String get colSheet;

  /// No description provided for @geneNamesAllKnown.
  ///
  /// In en, this message translates to:
  /// **'All gene names recognized'**
  String get geneNamesAllKnown;

  /// No description provided for @geneReviewAgain.
  ///
  /// In en, this message translates to:
  /// **'Review gene names again ({n})'**
  String geneReviewAgain(String n);

  /// No description provided for @geneReviewWhichTitle.
  ///
  /// In en, this message translates to:
  /// **'The review shows these genes again:'**
  String get geneReviewWhichTitle;

  /// No description provided for @geneReviewWhichUnknown.
  ///
  /// In en, this message translates to:
  /// **'Not found in SynGenes: {names}'**
  String geneReviewWhichUnknown(String names);

  /// No description provided for @geneReviewWhichConflicts.
  ///
  /// In en, this message translates to:
  /// **'/gene and /product disagree: {pairs}'**
  String geneReviewWhichConflicts(String pairs);

  /// No description provided for @clipkitCodonMode.
  ///
  /// In en, this message translates to:
  /// **'Coding markers: ClipKIT mode (--codon)'**
  String get clipkitCodonMode;

  /// No description provided for @clipkitCodonInfo.
  ///
  /// In en, this message translates to:
  /// **'To keep the reading frame, the coding markers are trimmed by ClipKIT with its --codon option: when one position of a codon should go, ClipKIT removes the whole codon (Steenwyk et al. 2020). Non-coding markers are trimmed by trimAl with the method above.'**
  String get clipkitCodonInfo;

  /// No description provided for @outgroupSuggesting.
  ///
  /// In en, this message translates to:
  /// **'Measuring the distances between the alignments\' taxa to suggest an outgroup…'**
  String get outgroupSuggesting;

  /// No description provided for @outgroupSuggestNone.
  ///
  /// In en, this message translates to:
  /// **'Align the markers to get an outgroup suggested from the distances.'**
  String get outgroupSuggestNone;

  /// No description provided for @outgroupSuggestTitle.
  ///
  /// In en, this message translates to:
  /// **'Suggested outgroup (alignment distances)'**
  String get outgroupSuggestTitle;

  /// No description provided for @outgroupSuggestWeak.
  ///
  /// In en, this message translates to:
  /// **'No clearly separate group: the most distant taxon is suggested'**
  String get outgroupSuggestWeak;

  /// No description provided for @outgroupSuggestInfo.
  ///
  /// In en, this message translates to:
  /// **'p-distances between every pair of taxa over the sites both have in each alignment, summed across markers; an average-linkage (UPGMA) tree of them; the smaller of the two groups joined last is the one farthest from the rest. The value beside each taxon is its mean distance to all the others. Check it against what you know of the group: a long branch (fast-evolving or contaminated sequence) can also look distant.'**
  String get outgroupSuggestInfo;

  /// No description provided for @outgroupSuggestApplied.
  ///
  /// In en, this message translates to:
  /// **'In use'**
  String get outgroupSuggestApplied;

  /// No description provided for @outgroupSuggestUse.
  ///
  /// In en, this message translates to:
  /// **'Use as outgroup'**
  String get outgroupSuggestUse;

  /// No description provided for @outgroupSuggestStats.
  ///
  /// In en, this message translates to:
  /// **'Mean p-distance to the ingroup {between}; within the ingroup {within} · {loci} alignments, {taxa} taxa.'**
  String outgroupSuggestStats(
    String between,
    String within,
    String loci,
    String taxa,
  );

  /// No description provided for @outgroupSuggestNext.
  ///
  /// In en, this message translates to:
  /// **'Next most distant'**
  String get outgroupSuggestNext;

  /// No description provided for @treeArrowShape.
  ///
  /// In en, this message translates to:
  /// **'Arrow shape'**
  String get treeArrowShape;

  /// No description provided for @treeArrowStraight.
  ///
  /// In en, this message translates to:
  /// **'Straight'**
  String get treeArrowStraight;

  /// No description provided for @treeArrowCurved.
  ///
  /// In en, this message translates to:
  /// **'Curved'**
  String get treeArrowCurved;

  /// No description provided for @treeArrowShapeInfo.
  ///
  /// In en, this message translates to:
  /// **'Values never overlap: when a value\'s place is taken it moves away and an arrow points to its node. The arrows can be straight or curved.'**
  String get treeArrowShapeInfo;

  /// No description provided for @treeSupportClasses.
  ///
  /// In en, this message translates to:
  /// **'Support classes'**
  String get treeSupportClasses;

  /// No description provided for @treeClassesAuto.
  ///
  /// In en, this message translates to:
  /// **'Automatic ({kind})'**
  String treeClassesAuto(String kind);

  /// No description provided for @treeClassesN.
  ///
  /// In en, this message translates to:
  /// **'{n} classes'**
  String treeClassesN(String n);

  /// No description provided for @treeClassesTertiles.
  ///
  /// In en, this message translates to:
  /// **'Tertiles of this tree'**
  String get treeClassesTertiles;

  /// No description provided for @treeClassHidden.
  ///
  /// In en, this message translates to:
  /// **'below \"show support from\"'**
  String get treeClassHidden;

  /// No description provided for @treeClassRemove.
  ///
  /// In en, this message translates to:
  /// **'Remove this class'**
  String get treeClassRemove;

  /// No description provided for @treeClassAdd.
  ///
  /// In en, this message translates to:
  /// **'Add a class'**
  String get treeClassAdd;

  /// No description provided for @tipFieldConservation.
  ///
  /// In en, this message translates to:
  /// **'Category (SALVE)'**
  String get tipFieldConservation;

  /// No description provided for @treeHeatStatus.
  ///
  /// In en, this message translates to:
  /// **'Category'**
  String get treeHeatStatus;

  /// No description provided for @treeHeatStatusLegend.
  ///
  /// In en, this message translates to:
  /// **'Extinction risk (SALVE / ICMBio)'**
  String get treeHeatStatusLegend;

  /// No description provided for @treeHeatYes.
  ///
  /// In en, this message translates to:
  /// **'yes'**
  String get treeHeatYes;

  /// No description provided for @treeHeatNo.
  ///
  /// In en, this message translates to:
  /// **'no'**
  String get treeHeatNo;

  /// No description provided for @treeHeatmap.
  ///
  /// In en, this message translates to:
  /// **'SALVE heatmap'**
  String get treeHeatmap;

  /// No description provided for @treeHeatmapInfo.
  ///
  /// In en, this message translates to:
  /// **'Coloured columns beside the labels from the species\' SALVE sheets: category in the IUCN colours, endemic to Brazil, and one column per biome. Rectangular layout only. The same information can go into the labels as text (label parts above).'**
  String get treeHeatmapInfo;

  /// No description provided for @satParamsTitle.
  ///
  /// In en, this message translates to:
  /// **'Test parameters'**
  String get satParamsTitle;

  /// No description provided for @satParamsBody.
  ///
  /// In en, this message translates to:
  /// **'How the saturation of each marker is measured and when it is called moderate or saturated. Changes apply on the next run.'**
  String get satParamsBody;

  /// No description provided for @satParamMethod.
  ///
  /// In en, this message translates to:
  /// **'Method'**
  String get satParamMethod;

  /// No description provided for @satParamMethodValue.
  ///
  /// In en, this message translates to:
  /// **'[[PhyKIT]] saturation, method of [[Philippe 2011|]]: slope of the uncorrected p-distance of every pair of taxa on their patristic distance in the tree'**
  String get satParamMethodValue;

  /// No description provided for @satParamDistance.
  ///
  /// In en, this message translates to:
  /// **'Trees'**
  String get satParamDistance;

  /// No description provided for @satParamDistanceValue.
  ///
  /// In en, this message translates to:
  /// **'IQ-TREE, one per gene (fast search) with the model below; the branch lengths of each codon subset are estimated on that topology (-te)'**
  String get satParamDistanceValue;

  /// No description provided for @satParamSubsets.
  ///
  /// In en, this message translates to:
  /// **'Sites tested'**
  String get satParamSubsets;

  /// No description provided for @satParamSubsetsValue.
  ///
  /// In en, this message translates to:
  /// **'1st, 2nd and 3rd codon positions, 1st+2nd and all sites (non-coding markers: all sites)'**
  String get satParamSubsetsValue;

  /// No description provided for @satParamCode.
  ///
  /// In en, this message translates to:
  /// **'Genetic code'**
  String get satParamCode;

  /// No description provided for @satParamFrame.
  ///
  /// In en, this message translates to:
  /// **'Reading frame'**
  String get satParamFrame;

  /// No description provided for @satParamFrameValue.
  ///
  /// In en, this message translates to:
  /// **'from the alignment (CDS start, or the frame with fewest internal stops)'**
  String get satParamFrameValue;

  /// No description provided for @satParamsThresholds.
  ///
  /// In en, this message translates to:
  /// **'Thresholds and settings'**
  String get satParamsThresholds;

  /// No description provided for @satParamSlopeModerate.
  ///
  /// In en, this message translates to:
  /// **'Saturation · moderate <'**
  String get satParamSlopeModerate;

  /// No description provided for @satParamSlopeSaturated.
  ///
  /// In en, this message translates to:
  /// **'Saturation · saturated <'**
  String get satParamSlopeSaturated;

  /// No description provided for @satParamSlopeInfo.
  ///
  /// In en, this message translates to:
  /// **'PhyKIT\'s saturation value: slope of the p-distance on the patristic distance. 1 = no saturation; the closer to 0, the more saturated. The defaults (0.5 and 0.3) are heuristic: the 2nd codon positions usually sit near 1 and saturated 3rd positions well below 0.5.'**
  String get satParamSlopeInfo;

  /// No description provided for @satParamIss.
  ///
  /// In en, this message translates to:
  /// **'Iss · saturated from'**
  String get satParamIss;

  /// No description provided for @satParamIssInfo.
  ///
  /// In en, this message translates to:
  /// **'Xia\'s index: mean site entropy over the entropy expected under full saturation. Near 1, the sites look random.'**
  String get satParamIssInfo;

  /// No description provided for @satParamUndefinedModerate.
  ///
  /// In en, this message translates to:
  /// **'No correction · moderate from'**
  String get satParamUndefinedModerate;

  /// No description provided for @satParamUndefinedSaturated.
  ///
  /// In en, this message translates to:
  /// **'No correction · saturated from'**
  String get satParamUndefinedSaturated;

  /// No description provided for @satParamUndefinedInfo.
  ///
  /// In en, this message translates to:
  /// **'Pairs whose corrected distance cannot be computed (beyond what the model corrects): a strong sign of saturation.'**
  String get satParamUndefinedInfo;

  /// No description provided for @satParamLowDivergence.
  ///
  /// In en, this message translates to:
  /// **'Low divergence · p below'**
  String get satParamLowDivergence;

  /// No description provided for @satParamLowDivergenceInfo.
  ///
  /// In en, this message translates to:
  /// **'With a mean p-distance this low, saturation does not matter and the subset is called OK.'**
  String get satParamLowDivergenceInfo;

  /// No description provided for @satParamMinSites.
  ///
  /// In en, this message translates to:
  /// **'Minimum sites'**
  String get satParamMinSites;

  /// No description provided for @satParamMinPairs.
  ///
  /// In en, this message translates to:
  /// **'Minimum pairs'**
  String get satParamMinPairs;

  /// No description provided for @satParamReplicates.
  ///
  /// In en, this message translates to:
  /// **'Iss · replicates'**
  String get satParamReplicates;

  /// No description provided for @satParamReplicatesInfo.
  ///
  /// In en, this message translates to:
  /// **'Random sites drawn to estimate the entropy expected under full saturation, for each number of sequences. More replicates, steadier Iss, longer run.'**
  String get satParamReplicatesInfo;

  /// No description provided for @satParamPlotPoints.
  ///
  /// In en, this message translates to:
  /// **'Plot points'**
  String get satParamPlotPoints;

  /// No description provided for @satParamPlotPointsInfo.
  ///
  /// In en, this message translates to:
  /// **'Pairs of sequences drawn in the plots (a sample when there are more); the statistics use every pair.'**
  String get satParamPlotPointsInfo;

  /// No description provided for @satParamsReset.
  ///
  /// In en, this message translates to:
  /// **'Restore defaults'**
  String get satParamsReset;

  /// No description provided for @satParamsStale.
  ///
  /// In en, this message translates to:
  /// **'The parameters changed since the last run: run the test again to update the results.'**
  String get satParamsStale;

  /// No description provided for @colTrend.
  ///
  /// In en, this message translates to:
  /// **'Population trend'**
  String get colTrend;

  /// No description provided for @colGeneration.
  ///
  /// In en, this message translates to:
  /// **'Generation time'**
  String get colGeneration;

  /// No description provided for @colMigratory.
  ///
  /// In en, this message translates to:
  /// **'Migratory'**
  String get colMigratory;

  /// No description provided for @colHabitats.
  ///
  /// In en, this message translates to:
  /// **'Habitats'**
  String get colHabitats;

  /// No description provided for @colThreats.
  ///
  /// In en, this message translates to:
  /// **'Threats'**
  String get colThreats;

  /// No description provided for @colOfficialList.
  ///
  /// In en, this message translates to:
  /// **'Official list'**
  String get colOfficialList;

  /// No description provided for @colUcs.
  ///
  /// In en, this message translates to:
  /// **'Protected areas'**
  String get colUcs;

  /// No description provided for @salveReferences.
  ///
  /// In en, this message translates to:
  /// **'References in the sheet'**
  String get salveReferences;

  /// No description provided for @salveHistory.
  ///
  /// In en, this message translates to:
  /// **'Assessment history'**
  String get salveHistory;

  /// No description provided for @salveLists.
  ///
  /// In en, this message translates to:
  /// **'Official lists and conventions'**
  String get salveLists;

  /// No description provided for @salveSecDistribution.
  ///
  /// In en, this message translates to:
  /// **'Distribution'**
  String get salveSecDistribution;

  /// No description provided for @salveNational.
  ///
  /// In en, this message translates to:
  /// **'In Brazil'**
  String get salveNational;

  /// No description provided for @salveGlobal.
  ///
  /// In en, this message translates to:
  /// **'Worldwide'**
  String get salveGlobal;

  /// No description provided for @salveSecNatural.
  ///
  /// In en, this message translates to:
  /// **'Natural history and population'**
  String get salveSecNatural;

  /// No description provided for @salveDiet.
  ///
  /// In en, this message translates to:
  /// **'Diet'**
  String get salveDiet;

  /// No description provided for @salveReproduction.
  ///
  /// In en, this message translates to:
  /// **'Reproduction'**
  String get salveReproduction;

  /// No description provided for @salveSexRatio.
  ///
  /// In en, this message translates to:
  /// **'Sex ratio'**
  String get salveSexRatio;

  /// No description provided for @salveNaturalHistory.
  ///
  /// In en, this message translates to:
  /// **'Natural history'**
  String get salveNaturalHistory;

  /// No description provided for @salveInteractions.
  ///
  /// In en, this message translates to:
  /// **'Interactions with other species'**
  String get salveInteractions;

  /// No description provided for @salvePopulation.
  ///
  /// In en, this message translates to:
  /// **'Population'**
  String get salvePopulation;

  /// No description provided for @salveGenetics.
  ///
  /// In en, this message translates to:
  /// **'Genetic characteristics'**
  String get salveGenetics;

  /// No description provided for @salveSecThreats.
  ///
  /// In en, this message translates to:
  /// **'Threats, uses, conservation and research'**
  String get salveSecThreats;

  /// No description provided for @salveUses.
  ///
  /// In en, this message translates to:
  /// **'Uses'**
  String get salveUses;

  /// No description provided for @salveActions.
  ///
  /// In en, this message translates to:
  /// **'Conservation actions'**
  String get salveActions;

  /// No description provided for @salveConservation.
  ///
  /// In en, this message translates to:
  /// **'Conservation'**
  String get salveConservation;

  /// No description provided for @salveResearch.
  ///
  /// In en, this message translates to:
  /// **'Research'**
  String get salveResearch;

  /// No description provided for @trendDeclining.
  ///
  /// In en, this message translates to:
  /// **'declining'**
  String get trendDeclining;

  /// No description provided for @trendStable.
  ///
  /// In en, this message translates to:
  /// **'stable'**
  String get trendStable;

  /// No description provided for @trendIncreasing.
  ///
  /// In en, this message translates to:
  /// **'increasing'**
  String get trendIncreasing;

  /// No description provided for @trendUnknown.
  ///
  /// In en, this message translates to:
  /// **'unknown'**
  String get trendUnknown;

  /// No description provided for @treeHeatFew.
  ///
  /// In en, this message translates to:
  /// **'few'**
  String get treeHeatFew;

  /// No description provided for @treeHeatMany.
  ///
  /// In en, this message translates to:
  /// **'many'**
  String get treeHeatMany;

  /// No description provided for @trimNothingRemoved.
  ///
  /// In en, this message translates to:
  /// **'{method} kept all {n} columns: this alignment has no gappy region to remove, so nothing is grey (the kpic modes also remove uninformative sites)'**
  String trimNothingRemoved(String method, String n);

  /// No description provided for @treeBranchWidth.
  ///
  /// In en, this message translates to:
  /// **'Branch width'**
  String get treeBranchWidth;

  /// No description provided for @treeBranchWidthInfo.
  ///
  /// In en, this message translates to:
  /// **'Width of the tree itself (branches), in cm. 0 = automatic: what is left beside labels, heatmap, bars and legend, never under 45 % of the figure — with many extras the figure grows wider instead of squeezing the tree. Set a width to fix it; the figure grows to fit the rest.'**
  String get treeBranchWidthInfo;

  /// No description provided for @treeHeatCell.
  ///
  /// In en, this message translates to:
  /// **'Heatmap cell'**
  String get treeHeatCell;

  /// No description provided for @treeHeatCellInfo.
  ///
  /// In en, this message translates to:
  /// **'Side of each heatmap cell, in pt. 0 = automatic (from the row height).'**
  String get treeHeatCellInfo;

  /// No description provided for @citeVersions.
  ///
  /// In en, this message translates to:
  /// **'Versions used'**
  String get citeVersions;

  /// No description provided for @satParamVersions.
  ///
  /// In en, this message translates to:
  /// **'Programs'**
  String get satParamVersions;

  /// No description provided for @satParamModel.
  ///
  /// In en, this message translates to:
  /// **'IQ-TREE model'**
  String get satParamModel;

  /// No description provided for @satParamModelInfo.
  ///
  /// In en, this message translates to:
  /// **'Model of the gene trees and of the branch lengths of each codon subset (IQ-TREE -m), e.g. GTR+G, HKY+G, GTR+F+I+G4.'**
  String get satParamModelInfo;

  /// No description provided for @satParamGaps.
  ///
  /// In en, this message translates to:
  /// **'Gaps'**
  String get satParamGaps;

  /// No description provided for @satParamGapsSwitch.
  ///
  /// In en, this message translates to:
  /// **'Ignore sites with gaps (-e)'**
  String get satParamGapsSwitch;

  /// No description provided for @satParamGapsInfo.
  ///
  /// In en, this message translates to:
  /// **'PhyKIT --exclude_gaps: sites with a gap in a pair are left out of its p-distance.'**
  String get satParamGapsInfo;

  /// No description provided for @phykitRequired.
  ///
  /// In en, this message translates to:
  /// **'The saturation test uses PhyKIT, which is not installed yet. Install it here (needs Python 3) or with: pip install phykit'**
  String get phykitRequired;

  /// No description provided for @installPhykit.
  ///
  /// In en, this message translates to:
  /// **'Install PhyKIT'**
  String get installPhykit;

  /// No description provided for @plotPhykit.
  ///
  /// In en, this message translates to:
  /// **'p × tree distance · PhyKIT {s}'**
  String plotPhykit(String s);

  /// No description provided for @plotPhykitBody.
  ///
  /// In en, this message translates to:
  /// **'Each point is a pair of taxa: uncorrected p-distance against their patristic distance in the IQ-TREE tree of these sites. Points far below the diagonal hide repeated substitutions.'**
  String get plotPhykitBody;

  /// No description provided for @plotPositionsPatristic.
  ///
  /// In en, this message translates to:
  /// **'Codon positions × tree distance'**
  String get plotPositionsPatristic;

  /// No description provided for @patristicDistance.
  ///
  /// In en, this message translates to:
  /// **'Patristic distance'**
  String get patristicDistance;

  /// No description provided for @metricMeanPatristic.
  ///
  /// In en, this message translates to:
  /// **'Mean patristic'**
  String get metricMeanPatristic;

  /// No description provided for @metricPhykit.
  ///
  /// In en, this message translates to:
  /// **'PhyKIT saturation'**
  String get metricPhykit;

  /// No description provided for @trimToolLabel.
  ///
  /// In en, this message translates to:
  /// **'Trimming tool'**
  String get trimToolLabel;

  /// No description provided for @trimToolTrimal.
  ///
  /// In en, this message translates to:
  /// **'trimAl · column by column'**
  String get trimToolTrimal;

  /// No description provided for @trimToolClipkit.
  ///
  /// In en, this message translates to:
  /// **'ClipKIT · whole codons'**
  String get trimToolClipkit;

  /// No description provided for @trimToolTrimalBody.
  ///
  /// In en, this message translates to:
  /// **'trimAl trims every marker column by column. In coding genes this can cut codons apart: their reading frame is lost, they become one partition without codon positions and stay out of the saturation test.'**
  String get trimToolTrimalBody;

  /// No description provided for @trimToolClipkitBody.
  ///
  /// In en, this message translates to:
  /// **'ClipKIT trims every marker. In the coding genes it runs with --codon: when one position of a codon should go, the whole codon goes, so the reading frame and the codon positions are kept. Non-coding markers are trimmed by ClipKIT without --codon. trimAl is not used.'**
  String get trimToolClipkitBody;

  /// No description provided for @mtShowBases.
  ///
  /// In en, this message translates to:
  /// **'Bases'**
  String get mtShowBases;

  /// No description provided for @mtShowPercent.
  ///
  /// In en, this message translates to:
  /// **'% of the marker'**
  String get mtShowPercent;

  /// No description provided for @mtOrderMatrix.
  ///
  /// In en, this message translates to:
  /// **'Matrix order'**
  String get mtOrderMatrix;

  /// No description provided for @mtOrderComplete.
  ///
  /// In en, this message translates to:
  /// **'Most complete first'**
  String get mtOrderComplete;

  /// No description provided for @mtScaleComplete.
  ///
  /// In en, this message translates to:
  /// **'Completeness'**
  String get mtScaleComplete;

  /// No description provided for @mtScaleMissing.
  ///
  /// In en, this message translates to:
  /// **'Missing / gaps'**
  String get mtScaleMissing;

  /// No description provided for @treeMarkShape.
  ///
  /// In en, this message translates to:
  /// **'Point shape'**
  String get treeMarkShape;

  /// No description provided for @treeShapeCircle.
  ///
  /// In en, this message translates to:
  /// **'Circle'**
  String get treeShapeCircle;

  /// No description provided for @treeShapeSquare.
  ///
  /// In en, this message translates to:
  /// **'Square'**
  String get treeShapeSquare;

  /// No description provided for @treeShapeDiamond.
  ///
  /// In en, this message translates to:
  /// **'Diamond'**
  String get treeShapeDiamond;

  /// No description provided for @treeShapeTriangle.
  ///
  /// In en, this message translates to:
  /// **'Triangle'**
  String get treeShapeTriangle;

  /// No description provided for @treeBarWidth.
  ///
  /// In en, this message translates to:
  /// **'Clade bar width'**
  String get treeBarWidth;

  /// No description provided for @treeBarWidthInfo.
  ///
  /// In en, this message translates to:
  /// **'Width of the bars beside clades, in pt. 0 = automatic (from the branch width).'**
  String get treeBarWidthInfo;

  /// No description provided for @treeBarFont.
  ///
  /// In en, this message translates to:
  /// **'Clade label size'**
  String get treeBarFont;

  /// No description provided for @treeBarFontInfo.
  ///
  /// In en, this message translates to:
  /// **'Size of the names beside the clade bars, in pt. 0 = the tip label size.'**
  String get treeBarFontInfo;

  /// No description provided for @treeHeatGap.
  ///
  /// In en, this message translates to:
  /// **'Heatmap column gap'**
  String get treeHeatGap;

  /// No description provided for @trimToolTrimalCodon.
  ///
  /// In en, this message translates to:
  /// **'trimAl · whole codons'**
  String get trimToolTrimalCodon;

  /// No description provided for @trimToolTrimalCodonBody.
  ///
  /// In en, this message translates to:
  /// **'trimAl trims every marker. Coding genes are read codon by codon and translated with the project\'s genetic code; trimAl trims that protein alignment and, with -backtrans, writes back the codons of the kept amino acids, so whole codons go or stay and the codon positions are kept. Non-coding markers are trimmed column by column.'**
  String get trimToolTrimalCodonBody;

  /// No description provided for @frameStepCodonTool.
  ///
  /// In en, this message translates to:
  /// **'{tool} takes the whole codon'**
  String frameStepCodonTool(String tool);

  /// No description provided for @trimCaseCodonsTool.
  ///
  /// In en, this message translates to:
  /// **'Keep the reading frame ({tool})'**
  String trimCaseCodonsTool(String tool);

  /// No description provided for @trimCasesCaptionTool.
  ///
  /// In en, this message translates to:
  /// **'Same alignment, two ways. Column by column, trimAl removes the columns it marks (×) one by one, so the codon positions (coloured ruler) fall out of order (!). With {tool} the whole codon of every marked column goes (+): the 1st, 2nd and 3rd positions stay in order for the partitions and the saturation test.'**
  String trimCasesCaptionTool(String tool);

  /// No description provided for @delimTitle.
  ///
  /// In en, this message translates to:
  /// **'Species delimitation'**
  String get delimTitle;

  /// No description provided for @delimIntro.
  ///
  /// In en, this message translates to:
  /// **'IQ-TREE and RAxML-NG infer the tree but do not delimit species, so SPLACE runs published delimitation programs: [[mPTP]] and [[GMYC]] on the maximum-likelihood trees, and [[ASAP]] and [[ABGD]] on the alignments (pairwise distances). Each result gets its own button on the tree (step 8). A delimitation is a hypothesis: compare the methods and confirm it with other evidence.'**
  String get delimIntro;

  /// No description provided for @delimTrees.
  ///
  /// In en, this message translates to:
  /// **'ML trees to delimit'**
  String get delimTrees;

  /// No description provided for @delimModels.
  ///
  /// In en, this message translates to:
  /// **'Models'**
  String get delimModels;

  /// No description provided for @delimModelMulti.
  ///
  /// In en, this message translates to:
  /// **'mPTP (one rate per species)'**
  String get delimModelMulti;

  /// No description provided for @delimModelSingle.
  ///
  /// In en, this message translates to:
  /// **'PTP (one rate for all)'**
  String get delimModelSingle;

  /// No description provided for @delimMinbr.
  ///
  /// In en, this message translates to:
  /// **'Minimum branch length from the data (--minbr_auto)'**
  String get delimMinbr;

  /// No description provided for @delimMinbrInfo.
  ///
  /// In en, this message translates to:
  /// **'Branches shorter than the smallest p-distance of the supermatrix carry no substitution and are ignored. Off: mPTP\'s default (0.0001).'**
  String get delimMinbrInfo;

  /// No description provided for @delimCrop.
  ///
  /// In en, this message translates to:
  /// **'Root on the outgroup and leave it out (--outgroup_crop)'**
  String get delimCrop;

  /// No description provided for @delimCropInfo.
  ///
  /// In en, this message translates to:
  /// **'The outgroups of step 6c root the tree and are not delimited, so their long branches do not change the result.'**
  String get delimCropInfo;

  /// No description provided for @delimRun.
  ///
  /// In en, this message translates to:
  /// **'Delimit species'**
  String get delimRun;

  /// No description provided for @delimRunning.
  ///
  /// In en, this message translates to:
  /// **'mPTP: species delimitation'**
  String get delimRunning;

  /// No description provided for @delimNone.
  ///
  /// In en, this message translates to:
  /// **'No delimitation yet.'**
  String get delimNone;

  /// No description provided for @delimSpecies.
  ///
  /// In en, this message translates to:
  /// **'{n} species'**
  String delimSpecies(String n);

  /// No description provided for @delimSingletons.
  ///
  /// In en, this message translates to:
  /// **'{n} with one tip'**
  String delimSingletons(String n);

  /// No description provided for @delimLrtFailed.
  ///
  /// In en, this message translates to:
  /// **'LRT failed: one species preferred'**
  String get delimLrtFailed;

  /// No description provided for @delimViewOnTree.
  ///
  /// In en, this message translates to:
  /// **'Show on the tree'**
  String get delimViewOnTree;

  /// No description provided for @delimMinbrCommand.
  ///
  /// In en, this message translates to:
  /// **'Minimum branch length'**
  String get delimMinbrCommand;

  /// No description provided for @delimSpeciesList.
  ///
  /// In en, this message translates to:
  /// **'Delimited species'**
  String get delimSpeciesList;

  /// No description provided for @treeDelimitation.
  ///
  /// In en, this message translates to:
  /// **'Species delimitation'**
  String get treeDelimitation;

  /// No description provided for @treeDelimitationInfo.
  ///
  /// In en, this message translates to:
  /// **'Each button adds a column of bars beside the labels: one bar per species delimited by that tool on that tree (consecutive species alternate colours). Tips are matched by name, so a delimitation can be compared on any tree. Run them in step 7e.'**
  String get treeDelimitationInfo;

  /// No description provided for @treeDelimitationCircular.
  ///
  /// In en, this message translates to:
  /// **'Shown in the rectangular layout only'**
  String get treeDelimitationCircular;

  /// No description provided for @citeFeatured.
  ///
  /// In en, this message translates to:
  /// **'Gene names and taxonomy'**
  String get citeFeatured;

  /// No description provided for @citeRoleHint.
  ///
  /// In en, this message translates to:
  /// **'Each reference says what the tool did in this project.'**
  String get citeRoleHint;

  /// No description provided for @citeRoleLabel.
  ///
  /// In en, this message translates to:
  /// **'What it did:'**
  String get citeRoleLabel;

  /// No description provided for @delimOnTrees.
  ///
  /// In en, this message translates to:
  /// **'On the trees (mPTP, GMYC)'**
  String get delimOnTrees;

  /// No description provided for @delimTreesIntro.
  ///
  /// In en, this message translates to:
  /// **'mPTP finds where the branching of the ML tree changes from speciation (between species) to coalescence (within species). The [[PTP]] model uses one coalescent rate for all species; mPTP gives each species its own rate, which copes better with uneven sampling. [[GMYC]] (single threshold, pGMYC) reads the ML tree made ultrametric by [[LSD2]] (relative dates) and finds where the waiting times between branchings change from speciation to coalescence.'**
  String get delimTreesIntro;

  /// No description provided for @delimOnAlignments.
  ///
  /// In en, this message translates to:
  /// **'On the alignments (ASAP, ABGD)'**
  String get delimOnAlignments;

  /// No description provided for @delimAlignmentsIntro.
  ///
  /// In en, this message translates to:
  /// **'Single-locus methods based on pairwise distances, usually run on a barcode marker (COI, 16S…). ASAP builds partitions by hierarchical clustering and ranks them by its score; ABGD looks for the barcode gap between intra- and interspecific distances over a range of priors.'**
  String get delimAlignmentsIntro;

  /// No description provided for @delimAlignments.
  ///
  /// In en, this message translates to:
  /// **'Alignments'**
  String get delimAlignments;

  /// No description provided for @delimNoAlignments.
  ///
  /// In en, this message translates to:
  /// **'No DNA alignment yet.'**
  String get delimNoAlignments;

  /// No description provided for @delimSupermatrix.
  ///
  /// In en, this message translates to:
  /// **'Supermatrix'**
  String get delimSupermatrix;

  /// No description provided for @delimDistance.
  ///
  /// In en, this message translates to:
  /// **'Distance'**
  String get delimDistance;

  /// No description provided for @delimDistanceInfo.
  ///
  /// In en, this message translates to:
  /// **'Pairwise distance of ASAP and ABGD (both default to Jukes-Cantor). ABGD offers JC69 and K2P only; with another distance ABGD uses JC69.'**
  String get delimDistanceInfo;

  /// No description provided for @delimPMin.
  ///
  /// In en, this message translates to:
  /// **'ABGD prior Pmin'**
  String get delimPMin;

  /// No description provided for @delimPMax.
  ///
  /// In en, this message translates to:
  /// **'ABGD prior Pmax'**
  String get delimPMax;

  /// No description provided for @delimSteps.
  ///
  /// In en, this message translates to:
  /// **'ABGD steps'**
  String get delimSteps;

  /// No description provided for @delimSlope.
  ///
  /// In en, this message translates to:
  /// **'ABGD gap width (X)'**
  String get delimSlope;

  /// No description provided for @delimPriorInfo.
  ///
  /// In en, this message translates to:
  /// **'ABGD tries priors of maximum intraspecific divergence from Pmin to Pmax (log scale) and gives one partition per prior. Defaults are ABGD\'s own (0.001 to 0.1, 10 steps).'**
  String get delimPriorInfo;

  /// No description provided for @delimSlopeInfo.
  ///
  /// In en, this message translates to:
  /// **'Relative gap width (ABGD -X, default 1.5): lower values find smaller gaps.'**
  String get delimSlopeInfo;

  /// No description provided for @delimPartition.
  ///
  /// In en, this message translates to:
  /// **'Partition shown'**
  String get delimPartition;

  /// No description provided for @delimPartitionAsapInfo.
  ///
  /// In en, this message translates to:
  /// **'ASAP proposes several partitions ranked by its score (lower is better); the best one is shown first. The one chosen here goes to the tree and to the methods text.'**
  String get delimPartitionAsapInfo;

  /// No description provided for @delimPartitionAbgdInfo.
  ///
  /// In en, this message translates to:
  /// **'ABGD gives an initial and a recursive partition for each prior P. Shown first: the initial partition at the prior nearest 0.01; look for partitions that stay the same over several priors. The one chosen here goes to the tree and to the methods text.'**
  String get delimPartitionAbgdInfo;

  /// No description provided for @delimAsapScore.
  ///
  /// In en, this message translates to:
  /// **'ASAP score {score}'**
  String delimAsapScore(String score);

  /// No description provided for @delimInitial.
  ///
  /// In en, this message translates to:
  /// **'initial'**
  String get delimInitial;

  /// No description provided for @delimRecursive.
  ///
  /// In en, this message translates to:
  /// **'recursive'**
  String get delimRecursive;

  /// No description provided for @delimModelGmyc.
  ///
  /// In en, this message translates to:
  /// **'GMYC (ultrametric tree by LSD2)'**
  String get delimModelGmyc;

  /// No description provided for @delimLsd2Command.
  ///
  /// In en, this message translates to:
  /// **'Ultrametric tree (LSD2)'**
  String get delimLsd2Command;

  /// No description provided for @exportGbButton.
  ///
  /// In en, this message translates to:
  /// **'Export corrected GenBank'**
  String get exportGbButton;

  /// No description provided for @exportGbChoose.
  ///
  /// In en, this message translates to:
  /// **'Save here'**
  String get exportGbChoose;

  /// No description provided for @exportGbInfo.
  ///
  /// In en, this message translates to:
  /// **'Writes every record as a GenBank file (.gb) with the species name as identified or edited, the lineage of the taxonomy and the gene names standardized by SynGenes (the original name is kept in a /note). Useful for a newly assembled genome before submission. A copy goes to the project folder (00_genbank_corrected).'**
  String get exportGbInfo;

  /// No description provided for @barcodingSwitch.
  ///
  /// In en, this message translates to:
  /// **'These sequences are DNA barcodes'**
  String get barcodingSwitch;

  /// No description provided for @barcodingSwitchInfo.
  ///
  /// In en, this message translates to:
  /// **'One specimen per sequence (COI, 16S, ITS…): identify each one against GenBank and BOLD and export the BOLD submission spreadsheet.'**
  String get barcodingSwitchInfo;

  /// No description provided for @barcodingLocus.
  ///
  /// In en, this message translates to:
  /// **'Barcode marker'**
  String get barcodingLocus;

  /// No description provided for @barcodingLocusInfo.
  ///
  /// In en, this message translates to:
  /// **'The locus of the FASTA files with the barcodes. COI is sent to BOLD and GenBank; other markers to GenBank only.'**
  String get barcodingLocusInfo;

  /// No description provided for @barcodingCount.
  ///
  /// In en, this message translates to:
  /// **'{n} specimens'**
  String barcodingCount(String n);

  /// No description provided for @barcodingSourcesCoi.
  ///
  /// In en, this message translates to:
  /// **'BOLD ID Engine + GenBank BLAST'**
  String get barcodingSourcesCoi;

  /// No description provided for @barcodingSourcesGenbank.
  ///
  /// In en, this message translates to:
  /// **'GenBank BLAST (BOLD identifies COI only)'**
  String get barcodingSourcesGenbank;

  /// No description provided for @barcodingIdentify.
  ///
  /// In en, this message translates to:
  /// **'Identify specimens'**
  String get barcodingIdentify;

  /// No description provided for @barcodingTaxonomy.
  ///
  /// In en, this message translates to:
  /// **'Taxonomy (GBIF)'**
  String get barcodingTaxonomy;

  /// No description provided for @barcodingTaxonomyInfo.
  ///
  /// In en, this message translates to:
  /// **'Fills phylum, class, order and family of the names kept, from GBIF (dataFishing method), for the BOLD sheet.'**
  String get barcodingTaxonomyInfo;

  /// No description provided for @barcodingExport.
  ///
  /// In en, this message translates to:
  /// **'Export BOLD spreadsheet (.xlsx)'**
  String get barcodingExport;

  /// No description provided for @barcodingExportInfo.
  ///
  /// In en, this message translates to:
  /// **'Writes the official BOLD5 single-page specimen template (sheet Submissions): sample id, taxonomy and identification filled; voucher and collection data left for you. A FASTA with the sample ids as headers is saved next to it for the sequence upload.'**
  String get barcodingExportInfo;

  /// No description provided for @barcodingExported.
  ///
  /// In en, this message translates to:
  /// **'BOLD spreadsheet saved with {n} specimens (and the sequences FASTA)'**
  String barcodingExported(String n);

  /// No description provided for @barcodingColSample.
  ///
  /// In en, this message translates to:
  /// **'Sample'**
  String get barcodingColSample;

  /// No description provided for @barcodingColMatches.
  ///
  /// In en, this message translates to:
  /// **'Best matches'**
  String get barcodingColMatches;

  /// No description provided for @barcodingColName.
  ///
  /// In en, this message translates to:
  /// **'Name kept'**
  String get barcodingColName;

  /// No description provided for @barcodingNoMatch.
  ///
  /// In en, this message translates to:
  /// **'no match'**
  String get barcodingNoMatch;

  /// No description provided for @barcodingNameHint.
  ///
  /// In en, this message translates to:
  /// **'Genus species'**
  String get barcodingNameHint;

  /// No description provided for @barcodingPick.
  ///
  /// In en, this message translates to:
  /// **'Choose a match'**
  String get barcodingPick;

  /// No description provided for @tutorialTitle.
  ///
  /// In en, this message translates to:
  /// **'How to use SPLACE'**
  String get tutorialTitle;

  /// No description provided for @tutorialSubtitle.
  ///
  /// In en, this message translates to:
  /// **'Each step with what it is for and what to do there. Close it at any time.'**
  String get tutorialSubtitle;

  /// No description provided for @tutorialInterface.
  ///
  /// In en, this message translates to:
  /// **'The interface'**
  String get tutorialInterface;

  /// No description provided for @tutorialWhatToDo.
  ///
  /// In en, this message translates to:
  /// **'What to do'**
  String get tutorialWhatToDo;

  /// No description provided for @tutorialPrev.
  ///
  /// In en, this message translates to:
  /// **'Previous'**
  String get tutorialPrev;

  /// No description provided for @tutorialNext.
  ///
  /// In en, this message translates to:
  /// **'Next'**
  String get tutorialNext;

  /// No description provided for @tutorialFinish.
  ///
  /// In en, this message translates to:
  /// **'Start working'**
  String get tutorialFinish;

  /// No description provided for @tutorialShowAtStart.
  ///
  /// In en, this message translates to:
  /// **'Show this tutorial after Start App'**
  String get tutorialShowAtStart;

  /// No description provided for @tutorialButton.
  ///
  /// In en, this message translates to:
  /// **'Tutorial'**
  String get tutorialButton;

  /// No description provided for @tutorialPageOf.
  ///
  /// In en, this message translates to:
  /// **'{n} of {total}'**
  String tutorialPageOf(String n, String total);

  /// No description provided for @newProject.
  ///
  /// In en, this message translates to:
  /// **'New project'**
  String get newProject;

  /// No description provided for @newProjectTitle.
  ///
  /// In en, this message translates to:
  /// **'Start a new project?'**
  String get newProjectTitle;

  /// No description provided for @newProjectBody.
  ///
  /// In en, this message translates to:
  /// **'SPLACE starts again with no data, in a new project folder, so you can analyse another group. Settings, the NCBI key and the tools found are kept.'**
  String get newProjectBody;

  /// No description provided for @newProjectKept.
  ///
  /// In en, this message translates to:
  /// **'The files of the current project stay in {dir}.'**
  String newProjectKept(String dir);

  /// No description provided for @newProjectConfirm.
  ///
  /// In en, this message translates to:
  /// **'New project'**
  String get newProjectConfirm;

  /// No description provided for @tutorialPlay.
  ///
  /// In en, this message translates to:
  /// **'Play'**
  String get tutorialPlay;

  /// No description provided for @tutorialPause.
  ///
  /// In en, this message translates to:
  /// **'Pause'**
  String get tutorialPause;

  /// No description provided for @tutorialItemOf.
  ///
  /// In en, this message translates to:
  /// **'Item {n} of {total}'**
  String tutorialItemOf(String n, String total);

  /// No description provided for @tutorialNextStep.
  ///
  /// In en, this message translates to:
  /// **'Next: {step}'**
  String tutorialNextStep(String step);

  /// No description provided for @taxSpeciesRecognized.
  ///
  /// In en, this message translates to:
  /// **'Species recognized:'**
  String get taxSpeciesRecognized;

  /// No description provided for @taxSynonyms.
  ///
  /// In en, this message translates to:
  /// **'Also search its synonyms (GBIF)'**
  String get taxSynonyms;

  /// No description provided for @taxSynonymsInfo.
  ///
  /// In en, this message translates to:
  /// **'Records may be filed under an older name (e.g. Bufo marinus for Rhinella marina). NCBI already finds the synonyms it knows; GBIF adds the accepted name and every synonym it lists, all in one NCBI query. The records are imported under their current NCBI name.'**
  String get taxSynonymsInfo;

  /// No description provided for @taxSynonymsLoading.
  ///
  /// In en, this message translates to:
  /// **'Looking up the synonyms in GBIF…'**
  String get taxSynonymsLoading;

  /// No description provided for @taxSynonymsError.
  ///
  /// In en, this message translates to:
  /// **'GBIF could not be reached: {error}'**
  String taxSynonymsError(String error);

  /// No description provided for @taxSynonymsNone.
  ///
  /// In en, this message translates to:
  /// **'GBIF has no synonyms for this name: it is searched alone.'**
  String get taxSynonymsNone;

  /// No description provided for @taxSynonymsFound.
  ///
  /// In en, this message translates to:
  /// **'{n} names searched together:'**
  String taxSynonymsFound(String n);

  /// No description provided for @examplesLabel.
  ///
  /// In en, this message translates to:
  /// **'Example data:'**
  String get examplesLabel;

  /// No description provided for @examplesEmpty.
  ///
  /// In en, this message translates to:
  /// **'No example files yet: put them in {folder} and rebuild the app.'**
  String examplesEmpty(String folder);

  /// No description provided for @examplesFiles.
  ///
  /// In en, this message translates to:
  /// **'Loads the {n} file(s) of {folder}.'**
  String examplesFiles(String n, String folder);

  /// No description provided for @exampleGenbank.
  ///
  /// In en, this message translates to:
  /// **'GenBank example'**
  String get exampleGenbank;

  /// No description provided for @exampleLoci.
  ///
  /// In en, this message translates to:
  /// **'Loci example (FASTA)'**
  String get exampleLoci;

  /// No description provided for @exampleBarcoding.
  ///
  /// In en, this message translates to:
  /// **'Barcoding example (FASTA)'**
  String get exampleBarcoding;

  /// No description provided for @reviewFiledAs.
  ///
  /// In en, this message translates to:
  /// **'filed as'**
  String get reviewFiledAs;

  /// No description provided for @markerBold.
  ///
  /// In en, this message translates to:
  /// **'Also search BOLD Systems'**
  String get markerBold;

  /// No description provided for @markerBoldInfo.
  ///
  /// In en, this message translates to:
  /// **'Public records of BOLD for {markers}. Each record shows where it came from: NCBI, BOLD, or NCBI + BOLD when BOLD mirrors a GenBank record (imported once, from GenBank). Any rank can be searched (species, genus, family…); synonyms BOLD knows are searched too.'**
  String markerBoldInfo(String markers);

  /// No description provided for @markerBoldNone.
  ///
  /// In en, this message translates to:
  /// **'BOLD holds {markers} for this organelle: select one to search BOLD too.'**
  String markerBoldNone(String markers);

  /// No description provided for @siFoundAs.
  ///
  /// In en, this message translates to:
  /// **'also as {names}'**
  String siFoundAs(String names);

  /// No description provided for @siVerdictAccepted.
  ///
  /// In en, this message translates to:
  /// **'The accepted name in GBIF is {name}: choose it to update the record.'**
  String siVerdictAccepted(String name);

  /// No description provided for @headerConflictUndetermined.
  ///
  /// In en, this message translates to:
  /// **'{n} name(s) repeat inside the same marker. Records of one species are joined into one sample, but records not identified to species (sp., cf., aff., a genus only) or a FASTA sample with the same name are different specimens. Name them in step 2 (species identification or the records table) or remove them: the alignment and every step after it stay locked until every name is unique.'**
  String headerConflictUndetermined(String n);

  /// No description provided for @tipFieldCountry.
  ///
  /// In en, this message translates to:
  /// **'Country'**
  String get tipFieldCountry;

  /// No description provided for @tipFieldLocality.
  ///
  /// In en, this message translates to:
  /// **'Locality'**
  String get tipFieldLocality;

  /// No description provided for @schemesOneCoding.
  ///
  /// In en, this message translates to:
  /// **'One marker: \"by locus\" would be a single partition, so it can be partitioned by codon position (its reading frame is known) or kept as one partition.'**
  String get schemesOneCoding;

  /// No description provided for @schemesOneNoCodons.
  ///
  /// In en, this message translates to:
  /// **'One marker without known codon positions (non-coding, protein, fragments without a confirmed reading frame, or frame lost in trimming): the matrix is a single partition.'**
  String get schemesOneNoCodons;

  /// No description provided for @schemesNoCodons.
  ///
  /// In en, this message translates to:
  /// **'No marker has known codon positions (non-coding, protein, fragments without a confirmed reading frame, or frame lost in trimming): the markers can only be partitioned by locus.'**
  String get schemesNoCodons;

  /// No description provided for @trimPlanTitle.
  ///
  /// In en, this message translates to:
  /// **'How each marker will be trimmed'**
  String get trimPlanTitle;

  /// No description provided for @trimWhyNoFrame.
  ///
  /// In en, this message translates to:
  /// **'partial sequences without a confirmed reading frame in the nucleotide alignment (align in Codon-aware mode to keep codons); one partition without codon positions'**
  String get trimWhyNoFrame;

  /// No description provided for @trimWhyNonCoding.
  ///
  /// In en, this message translates to:
  /// **'not protein-coding'**
  String get trimWhyNonCoding;

  /// No description provided for @trimWhyColumns.
  ///
  /// In en, this message translates to:
  /// **'trimAl by columns breaks codons: its codon positions are lost'**
  String get trimWhyColumns;

  /// No description provided for @stepMarkers.
  ///
  /// In en, this message translates to:
  /// **'Features and Markers'**
  String get stepMarkers;

  /// No description provided for @stepMarkersFasta.
  ///
  /// In en, this message translates to:
  /// **'Samples and Markers'**
  String get stepMarkersFasta;

  /// No description provided for @splashStepMarkers.
  ///
  /// In en, this message translates to:
  /// **'Choose which feature types (CDS, rRNA, tRNA) enter the analysis, or show SPLACE where FASTA headers keep the sample and locus names; then pick the markers for the matrix, reviewing missing and duplicated sequences in the presence and length heatmaps.'**
  String get splashStepMarkers;

  /// No description provided for @dropUndetermined.
  ///
  /// In en, this message translates to:
  /// **'Leave out the {n} record(s) not identified to species'**
  String dropUndetermined(String n);

  /// No description provided for @idAlnTitle.
  ///
  /// In en, this message translates to:
  /// **'Alignments of the identification'**
  String get idAlnTitle;

  /// No description provided for @idAlnOpen.
  ///
  /// In en, this message translates to:
  /// **'See the alignments (your sequence and the hits)'**
  String get idAlnOpen;

  /// No description provided for @idAlnAll.
  ///
  /// In en, this message translates to:
  /// **'All sources'**
  String get idAlnAll;

  /// No description provided for @idAlnQuery.
  ///
  /// In en, this message translates to:
  /// **'Query'**
  String get idAlnQuery;

  /// No description provided for @idAlnHit.
  ///
  /// In en, this message translates to:
  /// **'hit'**
  String get idAlnHit;

  /// No description provided for @idAlnCopy.
  ///
  /// In en, this message translates to:
  /// **'Copy as FASTA'**
  String get idAlnCopy;

  /// No description provided for @idAlnNone.
  ///
  /// In en, this message translates to:
  /// **'No hits for this marker.'**
  String get idAlnNone;

  /// No description provided for @idAlnIdentity.
  ///
  /// In en, this message translates to:
  /// **'identity'**
  String get idAlnIdentity;

  /// No description provided for @idAlnCoverage.
  ///
  /// In en, this message translates to:
  /// **'coverage'**
  String get idAlnCoverage;

  /// No description provided for @idAlnMismatches.
  ///
  /// In en, this message translates to:
  /// **'mismatches'**
  String get idAlnMismatches;

  /// No description provided for @idAlnGaps.
  ///
  /// In en, this message translates to:
  /// **'gaps'**
  String get idAlnGaps;

  /// No description provided for @idAlnLoading.
  ///
  /// In en, this message translates to:
  /// **'BOLD gives where the hit aligns, not the alignment: fetching its sequence from the BOLD portal and aligning it…'**
  String get idAlnLoading;

  /// No description provided for @idAlnLoadError.
  ///
  /// In en, this message translates to:
  /// **'The BOLD sequence could not be fetched: {error}'**
  String idAlnLoadError(String error);

  /// No description provided for @idpTitle.
  ///
  /// In en, this message translates to:
  /// **'Search parameters'**
  String get idpTitle;

  /// No description provided for @idpSubtitle.
  ///
  /// In en, this message translates to:
  /// **'GenBank BLAST and the BOLD ID engine. Changing a parameter sends the searches again; answers with the same parameters are reused.'**
  String get idpSubtitle;

  /// No description provided for @idpDatabase.
  ///
  /// In en, this message translates to:
  /// **'Database'**
  String get idpDatabase;

  /// No description provided for @idpBlastDbInfo.
  ///
  /// In en, this message translates to:
  /// **'core_nt: the default nucleotide database of web BLAST (core sequences, faster). nt: the complete nucleotide collection. refseq_select_rna / refseq_rna: RefSeq transcripts only.'**
  String get idpBlastDbInfo;

  /// No description provided for @idpProgram.
  ///
  /// In en, this message translates to:
  /// **'Program'**
  String get idpProgram;

  /// No description provided for @idpProgramInfo.
  ///
  /// In en, this message translates to:
  /// **'megablast: for the same or very close species (fast, long words). blastn: more sensitive, finds more distant relatives (slower).'**
  String get idpProgramInfo;

  /// No description provided for @idpHits.
  ///
  /// In en, this message translates to:
  /// **'Max results'**
  String get idpHits;

  /// No description provided for @idpHitsInfo.
  ///
  /// In en, this message translates to:
  /// **'How many target sequences BLAST returns per query (HITLIST_SIZE).'**
  String get idpHitsInfo;

  /// No description provided for @idpExpect.
  ///
  /// In en, this message translates to:
  /// **'E-value'**
  String get idpExpect;

  /// No description provided for @idpExpectInfo.
  ///
  /// In en, this message translates to:
  /// **'Highest expect value kept (EXPECT): lower is stricter.'**
  String get idpExpectInfo;

  /// No description provided for @idpWord.
  ///
  /// In en, this message translates to:
  /// **'Word size'**
  String get idpWord;

  /// No description provided for @idpWordInfo.
  ///
  /// In en, this message translates to:
  /// **'Length of the seed matches (WORD_SIZE). Empty: the program\'s default (28 for megablast, 11 for blastn).'**
  String get idpWordInfo;

  /// No description provided for @idpCoverage.
  ///
  /// In en, this message translates to:
  /// **'Min coverage (%)'**
  String get idpCoverage;

  /// No description provided for @idpCoverageInfo.
  ///
  /// In en, this message translates to:
  /// **'Share of your sequence a hit must cover to count (applied by SPLACE; changing it needs no new search).'**
  String get idpCoverageInfo;

  /// No description provided for @idpEntrez.
  ///
  /// In en, this message translates to:
  /// **'Restrict to (Entrez query)'**
  String get idpEntrez;

  /// No description provided for @idpEntrezInfo.
  ///
  /// In en, this message translates to:
  /// **'Optional: search only records matching an Entrez query, e.g. Bufonidae[Organism] or Anura[Organism] (ENTREZ_QUERY).'**
  String get idpEntrezInfo;

  /// No description provided for @idpBoldDbInfo.
  ///
  /// In en, this message translates to:
  /// **'public.tax-derep: public animal COI, one record per haplotype and taxon (recommended). species: records identified to species. all.tax-derep: public and private records (taxonomy only). all.animal-alt: animal records of other markers. public.plants: rbcL and matK. public.fungi: ITS.'**
  String get idpBoldDbInfo;

  /// No description provided for @idpMinIdentity.
  ///
  /// In en, this message translates to:
  /// **'Min identity (%)'**
  String get idpMinIdentity;

  /// No description provided for @idpMinIdentityInfo.
  ///
  /// In en, this message translates to:
  /// **'Hits below this identity are not returned (70–100 %).'**
  String get idpMinIdentityInfo;

  /// No description provided for @idpMinOverlap.
  ///
  /// In en, this message translates to:
  /// **'Min overlap (bp)'**
  String get idpMinOverlap;

  /// No description provided for @idpMinOverlapInfo.
  ///
  /// In en, this message translates to:
  /// **'Shortest alignment with your sequence that counts.'**
  String get idpMinOverlapInfo;

  /// No description provided for @idpMaxHits.
  ///
  /// In en, this message translates to:
  /// **'Max hits'**
  String get idpMaxHits;

  /// No description provided for @idpMaxHitsInfo.
  ///
  /// In en, this message translates to:
  /// **'How many hits BOLD returns per sequence.'**
  String get idpMaxHitsInfo;

  /// No description provided for @idpBoldMarkers.
  ///
  /// In en, this message translates to:
  /// **'This database identifies: {markers}.'**
  String idpBoldMarkers(String markers);

  /// No description provided for @speciesIdBoldBodyFor.
  ///
  /// In en, this message translates to:
  /// **'ID engine, database {db}: {markers}.'**
  String speciesIdBoldBodyFor(String db, String markers);

  /// No description provided for @speciesIdBoldOffFor.
  ///
  /// In en, this message translates to:
  /// **'Select {markers} to search BOLD too (the database can be changed in the parameters).'**
  String speciesIdBoldOffFor(String markers);

  /// No description provided for @exampleLociShort.
  ///
  /// In en, this message translates to:
  /// **'Loci example'**
  String get exampleLociShort;

  /// No description provided for @exampleBarcodingShort.
  ///
  /// In en, this message translates to:
  /// **'Barcoding example'**
  String get exampleBarcodingShort;

  /// No description provided for @speciesIdGenesHintFor.
  ///
  /// In en, this message translates to:
  /// **'Choose one or more genes. Every gene is searched in GenBank with BLAST; BOLD takes the markers of its database ({markers}).'**
  String speciesIdGenesHintFor(String markers);

  /// No description provided for @speciesIdGenbankBodyFor.
  ///
  /// In en, this message translates to:
  /// **'BLAST ({program}) against {db}, up to {hits} results per gene.'**
  String speciesIdGenbankBodyFor(String program, String db, String hits);

  /// No description provided for @trimPlanCodonsTitle.
  ///
  /// In en, this message translates to:
  /// **'Whole codons (reading frame known)'**
  String get trimPlanCodonsTitle;

  /// No description provided for @trimPlanColumnsTitle.
  ///
  /// In en, this message translates to:
  /// **'Column by column — {why}'**
  String trimPlanColumnsTitle(String why);

  /// No description provided for @openProject.
  ///
  /// In en, this message translates to:
  /// **'Open project'**
  String get openProject;

  /// No description provided for @saveProject.
  ///
  /// In en, this message translates to:
  /// **'Save project'**
  String get saveProject;

  /// No description provided for @saveProjectTooltip.
  ///
  /// In en, this message translates to:
  /// **'Save project (.splc): data, parameters and results, to continue later or share (Ctrl+S; Save as: Ctrl+Shift+S)'**
  String get saveProjectTooltip;

  /// No description provided for @openProjectTooltip.
  ///
  /// In en, this message translates to:
  /// **'Open a SPLACE project (.splc) (Ctrl+O)'**
  String get openProjectTooltip;

  /// No description provided for @openProjectTitle.
  ///
  /// In en, this message translates to:
  /// **'Open another project?'**
  String get openProjectTitle;

  /// No description provided for @openProjectBody.
  ///
  /// In en, this message translates to:
  /// **'{file} opens where it was saved. The current analysis is closed: save it first (Save project) to come back to it.'**
  String openProjectBody(String file);

  /// No description provided for @projectSaveFailed.
  ///
  /// In en, this message translates to:
  /// **'The project could not be saved: see the log.'**
  String get projectSaveFailed;

  /// No description provided for @splashOpenProject.
  ///
  /// In en, this message translates to:
  /// **'Open project'**
  String get splashOpenProject;

  /// No description provided for @windowTitleNew.
  ///
  /// In en, this message translates to:
  /// **'Unsaved project'**
  String get windowTitleNew;

  /// No description provided for @windowTitleUnsaved.
  ///
  /// In en, this message translates to:
  /// **'unsaved changes'**
  String get windowTitleUnsaved;

  /// No description provided for @windowTitleSaved.
  ///
  /// In en, this message translates to:
  /// **'saved'**
  String get windowTitleSaved;

  /// No description provided for @projectSavedTitle.
  ///
  /// In en, this message translates to:
  /// **'Project saved'**
  String get projectSavedTitle;

  /// No description provided for @exportGbDoneTitle.
  ///
  /// In en, this message translates to:
  /// **'{n} GenBank file(s) saved'**
  String exportGbDoneTitle(String n);

  /// No description provided for @delimSettings.
  ///
  /// In en, this message translates to:
  /// **'Settings'**
  String get delimSettings;

  /// No description provided for @delimLrtLevel.
  ///
  /// In en, this message translates to:
  /// **'LRT p-value · {tool}'**
  String delimLrtLevel(String tool);

  /// No description provided for @delimLrtInfoMptp.
  ///
  /// In en, this message translates to:
  /// **'Significance level of the likelihood-ratio test of mPTP (--pvalue, default 0.001). When the p-value of the test is above it, the one-species model is kept and every tip goes into one species.'**
  String get delimLrtInfoMptp;

  /// No description provided for @delimLrtInfoGmyc.
  ///
  /// In en, this message translates to:
  /// **'Significance level of the GMYC likelihood-ratio test (default 0.01). Above it GMYC does not reject a single species and joins every tip.'**
  String get delimLrtInfoGmyc;

  /// No description provided for @delimOnePerSpecies.
  ///
  /// In en, this message translates to:
  /// **'Every tip of the supermatrix is a different species (one sequence per species). PTP, mPTP and GMYC look for where the tree changes from branching between species to branching within them; without several individuals per species that change is not there, the test tends to fail and the whole tree comes out as one species. Use several sequences per species, or ASAP/ABGD on a marker.'**
  String get delimOnePerSpecies;

  /// No description provided for @delimCollapsed.
  ///
  /// In en, this message translates to:
  /// **'The likelihood-ratio test did not reject the one-species model (p = {p}, above the limit {alpha}), so {tool} put all {tips} tips in one species.'**
  String delimCollapsed(String tool, String p, String alpha, String tips);

  /// No description provided for @delimCollapsedBest.
  ///
  /// In en, this message translates to:
  /// **'Its best model had found {n} species, without statistical support.'**
  String delimCollapsedBest(String n);

  /// No description provided for @delimCollapsedOnePerSpecies.
  ///
  /// In en, this message translates to:
  /// **'Here every tip is a different species: these methods need several individuals per species. With one sequence each, use ASAP or ABGD on a marker, or compare the tree with the taxonomy.'**
  String get delimCollapsedOnePerSpecies;

  /// No description provided for @delimCollapsedTry.
  ///
  /// In en, this message translates to:
  /// **'Try a less strict limit (LRT p-value in the settings above), more sequences per species, or ASAP/ABGD on a marker.'**
  String get delimCollapsedTry;

  /// No description provided for @delimLrtNs.
  ///
  /// In en, this message translates to:
  /// **'LRT not significant'**
  String get delimLrtNs;

  /// No description provided for @treeIdentify.
  ///
  /// In en, this message translates to:
  /// **'Identify by BLAST…'**
  String get treeIdentify;

  /// No description provided for @treeIdTitle.
  ///
  /// In en, this message translates to:
  /// **'Identify by BLAST'**
  String get treeIdTitle;

  /// No description provided for @treeIdSubtitle.
  ///
  /// In en, this message translates to:
  /// **'{n} tip(s) of the tree'**
  String treeIdSubtitle(String n);

  /// No description provided for @treeIdIntro.
  ///
  /// In en, this message translates to:
  /// **'Each tip is sent with its sequence of the marker you choose to GenBank BLAST (and BOLD for COI, or rbcL, matK and ITS in its plant and fungal databases), with the options of the identification in step 2. Useful when a genus falls inside another one: the best match tells whether the sample is what its name says.'**
  String get treeIdIntro;

  /// No description provided for @treeIdMarker.
  ///
  /// In en, this message translates to:
  /// **'Marker'**
  String get treeIdMarker;

  /// No description provided for @treeIdMarkerTips.
  ///
  /// In en, this message translates to:
  /// **'{n} of {total} tips have it'**
  String treeIdMarkerTips(String n, String total);

  /// No description provided for @treeIdNoMarker.
  ///
  /// In en, this message translates to:
  /// **'None of these tips has a DNA marker sequence of at least 100 bases.'**
  String get treeIdNoMarker;

  /// No description provided for @treeIdOptions.
  ///
  /// In en, this message translates to:
  /// **'BLAST: {program} on {db}, hits covering at least {coverage} % of the sequence; same species from {identity} % identity. Answers already received are reused.'**
  String treeIdOptions(
    String program,
    String db,
    String coverage,
    String identity,
  );

  /// No description provided for @treeIdRun.
  ///
  /// In en, this message translates to:
  /// **'Identify'**
  String get treeIdRun;

  /// No description provided for @treeIdRunning.
  ///
  /// In en, this message translates to:
  /// **'Identifying {n} tip(s) with {marker}…'**
  String treeIdRunning(String marker, String n);

  /// No description provided for @treeIdOtherMarker.
  ///
  /// In en, this message translates to:
  /// **'Another marker'**
  String get treeIdOtherMarker;

  /// No description provided for @treeIdRename.
  ///
  /// In en, this message translates to:
  /// **'Rename {n} tip(s) in the figure'**
  String treeIdRename(String n);

  /// No description provided for @treeIdRenameOne.
  ///
  /// In en, this message translates to:
  /// **'Rename in the figure as {name}'**
  String treeIdRenameOne(String name);

  /// No description provided for @treeIdConfirmed.
  ///
  /// In en, this message translates to:
  /// **'Confirmed'**
  String get treeIdConfirmed;

  /// No description provided for @treeIdOtherSpecies.
  ///
  /// In en, this message translates to:
  /// **'Another species'**
  String get treeIdOtherSpecies;

  /// No description provided for @treeIdOtherGenus.
  ///
  /// In en, this message translates to:
  /// **'Another genus'**
  String get treeIdOtherGenus;

  /// No description provided for @treeIdAmbiguous.
  ///
  /// In en, this message translates to:
  /// **'Ambiguous (its species matches as well)'**
  String get treeIdAmbiguous;

  /// No description provided for @treeIdNoMatch.
  ///
  /// In en, this message translates to:
  /// **'No species at ≥ {identity} %'**
  String treeIdNoMatch(String identity);

  /// No description provided for @treeIdNone.
  ///
  /// In en, this message translates to:
  /// **'No result'**
  String get treeIdNone;

  /// No description provided for @treeIdCurrent.
  ///
  /// In en, this message translates to:
  /// **'Name:'**
  String get treeIdCurrent;

  /// No description provided for @treeIdCurrentAt.
  ///
  /// In en, this message translates to:
  /// **'its name at {pct}'**
  String treeIdCurrentAt(String pct);

  /// No description provided for @treeIdHits.
  ///
  /// In en, this message translates to:
  /// **'Alignments'**
  String get treeIdHits;
}

class _AppLocalizationsDelegate
    extends LocalizationsDelegate<AppLocalizations> {
  const _AppLocalizationsDelegate();

  @override
  Future<AppLocalizations> load(Locale locale) {
    return SynchronousFuture<AppLocalizations>(lookupAppLocalizations(locale));
  }

  @override
  bool isSupported(Locale locale) =>
      <String>['en', 'es', 'pt'].contains(locale.languageCode);

  @override
  bool shouldReload(_AppLocalizationsDelegate old) => false;
}

AppLocalizations lookupAppLocalizations(Locale locale) {
  // Lookup logic when only language code is specified.
  switch (locale.languageCode) {
    case 'en':
      return AppLocalizationsEn();
    case 'es':
      return AppLocalizationsEs();
    case 'pt':
      return AppLocalizationsPt();
  }

  throw FlutterError(
    'AppLocalizations.delegate failed to load unsupported locale "$locale". This is likely '
    'an issue with the localizations generation tool. Please file an issue '
    'on GitHub with a reproducible sample app and the gen-l10n configuration '
    'that was used.',
  );
}
