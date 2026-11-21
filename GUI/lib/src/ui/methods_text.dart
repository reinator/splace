import '../core/delimitation.dart';
import '../core/genetic_code.dart';
import '../services/tool_locator.dart' show ExternalTool;
import '../core/tool_commands.dart';
import '../core/saturation.dart';
import '../core/supermatrix.dart';
import '../core/syngenes.dart';
import '../core/taxonomy_sources.dart';
import '../state/project_controller.dart';
import 'common.dart';
import 'widgets/citations.dart';

/// Methods paragraphs written from what was done in the project, in
/// academic language, with the tools cited (author–year) and the list of
/// the works cited.
class MethodsText {
  MethodsText(this.paragraphs, this.cited, [this.roles = const {}]);
  final List<String> paragraphs;

  /// What each cited tool did in this project, by citation name.
  final Map<String, String> roles;

  /// Names of [citations] cited in the text, in order of first citation.
  final List<String> cited;

  bool get isEmpty => paragraphs.isEmpty;

  /// Paragraphs and the references, as plain text.
  String get plain {
    final refs = [for (final n in cited) citations.firstWhere((c) => c.name == n).text];
    return [...paragraphs, '', ...refs].join('\n\n').replaceAll('\n\n\n\n', '\n\n');
  }
}

/// Tools of our group highlighted in "How to cite" (name in bold, the
/// repository with the first citation).
const highlightedTools = ['SynGenes', 'dataFishing'];

/// Builds the methods text of [c] in [lang] (`en`, `pt` or `es`).
MethodsText buildMethodsText(ProjectController c, String lang) => _Writer(c, lang).build();

class _Writer {
  _Writer(this.c, String lang) : lang = const {'pt', 'es'}.contains(lang) ? lang : 'en';
  final ProjectController c;
  final String lang;
  final cited = <String>[];
  final roles = <String, String>{};

  /// Records what [name] did in this project (shown with its reference).
  void role(String name, String text) {
    final old = roles[name];
    if (old == null) {
      roles[name] = text;
    } else if (!old.contains(text)) {
      roles[name] = '$old $text';
    }
  }

  String tr(String en, String pt, String es) => switch (lang) { 'pt' => pt, 'es' => es, _ => en };

  /// "(Katoh & Standley 2013)" or "(A 2017; B 2016)".
  String cite(List<String> names) {
    final parts = <String>[];
    for (final n in names) {
      final ref = citations.firstWhere((x) => x.name == n);
      // SynGenes and dataFishing: their repository with the first
      // citation, so readers find the code.
      final first = !cited.contains(n);
      if (first) cited.add(n);
      parts.add('${shortAuthors(ref.authors)} ${ref.year}${first && highlightedTools.contains(n) && ref.url != null ? '; ${ref.url}' : ''}');
    }
    return '(${parts.join('; ')})';
  }

  /// "a, b and c" in the language of the text.
  String list(List<String> items) {
    if (items.length <= 1) return items.join();
    final and = tr('and', 'e', 'y');
    return '${items.sublist(0, items.length - 1).join(', ')} $and ${items.last}';
  }

  String plural(int n, String one, String many) => n == 1 ? one : many;

  /// Version of [tool] found on this computer ("7.526"), or null.
  String? version(ExternalTool tool) => c.tools[tool]?.version;

  /// "MAFFT v7.526" (or "MAFFT v[fallback]" when unknown).
  String named(String name, ExternalTool tool, [String fallback = '']) {
    final v = version(tool) ?? fallback;
    return v.isEmpty ? name : '$name v$v';
  }

  /// IQ-TREE 3 or 2, by the version that ran.
  bool get iqtree3 => !(version(ExternalTool.iqtree)?.startsWith('2') ?? false);
  String get iqtreeRef => iqtree3 ? 'IQ-TREE 3' : 'IQ-TREE 2';

  String ordinal(int p) => switch (lang) {
        'pt' => '$pª',
        'es' => '$pª',
        _ => switch (p) { 1 => '1st', 2 => '2nd', 3 => '3rd', _ => '${p}th' },
      };

  MethodsText build() {
    final out = <String>[
      ?_data(),
      ?_alignment(),
      ?_trimming(),
      ?_saturation(),
      ?_matrix(),
      ..._inference(),
      ?_delimitation(),
    ];
    if (out.isNotEmpty) {
      role('SPLACE', tr(
        'Ran the whole workflow: data retrieval, gene extraction, alignment, trimming, saturation test, supermatrix, tree inference and figures (v$appVersion).',
        'Executou todo o fluxo: obtenção dos dados, extração dos genes, alinhamento, trimagem, teste de saturação, supermatriz, inferência das árvores e figuras (v$appVersion).',
        'Ejecutó todo el flujo: obtención de los datos, extracción de los genes, alineamiento, recorte, prueba de saturación, supermatriz, inferencia de los árboles y figuras (v$appVersion).',
      ));
    }
    return MethodsText(out, cited, roles);
  }

  // -------------------------------------------------------------------------

  /// Individual markers searched with queries generated from SynGenes.
  String _markerSearch(int gb) {
    final searches = c.markerSearches.map((m) => m.search).toList();
    final bold = c.records.where((r) => r.source == 'bold').length;
    final fromNcbi = gb - bold;
    final genes = <String>{for (final m in searches) ...m.genes}.toList();
    final taxa = <String>{for (final m in searches) m.taxon.trim()}.toList();
    final species = {for (final r in c.records) r.organismName}.where((x) => x.isNotEmpty).length;
    final mt = searches.first.organelle == OrganelleType.mt;
    final short = searches.every((m) => !m.includeGenomes);
    final refseq = searches.every((m) => m.refseqOnly);
    final maxKb = searches.first.maxLength ~/ 1000;
    final out = StringBuffer(tr(
      'Sequences of the ${plural(genes.length, 'marker', 'markers')} ${list(genes)} of ${list(taxa)} were retrieved from the NCBI nucleotide database '
          '($fromNcbi GenBank ${plural(fromNcbi, 'record', 'records')}, $species ${plural(species, 'species', 'species')}) with Entrez queries generated by SynGenes ${cite(['SynGenes'])}, '
          'which gathers the synonymous names of each ${mt ? 'mitochondrial' : 'chloroplast'} gene into a single search (gene symbols in the gene field and product names in the title), '
          'restricted to ${mt ? 'mitochondrial' : 'chloroplast'} genomic DNA${short ? ' and to records of up to $maxKb kb' : ''}${refseq ? ' from RefSeq' : ''}; '
          'SynGenes also standardized the gene annotations of the records, and the queries used are saved with the project (00_ncbi/marker_queries.txt). ',
      'As sequências ${plural(genes.length, 'do marcador', 'dos marcadores')} ${list(genes)} de ${list(taxa)} foram obtidas na base de nucleotídeos do NCBI '
          '($fromNcbi ${plural(fromNcbi, 'registro', 'registros')} do GenBank, $species ${plural(species, 'espécie', 'espécies')}) com buscas Entrez geradas pelo SynGenes ${cite(['SynGenes'])}, '
          'que reúne os nomes sinônimos de cada gene ${mt ? 'mitocondrial' : 'de cloroplasto'} em uma única busca (símbolos do gene no campo gene e nomes do produto no título), '
          'restritas ao DNA genômico ${mt ? 'mitocondrial' : 'de cloroplasto'}${short ? ' e a registros de até $maxKb kb' : ''}${refseq ? ' do RefSeq' : ''}; '
          'o SynGenes também padronizou as anotações dos genes dos registros, e as buscas usadas ficam salvas no projeto (00_ncbi/marker_queries.txt). ',
      'Las secuencias ${plural(genes.length, 'del marcador', 'de los marcadores')} ${list(genes)} de ${list(taxa)} se obtuvieron de la base de nucleótidos del NCBI '
          '($fromNcbi ${plural(fromNcbi, 'registro', 'registros')} de GenBank, $species ${plural(species, 'especie', 'especies')}) con búsquedas Entrez generadas por SynGenes ${cite(['SynGenes'])}, '
          'que reúne los nombres sinónimos de cada gen ${mt ? 'mitocondrial' : 'de cloroplasto'} en una sola búsqueda (símbolos del gen en el campo gene y nombres del producto en el título), '
          'restringidas al ADN genómico ${mt ? 'mitocondrial' : 'de cloroplasto'}${short ? ' y a registros de hasta $maxKb kb' : ''}${refseq ? ' de RefSeq' : ''}; '
          'SynGenes también estandarizó las anotaciones de los genes de los registros, y las búsquedas usadas se guardan en el proyecto (00_ncbi/marker_queries.txt). ',
    ));
    role('SynGenes', tr(
      'Generated the NCBI queries of ${list(genes)} (every synonymous name of each gene in a single search) and standardized the gene names of the $gb records.',
      'Gerou as buscas no NCBI de ${list(genes)} (todos os nomes sinônimos de cada gene em uma única busca) e padronizou os nomes dos genes dos $gb registros.',
      'Generó las búsquedas en el NCBI de ${list(genes)} (todos los nombres sinónimos de cada gen en una sola búsqueda) y estandarizó los nombres de los genes de los $gb registros.',
    ));
    // Records found only in BOLD Systems.
    if (bold > 0) {
      final codes = <String>{for (final m in searches) ...m.boldCodes}.toList();
      out.write(tr(
        'The public records of BOLD Systems ${cite(['BOLD'])} were searched for the same taxa (${list(codes)}, through the BOLD v5 portal API); records that BOLD mirrors from GenBank were taken from GenBank, '
            'and $bold ${plural(bold, 'record found only in BOLD was', 'records found only in BOLD were')} added, each annotated with its BOLD process id, BIN and reading frame. ',
        'Os registros públicos do BOLD Systems ${cite(['BOLD'])} foram buscados para os mesmos táxons (${list(codes)}, pela API do portal BOLD v5); os registros que o BOLD espelha do GenBank foram obtidos do GenBank, '
            'e $bold ${plural(bold, 'registro encontrado só no BOLD foi acrescentado', 'registros encontrados só no BOLD foram acrescentados')}, cada um anotado com o process id, o BIN e a fase de leitura. ',
        'Los registros públicos de BOLD Systems ${cite(['BOLD'])} se buscaron para los mismos taxones (${list(codes)}, mediante la API del portal BOLD v5); los registros que BOLD replica de GenBank se tomaron de GenBank, '
            'y se añadieron $bold ${plural(bold, 'registro encontrado solo en BOLD', 'registros encontrados solo en BOLD')}, cada uno anotado con su process id, BIN y marco de lectura. ',
      ));
      role('BOLD', tr(
        'Searched for ${list(codes)} with NCBI; $bold ${plural(bold, 'record', 'records')} came only from BOLD.',
        'Buscado para ${list(codes)} junto com o NCBI; $bold ${plural(bold, 'registro veio', 'registros vieram')} só do BOLD.',
        'Buscado para ${list(codes)} junto con el NCBI; $bold ${plural(bold, 'registro vino', 'registros vinieron')} solo de BOLD.',
      ));
    }
    if (c.onePerTaxonMode) {
      out.write(tr(
        'As the markers come from separate records, the sequences of each species were joined into a single terminal, keeping the longest sequence when a species had more than one record of a marker; '
            'markers missing from a species were coded as missing data (?). ',
        'Como os marcadores vêm de registros separados, as sequências de cada espécie foram reunidas em um único terminal, mantendo a sequência mais longa quando uma espécie tinha mais de um registro do mesmo marcador; '
            'os marcadores ausentes em uma espécie foram codificados como dados faltantes (?). ',
        'Como los marcadores provienen de registros separados, las secuencias de cada especie se unieron en un único terminal, conservando la secuencia más larga cuando una especie tenía más de un registro del mismo marcador; '
            'los marcadores ausentes en una especie se codificaron como datos faltantes (?). ',
      ));
    }
    return out.toString();
  }

  String? _data() {
    final gb = c.records.length;
    final fastaFiles = c.inputs.length;
    final fastaSamples = c.organized?.samples.length ?? 0;
    if (gb == 0 && fastaFiles == 0) return null;
    final s = StringBuffer(tr(
      'Sequence data were obtained and processed in SPLACE v$appVersion ${cite(['SPLACE'])}. ',
      'Os dados de sequências foram obtidos e processados no SPLACE v$appVersion ${cite(['SPLACE'])}. ',
      'Los datos de secuencias se obtuvieron y procesaron en SPLACE v$appVersion ${cite(['SPLACE'])}. ',
    ));
    if (gb > 0 && c.markerSearches.isNotEmpty) {
      s.write(_markerSearch(gb));
    } else if (gb > 0) {
      final species = {for (final r in c.records) r.organismName}.where((x) => x.isNotEmpty).length;
      final mt = c.organelle == OrganelleType.mt;
      final genomes = mt ? tr('mitochondrial', 'mitocondriais', 'mitocondriales') : tr('chloroplast', 'de cloroplasto', 'de cloroplasto');
      s.write(tr(
        'A total of $gb $genomes GenBank ${plural(gb, 'record', 'records')} ($species ${plural(species, 'species', 'species')}) '
            'were retrieved from the NCBI nucleotide database, and their gene annotations were standardized with SynGenes ${cite(['SynGenes'])}, '
            'which unifies the synonymous names given to the same ${mt ? 'mitochondrial' : 'chloroplast'} gene. ',
        'Foram obtidos $gb ${plural(gb, 'registro', 'registros')} $genomes do GenBank ($species ${plural(species, 'espécie', 'espécies')}) '
            'na base de nucleotídeos do NCBI, e as anotações dos genes foram padronizadas com o SynGenes ${cite(['SynGenes'])}, '
            'que unifica os nomes sinônimos dados a um mesmo gene ${mt ? 'mitocondrial' : 'de cloroplasto'}. ',
        'Se obtuvieron $gb ${plural(gb, 'registro', 'registros')} $genomes de GenBank ($species ${plural(species, 'especie', 'especies')}) '
            'de la base de nucleótidos del NCBI, y las anotaciones de los genes se estandarizaron con SynGenes ${cite(['SynGenes'])}, '
            'que unifica los nombres sinónimos dados a un mismo gen ${mt ? 'mitocondrial' : 'de cloroplasto'}. ',
      ));
      role('SynGenes', tr(
        'Standardized the gene names of the $gb GenBank records, so that each gene has the same name in every genome.',
        'Padronizou os nomes dos genes dos $gb registros do GenBank, para que cada gene tenha o mesmo nome em todos os genomas.',
        'Estandarizó los nombres de los genes de los $gb registros de GenBank, para que cada gen tenga el mismo nombre en todos los genomas.',
      ));
    }
    if (fastaFiles > 0) {
      s.write(tr(
        '${gb > 0 ? 'In addition, ' : ''}$fastaFiles FASTA ${plural(fastaFiles, 'file', 'files')}'
            '${fastaSamples > 0 ? ' with $fastaSamples ${plural(fastaSamples, 'sample', 'samples')}' : ''} ${plural(fastaFiles, 'was', 'were')} provided by the user, '
            'and their locus names were matched to the standard gene names with SynGenes${gb > 0 ? '' : ' ${cite(['SynGenes'])}'}. ',
        '${gb > 0 ? 'Além disso, ' : ''}$fastaFiles ${plural(fastaFiles, 'arquivo', 'arquivos')} FASTA'
            '${fastaSamples > 0 ? ' com $fastaSamples ${plural(fastaSamples, 'amostra', 'amostras')}' : ''} ${plural(fastaFiles, 'foi fornecido', 'foram fornecidos')} pelo usuário, '
            'e os nomes dos loci foram associados aos nomes padronizados dos genes com o SynGenes${gb > 0 ? '' : ' ${cite(['SynGenes'])}'}. ',
        '${gb > 0 ? 'Además, ' : ''}$fastaFiles ${plural(fastaFiles, 'archivo', 'archivos')} FASTA'
            '${fastaSamples > 0 ? ' con $fastaSamples ${plural(fastaSamples, 'muestra', 'muestras')}' : ''} ${plural(fastaFiles, 'fue proporcionado', 'fueron proporcionados')} por el usuario, '
            'y los nombres de los loci se asociaron a los nombres estandarizados de los genes con SynGenes${gb > 0 ? '' : ' ${cite(['SynGenes'])}'}. ',
      ));
      role('SynGenes', tr(
        'Matched the locus names of the $fastaFiles FASTA ${plural(fastaFiles, 'file', 'files')} to the standard gene names.',
        'Associou os nomes dos loci ${plural(fastaFiles, 'do arquivo', 'dos $fastaFiles arquivos')} FASTA aos nomes padronizados dos genes.',
        'Asoció los nombres de los loci ${plural(fastaFiles, 'del archivo', 'de los $fastaFiles archivos')} FASTA a los nombres estandarizados de los genes.',
      ));
    }
    if (c.taxonomyResolved && gb > 0) {
      final source = switch (c.taxonomySource) {
        TaxonomySource.gbif => 'GBIF Backbone Taxonomy',
        TaxonomySource.ncbi => 'NCBI Taxonomy',
        TaxonomySource.salve => tr('SALVE (ICMBio, Brazil)', 'SALVE (ICMBio)', 'SALVE (ICMBio, Brasil)'),
      };
      s.write(tr(
        'The taxonomic classification of each record (kingdom to species) was retrieved from $source following the dataFishing approach ${cite(['dataFishing'])}. ',
        'A classificação taxonômica de cada registro (de reino a espécie) foi obtida do $source seguindo a abordagem do dataFishing ${cite(['dataFishing'])}. ',
        'La clasificación taxonómica de cada registro (de reino a especie) se obtuvo de $source siguiendo el enfoque de dataFishing ${cite(['dataFishing'])}. ',
      ));
      final species = {for (final r in c.records) r.organismName}.where((x) => x.isNotEmpty).length;
      role('dataFishing', tr(
        'Retrieved the taxonomic classification (kingdom to species) of the $species species from $source, used in the FASTA headers, the records table and the tree labels and colours.',
        'Obteve a classificação taxonômica (de reino a espécie) das $species espécies no $source, usada nos cabeçalhos FASTA, na tabela de registros e nos rótulos e cores da árvore.',
        'Obtuvo la clasificación taxonómica (de reino a especie) de las $species especies en $source, usada en los encabezados FASTA, la tabla de registros y las etiquetas y colores del árbol.',
      ));
    }
    if (c.speciesIds.isNotEmpty) {
      s.write(tr(
        'The identity of the sequences was checked by similarity searches against GenBank (BLAST) and BOLD Systems. ',
        'A identidade das sequências foi conferida por buscas de similaridade no GenBank (BLAST) e no BOLD Systems. ',
        'La identidad de las secuencias se verificó mediante búsquedas de similitud en GenBank (BLAST) y BOLD Systems. ',
      ));
    }
    return s.toString().trim();
  }

  List<LocusState> get _aligned => [for (final l in c.loci.values) if (l.included && l.aligned != null) l];

  String _strategy(String s) {
    final f = s.trim();
    final name = switch (f) {
      '--auto' => tr('automatic strategy selection', 'seleção automática de estratégia', 'selección automática de estrategia'),
      '--localpair' || '--localpair --maxiterate 1000' => 'L-INS-i',
      '--globalpair' || '--globalpair --maxiterate 1000' => 'G-INS-i',
      '--genafpair' || '--genafpair --maxiterate 1000' => 'E-INS-i',
      '--retree 2' => 'FFT-NS-2',
      '--retree 1' => 'FFT-NS-1',
      '--retree 2 --maxiterate 2' || '--maxiterate 2' => 'FFT-NS-i',
      _ => null,
    };
    return name == null ? f : '$name ($f)';
  }

  String? _alignment() {
    final loci = _aligned;
    if (loci.isEmpty) return null;
    final names = [for (final l in loci) l.name];
    final s = StringBuffer(tr(
      'A total of ${loci.length} ${plural(loci.length, 'marker was', 'markers were')} analysed (${list(names)}). ',
      '${plural(loci.length, 'Foi analisado ${loci.length} marcador', 'Foram analisados ${loci.length} marcadores')} (${list(names)}). ',
      '${plural(loci.length, 'Se analizó ${loci.length} marcador', 'Se analizaron ${loci.length} marcadores')} (${list(names)}). ',
    ));
    final m = c.mafft;
    final strategy = _strategy(m.strategy);
    final extras = [
      if (m.maxIterate > 0 && !m.strategy.contains('--maxiterate'))
        tr('up to ${m.maxIterate} refinement iterations', 'até ${m.maxIterate} iterações de refinamento', 'hasta ${m.maxIterate} iteraciones de refinamiento'),
      if (m.adjustDirection && c.seqMode == SeqMode.nucleotide)
        tr('automatic correction of reverse-complemented sequences (--adjustdirection)', 'correção automática de sequências em sentido reverso (--adjustdirection)',
            'corrección automática de secuencias en sentido inverso (--adjustdirection)'),
    ];
    final withExtras = extras.isEmpty ? '' : tr(', with ${list(extras)}', ', com ${list(extras)}', ', con ${list(extras)}');
    final name = geneticCodeName(c.geneticCode);
    final code = name.contains('NCBI') ? name : '$name (${tr('NCBI table', 'tabela NCBI', 'tabla NCBI')} ${c.geneticCode})';
    final coding = loci.where((l) => l.coding && !l.isProtein).length;
    final nonCoding = loci.length - coding;
    final mafft = cite(['MAFFT']);
    role('MAFFT', tr(
      'Aligned the ${loci.length} ${plural(loci.length, 'marker', 'markers')} ($strategy${c.seqMode == SeqMode.codon ? ', codon-aware: proteins aligned and back-translated' : c.seqMode == SeqMode.aminoAcid ? ', as amino acids' : ''}).',
      'Alinhou ${plural(loci.length, 'o marcador', 'os ${loci.length} marcadores')} ($strategy${c.seqMode == SeqMode.codon ? ', preservando os códons: proteínas alinhadas e retrotraduzidas' : c.seqMode == SeqMode.aminoAcid ? ', como aminoácidos' : ''}).',
      'Alineó ${plural(loci.length, 'el marcador', 'los ${loci.length} marcadores')} ($strategy${c.seqMode == SeqMode.codon ? ', preservando los codones: proteínas alineadas y retrotraducidas' : c.seqMode == SeqMode.aminoAcid ? ', como aminoácidos' : ''}).',
    ));
    switch (c.seqMode) {
      case SeqMode.codon:
        s.write(tr(
          'Protein-coding genes were aligned in a codon-aware way: their sequences were translated with the $code genetic code, aligned as amino acids '
              'with ${named('MAFFT', ExternalTool.mafft, '7')} $mafft using $strategy$withExtras, and back-translated to nucleotides, so that codons were never split by gaps. ',
          'Os genes codificadores de proteínas foram alinhados preservando os códons: as sequências foram traduzidas com o código genético $code, alinhadas '
              'como aminoácidos no ${named('MAFFT', ExternalTool.mafft, '7')} $mafft com $strategy$withExtras e retrotraduzidas para nucleotídeos, de modo que nenhum códon fosse partido por gaps. ',
          'Los genes codificantes de proteínas se alinearon preservando los codones: las secuencias se tradujeron con el código genético $code, se alinearon '
              'como aminoácidos en ${named('MAFFT', ExternalTool.mafft, '7')} $mafft con $strategy$withExtras y se retrotradujeron a nucleótidos, de modo que ningún codón quedara dividido por gaps. ',
        ));
        if (nonCoding > 0) {
          s.write(tr(
            'Non-coding markers ($nonCoding) were aligned as nucleotides with the same settings. ',
            'Os marcadores não codificadores ($nonCoding) foram alinhados como nucleotídeos com as mesmas configurações. ',
            'Los marcadores no codificantes ($nonCoding) se alinearon como nucleótidos con la misma configuración. ',
          ));
        }
      case SeqMode.aminoAcid:
        s.write(tr(
          'Protein-coding genes were translated with the $code genetic code and aligned as amino-acid sequences with ${named('MAFFT', ExternalTool.mafft, '7')} $mafft using $strategy$withExtras. ',
          'Os genes codificadores de proteínas foram traduzidos com o código genético $code e alinhados como sequências de aminoácidos no ${named('MAFFT', ExternalTool.mafft, '7')} $mafft com $strategy$withExtras. ',
          'Los genes codificantes de proteínas se tradujeron con el código genético $code y se alinearon como secuencias de aminoácidos en ${named('MAFFT', ExternalTool.mafft, '7')} $mafft con $strategy$withExtras. ',
        ));
      case SeqMode.nucleotide:
        s.write(tr(
          'Each marker was aligned separately as nucleotides with ${named('MAFFT', ExternalTool.mafft, '7')} $mafft using $strategy$withExtras. ',
          'Cada marcador foi alinhado separadamente como nucleotídeos no ${named('MAFFT', ExternalTool.mafft, '7')} $mafft com $strategy$withExtras. ',
          'Cada marcador se alineó por separado como nucleótidos en ${named('MAFFT', ExternalTool.mafft, '7')} $mafft con $strategy$withExtras. ',
        ));
    }
    return s.toString().trim();
  }

  String? _trimming() {
    final loci = _aligned;
    if (loci.isEmpty) return null;
    final trimmed = [for (final l in loci) if (l.trimmed != null) l];
    if (trimmed.isEmpty || !c.useTrimmed) {
      return c.trimSkipped || trimmed.isEmpty
          ? tr('Alignments were used without trimming.', 'Os alinhamentos foram usados sem trimagem.', 'Los alineamientos se usaron sin recorte.')
          : tr('Trimmed alignments were produced, but the untrimmed alignments were used in the analyses.',
              'Alinhamentos trimados foram gerados, mas os alinhamentos sem trimagem foram usados nas análises.',
              'Se generaron alineamientos recortados, pero en los análisis se usaron los alineamientos sin recorte.');
    }
    final clip = trimmed.where((l) => (l.trimMethod ?? '').startsWith('ClipKIT')).length;
    final trimal = trimmed.length - clip;
    final s = StringBuffer();
    if (trimal > 0) {
      final bt = trimmed.where((l) => (l.trimMethod ?? '').contains('-backtrans')).length;
      final ta = named('trimAl', ExternalTool.trimal);
      role('trimAl', tr(
        'Removed poorly aligned, gap-rich columns from $trimal ${plural(trimal, 'alignment', 'alignments')} (${c.trimal.label})${bt > 0 ? ', keeping whole codons in $bt coding ${plural(bt, 'gene', 'genes')} (-backtrans)' : ''}.',
        'Removeu colunas mal alinhadas e com muitos gaps de $trimal ${plural(trimal, 'alinhamento', 'alinhamentos')} (${c.trimal.label})${bt > 0 ? ', mantendo códons inteiros em $bt ${plural(bt, 'gene codificador', 'genes codificadores')} (-backtrans)' : ''}.',
        'Eliminó columnas mal alineadas y con muchos gaps de $trimal ${plural(trimal, 'alineamiento', 'alineamientos')} (${c.trimal.label})${bt > 0 ? ', conservando codones enteros en $bt ${plural(bt, 'gen codificante', 'genes codificantes')} (-backtrans)' : ''}.',
      ));
      s.write(tr(
        'Poorly aligned and gap-rich regions were removed with $ta ${cite(['trimAl'])} (${c.trimal.label})'
            '${bt > 0 ? '; the $bt coding ${plural(bt, 'gene was', 'genes were')} translated codon by codon with the ${c.geneticCode == 1 ? 'standard' : 'selected'} genetic code, trimmed as amino acids and back-translated by trimAl (-backtrans), so that whole codons were kept or removed and the reading frame was preserved' : ''}. ',
        'Regiões mal alinhadas e com muitos gaps foram removidas com o $ta ${cite(['trimAl'])} (${c.trimal.label})'
            '${bt > 0 ? '; ${plural(bt, 'o gene codificador foi traduzido', 'os $bt genes codificadores foram traduzidos')} códon a códon, ${plural(bt, 'trimado', 'trimados')} como aminoácidos e ${plural(bt, 'retrotraduzido', 'retrotraduzidos')} pelo trimAl (-backtrans), de modo que códons inteiros foram mantidos ou removidos e a fase de leitura foi preservada' : ''}. ',
        'Las regiones mal alineadas y con muchos gaps se eliminaron con $ta ${cite(['trimAl'])} (${c.trimal.label})'
            '${bt > 0 ? '; ${plural(bt, 'el gen codificante se tradujo', 'los $bt genes codificantes se tradujeron')} codón a codón, se ${plural(bt, 'recortó', 'recortaron')} como aminoácidos y trimAl ${plural(bt, 'lo retrotradujo', 'los retrotradujo')} (-backtrans), de modo que se conservaron o eliminaron codones enteros y se preservó el marco de lectura' : ''}. ',
      ));
    }
    if (clip > 0) {
      final codons = trimmed.where((l) => (l.trimMethod ?? '').contains('--codon')).length;
      final plain = clip - codons;
      final ck = named('ClipKIT', ExternalTool.clipkit);
      role('ClipKIT', tr(
        'Trimmed $clip ${plural(clip, 'alignment', 'alignments')} (${c.clipkitMode} mode)${codons > 0 ? ', removing whole codons in $codons coding ${plural(codons, 'gene', 'genes')} (--codon)' : ''}.',
        'Trimou $clip ${plural(clip, 'alinhamento', 'alinhamentos')} (modo ${c.clipkitMode})${codons > 0 ? ', removendo códons inteiros em $codons ${plural(codons, 'gene codificador', 'genes codificadores')} (--codon)' : ''}.',
        'Recortó $clip ${plural(clip, 'alineamiento', 'alineamientos')} (modo ${c.clipkitMode})${codons > 0 ? ', eliminando codones enteros en $codons ${plural(codons, 'gen codificante', 'genes codificantes')} (--codon)' : ''}.',
      ));
      s.write(tr(
        'Alignments were trimmed with $ck ${cite(['ClipKIT'])} in the ${c.clipkitMode} mode'
            '${codons > 0 ? '; the $codons coding ${plural(codons, 'gene was', 'genes were')} trimmed with its codon option (--codon), which removes whole codons so that the reading frame is preserved' : ''}'
            '${plain > 0 && codons > 0 ? ', and the $plain non-coding ${plural(plain, 'marker', 'markers')} without it' : ''}. ',
        'Os alinhamentos foram trimados com o $ck ${cite(['ClipKIT'])} no modo ${c.clipkitMode}'
            '${codons > 0 ? '; ${plural(codons, 'o gene codificador foi trimado', 'os $codons genes codificadores foram trimados')} com a opção de códons (--codon), que remove códons inteiros e preserva a fase de leitura' : ''}'
            '${plain > 0 && codons > 0 ? ', e ${plural(plain, 'o marcador não codificante', 'os $plain marcadores não codificantes')} sem ela' : ''}. ',
        'Los alineamientos se recortaron con $ck ${cite(['ClipKIT'])} en el modo ${c.clipkitMode}'
            '${codons > 0 ? '; ${plural(codons, 'el gen codificante se recortó', 'los $codons genes codificantes se recortaron')} con su opción de codones (--codon), que elimina codones completos y preserva el marco de lectura' : ''}'
            '${plain > 0 && codons > 0 ? ', y ${plural(plain, 'el marcador no codificante', 'los $plain marcadores no codificantes')} sin ella' : ''}. ',
      ));
    }
    return s.toString().trim();
  }

  String _verdict(SaturationVerdict v) => switch (v) {
        SaturationVerdict.ok => tr('little or no saturation', 'pouca ou nenhuma saturação', 'poca o ninguna saturación'),
        SaturationVerdict.moderate => tr('moderate saturation', 'saturação moderada', 'saturación moderada'),
        SaturationVerdict.saturated => tr('substantial saturation', 'saturação acentuada', 'saturación marcada'),
        SaturationVerdict.insufficient => tr('too little data to be assessed', 'dados insuficientes para avaliação', 'datos insuficientes para evaluarla'),
      };

  String? _saturation() {
    final pooled = c.pooledSaturation;
    if (pooled == null) return null;
    final model = c.satModel;
    final n = pooled.reports.length;
    final th = c.saturationThresholds;
    final iq = named(iqtreeRef, ExternalTool.iqtree);
    final pk = named('PhyKIT', ExternalTool.phykit);
    role('PhyKIT', tr(
      'Measured substitution saturation of ${plural(n, 'the coding gene', 'the $n coding genes')} at each codon position (regression of p-distances on tree distances).',
      'Mediu a saturação de substituições ${plural(n, 'do gene codificador', 'dos $n genes codificadores')} em cada posição do códon (regressão das distâncias p sobre as distâncias na árvore).',
      'Midió la saturación de sustituciones ${plural(n, 'del gen codificante', 'de los $n genes codificantes')} en cada posición del codón (regresión de las distancias p sobre las distancias en el árbol).',
    ));
    role(iqtreeRef, tr('Inferred the tree of each gene for the saturation test.', 'Inferiu a árvore de cada gene para o teste de saturação.',
        'Infirió el árbol de cada gen para la prueba de saturación.'));
    role('Philippe 2011', tr('Method of the saturation test (p-distance against tree distance).', 'Método do teste de saturação (distância p contra distância na árvore).',
        'Método de la prueba de saturación (distancia p contra distancia en el árbol).'));
    final s = StringBuffer(tr(
      'Substitution saturation was assessed for the 1st, 2nd and 3rd codon positions, the 1st + 2nd positions and all positions of '
          '${n == 1 ? 'the coding gene' : 'each of the $n coding genes separately'} with the saturation function of $pk ${cite(['PhyKIT'])}, '
          'which regresses the uncorrected p-distances between taxa on their patristic distances ${cite(['Philippe 2011'])}. '
          'The tree of each gene was inferred in $iq ${cite([iqtreeRef])} under $model (fast search), and the branch lengths of each codon subset '
          'were estimated on that topology. Slopes below ${th.slopeModerate} were taken as moderate and below ${th.slopeSaturated} as substantial saturation. ',
      'A saturação de substituições foi avaliada para a 1ª, a 2ª e a 3ª posições do códon, para as posições 1ª + 2ª e para todas as posições '
          '${n == 1 ? 'do gene codificador' : 'de cada um dos $n genes codificadores separadamente'} com a função de saturação do $pk ${cite(['PhyKIT'])}, '
          'que faz a regressão das distâncias p não corrigidas entre táxons sobre as suas distâncias patrísticas ${cite(['Philippe 2011'])}. '
          'A árvore de cada gene foi inferida no $iq ${cite([iqtreeRef])} sob $model (busca rápida), e os comprimentos de ramo de cada subconjunto de posições '
          'foram estimados nessa topologia. Inclinações abaixo de ${th.slopeModerate} foram consideradas saturação moderada e abaixo de ${th.slopeSaturated}, acentuada. ',
      'La saturación de sustituciones se evaluó para la 1ª, la 2ª y la 3ª posiciones del codón, para las posiciones 1ª + 2ª y para todas las posiciones '
          '${n == 1 ? 'del gen codificante' : 'de cada uno de los $n genes codificantes por separado'} con la función de saturación de $pk ${cite(['PhyKIT'])}, '
          'que hace la regresión de las distancias p no corregidas entre taxones sobre sus distancias patrísticas ${cite(['Philippe 2011'])}. '
          'El árbol de cada gen se infirió en $iq ${cite([iqtreeRef])} bajo $model (búsqueda rápida), y las longitudes de rama de cada subconjunto de posiciones '
          'se estimaron sobre esa topología. Pendientes por debajo de ${th.slopeModerate} se consideraron saturación moderada y por debajo de ${th.slopeSaturated}, marcada. ',
    ));
    // Results by codon position: one gene by verdict, several by count.
    if (n == 1) {
      final byVerdict = <SaturationVerdict, List<int>>{};
      for (final p in [1, 2, 3]) {
        final v = pooled.verdict(p);
        if (v != SaturationVerdict.insufficient) byVerdict.putIfAbsent(v, () => []).add(p);
      }
      if (byVerdict.isNotEmpty) {
        final parts = [
          for (final e in byVerdict.entries)
            tr(
              '${_verdict(e.key)} at the ${list([for (final p in e.value) ordinal(p)])} ${plural(e.value.length, 'position', 'positions')}',
              '${_verdict(e.key)} na${e.value.length == 1 ? '' : 's'} ${list([for (final p in e.value) ordinal(p)])} ${plural(e.value.length, 'posição', 'posições')}',
              '${_verdict(e.key)} en la${e.value.length == 1 ? '' : 's'} ${list([for (final p in e.value) ordinal(p)])} ${plural(e.value.length, 'posición', 'posiciones')}',
            ),
        ];
        final s3 = pooled.reports.single.subsets[CodonSubset.pos3]?.slope;
        s.write(tr(
          '${pooled.reports.single.name} showed ${list(parts)}${s3 == null ? '' : ' (saturation of the 3rd position = ${s3.toStringAsFixed(3)})'}. ',
          '${pooled.reports.single.name} apresentou ${list(parts)}${s3 == null ? '' : ' (saturação da 3ª posição = ${s3.toStringAsFixed(3)})'}. ',
          '${pooled.reports.single.name} mostró ${list(parts)}${s3 == null ? '' : ' (saturación de la 3ª posición = ${s3.toStringAsFixed(3)})'}. ',
        ));
      }
    } else {
      final parts = [
        for (final p in [1, 2, 3])
          if (pooled.assessed(p) > 0)
            tr(
              'the ${ordinal(p)} position was saturated in ${pooled.count(p, SaturationVerdict.saturated)} of ${pooled.assessed(p)} genes',
              'a ${ordinal(p)} posição apresentou saturação acentuada em ${pooled.count(p, SaturationVerdict.saturated)} de ${pooled.assessed(p)} genes',
              'la ${ordinal(p)} posición mostró saturación marcada en ${pooled.count(p, SaturationVerdict.saturated)} de ${pooled.assessed(p)} genes',
            ),
      ];
      if (parts.isNotEmpty) {
        final joined = list(parts);
        final cap = '${joined[0].toUpperCase()}${joined.substring(1)}';
        s.write(tr(
          '$cap; a position was considered saturated when it was saturated in more than half of the genes. ',
          '$cap; uma posição foi considerada saturada quando estava saturada em mais da metade dos genes. ',
          '$cap; una posición se consideró saturada cuando lo estaba en más de la mitad de los genes. ',
        ));
      }
    }
    final excluded = c.concat.excludePositions.toList()..sort();
    final perLocusExcl = c.concat.perLocusExclude.entries.where((e) => e.value.isNotEmpty).toList();
    if (excluded.isNotEmpty) {
      final pos = list([for (final p in excluded) ordinal(p)]);
      s.write(tr(
        'Based on these results, the $pos codon ${plural(excluded.length, 'position was', 'positions were')} excluded from the coding genes of the supermatrix.',
        'Com base nesses resultados, ${plural(excluded.length, 'a $pos posição do códon foi excluída', 'as $pos posições do códon foram excluídas')} dos genes codificadores da supermatriz.',
        'Con base en estos resultados, ${plural(excluded.length, 'la $pos posición del codón se excluyó', 'las $pos posiciones del codón se excluyeron')} de los genes codificantes de la supermatriz.',
      ));
    } else if (perLocusExcl.isNotEmpty) {
      final items = [for (final e in perLocusExcl) '${e.key} (${list([for (final p in (e.value.toList()..sort())) ordinal(p)])})'];
      s.write(tr(
        'Based on the results of each gene, saturated codon positions were excluded gene by gene: ${list(items)}.',
        'Com base nos resultados de cada gene, as posições saturadas do códon foram excluídas gene a gene: ${list(items)}.',
        'Con base en los resultados de cada gen, las posiciones saturadas del codón se excluyeron gen por gen: ${list(items)}.',
      ));
    } else {
      s.write(tr('All codon positions were retained.', 'Todas as posições do códon foram mantidas.', 'Se conservaron todas las posiciones del codón.'));
    }
    return s.toString().trim();
  }

  String? _matrix() {
    final m = c.matrix;
    if (m == null) return null;
    final loci = m.locusRanges.length;
    final sites = m.isProtein ? tr('amino-acid sites', 'sítios de aminoácidos', 'sitios de aminoácidos') : tr('nucleotide sites', 'sítios de nucleotídeos', 'sitios de nucleótidos');
    final s = StringBuffer(tr(
      'The $loci alignments were concatenated into a supermatrix of ${m.taxa.length} taxa and ${m.length} $sites; taxa lacking a gene were coded as missing data (?)'
          '${c.concat.replaceGapsWithMissing ? ', as were alignment gaps' : ''}. ',
      'Os $loci alinhamentos foram concatenados em uma supermatriz de ${m.taxa.length} táxons e ${m.length} $sites; os táxons sem um gene foram codificados como dados ausentes (?)'
          '${c.concat.replaceGapsWithMissing ? ', assim como os gaps do alinhamento' : ''}. ',
      'Los $loci alineamientos se concatenaron en una supermatriz de ${m.taxa.length} taxones y ${m.length} $sites; los taxones sin un gen se codificaron como datos faltantes (?)'
          '${c.concat.replaceGapsWithMissing ? ', al igual que los gaps del alineamiento' : ''}. ',
    ));
    final n = m.partitions.length;
    if (c.effectiveScheme == PartitionScheme.none || n == 0) {
      s.write(tr('The supermatrix was analysed without partitions.', 'A supermatriz foi analisada sem partições.', 'La supermatriz se analizó sin particiones.'));
    } else {
      final how = switch (c.effectiveScheme) {
        PartitionScheme.byLocus => tr('by gene', 'por gene', 'por gen'),
        PartitionScheme.byLocusAndCodon => tr('by gene and, in the coding genes, by codon position', 'por gene e, nos genes codificadores, por posição do códon',
            'por gen y, en los genes codificantes, por posición del codón'),
        PartitionScheme.byCodonPosition => tr('by codon position (across all coding genes)', 'por posição do códon (somando todos os genes codificadores)',
            'por posición del codón (sumando todos los genes codificantes)'),
        PartitionScheme.none => '',
      };
      s.write(tr(
        'The supermatrix was partitioned $how, giving $n ${plural(n, 'partition', 'partitions')}.',
        'A supermatriz foi particionada $how, resultando em $n ${plural(n, 'partição', 'partições')}.',
        'La supermatriz se particionó $how, lo que dio $n ${plural(n, 'partición', 'particiones')}.',
      ));
    }
    return s.toString().trim();
  }

  /// One paragraph per inference program (the latest finished run of each).
  List<String> _inference() {
    final out = <String>[];
    final done = c.treeRuns.where((r) => r.status == StepStatus.done);
    final iq = done.where((r) => r.engine == 'IQ-TREE' && r.iqtree != null).firstOrNull;
    final rx = done.where((r) => r.engine == 'RAxML-NG' && r.raxml != null).firstOrNull;
    final partitioned = c.matrix != null && c.matrix!.partitions.isNotEmpty && c.effectiveScheme != PartitionScheme.none;
    if (iq != null) out.add(_iqtree(iq, partitioned));
    if (rx != null) out.add(_raxml(rx, partitioned));
    if (out.isNotEmpty && c.concat.outgroups.isNotEmpty) {
      final info = c.tipInfo();
      final og = list([
        for (final o in c.concat.outgroups)
          if (info[o] case final i? when i.species.isNotEmpty) '${i.species}${i.accession.isNotEmpty ? ' (${i.accession})' : ''}' else o.replaceAll('_', ' '),
      ]);
      out[out.length - 1] = '${out.last} ${tr(
        'Trees were rooted with $og as the outgroup.',
        'As árvores foram enraizadas com $og como grupo externo.',
        'Los árboles se enraizaron con $og como grupo externo.',
      )}';
    }
    return out;
  }

  String _iqtree(TreeRun run, bool partitioned) {
    final o = run.iqtree!;
    final s = StringBuffer(tr(
      'Maximum-likelihood (ML) phylogenetic inference was performed in ${named(iqtreeRef, ExternalTool.iqtree)} ${cite([iqtreeRef])}. ',
      'A inferência filogenética por máxima verossimilhança (ML) foi feita no ${named(iqtreeRef, ExternalTool.iqtree)} ${cite([iqtreeRef])}. ',
      'La inferencia filogenética por máxima verosimilitud (ML) se realizó en ${named(iqtreeRef, ExternalTool.iqtree)} ${cite([iqtreeRef])}. ',
    ));
    final merge = o.model.contains('MERGE');
    final rcluster = run.command?.contains('--rcluster') ?? false;
    role(iqtreeRef, tr(
      'Inferred the maximum-likelihood tree of the supermatrix${run.bestModel == null ? '' : ' (${run.bestModel})'}.',
      'Inferiu a árvore de máxima verossimilhança da supermatriz${run.bestModel == null ? '' : ' (${run.bestModel})'}.',
      'Infirió el árbol de máxima verosimilitud de la supermatriz${run.bestModel == null ? '' : ' (${run.bestModel})'}.',
    ));
    if (o.model.startsWith('MF')) {
      role('ModelFinder', tr('Selected the substitution model${partitioned ? ' of each partition' : ''} (BIC).',
          'Escolheu o modelo de substituição${partitioned ? ' de cada partição' : ''} (BIC).', 'Eligió el modelo de sustitución${partitioned ? ' de cada partición' : ''} (BIC).'));
      if (merge && partitioned) {
        role('Partition merging', tr('Merged partitions with similar models.', 'Uniu as partições com modelos semelhantes.', 'Unió las particiones con modelos semejantes.'));
      }
      if (merge && partitioned && rcluster) {
        role('rcluster', tr('Relaxed clustering of the partition merging (top 10% of pairs).', 'Agrupamento relaxado da união das partições (10% dos pares).',
            'Agrupamiento relajado de la unión de particiones (10% de los pares).'));
      }
    }
    if (o.ufboot > 0) {
      role('UFBoot2', tr('Branch support: ${o.ufboot} ultrafast bootstrap replicates.', 'Suporte dos ramos: ${o.ufboot} réplicas de ultrafast bootstrap.',
          'Soporte de las ramas: ${o.ufboot} réplicas de ultrafast bootstrap.'));
    }
    if (o.alrt > 0) {
      role('SH-aLRT', tr('Branch support: ${o.alrt} SH-aLRT replicates.', 'Suporte dos ramos: ${o.alrt} réplicas do SH-aLRT.', 'Soporte de las ramas: ${o.alrt} réplicas de SH-aLRT.'));
    }
    if (o.model.startsWith('MF')) {
      final best = !partitioned && run.bestModel != null
          ? tr(' (best-fit model: ${run.bestModel})', ' (modelo escolhido: ${run.bestModel})', ' (modelo elegido: ${run.bestModel})')
          : '';
      if (partitioned) {
        s.write(tr(
          'The substitution model of each partition was selected with ModelFinder ${cite(['ModelFinder'])} by the Bayesian information criterion'
              '${merge ? ', merging partitions with similar models ${cite(['Partition merging'])}${rcluster ? ' with the relaxed clustering algorithm (top 10% of pairs) ${cite(['rcluster'])}' : ''}' : ''}. ',
          'O modelo de substituição de cada partição foi escolhido com o ModelFinder ${cite(['ModelFinder'])} pelo critério de informação bayesiano (BIC)'
              '${merge ? ', unindo as partições com modelos semelhantes ${cite(['Partition merging'])}${rcluster ? ' com o algoritmo de agrupamento relaxado (10% dos pares) ${cite(['rcluster'])}' : ''}' : ''}. ',
          'El modelo de sustitución de cada partición se eligió con ModelFinder ${cite(['ModelFinder'])} según el criterio de información bayesiano (BIC)'
              '${merge ? ', uniendo las particiones con modelos semejantes ${cite(['Partition merging'])}${rcluster ? ' con el algoritmo de agrupamiento relajado (10% de los pares) ${cite(['rcluster'])}' : ''}' : ''}. ',
        ));
      } else {
        s.write(tr(
          'The best-fit substitution model was selected with ModelFinder ${cite(['ModelFinder'])} by the Bayesian information criterion$best. ',
          'O modelo de substituição mais adequado foi escolhido com o ModelFinder ${cite(['ModelFinder'])} pelo critério de informação bayesiano (BIC)$best. ',
          'El modelo de sustitución más adecuado se eligió con ModelFinder ${cite(['ModelFinder'])} según el criterio de información bayesiano (BIC)$best. ',
        ));
      }
    } else {
      s.write(tr('The ${o.model} model was used. ', 'Foi usado o modelo ${o.model}. ', 'Se usó el modelo ${o.model}. '));
    }
    if (o.runs > 1) {
      s.write(tr(
        'The tree search was repeated in ${o.runs} independent runs and the tree with the highest likelihood was kept. ',
        'A busca da árvore foi repetida em ${o.runs} execuções independentes, e a árvore de maior verossimilhança foi mantida. ',
        'La búsqueda del árbol se repitió en ${o.runs} ejecuciones independientes y se conservó el árbol de mayor verosimilitud. ',
      ));
    }
    if (o.nstop != null || o.perturbation != null) {
      final parts = [
        if (o.nstop != null) tr('stopping after ${o.nstop} unsuccessful iterations', 'parando após ${o.nstop} iterações sem melhora', 'deteniéndose tras ${o.nstop} iteraciones sin mejora'),
        if (o.perturbation != null) tr('perturbation strength ${o.perturbation}', 'força de perturbação ${o.perturbation}', 'fuerza de perturbación ${o.perturbation}'),
      ];
      s.write(tr('Search settings: ${list(parts)}. ', 'Ajustes da busca: ${list(parts)}. ', 'Ajustes de la búsqueda: ${list(parts)}. '));
    }
    final support = [
      if (o.ufboot > 0)
        tr('${o.ufboot} ultrafast bootstrap replicates (UFBoot2) ${cite(['UFBoot2'])}${o.bnni ? ', each bootstrap tree optimised by NNI (--bnni)' : ''}',
            '${o.ufboot} réplicas de ultrafast bootstrap (UFBoot2) ${cite(['UFBoot2'])}${o.bnni ? ', cada árvore de bootstrap otimizada por NNI (--bnni)' : ''}',
            '${o.ufboot} réplicas de ultrafast bootstrap (UFBoot2) ${cite(['UFBoot2'])}${o.bnni ? ', cada árbol de bootstrap optimizado por NNI (--bnni)' : ''}'),
      if (o.alrt > 0)
        tr('${o.alrt} SH-aLRT replicates ${cite(['SH-aLRT'])}', '${o.alrt} réplicas do teste SH-aLRT ${cite(['SH-aLRT'])}', '${o.alrt} réplicas de la prueba SH-aLRT ${cite(['SH-aLRT'])}'),
      if (o.classicBootstrap > 0)
        tr('${o.classicBootstrap} standard non-parametric bootstrap replicates', '${o.classicBootstrap} réplicas de bootstrap não paramétrico padrão',
            '${o.classicBootstrap} réplicas de bootstrap no paramétrico estándar'),
    ];
    if (support.isNotEmpty) {
      s.write(tr(
        'Branch support was estimated with ${list(support)}.',
        'O suporte dos ramos foi estimado com ${list(support)}.',
        'El soporte de las ramas se estimó con ${list(support)}.',
      ));
      if (o.ufboot > 0 && o.alrt > 0) {
        s.write(tr(
          ' Clades with SH-aLRT ≥ 80% and UFBoot ≥ 95% were considered strongly supported.',
          ' Clados com SH-aLRT ≥ 80% e UFBoot ≥ 95% foram considerados fortemente sustentados.',
          ' Los clados con SH-aLRT ≥ 80% y UFBoot ≥ 95% se consideraron fuertemente sustentados.',
        ));
      } else if (o.ufboot > 0) {
        s.write(tr(
          ' Clades with UFBoot ≥ 95% were considered strongly supported.',
          ' Clados com UFBoot ≥ 95% foram considerados fortemente sustentados.',
          ' Los clados con UFBoot ≥ 95% se consideraron fuertemente sustentados.',
        ));
      }
    }
    if (run.outputs['consensus'] != null) {
      s.write(tr(
        'The majority-rule consensus tree of the bootstrap replicates${run.engine == 'IQ-TREE' ? '' : ' (RAxML-NG --consense MR)'} is presented. ',
        'É apresentada a árvore de consenso de maioria das réplicas de bootstrap${run.engine == 'IQ-TREE' ? '' : ' (RAxML-NG --consense MR)'}. ',
        'Se presenta el árbol de consenso de mayoría de las réplicas de bootstrap${run.engine == 'IQ-TREE' ? '' : ' (RAxML-NG --consense MR)'}. ',
      ));
    }
    return s.toString().trim();
  }

  String _raxml(TreeRun run, bool partitioned) {
    final o = run.raxml!;
    role('RAxML-NG', tr(
      'Inferred a second maximum-likelihood tree of the supermatrix${run.bestModel == null ? '' : ' (${run.bestModel})'}, with bootstrap support.',
      'Inferiu uma segunda árvore de máxima verossimilhança da supermatriz${run.bestModel == null ? '' : ' (${run.bestModel})'}, com suporte de bootstrap.',
      'Infirió un segundo árbol de máxima verosimilitud de la supermatriz${run.bestModel == null ? '' : ' (${run.bestModel})'}, con soporte de bootstrap.',
    ));
    role(o.fast ? 'RAxML-NG fast' : 'Adaptive RAxML-NG', tr(
      o.fast ? 'Fast tree search with early stopping (--fast).' : 'Adaptive tree search, its effort set by the difficulty of the alignment.',
      o.fast ? 'Busca rápida com parada antecipada (--fast).' : 'Busca adaptativa, com esforço ajustado à dificuldade do alinhamento.',
      o.fast ? 'Búsqueda rápida con parada anticipada (--fast).' : 'Búsqueda adaptativa, con esfuerzo ajustado a la dificultad del alineamiento.',
    ));
    for (final m in o.metrics) {
      switch (m) {
        case RaxmlSupport.tbe:
          role('TBE', tr('Branch support of RAxML-NG: transfer bootstrap expectation.', 'Suporte dos ramos no RAxML-NG: transfer bootstrap expectation.',
              'Soporte de las ramas en RAxML-NG: transfer bootstrap expectation.'));
        case RaxmlSupport.rbs:
          role('RBS', tr('Branch support of RAxML-NG: rapid bootstrap.', 'Suporte dos ramos no RAxML-NG: bootstrap rápido.', 'Soporte de las ramas en RAxML-NG: bootstrap rápido.'));
        case RaxmlSupport.ebg:
          role('EBG', tr('Branch support of RAxML-NG predicted by machine learning (EBG).', 'Suporte dos ramos no RAxML-NG previsto por aprendizado de máquina (EBG).',
              'Soporte de las ramas en RAxML-NG predicho por aprendizaje automático (EBG).'));
        case RaxmlSupport.sh:
          role('SH-aLRT', tr('Branch support of RAxML-NG: SH-like aLRT.', 'Suporte dos ramos no RAxML-NG: SH-aLRT.', 'Soporte de las ramas en RAxML-NG: SH-aLRT.'));
        case RaxmlSupport.fbp:
          break;
      }
    }
    if (run.models == RaxmlModels.modelFinder) {
      role('ModelFinder', tr('Selected the models used by RAxML-NG.', 'Escolheu os modelos usados no RAxML-NG.', 'Eligió los modelos usados en RAxML-NG.'));
    }
    final s = StringBuffer(tr(
      'ML inference was also performed with ${named('RAxML-NG', ExternalTool.raxmlng)} ${cite(['RAxML-NG'])}',
      'A inferência por ML também foi feita no ${named('RAxML-NG', ExternalTool.raxmlng)} ${cite(['RAxML-NG'])}',
      'La inferencia por ML también se realizó en ${named('RAxML-NG', ExternalTool.raxmlng)} ${cite(['RAxML-NG'])}',
    ));
    switch (run.models) {
      case RaxmlModels.moose:
        final chosen = run.mooseModels;
        final crit = o.criterion;
        s.write(chosen.length == 1
            ? tr(
                ' under the ${chosen.single} model, selected by its built-in MOOSE module (RAxML-NG v2) according to $crit. ',
                ' com o modelo ${chosen.single}, escolhido pelo módulo MOOSE do próprio programa (RAxML-NG v2) segundo o $crit. ',
                ' con el modelo ${chosen.single}, elegido por su módulo MOOSE (RAxML-NG v2) según el $crit. ',
              )
            : tr(
                ', with ${partitioned ? 'one substitution model per partition' : 'the substitution model'} selected by its built-in MOOSE module (RAxML-NG v2) according to $crit. ',
                ', com ${partitioned ? 'um modelo de substituição por partição' : 'o modelo de substituição'} escolhido pelo módulo MOOSE do próprio programa (RAxML-NG v2) segundo o $crit. ',
                ', con ${partitioned ? 'un modelo de sustitución por partición' : 'el modelo de sustitución'} elegido por su módulo MOOSE (RAxML-NG v2) según el $crit. ',
              ));
      case RaxmlModels.modelFinder:
        s.write(tr(
          ', with the ${partitioned ? 'partitioning scheme and models' : 'substitution model'} selected by ModelFinder ${cite(['ModelFinder'])} among the models available in RAxML-NG. ',
          ', com ${partitioned ? 'o esquema de partições e os modelos escolhidos' : 'o modelo de substituição escolhido'} pelo ModelFinder ${cite(['ModelFinder'])} entre os modelos disponíveis no RAxML-NG. ',
          ', con ${partitioned ? 'el esquema de particiones y los modelos elegidos' : 'el modelo de sustitución elegido'} por ModelFinder ${cite(['ModelFinder'])} entre los modelos disponibles en RAxML-NG. ',
        ));
      case RaxmlModels.fixed:
        final model = run.fixedModel ?? (c.matrix?.isProtein ?? false ? 'LG+G' : 'GTR+G');
        s.write(tr(' under the $model model. ', ' com o modelo $model. ', ' con el modelo $model. '));
    }
    final starts = RegExp(r'(pars|rand)\{(\d+)\}').allMatches(o.startTrees).map((m) {
      final n = m.group(2)!;
      return m.group(1) == 'pars' ? tr('$n parsimony', '$n de parcimônia', '$n de parsimonia') : tr('$n random', '$n aleatórias', '$n aleatorios');
    }).toList();
    if (o.fast) {
      s.write(tr(
        'The tree search used the fast mode (--fast), a simplified heuristic with early stopping ${cite(['RAxML-NG fast'])}. ',
        'A busca usou o modo rápido (--fast), uma heurística simplificada com parada antecipada ${cite(['RAxML-NG fast'])}. ',
        'La búsqueda usó el modo rápido (--fast), una heurística simplificada con parada anticipada ${cite(['RAxML-NG fast'])}. ',
      ));
    } else {
      s.write(tr(
        'The tree search used the adaptive heuristic, which sets its effort from the predicted difficulty of the alignment ${cite(['Adaptive RAxML-NG'])}. ',
        'A busca usou a heurística adaptativa, que ajusta o esforço à dificuldade prevista do alinhamento ${cite(['Adaptive RAxML-NG'])}. ',
        'La búsqueda usó la heurística adaptativa, que ajusta el esfuerzo a la dificultad prevista del alineamiento ${cite(['Adaptive RAxML-NG'])}. ',
      ));
    }
    if (starts.isNotEmpty) {
      s.write(tr(
        'It started from ${list(starts)} trees. ',
        'Ela partiu de ${list(starts)} árvores iniciais. ',
        'Partió de ${list(starts)} árboles iniciales. ',
      ));
    }
    final boots = o.bootstraps.trim();
    final auto = RegExp(r'autoMRE\{(\d+)\}').firstMatch(boots);
    final metrics = o.metrics.toSet();
    final summary = [
      if (metrics.contains(RaxmlSupport.fbp)) tr('Felsenstein bootstrap (FBP)', 'bootstrap de Felsenstein (FBP)', 'bootstrap de Felsenstein (FBP)'),
      if (metrics.contains(RaxmlSupport.tbe))
        tr('transfer bootstrap expectation (TBE) ${cite(['TBE'])}', 'transfer bootstrap expectation (TBE) ${cite(['TBE'])}', 'transfer bootstrap expectation (TBE) ${cite(['TBE'])}'),
      if (metrics.contains(RaxmlSupport.rbs)) tr('rapid bootstrap (RBS) ${cite(['RBS'])}', 'bootstrap rápido (RBS) ${cite(['RBS'])}', 'bootstrap rápido (RBS) ${cite(['RBS'])}'),
    ];
    final direct = [
      if (metrics.contains(RaxmlSupport.sh)) tr('the SH-like aLRT test ${cite(['SH-aLRT'])}', 'o teste SH-aLRT ${cite(['SH-aLRT'])}', 'la prueba SH-aLRT ${cite(['SH-aLRT'])}'),
      if (metrics.contains(RaxmlSupport.ebg))
        tr('the Educated Bootstrap Guesser (EBG), which predicts bootstrap values by machine learning ${cite(['EBG'])}',
            'o Educated Bootstrap Guesser (EBG), que prevê os valores de bootstrap por aprendizado de máquina ${cite(['EBG'])}',
            'el Educated Bootstrap Guesser (EBG), que predice los valores de bootstrap por aprendizaje automático ${cite(['EBG'])}'),
    ];
    if (o.needsReplicates && boots.isNotEmpty && boots != '0') {
      final reps = auto != null
          ? tr('up to ${auto.group(1)} bootstrap replicates, stopped by the autoMRE convergence criterion',
              'até ${auto.group(1)} réplicas de bootstrap, interrompidas pelo critério de convergência autoMRE',
              'hasta ${auto.group(1)} réplicas de bootstrap, detenidas por el criterio de convergencia autoMRE')
          : tr('$boots bootstrap replicates', '$boots réplicas de bootstrap', '$boots réplicas de bootstrap');
      s.write(tr(
        'Branch support was estimated with $reps${summary.isEmpty ? '' : ' and summarized as ${list(summary)}'}. ',
        'O suporte dos ramos foi estimado com $reps${summary.isEmpty ? '' : ' e resumido como ${list(summary)}'}. ',
        'El soporte de las ramas se estimó con $reps${summary.isEmpty ? '' : ' y se resumió como ${list(summary)}'}. ',
      ));
    }
    if (direct.isNotEmpty) {
      s.write(tr(
        'Branch support was ${o.needsReplicates ? 'also ' : ''}computed on the ML tree with ${list(direct)}.',
        'O suporte dos ramos ${o.needsReplicates ? 'também ' : ''}foi calculado na árvore de ML com ${list(direct)}.',
        'El soporte de las ramas ${o.needsReplicates ? 'también ' : ''}se calculó en el árbol de ML con ${list(direct)}.',
      ));
    }
    if (run.outputs['consensus'] != null) {
      s.write(tr(
        'The majority-rule consensus tree of the bootstrap replicates${run.engine == 'IQ-TREE' ? '' : ' (RAxML-NG --consense MR)'} is presented. ',
        'É apresentada a árvore de consenso de maioria das réplicas de bootstrap${run.engine == 'IQ-TREE' ? '' : ' (RAxML-NG --consense MR)'}. ',
        'Se presenta el árbol de consenso de mayoría de las réplicas de bootstrap${run.engine == 'IQ-TREE' ? '' : ' (RAxML-NG --consense MR)'}. ',
      ));
    }
    return s.toString().trim();
  }

  /// Species delimitation: mPTP on the ML trees, ASAP and ABGD on the
  /// alignments (the finished runs).
  String? _delimitation() {
    final all = c.delimitationsDone;
    if (all.isEmpty) return null;
    final s = StringBuffer();
    String n(int k) => '$k';
    String src(DelimitationRun d) => d.source == ProjectController.supermatrixAlignment ? tr('the supermatrix', 'a supermatriz', 'la supermatriz') : d.source;

    // On the trees.
    final done = [for (final d in all) if ((d.method == DelimMethod.ptp || d.method == DelimMethod.mptp) && d.run != null) d];
    final results = [
      for (final d in done)
        tr('${n(d.result!.species.length)} species with ${d.method.label} on the ${d.source} tree', '${n(d.result!.species.length)} espécies com o ${d.method.label} na árvore do ${d.source}',
            '${n(d.result!.species.length)} especies con ${d.method.label} en el árbol de ${d.source}'),
    ];
    if (done.isNotEmpty) {
      final methods = {for (final d in done) d.method};
      final engines = [for (final e in ['IQ-TREE', 'RAxML-NG']) if (done.any((d) => d.source == e)) e];
      final named0 = named('mPTP', ExternalTool.mptp);
      final modelText = [
        if (methods.contains(DelimMethod.mptp))
          tr('the multi-rate PTP model (mPTP), which gives each species its own coalescent rate', 'o modelo PTP de múltiplas taxas (mPTP), que dá a cada espécie a sua própria taxa de coalescência',
              'el modelo PTP de múltiples tasas (mPTP), que da a cada especie su propia tasa de coalescencia'),
        if (methods.contains(DelimMethod.ptp))
          tr('the single-rate PTP model ${cite(['PTP'])}', 'o modelo PTP de taxa única ${cite(['PTP'])}', 'el modelo PTP de tasa única ${cite(['PTP'])}'),
      ];
      final mptpCite = cite(['mPTP']);
      s.write(tr(
        'Species were delimited on the ML ${plural(engines.length, 'tree', 'trees')} of ${list(engines)} with $named0 $mptpCite under ${list(modelText)}, by maximum likelihood; '
            'these models find where the branching pattern of the tree changes from speciation (between species) to coalescence (within species). ',
        'As espécies foram delimitadas ${plural(engines.length, 'na árvore', 'nas árvores')} de ML do ${list(engines)} com o $named0 $mptpCite sob ${list(modelText)}, por máxima verossimilhança; '
            'esses modelos identificam onde o padrão de ramificação da árvore passa de especiação (entre espécies) para coalescência (dentro das espécies). ',
        'Las especies se delimitaron en ${plural(engines.length, 'el árbol', 'los árboles')} de ML de ${list(engines)} con $named0 $mptpCite bajo ${list(modelText)}, por máxima verosimilitud; '
            'estos modelos identifican dónde el patrón de ramificación del árbol pasa de especiación (entre especies) a coalescencia (dentro de las especies). ',
      ));
      final crop = done.any((d) => d.cropped);
      final minbr = done.map((d) => d.minbrAuto).nonNulls.firstOrNull;
      if (crop) {
        s.write(tr(
          'The trees were rooted on the outgroup, which was then removed from the delimitation (--outgroup_crop). ',
          'As árvores foram enraizadas no grupo externo, que depois foi retirado da delimitação (--outgroup_crop). ',
          'Los árboles se enraizaron en el grupo externo, que luego se retiró de la delimitación (--outgroup_crop). ',
        ));
      }
      if (minbr != null) {
        s.write(tr(
          'Branches shorter than ${minbr.toStringAsPrecision(3)} substitutions per site, the minimum length detected from the p-distances of the supermatrix (--minbr_auto), were ignored. ',
          'Ramos mais curtos que ${minbr.toStringAsPrecision(3)} substituições por sítio, o comprimento mínimo detectado a partir das distâncias p da supermatriz (--minbr_auto), foram ignorados. ',
          'Se ignoraron las ramas más cortas que ${minbr.toStringAsPrecision(3)} sustituciones por sitio, la longitud mínima detectada a partir de las distancias p de la supermatriz (--minbr_auto). ',
        ));
      }
      if (done.any((d) => d.result!.lrtPassed == false)) {
        final k = done.where((d) => d.result!.lrtPassed == false).length;
        s.write(tr(
          'In ${plural(k, 'one analysis', 'some analyses')} the likelihood-ratio test did not reject the null model of a single species. ',
          'Em ${plural(k, 'uma análise', 'algumas análises')} o teste da razão de verossimilhança não rejeitou o modelo nulo de uma única espécie. ',
          'En ${plural(k, 'un análisis', 'algunos análisis')} la prueba de razón de verosimilitud no rechazó el modelo nulo de una sola especie. ',
        ));
      }
      role('mPTP', tr('Delimited species on the ML trees: ${list(results)}.', 'Delimitou as espécies nas árvores de ML: ${list(results)}.',
          'Delimitó las especies en los árboles de ML: ${list(results)}.'));
      if (methods.contains(DelimMethod.ptp)) {
        role('PTP', tr('Single-rate model of the delimitation (one coalescent rate for all species).', 'Modelo de taxa única da delimitação (uma taxa de coalescência para todas as espécies).',
            'Modelo de tasa única de la delimitación (una tasa de coalescencia para todas las especies).'));
      }
    }

    // GMYC on the ultrametric trees.
    final gmyc = [for (final d in all) if (d.method == DelimMethod.gmyc && d.run != null) d];
    if (gmyc.isNotEmpty) {
      final engines = [for (final e in ['IQ-TREE', 'RAxML-NG']) if (gmyc.any((d) => d.source == e)) e];
      final parts = [
        for (final d in gmyc)
          tr('${n(d.result!.species.length)} species on the ${d.source} tree', '${n(d.result!.species.length)} espécies na árvore do ${d.source}',
              '${n(d.result!.species.length)} especies en el árbol de ${d.source}'),
      ];
      final crop = gmyc.any((d) => d.cropped);
      s.write(tr(
        '${s.isEmpty ? 'Species were' : 'Species were also'} delimited with the single-threshold GMYC model ${cite(['GMYC'])} as implemented in pGMYC ${cite(['PTP'])}. '
            'The ML ${plural(engines.length, 'tree', 'trees')} of ${list(engines)} ${plural(engines.length, 'was', 'were')} first made ultrametric with LSD2 ${cite(['LSD2'])}, estimating relative dates (root at 0, all tips at 1)'
            '${crop ? ' after rooting on the outgroup, which was then removed' : ' with the root position estimated by LSD2'}; '
            'GMYC delimited ${list(parts)} (likelihood-ratio test against the one-species model, p < 0.01). ',
        'As espécies ${s.isEmpty ? 'foram' : 'também foram'} delimitadas com o modelo GMYC de limiar único ${cite(['GMYC'])}, como implementado no pGMYC ${cite(['PTP'])}. '
            '${plural(engines.length, 'A árvore', 'As árvores')} de ML do ${list(engines)} ${plural(engines.length, 'foi tornada', 'foram tornadas')} ultramétrica${engines.length == 1 ? '' : 's'} com o LSD2 ${cite(['LSD2'])}, estimando datas relativas (raiz em 0, todos os terminais em 1)'
            '${crop ? ', após o enraizamento no grupo externo, que depois foi retirado' : ', com a posição da raiz estimada pelo LSD2'}; '
            'o GMYC delimitou ${list(parts)} (teste da razão de verossimilhança contra o modelo de uma espécie, p < 0,01). ',
        'Las especies ${s.isEmpty ? 'se' : 'también se'} delimitaron con el modelo GMYC de umbral único ${cite(['GMYC'])}, implementado en pGMYC ${cite(['PTP'])}. '
            '${plural(engines.length, 'El árbol', 'Los árboles')} de ML de ${list(engines)} se ${plural(engines.length, 'hizo', 'hicieron')} ultramétrico${engines.length == 1 ? '' : 's'} con LSD2 ${cite(['LSD2'])}, estimando fechas relativas (raíz en 0, todos los terminales en 1)'
            '${crop ? ', tras enraizar en el grupo externo, que luego se retiró' : ', con la posición de la raíz estimada por LSD2'}; '
            'GMYC delimitó ${list(parts)} (prueba de razón de verosimilitud contra el modelo de una especie, p < 0,01). ',
      ));
      role('GMYC', tr('Delimited species on the ultrametric trees: ${list(parts)}.', 'Delimitou as espécies nas árvores ultramétricas: ${list(parts)}.',
          'Delimitó las especies en los árboles ultramétricos: ${list(parts)}.'));
      role('LSD2', tr('Made the ML trees ultrametric (relative dates) for GMYC.', 'Tornou as árvores de ML ultramétricas (datas relativas) para o GMYC.',
          'Hizo ultramétricos los árboles de ML (fechas relativas) para GMYC.'));
      role('PTP', tr('Authors of pGMYC, the GMYC implementation used.', 'Autores do pGMYC, a implementação do GMYC usada.', 'Autores de pGMYC, la implementación de GMYC usada.'));
    }

    // On the alignments: ASAP and ABGD.
    final asap = [for (final d in all) if (d.method == DelimMethod.asap) d];
    final abgd = [for (final d in all) if (d.method == DelimMethod.abgd) d];
    if (asap.isNotEmpty) {
      final dist = asap.first.distance?.label ?? 'Jukes-Cantor (JC69)';
      final parts = [
        for (final d in asap)
          tr('${n(d.result!.species.length)} species in ${src(d)}${d.chosen == 0 ? ' (best ASAP score)' : ' (partition ${d.chosen + 1} by ASAP score)'}',
              '${n(d.result!.species.length)} espécies em ${src(d)}${d.chosen == 0 ? ' (melhor escore do ASAP)' : ' (partição ${d.chosen + 1} pelo escore do ASAP)'}',
              '${n(d.result!.species.length)} especies en ${src(d)}${d.chosen == 0 ? ' (mejor puntuación de ASAP)' : ' (partición ${d.chosen + 1} por la puntuación de ASAP)'}'),
      ];
      s.write(tr(
        '${s.isEmpty ? 'Species' : 'Species were also'} delimited from pairwise $dist distances with ASAP ${cite(['ASAP'])}, which builds partitions by hierarchical clustering and ranks them by the ASAP score; it proposed ${list(parts)}. ',
        'As espécies ${s.isEmpty ? 'foram' : 'também foram'} delimitadas a partir das distâncias par a par $dist com o ASAP ${cite(['ASAP'])}, que constrói partições por agrupamento hierárquico e as ordena pelo escore ASAP; ele propôs ${list(parts)}. ',
        'Las especies ${s.isEmpty ? 'se' : 'también se'} delimitaron a partir de las distancias por pares $dist con ASAP ${cite(['ASAP'])}, que construye particiones por agrupamiento jerárquico y las ordena por la puntuación ASAP; propuso ${list(parts)}. ',
      ));
      role('ASAP', tr('Delimited species from pairwise distances: ${list(parts)}.', 'Delimitou as espécies a partir das distâncias par a par: ${list(parts)}.',
          'Delimitó las especies a partir de las distancias por pares: ${list(parts)}.'));
    }
    if (abgd.isNotEmpty) {
      final dist = abgd.first.distance?.label ?? 'Jukes-Cantor (JC69)';
      String detail(DelimitationRun d) {
        final note = d.chosen < d.partitionNotes.length ? d.partitionNotes[d.chosen].split('|') : const ['init', ''];
        final prior = double.tryParse(note.last);
        final kind = note.first == 'rec' ? tr('recursive', 'recursiva', 'recursiva') : tr('initial', 'inicial', 'inicial');
        return tr('${n(d.result!.species.length)} groups in ${src(d)} (${prior == null ? '' : 'prior P = ${prior.toStringAsPrecision(3)}, '}$kind partition)',
            '${n(d.result!.species.length)} grupos em ${src(d)} (${prior == null ? '' : 'P a priori = ${prior.toStringAsPrecision(3)}, '}partição $kind)',
            '${n(d.result!.species.length)} grupos en ${src(d)} (${prior == null ? '' : 'P a priori = ${prior.toStringAsPrecision(3)}, '}partición $kind)');
      }

      final parts = [for (final d in abgd) detail(d)];
      s.write(tr(
        'The barcode gap was searched with ABGD ${cite(['ABGD'])} on $dist distances, over a range of prior maximum intraspecific divergences; it found ${list(parts)}.',
        'O barcode gap foi buscado com o ABGD ${cite(['ABGD'])} em distâncias $dist, numa faixa de divergências intraespecíficas máximas a priori; ele encontrou ${list(parts)}.',
        'El barcode gap se buscó con ABGD ${cite(['ABGD'])} en distancias $dist, en un intervalo de divergencias intraespecíficas máximas a priori; encontró ${list(parts)}.',
      ));
      role('ABGD', tr('Found the barcode gap and delimited: ${list(parts)}.', 'Encontrou o barcode gap e delimitou: ${list(parts)}.', 'Encontró el barcode gap y delimitó: ${list(parts)}.'));
    }
    return s.toString().trim();
  }
}

/// What each tool of SPLACE does, for the references of tools that did
/// not run in the project yet (in [lang]).
String defaultCitationRole(String name, String lang) {
  String tr(String en, String pt, String es) => switch (lang) { 'pt' => pt, 'es' => es, _ => en };
  return switch (name) {
    'SPLACE' => tr('Runs the whole workflow, from the GenBank or FASTA data to the tree figure.', 'Executa todo o fluxo, dos dados do GenBank ou FASTA até a figura da árvore.',
        'Ejecuta todo el flujo, de los datos de GenBank o FASTA hasta la figura del árbol.'),
    'SynGenes' => tr(
        'Standardizes the names of mitochondrial and chloroplast genes (COX1, COI, cox1… become one name) and builds the NCBI queries of individual markers with every synonym.',
        'Padroniza os nomes dos genes mitocondriais e de cloroplasto (COX1, COI, cox1… viram um só nome) e monta as buscas no NCBI dos marcadores individuais com todos os sinônimos.',
        'Estandariza los nombres de los genes mitocondriales y de cloroplasto (COX1, COI, cox1… pasan a ser un solo nombre) y arma las búsquedas en el NCBI de los marcadores individuales con todos los sinónimos.'),
    'dataFishing' => tr(
        'Retrieves the taxonomy of the records (GBIF, NCBI Taxonomy, ICMBio SALVE) and the biodiversity data used in the headers, the records table and the tree.',
        'Obtém a taxonomia dos registros (GBIF, NCBI Taxonomy, SALVE/ICMBio) e os dados de biodiversidade usados nos cabeçalhos, na tabela de registros e na árvore.',
        'Obtiene la taxonomía de los registros (GBIF, NCBI Taxonomy, SALVE/ICMBio) y los datos de biodiversidad usados en los encabezados, la tabla de registros y el árbol.'),
    'BOLD' => tr('Public barcode records (COI and plant barcodes) searched together with NCBI, and identification of specimens.',
        'Registros públicos de barcodes (COI e barcodes de plantas) buscados junto com o NCBI, e identificação de espécimes.',
        'Registros públicos de barcodes (COI y barcodes de plantas) buscados junto con el NCBI, e identificación de especímenes.'),
    'MAFFT' => tr('Aligns each marker.', 'Alinha cada marcador.', 'Alinea cada marcador.'),
    'trimAl' => tr('Removes poorly aligned regions (whole codons with -backtrans).', 'Remove regiões mal alinhadas (códons inteiros com -backtrans).',
        'Elimina regiones mal alineadas (codones enteros con -backtrans).'),
    'ClipKIT' => tr('Removes poorly aligned regions (whole codons with --codon).', 'Remove regiões mal alinhadas (códons inteiros com --codon).',
        'Elimina regiones mal alineadas (codones enteros con --codon).'),
    'PhyKIT' => tr('Tests substitution saturation per gene and codon position.', 'Testa a saturação de substituições por gene e posição do códon.',
        'Prueba la saturación de sustituciones por gen y posición del codón.'),
    'IQ-TREE 3' || 'IQ-TREE 2' => tr('Infers maximum-likelihood trees, with model selection and branch support.', 'Infere árvores de máxima verossimilhança, com escolha de modelo e suporte dos ramos.',
        'Infiere árboles de máxima verosimilitud, con elección de modelo y soporte de las ramas.'),
    'RAxML-NG' => tr('Infers maximum-likelihood trees with bootstrap support.', 'Infere árvores de máxima verossimilhança com suporte de bootstrap.',
        'Infiere árboles de máxima verosimilitud con soporte de bootstrap.'),
    'mPTP' => tr('Delimits species on the ML trees (PTP and mPTP models).', 'Delimita espécies nas árvores de ML (modelos PTP e mPTP).',
        'Delimita especies en los árboles de ML (modelos PTP y mPTP).'),
    'GMYC' => tr('Delimits species on an ultrametric tree (single-threshold GMYC, pGMYC).', 'Delimita espécies numa árvore ultramétrica (GMYC de limiar único, pGMYC).',
        'Delimita especies en un árbol ultramétrico (GMYC de umbral único, pGMYC).'),
    'LSD2' => tr('Dates the ML tree (relative dates) so GMYC can read it.', 'Data a árvore de ML (datas relativas) para o GMYC.', 'Data el árbol de ML (fechas relativas) para GMYC.'),
    'ASAP' => tr('Delimits species from the pairwise distances of an alignment (ranked partitions).', 'Delimita espécies a partir das distâncias par a par de um alinhamento (partições ordenadas por escore).',
        'Delimita especies a partir de las distancias por pares de un alineamiento (particiones ordenadas por puntuación).'),
    'ABGD' => tr('Finds the barcode gap in the pairwise distances of an alignment and delimits species.', 'Encontra o barcode gap nas distâncias par a par de um alinhamento e delimita espécies.',
        'Encuentra el barcode gap en las distancias por pares de un alineamiento y delimita especies.'),
    _ => '',
  };
}
