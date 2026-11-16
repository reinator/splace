// ignore: unused_import
import 'package:intl/intl.dart' as intl;

import 'app_localizations.dart';

// ignore_for_file: type=lint

/// The translations for Portuguese (`pt`).
class AppLocalizationsPt extends AppLocalizations {
  AppLocalizationsPt([String locale = 'pt']) : super(locale);

  @override
  String get navProgress => 'Progresso';

  @override
  String get navSteps => 'Etapas do pipeline';

  @override
  String get logTitle => 'Registro de atividade';

  @override
  String get logEmpty =>
      'Nenhuma atividade ainda. Os logs de busca, alinhamento e análise aparecerão aqui.';

  @override
  String get logCopy => 'Copiar tudo';

  @override
  String get logDownload => 'Baixar como .txt';

  @override
  String get logClear => 'Limpar';

  @override
  String get logCollapse => 'Ocultar painel de logs';

  @override
  String get logWelcome => 'SPLACE pronto. Comece carregando dados na etapa 1.';

  @override
  String get settingsTheme => 'Tema';

  @override
  String get themeLight => 'Claro';

  @override
  String get themeDark => 'Escuro';

  @override
  String get themeAuto => 'Auto';

  @override
  String get settingsDensity => 'Densidade';

  @override
  String get densityComfy => 'Confortável';

  @override
  String get densityCompact => 'Compacta';

  @override
  String get settingsFontSize => 'Tamanho do texto';

  @override
  String get fontSizeReset => 'Redefinir';

  @override
  String get minimapLocked => 'Conclua a etapa anterior para desbloquear';

  @override
  String get tooltipMenu => 'Mostrar navegador do pipeline';

  @override
  String get tooltipLogs => 'Mostrar registro de atividade';

  @override
  String get tooltipSettings => 'Exibição e idioma';

  @override
  String get tooltipHelp => 'Mostrar tela de boas-vindas';

  @override
  String get tooltipTools => 'Ferramentas externas e configurações de análise';

  @override
  String get cancel => 'Cancelar';

  @override
  String get done => 'Concluído';

  @override
  String andMore(String n) {
    return '…e mais $n';
  }

  @override
  String get splashSummary =>
      'Um kit completo para filogenômica: obtenha registros do GenBank ou organize loci FASTA, alinhe, teste a saturação, concatene e construa árvores.';

  @override
  String get splashStart => 'Iniciar';

  @override
  String get splashStep1 =>
      'Carregue genomas ou marcadores individuais (COI, 12S, CYTB…) do NCBI, ou arquivos GenBank e FASTA locais. O [[SynGenes]] monta as buscas por marcador com todos os sinônimos de cada gene e padroniza os nomes dos marcadores para as próximas etapas.';

  @override
  String get splashStep2 =>
      'Inspecione os registros importados, ajuste nomes de espécies e recupere a taxonomia do GBIF, NCBI ou SALVE com [[dataFishing]] antes da seleção de marcadores.';

  @override
  String get splashStep3 =>
      'Escolha quais tipos de feature (CDS, rRNA, tRNA) entram na análise, ou indique ao SPLACE onde os cabeçalhos FASTA guardam os nomes de amostra e locus.';

  @override
  String get splashStep4 =>
      'Escolha os marcadores da matriz e revise sequências ausentes e duplicadas nos mapas de presença e comprimento dos genes.';

  @override
  String get splashStep5 =>
      'Alinhe cada marcador com [[MAFFT]] e trime com [[trimAl]] ou [[ClipKIT]], mantendo códons inteiros nos loci codificantes.';

  @override
  String get splashStep6 =>
      'Teste a saturação de substituições por posição do códon com o índice de [[Iss (saturation)|Xia]] e decida quais posições manter.';

  @override
  String get splashStep7 =>
      'Concatene os loci alinhados numa supermatriz e exporte arquivos NEXUS, PHYLIP e de partições, opcionalmente divididos por posição do códon.';

  @override
  String get splashStep8 =>
      'Estime a filogenia de máxima verossimilhança com [[IQ-TREE 3|IQ-TREE]] ou [[RAxML-NG]], com suporte dos ramos e árvores por gene opcionais.';

  @override
  String get step1 => 'Carregar Dados';

  @override
  String get step2 => 'Registros Carregados';

  @override
  String get step3 => 'Tipos de Feature';

  @override
  String get step4 => 'Selecionar Marcadores';

  @override
  String get step5 => 'Alinhamento de Sequências';

  @override
  String get step6 => 'Teste de Saturação';

  @override
  String get step7 => 'Concatenar e NEXUS';

  @override
  String get step8 => 'Inferência Filogenética';

  @override
  String get step1Intro =>
      'Escolha como carregar dados no SPLACE: registros GenBank do NCBI ou de arquivos locais, ou seus próprios arquivos FASTA.';

  @override
  String get tipsKicker => 'Dicas';

  @override
  String get tipsTitle => 'Antes de começar';

  @override
  String get tipsIntro =>
      'Monte a **tabela de taxonomia** gradualmente e combine **fontes complementares** quando necessário.';

  @override
  String get tipsExamplePrefix => 'Cada opção abaixo tem um botão';

  @override
  String get tipsExampleSuffix => 'para ajudar você a começar rapidamente.';

  @override
  String get tipSampleTitle => 'Um arquivo por amostra';

  @override
  String get tipSampleBody =>
      'cada cabeçalho é um **locus** (ex.: NC_030596_Oryza.fasta com >psbA, >matK). O SPLACE reagrupa em **um arquivo por locus**.';

  @override
  String get tipLocusTitle => 'Um arquivo por locus';

  @override
  String get tipLocusBody =>
      'cada cabeçalho é uma **amostra** (ex.: COI.fasta). Já está no formato certo.';

  @override
  String get tipMixedTitle => 'Cabeçalhos mistos';

  @override
  String get tipMixedBody =>
      'amostra e locus no **mesmo cabeçalho** (>Bufo_bufo_COI) ou **tags do NCBI** ([gene=COX1]).';

  @override
  String get tipAlignedTitle => 'Já alinhado?';

  @override
  String get tipAlignedBody =>
      'Você pode **pular o MAFFT** na etapa 5 e ir direto ao trimming, saturação e árvores.';

  @override
  String optionN(String n) {
    return 'Opção $n';
  }

  @override
  String get filesTitle => 'Usar arquivos FASTA locais';

  @override
  String get filesInfo =>
      'Importe arquivos .fasta, .fas, .fa ou .fna que já estão no seu computador.';

  @override
  String get filesBody =>
      'Importe seus próprios arquivos .fasta, .fas, .fa ou .fna do disco.';

  @override
  String get filesDrop => 'Arraste e solte arquivos FASTA ou pastas aqui';

  @override
  String get filesBrowse => 'ou clique para procurar (.fasta, .fas, .fa, .fna)';

  @override
  String get folderTitle => 'Importar uma pasta inteira';

  @override
  String get folderInfo =>
      'Lê todos os arquivos FASTA da pasta e de suas subpastas.';

  @override
  String get folderBody => 'Útil para centenas de arquivos, um por amostra.';

  @override
  String get addFolder => 'Escolher pasta';

  @override
  String inputSummary(String files, String sequences) {
    return '$files arquivos · $sequences sequências carregadas';
  }

  @override
  String get clearRecords => 'Limpar registros';

  @override
  String get step2Intro =>
      'Todos os cabeçalhos FASTA estão listados abaixo com a amostra e o locus que o SPLACE lê de cada um. Ajuste as regras na etapa 3 até todos os registros terem ambos.';

  @override
  String get detectedLayout => 'Formato detectado';

  @override
  String get layoutSamplePerFile =>
      'um arquivo por amostra — cabeçalhos são loci';

  @override
  String get layoutLocusPerFile =>
      'um arquivo por locus — cabeçalhos são amostras';

  @override
  String get layoutMixedHeader => 'amostra e locus estão no cabeçalho';

  @override
  String get layoutTaggedHeader =>
      'locus lido das tags do cabeçalho ([gene=…])';

  @override
  String get layoutUnknown => 'não reconhecido — escolha os campos manualmente';

  @override
  String confidence(String percent) {
    return 'confiança $percent%';
  }

  @override
  String get organizing => 'Atualizando…';

  @override
  String get tokenColumns =>
      'Partes do cabeçalho (separadas pelos delimitadores do locus)';

  @override
  String tokenColumnInfo(String distinct, String geneLike) {
    return '$distinct distintos · $geneLike% parecem loci';
  }

  @override
  String get searchHeaders => 'Filtrar cabeçalhos ou arquivos';

  @override
  String get colFile => 'Arquivo';

  @override
  String get colHeader => 'Cabeçalho';

  @override
  String get colSample => 'Amostra';

  @override
  String get colLocus => 'Locus';

  @override
  String get colLength => 'Tamanho';

  @override
  String get step3Intro =>
      'Escolha de onde vêm o nome da amostra e o nome do locus. A prévia e a tabela da etapa 2 se atualizam enquanto você digita.';

  @override
  String fieldN(String n) {
    return 'Campo $n';
  }

  @override
  String get sampleField => 'Nome da amostra';

  @override
  String get sampleFieldInfo =>
      'Vira o nome da sequência em todos os arquivos de locus; precisa ser idêntico entre loci para a concatenação.';

  @override
  String get locusField => 'Nome do locus';

  @override
  String get locusFieldInfo =>
      'Define para qual arquivo cada sequência vai. Sinônimos como COX1/COI/CO1 podem ser unidos.';

  @override
  String get sourceFileName => 'Nome do arquivo';

  @override
  String get sourceHeader => 'Partes do cabeçalho';

  @override
  String get sourceTag => 'Tag';

  @override
  String get sourceRegex => 'Regex';

  @override
  String get delimiters => 'Delimitadores';

  @override
  String get tokenStart => 'Primeira parte';

  @override
  String get tokenEnd => 'Fim (vazio = última)';

  @override
  String get tagKey => 'Chave da tag';

  @override
  String get regexPattern => 'Padrão (usa o primeiro grupo)';

  @override
  String get regexOnFileName => 'Aplicar ao nome do arquivo';

  @override
  String get previewLabel => 'Prévia';

  @override
  String get optionsKicker => 'Opções';

  @override
  String get duplicates => 'Duplicatas e nomes';

  @override
  String get duplicatesInfo =>
      'O que fazer quando a mesma amostra tem o mesmo locus mais de uma vez (ex.: cópias da região invertida repetida em plastomas).';

  @override
  String get dupLongest => 'Manter a mais longa';

  @override
  String get dupMostComplete => 'Manter a com menos N/gaps';

  @override
  String get dupFirst => 'Manter a primeira';

  @override
  String get dupAll => 'Manter todas (renomear)';

  @override
  String get mergeSynonyms => 'Unir sinônimos de locus (COX1 = COI = CO1…)';

  @override
  String get sanitizeNames => 'Nomes seguros para programas de árvore';

  @override
  String get selectAll => 'Selecionar Todos';

  @override
  String get selectNone => 'Nenhum';

  @override
  String get sharedOnly => 'Apenas Compartilhados';

  @override
  String get legendSelectAll => 'seleciona todos os marcadores visíveis';

  @override
  String get legendNone => 'limpa todas as seleções';

  @override
  String get legendShared => 'marcadores presentes em todas as amostras';

  @override
  String get missingTitle => 'Marcadores ausentes';

  @override
  String missingIn(String n) {
    return 'ausente em $n amostra(s)';
  }

  @override
  String get duplicatesTitle => 'Marcadores duplicados';

  @override
  String duplicatedIn(String n) {
    return 'duplicado em $n amostra(s)';
  }

  @override
  String get colCopies => 'Cópias (tamanhos)';

  @override
  String get colUsed => 'Usada';

  @override
  String get unassignedTitle => 'Registros sem amostra ou locus';

  @override
  String get missingSample => 'sem amostra';

  @override
  String get missingLocus => 'sem locus';

  @override
  String useLoci(String n) {
    return 'Usar $n marcadores';
  }

  @override
  String get exportLoci => 'Exportar FASTA por locus…';

  @override
  String get heatmapTitle => 'Mapa de Presença de Genes';

  @override
  String get heatmapAbsent => 'Ausente';

  @override
  String get heatmapPresent => 'Presente';

  @override
  String matrixSummary(String loci, String samples, String occupancy) {
    return '$loci marcadores × $samples amostras · $occupancy% preenchida. Clique no nome de um marcador para incluí-lo ou excluí-lo.';
  }

  @override
  String get panelMafft => 'Configurações do MAFFT';

  @override
  String get panelMafftBody =>
      'Cada marcador é alinhado separadamente; vários marcadores rodam em paralelo.';

  @override
  String get mafftStrategy => 'Método de Alinhamento';

  @override
  String get directionLabel => 'Direção das sequências';

  @override
  String get adjustDirection => 'Corrigir sequências invertidas';

  @override
  String get panelTrimal => 'Trimming (trimAl / ClipKIT)';

  @override
  String get trimalMethod => 'Modo de trimagem';

  @override
  String get sequenceType => 'Tipo de Sequência';

  @override
  String get codonAwareTrimming => 'Manter o quadro de leitura';

  @override
  String get codonAwareHint =>
      'Remover colunas isoladas de alinhamentos codificantes desloca o quadro de leitura e invalida as partições por códon e o teste de saturação.';

  @override
  String get panelRun => 'Executar alinhamento';

  @override
  String panelRunBody(String dir) {
    return 'Alinhe todos os marcadores selecionados com o MAFFT. A saída é salva em $dir.';
  }

  @override
  String get runMafft => 'Alinhar (MAFFT)';

  @override
  String get runTrimal => 'Trimar (trimAl)';

  @override
  String get useAsAlignment => 'A entrada já está alinhada';

  @override
  String get colCoding => 'Codificante';

  @override
  String get colFrame => 'Quadro de leitura';

  @override
  String get colInclude => 'Usar este marcador';

  @override
  String get statRaw => 'Bruto';

  @override
  String get statAligned => 'Alinhado';

  @override
  String get statTrimmed => 'Trimado';

  @override
  String get statusError => 'Erro';

  @override
  String get saturationIntro =>
      'Cada marcador é testado com o PhyKIT (Philippe et al. 2011): o IQ-TREE infere a árvore do gene e, para a 1ª, 2ª e 3ª posições do códon (e 1ª+2ª, todos os sítios), o PhyKIT faz a regressão da distância p não corrigida de cada par de táxons sobre a distância entre eles na árvore. Um valor perto de 1 significa sem saturação; quanto mais perto de 0, mais as diferenças observadas escondem substituições repetidas.';

  @override
  String get saturationThresholds =>
      'Limiares heurísticos, definidos em 6a (parâmetros do teste): moderada abaixo do limite moderado do valor de saturação do PhyKIT, saturada abaixo do limite de saturação; divergência baixa demais (p médio) conta como OK.';

  @override
  String get runSaturation => 'Executar teste de saturação';

  @override
  String get recommendation => 'Recomendação';

  @override
  String get recKeepAll =>
      'Mantenha as três posições do códon. Na etapa 7b, particione a matriz por posição do códon para que cada uma tenha seu próprio modelo de substituição.';

  @override
  String get recAminoAcids =>
      'As três posições estão saturadas: nessa divergência os nucleotídeos dizem pouco sobre as relações. Alinhe como aminoácidos (etapa 5b) ou use marcadores mais lentos.';

  @override
  String get recUndetermined =>
      'Há poucos sítios ou pares de sequências para avaliar a saturação. Mantenha todas as posições.';

  @override
  String get applyPerLocus => 'Usar o resultado de cada gene';

  @override
  String get concatenated => 'Concatenado';

  @override
  String get saturationHeatmapTitle => 'Saturação por locus e posição do códon';

  @override
  String get selectLocusHint =>
      'Passe o mouse numa célula para ver os valores; clique numa linha para ver os gráficos.';

  @override
  String get verdictOk => 'OK';

  @override
  String get verdictModerate => 'Moderada';

  @override
  String get verdictSaturated => 'Saturada';

  @override
  String get verdictInsufficient => 'Poucos dados';

  @override
  String plotObservedVsCorrected(String model) {
    return 'p observada × distância $model';
  }

  @override
  String plotTsTv(String model) {
    return 'Transições e transversões × distância $model';
  }

  @override
  String get transitions => 'Transições';

  @override
  String get transversions => 'Transversões';

  @override
  String get metricSites => 'Sítios';

  @override
  String get metricPairs => 'Pares';

  @override
  String get metricMeanP => 'p média';

  @override
  String get metricSlope => 'Inclinação';

  @override
  String get metricTsTv => 'Ts/Tv';

  @override
  String get metricIss => 'Iss';

  @override
  String get metricUndefined => 'Sem correção';

  @override
  String get metricVerdict => 'Resultado';

  @override
  String frameLabel(String frame) {
    return 'códons começam na base $frame';
  }

  @override
  String get step7Intro =>
      'Monte uma supermatriz a partir dos genes alinhados e exporte para a análise filogenética.';

  @override
  String get srcTrimmed => 'Sequências trimadas';

  @override
  String get recommendedTag => '(recomendado)';

  @override
  String get srcTrimmedBody =>
      'Usa a saída do trimAl — remove posições mal alinhadas.';

  @override
  String get srcAligned => 'Sequências alinhadas';

  @override
  String get srcAlignedBody => 'Usa a saída do MAFFT sem trimming.';

  @override
  String get partitionScheme => 'Partições';

  @override
  String get schemeNone => 'Nenhuma';

  @override
  String get schemeByLocus => 'Por locus';

  @override
  String get schemeByLocusCodon => 'Por locus e posição do códon';

  @override
  String get schemeByCodon => 'Por posição do códon (todos os loci)';

  @override
  String get codonExclusion => 'Posições do códon';

  @override
  String get codonExclusionInfo =>
      'As posições excluídas são removidas só dos loci codificantes. O teste de saturação (etapa 6) pode preencher isto.';

  @override
  String excludePosition(String n) {
    return 'Excluir posição $n';
  }

  @override
  String perLocusExclusions(String list) {
    return 'Por locus: $list';
  }

  @override
  String get clearExclusions => 'Limpar exclusões por locus';

  @override
  String get buildMatrix => 'Gerar supermatriz';

  @override
  String get matrixReady => 'Supermatriz pronta';

  @override
  String matrixStats(
    String taxa,
    String sites,
    String partitions,
    String missing,
  ) {
    return '$taxa táxons × $sites sítios · $partitions partições · $missing% faltantes';
  }

  @override
  String get excludedLabel => 'Excluído';

  @override
  String get outputFiles => 'Arquivos gerados';

  @override
  String get openFolder => 'Abrir pasta';

  @override
  String get partitionPreview => 'Partições RAxML-NG';

  @override
  String get step8Intro =>
      'Estimativa de árvore por máxima verossimilhança a partir da supermatriz concatenada com IQ-TREE3 ou RAxML-NG.';

  @override
  String get outgroup => 'Grupo externo';

  @override
  String get outgroupInfo =>
      'As espécies selecionadas são gravadas no NEXUS e usadas para enraizar as árvores do IQ-TREE e do RAxML-NG.';

  @override
  String get outgroupHint => 'Buscar espécies…';

  @override
  String engineN(String n) {
    return 'Programa $n';
  }

  @override
  String get iqtreeInfo =>
      'O ModelFinder escolhe os modelos (e pode unir partições); bootstrap ultrarrápido e SH-aLRT dão os suportes.';

  @override
  String get raxmlInfo =>
      'O RAxML-NG 2 escolhe os modelos sozinho (MOOSE), ajusta a busca à dificuldade do alinhamento e oferece medidas de suporte rápidas. Assim como o IQ-TREE, é um programa de máxima verossimilhança (não bayesiano).';

  @override
  String get model => 'Modelo de Substituição';

  @override
  String get ufboot => 'Bootstrap ultrarrápido (-B)';

  @override
  String get alrt => 'SH-aLRT (--alrt)';

  @override
  String get runIqtree => 'Executar IQ-TREE3';

  @override
  String get raxmlModelFinder =>
      'ModelFinder (IQ-TREE, modelos compatíveis com RAxML)';

  @override
  String get bootstraps => 'Árvores de bootstrap';

  @override
  String get supportMetrics => 'Suporte';

  @override
  String get seed => 'Semente';

  @override
  String get runRaxml => 'Executar RAxML-NG';

  @override
  String get runs => 'Execuções';

  @override
  String get noRuns => 'Nenhuma árvore ainda.';

  @override
  String toolMissing(String tool) {
    return '$tool não encontrado';
  }

  @override
  String get toolsTitle => 'Ferramentas externas e configurações de análise';

  @override
  String get toolsBody =>
      'O SPLACE procura os binários embutidos, depois o PATH, depois o WSL. Defina um caminho se necessário.';

  @override
  String get tools => 'Ferramentas';

  @override
  String get redetect => 'Detectar novamente';

  @override
  String get notFound => 'Não encontrado';

  @override
  String get overridePath => 'Caminho personalizado (ou wsl:/caminho/no/wsl)';

  @override
  String get browse => 'Procurar…';

  @override
  String get analysisSettings => 'Análise';

  @override
  String threads(String n) {
    return 'Threads: $n';
  }

  @override
  String get geneticCode => 'Código genético';

  @override
  String get projectFolder => 'Pasta do projeto';

  @override
  String get change => 'Alterar…';

  @override
  String get tipComplete =>
      'Prefira **genomas completos** sempre que possível, porque eles costumam melhorar as análises seguintes.';

  @override
  String get tipMultipleSearches =>
      'Faça **várias buscas taxonômicas** para adicionar **grupos internos** e **grupos externos**. Os novos resultados são acumulados na mesma tabela de taxonomia.';

  @override
  String get tipCombine =>
      'Você pode **combinar** buscas taxonômicas, importação por acessos e arquivos locais **.gb, .gbk ou .genbank** no mesmo conjunto de dados.';

  @override
  String get example => 'Exemplo';

  @override
  String get taxTitle => 'Buscar por nome taxonômico';

  @override
  String get taxInfo =>
      'Pesquise no NCBI por nome taxonômico, tipo de genoma, banco e completude antes de importar os genomas escolhidos.';

  @override
  String get taxBody =>
      'Pesquise no NCBI e baixe automaticamente os registros GenBank correspondentes para a ferramenta.';

  @override
  String get taxExampleTooltip =>
      'Preencher a busca taxonômica com Carangidae, mitocôndria, genomas completos e somente RefSeq.';

  @override
  String get taxApiHint =>
      'Adicione uma chave de API do NCBI para acelerar as importações (botão da chave na barra superior).';

  @override
  String get taxName => 'Nome taxonômico';

  @override
  String get taxOrganelle => 'Tipo de genoma';

  @override
  String get organelleMito => 'Mitocôndria';

  @override
  String get organelleChloro => 'Cloroplasto';

  @override
  String get taxScope => 'Escopo do genoma';

  @override
  String get taxCompleteOnly => 'Apenas completos';

  @override
  String get taxScopeInfo =>
      'Apenas completos busca no NCBI os registros com “complete genome” no título. Desligado, também inclui os registros com “partial genome” no título.';

  @override
  String get scopeComplete => 'Genomas completos';

  @override
  String get scopePartial => 'Genomas parciais';

  @override
  String get taxDatabase => 'Banco';

  @override
  String get taxRefseqOnly => 'Somente RefSeq';

  @override
  String get taxQueryPreview => 'Prévia da query';

  @override
  String get search => 'Buscar';

  @override
  String get accTitle => 'Colar números de acesso';

  @override
  String get accInfo =>
      'Cole um ou mais acessos conhecidos e busque exatamente esses registros GenBank no NCBI.';

  @override
  String get accBody =>
      'Informe os números de acesso separados por vírgulas, espaços ou quebras de linha.';

  @override
  String get accExampleTooltip =>
      'Carregar acessos de exemplo de Balcázar et al. (2026)';

  @override
  String get accFetchTooltip =>
      'Buscar registros GenBank no NCBI para os acessos informados';

  @override
  String get fetch => 'Buscar';

  @override
  String get gbFilesTitle => 'Usar arquivos GenBank locais';

  @override
  String get gbFilesInfo =>
      'Importe arquivos locais .gb, .gbk ou .genbank quando os registros já estiverem no seu computador.';

  @override
  String get gbFilesBody =>
      'Importe seus próprios arquivos .gb, .gbk ou .genbank do disco.';

  @override
  String get gbFilesDrop => 'Arraste e solte arquivos ou pastas GenBank aqui';

  @override
  String get gbFilesBrowse => 'ou clique para procurar (.gb, .gbk, .genbank)';

  @override
  String gbSummary(String n) {
    return '$n registros GenBank carregados';
  }

  @override
  String get clearFasta => 'Limpar arquivos FASTA';

  @override
  String get close => 'Fechar';

  @override
  String get apiKeyTitle => 'API key do NCBI';

  @override
  String get apiKeyHelper =>
      'Usada nas Opções 1 e 2 para requisições ao NCBI. Sem API key, o SPLACE roda em cerca de 3 requisições por segundo. Com uma chave válida salva, o SPLACE usa até 10 requisições por segundo.';

  @override
  String get apiKeySave => 'Salvar chave';

  @override
  String get apiKeyClear => 'Apagar chave';

  @override
  String get apiKeyGetLink =>
      'Obter uma API key nas configurações da conta NCBI';

  @override
  String get apiKeySaved =>
      'API key salva. As Opções 1 e 2 agora usam até 10 requisições por segundo.';

  @override
  String get apiKeyInvalid =>
      'A API key parece inválida. Verifique antes de salvar.';

  @override
  String get apiKeyNone =>
      'Nenhuma API key salva. As Opções 1 e 2 estão usando até 3 requisições por segundo.';

  @override
  String get apiClearTitle => 'Remover a API key salva?';

  @override
  String get apiClearBody =>
      'As Opções 1 e 2 voltarão para a taxa padrão de requisições ao NCBI.';

  @override
  String get reviewTitle => 'Revisar genomas antes de importar';

  @override
  String get reviewSubtitle =>
      'Escolha quais genomas devem ser importados. Todos os switches começam ativados.';

  @override
  String get reviewDetailTaxon => 'Consulta taxonômica';

  @override
  String get reviewDetailGenomeType => 'Tipo de genoma';

  @override
  String get reviewDetailSearchTypes => 'Tipos de busca';

  @override
  String get reviewDetailRefseq => 'RefSeq';

  @override
  String get reviewRefseqOn => 'Usado';

  @override
  String get reviewRefseqOff => 'Não usado';

  @override
  String get reviewSelectAll => 'Ativar todos';

  @override
  String get reviewSelectNone => 'Desativar todos';

  @override
  String get reviewDisableUnverified => 'Desativar UNVERIFIED';

  @override
  String get reviewOnlyComplete => 'Somente completos';

  @override
  String get reviewOnlyPartial => 'Somente parciais';

  @override
  String get reviewColImport => 'Importar';

  @override
  String get reviewColAccession => 'Acesso';

  @override
  String get reviewColDescription => 'Descrição do genoma';

  @override
  String get reviewColScope => 'Classificação';

  @override
  String get reviewImport => 'Importar genomas selecionados';

  @override
  String get reviewTotal => 'Total';

  @override
  String get reviewSelected => 'Selecionados';

  @override
  String get reviewComplete => 'Completo';

  @override
  String get reviewPartial => 'Parcial';

  @override
  String get reviewMixed => 'Título misto';

  @override
  String get reviewUnknown => 'Não classificado';

  @override
  String get reviewUnverified => 'UNVERIFIED';

  @override
  String reviewUnverifiedCount(String count) {
    return 'Genomas UNVERIFIED encontrados: $count';
  }

  @override
  String get reviewUnverifiedAdvice =>
      'Recomendamos remover os genomas UNVERIFIED antes da importação.';

  @override
  String get reviewCancelled => 'Seleção de genomas cancelada.';

  @override
  String get reviewSearch => 'Filtrar por acesso, espécie ou descrição';

  @override
  String reviewShowing(String shown, String total) {
    return 'Mostrando $shown/$total genomas que correspondem à busca.';
  }

  @override
  String get gbRecordsHeading => 'Registros GenBank';

  @override
  String get fastaHeadersHeading => 'Cabeçalhos FASTA';

  @override
  String get clearAll => 'Limpar tudo';

  @override
  String get clearTitle => 'Limpar todos os registros?';

  @override
  String get clearBody =>
      'Esta ação não pode ser desfeita. Todos os registros e seleções serão removidos.';

  @override
  String get recordsReview =>
      'Revise e corrija os nomes das espécies abaixo usando o botão de edição.';

  @override
  String get recordsReviewMethod =>
      'Revise e corrija os nomes das espécies abaixo usando o botão de edição e depois obtenha a taxonomia com o método do [[dataFishing]].';

  @override
  String get recordsFetchPrefix => 'Depois clique em';

  @override
  String get recordsFetchSuffix =>
      'para obter os ranks taxonômicos do GBIF via';

  @override
  String get internetRequired => 'Internet necessária';

  @override
  String get fetchTaxonomy => 'Buscar taxonomia';

  @override
  String get skipTaxonomy => 'Pular taxonomia e continuar';

  @override
  String get skipTaxBody1 =>
      'Se nomes de espécies contiverem caracteres especiais ou houver entradas duplicadas, podem ocorrer erros nas próximas análises.';

  @override
  String get skipTaxBody2 =>
      'Recomenda-se corrigir os nomes e buscar a taxonomia para melhores resultados.';

  @override
  String get taxonomyResolved => 'Taxonomia pronta';

  @override
  String get taxonomyTree => 'Árvore taxonômica';

  @override
  String get taxonomyTreeBody =>
      'Resumo taxonômico indentado com os mesmos registros mostrados abaixo.';

  @override
  String get columnsLabel => 'Colunas:';

  @override
  String get colAccession => 'Acesso';

  @override
  String get colKingdom => 'Reino';

  @override
  String get colPhylum => 'Filo';

  @override
  String get colClass => 'Classe';

  @override
  String get colOrder => 'Ordem';

  @override
  String get colFamily => 'Família';

  @override
  String get colGenus => 'Gênero';

  @override
  String get colSpecies => 'Espécie';

  @override
  String get colAuthorship => 'Autoria';

  @override
  String get colSource => 'Origem';

  @override
  String get sourceFile => 'Arquivo';

  @override
  String get recordsSearch => 'Buscar acesso, espécie, família...';

  @override
  String editRecordTooltip(String species, String accession) {
    return 'Editar $species ($accession)';
  }

  @override
  String removeRecordTooltip(String species, String accession) {
    return 'Remover $species ($accession)';
  }

  @override
  String get removeTitle => 'Remover registro?';

  @override
  String get removeBody => 'Esta ação não pode ser desfeita.';

  @override
  String get remove => 'Remover';

  @override
  String get editTitle => 'Editar Registro';

  @override
  String get editHint =>
      'Revise o registro atual e atualize abaixo os campos de espécie e taxonomia.';

  @override
  String get editGbif => 'Buscar Taxonomia no GBIF';

  @override
  String get save => 'Salvar';

  @override
  String get autoLabel => '(automático)';

  @override
  String get fileSize => 'Tamanho';

  @override
  String get sequenceLabel => 'Sequência';

  @override
  String get featuresLabel => 'Features';

  @override
  String get lineageLabel => 'Linhagem';

  @override
  String get speciesRequired => 'O nome da espécie é obrigatório';

  @override
  String get speciesLetters => 'Apenas letras, espaços e hífens são permitidos';

  @override
  String get speciesBinomial => 'Use gênero e epíteto (ex.: Rhinella marina)';

  @override
  String get speciesGenusUpper => 'O gênero deve começar com letra maiúscula';

  @override
  String get familySuffix => 'Nomes de família terminam em -idae ou -aceae';

  @override
  String get gbifNoMatch => 'Nenhuma correspondência no GBIF para este nome.';

  @override
  String get gbifMatched => 'Correspondência encontrada no GBIF.';

  @override
  String get step3Fasta => 'Campos de Amostra e Locus';

  @override
  String get featureTypesIntro =>
      'Selecione as categorias de genes que deseja incluir na análise. Ative ou desative cada categoria com o botão.';

  @override
  String get standardizedPrefix =>
      'Os nomes dos genes foram padronizados usando';

  @override
  String get featureCdsInfo =>
      'Sequências codificantes de proteínas (genes traduzidos em proteínas)';

  @override
  String get featureRrnaInfo => 'Genes de RNA ribossomal (12S, 16S etc.)';

  @override
  String get featureTrnaInfo => 'Genes de RNA transportador';

  @override
  String get featureTypesHint =>
      'Ative pelo menos uma categoria para listar os marcadores na etapa 4.';

  @override
  String unknownGenesButton(String n) {
    return '$n nomes desconhecidos';
  }

  @override
  String get unknownTitle =>
      'Genes Não Encontrados no SynGenes (Rabelo et al. 2024)';

  @override
  String get unknownSubtitle =>
      'Selecione o nome padronizado correto para cada gene encontrado nos seus arquivos.';

  @override
  String unknownCount(String count) {
    return '$count gene(s) não puderam ser padronizados automaticamente.';
  }

  @override
  String get unknownHowto =>
      'Selecione o nome padronizado no menu suspenso da coluna Nome Correto. Escolha (ignorar) para excluir o gene da análise.';

  @override
  String get unknownSkip => '(ignorar)';

  @override
  String geneSuggested(String name) {
    return 'Sugestão: $name';
  }

  @override
  String geneFromProduct(String value) {
    return 'pelo /product “$value”';
  }

  @override
  String geneFromGene(String value) {
    return 'pelo /gene “$value”';
  }

  @override
  String geneFromSequence(String percent, String name, String n) {
    return 'sequência $percent% similar a $name em $n outros genomas';
  }

  @override
  String geneSequenceAgrees(String percent) {
    return 'a sequência confirma ($percent%)';
  }

  @override
  String geneSequenceDisagrees(String name, String percent) {
    return 'mas a sequência parece $name ($percent%)';
  }

  @override
  String get geneNoSuggestion =>
      'Sem correspondência no SynGenes nem nos outros genomas; confira a anotação.';

  @override
  String geneQualifierMismatch(String gene, String product) {
    return '/gene “$gene” e /product “$product” são diferentes';
  }

  @override
  String conflictsBanner(String n) {
    return '$n feature(s) têm /gene e /product indicando genes diferentes.';
  }

  @override
  String get conflictsReview => 'Revisar';

  @override
  String get conflictsTitle => '/gene e /product não conferem';

  @override
  String get conflictsSubtitle =>
      'O SPLACE usa o nome do /product. Cada sequência é comparada com o mesmo marcador nos outros genomas para mostrar qual qualificador está certo.';

  @override
  String conflictUsed(String name) {
    return 'Usado: $name';
  }

  @override
  String conflictSeqGene(String name, String percent) {
    return 'A sequência corresponde a $name ($percent%): o /gene parece certo';
  }

  @override
  String conflictSeqUsed(String name, String percent) {
    return 'A sequência corresponde a $name ($percent%): o nome usado parece certo';
  }

  @override
  String conflictSeqOther(String name, String percent) {
    return 'A sequência corresponde a $name ($percent%), nenhum dos dois';
  }

  @override
  String get conflictSeqNone => 'Sequência inconclusiva';

  @override
  String conflictUseGene(String name) {
    return 'Usar $name';
  }

  @override
  String get speciesIdRunning => 'Identificando espécies';

  @override
  String speciesIdGenesUsed(String genes) {
    return 'BLAST no GenBank e BOLD · $genes';
  }

  @override
  String get speciesIdFooter =>
      'O BLAST no NCBI costuma levar de 1 a 3 minutos; o BOLD responde em segundos.';

  @override
  String get speciesIdButton => 'Identificar espécies (GenBank + BOLD)';

  @override
  String get speciesIdShowResults => 'Resultado da identificação';

  @override
  String get speciesIdSetupTitle => 'Identificar espécies com GenBank e BOLD';

  @override
  String get speciesIdSetupSubtitle =>
      'Compare os genes de cada genoma com o GenBank (BLAST) e o BOLD para conferir os nomes das espécies antes de obter a taxonomia.';

  @override
  String get speciesIdGenesLabel => 'Genes para comparar';

  @override
  String get speciesIdGenesHint =>
      'Escolha um ou mais genes. Cada gene é buscado no GenBank com BLAST; o BOLD aceita apenas o COI.';

  @override
  String speciesIdInGenomes(String n) {
    return 'em $n genomas';
  }

  @override
  String get speciesIdSources => 'Bancos de dados';

  @override
  String get speciesIdGenbankBody =>
      'BLAST (megablast) contra o core_nt, para cada gene escolhido.';

  @override
  String get speciesIdBoldBody =>
      'ID engine do BOLD, biblioteca pública de espécies, apenas COI.';

  @override
  String get speciesIdBoldOff => 'Escolha o COI para buscar também no BOLD.';

  @override
  String get speciesIdRun => 'Identificar';

  @override
  String get speciesIdResultsTitle => 'Identificação das espécies';

  @override
  String speciesIdResultsSubtitle(String genes) {
    return 'BLAST no GenBank e BOLD · $genes. Escolha o nome a manter para cada genoma.';
  }

  @override
  String get speciesIdOnlyDisagree => 'Só possíveis identificações erradas';

  @override
  String speciesIdSummary(String disagree, String total) {
    return '$disagree de $total genomas podem estar com identificação errada';
  }

  @override
  String get speciesIdColGenome => 'Genoma';

  @override
  String get speciesIdColCandidates => 'Espécies encontradas (escolha uma)';

  @override
  String get speciesIdMismatch => 'Possível identificação errada';

  @override
  String get speciesIdMatch => 'Confirmada';

  @override
  String get speciesIdNoHits => 'Sem correspondência no GenBank ou no BOLD';

  @override
  String speciesIdKeepCurrent(String name) {
    return 'Manter $name';
  }

  @override
  String speciesIdHits(String n) {
    return '$n hits';
  }

  @override
  String speciesIdApply(String n) {
    return 'Aplicar $n alteração(ões)';
  }

  @override
  String get speciesIdNote =>
      'Hits do próprio genoma (ou do registro do GenBank de onde veio a cópia RefSeq NC_) aparecem como “este registro” e não contam como evidência.';

  @override
  String get geneReviewButton => 'Revisar nomes dos genes';

  @override
  String geneReviewButtonN(String n) {
    return 'Revisar nomes dos genes ($n)';
  }

  @override
  String get geneReviewDone => 'Nomes dos genes revisados';

  @override
  String get geneReviewTitle => 'Revisão dos nomes dos genes';

  @override
  String get geneReviewSubtitle =>
      'Nomes que o SynGenes não reconheceu e genes cujo /gene e /product não conferem. Revise antes da taxonomia.';

  @override
  String geneReviewUnknownSection(String n) {
    return 'Nomes não encontrados no SynGenes ($n)';
  }

  @override
  String geneReviewConflictSection(String n) {
    return '/gene e /product indicam genes diferentes ($n)';
  }

  @override
  String get geneReviewConflictHelp =>
      'Os dois qualificadores foram normalizados com o SynGenes. Escolha qual nome esta feature recebe; só este genoma muda.';

  @override
  String get geneReviewNothing =>
      'Todos os nomes de genes foram padronizados e cada /gene confere com seu /product.';

  @override
  String conflictPickGene(String name) {
    return '/gene → $name';
  }

  @override
  String conflictPickProduct(String name) {
    return '/product → $name';
  }

  @override
  String get grHowTitle => 'Como ler esta revisão';

  @override
  String get grHowSynGenes =>
      'O SynGenes reconhece uma tag quando ela é um nome padrão de gene ou um de seus sinônimos conhecidos.';

  @override
  String get grHowSequence =>
      'Similaridade é a fração da sequência (proteína para CDS, DNA para rRNA) que também aparece no mesmo gene dos outros genomas deste conjunto: acima de 50% é forte, de 35 a 50% moderada e abaixo de 35% fraca.';

  @override
  String get grSequence => 'Sequência';

  @override
  String grRecognized(String name) {
    return 'SynGenes: $name';
  }

  @override
  String get grNotRecognized => 'não está no SynGenes';

  @override
  String get grAbsent => 'ausente';

  @override
  String grInGenomes(String n) {
    return 'em $n genoma(s)';
  }

  @override
  String get grUseName => 'Nome a usar';

  @override
  String get grNoSequenceRef => 'nenhum outro genoma para comparar';

  @override
  String grWhyProduct(String value, String name) {
    return 'O /product “$value” é um nome de $name no SynGenes.';
  }

  @override
  String grWhyGene(String value, String name) {
    return 'O /gene “$value” é um nome de $name no SynGenes.';
  }

  @override
  String grWhySequence(String percent, String name, String n) {
    return 'Nenhuma das tags está no SynGenes, mas a sequência é $percent% similar a $name em $n outros genomas.';
  }

  @override
  String get grWhyNone =>
      'Nenhuma das tags está no SynGenes e a sequência não corresponde claramente a nenhum marcador. Confira a anotação ou ignore este nome.';

  @override
  String grSeqAgrees(String percent) {
    return 'A sequência confirma ($percent% de similaridade).';
  }

  @override
  String grSeqDisagrees(String name, String percent) {
    return 'Atenção: a sequência se parece mais com $name ($percent%).';
  }

  @override
  String grSays(String tag, String name) {
    return '$tag diz $name';
  }

  @override
  String grSimilarityTo(String name) {
    return 'Similaridade com $name nos outros genomas';
  }

  @override
  String get grSupported => 'Confirmado pela sequência';

  @override
  String get grNotSupported => 'Não confirmado pela sequência';

  @override
  String grConflictGeneRight(
    String gp,
    String gene,
    String pp,
    String product,
  ) {
    return 'A sequência é $gp% similar a $gene e $pp% a $product nos outros genomas; portanto, a tag /gene ($gene) provavelmente está certa e a tag /product está errada.';
  }

  @override
  String grConflictProductRight(
    String pp,
    String product,
    String gp,
    String gene,
  ) {
    return 'A sequência é $pp% similar a $product e $gp% a $gene nos outros genomas; portanto, a tag /product ($product) provavelmente está certa e a tag /gene está errada.';
  }

  @override
  String get grConflictUnclear =>
      'A sequência não favorece claramente nenhum dos dois genes. Confira a anotação deste genoma antes de escolher.';

  @override
  String get grConflictIntro =>
      'As tags /gene e /product desta feature indicam genes diferentes. Clique no nome a usar; só este genoma muda.';

  @override
  String get siHowTitle => 'Como ler a identificação';

  @override
  String get siHowBody =>
      'Os genes de cada genoma foram comparados com o GenBank (BLAST) e o BOLD. Identidade é a fração de bases idênticas com a sequência mais parecida encontrada.';

  @override
  String get siLegendSame => '≥ 98%: mesma espécie';

  @override
  String get siLegendClose => '95–98%: espécie próxima';

  @override
  String get siLegendOther => '< 95%: outra espécie';

  @override
  String get siColUse => 'Usar';

  @override
  String get siColSpecies => 'Espécie';

  @override
  String get siColGenbank => 'GenBank (BLAST)';

  @override
  String get siColBold => 'BOLD (COI)';

  @override
  String get siColSeqs => 'Sequências encontradas';

  @override
  String get siCurrentTag => 'nome atual';

  @override
  String get siNotFound => 'não encontrada';

  @override
  String siVerdictConfirmed(String identity) {
    return 'As correspondências confirmam o nome atual (melhor identidade $identity%).';
  }

  @override
  String siVerdictMismatch(
    String best,
    String identity,
    String current,
    String currentIdentity,
  ) {
    return 'A correspondência mais próxima é $best ($identity%), enquanto o nome atual $current chega a $currentIdentity%. O genoma pode estar identificado errado.';
  }

  @override
  String siVerdictMissing(String best, String identity, String current) {
    return 'A correspondência mais próxima é $best ($identity%) e $current não aparece entre as correspondências. O genoma pode estar identificado errado, ou o nome pode ser um sinônimo.';
  }

  @override
  String get siVerdictNone =>
      'Nenhuma sequência do GenBank ou do BOLD correspondeu aos genes deste genoma.';

  @override
  String noMarkersButton(String n) {
    return 'Excluir genomas sem marcadores ($n)';
  }

  @override
  String get noMarkersTitle => 'Genomas sem marcadores';

  @override
  String get noMarkersSubtitle =>
      'Estes genomas não têm nenhum gene codificante ou rRNA que o SPLACE possa usar, então não acrescentam nada à matriz. Confira a lista antes de excluir.';

  @override
  String get noMarkersReasonNone =>
      'Nenhum gene anotado (0 PCGs, 0 rRNAs, 0 tRNAs): o registro do GenBank não tem anotação de genes.';

  @override
  String noMarkersReasonTrna(String n) {
    return 'Só $n tRNA(s): nenhum PCG ou rRNA.';
  }

  @override
  String noMarkersReasonUnnamed(String pcgs, String rrnas) {
    return '$pcgs PCG(s) e $rrnas rRNA(s) anotados, mas nenhum tem nome em /gene ou /product para identificá-lo.';
  }

  @override
  String noMarkersRemove(String n) {
    return 'Excluir $n genoma(s)';
  }

  @override
  String get heatmapShow => 'Mostrar heatmap';

  @override
  String get heatmapHide => 'Ocultar heatmap';

  @override
  String get heatmapHoverHint =>
      'Aponte para uma célula, um marcador ou uma amostra para ver os detalhes.';

  @override
  String get headerPresets => 'Modelos prontos';

  @override
  String get headerDragHint =>
      'Arraste para reordenar. O exemplo vem do primeiro registro.';

  @override
  String get headerAddHint =>
      'Clique em um campo para adicioná-lo ao fim do cabeçalho.';

  @override
  String get headerFirstRecord => 'Primeiro registro';

  @override
  String get headerFirstRecords => 'Cabeçalhos dos primeiros registros';

  @override
  String siSourceFailed(String source, String error) {
    return '$source não respondeu, por isso a coluna dele está sem resultados: $error';
  }

  @override
  String get whenToUse => 'Quando usar';

  @override
  String get trimCodonHandling => 'Genes codificantes';

  @override
  String get trimCodonHandlingInfo =>
      'As posições do códon são usadas depois nas partições e no teste de saturação (passo 6), então genes codificantes devem manter códons inteiros. O trimAl faz isso com -backtrans (apara o alinhamento traduzido e devolve os códons; -ignorestopcodon porque ele confere stops com o código padrão) e o ClipKIT com --codon. Excluir uma posição do códon da supermatriz é decidido no passo 7b, não aqui.';

  @override
  String get alignRunning => 'Alinhando com MAFFT';

  @override
  String get trimRunning => 'Trimando os alinhamentos';

  @override
  String markersProgress(String done, String total, String mode) {
    return '$done de $total marcadores · $mode';
  }

  @override
  String get alignFooter =>
      'O MAFFT roda neste computador, vários marcadores ao mesmo tempo.';

  @override
  String get trimFooter =>
      'O trimAl/ClipKIT rodam neste computador, vários marcadores ao mesmo tempo.';

  @override
  String alignCurrent(String markers) {
    return 'Em andamento: $markers';
  }

  @override
  String get colReason => 'Motivo';

  @override
  String get unknownSkipAll => 'Ignorar Todos';

  @override
  String get unknownApply => 'Aplicar Correções';

  @override
  String get more => 'mais';

  @override
  String get defaultPhylogeny => 'Padrão para filogenia';

  @override
  String get legendDefault => 'conjunto recomendado para análise filogenética';

  @override
  String get tooltipAll => 'Selecionar todos os marcadores visíveis';

  @override
  String get tooltipNone => 'Limpar todos os marcadores selecionados';

  @override
  String get tooltipShared =>
      'Selecionar apenas genes presentes em todos os genomas (sem dados ausentes)';

  @override
  String get tooltipDefault =>
      'Selecionar o conjunto padrão de genes para análise filogenética';

  @override
  String get duplicateChoiceHelp =>
      'A supermatriz usa uma cópia por genoma: clique na cópia a usar. A mais longa começa selecionada; “mais próxima” marca a cópia com comprimento mais perto dos genomas com uma só cópia, em geral a certa quando as cópias diferem.';

  @override
  String get defaultLabel => 'padrão';

  @override
  String dupReference(
    String gene,
    String n,
    String avg,
    String min,
    String max,
  ) {
    return '$gene em $n genomas com uma única cópia: média de $avg pb (de $min a $max pb). Escolha a cópia mais próxima dela.';
  }

  @override
  String get dupNoReference =>
      'Nenhum outro genoma tem uma única cópia deste gene para comparar os tamanhos.';

  @override
  String get dupCopiesColumn => 'Cópias (clique na que será usada)';

  @override
  String dupVsAverage(String diff) {
    return '$diff pb em relação à média';
  }

  @override
  String get dupSameAsAverage => 'igual à média';

  @override
  String get dupClosest => 'mais próxima da média';

  @override
  String get dupLongestTag => 'mais longa';

  @override
  String get noMarkers =>
      'Nenhum marcador encontrado para os tipos selecionados.';

  @override
  String get featureLabel => 'Feature';

  @override
  String get recordsLabel => 'Registros';

  @override
  String heatmapSummary(String markers, String samples) {
    return '$markers marcadores × $samples amostras. Clique no nome de um marcador para selecioná-lo.';
  }

  @override
  String get lengthHeatmapTitle => 'Mapa de Comprimento dos Marcadores';

  @override
  String get lengthHeatmapBody =>
      'Cada marcador selecionado é escalado de forma independente pelo seu comprimento mínimo e máximo.';

  @override
  String get panelHeader => 'Configurar Cabeçalho FASTA';

  @override
  String get panelHeaderBody =>
      'Identificador das sequências nos arquivos de alinhamento. O padrão (acesso_família_gênero_espécie) é usado, a menos que você o altere no construtor de cabeçalho.';

  @override
  String get openHeaderBuilder => 'Abrir construtor de cabeçalho';

  @override
  String get headerBuilderTitle => 'Configurar Cabeçalho FASTA';

  @override
  String get headerBuilderSubtitle =>
      'Escolha, ordene e separe os campos de cada nome de sequência.';

  @override
  String get headerFields => 'Campos do cabeçalho (em ordem)';

  @override
  String get headerAvailable => 'Campos disponíveis';

  @override
  String get headerSeparator => 'Separador';

  @override
  String get headerEmpty => 'Adicione pelo menos um campo.';

  @override
  String get headerConfirm => 'Confirmar cabeçalho';

  @override
  String get headerConfirmed => 'Cabeçalho confirmado';

  @override
  String get headerUnique => 'Todos os cabeçalhos são únicos.';

  @override
  String headerDuplicates(String n) {
    return '$n cabeçalho(s) se repetem em mais de um registro; sufixos numéricos (_2, _3…) serão adicionados. Inclua o Acesso para torná-los únicos.';
  }

  @override
  String get moveLeft => 'Mover para antes';

  @override
  String get moveRight => 'Mover para depois';

  @override
  String get panelSeqType => 'Tipo de sequência';

  @override
  String get seqTypeBody =>
      'Escolha como alinhar as sequências. O modo Codon-aware alinha as traduções proteicas com MAFFT e depois retrotraduz para DNA, preservando o quadro de leitura.';

  @override
  String get seqNtTitle => 'Nucleotídeo';

  @override
  String get seqNtBody => 'Alinhar sequências de DNA/RNA diretamente';

  @override
  String get seqAaTitle => 'Aminoácido';

  @override
  String get seqAaBody =>
      'Traduz (codon_start, transl_table) e alinha as proteínas';

  @override
  String get seqCodonTitle => 'Codon-aware';

  @override
  String get seqCodonBody =>
      'Alinhamento de códons guiado por proteína (apenas CDS)';

  @override
  String get recNtTitle => 'Nucleotídeo, alinhamento direto do DNA';

  @override
  String get recNtBody =>
      'Melhor para rRNAs, tRNAs, regiões não codificantes e conjuntos exploratórios. O FASTA vai direto para o MAFFT.';

  @override
  String get recAaTitle => 'Aminoácido, alinhamento de proteínas';

  @override
  String get recAaBody =>
      'Melhor quando os genes codificantes são divergentes. O DNA é traduzido, o MAFFT alinha as proteínas e o resultado é analisado como proteína.';

  @override
  String get recCodonTitle => 'Códon-aware, filogenia de CDS';

  @override
  String get recCodonBody =>
      'Recomendado para matrizes de CDS com saída em DNA e códons preservados: valida o tamanho das CDS, remove stops terminais, alinha proteínas e retrotraduz.';

  @override
  String get geneticCodeInfo =>
      'Usado na tradução, na detecção do quadro de leitura e no teste de saturação. O qualificador /transl_table do GenBank tem prioridade.';

  @override
  String get translationInfo =>
      'As sequências são traduzidas antes do alinhamento. O MAFFT roda em modo --amino.';

  @override
  String get codonWarning =>
      'O modo Codon-aware processa apenas features CDS. Genes sem CDS anotado (rRNA, tRNA etc.) serão ignorados. O alinhamento resultante preserva o quadro de leitura e suporta partição por posição de códon.';

  @override
  String get codonWorkflow => 'Fluxo códon-aware';

  @override
  String get codonStep1 => 'Extrair CDS';

  @override
  String get codonStep1Body =>
      'Apenas CDS anotadas são usadas; rRNAs e tRNAs são ignorados.';

  @override
  String get codonStep2 => 'Validar quadro';

  @override
  String get codonStep2Body =>
      'Tamanho divisível por 3 após codon_start e tratamento do stop terminal.';

  @override
  String get codonStep3 => 'Traduzir';

  @override
  String get codonStep3Body => 'Com o código genético NCBI de cada CDS.';

  @override
  String get codonStep4 => 'MAFFT proteína';

  @override
  String get codonStep4Body => 'Proteínas alinhadas com mafft --amino.';

  @override
  String get codonStep5 => 'Retrotraduzir';

  @override
  String get codonStep5Body =>
      'O alinhamento de proteínas guia um alinhamento de DNA que preserva os códons.';

  @override
  String get allowPseudo => 'Permitir CDS com /pseudo e /pseudogene';

  @override
  String get allowPseudoBody =>
      'Desativado por padrão. Quando desligado, pseudo CDS são excluídos antes da tradução.';

  @override
  String get allowStops => 'Permitir stop codons internos';

  @override
  String get allowStopsBody =>
      'Desativado por padrão. Quando ligado, sequências com stops internos permanecem na análise e são reportadas com avisos.';

  @override
  String get allowBtWarnings =>
      'Continuar quando as verificações de aminoácido e códon gerarem avisos';

  @override
  String get allowBtWarningsBody =>
      'Desativado por padrão. Desligado, qualquer divergência aminoácido/códon após o MAFFT falha o marcador; ligado, os avisos ficam no log.';

  @override
  String get pythonTitle => 'ClipKIT (corte codon-aware)';

  @override
  String get pythonBody =>
      'O ClipKIT roda com Python 3. O SPLACE o usa nos alinhamentos codon-aware: ele corta o alinhamento de proteínas e o SPLACE retrotraduz o resultado para códons.';

  @override
  String get pythonNotChecked => 'O Python ainda não foi verificado.';

  @override
  String get pythonMissing => 'Python 3 não encontrado no PATH.';

  @override
  String pythonFound(String exe, String version) {
    return '$exe (Python $version)';
  }

  @override
  String clipkitFound(String version) {
    return 'ClipKIT $version disponível';
  }

  @override
  String get clipkitMissing => 'ClipKIT não instalado';

  @override
  String get checkPython => 'Verificar Python';

  @override
  String get installClipkit => 'Instalar ClipKIT';

  @override
  String get mafftAuto => '--auto (automático, recomendado)';

  @override
  String get mafftLocalpair => '--localpair (L-INS-i, mais preciso, lento)';

  @override
  String get mafftGlobalpair => '--globalpair (G-INS-i, global, lento)';

  @override
  String get mafftGenafpair => '--genafpair (E-INS-i, gaps longos, lento)';

  @override
  String get mafftRetree2 => '--retree 2 (FFT-NS-2, rápido)';

  @override
  String get mafftRetree1 => '--retree 1 (FFT-NS-1, mais rápido)';

  @override
  String get mafftMaxIterate => 'Máximo de iterações';

  @override
  String get mafftMaxIterateHelp => '0 = sem refinamento; 1000 = máximo';

  @override
  String get mafftThreads => 'Threads';

  @override
  String get mafftThreadsHelp => '(por job)';

  @override
  String mafftConcurrency(String cpus, String jobs, String threads) {
    return '$cpus CPUs lógicas: $jobs job(s) × $threads thread(s)';
  }

  @override
  String get mafftGapOpening => 'Abertura de gap';

  @override
  String get mafftGapExtension => 'Extensão de gap';

  @override
  String get mafftExtra => 'Parâmetros adicionais';

  @override
  String get adjustDirectionHelp => 'Inverte sequências da fita oposta';

  @override
  String get reorderHelp => 'Ordena a saída por similaridade';

  @override
  String get optional => '(opcional)';

  @override
  String get optionalLabel => 'Opcional';

  @override
  String get mafftAminoNote =>
      'Os modos aminoácido e códon-aware adicionam --amino automaticamente (--adjustdirection não se aplica a proteínas).';

  @override
  String runMafftMode(String mode) {
    return 'Executar Alinhamento MAFFT ($mode)';
  }

  @override
  String markersToAlign(String n) {
    return '$n marcador(es) selecionado(s)';
  }

  @override
  String get trimPromptTitle =>
      'Alinhamento concluído - trimar regiões mal alinhadas?';

  @override
  String get trimPromptBody =>
      'O trimAl remove colunas com gaps excessivos ou baixa conservação, melhorando a qualidade da inferência filogenética.';

  @override
  String get trimPromptCodon =>
      'No modo Codon-aware, o ClipKIT corta o alinhamento de proteínas e o SPLACE retrotraduz para DNA preservando os códons.';

  @override
  String get trimPromptNo => 'Não, pular trimagem';

  @override
  String get trimPromptYes => 'Sim, configurar trimagem';

  @override
  String get trimAutomated1 => '-automated1 (heurístico, recomendado)';

  @override
  String get trimStrict => '-strict (automatizado estrito)';

  @override
  String get trimStrictplus => '-strictplus (estrito+)';

  @override
  String get trimGappyout => '-gappyout (remove colunas com muitos gaps)';

  @override
  String get trimNogaps => '-nogaps (remove todas as colunas com gaps)';

  @override
  String get trimNoallgaps => '-noallgaps (somente colunas totalmente vazias)';

  @override
  String get trimManual => 'Limiares manuais';

  @override
  String get trimGapThreshold => 'Limite de gaps';

  @override
  String get trimSimilarity => 'Similaridade';

  @override
  String get trimConsistency => 'Consistência';

  @override
  String get trimWindow => 'Tamanho da janela';

  @override
  String get trimExtra => 'Parâmetros adicionais';

  @override
  String get trimCodonNote =>
      'Nas execuções códon-aware o MAFFT alinhou proteínas: o trimmer corta o alinhamento de proteínas, o SPLACE retrotraduz para DNA preservando códons e valida o quadro.';

  @override
  String get codonTrimBackend => 'Ferramenta de trimming de códons';

  @override
  String get codonTrimAuto => 'Auto: trimAl e, se necessário, ClipKIT';

  @override
  String get codonTrimTrimal =>
      'trimAl, alinhamento de proteínas com retrotradução';

  @override
  String get codonTrimClipkit => 'ClipKIT (requer Python)';

  @override
  String get clipkitMode => 'Modo do ClipKIT';

  @override
  String get codonTrimmingTitle => 'Trimming de códons';

  @override
  String get codonTrimmingBody =>
      'Remove blocos de códons mal alinhados antes da concatenação e de qualquer filtragem por posição do códon.';

  @override
  String get trimSkippedTag => 'Trimming ignorado';

  @override
  String get trimming => 'Trimming';

  @override
  String get resultsTitle => 'Resultados do alinhamento';

  @override
  String get noTrimming => 'Sem trimagem';

  @override
  String get viewAlignment => 'Mostrar alinhamento';

  @override
  String get tabProtein => 'Proteína';

  @override
  String get viewDna => 'Visualização DNA';

  @override
  String get viewTranslate => 'Traduzir para aminoácidos';

  @override
  String possibleStops(String n) {
    return '$n possível(is) stop(s) interno(s)';
  }

  @override
  String get conservation => 'CONSERVAÇÃO';

  @override
  String get stopLegend => 'Stop';

  @override
  String get gapLegend => 'Gap';

  @override
  String get replaceGaps => 'Substituir gaps do alinhamento por ?';

  @override
  String get replaceGapsBody =>
      'Use quando o programa seguinte deve tratar gaps como dados faltantes. As coordenadas das partições são preservadas.';

  @override
  String get selectedLabel => 'Selecionadas';

  @override
  String get saveNexus => 'Salvar NEXUS';

  @override
  String get savePartition => 'Salvar tabela de partições';

  @override
  String get saveFasta => 'Salvar FASTA';

  @override
  String get savePhylip => 'Salvar PHYLIP';

  @override
  String get partitionName => 'Partição';

  @override
  String get partitionRange => 'Intervalo';

  @override
  String get partitionSites => 'Sítios';

  @override
  String get partitionType => 'Tipo';

  @override
  String get noneSelected => 'nenhum selecionado';

  @override
  String get outgroupFromStep7 => 'escolhido na etapa 7';

  @override
  String get substitutionModel => 'Modelo de Substituição';

  @override
  String get modelTest => 'ModelFinder: melhor modelo automático (TEST)';

  @override
  String get modelMfp => 'MFP: ModelFinder Plus';

  @override
  String get modelMfpMerge =>
      'MFP+MERGE: busca de modelo por partição (recomendado)';

  @override
  String get modelCustom => 'Personalizado…';

  @override
  String get bootstrapSupport => 'Suporte de Bootstrap';

  @override
  String get bootUf1000 => 'UFBoot2: 1000 réplicas (recomendado)';

  @override
  String get bootUf5000 => 'UFBoot2: 5000 réplicas';

  @override
  String get bootClassic => 'Bootstrap clássico: 100 réplicas (lento)';

  @override
  String get bootNone => 'Sem bootstrap';

  @override
  String get threadsLabel => 'Threads';

  @override
  String get extraFlags => 'Parâmetros extras';

  @override
  String get spaceSeparated => '(separados por espaço)';

  @override
  String get perGeneTrees => 'Inferir também árvores por gene';

  @override
  String get perGeneTreesBody =>
      'Executa uma análise IQ-TREE separada para cada marcador (salva em gene_trees/).';

  @override
  String get commandPreview => 'Prévia do comando';

  @override
  String get saveZip => 'Salvar resultados em ZIP';

  @override
  String get analysisComplete => 'Análise concluída';

  @override
  String get outputDirectory => 'Pasta de saída';

  @override
  String get geneTreesLabel => 'Árvores por gene';

  @override
  String get howToCite => 'Como citar';

  @override
  String get howToCiteHint =>
      'Se usar estes resultados, cite o SPLACE e as ferramentas que ele executou';

  @override
  String get copy => 'Copiar';

  @override
  String metaGenomes(String n) {
    return '$n genomas';
  }

  @override
  String metaFastaFiles(String n) {
    return '$n arquivos FASTA';
  }

  @override
  String metaSamples(String n) {
    return '$n amostras';
  }

  @override
  String metaSelected(String selected, String total) {
    return '$selected de $total selecionados';
  }

  @override
  String metaAligned(String done, String total) {
    return '$done/$total alinhados';
  }

  @override
  String metaTrimmed(String n) {
    return '$n trimados';
  }

  @override
  String metaMatrix(String taxa, String sites) {
    return '$taxa táxons × $sites sítios';
  }

  @override
  String reviewSelectedCount(String selected, String total) {
    return '$selected de $total genomas ativados';
  }

  @override
  String get progressSearchTitle => 'Buscando no NCBI por nome taxonômico';

  @override
  String progressSearchingFor(String name) {
    return 'Buscando por $name';
  }

  @override
  String get progressFetchTitle => 'Importando registros GenBank';

  @override
  String progressCounter(String current, String total) {
    return '$current de $total';
  }

  @override
  String progressSpecies(String current, String total) {
    return '$current de $total espécies';
  }

  @override
  String get progressRateLabel => 'Taxa de requisições';

  @override
  String get progressRateWithKey => 'Até 10 requisições/s';

  @override
  String get progressRateDefault => 'Até 3 requisições/s';

  @override
  String get progressApiEnabled => 'Salva e ativa';

  @override
  String get progressApiDisabled => 'Não informada';

  @override
  String get progressCurrentLabel => 'Item atual';

  @override
  String get progressCurrentSpecies => 'Espécie atual';

  @override
  String get progressStatusIdle => 'Aguardando início';

  @override
  String get progressFooter =>
      'O SPLACE está consultando o NCBI e atualizando o progresso continuamente.';

  @override
  String progressTaxonomyTitle(String source) {
    return 'Buscando taxonomia em $source';
  }

  @override
  String progressTaxonomyFooter(String source) {
    return 'O SPLACE está consultando $source e atualizando os registros conforme encontra correspondências.';
  }

  @override
  String get progressWorking => 'Processando…';

  @override
  String get colStatus => 'Categoria (SALVE)';

  @override
  String fetchTaxonomyFrom(String source) {
    return 'Buscar taxonomia ($source)';
  }

  @override
  String taxonomyMatched(String source) {
    return 'Correspondência encontrada em $source.';
  }

  @override
  String taxonomyNoMatch(String source) {
    return 'Nenhuma correspondência em $source para este nome.';
  }

  @override
  String get taxonomyCardTitle => 'Taxonomia dos registros';

  @override
  String get taxonomySourceLabel => 'Fonte da taxonomia';

  @override
  String get sourceGbifBody =>
      'Nomes aceitos e autoria do backbone do Global Biodiversity Information Facility.';

  @override
  String sourceNcbiBody(String withId, String total) {
    return 'Linhagem do taxon id de cada registro GenBank ($withId de $total registros têm; os demais são buscados pelo nome).';
  }

  @override
  String get sourceSalveBody =>
      'Só fauna brasileira. Avaliação nacional do ICMBio: taxonomia, categoria de risco de extinção e a ficha da espécie, pela API pública do SALVE.';

  @override
  String get badgeOnline => 'online';

  @override
  String get salveHint =>
      '**Somente espécies da fauna brasileira.** O SALVE reúne as avaliações nacionais de risco de extinção do ICMBio; espécies de outros lugares não têm ficha e aparecem como não encontradas (use GBIF ou NCBI para elas). Todas as seções de cada ficha são lidas online: ranques, categoria nas cores da IUCN e critérios, nomes comuns, distribuição (endemismo, biomas, estados, bacias), história natural (habitats, alimentação, reprodução, migração), população (tendência, tempo de geração), ameaças, usos, listas oficiais, histórico de avaliações, ações de conservação, unidades de conservação, pesquisa, avaliadores, citação e DOI.';

  @override
  String get poweredByDataFishing => 'Método do';

  @override
  String get metaSatKeep => 'manter todas as posições';

  @override
  String metaSatExclude(String positions) {
    return 'remover posição $positions';
  }

  @override
  String get metaSatAmino => 'usar aminoácidos';

  @override
  String get apiKeyActive => 'Ativa';

  @override
  String get apiKeyInactive => 'Não configurada';

  @override
  String get apiKeyWithout => 'Sem chave';

  @override
  String get apiKeyWith => 'Com chave';

  @override
  String get apiKeyUsedFor =>
      'Usada na busca do NCBI, no download por acessos e no NCBI Taxonomy.';

  @override
  String get apiKeyHowTitle => 'Como obter uma chave';

  @override
  String get apiKeyHow1 => 'Entre na sua conta do NCBI (gratuita).';

  @override
  String get apiKeyHow2 =>
      'Abra Account settings → API Key Management e clique em “Create an API Key”.';

  @override
  String get apiKeyHow3 => 'Copie a chave e cole no campo abaixo.';

  @override
  String get apiKeyOpenNcbi => 'Abrir conta do NCBI';

  @override
  String get apiKeyPaste => 'Colar';

  @override
  String get apiKeyFieldLabel => 'Sua API key';

  @override
  String get apiKeyLooksValid => 'Formato válido';

  @override
  String get apiKeyStored =>
      'Guardada só neste computador (Documentos/SPLACE/settings.json).';

  @override
  String get apiKeyEstimateTitle => 'Tempo estimado de download';

  @override
  String apiKeyEstimateGenomes(String n) {
    return '$n genomas';
  }

  @override
  String apiKeyEstimateDetail(String batch, String rate) {
    return '$batch genomas por requisição · até $rate requisições/s';
  }

  @override
  String apiKeyEstimateGain(String factor, String saved) {
    return '$factor× mais rápido com chave · $saved a menos';
  }

  @override
  String get apiKeyEstimateNote =>
      'Download de genomas mitocondriais que ainda não estão no cache local, incluindo a leitura, a partir de tempos medidos com o app em 08/11/2026 (Testudines, cerca de 75 kB por registro): 50, 200 e 700 genomas levaram 1,6, 4,4 e 10,1 s com chave e 5,5, 16,9 e 64,6 s sem chave. A busca e a lista de resultados somam alguns segundos; registros já no cache carregam na hora. O tempo real depende da sua conexão e da carga do NCBI.';

  @override
  String get apiKeyCurrent => 'Atual';

  @override
  String reviewFound(String n, String query) {
    return '$n genomas encontrados · $query';
  }

  @override
  String get reviewColOrganism => 'Organismo';

  @override
  String get reviewColStatus => 'Situação';

  @override
  String reviewUnverifiedNote(String n) {
    return '$n genomas estão marcados como UNVERIFIED pelo NCBI. Recomendamos não importá-los.';
  }

  @override
  String get reviewDeselectUnverified => 'Desmarcar UNVERIFIED';

  @override
  String get reviewSelectMenu => 'Selecionar';

  @override
  String reviewImportN(String n) {
    return 'Importar $n';
  }

  @override
  String get reviewNoMatch => 'Nenhum genoma corresponde ao filtro.';

  @override
  String reviewSelectedOf(String selected, String total) {
    return '$selected de $total selecionados';
  }

  @override
  String metaTrees(String n) {
    return '$n árvore(s)';
  }

  @override
  String get itemQueued => 'Na fila';

  @override
  String itemCached(String note) {
    return 'Do cache · $note';
  }

  @override
  String get itemSkipped => 'Ignorado';

  @override
  String itemsRunning(String n) {
    return '$n em andamento';
  }

  @override
  String itemsQueued(String n) {
    return '$n na fila';
  }

  @override
  String itemsDone(String n) {
    return '$n concluído(s)';
  }

  @override
  String itemsCached(String n) {
    return '$n do cache (sem mudanças)';
  }

  @override
  String itemsSkipped(String n) {
    return '$n ignorado(s)';
  }

  @override
  String itemsFailed(String n) {
    return '$n com erro';
  }

  @override
  String get trimToolNt =>
      'Escolha uma ferramenta: trimAl ou ClipKIT mantendo códons inteiros nos genes codificantes (recomendado), ou trimAl coluna a coluna (parte os códons).';

  @override
  String get trimToolAa => 'O trimAl corta o alinhamento de proteínas.';

  @override
  String get trimToolCodon =>
      'O alinhamento de códons é aparado mantendo códons inteiros: pelo trimAl (-backtrans) ou pelo ClipKIT (--codon).';

  @override
  String runTrimWith(String tool) {
    return 'Cortar com $tool';
  }

  @override
  String get trimDescAutomated1 =>
      'O trimAl analisa o alinhamento e escolhe sozinho entre -gappyout e -strict. Recomendado para árvores de máxima verossimilhança.';

  @override
  String get trimDescStrict =>
      'Remove colunas com muitos gaps e também colunas com baixa similaridade entre as sequências.';

  @override
  String get trimDescStrictplus =>
      'Como o -strict, mas ajustado para árvores de neighbour-joining; remove mais.';

  @override
  String get trimDescGappyout =>
      'Remove só as colunas com muitos gaps, com um limite calculado pela distribuição de gaps do alinhamento.';

  @override
  String get trimDescNogaps =>
      'Remove toda coluna que tenha pelo menos um gap. Pode remover muito.';

  @override
  String get trimDescNoallgaps =>
      'Remove apenas colunas formadas só por gaps. A opção mais leve.';

  @override
  String get trimDescManual =>
      'Você define os limites de gaps, similaridade e consistência.';

  @override
  String get clipkitDescGappy =>
      'Remove colunas com mais gaps que o limite (90 %). Simples e conservador.';

  @override
  String get clipkitDescKpic =>
      'Mantém apenas colunas informativas para parcimônia e constantes; remove o resto.';

  @override
  String get clipkitDescKpicGappy =>
      'kpic, e também remove colunas com muitos gaps.';

  @override
  String get clipkitDescSmartGap =>
      'Escolhe o limite de gaps a partir do próprio alinhamento (padrão do ClipKIT).';

  @override
  String get clipkitRequired =>
      'O corte codon-aware usa o ClipKIT, que ainda não está instalado. Instale abaixo (precisa do Python 3).';

  @override
  String get frameInfo =>
      'Onde começa o primeiro códon completo deste alinhamento. Base 1 (+1): o alinhamento começa na 1ª posição de um códon, como na maioria dos CDS. +2 ou +3: começa 1 ou 2 bases dentro de um códon. O SPLACE escolhe o quadro com menos códons de parada; ele é usado para traduzir, cortar códons inteiros e separar a 1ª, 2ª e 3ª posições do códon (teste de saturação e partições).';

  @override
  String frameOption(String n) {
    return 'Base $n (+$n)';
  }

  @override
  String frameTag(String n) {
    return 'Códons começam na base $n';
  }

  @override
  String get viewerHint =>
      'Passe o mouse sobre uma base para ver a posição · arraste com o botão direito para rolar · arraste a borda dos nomes para mudar a largura, clique duas vezes para ajustar.';

  @override
  String viewerColumn(String n) {
    return 'Coluna $n';
  }

  @override
  String viewerSeqPos(String n) {
    return 'posição $n na sequência (sem contar gaps)';
  }

  @override
  String get viewerGap => 'gap';

  @override
  String viewerIdentical(String pct) {
    return '$pct% das sequências têm esta base';
  }

  @override
  String viewerCodonPos(String n) {
    return '$nª posição do códon';
  }

  @override
  String get viewerNames => 'Sequências';

  @override
  String siCacheInfo(String cached, String total, String todo) {
    return '$cached de $total buscas já foram feitas antes e serão reaproveitadas; apenas $todo serão enviadas.';
  }

  @override
  String siCacheAll(String total) {
    return 'Todas as $total buscas já foram feitas antes: os resultados aparecem na hora, nada é enviado.';
  }

  @override
  String get siCacheClear => 'Limpar resultados salvos';

  @override
  String get trimFrameWarning =>
      '**O trimAl remove colunas individuais, não códons.** Num gene codificante alinhado como DNA, remover uma ou duas bases de um códon desloca o quadro de leitura de tudo o que vem depois: as 1ª, 2ª e 3ª posições se misturam. Sem \"Manter o quadro de leitura\", esses genes viram uma partição única, sem posições do códon, e ficam fora do teste de saturação.';

  @override
  String get frameLostTag => 'Quadro de leitura perdido no corte';

  @override
  String get frameLostInfo =>
      'O trimAl removeu colunas individuais deste gene codificante, então o tamanho pode não ser múltiplo de 3 e as posições do códon são desconhecidas. Ele é concatenado como uma partição única e fica fora do teste de saturação. Ative \"Manter o quadro de leitura\" em 5e e corte de novo para manter as partições por códon.';

  @override
  String frameShiftTag(String n) {
    return '$n sequência(s) com quadro deslocado';
  }

  @override
  String get frameShiftInfo =>
      'O MAFFT colocou gaps que não são múltiplos de 3 dentro destas sequências: a partir desse gap, as bases ficam na posição do códon errada da coluna. Isso é comum ao alinhar genes codificantes como nucleotídeos. As partições por códon e o teste de saturação ficam aproximados para elas; o modo Codon-aware (5b) evita isso.';

  @override
  String get saturationRunning => 'Teste de saturação';

  @override
  String saturationProgress(String done, String total) {
    return '$done de $total alinhamentos';
  }

  @override
  String get matrixRunning => 'Montando a supermatriz';

  @override
  String get matrixRunningBody => 'Concatenando os alinhamentos';

  @override
  String get iqtreeRunning => 'IQ-TREE: árvore de máxima verossimilhança';

  @override
  String get raxmlRunning => 'RAxML-NG: árvore de máxima verossimilhança';

  @override
  String treeRunSubtitle(String model, String elapsed) {
    return 'Modelo $model · tempo $elapsed';
  }

  @override
  String get geneTreesProgress => 'Árvores por gene';

  @override
  String get localFooter => 'Roda neste computador.';

  @override
  String get treeFooter =>
      'Roda neste computador. Análises longas podem levar minutos ou horas; acompanhe cada etapa aqui.';

  @override
  String get treeDetailPrepare => 'Lendo o alinhamento e as partições';

  @override
  String get treeDetailModelStart => 'Testando modelos de substituição';

  @override
  String treeDetailModels(String done, String total) {
    return 'Modelo $done de até $total testados';
  }

  @override
  String treeDetailStarts(String done, String total) {
    return 'Busca ML a partir da árvore inicial $done de $total';
  }

  @override
  String get treeDetailSearchStart => 'Construindo as árvores iniciais';

  @override
  String treeDetailIteration(String n) {
    return 'Busca da árvore: iteração $n';
  }

  @override
  String treeDetailLeft(String time) {
    return '~$time restantes';
  }

  @override
  String treeDetailBoot(String done, String total) {
    return 'Réplica de bootstrap $done de $total';
  }

  @override
  String get treeDetailSupport => 'Calculando o suporte dos ramos';

  @override
  String get treeDetailFinish => 'Gravando os resultados';

  @override
  String treeDetailBestModel(String model) {
    return 'Melhor modelo: $model';
  }

  @override
  String treeDetailUfboot(String n) {
    return '$n réplicas UFBoot';
  }

  @override
  String treeDetailConverged(String n) {
    return 'Bootstrap convergiu após $n';
  }

  @override
  String get treeStagePrepare => 'Alinhamento';

  @override
  String get treeStageModel => 'Seleção do modelo';

  @override
  String get treeStageSearch => 'Busca da árvore ML';

  @override
  String get treeStageSupport => 'Suporte dos ramos';

  @override
  String get treeStageBootstrap => 'Bootstrap';

  @override
  String get treeStageFinish => 'Resultados';

  @override
  String get step9 => 'Visualizar Árvore';

  @override
  String metaViewable(String n) {
    return '$n árvore(s)';
  }

  @override
  String get viewTree => 'Ver árvore';

  @override
  String get treeViewIntro =>
      'Veja as árvores inferidas. Enraíze pelo grupo externo, pelo ponto médio ou em qualquer ramo, confira o suporte de cada clado, encontre espécies e exporte a figura.';

  @override
  String get treeChoose => 'Árvore';

  @override
  String get treeMissing => 'O arquivo da árvore não foi encontrado.';

  @override
  String get treeKindMl => 'Árvore ML com suporte';

  @override
  String get treeKindConsensus => 'consenso do bootstrap';

  @override
  String get treeKindFbp => 'árvore ML, suporte bootstrap (FBP)';

  @override
  String get treeKindTbe => 'árvore ML, transfer bootstrap (TBE)';

  @override
  String get treeKindBest => 'melhor árvore ML (sem suporte)';

  @override
  String treeKindGene(String gene) {
    return 'árvore do gene $gene';
  }

  @override
  String treeParseError(String error) {
    return 'Não foi possível ler esta árvore: $error';
  }

  @override
  String get treeRootAsIs => 'Como inferida';

  @override
  String get treeRootOutgroup => 'Grupo externo';

  @override
  String get treeRootMidpoint => 'Ponto médio';

  @override
  String get treeRootManual => 'Enraizada no ramo escolhido';

  @override
  String get treePhylogram => 'Filograma';

  @override
  String get treeCladogram => 'Cladograma';

  @override
  String get treeShowSupport => 'Suporte';

  @override
  String get treeItalic => 'Itálico';

  @override
  String get treeSearch => 'Buscar espécie';

  @override
  String get treeFont => 'Tamanho da fonte';

  @override
  String treeStats(String tips, String length) {
    return '$tips terminais · comprimento total $length';
  }

  @override
  String get treeHint =>
      'Clique num ramo ou nome para selecioná-lo: enraizar ali, girar, colapsar ou encurtar. Ctrl+Z desfaz.';

  @override
  String treeClade(String n) {
    return 'Clado com $n terminais';
  }

  @override
  String treeSupportValue(String value) {
    return 'suporte $value';
  }

  @override
  String get treeRootHere => 'Enraizar aqui';

  @override
  String get treeCollapse => 'Recolher';

  @override
  String get treeExpand => 'Expandir';

  @override
  String get treeClearSelection => 'Limpar';

  @override
  String get treeSupportLegend => 'Suporte dos ramos:';

  @override
  String get treeSupportInfoUfboot =>
      'Os rótulos do IQ-TREE são SH-aLRT/UFBoot (ou só UFBoot); as cores usam o último valor. UFBoot ≥ 95 é suporte forte (Hoang et al. 2018); SH-aLRT ≥ 80 junto com ele deixa o clado bem suportado.';

  @override
  String get treeSupportInfoBootstrap =>
      'Bootstrap padrão (não paramétrico): ≥ 70 costuma ser considerado bom suporte (Hillis & Bull 1993).';

  @override
  String get siSelfTag => 'este registro';

  @override
  String siSelfInfo(String accessions) {
    return 'O BLAST encontrou a própria sequência deste genoma no GenBank ($accessions). O dado está lá, mas um registro não confirma o próprio nome, então não conta como evidência. Outros registros da espécie contariam.';
  }

  @override
  String siVerdictSelfOnly(String accessions, String current) {
    return 'O GenBank tem a sequência deste genoma ($accessions), mas nenhum outro registro de $current, então o GenBank não consegue confirmar o nome.';
  }

  @override
  String siVerdictNoReference(String current, String best, String identity) {
    return 'Nenhuma outra sequência de $current foi encontrada. A mais próxima é $best ($identity%), abaixo de 98% — distância comum entre espécies —, então o nome atual não é contestado.';
  }

  @override
  String get siNoReference => 'Sem referência';

  @override
  String get apiKeyEstimateColSize => 'Download';

  @override
  String get apiKeyEstimateColGain => 'Ganho';

  @override
  String apiKeyEstimateSaved(String n, String saved) {
    return 'Para $n genomas, a chave economiza cerca de $saved.';
  }

  @override
  String missTypeNotSelected(String type) {
    return 'anotado como $type';
  }

  @override
  String missTypeNotSelectedInfo(String gene, String type) {
    return '$gene está anotado neste registro como $type, um tipo de feature não selecionado na etapa 3. Selecione $type lá para incluí-lo.';
  }

  @override
  String get missGeneOnly => 'só feature \'gene\'';

  @override
  String missGeneOnlyInfo(String gene, String location) {
    return 'O registro tem uma feature \'gene\' para $gene ($location), mas não a feature CDS/rRNA, que é a que o SPLACE lê. O gene está na sequência; a anotação é que está incompleta. Use um registro mais bem anotado ou anote o genoma (ex.: MitoS ou MitoFinder).';
  }

  @override
  String missNoAnnotation(String types) {
    return 'sem $types anotados';
  }

  @override
  String missNoAnnotationInfo(String types) {
    return 'Este registro não tem nenhum $types anotado, então nenhum gene pode ser extraído dele. Provavelmente é uma submissão sem anotação; anote o genoma (ex.: MitoS ou MitoFinder) ou remova-o.';
  }

  @override
  String get missUnknownName => 'nome não reconhecido';

  @override
  String missUnknownNameInfo(String names, String gene) {
    return 'O registro tem genes com nomes que o SPLACE não reconheceu ($names). Um deles pode ser $gene: corrija o nome na etapa 3 (nomes desconhecidos).';
  }

  @override
  String get missPartial => 'genoma parcial';

  @override
  String missPartialInfo(String length, String typical, String gene) {
    return 'O registro tem $length pb, bem menos que o genoma típico deste conjunto ($typical pb): é um genoma parcial e não cobre $gene.';
  }

  @override
  String get missNotAnnotated => 'não anotado';

  @override
  String missNotAnnotatedInfo(String gene) {
    return '$gene não está entre os genes anotados deste registro: está ausente da sequência ou os autores não o anotaram.';
  }

  @override
  String get missNotInFiles => 'ausente nos arquivos';

  @override
  String missNotInFilesInfo(String gene) {
    return 'Nenhum arquivo FASTA traz sequência de $gene para esta amostra.';
  }

  @override
  String missPresentIn(String present, String total, String missing) {
    return 'em $present de $total amostras · falta em $missing';
  }

  @override
  String missSampleCount(String missing, String total) {
    return 'faltam $missing de $total marcadores';
  }

  @override
  String missSummary(
    String markers,
    String selected,
    String samples,
    String total,
  ) {
    return '$markers dos $selected marcadores selecionados faltam em pelo menos uma amostra; $samples de $total amostras não têm pelo menos um marcador.';
  }

  @override
  String get missMatrixNote =>
      'Na supermatriz (etapa 7), o marcador ausente vira \'?\' para a amostra quando “Permitir sequências ausentes” está ligado; senão, os marcadores que faltam em algum táxon ficam de fora. Algumas lacunas são normais; uma amostra sem a maioria dos marcadores contribui pouco e pode valer a pena removê-la.';

  @override
  String missFixableNote(String n) {
    return '$n lacuna(s) podem ser corrigidas no SPLACE (etiquetas laranja). Passe o mouse na etiqueta para ver como.';
  }

  @override
  String get missByMarker => 'Por marcador';

  @override
  String get missBySample => 'Por amostra';

  @override
  String get missHoverHint =>
      'Cada etiqueta diz por que o marcador falta; passe o mouse para detalhes.';

  @override
  String dupSummary(String genes, String genomes) {
    return '$genes gene(s) têm mais de uma cópia anotada em $genomes genoma(s).';
  }

  @override
  String get dupFastaHelp =>
      'Estas amostras têm mais de uma sequência do mesmo locus nos arquivos FASTA. A coluna “Usada” mostra o que entra no alinhamento, conforme a opção “Duplicatas e nomes”.';

  @override
  String get dupIdentical => 'cópias idênticas';

  @override
  String get dupIdenticalInfo =>
      'Mesma sequência (ex.: cópias da repetição invertida em plastomas): qualquer cópia dá o mesmo resultado.';

  @override
  String get dupContained => 'fragmento de outra cópia';

  @override
  String get dupContainedInfo =>
      'A cópia menor é parte da maior: em geral uma anotação duplicada ou truncada. Use a cópia completa.';

  @override
  String get dupDifferent => 'cópias diferentes';

  @override
  String get dupDifferentInfo =>
      'As cópias diferem na sequência (parálogos, cópias tipo NUMT ou erros de anotação). A escolha muda o alinhamento: prefira a cópia mais próxima dos outros genomas.';

  @override
  String get siCardTitle => 'Identificação das espécies';

  @override
  String get siCardBody =>
      'Compara um gene de barcode (COI, 16S…) de cada genoma com o GenBank (BLAST) e o BOLD para encontrar registros com nome errado antes de obter a taxonomia, que depende dos nomes.';

  @override
  String get siCardOptional => 'Opcional';

  @override
  String get siCardRerun => 'Identificar novamente';

  @override
  String siSumChecked(String n) {
    return '$n genomas verificados';
  }

  @override
  String siSumConfirmed(String n) {
    return '$n nomes confirmados';
  }

  @override
  String siSumFlagged(String n) {
    return '$n possíveis erros de identificação';
  }

  @override
  String siSumNoReference(String n) {
    return '$n sem referência';
  }

  @override
  String get headerDefault => 'Cabeçalho padrão · usado ao alinhar';

  @override
  String viewerMafftColumn(String n) {
    return 'coluna $n do alinhamento MAFFT';
  }

  @override
  String get satSourceMafft =>
      'Testado no alinhamento do MAFFT (sem trimagem).';

  @override
  String get satSourceClipkit =>
      'Testado no alinhamento trimado (ClipKIT nas proteínas, retrotraduzido): só códons inteiros foram removidos, então cada base mantém sua posição no códon.';

  @override
  String satSourceTrimmed(String kept, String total) {
    return 'Testado no alinhamento trimado: o trimAl manteve $kept dos $total códons do MAFFT, e o SPLACE removeu só códons inteiros. Assim, cada coluna mantém a posição no códon que tinha no alinhamento do MAFFT (passe o mouse numa coluna trimada no passo 5 para ver de qual coluna do MAFFT ela veio).';
  }

  @override
  String satSourceFrame(String frame) {
    return 'Os códons agora começam na base 1 porque o códon incompleto antes da base $frame do alinhamento do MAFFT foi removido; as posições são as mesmas.';
  }

  @override
  String get preserveCaseHelp => 'Mantém maiúsculas e minúsculas';

  @override
  String get configureTrimming => 'Configurar trimagem';

  @override
  String panelTrimParams(String tool) {
    return 'Parâmetros do $tool';
  }

  @override
  String trimmedCount(String n, String total) {
    return '$n de $total marcadores trimados';
  }

  @override
  String get showRemovedColumns => 'Mostrar regiões removidas';

  @override
  String get nextStop => 'Próximo stop';

  @override
  String get translationTable => 'Tabela de tradução';

  @override
  String translationTableTitle(String code, String name) {
    return 'Código genético $code · $name';
  }

  @override
  String get translationTableBody =>
      'Códons do código genético escolhido acima para a visão em aminoácidos (ele nunca é deduzido dos dados). Os stops estão em vermelho; o selo conta quantos stops internos deste alinhamento usam cada um.';

  @override
  String codonSecondBase(String base) {
    return '2ª $base';
  }

  @override
  String codonDiffers(String codon, String aa, String standard) {
    return '$codon → $aa neste código ($standard no código padrão)';
  }

  @override
  String get codonDiffersLegend => 'Difere do código padrão';

  @override
  String stopCodonsFound(String list) {
    return 'Stops internos neste alinhamento: $list';
  }

  @override
  String removedColumns(String n, String total) {
    return '$n de $total colunas removidas';
  }

  @override
  String get removedColumnsInfo =>
      'O alinhamento MAFFT inteiro é mostrado; as colunas em cinza foram removidas pela trimagem. Desligue \"Mostrar regiões removidas\" para ver só as colunas mantidas.';

  @override
  String get removedLegend => 'Removido';

  @override
  String get viewerRemoved => 'removida pela trimagem';

  @override
  String viewerCodon(String codon) {
    return 'códon $codon';
  }

  @override
  String get plotObservedVsCorrectedBody =>
      'Cada ponto é um par de sequências. Pontos sobre a linha tracejada: todas as substituições são visíveis; pontos que se achatam abaixo dela: substituições múltiplas escondidas (saturação).';

  @override
  String get plotTsTvBody =>
      'As transições (vermelho) costumam saturar primeiro: quando param de crescer enquanto as transversões (azul) continuam subindo, o conjunto está saturado.';

  @override
  String plotPositions(String model) {
    return 'p de cada posição do códon × distância $model';
  }

  @override
  String get plotPositionsBody =>
      'Mesmos pares, uma cor por posição do códon, contra a distância do gene inteiro. Uma posição cujos pontos se achatam enquanto as outras continuam subindo está saturada (geralmente a 3ª).';

  @override
  String get allPositions => 'todas as posições';

  @override
  String get compositionTitle => 'Composição e variabilidade';

  @override
  String get compositionInfo =>
      'Teste qui-quadrado de homogeneidade das frequências de bases entre os táxons (teste de contingência, como no PAUP*) e o número de táxons cuja composição difere da do conjunto (gl = 3, p < 0,05, como o IQ-TREE informa). Posições saturadas costumam variar na composição (ex.: GC3), o que também pode enviesar a árvore.';

  @override
  String get metricVariable => 'Variáveis';

  @override
  String get metricInformative => 'Informativos (parcimônia)';

  @override
  String get metricTaxaFailing => 'Táxons reprovados';

  @override
  String get uncorrectablePairs => 'pares sem correção possível';

  @override
  String get heatmapIntensity =>
      'O tom dentro de cada cor segue a inclinação e os pares sem correção.';

  @override
  String get treeRootLabel => 'Raiz';

  @override
  String get treeShapeLabel => 'Forma';

  @override
  String get treeSupportValues => 'Valores';

  @override
  String get treeSupportDots => 'Pontos';

  @override
  String get treeSupportHidden => 'Ocultar';

  @override
  String get treeDisplay => 'Exibir';

  @override
  String get treeAlignLabels => 'Alinhados à direita';

  @override
  String get treeExport => 'Exportar';

  @override
  String treeMatches(String n) {
    return '$n encontrado(s)';
  }

  @override
  String treeHoverTip(String name, String length) {
    return '$name · comprimento do ramo $length';
  }

  @override
  String treeHoverClade(String n, String support, String length) {
    return 'Clado com $n terminais · suporte $support · comprimento do ramo $length — clique para selecionar';
  }

  @override
  String get treeOutgroupLegend => 'Grupo externo';

  @override
  String get ord1 => '1ª';

  @override
  String get ord2 => '2ª';

  @override
  String get ord3 => '3ª';

  @override
  String get listAnd => 'e';

  @override
  String codonPositionN(String n) {
    return '$n posição';
  }

  @override
  String reasonInsufficient(String sites, String pairs) {
    return 'Poucos dados para avaliar ($sites sítios, $pairs pares de sequências).';
  }

  @override
  String reasonUndefined(String pct) {
    return '$pct% dos pares de sequências são tão diferentes que a distância não pode ser corrigida.';
  }

  @override
  String reasonLowDivergence(String p) {
    return 'As sequências diferem muito pouco (p médio = $p), então a saturação não é um problema.';
  }

  @override
  String reasonSlopeOk(String slope) {
    return 'Saturação do PhyKIT $slope: as diferenças observadas ainda crescem junto com a distância na árvore.';
  }

  @override
  String reasonSlopeModerate(String slope, String limit) {
    return 'Saturação do PhyKIT $slope (abaixo de $limit): as diferenças observadas crescem mais devagar que a distância na árvore.';
  }

  @override
  String reasonSlopeSaturated(String slope, String limit) {
    return 'Saturação do PhyKIT $slope (abaixo de $limit): as diferenças observadas quase pararam de crescer com a distância na árvore; muitas mudanças foram apagadas por outras posteriores.';
  }

  @override
  String reasonIss(String iss) {
    return 'Iss $iss (0,8 ou mais): os sítios parecem quase aleatórios.';
  }

  @override
  String get recTitle => 'O que fazer com as posições do códon';

  @override
  String recBasedOn(String taxa, String length) {
    return 'Com base na matriz concatenada ($taxa táxons, $length pb). Cada posição é testada separadamente.';
  }

  @override
  String get recWhatToDo => 'Sugestão';

  @override
  String recExcludeOne(String positions) {
    return 'Remova a $positions posição do códon da matriz (etapa 7b) e mantenha as outras duas, particionadas por posição do códon.';
  }

  @override
  String recExcludeMany(String positions) {
    return 'Remova as posições $positions do códon da matriz (etapa 7b) e mantenha a restante.';
  }

  @override
  String get recWhy =>
      'Por quê: uma posição saturada mudou tantas vezes que suas diferenças já não mostram há quanto tempo as linhagens se separaram. Ela acrescenta principalmente ruído (homoplasia) e pode agrupar ramos longos.';

  @override
  String recCautionOne(String positions) {
    return 'A $positions posição tem saturação moderada: mantenha-a numa partição própria e veja se a árvore muda sem ela.';
  }

  @override
  String recCautionMany(String positions) {
    return 'As posições $positions têm saturação moderada: mantenha-as em partições próprias e veja se a árvore muda sem elas.';
  }

  @override
  String applyExclude(String positions) {
    return 'Marcar $positions na etapa 7b';
  }

  @override
  String get applyKeepAll =>
      'Deixar todos os \"Excluir posição\" desmarcados na etapa 7b';

  @override
  String applyPerLocusInfo(String loci) {
    return 'Cada gene exclui só as posições saturadas nele. Genes com resultado diferente da sugestão geral: $loci.';
  }

  @override
  String recLociDiffer(String n) {
    return '$n gene(s) têm resultado diferente da sugestão geral; veja as linhas deles no heatmap.';
  }

  @override
  String get recStatusKeep => 'A etapa 7b agora mantém todas as posições.';

  @override
  String recStatusExclude(String positions) {
    return 'A etapa 7b agora exclui: $positions posição(ões).';
  }

  @override
  String recStatusPerLocus(String n) {
    return 'A etapa 7b agora exclui posições por gene ($n genes).';
  }

  @override
  String get saturationHeatmapHelp =>
      'Cada célula mostra o valor de saturação do PhyKIT (inclinação de p sobre a distância na árvore) — 1 significa sem saturação, perto de 0 saturada. Azul é OK, vermelho é saturada.';

  @override
  String get rankSpecies => 'Espécie';

  @override
  String get rankGenus => 'Gênero';

  @override
  String get rankFamily => 'Família';

  @override
  String get rankOrder => 'Ordem';

  @override
  String get rankClass => 'Classe';

  @override
  String get rankPhylum => 'Filo';

  @override
  String get rankKingdom => 'Reino';

  @override
  String get tipFieldOriginal => 'Nome como na árvore';

  @override
  String get tipFieldAccession => 'Acesso';

  @override
  String get tipFieldVoucher => 'Voucher';

  @override
  String get treeBold => 'Negrito';

  @override
  String get treeParens => 'Entre parênteses';

  @override
  String get treeRemovePart => 'Remover esta parte';

  @override
  String get treeAddPart => 'Adicionar campo';

  @override
  String get treeMovePart => 'Mover para a esquerda';

  @override
  String get treeLabelFormat => 'Formato dos rótulos';

  @override
  String get treeLabelParts => 'Partes do rótulo';

  @override
  String get treeLabelPartsInfo =>
      'Cada rótulo é montado com estes campos, nesta ordem. I = itálico, B = negrito, ( ) = entre parênteses. Os valores vêm dos registros do GenBank e da taxonomia obtida na etapa 2; um campo que a ponta não tem é omitido.';

  @override
  String get treeTabTree => 'Árvore';

  @override
  String get treeTabBranches => 'Ramos e suporte';

  @override
  String get treeTabLabels => 'Rótulos';

  @override
  String get treeTabColours => 'Cores';

  @override
  String get treeOrderLabel => 'Ordem e rotação';

  @override
  String get treeLadderUp => 'Escada ↑';

  @override
  String get treeLadderDown => 'Escada ↓';

  @override
  String get treeFlip => 'Inverter árvore';

  @override
  String get treeRotateClade => 'Girar clado';

  @override
  String get treeRotateInfo =>
      'Girar um clado troca a ordem dos seus ramos; as relações não mudam. Selecione um clado clicando no seu ramo e depois gire. Inverter árvore espelha a árvore toda.';

  @override
  String get treeLineWidth => 'Espessura das linhas';

  @override
  String get treePointSize => 'Tamanho dos pontos';

  @override
  String get treePoints => 'Pontos';

  @override
  String get treeTipDots => 'Pontos nas pontas';

  @override
  String get treeSupportMin => 'Mostrar suporte a partir de';

  @override
  String get treeSupportMinInfo =>
      'Só os nós com suporte igual ou acima deste valor recebem número ou ponto (por exemplo, pontos só nos nós com UFBoot ≥ 95).';

  @override
  String get treeSupportAll => 'Todos';

  @override
  String get treeLongBranches => 'Ramos longos';

  @override
  String get treeLongBranchesInfo =>
      'Um ramo encurtado é desenhado mais curto do que é e marcado com //, para sobrar espaço para o resto da árvore. Contam como longos os ramos com mais de 5× o comprimento mediano e mais de 30% da altura da árvore. Também dá para selecionar qualquer ramo e encurtá-lo. O comprimento real continua na barra de informação e no Newick exportado.';

  @override
  String treeShortenLong(String n) {
    return 'Encurtar ramos longos ($n)';
  }

  @override
  String get treeRestoreLengths => 'Restaurar comprimentos';

  @override
  String get treeShortenBranch => 'Encurtar ramo //';

  @override
  String get treeRestoreBranch => 'Restaurar comprimento';

  @override
  String get treeShortenedLegend => 'ramo encurtado (fora de escala)';

  @override
  String get treeColorRows => 'Colorir linhas por';

  @override
  String get treeColorRowsInfo =>
      'Cada táxon do nível escolhido ganha uma cor de fundo na linha inteira, com o nome à direita. Os níveis vêm da taxonomia obtida na etapa 2 (ou da espécie e do gênero nos nomes das pontas).';

  @override
  String get treeColorNone => 'Nenhum';

  @override
  String get treeNoTaxonomy =>
      'Sem taxonomia para estas pontas. Obtenha-a na etapa 2.';

  @override
  String get treeExportInfo =>
      'SVG e PNG são a figura como você a vê (rótulos, cores, legenda, ramos encurtados), sem as marcas de seleção e de busca. O Newick tem a topologia e os comprimentos reais.';

  @override
  String get treeUndo => 'Desfazer';

  @override
  String get treeRedo => 'Refazer';

  @override
  String get treeReset => 'Recomeçar';

  @override
  String get zipTitle => 'Salvar tudo';

  @override
  String get zipBody =>
      'Um ZIP com as sequências, alinhamentos, tabelas de saturação, matrizes, partições, árvores e logs deste projeto.';

  @override
  String get concatHowTitle => 'Como a supermatriz é montada';

  @override
  String get concatHow1 =>
      'O alinhamento de cada locus é colocado lado a lado, uma linha por táxon (casada pelo cabeçalho FASTA).';

  @override
  String get concatHow2 =>
      'Um táxon que não tem um locus recebe ? nesse trecho (dado faltante), e todas as linhas ficam com o mesmo comprimento.';

  @override
  String get concatHow3 =>
      'O arquivo de partições diz ao IQ-TREE / RAxML-NG onde cada locus (e cada posição do códon) começa e termina, para que cada parte tenha seu próprio modelo de substituição.';

  @override
  String get concatPictureCaption => '? = táxon sem aquele locus';

  @override
  String get concatBlockSource => 'Quais alinhamentos juntar';

  @override
  String get concatBlockSourceBody =>
      'Os alinhamentos trimados da etapa 5 (sem as colunas mal alinhadas) ou a saída do MAFFT como está.';

  @override
  String get concatBlockPartitions => 'Como dividir a matriz em partições';

  @override
  String get concatBlockPartitionsBody =>
      'Cada partição recebe seu próprio modelo de substituição na busca da árvore. Aqui também dá para remover posições do códon.';

  @override
  String get schemeByLocusCodonBody =>
      'Cada locus é uma partição e, nos loci codificantes, cada posição do códon (1ª, 2ª, 3ª) tem a sua. As posições evoluem em taxas muito diferentes; com MFP+MERGE o IQ-TREE junta as partições que se comportam igual.';

  @override
  String get schemeByLocusBody =>
      'Uma partição por locus: cada gene tem seu modelo, e as três posições do códon o compartilham. Bom para rRNAs, loci não codificantes ou com poucos sítios.';

  @override
  String get schemeByCodonBody =>
      'Três partições compartilhadas por todos os loci codificantes (todas as 1ªs posições juntas, todas as 2ªs, todas as 3ªs); cada locus não codificante continua sendo uma partição. Menos modelos que por locus e posição.';

  @override
  String get schemeNoneBody =>
      'Um modelo para a matriz inteira (sem arquivo de partições). Mais simples, mas costuma se ajustar pior quando genes ou posições evoluem de formas diferentes.';

  @override
  String get codonExclusionHelp =>
      'As posições removidas saem só dos loci codificantes. A etapa 6 (saturação) pode preencher isto.';

  @override
  String get concatSummaryTitle => 'O que vai ser gerado';

  @override
  String get concatSummaryLoci => 'Loci';

  @override
  String get concatSummarySource => 'Alinhamentos';

  @override
  String get concatNoOutgroup => 'nenhum (ponto médio)';

  @override
  String get concatNothingExcluded => 'nada';

  @override
  String get treeXScale => 'Escala horizontal';

  @override
  String get treeXScaleInfo =>
      'Estica (acima de 1) ou encolhe (abaixo de 1) os ramos no eixo x. Acima de 1 a figura fica mais larga que a largura definida.';

  @override
  String get treeYScale => 'Escala vertical';

  @override
  String get treeYScaleInfo =>
      'Espaço entre as linhas no eixo y (1 = 1,7 vez o tamanho do rótulo).';

  @override
  String get treeBranchLengths => 'Comprimento dos ramos';

  @override
  String get treeShowLengths => 'Mostrar como rótulo';

  @override
  String get treeSupportFont => 'Tamanho do texto do suporte';

  @override
  String get treeProportional => 'Tamanho pelo valor do suporte';

  @override
  String get treeShowLabels => 'Mostrar rótulos das pontas';

  @override
  String get treeTabFigure => 'Figura e exportação';

  @override
  String get treeScaleBar => 'Barra de escala';

  @override
  String get treeShowScaleBar => 'Mostrar';

  @override
  String get treeScaleBarLength => 'Comprimento da barra de escala';

  @override
  String get treeScaleBarLengthInfo =>
      'Em substituições por sítio. 0 escolhe um valor redondo de cerca de um quinto da árvore.';

  @override
  String get treeFigureLegend => 'Legenda';

  @override
  String get treeLegendInFigure => 'Desenhar na figura';

  @override
  String get treeFigureLegendInfo =>
      'As cores do suporte, o grupo externo, os ramos encurtados e as cores das linhas são desenhados abaixo da árvore e vão junto no PNG e no SVG.';

  @override
  String get treeAtTips => 'Nas pontas';

  @override
  String get treeLegendMarks => 'Marcas';

  @override
  String get treeBranchColor => 'Cor dos ramos';

  @override
  String get treeSupportPosition => 'Posição do suporte';

  @override
  String get treeAboveBranch => 'Acima do ramo';

  @override
  String get treeBelowBranch => 'Abaixo do ramo';

  @override
  String get treeSupportColor => 'Cor';

  @override
  String get treeFontFamily => 'Fonte';

  @override
  String get treeLabelColor => 'Cor dos rótulos';

  @override
  String get treeLabelPosition => 'Posição dos rótulos';

  @override
  String get treeLabelGap => 'Distância da ponta';

  @override
  String get treeLabelGapInfo =>
      'Espaço entre o fim do ramo (ou a coluna alinhada) e o rótulo.';

  @override
  String get treeColorApply => 'Aplicar as cores em';

  @override
  String get treeColorApplyRows => 'Fundo da linha';

  @override
  String get treeColorApplyBranches => 'Ramos';

  @override
  String get treeColorApplyLabels => 'Rótulos';

  @override
  String get treeCladeColors => 'Cores de clados';

  @override
  String get treeCladeColorsInfo =>
      'Selecione um clado clicando no seu ramo e dê uma cor a ele: os ramos do clado passam a usá-la. O desfazer também vale aqui.';

  @override
  String get treeColorClade => 'Colorir clado';

  @override
  String get treeRemoveColor => 'Remover cor';

  @override
  String treeClearCladeColors(String n) {
    return 'Remover todas ($n)';
  }

  @override
  String get treeBackground => 'Fundo da figura';

  @override
  String get treeLegendPosition => 'Posição';

  @override
  String get treeLegendBelow => 'Abaixo da árvore';

  @override
  String get treeLegendRight => 'À direita da árvore';

  @override
  String get treeLegendStyle => 'Estilo';

  @override
  String get treeLegendBox => 'Em uma caixa';

  @override
  String get treeLegendTitle => 'Título da legenda';

  @override
  String get treeLegendTitleHint => 'ex.: Legenda';

  @override
  String get treeLegendFont => 'Tamanho do texto da legenda';

  @override
  String get treeFigureTitle => 'Título da figura';

  @override
  String get treeFigureTitleHint => 'Opcional, desenhado acima da árvore';

  @override
  String get treeDpi => 'Resolução do PNG';

  @override
  String get treeDpiInfo =>
      'Pontos por polegada do PNG (300 para a maioria das revistas, 600 para desenhos de linha). O DPI fica gravado no arquivo. Figuras muito grandes são limitadas a 16 000 px no lado maior.';

  @override
  String get treePngBackground => 'Fundo do PNG';

  @override
  String get treeTransparent => 'Transparente';

  @override
  String get treePickColor => 'Escolha uma cor';

  @override
  String get concatBlockOutgroup =>
      'Quais táxons enraízam a árvore (grupo externo)';

  @override
  String get concatBlockOutgroupBody =>
      'Escolha um ou mais táxons fora do grupo estudado; as árvores das etapas 8 e 9 são enraizadas neles. Sem grupo externo, o visualizador enraíza no ponto médio.';

  @override
  String get schemeOneModel => 'um modelo para tudo';

  @override
  String schemePartitionsCount(String n) {
    return '$n partição(ões)';
  }

  @override
  String get resCpu => 'CPU';

  @override
  String get resCpuInfo =>
      'Uso de CPU da máquina, incluindo MAFFT, trimAl, IQ-TREE e as outras ferramentas que o SPLACE executa.';

  @override
  String get resMemory => 'RAM';

  @override
  String resMemoryInfo(String app, String used, String total) {
    return 'O SPLACE usa $app. Máquina: $used de $total em uso (a barra).';
  }

  @override
  String resMemoryInfoApp(String app) {
    return 'O SPLACE usa $app.';
  }

  @override
  String get resDisk => 'Disco';

  @override
  String resFree(String v) {
    return '$v livres';
  }

  @override
  String resDiskInfo(String free, String total) {
    return '$free livres de $total no disco da pasta do projeto (a barra é o espaço em uso).';
  }

  @override
  String get resDiskInfoNone => 'Espaço no disco da pasta do projeto.';

  @override
  String get resNetwork => 'Rede';

  @override
  String resNetworkInfo(String total) {
    return 'Dados que o SPLACE troca com NCBI, GBIF e BLAST, por segundo. Desde que o app abriu: $total.';
  }

  @override
  String get methodsTitle => 'Texto de métodos';

  @override
  String get methodsCopy => 'Copiar texto e referências';

  @override
  String get methodsHelp =>
      'Escrito a partir do que foi feito neste projeto (dados, alinhamento, trimagem, saturação, supermatriz e árvores). Revise e adapte ao seu manuscrito.';

  @override
  String get citeReferences => 'Referências citadas no texto';

  @override
  String get citeTools => 'Ferramentas usadas pelo SPLACE';

  @override
  String get treeBranchStyle => 'Estilo dos ramos';

  @override
  String get treeBranchSquare => 'Ângulos retos';

  @override
  String get treeBranchRounded => 'Cantos arredondados';

  @override
  String get treeBranchSlanted => 'Diagonais';

  @override
  String get treeZoom => 'Zoom (só na tela)';

  @override
  String get treeZoomFit => 'Ajustar à janela';

  @override
  String get treeZoomInfo =>
      'Aumenta ou reduz a visualização da figura nesta tela. A figura exportada mantém a largura definida em Figura e exportação, qualquer que seja o monitor.';

  @override
  String get treeTipDotsLabel => 'Pontas';

  @override
  String get treeLabelCustom => 'Personalizado';

  @override
  String get treeMarkOutgroup => 'Caixa cinza atrás';

  @override
  String get treeRenamed => 'Rótulos editados';

  @override
  String get treeRenamedInfo =>
      'Clique numa ponta e escolha Editar rótulo para escrever um texto próprio (com itálico e negrito). Apague o texto para voltar ao rótulo montado.';

  @override
  String get treeRenameTip => 'Editar rótulo';

  @override
  String treeRestoreNames(String n) {
    return 'Restaurar rótulos editados ($n)';
  }

  @override
  String get treeBarsByRank => 'Barras com nomes de táxons';

  @override
  String get treeBarsByRankInfo =>
      'Desenha uma barra ao lado de cada grupo de pontas que compartilham o táxon deste nível (por exemplo, uma barra por família), com o nome dele. Táxons separados na árvore ganham uma barra por parte.';

  @override
  String get treeHighlightOpacity => 'Opacidade das caixas';

  @override
  String get treeHighlightClade => 'Destacar clado';

  @override
  String get treeRemoveHighlight => 'Remover destaque';

  @override
  String get treeCladeLabel => 'Nomear clado';

  @override
  String treeClearHighlights(String n) {
    return 'Remover destaques ($n)';
  }

  @override
  String treeClearCladeLabels(String n) {
    return 'Remover nomes de clados ($n)';
  }

  @override
  String get treeRemoveFromFigure => 'Remover da árvore';

  @override
  String get treeFigureWidth => 'Largura da figura';

  @override
  String get treeFigureWidthInfo =>
      'Largura real da figura exportada. Os tamanhos de fonte são em pontos, então um rótulo de 9 pt sai impresso com 9 pt nessa largura. A altura vem do número de pontas e da escala vertical.';

  @override
  String get treeWidthPreset => 'Largura de revista';

  @override
  String get treeWidthOneCol => '1 coluna (8,5 cm)';

  @override
  String get treeWidthOneHalf => '1,5 coluna (12 cm)';

  @override
  String get treeWidthTwoCol => 'Página inteira (17,5 cm)';

  @override
  String get treeWidthCustom => 'Personalizada';

  @override
  String get treeMargin => 'Margem';

  @override
  String get treeTitleFont => 'Tamanho do título';

  @override
  String treeFigureSize(String w, String h, String pw, String ph, String dpi) {
    return 'Figura $w × $h cm · PNG $pw × $ph px a $dpi dpi';
  }

  @override
  String get treeTextHint => 'Texto';

  @override
  String get treeTextEmptyHint => 'Deixe vazio para remover.';

  @override
  String get concatHowOpen =>
      'Clique para ver como os loci são colocados lado a lado';

  @override
  String get treeOpenEditor => 'Abrir o editor da árvore';

  @override
  String get treeEditorTitle => 'Editor da árvore';

  @override
  String get treeEditorHint =>
      'Prévia da figura no tamanho real. Clique nela, ou no botão acima, para editá-la na janela inteira.';

  @override
  String get treeCompare => 'Comparar árvores';

  @override
  String get treeModeEdit => 'Editar figura';

  @override
  String get treeDone => 'Concluir';

  @override
  String get treeExportMenu => 'Exportar';

  @override
  String get treeZoomIn => 'Aproximar (Ctrl + roda)';

  @override
  String get treeZoomOut => 'Afastar (Ctrl + roda)';

  @override
  String get treeFitWidth => 'Ajustar à largura';

  @override
  String get treeFitPage => 'Figura inteira';

  @override
  String get treeCanvasHint =>
      'Arraste para mover · a roda rola, Ctrl + roda dá zoom · clique seleciona, Shift + clique seleciona o clado que une duas pontas · botão direito para opções · Del remove · Esc limpa';

  @override
  String get treeTabCompare => 'Comparar';

  @override
  String get treeCompareWith => 'Comparar com';

  @override
  String get treeUntangle => 'Desembaraçar';

  @override
  String get treeUntangleText => 'Girar a outra árvore';

  @override
  String get treeUntangleInfo =>
      'Gira os clados da outra árvore (a topologia não muda) para cada ponta ficar de frente para o seu par na árvore em edição, com o mínimo de linhas cruzadas.';

  @override
  String get treeConflicts => 'Clados diferentes';

  @override
  String get treeConflictsText => 'Marcar';

  @override
  String get treeConflictsInfo =>
      'Os ramos que levam a um clado que a outra árvore não tem são desenhados mais grossos, nesta cor.';

  @override
  String get treeConflictColor => 'Cor dos clados diferentes';

  @override
  String get treeConnectors => 'Linhas entre as árvores';

  @override
  String get treeConnectorsGroup => 'Cores das linhas';

  @override
  String get treeConnectorsInfo =>
      'Com as linhas coloridas por uma categoria (aba Cores), cada ligação fica com a cor do seu táxon; senão, fica cinza.';

  @override
  String get treeCompareGap => 'Espaço entre as árvores';

  @override
  String get treeCompareResult => 'Resultado';

  @override
  String treeCompareCommon(String n) {
    return '$n pontas nas duas árvores';
  }

  @override
  String treeCompareOnly(String a, String b) {
    return '$a só na esquerda · $b só na direita';
  }

  @override
  String treeCompareRf(String rf, String nrf) {
    return 'Distância de Robinson–Foulds $rf (normalizada $nrf)';
  }

  @override
  String treeCompareShared(String n, String a, String b) {
    return '$n clados em comum · $a só na esquerda · $b só na direita';
  }

  @override
  String get treeCompareRfInfo =>
      'Conta os clados que estão numa árvore e não na outra (sem raiz, sobre as pontas que as duas têm). 0 indica a mesma topologia; o valor normalizado divide pela maior distância possível.';

  @override
  String get treeCompareLegend => 'Clado ausente na outra árvore';

  @override
  String get treeCompareTooFew =>
      'As árvores têm menos de 4 pontas em comum: os clados não podem ser comparados.';

  @override
  String treeCompareError(String e) {
    return 'Não foi possível ler a outra árvore: $e';
  }

  @override
  String get treeCompareIntro =>
      'A árvore em edição e a escolhida aqui são desenhadas frente a frente, uma espelhada, com uma linha ligando cada ponta ao seu par. Os estilos (fontes, ramos, suporte, cores) vêm das outras abas.';

  @override
  String get treeSelectParent => 'Clado acima';

  @override
  String get treeLeft => 'Esquerda';

  @override
  String treeDetailMoose(String n) {
    return 'O MOOSE está escolhendo os modelos ($n partição(ões))';
  }

  @override
  String get treeKindRbs => 'árvore ML, bootstrap rápido (RBS)';

  @override
  String get treeKindSh => 'árvore ML, SH-aLRT';

  @override
  String get treeKindEbg => 'árvore ML, bootstrap previsto (EBG)';

  @override
  String raxmlOld(String version) {
    return 'Foi encontrado o RAxML-NG $version. O MOOSE, a busca rápida e os suportes SH-aLRT, EBG e RBS precisam da versão 2.0 ou mais nova (o SPLACE já traz a 2.0.4).';
  }

  @override
  String get raxmlModelsHeading => 'Modelos de substituição';

  @override
  String get raxmlMoose => 'MOOSE (escolha de modelos do próprio RAxML-NG)';

  @override
  String get raxmlMooseBody =>
      'Um modelo por partição, escolhido na mesma execução. Recomendado.';

  @override
  String get raxmlModelFinderBody =>
      'O IQ-TREE escolhe os modelos e pode juntar partições; depois o RAxML-NG faz a árvore.';

  @override
  String get raxmlFixedModel => 'Um modelo fixo';

  @override
  String get raxmlFixedBody =>
      'O mesmo modelo em todas as partições (por exemplo GTR+G).';

  @override
  String get raxmlCriterion => 'Critério';

  @override
  String get raxmlSearch => 'Busca da árvore';

  @override
  String get raxmlSearchAdaptive => 'Adaptativa (padrão)';

  @override
  String get raxmlSearchFast => 'Rápida (--fast)';

  @override
  String get raxmlSupportInfo =>
      'FBP, TBE e RBS vêm de árvores de réplicas de bootstrap. SH-aLRT e EBG são calculados só na árvore de ML e levam segundos. FBP e RBS não podem ser usados juntos.';

  @override
  String get treeLayoutLabel => 'Formato';

  @override
  String get treeLayoutRect => 'Retangular';

  @override
  String get treeLayoutCircle => 'Circular';

  @override
  String get treeLayoutInfo =>
      'Circular: a raiz no centro e as pontas em volta, bom para árvores com muitas pontas. Os nomes seguem o raio; as cores de táxon viram um anel por fora dos nomes e as barras de clado viram arcos.';

  @override
  String get treeArc => 'Abertura';

  @override
  String get treeArcInfo =>
      'Quanto do círculo as pontas ocupam: 360° fecha o círculo, 180° faz um semicírculo (leque).';

  @override
  String get treeRotation => 'Rotação';

  @override
  String get treeSwapSides => 'Trocar lados';

  @override
  String get treeSwapSidesText => 'Esta árvore à direita';

  @override
  String get treeSwapSidesInfo =>
      'Coloca a árvore em edição à direita e a outra à esquerda.';

  @override
  String get iqAlrt => 'Teste SH-aLRT';

  @override
  String get iqAlrtInfo =>
      'Um segundo valor de suporte, rápido, para cada ramo (mostrado como SH-aLRT/UFBoot). Ramos com SH-aLRT ≥ 80 e UFBoot ≥ 95 são bem sustentados.';

  @override
  String get iqAlrt1000 => '1000 réplicas (recomendado)';

  @override
  String get iqSearchHeading => 'Busca da árvore';

  @override
  String get iqSearchInfo =>
      'Os valores padrão servem para a maioria dos dados. Mais execuções ou mais iterações demoram mais, mas exploram mais árvores; útil quando o resultado muda entre execuções.';

  @override
  String get iqRuns => 'Execuções independentes';

  @override
  String get iqRunsInfo =>
      '--runs: repete a busca inteira a partir de árvores iniciais diferentes e guarda a melhor árvore.';

  @override
  String get iqNstop => 'Parar após (iterações)';

  @override
  String get iqNstopInfo =>
      '-nstop: a busca para depois deste número de iterações sem uma árvore melhor (padrão 100). Use mais para dados grandes ou difíceis.';

  @override
  String get iqPers => 'Força de perturbação';

  @override
  String get iqPersInfo =>
      '-pers: quanto a árvore atual é alterada a cada iteração (padrão 0,5). Valores menores ajudam em alinhamentos curtos.';

  @override
  String get iqSeedInfo =>
      'A mesma semente repete o mesmo resultado. Vazio: semente aleatória (anotada no log).';

  @override
  String get iqSeedRandom => 'aleatória';

  @override
  String get iqBnni => 'Refinar as árvores do UFBoot (--bnni)';

  @override
  String get iqBnniBody =>
      'Otimiza cada árvore de bootstrap com movimentos NNI. Recomendado: reduz suportes superestimados quando o modelo não se ajusta perfeitamente aos dados.';

  @override
  String get bothTitle => 'Rodar o IQ-TREE e o RAxML-NG juntos';

  @override
  String get bothInfo =>
      'Os dois são métodos de máxima verossimilhança. Comparar as árvores (etapa 9, Comparar árvores) mostra quais clados não dependem do programa.';

  @override
  String bothBody(String iq, String rx) {
    return 'As duas análises começam ao mesmo tempo, com as opções escolhidas nos dois cards acima. As threads são divididas para o computador continuar usável: o IQ-TREE fica com $iq e o RAxML-NG com $rx. O total se muda nas configurações (barra superior).';
  }

  @override
  String bothThreads(String n) {
    return '$n threads';
  }

  @override
  String get runBoth => 'Rodar os dois';

  @override
  String get bothRunning => 'IQ-TREE e RAxML-NG';

  @override
  String bothRunningSub(String threads, String time) {
    return 'Lado a lado · threads $threads · $time';
  }

  @override
  String get satUnlocksStep7 =>
      'Rode o teste para continuar: a sugestão dele (quais posições do códon manter) é usada para montar a matriz na etapa 7, que se abre em seguida.';

  @override
  String get satAllGenes => 'Todos os genes juntos';

  @override
  String get satAllGenesInfo =>
      '\"Todos os genes juntos\" une os genes codificantes só para este teste (a matriz em si é montada na etapa 7). É essa linha que decide a sugestão de quais posições do códon tirar da matriz inteira; a linha de cada gene permite escolher gene a gene.';

  @override
  String get exampleLabel => 'Exemplo';

  @override
  String get exampleNt =>
      'O mesmo alinhamento nos três modos. Aqui o gap da espécie B cai dentro dos códons 3 e 4 e desloca o quadro de leitura.';

  @override
  String get exampleAa =>
      'Os mesmos genes traduzidos: uma letra por códon, logo menos sítios e mais conservados.';

  @override
  String get exampleCodon =>
      'O gap da espécie B é o códon 3 inteiro: o quadro de leitura é preservado.';

  @override
  String exampleTrim(String n, String total) {
    return 'A trimagem remove colunas pouco confiáveis (quase só gaps, ou sem sinal filogenético) antes de construir a árvore. Aqui $n de $total colunas são removidas (×).';
  }

  @override
  String get exampleCodonKeepAll => 'As três posições do códon são mantidas.';

  @override
  String exampleCodonKept(String positions, String pct) {
    return 'Mantidas: posições $positions ($pct% dos sítios codificantes).';
  }

  @override
  String get exampleSatLabel => 'Como ler os gráficos';

  @override
  String get exampleSatOk => 'Sem saturação';

  @override
  String get exampleSatOkBody =>
      'As diferenças observadas (p) continuam subindo com a distância (d), perto da linha tracejada y = x.';

  @override
  String get exampleSatBad => 'Saturado';

  @override
  String get exampleSatBadBody =>
      'p se estabiliza enquanto d continua crescendo: mudanças novas apagam as antigas e o sinal se perde.';

  @override
  String get plotTrend => 'tendência';

  @override
  String exampleSpecies(String letter) {
    return 'Espécie $letter';
  }

  @override
  String get examplePlayAgain => 'Ver de novo';

  @override
  String get trimStepOriginal => 'Alinhamento';

  @override
  String get trimStepMarked => 'Colunas removidas';

  @override
  String get trimStepResult => 'Alinhamento trimado';

  @override
  String get treeSupportArrowLabel => 'Setas';

  @override
  String get treeSupportArrow => 'Valor com seta para o ramo';

  @override
  String get treeSupportArrowInfo =>
      'Afasta cada valor de suporte um pouco dos ramos (para cima ou para baixo, como em Posição do suporte) e o liga ao seu nó com uma pequena seta. Útil quando os valores caem sobre ramos curtos ou em cima de outras linhas.';

  @override
  String get orgOpen => 'Amostra e marcador';

  @override
  String get orgTitle => 'Amostra e marcador nos cabeçalhos';

  @override
  String orgSubtitle(String files, String seqs) {
    return '$files arquivo(s) FASTA · $seqs sequências';
  }

  @override
  String orgSummary(String samples, String markers, String pct) {
    return '$samples amostras × $markers marcadores · $pct% preenchido';
  }

  @override
  String orgApply(String n) {
    return 'Gerar $n FASTAs de marcador';
  }

  @override
  String get orgIntro =>
      'Cada cabeçalho precisa dizer a qual amostra (ID) e a qual marcador a sequência pertence. Escolha os caracteres que separam as partes do cabeçalho e clique nas partes abaixo: uma vez para amostra, duas para marcador, uma terceira para ignorar. O SPLACE então gera um FASTA por marcador com todas as amostras, pronto para alinhar.';

  @override
  String get orgSampleFrom => 'Amostra (ID) vem de';

  @override
  String get orgSampleFromInfo =>
      'Uma parte do cabeçalho, ou o nome do arquivo quando cada arquivo tem uma amostra e os cabeçalhos só dizem o marcador.';

  @override
  String get orgSampleHeader => 'Cabeçalho';

  @override
  String get orgSampleFile => 'Nome do arquivo';

  @override
  String get orgDelimitersInfo =>
      'Cada caractere ligado separa o cabeçalho. Vários podem ficar ligados quando os arquivos usam separadores diferentes (ITV01-L1 e ITV02_L1).';

  @override
  String get orgSpace => 'espaço';

  @override
  String get orgOther => 'outros';

  @override
  String get orgParts =>
      'Partes do cabeçalho (exemplos de arquivos diferentes)';

  @override
  String get orgPartsInfo =>
      'As posições valem para todos os cabeçalhos. Uma parte que vai até o fim do cabeçalho é contada a partir do fim, então cabeçalhos com mais ou menos partes também funcionam.';

  @override
  String get orgLegendSample => 'Amostra (ID)';

  @override
  String get orgLegendMarker => 'Marcador';

  @override
  String get orgLegendIgnored => 'Ignorado';

  @override
  String get orgAdvanced => 'Avançado: tags e expressões regulares';

  @override
  String get orgAdvancedActive =>
      'Uma tag ou expressão regular está em uso (veja Avançado); o clique nas partes fica desligado.';

  @override
  String orgImplausible(String markers, String seqs) {
    return '$markers marcadores para $seqs sequências';
  }

  @override
  String get orgImplausibleBody =>
      'Quase cada sequência virou um marcador próprio, então a parte do marcador provavelmente está errada (por exemplo, o cabeçalho inteiro). Escolha o delimitador e clique na parte que nomeia o marcador.';

  @override
  String orgUnassigned(String n) {
    return '$n sequência(s) sem amostra ou sem marcador (ficam de fora)';
  }

  @override
  String get orgUnassignedLabel => 'Sem amostra ou marcador';

  @override
  String get orgVariations => 'Variações de nome';

  @override
  String get orgVariationsInfo =>
      'A mesma amostra ou marcador pode estar escrita de jeitos um pouco diferentes. Estas opções juntam as variações; você também pode renomear qualquer amostra na lista abaixo.';

  @override
  String get orgIgnoreCase => 'Ignorar maiúsculas e minúsculas';

  @override
  String get orgLeadingZeros => 'Ignorar zeros à esquerda nos números';

  @override
  String orgFileMismatch(String n) {
    return '$n arquivo(s) com nome diferente da amostra nos cabeçalhos';
  }

  @override
  String orgRenameTo(String name) {
    return 'Chamar de $name';
  }

  @override
  String get orgRenameAll => 'Usar o nome do arquivo nestas';

  @override
  String get orgUseFileNames => 'Tirar todas as amostras do nome do arquivo';

  @override
  String orgMerged(String n) {
    return '$n amostra(s) escrita(s) de mais de um jeito, unidas:';
  }

  @override
  String orgSamples(String n) {
    return 'Amostras ($n)';
  }

  @override
  String get orgSamplesInfo =>
      'Clique em um nome para renomear (Enter aplica). Dar o mesmo nome a duas amostras junta as duas.';

  @override
  String get orgFindSample => 'Buscar amostra';

  @override
  String orgSampleMarkers(String n, String total) {
    return '$n de $total marcadores';
  }

  @override
  String orgMarkers(String n) {
    return 'Marcadores ($n)';
  }

  @override
  String get orgMarkersInfo =>
      'Quantas amostras cada marcador tem, e o que fazer quando uma amostra tem mais de uma sequência do mesmo marcador. Todos os marcadores vão para a etapa 4, onde você escolhe quais usar.';

  @override
  String orgDuplicates(String n) {
    return '$n caso(s) de amostra com mais de uma sequência do mesmo marcador';
  }

  @override
  String get orgCoverage => 'Marcadores por número de amostras';

  @override
  String orgCoverageBar(String markers, String samples) {
    return '$markers marcadores em $samples amostras';
  }

  @override
  String get orgSamplesAxis => 'amostras';

  @override
  String orgCoverageSummary(String complete, String half) {
    return '$complete em todas as amostras · $half em pelo menos metade';
  }

  @override
  String get orgPending =>
      'Diga ao SPLACE qual parte dos cabeçalhos é a amostra e qual é o marcador';

  @override
  String orgProblem(String problem) {
    return 'A divisão não foi usada: $problem';
  }

  @override
  String get orgPendingBody =>
      'Depois é gerado um FASTA por marcador com todas as amostras, e os marcadores aparecem na etapa 4.';

  @override
  String get orgMatrix => 'Matriz';

  @override
  String get orgWritten => 'FASTA por marcador';

  @override
  String get findMarker => 'Buscar marcador';

  @override
  String resultsShowing(String n, String total) {
    return 'Mostrando $n de $total';
  }

  @override
  String showMoreN(String n) {
    return 'Mostrar mais $n';
  }

  @override
  String get showAll => 'Mostrar todos';

  @override
  String get noMatches => 'Nada encontrado.';

  @override
  String orgFormats(String n) {
    return '$n formatos de cabeçalho';
  }

  @override
  String perLocusCount(String positions, String n) {
    return '$positions em $n marcadores';
  }

  @override
  String get satOptionalLoci =>
      'Opcional para loci: a etapa 7 já está liberada. Os marcadores marcados como codificantes são testados por posição do códon; os demais, em todos os sítios.';

  @override
  String get satOptionalLociNoCoding =>
      'Opcional para loci: a etapa 7 já está liberada. Nenhum marcador é codificante, então cada locus é testado em todos os sítios (não há posições de códon): o teste mostra quais loci estão saturados.';

  @override
  String exampleTrimNone(String total) {
    return 'Nenhuma das $total colunas é removida: um alinhamento do MAFFT não tem colunas formadas só por gaps. Esta opção só faz diferença em alinhamentos editados à mão ou em subconjuntos de sequências.';
  }

  @override
  String get exampleGapsRow => 'gaps';

  @override
  String get runFailed => 'A análise falhou.';

  @override
  String get runLastLines => 'Últimas linhas do programa';

  @override
  String get runOpenLog => 'Abrir o log';

  @override
  String orgHeaderFormats(String n) {
    return 'Formatos de cabeçalho ($n)';
  }

  @override
  String get orgHeaderFormatsInfo =>
      'Os cabeçalhos estão escritos de mais de um jeito (outros delimitadores, mais ou menos partes). Cada formato mostra as suas sequências e pode ter as próprias partes de amostra e marcador: clique nas partes do exemplo dele. No fim, todas as sequências ficam com o mesmo padrão: a amostra no cabeçalho e um FASTA por marcador.';

  @override
  String orgFormatCount(String n, String files) {
    return '$n sequências em $files arquivos';
  }

  @override
  String get orgFormatOwnRule => 'Partes próprias';

  @override
  String get orgFormatSkipped => 'Deixado de fora';

  @override
  String orgFormatShow(String n) {
    return 'Mostrar as $n sequências';
  }

  @override
  String get orgFormatHide => 'Ocultar as sequências';

  @override
  String get orgFormatUseMain => 'Mesmas partes da regra principal';

  @override
  String get orgFormatSkip => 'Deixar estas sequências de fora';

  @override
  String get orgFormatUse => 'Usar estas sequências';

  @override
  String orgNameStructures(String n) {
    return 'Estruturas dos nomes de amostra ($n)';
  }

  @override
  String get orgNameStructuresInfo =>
      'Nomes de amostra agrupados pelo jeito como estão escritos: letras como estão, cada dígito como #. Nomes de uma coleção escritos com mais ou menos dígitos (ITV1054 e ITV01054) podem ser padronizados para a estrutura que você escolher: zeros são adicionados ou removidos à esquerda. Passe o mouse no botão para ver as renomeações.';

  @override
  String orgStructureCount(int n) {
    String _temp0 = intl.Intl.pluralLogic(
      n,
      locale: localeName,
      other: '$n amostras',
      one: '1 amostra',
    );
    return '$_temp0';
  }

  @override
  String orgStandardize(String structure) {
    return 'Padronizar como $structure';
  }

  @override
  String get satAllSitesNonCoding => 'todos os sítios (não codificante)';

  @override
  String schemeMixed(String coding, String other) {
    return '$coding marcadores codificantes divididos por códon; $other sem códons ficam inteiros';
  }

  @override
  String get exampleCodonGene => 'Gene';

  @override
  String get exampleCodonMatrix => 'Na matriz';

  @override
  String get mtTitle => 'Táxons × marcadores';

  @override
  String get mtIntro =>
      'Como no SequenceMatrix: cada célula é o número de bases que o táxon tem naquele marcador (verde mais escuro = mais completo), \"—\" quando não há sequência. O cabeçalho mostra o tamanho de cada alinhamento, quantos táxons o têm e os gaps. Passe o mouse numa célula ou cabeçalho para ver detalhes; a mesma tabela é salva em matrix_overview.tsv.';

  @override
  String get mtTaxon => 'Táxon';

  @override
  String get mtLoci => 'Marcadores';

  @override
  String get mtBases => 'Bases';

  @override
  String get mtMissing => 'Faltando';

  @override
  String get mtPartitions => 'Partições';

  @override
  String get mtLegend =>
      'Célula: bases do táxon no marcador · — sem sequência · cabeçalho: tamanho · táxons com o marcador · % de gaps';

  @override
  String mtLocusInfo(
    String locus,
    String length,
    String n,
    String total,
    String gaps,
  ) {
    return '$locus: $length pb, em $n de $total táxons, $gaps% de gaps';
  }

  @override
  String get mtMissingCell => 'sem sequência (preenchido com ?)';

  @override
  String mtCellInfo(String bases, String length, String gaps) {
    return '$bases de $length pb, $gaps gaps';
  }

  @override
  String get searchGenomes => 'Genomas';

  @override
  String get searchMarkers => 'Marcadores individuais';

  @override
  String get markerIncludeGenomes => 'Incluir genomas';

  @override
  String get markerGenomesInfo =>
      'Desligado: só registros de até 10 kb (sequências de genes, como barcodes). Ligado: também genomas completos e fragmentos longos, que também contêm os marcadores (e são bem maiores para baixar).';

  @override
  String markerPick(String n) {
    return 'Marcadores (SynGenes) · $n selecionados';
  }

  @override
  String get markerPickInfo =>
      'Genes codificantes e rRNAs das tabelas do SynGenes. Cada marcador é buscado com todos os nomes que o SynGenes conhece (símbolos e nomes do produto, ex.: COI, COX1, cytochrome c oxidase subunit I), então registros anotados de jeitos diferentes são encontrados. Registros da mesma espécie viram uma amostra só, e espécies sem um marcador ficam com dado faltante nele.';

  @override
  String reviewMarkerCount(
    String marker,
    String selected,
    String found,
    String species,
  ) {
    return '$marker: $selected/$found registros · $species espécies';
  }

  @override
  String get reviewOnePerSpecies => 'Um por espécie e marcador';

  @override
  String get reviewOnePerSpeciesInfo =>
      'Mantém, para cada marcador, o registro mais longo de cada espécie (registros verificados primeiro). É a seleção com que a revisão começa.';

  @override
  String get reviewColMarkers => 'Marcadores';

  @override
  String get onePerTaxon => 'Uma amostra por espécie (marcadores individuais)';

  @override
  String get onePerTaxonBody =>
      'Junta registros da mesma espécie (um COI de um voucher, um CYTB de outro) numa linha da matriz: o cabeçalho é família_espécie e, para cada marcador, fica a sequência mais longa da espécie. Espécies sem um marcador ficam com dado faltante (?).';

  @override
  String onePerTaxonJoined(String n) {
    return '$n espécies têm vários registros: cada uma vira uma amostra só.';
  }

  @override
  String satOptionalPartial(String markers) {
    return 'Opcional aqui: a etapa 7 já está liberada. $markers vieram de registros de marcadores individuais com sequências parciais (fragmentos como barcodes, marcados com < ou > no GenBank), então o teste serve de orientação, não é obrigatório.';
  }

  @override
  String satPartialMixed(String markers) {
    return '$markers têm sequências parciais: para eles o teste é opcional. A etapa 7 espera só pelos marcadores completos.';
  }

  @override
  String metaRecords(String n) {
    return '$n registros';
  }

  @override
  String headerConflictBlock(String n) {
    return '$n cabeçalho(s) se repetem dentro do mesmo marcador, então essas sequências não poderiam ser diferenciadas no alinhamento nem na árvore. Inclua um campo que as diferencie (Acesso ou Voucher) ou ligue \"Uma amostra por espécie\" para juntá-las. O alinhamento fica bloqueado até todos os cabeçalhos serem únicos.';
  }

  @override
  String get headerAddAccession => 'Incluir o Acesso';

  @override
  String get alignBlockedHeader => 'Bloqueado: cabeçalhos FASTA repetidos (5a)';

  @override
  String recGenesCount(String n, String total) {
    return 'Saturada em $n de $total genes.';
  }

  @override
  String recBasedOnGene(String gene, String taxa, String length) {
    return 'Resultado do $gene ($taxa táxons, $length pb). Cada posição é testada separadamente.';
  }

  @override
  String recBasedOnGenes(String n) {
    return '$n genes codificantes, cada um testado sozinho (nada é concatenado para o teste). Uma posição conta como saturada quando está saturada em mais da metade dos genes.';
  }

  @override
  String rec7bPreview(String section) {
    return 'Como \"$section\" vai ficar na etapa 7b:';
  }

  @override
  String get frameExampleLabel => 'Exemplo · manter o quadro de leitura';

  @override
  String get frameExampleCaption =>
      'O trimAl remove colunas soltas. Com esta opção, quando o trimAl remove uma coluna (×), o SPLACE remove o códon inteiro a que ela pertence (+): o gene continua com códons inteiros e as posições 1, 2 e 3 seguem na ordem (a régua colorida).';

  @override
  String get frameStepOriginal => 'Códons 1-2-3';

  @override
  String get frameStepColumn => 'O trimAl remove uma coluna';

  @override
  String get frameStepCodon => 'O ClipKIT --codon tira o códon inteiro';

  @override
  String get frameStepResult => 'Quadro mantido';

  @override
  String get trimCasesLabel =>
      'Exemplo · aparar um gene codificante, de duas formas';

  @override
  String get trimCasesStepMark => 'Colunas marcadas';

  @override
  String get trimCaseColumns => 'Coluna a coluna (trimAl)';

  @override
  String get trimCaseCodons => 'Manter o quadro de leitura (ClipKIT --codon)';

  @override
  String get trimCaseSelected => 'selecionado';

  @override
  String trimCaseFrameLost(String n) {
    return 'Quadro perdido: $n colunas fora da posição do códon';
  }

  @override
  String trimCaseFrameKept(String k, String total) {
    return 'Quadro mantido: $k de $total códons, posições 1-2-3 em ordem';
  }

  @override
  String get trimCaseFrameKeptCols =>
      'Quadro mantido: nenhuma coluna fora de posição';

  @override
  String get trimCasesCaption =>
      'Mesmo alinhamento, duas ferramentas. O trimAl remove uma a uma as colunas que marca (×), então as posições do códon (régua colorida) saem de ordem (!). O ClipKIT com --codon remove o códon inteiro de cada coluna que marca (+): as 1ª, 2ª e 3ª posições continuam em ordem para as partições e o teste de saturação.';

  @override
  String get subRecordsTable => 'Tabela de registros';

  @override
  String get subRecordsTableBody =>
      'Cada registro com sua taxonomia e contagem de genes; edite ou remova registros aqui.';

  @override
  String get subFeatureTypes => 'Anotações usadas (CDS, rRNA, tRNA)';

  @override
  String get subChooseMarkers => 'Escolher os marcadores';

  @override
  String get subRunSaturation => 'Rodar o teste';

  @override
  String get subBuildMatrix => 'Montar a matriz';

  @override
  String get subTree => 'Árvore e figura';

  @override
  String get viewCodons => 'Códons';

  @override
  String get viewDots => 'Pontos';

  @override
  String get viewDotsInfo =>
      'Mostra um ponto no lugar das bases nas colunas em que todas as sequências têm a mesma base (a primeira sequência mantém as letras), para destacar as colunas variáveis, como no MEGA e em outros visualizadores. Na visão por códons, os códons iguais em todas as sequências viram pontos.';

  @override
  String codonsVarying(String n, String total, String nonsyn) {
    return '$n de $total códons variam · $nonsyn mudam o aminoácido';
  }

  @override
  String get codonVariation => 'Variação do códon';

  @override
  String viewerCodonN(String n) {
    return 'códon $n';
  }

  @override
  String codonColumnInfo(String codons, String aa) {
    return '$codons códon(s) diferente(s) aqui, $aa aminoácido(s)';
  }

  @override
  String get codonLegendAa => 'cor = aminoácido';

  @override
  String get codonLegendSyn => 'mudança sinônima';

  @override
  String get codonLegendNonSyn => 'muda o aminoácido';

  @override
  String get mtResizeNames =>
      'Arraste para alargar os nomes; clique duas vezes para caber o mais longo';

  @override
  String treeLegendUnmarked(String min) {
    return 'sem marca: suporte < $min';
  }

  @override
  String get treeStrongFrom => 'Suporte forte a partir de';

  @override
  String get treeModerateFrom => 'Suporte moderado a partir de';

  @override
  String get treeSupportClassesInfo =>
      'Cada classe colore os pontos (e valores) do seu limite até o próximo. Escolha um modelo — automático para o tipo de suporte, de 2 a 5 classes, ou cortado nos tercis dos valores desta árvore — ou defina você mesmo os limites e as cores, adicionando ou removendo classes.';

  @override
  String get treeLegendSupport => 'Suporte dos ramos';

  @override
  String get treeLegendSupportShow => 'Na legenda';

  @override
  String get treeLegendSupportTitle => 'Título do suporte';

  @override
  String get badgeBrazil => 'só Brasil';

  @override
  String get salvePossiblyExtinct => 'Possivelmente extinta';

  @override
  String get salveSheetSubtitle =>
      'Ficha do SALVE (ICMBio) · avaliação nacional do risco de extinção';

  @override
  String get salveOpenSheet => 'Abrir a ficha no SALVE';

  @override
  String get salveWatersheds => 'Bacias';

  @override
  String get salvePublished => 'Publicada em';

  @override
  String get salveSituation => 'Situação da ficha';

  @override
  String get salveSynonyms => 'Nomes antigos';

  @override
  String get salveJustification => 'Justificativa';

  @override
  String get salveAssessors => 'Avaliadores';

  @override
  String get salveCitation => 'Como citar';

  @override
  String get colCriteria => 'Critérios';

  @override
  String get colCommonNames => 'Nomes comuns';

  @override
  String get colSalveGroup => 'Grupo SALVE';

  @override
  String get colEndemic => 'Endêmica do Brasil';

  @override
  String get colBiomes => 'Biomas';

  @override
  String get colStates => 'Estados';

  @override
  String get colAssessed => 'Avaliada em';

  @override
  String get colSheet => 'Ficha';

  @override
  String get geneNamesAllKnown => 'Todos os nomes de genes reconhecidos';

  @override
  String geneReviewAgain(String n) {
    return 'Revisar de novo ($n)';
  }

  @override
  String get clipkitCodonMode =>
      'Marcadores codificantes: modo do ClipKIT (--codon)';

  @override
  String get clipkitCodonInfo =>
      'Para manter o quadro de leitura, os marcadores codificantes são aparados pelo ClipKIT com a opção --codon: quando uma posição de um códon deve sair, o ClipKIT remove o códon inteiro (Steenwyk et al. 2020). Os marcadores não codificantes são aparados pelo trimAl com o método acima.';

  @override
  String get outgroupSuggesting =>
      'Medindo as distâncias entre os táxons dos alinhamentos para sugerir um grupo externo…';

  @override
  String get outgroupSuggestNone =>
      'Alinhe os marcadores para receber uma sugestão de grupo externo pelas distâncias.';

  @override
  String get outgroupSuggestTitle =>
      'Grupo externo sugerido (distâncias dos alinhamentos)';

  @override
  String get outgroupSuggestWeak =>
      'Nenhum grupo claramente separado: o táxon mais distante é sugerido';

  @override
  String get outgroupSuggestInfo =>
      'Distâncias p entre cada par de táxons nos sítios que ambos têm em cada alinhamento, somadas entre os marcadores; uma árvore de ligação média (UPGMA) dessas distâncias; o menor dos dois grupos unidos por último é o mais distante do resto. O valor ao lado de cada táxon é a sua distância média a todos os outros. Confira com o que você sabe do grupo: um ramo longo (sequência de evolução rápida ou contaminada) também pode parecer distante.';

  @override
  String get outgroupSuggestApplied => 'Em uso';

  @override
  String get outgroupSuggestUse => 'Usar como grupo externo';

  @override
  String outgroupSuggestStats(
    String between,
    String within,
    String loci,
    String taxa,
  ) {
    return 'Distância p média ao grupo interno $between; dentro do grupo interno $within · $loci alinhamentos, $taxa táxons.';
  }

  @override
  String get outgroupSuggestNext => 'Próximos mais distantes';

  @override
  String get treeArrowShape => 'Forma da seta';

  @override
  String get treeArrowStraight => 'Reta';

  @override
  String get treeArrowCurved => 'Curva';

  @override
  String get treeArrowShapeInfo =>
      'Os valores nunca se sobrepõem: quando o lugar de um valor está ocupado ele se afasta e uma seta aponta para o seu nó. As setas podem ser retas ou curvas.';

  @override
  String get treeSupportClasses => 'Classes de suporte';

  @override
  String treeClassesAuto(String kind) {
    return 'Automático ($kind)';
  }

  @override
  String treeClassesN(String n) {
    return '$n classes';
  }

  @override
  String get treeClassesTertiles => 'Tercis desta árvore';

  @override
  String get treeClassHidden => 'abaixo de \"mostrar suporte a partir de\"';

  @override
  String get treeClassRemove => 'Remover esta classe';

  @override
  String get treeClassAdd => 'Adicionar classe';

  @override
  String get tipFieldConservation => 'Categoria (SALVE)';

  @override
  String get treeHeatStatus => 'Categoria';

  @override
  String get treeHeatStatusLegend => 'Risco de extinção (SALVE / ICMBio)';

  @override
  String get treeHeatYes => 'sim';

  @override
  String get treeHeatNo => 'não';

  @override
  String get treeHeatmap => 'Heatmap do SALVE';

  @override
  String get treeHeatmapInfo =>
      'Colunas coloridas ao lado dos rótulos a partir das fichas do SALVE das espécies: categoria nas cores da IUCN, endemismo no Brasil e uma coluna por bioma. Só no layout retangular. A mesma informação pode entrar nos rótulos como texto (partes do rótulo acima).';

  @override
  String get satParamsTitle => 'Parâmetros do teste';

  @override
  String get satParamsBody =>
      'Como a saturação de cada marcador é medida e quando ela é considerada moderada ou saturada. As mudanças valem na próxima execução.';

  @override
  String get satParamMethod => 'Método';

  @override
  String get satParamMethodValue =>
      'saturação do PhyKIT (Steenwyk et al. 2021; Philippe et al. 2011): inclinação da distância p não corrigida de cada par de táxons sobre a distância patrística na árvore';

  @override
  String get satParamDistance => 'Árvores';

  @override
  String get satParamDistanceValue =>
      'IQ-TREE, uma por gene (busca rápida) com o modelo abaixo; os comprimentos de ramo de cada subconjunto de posições são estimados nessa topologia (-te)';

  @override
  String get satParamSubsets => 'Sítios testados';

  @override
  String get satParamSubsetsValue =>
      '1ª, 2ª e 3ª posições do códon, 1ª+2ª e todos os sítios (marcadores não codificantes: todos os sítios)';

  @override
  String get satParamCode => 'Código genético';

  @override
  String get satParamFrame => 'Quadro de leitura';

  @override
  String get satParamFrameValue =>
      'do alinhamento (início do CDS, ou o quadro com menos stops internos)';

  @override
  String get satParamsThresholds => 'Limiares e ajustes';

  @override
  String get satParamSlopeModerate => 'Saturação · moderada <';

  @override
  String get satParamSlopeSaturated => 'Saturação · saturada <';

  @override
  String get satParamSlopeInfo =>
      'Valor de saturação do PhyKIT: inclinação da distância p sobre a distância patrística. 1 = sem saturação; quanto mais perto de 0, mais saturada. Os padrões (0,5 e 0,3) são heurísticos: as 2ª posições do códon costumam ficar perto de 1 e 3ª posições saturadas bem abaixo de 0,5.';

  @override
  String get satParamIss => 'Iss · saturada ≥';

  @override
  String get satParamIssInfo =>
      'Índice de Xia: entropia média dos sítios sobre a entropia esperada sob saturação total. Perto de 1, os sítios parecem aleatórios.';

  @override
  String get satParamUndefinedModerate => 'Sem correção · moderada ≥';

  @override
  String get satParamUndefinedSaturated => 'Sem correção · saturada ≥';

  @override
  String get satParamUndefinedInfo =>
      'Pares cuja distância corrigida não pode ser calculada (além do que o modelo corrige): forte sinal de saturação.';

  @override
  String get satParamLowDivergence => 'Divergência baixa · p <';

  @override
  String get satParamLowDivergenceInfo =>
      'Com uma distância p média tão baixa, a saturação não importa e o subconjunto é considerado OK.';

  @override
  String get satParamMinSites => 'Mínimo de sítios';

  @override
  String get satParamMinPairs => 'Mínimo de pares';

  @override
  String get satParamReplicates => 'Iss · réplicas';

  @override
  String get satParamReplicatesInfo =>
      'Sítios aleatórios sorteados para estimar a entropia esperada sob saturação total, para cada número de sequências. Mais réplicas, Iss mais estável e execução mais longa.';

  @override
  String get satParamPlotPoints => 'Pontos no gráfico';

  @override
  String get satParamPlotPointsInfo =>
      'Pares de sequências desenhados nos gráficos (uma amostra quando há mais); as estatísticas usam todos os pares.';

  @override
  String get satParamsReset => 'Restaurar padrões';

  @override
  String get satParamsStale =>
      'Os parâmetros mudaram desde a última execução: rode o teste de novo para atualizar os resultados.';

  @override
  String get colTrend => 'Tendência populacional';

  @override
  String get colGeneration => 'Tempo de geração';

  @override
  String get colMigratory => 'Migratória';

  @override
  String get colHabitats => 'Habitats';

  @override
  String get colThreats => 'Ameaças';

  @override
  String get colOfficialList => 'Lista oficial';

  @override
  String get colUcs => 'Unidades de conservação';

  @override
  String get salveReferences => 'Referências na ficha';

  @override
  String get salveHistory => 'Histórico de avaliações';

  @override
  String get salveLists => 'Listas oficiais e convenções';

  @override
  String get salveSecDistribution => 'Distribuição';

  @override
  String get salveNational => 'No Brasil';

  @override
  String get salveGlobal => 'No mundo';

  @override
  String get salveSecNatural => 'História natural e população';

  @override
  String get salveDiet => 'Hábito alimentar';

  @override
  String get salveReproduction => 'Reprodução';

  @override
  String get salveSexRatio => 'Razão sexual';

  @override
  String get salveNaturalHistory => 'História natural';

  @override
  String get salveInteractions => 'Interações com outras espécies';

  @override
  String get salvePopulation => 'População';

  @override
  String get salveGenetics => 'Características genéticas';

  @override
  String get salveSecThreats => 'Ameaças, usos, conservação e pesquisa';

  @override
  String get salveUses => 'Usos';

  @override
  String get salveActions => 'Ações de conservação';

  @override
  String get salveConservation => 'Conservação';

  @override
  String get salveResearch => 'Pesquisa';

  @override
  String get trendDeclining => 'declinando';

  @override
  String get trendStable => 'estável';

  @override
  String get trendIncreasing => 'aumentando';

  @override
  String get trendUnknown => 'desconhecida';

  @override
  String get treeHeatFew => 'poucas';

  @override
  String get treeHeatMany => 'muitas';

  @override
  String trimNothingRemoved(String method, String n) {
    return '$method manteve todas as $n colunas: este alinhamento não tem região com gaps para remover, então nada fica em cinza (os modos kpic também removem sítios não informativos)';
  }

  @override
  String get treeBranchWidth => 'Largura dos ramos';

  @override
  String get treeBranchWidthInfo =>
      'Largura da árvore em si (ramos), em cm. 0 = automático: o que sobra ao lado de rótulos, heatmap, barras e legenda, nunca menos que 45 % da figura — com muitos elementos a figura fica mais larga em vez de espremer a árvore. Defina uma largura para fixá-la; a figura cresce para caber o resto.';

  @override
  String get treeHeatCell => 'Célula do heatmap';

  @override
  String get treeHeatCellInfo =>
      'Lado de cada célula do heatmap, em pt. 0 = automático (pela altura da linha).';

  @override
  String get citeVersions => 'Versões usadas';

  @override
  String get satParamVersions => 'Programas';

  @override
  String get satParamModel => 'Modelo do IQ-TREE';

  @override
  String get satParamModelInfo =>
      'Modelo das árvores dos genes e dos comprimentos de ramo de cada subconjunto (IQ-TREE -m), ex.: GTR+G, HKY+G, GTR+F+I+G4.';

  @override
  String get satParamGaps => 'Gaps';

  @override
  String get satParamGapsSwitch => 'Ignorar sítios com gap (-e)';

  @override
  String get satParamGapsInfo =>
      'PhyKIT --exclude_gaps: sítios com gap em um par ficam fora da sua distância p.';

  @override
  String get phykitRequired =>
      'O teste de saturação usa o PhyKIT, que ainda não está instalado. Instale aqui (precisa do Python 3) ou com: pip install phykit';

  @override
  String get installPhykit => 'Instalar PhyKIT';

  @override
  String plotPhykit(String s) {
    return 'p × distância na árvore · PhyKIT $s';
  }

  @override
  String get plotPhykitBody =>
      'Cada ponto é um par de táxons: distância p não corrigida contra a distância patrística na árvore do IQ-TREE destes sítios. Pontos muito abaixo da diagonal escondem substituições repetidas.';

  @override
  String get plotPositionsPatristic =>
      'Posições do códon × distância na árvore';

  @override
  String get patristicDistance => 'Distância patrística';

  @override
  String get metricMeanPatristic => 'Patrística média';

  @override
  String get metricPhykit => 'Saturação PhyKIT';

  @override
  String get trimToolLabel => 'Ferramenta de trimagem';

  @override
  String get trimToolTrimal => 'trimAl · coluna a coluna';

  @override
  String get trimToolClipkit => 'ClipKIT · códons inteiros';

  @override
  String get trimToolTrimalBody =>
      'O trimAl apara todos os marcadores coluna a coluna. Em genes codificantes isso pode partir códons: o quadro de leitura se perde, eles viram uma partição sem posições do códon e ficam fora do teste de saturação.';

  @override
  String get trimToolClipkitBody =>
      'O ClipKIT apara todos os marcadores. Nos genes codificantes ele roda com --codon: quando uma posição de um códon deve sair, o códon inteiro sai, então o quadro de leitura e as posições do códon são mantidos. Marcadores não codificantes são aparados pelo ClipKIT sem --codon. O trimAl não é usado.';

  @override
  String get mtShowBases => 'Bases';

  @override
  String get mtShowPercent => '% do marcador';

  @override
  String get mtOrderMatrix => 'Ordem da matriz';

  @override
  String get mtOrderComplete => 'Mais completos primeiro';

  @override
  String get mtScaleComplete => 'Completude';

  @override
  String get mtScaleMissing => 'Faltantes / gaps';

  @override
  String get treeMarkShape => 'Forma dos pontos';

  @override
  String get treeShapeCircle => 'Círculo';

  @override
  String get treeShapeSquare => 'Quadrado';

  @override
  String get treeShapeDiamond => 'Losango';

  @override
  String get treeShapeTriangle => 'Triângulo';

  @override
  String get treeBarWidth => 'Espessura das barras de clado';

  @override
  String get treeBarWidthInfo =>
      'Espessura das barras ao lado dos clados, em pt. 0 = automático (pela espessura dos ramos).';

  @override
  String get treeBarFont => 'Tamanho dos nomes de clado';

  @override
  String get treeBarFontInfo =>
      'Tamanho dos nomes ao lado das barras de clado, em pt. 0 = o tamanho dos rótulos das pontas.';

  @override
  String get treeHeatGap => 'Espaço entre colunas do heatmap';

  @override
  String get trimToolTrimalCodon => 'trimAl · códons inteiros';

  @override
  String get trimToolTrimalCodonBody =>
      'O trimAl apara todos os marcadores. Os genes codificantes são lidos códon a códon e traduzidos com o código genético do projeto; o trimAl apara esse alinhamento de proteínas e, com -backtrans, devolve os códons dos aminoácidos mantidos, então códons inteiros saem ou ficam e as posições do códon são mantidas. Marcadores não codificantes são aparados coluna a coluna.';

  @override
  String frameStepCodonTool(String tool) {
    return 'O $tool tira o códon inteiro';
  }

  @override
  String trimCaseCodonsTool(String tool) {
    return 'Manter o quadro de leitura ($tool)';
  }

  @override
  String trimCasesCaptionTool(String tool) {
    return 'Mesmo alinhamento, duas formas. Coluna a coluna, o trimAl remove uma a uma as colunas que marca (×), então as posições do códon (régua colorida) saem de ordem (!). Com o $tool, sai o códon inteiro de cada coluna marcada (+): as 1ª, 2ª e 3ª posições continuam em ordem para as partições e o teste de saturação.';
  }
}
