// ignore: unused_import
import 'package:intl/intl.dart' as intl;

import 'app_localizations.dart';

// ignore_for_file: type=lint

/// The translations for Spanish Castilian (`es`).
class AppLocalizationsEs extends AppLocalizations {
  AppLocalizationsEs([String locale = 'es']) : super(locale);

  @override
  String get navProgress => 'Progreso';

  @override
  String get navSteps => 'Etapas del flujo';

  @override
  String get logTitle => 'Registro de actividad';

  @override
  String get logEmpty =>
      'Aún no hay actividad. Aquí aparecerán los registros de búsqueda, alineamiento y análisis.';

  @override
  String get logCopy => 'Copiar todo';

  @override
  String get logDownload => 'Descargar como .txt';

  @override
  String get logClear => 'Limpiar';

  @override
  String get logCollapse => 'Ocultar panel de registros';

  @override
  String get logWelcome =>
      'SPLACE está listo. Comience cargando datos en la etapa 1.';

  @override
  String get settingsTheme => 'Tema';

  @override
  String get themeLight => 'Claro';

  @override
  String get themeDark => 'Oscuro';

  @override
  String get themeAuto => 'Auto';

  @override
  String get settingsDensity => 'Densidad';

  @override
  String get densityComfy => 'Cómoda';

  @override
  String get densityCompact => 'Compacta';

  @override
  String get settingsFontSize => 'Tamaño del texto';

  @override
  String get fontSizeReset => 'Restablecer';

  @override
  String get minimapLocked => 'Complete la etapa anterior para desbloquear';

  @override
  String get tooltipMenu => 'Mostrar navegador del flujo';

  @override
  String get tooltipLogs => 'Mostrar registro de actividad';

  @override
  String get tooltipSettings => 'Visualización e idioma';

  @override
  String get tooltipHelp => 'Mostrar pantalla de bienvenida';

  @override
  String get tooltipTools =>
      'Herramientas externas y configuración de análisis';

  @override
  String get cancel => 'Cancelar';

  @override
  String get done => 'Listo';

  @override
  String andMore(String n) {
    return '…y $n más';
  }

  @override
  String get splashSummary =>
      'Un kit completo para filogenómica: obtenga registros de GenBank u organice loci FASTA, alinee, pruebe la saturación, concatene y construya árboles.';

  @override
  String get splashStart => 'Iniciar';

  @override
  String get splashStep1 =>
      'Cargue genomas o marcadores individuales (COI, 12S, CYTB…) del NCBI, o archivos GenBank y FASTA locales. [[SynGenes]] arma las búsquedas por marcador con todos los sinónimos de cada gen y estandariza los nombres de los marcadores para los pasos siguientes.';

  @override
  String get splashStep2 =>
      'Inspeccione los registros importados, ajuste los nombres de especies y recupere la taxonomía de GBIF, NCBI o SALVE con [[dataFishing]] antes de seleccionar marcadores.';

  @override
  String get splashStep3 =>
      'Elija qué tipos de feature (CDS, rRNA, tRNA) entran en el análisis, o indique a SPLACE dónde los encabezados FASTA guardan los nombres de muestra y locus.';

  @override
  String get splashStep4 =>
      'Elija los marcadores de la matriz y revise las secuencias ausentes y duplicadas en los mapas de presencia y longitud de genes.';

  @override
  String get splashStep5 =>
      'Alinee cada marcador con [[MAFFT]] y recórtelo con [[trimAl]] o [[ClipKIT]], manteniendo codones completos en los loci codificantes.';

  @override
  String get splashStep6 =>
      'Pruebe la saturación de sustituciones por posición del codón con el índice de [[Iss (saturation)|Xia]] y decida qué posiciones conservar.';

  @override
  String get splashStep7 =>
      'Concatene los loci alineados en una supermatriz y exporte archivos NEXUS, PHYLIP y de particiones, opcionalmente divididos por posición del codón.';

  @override
  String get splashStep8 =>
      'Estime la filogenia de máxima verosimilitud con [[IQ-TREE 3|IQ-TREE]] o [[RAxML-NG]], con soporte de ramas y árboles por gen opcionales.';

  @override
  String get step1 => 'Cargar datos';

  @override
  String get step2 => 'Registros cargados';

  @override
  String get step3 => 'Tipos de feature';

  @override
  String get step4 => 'Seleccionar marcadores';

  @override
  String get step5 => 'Alineamiento de secuencias';

  @override
  String get step6 => 'Prueba de saturación';

  @override
  String get step7 => 'Concatenar y NEXUS';

  @override
  String get step8 => 'Inferencia filogenética';

  @override
  String get step1Intro =>
      'Elija cómo cargar datos en SPLACE: registros GenBank de NCBI o de archivos locales, o sus propios archivos FASTA.';

  @override
  String get tipsKicker => 'Consejos';

  @override
  String get tipsTitle => 'Antes de empezar';

  @override
  String get tipsIntro =>
      'Construya la **tabla taxonómica** gradualmente y combine **fuentes complementarias** cuando sea necesario.';

  @override
  String get tipsExamplePrefix => 'Cada opción de abajo incluye un botón';

  @override
  String get tipsExampleSuffix => 'para ayudarle a empezar rápidamente.';

  @override
  String get tipSampleTitle => 'Un archivo por muestra';

  @override
  String get tipSampleBody =>
      'cada encabezado es un **locus** (p. ej. NC_030596_Oryza.fasta con >psbA, >matK). SPLACE los reagrupa en **un archivo por locus**.';

  @override
  String get tipLocusTitle => 'Un archivo por locus';

  @override
  String get tipLocusBody =>
      'cada encabezado es una **muestra** (p. ej. COI.fasta). Ya tiene el formato correcto.';

  @override
  String get tipMixedTitle => 'Encabezados mixtos';

  @override
  String get tipMixedBody =>
      'muestra y locus en el **mismo encabezado** (>Bufo_bufo_COI) o **etiquetas del NCBI** ([gene=COX1]).';

  @override
  String get tipAlignedTitle => '¿Ya alineado?';

  @override
  String get tipAlignedBody =>
      'Puede **omitir MAFFT** en la etapa 5 e ir directo al recorte, la saturación y los árboles.';

  @override
  String optionN(String n) {
    return 'Opción $n';
  }

  @override
  String get filesTitle => 'Usar archivos FASTA locales';

  @override
  String get filesInfo =>
      'Importe archivos .fasta, .fas, .fa o .fna que ya están en su computadora.';

  @override
  String get filesBody =>
      'Importe sus propios archivos .fasta, .fas, .fa o .fna desde el disco.';

  @override
  String get filesDrop => 'Arrastre y suelte archivos FASTA o carpetas aquí';

  @override
  String get filesBrowse => 'o haga clic para buscar (.fasta, .fas, .fa, .fna)';

  @override
  String get folderTitle => 'Importar una carpeta completa';

  @override
  String get folderInfo =>
      'Lee todos los archivos FASTA de la carpeta y sus subcarpetas.';

  @override
  String get folderBody => 'Útil para cientos de archivos, uno por muestra.';

  @override
  String get addFolder => 'Elegir carpeta';

  @override
  String inputSummary(String files, String sequences) {
    return '$files archivos · $sequences secuencias cargadas';
  }

  @override
  String get clearRecords => 'Limpiar registros';

  @override
  String get step2Intro =>
      'Todos los encabezados FASTA se listan abajo con la muestra y el locus que SPLACE lee de cada uno. Ajuste las reglas en la etapa 3 hasta que todos los registros tengan ambos.';

  @override
  String get detectedLayout => 'Formato detectado';

  @override
  String get layoutSamplePerFile =>
      'un archivo por muestra — los encabezados son loci';

  @override
  String get layoutLocusPerFile =>
      'un archivo por locus — los encabezados son muestras';

  @override
  String get layoutMixedHeader => 'muestra y locus están en el encabezado';

  @override
  String get layoutTaggedHeader =>
      'locus leído de las etiquetas del encabezado ([gene=…])';

  @override
  String get layoutUnknown => 'no reconocido — elija los campos manualmente';

  @override
  String confidence(String percent) {
    return 'confianza $percent%';
  }

  @override
  String get organizing => 'Actualizando…';

  @override
  String get tokenColumns =>
      'Partes del encabezado (separadas por los delimitadores del locus)';

  @override
  String tokenColumnInfo(String distinct, String geneLike) {
    return '$distinct distintos · $geneLike% parecen loci';
  }

  @override
  String get searchHeaders => 'Filtrar encabezados o archivos';

  @override
  String get colFile => 'Archivo';

  @override
  String get colHeader => 'Encabezado';

  @override
  String get colSample => 'Muestra';

  @override
  String get colLocus => 'Locus';

  @override
  String get colLength => 'Longitud';

  @override
  String get step3Intro =>
      'Elija de dónde vienen el nombre de la muestra y el del locus. La vista previa y la tabla de la etapa 2 se actualizan mientras escribe.';

  @override
  String fieldN(String n) {
    return 'Campo $n';
  }

  @override
  String get sampleField => 'Nombre de la muestra';

  @override
  String get sampleFieldInfo =>
      'Se convierte en el nombre de la secuencia en cada archivo de locus; debe ser idéntico entre loci para la concatenación.';

  @override
  String get locusField => 'Nombre del locus';

  @override
  String get locusFieldInfo =>
      'Define a qué archivo va cada secuencia. Sinónimos como COX1/COI/CO1 pueden unirse.';

  @override
  String get sourceFileName => 'Nombre del archivo';

  @override
  String get sourceHeader => 'Partes del encabezado';

  @override
  String get sourceTag => 'Etiqueta';

  @override
  String get sourceRegex => 'Regex';

  @override
  String get delimiters => 'Delimitadores';

  @override
  String get tokenStart => 'Primera parte';

  @override
  String get tokenEnd => 'Fin (vacío = última)';

  @override
  String get tagKey => 'Clave de la etiqueta';

  @override
  String get regexPattern => 'Patrón (se usa el primer grupo)';

  @override
  String get regexOnFileName => 'Aplicar al nombre del archivo';

  @override
  String get previewLabel => 'Vista previa';

  @override
  String get optionsKicker => 'Opciones';

  @override
  String get duplicates => 'Duplicados y nombres';

  @override
  String get duplicatesInfo =>
      'Qué hacer cuando la misma muestra tiene el mismo locus más de una vez (p. ej. copias de la repetición invertida en plastomas).';

  @override
  String get dupLongest => 'Mantener la más larga';

  @override
  String get dupMostComplete => 'Mantener la con menos N/gaps';

  @override
  String get dupFirst => 'Mantener la primera';

  @override
  String get dupAll => 'Mantener todas (renombrar)';

  @override
  String get mergeSynonyms => 'Unir sinónimos de locus (COX1 = COI = CO1…)';

  @override
  String get sanitizeNames => 'Nombres seguros para programas de árboles';

  @override
  String get selectAll => 'Seleccionar todos';

  @override
  String get selectNone => 'Ninguno';

  @override
  String get sharedOnly => 'Solo compartidos';

  @override
  String get legendSelectAll => 'selecciona todos los marcadores visibles';

  @override
  String get legendNone => 'borra todas las selecciones';

  @override
  String get legendShared => 'marcadores presentes en todas las muestras';

  @override
  String get missingTitle => 'Marcadores ausentes';

  @override
  String missingIn(String n) {
    return 'ausente en $n muestra(s)';
  }

  @override
  String get duplicatesTitle => 'Marcadores duplicados';

  @override
  String duplicatedIn(String n) {
    return 'duplicado en $n muestra(s)';
  }

  @override
  String get colCopies => 'Copias (tamaños)';

  @override
  String get colUsed => 'Usada';

  @override
  String get unassignedTitle => 'Registros sin muestra o locus';

  @override
  String get missingSample => 'sin muestra';

  @override
  String get missingLocus => 'sin locus';

  @override
  String useLoci(String n) {
    return 'Usar $n marcadores';
  }

  @override
  String get exportLoci => 'Exportar FASTA por locus…';

  @override
  String get heatmapTitle => 'Mapa de presencia de genes';

  @override
  String get heatmapAbsent => 'Ausente';

  @override
  String get heatmapPresent => 'Presente';

  @override
  String matrixSummary(String loci, String samples, String occupancy) {
    return '$loci marcadores × $samples muestras · $occupancy% completa. Haga clic en el nombre de un marcador para incluirlo o excluirlo.';
  }

  @override
  String get panelMafft => 'Configuración de MAFFT';

  @override
  String get panelMafftBody =>
      'Cada marcador se alinea por separado; varios marcadores se ejecutan en paralelo.';

  @override
  String get mafftStrategy => 'Método de alineamiento';

  @override
  String get directionLabel => 'Dirección de las secuencias';

  @override
  String get adjustDirection => 'Corregir secuencias invertidas';

  @override
  String get panelTrimal => 'Recorte (trimAl / ClipKIT)';

  @override
  String get trimalMethod => 'Modo de recorte';

  @override
  String get sequenceType => 'Tipo de secuencia';

  @override
  String get codonAwareTrimming => 'Mantener el marco de lectura';

  @override
  String get codonAwareHint =>
      'Eliminar columnas sueltas de alineamientos codificantes desplaza el marco de lectura e invalida las particiones por codón y la prueba de saturación.';

  @override
  String get panelRun => 'Ejecutar alineamiento';

  @override
  String panelRunBody(String dir) {
    return 'Alinee todos los marcadores seleccionados con MAFFT. La salida se guarda en $dir.';
  }

  @override
  String get runMafft => 'Alinear (MAFFT)';

  @override
  String get runTrimal => 'Recortar (trimAl)';

  @override
  String get useAsAlignment => 'La entrada ya está alineada';

  @override
  String get colCoding => 'Codificante';

  @override
  String get colFrame => 'Marco de lectura';

  @override
  String get colInclude => 'Usar este marcador';

  @override
  String get statRaw => 'Crudo';

  @override
  String get statAligned => 'Alineado';

  @override
  String get statTrimmed => 'Recortado';

  @override
  String get statusError => 'Error';

  @override
  String get saturationIntro =>
      'Cada marcador se prueba con PhyKIT (Philippe et al. 2011): IQ-TREE infiere el árbol del gen y, para la 1.ª, 2.ª y 3.ª posiciones del codón (y 1.ª+2.ª, todos los sitios), PhyKIT hace la regresión de la distancia p no corregida de cada par de taxones sobre su distancia en el árbol. Un valor cerca de 1 significa sin saturación; cuanto más cerca de 0, más esconden las diferencias observadas sustituciones repetidas.';

  @override
  String get saturationThresholds =>
      'Umbrales heurísticos, definidos en 6a (parámetros de la prueba): moderada por debajo del límite moderado del valor de saturación de PhyKIT, saturada por debajo del límite de saturación; una divergencia demasiado baja (p media) cuenta como OK.';

  @override
  String get runSaturation => 'Ejecutar prueba de saturación';

  @override
  String get recommendation => 'Recomendación';

  @override
  String get recKeepAll =>
      'Mantenga las tres posiciones del codón. En el paso 7b, particione la matriz por posición del codón para que cada una tenga su propio modelo de sustitución.';

  @override
  String get recAminoAcids =>
      'Las tres posiciones están saturadas: con esta divergencia los nucleótidos dicen poco sobre las relaciones. Alinee como aminoácidos (paso 5b) o use marcadores más lentos.';

  @override
  String get recUndetermined =>
      'Hay muy pocos sitios o pares de secuencias para evaluar la saturación. Mantenga todas las posiciones.';

  @override
  String get applyPerLocus => 'Usar el resultado de cada gen';

  @override
  String get concatenated => 'Concatenado';

  @override
  String get saturationHeatmapTitle =>
      'Saturación por locus y posición del codón';

  @override
  String get selectLocusHint =>
      'Pase el ratón por una celda para ver sus valores; haga clic en una fila para ver sus gráficos.';

  @override
  String get verdictOk => 'OK';

  @override
  String get verdictModerate => 'Moderada';

  @override
  String get verdictSaturated => 'Saturada';

  @override
  String get verdictInsufficient => 'Pocos datos';

  @override
  String plotObservedVsCorrected(String model) {
    return 'p observada vs distancia $model';
  }

  @override
  String plotTsTv(String model) {
    return 'Transiciones y transversiones vs distancia $model';
  }

  @override
  String get transitions => 'Transiciones';

  @override
  String get transversions => 'Transversiones';

  @override
  String get metricSites => 'Sitios';

  @override
  String get metricPairs => 'Pares';

  @override
  String get metricMeanP => 'p media';

  @override
  String get metricSlope => 'Pendiente';

  @override
  String get metricTsTv => 'Ts/Tv';

  @override
  String get metricIss => 'Iss';

  @override
  String get metricUndefined => 'Sin corrección';

  @override
  String get metricVerdict => 'Resultado';

  @override
  String frameLabel(String frame) {
    return 'codones empiezan en la base $frame';
  }

  @override
  String get step7Intro =>
      'Construya una supermatriz a partir de los genes alineados y exporte para el análisis filogenético.';

  @override
  String get srcTrimmed => 'Secuencias recortadas';

  @override
  String get recommendedTag => '(recomendado)';

  @override
  String get srcTrimmedBody =>
      'Usa la salida de trimAl — elimina posiciones mal alineadas.';

  @override
  String get srcAligned => 'Secuencias alineadas';

  @override
  String get srcAlignedBody => 'Usa la salida de MAFFT sin recorte.';

  @override
  String get partitionScheme => 'Particiones';

  @override
  String get schemeNone => 'Ninguna';

  @override
  String get schemeByLocus => 'Por locus';

  @override
  String get schemeByLocusCodon => 'Por locus y posición del codón';

  @override
  String get schemeByCodon => 'Por posición del codón (todos los loci)';

  @override
  String get codonExclusion => 'Posiciones del codón';

  @override
  String get codonExclusionInfo =>
      'Las posiciones excluidas se eliminan solo de los loci codificantes. La prueba de saturación (etapa 6) puede completarlo.';

  @override
  String excludePosition(String n) {
    return 'Excluir posición $n';
  }

  @override
  String perLocusExclusions(String list) {
    return 'Por locus: $list';
  }

  @override
  String get clearExclusions => 'Borrar exclusiones por locus';

  @override
  String get buildMatrix => 'Generar supermatriz';

  @override
  String get matrixReady => 'Supermatriz lista';

  @override
  String matrixStats(
    String taxa,
    String sites,
    String partitions,
    String missing,
  ) {
    return '$taxa taxones × $sites sitios · $partitions particiones · $missing% faltantes';
  }

  @override
  String get excludedLabel => 'Excluido';

  @override
  String get outputFiles => 'Archivos generados';

  @override
  String get openFolder => 'Abrir carpeta';

  @override
  String get partitionPreview => 'Particiones RAxML-NG';

  @override
  String get step8Intro =>
      'Estimación de árbol por máxima verosimilitud a partir de la supermatriz concatenada con IQ-TREE3 o RAxML-NG.';

  @override
  String get outgroup => 'Grupo externo';

  @override
  String get outgroupInfo =>
      'Las especies seleccionadas se escriben en el NEXUS y se usan para enraizar los árboles de IQ-TREE y RAxML-NG.';

  @override
  String get outgroupHint => 'Buscar especies…';

  @override
  String engineN(String n) {
    return 'Programa $n';
  }

  @override
  String get iqtreeInfo =>
      'ModelFinder elige los modelos (y puede unir particiones); el bootstrap ultrarrápido y SH-aLRT dan los soportes.';

  @override
  String get raxmlInfo =>
      'RAxML-NG 2 elige los modelos por sí mismo (MOOSE), ajusta la búsqueda a la dificultad del alineamiento y ofrece medidas de soporte rápidas. Como IQ-TREE, es un programa de máxima verosimilitud (no bayesiano).';

  @override
  String get model => 'Modelo de sustitución';

  @override
  String get ufboot => 'Bootstrap ultrarrápido (-B)';

  @override
  String get alrt => 'SH-aLRT (--alrt)';

  @override
  String get runIqtree => 'Ejecutar IQ-TREE3';

  @override
  String get raxmlModelFinder =>
      'ModelFinder (IQ-TREE, modelos compatibles con RAxML)';

  @override
  String get bootstraps => 'Árboles de bootstrap';

  @override
  String get supportMetrics => 'Soporte';

  @override
  String get seed => 'Semilla';

  @override
  String get runRaxml => 'Ejecutar RAxML-NG';

  @override
  String get runs => 'Ejecuciones';

  @override
  String get noRuns => 'Aún no hay árboles.';

  @override
  String toolMissing(String tool) {
    return '$tool no encontrado';
  }

  @override
  String get toolsTitle => 'Herramientas externas y configuración de análisis';

  @override
  String get toolsBody =>
      'SPLACE busca los binarios incluidos, luego el PATH y luego WSL. Defina una ruta si es necesario.';

  @override
  String get tools => 'Herramientas';

  @override
  String get redetect => 'Detectar de nuevo';

  @override
  String get notFound => 'No encontrado';

  @override
  String get overridePath => 'Ruta personalizada (o wsl:/ruta/en/wsl)';

  @override
  String get browse => 'Buscar…';

  @override
  String get analysisSettings => 'Análisis';

  @override
  String threads(String n) {
    return 'Hilos: $n';
  }

  @override
  String get geneticCode => 'Código genético';

  @override
  String get projectFolder => 'Carpeta del proyecto';

  @override
  String get change => 'Cambiar…';

  @override
  String get tipComplete =>
      'Prefiera **genomas completos** siempre que sea posible, porque suelen mejorar los análisis posteriores.';

  @override
  String get tipMultipleSearches =>
      'Realice **varias búsquedas taxonómicas** para agregar **grupos internos** y **externos**. Los nuevos resultados se añaden a la misma tabla taxonómica.';

  @override
  String get tipCombine =>
      'Puede **combinar** búsquedas taxonómicas, importación por accesos y archivos locales **.gb, .gbk o .genbank** en el mismo conjunto de datos.';

  @override
  String get example => 'Ejemplo';

  @override
  String get taxTitle => 'Buscar por nombre taxonómico';

  @override
  String get taxInfo =>
      'Busque en NCBI por nombre taxonómico, tipo de genoma, base de datos y completitud antes de importar los genomas elegidos.';

  @override
  String get taxBody =>
      'Busque en NCBI y descargue automáticamente los registros GenBank coincidentes en la herramienta.';

  @override
  String get taxExampleTooltip =>
      'Rellenar la búsqueda taxonómica con Carangidae, mitocondria, genomas completos y solo RefSeq.';

  @override
  String get taxApiHint =>
      'Agregue una clave de API de NCBI para acelerar las importaciones (botón de la llave en la barra superior).';

  @override
  String get taxName => 'Nombre taxonómico';

  @override
  String get taxOrganelle => 'Tipo de genoma';

  @override
  String get organelleMito => 'Mitocondria';

  @override
  String get organelleChloro => 'Cloroplasto';

  @override
  String get taxScope => 'Alcance del genoma';

  @override
  String get taxCompleteOnly => 'Solo completos';

  @override
  String get taxScopeInfo =>
      'Solo completos busca en NCBI los registros con “complete genome” en el título. Desactivado, también incluye los registros con “partial genome” en el título.';

  @override
  String get scopeComplete => 'Genomas completos';

  @override
  String get scopePartial => 'Genomas parciales';

  @override
  String get taxDatabase => 'Base de datos';

  @override
  String get taxRefseqOnly => 'Solo RefSeq';

  @override
  String get taxQueryPreview => 'Vista previa de la query';

  @override
  String get search => 'Buscar';

  @override
  String get accTitle => 'Pegar números de acceso';

  @override
  String get accInfo =>
      'Pegue uno o más accesos conocidos y recupere exactamente esos registros GenBank desde NCBI.';

  @override
  String get accBody =>
      'Introduzca los números de acceso separados por comas, espacios o saltos de línea.';

  @override
  String get accExampleTooltip =>
      'Cargar accesos de ejemplo de Balcázar et al. (2026)';

  @override
  String get accFetchTooltip =>
      'Buscar registros GenBank en NCBI para los accesos introducidos';

  @override
  String get fetch => 'Buscar';

  @override
  String get gbFilesTitle => 'Usar archivos GenBank locales';

  @override
  String get gbFilesInfo =>
      'Importe archivos locales .gb, .gbk o .genbank cuando los registros ya estén en su ordenador.';

  @override
  String get gbFilesBody =>
      'Importe sus propios archivos .gb, .gbk o .genbank desde el disco.';

  @override
  String get gbFilesDrop =>
      'Arrastre y suelte archivos o carpetas GenBank aquí';

  @override
  String get gbFilesBrowse => 'o haga clic para explorar (.gb, .gbk, .genbank)';

  @override
  String gbSummary(String n) {
    return '$n registros GenBank cargados';
  }

  @override
  String get clearFasta => 'Limpiar archivos FASTA';

  @override
  String get close => 'Cerrar';

  @override
  String get apiKeyTitle => 'API key de NCBI';

  @override
  String get apiKeyHelper =>
      'Se usa en la Opción 1 y la Opción 2 para las solicitudes a NCBI. Sin API key, SPLACE funciona a unas 3 solicitudes por segundo. Con una clave válida guardada, SPLACE usa hasta 10 solicitudes por segundo.';

  @override
  String get apiKeySave => 'Guardar clave';

  @override
  String get apiKeyClear => 'Borrar clave';

  @override
  String get apiKeyGetLink =>
      'Obtener una API key en la configuración de la cuenta de NCBI';

  @override
  String get apiKeySaved =>
      'API key guardada. La Opción 1 y la Opción 2 ahora usan hasta 10 solicitudes por segundo.';

  @override
  String get apiKeyInvalid =>
      'La API key parece inválida. Revísela antes de guardarla.';

  @override
  String get apiKeyNone =>
      'No hay API key guardada. La Opción 1 y la Opción 2 están usando hasta 3 solicitudes por segundo.';

  @override
  String get apiClearTitle => '¿Eliminar la API key guardada?';

  @override
  String get apiClearBody =>
      'La Opción 1 y la Opción 2 volverán a la tasa predeterminada de solicitudes a NCBI.';

  @override
  String get reviewTitle => 'Revisar genomas antes de importar';

  @override
  String get reviewSubtitle =>
      'Elija qué genomas deben importarse. Todos los interruptores empiezan activados.';

  @override
  String get reviewDetailTaxon => 'Consulta taxonómica';

  @override
  String get reviewDetailGenomeType => 'Tipo de genoma';

  @override
  String get reviewDetailSearchTypes => 'Tipos de búsqueda';

  @override
  String get reviewDetailRefseq => 'RefSeq';

  @override
  String get reviewRefseqOn => 'Usado';

  @override
  String get reviewRefseqOff => 'No usado';

  @override
  String get reviewSelectAll => 'Activar todos';

  @override
  String get reviewSelectNone => 'Desactivar todos';

  @override
  String get reviewDisableUnverified => 'Desactivar UNVERIFIED';

  @override
  String get reviewOnlyComplete => 'Solo completos';

  @override
  String get reviewOnlyPartial => 'Solo parciales';

  @override
  String get reviewColImport => 'Importar';

  @override
  String get reviewColAccession => 'Acceso';

  @override
  String get reviewColDescription => 'Descripción del genoma';

  @override
  String get reviewColScope => 'Clasificación';

  @override
  String get reviewImport => 'Importar genomas seleccionados';

  @override
  String get reviewTotal => 'Total';

  @override
  String get reviewSelected => 'Seleccionados';

  @override
  String get reviewComplete => 'Completo';

  @override
  String get reviewPartial => 'Parcial';

  @override
  String get reviewMixed => 'Título mixto';

  @override
  String get reviewUnknown => 'Sin clasificar';

  @override
  String get reviewUnverified => 'UNVERIFIED';

  @override
  String reviewUnverifiedCount(String count) {
    return 'Genomas UNVERIFIED encontrados: $count';
  }

  @override
  String get reviewUnverifiedAdvice =>
      'Se recomienda eliminar los genomas UNVERIFIED antes de importar.';

  @override
  String get reviewCancelled => 'Selección de genomas cancelada.';

  @override
  String get reviewSearch => 'Filtrar por accesión, especie o descripción';

  @override
  String reviewShowing(String shown, String total) {
    return 'Mostrando $shown/$total genomas que coinciden con la búsqueda.';
  }

  @override
  String get gbRecordsHeading => 'Registros GenBank';

  @override
  String get fastaHeadersHeading => 'Encabezados FASTA';

  @override
  String get clearAll => 'Limpiar todo';

  @override
  String get clearTitle => '¿Limpiar todos los registros?';

  @override
  String get clearBody =>
      'Esta acción no se puede deshacer. Se eliminarán todos los registros cargados y las selecciones.';

  @override
  String get recordsReview =>
      'Revise y corrija los nombres de las especies usando el botón de edición.';

  @override
  String get recordsReviewMethod =>
      'Revise y corrija los nombres de las especies abajo con el botón de edición y luego obtenga la taxonomía con el método de [[dataFishing]].';

  @override
  String get recordsFetchPrefix => 'Luego haga clic en';

  @override
  String get recordsFetchSuffix =>
      'para obtener los rangos taxonómicos desde GBIF mediante';

  @override
  String get internetRequired => 'Internet requerida';

  @override
  String get fetchTaxonomy => 'Obtener taxonomía';

  @override
  String get skipTaxonomy => 'Omitir taxonomía y continuar';

  @override
  String get skipTaxBody1 =>
      'Si los nombres de especies contienen caracteres especiales o hay entradas duplicadas, pueden producirse errores en análisis posteriores.';

  @override
  String get skipTaxBody2 =>
      'Se recomienda corregir los nombres y recuperar la taxonomía para obtener mejores resultados.';

  @override
  String get taxonomyResolved => 'Taxonomía lista';

  @override
  String get taxonomyTree => 'Árbol taxonómico';

  @override
  String get taxonomyTreeBody =>
      'Resumen taxonómico indentado con los mismos registros que se muestran abajo.';

  @override
  String get columnsLabel => 'Columnas:';

  @override
  String get colAccession => 'Acceso';

  @override
  String get colKingdom => 'Reino';

  @override
  String get colPhylum => 'Filo';

  @override
  String get colClass => 'Clase';

  @override
  String get colOrder => 'Orden';

  @override
  String get colFamily => 'Familia';

  @override
  String get colGenus => 'Género';

  @override
  String get colSpecies => 'Especie';

  @override
  String get colAuthorship => 'Autoría';

  @override
  String get colSource => 'Origen';

  @override
  String get sourceFile => 'Archivo';

  @override
  String get recordsSearch => 'Buscar acceso, especie, familia...';

  @override
  String editRecordTooltip(String species, String accession) {
    return 'Editar $species ($accession)';
  }

  @override
  String removeRecordTooltip(String species, String accession) {
    return 'Eliminar $species ($accession)';
  }

  @override
  String get removeTitle => '¿Eliminar registro?';

  @override
  String get removeBody => 'Esta acción no se puede deshacer.';

  @override
  String get remove => 'Eliminar';

  @override
  String get editTitle => 'Editar registro';

  @override
  String get editHint =>
      'Revise el registro actual y actualice abajo los campos de especie y taxonomía.';

  @override
  String get editGbif => 'Obtener taxonomía desde GBIF';

  @override
  String get save => 'Guardar';

  @override
  String get autoLabel => '(automático)';

  @override
  String get fileSize => 'Tamaño';

  @override
  String get sequenceLabel => 'Secuencia';

  @override
  String get featuresLabel => 'Features';

  @override
  String get lineageLabel => 'Linaje';

  @override
  String get speciesRequired => 'El nombre de la especie es obligatorio';

  @override
  String get speciesLetters => 'Solo se permiten letras, espacios y guiones';

  @override
  String get speciesBinomial => 'Use género y epíteto (p. ej. Rhinella marina)';

  @override
  String get speciesGenusUpper => 'El género debe empezar con mayúscula';

  @override
  String get familySuffix =>
      'Los nombres de familia terminan en -idae o -aceae';

  @override
  String get gbifNoMatch => 'Sin coincidencias en GBIF para este nombre.';

  @override
  String get gbifMatched => 'Coincidencia encontrada en GBIF.';

  @override
  String get step3Fasta => 'Campos de muestra y locus';

  @override
  String get featureTypesIntro =>
      'Seleccione las categorías de genes que desea incluir en su análisis. Active o desactive cada categoría.';

  @override
  String get standardizedPrefix =>
      'Los nombres de los genes se han estandarizado usando';

  @override
  String get featureCdsInfo =>
      'Secuencias codificantes de proteínas (genes traducidos a proteínas)';

  @override
  String get featureRrnaInfo => 'Genes de ARN ribosómico (12S, 16S, etc.)';

  @override
  String get featureTrnaInfo => 'Genes de ARN de transferencia';

  @override
  String get featureTypesHint =>
      'Active al menos una categoría para listar los marcadores en la etapa 4.';

  @override
  String unknownGenesButton(String n) {
    return '$n nombres desconocidos';
  }

  @override
  String get unknownTitle =>
      'Genes No Encontrados en SynGenes (Rabelo et al. 2024)';

  @override
  String get unknownSubtitle =>
      'Seleccione el nombre estandarizado correcto para cada gen encontrado en sus archivos.';

  @override
  String unknownCount(String count) {
    return '$count gen(es) no pudieron ser estandarizados automáticamente.';
  }

  @override
  String get unknownHowto =>
      'Seleccione el nombre estandarizado en el menú desplegable de la columna Nombre Correcto. Elija (omitir) para excluir el gen del análisis.';

  @override
  String get unknownSkip => '(omitir)';

  @override
  String geneSuggested(String name) {
    return 'Sugerencia: $name';
  }

  @override
  String geneFromProduct(String value) {
    return 'por el /product “$value”';
  }

  @override
  String geneFromGene(String value) {
    return 'por el /gene “$value”';
  }

  @override
  String geneFromSequence(String percent, String name, String n) {
    return 'secuencia $percent% similar a $name en $n otros genomas';
  }

  @override
  String geneSequenceAgrees(String percent) {
    return 'la secuencia lo confirma ($percent%)';
  }

  @override
  String geneSequenceDisagrees(String name, String percent) {
    return 'pero la secuencia parece $name ($percent%)';
  }

  @override
  String get geneNoSuggestion =>
      'Sin coincidencia en SynGenes ni en los otros genomas; revise la anotación.';

  @override
  String geneQualifierMismatch(String gene, String product) {
    return '/gene “$gene” y /product “$product” son diferentes';
  }

  @override
  String conflictsBanner(String n) {
    return '$n feature(s) tienen /gene y /product que indican genes diferentes.';
  }

  @override
  String get conflictsReview => 'Revisar';

  @override
  String get conflictsTitle => '/gene y /product no coinciden';

  @override
  String get conflictsSubtitle =>
      'SPLACE usa el nombre del /product. Cada secuencia se compara con el mismo marcador en los otros genomas para mostrar qué calificador es correcto.';

  @override
  String conflictUsed(String name) {
    return 'Usado: $name';
  }

  @override
  String conflictSeqGene(String name, String percent) {
    return 'La secuencia corresponde a $name ($percent%): el /gene parece correcto';
  }

  @override
  String conflictSeqUsed(String name, String percent) {
    return 'La secuencia corresponde a $name ($percent%): el nombre usado parece correcto';
  }

  @override
  String conflictSeqOther(String name, String percent) {
    return 'La secuencia corresponde a $name ($percent%), ninguno de los dos';
  }

  @override
  String get conflictSeqNone => 'Secuencia no concluyente';

  @override
  String conflictUseGene(String name) {
    return 'Usar $name';
  }

  @override
  String get speciesIdRunning => 'Identificando especies';

  @override
  String speciesIdGenesUsed(String genes) {
    return 'BLAST en GenBank y BOLD · $genes';
  }

  @override
  String get speciesIdFooter =>
      'El BLAST en NCBI suele tardar de 1 a 3 minutos; BOLD responde en segundos.';

  @override
  String get speciesIdButton => 'Identificar especies (GenBank + BOLD)';

  @override
  String get speciesIdShowResults => 'Resultado de la identificación';

  @override
  String get speciesIdSetupTitle => 'Identificar especies con GenBank y BOLD';

  @override
  String get speciesIdSetupSubtitle =>
      'Compare los genes de cada genoma con GenBank (BLAST) y BOLD para verificar los nombres de las especies antes de obtener la taxonomía.';

  @override
  String get speciesIdGenesLabel => 'Genes para comparar';

  @override
  String get speciesIdGenesHint =>
      'Elija uno o más genes. Cada gen se busca en GenBank con BLAST; BOLD acepta solo COI.';

  @override
  String speciesIdInGenomes(String n) {
    return 'en $n genomas';
  }

  @override
  String get speciesIdSources => 'Bases de datos';

  @override
  String get speciesIdGenbankBody =>
      'BLAST (megablast) contra core_nt, para cada gen elegido.';

  @override
  String get speciesIdBoldBody =>
      'ID engine de BOLD, biblioteca pública de especies, solo COI.';

  @override
  String get speciesIdBoldOff => 'Elija COI para buscar también en BOLD.';

  @override
  String get speciesIdRun => 'Identificar';

  @override
  String get speciesIdResultsTitle => 'Identificación de las especies';

  @override
  String speciesIdResultsSubtitle(String genes) {
    return 'BLAST en GenBank y BOLD · $genes. Elija el nombre que se conserva para cada genoma.';
  }

  @override
  String get speciesIdOnlyDisagree => 'Solo posibles identificaciones erróneas';

  @override
  String speciesIdSummary(String disagree, String total) {
    return '$disagree de $total genomas pueden estar mal identificados';
  }

  @override
  String get speciesIdColGenome => 'Genoma';

  @override
  String get speciesIdColCandidates => 'Especies encontradas (elija una)';

  @override
  String get speciesIdMismatch => 'Posible identificación errónea';

  @override
  String get speciesIdMatch => 'Confirmada';

  @override
  String get speciesIdNoHits => 'Sin coincidencias en GenBank ni en BOLD';

  @override
  String speciesIdKeepCurrent(String name) {
    return 'Mantener $name';
  }

  @override
  String speciesIdHits(String n) {
    return '$n hits';
  }

  @override
  String speciesIdApply(String n) {
    return 'Aplicar $n cambio(s)';
  }

  @override
  String get speciesIdNote =>
      'Los hits del propio genoma (o del registro de GenBank del que proviene la copia RefSeq NC_) aparecen como “este registro” y no cuentan como evidencia.';

  @override
  String get geneReviewButton => 'Revisar nombres de genes';

  @override
  String geneReviewButtonN(String n) {
    return 'Revisar nombres de genes ($n)';
  }

  @override
  String get geneReviewDone => 'Nombres de genes revisados';

  @override
  String get geneReviewTitle => 'Revisión de los nombres de genes';

  @override
  String get geneReviewSubtitle =>
      'Nombres que SynGenes no reconoció y genes cuyo /gene y /product no coinciden. Revíselos antes de la taxonomía.';

  @override
  String geneReviewUnknownSection(String n) {
    return 'Nombres no encontrados en SynGenes ($n)';
  }

  @override
  String geneReviewConflictSection(String n) {
    return '/gene y /product indican genes diferentes ($n)';
  }

  @override
  String get geneReviewConflictHelp =>
      'Ambos calificadores se normalizaron con SynGenes. Elija qué nombre recibe esta feature; solo cambia este genoma.';

  @override
  String get geneReviewNothing =>
      'Todos los nombres de genes se estandarizaron y cada /gene coincide con su /product.';

  @override
  String conflictPickGene(String name) {
    return '/gene → $name';
  }

  @override
  String conflictPickProduct(String name) {
    return '/product → $name';
  }

  @override
  String get grHowTitle => 'Cómo leer esta revisión';

  @override
  String get grHowSynGenes =>
      'SynGenes reconoce una etiqueta cuando es un nombre estándar de gen o uno de sus sinónimos conocidos.';

  @override
  String get grHowSequence =>
      'La similitud es la fracción de la secuencia (proteína para CDS, ADN para rRNA) que también aparece en el mismo gen de los otros genomas de este conjunto: más de 50% es fuerte, de 35 a 50% moderada y menos de 35% débil.';

  @override
  String get grSequence => 'Secuencia';

  @override
  String grRecognized(String name) {
    return 'SynGenes: $name';
  }

  @override
  String get grNotRecognized => 'no está en SynGenes';

  @override
  String get grAbsent => 'ausente';

  @override
  String grInGenomes(String n) {
    return 'en $n genoma(s)';
  }

  @override
  String get grUseName => 'Nombre a usar';

  @override
  String get grNoSequenceRef => 'ningún otro genoma para comparar';

  @override
  String grWhyProduct(String value, String name) {
    return 'El /product “$value” es un nombre de $name en SynGenes.';
  }

  @override
  String grWhyGene(String value, String name) {
    return 'El /gene “$value” es un nombre de $name en SynGenes.';
  }

  @override
  String grWhySequence(String percent, String name, String n) {
    return 'Ninguna etiqueta está en SynGenes, pero la secuencia es $percent% similar a $name en $n otros genomas.';
  }

  @override
  String get grWhyNone =>
      'Ninguna etiqueta está en SynGenes y la secuencia no corresponde claramente a ningún marcador. Revise la anotación u omita este nombre.';

  @override
  String grSeqAgrees(String percent) {
    return 'La secuencia lo confirma ($percent% de similitud).';
  }

  @override
  String grSeqDisagrees(String name, String percent) {
    return 'Atención: la secuencia se parece más a $name ($percent%).';
  }

  @override
  String grSays(String tag, String name) {
    return '$tag dice $name';
  }

  @override
  String grSimilarityTo(String name) {
    return 'Similitud con $name en los otros genomas';
  }

  @override
  String get grSupported => 'Confirmado por la secuencia';

  @override
  String get grNotSupported => 'No confirmado por la secuencia';

  @override
  String grConflictGeneRight(
    String gp,
    String gene,
    String pp,
    String product,
  ) {
    return 'La secuencia es $gp% similar a $gene y $pp% a $product en los otros genomas; por lo tanto, la etiqueta /gene ($gene) probablemente es correcta y la etiqueta /product está equivocada.';
  }

  @override
  String grConflictProductRight(
    String pp,
    String product,
    String gp,
    String gene,
  ) {
    return 'La secuencia es $pp% similar a $product y $gp% a $gene en los otros genomas; por lo tanto, la etiqueta /product ($product) probablemente es correcta y la etiqueta /gene está equivocada.';
  }

  @override
  String get grConflictUnclear =>
      'La secuencia no favorece claramente a ninguno de los dos genes. Revise la anotación de este genoma antes de elegir.';

  @override
  String get grConflictIntro =>
      'Las etiquetas /gene y /product de esta feature indican genes diferentes. Haga clic en el nombre a usar; solo cambia este genoma.';

  @override
  String get siHowTitle => 'Cómo leer la identificación';

  @override
  String get siHowBody =>
      'Los genes de cada genoma se compararon con GenBank (BLAST) y BOLD. La identidad es la fracción de bases idénticas con la secuencia más parecida encontrada.';

  @override
  String get siLegendSame => '≥ 98%: misma especie';

  @override
  String get siLegendClose => '95–98%: especie cercana';

  @override
  String get siLegendOther => '< 95%: otra especie';

  @override
  String get siColUse => 'Usar';

  @override
  String get siColSpecies => 'Especie';

  @override
  String get siColGenbank => 'GenBank (BLAST)';

  @override
  String get siColBold => 'BOLD (COI)';

  @override
  String get siColSeqs => 'Secuencias encontradas';

  @override
  String get siCurrentTag => 'nombre actual';

  @override
  String get siNotFound => 'no encontrada';

  @override
  String siVerdictConfirmed(String identity) {
    return 'Las coincidencias confirman el nombre actual (mejor identidad $identity%).';
  }

  @override
  String siVerdictMismatch(
    String best,
    String identity,
    String current,
    String currentIdentity,
  ) {
    return 'La coincidencia más cercana es $best ($identity%), mientras que el nombre actual $current llega a $currentIdentity%. El genoma puede estar mal identificado.';
  }

  @override
  String siVerdictMissing(String best, String identity, String current) {
    return 'La coincidencia más cercana es $best ($identity%) y $current no aparece entre las coincidencias. El genoma puede estar mal identificado, o el nombre puede ser un sinónimo.';
  }

  @override
  String get siVerdictNone =>
      'Ninguna secuencia de GenBank o BOLD coincidió con los genes de este genoma.';

  @override
  String noMarkersButton(String n) {
    return 'Eliminar genomas sin marcadores ($n)';
  }

  @override
  String get noMarkersTitle => 'Genomas sin marcadores';

  @override
  String get noMarkersSubtitle =>
      'Estos genomas no tienen ningún gen codificante ni rRNA que SPLACE pueda usar, así que no aportan nada a la matriz. Revise la lista antes de eliminarlos.';

  @override
  String get noMarkersReasonNone =>
      'Ningún gen anotado (0 PCGs, 0 rRNAs, 0 tRNAs): el registro de GenBank no tiene anotación de genes.';

  @override
  String noMarkersReasonTrna(String n) {
    return 'Solo $n tRNA(s): ningún PCG ni rRNA.';
  }

  @override
  String noMarkersReasonUnnamed(String pcgs, String rrnas) {
    return '$pcgs PCG(s) y $rrnas rRNA(s) anotados, pero ninguno tiene nombre en /gene o /product para identificarlo.';
  }

  @override
  String noMarkersRemove(String n) {
    return 'Eliminar $n genoma(s)';
  }

  @override
  String get heatmapShow => 'Mostrar heatmap';

  @override
  String get heatmapHide => 'Ocultar heatmap';

  @override
  String get heatmapHoverHint =>
      'Apunte a una celda, un marcador o una muestra para ver los detalles.';

  @override
  String get headerPresets => 'Plantillas';

  @override
  String get headerDragHint =>
      'Arrastre para reordenar. El ejemplo viene del primer registro.';

  @override
  String get headerAddHint =>
      'Haga clic en un campo para añadirlo al final del encabezado.';

  @override
  String get headerFirstRecord => 'Primer registro';

  @override
  String get headerFirstRecords => 'Encabezados de los primeros registros';

  @override
  String siSourceFailed(String source, String error) {
    return '$source no respondió, por eso su columna no tiene resultados: $error';
  }

  @override
  String get whenToUse => 'Cuándo usar';

  @override
  String get trimCodonHandling => 'Genes codificantes';

  @override
  String get trimCodonHandlingInfo =>
      'Las posiciones del codón se usan después en las particiones y la prueba de saturación (paso 6), así que los genes codificantes deben conservar codones enteros. trimAl lo hace con -backtrans (recorta el alineamiento traducido y devuelve los codones; -ignorestopcodon porque comprueba los stops con el código estándar) y ClipKIT con --codon. Excluir una posición del codón de la supermatriz se decide en el paso 7b, no aquí.';

  @override
  String get alignRunning => 'Alineando con MAFFT';

  @override
  String get trimRunning => 'Recortando los alineamientos';

  @override
  String markersProgress(String done, String total, String mode) {
    return '$done de $total marcadores · $mode';
  }

  @override
  String get alignFooter =>
      'MAFFT se ejecuta en este equipo, varios marcadores a la vez.';

  @override
  String get trimFooter =>
      'trimAl/ClipKIT se ejecutan en este equipo, varios marcadores a la vez.';

  @override
  String alignCurrent(String markers) {
    return 'En curso: $markers';
  }

  @override
  String get colReason => 'Motivo';

  @override
  String get unknownSkipAll => 'Omitir Todos';

  @override
  String get unknownApply => 'Aplicar Correcciones';

  @override
  String get more => 'más';

  @override
  String get defaultPhylogeny => 'Por defecto para filogenia';

  @override
  String get legendDefault => 'conjunto recomendado para análisis filogenético';

  @override
  String get tooltipAll => 'Seleccionar todos los marcadores visibles';

  @override
  String get tooltipNone => 'Limpiar todos los marcadores seleccionados';

  @override
  String get tooltipShared =>
      'Seleccionar solo los genes presentes en todos los genomas (sin datos faltantes)';

  @override
  String get tooltipDefault =>
      'Seleccionar el conjunto predeterminado de genes para análisis filogenético';

  @override
  String get duplicateChoiceHelp =>
      'La supermatriz usa una copia por genoma: haga clic en la copia a usar. La más larga empieza seleccionada; “más cercana” marca la copia cuya longitud se acerca más a los genomas con una sola copia, normalmente la correcta cuando las copias difieren.';

  @override
  String get defaultLabel => 'predeterminado';

  @override
  String dupReference(
    String gene,
    String n,
    String avg,
    String min,
    String max,
  ) {
    return '$gene en $n genomas con una sola copia: promedio de $avg pb (de $min a $max pb). Elija la copia más cercana.';
  }

  @override
  String get dupNoReference =>
      'Ningún otro genoma tiene una sola copia de este gen para comparar las longitudes.';

  @override
  String get dupCopiesColumn => 'Copias (haga clic en la que se usará)';

  @override
  String dupVsAverage(String diff) {
    return '$diff pb respecto al promedio';
  }

  @override
  String get dupSameAsAverage => 'igual al promedio';

  @override
  String get dupClosest => 'más cercana al promedio';

  @override
  String get dupLongestTag => 'más larga';

  @override
  String get noMarkers =>
      'No se encontraron marcadores para los tipos seleccionados.';

  @override
  String get featureLabel => 'Feature';

  @override
  String get recordsLabel => 'Registros';

  @override
  String heatmapSummary(String markers, String samples) {
    return '$markers marcadores × $samples muestras. Haga clic en el nombre de un marcador para seleccionarlo.';
  }

  @override
  String get lengthHeatmapTitle => 'Mapa de longitud de marcadores';

  @override
  String get lengthHeatmapBody =>
      'Cada marcador seleccionado se escala de forma independiente según su longitud mínima y máxima.';

  @override
  String get panelHeader => 'Configurar encabezado FASTA';

  @override
  String get panelHeaderBody =>
      'Identificador de las secuencias en los archivos de alineamiento. Se usa el predeterminado (acceso_familia_género_especie) salvo que lo cambie en el constructor de encabezados.';

  @override
  String get openHeaderBuilder => 'Abrir constructor de encabezado';

  @override
  String get headerBuilderTitle => 'Constructor de encabezado FASTA';

  @override
  String get headerBuilderSubtitle =>
      'Elija, ordene y separe los campos de cada nombre de secuencia.';

  @override
  String get headerFields => 'Campos del encabezado (en orden)';

  @override
  String get headerAvailable => 'Campos disponibles';

  @override
  String get headerSeparator => 'Separador';

  @override
  String get headerEmpty => 'Agregue al menos un campo.';

  @override
  String get headerConfirm => 'Confirmar encabezado';

  @override
  String get headerConfirmed => 'Encabezado confirmado';

  @override
  String get headerUnique => 'Todos los encabezados son únicos.';

  @override
  String headerDuplicates(String n) {
    return '$n encabezado(s) se repiten en más de un registro; se agregarán sufijos (_2, _3…). Incluya la Accesión para hacerlos únicos.';
  }

  @override
  String get moveLeft => 'Mover antes';

  @override
  String get moveRight => 'Mover después';

  @override
  String get panelSeqType => 'Tipo de secuencia';

  @override
  String get seqTypeBody =>
      'Elija cómo alinear las secuencias. El modo Codon-aware alinea las traducciones proteicas con MAFFT y luego retrotraduce a DNA, preservando el marco de lectura.';

  @override
  String get seqNtTitle => 'Nucleótido';

  @override
  String get seqNtBody => 'Alinear secuencias de DNA/RNA directamente';

  @override
  String get seqAaTitle => 'Aminoácido';

  @override
  String get seqAaBody =>
      'Traduce (codon_start, transl_table) y alinea las proteínas';

  @override
  String get seqCodonTitle => 'Codon-aware';

  @override
  String get seqCodonBody =>
      'Alineamiento de codones guiado por proteína (solo CDS)';

  @override
  String get recNtTitle => 'Nucleótido, alineamiento directo del ADN';

  @override
  String get recNtBody =>
      'Mejor para rRNA, tRNA, regiones no codificantes y conjuntos exploratorios. El FASTA va directo a MAFFT.';

  @override
  String get recAaTitle => 'Aminoácido, alineamiento de proteínas';

  @override
  String get recAaBody =>
      'Mejor cuando los genes codificantes son divergentes. El ADN se traduce, MAFFT alinea las proteínas y el resultado se analiza como proteína.';

  @override
  String get recCodonTitle => 'Codón-aware, filogenia de CDS';

  @override
  String get recCodonBody =>
      'Recomendado para matrices de CDS con salida en ADN y codones preservados: valida la longitud de las CDS, elimina stops terminales, alinea proteínas y retrotraduce.';

  @override
  String get geneticCodeInfo =>
      'Se usa en la traducción, la detección del marco de lectura y la prueba de saturación. El calificador /transl_table de GenBank tiene prioridad.';

  @override
  String get translationInfo =>
      'Las secuencias se traducen antes del alineamiento. MAFFT se ejecuta en modo --amino.';

  @override
  String get codonWarning =>
      'El modo Codon-aware solo procesa features CDS. Los genes sin CDS anotado (rRNA, tRNA, etc.) se omitirán. El alineamiento resultante preserva el marco de lectura y admite partición por posición de codón.';

  @override
  String get codonWorkflow => 'Flujo codón-aware';

  @override
  String get codonStep1 => 'Extraer CDS';

  @override
  String get codonStep1Body =>
      'Solo se usan CDS anotadas; se omiten rRNA y tRNA.';

  @override
  String get codonStep2 => 'Validar marco';

  @override
  String get codonStep2Body =>
      'Longitud divisible por 3 tras codon_start y el tratamiento del stop terminal.';

  @override
  String get codonStep3 => 'Traducir';

  @override
  String get codonStep3Body => 'Con el código genético NCBI de cada CDS.';

  @override
  String get codonStep4 => 'MAFFT proteína';

  @override
  String get codonStep4Body => 'Proteínas alineadas con mafft --amino.';

  @override
  String get codonStep5 => 'Retrotraducir';

  @override
  String get codonStep5Body =>
      'El alineamiento de proteínas guía un alineamiento de ADN que preserva los codones.';

  @override
  String get allowPseudo => 'Permitir CDS con /pseudo y /pseudogene';

  @override
  String get allowPseudoBody =>
      'Desactivado por defecto. Cuando está apagado, los pseudo CDS se excluyen antes de la traducción.';

  @override
  String get allowStops => 'Permitir stop codons internos';

  @override
  String get allowStopsBody =>
      'Desactivado por defecto. Cuando está activado, las secuencias con stops internos permanecen en la corrida y se informan con advertencias.';

  @override
  String get allowBtWarnings =>
      'Continuar cuando las verificaciones de aminoácido y codón generen avisos';

  @override
  String get allowBtWarningsBody =>
      'Desactivado por defecto. Apagado, cualquier discrepancia aminoácido/codón tras MAFFT falla el marcador; encendido, los avisos quedan en el registro.';

  @override
  String get pythonTitle => 'ClipKIT (recorte codon-aware)';

  @override
  String get pythonBody =>
      'ClipKIT funciona con Python 3. SPLACE lo usa en los alineamientos codon-aware: recorta el alineamiento de proteínas y SPLACE retrotraduce el resultado a codones.';

  @override
  String get pythonNotChecked => 'Python aún no fue verificado.';

  @override
  String get pythonMissing => 'No se encontró Python 3 en el PATH.';

  @override
  String pythonFound(String exe, String version) {
    return '$exe (Python $version)';
  }

  @override
  String clipkitFound(String version) {
    return 'ClipKIT $version disponible';
  }

  @override
  String get clipkitMissing => 'ClipKIT no instalado';

  @override
  String get checkPython => 'Verificar Python';

  @override
  String get installClipkit => 'Instalar ClipKIT';

  @override
  String get mafftAuto => '--auto (automático, recomendado)';

  @override
  String get mafftLocalpair => '--localpair (L-INS-i, más preciso, lento)';

  @override
  String get mafftGlobalpair => '--globalpair (G-INS-i, global, lento)';

  @override
  String get mafftGenafpair => '--genafpair (E-INS-i, gaps largos, lento)';

  @override
  String get mafftRetree2 => '--retree 2 (FFT-NS-2, rápido)';

  @override
  String get mafftRetree1 => '--retree 1 (FFT-NS-1, más rápido)';

  @override
  String get mafftMaxIterate => 'Máximo de iteraciones';

  @override
  String get mafftMaxIterateHelp => '0 = sin refinamiento; 1000 = máximo';

  @override
  String get mafftThreads => 'Hilos';

  @override
  String get mafftThreadsHelp => '(por job)';

  @override
  String mafftConcurrency(String cpus, String jobs, String threads) {
    return '$cpus CPU lógicas: $jobs trabajo(s) × $threads hilo(s)';
  }

  @override
  String get mafftGapOpening => 'Apertura de gap';

  @override
  String get mafftGapExtension => 'Extensión de gap';

  @override
  String get mafftExtra => 'Parámetros adicionales';

  @override
  String get adjustDirectionHelp => 'Invierte secuencias de la hebra opuesta';

  @override
  String get reorderHelp => 'Ordena la salida por similitud';

  @override
  String get optional => '(opcional)';

  @override
  String get optionalLabel => 'Opcional';

  @override
  String get mafftAminoNote =>
      'Los modos aminoácido y codón-aware agregan --amino automáticamente (--adjustdirection no se aplica a proteínas).';

  @override
  String runMafftMode(String mode) {
    return 'Ejecutar alineamiento MAFFT ($mode)';
  }

  @override
  String markersToAlign(String n) {
    return '$n marcador(es) seleccionado(s)';
  }

  @override
  String get trimPromptTitle =>
      'Alineamiento completo - ¿recortar regiones mal alineadas?';

  @override
  String get trimPromptBody =>
      'trimAl elimina columnas con gaps excesivos o baja conservación, mejorando la calidad de la inferencia filogenética.';

  @override
  String get trimPromptCodon =>
      'En modo Codon-aware, ClipKIT recorta el alineamiento de proteínas y SPLACE lo retrotraduce a ADN preservando los codones.';

  @override
  String get trimPromptNo => 'No, omitir recorte';

  @override
  String get trimPromptYes => 'Sí, configurar el recorte';

  @override
  String get trimAutomated1 => '-automated1 (heurístico, recomendado)';

  @override
  String get trimStrict => '-strict (automático estricto)';

  @override
  String get trimStrictplus => '-strictplus (estricto+)';

  @override
  String get trimGappyout => '-gappyout (elimina columnas con muchos gaps)';

  @override
  String get trimNogaps => '-nogaps (elimina todas las columnas con gaps)';

  @override
  String get trimNoallgaps => '-noallgaps (solo columnas totalmente vacías)';

  @override
  String get trimManual => 'Umbrales manuales';

  @override
  String get trimGapThreshold => 'Umbral de gaps';

  @override
  String get trimSimilarity => 'Similitud';

  @override
  String get trimConsistency => 'Consistencia';

  @override
  String get trimWindow => 'Tamaño de ventana';

  @override
  String get trimExtra => 'Parámetros adicionales';

  @override
  String get trimCodonNote =>
      'En las ejecuciones codón-aware MAFFT alineó proteínas: el recortador recorta el alineamiento de proteínas, SPLACE lo retrotraduce a ADN preservando codones y valida el marco.';

  @override
  String get codonTrimBackend => 'Herramienta de recorte de codones';

  @override
  String get codonTrimAuto => 'Auto: trimAl y, si hace falta, ClipKIT';

  @override
  String get codonTrimTrimal =>
      'trimAl, alineamiento de proteínas con retrotraducción';

  @override
  String get codonTrimClipkit => 'ClipKIT (requiere Python)';

  @override
  String get clipkitMode => 'Modo de ClipKIT';

  @override
  String get codonTrimmingTitle => 'Recorte de codones';

  @override
  String get codonTrimmingBody =>
      'Elimina bloques de codones mal alineados antes de la concatenación y de cualquier filtrado por posición del codón.';

  @override
  String get trimSkippedTag => 'Recorte omitido';

  @override
  String get trimming => 'Recorte';

  @override
  String get resultsTitle => 'Resultados del alineamiento';

  @override
  String get noTrimming => 'Sin recorte';

  @override
  String get viewAlignment => 'Mostrar alineamiento';

  @override
  String get tabProtein => 'Proteína';

  @override
  String get viewDna => 'Vista DNA';

  @override
  String get viewTranslate => 'Traducir a aminoácidos';

  @override
  String possibleStops(String n) {
    return '$n posible(s) stop(s) interno(s)';
  }

  @override
  String get conservation => 'CONSERVACIÓN';

  @override
  String get stopLegend => 'Stop';

  @override
  String get gapLegend => 'Gap';

  @override
  String get replaceGaps => 'Reemplazar gaps del alineamiento por ?';

  @override
  String get replaceGapsBody =>
      'Úselo cuando el programa siguiente deba tratar los gaps como datos faltantes. Las coordenadas de las particiones se conservan.';

  @override
  String get selectedLabel => 'Seleccionadas';

  @override
  String get saveNexus => 'Guardar NEXUS';

  @override
  String get savePartition => 'Guardar tabla de particiones';

  @override
  String get saveFasta => 'Guardar FASTA';

  @override
  String get savePhylip => 'Guardar PHYLIP';

  @override
  String get partitionName => 'Partición';

  @override
  String get partitionRange => 'Rango';

  @override
  String get partitionSites => 'Sitios';

  @override
  String get partitionType => 'Tipo';

  @override
  String get noneSelected => 'ninguno seleccionado';

  @override
  String get outgroupFromStep7 => 'elegido en la etapa 7';

  @override
  String get substitutionModel => 'Modelo de sustitución';

  @override
  String get modelTest => 'ModelFinder: mejor modelo automático (TEST)';

  @override
  String get modelMfp => 'MFP: ModelFinder Plus';

  @override
  String get modelMfpMerge =>
      'MFP+MERGE: búsqueda de modelo por partición (recomendado)';

  @override
  String get modelCustom => 'Personalizado…';

  @override
  String get bootstrapSupport => 'Soporte bootstrap';

  @override
  String get bootUf1000 => 'UFBoot2: 1000 réplicas (recomendado)';

  @override
  String get bootUf5000 => 'UFBoot2: 5000 réplicas';

  @override
  String get bootClassic => 'Bootstrap clásico: 100 réplicas (lento)';

  @override
  String get bootNone => 'Sin bootstrap';

  @override
  String get threadsLabel => 'Hilos';

  @override
  String get extraFlags => 'Parámetros extra';

  @override
  String get spaceSeparated => '(separados por espacio)';

  @override
  String get perGeneTrees => 'Inferir también árboles por gen';

  @override
  String get perGeneTreesBody =>
      'Ejecuta un análisis IQ-TREE separado para cada marcador (guardado en gene_trees/).';

  @override
  String get commandPreview => 'Vista previa del comando';

  @override
  String get saveZip => 'Guardar resultados en ZIP';

  @override
  String get analysisComplete => 'Análisis completado';

  @override
  String get outputDirectory => 'Carpeta de salida';

  @override
  String get geneTreesLabel => 'Árboles por gen';

  @override
  String get howToCite => 'Cómo citar';

  @override
  String get howToCiteHint =>
      'Si usa estos resultados, cite SPLACE y las herramientas que ejecutó';

  @override
  String get copy => 'Copiar';

  @override
  String metaGenomes(String n) {
    return '$n genomas';
  }

  @override
  String metaFastaFiles(String n) {
    return '$n archivos FASTA';
  }

  @override
  String metaSamples(String n) {
    return '$n muestras';
  }

  @override
  String metaSelected(String selected, String total) {
    return '$selected de $total seleccionados';
  }

  @override
  String metaAligned(String done, String total) {
    return '$done/$total alineados';
  }

  @override
  String metaTrimmed(String n) {
    return '$n recortados';
  }

  @override
  String metaMatrix(String taxa, String sites) {
    return '$taxa taxones × $sites sitios';
  }

  @override
  String reviewSelectedCount(String selected, String total) {
    return '$selected de $total genomas activados';
  }

  @override
  String get progressSearchTitle => 'Buscando en NCBI por nombre taxonómico';

  @override
  String progressSearchingFor(String name) {
    return 'Buscando $name';
  }

  @override
  String get progressFetchTitle => 'Importando registros GenBank';

  @override
  String progressCounter(String current, String total) {
    return '$current de $total';
  }

  @override
  String progressSpecies(String current, String total) {
    return '$current de $total especies';
  }

  @override
  String get progressRateLabel => 'Tasa de solicitudes';

  @override
  String get progressRateWithKey => 'Hasta 10 solicitudes/s';

  @override
  String get progressRateDefault => 'Hasta 3 solicitudes/s';

  @override
  String get progressApiEnabled => 'Guardada y activa';

  @override
  String get progressApiDisabled => 'No proporcionada';

  @override
  String get progressCurrentLabel => 'Elemento actual';

  @override
  String get progressCurrentSpecies => 'Especie actual';

  @override
  String get progressStatusIdle => 'Esperando para empezar';

  @override
  String get progressFooter =>
      'SPLACE está consultando NCBI y actualizando el progreso continuamente.';

  @override
  String progressTaxonomyTitle(String source) {
    return 'Obteniendo taxonomía de $source';
  }

  @override
  String progressTaxonomyFooter(String source) {
    return 'SPLACE está consultando $source y actualizando los registros a medida que encuentra coincidencias.';
  }

  @override
  String get progressWorking => 'Procesando…';

  @override
  String get colStatus => 'Categoría (SALVE)';

  @override
  String fetchTaxonomyFrom(String source) {
    return 'Obtener taxonomía ($source)';
  }

  @override
  String taxonomyMatched(String source) {
    return 'Coincidencia encontrada en $source.';
  }

  @override
  String taxonomyNoMatch(String source) {
    return 'Sin coincidencias en $source para este nombre.';
  }

  @override
  String get taxonomyCardTitle => 'Taxonomía de los registros';

  @override
  String get taxonomySourceLabel => 'Fuente de la taxonomía';

  @override
  String get sourceGbifBody =>
      'Nombres aceptados y autoría del backbone del Global Biodiversity Information Facility.';

  @override
  String sourceNcbiBody(String withId, String total) {
    return 'Linaje del taxon id de cada registro GenBank ($withId de $total registros lo tienen; los demás se buscan por nombre).';
  }

  @override
  String get sourceSalveBody =>
      'Solo fauna brasileña. Evaluación nacional del ICMBio: taxonomía, categoría de riesgo de extinción y la ficha de la especie, por la API pública de SALVE.';

  @override
  String get badgeOnline => 'en línea';

  @override
  String get salveHint =>
      '**Solo especies de la fauna brasileña.** SALVE reúne las evaluaciones nacionales de riesgo de extinción del ICMBio; las especies de otros lugares no tienen ficha y aparecen como no encontradas (use GBIF o NCBI para ellas). Se leen en línea todas las secciones de cada ficha: rangos, categoría en los colores de la UICN y criterios, nombres comunes, distribución (endemismo, biomas, estados, cuencas), historia natural (hábitats, alimentación, reproducción, migración), población (tendencia, tiempo generacional), amenazas, usos, listas oficiales, historial de evaluaciones, acciones de conservación, áreas protegidas, investigación, evaluadores, cita y DOI.';

  @override
  String get poweredByDataFishing => 'Método de';

  @override
  String get metaSatKeep => 'mantener todas las posiciones';

  @override
  String metaSatExclude(String positions) {
    return 'eliminar posición $positions';
  }

  @override
  String get metaSatAmino => 'usar aminoácidos';

  @override
  String get apiKeyActive => 'Activa';

  @override
  String get apiKeyInactive => 'No configurada';

  @override
  String get apiKeyWithout => 'Sin clave';

  @override
  String get apiKeyWith => 'Con clave';

  @override
  String get apiKeyUsedFor =>
      'Usada en la búsqueda de NCBI, la descarga por accesiones y NCBI Taxonomy.';

  @override
  String get apiKeyHowTitle => 'Cómo obtener una clave';

  @override
  String get apiKeyHow1 => 'Inicie sesión en su cuenta de NCBI (gratuita).';

  @override
  String get apiKeyHow2 =>
      'Abra Account settings → API Key Management y pulse “Create an API Key”.';

  @override
  String get apiKeyHow3 => 'Copie la clave y péguela en el campo de abajo.';

  @override
  String get apiKeyOpenNcbi => 'Abrir cuenta de NCBI';

  @override
  String get apiKeyPaste => 'Pegar';

  @override
  String get apiKeyFieldLabel => 'Su API key';

  @override
  String get apiKeyLooksValid => 'Formato válido';

  @override
  String get apiKeyStored =>
      'Guardada solo en este equipo (Documentos/SPLACE/settings.json).';

  @override
  String get apiKeyEstimateTitle => 'Tiempo estimado de descarga';

  @override
  String apiKeyEstimateGenomes(String n) {
    return '$n genomas';
  }

  @override
  String apiKeyEstimateDetail(String batch, String rate) {
    return '$batch genomas por solicitud · hasta $rate solicitudes/s';
  }

  @override
  String apiKeyEstimateGain(String factor, String saved) {
    return '$factor× más rápido con clave · $saved menos';
  }

  @override
  String get apiKeyEstimateNote =>
      'Descarga de genomas mitocondriales que aún no están en la caché local, lectura incluida, a partir de tiempos medidos con la app el 08/11/2026 (Testudines, unos 75 kB por registro): 50, 200 y 700 genomas tardaron 1,6, 4,4 y 10,1 s con clave y 5,5, 16,9 y 64,6 s sin ella. La búsqueda y la lista de resultados suman unos segundos; los registros ya en caché cargan al instante. El tiempo real depende de su conexión y de la carga del NCBI.';

  @override
  String get apiKeyCurrent => 'Actual';

  @override
  String reviewFound(String n, String query) {
    return '$n genomas encontrados · $query';
  }

  @override
  String get reviewColOrganism => 'Organismo';

  @override
  String get reviewColStatus => 'Estado';

  @override
  String reviewUnverifiedNote(String n) {
    return '$n genomas están marcados como UNVERIFIED por NCBI. Recomendamos no importarlos.';
  }

  @override
  String get reviewDeselectUnverified => 'Desmarcar UNVERIFIED';

  @override
  String get reviewSelectMenu => 'Seleccionar';

  @override
  String reviewImportN(String n) {
    return 'Importar $n';
  }

  @override
  String get reviewNoMatch => 'Ningún genoma coincide con el filtro.';

  @override
  String reviewSelectedOf(String selected, String total) {
    return '$selected de $total seleccionados';
  }

  @override
  String metaTrees(String n) {
    return '$n árbol(es)';
  }

  @override
  String get itemQueued => 'En cola';

  @override
  String itemCached(String note) {
    return 'De la caché · $note';
  }

  @override
  String get itemSkipped => 'Omitido';

  @override
  String itemsRunning(String n) {
    return '$n en curso';
  }

  @override
  String itemsQueued(String n) {
    return '$n en cola';
  }

  @override
  String itemsDone(String n) {
    return '$n terminado(s)';
  }

  @override
  String itemsCached(String n) {
    return '$n de la caché (sin cambios)';
  }

  @override
  String itemsSkipped(String n) {
    return '$n omitido(s)';
  }

  @override
  String itemsFailed(String n) {
    return '$n con error';
  }

  @override
  String get trimToolNt =>
      'Elija una herramienta: trimAl o ClipKIT conservando codones enteros en los genes codificantes (recomendado), o trimAl columna a columna (parte los codones).';

  @override
  String get trimToolAa => 'trimAl recorta el alineamiento de proteínas.';

  @override
  String get trimToolCodon =>
      'El alineamiento de codones se recorta conservando codones enteros: con trimAl (-backtrans) o ClipKIT (--codon).';

  @override
  String runTrimWith(String tool) {
    return 'Recortar con $tool';
  }

  @override
  String get trimDescAutomated1 =>
      'trimAl analiza el alineamiento y elige por sí solo entre -gappyout y -strict. Recomendado para árboles de máxima verosimilitud.';

  @override
  String get trimDescStrict =>
      'Elimina columnas con muchos gaps y también columnas con baja similitud entre secuencias.';

  @override
  String get trimDescStrictplus =>
      'Como -strict, pero ajustado para árboles neighbour-joining; elimina más.';

  @override
  String get trimDescGappyout =>
      'Elimina solo las columnas con muchos gaps, con un umbral calculado a partir de la distribución de gaps del alineamiento.';

  @override
  String get trimDescNogaps =>
      'Elimina toda columna que tenga al menos un gap. Puede eliminar mucho.';

  @override
  String get trimDescNoallgaps =>
      'Elimina solo columnas formadas únicamente por gaps. La opción más ligera.';

  @override
  String get trimDescManual =>
      'Usted define los umbrales de gaps, similitud y consistencia.';

  @override
  String get clipkitDescGappy =>
      'Elimina columnas con más gaps que el umbral (90 %). Simple y conservador.';

  @override
  String get clipkitDescKpic =>
      'Mantiene solo columnas informativas para parsimonia y constantes; elimina el resto.';

  @override
  String get clipkitDescKpicGappy =>
      'kpic, y además elimina columnas con muchos gaps.';

  @override
  String get clipkitDescSmartGap =>
      'Elige el umbral de gaps a partir del propio alineamiento (predeterminado de ClipKIT).';

  @override
  String get clipkitRequired =>
      'El recorte codon-aware usa ClipKIT, que aún no está instalado. Instálelo abajo (requiere Python 3).';

  @override
  String get frameInfo =>
      'Dónde empieza el primer codón completo de este alineamiento. Base 1 (+1): el alineamiento empieza en la 1.ª posición de un codón, como en la mayoría de los CDS. +2 o +3: empieza 1 o 2 bases dentro de un codón. SPLACE elige el marco con menos codones de parada; se usa para traducir, recortar codones completos y separar la 1.ª, 2.ª y 3.ª posiciones del codón (prueba de saturación y particiones).';

  @override
  String frameOption(String n) {
    return 'Base $n (+$n)';
  }

  @override
  String frameTag(String n) {
    return 'Codones empiezan en la base $n';
  }

  @override
  String get viewerHint =>
      'Pase el ratón sobre una base para ver su posición · arrastre con el botón derecho para desplazarse · arrastre el borde de los nombres para cambiar el ancho, doble clic para ajustar.';

  @override
  String viewerColumn(String n) {
    return 'Columna $n';
  }

  @override
  String viewerSeqPos(String n) {
    return 'posición $n en la secuencia (sin contar gaps)';
  }

  @override
  String get viewerGap => 'gap';

  @override
  String viewerIdentical(String pct) {
    return '$pct% de las secuencias tienen esta base';
  }

  @override
  String viewerCodonPos(String n) {
    return '$n.ª posición del codón';
  }

  @override
  String get viewerNames => 'Secuencias';

  @override
  String siCacheInfo(String cached, String total, String todo) {
    return '$cached de $total búsquedas ya se hicieron antes y se reutilizarán; solo se enviarán $todo.';
  }

  @override
  String siCacheAll(String total) {
    return 'Las $total búsquedas ya se hicieron antes: los resultados aparecen al instante, no se envía nada.';
  }

  @override
  String get siCacheClear => 'Borrar resultados guardados';

  @override
  String get trimFrameWarning =>
      '**trimAl elimina columnas sueltas, no codones.** En un gen codificante alineado como ADN, eliminar una o dos bases de un codón desplaza el marco de lectura de todo lo que sigue: las posiciones 1.ª, 2.ª y 3.ª se mezclan. Sin \"Mantener el marco de lectura\", estos genes pasan a ser una partición única sin posiciones del codón y quedan fuera de la prueba de saturación.';

  @override
  String get frameLostTag => 'Marco de lectura perdido en el recorte';

  @override
  String get frameLostInfo =>
      'trimAl eliminó columnas sueltas de este gen codificante, así que su longitud puede no ser múltiplo de 3 y las posiciones del codón son desconocidas. Se concatena como una partición única y queda fuera de la prueba de saturación. Active \"Mantener el marco de lectura\" en 5e y recorte de nuevo para conservar las particiones por codón.';

  @override
  String frameShiftTag(String n) {
    return '$n secuencia(s) con marco desplazado';
  }

  @override
  String get frameShiftInfo =>
      'MAFFT colocó gaps que no son múltiplos de 3 dentro de estas secuencias: desde ese gap, sus bases quedan en la posición del codón equivocada de la columna. Es común al alinear genes codificantes como nucleótidos. Las particiones por codón y la prueba de saturación son aproximadas para ellas; el modo Codon-aware (5b) lo evita.';

  @override
  String get saturationRunning => 'Prueba de saturación';

  @override
  String saturationProgress(String done, String total) {
    return '$done de $total alineamientos';
  }

  @override
  String get matrixRunning => 'Construyendo la supermatriz';

  @override
  String get matrixRunningBody => 'Concatenando los alineamientos';

  @override
  String get iqtreeRunning => 'IQ-TREE: árbol de máxima verosimilitud';

  @override
  String get raxmlRunning => 'RAxML-NG: árbol de máxima verosimilitud';

  @override
  String treeRunSubtitle(String model, String elapsed) {
    return 'Modelo $model · tiempo $elapsed';
  }

  @override
  String get geneTreesProgress => 'Árboles por gen';

  @override
  String get localFooter => 'Se ejecuta en este equipo.';

  @override
  String get treeFooter =>
      'Se ejecuta en este equipo. Los análisis largos pueden tardar minutos u horas; siga cada etapa aquí.';

  @override
  String get treeDetailPrepare => 'Leyendo el alineamiento y las particiones';

  @override
  String get treeDetailModelStart => 'Probando modelos de sustitución';

  @override
  String treeDetailModels(String done, String total) {
    return 'Modelo $done de hasta $total probados';
  }

  @override
  String treeDetailStarts(String done, String total) {
    return 'Búsqueda ML desde el árbol inicial $done de $total';
  }

  @override
  String get treeDetailSearchStart => 'Construyendo los árboles iniciales';

  @override
  String treeDetailIteration(String n) {
    return 'Búsqueda del árbol: iteración $n';
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
  String get treeDetailSupport => 'Calculando el soporte de las ramas';

  @override
  String get treeDetailFinish => 'Guardando los resultados';

  @override
  String treeDetailBestModel(String model) {
    return 'Mejor modelo: $model';
  }

  @override
  String treeDetailUfboot(String n) {
    return '$n réplicas UFBoot';
  }

  @override
  String treeDetailConverged(String n) {
    return 'Bootstrap convergió tras $n';
  }

  @override
  String get treeStagePrepare => 'Alineamiento';

  @override
  String get treeStageModel => 'Selección del modelo';

  @override
  String get treeStageSearch => 'Búsqueda del árbol ML';

  @override
  String get treeStageSupport => 'Soporte de ramas';

  @override
  String get treeStageBootstrap => 'Bootstrap';

  @override
  String get treeStageFinish => 'Resultados';

  @override
  String get step9 => 'Visualizar Árbol';

  @override
  String metaViewable(String n) {
    return '$n árbol(es)';
  }

  @override
  String get viewTree => 'Ver árbol';

  @override
  String get treeViewIntro =>
      'Vea los árboles inferidos. Enraícelos con el grupo externo, en el punto medio o en cualquier rama, revise el soporte de cada clado, busque especies y exporte la figura.';

  @override
  String get treeChoose => 'Árbol';

  @override
  String get treeMissing => 'No se encontró el archivo del árbol.';

  @override
  String get treeKindMl => 'Árbol ML con soporte';

  @override
  String get treeKindConsensus => 'consenso del bootstrap';

  @override
  String get treeKindFbp => 'árbol ML, soporte bootstrap (FBP)';

  @override
  String get treeKindTbe => 'árbol ML, transfer bootstrap (TBE)';

  @override
  String get treeKindBest => 'mejor árbol ML (sin soporte)';

  @override
  String treeKindGene(String gene) {
    return 'árbol del gen $gene';
  }

  @override
  String treeParseError(String error) {
    return 'No se pudo leer este árbol: $error';
  }

  @override
  String get treeRootAsIs => 'Como se infirió';

  @override
  String get treeRootOutgroup => 'Grupo externo';

  @override
  String get treeRootMidpoint => 'Punto medio';

  @override
  String get treeRootManual => 'Enraizado en la rama elegida';

  @override
  String get treePhylogram => 'Filograma';

  @override
  String get treeCladogram => 'Cladograma';

  @override
  String get treeShowSupport => 'Soporte';

  @override
  String get treeItalic => 'Cursiva';

  @override
  String get treeSearch => 'Buscar especie';

  @override
  String get treeFont => 'Tamaño de fuente';

  @override
  String treeStats(String tips, String length) {
    return '$tips terminales · longitud total $length';
  }

  @override
  String get treeHint =>
      'Haga clic en una rama o nombre para seleccionarlo: enraizar allí, girar, colapsar o acortar. Ctrl+Z deshace.';

  @override
  String treeClade(String n) {
    return 'Clado con $n terminales';
  }

  @override
  String treeSupportValue(String value) {
    return 'soporte $value';
  }

  @override
  String get treeRootHere => 'Enraizar aquí';

  @override
  String get treeCollapse => 'Contraer';

  @override
  String get treeExpand => 'Expandir';

  @override
  String get treeClearSelection => 'Limpiar';

  @override
  String get treeSupportLegend => 'Soporte de ramas:';

  @override
  String get treeSupportInfoUfboot =>
      'Las etiquetas de IQ-TREE son SH-aLRT/UFBoot (o solo UFBoot); los colores usan el último valor. UFBoot ≥ 95 es soporte fuerte (Hoang et al. 2018); SH-aLRT ≥ 80 junto con él deja el clado bien soportado.';

  @override
  String get treeSupportInfoBootstrap =>
      'Bootstrap estándar (no paramétrico): ≥ 70 suele considerarse buen soporte (Hillis & Bull 1993).';

  @override
  String get siSelfTag => 'este registro';

  @override
  String siSelfInfo(String accessions) {
    return 'BLAST encontró la propia secuencia de este genoma en GenBank ($accessions). El dato está ahí, pero un registro no confirma su propio nombre, así que no cuenta como evidencia. Otros registros de la especie sí contarían.';
  }

  @override
  String siVerdictSelfOnly(String accessions, String current) {
    return 'GenBank tiene la secuencia de este genoma ($accessions), pero ningún otro registro de $current, así que GenBank no puede confirmar el nombre.';
  }

  @override
  String siVerdictNoReference(String current, String best, String identity) {
    return 'No se encontró otra secuencia de $current. La más cercana es $best ($identity%), por debajo del 98% —la distancia habitual entre especies—, así que el nombre actual no se cuestiona.';
  }

  @override
  String get siNoReference => 'Sin referencia';

  @override
  String get apiKeyEstimateColSize => 'Descarga';

  @override
  String get apiKeyEstimateColGain => 'Ganancia';

  @override
  String apiKeyEstimateSaved(String n, String saved) {
    return 'Para $n genomas, la clave ahorra unos $saved.';
  }

  @override
  String missTypeNotSelected(String type) {
    return 'anotado como $type';
  }

  @override
  String missTypeNotSelectedInfo(String gene, String type) {
    return '$gene está anotado en este registro como $type, un tipo de feature no seleccionado en el paso 3. Seleccione $type allí para incluirlo.';
  }

  @override
  String get missGeneOnly => 'solo feature \'gene\'';

  @override
  String missGeneOnlyInfo(String gene, String location) {
    return 'El registro tiene una feature \'gene\' para $gene ($location), pero no la feature CDS/rRNA, que es la que SPLACE lee. El gen está en la secuencia; la anotación está incompleta. Use un registro mejor anotado o anote el genoma (p. ej., con MitoS o MitoFinder).';
  }

  @override
  String missNoAnnotation(String types) {
    return 'sin $types anotados';
  }

  @override
  String missNoAnnotationInfo(String types) {
    return 'Este registro no tiene ningún $types anotado, así que no se puede extraer ningún gen. Probablemente es un envío sin anotar; anote el genoma (p. ej., con MitoS o MitoFinder) o elimínelo.';
  }

  @override
  String get missUnknownName => 'nombre no reconocido';

  @override
  String missUnknownNameInfo(String names, String gene) {
    return 'El registro tiene genes con nombres que SPLACE no reconoció ($names). Uno de ellos puede ser $gene: corrija el nombre en el paso 3 (nombres desconocidos).';
  }

  @override
  String get missPartial => 'genoma parcial';

  @override
  String missPartialInfo(String length, String typical, String gene) {
    return 'El registro tiene $length pb, muy por debajo del genoma típico de este conjunto ($typical pb): es un genoma parcial y no cubre $gene.';
  }

  @override
  String get missNotAnnotated => 'no anotado';

  @override
  String missNotAnnotatedInfo(String gene) {
    return '$gene no está entre los genes anotados de este registro: falta en la secuencia o los autores no lo anotaron.';
  }

  @override
  String get missNotInFiles => 'ausente en los archivos';

  @override
  String missNotInFilesInfo(String gene) {
    return 'Ningún archivo FASTA trae una secuencia de $gene para esta muestra.';
  }

  @override
  String missPresentIn(String present, String total, String missing) {
    return 'en $present de $total muestras · falta en $missing';
  }

  @override
  String missSampleCount(String missing, String total) {
    return 'faltan $missing de $total marcadores';
  }

  @override
  String missSummary(
    String markers,
    String selected,
    String samples,
    String total,
  ) {
    return '$markers de los $selected marcadores seleccionados faltan en al menos una muestra; a $samples de $total muestras les falta al menos un marcador.';
  }

  @override
  String get missMatrixNote =>
      'En la supermatriz (paso 7), el marcador ausente pasa a \'?\' para esa muestra cuando “Permitir secuencias ausentes” está activado; si no, los marcadores ausentes en algún taxón se excluyen. Algunos huecos son normales; una muestra sin la mayoría de los marcadores aporta poco y puede convenir eliminarla.';

  @override
  String missFixableNote(String n) {
    return '$n hueco(s) se pueden corregir en SPLACE (etiquetas naranjas). Pase el ratón por la etiqueta para ver cómo.';
  }

  @override
  String get missByMarker => 'Por marcador';

  @override
  String get missBySample => 'Por muestra';

  @override
  String get missHoverHint =>
      'Cada etiqueta indica por qué falta el marcador; pase el ratón para ver detalles.';

  @override
  String dupSummary(String genes, String genomes) {
    return '$genes gen(es) tienen más de una copia anotada en $genomes genoma(s).';
  }

  @override
  String get dupFastaHelp =>
      'Estas muestras tienen más de una secuencia del mismo locus en los archivos FASTA. La columna “Usada” muestra lo que entra en el alineamiento, según la opción “Duplicados y nombres”.';

  @override
  String get dupIdentical => 'copias idénticas';

  @override
  String get dupIdenticalInfo =>
      'Misma secuencia (p. ej., copias de la repetición invertida en plastomas): cualquier copia da el mismo resultado.';

  @override
  String get dupContained => 'fragmento de otra copia';

  @override
  String get dupContainedInfo =>
      'La copia más corta es parte de la más larga: normalmente una anotación duplicada o truncada. Use la copia completa.';

  @override
  String get dupDifferent => 'copias diferentes';

  @override
  String get dupDifferentInfo =>
      'Las copias difieren en secuencia (parálogos, copias tipo NUMT o errores de anotación). La elección cambia el alineamiento: prefiera la copia más cercana a los otros genomas.';

  @override
  String get siCardTitle => 'Identificación de especies';

  @override
  String get siCardBody =>
      'Compara un gen de código de barras (COI, 16S…) de cada genoma con GenBank (BLAST) y BOLD para detectar registros mal identificados antes de obtener la taxonomía, que depende de los nombres.';

  @override
  String get siCardOptional => 'Opcional';

  @override
  String get siCardRerun => 'Identificar de nuevo';

  @override
  String siSumChecked(String n) {
    return '$n genomas verificados';
  }

  @override
  String siSumConfirmed(String n) {
    return '$n nombres confirmados';
  }

  @override
  String siSumFlagged(String n) {
    return '$n posibles errores de identificación';
  }

  @override
  String siSumNoReference(String n) {
    return '$n sin referencia';
  }

  @override
  String get headerDefault => 'Encabezado predeterminado · se usa al alinear';

  @override
  String viewerMafftColumn(String n) {
    return 'columna $n del alineamiento MAFFT';
  }

  @override
  String get satSourceMafft =>
      'Probado en el alineamiento de MAFFT (sin recorte).';

  @override
  String get satSourceClipkit =>
      'Probado en el alineamiento recortado (ClipKIT en las proteínas, retrotraducido): solo se quitaron codones enteros, así que cada base conserva su posición en el codón.';

  @override
  String satSourceTrimmed(String kept, String total) {
    return 'Probado en el alineamiento recortado: trimAl conservó $kept de los $total codones de MAFFT, y SPLACE quitó solo codones enteros, así que cada columna conserva la posición en el codón que tenía en el alineamiento de MAFFT (pase el ratón sobre una columna recortada en el paso 5 para ver de qué columna de MAFFT viene).';
  }

  @override
  String satSourceFrame(String frame) {
    return 'Los codones ahora empiezan en la base 1 porque se quitó el codón incompleto antes de la base $frame del alineamiento de MAFFT; las posiciones son las mismas.';
  }

  @override
  String get preserveCaseHelp => 'Mantiene mayúsculas y minúsculas';

  @override
  String get configureTrimming => 'Configurar el recorte';

  @override
  String panelTrimParams(String tool) {
    return 'Parámetros de $tool';
  }

  @override
  String trimmedCount(String n, String total) {
    return '$n de $total marcadores recortados';
  }

  @override
  String get showRemovedColumns => 'Mostrar regiones eliminadas';

  @override
  String get nextStop => 'Siguiente stop';

  @override
  String get translationTable => 'Tabla de traducción';

  @override
  String translationTableTitle(String code, String name) {
    return 'Código genético $code · $name';
  }

  @override
  String get translationTableBody =>
      'Codones del código genético elegido arriba para la vista de aminoácidos (nunca se deduce de los datos). Los stops están en rojo; la insignia cuenta cuántos stops internos de este alineamiento usan cada uno.';

  @override
  String codonSecondBase(String base) {
    return '2.ª $base';
  }

  @override
  String codonDiffers(String codon, String aa, String standard) {
    return '$codon → $aa en este código ($standard en el código estándar)';
  }

  @override
  String get codonDiffersLegend => 'Difiere del código estándar';

  @override
  String stopCodonsFound(String list) {
    return 'Stops internos en este alineamiento: $list';
  }

  @override
  String removedColumns(String n, String total) {
    return '$n de $total columnas eliminadas';
  }

  @override
  String get removedColumnsInfo =>
      'Se muestra todo el alineamiento MAFFT; las columnas en gris fueron eliminadas por el recorte. Desactive \"Mostrar regiones eliminadas\" para ver solo las columnas conservadas.';

  @override
  String get removedLegend => 'Eliminado';

  @override
  String get viewerRemoved => 'eliminada por el recorte';

  @override
  String viewerCodon(String codon) {
    return 'codón $codon';
  }

  @override
  String get plotObservedVsCorrectedBody =>
      'Cada punto es un par de secuencias. Puntos sobre la línea discontinua: todas las sustituciones son visibles; puntos que se aplanan por debajo: sustituciones múltiples ocultas (saturación).';

  @override
  String get plotTsTvBody =>
      'Las transiciones (rojo) suelen saturarse primero: cuando dejan de crecer mientras las transversiones (azul) siguen subiendo, el conjunto está saturado.';

  @override
  String plotPositions(String model) {
    return 'p de cada posición del codón × distancia $model';
  }

  @override
  String get plotPositionsBody =>
      'Mismos pares, un color por posición del codón, frente a la distancia del gen completo. Una posición cuyos puntos se aplanan mientras las otras siguen subiendo está saturada (normalmente la 3.ª).';

  @override
  String get allPositions => 'todas las posiciones';

  @override
  String get compositionTitle => 'Composición y variabilidad';

  @override
  String get compositionInfo =>
      'Prueba chi-cuadrado de homogeneidad de las frecuencias de bases entre taxones (prueba de contingencia, como en PAUP*) y el número de taxones cuya composición difiere de la del conjunto (gl = 3, p < 0,05, como informa IQ-TREE). Las posiciones saturadas suelen variar en composición (p. ej. GC3), lo que también puede sesgar el árbol.';

  @override
  String get metricVariable => 'Variables';

  @override
  String get metricInformative => 'Informativos (parsimonia)';

  @override
  String get metricTaxaFailing => 'Taxones rechazados';

  @override
  String get uncorrectablePairs => 'pares sin corrección posible';

  @override
  String get heatmapIntensity =>
      'El tono dentro de cada color sigue la pendiente y los pares sin corrección.';

  @override
  String get treeRootLabel => 'Raíz';

  @override
  String get treeShapeLabel => 'Forma';

  @override
  String get treeSupportValues => 'Valores';

  @override
  String get treeSupportDots => 'Puntos';

  @override
  String get treeSupportHidden => 'Ocultar';

  @override
  String get treeDisplay => 'Mostrar';

  @override
  String get treeAlignLabels => 'Alineadas a la derecha';

  @override
  String get treeExport => 'Exportar';

  @override
  String treeMatches(String n) {
    return '$n encontrado(s)';
  }

  @override
  String treeHoverTip(String name, String length) {
    return '$name · longitud de rama $length';
  }

  @override
  String treeHoverClade(String n, String support, String length) {
    return 'Clado con $n terminales · soporte $support · longitud de rama $length — haga clic para seleccionarlo';
  }

  @override
  String get treeOutgroupLegend => 'Grupo externo';

  @override
  String get ord1 => '1.ª';

  @override
  String get ord2 => '2.ª';

  @override
  String get ord3 => '3.ª';

  @override
  String get listAnd => 'y';

  @override
  String codonPositionN(String n) {
    return '$n posición';
  }

  @override
  String reasonInsufficient(String sites, String pairs) {
    return 'Muy pocos datos para evaluar ($sites sitios, $pairs pares de secuencias).';
  }

  @override
  String reasonUndefined(String pct) {
    return 'El $pct% de los pares de secuencias difieren tanto que la distancia no puede corregirse.';
  }

  @override
  String reasonLowDivergence(String p) {
    return 'Las secuencias difieren muy poco (p media = $p), así que la saturación no es un problema.';
  }

  @override
  String reasonSlopeOk(String slope) {
    return 'Saturación de PhyKIT $slope: las diferencias observadas aún crecen junto con la distancia en el árbol.';
  }

  @override
  String reasonSlopeModerate(String slope, String limit) {
    return 'Saturación de PhyKIT $slope (por debajo de $limit): las diferencias observadas crecen más despacio que la distancia en el árbol.';
  }

  @override
  String reasonSlopeSaturated(String slope, String limit) {
    return 'Saturación de PhyKIT $slope (por debajo de $limit): las diferencias observadas casi dejaron de crecer con la distancia en el árbol; muchos cambios fueron borrados por otros posteriores.';
  }

  @override
  String reasonIss(String iss) {
    return 'Iss $iss (0,8 o más): los sitios parecen casi aleatorios.';
  }

  @override
  String get recTitle => 'Qué hacer con las posiciones del codón';

  @override
  String recBasedOn(String taxa, String length) {
    return 'Según la matriz concatenada ($taxa taxones, $length pb). Cada posición se prueba por separado.';
  }

  @override
  String get recWhatToDo => 'Sugerencia';

  @override
  String recExcludeOne(String positions) {
    return 'Elimine la $positions posición del codón de la matriz (paso 7b) y mantenga las otras dos, particionadas por posición del codón.';
  }

  @override
  String recExcludeMany(String positions) {
    return 'Elimine las posiciones $positions del codón de la matriz (paso 7b) y mantenga la restante.';
  }

  @override
  String get recWhy =>
      'Por qué: una posición saturada cambió tantas veces que sus diferencias ya no muestran hace cuánto se separaron los linajes. Aporta sobre todo ruido (homoplasia) y puede agrupar ramas largas.';

  @override
  String recCautionOne(String positions) {
    return 'La $positions posición muestra saturación moderada: manténgala en su propia partición y compruebe si el árbol cambia sin ella.';
  }

  @override
  String recCautionMany(String positions) {
    return 'Las posiciones $positions muestran saturación moderada: manténgalas en sus propias particiones y compruebe si el árbol cambia sin ellas.';
  }

  @override
  String applyExclude(String positions) {
    return 'Marcar $positions en el paso 7b';
  }

  @override
  String get applyKeepAll =>
      'Dejar todos los \"Excluir posición\" sin marcar en el paso 7b';

  @override
  String applyPerLocusInfo(String loci) {
    return 'Cada gen excluye solo las posiciones saturadas en él. Genes con un resultado distinto de la sugerencia general: $loci.';
  }

  @override
  String recLociDiffer(String n) {
    return '$n gen(es) tienen un resultado distinto de la sugerencia general; vea sus filas en el mapa de calor.';
  }

  @override
  String get recStatusKeep => 'El paso 7b ahora mantiene todas las posiciones.';

  @override
  String recStatusExclude(String positions) {
    return 'El paso 7b ahora excluye: $positions posición(es).';
  }

  @override
  String recStatusPerLocus(String n) {
    return 'El paso 7b ahora excluye posiciones por gen ($n genes).';
  }

  @override
  String get saturationHeatmapHelp =>
      'Cada celda muestra el valor de saturación de PhyKIT (pendiente de p sobre la distancia en el árbol) — 1 significa sin saturación, cerca de 0 saturada. Azul es OK, rojo es saturada.';

  @override
  String get rankSpecies => 'Especie';

  @override
  String get rankGenus => 'Género';

  @override
  String get rankFamily => 'Familia';

  @override
  String get rankOrder => 'Orden';

  @override
  String get rankClass => 'Clase';

  @override
  String get rankPhylum => 'Filo';

  @override
  String get rankKingdom => 'Reino';

  @override
  String get tipFieldOriginal => 'Nombre como en el árbol';

  @override
  String get tipFieldAccession => 'Acceso';

  @override
  String get tipFieldVoucher => 'Voucher';

  @override
  String get treeBold => 'Negrita';

  @override
  String get treeParens => 'Entre paréntesis';

  @override
  String get treeRemovePart => 'Quitar esta parte';

  @override
  String get treeAddPart => 'Añadir campo';

  @override
  String get treeMovePart => 'Mover a la izquierda';

  @override
  String get treeLabelFormat => 'Formato de las etiquetas';

  @override
  String get treeLabelParts => 'Partes de la etiqueta';

  @override
  String get treeLabelPartsInfo =>
      'Cada etiqueta se arma con estos campos, en este orden. I = cursiva, B = negrita, ( ) = entre paréntesis. Los valores vienen de los registros de GenBank y de la taxonomía obtenida en el paso 2; un campo que la punta no tiene se omite.';

  @override
  String get treeTabTree => 'Árbol';

  @override
  String get treeTabBranches => 'Ramas y soporte';

  @override
  String get treeTabLabels => 'Etiquetas';

  @override
  String get treeTabColours => 'Colores';

  @override
  String get treeOrderLabel => 'Orden y rotación';

  @override
  String get treeLadderUp => 'Escalera ↑';

  @override
  String get treeLadderDown => 'Escalera ↓';

  @override
  String get treeFlip => 'Invertir árbol';

  @override
  String get treeRotateClade => 'Girar clado';

  @override
  String get treeRotateInfo =>
      'Girar un clado intercambia el orden de sus ramas; las relaciones no cambian. Seleccione un clado haciendo clic en su rama y luego gírelo. Invertir árbol refleja el árbol entero.';

  @override
  String get treeLineWidth => 'Grosor de las líneas';

  @override
  String get treePointSize => 'Tamaño de los puntos';

  @override
  String get treePoints => 'Puntos';

  @override
  String get treeTipDots => 'Puntos en las puntas';

  @override
  String get treeSupportMin => 'Mostrar soporte desde';

  @override
  String get treeSupportMinInfo =>
      'Solo los nodos con soporte igual o mayor a este valor reciben número o punto (por ejemplo, puntos solo en los nodos con UFBoot ≥ 95).';

  @override
  String get treeSupportAll => 'Todos';

  @override
  String get treeLongBranches => 'Ramas largas';

  @override
  String get treeLongBranchesInfo =>
      'Una rama acortada se dibuja más corta de lo que es y se marca con //, para dejar espacio al resto del árbol. Cuentan como largas las ramas de más de 5× la longitud mediana y más del 30% de la altura del árbol. También puede seleccionar cualquier rama y acortarla. La longitud real queda en la barra de información y en el Newick exportado.';

  @override
  String treeShortenLong(String n) {
    return 'Acortar ramas largas ($n)';
  }

  @override
  String get treeRestoreLengths => 'Restaurar longitudes';

  @override
  String get treeShortenBranch => 'Acortar rama //';

  @override
  String get treeRestoreBranch => 'Restaurar longitud';

  @override
  String get treeShortenedLegend => 'rama acortada (fuera de escala)';

  @override
  String get treeColorRows => 'Colorear filas por';

  @override
  String get treeColorRowsInfo =>
      'Cada taxón del rango elegido recibe un color de fondo en toda la fila, con su nombre a la derecha. Los rangos vienen de la taxonomía obtenida en el paso 2 (o de la especie y el género en los nombres de las puntas).';

  @override
  String get treeColorNone => 'Ninguno';

  @override
  String get treeNoTaxonomy =>
      'No hay taxonomía para estas puntas. Obténgala en el paso 2.';

  @override
  String get treeExportInfo =>
      'SVG y PNG son la figura tal como la ve (etiquetas, colores, leyenda, ramas acortadas), sin las marcas de selección ni de búsqueda. El Newick tiene la topología y las longitudes reales.';

  @override
  String get treeUndo => 'Deshacer';

  @override
  String get treeRedo => 'Rehacer';

  @override
  String get treeReset => 'Empezar de nuevo';

  @override
  String get zipTitle => 'Guardar todo';

  @override
  String get zipBody =>
      'Un ZIP con las secuencias, alineamientos, tablas de saturación, matrices, particiones, árboles y registros de este proyecto.';

  @override
  String get concatHowTitle => 'Cómo se arma la supermatriz';

  @override
  String get concatHow1 =>
      'El alineamiento de cada locus se coloca lado a lado, una fila por taxón (emparejada por el encabezado FASTA).';

  @override
  String get concatHow2 =>
      'Un taxón sin un locus recibe ? en ese tramo (dato faltante), así todas las filas tienen la misma longitud.';

  @override
  String get concatHow3 =>
      'El archivo de particiones indica a IQ-TREE / RAxML-NG dónde empieza y termina cada locus (y cada posición del codón), para que cada parte tenga su propio modelo de sustitución.';

  @override
  String get concatPictureCaption => '? = taxón sin ese locus';

  @override
  String get concatBlockSource => 'Qué alineamientos unir';

  @override
  String get concatBlockSourceBody =>
      'Los alineamientos recortados del paso 5 (sin las columnas mal alineadas) o la salida de MAFFT tal cual.';

  @override
  String get concatBlockPartitions => 'Cómo dividir la matriz en particiones';

  @override
  String get concatBlockPartitionsBody =>
      'Cada partición recibe su propio modelo de sustitución en la búsqueda del árbol. Aquí también se pueden quitar posiciones del codón.';

  @override
  String get schemeByLocusCodonBody =>
      'Cada locus es una partición y, en los loci codificantes, cada posición del codón (1.ª, 2.ª, 3.ª) tiene la suya. Las posiciones evolucionan a tasas muy distintas; con MFP+MERGE, IQ-TREE une las particiones que se comportan igual.';

  @override
  String get schemeByLocusBody =>
      'Una partición por locus: cada gen tiene su modelo y sus tres posiciones del codón lo comparten. Bueno para rRNAs, loci no codificantes o con pocos sitios.';

  @override
  String get schemeByCodonBody =>
      'Tres particiones compartidas por todos los loci codificantes (todas las 1.as posiciones juntas, todas las 2.as, todas las 3.as); cada locus no codificante sigue siendo una partición. Menos modelos que por locus y posición.';

  @override
  String get schemeNoneBody =>
      'Un modelo para toda la matriz (sin archivo de particiones). Más simple, pero suele ajustar peor cuando genes o posiciones evolucionan distinto.';

  @override
  String get codonExclusionHelp =>
      'Las posiciones quitadas salen solo de los loci codificantes. El paso 6 (saturación) puede completar esto.';

  @override
  String get concatSummaryTitle => 'Qué se va a generar';

  @override
  String get concatSummaryLoci => 'Loci';

  @override
  String get concatSummarySource => 'Alineamientos';

  @override
  String get concatNoOutgroup => 'ninguno (punto medio)';

  @override
  String get concatNothingExcluded => 'nada';

  @override
  String get treeXScale => 'Escala horizontal';

  @override
  String get treeXScaleInfo =>
      'Estira (por encima de 1) o encoge (por debajo de 1) las ramas en el eje x. Por encima de 1 la figura queda más ancha que el ancho definido.';

  @override
  String get treeYScale => 'Escala vertical';

  @override
  String get treeYScaleInfo =>
      'Espacio entre las filas en el eje y (1 = 1,7 veces el tamaño de la etiqueta).';

  @override
  String get treeBranchLengths => 'Longitud de las ramas';

  @override
  String get treeShowLengths => 'Mostrar como etiqueta';

  @override
  String get treeSupportFont => 'Tamaño del texto del soporte';

  @override
  String get treeProportional => 'Tamaño según el soporte';

  @override
  String get treeShowLabels => 'Mostrar etiquetas de las puntas';

  @override
  String get treeTabFigure => 'Figura y exportación';

  @override
  String get treeScaleBar => 'Barra de escala';

  @override
  String get treeShowScaleBar => 'Mostrar';

  @override
  String get treeScaleBarLength => 'Longitud de la barra de escala';

  @override
  String get treeScaleBarLengthInfo =>
      'En sustituciones por sitio. 0 elige un valor redondo de cerca de un quinto del árbol.';

  @override
  String get treeFigureLegend => 'Leyenda';

  @override
  String get treeLegendInFigure => 'Dibujar en la figura';

  @override
  String get treeFigureLegendInfo =>
      'Los colores del soporte, el grupo externo, las ramas acortadas y los colores de las filas se dibujan debajo del árbol y van con el PNG y el SVG.';

  @override
  String get treeAtTips => 'En las puntas';

  @override
  String get treeLegendMarks => 'Marcas';

  @override
  String get treeBranchColor => 'Color de las ramas';

  @override
  String get treeSupportPosition => 'Posición del soporte';

  @override
  String get treeAboveBranch => 'Sobre la rama';

  @override
  String get treeBelowBranch => 'Bajo la rama';

  @override
  String get treeSupportColor => 'Color';

  @override
  String get treeFontFamily => 'Fuente';

  @override
  String get treeLabelColor => 'Color de las etiquetas';

  @override
  String get treeLabelPosition => 'Posición de las etiquetas';

  @override
  String get treeLabelGap => 'Distancia a la punta';

  @override
  String get treeLabelGapInfo =>
      'Espacio entre el final de la rama (o la columna alineada) y la etiqueta.';

  @override
  String get treeColorApply => 'Aplicar los colores a';

  @override
  String get treeColorApplyRows => 'Fondo de la fila';

  @override
  String get treeColorApplyBranches => 'Ramas';

  @override
  String get treeColorApplyLabels => 'Etiquetas';

  @override
  String get treeCladeColors => 'Colores de clados';

  @override
  String get treeCladeColorsInfo =>
      'Seleccione un clado haciendo clic en su rama y asígnele un color: sus ramas lo toman. Deshacer también funciona aquí.';

  @override
  String get treeColorClade => 'Colorear clado';

  @override
  String get treeRemoveColor => 'Quitar color';

  @override
  String treeClearCladeColors(String n) {
    return 'Quitar todos ($n)';
  }

  @override
  String get treeBackground => 'Fondo de la figura';

  @override
  String get treeLegendPosition => 'Posición';

  @override
  String get treeLegendBelow => 'Bajo el árbol';

  @override
  String get treeLegendRight => 'A la derecha del árbol';

  @override
  String get treeLegendStyle => 'Estilo';

  @override
  String get treeLegendBox => 'En un recuadro';

  @override
  String get treeLegendTitle => 'Título de la leyenda';

  @override
  String get treeLegendTitleHint => 'p. ej. Leyenda';

  @override
  String get treeLegendFont => 'Tamaño del texto de la leyenda';

  @override
  String get treeFigureTitle => 'Título de la figura';

  @override
  String get treeFigureTitleHint => 'Opcional, se dibuja sobre el árbol';

  @override
  String get treeDpi => 'Resolución del PNG';

  @override
  String get treeDpiInfo =>
      'Puntos por pulgada del PNG (300 para la mayoría de las revistas, 600 para dibujos de línea). El DPI queda grabado en el archivo. Las figuras muy grandes se limitan a 16 000 px en el lado mayor.';

  @override
  String get treePngBackground => 'Fondo del PNG';

  @override
  String get treeTransparent => 'Transparente';

  @override
  String get treePickColor => 'Elija un color';

  @override
  String get concatBlockOutgroup =>
      'Qué taxones enraízan el árbol (grupo externo)';

  @override
  String get concatBlockOutgroupBody =>
      'Elija uno o más taxones fuera del grupo estudiado; los árboles de los pasos 8 y 9 se enraízan en ellos. Sin grupo externo, el visor enraíza en el punto medio.';

  @override
  String get schemeOneModel => 'un modelo para todo';

  @override
  String schemePartitionsCount(String n) {
    return '$n partición(es)';
  }

  @override
  String get resCpu => 'CPU';

  @override
  String get resCpuInfo =>
      'Uso de CPU de la máquina, incluidos MAFFT, trimAl, IQ-TREE y las demás herramientas que ejecuta SPLACE.';

  @override
  String get resMemory => 'RAM';

  @override
  String resMemoryInfo(String app, String used, String total) {
    return 'SPLACE usa $app. Máquina: $used de $total en uso (la barra).';
  }

  @override
  String resMemoryInfoApp(String app) {
    return 'SPLACE usa $app.';
  }

  @override
  String get resDisk => 'Disco';

  @override
  String resFree(String v) {
    return '$v libres';
  }

  @override
  String resDiskInfo(String free, String total) {
    return '$free libres de $total en el disco de la carpeta del proyecto (la barra es el espacio en uso).';
  }

  @override
  String get resDiskInfoNone =>
      'Espacio en el disco de la carpeta del proyecto.';

  @override
  String get resNetwork => 'Red';

  @override
  String resNetworkInfo(String total) {
    return 'Datos que SPLACE intercambia con NCBI, GBIF y BLAST, por segundo. Desde que se abrió la app: $total.';
  }

  @override
  String get methodsTitle => 'Texto de métodos';

  @override
  String get methodsCopy => 'Copiar texto y referencias';

  @override
  String get methodsHelp =>
      'Escrito a partir de lo que se hizo en este proyecto (datos, alineamiento, recorte, saturación, supermatriz y árboles). Revíselo y adáptelo a su manuscrito.';

  @override
  String get citeReferences => 'Referencias citadas en el texto';

  @override
  String get citeTools => 'Herramientas usadas por SPLACE';

  @override
  String get treeBranchStyle => 'Estilo de las ramas';

  @override
  String get treeBranchSquare => 'Ángulos rectos';

  @override
  String get treeBranchRounded => 'Esquinas redondeadas';

  @override
  String get treeBranchSlanted => 'Diagonales';

  @override
  String get treeZoom => 'Zoom (solo en pantalla)';

  @override
  String get treeZoomFit => 'Ajustar a la ventana';

  @override
  String get treeZoomInfo =>
      'Amplía o reduce la vista de la figura en esta pantalla. La figura exportada mantiene el ancho definido en Figura y exportación, sea cual sea el monitor.';

  @override
  String get treeTipDotsLabel => 'Puntas';

  @override
  String get treeLabelCustom => 'Personalizado';

  @override
  String get treeMarkOutgroup => 'Caja gris detrás';

  @override
  String get treeRenamed => 'Etiquetas editadas';

  @override
  String get treeRenamedInfo =>
      'Haga clic en una punta y elija Editar etiqueta para escribir un texto propio (con cursiva y negrita). Borre el texto para volver a la etiqueta construida.';

  @override
  String get treeRenameTip => 'Editar etiqueta';

  @override
  String treeRestoreNames(String n) {
    return 'Restaurar etiquetas editadas ($n)';
  }

  @override
  String get treeBarsByRank => 'Barras con nombres de taxones';

  @override
  String get treeBarsByRankInfo =>
      'Dibuja una barra junto a cada grupo de puntas que comparten el taxón de este rango (por ejemplo, una barra por familia), con su nombre. Los taxones separados en el árbol reciben una barra por parte.';

  @override
  String get treeHighlightOpacity => 'Opacidad de las cajas';

  @override
  String get treeHighlightClade => 'Destacar clado';

  @override
  String get treeRemoveHighlight => 'Quitar destaque';

  @override
  String get treeCladeLabel => 'Nombrar clado';

  @override
  String treeClearHighlights(String n) {
    return 'Quitar destaques ($n)';
  }

  @override
  String treeClearCladeLabels(String n) {
    return 'Quitar nombres de clados ($n)';
  }

  @override
  String get treeRemoveFromFigure => 'Quitar del árbol';

  @override
  String get treeFigureWidth => 'Ancho de la figura';

  @override
  String get treeFigureWidthInfo =>
      'Ancho real de la figura exportada. Los tamaños de fuente están en puntos, así que una etiqueta de 9 pt se imprime con 9 pt en este ancho. La altura depende del número de puntas y de la escala vertical.';

  @override
  String get treeWidthPreset => 'Ancho de revista';

  @override
  String get treeWidthOneCol => '1 columna (8,5 cm)';

  @override
  String get treeWidthOneHalf => '1,5 columnas (12 cm)';

  @override
  String get treeWidthTwoCol => 'Página completa (17,5 cm)';

  @override
  String get treeWidthCustom => 'Personalizado';

  @override
  String get treeMargin => 'Margen';

  @override
  String get treeTitleFont => 'Tamaño del título';

  @override
  String treeFigureSize(String w, String h, String pw, String ph, String dpi) {
    return 'Figura $w × $h cm · PNG $pw × $ph px a $dpi dpi';
  }

  @override
  String get treeTextHint => 'Texto';

  @override
  String get treeTextEmptyHint => 'Déjelo vacío para quitarlo.';

  @override
  String get concatHowOpen =>
      'Haga clic para ver cómo se colocan los loci uno al lado del otro';

  @override
  String get treeOpenEditor => 'Abrir el editor del árbol';

  @override
  String get treeEditorTitle => 'Editor del árbol';

  @override
  String get treeEditorHint =>
      'Vista previa de la figura a su tamaño real. Haga clic en ella, o en el botón de arriba, para editarla en la ventana completa.';

  @override
  String get treeCompare => 'Comparar árboles';

  @override
  String get treeModeEdit => 'Editar figura';

  @override
  String get treeDone => 'Listo';

  @override
  String get treeExportMenu => 'Exportar';

  @override
  String get treeZoomIn => 'Acercar (Ctrl + rueda)';

  @override
  String get treeZoomOut => 'Alejar (Ctrl + rueda)';

  @override
  String get treeFitWidth => 'Ajustar al ancho';

  @override
  String get treeFitPage => 'Figura completa';

  @override
  String get treeCanvasHint =>
      'Arrastre para mover · la rueda desplaza, Ctrl + rueda hace zoom · clic selecciona, Shift + clic selecciona el clado que une dos puntas · clic derecho para opciones · Supr elimina · Esc limpia';

  @override
  String get treeTabCompare => 'Comparar';

  @override
  String get treeCompareWith => 'Comparar con';

  @override
  String get treeUntangle => 'Desenredar';

  @override
  String get treeUntangleText => 'Girar el otro árbol';

  @override
  String get treeUntangleInfo =>
      'Gira los clados del otro árbol (su topología no cambia) para que cada punta quede frente a su par en el árbol en edición, con el mínimo de líneas cruzadas.';

  @override
  String get treeConflicts => 'Clados diferentes';

  @override
  String get treeConflictsText => 'Marcarlos';

  @override
  String get treeConflictsInfo =>
      'Las ramas que llevan a un clado que el otro árbol no tiene se dibujan más gruesas, en este color.';

  @override
  String get treeConflictColor => 'Color de los clados diferentes';

  @override
  String get treeConnectors => 'Líneas entre los árboles';

  @override
  String get treeConnectorsGroup => 'Colores de las filas';

  @override
  String get treeConnectorsInfo =>
      'Con las filas coloreadas por un rango (pestaña Colores), cada línea toma el color de su taxón; si no, es gris.';

  @override
  String get treeCompareGap => 'Espacio entre los árboles';

  @override
  String get treeCompareResult => 'Resultado';

  @override
  String treeCompareCommon(String n) {
    return '$n puntas en ambos árboles';
  }

  @override
  String treeCompareOnly(String a, String b) {
    return '$a solo a la izquierda · $b solo a la derecha';
  }

  @override
  String treeCompareRf(String rf, String nrf) {
    return 'Distancia de Robinson–Foulds $rf (normalizada $nrf)';
  }

  @override
  String treeCompareShared(String n, String a, String b) {
    return '$n clados en común · $a solo a la izquierda · $b solo a la derecha';
  }

  @override
  String get treeCompareRfInfo =>
      'Cuenta los clados que están en un árbol y no en el otro (sin raíz, sobre las puntas que ambos tienen). 0 indica la misma topología; el valor normalizado divide por la mayor distancia posible.';

  @override
  String get treeCompareLegend => 'Clado ausente en el otro árbol';

  @override
  String get treeCompareTooFew =>
      'Los árboles comparten menos de 4 puntas: sus clados no se pueden comparar.';

  @override
  String treeCompareError(String e) {
    return 'No se pudo leer el otro árbol: $e';
  }

  @override
  String get treeCompareIntro =>
      'El árbol en edición y el elegido aquí se dibujan frente a frente, uno en espejo, con una línea que une cada punta con su par. Los estilos (fuentes, ramas, soporte, colores) vienen de las otras pestañas.';

  @override
  String get treeSelectParent => 'Clado superior';

  @override
  String get treeLeft => 'Izquierda';

  @override
  String treeDetailMoose(String n) {
    return 'MOOSE está eligiendo los modelos ($n partición(es))';
  }

  @override
  String get treeKindRbs => 'árbol ML, bootstrap rápido (RBS)';

  @override
  String get treeKindSh => 'árbol ML, SH-aLRT';

  @override
  String get treeKindEbg => 'árbol ML, bootstrap predicho (EBG)';

  @override
  String raxmlOld(String version) {
    return 'Se encontró RAxML-NG $version. MOOSE, la búsqueda rápida y los soportes SH-aLRT, EBG y RBS necesitan la versión 2.0 o más reciente (SPLACE incluye la 2.0.4).';
  }

  @override
  String get raxmlModelsHeading => 'Modelos de sustitución';

  @override
  String get raxmlMoose => 'MOOSE (selección de modelos del propio RAxML-NG)';

  @override
  String get raxmlMooseBody =>
      'Un modelo por partición, elegido en la misma ejecución. Recomendado.';

  @override
  String get raxmlModelFinderBody =>
      'IQ-TREE elige los modelos y puede unir particiones; luego RAxML-NG construye el árbol.';

  @override
  String get raxmlFixedModel => 'Un modelo fijo';

  @override
  String get raxmlFixedBody =>
      'El mismo modelo en todas las particiones (por ejemplo GTR+G).';

  @override
  String get raxmlCriterion => 'Criterio';

  @override
  String get raxmlSearch => 'Búsqueda del árbol';

  @override
  String get raxmlSearchAdaptive => 'Adaptativa (por defecto)';

  @override
  String get raxmlSearchFast => 'Rápida (--fast)';

  @override
  String get raxmlSupportInfo =>
      'FBP, TBE y RBS provienen de árboles de réplicas de bootstrap. SH-aLRT y EBG se calculan solo sobre el árbol de ML y tardan segundos. FBP y RBS no se pueden usar juntos.';

  @override
  String get treeLayoutLabel => 'Disposición';

  @override
  String get treeLayoutRect => 'Rectangular';

  @override
  String get treeLayoutCircle => 'Circular';

  @override
  String get treeLayoutInfo =>
      'Circular: la raíz en el centro y las puntas alrededor, útil para árboles con muchas puntas. Los nombres siguen el radio; los colores de taxón forman un anillo fuera de los nombres y las barras de clado se vuelven arcos.';

  @override
  String get treeArc => 'Apertura';

  @override
  String get treeArcInfo =>
      'Cuánto del círculo ocupan las puntas: 360° lo cierra, 180° forma un semicírculo (abanico).';

  @override
  String get treeRotation => 'Rotación';

  @override
  String get treeSwapSides => 'Cambiar lados';

  @override
  String get treeSwapSidesText => 'Este árbol a la derecha';

  @override
  String get treeSwapSidesInfo =>
      'Pone el árbol en edición a la derecha y el otro a la izquierda.';

  @override
  String get iqAlrt => 'Prueba SH-aLRT';

  @override
  String get iqAlrtInfo =>
      'Un segundo valor de soporte, rápido, para cada rama (se muestra como SH-aLRT/UFBoot). Las ramas con SH-aLRT ≥ 80 y UFBoot ≥ 95 están bien sustentadas.';

  @override
  String get iqAlrt1000 => '1000 réplicas (recomendado)';

  @override
  String get iqSearchHeading => 'Búsqueda del árbol';

  @override
  String get iqSearchInfo =>
      'Los valores por defecto sirven para la mayoría de los datos. Más ejecuciones o más iteraciones tardan más, pero exploran más árboles; útil cuando el resultado cambia entre ejecuciones.';

  @override
  String get iqRuns => 'Ejecuciones independientes';

  @override
  String get iqRunsInfo =>
      '--runs: repite toda la búsqueda desde árboles iniciales distintos y conserva el mejor árbol.';

  @override
  String get iqNstop => 'Parar tras (iteraciones)';

  @override
  String get iqNstopInfo =>
      '-nstop: la búsqueda se detiene tras este número de iteraciones sin un árbol mejor (por defecto 100). Use más para datos grandes o difíciles.';

  @override
  String get iqPers => 'Fuerza de perturbación';

  @override
  String get iqPersInfo =>
      '-pers: cuánto se modifica el árbol actual en cada iteración (por defecto 0,5). Valores menores ayudan con alineamientos cortos.';

  @override
  String get iqSeedInfo =>
      'La misma semilla repite el mismo resultado. Vacío: semilla aleatoria (anotada en el registro).';

  @override
  String get iqSeedRandom => 'aleatoria';

  @override
  String get iqBnni => 'Refinar los árboles del UFBoot (--bnni)';

  @override
  String get iqBnniBody =>
      'Optimiza cada árbol de bootstrap con movimientos NNI. Recomendado: reduce soportes sobreestimados cuando el modelo no se ajusta perfectamente a los datos.';

  @override
  String get bothTitle => 'Ejecutar IQ-TREE y RAxML-NG juntos';

  @override
  String get bothInfo =>
      'Ambos son métodos de máxima verosimilitud. Comparar los árboles (paso 9, Comparar árboles) muestra qué clados no dependen del programa.';

  @override
  String bothBody(String iq, String rx) {
    return 'Ambos análisis empiezan a la vez con las opciones elegidas en las dos tarjetas de arriba. Los hilos se dividen para que el equipo siga usable: IQ-TREE recibe $iq y RAxML-NG $rx. El total se cambia en la configuración (barra superior).';
  }

  @override
  String bothThreads(String n) {
    return '$n hilos';
  }

  @override
  String get runBoth => 'Ejecutar ambos';

  @override
  String get bothRunning => 'IQ-TREE y RAxML-NG';

  @override
  String bothRunningSub(String threads, String time) {
    return 'En paralelo · hilos $threads · $time';
  }

  @override
  String get satUnlocksStep7 =>
      'Ejecute la prueba para continuar: su sugerencia (qué posiciones del codón conservar) se usa para construir la matriz en el paso 7, que se abre a continuación.';

  @override
  String get satAllGenes => 'Todos los genes juntos';

  @override
  String get satAllGenesInfo =>
      '\"Todos los genes juntos\" une los genes codificantes solo para esta prueba (la matriz se construye en el paso 7). Es la fila que decide la sugerencia de qué posiciones del codón quitar de toda la matriz; la fila de cada gen permite elegir gen por gen.';

  @override
  String get exampleLabel => 'Ejemplo';

  @override
  String get exampleNt =>
      'El mismo alineamiento en los tres modos. Aquí el gap de la especie B cae dentro de los codones 3 y 4 y desplaza el marco de lectura.';

  @override
  String get exampleAa =>
      'Los mismos genes traducidos: una letra por codón, así que menos sitios y más conservados.';

  @override
  String get exampleCodon =>
      'El gap de la especie B es el codón 3 entero: se conserva el marco de lectura.';

  @override
  String exampleTrim(String n, String total) {
    return 'El recorte elimina columnas poco fiables (casi solo gaps, o sin señal filogenética) antes de construir el árbol. Aquí se eliminan $n de $total columnas (×).';
  }

  @override
  String get exampleCodonKeepAll =>
      'Se conservan las tres posiciones del codón.';

  @override
  String exampleCodonKept(String positions, String pct) {
    return 'Se conservan: posiciones $positions ($pct% de los sitios codificantes).';
  }

  @override
  String get exampleSatLabel => 'Cómo leer los gráficos';

  @override
  String get exampleSatOk => 'Sin saturación';

  @override
  String get exampleSatOkBody =>
      'Las diferencias observadas (p) siguen subiendo con la distancia (d), cerca de la línea discontinua y = x.';

  @override
  String get exampleSatBad => 'Saturado';

  @override
  String get exampleSatBadBody =>
      'p se estabiliza mientras d sigue creciendo: los cambios nuevos borran los antiguos y la señal se pierde.';

  @override
  String get plotTrend => 'tendencia';

  @override
  String exampleSpecies(String letter) {
    return 'Especie $letter';
  }

  @override
  String get examplePlayAgain => 'Ver de nuevo';

  @override
  String get trimStepOriginal => 'Alineamiento';

  @override
  String get trimStepMarked => 'Columnas eliminadas';

  @override
  String get trimStepResult => 'Alineamiento recortado';

  @override
  String get treeSupportArrowLabel => 'Flechas';

  @override
  String get treeSupportArrow => 'Valor con flecha a la rama';

  @override
  String get treeSupportArrowInfo =>
      'Aleja cada valor de soporte un poco de las ramas (arriba o abajo, como en Posición del soporte) y lo une a su nodo con una pequeña flecha. Útil cuando los valores caen sobre ramas cortas o encima de otras líneas.';

  @override
  String get orgOpen => 'Muestra y marcador';

  @override
  String get orgTitle => 'Muestra y marcador en los encabezados';

  @override
  String orgSubtitle(String files, String seqs) {
    return '$files archivo(s) FASTA · $seqs secuencias';
  }

  @override
  String orgSummary(String samples, String markers, String pct) {
    return '$samples muestras × $markers marcadores · $pct% completo';
  }

  @override
  String orgApply(String n) {
    return 'Generar $n FASTA de marcador';
  }

  @override
  String get orgIntro =>
      'Cada encabezado debe indicar a qué muestra (ID) y a qué marcador pertenece la secuencia. Elija los caracteres que separan las partes del encabezado y haga clic en las partes de abajo: una vez para muestra, dos para marcador, una tercera para ignorarla. SPLACE escribe entonces un FASTA por marcador con todas las muestras, listo para alinear.';

  @override
  String get orgSampleFrom => 'Muestra (ID) desde';

  @override
  String get orgSampleFromInfo =>
      'Una parte del encabezado, o el nombre del archivo cuando cada archivo tiene una muestra y los encabezados solo nombran el marcador.';

  @override
  String get orgSampleHeader => 'Encabezado';

  @override
  String get orgSampleFile => 'Nombre del archivo';

  @override
  String get orgDelimitersInfo =>
      'Cada carácter activo separa el encabezado. Puede haber varios a la vez cuando los archivos usan separadores distintos (ITV01-L1 e ITV02_L1).';

  @override
  String get orgSpace => 'espacio';

  @override
  String get orgOther => 'otros';

  @override
  String get orgParts =>
      'Partes del encabezado (ejemplos de archivos distintos)';

  @override
  String get orgPartsInfo =>
      'Las posiciones valen para todos los encabezados. Una parte que llega al final del encabezado se cuenta desde el final, así que los encabezados con más o menos partes también funcionan.';

  @override
  String get orgLegendSample => 'Muestra (ID)';

  @override
  String get orgLegendMarker => 'Marcador';

  @override
  String get orgLegendIgnored => 'Ignorado';

  @override
  String get orgAdvanced => 'Avanzado: etiquetas y expresiones regulares';

  @override
  String get orgAdvancedActive =>
      'Hay una etiqueta o expresión regular en uso (vea Avanzado); el clic en las partes está desactivado.';

  @override
  String orgImplausible(String markers, String seqs) {
    return '$markers marcadores para $seqs secuencias';
  }

  @override
  String get orgImplausibleBody =>
      'Casi cada secuencia se volvió su propio marcador, así que la parte del marcador probablemente está mal (por ejemplo, el encabezado entero). Elija el delimitador y haga clic en la parte que nombra el marcador.';

  @override
  String orgUnassigned(String n) {
    return '$n secuencia(s) sin muestra o sin marcador (quedan fuera)';
  }

  @override
  String get orgUnassignedLabel => 'Sin muestra o marcador';

  @override
  String get orgVariations => 'Variaciones de nombre';

  @override
  String get orgVariationsInfo =>
      'La misma muestra o marcador puede estar escrita de formas un poco distintas. Estas opciones las unen; también puede renombrar cualquier muestra en la lista de abajo.';

  @override
  String get orgIgnoreCase => 'Ignorar mayúsculas y minúsculas';

  @override
  String get orgLeadingZeros => 'Ignorar ceros a la izquierda en los números';

  @override
  String orgFileMismatch(String n) {
    return '$n archivo(s) cuyo nombre difiere de la muestra en sus encabezados';
  }

  @override
  String orgRenameTo(String name) {
    return 'Llamarla $name';
  }

  @override
  String get orgRenameAll => 'Usar el nombre del archivo en estas';

  @override
  String get orgUseFileNames =>
      'Tomar todas las muestras del nombre del archivo';

  @override
  String orgMerged(String n) {
    return '$n muestra(s) escrita(s) de más de una forma, unidas:';
  }

  @override
  String orgSamples(String n) {
    return 'Muestras ($n)';
  }

  @override
  String get orgSamplesInfo =>
      'Haga clic en un nombre para renombrarlo (Enter aplica). Dar el mismo nombre a dos muestras las une.';

  @override
  String get orgFindSample => 'Buscar muestra';

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
      'Cuántas muestras tiene cada marcador, y qué hacer cuando una muestra tiene más de una secuencia del mismo marcador. Todos los marcadores pasan al paso 4, donde elige cuáles usar.';

  @override
  String orgDuplicates(String n) {
    return '$n caso(s) de muestra con más de una secuencia del mismo marcador';
  }

  @override
  String get orgCoverage => 'Marcadores por número de muestras';

  @override
  String orgCoverageBar(String markers, String samples) {
    return '$markers marcadores en $samples muestras';
  }

  @override
  String get orgSamplesAxis => 'muestras';

  @override
  String orgCoverageSummary(String complete, String half) {
    return '$complete en todas las muestras · $half en al menos la mitad';
  }

  @override
  String get orgPending =>
      'Indique a SPLACE qué parte de los encabezados es la muestra y cuál el marcador';

  @override
  String orgProblem(String problem) {
    return 'La división no se usa: $problem';
  }

  @override
  String get orgPendingBody =>
      'Después se escribe un FASTA por marcador con todas las muestras, y los marcadores aparecen en el paso 4.';

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
    return 'Mostrar $n más';
  }

  @override
  String get showAll => 'Mostrar todos';

  @override
  String get noMatches => 'Nada coincide.';

  @override
  String orgFormats(String n) {
    return '$n formatos de encabezado';
  }

  @override
  String perLocusCount(String positions, String n) {
    return '$positions en $n marcadores';
  }

  @override
  String get satOptionalLoci =>
      'Opcional para loci: el paso 7 ya está disponible. Los marcadores marcados como codificantes se prueban por posición del codón; los demás, en todos sus sitios.';

  @override
  String get satOptionalLociNoCoding =>
      'Opcional para loci: el paso 7 ya está disponible. Ningún marcador es codificante, así que cada locus se prueba en todos sus sitios (no hay posiciones de codón): la prueba muestra qué loci están saturados.';

  @override
  String exampleTrimNone(String total) {
    return 'No se elimina ninguna de las $total columnas: un alineamiento de MAFFT no tiene columnas formadas solo por gaps. Esta opción solo importa en alineamientos editados a mano o en subconjuntos de secuencias.';
  }

  @override
  String get exampleGapsRow => 'gaps';

  @override
  String get runFailed => 'El análisis falló.';

  @override
  String get runLastLines => 'Últimas líneas del programa';

  @override
  String get runOpenLog => 'Abrir el log';

  @override
  String orgHeaderFormats(String n) {
    return 'Formatos de encabezado ($n)';
  }

  @override
  String get orgHeaderFormatsInfo =>
      'Los encabezados están escritos de más de una forma (otros delimitadores, más o menos partes). Cada formato muestra sus secuencias y puede tener sus propias partes de muestra y marcador: haga clic en las partes de su ejemplo. Al final, todas las secuencias quedan con el mismo patrón: la muestra en el encabezado y un FASTA por marcador.';

  @override
  String orgFormatCount(String n, String files) {
    return '$n secuencias en $files archivos';
  }

  @override
  String get orgFormatOwnRule => 'Partes propias';

  @override
  String get orgFormatSkipped => 'Excluido';

  @override
  String orgFormatShow(String n) {
    return 'Mostrar las $n secuencias';
  }

  @override
  String get orgFormatHide => 'Ocultar las secuencias';

  @override
  String get orgFormatUseMain => 'Mismas partes que la regla principal';

  @override
  String get orgFormatSkip => 'Excluir estas secuencias';

  @override
  String get orgFormatUse => 'Usar estas secuencias';

  @override
  String orgNameStructures(String n) {
    return 'Estructuras de los nombres de muestra ($n)';
  }

  @override
  String get orgNameStructuresInfo =>
      'Nombres de muestra agrupados según cómo están escritos: letras tal cual, cada dígito como #. Los nombres de una colección escritos con más o menos dígitos (ITV1054 e ITV01054) pueden estandarizarse a la estructura que elija: se agregan o quitan ceros a la izquierda. Pase el ratón sobre el botón para ver los cambios de nombre.';

  @override
  String orgStructureCount(int n) {
    String _temp0 = intl.Intl.pluralLogic(
      n,
      locale: localeName,
      other: '$n muestras',
      one: '1 muestra',
    );
    return '$_temp0';
  }

  @override
  String orgStandardize(String structure) {
    return 'Estandarizar como $structure';
  }

  @override
  String get satAllSitesNonCoding => 'todos los sitios (no codificante)';

  @override
  String schemeMixed(String coding, String other) {
    return '$coding marcadores codificantes divididos por codón; $other sin codones quedan enteros';
  }

  @override
  String get exampleCodonGene => 'Gen';

  @override
  String get exampleCodonMatrix => 'En la matriz';

  @override
  String get mtTitle => 'Taxones × marcadores';

  @override
  String get mtIntro =>
      'Como en SequenceMatrix: cada celda es el número de bases que el taxón tiene en ese marcador (verde más oscuro = más completo), \"—\" cuando no tiene secuencia. El encabezado muestra la longitud de cada alineamiento, cuántos taxones lo tienen y sus gaps. Pase el ratón sobre una celda o encabezado para ver detalles; la misma tabla se guarda en matrix_overview.tsv.';

  @override
  String get mtTaxon => 'Taxón';

  @override
  String get mtLoci => 'Marcadores';

  @override
  String get mtBases => 'Bases';

  @override
  String get mtMissing => 'Faltante';

  @override
  String get mtPartitions => 'Particiones';

  @override
  String get mtLegend =>
      'Celda: bases del taxón en el marcador · — sin secuencia · encabezado: longitud · taxones con el marcador · % de gaps';

  @override
  String mtLocusInfo(
    String locus,
    String length,
    String n,
    String total,
    String gaps,
  ) {
    return '$locus: $length pb, en $n de $total taxones, $gaps% de gaps';
  }

  @override
  String get mtMissingCell => 'sin secuencia (rellenado con ?)';

  @override
  String mtCellInfo(String bases, String length, String gaps) {
    return '$bases de $length pb, $gaps gaps';
  }

  @override
  String get searchGenomes => 'Genomas';

  @override
  String get searchMarkers => 'Marcadores individuales';

  @override
  String get markerIncludeGenomes => 'Incluir genomas';

  @override
  String get markerGenomesInfo =>
      'Apagado: solo registros de hasta 10 kb (secuencias de genes, como barcodes). Encendido: también genomas completos y fragmentos largos, que también contienen los marcadores (y son mucho más grandes de descargar).';

  @override
  String markerPick(String n) {
    return 'Marcadores (SynGenes) · $n seleccionados';
  }

  @override
  String get markerPickInfo =>
      'Genes codificantes y rRNAs de las tablas de SynGenes. Cada marcador se busca con todos los nombres que SynGenes conoce (símbolos y nombres del producto, p. ej. COI, COX1, cytochrome c oxidase subunit I), así que se encuentran registros anotados de formas distintas. Los registros de la misma especie forman una sola muestra, y las especies sin un marcador quedan con datos faltantes en él.';

  @override
  String reviewMarkerCount(
    String marker,
    String selected,
    String found,
    String species,
  ) {
    return '$marker: $selected/$found registros · $species especies';
  }

  @override
  String get reviewOnePerSpecies => 'Uno por especie y marcador';

  @override
  String get reviewOnePerSpeciesInfo =>
      'Conserva, para cada marcador, el registro más largo de cada especie (registros verificados primero). Es la selección con la que empieza la revisión.';

  @override
  String get reviewColMarkers => 'Marcadores';

  @override
  String get onePerTaxon => 'Una muestra por especie (marcadores individuales)';

  @override
  String get onePerTaxonBody =>
      'Une los registros de la misma especie (un COI de un voucher, un CYTB de otro) en una fila de la matriz: el encabezado es familia_especie y, para cada marcador, se conserva la secuencia más larga de la especie. Las especies sin un marcador quedan con datos faltantes (?).';

  @override
  String onePerTaxonJoined(String n) {
    return '$n especies tienen varios registros: cada una forma una sola muestra.';
  }

  @override
  String satOptionalPartial(String markers) {
    return 'Opcional aquí: el paso 7 ya está disponible. $markers vienen de registros de marcadores individuales con secuencias parciales (fragmentos como barcodes, marcados con < o > en GenBank), así que la prueba es una guía, no un requisito.';
  }

  @override
  String satPartialMixed(String markers) {
    return '$markers tienen secuencias parciales: para ellos la prueba es opcional. El paso 7 espera solo a los marcadores completos.';
  }

  @override
  String metaRecords(String n) {
    return '$n registros';
  }

  @override
  String headerConflictBlock(String n) {
    return '$n encabezado(s) se repiten dentro del mismo marcador, así que esas secuencias no podrían distinguirse en el alineamiento ni en el árbol. Agregue un campo que las distinga (Accesión o Voucher) o active \"Una muestra por especie\" para unirlas. El alineamiento queda bloqueado hasta que todos los encabezados sean únicos.';
  }

  @override
  String get headerAddAccession => 'Agregar la Accesión';

  @override
  String get alignBlockedHeader =>
      'Bloqueado: encabezados FASTA repetidos (5a)';

  @override
  String recGenesCount(String n, String total) {
    return 'Saturada en $n de $total genes.';
  }

  @override
  String recBasedOnGene(String gene, String taxa, String length) {
    return 'Resultado de $gene ($taxa taxones, $length pb). Cada posición se prueba por separado.';
  }

  @override
  String recBasedOnGenes(String n) {
    return '$n genes codificantes, cada uno probado por separado (nada se concatena para la prueba). Una posición cuenta como saturada cuando está saturada en más de la mitad de los genes.';
  }

  @override
  String rec7bPreview(String section) {
    return 'Cómo quedará \"$section\" en el paso 7b:';
  }

  @override
  String get frameExampleLabel => 'Ejemplo · mantener el marco de lectura';

  @override
  String get frameExampleCaption =>
      'trimAl elimina columnas sueltas. Con esta opción, cuando trimAl elimina una columna (×), SPLACE elimina el codón entero al que pertenece (+): el gen conserva codones enteros y las posiciones 1, 2 y 3 siguen en orden (la regla de colores).';

  @override
  String get frameStepOriginal => 'Codones 1-2-3';

  @override
  String get frameStepColumn => 'trimAl elimina una columna';

  @override
  String get frameStepCodon => 'ClipKIT --codon quita el codón entero';

  @override
  String get frameStepResult => 'Marco mantenido';

  @override
  String get trimCasesLabel =>
      'Ejemplo · recortar un gen codificante, de dos formas';

  @override
  String get trimCasesStepMark => 'Columnas marcadas';

  @override
  String get trimCaseColumns => 'Columna a columna (trimAl)';

  @override
  String get trimCaseCodons => 'Mantener el marco de lectura (ClipKIT --codon)';

  @override
  String get trimCaseSelected => 'seleccionado';

  @override
  String trimCaseFrameLost(String n) {
    return 'Marco perdido: $n columnas fuera de la posición del codón';
  }

  @override
  String trimCaseFrameKept(String k, String total) {
    return 'Marco mantenido: $k de $total codones, posiciones 1-2-3 en orden';
  }

  @override
  String get trimCaseFrameKeptCols =>
      'Marco mantenido: ninguna columna fuera de posición';

  @override
  String get trimCasesCaption =>
      'Mismo alineamiento, dos herramientas. trimAl elimina una a una las columnas que marca (×), así que las posiciones del codón (regla de colores) quedan fuera de orden (!). ClipKIT con --codon elimina el codón entero de cada columna que marca (+): las posiciones 1.ª, 2.ª y 3.ª siguen en orden para las particiones y la prueba de saturación.';

  @override
  String get subRecordsTable => 'Tabla de registros';

  @override
  String get subRecordsTableBody =>
      'Cada registro con su taxonomía y conteo de genes; edite o elimine registros aquí.';

  @override
  String get subFeatureTypes => 'Anotaciones usadas (CDS, rRNA, tRNA)';

  @override
  String get subChooseMarkers => 'Elegir los marcadores';

  @override
  String get subRunSaturation => 'Ejecutar la prueba';

  @override
  String get subBuildMatrix => 'Construir la matriz';

  @override
  String get subTree => 'Árbol y figura';

  @override
  String get viewCodons => 'Codones';

  @override
  String get viewDots => 'Puntos';

  @override
  String get viewDotsInfo =>
      'Muestra un punto en lugar de las bases en las columnas donde todas las secuencias tienen la misma base (la primera secuencia conserva sus letras), para destacar las columnas variables, como en MEGA y otros visores. En la vista por codones, los codones iguales en todas las secuencias se vuelven puntos.';

  @override
  String codonsVarying(String n, String total, String nonsyn) {
    return '$n de $total codones varían · $nonsyn cambian el aminoácido';
  }

  @override
  String get codonVariation => 'Variación del codón';

  @override
  String viewerCodonN(String n) {
    return 'codón $n';
  }

  @override
  String codonColumnInfo(String codons, String aa) {
    return '$codons codón(es) distinto(s) aquí, $aa aminoácido(s)';
  }

  @override
  String get codonLegendAa => 'color = aminoácido';

  @override
  String get codonLegendSyn => 'cambio sinónimo';

  @override
  String get codonLegendNonSyn => 'cambia el aminoácido';

  @override
  String get mtResizeNames =>
      'Arrastre para ensanchar los nombres; doble clic para ajustar al más largo';

  @override
  String treeLegendUnmarked(String min) {
    return 'sin marca: soporte < $min';
  }

  @override
  String get treeStrongFrom => 'Soporte fuerte desde';

  @override
  String get treeModerateFrom => 'Soporte moderado desde';

  @override
  String get treeSupportClassesInfo =>
      'Cada clase colorea los puntos (y valores) desde su umbral hasta el siguiente. Elija un modelo — automático según el tipo de soporte, de 2 a 5 clases, o cortado en los terciles de los valores de este árbol — o defina usted mismo los umbrales y colores, añadiendo o quitando clases.';

  @override
  String get treeLegendSupport => 'Soporte de las ramas';

  @override
  String get treeLegendSupportShow => 'En la leyenda';

  @override
  String get treeLegendSupportTitle => 'Título del soporte';

  @override
  String get badgeBrazil => 'solo Brasil';

  @override
  String get salvePossiblyExtinct => 'Posiblemente extinta';

  @override
  String get salveSheetSubtitle =>
      'Ficha de SALVE (ICMBio) · evaluación nacional del riesgo de extinción';

  @override
  String get salveOpenSheet => 'Abrir la ficha en SALVE';

  @override
  String get salveWatersheds => 'Cuencas';

  @override
  String get salvePublished => 'Publicada';

  @override
  String get salveSituation => 'Estado de la ficha';

  @override
  String get salveSynonyms => 'Nombres anteriores';

  @override
  String get salveJustification => 'Justificación';

  @override
  String get salveAssessors => 'Evaluadores';

  @override
  String get salveCitation => 'Cómo citar';

  @override
  String get colCriteria => 'Criterios';

  @override
  String get colCommonNames => 'Nombres comunes';

  @override
  String get colSalveGroup => 'Grupo SALVE';

  @override
  String get colEndemic => 'Endémica de Brasil';

  @override
  String get colBiomes => 'Biomas';

  @override
  String get colStates => 'Estados';

  @override
  String get colAssessed => 'Evaluada';

  @override
  String get colSheet => 'Ficha';

  @override
  String get geneNamesAllKnown => 'Todos los nombres de genes reconocidos';

  @override
  String geneReviewAgain(String n) {
    return 'Revisar de nuevo ($n)';
  }

  @override
  String get clipkitCodonMode =>
      'Marcadores codificantes: modo de ClipKIT (--codon)';

  @override
  String get clipkitCodonInfo =>
      'Para mantener el marco de lectura, los marcadores codificantes se recortan con ClipKIT y su opción --codon: cuando una posición de un codón debe salir, ClipKIT elimina el codón entero (Steenwyk et al. 2020). Los marcadores no codificantes se recortan con trimAl con el método de arriba.';

  @override
  String get outgroupSuggesting =>
      'Midiendo las distancias entre los taxones de los alineamientos para sugerir un grupo externo…';

  @override
  String get outgroupSuggestNone =>
      'Alinee los marcadores para recibir una sugerencia de grupo externo por las distancias.';

  @override
  String get outgroupSuggestTitle =>
      'Grupo externo sugerido (distancias de los alineamientos)';

  @override
  String get outgroupSuggestWeak =>
      'Ningún grupo claramente separado: se sugiere el taxón más distante';

  @override
  String get outgroupSuggestInfo =>
      'Distancias p entre cada par de taxones en los sitios que ambos tienen en cada alineamiento, sumadas entre marcadores; un árbol de ligamiento promedio (UPGMA) de ellas; el menor de los dos grupos unidos al final es el más alejado del resto. El valor junto a cada taxón es su distancia media a todos los demás. Contrástelo con lo que sabe del grupo: una rama larga (secuencia de evolución rápida o contaminada) también puede parecer distante.';

  @override
  String get outgroupSuggestApplied => 'En uso';

  @override
  String get outgroupSuggestUse => 'Usar como grupo externo';

  @override
  String outgroupSuggestStats(
    String between,
    String within,
    String loci,
    String taxa,
  ) {
    return 'Distancia p media al grupo interno $between; dentro del grupo interno $within · $loci alineamientos, $taxa taxones.';
  }

  @override
  String get outgroupSuggestNext => 'Siguientes más distantes';

  @override
  String get treeArrowShape => 'Forma de la flecha';

  @override
  String get treeArrowStraight => 'Recta';

  @override
  String get treeArrowCurved => 'Curva';

  @override
  String get treeArrowShapeInfo =>
      'Los valores nunca se superponen: cuando el lugar de un valor está ocupado se aleja y una flecha apunta a su nodo. Las flechas pueden ser rectas o curvas.';

  @override
  String get treeSupportClasses => 'Clases de soporte';

  @override
  String treeClassesAuto(String kind) {
    return 'Automático ($kind)';
  }

  @override
  String treeClassesN(String n) {
    return '$n clases';
  }

  @override
  String get treeClassesTertiles => 'Terciles de este árbol';

  @override
  String get treeClassHidden => 'por debajo de \"mostrar soporte desde\"';

  @override
  String get treeClassRemove => 'Quitar esta clase';

  @override
  String get treeClassAdd => 'Añadir clase';

  @override
  String get tipFieldConservation => 'Categoría (SALVE)';

  @override
  String get treeHeatStatus => 'Categoría';

  @override
  String get treeHeatStatusLegend => 'Riesgo de extinción (SALVE / ICMBio)';

  @override
  String get treeHeatYes => 'sí';

  @override
  String get treeHeatNo => 'no';

  @override
  String get treeHeatmap => 'Heatmap de SALVE';

  @override
  String get treeHeatmapInfo =>
      'Columnas de colores junto a las etiquetas a partir de las fichas de SALVE de las especies: categoría en los colores de la UICN, endemismo en Brasil y una columna por bioma. Solo en el diseño rectangular. La misma información puede ir en las etiquetas como texto (partes de la etiqueta arriba).';

  @override
  String get satParamsTitle => 'Parámetros de la prueba';

  @override
  String get satParamsBody =>
      'Cómo se mide la saturación de cada marcador y cuándo se considera moderada o saturada. Los cambios valen en la próxima ejecución.';

  @override
  String get satParamMethod => 'Método';

  @override
  String get satParamMethodValue =>
      'saturación de PhyKIT (Steenwyk et al. 2021; Philippe et al. 2011): pendiente de la distancia p no corregida de cada par de taxones sobre su distancia patrística en el árbol';

  @override
  String get satParamDistance => 'Árboles';

  @override
  String get satParamDistanceValue =>
      'IQ-TREE, uno por gen (búsqueda rápida) con el modelo de abajo; las longitudes de rama de cada subconjunto de posiciones se estiman sobre esa topología (-te)';

  @override
  String get satParamSubsets => 'Sitios probados';

  @override
  String get satParamSubsetsValue =>
      '1.ª, 2.ª y 3.ª posiciones del codón, 1.ª+2.ª y todos los sitios (marcadores no codificantes: todos los sitios)';

  @override
  String get satParamCode => 'Código genético';

  @override
  String get satParamFrame => 'Marco de lectura';

  @override
  String get satParamFrameValue =>
      'del alineamiento (inicio del CDS, o el marco con menos stops internos)';

  @override
  String get satParamsThresholds => 'Umbrales y ajustes';

  @override
  String get satParamSlopeModerate => 'Saturación · moderada <';

  @override
  String get satParamSlopeSaturated => 'Saturación · saturada <';

  @override
  String get satParamSlopeInfo =>
      'Valor de saturación de PhyKIT: pendiente de la distancia p sobre la distancia patrística. 1 = sin saturación; cuanto más cerca de 0, más saturada. Los valores por defecto (0,5 y 0,3) son heurísticos: las 2.ª posiciones del codón suelen quedar cerca de 1 y las 3.ª saturadas muy por debajo de 0,5.';

  @override
  String get satParamIss => 'Iss · saturada ≥';

  @override
  String get satParamIssInfo =>
      'Índice de Xia: entropía media de los sitios sobre la entropía esperada bajo saturación total. Cerca de 1, los sitios parecen aleatorios.';

  @override
  String get satParamUndefinedModerate => 'Sin corrección · moderada ≥';

  @override
  String get satParamUndefinedSaturated => 'Sin corrección · saturada ≥';

  @override
  String get satParamUndefinedInfo =>
      'Pares cuya distancia corregida no se puede calcular (más allá de lo que corrige el modelo): fuerte señal de saturación.';

  @override
  String get satParamLowDivergence => 'Divergencia baja · p <';

  @override
  String get satParamLowDivergenceInfo =>
      'Con una distancia p media tan baja, la saturación no importa y el subconjunto se considera OK.';

  @override
  String get satParamMinSites => 'Mínimo de sitios';

  @override
  String get satParamMinPairs => 'Mínimo de pares';

  @override
  String get satParamReplicates => 'Iss · réplicas';

  @override
  String get satParamReplicatesInfo =>
      'Sitios aleatorios sorteados para estimar la entropía esperada bajo saturación total, para cada número de secuencias. Más réplicas, Iss más estable y ejecución más larga.';

  @override
  String get satParamPlotPoints => 'Puntos del gráfico';

  @override
  String get satParamPlotPointsInfo =>
      'Pares de secuencias dibujados en los gráficos (una muestra cuando hay más); las estadísticas usan todos los pares.';

  @override
  String get satParamsReset => 'Restaurar valores por defecto';

  @override
  String get satParamsStale =>
      'Los parámetros cambiaron desde la última ejecución: vuelva a ejecutar la prueba para actualizar los resultados.';

  @override
  String get colTrend => 'Tendencia poblacional';

  @override
  String get colGeneration => 'Tiempo generacional';

  @override
  String get colMigratory => 'Migratoria';

  @override
  String get colHabitats => 'Hábitats';

  @override
  String get colThreats => 'Amenazas';

  @override
  String get colOfficialList => 'Lista oficial';

  @override
  String get colUcs => 'Áreas protegidas';

  @override
  String get salveReferences => 'Referencias en la ficha';

  @override
  String get salveHistory => 'Historial de evaluaciones';

  @override
  String get salveLists => 'Listas oficiales y convenciones';

  @override
  String get salveSecDistribution => 'Distribución';

  @override
  String get salveNational => 'En Brasil';

  @override
  String get salveGlobal => 'En el mundo';

  @override
  String get salveSecNatural => 'Historia natural y población';

  @override
  String get salveDiet => 'Hábito alimentario';

  @override
  String get salveReproduction => 'Reproducción';

  @override
  String get salveSexRatio => 'Razón de sexos';

  @override
  String get salveNaturalHistory => 'Historia natural';

  @override
  String get salveInteractions => 'Interacciones con otras especies';

  @override
  String get salvePopulation => 'Población';

  @override
  String get salveGenetics => 'Características genéticas';

  @override
  String get salveSecThreats => 'Amenazas, usos, conservación e investigación';

  @override
  String get salveUses => 'Usos';

  @override
  String get salveActions => 'Acciones de conservación';

  @override
  String get salveConservation => 'Conservación';

  @override
  String get salveResearch => 'Investigación';

  @override
  String get trendDeclining => 'en declive';

  @override
  String get trendStable => 'estable';

  @override
  String get trendIncreasing => 'en aumento';

  @override
  String get trendUnknown => 'desconocida';

  @override
  String get treeHeatFew => 'pocas';

  @override
  String get treeHeatMany => 'muchas';

  @override
  String trimNothingRemoved(String method, String n) {
    return '$method mantuvo las $n columnas: este alineamiento no tiene regiones con gaps para quitar, así que nada queda en gris (los modos kpic también quitan sitios no informativos)';
  }

  @override
  String get treeBranchWidth => 'Ancho de las ramas';

  @override
  String get treeBranchWidthInfo =>
      'Ancho del árbol en sí (ramas), en cm. 0 = automático: lo que queda junto a etiquetas, heatmap, barras y leyenda, nunca menos del 45 % de la figura — con muchos elementos la figura se ensancha en lugar de comprimir el árbol. Fije un ancho para mantenerlo; la figura crece para que quepa el resto.';

  @override
  String get treeHeatCell => 'Celda del heatmap';

  @override
  String get treeHeatCellInfo =>
      'Lado de cada celda del heatmap, en pt. 0 = automático (según la altura de fila).';

  @override
  String get citeVersions => 'Versiones usadas';

  @override
  String get satParamVersions => 'Programas';

  @override
  String get satParamModel => 'Modelo de IQ-TREE';

  @override
  String get satParamModelInfo =>
      'Modelo de los árboles de los genes y de las longitudes de rama de cada subconjunto (IQ-TREE -m), p. ej. GTR+G, HKY+G, GTR+F+I+G4.';

  @override
  String get satParamGaps => 'Gaps';

  @override
  String get satParamGapsSwitch => 'Ignorar sitios con gaps (-e)';

  @override
  String get satParamGapsInfo =>
      'PhyKIT --exclude_gaps: los sitios con un gap en un par quedan fuera de su distancia p.';

  @override
  String get phykitRequired =>
      'La prueba de saturación usa PhyKIT, que aún no está instalado. Instálelo aquí (requiere Python 3) o con: pip install phykit';

  @override
  String get installPhykit => 'Instalar PhyKIT';

  @override
  String plotPhykit(String s) {
    return 'p × distancia en el árbol · PhyKIT $s';
  }

  @override
  String get plotPhykitBody =>
      'Cada punto es un par de taxones: distancia p no corregida frente a su distancia patrística en el árbol de IQ-TREE de estos sitios. Los puntos muy por debajo de la diagonal esconden sustituciones repetidas.';

  @override
  String get plotPositionsPatristic =>
      'Posiciones del codón × distancia en el árbol';

  @override
  String get patristicDistance => 'Distancia patrística';

  @override
  String get metricMeanPatristic => 'Patrística media';

  @override
  String get metricPhykit => 'Saturación PhyKIT';

  @override
  String get trimToolLabel => 'Herramienta de recorte';

  @override
  String get trimToolTrimal => 'trimAl · columna a columna';

  @override
  String get trimToolClipkit => 'ClipKIT · codones enteros';

  @override
  String get trimToolTrimalBody =>
      'trimAl recorta todos los marcadores columna a columna. En genes codificantes esto puede partir codones: se pierde el marco de lectura, pasan a ser una partición sin posiciones del codón y quedan fuera de la prueba de saturación.';

  @override
  String get trimToolClipkitBody =>
      'ClipKIT recorta todos los marcadores. En los genes codificantes se ejecuta con --codon: cuando una posición de un codón debe salir, sale el codón entero, así se mantienen el marco de lectura y las posiciones del codón. Los marcadores no codificantes se recortan con ClipKIT sin --codon. trimAl no se usa.';

  @override
  String get mtShowBases => 'Bases';

  @override
  String get mtShowPercent => '% del marcador';

  @override
  String get mtOrderMatrix => 'Orden de la matriz';

  @override
  String get mtOrderComplete => 'Más completos primero';

  @override
  String get mtScaleComplete => 'Completitud';

  @override
  String get mtScaleMissing => 'Faltantes / gaps';

  @override
  String get treeMarkShape => 'Forma de los puntos';

  @override
  String get treeShapeCircle => 'Círculo';

  @override
  String get treeShapeSquare => 'Cuadrado';

  @override
  String get treeShapeDiamond => 'Rombo';

  @override
  String get treeShapeTriangle => 'Triángulo';

  @override
  String get treeBarWidth => 'Grosor de las barras de clado';

  @override
  String get treeBarWidthInfo =>
      'Grosor de las barras junto a los clados, en pt. 0 = automático (según el grosor de las ramas).';

  @override
  String get treeBarFont => 'Tamaño de los nombres de clado';

  @override
  String get treeBarFontInfo =>
      'Tamaño de los nombres junto a las barras de clado, en pt. 0 = el tamaño de las etiquetas de las puntas.';

  @override
  String get treeHeatGap => 'Espacio entre columnas del heatmap';

  @override
  String get trimToolTrimalCodon => 'trimAl · codones enteros';

  @override
  String get trimToolTrimalCodonBody =>
      'trimAl recorta todos los marcadores. Los genes codificantes se leen codón a codón y se traducen con el código genético del proyecto; trimAl recorta ese alineamiento de proteínas y, con -backtrans, devuelve los codones de los aminoácidos conservados, así que salen o quedan codones enteros y se mantienen las posiciones del codón. Los marcadores no codificantes se recortan columna a columna.';

  @override
  String frameStepCodonTool(String tool) {
    return '$tool quita el codón entero';
  }

  @override
  String trimCaseCodonsTool(String tool) {
    return 'Mantener el marco de lectura ($tool)';
  }

  @override
  String trimCasesCaptionTool(String tool) {
    return 'Mismo alineamiento, dos formas. Columna a columna, trimAl elimina una a una las columnas que marca (×), así que las posiciones del codón (regla de colores) quedan fuera de orden (!). Con $tool sale el codón entero de cada columna marcada (+): las posiciones 1.ª, 2.ª y 3.ª siguen en orden para las particiones y la prueba de saturación.';
  }
}
