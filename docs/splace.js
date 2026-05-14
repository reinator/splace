// ========================================================================
// State
// ========================================================================
const state = {
    records: [],
    selectedGenes: new Set(),
    selectedFeatureTypes: new Set(),
    detectedDataType: "mt",
    pendingRemoveIndex: null,
    pendingEditIndex: null,
    headerTemplate: ["accession", "family", "genus", "species"],
    hiddenColumns: new Set(["source", "accession", "kingdom", "phylum"]),
    taxonomyImportCandidates: [],
    taxonomyImportResolver: null,
    duplicateGeneChoices: new Map(),
    progressUi: {
        mode: "generic",
        showMeta: false,
        showCurrent: true,
        currentLabel: "",
        footerText: "",
        idleText: "",
        successText: "",
        counterFormatter: null,
    },
    alignmentResultsData: null,
};

const MT_DEFAULT_GENES = ["COI", "COII", "COIII", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "ATP6", "ATP8"];
const CP_DEFAULT_GENES = ["rbcL", "matK", "ndhF", "atpB", "psaA", "psbA", "psbB", "psbC", "psbD", "psbE", "psbF", "psbH", "psbI", "psbJ", "psbK", "psbL", "psbM", "psbN", "psbT"];

const EXAMPLE_ACCESSIONS = "NC_024026.1, KJ184305.1, NC_024184.1, KJ556976.1, OP056886.1, PV818392.1, PX528615.1, NC_022707.1, KF356397.1, KJ642220.1, NC_057648.1, NC_082556.1, NC_023954.1";

const FEATURE_TYPE_INFO = {
    CDS: "Protein-coding sequences (genes translated into proteins)",
    rRNA: "Ribosomal RNA genes (12S, 16S, etc.)",
    tRNA: "Transfer RNA genes",
};

const FEATURE_TYPE_LABELS = {
    CDS: "PCGs",
    rRNA: "rRNAs",
    tRNA: "tRNAs",
};

function i18nText(key, fallback, params = {}) {
    const source = window.SPLACE_I18N;
    const raw = source && typeof source.t === "function" ? source.t(key) : null;
    const template = typeof raw === "string" && raw !== key ? raw : fallback;
    return Object.entries(params).reduce((text, [paramKey, value]) => {
        return text.replaceAll(`{${paramKey}}`, String(value));
    }, template);
}

// ========================================================================
// Gene Name Standardization Dictionaries
// ========================================================================
const MitochondrialGenes = {
    "12S": ["small subunit ribosomal RNA", "s-rRNA", "12S ribosormal RNA", "small ribosomal RNA subunit RNA", "12SrRNA", "12 ribosomal RNA", "rrnS", "12S ribosomal RNA subunit", "12S", "small ribosomal RNA", "small subunit ribosormal RNA", "12 rRNA", "12 S ribosomal RNA", "12S small subunit ribosomal RNA", "trnS", "Product small subunit ribosomal RNA", "12S-rRNA", "rRNA-12S", "12S ribosonal RNA", "12Srrn", "12S ribosome RNA", "12S ribsomal RNA", "12S rRNA", "12S ribosomal RNA", "12S ribosomal ribonucleic acid"],
    "16S": ["large subunit ribosomal RNA", "l-rRNA", "16S bibosomal RNA", "large ribosomal RNA subunit RNA", "rrnL", "16S ribosomal RNA subunit", "16S rivbosomal RNA", "l6S ribosomal RNA", "16S", "16S ribosamal RNA", "large ribosomal RNA", "16S rRNA", "16 S ribosomal RNA", "16S large subunit ribosomal RNA", "l-RNA", "16S-rRNA", "16Srrn", "16S ribosommal RNA", "16S ribosomal ribonucleic acid", "16S ribosomal RNA", "16S ribosome RNA", "16S recombinant RNA", "16S ribosomal RNA gene", "16S recombinant ribonucleic acid"],
    "ATP6": ["ATPase F0 subunit 6", "ATP synthase F0 subunit 6", "ATP synthase subunit 6", "ATPase 6", "ATP6", "ATP synthase FO subunit 6", "ATP synthase protein 6", "ATPase subunits 6", "ATP subunit 6", "ATPase subunit-6", "ATP synthase subunit F0 6", "adenosine triphosphatase subunit 6", "ATP synthase subunit-6", "ATP sythase F0 subunit 6", "ATPase 6 protein", "ATP synthase 6", "ATP synthetase F0 Subunit 6", "ATPase subuint 6", "ATPase sununit 6", "ATPase6 protein", "ATP Synthase Membrane Subunit 6", "TP synthase F0 subunit 6", "ATP sybthase F0 subunit 6", "ATP6ase", "ATP synthase A chain protein 6", "ATP synthetase subunit 6", "F0/F1 ATP synthase subunit 6", "disrupted ATP synthase F0 subunit 6", "ATPsynthase F0 subunit 6", "F1 ATPase subunit 6", "ATP sythase subunit 6", "adenine triphosphatase subunit 6", "F0-ATP synthase subunit 6", "F0-ATP synthase subunit6", "F0-ATPase subunit6", "ATP 6synthase 6", "adenosine triphosphate subunit 6", "ATPase subunit 6", "ATPase6", "adenosine triphosphate synthase-6"],
    "ATP8": ["ATPase F0 subunit 8", "ATP synthase F0 subunit 8", "ATP synthase subunit 8", "ATPase 8", "ATP8", "ATP synthase FO submit 8", "ATPase8", "ATP synthase protein 8", "ATPase subunits 8", "ATP subunit 8", "ATPase subunit-8", "ATP synthase subunit F0 8", "adenosine triphosphatase subunit 8", "ATP synthase subunit-8", "ATP sythase F0 subunit 8", "ATP synthase FO subunit 8", "ATPase 8 protein", "adenosine triphoshatase subunit 8", "ATP synthase 8", "ATPase subunit8", "ATP synthetase F0 Subunit 8", "adenosine triphosphate subunit 8", "ATPase8 protein", "ATP Synthase Membrane Subunit 8", "TP synthase F0 subunit 8", "product ATP synthase F0 subunit 8", "ATP sybthase F0 subunit 8", "ATP8ase", "ATP synthetase subunit 8", "F0/F1 ATP synthase subunit 8", "ATPsynthase F0 subunit 8", "F1 ATPase subunit 8", "ATP sythase subunit 8", "adenine triphosphatase subunit 8", "ATP syntahse F0 subunit 8", "F0-ATP synthase subunit 8", "F0-ATP synthase subunit8", "ATPase subunit 8", "F0-ATPase subunit8", "ATP synthase F0 subunit"],
    "COI": ["COX1", "cytochrome c oxidase subunit 1", "cytochrome oxidase gene", "coxidase subunit I", "COX", "c oxidase subunits I", "cytochrome-c-oxidase I", "Cytochrome Oxidase subunit I region", "c-oxidase subunit I", "c oxi- dase I", "c oxidase subunit-I", "cytochrome oxidase I region", "cytochrome oxydase I", "c oxydase I", "cytochrome-oxidase 1", "C Oxidase Gene Subunit-I", "C Oxidation I", "c oxidase I subunit", "cytochrome c oxidase I", "Cythocrome Oxidase I", "cytochrome oxidase I subunit", "Citochrome Oxidase I", "cytochromec oxidase I", "c oxidase submit I", "c oxidase unit I", "c oxidate subunit I", "cytochrome I", "cytochome oxidase subunit I", "cytochrome-c oxidase, subunit I", "cytochrome b oxidase subunit I", "cytochrome subunit I", "cytochrome-oxidase I", "COX-1", "cytochromoxidase I", "cytochrome oxidase 1", "cytochrome oxidase subunit 1", "cytochrome oxidase I", "C Oxidase type I", "cytochrome oxidase subunit I", "cytochrome oxidase I locus", "c oxidase subunit I sequences", "coxidase I", "c oxidase subunit I locus", "Cytochrome Oxidase unit I", "cytochrome oxidase subunits I", "cytochrome oxidase subunit I mtDNA", "cytochrome C oxidase subunit I", "Markers-Cytochrome Oxidase Subunit I", "C Oxidase, Subunit I", "chytochrome c oxidase subunit I", "I", "cytochrome oxidase subunit-1", "Cytochrome oxydase subunit 1", "cytochrome c oxidase subunit idase subunit I", "cytochrome c-oxidase subunit I", "cytochrome c oxidase subunits I", "cytchrome c oxidase subunit I", "subunit 1 of the cytochrome c oxidase", "cytochorome oxidase subunit I", "COI", "cytochrome c oxydase subunit 1", "cytochrome oxidase1", "COI protein", "cyt oxidase subunit 1", "cytochrome oxidase c subunit 1", "cytochrome oxidase c subunit I", "cytochrome oxydase subunit I", "cytochrome c oxidase polypeptide I", "cytochrome coxidase subunit I", "cytochrome c-oxidase subunit 1", "cytochrome c oxidase polypeptide 1", "CO I", "product: cytochrome c oxidase subunit I", "cytochome c oxidase subunit 1", "Cytochrome c oxidase subunit1", "cytochrome coxidase subunit 1", "cytochrome c oxidase subunit I"],
    "COII": ["COX-2", "c oxidase subunit II gene", "cytochrome oxidase II gene", "cytochrome oxidase subunit II gene", "cytochrome oxidase subunit II region", "coxidase subunit II", "c oxidase II gene", "cytochrome oxidase-II", "cytochrome c oxidase subunit 2", "chytochrome c oxidase subunit II", "II", "cytochrome c oxidase subunit II", "cytochrome oxidase subunit II", "cytochrome oxidase subunit 2", "cytohrome oxidase subunit II", "cytochrome oxidase subunit-2", "Cytochrome oxydase subunit 2", "cytochrome c oxidase subunit idase subunit II", "cytochrome c-oxidase subunit II", "cytochrome c oxidase subunits II", "cytochrome c oxidase II", "cytochrome oxidase II", "cytchrome c oxidase subunit II", "subunit 2 of the cytochrome c oxidase", "COX2", "cytochorom oxidase subunit II", "COII", "cytochrome c oxydase subunit 2", "CO2", "cytochrome oxidase subunit2", "COII protein", "cyt oxidase subunit 2", "cytochrome oxidase c subunit 2", "cytochrome oxidase c subunit II", "cytochrome oxydase subunit II", "cytochrome c oxidase polypeptide II", "cytochrome coxidase subunit II", "cytochome oxidase II", "cytochrome c-oxidase subunit 2", "cytochrome c oxidase polypeptide 2", "CO II", "cytochome c oxidase subunit 2", "Cytochrome c oxidase subunit2", "cytochrome coxidase subunit 2", "cytochome oxidase subunit 2", "C-terminal domain of Cytochrome c Oxidase subunit II"],
    "COIII": ["cytochrome oxidase subunit III gene", "c oxidase subunit III gene", "COX3", "COX 3", "COX-3", "c oxidase mitochondrial subunit III", "cytochrome c oxidase subunit 3", "cytochrome c oxidase subunit III", "chytochrome c oxidase subunit III", "cytochrome oxidase subunit III", "cytochrome oxidase subunit 3", "cytohrome oxidase subunit III", "cyctochrome c oxidase subunit III", "cutochrome oxidase subunit 3", "cytochrome oxidase subunit-3", "Cytochrome oxydase subunit 3", "cytochrome c oxidase subunit idase subunit III", "cytochrome c-oxidase subunit III", "cytochrome c oxidase subunits III", "cytochrome c oxidase III", "cytochrome oxidase III", "cytchrome c oxidase subunit III", "subunit 3 of the cytochrome c oxidase", "cytochorome oxidase subunit III", "COIII", "cytochrome c oxydase subunit 3", "CO3", "cytochrome oxidase subunit3", "cytochrome coxidase subunit III", "COIII protein", "cyt oxidase subunit 3", "cytochrome oxidase c subunit 3", "cytochrome oxidase c subunit III", "cytochrome oxydase subunit III", "cytochrome co oxidase subunit III", "cytochrome c oxidase subnunit III", "cytochrome oxidase sununit 3", "Cytochrome c oxidase polypeptide III", "cytochrome C oxidase asubunit 3", "cytochrome c oxidase sununit III", "cytochrome c-oxidase subunit 3", "cytochrome c oxidase polypeptide 3", "CO III", "cytochome c oxidase subunit 3", "cytochrome c oxidase subunit3", "cytochrome coxidase subunit 3"],
    "CYTB": ["Cytochrome b apoenzyme", "apoenzyme", "cytohrome b", "cytochome b", "cytochorome b", "cytchrome b", "cob", "cytochrome b protein", "Cythocrome b", "Cytb protein", "cytchorome b", "CYTB", "Cytochrome-b", "cytbochrome b", "apocytochrome b", "cytochromeb", "ctyb", "Cyt b", "apocytochome b", "cytochrome b", "cytochrome bc1"],
    "ND1": ["NADH dehydrogenase subunit 1", "NAD dehydrogenase subunit 1", "NADH dehydrogenase subunit-1", "NADH denydrogenase subunit 1", "NADH dehydrogenase 1", "NADH dehydrogenase subunits 1", "NADH dehydogenase subunit 1", "NADH dehydrogenase subunit1", "NADH dehydrogenase subunit #1", "Subunit 1 of the NADH ubiquinone oxidoreductase complex", "NADHsubunit 1", "NADH dehydrogenase subnit 1", "ND1", "NADH dehydrogenase subunit I", "NADH1", "NADH1 protein", "NADH subunit 1", "NADH-ubiquinone oxidoreductase chain 1", "NADH dehydrogenase, subunit 1", "NADH-ubiquinone oxidoreductase subunit I", "NADH 1", "NADH dehydrogenase subumit 1", "NADH-ubiquinone oxidoreductase subunit 1", "NADH ubiquinone oxidoreductase subunit 1", "nicotinamide adenine dinucleotide dehydrogenase subunit 1", "truncated NADH dehydrogenase subunit 1", "NADH-1", "NADH dehdrogenase subunit 1", "NaD1", "NADH dehydrogynase subunit 1"],
    "ND2": ["NADH dehydrogenase subunit 2", "NADH dehydrogenase subunit-2", "#NADH dehydrogenase subunit 2", "NADH denydrogenase subunit 2", "NADH dehydrogenase 2", "NADH dehydrogenase subunits 2", "NADH dehydrogenase subunit #2", "subunit 2 of the NADH ubiquinone oxidoreductase complex", "NADH deydrogenase subunit 2", "NADH dehydrogenase subnit 2", "ND2", "NADH dehydrogenase subunit II", "NADH2", "NADH2 protein", "NADH subunit 2", "NADH-ubiquinone oxidoreductase chain 2", "NADH dehdrogenase subnuit 2", "NADH-ubiquinone oxidoreductase subunit II", "NADH 2", "NADH dehydrogenase subumit 2", "NADH-ubiquinone oxidoreductase subunit 2", "NADH ubiquinone oxidoreductase subunit 2", "NADH dehydrognase subunit II", "nicotinamide adenine dinucleotide dehydrogenase subunit 2", "NADH dehydrogenase subunit2"],
    "ND3": ["NADH dehydrogenase subunit 3", "NAD dehydrogenase subunit 3", "NADH dehydrogenase subunit-3", "NADH denydrogenase subunit 3", "NADH dehydrogenase 3", "NADH dehydrogenase subunits 3", "subunit 3 of the NADH ubiquinone oxidoreductase complex", "NADH deydrogenase subunit 3", "NADH dehydrogenase subnit 3", "ND3", "NADH3", "NADH dehydrogenase subunit III", "NADH3 protein", "NADH subunit 3", "NADH-ubiquinone oxidoreductase chain 3", "ND3 NADH dehydrogenase subunit 3", "NADH dehydrogenasesubunit 3", "NADH dehydrogenase, subunit 3", "NADH-ubiquinone oxidoreductase subunit III", "truncated NADH dehydrogenase subunit 3", "NADH 3", "NADH dehydrogenase subumit 3", "NADH-ubiquinone oxidoreductase subunit 3", "NADH ubiquinone oxidoreductase subunit 3", "NADH dehydrogenase subunit3"],
    "ND4": ["NADH dehydrogenase subunit 4", "NAD dehydrogenase subunit 4", "NADH hehydrogenase subunit 4", "NADH dehrogenase subunit 4", "NADH dehydrogenase subunit-4", "NADH denydrogenase subunit 4", "NADH dehydrogenase 4", "NADH dehydrogenase subunits 4", "NADH dehydrosenase subunit 4", "subunit 4 of the NADH ubiquinone oxidoreductase complex", "NADH deydrogenase subunit 4", "NADH dehydrogenase subnit 4", "ND4", "NADH dehyodrogenase subunit 4", "NADH4", "NADH dehydrogenase subunit4", "NADH dehydrogenase subunit IV", "NADH4 protein", "NADH subunit 4", "NADH dehydrogenase sunbunit 4", "NADH-ubiquinone oxidoreductase chain 4", "NADH-ubiquinone oxidoreductase subunit IV", "NADH 4", "NADH dehydrogenase subumit 4", "NADH-ubiquinone oxidoreductase subunit 4", "NADH ubiquinone oxidoreductase subunit 4", "nicotinamide adenine dinucleotide dehydrogenase subunit 4", "truncated NADH dehydrogenase subunit 4"],
    "ND4L": ["NADH dehydrogenase subunit 4L", "NADH dehydrogenase subunit-4L", "NADH denydrogenase subunit 4L", "NADH dehydrogenase 4L", "ND4L", "NADH dehydrogenase subunits 4L", "subunit ND4L of the NADH ubiquinone oxidoreductase complex", "NADH4L protein", "NADH dehydrogenase subnit 4L", "NADH4L", "NADH dehydrogenase subunit 4 L", "NADH dehydrogenase subunit IV L", "NADH subunit 4L", "NADH-ubiquinone oxidoreductase chain 4L", "NADH-ubiquinone oxidoreductase subunit 4L", "NADH dehydrogenase, subunit 4L (complex I)", "NADH 4L", "NADH dehydrogenase subumit 4L", "HADH dehydrogenase 4L", "NADH dehydrogenase subujnit 4L", "NADH ubiquinone oxidoreductase subunit 4L", "nicotinamide adenine dinucleotide dehydrogenase subunit 4L", "NADH dehydrogenase subunit4L"],
    "ND5": ["NADH dehydrogenase subunit 5", "NAD dehydrogenase subunit 5", "NADH dehrogenase subunit 5", "NADH dehydrogenase subunit-5", "NADH hehydrogenase subunit 5", "NADH denydrogenase subunit 5", "NADH dehydrogenase 5", "ND5", "NADH dehydrogenase subunits 5", "NADH hydrogenase subunit 5", "subunit 5 of the NADH ubiquinone oxidoreductase complex", "NADH deydrogenase subunit 5", "NADH dehydrogenase subnit 5", "NADH dehydrodenase subunit 5", "NADH5", "HADA dehydrogenase subunit 5", "NADH dehydrogenase subunit V", "NADH5 protein", "NADH subunit 5", "NADH dehydrogenase subunit 5-0", "NADH-ubiquinone oxidoreductase chain 5", "NADH-ubiquinone oxidoreductase subunit V", "NADH dehydrogenase, subunit 5", "NADH dehydrogenase, subunit 5 complex I", "NADH 5", "NADH dehydrogenase subumit 5", "NADH-ubiquinone oxidoreductase subunit 5", "NADH ubiquinone oxidoreductase subunit 5", "nicotinamide adenine dinucleotide dehydrogenase subunit 5", "truncated NADH dehydrogenase subunit 5", "NADH dehydrogenase subunit5", "NADH dehydroghenase subunit 5"],
    "ND6": ["NADH dehydrogenase subunit 6", "NAD dehydrogenase subunit 6", "NADH dehydrogenase subunit-6", "NADH denydrogenase subunit 6", "NADH dehydrogenase 6", "ND6", "NADH dehydrogenase subunits 6", "subunit 6 of the NADH ubiquinone oxidoreductase complex", "NADH deydrogenase subunit 6", "NADH dehydrogenase subnit 6", "NADH6", "NADH dehydrogenase subunit VI", "NADH6 protein", "NADH subunit 6", "truncated NADH dehydrogenase subunit 6", "NADH dehygrogenase subunit 6", "NADH dehydrogenease subunit 6", "NADH-ubiquinone oxidoreductase chain 6", "NADH dehydrogenase, subunit 6", "NADH dsehydrogenase subunit 6", "NADH-ubiquinone oxidoreductase subunit VI", "NADH 6", "NADH dehydrogenase subumit 6", "NADH-ubiquinone oxidoreductase subunit 6", "NADH ubiquinone oxidoreductase subunit 6", "nicotinamide adenine dinucleotide dehydrogenase subunit 6", "NADH dehydrogenase subunit6"],
    "CR": ["Control region", "non-coding region", "putative control region", "control region 1", "control region ii", "control region i", "control region 2", "noncoding region", "pseudo control region", "cr", "control region (d-loop)", "d-loop control region", "similar to control region", "non coding region", "conrol region", "region: control region", "a+t-rich region", "putative control region 2", "d-loop region (= control reagion)", "pseudo contorl region", "a+t rich region", "c-rich region", "largest non-coding region", "d-loop", "d loop", "control region cr", "the control region", "control region c-rich sequence", "control region coretas sequence", "putative d-loop/control region", "d-loop containing region", "a+t rich", "at-rich region"],
};

const ChloroplastGenes = {
    "accD": ["acetyl-CoA carboxylase beta subunit", "acetyl-CoA carboxylase carboxyl transferase subunit beta", "beta-carboxyltransferase", "AccD", "acetyl-CoA carboxylase carboxyltransferase beta subunit", "accD protein", "Acetyl-CoA carboxylase, biotin carboxylase", "carboxyl transferase subunit beta", "acetyl-CoAcarboxylase beta subunit", "accD gene product", "Accda; acetyl-CoA carboxylase beta subunit"],
    "atpA": ["ATP synthase CF1 alpha subunit", "ATP synthase alpha subunit", "CF1 alpha subunit", "ATP synthase CF1 alpha chain", "AtpA", "atpA protein", "ATP synthase subunit alpha", "F1-ATPase alpha subunit", "ATP synthase alpha chain", "CF1-alpha", "Atpa; ATP synthase CF1 alpha subunit"],
    "atpB": ["ATP synthase CF1 beta subunit", "ATP synthase beta subunit", "CF1 beta subunit", "ATP synthase CF1 beta chain", "AtpB", "atpB protein", "ATP synthase subunit beta", "F1-ATPase beta subunit", "ATP synthase beta chain", "CF1-beta", "Atpb; ATP synthase CF1 beta subunit"],
    "atpE": ["ATP synthase CF1 epsilon subunit", "ATP synthase epsilon subunit", "CF1 epsilon subunit", "AtpE", "atpE protein", "ATP synthase subunit epsilon", "Atpe; ATP synthase CF1 epsilon subunit"],
    "atpF": ["ATP synthase CF0 subunit I", "ATP synthase CF0 I subunit", "AtpF", "atpF protein", "ATP synthase subunit b", "CF0 subunit I", "Atpf; ATP synthase CF0 subunit I"],
    "atpH": ["ATP synthase CF0 subunit III", "ATP synthase CF0 III subunit", "AtpH", "atpH protein", "CF0 subunit III", "ATP synthase subunit c", "Atph; ATP synthase CF0 subunit III"],
    "atpI": ["ATP synthase CF0 subunit IV", "ATP synthase CF0 IV subunit", "AtpI", "atpI protein", "CF0 subunit IV", "ATP synthase subunit a", "Atpi; ATP synthase CF0 subunit IV"],
    "ccsA": ["cytochrome c biogenesis protein", "CcsA", "cytochrome c heme attachment protein", "ccsA protein", "cytochrome c assembly protein"],
    "cemA": ["envelope membrane protein", "CemA", "cemA protein", "inner envelope membrane protein", "chloroplast envelope membrane protein"],
    "chlB": ["light-independent protochlorophyllide reductase subunit B", "ChlB", "protochlorophyllide reductase subunit B", "chlB protein"],
    "chlL": ["light-independent protochlorophyllide reductase subunit L", "ChlL", "protochlorophyllide reductase subunit L", "chlL protein"],
    "chlN": ["light-independent protochlorophyllide reductase subunit N", "ChlN", "protochlorophyllide reductase subunit N", "chlN protein"],
    "clpP": ["ATP-dependent Clp protease proteolytic subunit", "ClpP", "clpP protein", "Clp protease", "ATP-dependent Clp protease", "proteolytic subunit of Clp protease", "clpP1 protein"],
    "clpP1": ["ATP-dependent Clp protease proteolytic subunit 1", "ClpP1"],
    "infA": ["translation initiation factor 1", "IF1", "InfA", "infA protein", "translational initiation factor 1", "Infa; translation initiation factor 1", "initiation factor 1"],
    "matK": ["maturase K", "MatK", "maturase", "matK protein", "intron maturase", "maturaseK", "mat K", "maturase k protein", "maturaseK protein", "maturase enzyme", "MatK; maturase K", "Maturase K protein", "maturase type II intron", "maturase-K", "maturase protein K", "maturase-like protein", "maturase K-like protein", "intron-encoded maturase", "Group II intron maturase", "maturase protein"],
    "ndhA": ["NADH dehydrogenase subunit A", "NADH-plastoquinone oxidoreductase subunit A", "NdhA", "NAD(P)H dehydrogenase subunit A", "NADH dehydogenase subunit A", "ndhA protein", "NADH dehydrogenase ND1", "Ndha; NADH dehydrogenase subunit A"],
    "ndhB": ["NADH dehydrogenase subunit B", "NADH-plastoquinone oxidoreductase subunit B", "NdhB", "NAD(P)H dehydrogenase subunit B", "ndhB protein", "NADH dehydrogenase ND2", "Ndhb; NADH dehydrogenase subunit B"],
    "ndhC": ["NADH dehydrogenase subunit C", "NADH-plastoquinone oxidoreductase subunit C", "NdhC", "NAD(P)H dehydrogenase subunit C", "ndhC protein", "NADH dehydrogenase ND3", "Ndhc; NADH dehydrogenase subunit C"],
    "ndhD": ["NADH dehydrogenase subunit D", "NADH-plastoquinone oxidoreductase subunit D", "NdhD", "NAD(P)H dehydrogenase subunit D", "ndhD protein", "NADH dehydrogenase ND4", "Ndhd; NADH dehydrogenase subunit D"],
    "ndhE": ["NADH dehydrogenase subunit E", "NADH-plastoquinone oxidoreductase subunit E", "NdhE", "NAD(P)H dehydrogenase subunit E", "ndhE protein", "NADH dehydrogenase ND4L", "Ndhe; NADH dehydrogenase subunit E"],
    "ndhF": ["NADH dehydrogenase subunit F", "NADH-plastoquinone oxidoreductase subunit F", "NdhF", "NAD(P)H dehydrogenase subunit F", "ndhF protein", "NADH dehydrogenase ND5", "Ndhf; NADH dehydrogenase subunit F", "NADH dehydrogenase subunit 5"],
    "ndhG": ["NADH dehydrogenase subunit G", "NADH-plastoquinone oxidoreductase subunit G", "NdhG", "NAD(P)H dehydrogenase subunit G", "ndhG protein", "NADH dehydrogenase ND6", "Ndhg; NADH dehydrogenase subunit G"],
    "ndhH": ["NADH dehydrogenase subunit H", "NADH-plastoquinone oxidoreductase subunit H", "NdhH", "NAD(P)H dehydrogenase subunit H", "ndhH protein", "Ndhh; NADH dehydrogenase subunit H"],
    "ndhI": ["NADH dehydrogenase subunit I", "NADH-plastoquinone oxidoreductase subunit I", "NdhI", "NAD(P)H dehydrogenase subunit I", "ndhI protein", "Ndhi; NADH dehydrogenase subunit I"],
    "ndhJ": ["NADH dehydrogenase subunit J", "NADH-plastoquinone oxidoreductase subunit J", "NdhJ", "NAD(P)H dehydrogenase subunit J", "ndhJ protein", "Ndhj; NADH dehydrogenase subunit J"],
    "ndhK": ["NADH dehydrogenase subunit K", "NADH-plastoquinone oxidoreductase subunit K", "NdhK", "NAD(P)H dehydrogenase subunit K", "ndhK protein", "Ndhk; NADH dehydrogenase subunit K"],
    "petA": ["cytochrome f", "apocytochrome f", "PetA", "cytochrome f apoprotein", "petA protein", "Peta; cytochrome f"],
    "petB": ["cytochrome b6", "apocytochrome b6", "PetB", "cytochrome b6 apoprotein", "petB protein", "Petb; cytochrome b6"],
    "petD": ["cytochrome b6/f complex subunit IV", "PetD", "subunit IV of cytochrome b6/f complex", "petD protein", "Petd; cytochrome b6/f complex subunit IV"],
    "petG": ["cytochrome b6/f complex subunit V", "PetG", "petG protein", "Petg; cytochrome b6/f complex subunit V"],
    "petL": ["cytochrome b6/f complex subunit VI", "PetL", "petL protein", "Petl; cytochrome b6/f complex subunit VI"],
    "petN": ["cytochrome b6/f complex subunit VIII", "PetN", "petN protein", "Petn; cytochrome b6/f complex subunit VIII"],
    "psaA": ["photosystem I P700 chlorophyll a apoprotein A1", "PsaA", "PSI-A", "photosystem I subunit A", "psaA protein", "P700 apoprotein A1", "photosystem I reaction center subunit II", "Psaa; photosystem I P700 chlorophyll a apoprotein A1"],
    "psaB": ["photosystem I P700 chlorophyll a apoprotein A2", "PsaB", "PSI-B", "photosystem I subunit B", "psaB protein", "P700 apoprotein A2", "photosystem I reaction center subunit I", "Psab; photosystem I P700 chlorophyll a apoprotein A2"],
    "psaC": ["photosystem I subunit VII", "PsaC", "PSI-C", "psaC protein", "photosystem I iron-sulfur center", "Psac; photosystem I subunit VII"],
    "psaI": ["photosystem I subunit VIII", "PsaI", "PSI-I", "psaI protein", "Psai; photosystem I subunit VIII"],
    "psaJ": ["photosystem I subunit IX", "PsaJ", "PSI-J", "psaJ protein", "Psaj; photosystem I subunit IX"],
    "psbA": ["photosystem II protein D1", "PsbA", "PSII-D1", "D1 protein", "photosystem II D1 protein", "psbA protein", "photosystem II reaction center protein D1", "Psba; photosystem II protein D1", "photosystem II 32 kDa protein"],
    "psbB": ["photosystem II CP47 chlorophyll apoprotein", "PsbB", "PSII-B", "CP47", "photosystem II 47 kDa protein", "psbB protein", "photosystem II CP47 protein", "Psbb; photosystem II CP47 chlorophyll apoprotein"],
    "psbC": ["photosystem II CP43 chlorophyll apoprotein", "PsbC", "PSII-C", "CP43", "photosystem II 43 kDa protein", "psbC protein", "photosystem II CP43 protein", "Psbc; photosystem II CP43 chlorophyll apoprotein"],
    "psbD": ["photosystem II protein D2", "PsbD", "PSII-D2", "D2 protein", "photosystem II D2 protein", "psbD protein", "photosystem II reaction center protein D2", "Psbd; photosystem II protein D2"],
    "psbE": ["photosystem II cytochrome b559 alpha subunit", "PsbE", "cytochrome b559 alpha subunit", "psbE protein", "Psbe; photosystem II cytochrome b559 alpha subunit"],
    "psbF": ["photosystem II cytochrome b559 beta subunit", "PsbF", "cytochrome b559 beta subunit", "psbF protein", "Psbf; photosystem II cytochrome b559 beta subunit"],
    "psbH": ["photosystem II phosphoprotein", "PsbH", "photosystem II 10 kDa phosphoprotein", "psbH protein", "Psbh; photosystem II phosphoprotein"],
    "psbI": ["photosystem II protein I", "PsbI", "photosystem II reaction center subunit I", "psbI protein", "Psbi; photosystem II protein I"],
    "psbJ": ["photosystem II protein J", "PsbJ", "psbJ protein", "Psbj; photosystem II protein J"],
    "psbK": ["photosystem II protein K", "PsbK", "psbK protein", "Psbk; photosystem II protein K"],
    "psbL": ["photosystem II protein L", "PsbL", "psbL protein", "Psbl; photosystem II protein L"],
    "psbM": ["photosystem II protein M", "PsbM", "psbM protein", "Psbm; photosystem II protein M"],
    "psbN": ["photosystem II protein N", "PsbN", "psbN protein", "Psbn; photosystem II protein N"],
    "psbT": ["photosystem II protein T", "PsbT", "psbT protein", "Psbt; photosystem II protein T"],
    "psbZ": ["photosystem II protein Z", "PsbZ", "psbZ protein", "Psbz; photosystem II protein Z"],
    "rbcL": ["ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit", "ribulose 1,5-bisphosphate carboxylase/oxygenase large subunit", "ribulose bisophosphate carboxylase", "ribulose bisphosphate carboxylase large subunit", "Rubisco", "ribulose-1,5-bisphosphate carboxylase/oxygenase", "large subunit of riblose-1,5-bisphosphate carboxylase/oxygenase", "ribulose bisphosphate carboxylase large chain", "ribulose-1,5-bisphosphate carboxylase", "ribulose-bisphosphate carboxylase large chain", "ribulose-bisphosphate carboxylase large subunit", "RuBisCO large chain", "ribulose bisphosphate carboxylase/oxygenase", "ribulose-bisphosphate carboxylase/oxygenase", "ribulose-bisphosphate carboxylase", "ribulose 1,5-bisphosphate carboxylase/oxygenase", "ribulose 1,5-bisphosphate carboxylase large subunit", "ribulose 1,5-bisphosphate carboxylase", "Rubisco large subunit", "large subunit of RuBisCO", "rbcL protein", "RuBisCO,Ribulose Bisphosphate Carboxylase/Oxygenase large subunit", "ribulose 1,5-bisphosphate carboxylase oxygenase large subunit", "ribulose-1,5-bisphosphate carboxylase oxygenase large subunit", "ribulose-1,5-bisphosphate carboxylase oxygenase", "ribulose bisphosphate carboxylase large subunit precursor", "Ribulose-1,5-Bisphosphate Carboxylase-Oxygenase"],
    "rpl2": ["ribosomal protein L2", "50S ribosomal protein L2", "rpl2 protein", "large subunit ribosomal protein L2"],
    "rpl14": ["ribosomal protein L14", "50S ribosomal protein L14", "rpl14 protein", "large subunit ribosomal protein L14"],
    "rpl16": ["ribosomal protein L16", "50S ribosomal protein L16", "rpl16 protein", "large subunit ribosomal protein L16"],
    "rpl20": ["ribosomal protein L20", "50S ribosomal protein L20", "rpl20 protein", "large subunit ribosomal protein L20"],
    "rpl22": ["ribosomal protein L22", "50S ribosomal protein L22", "rpl22 protein"],
    "rpl23": ["ribosomal protein L23", "50S ribosomal protein L23", "rpl23 protein", "large subunit ribosomal protein L23"],
    "rpl32": ["ribosomal protein L32", "50S ribosomal protein L32", "rpl32 protein"],
    "rpl33": ["ribosomal protein L33", "50S ribosomal protein L33", "rpl33 protein", "large subunit ribosomal protein L33"],
    "rpl36": ["ribosomal protein L36", "50S ribosomal protein L36", "rpl36 protein"],
    "rpoA": ["RNA polymerase alpha subunit", "RpoA", "DNA-directed RNA polymerase alpha subunit", "RNA polymerase subunit alpha", "rpoA protein"],
    "rpoB": ["RNA polymerase beta subunit", "RpoB", "DNA-directed RNA polymerase beta subunit", "RNA polymerase subunit beta", "rpoB protein"],
    "rpoC1": ["RNA polymerase beta' subunit", "RpoC1", "DNA-directed RNA polymerase subunit gamma", "RNA polymerase subunit C1", "rpoC1 protein"],
    "rpoC2": ["RNA polymerase beta'' subunit", "RpoC2", "DNA-directed RNA polymerase beta'' subunit", "RNA polymerase subunit beta''", "rpoC2 protein"],
    "rps2": ["ribosomal protein S2", "30S ribosomal protein S2", "rps2 protein", "small subunit ribosomal protein S2"],
    "rps3": ["ribosomal protein S3", "30S ribosomal protein S3", "rps3 protein"],
    "rps4": ["ribosomal protein S4", "30S ribosomal protein S4", "rps4 protein", "small subunit ribosomal protein S4"],
    "rps7": ["ribosomal protein S7", "30S ribosomal protein S7", "rps7 protein", "small subunit ribosomal protein S7"],
    "rps8": ["ribosomal protein S8", "30S ribosomal protein S8", "rps8 protein"],
    "rps11": ["ribosomal protein S11", "30S ribosomal protein S11", "rps11 protein", "small subunit ribosomal protein S11"],
    "rps12": ["ribosomal protein S12", "30S ribosomal protein S12", "rps12 protein", "small subunit ribosomal protein S12"],
    "rps14": ["ribosomal protein S14", "30S ribosomal protein S14", "rps14 protein"],
    "rps15": ["ribosomal protein S15", "30S ribosomal protein S15", "rps15 protein", "small subunit ribosomal protein S15"],
    "rps16": ["ribosomal protein S16", "30S ribosomal protein S16", "rps16"],
    "rps18": ["ribosomal protein S18", "30S ribosomal protein S18", "rps18 protein"],
    "rps19": ["ribosomal protein S19", "30S ribosomal protein S19", "rps19 protein"],
    "rrn16S": ["16S ribosomal RNA", "small subunit ribosomal RNA", "rrn16S", "ribosomal RNA"],
    "rrn23S": ["23S ribosomal RNA", "large subunit ribosomal RNA", "ribosomal 23S RNA"],
    "rrn4.5S": ["4.5S ribosomal RNA", "ribosomal 4.5S RNA"],
    "rrn5S": ["5S ribosomal RNA", "ribosomal 5S RNA"],
    "ycf1": ["Ycf1 protein", "ycf1"],
    "ycf2": ["Ycf2 protein", "ycf2", "Ycf2"],
};

// ========================================================================
// Gene Name Standardization
// ========================================================================
function standardizeGeneName(rawName, product, dataType) {
    const dict = dataType === "cp" ? ChloroplastGenes : MitochondrialGenes;
    const altDict = dataType === "cp" ? MitochondrialGenes : ChloroplastGenes;

    // First try matching product against synonym lists (case-insensitive)
    if (product) {
        const productLower = product.toLowerCase().trim();
        for (const [key, synonyms] of Object.entries(dict)) {
            for (const syn of synonyms) {
                if (syn.toLowerCase() === productLower) {
                    return key;
                }
            }
        }
        // Try alt dictionary
        for (const [key, synonyms] of Object.entries(altDict)) {
            for (const syn of synonyms) {
                if (syn.toLowerCase() === productLower) {
                    return key;
                }
            }
        }
    }

    // Then try rawName (gene qualifier) against dictionary keys directly
    if (rawName) {
        for (const key of Object.keys(dict)) {
            if (key.toUpperCase() === rawName.toUpperCase()) return key;
        }
        for (const key of Object.keys(altDict)) {
            if (key.toUpperCase() === rawName.toUpperCase()) return key;
        }
        // Try rawName against synonyms too
        const rawLower = rawName.toLowerCase().trim();
        for (const [key, synonyms] of Object.entries(dict)) {
            for (const syn of synonyms) {
                if (syn.toLowerCase() === rawLower) {
                    return key;
                }
            }
        }
        for (const [key, synonyms] of Object.entries(altDict)) {
            for (const syn of synonyms) {
                if (syn.toLowerCase() === rawLower) {
                    return key;
                }
            }
        }
    }

    // Fallback: return raw name
    return rawName || product || null;
}

// ========================================================================
// tRNA Name Standardization (Leu1/Leu2, Ser1/Ser2)
// ========================================================================
function standardizeTrnaName(feat) {
    const product = feat.qualifiers.product || "";
    const gene = feat.qualifiers.gene || "";
    const note = feat.qualifiers.note || "";
    const anticodon = feat.qualifiers.anticodon || "";
    const raw = product || gene;

    if (!raw) {
        // Try to infer from note or anticodon
        if (note) {
            const noteAa = note.match(/tRNA-(\w+)/i) || note.match(/transfer\s+RNA[- ]+(\w+)/i);
            if (noteAa) return `tRNA-${noteAa[1].charAt(0).toUpperCase()}${noteAa[1].slice(1).toLowerCase()}`;
        }
        return "tRNA-Unknown";
    }

    const rawLower = raw.toLowerCase();
    const geneLower = gene.toLowerCase();
    const combined = (anticodon + " " + product + " " + note + " " + gene).toLowerCase();

    // Detect Leucine tRNAs (product "tRNA-Leu" or gene "trnL*")
    if (rawLower.includes("leu") || geneLower.startsWith("trnl")) {
        // Check gene qualifier for direct numbering (trnL1, trnL2, trnL(UAA), trnL(UAG))
        if (geneLower === "trnl1" || geneLower.includes("trnl(uaa)") || geneLower.includes("trnl(taa)")) {
            return "tRNA-Leu1";
        }
        if (geneLower === "trnl2" || geneLower.includes("trnl(uag)") || geneLower.includes("trnl(tag)")) {
            return "tRNA-Leu2";
        }
        // Check anticodon and other qualifiers
        // Leu1 = trnL(UAA) = anticodon TAA
        if (combined.includes("taa") || combined.includes("uaa")) {
            return "tRNA-Leu1";
        }
        // Leu2 = trnL(UAG) = anticodon TAG, recognizes CUN codons
        if (combined.includes("tag") || combined.includes("uag") || combined.includes("cun")) {
            return "tRNA-Leu2";
        }
        return "tRNA-Leu1"; // Default to Leu1 when anticodon cannot be determined
    }

    // Detect Serine tRNAs (product "tRNA-Ser" or gene "trnS*")
    if (rawLower.includes("ser") || geneLower.startsWith("trns")) {
        // Check gene qualifier for direct numbering (trnS1, trnS2, trnS(UGA), trnS(GCU))
        if (geneLower === "trns1" || geneLower.includes("trns(uga)") || geneLower.includes("trns(tga)")) {
            return "tRNA-Ser1";
        }
        if (geneLower === "trns2" || geneLower.includes("trns(gcu)") || geneLower.includes("trns(gct)")) {
            return "tRNA-Ser2";
        }
        // Check anticodon and other qualifiers
        // Ser1 = trnS(UGA) = anticodon TGA (recognizes UCN)
        if (combined.includes("tga") || combined.includes("uga") || combined.includes("ucn")) {
            return "tRNA-Ser1";
        }
        // Ser2 = trnS(GCU) = anticodon GCT (recognizes AGN/AGY)
        if (combined.includes("gct") || combined.includes("gcu") || combined.includes("agn") || combined.includes("agy")) {
            return "tRNA-Ser2";
        }
        return "tRNA-Ser1"; // Default to Ser1 when anticodon cannot be determined
    }

    // For other tRNAs, return a clean standardized name
    const aaMatch = raw.match(/tRNA-(\w+)/i);
    if (aaMatch) {
        return `tRNA-${aaMatch[1].charAt(0).toUpperCase()}${aaMatch[1].slice(1).toLowerCase()}`;
    }

    // Handle gene-qualifier format: trnX → tRNA-Xxx
    const trnMatch = raw.match(/^trn([A-Za-z])$/i);
    if (trnMatch) {
        const letter = trnMatch[1].toUpperCase();
        const aaNames = {
            A: "Ala", R: "Arg", N: "Asn", D: "Asp", C: "Cys", E: "Glu", Q: "Gln",
            G: "Gly", H: "His", I: "Ile", L: "Leu", K: "Lys", M: "Met", F: "Phe", P: "Pro",
            S: "Ser", T: "Thr", W: "Trp", Y: "Tyr", V: "Val"
        };
        if (aaNames[letter]) return `tRNA-${aaNames[letter]}`;
        return `tRNA-${letter}`;
    }

    // Handle "transfer RNA-Xxx" format
    const transferMatch = raw.match(/transfer\s+RNA[- ]+(\w+)/i);
    if (transferMatch) {
        return `tRNA-${transferMatch[1].charAt(0).toUpperCase()}${transferMatch[1].slice(1).toLowerCase()}`;
    }

    return raw;
}

// ========================================================================
// Auto Data Type Detection
// ========================================================================
function detectDataType() {
    let mtScore = 0;
    let cpScore = 0;

    for (const record of state.records) {
        for (const feat of record.features) {
            const rawGene = feat.qualifiers.gene || "";
            const rawProduct = feat.qualifiers.product || "";
            const combined = (rawGene + " " + rawProduct).toLowerCase();

            // Check for mt-specific markers
            for (const key of Object.keys(MitochondrialGenes)) {
                if (key.toLowerCase() === rawGene.toLowerCase() ||
                    key.toLowerCase() === rawProduct.toLowerCase()) {
                    mtScore++;
                    break;
                }
            }

            // Check for cp-specific markers
            for (const key of Object.keys(ChloroplastGenes)) {
                if (key.toLowerCase() === rawGene.toLowerCase() ||
                    key.toLowerCase() === rawProduct.toLowerCase()) {
                    cpScore++;
                    break;
                }
            }
        }
    }

    state.detectedDataType = cpScore > mtScore ? "cp" : "mt";
    return state.detectedDataType;
}

// ========================================================================
// GenBank Parser
// ========================================================================
function parseGenBank(text) {
    // Normalize line endings (Windows \r\n → \n) to avoid regex $ failures
    text = text.replace(/\r\n/g, "\n").replace(/\r/g, "\n");

    const result = { accession: "", organism: "", taxonomy: "", features: [], sequence: "", taxonomyRanks: { kingdom: "", phylum: "", class: "", order: "", family: "", authorship: "" } };

    const locusMatch = text.match(/^LOCUS\s+(\S+)/m);
    if (locusMatch) result.accession = locusMatch[1];

    const versionMatch = text.match(/^VERSION\s+(\S+)/m);
    if (versionMatch) result.accession = versionMatch[1];

    const orgMatch = text.match(/^\s+ORGANISM\s+(.+)/m);
    if (orgMatch) result.organism = orgMatch[1].trim();

    // Parse taxonomy lineage (lines after ORGANISM until next section)
    const orgIdx = text.indexOf("ORGANISM");
    if (orgIdx !== -1) {
        const afterOrg = text.substring(orgIdx);
        const orgLines = afterOrg.split("\n");
        const taxLines = [];
        for (let i = 1; i < orgLines.length; i++) {
            const line = orgLines[i];
            // Taxonomy lines are indented and contain semicolons or end with period
            if (/^\s{10,}/.test(line) && (line.includes(";") || line.trim().endsWith("."))) {
                taxLines.push(line.trim().replace(/\.$/, ""));
            } else {
                break;
            }
        }
        result.taxonomy = taxLines.join(" ");
    }

    const featStart = text.indexOf("FEATURES");
    const originStart = text.indexOf("ORIGIN");
    if (featStart !== -1 && originStart !== -1) {
        const featText = text.substring(featStart, originStart);
        result.features = parseFeatures(featText);
    }

    if (originStart !== -1) {
        const originText = text.substring(originStart);
        const seqLines = originText.split("\n").slice(1);
        const seqParts = [];
        for (const line of seqLines) {
            if (line.startsWith("//")) break;
            seqParts.push(line.replace(/[\s\d]/g, ""));
        }
        result.sequence = seqParts.join("").toLowerCase();
    }

    return result;
}

function extractTaxonomy(record) {
    const organism = record.editedOrganism || record.organism || "";
    const parts = organism.split(/\s+/);
    const genus = parts[0] || "";
    const species = parts.slice(1).join(" ") || "";

    const ranks = record.taxonomyRanks || {};

    const kingdom = ranks.kingdom || "";
    const phylum = ranks.phylum || "";
    const klass = ranks.class || "";
    const order = ranks.order || "";

    // Family: prefer taxonomyRanks, then editedFamily, then lineage parse
    let family = ranks.family || record.editedFamily || "";
    if (!family && record.taxonomy) {
        const taxParts = record.taxonomy.split(/;\s*/);
        for (const part of taxParts) {
            const trimmed = part.trim();
            if (trimmed.match(/idae$|aceae$|ales$/i)) {
                if (trimmed.match(/idae$|aceae$/i)) {
                    family = trimmed;
                } else if (!family) {
                    family = trimmed;
                }
            }
        }
    }

    return { kingdom, phylum, class: klass, order, family, genus, species, authorship: ranks.authorship || "" };
}

function countFeatureTypes(record) {
    let pcgs = 0, rrnas = 0, trnas = 0;
    for (const feat of record.features) {
        if (feat.type === "CDS") pcgs++;
        else if (feat.type === "rRNA") rrnas++;
        else if (feat.type === "tRNA") trnas++;
    }
    return { pcgs, rrnas, trnas };
}

function parseFeatures(featText) {
    const features = [];
    const lines = featText.split("\n");
    let current = null;
    let lastQualKey = null;

    for (let i = 1; i < lines.length; i++) {
        const line = lines[i];
        if (!line || line.trim() === "") continue;

        const featureMatch = line.match(/^     (\S+)\s+([\w<>().,:]+)/);
        if (featureMatch) {
            if (current) features.push(current);
            current = {
                type: featureMatch[1],
                locationStr: featureMatch[2],
                qualifiers: {},
            };
            lastQualKey = null;

            let j = i + 1;
            while (j < lines.length && lines[j].match(/^\s{21}[^/]/) && !lines[j].match(/^\s{21}\//)) {
                current.locationStr += lines[j].trim();
                j++;
            }
            continue;
        }

        if (current) {
            const qualMatch = line.match(/^\s{21}\/(\w+)="?([^"]*)"?$/);
            const qualMatchStart = line.match(/^\s{21}\/(\w+)="([^"]*)$/);
            if (qualMatch) {
                current.qualifiers[qualMatch[1]] = qualMatch[2];
                lastQualKey = qualMatch[1];
            } else if (qualMatchStart) {
                current.qualifiers[qualMatchStart[1]] = qualMatchStart[2];
                lastQualKey = qualMatchStart[1];
            } else if (lastQualKey && line.match(/^\s{21}[^/]/)) {
                const val = line.trim().replace(/"$/, "");
                current.qualifiers[lastQualKey] = (current.qualifiers[lastQualKey] || "") + " " + val;
            }
        }
    }
    if (current) features.push(current);
    return features;
}

function parseLocation(locStr, seqLen) {
    const segments = [];
    const isComplement = locStr.startsWith("complement(");
    let inner = locStr;
    if (isComplement) inner = inner.replace(/^complement\(/, "").replace(/\)$/, "");

    const isJoin = inner.startsWith("join(");
    if (isJoin) inner = inner.replace(/^join\(/, "").replace(/\)$/, "");

    const parts = inner.split(",");
    for (const part of parts) {
        const rangeMatch = part.trim().match(/<?\s*(\d+)\s*\.\.\s*>?\s*(\d+)/);
        if (rangeMatch) {
            segments.push({
                start: parseInt(rangeMatch[1]) - 1,
                end: parseInt(rangeMatch[2]),
                complement: isComplement,
            });
        }
    }
    return segments;
}

function extractSequence(fullSeq, locationStr) {
    const segments = parseLocation(locationStr, fullSeq.length);
    let seq = "";
    for (const seg of segments) {
        seq += fullSeq.substring(seg.start, seg.end);
    }
    if (segments.length > 0 && segments[0].complement) {
        seq = reverseComplement(seq);
    }
    return seq.toUpperCase();
}

function reverseComplement(seq) {
    const comp = {
        a: "t", t: "a", c: "g", g: "c", n: "n",
        A: "T", T: "A", C: "G", G: "C", N: "N"
    };
    return seq.split("").reverse().map(c => comp[c] || c).join("");
}

// ========================================================================
// NCBI Fetcher
// ========================================================================
const NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi";
const NCBI_SEARCH_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi";
const NCBI_SUMMARY_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi";
const NCBI_API_KEY_STORAGE = "splace.ncbi.apiKey";
const FETCH_DELAY_DEFAULT = 334;
const FETCH_DELAY_WITH_KEY = 100;

function getSavedNcbiApiKey() {
    try {
        return (localStorage.getItem(NCBI_API_KEY_STORAGE) || "").trim();
    } catch {
        return "";
    }
}

function looksLikeValidNcbiApiKey(value) {
    return /^[A-Za-z0-9_-]{20,}$/.test((value || "").trim());
}

function hasUsableNcbiApiKey() {
    return looksLikeValidNcbiApiKey(getSavedNcbiApiKey());
}

function getNcbiFetchDelay() {
    return hasUsableNcbiApiKey() ? FETCH_DELAY_WITH_KEY : FETCH_DELAY_DEFAULT;
}

function formatProgressCounter(current, total) {
    return i18nText("step1.progress.counter", { current, total });
}

function appendNcbiApiKey(params) {
    const apiKey = getSavedNcbiApiKey();
    if (looksLikeValidNcbiApiKey(apiKey)) {
        params.set("api_key", apiKey);
    }
    return params;
}

function renderTaxonomySearchButton(stateKey = "step1.action.search", spinning = false) {
    const btn = document.getElementById("taxonomySearchBtn");
    if (!btn) return;
    btn.innerHTML = `${spinning ? '<i class="fa-solid fa-spinner fa-spin mr-1"></i>' : '<i class="fa-solid fa-magnifying-glass mr-1"></i>'}<span data-i18n="${stateKey}">${i18nText(stateKey)}</span>`;
}

function renderFetchButton(stateKey = "step1.action.fetch", spinning = false) {
    const btn = document.getElementById("fetchBtn");
    if (!btn) return;
    btn.innerHTML = `${spinning ? '<i class="fa-solid fa-spinner fa-spin mr-1"></i>' : '<i class="fa-solid fa-magnifying-glass mr-1"></i>'}<span data-i18n="${stateKey}">${i18nText(stateKey)}</span>`;
}

function syncNcbiApiUi(messageKey) {
    const input = document.getElementById("ncbiApiKeyInput");
    const status = document.getElementById("ncbiApiKeyStatus");
    if (!input || !status) return;
    input.value = getSavedNcbiApiKey();
    const key = messageKey || (hasUsableNcbiApiKey() ? "step1.api.saved" : "step1.api.none");
    status.classList.toggle("saved", key === "step1.api.saved");
    status.classList.toggle("invalid", key === "step1.api.invalid");
    status.textContent = i18nText(key);
}

function openNcbiApiKeyModal() {
    syncNcbiApiUi();
    document.getElementById("ncbiApiKeyModal")?.classList.remove("hidden");
}

function closeNcbiApiKeyModal() {
    document.getElementById("ncbiApiKeyModal")?.classList.add("hidden");
}

function openClearNcbiApiKeyModal() {
    document.getElementById("clearNcbiApiKeyModal")?.classList.remove("hidden");
}

function closeClearNcbiApiKeyModal() {
    document.getElementById("clearNcbiApiKeyModal")?.classList.add("hidden");
}

function resetStep1Inputs() {
    const taxonomyInput = document.getElementById("taxonomySearchInput");
    const organelleSelect = document.getElementById("taxonomyOrganelleSelect");
    const completeCheckbox = document.getElementById("taxonomyCompleteCheckbox");
    const partialCheckbox = document.getElementById("taxonomyPartialCheckbox");
    const refseqOnly = document.getElementById("taxonomyRefseqOnly");
    const accessionInput = document.getElementById("accessionInput");
    const fetchStatus = document.getElementById("fetchStatus");
    const taxonomyFetchStatus = document.getElementById("taxonomyFetchStatus");
    const queryPreview = document.getElementById("taxonomyQueryPreview");

    if (taxonomyInput) taxonomyInput.value = "";
    if (organelleSelect) organelleSelect.value = "mitochondrion";
    if (completeCheckbox) completeCheckbox.checked = false;
    if (partialCheckbox) partialCheckbox.checked = false;
    if (refseqOnly) refseqOnly.checked = false;
    if (accessionInput) accessionInput.value = "";
    if (fetchStatus) {
        fetchStatus.textContent = "";
        fetchStatus.classList.add("hidden");
    }
    if (taxonomyFetchStatus) {
        taxonomyFetchStatus.textContent = "";
        taxonomyFetchStatus.classList.add("hidden");
    }
    if (queryPreview) queryPreview.textContent = "";
}

function clearNcbiApiKey() {
    const input = document.getElementById("ncbiApiKeyInput");
    if (input) input.value = "";
    try {
        localStorage.removeItem(NCBI_API_KEY_STORAGE);
    } catch {
        // Storage unavailable, ignored
    }
    syncNcbiApiUi("step1.api.none");
    closeClearNcbiApiKeyModal();
}

function getDuplicateChoiceKey(recordKey, geneName) {
    return `${recordKey}::${geneName}`;
}

window.setDuplicateGeneChoice = function (recordKey, geneName, selectedIndex) {
    state.duplicateGeneChoices.set(getDuplicateChoiceKey(recordKey, geneName), Number(selectedIndex));
    renderGeneSelection();
};

window.toggleDuplicateGeneChoice = function (recordKey, geneName, selectedIndex) {
    const choiceKey = getDuplicateChoiceKey(recordKey, geneName);
    state.duplicateGeneChoices.set(choiceKey, Number(selectedIndex));
    renderGeneSelection();
};

function saveNcbiApiKeyFromInput() {
    const input = document.getElementById("ncbiApiKeyInput");
    if (!input) return;
    const value = input.value.trim();
    if (value && !looksLikeValidNcbiApiKey(value)) {
        syncNcbiApiUi("step1.api.invalid");
        return;
    }
    try {
        if (value) localStorage.setItem(NCBI_API_KEY_STORAGE, value);
        else localStorage.removeItem(NCBI_API_KEY_STORAGE);
    } catch {
        // Storage unavailable, ignored
    }
    syncNcbiApiUi(value ? "step1.api.saved" : "step1.api.none");
}

window.refreshStep1UiTranslations = function () {
    renderTaxonomySearchButton();
    renderFetchButton();
    syncNcbiApiUi();
    updateProgressMeta();
};

window.refreshRecordsUiTranslations = function () {
    if (state.records.length > 0) {
        renderRecords();
    }
};

async function fetchAccession(accession) {
    const cached = await cacheGet(accession);
    if (cached) {
        return cached;
    }

    const url = `${NCBI_BASE}?db=nucleotide&id=${encodeURIComponent(accession)}&rettype=gb&retmode=text`;
    const params = appendNcbiApiKey(new URLSearchParams({
        db: "nucleotide",
        id: accession,
        rettype: "gb",
        retmode: "text",
    }));

    const resp = await fetch(`${NCBI_BASE}?${params.toString()}`);
    if (!resp.ok) throw new Error(`NCBI returned ${resp.status} for ${accession}`);

    const text = await resp.text();
    if (text.includes("Error") && text.length < 500) throw new Error(`NCBI error for ${accession}: ${text.trim()}`);

    const parsed = parseGenBank(text);
    parsed.source = "ncbi";

    await cachePut(accession, parsed, text);
    return parsed;
}

async function fetchMultiple(accessions, options = {}) {
    const results = [];
    const statusEl = document.getElementById("fetchStatus");
    const localStatusEl = options.statusElementId ? document.getElementById(options.statusElementId) : statusEl;
    localStatusEl.classList.remove("hidden");

    showProgress(
        options.progressTitle || i18nText("step1.progress.fetchGenericTitle"),
        formatProgressCounter(0, accessions.length),
        options.progressDetail || "",
        {
            mode: "ncbi",
            showMeta: true,
            currentLabel: i18nText("step1.progress.currentLabel"),
            footerText: i18nText("step1.progress.footer"),
            idleText: i18nText("step1.progress.statusIdle"),
            successText: i18nText("step1.progress.complete"),
            counterFormatter: formatProgressCounter,
        }
    );

    for (let i = 0; i < accessions.length; i++) {
        const acc = accessions[i].trim();
        if (!acc) continue;
        const currentLabel = options.progressItems?.[acc] || acc;
        localStatusEl.textContent = `Fetching ${i + 1}/${accessions.length}: ${acc}...`;
        updateProgress(i + 1, accessions.length, currentLabel);
        try {
            const record = await fetchAccession(acc);
            results.push(record);
            await new Promise(r => setTimeout(r, getNcbiFetchDelay()));
        } catch (e) {
            // Error fetching record, skipped
        }
    }
    localStatusEl.textContent = `Fetched ${results.length}/${accessions.length} records`;
    setTimeout(() => localStatusEl.classList.add("hidden"), 3000);

    hideProgress(`Loaded ${results.length} record${results.length !== 1 ? 's' : ''}`);

    return results;
}

async function searchNcbiGenomeIds(options) {
    const config = buildNcbiGenomeSearchConfig(options);
    if (!config) return [];

    const params = appendNcbiApiKey(new URLSearchParams({
        db: config.db,
        retmode: "json",
        retmax: "200",
        term: config.term,
    }));

    const resp = await fetch(`${NCBI_SEARCH_BASE}?${params.toString()}`);
    if (!resp.ok) throw new Error(`NCBI search returned ${resp.status}`);

    const data = await resp.json();
    return data?.esearchresult?.idlist || [];
}

function getSelectedGenomeScopes() {
    const scopes = [];
    if (document.getElementById("taxonomyCompleteCheckbox")?.checked) scopes.push("complete genome");
    if (document.getElementById("taxonomyPartialCheckbox")?.checked) scopes.push("partial genome");
    return scopes;
}

function buildTaxonomyReviewDetails(options) {
    const organelleKey = options?.organelle === "chloroplast"
        ? "step1.taxonomy.organelle.chloro"
        : "step1.taxonomy.organelle.mito";

    const scopeLabels = (Array.isArray(options?.completeness) ? options.completeness : []).map((scope) => {
        if (scope === "complete genome") return i18nText("step1.taxonomy.completeness.complete");
        if (scope === "partial genome") return i18nText("step1.taxonomy.completeness.partial");
        return scope;
    });

    return [
        {
            label: i18nText("step1.review.detailTaxon"),
            value: options?.taxon || "-",
        },
        {
            label: i18nText("step1.review.detailGenomeType"),
            value: i18nText(organelleKey),
        },
        {
            label: i18nText("step1.review.detailSearchTypes"),
            value: scopeLabels.join(" + "),
        },
        {
            label: i18nText("step1.review.detailRefseq"),
            value: i18nText(options?.refseqOnly ? "step1.review.refseq.on" : "step1.review.refseq.off"),
        },
    ];
}

function buildNcbiGenomeSearchConfig(options) {
    const taxon = (options.taxon || "").trim();
    if (!taxon) return null;

    const organelle = options.organelle === "chloroplast" ? "chloroplast" : "mitochondrion";
    const completeness = Array.isArray(options.completeness) ? options.completeness : [];
    if (completeness.length === 0) return null;

    const parts = [
        `\"${taxon}\"[Organism]`,
        completeness.length === 1
            ? `\"${completeness[0]}\"[Title]`
            : `(${completeness.map((entry) => `\"${entry}\"[Title]`).join(" OR ")})`,
        organelle === "chloroplast"
            ? `(chloroplast[Title] OR plastid[Title] OR chloroplast[Filter])`
            : `(mitochondrion[Title] OR mitochondrial[Title] OR mitochondrion[Filter])`,
        "biomol_genomic[PROP]",
    ];

    if (options.refseqOnly) {
        parts.push("srcdb_refseq[PROP]");
    }

    return {
        db: "nucleotide",
        term: parts.join(" AND "),
    };
}

async function fetchNcbiSummaries(ids) {
    if (!ids.length) return [];
    const params = appendNcbiApiKey(new URLSearchParams({
        db: "nucleotide",
        retmode: "json",
        id: ids.join(","),
    }));
    const resp = await fetch(`${NCBI_SUMMARY_BASE}?${params.toString()}`);
    if (!resp.ok) throw new Error(`NCBI summary returned ${resp.status}`);
    const data = await resp.json();
    const result = data?.result || {};
    return ids.map((id) => result[id]).filter(Boolean);
}

function classifyGenomeTitle(title) {
    const normalized = (title || "").toLowerCase();
    const isComplete = normalized.includes("complete");
    const isPartial = normalized.includes("partial") || normalized.includes("incomplete");
    return { isComplete, isPartial };
}

function summarizeGenomeTitles(records) {
    return records.reduce((summary, record) => {
        const title = record?.title || "";
        const { isComplete, isPartial } = classifyGenomeTitle(title);
        summary.total += 1;
        if (isComplete) summary.complete += 1;
        if (isPartial) summary.partial += 1;
        return summary;
    }, { total: 0, complete: 0, partial: 0 });
}

function getGenomeClassificationKey(title) {
    const { isComplete, isPartial } = classifyGenomeTitle(title);
    if (isComplete && isPartial) return "mixed";
    if (isComplete) return "complete";
    if (isPartial) return "partial";
    return "unknown";
}

function isUnverifiedGenome(title) {
    return /\bunverified\b/i.test(title || "");
}

function buildTaxonomyImportCandidates(ids, summaries) {
    return ids.map((id, index) => {
        const summary = summaries[index] || {};
        const title = summary.title || "";
        return {
            id,
            accession: summary.caption || summary.accessionversion || summary.uid || id,
            title,
            classification: getGenomeClassificationKey(title),
            unverified: isUnverifiedGenome(title),
            selected: true,
        };
    });
}

function summarizeUnverifiedCandidates(candidates) {
    return candidates.reduce((count, candidate) => count + (candidate.unverified ? 1 : 0), 0);
}

function renderTaxonomyImportCandidates() {
    const tbody = document.getElementById("taxonomyImportTableBody");
    const countEl = document.getElementById("taxonomyImportSelectionCount");
    const confirmBtn = document.getElementById("taxonomyImportConfirm");
    if (!tbody || !countEl || !confirmBtn) return;

    const total = state.taxonomyImportCandidates.length;
    const selectedCandidates = state.taxonomyImportCandidates.filter((candidate) => candidate.selected);
    const selected = selectedCandidates.length;

    // Calculate proportions for all genomes
    const completeCount = state.taxonomyImportCandidates.filter((c) => c.classification === "complete").length;
    const partialCount = state.taxonomyImportCandidates.filter((c) => c.classification === "partial").length;
    const unverifiedCount = state.taxonomyImportCandidates.filter((c) => c.unverified).length;

    const completePercent = total > 0 ? (completeCount / total) * 100 : 0;
    const partialPercent = total > 0 ? (partialCount / total) * 100 : 0;
    const unverifiedPercent = total > 0 ? (unverifiedCount / total) * 100 : 0;

    // Calculate proportions for selected genomes
    const selectedCompleteCount = selectedCandidates.filter((c) => c.classification === "complete").length;
    const selectedPartialCount = selectedCandidates.filter((c) => c.classification === "partial").length;
    const selectedUnverifiedCount = selectedCandidates.filter((c) => c.unverified).length;

    const selectedCompletePercent = selected > 0 ? (selectedCompleteCount / selected) * 100 : 0;
    const selectedPartialPercent = selected > 0 ? (selectedPartialCount / selected) * 100 : 0;
    const selectedUnverifiedPercent = selected > 0 ? (selectedUnverifiedCount / selected) * 100 : 0;

    const stackedBarHtml = `
        <div class="taxonomy-stacked-bar-container">
            <div class="taxonomy-stacked-bar-group">
                <div class="taxonomy-stacked-bar-label">${i18nText("step1.review.total")} (${total})</div>
                <div class="taxonomy-stacked-bar">
                    <div class="taxonomy-stacked-bar-segment complete" style="width: ${completePercent}%" title="${completeCount} complete">${completePercent > 8 ? Math.round(completePercent) + '%' : ''}</div>
                    <div class="taxonomy-stacked-bar-segment partial" style="width: ${partialPercent}%" title="${partialCount} partial">${partialPercent > 8 ? Math.round(partialPercent) + '%' : ''}</div>
                    <div class="taxonomy-stacked-bar-segment unverified" style="width: ${unverifiedPercent}%" title="${unverifiedCount} unverified">${unverifiedPercent > 0 ? Math.round(unverifiedPercent) + '%' : ''}</div>
                </div>
            </div>
            <div class="taxonomy-stacked-bar-group">
                <div class="taxonomy-stacked-bar-label">${i18nText("step1.review.selected")} (${selected})</div>
                <div class="taxonomy-stacked-bar">
                    <div class="taxonomy-stacked-bar-segment complete" style="width: ${selectedCompletePercent}%" title="${selectedCompleteCount} complete">${selectedCompletePercent > 8 ? Math.round(selectedCompletePercent) + '%' : ''}</div>
                    <div class="taxonomy-stacked-bar-segment partial" style="width: ${selectedPartialPercent}%" title="${selectedPartialCount} partial">${selectedPartialPercent > 8 ? Math.round(selectedPartialPercent) + '%' : ''}</div>
                    <div class="taxonomy-stacked-bar-segment unverified" style="width: ${selectedUnverifiedPercent}%" title="${selectedUnverifiedCount} unverified">${selectedUnverifiedPercent > 0 ? Math.round(selectedUnverifiedPercent) + '%' : ''}</div>
                </div>
            </div>
            <div class="taxonomy-stacked-bar-legend">
                <span class="taxonomy-stacked-bar-legend-item"><span class="taxonomy-stacked-bar-legend-dot complete"></span> ${i18nText("step1.review.class.complete")}</span>
                <span class="taxonomy-stacked-bar-legend-item"><span class="taxonomy-stacked-bar-legend-dot partial"></span> ${i18nText("step1.review.class.partial")}</span>
                ${unverifiedCount > 0 ? `<span class="taxonomy-stacked-bar-legend-item"><span class="taxonomy-stacked-bar-legend-dot unverified"></span> ${i18nText("step1.review.class.unverified")}</span>` : ""}
            </div>
        </div>
    `;

    countEl.innerHTML = stackedBarHtml;
    confirmBtn.disabled = selected === 0;

    tbody.innerHTML = state.taxonomyImportCandidates.map((candidate) => {
        const isActive = candidate.selected;
        const label = i18nText(`step1.review.class.${candidate.classification}`);
        const unverifiedLabel = i18nText("step1.review.class.unverified");
        const classificationHtml = candidate.unverified
            ? `<span class="taxonomy-review-badge unverified">${escHtml(unverifiedLabel)}</span>`
            : `<span class="taxonomy-review-badge ${candidate.classification}">${escHtml(label)}</span>`;
        return `
            <tr class="${candidate.unverified ? "unverified" : ""}">
                <td class="px-4 py-3 text-center">
                    <button
                        type="button"
                        class="taxonomy-review-switch ${isActive ? "active" : "inactive"}"
                        onclick="toggleTaxonomyImportCandidate('${escHtml(String(candidate.id))}', ${!isActive}); event.preventDefault();"
                        aria-pressed="${isActive ? "true" : "false"}"
                    >
                        <span class="switch-track ${isActive ? "active" : "inactive"}"><span class="switch-knob"></span></span>
                        <span class="taxonomy-review-switch-label ${isActive ? "active" : "inactive"}">${escHtml(i18nText(isActive ? "step1.review.switch.on" : "step1.review.switch.off"))}</span>
                    </button>
                </td>
                <td class="px-4 py-3 text-gray-700 text-center">${escHtml(String(candidate.accession))}</td>
                <td class="px-4 py-3 taxonomy-review-description">
                    <strong>${escHtml(String(candidate.accession))}</strong>
                    <p class="${candidate.unverified ? "unverified" : ""}">${escHtml(candidate.title || "-")}</p>
                </td>
                <td class="px-4 py-3 text-center">${classificationHtml}</td>
            </tr>
        `;
    }).join("");
}

function toggleTaxonomyImportCandidate(id, checked) {
    state.taxonomyImportCandidates = state.taxonomyImportCandidates.map((candidate) => (
        String(candidate.id) === String(id)
            ? { ...candidate, selected: checked }
            : candidate
    ));
    renderTaxonomyImportCandidates();
}

function setTaxonomyImportSelection(mode) {
    state.taxonomyImportCandidates = state.taxonomyImportCandidates.map((candidate) => {
        if (mode === "all") return { ...candidate, selected: true };
        if (mode === "none") return { ...candidate, selected: false };
        if (mode === "disable-unverified") return { ...candidate, selected: candidate.unverified ? false : candidate.selected };
        if (mode === "complete") return { ...candidate, selected: candidate.classification === "complete" && !candidate.unverified };
        if (mode === "partial") return { ...candidate, selected: candidate.classification === "partial" && !candidate.unverified };
        return candidate;
    });
    renderTaxonomyImportCandidates();
}

function closeTaxonomyImportModal(selectedIds = null) {
    document.getElementById("taxonomyImportModal").classList.add("hidden");
    const resolver = state.taxonomyImportResolver;
    state.taxonomyImportResolver = null;
    if (typeof resolver === "function") {
        resolver(selectedIds);
    }
}

function openTaxonomyImportModal({ taxon, summaryHtml, candidates, searchDetails = [] }) {
    document.getElementById("taxonomyImportTitle").textContent = i18nText("step1.review.title");
    const searchDetailsHtml = searchDetails.length > 0
        ? `<div class="taxonomy-review-meta">${searchDetails.map((detail) => `<span class="step1-filter-chip"><strong>${escHtml(detail.label)}:</strong> ${escHtml(detail.value)}</span>`).join("")}</div>`
        : "";
    document.getElementById("taxonomyImportSubtitle").innerHTML = `${i18nText("step1.review.subtitle")}${searchDetailsHtml}`;
    const unverifiedCount = summarizeUnverifiedCandidates(candidates);
    const unverifiedHtml = unverifiedCount > 0
        ? `<div class="mt-2"><span class="taxonomy-review-badge unverified">${i18nText("step1.review.class.unverified")}</span> ${i18nText("step1.review.unverifiedCount", { count: unverifiedCount })}<div class="mt-1 text-sm text-splace-red-700">${i18nText("step1.review.unverifiedAdvice")}</div></div>`
        : "";
    document.getElementById("taxonomyImportSummary").innerHTML = `${summaryHtml}${unverifiedHtml}`;
    state.taxonomyImportCandidates = candidates;
    renderTaxonomyImportCandidates();
    document.getElementById("taxonomyImportModal").classList.remove("hidden");
    return new Promise((resolve) => {
        state.taxonomyImportResolver = resolve;
    });
}

function updateProgressMeta() {
    const rateLabel = document.getElementById("progressRateLabel");
    const rateValue = document.getElementById("progressRateValue");
    const apiLabel = document.getElementById("progressApiLabel");
    const apiValue = document.getElementById("progressApiValue");
    const currentLabel = document.getElementById("progressCurrentLabel");
    const currentItem = document.getElementById("progressCurrentItem");
    const footerNote = document.getElementById("progressFooterNote");
    const metaGrid = document.getElementById("progressMetaGrid");
    const currentPanel = document.getElementById("progressCurrentPanel");
    if (!rateLabel || !rateValue || !apiLabel || !apiValue || !currentLabel || !currentItem || !metaGrid || !currentPanel) return;

    metaGrid.classList.toggle("hidden", !state.progressUi.showMeta);
    currentPanel.classList.toggle("hidden", state.progressUi.showCurrent === false);

    if (!state.progressUi.showMeta) {
        currentLabel.textContent = state.progressUi.currentLabel || "";
        if (!currentItem.textContent.trim() && state.progressUi.idleText) {
            currentItem.textContent = state.progressUi.idleText;
        }
        if (footerNote) footerNote.textContent = state.progressUi.footerText || "";
        return;
    }

    rateLabel.textContent = i18nText("step1.progress.rateLabel");
    rateValue.textContent = i18nText(hasUsableNcbiApiKey() ? "step1.progress.rateWithKey" : "step1.progress.rateDefault");
    apiLabel.textContent = i18nText("step1.progress.apiLabel");
    apiValue.textContent = i18nText(hasUsableNcbiApiKey() ? "step1.progress.apiEnabled" : "step1.progress.apiDisabled");
    currentLabel.textContent = i18nText("step1.progress.currentLabel");
    if (!currentItem.textContent.trim()) {
        currentItem.textContent = i18nText("step1.progress.statusIdle");
    }
    if (footerNote) footerNote.textContent = i18nText("step1.progress.footer");
}

function setProgressCurrentItem(text, options = {}) {
    const currentItem = document.getElementById("progressCurrentItem");
    if (!currentItem) return;
    currentItem.textContent = text || "";
    currentItem.classList.toggle("species-italic", !!options.italic);
}

// ========================================================================
// IndexedDB Cache
// ========================================================================
const DB_NAME = "splace_cache";
const DB_VERSION = 1;
const STORE_NAME = "genbank_records";
const CACHE_TTL = 7 * 24 * 60 * 60 * 1000;

function openDB() {
    return new Promise((resolve, reject) => {
        const req = indexedDB.open(DB_NAME, DB_VERSION);
        req.onupgradeneeded = () => {
            const db = req.result;
            if (!db.objectStoreNames.contains(STORE_NAME)) {
                db.createObjectStore(STORE_NAME, { keyPath: "accession" });
            }
        };
        req.onsuccess = () => resolve(req.result);
        req.onerror = () => reject(req.error);
    });
}

async function cacheGet(accession) {
    try {
        const db = await openDB();
        return new Promise((resolve) => {
            const tx = db.transaction(STORE_NAME, "readonly");
            const store = tx.objectStore(STORE_NAME);
            const req = store.get(accession);
            req.onsuccess = () => {
                const record = req.result;
                if (record && (Date.now() - record.timestamp < CACHE_TTL)) {
                    resolve(record.data);
                } else {
                    resolve(null);
                }
            };
            req.onerror = () => resolve(null);
        });
    } catch { return null; }
}

async function cachePut(accession, data, rawText) {
    try {
        const db = await openDB();
        const tx = db.transaction(STORE_NAME, "readwrite");
        tx.objectStore(STORE_NAME).put({
            accession,
            data,
            rawText,
            timestamp: Date.now(),
        });
    } catch {
        // Cache write error, ignored
    }
}

// ========================================================================
// Step Badge System
// ========================================================================

function setBadge(id, status) {
    const el = document.getElementById(id);
    if (!el) return;
    el.classList.remove("badge-required", "badge-optional", "badge-done");
    el.classList.add("badge-" + status);
}

function updateStepBadges() {
    const hasRecords = state.records.length > 0;
    const hasGenes = state.selectedGenes.size > 0;

    // Step 1: Load Data
    setBadge("step1-badge", hasRecords ? "done" : "required");

    // Step 2: Loaded Records (informational — green as soon as visible)
    setBadge("step2-badge", "done");

    // Step 3: Feature Types (optional — amber, user can keep defaults)
    setBadge("step3-badge", "optional");

    // Step 4: Select Markers
    if (hasRecords) setBadge("step4-badge", hasGenes ? "done" : "required");

    // Step 5: Heatmap — no badge (merged with step 4)

    // Step 5 (was 6): Export / Align
    if (window.isElectron) {
        const headerConfirmed = !!document.getElementById("headerPreviewInline")?.textContent?.trim();
        setBadge("step6-desktop-badge", headerConfirmed ? "done" : "required");
        setBadge("step6a-badge", headerConfirmed ? "done" : "required");
        setBadge("step6b-badge", "optional");
        setBadge("step6c-badge", headerConfirmed ? "optional" : "required");
    } else {
        setBadge("step6-web-badge", "optional");
    }
}

// ========================================================================
// UI Rendering
// ========================================================================
function renderRecords() {
    const section = document.getElementById("recordsSection");
    const tbody = document.getElementById("recordsTableBody");
    const count = document.getElementById("recordCount");

    if (state.records.length === 0) {
        section.classList.add("hidden");
        document.getElementById("featureTypesSection").classList.add("hidden");
        document.getElementById("genesSection").classList.add("hidden");
        document.getElementById("heatmapSection").classList.add("hidden");
        document.getElementById("downloadSection").classList.add("hidden");
        return;
    }

    // Auto-detect data type
    detectDataType();

    section.classList.remove("hidden");
    count.textContent = `(${state.records.length})`;

    // Sort records alphabetically by organism name
    const sortedIndices = state.records
        .map((r, i) => ({ record: r, originalIndex: i }))
        .sort((a, b) => (a.record.organism || "").localeCompare(b.record.organism || ""));

    tbody.innerHTML = sortedIndices.map(({ record: r, originalIndex: i }) => {
        const tax = extractTaxonomy(r);
        const counts = countFeatureTypes(r);
        const sourceClass = r.source === 'ncbi' ? 'bg-splace-blue-50 text-splace-blue-600' : 'bg-gray-100 text-gray-600';
        const sourceLabel = r.source === 'ncbi' ? i18nText('step2.record.source.ncbi') : i18nText('step2.record.source.file');
        const abbrevSpecies = tax.genus ? `${tax.genus.charAt(0)}.${tax.species ? ' ' + tax.species : ''}` : '';
        const fullSpecies = `${tax.genus}${tax.species ? ' ' + tax.species : ''}`;
        const displayName = fullSpecies || (r.editedOrganism || r.organism || r.accession || "Unknown");

        return `
            <tr class="group">
                <td class="px-3 py-2 font-mono text-gray-600" data-col="accession">${r.accession}</td>
                <td class="px-3 py-2 text-gray-500" data-col="kingdom">${tax.kingdom || "—"}</td>
                <td class="px-3 py-2 text-gray-500" data-col="phylum">${tax.phylum || "—"}</td>
                <td class="px-3 py-2 text-gray-500" data-col="class">${tax.class || "—"}</td>
                <td class="px-3 py-2 text-gray-500" data-col="order">${tax.order || "—"}</td>
                <td class="px-3 py-2 text-gray-500" data-col="family">${tax.family || "—"}</td>
                <td class="px-3 py-2 text-gray-500" data-col="genus"><span class="genus-italic">${tax.genus}</span></td>
                <td class="px-3 py-2 text-gray-500" data-col="species"><span class="species-italic">${abbrevSpecies}</span></td>
                <td class="px-3 py-2 text-gray-500" data-col="authorship">${tax.authorship || "—"}</td>
                <td class="px-3 py-2 text-center text-gray-600" data-col="pcgs">${counts.pcgs}</td>
                <td class="px-3 py-2 text-center text-gray-600" data-col="rrnas">${counts.rrnas}</td>
                <td class="px-3 py-2 text-center text-gray-600" data-col="trnas">${counts.trnas}</td>
                <td class="px-3 py-2 text-center" data-col="source"><span class="text-xs px-1.5 py-0.5 rounded ${sourceClass}">${sourceLabel}</span></td>
                <td class="px-3 py-2 text-right whitespace-nowrap">
                    <button onclick="editRecord(${i})" class="record-action-edit transition-colors text-sm leading-none mr-2" data-tippy-content="${i18nText('step2.record.tooltip.edit', { species: `<em>${escHtml(displayName)}</em>`, accession: escHtml(r.accession) })}"><i class="fa-solid fa-pen-to-square"></i></button>
                    <button onclick="removeRecord(${i})" class="record-action-remove transition-colors text-lg leading-none" data-tippy-content="${i18nText('step2.record.tooltip.remove', { species: `<em>${escHtml(displayName)}</em>`, accession: escHtml(r.accession) })}"><i class="fa-solid fa-xmark"></i></button>
                </td>
            </tr>
        `;
    }).join("");

    tbody.querySelectorAll("[data-tippy-content]").forEach((el) => {
        const htmlContent = el.getAttribute("data-tippy-content") || "";
        const titleText = htmlContent.replace(/<[^>]+>/g, " ").replace(/\s+/g, " ").trim();
        if (titleText) {
            el.setAttribute("title", titleText);
            el.setAttribute("aria-label", titleText);
        }
        window.createTippy && window.createTippy(el, el.getAttribute("data-tippy-content"));
    });
    window.initTooltips && window.initTooltips();

    // Only update gene selection if steps 3/4 are already visible (user already proceeded)
    if (!document.getElementById("featureTypesSection").classList.contains("hidden")) {
        renderGeneSelection();
    }
    applyColumnVisibility();
    updateStepBadges();
    filterRecords();
}

function filterRecords() {
    const input = document.getElementById("recordSearch");
    const clearBtn = document.getElementById("recordSearchClear");
    const q = (input?.value || "").toLowerCase().trim();

    if (clearBtn) clearBtn.classList.toggle("hidden", !q);

    const tbody = document.getElementById("recordsTableBody");
    if (!tbody) return;

    const rows = tbody.querySelectorAll("tr");
    let visible = 0;
    rows.forEach(row => {
        if (row.id === "recordSearchNoResults") return;

        // Save original cell HTML once per render cycle (cleared when tbody is replaced)
        if (!row.dataset.originals) {
            const cells = [...row.querySelectorAll("td")];
            row.dataset.originals = JSON.stringify(cells.map(td => td.innerHTML));
        }

        // Always restore clean HTML before re-applying highlight
        const originals = JSON.parse(row.dataset.originals);
        row.querySelectorAll("td").forEach((td, i) => { td.innerHTML = originals[i]; });

        const match = !q || row.textContent.toLowerCase().includes(q);
        row.style.display = match ? "" : "none";

        if (q && match) {
            row.querySelectorAll("td").forEach(td => _highlightTextNodes(td, q));
        }

        if (match) visible++;
    });

    // No-results message
    let noRow = document.getElementById("recordSearchNoResults");
    if (!q || visible > 0) {
        if (noRow) noRow.remove();
    } else {
        if (!noRow) {
            noRow = document.createElement("tr");
            noRow.id = "recordSearchNoResults";
            noRow.innerHTML = `<td colspan="15" class="px-4 py-6 text-center text-gray-400 italic">No records match "<span class="text-gray-600 not-italic font-medium"></span>"</td>`;
            tbody.appendChild(noRow);
        }
        noRow.querySelector("span").textContent = input.value;
    }

    tbody.querySelectorAll("[data-tippy-content]").forEach((el) => {
        const htmlContent = el.getAttribute("data-tippy-content") || "";
        const titleText = htmlContent.replace(/<[^>]+>/g, " ").replace(/\s+/g, " ").trim();
        if (titleText) {
            el.setAttribute("title", titleText);
            el.setAttribute("aria-label", titleText);
        }
        window.createTippy && window.createTippy(el, htmlContent);
    });
}

// Walk text nodes inside el and wrap occurrences of query with a <mark>
function _highlightTextNodes(el, query) {
    const walker = document.createTreeWalker(el, NodeFilter.SHOW_TEXT, null, false);
    const hits = [];
    let node;
    while ((node = walker.nextNode())) {
        if (node.nodeValue.toLowerCase().includes(query)) hits.push(node);
    }
    hits.forEach(textNode => {
        const text = textNode.nodeValue;
        const lower = text.toLowerCase();
        const frag = document.createDocumentFragment();
        let last = 0, idx;
        while ((idx = lower.indexOf(query, last)) !== -1) {
            if (idx > last) frag.appendChild(document.createTextNode(text.slice(last, idx)));
            const mark = document.createElement("mark");
            mark.className = "bg-yellow-200 text-yellow-900 rounded-sm px-0.5";
            mark.textContent = text.slice(idx, idx + query.length);
            frag.appendChild(mark);
            last = idx + query.length;
        }
        if (last < text.length) frag.appendChild(document.createTextNode(text.slice(last)));
        textNode.parentNode.replaceChild(frag, textNode);
    });
}

function applyColumnVisibility() {
    const table = document.querySelector(".records-table");
    if (!table) return;
    const cells = table.querySelectorAll("[data-col]");
    cells.forEach(cell => {
        const col = cell.dataset.col;
        cell.style.display = state.hiddenColumns.has(col) ? "none" : "";
    });
    // Sync eye icons
    document.querySelectorAll(".col-eye-toggle[data-col-toggle]").forEach(btn => {
        const col = btn.dataset.colToggle;
        const icon = btn.querySelector("i");
        if (state.hiddenColumns.has(col)) {
            icon.className = "fa-solid fa-eye-slash";
            btn.classList.add("text-gray-300");
            btn.classList.remove("text-gray-500");
        } else {
            icon.className = "fa-solid fa-eye";
            btn.classList.remove("text-gray-300");
            btn.classList.add("text-gray-500");
        }
    });
}

function toggleColumn(colKey, visible) {
    if (visible) {
        state.hiddenColumns.delete(colKey);
    } else {
        state.hiddenColumns.add(colKey);
    }
    applyColumnVisibility();
}

function renderGeneSelection() {
    if (state.records.length === 0) return;

    // Step 3: Feature Types — always visible once user proceeds
    document.getElementById("featureTypesSection").classList.remove("hidden");

    // Step 4: Select Markers — only visible when at least one feature type is selected
    const hasFeatureType = state.selectedFeatureTypes.size > 0;
    document.getElementById("genesSection").classList.toggle("hidden", !hasFeatureType);

    // Download / alignment section appear only once genes are selected
    const hasGenes = state.selectedGenes.size > 0;
    if (window.isElectron) {
        document.getElementById("downloadSection").classList.add("hidden");
        document.getElementById("nextStepSection").classList.toggle("hidden", !hasGenes);
    } else {
        document.getElementById("downloadSection").classList.toggle("hidden", !hasGenes);
    }

    const dataType = state.detectedDataType;

    // Collect all feature types and genes
    const featureTypes = new Map();
    const genesByType = new Map();

    for (let ri = 0; ri < state.records.length; ri++) {
        const record = state.records[ri];
        const assignedTrnas = new Set();
        for (const feat of record.features) {
            featureTypes.set(feat.type, (featureTypes.get(feat.type) || 0) + 1);

            let geneName = null;

            if (feat.type === "tRNA") {
                geneName = standardizeTrnaName(feat);
                // Dedup: if same tRNA name already seen in this record, use alternate
                if (geneName && assignedTrnas.has(geneName)) {
                    if (geneName === "tRNA-Ser1") geneName = "tRNA-Ser2";
                    else if (geneName === "tRNA-Leu1") geneName = "tRNA-Leu2";
                }
                if (geneName) assignedTrnas.add(geneName);
            } else {
                const rawGene = feat.qualifiers.gene || null;
                const rawProduct = feat.qualifiers.product || null;
                if (!rawGene && !rawProduct) continue;
                geneName = standardizeGeneName(rawGene, rawProduct, dataType);
            }

            if (!geneName) continue;

            if (!genesByType.has(feat.type)) genesByType.set(feat.type, new Map());
            const geneMap = genesByType.get(feat.type);
            const rawProduct = feat.qualifiers.product || geneName;

            // Compute sequence length for this feature
            let seqLen = 0;
            try {
                const seq = extractSequence(record.sequence, feat.locationStr);
                seqLen = seq.length;
            } catch (e) { /* ignore */ }

            const existing = geneMap.get(geneName) || { product: rawProduct, count: 0, recordIndices: new Set(), seqLengths: [], recordLengths: new Map() };
            existing.count++;
            existing.recordIndices.add(ri);
            if (seqLen > 0) {
                existing.seqLengths.push(seqLen);
                if (!existing.recordLengths.has(ri)) existing.recordLengths.set(ri, []);
                existing.recordLengths.get(ri).push(seqLen);
            }
            geneMap.set(geneName, existing);
        }
    }

    // Render feature type switches
    const ftList = document.getElementById("featureTypesList");
    const relevantTypes = ["CDS", "rRNA", "tRNA"];
    const availableTypes = relevantTypes.filter(t => featureTypes.has(t));

    ftList.innerHTML = availableTypes.map(type => {
        const count = featureTypes.get(type) || 0;
        const active = state.selectedFeatureTypes.has(type);
        const info = FEATURE_TYPE_INFO[type] || `Feature type: ${type}`;
        const displayLabel = FEATURE_TYPE_LABELS[type] || type;
        return `
            <label class="feature-switch inline-flex items-center gap-2 cursor-pointer select-none" data-tippy-content="${info}">
                <span class="switch-track ${active ? 'active' : ''}" onclick="toggleFeatureType('${type}', !${active}); event.preventDefault();">
                    <span class="switch-knob"></span>
                </span>
                <span class="text-sm font-medium ${active ? 'text-gray-800' : 'text-gray-500'}">${displayLabel}</span>
                <span class="text-xs text-gray-400">(${count})</span>
            </label>
        `;
    }).join("");

    renderGenes(genesByType);

    // Init Tippy on feature switches and record buttons
    window.initTooltips && window.initTooltips();
}

function renderGenes(genesByType) {
    if (!genesByType) {
        genesByType = new Map();
        const dataType = state.detectedDataType;
        for (let ri = 0; ri < state.records.length; ri++) {
            const record = state.records[ri];
            const assignedTrnas = new Set();
            for (const feat of record.features) {
                let geneName = null;

                if (feat.type === "tRNA") {
                    geneName = standardizeTrnaName(feat);
                    if (geneName && assignedTrnas.has(geneName)) {
                        if (geneName === "tRNA-Ser1") geneName = "tRNA-Ser2";
                        else if (geneName === "tRNA-Leu1") geneName = "tRNA-Leu2";
                    }
                    if (geneName) assignedTrnas.add(geneName);
                } else {
                    const rawGene = feat.qualifiers.gene || null;
                    const rawProduct = feat.qualifiers.product || null;
                    if (!rawGene && !rawProduct) continue;
                    geneName = standardizeGeneName(rawGene, rawProduct, dataType);
                }

                if (!geneName) continue;
                if (!genesByType.has(feat.type)) genesByType.set(feat.type, new Map());
                const geneMap = genesByType.get(feat.type);
                const rawProduct = feat.qualifiers.product || geneName;

                let seqLen = 0;
                try {
                    const seq = extractSequence(record.sequence, feat.locationStr);
                    seqLen = seq.length;
                } catch (e) { /* ignore */ }

                const existing = geneMap.get(geneName) || { product: rawProduct, count: 0, recordIndices: new Set(), seqLengths: [], recordLengths: new Map() };
                existing.count++;
                existing.recordIndices.add(ri);
                if (seqLen > 0) {
                    existing.seqLengths.push(seqLen);
                    if (!existing.recordLengths.has(ri)) existing.recordLengths.set(ri, []);
                    existing.recordLengths.get(ri).push(seqLen);
                }
                geneMap.set(geneName, existing);
            }
        }
    }

    const genesList = document.getElementById("genesList");
    const visibleGenes = new Map();

    for (const [type, geneMap] of genesByType) {
        if (!state.selectedFeatureTypes.has(type)) continue;
        for (const [name, info] of geneMap) {
            const existing = visibleGenes.get(name) || { product: info.product, count: 0, recordIndices: new Set(), types: [], seqLengths: [] };
            existing.count += info.count;
            if (info.recordIndices) {
                for (const idx of info.recordIndices) existing.recordIndices.add(idx);
            }
            existing.types.push(type);
            existing.seqLengths.push(...(info.seqLengths || []));
            visibleGenes.set(name, existing);
        }
    }

    const sortedGenes = [...visibleGenes.entries()].sort((a, b) => a[0].localeCompare(b[0]));

    genesList.innerHTML = sortedGenes.map(([name, info]) => {
        const active = state.selectedGenes.has(name);
        const lengths = info.seqLengths;
        const uniqueRecords = info.recordIndices ? info.recordIndices.size : info.count;
        let avgLen = 0, minLen = 0, maxLen = 0;
        if (lengths.length > 0) {
            avgLen = Math.round(lengths.reduce((a, b) => a + b, 0) / lengths.length);
            minLen = Math.min(...lengths);
            maxLen = Math.max(...lengths);
        }

        const escapedName = name.replace(/'/g, "\\'");

        return `
            <div class="gene-chip ${active ? 'active' : ''}" data-gene="${name}"
                onclick="toggleGene('${escapedName}', !state.selectedGenes.has('${escapedName}'))">
                <span class="switch-track-sm ${active ? 'active' : ''}">
                    <span class="switch-knob-sm"></span>
                </span>
                <span class="gene-name">${name}</span>
                <span class="gene-count">${uniqueRecords}/${state.records.length}</span>
                ${lengths.length > 0 ? `<span class="gene-len">~${avgLen} bp</span>` : ''}
            </div>
        `;
    }).join("");

    if (sortedGenes.length === 0) {
        genesList.innerHTML = '<p class="text-sm text-gray-400 italic">No genes found for selected feature types</p>';
    }

    // Init Tippy tooltips for gene chips
    if (typeof tippy !== "undefined") {
        // Destroy old instances
        document.querySelectorAll('.gene-chip').forEach(el => {
            if (el._tippy) el._tippy.destroy();
        });

        sortedGenes.forEach(([name, info]) => {
            const lengths = info.seqLengths;
            const uniqueRecords = info.recordIndices ? info.recordIndices.size : info.count;
            let avgLen = 0, minLen = 0, maxLen = 0;
            if (lengths.length > 0) {
                avgLen = Math.round(lengths.reduce((a, b) => a + b, 0) / lengths.length);
                minLen = Math.min(...lengths);
                maxLen = Math.max(...lengths);
            }

            const duplicateNote = info.count > uniqueRecords
                ? `<br><em>${info.count} features in ${uniqueRecords} records (some records have duplicates)</em>`
                : "";

            const tooltipHtml = `<strong>${name}</strong><br>` +
                `Feature: ${info.types.join(", ")}<br>` +
                `Records: ${uniqueRecords}/${state.records.length}` +
                duplicateNote +
                (lengths.length > 0
                    ? `<br>Avg: ${avgLen.toLocaleString()} bp<br>Min: ${minLen.toLocaleString()} bp | Max: ${maxLen.toLocaleString()} bp`
                    : "<br>No sequence data");

            const el = document.querySelector(`.gene-chip[data-gene="${name}"]`);
            if (el) {
                window.createTippy(el, tooltipHtml);
            }
        });
    }

    // Render missing gene warnings for selected genes
    const warningsDiv = document.getElementById("geneMissingWarnings");
    if (warningsDiv) {
        const selectedVisible = sortedGenes.filter(([name]) => state.selectedGenes.has(name));
        const warningItems = [];
        const missingSummaryTemplate = i18nText("step4.warning.missing.summary", "missing in {count} genome{suffix}");
        const missingTitle = i18nText("step4.warning.missing.title", "Missing Genes");
        const speciesLabel = i18nText("step4.warning.species", "Species");
        const accessionLabel = i18nText("step4.warning.accession", "Accession");
        const fileLabel = i18nText("step4.warning.file", "File");

        for (const [geneName, info] of selectedVisible) {
            const presentIndices = info.recordIndices || new Set();
            if (presentIndices.size >= state.records.length) continue;

            const missingRecords = [];
            for (let i = 0; i < state.records.length; i++) {
                if (!presentIndices.has(i)) {
                    const r = state.records[i];
                    const tax = extractTaxonomy(r);
                    const speciesName = `${tax.genus} ${tax.species}`.trim() || r.organism || "Unknown";
                    const accession = r.accession || 'N/A';
                    const ncbiUrl = r.accession ? `https://www.ncbi.nlm.nih.gov/nuccore/${r.accession}` : null;
                    const accHtml = ncbiUrl
                        ? `<a href="${ncbiUrl}" target="_blank" rel="noopener" class="text-splace-blue-600 hover:underline">${accession}</a>`
                        : `<span class="text-gray-400">${accession}</span>`;
                    const fileName = r.fileName || '—';
                    missingRecords.push({ speciesName, accHtml, fileName });
                }
            }

            if (missingRecords.length > 0) {
                let tableRows = missingRecords.map((rec, idx) => {
                    const bg = idx % 2 === 0 ? 'bg-amber-50' : 'bg-white';
                    return `<tr class="${bg}"><td class="px-3 py-1.5 text-sm italic">${rec.speciesName}</td><td class="px-3 py-1.5 text-sm">${rec.accHtml}</td><td class="px-3 py-1.5 text-sm text-gray-500">${rec.fileName}</td></tr>`;
                }).join("");
                warningItems.push(
                    `<div class="mb-3">` +
                    `<div class="flex items-center gap-2 mb-1">` +
                    `<i class="fa-solid fa-triangle-exclamation text-amber-500 flex-shrink-0"></i>` +
                    `<span class="text-sm font-semibold text-amber-800">${geneName}</span>` +
                    `<span class="text-xs text-amber-600">— ${i18nText("step4.warning.missing.summary", missingSummaryTemplate, { count: missingRecords.length, suffix: missingRecords.length > 1 ? "s" : "" })}</span>` +
                    `</div>` +
                    `<table class="w-full border border-amber-200 rounded-md overflow-hidden text-left">` +
                    `<thead><tr class="bg-amber-100"><th class="px-3 py-1 text-xs font-semibold text-amber-800 uppercase">${speciesLabel}</th><th class="px-3 py-1 text-xs font-semibold text-amber-800 uppercase">${accessionLabel}</th><th class="px-3 py-1 text-xs font-semibold text-amber-800 uppercase">${fileLabel}</th></tr></thead>` +
                    `<tbody>${tableRows}</tbody></table>` +
                    `</div>`
                );
            }
        }

        if (warningItems.length > 0) {
            warningsDiv.innerHTML =
                `<div class="bg-amber-50 border border-amber-200 rounded-lg p-4">` +
                `<div class="text-xs font-semibold text-amber-700 uppercase tracking-wider mb-3"><i class="fa-solid fa-circle-exclamation mr-1"></i>${missingTitle}</div>` +
                warningItems.join("") +
                `</div>`;
            warningsDiv.classList.remove("hidden");
        } else {
            warningsDiv.innerHTML = "";
            warningsDiv.classList.add("hidden");
        }
    }

    // Render duplicate gene warnings for PCGs and rRNAs
    const dupWarningsDiv = document.getElementById("geneDuplicateWarnings");
    if (dupWarningsDiv && genesByType) {
        const dupItems = [];
        const dupFeatureTypes = ["CDS", "rRNA"];
        const typeLabels = { CDS: "PCG", rRNA: "rRNA" };
        const duplicateTitle = i18nText("step4.warning.duplicate.title", "Duplicate Genes");
        const duplicateSummaryTemplate = i18nText("step4.warning.duplicate.summary", "duplicated in {count} genome{suffix}");
        const speciesLabel = i18nText("step4.warning.species", "Species");
        const accessionLabel = i18nText("step4.warning.accession", "Accession");
        const copiesLabel = i18nText("step4.warning.copies", "Copies (sizes)");
        const usedLabel = i18nText("step4.warning.used", "Used");
        const choiceHelp = i18nText("step4.warning.choiceHelp", "One copy is used per genome. The longest starts selected.");
        const defaultLabel = i18nText("step4.warning.default", "default");

        for (const ftype of dupFeatureTypes) {
            const fMap = genesByType.get(ftype);
            if (!fMap) continue;

            for (const [geneName, info] of fMap) {
                if (!state.selectedGenes.has(geneName)) continue;
                if (!info.recordLengths) continue;

                const dupRecords = [];
                for (const [ri, lengths] of info.recordLengths) {
                    if (lengths.length < 2) continue;
                    const r = state.records[ri];
                    const tax = extractTaxonomy(r);
                    const speciesName = `${tax.genus} ${tax.species}`.trim() || r.organism || "Unknown";
                    const ncbiUrl = r.accession ? `https://www.ncbi.nlm.nih.gov/nuccore/${r.accession}` : null;
                    const accHtml = ncbiUrl
                        ? `<a href="${ncbiUrl}" target="_blank" rel="noopener" class="text-splace-blue-600 hover:underline">${r.accession}</a>`
                        : `<span class="text-gray-400">${r.accession || 'N/A'}</span>`;
                    const recordKey = r.accession || String(ri);
                    const maxLen = Math.max(...lengths);
                    const lensList = lengths.map(l => `${l.toLocaleString()} bp`).join(", ");
                    const otherLengths = [];
                    for (const [otherRi, otherGeneLengths] of info.recordLengths) {
                        if (otherRi === ri) continue;
                        otherLengths.push(...otherGeneLengths);
                    }
                    const otherMean = otherLengths.length > 0
                        ? Math.round(otherLengths.reduce((sum, value) => sum + value, 0) / otherLengths.length)
                        : null;
                    const defaultChoiceIndex = lengths.findIndex((value) => value === maxLen);
                    const choiceKey = getDuplicateChoiceKey(recordKey, geneName);
                    const selectedIndex = state.duplicateGeneChoices.has(choiceKey)
                        ? state.duplicateGeneChoices.get(choiceKey)
                        : defaultChoiceIndex;
                    if (!state.duplicateGeneChoices.has(choiceKey)) {
                        state.duplicateGeneChoices.set(choiceKey, defaultChoiceIndex);
                    }
                    dupRecords.push({ speciesName, accHtml, maxLen, lensList, otherMean, lengths, selectedIndex, recordKey, geneName });
                }

                if (dupRecords.length === 0) continue;

                const typeLabel = typeLabels[ftype] || ftype;
                const tableRows = dupRecords.map((rec, idx) => {
                    const bg = idx % 2 === 0 ? "bg-orange-50" : "bg-white";
                    const choiceOptions = rec.lengths.map((value, optionIndex) => {
                        const checked = optionIndex === rec.selectedIndex ? " checked" : "";
                        const hint = optionIndex === rec.lengths.findIndex((len) => len === rec.maxLen) ? ` <span class="duplicate-choice-hint">${defaultLabel}</span>` : "";
                        return `<label class="duplicate-choice-item">` +
                            `<input type="checkbox"${checked} onchange="toggleDuplicateGeneChoice('${escHtml(String(rec.recordKey))}', '${escHtml(String(rec.geneName))}', ${optionIndex})">` +
                            `<span>${value.toLocaleString()} bp</span>${hint}` +
                            `</label>`;
                    }).join("");
                    return `<tr class="${bg}">` +
                        `<td class="px-3 py-1.5 text-sm italic">${rec.speciesName}</td>` +
                        `<td class="px-3 py-1.5 text-sm">${rec.accHtml}</td>` +
                        `<td class="px-3 py-1.5 text-sm text-gray-600">${rec.lensList}</td>` +
                        `<td class="px-3 py-1.5 text-sm">` +
                        `<div class="duplicate-choice-list">${choiceOptions}</div>` +
                        ` <div class="text-xs text-gray-400 mt-1">${choiceHelp}</div>` +
                        `${rec.otherMean ? ` <div class="text-xs text-gray-500 mt-0.5">${i18nText("step4.warning.meanOther", "Mean in other species: {value} bp", { value: rec.otherMean.toLocaleString() })}</div>` : ""}</td>` +
                        `</tr>`;
                }).join("");

                dupItems.push(
                    `<div class="mb-3">` +
                    `<div class="flex items-center gap-2 mb-1.5">` +
                    `<i class="fa-solid fa-copy text-orange-500 flex-shrink-0"></i>` +
                    `<span class="font-semibold text-orange-800">${geneName}</span>` +
                    `<span class="text-xs font-semibold text-orange-500 uppercase tracking-wide">${typeLabel}</span>` +
                    `<span class="text-sm text-orange-600">— ${i18nText("step4.warning.duplicate.summary", duplicateSummaryTemplate, { count: dupRecords.length, suffix: dupRecords.length > 1 ? "s" : "" })}</span>` +
                    `</div>` +
                    `<table class="w-full border border-orange-200 rounded-md overflow-hidden text-left">` +
                    `<thead><tr class="bg-orange-100">` +
                    `<th class="px-3 py-1 text-xs font-semibold text-orange-800 uppercase tracking-wide">${speciesLabel}</th>` +
                    `<th class="px-3 py-1 text-xs font-semibold text-orange-800 uppercase tracking-wide">${accessionLabel}</th>` +
                    `<th class="px-3 py-1 text-xs font-semibold text-orange-800 uppercase tracking-wide">${copiesLabel}</th>` +
                    `<th class="px-3 py-1 text-xs font-semibold text-orange-800 uppercase tracking-wide">${usedLabel}</th>` +
                    `</tr></thead>` +
                    `<tbody>${tableRows}</tbody></table>` +
                    `</div>`
                );
            }
        }

        if (dupItems.length > 0) {
            dupWarningsDiv.innerHTML =
                `<div class="bg-orange-50 border border-orange-200 rounded-lg p-4">` +
                `<div class="font-semibold text-orange-700 uppercase tracking-wider mb-3 text-sm">` +
                `<i class="fa-solid fa-copy mr-1"></i>${duplicateTitle}</div>` +
                dupItems.join("") +
                `</div>`;
            dupWarningsDiv.classList.remove("hidden");
        } else {
            dupWarningsDiv.innerHTML = "";
            dupWarningsDiv.classList.add("hidden");
        }
    }

    // Render heatmap
    renderHeatmap();
}

// ========================================================================
// Gene Presence Heatmap
// ========================================================================
function renderHeatmap() {
    const section = document.getElementById("heatmapSection");
    const container = document.getElementById("heatmapContainer");
    if (!section || !container) return;

    const selectedGenes = [...state.selectedGenes].sort();
    if (selectedGenes.length === 0 || state.records.length === 0) {
        section.classList.add("hidden");
        updateStepBadges();
        return;
    }

    // Build presence map: gene -> Set of record indices that have it
    const presenceMap = new Map();
    for (const gene of selectedGenes) {
        presenceMap.set(gene, new Set());
    }

    const dataType = state.detectedDataType;
    for (let ri = 0; ri < state.records.length; ri++) {
        const record = state.records[ri];
        const assignedTrnas = new Set();
        for (const feat of record.features) {
            if (!state.selectedFeatureTypes.has(feat.type)) continue;

            let geneName = null;
            if (feat.type === "tRNA") {
                geneName = standardizeTrnaName(feat);
                if (geneName && assignedTrnas.has(geneName)) {
                    if (geneName === "tRNA-Ser1") geneName = "tRNA-Ser2";
                    else if (geneName === "tRNA-Leu1") geneName = "tRNA-Leu2";
                }
                if (geneName) assignedTrnas.add(geneName);
            } else {
                const rawGene = feat.qualifiers.gene || null;
                const rawProduct = feat.qualifiers.product || null;
                if (!rawGene && !rawProduct) continue;
                geneName = standardizeGeneName(rawGene, rawProduct, dataType);
            }

            if (geneName && presenceMap.has(geneName)) {
                presenceMap.get(geneName).add(ri);
            }
        }
    }

    // Gene color by feature type
    const geneTypeColor = {};
    for (const gene of selectedGenes) {
        geneTypeColor[gene] = "#323795"; // default blue
    }

    // Build table HTML
    let html = '<table style="border-spacing:0;border-collapse:separate;width:100%;min-width:max-content;">';

    // Header row with gene names (rotated) — sticky top + sticky col for corner
    html += '<thead><tr class="heatmap-sticky-header"><td class="heatmap-sticky-col" style="min-width:160px;"></td>';
    for (const gene of selectedGenes) {
        const color = geneTypeColor[gene];
        html += `<td class="px-0.5" style="vertical-align:bottom;text-align:center;">` +
            `<div class="heatmap-gene-label" style="color:${color};">${gene}</div></td>`;
    }
    html += '</tr></thead><tbody>';

    // Build sorted index array (alphabetical by species name)
    const sortedRows = state.records.map((r, ri) => {
        const tax = typeof extractTaxonomy === "function" ? extractTaxonomy(r) : {};
        const speciesName = (tax.genus && tax.species)
            ? `${tax.genus} ${tax.species}`
            : (r.organism || r.accession || "Record " + (ri + 1));
        return { ri, speciesName, r };
    });
    sortedRows.sort((a, b) => a.speciesName.localeCompare(b.speciesName));

    // Data rows
    sortedRows.forEach(({ ri, speciesName, r }, ridx) => {
        const ncbiUrl = r.accession ? `https://www.ncbi.nlm.nih.gov/nuccore/${r.accession}` : null;
        const accessionHtml = ncbiUrl
            ? `<a href="${ncbiUrl}" target="_blank" rel="noopener" class="text-splace-blue-600 hover:underline font-semibold">(${r.accession})</a>`
            : `<strong>(${r.accession || 'N/A'})</strong>`;
        html += `<tr data-ridx="${ridx}">`;
        html += `<td class="heatmap-sticky-col pr-2 text-gray-600 whitespace-nowrap" style="min-width:160px;">` +
            `<span class="italic">${speciesName}</span> ${accessionHtml}</td>`;

        selectedGenes.forEach((gene, ci) => {
            const present = presenceMap.get(gene).has(ri);
            html += `<td class="px-0.5 py-0.5">` +
                `<div class="heatmap-cell ${present ? 'present' : 'absent'}" data-ridx="${ridx}" data-col="${ci}">` +
                `</div></td>`;
        });
        html += '</tr>';
    });

    html += '</tbody></table>';
    container.innerHTML = html;
    section.classList.remove("hidden");

    // Drag-to-scroll (only on heatmap cells, not species column)
    let isDown = false, startX, startY, scrollLeft, scrollTop, didDrag = false;
    container.addEventListener('mousedown', e => {
        // Only start drag from heatmap cells or gene headers, not species names or links
        if (e.target.closest('.heatmap-sticky-col') || e.target.closest('a')) return;
        isDown = true;
        didDrag = false;
        container.style.cursor = 'grabbing';
        container.style.userSelect = 'none';
        startX = e.pageX - container.offsetLeft;
        startY = e.pageY - container.offsetTop;
        scrollLeft = container.scrollLeft;
        scrollTop = container.scrollTop;
    });
    container.addEventListener('mouseleave', () => { isDown = false; container.style.cursor = 'grab'; container.style.userSelect = ''; });
    container.addEventListener('mouseup', () => { isDown = false; container.style.cursor = 'grab'; container.style.userSelect = ''; });
    container.addEventListener('mousemove', e => {
        if (!isDown) return;
        e.preventDefault();
        didDrag = true;
        container.scrollLeft = scrollLeft - (e.pageX - container.offsetLeft - startX);
        container.scrollTop = scrollTop - (e.pageY - container.offsetTop - startY);
    });

    // Crosshair: click to highlight, mouseleave clicked cell to reset
    const table = container.querySelector('table');
    const allCells = () => table ? table.querySelectorAll('tbody .heatmap-cell') : [];
    const allBodyRows = () => table ? table.querySelectorAll('tbody tr') : [];

    let crosshairActive = false;
    let lockedRidx = -1;
    let lockedCol = -1;
    let lockedCell = null;

    function activateCrosshair(ridx, col) {
        crosshairActive = true;
        // Dim all cells to 25%
        allCells().forEach(c => { c.style.opacity = '0.25'; });
        const rows = allBodyRows();
        rows.forEach((tr, rowIdx) => {
            const cells = tr.querySelectorAll('.heatmap-cell');
            if (rowIdx === ridx) {
                // Hovered row: only up to the selected column (inclusive)
                cells.forEach((c, ci) => {
                    c.style.opacity = ci <= col ? '1' : '0.25';
                });
            } else if (rowIdx < ridx) {
                // Column cells ABOVE the hovered row at 100%
                if (cells[col]) cells[col].style.opacity = '1';
            }
        });
        // Header gene names
        const hdrCells = table.querySelectorAll('thead .heatmap-gene-label');
        hdrCells.forEach((lbl, i) => {
            lbl.style.opacity = i === col ? '1' : '0.25';
        });
        // Species text + link
        rows.forEach((tr, rowIdx) => {
            const speciesTd = tr.querySelector('.heatmap-sticky-col');
            if (speciesTd) {
                const dimmed = rowIdx !== ridx;
                speciesTd.style.color = dimmed ? 'rgba(107,114,128,0.25)' : '';
                const link = speciesTd.querySelector('a');
                if (link) link.style.color = dimmed ? 'rgba(107,114,128,0.25)' : '';
            }
        });
    }

    function resetCrosshair() {
        crosshairActive = false;
        lockedRidx = -1;
        lockedCol = -1;
        lockedCell = null;
        allCells().forEach(c => { c.style.opacity = ''; });
        const hdrCells = table.querySelectorAll('thead .heatmap-gene-label');
        hdrCells.forEach(lbl => { lbl.style.opacity = ''; });
        allBodyRows().forEach(tr => {
            const speciesTd = tr.querySelector('.heatmap-sticky-col');
            if (speciesTd) {
                speciesTd.style.color = '';
                const link = speciesTd.querySelector('a');
                if (link) link.style.color = '';
            }
        });
    }

    // Click to activate crosshair (ignore if it was a drag)
    container.addEventListener('click', e => {
        if (didDrag) { didDrag = false; return; }
        const cell = e.target.closest('.heatmap-cell');
        if (!cell) return;
        const ridx = parseInt(cell.dataset.ridx, 10);
        const col = parseInt(cell.dataset.col, 10);
        if (crosshairActive && lockedRidx === ridx && lockedCol === col) {
            // Click same cell again: toggle off
            resetCrosshair();
        } else {
            lockedRidx = ridx;
            lockedCol = col;
            lockedCell = cell;
            activateCrosshair(ridx, col);
        }
    });

    // When mouse leaves the clicked cell, reset
    container.addEventListener('mouseover', e => {
        const cell = e.target.closest('.heatmap-cell');
        // If crosshair is active and mouse left the locked cell, reset
        if (crosshairActive && lockedCell && cell !== lockedCell) {
            resetCrosshair();
        }
    });

    container.addEventListener('mouseleave', () => {
        if (crosshairActive) resetCrosshair();
    });

    updateStepBadges();
}

// ========================================================================
// UI Event Handlers
// ========================================================================
function toggleFeatureType(type, checked) {
    if (checked) state.selectedFeatureTypes.add(type);
    else state.selectedFeatureTypes.delete(type);
    renderGeneSelection();
}

function toggleGene(name, checked) {
    if (checked) state.selectedGenes.add(name);
    else state.selectedGenes.delete(name);
    renderGenes();
    renderGeneSelection(); // updates download/alignment section visibility + heatmap
    updateStepBadges();
}

window.removeRecord = function (index) {
    state.pendingRemoveIndex = index;
    const r = state.records[index];
    const tax = extractTaxonomy(r);
    document.getElementById("removeModalText").innerHTML =
        `Remove <em>${tax.genus} ${tax.species}</em> (${r.accession})? This action cannot be undone.`;
    document.getElementById("removeModal").classList.remove("hidden");
};

window.toggleFeatureType = toggleFeatureType;
window.toggleGene = toggleGene;

function selectAllGenes() {
    document.querySelectorAll(".gene-chip").forEach(el => {
        const name = el.dataset.gene;
        if (name) state.selectedGenes.add(name);
    });
    renderGenes();
    renderGeneSelection();
    updateStepBadges();
}

function selectNoneGenes() {
    state.selectedGenes.clear();
    renderGenes();
    renderGeneSelection();
    updateStepBadges();
}

function selectCompleteGenes() {
    state.selectedGenes.clear();
    document.querySelectorAll(".gene-chip").forEach(el => {
        const name = el.dataset.gene;
        if (!name) return;
        const countText = el.querySelector('.gene-count');
        if (countText) {
            const parts = countText.textContent.split('/');
            if (parts.length === 2 && parts[0].trim() === parts[1].trim()) {
                state.selectedGenes.add(name);
            }
        }
    });
    renderGenes();
    renderGeneSelection();
    updateStepBadges();
}

function selectDefaultGenes() {
    state.selectedGenes.clear();
    // Ensure CDS feature type is enabled (default genes are PCGs)
    if (!state.selectedFeatureTypes.has("CDS")) {
        state.selectedFeatureTypes.add("CDS");
    }
    const defaults = state.detectedDataType === "mt" ? MT_DEFAULT_GENES : CP_DEFAULT_GENES;
    defaults.forEach(g => state.selectedGenes.add(g));
    renderGeneSelection();
    updateStepBadges();
}

function proceedToAnalysis() {
    renderGeneSelection();
    document.getElementById("featureTypesSection").scrollIntoView({ behavior: "smooth", block: "start" });
}

function fillExampleAccessions() {
    document.getElementById("accessionInput").value = EXAMPLE_ACCESSIONS;
    validateAccessionInput();
}

function fillTaxonomySearchExample() {
    document.getElementById("taxonomySearchInput").value = "Bufonidae";
    document.getElementById("taxonomyOrganelleSelect").value = "mitochondrion";
    document.getElementById("taxonomyRefseqOnly").checked = true;
    document.getElementById("taxonomyCompleteCheckbox").checked = true;
    document.getElementById("taxonomyPartialCheckbox").checked = false;
    validateTaxonomySearchInput();
}

// NCBI accession numbers: 1-6 letters (optional underscore) then 5+ digits, optionally .version
const ACCESSION_RE = /\b[A-Za-z]{1,6}_?\d{4,}(?:\.\d+)?\b/;

function validateAccessionInput() {
    const raw = document.getElementById("accessionInput").value;
    const btn = document.getElementById("fetchBtn");
    btn.disabled = !ACCESSION_RE.test(raw);
}

function validateTaxonomySearchInput() {
    const raw = document.getElementById("taxonomySearchInput")?.value || "";
    const btn = document.getElementById("taxonomySearchBtn");
    const scopes = getSelectedGenomeScopes();
    if (btn) btn.disabled = raw.trim().length < 3 || scopes.length === 0;

    const previewEl = document.getElementById("taxonomyQueryPreview");
    if (!previewEl) return;

    const config = buildNcbiGenomeSearchConfig({
        taxon: raw,
        organelle: document.getElementById("taxonomyOrganelleSelect")?.value,
        completeness: scopes,
        refseqOnly: document.getElementById("taxonomyRefseqOnly")?.checked,
    });
    previewEl.textContent = config ? config.term : "";
}

function mergeRecords(records) {
    const seen = new Set(state.records.map((record) => (record.accession || "").trim()).filter(Boolean));
    for (const record of records) {
        const accession = (record.accession || "").trim();
        if (accession && seen.has(accession)) continue;
        state.records.push(record);
        if (accession) seen.add(accession);
    }
}

function i18nText(key, fallbackOrParams, maybeParams) {
    const translate = window.SPLACE_I18N?.t;
    const hasFallback = typeof fallbackOrParams === "string";
    const params = hasFallback ? maybeParams : fallbackOrParams;
    let value = typeof translate === "function" ? translate(key, params) : key;

    if (value === key && hasFallback) {
        value = fallbackOrParams;
    }

    if (params && typeof params === "object") {
        value = value.replace(/\{(\w+)\}/g, (_, name) => {
            return params[name] == null ? `{${name}}` : String(params[name]);
        });
    }

    return value;
}

// ========================================================================
// Progress Modal
// ========================================================================
// ========================================================================
// Progress Modal Aprimorado
// ========================================================================
function showProgress(title, subtitle, sourceDetail = "", options = {}) {
    state.progressUi = {
        mode: options.mode || "generic",
        showMeta: !!options.showMeta,
        showCurrent: options.showCurrent !== false,
        currentLabel: options.currentLabel || "",
        footerText: options.footerText || "",
        idleText: options.idleText || "",
        successText: options.successText || "",
        counterFormatter: options.counterFormatter || null,
    };
    const progressCard = document.querySelector(".progress-modal-card");
    if (progressCard) {
        progressCard.classList.toggle("progress-modal-gbif", state.progressUi.mode === "gbif");
        progressCard.classList.toggle("progress-modal-ncbi", state.progressUi.mode === "ncbi");
    }
    document.getElementById("progressTitle").innerHTML = title;
    document.getElementById("progressSubtitle").textContent = subtitle || "";
    document.getElementById("progressBar").style.width = "0%";
    document.getElementById("progressPercentage").textContent = "0%";
    document.getElementById("progressDetail").innerHTML = sourceDetail;
    document.getElementById("progressStatusText").textContent = subtitle || "";
    setProgressCurrentItem(state.progressUi.idleText || "");
    document.getElementById("progressIcon").className = "fa-solid fa-spinner fa-spin text-splace-blue-600 text-2xl";
    updateProgressMeta();
    syncProgressPanels();
    document.getElementById("progressModal").classList.remove("hidden");
}

function updateProgress(current, total, detail, options = {}) {
    const pct = total > 0 ? Math.round((current / total) * 100) : 100;
    const counterText = typeof state.progressUi.counterFormatter === "function"
        ? state.progressUi.counterFormatter(current, total)
        : `${current} / ${total}`;
    document.getElementById("progressBar").style.width = pct + "%";
    document.getElementById("progressPercentage").textContent = pct + "%";
    document.getElementById("progressSubtitle").textContent = counterText;
    document.getElementById("progressStatusText").textContent = counterText;

    if (detail) {
        setProgressCurrentItem(detail, { italic: !!options.italicCurrentItem });
    }
    syncProgressPanels();
}

function hideProgress(successMsg) {
    const progressCard = document.querySelector(".progress-modal-card");
    if (progressCard) {
        progressCard.classList.remove("progress-modal-gbif", "progress-modal-ncbi");
    }
    document.getElementById("progressBar").style.width = "100%";
    document.getElementById("progressPercentage").textContent = "100%";
    document.getElementById("progressIcon").className = "fa-solid fa-circle-check text-green-500 text-2xl";
    if (successMsg) document.getElementById("progressTitle").innerHTML = successMsg;
    document.getElementById("progressStatusText").textContent = state.progressUi.successText || i18nText("step1.progress.complete");
    if (state.progressUi.mode !== "ncbi" || !document.getElementById("progressCurrentItem")?.textContent.trim()) {
        setProgressCurrentItem(successMsg || state.progressUi.successText || i18nText("step1.progress.complete"));
    }
    updateProgressMeta();
    syncProgressPanels();

    setTimeout(() => {
        document.getElementById("progressModal").classList.add("hidden");
    }, 2000); // Aumentei um pouco o tempo para o usuário conseguir ler o final
}

function syncProgressPanels() {
    const detailEl = document.getElementById("progressDetail");
    const metaGrid = document.getElementById("progressMetaGrid");
    const currentPanel = document.getElementById("progressCurrentPanel");
    const currentItem = document.getElementById("progressCurrentItem");
    const currentLabel = document.getElementById("progressCurrentLabel");
    if (!detailEl || !metaGrid || !currentPanel || !currentItem || !currentLabel) return;

    const detailText = detailEl.textContent.replace(/\s+/g, " ").trim();
    detailEl.classList.toggle("hidden", !detailText);

    const hasCurrentText = currentItem.textContent.replace(/\s+/g, " ").trim();
    const hasCurrentLabel = currentLabel.textContent.replace(/\s+/g, " ").trim();
    const showCurrent = state.progressUi.showCurrent !== false && !!(hasCurrentText || hasCurrentLabel);
    currentPanel.classList.toggle("hidden", !showCurrent);

    const metaVisible = !metaGrid.classList.contains("hidden");
    metaGrid.classList.toggle("progress-modal-meta-grid--standalone", metaVisible && !showCurrent && !detailText);
}

// ========================================================================
// Toggle Visual de Colunas (Chips)
// ========================================================================
function applyColumnVisibility() {
    const table = document.querySelector(".records-table");
    if (!table) return;
    const cells = table.querySelectorAll("[data-col]");
    cells.forEach(cell => {
        const col = cell.dataset.col;
        cell.style.display = state.hiddenColumns.has(col) ? "none" : "";
    });

    // Sincronizar ícones e cores dos chips
    document.querySelectorAll(".col-eye-toggle[data-col-toggle]").forEach(btn => {
        const col = btn.dataset.colToggle;
        const icon = btn.querySelector("i");

        if (state.hiddenColumns.has(col)) {
            icon.className = "fa-solid fa-eye-slash text-gray-400";
            btn.classList.replace("bg-splace-blue-50", "bg-white");
            btn.classList.replace("border-splace-blue-200", "border-gray-200");
            btn.classList.replace("text-splace-blue-700", "text-gray-600");
        } else {
            icon.className = "fa-solid fa-eye text-splace-blue-600";
            btn.classList.replace("bg-white", "bg-splace-blue-50");
            btn.classList.replace("border-gray-200", "border-splace-blue-200");
            btn.classList.replace("text-gray-600", "text-splace-blue-700");
        }
    });
}

// ========================================================================
// GBIF Taxonomy Lookup (dataFishing - Rabelo et al. 2025)
// ========================================================================
async function fetchTaxonomyFromGBIF(speciesName) {
    const url = `https://api.gbif.org/v1/species?name=${encodeURIComponent(speciesName)}`;
    const resp = await fetch(url, {
        method: "GET",
        headers: { "Accept": "application/json" }
    });
    if (!resp.ok) throw new Error(`GBIF returned ${resp.status}`);
    const data = await resp.json();
    if (data && data.results && data.results.length > 0) {
        const accepted = data.results.find(r => r.taxonomicStatus === "ACCEPTED") || data.results[0];
        return {
            kingdom: accepted.kingdom || "",
            phylum: accepted.phylum || "",
            class: accepted.class || "",
            order: accepted.order || "",
            family: accepted.family || "",
            authorship: accepted.authorship || "",
        };
    }
    return null;
}

async function fetchAllTaxonomy() {
    if (state.records.length === 0) return;

    const btn = document.getElementById("fetchTaxonomyBtn");
    btn.disabled = true;

    // Group records by organism name to avoid duplicate GBIF queries
    const speciesMap = new Map(); // organism name → [record indices]
    for (let i = 0; i < state.records.length; i++) {
        const organism = (state.records[i].editedOrganism || state.records[i].organism || "").trim();
        if (!organism) continue;
        if (!speciesMap.has(organism)) speciesMap.set(organism, []);
        speciesMap.get(organism).push(i);
    }

    const uniqueSpecies = [...speciesMap.keys()];
    const totalUnique = uniqueSpecies.length;
    const totalRecords = state.records.length;

    showProgress(
        i18nText("step2.progress.title"),
        i18nText("step2.progress.subtitle", { current: 0, total: totalUnique }),
        i18nText("step2.progress.detail"),
        {
            mode: "gbif",
            showMeta: false,
            currentLabel: i18nText("step2.progress.currentLabel"),
            footerText: i18nText("step2.progress.footer"),
            idleText: i18nText("step2.progress.subtitle", { current: 0, total: totalUnique }),
            successText: i18nText("step2.action.fetchTaxonomy"),
            counterFormatter: (current, total) => i18nText("step2.progress.subtitle", { current, total }),
        }
    );

    let successCount = 0;
    let errorCount = 0;
    let completed = 0;
    const BATCH_SIZE = 5;

    for (let i = 0; i < uniqueSpecies.length; i += BATCH_SIZE) {
        const batch = uniqueSpecies.slice(i, i + BATCH_SIZE);
        const promises = batch.map(async (organism) => {
            try {
                const taxonomy = await fetchTaxonomyFromGBIF(organism);
                if (taxonomy) {
                    // Apply taxonomy data to ALL records with this organism name
                    for (const idx of speciesMap.get(organism)) {
                        const record = state.records[idx];
                        record.taxonomyRanks = {
                            kingdom: taxonomy.kingdom,
                            phylum: taxonomy.phylum,
                            class: taxonomy.class,
                            order: taxonomy.order,
                            family: taxonomy.family,
                            authorship: taxonomy.authorship,
                        };
                        if (taxonomy.family && !record.editedFamily) {
                            record.editedFamily = taxonomy.family;
                        }
                    }
                    successCount++;
                } else {
                    errorCount++;
                }
            } catch {
                errorCount++;
            }

            completed++;
            updateProgress(completed, totalUnique, organism, { italicCurrentItem: true });
        });

        await Promise.all(promises);
        if (i + BATCH_SIZE < uniqueSpecies.length) {
            await new Promise(r => setTimeout(r, 200));
        }
    }

    const errorSuffix = errorCount > 0 ? i18nText("step2.progress.doneErrors", { count: errorCount }) : "";
    const uniqueSuffix = totalUnique < totalRecords ? i18nText("step2.progress.doneUnique", { unique: totalUnique, records: totalRecords }) : "";
    hideProgress(i18nText("step2.progress.done", { success: successCount, errors: errorSuffix, unique: uniqueSuffix }));
    btn.disabled = false;
    renderRecords();
    renderGeneSelection(); // reveal steps 3 & 4 now that taxonomy is done
}

// ========================================================================
// Edit Record
// ========================================================================
function validateSpeciesName(name) {
    if (!name || !name.trim()) return "Species name is required";
    if (!/^[A-Za-z\s\-]+$/.test(name.trim())) return "Only letters, spaces, and hyphens are allowed";
    const parts = name.trim().split(/\s+/);
    if (parts.length < 2) return "Species name must include genus and epithet (e.g., Rhinella marina)";
    if (parts[0][0] !== parts[0][0].toUpperCase()) return "Genus must start with an uppercase letter";
    return null;
}

function validateFamilyName(name) {
    if (!name || !name.trim()) return null; // family is optional
    if (!name.trim().match(/idae$|aceae$/i)) return "Family name must end with -idae or -aceae";
    return null;
}

function formatFileSize(bytes) {
    if (bytes < 1024) return bytes + " B";
    if (bytes < 1024 * 1024) return (bytes / 1024).toFixed(1) + " KB";
    return (bytes / (1024 * 1024)).toFixed(1) + " MB";
}

window.editRecord = function (index) {
    state.pendingEditIndex = index;
    const r = state.records[index];
    const tax = extractTaxonomy(r);
    const counts = countFeatureTypes(r);

    // Build info section
    let info = '<div class="space-y-1">';
    info += `<div><span class="font-medium">Accession:</span> ${r.accession}</div>`;
    info += `<div><span class="font-medium">Source:</span> ${r.source === "ncbi" ? "NCBI" : "File"}</div>`;
    if (r.fileName) info += `<div><span class="font-medium">File:</span> ${r.fileName}</div>`;
    if (r.fileSize) info += `<div><span class="font-medium">Size:</span> ${formatFileSize(r.fileSize)}</div>`;
    if (r.fileDate) info += `<div><span class="font-medium">Date:</span> ${r.fileDate.toLocaleDateString()}</div>`;
    info += `<div><span class="font-medium">Sequence:</span> ${r.sequence ? r.sequence.length.toLocaleString() + " bp" : "N/A"}</div>`;
    info += `<div><span class="font-medium">Features:</span> ${counts.pcgs} PCGs, ${counts.rrnas} rRNAs, ${counts.trnas} tRNAs</div>`;
    info += "</div>";
    document.getElementById("editModalInfo").innerHTML = info;

    // Populate editable fields
    const currentOrganism = r.editedOrganism || r.organism || "";
    document.getElementById("editSpeciesInput").value = currentOrganism;
    document.getElementById("editGenusDisplay").value = currentOrganism.split(/\s+/)[0] || "";
    document.getElementById("editKingdomInput").value = tax.kingdom || "";
    document.getElementById("editPhylumInput").value = tax.phylum || "";
    document.getElementById("editClassInput").value = tax.class || "";
    document.getElementById("editOrderInput").value = tax.order || "";
    document.getElementById("editFamilyInput").value = tax.family || "";

    // Clear errors/status
    document.getElementById("editSpeciesError").classList.add("hidden");
    document.getElementById("editFamilyError").classList.add("hidden");
    document.getElementById("editGbifStatus").classList.add("hidden");

    // Set GBIF button state based on current species
    const speciesErr = validateSpeciesName(currentOrganism);
    document.getElementById("editGbifBtn").disabled = !!speciesErr;
    document.getElementById("editGbifBtn").dataset.fetchedAuthorship = "";

    document.getElementById("editModal").classList.remove("hidden");
};

// ========================================================================
// Header Builder
// ========================================================================
const HEADER_FIELDS = [
    { key: "species", label: "Species", getter: (r, tax) => (r.editedOrganism || r.organism || "Unknown").replace(/\s+/g, "_") },
    { key: "accession", label: "Accession", getter: (r, tax) => r.accession || "NA" },
    { key: "kingdom", label: "Kingdom", getter: (r, tax) => tax.kingdom || "NA" },
    { key: "phylum", label: "Phylum", getter: (r, tax) => tax.phylum || "NA" },
    { key: "class", label: "Class", getter: (r, tax) => tax.class || "NA" },
    { key: "order", label: "Order", getter: (r, tax) => tax.order || "NA" },
    { key: "family", label: "Family", getter: (r, tax) => tax.family || "NA" },
    { key: "genus", label: "Genus", getter: (r, tax) => tax.genus || "NA" },
    { key: "authorship", label: "Authorship", getter: (r, tax) => tax.authorship || "NA" },
];

function getHeaderSeparator() {
    const checked = document.querySelector('input[name="headerSeparator"]:checked');
    return checked ? checked.value : "_";
}

function buildHeader(record, tax) {
    const sep = getHeaderSeparator();
    return state.headerTemplate
        .map(key => {
            const field = HEADER_FIELDS.find(f => f.key === key);
            return field ? field.getter(record, tax) : key;
        })
        .join(sep);
}

function checkHeaderDuplicates() {
    const headers = new Map();
    for (const record of state.records) {
        const tax = extractTaxonomy(record);
        const header = buildHeader(record, tax);
        const existing = headers.get(header) || [];
        existing.push(record.accession);
        headers.set(header, existing);
    }
    const duplicates = [];
    for (const [header, accessions] of headers) {
        if (accessions.length > 1) {
            duplicates.push({ header, accessions });
        }
    }
    return duplicates;
}

function renderHeaderBuilder() {
    const availableContainer = document.getElementById("headerAvailableFields");
    const selectedContainer = document.getElementById("headerSelectedFields");
    const previewEl = document.getElementById("headerPreview");
    const warningEl = document.getElementById("headerDuplicateWarning");
    const warningTextEl = document.getElementById("headerDuplicateText");
    const confirmBtn = document.getElementById("headerModalConfirm");

    // Available fields = all fields NOT in template
    const available = HEADER_FIELDS.filter(f => !state.headerTemplate.includes(f.key));
    availableContainer.innerHTML = available.map(f =>
        `<span class="header-chip available" onclick="addHeaderField('${f.key}')">${f.label}</span>`
    ).join("");

    // Selected fields = fields in template order
    selectedContainer.innerHTML = state.headerTemplate.map((key, idx) => {
        const field = HEADER_FIELDS.find(f => f.key === key);
        const label = field ? field.label : key;
        return `<span class="header-chip selected" draggable="true" data-idx="${idx}" onclick="removeHeaderField(${idx})">${label} <i class="fa-solid fa-xmark text-xs opacity-70"></i></span>`;
    }).join("");

    // Live preview from first record
    const sep = getHeaderSeparator();
    if (state.records.length > 0 && state.headerTemplate.length > 0) {
        const r = state.records[0];
        const tax = extractTaxonomy(r);
        previewEl.textContent = ">" + buildHeader(r, tax);
    } else if (state.headerTemplate.length === 0) {
        previewEl.textContent = i18nText("modal.header.empty", "(no fields selected)");
    } else {
        previewEl.innerHTML = "&gt;Species" + sep + "Accession";
    }

    // Duplicate check
    const duplicates = checkHeaderDuplicates();
    if (duplicates.length > 0 && state.headerTemplate.length > 0) {
        const totalDupes = duplicates.reduce((sum, d) => sum + d.accessions.length, 0);
        warningTextEl.textContent = i18nText("modal.header.duplicate.body", "{records} records produce {headers} duplicated header(s). Add more fields (e.g., Accession) to make headers unique.", {
            records: totalDupes,
            headers: duplicates.length,
        });
        warningEl.classList.remove("hidden");
        confirmBtn.disabled = true;
    } else if (state.headerTemplate.length === 0) {
        warningEl.classList.add("hidden");
        confirmBtn.disabled = true;
    } else {
        warningEl.classList.add("hidden");
        confirmBtn.disabled = false;
    }

    // Init drag-drop for selected chips
    initHeaderDragDrop();
}

function addHeaderField(key) {
    if (!state.headerTemplate.includes(key)) {
        state.headerTemplate.push(key);
        renderHeaderBuilder();
    }
}

function removeHeaderField(idx) {
    state.headerTemplate.splice(idx, 1);
    renderHeaderBuilder();
}

window.addHeaderField = addHeaderField;
window.removeHeaderField = removeHeaderField;

function initHeaderDragDrop() {
    const container = document.getElementById("headerSelectedFields");
    const chips = container.querySelectorAll(".header-chip.selected");
    let dragIdx = null;

    chips.forEach(chip => {
        chip.addEventListener("dragstart", (e) => {
            dragIdx = parseInt(chip.dataset.idx);
            chip.classList.add("dragging");
            e.dataTransfer.effectAllowed = "move";
        });

        chip.addEventListener("dragend", () => {
            chip.classList.remove("dragging");
            dragIdx = null;
        });

        chip.addEventListener("dragover", (e) => {
            e.preventDefault();
            e.dataTransfer.dropEffect = "move";
        });

        chip.addEventListener("drop", (e) => {
            e.preventDefault();
            const dropIdx = parseInt(chip.dataset.idx);
            if (dragIdx !== null && dragIdx !== dropIdx) {
                const item = state.headerTemplate.splice(dragIdx, 1)[0];
                state.headerTemplate.splice(dropIdx, 0, item);
                renderHeaderBuilder();
            }
        });
    });
}

function openHeaderBuilder() {
    renderHeaderBuilder();
    document.getElementById("headerModal").classList.remove("hidden");
}

// Update the inline header preview shown in the "Next Step" section
function updateHeaderPreviewInline() {
    const el = document.getElementById("headerPreviewInline");
    if (!el) return;
    const preview = document.getElementById("headerPreview");
    if (preview) {
        el.textContent = preview.textContent || preview.innerText;
        el.classList.remove("hidden");
    }
}

// Enable the Run Alignment button once the header has been confirmed
function enableRunButton() {
    const btn = document.getElementById("runAlignmentBtn");
    if (!btn) return;
    btn.disabled = false;
}

// ========================================================================
// FASTA Generation
// ========================================================================
function generateFastaFiles() {
    const fastaMap = new Map();
    const dataType = state.detectedDataType;

    for (const record of state.records) {
        const candidatesPerGene = new Map();
        const assignedTrnas = new Set();
        const tax = extractTaxonomy(record);
        const header = ">" + buildHeader(record, tax);
        const recordKey = record.accession || buildHeader(record, tax);

        for (const feat of record.features) {
            if (!state.selectedFeatureTypes.has(feat.type)) continue;

            let geneName = null;
            if (feat.type === "tRNA") {
                geneName = standardizeTrnaName(feat);
                if (geneName && assignedTrnas.has(geneName)) {
                    if (geneName === "tRNA-Ser1") geneName = "tRNA-Ser2";
                    else if (geneName === "tRNA-Leu1") geneName = "tRNA-Leu2";
                }
                if (geneName) assignedTrnas.add(geneName);
            } else {
                const rawGene = feat.qualifiers.gene || null;
                const rawProduct = feat.qualifiers.product || null;
                geneName = standardizeGeneName(rawGene, rawProduct, dataType);
            }

            if (!geneName || !state.selectedGenes.has(geneName)) continue;

            const seq = extractSequence(record.sequence, feat.locationStr);
            if (!seq) continue;

            const existing = candidatesPerGene.get(geneName) || [];
            existing.push({ header, seq, len: seq.length });
            candidatesPerGene.set(geneName, existing);
        }

        for (const [geneName, candidates] of candidatesPerGene) {
            const choiceKey = getDuplicateChoiceKey(recordKey, geneName);
            const longestIndex = candidates.reduce((bestIndex, candidate, index, all) => (
                candidate.len > all[bestIndex].len ? index : bestIndex
            ), 0);
            const selectedIndex = state.duplicateGeneChoices.has(choiceKey)
                ? state.duplicateGeneChoices.get(choiceKey)
                : longestIndex;
            const chosen = candidates[selectedIndex] || candidates[longestIndex];
            if (!chosen) continue;
            const existing = fastaMap.get(geneName) || "";
            fastaMap.set(geneName, existing + chosen.header + "\n" + chosen.seq + "\n");
        }
    }

    return fastaMap;
}

function downloadIndividualFasta() {
    const files = generateFastaFiles();
    if (files.size === 0) {
        alert("No sequences to download. Select genes and ensure records are loaded.");
        return;
    }
    for (const [name, content] of files) {
        const blob = new Blob([content], { type: "text/plain" });
        saveAs(blob, `${name}.fasta`);
    }
}

async function downloadZip() {
    const files = generateFastaFiles();
    if (files.size === 0) {
        alert("No sequences to download. Select genes and ensure records are loaded.");
        return;
    }
    const zip = new JSZip();
    for (const [name, content] of files) {
        zip.file(`${name}.fasta`, content);
    }
    const blob = await zip.generateAsync({ type: "blob" });
    saveAs(blob, `splace_${state.detectedDataType || "genes"}_genes.zip`);
}

// ========================================================================
// File Handling
// ========================================================================
const GENBANK_EXTENSIONS = /\.(gb|gbk|genbank)$/i;

function handleFiles(fileList) {
    const validFiles = [];
    for (const file of fileList) {
        if (!GENBANK_EXTENSIONS.test(file.name)) {
            continue;
        }
        validFiles.push(file);
    }

    if (validFiles.length === 0) return;

    showProgress("Loading files", `0 of ${validFiles.length}`);
    let loaded = 0;

    for (const file of validFiles) {
        const reader = new FileReader();
        reader.onload = (e) => {
            const text = e.target.result;
            try {
                const record = parseGenBank(text);
                record.source = "file";
                record.fileName = file.name;
                record.fileSize = file.size;
                record.fileDate = new Date(file.lastModified);
                if (!record.accession) record.accession = file.name.replace(GENBANK_EXTENSIONS, "");
                mergeRecords([record]);
            } catch {
                // Error parsing file, skipped
            }

            loaded++;
            updateProgress(loaded, validFiles.length, file.name);

            if (loaded === validFiles.length) {
                renderRecords();
                hideProgress(`Loaded ${validFiles.length} file${validFiles.length !== 1 ? 's' : ''}`);
            }
        };
        reader.readAsText(file);
    }
}

function handleElectronSelectedFiles(fileEntries) {
    const preparedFiles = (fileEntries || []).map((entry) => {
        if (entry instanceof File) {
            return entry;
        }
        if (!entry?.name || typeof entry.text !== "string") {
            return null;
        }
        return new File(
            [entry.text],
            entry.name,
            { type: "text/plain", lastModified: entry.lastModified || Date.now() }
        );
    }).filter(Boolean);

    if (!preparedFiles.length) {
        return;
    }

    handleFiles(preparedFiles);
}

function readEntryAsFile(fileEntry) {
    return new Promise((resolve, reject) => {
        fileEntry.file(resolve, reject);
    });
}

function readDirectoryEntries(dirReader) {
    return new Promise((resolve, reject) => {
        dirReader.readEntries(resolve, reject);
    });
}

async function traverseEntry(entry) {
    const files = [];
    if (entry.isFile) {
        if (GENBANK_EXTENSIONS.test(entry.name)) {
            try {
                const file = await readEntryAsFile(entry);
                files.push(file);
            } catch {
                // Error reading file entry, skipped
            }
        }
    } else if (entry.isDirectory) {
        const reader = entry.createReader();
        let batch;
        do {
            batch = await readDirectoryEntries(reader);
            for (const child of batch) {
                const childFiles = await traverseEntry(child);
                files.push(...childFiles);
            }
        } while (batch.length > 0);
    }
    return files;
}

async function handleDroppedItems(dataTransfer) {
    const items = dataTransfer.items;
    const allFiles = [];

    if (items) {
        const entries = [];
        for (let i = 0; i < items.length; i++) {
            const entry = items[i].webkitGetAsEntry ? items[i].webkitGetAsEntry() : null;
            if (entry) {
                entries.push(entry);
            }
        }

        if (entries.length > 0) {
            showProgress("Scanning folders", "Looking for GenBank files...");
            for (const entry of entries) {
                const files = await traverseEntry(entry);
                allFiles.push(...files);
            }
            if (allFiles.length > 0) {
                document.getElementById("progressModal").classList.add("hidden");
                handleFiles(allFiles);
            } else {
                hideProgress("No GenBank files found");
            }
            return;
        }
    }

    // Fallback: plain file drop
    handleFiles(dataTransfer.files);
}

// ========================================================================
// Event Listeners
// ========================================================================
document.addEventListener("DOMContentLoaded", () => {
    const dropZone = document.getElementById("dropZone");
    const fileInput = document.getElementById("fileInput");

    window.addEventListener("dragover", (e) => {
        e.preventDefault();
    });
    window.addEventListener("drop", (e) => {
        e.preventDefault();
    });

    dropZone.addEventListener("click", async () => {
        if (window.electronAPI?.selectGenbankInputs) {
            try {
                const selectedEntries = await window.electronAPI.selectGenbankInputs();
                if (Array.isArray(selectedEntries) && selectedEntries.length > 0) {
                    handleElectronSelectedFiles(selectedEntries);
                    return;
                }
            } catch (error) {
                console.error("Failed to open native GenBank picker", error);
            }
        }
        fileInput.value = "";
        fileInput.click();
    });
    dropZone.addEventListener("dragover", (e) => { e.preventDefault(); dropZone.classList.add("drop-active"); });
    dropZone.addEventListener("dragleave", () => dropZone.classList.remove("drop-active"));
    dropZone.addEventListener("drop", (e) => {
        e.preventDefault();
        dropZone.classList.remove("drop-active");
        handleDroppedItems(e.dataTransfer);
    });
    fileInput.addEventListener("change", (e) => {
        handleFiles(e.target.files);
        e.target.value = "";
    });

    resetStep1Inputs();

    document.getElementById("exampleBtn").addEventListener("click", fillExampleAccessions);
    document.getElementById("taxonomySearchExampleBtn").addEventListener("click", fillTaxonomySearchExample);
    document.getElementById("ncbiApiKeyBtn").addEventListener("click", openNcbiApiKeyModal);
    document.getElementById("ncbiApiKeyModalOverlay").addEventListener("click", closeNcbiApiKeyModal);
    document.getElementById("ncbiApiKeyModalClose").addEventListener("click", closeNcbiApiKeyModal);
    document.getElementById("saveNcbiApiKeyBtn").addEventListener("click", saveNcbiApiKeyFromInput);
    document.getElementById("clearNcbiApiKeyBtn").addEventListener("click", openClearNcbiApiKeyModal);
    document.getElementById("clearNcbiApiKeyModalOverlay").addEventListener("click", closeClearNcbiApiKeyModal);
    document.getElementById("clearNcbiApiKeyModalCancel").addEventListener("click", closeClearNcbiApiKeyModal);
    document.getElementById("clearNcbiApiKeyModalConfirm").addEventListener("click", clearNcbiApiKey);
    document.getElementById("taxonomyImportOverlay").addEventListener("click", () => closeTaxonomyImportModal(null));
    document.getElementById("taxonomyImportCancel").addEventListener("click", () => closeTaxonomyImportModal(null));
    document.getElementById("taxonomyImportSelectAll").addEventListener("click", () => setTaxonomyImportSelection("all"));
    document.getElementById("taxonomyImportSelectNone").addEventListener("click", () => setTaxonomyImportSelection("none"));
    document.getElementById("taxonomyImportDisableUnverified").addEventListener("click", () => setTaxonomyImportSelection("disable-unverified"));
    document.getElementById("taxonomyImportSelectComplete").addEventListener("click", () => setTaxonomyImportSelection("complete"));
    document.getElementById("taxonomyImportSelectPartial").addEventListener("click", () => setTaxonomyImportSelection("partial"));
    document.getElementById("taxonomyImportConfirm").addEventListener("click", () => {
        const selectedIds = state.taxonomyImportCandidates
            .filter((candidate) => candidate.selected)
            .map((candidate) => candidate.id);
        closeTaxonomyImportModal(selectedIds);
    });

    document.getElementById("accessionInput").addEventListener("input", validateAccessionInput);
    document.getElementById("taxonomySearchInput").addEventListener("input", validateTaxonomySearchInput);
    document.getElementById("taxonomyOrganelleSelect").addEventListener("change", validateTaxonomySearchInput);
    document.getElementById("taxonomyRefseqOnly").addEventListener("change", validateTaxonomySearchInput);
    document.getElementById("taxonomyCompleteCheckbox").addEventListener("change", validateTaxonomySearchInput);
    document.getElementById("taxonomyPartialCheckbox").addEventListener("change", validateTaxonomySearchInput);

    document.getElementById("fetchBtn").addEventListener("click", async () => {
        const raw = document.getElementById("accessionInput").value;
        const accessions = raw.split(/[\s,;\n]+/).filter(a => a.trim());
        if (accessions.length === 0) return;

        const btn = document.getElementById("fetchBtn");
        btn.disabled = true;
        renderFetchButton("step1.action.fetching", true);
        try {
            const records = await fetchMultiple(accessions);
            mergeRecords(records);
            renderRecords();
        } catch {
            // Fetch error, handled by progress UI
        }
        btn.disabled = false;
        renderFetchButton();
        validateAccessionInput();
    });

    document.getElementById("taxonomySearchBtn").addEventListener("click", async () => {
        const taxon = document.getElementById("taxonomySearchInput").value.trim();
        if (!taxon) return;

        const btn = document.getElementById("taxonomySearchBtn");
        const statusEl = document.getElementById("taxonomyFetchStatus");
        const options = {
            taxon,
            organelle: document.getElementById("taxonomyOrganelleSelect").value,
            refseqOnly: document.getElementById("taxonomyRefseqOnly").checked,
            completeness: getSelectedGenomeScopes(),
        };

        btn.disabled = true;
        renderTaxonomySearchButton("step1.action.searching", true);
        statusEl.classList.remove("hidden");
        statusEl.textContent = `Searching NCBI for ${taxon}...`;
        const t = window.SPLACE_I18N?.t || ((key) => key);
        showProgress(
            t("step1.progress.searchTitle"),
            t("step1.progress.searchSubtitle"),
            t("step1.progress.searchingFor", { name: taxon })
        );

        try {
            const ids = await searchNcbiGenomeIds(options);
            if (ids.length === 0) {
                statusEl.textContent = t("step1.status.noRecords", { name: taxon });
                hideProgress(t("step1.status.noRecords", { name: taxon }));
                return;
            }

            const summaries = await fetchNcbiSummaries(ids);
            const titleSummary = summarizeGenomeTitles(summaries);
            const reviewSummaryHtml = [
                t("step1.progress.foundSummary", titleSummary),
                `<div class="mt-1 text-xs text-gray-500">${t("step1.progress.fetchDetail")}</div>`
            ].join("");
            document.getElementById("progressDetail").innerHTML = reviewSummaryHtml;

            statusEl.textContent = t("step1.status.reviewing", { count: ids.length });
            document.getElementById("progressModal").classList.add("hidden");

            const candidates = buildTaxonomyImportCandidates(ids, summaries);
            const selectedIds = await openTaxonomyImportModal({
                taxon,
                summaryHtml: reviewSummaryHtml,
                candidates,
                searchDetails: buildTaxonomyReviewDetails(options),
            });

            if (!selectedIds || selectedIds.length === 0) {
                statusEl.textContent = t("step1.status.reviewCancelled");
                return;
            }

            statusEl.textContent = selectedIds.length === 1
                ? t("step1.status.foundImporting.one")
                : t("step1.status.foundImporting", { count: selectedIds.length });
            const progressItems = Object.fromEntries(candidates.map((candidate) => [String(candidate.id), candidate.title || candidate.accession || candidate.id]));
            const records = await fetchMultiple(selectedIds, {
                statusElementId: "taxonomyFetchStatus",
                progressTitle: t("step1.progress.fetchTitle"),
                progressDetail: reviewSummaryHtml,
                progressItems,
            });
            mergeRecords(records);
            renderRecords();
            statusEl.textContent = records.length === 1
                ? t("step1.status.imported.one", { name: taxon })
                : t("step1.status.imported", { count: records.length, name: taxon });
        } catch {
            statusEl.textContent = t("step1.status.searchFailed", { name: taxon });
        } finally {
            renderTaxonomySearchButton();
            validateTaxonomySearchInput();
            setTimeout(() => statusEl.classList.add("hidden"), 5000);
        }
    });

    renderTaxonomySearchButton();
    renderFetchButton();
    syncNcbiApiUi();
    validateAccessionInput();
    validateTaxonomySearchInput();
    updateProgressMeta();

    document.getElementById("clearRecords").addEventListener("click", () => {
        document.getElementById("clearModal").classList.remove("hidden");
    });

    document.getElementById("clearModalCancel").addEventListener("click", () => {
        document.getElementById("clearModal").classList.add("hidden");
    });

    document.getElementById("clearModalOverlay").addEventListener("click", () => {
        document.getElementById("clearModal").classList.add("hidden");
    });

    document.getElementById("clearModalConfirm").addEventListener("click", () => {
        document.getElementById("clearModal").classList.add("hidden");
        state.records = [];
        state.selectedGenes.clear();
        state.duplicateGeneChoices.clear();
        renderRecords();
    });

    document.getElementById("removeModalCancel").addEventListener("click", () => {
        document.getElementById("removeModal").classList.add("hidden");
    });

    document.getElementById("removeModalOverlay").addEventListener("click", () => {
        document.getElementById("removeModal").classList.add("hidden");
    });

    document.getElementById("removeModalConfirm").addEventListener("click", () => {
        document.getElementById("removeModal").classList.add("hidden");
        if (state.pendingRemoveIndex !== null && state.pendingRemoveIndex !== undefined) {
            state.records.splice(state.pendingRemoveIndex, 1);
            state.pendingRemoveIndex = null;
            renderRecords();
        }
    });

    // Edit modal listeners
    document.getElementById("editModalCancel").addEventListener("click", () => {
        document.getElementById("editModal").classList.add("hidden");
    });

    document.getElementById("editModalOverlay").addEventListener("click", () => {
        document.getElementById("editModal").classList.add("hidden");
    });

    document.getElementById("editSpeciesInput").addEventListener("input", function () {
        const parts = this.value.trim().split(/\s+/);
        document.getElementById("editGenusDisplay").value = parts[0] || "";
        // Enable/disable GBIF button based on species validation
        const err = validateSpeciesName(this.value);
        document.getElementById("editGbifBtn").disabled = !!err;
    });

    document.getElementById("editGbifBtn").addEventListener("click", async function () {
        const species = document.getElementById("editSpeciesInput").value.trim();
        if (!species) return;
        const statusEl = document.getElementById("editGbifStatus");
        statusEl.textContent = "Searching GBIF...";
        statusEl.className = "text-xs text-gray-500 mt-1";
        statusEl.classList.remove("hidden");
        this.disabled = true;
        try {
            const taxonomy = await fetchTaxonomyFromGBIF(species);
            if (taxonomy) {
                document.getElementById("editKingdomInput").value = taxonomy.kingdom || "";
                document.getElementById("editPhylumInput").value = taxonomy.phylum || "";
                document.getElementById("editClassInput").value = taxonomy.class || "";
                document.getElementById("editOrderInput").value = taxonomy.order || "";
                document.getElementById("editFamilyInput").value = taxonomy.family || "";
                // Store authorship for saving later
                this.dataset.fetchedAuthorship = taxonomy.authorship || "";
                statusEl.innerHTML = `Found via <span class="citation-ref" data-citation-key="datafishing"><em>dataFishing</em> (Rabelo et al. 2025)</span>`;
                window.initTooltips && window.initTooltips();
                statusEl.className = "text-xs text-green-600 mt-1";
                document.getElementById("editFamilyError").classList.add("hidden");
            } else {
                statusEl.textContent = "No results found in GBIF";
                statusEl.className = "text-xs text-orange-500 mt-1";
            }
        } catch (e) {
            statusEl.textContent = `Error: ${e.message}`;
            statusEl.className = "text-xs text-red-500 mt-1";
        }
        this.disabled = false;
    });

    document.getElementById("editModalSave").addEventListener("click", () => {
        const species = document.getElementById("editSpeciesInput").value.trim();
        const family = document.getElementById("editFamilyInput").value.trim();
        const kingdom = document.getElementById("editKingdomInput").value.trim();
        const phylum = document.getElementById("editPhylumInput").value.trim();
        const klass = document.getElementById("editClassInput").value.trim();
        const order = document.getElementById("editOrderInput").value.trim();

        const speciesError = validateSpeciesName(species);
        const familyError = validateFamilyName(family);

        const speciesErrEl = document.getElementById("editSpeciesError");
        const familyErrEl = document.getElementById("editFamilyError");

        if (speciesError) {
            speciesErrEl.textContent = speciesError;
            speciesErrEl.classList.remove("hidden");
        } else {
            speciesErrEl.classList.add("hidden");
        }

        if (familyError) {
            familyErrEl.textContent = familyError;
            familyErrEl.classList.remove("hidden");
        } else {
            familyErrEl.classList.add("hidden");
        }

        if (speciesError || familyError) return;

        const r = state.records[state.pendingEditIndex];
        r.editedOrganism = species;
        r.organism = species;
        if (family) r.editedFamily = family;

        // Save taxonomy ranks
        const fetchedAuthorship = document.getElementById("editGbifBtn").dataset.fetchedAuthorship || "";
        r.taxonomyRanks = {
            kingdom: kingdom,
            phylum: phylum,
            class: klass,
            order: order,
            family: family || (r.taxonomyRanks && r.taxonomyRanks.family) || "",
            authorship: fetchedAuthorship || (r.taxonomyRanks && r.taxonomyRanks.authorship) || "",
        };

        document.getElementById("editModal").classList.add("hidden");
        renderRecords();
    });

    document.getElementById("selectAllGenes").addEventListener("click", selectAllGenes);
    document.getElementById("selectNoneGenes").addEventListener("click", selectNoneGenes);
    document.getElementById("selectCompleteGenes").addEventListener("click", selectCompleteGenes);
    document.getElementById("selectDefaultGenes").addEventListener("click", selectDefaultGenes);
    document.getElementById("downloadBtn").addEventListener("click", () => openHeaderBuilder());

    // Fetch Taxonomy button
    document.getElementById("fetchTaxonomyBtn").addEventListener("click", fetchAllTaxonomy);

    // Column visibility eye toggles
    document.querySelectorAll(".col-eye-toggle[data-col-toggle]").forEach(btn => {
        btn.addEventListener("click", () => {
            const col = btn.dataset.colToggle;
            const isHidden = state.hiddenColumns.has(col);
            toggleColumn(col, isHidden);
        });
    });

    // Initialize column visibility checkboxes
    applyColumnVisibility();

    // Header builder modal listeners
    document.getElementById("headerModalCancel").addEventListener("click", () => {
        document.getElementById("headerModal").classList.add("hidden");
    });

    document.getElementById("headerModalOverlay").addEventListener("click", () => {
        document.getElementById("headerModal").classList.add("hidden");
    });

    document.getElementById("headerModalConfirm").addEventListener("click", () => {
        document.getElementById("headerModal").classList.add("hidden");
        if (window.isElectron) {
            updateHeaderPreviewInline();
            document.getElementById("panel6b").classList.remove("hidden");
            document.getElementById("panel6c").classList.remove("hidden");
            enableRunButton();
            updateStepBadges();
            return;
        }
        const format = document.querySelector('input[name="downloadFormat"]:checked').value;
        if (format === "individual") {
            downloadIndividualFasta();
        } else {
            downloadZip();
        }
    });

    // Re-render header builder when separator changes
    document.querySelectorAll('input[name="headerSeparator"]').forEach(radio => {
        radio.addEventListener("change", () => renderHeaderBuilder());
    });

    // Initialize Tippy on static elements
    window.initTooltips && window.initTooltips();

    // ----------------------------------------------------------------
    // Electron desktop mode
    // ----------------------------------------------------------------
    if (window.isElectron) {

        // In header modal: rename "Download" → "Confirm Header", hide format selector
        const headerModalConfirmBtn = document.getElementById("headerModalConfirm");
        headerModalConfirmBtn.innerHTML = '<i class="fa-solid fa-check mr-1"></i> Confirm Header';
        const dlFormatDiv = headerModalConfirmBtn.closest('.flex')?.previousElementSibling;
        if (dlFormatDiv && dlFormatDiv.querySelector('input[name="downloadFormat"]')) {
            dlFormatDiv.classList.add("hidden");
        }

        // Open header builder
        document.getElementById("openHeaderBuilderBtn").addEventListener("click", openHeaderBuilder);

        // Ask main process for CPU count; set default threads to all CPUs
        window.electronAPI.getCpuCount().then(cpus => {
            const threads = cpus;
            document.getElementById("mafftThreads").value = threads;
            const concurrent = Math.max(1, Math.floor(cpus / threads));
            document.getElementById("mafftConcurrencyInfo").textContent =
                `${cpus} logical CPUs — using ${threads} threads (${concurrent} job${concurrent !== 1 ? 's' : ''} in parallel)`;
        });
        document.getElementById("mafftThreads").addEventListener("input", () => {
            window.electronAPI.getCpuCount().then(cpus => {
                const threads = Math.max(1, parseInt(document.getElementById("mafftThreads").value) || 1);
                const concurrent = Math.max(1, Math.floor(cpus / threads));
                document.getElementById("mafftConcurrencyInfo").textContent =
                    `${cpus} logical CPUs — ${concurrent} job${concurrent !== 1 ? 's' : ''} × ${threads} thread${threads !== 1 ? 's' : ''}`;
            });
        });

        // trimAl mode toggle
        window.toggleTrimalManual = function (val) {
            document.getElementById("trimalManualFields").classList.toggle("hidden", val !== "manual");
        };

        // trimAl prompt buttons
        document.getElementById("trimalNoBtn").addEventListener("click", () => {
            document.getElementById("trimalPromptCard").classList.add("hidden");
            document.getElementById("alignmentResults").classList.remove("hidden");
            document.getElementById("showTrimalBtn").classList.remove("hidden");
            document.getElementById("alignmentResultsAccordion").scrollIntoView({ behavior: "smooth", block: "start" });
        });
        document.getElementById("trimalYesBtn").addEventListener("click", () => {
            document.getElementById("trimalPromptCard").classList.add("hidden");
            document.getElementById("panel6d").classList.remove("hidden");
            document.getElementById("panel6d").scrollIntoView({ behavior: "smooth", block: "start" });
        });
        document.getElementById("showTrimalBtn").addEventListener("click", () => {
            document.getElementById("showTrimalBtn").classList.add("hidden");
            document.getElementById("panel6d").classList.remove("hidden");
            document.getElementById("panel6d").scrollIntoView({ behavior: "smooth", block: "start" });
        });

        // Run trimAl
        document.getElementById("runTrimalBtn").addEventListener("click", () => {
            if (state.selectedGenes.size === 0) return;
            const mode = document.getElementById("trimalMode").value;
            const extra = document.getElementById("trimalExtra").value.trim();
            const params = [];
            if (mode === "manual") {
                const gt = document.getElementById("trimalGt").value;
                const st = document.getElementById("trimalSt").value;
                const cons = document.getElementById("trimalCons").value;
                const w = document.getElementById("trimalW").value;
                if (gt) params.push("-gt", gt);
                if (st) params.push("-st", st);
                if (cons) params.push("-cons", cons);
                if (w) params.push("-w", w);
            } else {
                params.push(mode);
            }
            if (extra) params.push(...extra.split(/\s+/).filter(Boolean));

            const markers = [...state.selectedGenes];
            openAnalysisModal("trimAl", markers);
            window.electronAPI.runTrimal({ markers, params });
        });

        // Analysis progress events
        window.electronAPI.onAnalysisProgress((data) => {
            updateAnalysisModal(data);
        });

        window.electronAPI.onAnalysisDone((result) => {
            finalizeAnalysisModal(result);
        });

        // Close modal
        document.getElementById("analysisModalClose").addEventListener("click", () => {
            document.getElementById("analysisModal").classList.add("hidden");
        });

        // Run alignment button
        document.getElementById("runAlignmentBtn").addEventListener("click", () => {
            if (state.selectedGenes.size === 0) { alert("No markers selected."); return; }
            const fastaFiles = generateFastaFiles();
            if (fastaFiles.size === 0) { alert("No sequences to align. Check records and gene selection."); return; }

            const method = document.getElementById("mafftMethod").value;
            const threads = parseInt(document.getElementById("mafftThreads").value) || 4;
            const maxIterate = parseInt(document.getElementById("mafftMaxIterate").value) || 0;
            const op = parseFloat(document.getElementById("mafftOp").value);
            const ep = parseFloat(document.getElementById("mafftEp").value);
            const adjustDir = document.getElementById("mafftAdjustDir").checked;
            const reorder = document.getElementById("mafftReorder").checked;
            const preserveCase = document.getElementById("mafftPreserveCase").checked;
            const extra = document.getElementById("mafftExtra").value.trim();

            const params = [...method.split(/\s+/).filter(Boolean)];
            if (maxIterate > 0) params.push("--maxiterate", String(maxIterate));
            if (!isNaN(op)) params.push("--op", String(op));
            if (!isNaN(ep)) params.push("--ep", String(ep));
            if (adjustDir) params.push("--adjustdirection");
            if (reorder) params.push("--reorder");
            if (preserveCase) params.push("--preservecase");
            if (extra) params.push(...extra.split(/\s+/).filter(Boolean));

            const files = {};
            for (const [name, content] of fastaFiles) files[name] = content;

            const markers = Object.keys(files);
            openAnalysisModal("mafft", markers);
            window.electronAPI.runAnalysis({ files, params, threads });
        });
    }
});

// ========================================================================
// Electron: Analysis Modal Helpers
// ========================================================================

let _analysisMarkers = [];

function openAnalysisModal(phase, markers) {
    _analysisMarkers = markers;
    const modal = document.getElementById("analysisModal");
    const icon = document.getElementById("analysisModalIcon");
    const title = document.getElementById("analysisModalTitle");
    const sub = document.getElementById("analysisModalSubtitle");
    const bar = document.getElementById("analysisModalProgressBar");
    const lbl = document.getElementById("analysisModalProgressLabel");
    const pct = document.getElementById("analysisModalProgressPct");
    const log = document.getElementById("analysisModalLog");
    const list = document.getElementById("analysisModalMarkerList");
    const close = document.getElementById("analysisModalClose");

    icon.className = "fa-solid fa-spinner fa-spin text-splace-blue-600 text-lg";
    title.textContent = phase === "mafft" ? "Running MAFFT Alignment…" : "Running trimAl…";
    sub.textContent = `${markers.length} marker${markers.length !== 1 ? 's' : ''} queued`;
    bar.style.width = "0%";
    lbl.textContent = `0 / ${markers.length} markers`;
    pct.textContent = "0%";
    log.textContent = "";
    close.classList.add("hidden");

    list.innerHTML = markers.map(m =>
        `<div id="modal-marker-${CSS.escape(m)}" class="flex items-center gap-2 px-2 py-1 rounded text-xs text-gray-600">
            <i class="fa-solid fa-clock text-gray-300 w-3"></i>
            <span class="font-mono">${m}</span>
        </div>`
    ).join("");

    modal.classList.remove("hidden");
}

window.toggleTaxonomyImportCandidate = toggleTaxonomyImportCandidate;

function updateAnalysisModal(data) {
    const log = document.getElementById("analysisModalLog");
    const bar = document.getElementById("analysisModalProgressBar");
    const lbl = document.getElementById("analysisModalProgressLabel");
    const pct = document.getElementById("analysisModalProgressPct");
    const sub = document.getElementById("analysisModalSubtitle");

    if (data.message) {
        log.textContent += data.message + "\n";
        log.scrollTop = log.scrollHeight;
    }

    if (data.total > 0) {
        const p = Math.round((data.done / data.total) * 100);
        bar.style.width = p + "%";
        lbl.textContent = `${data.done} / ${data.total} markers`;
        pct.textContent = p + "%";
    }

    // Update individual marker row
    if (data.marker) {
        const row = document.getElementById("modal-marker-" + CSS.escape(data.marker));
        if (row) {
            const icon = row.querySelector("i");
            if (data.status === "running") {
                icon.className = "fa-solid fa-spinner fa-spin text-splace-blue-500 w-3";
            } else if (data.status === "done") {
                icon.className = "fa-solid fa-check text-green-500 w-3";
                if (data.info) {
                    const span = row.querySelector("span");
                    span.innerHTML += ` <span class="text-gray-400">${data.info}</span>`;
                }
            } else if (data.status === "error") {
                icon.className = "fa-solid fa-xmark text-red-500 w-3";
            }
        }
    }

    if (data.phase) {
        sub.textContent = data.phase;
    }
}

function finalizeAnalysisModal(result) {
    const icon = document.getElementById("analysisModalIcon");
    const title = document.getElementById("analysisModalTitle");
    const close = document.getElementById("analysisModalClose");
    const bar = document.getElementById("analysisModalProgressBar");
    const lbl = document.getElementById("analysisModalProgressLabel");
    const pct = document.getElementById("analysisModalProgressPct");

    bar.style.width = "100%";
    lbl.textContent = `${result.aligned} / ${result.total} markers`;
    pct.textContent = "100%";

    if (result.phase === "mafft") {
        icon.className = "fa-solid fa-check text-green-600 text-lg";
        title.textContent = `Alignment complete — ${result.aligned}/${result.total} markers aligned`;
        close.classList.remove("hidden");
        // After modal is closed, show trimAl prompt and partial results
        close.onclick = () => {
            document.getElementById("analysisModal").classList.add("hidden");
            if (result.markerResults) {
                window._alignmentOutputDir = result.outputDir;
                renderAlignmentResults(result.markerResults, false);
                renderConcatSection(result.markerResults, false);
            }
            document.getElementById("trimalPromptCard").classList.remove("hidden");
            document.getElementById("trimalPromptCard").scrollIntoView({ behavior: "smooth", block: "start" });
            setBadge("step6c-badge", "done");
            updateStepBadges();
        };
    } else if (result.phase === "trimal") {
        icon.className = "fa-solid fa-check text-green-600 text-lg";
        title.textContent = `Trimming complete — ${result.aligned}/${result.total} markers trimmed`;
        close.classList.remove("hidden");
        close.onclick = () => {
            document.getElementById("analysisModal").classList.add("hidden");
            if (result.markerResults) {
                renderAlignmentResults(result.markerResults, true);
                renderConcatSection(result.markerResults, true);
            }
            document.getElementById("alignmentResults").classList.remove("hidden");
            document.getElementById("alignmentResults").scrollIntoView({ behavior: "smooth", block: "start" });
            setBadge("step6d-badge", "done");
            updateStepBadges();
        };
    }
}

// ========================================================================
// Electron: Alignment Visualization
// ========================================================================

function parseFasta(text) {
    const seqs = [];
    let cur = null;
    for (const line of text.split(/\r?\n/)) {
        if (line.startsWith(">")) {
            if (cur) seqs.push(cur);
            cur = { name: line.slice(1).trim(), seq: "" };
        } else if (cur) {
            cur.seq += line.trim();
        }
    }
    if (cur) seqs.push(cur);
    return seqs;
}

const NT_COLORS = { A: "#1AC253", T: "#E11218", U: "#E11218", G: "#ECA918", C: "#4272EE", "-": "#d1d5db" };
const AA_COLORS = {
    A: "#8dd3c7", C: "#ffffb3", D: "#fb8072", E: "#fb8072", F: "#80b1d3",
    G: "#b3de69", H: "#bebada", I: "#80b1d3", K: "#fdb462", L: "#80b1d3",
    M: "#8dd3c7", N: "#fccde5", P: "#d9d9d9", Q: "#fccde5", R: "#fdb462",
    S: "#b3de69", T: "#b3de69", V: "#80b1d3", W: "#bc80bd", Y: "#bc80bd",
    X: "#cbd5e1", "*": "#dc2626", "-": "#e2e8f0", "?": "#cbd5e1"
};
const CODON_TABLE = {
    TTT: "F", TTC: "F", TTA: "L", TTG: "L", TCT: "S", TCC: "S", TCA: "S", TCG: "S",
    TAT: "Y", TAC: "Y", TAA: "*", TAG: "*", TGT: "C", TGC: "C", TGA: "*", TGG: "W",
    CTT: "L", CTC: "L", CTA: "L", CTG: "L", CCT: "P", CCC: "P", CCA: "P", CCG: "P",
    CAT: "H", CAC: "H", CAA: "Q", CAG: "Q", CGT: "R", CGC: "R", CGA: "R", CGG: "R",
    ATT: "I", ATC: "I", ATA: "I", ATG: "M", ACT: "T", ACC: "T", ACA: "T", ACG: "T",
    AAT: "N", AAC: "N", AAA: "K", AAG: "K", AGT: "S", AGC: "S", AGA: "R", AGG: "R",
    GTT: "V", GTC: "V", GTA: "V", GTG: "V", GCT: "A", GCC: "A", GCA: "A", GCG: "A",
    GAT: "D", GAC: "D", GAA: "E", GAG: "E", GGT: "G", GGC: "G", GGA: "G", GGG: "G"
};
const GENETIC_CODE_TABLES = {
    "1": CODON_TABLE,
    "2": { ...CODON_TABLE, ATA: "M", TGA: "W", AGA: "*", AGG: "*" },
    "5": { ...CODON_TABLE, ATA: "M", TGA: "W", AGA: "S", AGG: "S" },
    "11": CODON_TABLE,
};

function hexToRgb(hex) {
    return {
        r: parseInt(hex.slice(1, 3), 16),
        g: parseInt(hex.slice(3, 5), 16),
        b: parseInt(hex.slice(5, 7), 16)
    };
}

function translateAlignedCodons(sequence, geneticCode = "1") {
    const normalized = (sequence || "").toUpperCase().replace(/U/g, "T");
    const padded = normalized + "-".repeat((3 - (normalized.length % 3 || 3)) % 3);
    const codonTable = GENETIC_CODE_TABLES[geneticCode] || CODON_TABLE;
    let protein = "";
    let stopCount = 0;
    for (let idx = 0; idx < padded.length; idx += 3) {
        const codon = padded.slice(idx, idx + 3);
        if (codon === "---") {
            protein += "-";
            continue;
        }
        if (/[^ACGT?\-]/.test(codon) || codon.includes("-") || codon.includes("?")) {
            protein += "X";
            continue;
        }
        const aa = codonTable[codon] || "X";
        if (aa === "*") stopCount++;
        protein += aa;
    }
    return { protein, stopCount };
}

function translateAlignmentSequences(sequences, geneticCode = "1") {
    return (sequences || []).map((seq) => {
        const translated = translateAlignedCodons(seq.seq, geneticCode);
        return {
            name: seq.name,
            seq: translated.protein,
            stopCount: translated.stopCount,
        };
    });
}

function buildSequenceSnippet(sequence, palette, maxChars = 72) {
    const slice = (sequence || "").slice(0, maxChars);
    const html = [];
    for (const residue of slice) {
        const color = palette[residue] || "#cbd5e1";
        const fg = residue === "*" ? "#fff" : "rgba(15,23,42,0.78)";
        html.push(`<span class="alignment-residue" style="background:${color};color:${fg};">${residue}</span>`);
    }
    if ((sequence || "").length > maxChars) html.push('<span class="text-gray-400">...</span>');
    return html.join("");
}

function renderAlignmentPreview(container, sequences, mode) {
    if (!container) return;
    const rows = sequences || [];
    if (!rows.length) {
        container.innerHTML = `<div class="alignment-preview"><div class="alignment-preview-empty">${i18nText("step5.results.preview.empty", "No sequences available for this view.")}</div></div>`;
        return;
    }

    const palette = mode === "aa" ? AA_COLORS : NT_COLORS;
    const totalStops = mode === "aa" ? rows.reduce((acc, seq) => acc + (seq.stopCount || 0), 0) : 0;
    const summaryKey = mode === "aa" ? "step5.results.preview.summaryStops" : "step5.results.preview.summary";
    const summaryFallback = mode === "aa"
        ? "{seqs} sequences rendered; {stops} possible stop codon(s)"
        : "{seqs} sequences rendered";
    const summary = i18nText(summaryKey, summaryFallback, { seqs: rows.length, stops: totalStops });
    const previewRows = rows.slice(0, 12).map((seq) => {
        const length = seq.seq.length;
        const stopDisplay = mode === "aa"
            ? String(seq.stopCount || 0)
            : i18nText("step5.results.preview.none", "none");
        return `
            <div class="alignment-preview-row">
                <div class="alignment-preview-name" title="${seq.name}">${seq.name}</div>
                <div class="alignment-preview-metric">${length}</div>
                <div class="alignment-preview-metric">${stopDisplay}</div>
                <div class="alignment-preview-seq">${buildSequenceSnippet(seq.seq, palette)}</div>
            </div>`;
    }).join("");

    container.innerHTML = `
        <div class="alignment-preview">
            <div class="alignment-preview-summary">
                <strong>${summary}</strong>
            </div>
            <div class="alignment-preview-table">
                <div class="alignment-preview-row alignment-preview-head">
                    <div>${i18nText("step5.results.preview.name", "Sequence")}</div>
                    <div>${i18nText("step5.results.preview.length", "Length")}</div>
                    <div>${i18nText("step5.results.preview.stops", "Stops")}</div>
                    <div>${i18nText("step5.results.preview.snippet", "Preview")}</div>
                </div>
                ${previewRows}
            </div>
        </div>`;
}

function getAlignmentViewOptions(mode) {
    if (mode === "aa") {
        return {
            colors: AA_COLORS,
            labels: {
                zoom: i18nText("step5.results.zoom", "ZOOM"),
                consensus: i18nText("step5.results.consensus", "CONSENSUS"),
                conservation: i18nText("step5.results.conservation", "CONSERVATION"),
            },
            legend: [
                [i18nText("step5.results.aaLegend", "AA"), "#80b1d3"],
                [i18nText("step5.results.stopLegend", "Stop"), "#dc2626"],
                [i18nText("step5.results.gapLegend", "Gap"), "#e2e8f0"],
                [i18nText("step5.results.unknownLegend", "Unknown"), "#cbd5e1"],
            ],
        };
    }
    return {
        colors: NT_COLORS,
        labels: {
            zoom: i18nText("step5.results.zoom", "ZOOM"),
            consensus: i18nText("step5.results.consensus", "CONSENSUS"),
            conservation: i18nText("step5.results.conservation", "CONSERVATION"),
        },
        legend: [["A", "#22c55e"], ["T/U", "#ef4444"], ["G", "#f59e0b"], ["C", "#3b82f6"], [i18nText("step5.results.gapLegend", "Gap"), "#e2e8f0"]],
    };
}

function renderAlignmentTab(markerData, tab, mode) {
    const contentMap = {
        raw: markerData.rawContent,
        aligned: markerData.alignedContent,
        trimmed: markerData.trimmedContent,
    };
    const content = contentMap[tab];
    const wrap = document.getElementById(`vizwrap-${markerData.id}-${tab}`);
    const preview = document.getElementById(`preview-${markerData.id}-${tab}`);
    if (!wrap || !preview || !content) return;
    const sequences = parseFasta(content);
    const geneticCode = document.getElementById(`genetic-code-${markerData.id}-${tab}`)?.value || "1";
    const displaySequences = mode === "aa" ? translateAlignmentSequences(sequences, geneticCode) : sequences;
    renderAlignmentPlotly(wrap, displaySequences, getAlignmentViewOptions(mode));
    renderAlignmentPreview(preview, displaySequences, mode);
    wrap.dataset.rendered = "1";
    wrap.dataset.view = mode;
}

function renderAlignmentCanvas(canvas, sequences) {
    if (!sequences || sequences.length === 0) return;
    const nSeq = sequences.length;
    const seqLen = sequences[0].seq.length || 1;

    const container = canvas.parentElement;
    const W = Math.min(seqLen, container ? container.clientWidth || 560 : 560);
    const cellH = Math.max(2, Math.min(8, Math.floor(160 / nSeq)));
    const H = nSeq * cellH;

    canvas.width = W;
    canvas.height = H;
    canvas.style.width = "100%";
    canvas.style.height = H + "px";
    canvas.style.imageRendering = "pixelated";

    const ctx = canvas.getContext("2d");
    const img = ctx.createImageData(W, H);
    const data = img.data;

    const colorCache = {};
    const getColor = (ch) => {
        if (!colorCache[ch]) colorCache[ch] = hexToRgb(NT_COLORS[ch] || "#9ca3af");
        return colorCache[ch];
    };

    for (let si = 0; si < nSeq; si++) {
        const seq = sequences[si].seq.toUpperCase();
        const rowY = si * cellH;
        for (let px = 0; px < W; px++) {
            const pi = Math.floor(px * seqLen / W);
            const c = getColor(seq[pi] || "");
            for (let dy = 0; dy < cellH; dy++) {
                const idx = ((rowY + dy) * W + px) * 4;
                data[idx] = c.r;
                data[idx + 1] = c.g;
                data[idx + 2] = c.b;
                data[idx + 3] = 255;
            }
        }
    }
    ctx.putImageData(img, 0, 0);
}

// MSA alignment viewer — interactive, Dash Bio-style with zoom/pan/consensus
function renderAlignmentPlotly(wrapper, sequences, options = {}) {
    if (!sequences || sequences.length === 0) return;
    wrapper.innerHTML = "";

    const nSeq = sequences.length;
    const seqLen = sequences[0].seq.length || 1;

    // ── constants ──────────────────────────────────────────────────────────
    const RULER_H = 24;
    const CONS_SEQ_H = Math.max(14, Math.min(22, Math.floor(480 / nSeq)));
    const CONS_BAR_H = 52;
    const NAME_W = 190;
    const cellH = Math.max(14, Math.min(22, Math.floor(480 / nSeq)));

    const palette = options.colors || {
        A: "#22c55e", T: "#ef4444", U: "#ef4444",
        G: "#f59e0b", C: "#3b82f6",
        "-": "#e2e8f0", "?": "#e2e8f0", N: "#cbd5e1",
    };
    const legendItems = options.legend || [["A", "#22c55e"], ["T/U", "#ef4444"], ["G", "#f59e0b"], ["C", "#3b82f6"], ["Gap", "#e2e8f0"]];
    const labels = options.labels || { zoom: "ZOOM", consensus: "CONSENSUS", conservation: "CONSERVATION" };
    const DEFAULT_COLOR = "#cbd5e1";

    // ── precompute consensus + conservation (independent of zoom) ──────────
    const consensusSeq = [];
    const consScores = [];
    for (let pi = 0; pi < seqLen; pi++) {
        const counts = {};
        let total = 0;
        for (let si = 0; si < nSeq; si++) {
            const ch = (sequences[si].seq[pi] || "-").toUpperCase();
            if (ch !== "-" && ch !== "?") { counts[ch] = (counts[ch] || 0) + 1; total++; }
        }
        const entries = Object.entries(counts);
        consensusSeq.push(entries.length ? entries.sort((a, b) => b[1] - a[1])[0][0] : "-");
        consScores.push(total > 0 ? Math.max(...Object.values(counts)) / nSeq : 0);
    }

    // ── mutable zoom state ─────────────────────────────────────────────────
    let cellW = seqLen <= 100 ? 14 : seqLen <= 400 ? 10 : seqLen <= 1000 ? 6 : seqLen <= 3000 ? 3 : 2;

    // ── outer container ────────────────────────────────────────────────────
    const container = document.createElement("div");
    container.style.cssText = "background:#f9fafb;border-radius:8px;overflow:hidden;display:flex;flex-direction:column;max-height:660px;user-select:none;font-family:system-ui,sans-serif;";
    wrapper.appendChild(container);

    // ── toolbar ────────────────────────────────────────────────────────────
    const toolbar = document.createElement("div");
    toolbar.style.cssText = "display:flex;align-items:center;gap:6px;padding:5px 10px;background:#f3f4f6;border-bottom:1px solid #e5e7eb;flex-shrink:0;";
    container.appendChild(toolbar);

    const zoomLabel = document.createElement("span");
    zoomLabel.style.cssText = "font-size:10px;color:#6b7280;font-weight:700;letter-spacing:.05em;";
    zoomLabel.textContent = labels.zoom;
    toolbar.appendChild(zoomLabel);

    function makeBtn(html, title, onClick) {
        const btn = document.createElement("button");
        btn.innerHTML = html; btn.title = title;
        btn.style.cssText = "padding:1px 8px;border:1px solid #d1d5db;border-radius:5px;background:#fff;color:#374151;font-size:14px;cursor:pointer;font-weight:700;line-height:1.5;transition:background .1s;";
        btn.onmouseenter = () => { btn.style.background = "#f9fafb"; };
        btn.onmouseleave = () => { btn.style.background = "#fff"; };
        btn.onclick = onClick;
        return btn;
    }

    const zoomInfo = document.createElement("span");
    zoomInfo.style.cssText = "font-size:11px;color:#9ca3af;min-width:52px;";

    const zoomOutBtn = makeBtn("−", "Zoom out (or Ctrl+scroll)", () => { if (cellW > 1) { cellW = Math.max(1, Math.floor(cellW * 0.65)); redraw(); } });
    const zoomInBtn = makeBtn("+", "Zoom in (or Ctrl+scroll)", () => { cellW = Math.min(28, Math.ceil(cellW * 1.55)); redraw(); });
    const zoomFitBtn = makeBtn("Fit", "Fit to view", () => {
        const avail = seqOuter.clientWidth || container.clientWidth - NAME_W;
        cellW = Math.max(1, Math.floor(avail / seqLen));
        redraw();
    });

    toolbar.appendChild(zoomOutBtn);
    toolbar.appendChild(zoomInBtn);
    toolbar.appendChild(zoomFitBtn);
    toolbar.appendChild(zoomInfo);

    // Legend
    const legend = document.createElement("div");
    legend.style.cssText = "display:flex;align-items:center;gap:10px;margin-left:auto;";
    legendItems.forEach(([lbl, col]) => {
        const chip = document.createElement("span");
        chip.style.cssText = "display:flex;align-items:center;gap:3px;font-size:11px;color:#4b5563;";
        chip.innerHTML = `<span style="display:inline-block;width:10px;height:10px;border-radius:2px;background:${col};border:1px solid rgba(0,0,0,.12);"></span>${lbl}`;
        legend.appendChild(chip);
    });
    toolbar.appendChild(legend);

    // ── ruler row ──────────────────────────────────────────────────────────
    const rulerRow = document.createElement("div");
    rulerRow.style.cssText = "display:flex;flex-shrink:0;";
    container.appendChild(rulerRow);

    const rulerSpacer = document.createElement("div");
    rulerSpacer.style.cssText = `width:${NAME_W}px;flex-shrink:0;height:${RULER_H}px;background:#f3f4f6;border-right:2px solid #d1d5db;`;
    rulerRow.appendChild(rulerSpacer);

    const rulerOuter = document.createElement("div");
    rulerOuter.style.cssText = "overflow:hidden;flex:1;";
    rulerRow.appendChild(rulerOuter);

    const rulerCanvas = document.createElement("canvas");
    rulerCanvas.height = RULER_H;
    rulerCanvas.style.cssText = "display:block;";
    rulerOuter.appendChild(rulerCanvas);

    // ── body row ───────────────────────────────────────────────────────────
    const bodyRow = document.createElement("div");
    bodyRow.style.cssText = "display:flex;flex:1;overflow-y:auto;min-height:0;";
    container.appendChild(bodyRow);

    const nameDiv = document.createElement("div");
    nameDiv.style.cssText = `width:${NAME_W}px;flex-shrink:0;overflow:hidden;background:#f3f4f6;border-right:2px solid #d1d5db;`;
    bodyRow.appendChild(nameDiv);

    const nameCanvas = document.createElement("canvas");
    nameCanvas.width = NAME_W;
    nameCanvas.style.cssText = "display:block;";
    nameDiv.appendChild(nameCanvas);

    const seqOuter = document.createElement("div");
    seqOuter.style.cssText = "overflow-x:auto;overflow-y:hidden;flex:1;cursor:grab;";
    bodyRow.appendChild(seqOuter);

    const mainCanvas = document.createElement("canvas");
    mainCanvas.style.cssText = "display:block;";
    seqOuter.appendChild(mainCanvas);

    // ── consensus row ──────────────────────────────────────────────────────
    const consSeqRow = document.createElement("div");
    consSeqRow.style.cssText = "display:flex;flex-shrink:0;border-top:2px solid #cbd5e1;";
    container.appendChild(consSeqRow);

    const consSeqNameDiv = document.createElement("div");
    consSeqNameDiv.style.cssText = `width:${NAME_W}px;flex-shrink:0;height:${CONS_SEQ_H}px;background:#dbeafe;border-right:2px solid #d1d5db;display:flex;align-items:center;justify-content:flex-end;padding-right:8px;`;
    consSeqNameDiv.innerHTML = `<span style="font-size:10px;color:#1e40af;font-weight:700;letter-spacing:.05em;">${labels.consensus}</span>`;
    consSeqRow.appendChild(consSeqNameDiv);

    const consSeqOuter = document.createElement("div");
    consSeqOuter.style.cssText = "overflow:hidden;flex:1;";
    consSeqRow.appendChild(consSeqOuter);

    const consSeqCanvas = document.createElement("canvas");
    consSeqCanvas.height = CONS_SEQ_H;
    consSeqCanvas.style.cssText = "display:block;";
    consSeqOuter.appendChild(consSeqCanvas);

    // ── conservation bar row ───────────────────────────────────────────────
    const consBarRow = document.createElement("div");
    consBarRow.style.cssText = "display:flex;flex-shrink:0;border-top:1px solid #e5e7eb;";
    container.appendChild(consBarRow);

    const consBarNameDiv = document.createElement("div");
    consBarNameDiv.style.cssText = `width:${NAME_W}px;flex-shrink:0;height:${CONS_BAR_H}px;background:#f3f4f6;border-right:2px solid #d1d5db;display:flex;align-items:center;justify-content:flex-end;padding-right:8px;`;
    consBarNameDiv.innerHTML = `<span style="font-size:10px;color:#6b7280;font-weight:600;letter-spacing:.04em;">${labels.conservation}</span>`;
    consBarRow.appendChild(consBarNameDiv);

    const consBarOuter = document.createElement("div");
    consBarOuter.style.cssText = "overflow:hidden;flex:1;";
    consBarRow.appendChild(consBarOuter);

    const consBarCanvas = document.createElement("canvas");
    consBarCanvas.height = CONS_BAR_H;
    consBarCanvas.style.cssText = "display:block;";
    consBarOuter.appendChild(consBarCanvas);

    // ── redraw function ────────────────────────────────────────────────────
    function redraw() {
        const W = seqLen * cellW;
        const H = nSeq * cellH;
        const gapX = cellW > 3 ? 1 : 0;
        const gapY = 1;
        const showLetter = cellW >= 6;
        const fontSize = showLetter ? Math.min(cellW - 1, cellH - 2, 12) : 0;

        zoomInfo.textContent = `${cellW}px/nt`;

        // --- Ruler ---
        rulerCanvas.width = W;
        const rCtx = rulerCanvas.getContext("2d");
        rCtx.fillStyle = "#f3f4f6";
        rCtx.fillRect(0, 0, W, RULER_H);
        rCtx.font = "10px monospace";
        // choose tick interval so labels don't overlap (each label ~30px)
        const minTickPx = 35;
        const rawStep = Math.ceil(minTickPx / cellW);
        const niceSteps = [1, 2, 5, 10, 25, 50, 100, 200, 500, 1000, 2000, 5000];
        const tickStep = niceSteps.find(s => s >= rawStep) || rawStep;
        for (let p = 0; p < seqLen; p++) {
            if (p === 0 || (p + 1) % tickStep === 0) {
                const x = p * cellW;
                rCtx.fillStyle = "#9ca3af";
                rCtx.fillRect(x, RULER_H - 6, 1, 6);
                rCtx.fillStyle = "#374151";
                rCtx.fillText(p + 1, x + 2, RULER_H - 9);
            }
        }

        // --- Names ---
        nameCanvas.height = H;
        const nCtx = nameCanvas.getContext("2d");
        nCtx.fillStyle = "#f3f4f6";
        nCtx.fillRect(0, 0, NAME_W, H);
        const nameFontSize = Math.min(12, cellH - 2);
        nCtx.font = `${nameFontSize}px system-ui,sans-serif`;
        nCtx.textBaseline = "middle";
        nCtx.textAlign = "right";
        for (let si = 0; si < nSeq; si++) {
            const y = si * cellH;
            if (si % 2 === 0) { nCtx.fillStyle = "#eef0f5"; nCtx.fillRect(0, y, NAME_W, cellH); }
            nCtx.fillStyle = "#374151";
            const name = sequences[si].name;
            const maxChars = Math.floor((NAME_W - 10) / (nameFontSize * 0.58));
            nCtx.fillText(name.length > maxChars ? name.slice(0, maxChars - 1) + "…" : name, NAME_W - 6, y + cellH / 2);
        }

        // --- Sequences ---
        mainCanvas.width = W;
        mainCanvas.height = H;
        const mCtx = mainCanvas.getContext("2d");
        if (showLetter) {
            mCtx.font = `bold ${fontSize}px monospace`;
            mCtx.textAlign = "center";
            mCtx.textBaseline = "middle";
        }
        for (let si = 0; si < nSeq; si++) {
            const y = si * cellH;
            if (si % 2 === 0) {
                mCtx.fillStyle = "rgba(0,0,0,0.025)";
                mCtx.fillRect(0, y, W, cellH);
            }
            for (let pi = 0; pi < seqLen; pi++) {
                const ch = (sequences[si].seq[pi] || "-").toUpperCase();
                mCtx.fillStyle = palette[ch] || DEFAULT_COLOR;
                mCtx.fillRect(pi * cellW, y + gapY, cellW - gapX, cellH - gapY * 2);
                if (showLetter) {
                    mCtx.fillStyle = (ch === "-" || ch === "?") ? "#94a3b8" : "rgba(0,0,0,0.65)";
                    mCtx.fillText(ch, pi * cellW + cellW / 2, y + cellH / 2);
                }
            }
        }

        // --- Consensus ---
        consSeqCanvas.width = W;
        const csCtx = consSeqCanvas.getContext("2d");
        csCtx.fillStyle = "#dbeafe";
        csCtx.fillRect(0, 0, W, CONS_SEQ_H);
        if (showLetter) {
            csCtx.font = `bold ${fontSize}px monospace`;
            csCtx.textAlign = "center";
            csCtx.textBaseline = "middle";
        }
        for (let pi = 0; pi < seqLen; pi++) {
            const ch = consensusSeq[pi];
            csCtx.fillStyle = palette[ch] || DEFAULT_COLOR;
            csCtx.fillRect(pi * cellW, gapY, cellW - gapX, CONS_SEQ_H - gapY * 2);
            if (showLetter) {
                csCtx.fillStyle = (ch === "-" || ch === "?") ? "#94a3b8" : "rgba(0,0,0,0.65)";
                csCtx.fillText(ch, pi * cellW + cellW / 2, CONS_SEQ_H / 2);
            }
        }

        // --- Conservation bar ---
        consBarCanvas.width = W;
        const cbCtx = consBarCanvas.getContext("2d");
        cbCtx.fillStyle = "#f9fafb";
        cbCtx.fillRect(0, 0, W, CONS_BAR_H);
        for (let pi = 0; pi < seqLen; pi++) {
            const score = consScores[pi];
            const barH = Math.round(score * (CONS_BAR_H - 6));
            cbCtx.fillStyle = score >= 0.8 ? "#16a34a" : score >= 0.5 ? "#d97706" : "#dc2626";
            cbCtx.fillRect(pi * cellW, CONS_BAR_H - barH - 2, Math.max(1, cellW - gapX), barH);
        }
    }

    // Initial draw
    redraw();

    // ── scroll sync ────────────────────────────────────────────────────────
    seqOuter.addEventListener("scroll", () => {
        rulerOuter.scrollLeft = seqOuter.scrollLeft;
        consSeqOuter.scrollLeft = seqOuter.scrollLeft;
        consBarOuter.scrollLeft = seqOuter.scrollLeft;
    });
    bodyRow.addEventListener("scroll", () => {
        nameDiv.scrollTop = bodyRow.scrollTop;
    });

    // ── Ctrl+wheel zoom ────────────────────────────────────────────────────
    seqOuter.addEventListener("wheel", (e) => {
        if (!e.ctrlKey && !e.metaKey) return;
        e.preventDefault();
        const prevW = cellW;
        if (e.deltaY < 0) cellW = Math.min(28, Math.ceil(cellW * 1.4));
        else cellW = Math.max(1, Math.floor(cellW * 0.72));
        if (cellW !== prevW) {
            const mouseRatio = (seqOuter.scrollLeft + e.offsetX) / (seqLen * prevW);
            redraw();
            seqOuter.scrollLeft = mouseRatio * seqLen * cellW - e.offsetX;
        }
    }, { passive: false });

    // ── click-drag pan ─────────────────────────────────────────────────────
    let dragging = false, dragX = 0, dragScroll = 0;
    mainCanvas.addEventListener("mousedown", (e) => {
        dragging = true; dragX = e.clientX; dragScroll = seqOuter.scrollLeft;
        mainCanvas.style.cursor = "grabbing";
        e.preventDefault();
    });
    document.addEventListener("mousemove", (e) => {
        if (!dragging) return;
        seqOuter.scrollLeft = dragScroll - (e.clientX - dragX);
    });
    document.addEventListener("mouseup", () => {
        if (!dragging) return;
        dragging = false;
        mainCanvas.style.cursor = "grab";
    });

    // ── tooltip on hover ───────────────────────────────────────────────────
    const tip = document.createElement("div");
    tip.style.cssText = "position:fixed;background:#1e293b;color:#f1f5f9;font:11px monospace;padding:4px 8px;border-radius:5px;pointer-events:none;z-index:9999;white-space:nowrap;display:none;box-shadow:0 2px 8px rgba(0,0,0,.4);";
    document.body.appendChild(tip);
    mainCanvas.addEventListener("mousemove", (e) => {
        if (dragging) { tip.style.display = "none"; return; }
        const rect = mainCanvas.getBoundingClientRect();
        const pi = Math.floor((e.clientX - rect.left + seqOuter.scrollLeft) / cellW);
        const si = Math.floor((e.clientY - rect.top + bodyRow.scrollTop) / cellH);
        if (pi >= 0 && pi < seqLen && si >= 0 && si < nSeq) {
            const ch = (sequences[si].seq[pi] || "-").toUpperCase();
            tip.textContent = `${sequences[si].name} | pos ${pi + 1} | ${ch}`;
            tip.style.display = "block";
            tip.style.left = (e.clientX + 14) + "px";
            tip.style.top = (e.clientY - 22) + "px";
        } else {
            tip.style.display = "none";
        }
    });
    mainCanvas.addEventListener("mouseleave", () => { tip.style.display = "none"; });
    wrapper.addEventListener("remove", () => tip.remove()); // cleanup
    // cleanup tip when tab changes / re-renders
    const obs = new MutationObserver(() => { if (!document.body.contains(wrapper)) tip.remove(); });
    obs.observe(document.body, { childList: true, subtree: true });
}

// Canvas-based alignment viewer (kept as fallback)
function renderAlignmentCanvasEnhanced(wrapper, sequences) {
    if (!sequences || sequences.length === 0) return;
    wrapper.innerHTML = "";

    const nSeq = sequences.length;
    const seqLen = sequences[0].seq.length || 1;
    const W = Math.min(seqLen, 2000);
    const cellH = Math.max(4, Math.min(10, Math.floor(300 / nSeq)));
    const RULER_H = 28;
    const CONS_H = 64;
    const LABEL_W = 0; // no per-row labels (overview mode)

    // --- scrollable outer container ---
    const outer = document.createElement("div");
    outer.style.cssText = "overflow-x:auto;overflow-y:auto;max-height:480px;background:#f9fafb;border-radius:8px;";
    wrapper.appendChild(outer);

    // --- inner container (natural width) ---
    const inner = document.createElement("div");
    inner.style.cssText = `width:${W}px;position:relative;`;
    outer.appendChild(inner);

    const colorCache = {};
    const getColor = (ch) => {
        if (!colorCache[ch]) colorCache[ch] = hexToRgb(NT_COLORS[ch.toUpperCase()] || "#9ca3af");
        return colorCache[ch];
    };

    // === RULER CANVAS ===
    const ruler = document.createElement("canvas");
    ruler.width = W;
    ruler.height = RULER_H;
    ruler.style.cssText = `display:block;width:${W}px;height:${RULER_H}px;position:sticky;top:0;z-index:2;background:#f3f4f6;`;
    inner.appendChild(ruler);

    const rctx = ruler.getContext("2d");
    rctx.fillStyle = "#f3f4f6";
    rctx.fillRect(0, 0, W, RULER_H);
    rctx.strokeStyle = "#9ca3af";
    rctx.fillStyle = "#6b7280";
    rctx.font = "9px monospace";
    rctx.textAlign = "left";

    // Aim for ~8 ticks visible
    const rawStep = seqLen / 8;
    const mag = Math.pow(10, Math.floor(Math.log10(rawStep)));
    const tick = Math.ceil(rawStep / mag) * mag;

    for (let pos = 0; pos <= seqLen; pos += tick) {
        const px = Math.floor(pos * W / seqLen);
        rctx.strokeStyle = "#9ca3af";
        rctx.lineWidth = 1;
        rctx.beginPath();
        rctx.moveTo(px + 0.5, RULER_H - 10);
        rctx.lineTo(px + 0.5, RULER_H - 2);
        rctx.stroke();
        if (pos > 0) {
            rctx.fillStyle = "#6b7280";
            rctx.fillText(pos.toLocaleString(), px + 2, RULER_H - 12);
        }
    }

    // === MAIN ALIGNMENT CANVAS ===
    const mainH = nSeq * cellH;
    const main = document.createElement("canvas");
    main.width = W;
    main.height = mainH;
    main.style.cssText = `display:block;width:${W}px;height:${mainH}px;image-rendering:pixelated;`;
    inner.appendChild(main);

    const mctx = main.getContext("2d");
    const mimg = mctx.createImageData(W, mainH);
    const mdata = mimg.data;

    for (let si = 0; si < nSeq; si++) {
        const seq = sequences[si].seq.toUpperCase();
        const rowY = si * cellH;
        for (let px = 0; px < W; px++) {
            const pi = Math.floor(px * seqLen / W);
            const c = getColor(seq[pi] || "");
            for (let dy = 0; dy < cellH; dy++) {
                const idx = ((rowY + dy) * W + px) * 4;
                mdata[idx] = c.r;
                mdata[idx + 1] = c.g;
                mdata[idx + 2] = c.b;
                mdata[idx + 3] = 255;
            }
        }
    }
    mctx.putImageData(mimg, 0, 0);

    // === CONSERVATION CANVAS ===
    const cons = document.createElement("canvas");
    cons.width = W;
    cons.height = CONS_H;
    cons.style.cssText = `display:block;width:${W}px;height:${CONS_H}px;background:#f9fafb;border-top:1px solid #e5e7eb;`;
    inner.appendChild(cons);

    const cctx = cons.getContext("2d");
    cctx.fillStyle = "#f9fafb";
    cctx.fillRect(0, 0, W, CONS_H);

    // Label
    cctx.font = "8px sans-serif";
    cctx.fillStyle = "#9ca3af";
    cctx.fillText("Conservation", 2, 9);

    const BAR_AREA = CONS_H - 14;

    for (let px = 0; px < W; px++) {
        const pi = Math.floor(px * seqLen / W);
        // Count bases at this position
        const counts = {};
        let total = 0;
        for (let si = 0; si < nSeq; si++) {
            const ch = sequences[si].seq[pi]?.toUpperCase();
            if (ch && ch !== "-" && ch !== "?") {
                counts[ch] = (counts[ch] || 0) + 1;
                total++;
            }
        }
        if (total === 0) continue;
        const maxCount = Math.max(...Object.values(counts));
        const score = maxCount / nSeq; // 0..1

        let color;
        if (score >= 0.8) color = "#16a34a"; // green
        else if (score >= 0.5) color = "#d97706"; // amber
        else color = "#dc2626"; // red

        const barH = Math.round(score * BAR_AREA);
        cctx.fillStyle = color;
        cctx.fillRect(px, CONS_H - 2 - barH, 1, barH);
    }

    // baseline
    cctx.strokeStyle = "#d1d5db";
    cctx.lineWidth = 1;
    cctx.beginPath();
    cctx.moveTo(0, CONS_H - 2);
    cctx.lineTo(W, CONS_H - 2);
    cctx.stroke();

    // === INTERACTIVE TOOLTIP on main canvas ===
    const tooltip = document.createElement("div");
    const tooltipTheme = document.documentElement.getAttribute("data-theme") === "dark"
        ? {
            background: "rgba(248, 250, 252, 0.98)",
            color: "#0f172a",
            border: "1px solid rgba(99, 102, 241, 0.24)",
            boxShadow: "0 8px 18px rgba(2, 6, 23, 0.28)",
        }
        : {
            background: "rgba(50, 55, 149, 0.96)",
            color: "#ffffff",
            border: "1px solid rgba(255, 255, 255, 0.12)",
            boxShadow: "0 6px 14px rgba(50, 55, 149, 0.24)",
        };
    tooltip.style.cssText = "position:absolute;pointer-events:none;font-size:11px;padding:3px 8px;border-radius:5px;white-space:nowrap;display:none;z-index:10;max-width:320px;overflow:hidden;text-overflow:ellipsis;";
    tooltip.style.background = tooltipTheme.background;
    tooltip.style.color = tooltipTheme.color;
    tooltip.style.border = tooltipTheme.border;
    tooltip.style.boxShadow = tooltipTheme.boxShadow;
    inner.style.position = "relative";
    inner.appendChild(tooltip);

    main.style.cursor = "crosshair";
    main.addEventListener("mousemove", (e) => {
        const rect = main.getBoundingClientRect();
        const canvasY = e.clientY - rect.top;
        // Adjust for devicePixelRatio if needed (canvas CSS vs pixel)
        const scaleY = main.height / rect.height;
        const row = Math.floor((canvasY * scaleY) / cellH);
        if (row >= 0 && row < nSeq) {
            tooltip.style.display = "block";
            tooltip.textContent = sequences[row].name;
            // Position relative to inner div
            const innerRect = inner.getBoundingClientRect();
            const tx = Math.min(e.clientX - innerRect.left + 10, W - 5);
            const ty = e.clientY - innerRect.top - 24;
            tooltip.style.left = tx + "px";
            tooltip.style.top = ty + "px";
        }
    });
    main.addEventListener("mouseleave", () => { tooltip.style.display = "none"; });
}

function renderAlignmentResults(markerResults, hasTrimal) {
    state.alignmentResultsData = { markerResults, hasTrimal };
    const accordion = document.getElementById("alignmentResultsAccordion");
    if (!accordion) return;
    accordion.innerHTML = "";
    document.getElementById("alignmentResults").classList.remove("hidden");

    for (const r of markerResults) {
        const id = "marker-result-" + r.marker.replace(/[^a-zA-Z0-9]/g, "_");
        const trimStat = r.trimmedStats
            ? `<div class="text-center"><div class="text-xs text-gray-400 mb-0.5">${i18nText("step5.results.stats.trimmed", "Trimmed")}</div><div class="text-xs font-semibold text-gray-700">${r.trimmedStats.numSeqs} seq · ${r.trimmedStats.length} bp</div><div class="text-xs text-gray-400">${r.trimmedStats.gapPct}% ${i18nText("step5.results.stats.gaps", "gaps")}</div></div>`
            : `<div class="text-center text-xs text-gray-300 italic">${i18nText("step5.results.stats.noTrimming", "No trimming")}</div>`;

        const tabs = [
            { key: "raw", label: i18nText("step5.results.tabs.raw", "Raw"), icon: "fa-dna" },
            { key: "aligned", label: i18nText("step5.results.tabs.aligned", "Aligned"), icon: "fa-align-left" },
            ...(r.trimmedStats ? [{ key: "trimmed", label: i18nText("step5.results.tabs.trimmed", "Trimmed"), icon: "fa-scissors" }] : []),
        ];

        const tabBtns = tabs.map((t, i) =>
            `<button class="result-tab-btn px-3 py-1 rounded-lg text-xs font-medium transition-colors ${i === 0 ? 'bg-splace-blue-600 text-white' : 'text-gray-500 hover:bg-gray-100'}"
                data-marker="${r.marker}" data-tab="${t.key}">
                <i class="fa-solid ${t.icon} mr-1"></i>${t.label}
             </button>`
        ).join("");

        const panels = tabs.map((t, i) => {
            return `<div id="${id}-${t.key}" class="result-tab-panel ${i !== 0 ? 'hidden' : ''}">
                <div class="alignment-view-shell">
                    <div class="alignment-view-toolbar">
                        <div class="alignment-view-toggle">
                            <button class="alignment-view-btn active" type="button" data-marker="${r.marker}" data-tab="${t.key}" data-view="nt">${i18nText("step5.results.view.nt", "DNA view")}</button>
                            <button class="alignment-view-btn" type="button" data-marker="${r.marker}" data-tab="${t.key}" data-view="aa">${i18nText("step5.results.view.aa", "Translate to amino acids")}</button>
                        </div>
                        <div class="alignment-toolbar-meta">
                            <span class="alignment-view-hint">${i18nText("step5.results.view.hint", "Use translated mode to spot possible stop codons in the alignment.")}</span>
                            <div class="alignment-genetic-code">
                                <label for="genetic-code-${id}-${t.key}">${i18nText("step5.results.view.geneticCode", "Genetic code")}</label>
                                <select id="genetic-code-${id}-${t.key}" class="alignment-genetic-code-select" data-marker="${r.marker}" data-tab="${t.key}">
                                    <option value="1">${i18nText("step5.results.view.code.standard", "Standard (NCBI 1)")}</option>
                                    <option value="2">${i18nText("step5.results.view.code.vertMito", "Vertebrate mitochondrial (NCBI 2)")}</option>
                                    <option value="5">${i18nText("step5.results.view.code.invMito", "Invertebrate mitochondrial (NCBI 5)")}</option>
                                    <option value="11">${i18nText("step5.results.view.code.bacterial", "Bacterial, archaeal and plant plastid (NCBI 11)")}</option>
                                </select>
                            </div>
                        </div>
                    </div>
                    <div id="vizwrap-${id}-${t.key}"></div>
                    <div id="preview-${id}-${t.key}"></div>
                </div>
            </div>`;
        }).join("");

        const panel = document.createElement("div");
        panel.className = "bg-white rounded-xl border border-gray-200 overflow-hidden";
        panel.innerHTML = `
            <button class="w-full px-5 py-3 flex items-center justify-between text-left hover:bg-gray-50 transition-colors" onclick="toggleMarkerResult('${id}')">
                <div class="flex items-center gap-3">
                    <span class="step-badge badge-done" style="font-size:0.55rem;width:1rem;height:1rem;">✓</span>
                    <span class="font-semibold text-sm text-gray-800">${r.marker}</span>
                    <span class="text-xs text-gray-400">${r.alignedStats?.numSeqs ?? "?"} seq · ${r.alignedStats?.length ?? "?"} bp ${i18nText("step5.results.stats.alignedSuffix", "aligned")}</span>
                </div>
                <i id="${id}-chevron" class="fa-solid fa-chevron-down text-gray-400 text-xs transition-transform"></i>
            </button>
            <div id="${id}-body" class="hidden px-5 pb-5">
                <div class="grid grid-cols-3 gap-3 mb-4 pt-2">
                    <div class="text-center"><div class="text-xs text-gray-400 mb-0.5">${i18nText("step5.results.stats.raw", "Raw")}</div><div class="text-xs font-semibold text-gray-700">${r.rawStats?.numSeqs ?? "?"} seq · ${r.rawStats?.avgLen ?? "?"} bp ${i18nText("step5.results.stats.avg", "avg")}</div></div>
                    <div class="text-center"><div class="text-xs text-gray-400 mb-0.5">${i18nText("step5.results.stats.aligned", "Aligned")}</div><div class="text-xs font-semibold text-gray-700">${r.alignedStats?.numSeqs ?? "?"} seq · ${r.alignedStats?.length ?? "?"} bp</div><div class="text-xs text-gray-400">${r.alignedStats?.gapPct ?? "?"}% ${i18nText("step5.results.stats.gaps", "gaps")}</div></div>
                    ${trimStat}
                </div>
                <div class="flex gap-2 mb-3">${tabBtns}</div>
                <div class="bg-gray-50 rounded-lg p-2">${panels}</div>
            </div>`;
        accordion.appendChild(panel);

        // Store sequences on wrappers for lazy rendering
        // Render raw tab immediately (it's visible by default)
        requestAnimationFrame(() => {
            if (r.rawContent) renderAlignmentTab({ id, rawContent: r.rawContent, alignedContent: r.alignedContent, trimmedContent: r.trimmedContent }, "raw", "nt");
        });

        // Store content references for lazy rendering on tab switch / accordion open
        panel._markerData = {
            id, rawContent: r.rawContent,
            alignedContent: r.alignedContent,
            trimmedContent: r.trimmedContent,
        };
    }

    // Tab switching — lazy render enhanced tabs when first shown
    accordion.querySelectorAll(".result-tab-btn").forEach(btn => {
        btn.addEventListener("click", () => {
            const marker = btn.dataset.marker;
            const tab = btn.dataset.tab;
            const id2 = "marker-result-" + marker.replace(/[^a-zA-Z0-9]/g, "_");
            btn.closest(".flex").querySelectorAll(".result-tab-btn").forEach(b => {
                b.className = "result-tab-btn px-3 py-1 rounded-lg text-xs font-medium transition-colors text-gray-500 hover:bg-gray-100";
            });
            btn.className = "result-tab-btn px-3 py-1 rounded-lg text-xs font-medium transition-colors bg-splace-blue-600 text-white";
            btn.closest(".flex").nextElementSibling.querySelectorAll(".result-tab-panel").forEach(p => p.classList.add("hidden"));
            const panelEl = document.getElementById(`${id2}-${tab}`);
            panelEl?.classList.remove("hidden");

            // Lazy render on first tab switch
            const panelDiv = btn.closest(".bg-white");
            const data = panelDiv?._markerData;
            if (!data) return;
            const contentMap = { raw: data.rawContent, aligned: data.alignedContent, trimmed: data.trimmedContent };
            const content = contentMap[tab];
            if (content) {
                const wrap = document.getElementById(`vizwrap-${data.id}-${tab}`);
                if (wrap && !wrap.dataset.rendered) {
                    requestAnimationFrame(() => renderAlignmentTab(data, tab, wrap.dataset.view || "nt"));
                }
            }
        });
    });

    accordion.querySelectorAll(".alignment-view-btn").forEach((btn) => {
        btn.addEventListener("click", () => {
            const marker = btn.dataset.marker;
            const tab = btn.dataset.tab;
            const mode = btn.dataset.view;
            const id2 = "marker-result-" + marker.replace(/[^a-zA-Z0-9]/g, "_");
            const panelEl = document.getElementById(`${id2}-${tab}`);
            panelEl?.querySelectorAll(".alignment-view-btn").forEach((candidate) => {
                candidate.classList.toggle("active", candidate === btn);
            });
            const panelDiv = btn.closest(".bg-white");
            const data = panelDiv?._markerData;
            if (!data) return;
            renderAlignmentTab(data, tab, mode);
        });
    });

    accordion.querySelectorAll(".alignment-genetic-code-select").forEach((select) => {
        select.addEventListener("change", () => {
            const marker = select.dataset.marker;
            const tab = select.dataset.tab;
            const id2 = "marker-result-" + marker.replace(/[^a-zA-Z0-9]/g, "_");
            const panelEl = document.getElementById(`${id2}-${tab}`);
            const activeBtn = panelEl?.querySelector(".alignment-view-btn.active");
            const mode = activeBtn?.dataset.view || "nt";
            const panelDiv = select.closest(".bg-white");
            const data = panelDiv?._markerData;
            if (!data) return;
            renderAlignmentTab(data, tab, mode);
        });
    });
}

window.refreshAlignmentUiTranslations = function () {
    if (state.alignmentResultsData) {
        renderAlignmentResults(state.alignmentResultsData.markerResults, state.alignmentResultsData.hasTrimal);
    }
};

window.toggleMarkerResult = function (id) {
    const body = document.getElementById(id + "-body");
    const chevron = document.getElementById(id + "-chevron");
    const open = body.classList.toggle("hidden");
    chevron.style.transform = open ? "" : "rotate(180deg)";
};

// ========================================================================
// NEXUS Concatenation
// ========================================================================

// State for concatenation
let _lastMarkerResultsUI = []; // mirror of backend results, kept in JS

function buildConcat(markerResults, useTrimmed, allowMissing) {
    // Collect all species and their sequences per gene
    const geneOrder = markerResults.map(r => r.marker);
    const geneSeqs = {};   // gene → {speciesName → sequence}
    const geneLens = {};   // gene → alignment length

    for (const r of markerResults) {
        const content = useTrimmed && r.trimmedContent ? r.trimmedContent : r.alignedContent;
        if (!content) continue;
        const seqs = parseFasta(content);
        if (!seqs.length) continue;
        geneLens[r.marker] = seqs[0].seq.length;
        geneSeqs[r.marker] = {};
        for (const s of seqs) {
            geneSeqs[r.marker][s.name] = s.seq;
        }
    }

    // Union of all species
    const allSpecies = [...new Set(
        Object.values(geneSeqs).flatMap(m => Object.keys(m))
    )].sort();

    // Check for missing
    const missingReport = []; // [{species, genes[]}]
    for (const sp of allSpecies) {
        const missing = geneOrder.filter(g => geneSeqs[g] && !geneSeqs[g][sp]);
        if (missing.length) missingReport.push({ species: sp, genes: missing });
    }

    // Determine final species list
    let finalSpecies;
    if (allowMissing) {
        finalSpecies = allSpecies; // keep all, fill missing with ?
    } else {
        const removedSet = new Set(missingReport.map(m => m.species));
        finalSpecies = allSpecies.filter(sp => !removedSet.has(sp));
    }

    // Concatenate
    const alignedSeqs = [];
    const partitions = [];
    let pos = 1;

    for (const sp of finalSpecies) {
        let concat = "";
        for (const gene of geneOrder) {
            if (!geneSeqs[gene]) continue;
            const len = geneLens[gene];
            concat += (geneSeqs[gene][sp] || "?".repeat(len));
        }
        alignedSeqs.push({ name: sp, seq: concat });
    }

    for (const gene of geneOrder) {
        if (!geneSeqs[gene]) continue;
        const len = geneLens[gene];
        partitions.push({ gene, start: pos, end: pos + len - 1, len });
        pos += len;
    }

    return { alignedSeqs, partitions, totalLen: pos - 1, species: finalSpecies, missingReport, removedCount: allowMissing ? 0 : missingReport.length };
}

function buildNexusString(alignedSeqs, partitions, outgroups) {
    const ntax = alignedSeqs.length;
    const nchar = alignedSeqs[0]?.seq.length || 0;
    // Pad names for alignment
    const maxNameLen = Math.max(...alignedSeqs.map(s => s.name.length));
    const matrix = alignedSeqs.map(s =>
        `    ${s.name.padEnd(maxNameLen + 2)}${s.seq}`
    ).join("\n");

    const charsets = partitions.map(p =>
        `    charset ${p.gene} = ${p.start}-${p.end};`
    ).join("\n");

    const outgroupBlock = outgroups && outgroups.length
        ? `\nbegin assumptions;\n    outgroup ${outgroups.join(" ")};\nend;\n`
        : "";

    return `#NEXUS

begin data;
    dimensions ntax=${ntax} nchar=${nchar};
    format datatype=dna missing=? gap=-;
    matrix
${matrix}
    ;
end;

begin sets;
${charsets}
end;
${outgroupBlock}`;
}

function buildPartitionString(partitions) {
    return partitions.map(p => `DNA, ${p.gene} = ${p.start}-${p.end}`).join("\n");
}

function renderConcatSection(markerResults, hasTrimal) {
    _lastMarkerResultsUI = markerResults;
    const sec = document.getElementById("concatSection");
    if (!sec) return;
    sec.classList.remove("hidden");
    sec.scrollIntoView({ behavior: "smooth", block: "start" });

    // Show trimmed option only if trimAl was run
    const trimmedOption = document.getElementById("concatTrimmedOption");
    if (trimmedOption) trimmedOption.classList.toggle("hidden", !hasTrimal);
    // If no trimAl, force aligned
    if (!hasTrimal) {
        const noRadio = document.getElementById("concatUseTrimmedNo");
        if (noRadio) noRadio.checked = true;
    }

    // Populate outgroup list from all species across markers
    const speciesSet = new Set();
    for (const r of markerResults) {
        const content = r.trimmedContent || r.alignedContent;
        if (content) parseFasta(content).forEach(s => speciesSet.add(s.name));
    }
    const ogList = document.getElementById("outgroupList");
    if (ogList) {
        ogList.innerHTML = [...speciesSet].sort().map(sp =>
            `<label class="flex items-center gap-2 text-xs text-gray-700 cursor-pointer hover:bg-gray-50 px-2 py-1 rounded">
                <input type="checkbox" class="outgroup-cb rounded text-splace-blue-600" value="${sp.replace(/"/g, '&quot;')}">
                <span class="truncate">${sp}</span>
            </label>`
        ).join("");
    }

    // Filter outgroup search — use assignment to avoid stacking listeners on re-render
    const ogSearch = document.getElementById("outgroupSearch");
    if (ogSearch) {
        ogSearch.value = "";        // reset on re-render
        ogSearch.oninput = filterOutgroups;
    }
}

function filterOutgroups() {
    const input = document.getElementById("outgroupSearch");
    const q = (input?.value || "").toLowerCase().trim();
    const ogList = document.getElementById("outgroupList");
    if (!ogList) return;
    ogList.querySelectorAll("label").forEach(label => {
        // Restore original text before re-highlighting
        const span = label.querySelector("span.truncate");
        if (!span) return;
        if (span.dataset.orig === undefined) span.dataset.orig = span.textContent;
        const orig = span.dataset.orig;
        const matches = q !== "" && orig.toLowerCase().includes(q);
        label.classList.toggle("hidden", q !== "" && !matches);
        if (matches) {
            // Highlight the match
            const idx = orig.toLowerCase().indexOf(q);
            span.innerHTML =
                escHtml(orig.slice(0, idx)) +
                `<mark class="bg-yellow-200 text-yellow-900 rounded-sm px-0.5">${escHtml(orig.slice(idx, idx + q.length))}</mark>` +
                escHtml(orig.slice(idx + q.length));
        } else {
            span.innerHTML = escHtml(orig);
        }
    });
}

function escHtml(s) {
    return s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
}

function buildIqtreeCommand() {
    const model = document.getElementById("iqtreeModel")?.value || "TEST";
    const customModel = document.getElementById("iqtreeModelCustom")?.value || "GTR+G";
    const bootstrap = document.getElementById("iqtreeBootstrap")?.value || "-B 1000";
    const threads = document.getElementById("iqtreeThreads")?.value || "4";
    const extra = (document.getElementById("iqtreeExtra")?.value || "").trim();
    const outgroups = window._nexusData?.outgroups || [];
    const dir = "~/Documents/SPLACE/iqtree";

    const parts = [
        "iqtree3",
        "-s", `${dir}/concat.nex`,
        "-p", `${dir}/partition.txt`,
        "--prefix", `${dir}/concat`,
        "-T", threads,
        "-m", model === "custom" ? customModel : model,
        ...(bootstrap ? bootstrap.split(" ").filter(Boolean) : []),
        ...(outgroups.length ? ["-o", outgroups.join(",")] : []),
        ...(extra ? extra.split(/\s+/).filter(Boolean) : []),
        "--redo",
    ];
    return parts.join(" ");
}

window.updateIqtreeCommand = function () {
    const el = document.getElementById("iqtreeCommandPreview");
    if (el) el.textContent = buildIqtreeCommand();
};

window.generateNexus = function () {
    const useTrimmed = document.getElementById("concatUseTrimmedYes")?.checked ?? true;
    const allowMissing = document.getElementById("concatAllowMissing")?.checked ?? false;
    const outgroups = [...document.querySelectorAll(".outgroup-cb:checked")].map(c => c.value);

    const result = buildConcat(_lastMarkerResultsUI, useTrimmed, allowMissing);

    const warn = document.getElementById("concatMissingWarning");
    warn.classList.add("hidden");
    warn.className = "text-xs rounded-lg px-3 py-2 border";

    if (result.missingReport.length) {
        if (allowMissing) {
            warn.className += " text-amber-700 bg-amber-50 border-amber-200";
            warn.textContent = `${result.missingReport.length} species have missing genes — filled with "?" in the matrix.`;
        } else {
            warn.className += " text-amber-700 bg-amber-50 border-amber-200";
            const names = result.missingReport.map(m => m.species).join(", ");
            warn.textContent = `${result.removedCount} species removed from matrix due to missing genes: ${names}. Enable "Allow missing sequences" to keep them with "?" characters.`;
        }
        warn.classList.remove("hidden");
    }

    const nexus = buildNexusString(result.alignedSeqs, result.partitions, outgroups);
    const partition = buildPartitionString(result.partitions);

    // Store for IQ-TREE section
    window._nexusData = { nexus, partition, result, outgroups };

    // Show partition table
    const tbl = document.getElementById("partitionTable");
    if (tbl) {
        tbl.innerHTML = `
            <table class="w-full text-xs border-collapse">
                <thead><tr class="bg-gray-100">
                    <th class="text-left px-3 py-2 font-semibold text-gray-600">Gene</th>
                    <th class="text-right px-3 py-2 font-semibold text-gray-600">Start</th>
                    <th class="text-right px-3 py-2 font-semibold text-gray-600">End</th>
                    <th class="text-right px-3 py-2 font-semibold text-gray-600">Length (bp)</th>
                </tr></thead>
                <tbody>
                    ${result.partitions.map((p, i) => `
                        <tr class="${i % 2 === 0 ? '' : 'bg-gray-50'}">
                            <td class="px-3 py-1.5 font-medium text-splace-blue-700">${p.gene}</td>
                            <td class="px-3 py-1.5 text-right text-gray-600">${p.start}</td>
                            <td class="px-3 py-1.5 text-right text-gray-600">${p.end}</td>
                            <td class="px-3 py-1.5 text-right text-gray-600">${p.len}</td>
                        </tr>`).join("")}
                    <tr class="border-t border-gray-200 font-semibold">
                        <td class="px-3 py-1.5 text-gray-700">Total</td>
                        <td class="px-3 py-1.5 text-right text-gray-500">—</td>
                        <td class="px-3 py-1.5 text-right text-gray-500">—</td>
                        <td class="px-3 py-1.5 text-right text-gray-700">${result.totalLen}</td>
                    </tr>
                </tbody>
            </table>`;
        tbl.classList.remove("hidden");
    }

    // Show summary
    const summary = document.getElementById("concatSummary");
    if (summary) {
        summary.textContent = `${result.alignedSeqs.length} species · ${result.partitions.length} genes · ${result.totalLen} bp`;
        summary.classList.remove("hidden");
    }

    // Update step badge for concat
    setBadge("step-concat-badge", "done");

    // Show IQ-TREE section
    renderIqtreeSection(outgroups);
};

window.downloadNexus = function () {
    if (!window._nexusData) return;
    const blob = new Blob([window._nexusData.nexus], { type: "text/plain" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = "concatenated.nex";
    a.click();
};

window.downloadPartition = function () {
    if (!window._nexusData) return;
    const blob = new Blob([window._nexusData.partition], { type: "text/plain" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = "partitions.txt";
    a.click();
};

// ========================================================================
// IQ-TREE3
// ========================================================================

function renderIqtreeSection(outgroups) {
    const sec = document.getElementById("iqtreeSection");
    if (!sec) return;
    sec.classList.remove("hidden");
    sec.scrollIntoView({ behavior: "smooth", block: "start" });

    // Show selected outgroups
    const ogLabel = document.getElementById("iqtreeOutgroupLabel");
    if (ogLabel) {
        ogLabel.textContent = outgroups && outgroups.length
            ? outgroups.join(", ")
            : "None selected";
    }
    // Store outgroups for command preview
    if (window._nexusData) window._nexusData.outgroups = outgroups || [];

    // Pre-fill thread count, then update command preview
    if (window.electronAPI?.getCpuCount) {
        window.electronAPI.getCpuCount().then(n => {
            const t = document.getElementById("iqtreeThreads");
            if (t && !t.value) t.value = Math.max(1, n - 2);
            updateIqtreeCommand();
        });
    } else {
        updateIqtreeCommand();
    }
}

window.runIqtree = async function () {
    if (!window._nexusData || !window.electronAPI) return;

    const model = document.getElementById("iqtreeModel")?.value || "TEST";
    const bootstrap = document.getElementById("iqtreeBootstrap")?.value || "-B 1000";
    const threads = parseInt(document.getElementById("iqtreeThreads")?.value) || 4;
    const extra = document.getElementById("iqtreeExtra")?.value || "";
    const perGene = document.getElementById("iqtreePerGene")?.checked ?? false;
    const outgroup = window._nexusData.outgroups;

    const params = [
        ...(model === "custom"
            ? ["-m", document.getElementById("iqtreeModelCustom")?.value || "GTR+G"]
            : ["-m", model]),
        ...bootstrap.split(" ").filter(Boolean),
        ...(extra ? extra.split(/\s+/).filter(Boolean) : []),
    ];

    if (!window.electronAPI.runIqtree) {
        alert("IQ-TREE IPC not available. Please restart the app.");
        return;
    }

    // Build per-gene FASTA map (trimmed preferred)
    let geneFiles = null;
    if (perGene) {
        geneFiles = {};
        const useTrimmed = document.getElementById("concatUseTrimmedYes")?.checked ?? true;
        for (const r of _lastMarkerResultsUI) {
            const content = (useTrimmed && r.trimmedContent) ? r.trimmedContent : r.alignedContent;
            if (content) geneFiles[r.marker] = content;
        }
    }

    openIqtreeModal();

    window.electronAPI.runIqtree({
        nexus: window._nexusData.nexus,
        partition: window._nexusData.partition,
        params,
        threads,
        outgroup,
        perGene,
        geneFiles,
    });
};

function openIqtreeModal() {
    const modal = document.getElementById("iqtreeModal");
    if (!modal) return;
    modal.classList.remove("hidden");
    document.getElementById("iqtreeModalLog").textContent = "";
    document.getElementById("iqtreeModalTitle").textContent = "Running IQ-TREE…";
    document.getElementById("iqtreeModalClose").classList.add("hidden");
    document.getElementById("iqtreeModalIcon").className = "fa-solid fa-spinner fa-spin text-splace-blue-600 text-lg";
}

if (window.electronAPI?.onIqtreeProgress) {
    window.electronAPI.onIqtreeProgress((data) => {
        const log = document.getElementById("iqtreeModalLog");
        if (log && data.message) {
            log.textContent += data.message + "\n";
            // scroll the overflow container (parent of <pre>), not the <pre> itself
            const logContainer = log.parentElement;
            if (logContainer) logContainer.scrollTop = logContainer.scrollHeight;
        }
    });
}

if (window.electronAPI?.onIqtreeDone) {
    window.electronAPI.onIqtreeDone((result) => {
        const icon = document.getElementById("iqtreeModalIcon");
        const title = document.getElementById("iqtreeModalTitle");
        const close = document.getElementById("iqtreeModalClose");

        if (result.success) {
            icon.className = "fa-solid fa-check text-green-600 text-lg";
            title.textContent = "IQ-TREE analysis complete";
        } else {
            icon.className = "fa-solid fa-xmark text-red-600 text-lg";
            title.textContent = "IQ-TREE failed";
        }
        close.classList.remove("hidden");
        close.onclick = () => {
            document.getElementById("iqtreeModal").classList.add("hidden");
            if (result.success) {
                setBadge("step-iqtree-badge", "done");
                showIqtreeResults(result);
            }
        };
    });
}

function showIqtreeResults(result) {
    const res = document.getElementById("iqtreeResults");
    if (!res) return;
    res.classList.remove("hidden");

    const fileList = document.getElementById("iqtreeFileList");
    if (fileList && result.files) {
        fileList.innerHTML = result.files.map(f =>
            `<div class="flex items-center gap-2 text-xs text-gray-700 py-1 border-b border-gray-100 last:border-0">
                <i class="fa-solid fa-file-code text-splace-blue-400 w-4"></i>
                <span class="font-mono truncate">${f}</span>
            </div>`
        ).join("");
    }

    // Per-gene summary
    if (result.geneResults && result.geneResults.length) {
        const geneDiv = document.getElementById("iqtreeGeneResults");
        if (geneDiv) {
            geneDiv.classList.remove("hidden");
            geneDiv.innerHTML = `
                <p class="text-xs font-semibold text-gray-700 mb-2">Per-gene trees</p>
                <div class="grid grid-cols-3 gap-1.5">
                    ${result.geneResults.map(r =>
                `<div class="flex items-center gap-1.5 text-xs px-2 py-1 rounded-lg ${r.success ? 'bg-green-50 text-green-700' : 'bg-red-50 text-red-700'}">
                            <i class="fa-solid ${r.success ? 'fa-check' : 'fa-xmark'} w-3"></i>
                            <span>${r.gene}</span>
                        </div>`
            ).join("")}
                </div>`;
        }
    }

    document.getElementById("iqtreeOutputDir").textContent = result.outputDir || "";
    window._iqtreeOutputDir = result.outputDir;
}

window.saveIqtreeZip = async function () {
    if (!window.electronAPI?.saveZip) return;
    if (!window._alignmentOutputDir && !window._iqtreeOutputDir) return;

    // Collect raw FASTAs from stored marker results
    const rawFastas = {};
    for (const r of _lastMarkerResultsUI) {
        if (r.rawContent) rawFastas[r.marker] = r.rawContent;
    }

    const result = await window.electronAPI.saveZip({
        alignmentDir: window._alignmentOutputDir || null,
        iqtreeDir: window._iqtreeOutputDir || null,
        rawFastas,
        suggestedName: "SPLACE_results.zip",
    });
    if (result.cancelled) return;
    if (result.success) {
        alert(`Saved: ${result.filePath}`);
    } else {
        alert(`Error creating ZIP: ${result.error}`);
    }
};

