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
    geneNameOverrides: new Map(),
    unknownGenesChecked: false,
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
function isKnownGeneName(name) {
    if (!name) return false;
    return (name in MitochondrialGenes) || (name in ChloroplastGenes);
}

function standardizeGeneName(rawName, product, dataType) {
    // Check user-supplied overrides first
    const overrideKey = rawName || product;
    if (overrideKey && state.geneNameOverrides.has(overrideKey)) {
        return state.geneNameOverrides.get(overrideKey);
    }

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
// Unknown Gene Names — collect & modal
// ========================================================================
function collectUnknownGenes() {
    const dataType = state.detectedDataType;
    // Map from rawKey → { rawKey, featureType, occurrences: [{species, accession}] }
    const unknownMap = new Map();

    for (const record of state.records) {
        for (const feat of record.features) {
            if (feat.type !== "CDS" && feat.type !== "rRNA") continue;
            const rawGene = feat.qualifiers.gene || null;
            const rawProduct = feat.qualifiers.product || null;
            if (!rawGene && !rawProduct) continue;

            const standardized = standardizeGeneName(rawGene, rawProduct, dataType);
            if (!standardized || isKnownGeneName(standardized)) continue;

            // standardizeGeneName fell back to the raw name — not in the dictionary
            const rawKey = standardized;
            if (!unknownMap.has(rawKey)) {
                unknownMap.set(rawKey, { rawKey, featureType: feat.type, occurrences: [] });
            }
            unknownMap.get(rawKey).occurrences.push({
                species: record.organism || "—",
                accession: record.accession || "—",
                source: record.source || "ncbi",
                fileName: record.fileName || null,
            });
        }
    }
    return Array.from(unknownMap.values());
}

function buildGeneOptions(dataType) {
    const primaryDict = dataType === "cp" ? ChloroplastGenes : MitochondrialGenes;
    const secondaryDict = dataType === "cp" ? MitochondrialGenes : ChloroplastGenes;
    const primaryLabel = dataType === "cp" ? "Chloroplast" : "Mitochondrial";
    const secondaryLabel = dataType === "cp" ? "Mitochondrial" : "Chloroplast";
    const skipLabel = i18nText("modal.unknowngenes.select.skip", "(skip)");
    const opts = [{ value: "", label: skipLabel, group: null }];
    for (const k of Object.keys(primaryDict).sort()) opts.push({ value: k, label: k, group: primaryLabel });
    for (const k of Object.keys(secondaryDict).sort()) opts.push({ value: k, label: k, group: secondaryLabel });
    return opts;
}

function buildSearchableSelect(geneOptions) {
    const skipLabel = geneOptions.find(o => o.value === "")?.label || "(skip)";
    const searchPlaceholder = i18nText("modal.unknowngenes.search", "Search gene…");
    const optItems = geneOptions.map(o => {
        if (!o.value) {
            return `<div class="gene-opt px-3 py-1.5 cursor-pointer hover:bg-gray-50 text-sm text-gray-400 italic border-b border-gray-100" data-value="" data-label="${escHtml(o.label)}">${escHtml(o.label)}</div>`;
        }
        return `<div class="gene-opt px-3 py-1.5 cursor-pointer hover:bg-splace-blue-50 text-sm text-gray-800 font-normal" data-value="${escHtml(o.value)}" data-label="${escHtml(o.label)}" data-group="${escHtml(o.group || "")}">${escHtml(o.label)}</div>`;
    }).join("");
    return `<div class="gene-combobox relative">
        <input type="hidden" class="gene-value" value="">
        <button type="button" class="gene-trigger w-full flex items-center justify-between gap-2 border border-gray-200 rounded-lg px-3 py-2 bg-white hover:border-splace-blue-400 transition-colors focus:outline-none focus:ring-2 focus:ring-splace-blue-200 text-left">
            <span class="gene-display text-sm text-gray-400 truncate flex-1">${escHtml(skipLabel)}</span>
            <i class="fa-solid fa-chevron-down text-xs text-gray-400 flex-shrink-0 transition-transform duration-150"></i>
        </button>
        <div class="gene-dropdown hidden absolute left-0 right-0 top-full mt-1 bg-white border border-gray-200 rounded-xl shadow-xl z-50 overflow-hidden" style="min-width:180px;">
            <div class="p-2 border-b border-gray-100">
                <input type="text" class="gene-search w-full text-sm px-2.5 py-1.5 rounded-lg border border-gray-200 focus:outline-none focus:border-splace-blue-400 bg-white placeholder-gray-400 transition-colors" placeholder="${escHtml(searchPlaceholder)}">
            </div>
            <div class="gene-opts-list overflow-y-auto" style="max-height:200px;">${optItems}</div>
        </div>
    </div>`;
}

function openUnknownGenesModal(entries) {
    const modal = document.getElementById("unknownGenesModal");
    if (!modal) return;

    const dataType = state.detectedDataType;
    const geneOptions = buildGeneOptions(dataType);

    const colSpecies = i18nText("modal.unknowngenes.col.species", "Species");
    const colSource = i18nText("modal.unknowngenes.col.source", "File / Accession");
    const colFound = i18nText("modal.unknowngenes.col.found", "Gene Not Found");
    const colCorrect = i18nText("modal.unknowngenes.col.correct", "Correct Name");
    const moreLabel = i18nText("modal.unknowngenes.more", "more");

    const tooltipData = [];
    const rows = entries.map((entry, idx) => {
        const first = entry.occurrences[0];
        const extra = entry.occurrences.length - 1;
        const speciesCell = `<span class="italic text-gray-700 text-sm">${escHtml(first.species)}</span>` +
            (extra > 0 ? ` <span class="text-sm text-gray-400 whitespace-nowrap">+${extra} ${moreLabel}</span>` : "");
        const isFile = first.source === "file" && first.fileName;
        const sourceCell = isFile
            ? `<span class="inline-flex items-center gap-1 text-sm text-gray-600"><i class="fa-regular fa-file-lines text-gray-400 flex-shrink-0"></i><span class="truncate">${escHtml(first.fileName)}</span></span>`
            : `<span class="font-mono text-sm text-gray-500">${escHtml(first.accession)}</span>`;
        const tipLines = entry.occurrences.map(o =>
            o.source === "file" && o.fileName
                ? `<em>${escHtml(o.species)}</em> &mdash; <code>${escHtml(o.fileName)}</code>`
                : `<em>${escHtml(o.species)}</em> &mdash; <code>${escHtml(o.accession)}</code>`
        ).join("<br>");
        tooltipData.push(tipLines);
        const stripeClass = idx % 2 === 1 ? "bg-gray-50/60" : "";
        return `<tr data-rawkey="${escHtml(entry.rawKey)}" class="border-t border-gray-100 hover:bg-amber-50/40 transition-colors ${stripeClass}">
            <td class="px-4 py-3">${speciesCell}</td>
            <td class="px-4 py-3 max-w-[200px] cursor-help" data-occ-idx="${idx}">${sourceCell}</td>
            <td class="px-4 py-3">
                <span class="inline-flex items-center gap-1 rounded-full bg-amber-100 text-amber-800 text-sm font-semibold px-2.5 py-1 leading-none">
                    <i class="fa-solid fa-circle-exclamation text-amber-500 text-xs"></i>
                    ${escHtml(entry.rawKey)}
                </span>
            </td>
            <td class="px-4 py-3" style="min-width:240px;">${buildSearchableSelect(geneOptions)}</td>
        </tr>`;
    }).join("");

    document.getElementById("unknownGenesTableHead").innerHTML =
        `<tr class="bg-amber-50 border-b border-amber-100">
            <th class="px-4 py-2.5 text-left text-xs font-semibold text-amber-800 uppercase tracking-wider">${escHtml(colSpecies)}</th>
            <th class="px-4 py-2.5 text-left text-xs font-semibold text-amber-800 uppercase tracking-wider">${escHtml(colSource)}</th>
            <th class="px-4 py-2.5 text-left text-xs font-semibold text-amber-800 uppercase tracking-wider">${escHtml(colFound)}</th>
            <th class="px-4 py-2.5 text-left text-xs font-semibold text-amber-800 uppercase tracking-wider" style="min-width:240px;">${escHtml(colCorrect)}</th>
        </tr>`;

    document.getElementById("unknownGenesTableBody").innerHTML = rows;

    if (typeof tippy !== "undefined") {
        document.querySelectorAll("#unknownGenesTableBody td[data-occ-idx]").forEach(td => {
            const occIdx = parseInt(td.dataset.occIdx);
            if (tooltipData[occIdx]) {
                tippy(td, { ...TIPPY_DEFAULTS, content: tooltipData[occIdx], placement: "right" });
            }
        });
    }

    const count = entries.length;
    document.getElementById("unknownGenesCount").textContent =
        i18nText("modal.unknowngenes.count", `${count} gene(s) not found in dictionary`, { count });

    modal.classList.remove("hidden");
    initGeneComboboxes();
}

window.closeUnknownGenesModal = function () {
    document.getElementById("unknownGenesModal").classList.add("hidden");
};

let _geneComboboxCloseHandler = null;

function initGeneComboboxes() {
    const tbody = document.getElementById("unknownGenesTableBody");
    if (!tbody) return;

    tbody.querySelectorAll(".gene-combobox").forEach(combobox => {
        const trigger = combobox.querySelector(".gene-trigger");
        const dropdown = combobox.querySelector(".gene-dropdown");
        const search = combobox.querySelector(".gene-search");
        const hidden = combobox.querySelector(".gene-value");
        const display = combobox.querySelector(".gene-display");
        const chevron = combobox.querySelector(".fa-chevron-down");

        function openDropdown() {
            document.querySelectorAll(".gene-dropdown:not(.hidden)").forEach(d => {
                if (d !== dropdown) {
                    d.classList.add("hidden");
                    d.style.cssText = "";
                    d.closest(".gene-combobox")?.querySelector(".fa-chevron-down")?.classList.remove("rotate-180");
                }
            });
            const rect = trigger.getBoundingClientRect();
            Object.assign(dropdown.style, {
                position: "fixed",
                top: (rect.bottom + 4) + "px",
                left: rect.left + "px",
                width: Math.max(rect.width, 240) + "px",
                zIndex: "9999",
            });
            dropdown.classList.remove("hidden");
            chevron.classList.add("rotate-180");
            search.value = "";
            combobox.querySelectorAll(".gene-opt").forEach(o => {
                o.style.display = "";
                o.textContent = o.dataset.label;
            });
            search.focus();
        }

        function closeDropdown() {
            dropdown.classList.add("hidden");
            dropdown.style.cssText = "";
            chevron.classList.remove("rotate-180");
        }

        trigger.addEventListener("click", e => {
            e.stopPropagation();
            dropdown.classList.contains("hidden") ? openDropdown() : closeDropdown();
        });

        search.addEventListener("input", () => {
            const query = search.value.trim().toLowerCase();
            combobox.querySelectorAll(".gene-opt").forEach((option) => {
                const label = option.dataset.label || option.textContent || "";
                const match = !query || label.toLowerCase().includes(query);
                option.style.display = match ? "" : "none";
                if (!match || !query) {
                    option.textContent = label;
                    return;
                }
                const idx = label.toLowerCase().indexOf(query);
                option.innerHTML =
                    escHtml(label.slice(0, idx)) +
                    `<mark class="bg-yellow-200 text-yellow-900 rounded-sm px-0.5">${escHtml(label.slice(idx, idx + query.length))}</mark>` +
                    escHtml(label.slice(idx + query.length));
            });
        });

        dropdown.querySelectorAll(".gene-opt").forEach(option => {
            option.addEventListener("click", () => {
                hidden.value = option.dataset.value || option.dataset.label || option.textContent.trim();
                display.textContent = option.dataset.label || option.textContent.trim();
                closeDropdown();
            });
        });
    });

    const scrollEl = document.querySelector("#unknownGenesModal .overflow-y-auto");
    if (scrollEl && !scrollEl._geneScrollHandler) {
        scrollEl._geneScrollHandler = () => {
            document.querySelectorAll(".gene-dropdown:not(.hidden)").forEach(d => {
                d.classList.add("hidden");
                d.style.cssText = "";
                d.closest(".gene-combobox")?.querySelector(".fa-chevron-down")?.classList.remove("rotate-180");
            });
        };
        scrollEl.addEventListener("scroll", scrollEl._geneScrollHandler, { passive: true });
    }
}

window.applyGeneNameOverrides = function () {
    const rows = document.querySelectorAll("#unknownGenesTableBody tr[data-rawkey]");
    rows.forEach(row => {
        const rawKey = row.dataset.rawkey;
        const valInput = row.querySelector(".gene-value");
        if (valInput && valInput.value) {
            state.geneNameOverrides.set(rawKey, valInput.value);
        }
    });
    window.closeUnknownGenesModal();
    renderGeneSelection();
};

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


const MARKER_KEYS = ["pcgs", "trnas", "rrnas"];
const MARKER_LABELS = { pcgs: "PCGs", trnas: "tRNAs", rrnas: "rRNAs" };
const SPLACE_HEATMAP_COLORS = {
    lossStart: [254, 242, 242],
    lossEnd: [191, 0, 48],
    gainStart: [238, 240, 251],
    gainEnd: [58, 71, 171],
};

function getModeCount(values) {
    const counts = new Map();
    values.forEach((value) => counts.set(value, (counts.get(value) || 0) + 1));
    return [...counts.entries()].sort((a, b) => b[1] - a[1] || a[0] - b[0])[0]?.[0] ?? 0;
}

function buildMarkerHeatmapStats(records) {
    const stats = {};
    for (const key of MARKER_KEYS) {
        const values = records.map((record) => countFeatureTypes(record)[key]);
        stats[key] = {
            min: Math.min(...values),
            max: Math.max(...values),
            mode: getModeCount(values),
        };
    }
    return stats;
}

function interpolateRgb(start, end, t) {
    const clamped = Math.max(0, Math.min(1, t));
    return start.map((value, index) => Math.round(value + (end[index] - value) * clamped));
}

function rgbToCss(rgb) {
    return `rgb(${rgb[0]}, ${rgb[1]}, ${rgb[2]})`;
}

function markerHeatmapCellData(value, key, stats) {
    const label = MARKER_LABELS[key] || key;
    const stat = stats?.[key] || { min: value, max: value, mode: value };
    const expected = stat.mode;
    let relation = "normal";
    let intensity = 0;
    let bg = "#f9fafb";
    let color = "#374151";
    let border = "rgba(209, 213, 219, 0.9)";
    let note = "matches the most common marker count in this dataset";

    if (value < expected) {
        relation = "loss";
        intensity = Math.max(0.28, (expected - value) / Math.max(1, expected - stat.min));
        bg = rgbToCss(interpolateRgb(SPLACE_HEATMAP_COLORS.lossStart, SPLACE_HEATMAP_COLORS.lossEnd, intensity));
        border = "rgba(191, 0, 48, 0.55)";
        color = intensity > 0.62 ? "#ffffff" : "#8a001f";
        note = "below the most common count, possible gene loss or missing annotation";
    } else if (value > expected) {
        relation = "duplication";
        intensity = Math.max(0.28, (value - expected) / Math.max(1, stat.max - expected));
        bg = rgbToCss(interpolateRgb(SPLACE_HEATMAP_COLORS.gainStart, SPLACE_HEATMAP_COLORS.gainEnd, intensity));
        border = "rgba(58, 71, 171, 0.55)";
        color = intensity > 0.62 ? "#ffffff" : "#191b4d";
        note = "above the most common count, possible duplicated marker or duplicated annotation";
    }

    return {
        value: escHtml(String(value)),
        relation,
        style: `background:${bg};color:${color};border-color:${border};`,
        tooltip: `${label}: ${value}. Common count in this dataset: ${expected}. This value ${note}.`,
    };
}

function markerHeatmapCell(value, key, stats) {
    const cell = markerHeatmapCellData(value, key, stats);
    return `<span class="marker-heatmap-cell marker-${cell.relation}" style="${cell.style}" data-tippy-content="${escHtml(cell.tooltip)}">${cell.value}</span>`;
}

function formatAverage(value) {
    return Number.isInteger(value) ? String(value) : value.toFixed(1);
}

function getFeatureLengthBp(feature, sequenceLength = 0) {
    const segments = parseLocation(feature?.locationStr || "", sequenceLength || 0);
    if (!segments.length) return 0;
    return segments.reduce((sum, segment) => sum + Math.max(0, segment.end - segment.start), 0);
}

function markerLengthTotals(record) {
    const totals = { pcgs: 0, trnas: 0, rrnas: 0 };
    const sequenceLength = record?.sequence?.length || 0;
    for (const feature of record?.features || []) {
        let key = null;
        if (feature.type === "CDS") key = "pcgs";
        else if (feature.type === "tRNA") key = "trnas";
        else if (feature.type === "rRNA") key = "rrnas";
        if (key) totals[key] += getFeatureLengthBp(feature, sequenceLength);
    }
    return totals;
}

function formatBp(value) {
    const rounded = Math.round(Number(value) || 0);
    return `${rounded.toLocaleString()} bp`;
}

function meanMarkerLength(records, key) {
    if (!records.length) return 0;
    const total = records.reduce((sum, record) => sum + markerLengthTotals(record)[key], 0);
    return total / records.length;
}

function markerMetricsHtml(records, useRealValues = false) {
    return MARKER_KEYS.map((key) => {
        const label = MARKER_LABELS[key] || key;
        if (useRealValues) {
            const values = records.map((record) => markerLengthTotals(record)[key]);
            const uniqueValues = [...new Set(values)].sort((a, b) => a - b);
            const text = uniqueValues.length === 1
                ? formatBp(uniqueValues[0])
                : uniqueValues.map(formatBp).join("/");
            return `<span class="taxonomy-tree-chip real">${label} ${text}</span>`;
        }
        return `<span class="taxonomy-tree-chip">${label} ~${formatBp(meanMarkerLength(records, key))}</span>`;
    }).join("");
}

function hasFetchedTaxonomy(record) {
    const ranks = record.taxonomyRanks || {};
    return Boolean(ranks.kingdom || ranks.phylum || ranks.class || ranks.order || ranks.family || ranks.authorship);
}

function addTaxonomyTreeNode(parent, rankKey, rankLabel, name, record) {
    const cleanName = String(name || "").trim();
    if (!cleanName) return parent;
    const mapKey = `${rankKey}:${cleanName}`;
    if (!parent.children.has(mapKey)) {
        parent.children.set(mapKey, {
            rankKey,
            rankLabel,
            name: cleanName,
            records: [],
            children: new Map(),
        });
    }
    const node = parent.children.get(mapKey);
    node.records.push(record);
    return node;
}

function buildTaxonomyTree(records) {
    const root = { children: new Map(), records };
    for (const record of records) {
        const tax = extractTaxonomy(record);
        const speciesName = `${tax.genus}${tax.species ? " " + tax.species : ""}`.trim() || record.editedOrganism || record.organism || record.accession || "Unknown species";
        const lineage = [
            ["kingdom", "Kingdom", tax.kingdom],
            ["phylum", "Phylum", tax.phylum],
            ["class", "Class", tax.class],
            ["order", "Order", tax.order],
            ["family", "Family", tax.family],
            ["genus", "Genus", tax.genus],
            ["species", "Species", speciesName],
        ];
        let current = root;
        lineage.forEach(([rankKey, rankLabel, name]) => {
            current = addTaxonomyTreeNode(current, rankKey, rankLabel, name, record);
        });
    }
    return root;
}

function renderTaxonomyTreeNode(node) {
    const childNodes = [...node.children.values()].sort((a, b) => a.name.localeCompare(b.name));
    const isSpecies = node.rankKey === "species";
    const isItalicTaxon = node.rankKey === "genus" || node.rankKey === "species";
    const metrics = markerMetricsHtml(node.records, isSpecies);
    const displayName = isItalicTaxon ? `<em>${escHtml(node.name)}</em>` : escHtml(node.name);
    const row = `
        <span class="taxonomy-tree-row">
            <span class="taxonomy-tree-identity"><span class="taxonomy-tree-rank">${escHtml(node.rankLabel)}</span><span class="taxonomy-tree-name">${displayName}</span> <span class="text-gray-400 text-xs">(${node.records.length})</span></span>
            <span class="taxonomy-tree-metrics">${metrics}</span>
        </span>
    `;

    if (!childNodes.length) {
        return `<li class="taxonomy-tree-node">${row}</li>`;
    }

    return `
        <li class="taxonomy-tree-node">
            <details class="taxonomy-tree-details">
                <summary>${row}</summary>
                <ul class="taxonomy-tree-children">${childNodes.map(renderTaxonomyTreeNode).join("")}</ul>
            </details>
        </li>
    `;
}

function renderTaxonomyExplorer() {
    const panel = document.getElementById("taxonomyTreePanel");
    const body = document.getElementById("taxonomyTreeBody");
    if (!panel || !body) return;

    const recordsWithTaxonomy = state.records.filter(hasFetchedTaxonomy);
    if (recordsWithTaxonomy.length === 0) {
        panel.classList.add("hidden");
        body.innerHTML = "";
        return;
    }

    const tree = buildTaxonomyTree(recordsWithTaxonomy);
    const topNodes = [...tree.children.values()].sort((a, b) => a.name.localeCompare(b.name));
    body.innerHTML = `<ul class="taxonomy-tree-list">${topNodes.map(renderTaxonomyTreeNode).join("")}</ul>`;
    panel.classList.remove("hidden");

    // v12: keep the Taxonomic tree panel visible with only top-level nodes shown.
    // The user expands each first-level node to inspect the nested taxonomy.
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
const TAXONOMY_API_HINT_DISMISSED_STORAGE = "splace.taxonomyApiHintDismissed";
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

function syncTaxonomyApiHint() {
    const hint = document.getElementById("taxonomyApiHint");
    if (!hint) return;
    let dismissed = false;
    try {
        dismissed = localStorage.getItem(TAXONOMY_API_HINT_DISMISSED_STORAGE) === "1";
    } catch {
        dismissed = false;
    }
    hint.classList.toggle("hidden", dismissed);
}

function dismissTaxonomyApiHint() {
    try {
        localStorage.setItem(TAXONOMY_API_HINT_DISMISSED_STORAGE, "1");
    } catch {
        // Storage unavailable, ignored
    }
    syncTaxonomyApiHint();
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
    const accessionEntries = accessions
        .map((accession, originalIndex) => ({ accession: accession.trim(), originalIndex }))
        .filter((entry) => entry.accession);
    const resultsByInput = new Array(accessionEntries.length);
    const statusEl = document.getElementById("fetchStatus");
    const localStatusEl = options.statusElementId ? document.getElementById(options.statusElementId) : statusEl;
    localStatusEl.classList.remove("hidden");

    showProgress(
        options.progressTitle || i18nText("step1.progress.fetchGenericTitle"),
        formatProgressCounter(0, accessionEntries.length),
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

    const batchSize = hasUsableNcbiApiKey() ? 10 : 3;
    const batchDelayMs = 1000;
    let completed = 0;

    for (let start = 0; start < accessionEntries.length; start += batchSize) {
        const batch = accessionEntries.slice(start, start + batchSize);
        const batchEnd = start + batch.length;
        localStatusEl.textContent = `Fetching ${start + 1}-${batchEnd}/${accessionEntries.length} records...`;

        await Promise.all(batch.map(async (entry, batchIndex) => {
            const acc = entry.accession;
            const currentLabel = options.progressItems?.[acc] || acc;
            try {
                const record = await fetchAccession(acc);
                resultsByInput[start + batchIndex] = record;
                console.log(`Fetched: ${record.organism || acc} (${record.accession || acc})`);
            } catch (e) {
                console.error(`Failed to fetch ${acc}: ${e.message}`);
            } finally {
                completed++;
                updateProgress(completed, accessionEntries.length, currentLabel);
            }
        }));

        if (start + batchSize < accessionEntries.length) {
            await new Promise((resolve) => setTimeout(resolve, batchDelayMs));
        }
    }

    const results = resultsByInput.filter(Boolean);
    localStatusEl.textContent = `Fetched ${results.length}/${accessionEntries.length} records`;
    setTimeout(() => localStatusEl.classList.add("hidden"), 3000);

    hideProgress(`Loaded ${results.length} record${results.length !== 1 ? 's' : ''}`);
    console.log(`NCBI fetch complete: ${results.length}/${accessionEntries.length} records using batches of ${batchSize}`);

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
    const completeOnly = document.getElementById("taxonomyCompleteCheckbox")?.checked;
    return completeOnly
        ? ["complete genome"]
        : ["complete genome", "partial genome"];
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


function getAccessionVersionFromSummary(summary, fallbackId) {
    return summary?.accessionversion
        || summary?.caption
        || summary?.extra?.gi
        || summary?.uid
        || fallbackId;
}

function extractSpeciesNameFromGenomeTitle(title) {
    const clean = String(title || "")
        .replace(/\b(?:complete|partial|incomplete|unverified)\b.*$/i, "")
        .replace(/\b(?:mitochondrion|mitochondrial|chloroplast|plastid)\b.*$/i, "")
        .trim();

    const match = clean.match(/^([A-Z][a-zA-Z.-]+\s+(?:cf\.\s+|aff\.\s+|sp\.\s+)?[a-z][a-zA-Z.-]+)/);
    return match ? match[1].replace(/\s+/g, " ").trim() : "";
}

function cleanGenomeDescription(title, accession, speciesName) {
    let description = String(title || "-");
    const baseAccession = String(accession || "").replace(/\.\d+$/, "");

    [accession, baseAccession].filter(Boolean).forEach((value) => {
        const safe = value.replace(/[.*+?^${}()|[\]\\]/g, "\\$&");
        description = description.replace(new RegExp(`\\b${safe}\\b`, "g"), "");
    });

    description = description
        .replace(/\b(voucher|isolate|strain|specimen)\s+[^,;:]+/ig, "")
        .replace(/\s{2,}/g, " ")
        .replace(/\s+([,;:])/g, "$1")
        .replace(/^[,;:\s]+|[,;:\s]+$/g, "")
        .trim();

    if (speciesName && !description.toLowerCase().startsWith(speciesName.toLowerCase())) {
        description = `${speciesName} ${description}`.trim();
    }

    return description || speciesName || "-";
}

function normalizeTaxonomyImportFilter(value) {
    return String(value || "")
        .normalize("NFD")
        .replace(/[\u0300-\u036f]/g, "")
        .toLowerCase()
        .trim();
}


function escapeRegExp(value) {
    return String(value || "").replace(/[.*+?^${}()|[\]\\]/g, "\\$&");
}

function highlightTaxonomyImportText(value, query) {
    const text = String(value ?? "");
    const cleanQuery = String(query || "").trim();
    if (!cleanQuery) return escHtml(text);

    const pattern = new RegExp(escapeRegExp(cleanQuery), "ig");
    let lastIndex = 0;
    let html = "";
    for (const match of text.matchAll(pattern)) {
        const start = match.index || 0;
        html += escHtml(text.slice(lastIndex, start));
        html += `<mark class="taxonomy-search-mark">${escHtml(match[0])}</mark>`;
        lastIndex = start + match[0].length;
    }
    html += escHtml(text.slice(lastIndex));
    return html;
}

function getVisibleTaxonomyImportCandidates() {
    const filter = normalizeTaxonomyImportFilter(document.getElementById("taxonomyImportSearchInput")?.value || "");
    if (!filter) return state.taxonomyImportCandidates;

    return state.taxonomyImportCandidates.filter((candidate) => {
        const haystack = normalizeTaxonomyImportFilter([
            candidate.accession,
            candidate.id,
            candidate.speciesName,
            candidate.description,
            candidate.title,
        ].filter(Boolean).join(" "));
        return haystack.includes(filter);
    });
}

function buildTaxonomyImportCandidates(ids, summaries) {
    return ids.map((id, index) => {
        const summary = summaries[index] || {};
        const title = summary.title || "";
        const accession = getAccessionVersionFromSummary(summary, id);
        const speciesName = extractSpeciesNameFromGenomeTitle(title);
        return {
            id,
            accession,
            speciesName,
            description: cleanGenomeDescription(title, accession, speciesName),
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

    const visibleCandidates = getVisibleTaxonomyImportCandidates();
    const rawFilter = (document.getElementById("taxonomyImportSearchInput")?.value || "").trim();
    const filter = normalizeTaxonomyImportFilter(rawFilter);
    const showingHtml = filter
        ? `<div class="mt-2 text-xs text-gray-500">Showing ${visibleCandidates.length}/${total} genomes matching the search.</div>`
        : "";
    countEl.insertAdjacentHTML("beforeend", showingHtml);

    if (visibleCandidates.length === 0) {
        tbody.innerHTML = `<tr><td colspan="4" class="taxonomy-review-empty">No genomes match the current search.</td></tr>`;
        return;
    }

    tbody.innerHTML = visibleCandidates.map((candidate) => {
        const isActive = candidate.selected;
        const label = i18nText(`step1.review.class.${candidate.classification}`);
        const unverifiedLabel = i18nText("step1.review.class.unverified");
        const classificationHtml = candidate.unverified
            ? `<span class="taxonomy-review-badge unverified">${escHtml(unverifiedLabel)}</span>`
            : `<span class="taxonomy-review-badge ${candidate.classification}">${escHtml(label)}</span>`;
        const speciesHtml = candidate.speciesName
            ? `<strong class="taxonomy-review-species"><em>${highlightTaxonomyImportText(candidate.speciesName, rawFilter)}</em></strong>`
            : `<strong class="taxonomy-review-species">${escHtml(i18nText("step1.review.unknownSpecies", "Unknown species"))}</strong>`;
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
                <td class="px-4 py-3 text-gray-700 text-center font-mono">${highlightTaxonomyImportText(candidate.accession, rawFilter)}</td>
                <td class="px-4 py-3 taxonomy-review-description">
                    ${speciesHtml}
                    <p class="${candidate.unverified ? "unverified" : ""}">${highlightTaxonomyImportText(candidate.description || candidate.title || "-", rawFilter)}</p>
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
        ? `<div class="taxonomy-review-meta">${searchDetails.map((detail) => `<span class="step1-filter-chip">${escHtml(detail.label)}: <strong>${escHtml(detail.value)}</strong></span>`).join("")}</div>`
        : "";
    document.getElementById("taxonomyImportSubtitle").innerHTML = `${i18nText("step1.review.subtitle")}${searchDetailsHtml}`;
    const unverifiedCount = summarizeUnverifiedCandidates(candidates);
    const unverifiedHtml = unverifiedCount > 0
        ? `<div class="mt-2"><span class="taxonomy-review-badge unverified">${i18nText("step1.review.class.unverified")}</span> ${i18nText("step1.review.unverifiedCount", { count: unverifiedCount })}<div class="mt-1 text-sm text-splace-red-700">${i18nText("step1.review.unverifiedAdvice")}</div></div>`
        : "";
    document.getElementById("taxonomyImportSummary").innerHTML = `${summaryHtml}${unverifiedHtml}`;
    const searchInput = document.getElementById("taxonomyImportSearchInput");
    if (searchInput) searchInput.value = "";
    state.taxonomyImportCandidates = candidates;
    renderTaxonomyImportCandidates();
    if (typeof window.initTooltips === "function") window.initTooltips();
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
        document.getElementById("taxonomyTreePanel")?.classList.add("hidden");
        document.getElementById("featureTypesSection").classList.add("hidden");
        document.getElementById("genesSection").classList.add("hidden");
        document.getElementById("heatmapSection").classList.add("hidden");
        document.getElementById("markerLengthHeatmapSection")?.classList.add("hidden");
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
    const markerStats = buildMarkerHeatmapStats(sortedIndices.map(({ record }) => record));

    renderTaxonomyExplorer();

    tbody.innerHTML = sortedIndices.map(({ record: r, originalIndex: i }) => {
        const tax = extractTaxonomy(r);
        const counts = countFeatureTypes(r);
        const sourceClass = r.source === 'ncbi' ? 'bg-splace-blue-50 text-splace-blue-600' : 'bg-gray-100 text-gray-600';
        const sourceLabel = r.source === 'ncbi' ? i18nText('step2.record.source.ncbi') : i18nText('step2.record.source.file');
        const abbrevSpecies = tax.genus ? `${tax.genus.charAt(0)}.${tax.species ? ' ' + tax.species : ''}` : '';
        const fullSpecies = `${tax.genus}${tax.species ? ' ' + tax.species : ''}`;
        const displayName = fullSpecies || (r.editedOrganism || r.organism || r.accession || "Unknown");

        const pcgHeat = markerHeatmapCell(counts.pcgs, "pcgs", markerStats);
        const rrnaHeat = markerHeatmapCell(counts.rrnas, "rrnas", markerStats);
        const trnaHeat = markerHeatmapCell(counts.trnas, "trnas", markerStats);

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
                <td class="px-3 py-2 text-center" data-col="pcgs">${pcgHeat}</td>
                <td class="px-3 py-2 text-center" data-col="rrnas">${rrnaHeat}</td>
                <td class="px-3 py-2 text-center" data-col="trnas">${trnaHeat}</td>
                <td class="px-3 py-2 text-center" data-col="source"><span class="text-xs px-1.5 py-0.5 rounded ${sourceClass}">${sourceLabel}</span></td>
                <td class="px-3 py-2 text-right whitespace-nowrap">
                    <button onclick="editRecord(${i})" class="record-action-edit transition-colors leading-none mr-2" data-tippy-content="${i18nText('step2.record.tooltip.edit', { species: `<em>${escHtml(displayName)}</em>`, accession: escHtml(r.accession) })}"><i class="fa-solid fa-pen-to-square"></i></button>
                    <button onclick="removeRecord(${i})" class="record-action-remove transition-colors leading-none" data-tippy-content="${i18nText('step2.record.tooltip.remove', { species: `<em>${escHtml(displayName)}</em>`, accession: escHtml(r.accession) })}"><i class="fa-solid fa-xmark"></i></button>
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

    syncCodonUiDefaults();

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

    // Check for unrecognized gene names (only once per load cycle)
    if (!state.unknownGenesChecked) {
        state.unknownGenesChecked = true;
        const unknowns = collectUnknownGenes();
        if (unknowns.length > 0) {
            openUnknownGenesModal(unknowns);
        }
    }
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

    const countEl = document.getElementById("selectedMarkersCount");
    if (countEl) {
        const n = state.selectedGenes.size;
        if (n > 0) { countEl.textContent = n; countEl.classList.remove("hidden"); }
        else { countEl.classList.add("hidden"); }
    }

    // Render the original gene-presence heatmap and the desktop-only length heatmap.
    renderHeatmap();
    renderMarkerLengthHeatmap();
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

function markerHeatmapSpeciesLabel(record, index) {
    const tax = typeof extractTaxonomy === "function" ? extractTaxonomy(record) : {};
    const speciesName = (tax.genus && tax.species)
        ? `${tax.genus} ${tax.species}`
        : (record.editedOrganism || record.organism || `Record ${index + 1}`);
    const accession = record.accession || "N/A";
    return { speciesName, accession };
}

function collectSelectedMarkerLengthMatrix(selectedGenes) {
    const selectedGeneSet = new Set(selectedGenes);
    const dataType = state.detectedDataType;

    const rows = state.records.map((record, ri) => {
        const label = markerHeatmapSpeciesLabel(record, ri);
        const lengthsByGene = new Map();
        const candidatesPerGene = new Map();
        const assignedTrnas = new Set();
        const recordKey = record.accession || `${label.speciesName}_${ri}`;

        for (const feat of record.features || []) {
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
                geneName = resolveSelectedFeatureGeneName(feat, dataType);
            }

            if (!geneName || !selectedGeneSet.has(geneName)) continue;

            let len = 0;
            try {
                const seq = extractSequence(record.sequence, feat.locationStr);
                len = seq ? seq.length : 0;
            } catch {
                len = getFeatureLengthBp(feat, record.sequence?.length || 0);
            }
            if (!len) continue;

            const list = candidatesPerGene.get(geneName) || [];
            list.push({ len });
            candidatesPerGene.set(geneName, list);
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
            if (chosen) lengthsByGene.set(geneName, chosen.len);
        }

        return { ...label, originalIndex: ri, lengthsByGene };
    }).sort((a, b) => a.speciesName.localeCompare(b.speciesName));

    return rows;
}

function buildMarkerLengthColumnStats(selectedGenes, rows) {
    const stats = new Map();
    selectedGenes.forEach((gene) => {
        const values = rows
            .map((row) => row.lengthsByGene.get(gene))
            .filter((value) => Number.isFinite(value) && value > 0);
        stats.set(gene, {
            min: values.length ? Math.min(...values) : 0,
            max: values.length ? Math.max(...values) : 1,
        });
    });
    return stats;
}

function markerLengthHeatColor(value, stat) {
    if (!Number.isFinite(value) || value <= 0) {
        return { bg: "#f3f4f6", color: "#6b7280", ratio: 0, equalColumn: false };
    }

    const min = stat?.min ?? value;
    const max = stat?.max ?? value;
    if (max === min) {
        return { bg: "#ffffff", color: "#111827", ratio: 0.5, equalColumn: true };
    }

    const ratio = (value - min) / Math.max(1, max - min);
    const rgb = interpolateRgb([191, 0, 48], [58, 71, 171], ratio);
    return {
        bg: rgbToCss(rgb),
        color: "#ffffff",
        ratio,
        equalColumn: false,
    };
}

function renderHighchartsMarkerLengthHeatmap(container, selectedGenes, rows) {
    if (!container || selectedGenes.length === 0 || rows.length === 0) return;
    const columnStats = buildMarkerLengthColumnStats(selectedGenes, rows);
    const data = [];

    rows.forEach((row, y) => {
        selectedGenes.forEach((gene, x) => {
            const value = row.lengthsByGene.has(gene) ? row.lengthsByGene.get(gene) : null;
            const stat = columnStats.get(gene);
            const colorInfo = markerLengthHeatColor(value, stat);
            data.push({
                x,
                y,
                value,
                color: colorInfo.bg,
                dataLabels: { style: { color: colorInfo.color, textOutline: "none" } },
                equalColumn: colorInfo.equalColumn === true,
                columnMin: stat?.min ?? 0,
                columnMax: stat?.max ?? 0,
                speciesName: row.speciesName,
                accession: row.accession,
                gene,
            });
        });
    });

    const manyCells = rows.length * selectedGenes.length > 180;
    const minWidth = Math.max(760, 220 + selectedGenes.length * 92);
    const minHeight = Math.max(360, 120 + rows.length * 28);

    if (container.clientWidth < 10) {
        requestAnimationFrame(() => renderHighchartsMarkerLengthHeatmap(container, selectedGenes, rows));
        return;
    }

    const yAxisSpeciesLabels = rows.map((row) => row.speciesName);
    const yAxisAccessions = rows.map((row) => row.accession);

    Highcharts.chart(container, {
        chart: {
            type: "heatmap",
            height: Math.min(Math.max(360, minHeight), 1200),
            scrollablePlotArea: {
                minWidth,
                minHeight,
            },
            backgroundColor: "#ffffff",
            spacingTop: 18,
            spacingRight: 24,
            spacingBottom: 24,
            spacingLeft: 12,
        },
        title: { text: null },
        subtitle: {
            text: "Each marker column uses its own min and max length scale.",
            align: "left",
            style: { color: "#6b7280", fontSize: "12px" },
        },
        credits: { enabled: false },
        exporting: { enabled: false },
        xAxis: {
            categories: selectedGenes,
            title: { text: null },
            labels: {
                useHTML: false,
                rotation: 0,
                formatter: function () {
                    const gene = selectedGenes[this.pos] || this.value;
                    const stat = columnStats.get(gene) || { min: 0, max: 0 };
                    const range = stat.min && stat.max
                        ? `${Number(stat.min).toLocaleString()}–${Number(stat.max).toLocaleString()} bp`
                        : "no values";
                    return `${gene} (${range})`;
                },
                style: { color: "#191b4d", fontWeight: "600", fontSize: "10px", textOverflow: "ellipsis" },
            },
            opposite: true,
        },
        yAxis: {
            categories: rows.map((row, index) => String(index)),
            title: { text: null },
            reversed: true,
            labels: {
                useHTML: false,
                formatter: function () {
                    const speciesName = yAxisSpeciesLabels[this.pos] || this.value || "";
                    const accession = yAxisAccessions[this.pos] || "";
                    return `<tspan style="font-style:italic">${escHtml(speciesName)}</tspan> <tspan style="fill:#9ca3af">(${escHtml(accession)})</tspan>`;
                },
                style: { color: "#374151", fontSize: "12px", whiteSpace: "nowrap" },
            },
        },
        legend: { enabled: false },
        tooltip: {
            useHTML: true,
            formatter: function () {
                const point = this.point;
                const value = point.value === null || typeof point.value === "undefined"
                    ? "missing"
                    : `${Number(point.value).toLocaleString()} bp`;
                const range = Number(point.columnMin) && Number(point.columnMax)
                    ? `${Number(point.columnMin).toLocaleString()}–${Number(point.columnMax).toLocaleString()} bp`
                    : "no comparable values";
                return `<strong>${escHtml(point.gene)}</strong><br>` +
                    `<em>${escHtml(point.speciesName)}</em><br>` +
                    `${escHtml(point.accession)}<br>` +
                    `<span style="color:#323795;font-weight:700">${value}</span><br>` +
                    `<span style="color:#6b7280">Column min–max: ${range}</span>`;
            },
        },
        plotOptions: {
            series: {
                borderWidth: 1,
                borderColor: "#ffffff",
                nullColor: "#f3f4f6",
                dataLabels: {
                    enabled: !manyCells,
                    formatter: function () {
                        return this.point.value ? Number(this.point.value).toLocaleString() : "";
                    },
                    style: { textOutline: "none", fontSize: "10px", fontWeight: "700" },
                },
            },
        },
        series: [{
            name: "Marker length",
            data,
            turboThreshold: 0,
        }],
        accessibility: {
            enabled: false,
        },
    });
}

function renderFallbackMarkerLengthHeatmap(container, selectedGenes, rows) {
    const columnStats = buildMarkerLengthColumnStats(selectedGenes, rows);

    function valueStyle(value, gene) {
        const stat = columnStats.get(gene);
        const colorInfo = markerLengthHeatColor(value, stat);
        return `background:${colorInfo.bg};color:${colorInfo.color};`;
    }

    let html = `<div class="text-xs text-gray-500 mb-2">Each marker column uses its own min and max length scale.</div>`;
    html += `<div class="marker-length-heatmap-wrap">`;
    html += `<table class="marker-length-heatmap-table" style="min-width:${Math.max(760, 220 + selectedGenes.length * 88)}px;width:100%;">`;
    html += `<thead><tr><th class="marker-length-heatmap-species px-2 py-2 text-left">Species</th>`;
    selectedGenes.forEach((gene) => {
        const stat = columnStats.get(gene) || { min: 0, max: 0 };
        const range = stat.min && stat.max ? `${stat.min.toLocaleString()}–${stat.max.toLocaleString()} bp` : "no values";
        html += `<th class="px-2 py-2 text-center text-xs text-splace-blue-900" data-tippy-content="Column min-max: ${escHtml(range)}"><div class="font-bold">${escHtml(gene)}</div><div class="font-normal text-[10px] text-gray-500">${escHtml(range)}</div></th>`;
    });
    html += `</tr></thead><tbody>`;
    rows.forEach((row) => {
        html += `<tr><td class="marker-length-heatmap-species px-2 py-1.5 whitespace-nowrap"><em>${escHtml(row.speciesName)}</em> <span class="text-gray-400">(${escHtml(row.accession)})</span></td>`;
        selectedGenes.forEach((gene) => {
            const value = row.lengthsByGene.get(gene) || 0;
            const stat = columnStats.get(gene) || { min: 0, max: 0 };
            const range = stat.min && stat.max ? `${stat.min.toLocaleString()}–${stat.max.toLocaleString()} bp` : "no values";
            const colorInfo = markerLengthHeatColor(value, stat);
            const equalClass = colorInfo.equalColumn ? " marker-length-equal" : "";
            html += `<td class="marker-length-heatmap-cell${equalClass} px-2 py-1.5 text-center font-semibold" style="${valueStyle(value, gene)}" data-tippy-content="${escHtml(gene)}: ${value ? value.toLocaleString() + ' bp' : 'missing'}. Column min-max: ${range}">${value ? value.toLocaleString() : "—"}</td>`;
        });
        html += `</tr>`;
    });
    html += `</tbody></table></div>`;
    container.innerHTML = html;
    container.querySelectorAll("[data-tippy-content]").forEach((el) => {
        const htmlContent = el.getAttribute("data-tippy-content") || "";
        const titleText = htmlContent.replace(/<[^>]+>/g, " ").replace(/\s+/g, " ").trim();
        if (titleText) {
            el.setAttribute("title", titleText);
            el.setAttribute("aria-label", titleText);
        }
        window.createTippy && window.createTippy(el, htmlContent);
    });
}



let highchartsAssetLoadPromise = null;
let markerLengthHeatmapRenderToken = 0;

function isDesktopRuntime() {
    return Boolean(window.isElectron && window.electronAPI);
}

function syncRuntimeAttribute() {
    document.documentElement.dataset.runtime = isDesktopRuntime() ? "desktop" : "web";
}

function syncPythonRuntimeCard() {
    const card = document.getElementById("pythonRuntimeCard");
    if (!card) return;
    const isCodon = document.getElementById("alignSeqCodon")?.checked === true;
    card.classList.toggle("hidden", !(isDesktopRuntime() && isCodon));
}

function syncCodonTrimToolUi() {
    const section = document.getElementById("codonTrimToolSection");
    const modeRow = document.getElementById("clipkitModeRow");
    const tool = document.getElementById("codonTrimTool")?.value || "auto";
    const isCodon = window._codonMode === true || document.getElementById("alignSeqCodon")?.checked === true;
    if (section) section.classList.toggle("hidden", !isCodon);
    if (modeRow) modeRow.classList.toggle("hidden", !(isCodon && tool === "clipkit"));
}

async function checkPythonRuntime() {
    const statusEl = document.getElementById("pythonRuntimeStatus");
    if (!statusEl) return;
    if (!window.electronAPI?.checkPythonEnvironment) {
        statusEl.innerHTML = `<span class="font-semibold text-splace-red-700">Desktop Python bridge is unavailable.</span><br><span class="text-gray-500">Replace both <code>electron/main.js</code> and <code>electron/preload.js</code>, then restart SPLACE.</span>`;
        return;
    }
    statusEl.textContent = "Checking Python and ClipKIT…";
    try {
        const result = await window.electronAPI.checkPythonEnvironment();
        if (result?.ipcError) {
            statusEl.innerHTML = `<span class="font-semibold text-splace-red-700">Desktop Python bridge is not registered.</span><br><span class="text-gray-500">${escHtml(result.error || "Replace electron/main.js and restart SPLACE.")}</span>`;
            return;
        }
        const pythonText = result?.pythonFound
            ? `Python ${result.version || "detected"} at ${result.executable || "system PATH"}`
            : "Python was not found.";
        const clipkitText = result?.clipkitFound
            ? `ClipKIT ${result.clipkitVersion || "detected"}`
            : "ClipKIT is not installed.";
        statusEl.innerHTML = `<span class="font-semibold ${result?.pythonFound ? "text-splace-blue-700" : "text-splace-red-700"}">${escHtml(pythonText)}</span><br><span class="${result?.clipkitFound ? "text-splace-blue-700" : "text-amber-700"}">${escHtml(clipkitText)}</span>`;
    } catch (error) {
        statusEl.innerHTML = `<span class="font-semibold text-splace-red-700">Python check failed.</span><br><span class="text-gray-500">${escHtml(error?.message || String(error))}</span><br><span class="text-gray-500">Make sure the updated <code>electron/main.js</code> is in place and restart the Electron app.</span>`;
    }
}

async function installPythonPackagesForSplace() {
    const statusEl = document.getElementById("pythonRuntimeStatus");
    if (!statusEl) return;
    if (!window.electronAPI?.installPythonPackages) {
        statusEl.innerHTML = `<span class="font-semibold text-splace-red-700">Desktop Python bridge is unavailable.</span><br><span class="text-gray-500">Replace both <code>electron/main.js</code> and <code>electron/preload.js</code>, then restart SPLACE.</span>`;
        return;
    }
    statusEl.textContent = "Installing ClipKIT with pip…";
    try {
        const result = await window.electronAPI.installPythonPackages({ packages: ["clipkit"] });
        if (result?.success) {
            statusEl.innerHTML = `<span class="font-semibold text-splace-blue-700">ClipKIT installed or upgraded successfully.</span><br><span class="text-gray-500">${escHtml(result.output || "")}</span>`;
            await checkPythonRuntime();
        } else {
            statusEl.innerHTML = `<span class="font-semibold text-splace-red-700">Install failed.</span><br><span class="text-gray-500 whitespace-pre-wrap">${escHtml(result?.error || result?.output || "Unknown error")}</span>`;
        }
    } catch (error) {
        statusEl.textContent = `Install failed: ${error?.message || error}`;
    }
}

async function ensureClipkitIfRequested(trimTool) {
    if (trimTool !== "clipkit") return true;
    if (!window.electronAPI?.checkPythonEnvironment) return true;
    try {
        const result = await window.electronAPI.checkPythonEnvironment();
        if (result?.pythonFound && result?.clipkitFound) return true;
    } catch (error) {
        alert(`ClipKIT check failed: ${error?.message || error}. Replace electron/main.js and electron/preload.js, then restart SPLACE.`);
        return false;
    }
    alert("ClipKIT requires Python with the clipkit package installed. Use the Check Python and Install ClipKIT buttons in Sequence Type first, or choose trimAl/Auto.");
    return false;
}

function loadScriptOnce(src, globalCheck) {
    if (globalCheck && globalCheck()) return Promise.resolve(true);
    return new Promise((resolve) => {
        if (!src) return resolve(false);
        const existing = document.querySelector(`script[data-splace-src="${CSS.escape(src)}"]`);
        if (existing) {
            existing.addEventListener("load", () => resolve(true), { once: true });
            existing.addEventListener("error", () => resolve(false), { once: true });
            if (globalCheck && globalCheck()) resolve(true);
            return;
        }
        const script = document.createElement("script");
        script.src = src;
        script.dataset.splaceSrc = src;
        script.onload = () => resolve(true);
        script.onerror = () => resolve(false);
        document.head.appendChild(script);
    });
}

async function ensureDesktopHighchartsLoaded() {
    if (window.Highcharts?.chart) return true;
    if (!isDesktopRuntime() || !window.electronAPI) return false;

    if (!highchartsAssetLoadPromise) {
        highchartsAssetLoadPromise = (async () => {
            let lastError = null;

            async function tryAssetUrls() {
                if (!window.electronAPI.getHighchartsAssetUrls) return false;
                const urls = await window.electronAPI.getHighchartsAssetUrls();
                const loadedCore = await loadScriptOnce(urls?.highcharts, () => Boolean(window.Highcharts?.chart));
                if (!loadedCore) return false;
                await loadScriptOnce(urls?.heatmap, () => Boolean(window.Highcharts?.seriesTypes?.heatmap));
                await loadScriptOnce(urls?.accessibility, () => Boolean(window.Highcharts?.AST || window.Highcharts?.A11yChartUtilities));
                return Boolean(window.Highcharts?.chart && window.Highcharts?.seriesTypes?.heatmap);
            }

            async function tryInlineAssets() {
                if (!window.electronAPI.loadHighchartsAssets) return false;
                const assets = await window.electronAPI.loadHighchartsAssets();
                const scripts = [assets?.highcharts, assets?.heatmap, assets?.accessibility];
                for (const code of scripts) {
                    if (!code) continue;
                    const blob = new Blob([code], { type: "text/javascript" });
                    const url = URL.createObjectURL(blob);
                    await loadScriptOnce(url, null);
                    URL.revokeObjectURL(url);
                }
                return Boolean(window.Highcharts?.chart && window.Highcharts?.seriesTypes?.heatmap);
            }

            try {
                const ok = await tryAssetUrls();
                if (ok) return true;
            } catch (error) {
                lastError = error;
                console.warn("Highcharts asset URL loading failed; trying inline desktop assets.", error);
            }

            try {
                const ok = await tryInlineAssets();
                if (ok) return true;
            } catch (error) {
                lastError = error;
                console.warn("Highcharts inline loading failed.", error);
            }

            if (lastError) {
                console.warn("Highcharts could not be loaded from the desktop package.", lastError);
            }
            return false;
        })();
    }
    return highchartsAssetLoadPromise;
}

async function renderMarkerLengthHeatmap() {
    syncRuntimeAttribute();
    const section = document.getElementById("markerLengthHeatmapSection");
    const container = document.getElementById("markerLengthHeatmapContainer");
    if (!section || !container) return;

    const selectedGenes = [...state.selectedGenes].sort();
    if (!isDesktopRuntime() || selectedGenes.length === 0 || state.records.length === 0) {
        section.classList.add("hidden");
        container.innerHTML = "";
        return;
    }

    const renderToken = ++markerLengthHeatmapRenderToken;
    const rows = collectSelectedMarkerLengthMatrix(selectedGenes);
    if (rows.length === 0) {
        section.classList.add("hidden");
        container.innerHTML = "";
        return;
    }
    section.classList.remove("hidden");
    container.innerHTML = `<div class="text-sm text-gray-500 p-4">Loading desktop Highcharts heatmap…</div>`;

    const highchartsReady = await ensureDesktopHighchartsLoaded();
    if (renderToken !== markerLengthHeatmapRenderToken) return;

    if (highchartsReady) {
        renderHighchartsMarkerLengthHeatmap(container, selectedGenes, rows);
    } else {
        container.innerHTML = `<div class="text-sm text-splace-red-700 bg-splace-red-50 border border-splace-red-100 rounded-xl p-3 mb-3">Highcharts is not available in the desktop package. Run <code>npm install highcharts --save</code> in the Electron folder and restart SPLACE.</div>`;
        renderFallbackMarkerLengthHeatmap(container, selectedGenes, rows);
    }
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
        `Remove <strong><em>${tax.genus} ${tax.species}</em></strong> (${r.accession})? This action cannot be undone.`;
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
    if (btn) btn.disabled = raw.trim().length < 3;

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
    let added = 0;
    for (const record of records) {
        const accession = (record.accession || "").trim();
        if (accession && seen.has(accession)) continue;
        state.records.push(record);
        if (accession) seen.add(accession);
        added++;
    }
    if (added > 0) {
        // New records may introduce unknown gene names — re-trigger modal check
        state.unknownGenesChecked = false;
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

function getDefaultCodonFallbackGeneticCode() {
    return state.detectedDataType === "cp" ? "11" : "2";
}

function syncCodonUiDefaults(force = false) {
    const geneticCodeSelect = document.getElementById("alignGeneticCode");
    if (!geneticCodeSelect) return;

    if (force || geneticCodeSelect.dataset.userSelected !== "true") {
        geneticCodeSelect.value = getDefaultCodonFallbackGeneticCode();
    }
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
    console.log(`Taxonomy fetched: ${successCount}/${totalUnique} species${errorCount > 0 ? `, ${errorCount} failed` : ""}`);
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
    try { renderHeaderBuilder(); } catch (e) { console.error("Header builder render error:", e); }
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

// -----------------------------------------------------------------------
// Codon-aware helpers: CDS validation + JS translation
// -----------------------------------------------------------------------

function getGeneticCodeName(code) {
    const select = document.getElementById("alignGeneticCode");
    const option = select?.querySelector(`option[value="${code}"]`);
    return option ? option.textContent.trim() : `NCBI ${code}`;
}

function getCodonRunOptions() {
    const fallbackGeneticCode = String(
        parseInt(document.getElementById("alignGeneticCode")?.value, 10) ||
        parseInt(getDefaultCodonFallbackGeneticCode(), 10)
    );

    return {
        fallbackGeneticCode,
        fallbackGeneticCodeName: getGeneticCodeName(fallbackGeneticCode),
        allowPseudo: document.getElementById("codonAllowPseudo")?.checked === true,
        allowBacktranslationWarnings: document.getElementById("codonAllowBacktranslationWarnings")?.checked === true,
        internalStopPolicy: document.getElementById("codonAllowInternalStops")?.checked === true
            ? "allow-with-warning"
            : "exclude",
    };
}

function normalizeQualifierText(value) {
    return typeof value === "string" ? value.replace(/\s+/g, " ").trim() : "";
}

function resolveCodonGeneName(feat, dataType) {
    const rawGene = normalizeQualifierText(feat.qualifiers.gene);
    const rawProduct = normalizeQualifierText(feat.qualifiers.product);
    const directMatch = standardizeGeneName(rawGene || null, rawProduct || null, dataType);
    if (directMatch) return directMatch;

    const fallbackFields = [
        normalizeQualifierText(feat.qualifiers.note),
        normalizeQualifierText(feat.qualifiers.locus_tag),
        normalizeQualifierText(feat.qualifiers.gene_synonym),
    ].filter(Boolean);

    for (const fallbackField of fallbackFields) {
        const fallbackMatch = standardizeGeneName(rawGene || fallbackField, fallbackField, dataType);
        if (fallbackMatch) return fallbackMatch;
    }

    return rawGene || rawProduct || fallbackFields[0] || null;
}

function resolveSelectedFeatureGeneName(feat, dataType) {
    if (feat.type === "tRNA") return standardizeTrnaName(feat);
    return resolveCodonGeneName(feat, dataType);
}

function isPseudoCdsFeature(feat) {
    const qualifiers = feat?.qualifiers || {};
    return Object.prototype.hasOwnProperty.call(qualifiers, "pseudo") ||
        Object.prototype.hasOwnProperty.call(qualifiers, "pseudogene");
}

function resolveCodonFeatureGeneticCode(feat, fallbackGeneticCode) {
    const featureCode = parseInt(feat.qualifiers.transl_table, 10);
    if (featureCode) return featureCode;

    const parsedFallback = parseInt(fallbackGeneticCode, 10);
    if (parsedFallback) return parsedFallback;

    return parseInt(getDefaultCodonFallbackGeneticCode(), 10);
}

/**
 * Validate a single CDS DNA sequence and translate it to protein.
 * Returns { validDna, protein } where validDna has the terminal stop codon
 * removed (if present), so protein length * 3 === validDna.length.
 * Returns null if validation fails.
 */
function validateAndTranslateCdsForCodon(geneName, seqHeader, dna, translTable, options, logFn) {
    if (!dna || dna.length === 0) {
        logFn(`[CDS validation] [${geneName}] SKIP "${seqHeader}": empty sequence`);
        return null;
    }

    const fallbackCode = String(parseInt(options?.fallbackGeneticCode, 10) || parseInt(getDefaultCodonFallbackGeneticCode(), 10));
    const internalStopPolicy = options?.internalStopPolicy || "exclude";

    let tableKey = String(translTable || fallbackCode);
    let table = GENETIC_CODE_TABLES[tableKey] || GENETIC_CODE_TABLES[fallbackCode] || GENETIC_CODE_TABLES["1"];
    if (!GENETIC_CODE_TABLES[tableKey]) {
        logFn(`[Translation] [${geneName}] WARN "${seqHeader}": unknown transl_table ${translTable}, using fallback ${getGeneticCodeName(fallbackCode)}`);
        tableKey = fallbackCode;
        table = GENETIC_CODE_TABLES[tableKey] || GENETIC_CODE_TABLES["1"];
    }

    let workDna = dna.toUpperCase().replace(/U/g, "T");
    let normalized = false;

    // Normalize incomplete terminal stop codon fragments common in mitochondrial
    // genes where the stop is completed by polyadenylation: T→TAA, TA→TAA.
    if (workDna.length % 3 !== 0) {
        const remainder = workDna.length % 3;
        const terminal = workDna.slice(-remainder);
        if (remainder === 1 && terminal === "T") {
            logFn(`[CDS validation] [${geneName}] "${seqHeader}": removed terminal incomplete stop codon fragment T (${dna.length} bp → ${dna.length - 1} bp)`);
            workDna = workDna.slice(0, -1);
            normalized = true;
        } else if (remainder === 2 && terminal === "TA") {
            logFn(`[CDS validation] [${geneName}] "${seqHeader}": removed terminal incomplete stop codon fragment TA (${dna.length} bp → ${dna.length - 2} bp)`);
            workDna = workDna.slice(0, -2);
            normalized = true;
        } else {
            logFn(`[CDS validation] [${geneName}] SKIP "${seqHeader}": length ${dna.length} bp is not a multiple of 3 (terminal "${terminal}" is not a recognizable incomplete stop fragment)`);
            return null;
        }
    }

    if (workDna.length === 0) {
        logFn(`[CDS validation] [${geneName}] SKIP "${seqHeader}": sequence empty after normalization`);
        return null;
    }

    let protein = "";
    let internalStops = 0;

    // Translate all codons except the last (check for internal stops)
    for (let i = 0; i < workDna.length - 3; i += 3) {
        const codon = workDna.slice(i, i + 3);
        const aa = table[codon] || "X";
        if (aa === "*") internalStops++;
        protein += aa;
    }

    // Last codon — check if it's a terminal stop and remove it from the DNA
    const lastCodon = workDna.slice(-3);
    const lastAa = table[lastCodon] || "X";
    let validDna;

    if (lastAa === "*") {
        validDna = workDna.slice(0, -3);
        logFn(`[CDS validation] [${geneName}] "${seqHeader}": removed terminal stop codon ${lastCodon} (${workDna.length / 3} codons, code ${tableKey})`);
        // protein already excludes the last codon
    } else {
        protein += lastAa;
        validDna = workDna;
        logFn(`[CDS validation] [${geneName}] "${seqHeader}": ${workDna.length / 3} codons, no terminal stop (code ${tableKey})`);
    }

    if (!validDna || validDna.length === 0) {
        logFn(`[CDS validation] [${geneName}] SKIP "${seqHeader}": no validated CDS remains after terminal stop removal`);
        return null;
    }

    if (validDna.length % 3 !== 0) {
        logFn(`[CDS validation] [${geneName}] SKIP "${seqHeader}": validated CDS length ${validDna.length} bp is not divisible by 3`);
        return null;
    }

    if (internalStops > 0) {
        if (internalStopPolicy === "exclude") {
            logFn(`[CDS validation] [${geneName}] SKIP "${seqHeader}": ${internalStops} internal stop codon(s) detected (policy: exclude)`);
            return null;
        }
        logFn(`[Translation] [${geneName}] WARN "${seqHeader}": ${internalStops} internal stop codon(s) detected (policy: allow with warning)`);
    }

    return {
        validDna,
        protein,
        normalized,
        terminalStopRemoved: lastAa === "*",
        internalStops,
        internalStopPolicy,
        translTableUsed: tableKey,
    };
}

// Build perGeneCandidates Map from state.records (first phase of codon validation).
function buildCodonCandidates(options, logFn) {
    const dataType = state.detectedDataType;
    const perGeneCandidates = new Map();
    const selectedGeneTypes = new Map();

    for (const record of state.records) {
        const candidatesPerGene = new Map();
        const tax = extractTaxonomy(record);
        const headerName = buildHeader(record, tax);
        const recordKey = record.accession || headerName;

        for (const feat of record.features) {
            const selectedGeneName = resolveSelectedFeatureGeneName(feat, dataType);
            if (selectedGeneName && state.selectedGenes.has(selectedGeneName)) {
                const types = selectedGeneTypes.get(selectedGeneName) || new Set();
                types.add(feat.type);
                selectedGeneTypes.set(selectedGeneName, types);
            }

            if (feat.type !== "CDS") continue;

            const geneName = resolveCodonGeneName(feat, dataType);
            if (!geneName || !state.selectedGenes.has(geneName)) continue;

            if (isPseudoCdsFeature(feat) && !options.allowPseudo) {
                logFn(`[CDS extraction] [${geneName}] SKIP "${headerName}": /pseudo or /pseudogene CDS excluded by policy`);
                continue;
            }

            const rawSeq = extractSequence(record.sequence, feat.locationStr);
            if (!rawSeq) {
                logFn(`[CDS extraction] [${geneName}] SKIP "${headerName}": unable to extract CDS sequence from GenBank location`);
                continue;
            }

            const codonStart = parseInt(feat.qualifiers.codon_start, 10) || 1;
            const offset = Math.max(0, codonStart - 1);
            const seq = offset > 0 ? rawSeq.slice(offset) : rawSeq;
            if (offset > 0) {
                logFn(`[CDS extraction] [${geneName}] "${headerName}": /codon_start=${codonStart}, trimmed ${offset} bp`);
            }

            const translTable = resolveCodonFeatureGeneticCode(feat, options.fallbackGeneticCode);
            const translTableSource = parseInt(feat.qualifiers.transl_table, 10)
                ? "/transl_table"
                : `fallback ${getGeneticCodeName(options.fallbackGeneticCode)}`;
            const existing = candidatesPerGene.get(geneName) || [];
            existing.push({ header: headerName, seq, translTable, translTableSource, len: seq.length });
            candidatesPerGene.set(geneName, existing);
        }

        for (const [geneName, candidates] of candidatesPerGene) {
            const choiceKey = getDuplicateChoiceKey(recordKey, geneName);
            const longestIndex = candidates.reduce((bi, c, i, all) => c.len > all[bi].len ? i : bi, 0);
            const selectedIndex = state.duplicateGeneChoices.has(choiceKey)
                ? state.duplicateGeneChoices.get(choiceKey)
                : longestIndex;
            const chosen = candidates[selectedIndex] || candidates[longestIndex];
            if (!chosen) continue;

            const byGene = perGeneCandidates.get(geneName) || [];
            byGene.push(chosen);
            perGeneCandidates.set(geneName, byGene);
        }
    }

    const cdsSelectedMarkers = [...state.selectedGenes]
        .filter((geneName) => selectedGeneTypes.get(geneName)?.has("CDS"))
        .sort();
    const nonCdsMarkers = [...state.selectedGenes]
        .filter((geneName) => {
            const featureTypes = selectedGeneTypes.get(geneName);
            return featureTypes && !featureTypes.has("CDS");
        })
        .sort();

    return {
        perGeneCandidates,
        summary: {
            cdsSelectedMarkers,
            nonCdsMarkers,
        },
    };
}

// Validate and translate one gene entry from perGeneCandidates.
// Returns { cdsDnaFasta, proteinFasta, stats } or null if not enough valid seqs.
function processCodonGeneEntry(geneName, seqs, options, logFn) {
    const headerCounts = {};
    for (const s of seqs) headerCounts[s.header] = (headerCounts[s.header] || 0) + 1;
    const dupHeaders = Object.entries(headerCounts).filter(([, c]) => c > 1).map(([h]) => h);
    if (dupHeaders.length) logFn(`[CDS extraction] [${geneName}] WARN: duplicate headers detected — ${dupHeaders.join("; ")}`);

    const codesUsed = [...new Set(seqs.map((s) => String(s.translTable)))];
    const codeNames = codesUsed.map((code) => getGeneticCodeName(code));
    logFn(`[CDS extraction] [${geneName}] genetic code(s): ${codeNames.join(", ")} (${seqs.length} seq)`);

    let cdsDnaFasta = "", proteinFasta = "";
    const sequenceMeta = {};
    let validCount = 0, skipCount = 0, normalizedCount = 0;
    let terminalStopsRemoved = 0, incompleteStopFragmentsRemoved = 0, internalStopsDetected = 0;

    for (const s of seqs) {
        const validated = validateAndTranslateCdsForCodon(geneName, s.header, s.seq, s.translTable, options, logFn);
        if (!validated) { skipCount++; continue; }
        if (validated.normalized) { normalizedCount++; incompleteStopFragmentsRemoved++; }
        if (validated.terminalStopRemoved) terminalStopsRemoved++;
        internalStopsDetected += validated.internalStops || 0;
        cdsDnaFasta += `>${s.header}\n${validated.validDna}\n`;
        proteinFasta += `>${s.header}\n${validated.protein}\n`;
        sequenceMeta[s.header] = {
            geneticCode: validated.translTableUsed,
            geneticCodeName: getGeneticCodeName(validated.translTableUsed),
            protein: validated.protein,
            validDnaLength: validated.validDna.length,
        };
        validCount++;
    }

    const rawCount = seqs.length;
    const stats = {
        rawCount,
        validCount,
        skippedCount: skipCount,
        normalizedCount,
        terminalStopsRemoved,
        incompleteStopFragmentsRemoved,
        internalStopsDetected,
        internalStopPolicy: options.internalStopPolicy,
        geneticCodesUsed: codesUsed,
        geneticCode: codesUsed.length === 1 ? codesUsed[0] : String(options.fallbackGeneticCode),
        geneticCodeName: codesUsed.length === 1
            ? getGeneticCodeName(codesUsed[0])
            : codeNames.join(", "),
    };

    if (validCount < 2) {
        logFn(`[CDS validation] [${geneName}] SKIP: ${rawCount} raw | ${validCount} valid | ${skipCount} skipped | ${normalizedCount} normalized (minimum 2 required)`);
        return { eligible: false, cdsDnaFasta, proteinFasta, stats, sequenceMeta };
    }
    logFn(`[CDS validation] [${geneName}] OK: ${rawCount} raw | ${validCount} valid | ${skipCount} skipped | ${normalizedCount} normalized`);
    return { eligible: true, cdsDnaFasta, proteinFasta, stats, sequenceMeta };
}

function generateCodonFastaPair(defaultGeneticCode, logFn) {
    const options = typeof defaultGeneticCode === "object"
        ? defaultGeneticCode
        : {
            ...getCodonRunOptions(),
            fallbackGeneticCode: String(defaultGeneticCode || getDefaultCodonFallbackGeneticCode()),
            fallbackGeneticCodeName: getGeneticCodeName(String(defaultGeneticCode || getDefaultCodonFallbackGeneticCode())),
        };
    const candidateContext = buildCodonCandidates(options, logFn);
    const result = new Map();
    for (const [geneName, seqs] of candidateContext.perGeneCandidates) {
        const pair = processCodonGeneEntry(geneName, seqs, options, logFn);
        if (pair?.eligible) result.set(geneName, pair);
    }
    return result;
}

window.runCodonAwareSelfTest = function () {
    const silentLog = () => { };
    const mtOptions = { fallbackGeneticCode: "2", internalStopPolicy: "exclude" };
    const cpOptions = { fallbackGeneticCode: "11", internalStopPolicy: "exclude" };

    const terminalStop = validateAndTranslateCdsForCodon("TEST", "terminal-stop", "ATGTAG", 2, mtOptions, silentLog);
    const incompleteStop = validateAndTranslateCdsForCodon("TEST", "incomplete-stop", "ATGTA", 2, mtOptions, silentLog);
    const mtTranslation = validateAndTranslateCdsForCodon("TEST", "mt-code", "ATATGATAA", 2, mtOptions, silentLog);
    const cpTranslation = validateAndTranslateCdsForCodon("TEST", "cp-code", "ATATTTTAA", 11, cpOptions, silentLog);
    const backtranslated = backTranslateProteinAlignment("M-F", "ATGTTT");
    const integrity = getCodonIntegrityChecks([
        { name: "tax1", seq: backtranslated },
        { name: "tax2", seq: "ATA---TTT" },
    ]);
    const partitions = [{ gene: "geneA", start: 1, end: 9, len: 9 }];
    const noThirdSeq = removeThirdCodonPositionsFromSequence("ATGCCCTTT");
    const concatFixture = {
        alignedSeqs: [
            { name: "tax1", seq: "ATGCCCTTT" },
            { name: "tax2", seq: "ATA---TTT" },
        ],
        partitions,
        geneOrder: ["geneA"],
        geneSeqs: {
            geneA: {
                tax1: "ATGCCCTTT",
                tax2: "ATA---TTT",
            },
        },
        geneLens: { geneA: 9 },
        allSpecies: ["tax1", "tax2"],
        totalLen: 9,
        missingReport: [],
    };
    const noThirdMatrix = buildNoThirdCodonMatrix(concatFixture);
    const nexus = buildNexusString([
        { name: "tax1", seq: "ATG---TTC" },
        { name: "tax2", seq: "ATA---TTT" },
    ], partitions, [], "dna", true);
    const partitionText = buildPartitionString(partitions, true, "dna");
    const positions12Nexus = buildNexusString(noThirdMatrix.alignedSeqs, noThirdMatrix.partitions, [], "dna", "positions12");
    const positions12Partition = buildPartitionString(noThirdMatrix.partitions, "positions12", "dna");

    let modalFooterVisible = false;
    const analysisModal = document.getElementById("analysisModal");
    const modalShell = analysisModal?.querySelector(".analysis-modal-shell");
    const modalMain = analysisModal?.querySelector(".analysis-modal-main");
    const modalFooter = document.getElementById("analysisModalFooter");
    if (analysisModal && modalShell && modalMain && modalFooter) {
        const previousHidden = analysisModal.classList.contains("hidden");
        const previousHeight = modalShell.style.height;
        const previousMaxHeight = modalShell.style.maxHeight;
        analysisModal.classList.remove("hidden");
        modalShell.style.height = "240px";
        modalShell.style.maxHeight = "240px";
        refreshManagedModalState("analysisModal");
        const shellRect = modalShell.getBoundingClientRect();
        const footerRect = modalFooter.getBoundingClientRect();
        modalFooterVisible = !modalMain.contains(modalFooter)
            && footerRect.top >= shellRect.top - 1
            && footerRect.bottom <= shellRect.bottom + 1;
        modalShell.style.height = previousHeight;
        modalShell.style.maxHeight = previousMaxHeight;
        if (previousHidden) analysisModal.classList.add("hidden");
        refreshManagedModalState("analysisModal");
    }

    const results = [
        { name: "terminal stop removal", pass: terminalStop?.validDna === "ATG" && terminalStop?.protein === "M" },
        { name: "incomplete stop fragment removal", pass: incompleteStop?.validDna === "ATG" && incompleteStop?.protein === "M" },
        { name: "NCBI 2 translation", pass: mtTranslation?.protein === "MW" },
        { name: "NCBI 11 translation", pass: cpTranslation?.protein === "IF" },
        { name: "back-translation length divisible by 3", pass: backtranslated.length % 3 === 0 },
        { name: "back-translation gaps in triplets", pass: integrity.gapBlocksAreTriplets && integrity.equalLengths },
        { name: "codon-position charsets", pass: nexus.includes("charset geneA_pos1 = 1-9\\3;") && partitionText.includes("DNA, geneA_pos3 = 3-9\\3") },
        { name: "remove third codon positions", pass: noThirdSeq === "ATCCTT" },
        { name: "no-third matrix is two thirds of original length", pass: noThirdMatrix.originalLength === 9 && noThirdMatrix.noThirdLength === 6 && noThirdMatrix.noThirdLength * 3 === noThirdMatrix.originalLength * 2 },
        { name: "positions12 partition output", pass: positions12Nexus.includes("charset geneA_pos1 = 1-6\\2;") && positions12Nexus.includes("charset geneA_pos2 = 2-6\\2;") && !positions12Nexus.includes("pos3") && positions12Partition.includes("DNA, geneA_pos2 = 2-6\\2") && !positions12Partition.includes("pos3") },
        { name: "analysis modal footer stays visible in compact height", pass: modalFooterVisible },
    ];

    const summary = {
        pass: results.every((result) => result.pass),
        results,
        artifacts: {
            terminalStopProtein: terminalStop?.protein || null,
            mtTranslationProtein: mtTranslation?.protein || null,
            cpTranslationProtein: cpTranslation?.protein || null,
            backtranslated,
            noThirdSeq,
            noThirdMatrix,
            nexus,
            partitionText,
            positions12Nexus,
            positions12Partition,
            modalFooterVisible,
        },
    };

    console.table(results);
    console.log("[codon-self-test]", summary);
    return summary;
};

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
    console.log(`Downloaded ZIP: ${files.size} gene${files.size !== 1 ? 's' : ''}`);
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
                console.log(`Parsed: ${record.organism || record.accession || file.name} (${file.name})`);
            } catch {
                console.error(`Failed to parse ${file.name}`);
            }

            loaded++;
            updateProgress(loaded, validFiles.length, file.name);
            if (loaded === validFiles.length) {
                renderRecords();
                hideProgress(`Loaded ${validFiles.length} file${validFiles.length !== 1 ? 's' : ''}`);
                console.log(`${validFiles.length} file${validFiles.length !== 1 ? 's' : ''} loaded`);
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

    registerManagedModal({
        modalId: "analysisModal",
        scrollRegionId: "analysisModalBody",
        footerId: "analysisModalFooter",
        scrollCueId: "analysisModalScrollCue",
        closeButtonId: "analysisModalClose",
    });
    registerManagedModal({
        modalId: "iqtreeModal",
        scrollRegionId: "iqtreeModalBody",
        footerId: "iqtreeModalFooter",
        scrollCueId: "iqtreeModalScrollCue",
        closeButtonId: "iqtreeModalClose",
    });
    document.addEventListener("keydown", (event) => {
        if (event.key !== "Escape") return;
        if (closeManagedModal("iqtreeModal")) {
            event.preventDefault();
            return;
        }
        if (closeManagedModal("analysisModal")) {
            event.preventDefault();
        }
    });

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
    document.getElementById("taxonomyApiHintClose")?.addEventListener("click", dismissTaxonomyApiHint);
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
    document.getElementById("taxonomyImportSearchInput")?.addEventListener("input", renderTaxonomyImportCandidates);
    document.getElementById("taxonomyImportConfirm").addEventListener("click", () => {
        const selectedIds = state.taxonomyImportCandidates
            .filter((candidate) => candidate.selected)
            .map((candidate) => candidate.id);
        closeTaxonomyImportModal(selectedIds);
    });
    document.getElementById("checkPythonBtn")?.addEventListener("click", checkPythonRuntime);
    document.getElementById("installPythonPackagesBtn")?.addEventListener("click", installPythonPackagesForSplace);
    document.getElementById("codonTrimTool")?.addEventListener("change", syncCodonTrimToolUi);
    document.getElementById("taxonomyTreePanel")?.addEventListener("toggle", () => {
        // The panel remains user-controlled after opening; new taxonomy renders reset it closed.
    });

    document.getElementById("accessionInput").addEventListener("input", validateAccessionInput);
    document.getElementById("taxonomySearchInput").addEventListener("input", validateTaxonomySearchInput);
    document.getElementById("taxonomyOrganelleSelect").addEventListener("change", validateTaxonomySearchInput);
    document.getElementById("taxonomyRefseqOnly").addEventListener("change", validateTaxonomySearchInput);
    document.getElementById("taxonomyCompleteCheckbox").addEventListener("change", validateTaxonomySearchInput);
    document.getElementById("taxonomyPartialCheckbox")?.addEventListener("change", validateTaxonomySearchInput);

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
    syncTaxonomyApiHint();
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
        state.geneNameOverrides.clear();
        state.unknownGenesChecked = false;
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
            document.getElementById("panel6seq").classList.remove("hidden");
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

    // Sequence-type radio toggle (Nucleotide / Amino Acid / Codon Aware)
    document.querySelectorAll('input[name="alignSeqType"]').forEach(radio => {
        radio.addEventListener("change", () => {
            const isAa = document.getElementById("alignSeqAa")?.checked;
            const isCodon = document.getElementById("alignSeqCodon")?.checked;
            const isNt = !isAa && !isCodon;
            const codeSection = document.getElementById("alignGeneticCodeSection");
            const codonWarning = document.getElementById("codonModeWarning");
            const codonPolicy = document.getElementById("codonPolicySection");
            if (codeSection) codeSection.classList.toggle("hidden", !(isAa || isCodon));
            if (codonWarning) codonWarning.classList.toggle("hidden", !isCodon);
            if (codonPolicy) codonPolicy.classList.toggle("hidden", !isCodon);
            if (isCodon) syncCodonUiDefaults();
            syncPythonRuntimeCard();
            syncCodonTrimToolUi();
            // Visual highlight on the label cards
            const ntLabel = document.getElementById("alignSeqNtLabel");
            const aaLabel = document.getElementById("alignSeqAaLabel");
            const codonLabel = document.getElementById("alignSeqCodonLabel");
            [[ntLabel, isNt], [aaLabel, isAa], [codonLabel, isCodon]].forEach(([el, active]) => {
                if (!el) return;
                el.classList.toggle("border-splace-blue-400", !!active);
                el.classList.toggle("bg-splace-blue-50", !!active);
                el.classList.toggle("border-gray-200", !active);
            });
        });
    });

    const geneticCodeSelect = document.getElementById("alignGeneticCode");
    if (geneticCodeSelect) {
        syncCodonUiDefaults();
        geneticCodeSelect.addEventListener("change", () => {
            geneticCodeSelect.dataset.userSelected = "true";
        });
    }

    // Re-render header builder when separator changes
    document.querySelectorAll('input[name="headerSeparator"]').forEach(radio => {
        radio.addEventListener("change", () => renderHeaderBuilder());
    });

    // Initialize Tippy on static elements
    window.initTooltips && window.initTooltips();
    syncRuntimeAttribute();
    syncPythonRuntimeCard();
    syncCodonTrimToolUi();

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

        // Ask main process for CPU count; set default threads to 4
        window.electronAPI.getCpuCount().then(cpus => {
            const threads = 4;
            document.getElementById("mafftThreads").value = threads;
            const concurrent = Math.max(1, Math.floor((cpus - 2) / threads));
            document.getElementById("mafftConcurrencyInfo").textContent =
                `${cpus} logical CPUs — using ${threads} threads (${concurrent} job${concurrent !== 1 ? 's' : ''} in parallel)`;
        });
        document.getElementById("mafftThreads").addEventListener("input", () => {
            window.electronAPI.getCpuCount().then(cpus => {
                const threads = Math.max(1, parseInt(document.getElementById("mafftThreads").value) || 1);
                const concurrent = Math.max(1, Math.floor((cpus - 2) / threads));
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
            syncCodonTrimToolUi();
            document.getElementById("panel6d").classList.remove("hidden");
            document.getElementById("panel6d").scrollIntoView({ behavior: "smooth", block: "start" });
        });
        document.getElementById("showTrimalBtn").addEventListener("click", () => {
            document.getElementById("showTrimalBtn").classList.add("hidden");
            syncCodonTrimToolUi();
            document.getElementById("panel6d").classList.remove("hidden");
            document.getElementById("panel6d").scrollIntoView({ behavior: "smooth", block: "start" });
        });

        // Run trimAl
        document.getElementById("runTrimalBtn").addEventListener("click", async () => {
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

            const isCodonMode = window._codonMode === true;
            const trimTool = isCodonMode ? (document.getElementById("codonTrimTool")?.value || "auto") : "trimal";
            const clipkitMode = document.getElementById("clipkitMode")?.value || "gappy";
            if (!(await ensureClipkitIfRequested(trimTool))) return;
            const markers = isCodonMode
                ? Object.keys(window._codonOriginalDna || {})
                : [...state.selectedGenes];
            openAnalysisModal("trimAl", markers, {
                alignMode: window._alignmentMode || (isCodonMode ? "codon" : "nt"),
                isTranslated: isCodonMode || window._alignmentMode === "aa",
                mafftMode: null,
                postProcess: isCodonMode ? (trimTool === "clipkit" ? "ClipKIT protein trim + codon back-translation" : "trimAl/ClipKIT codon back-translation") : "trimAl",
            });
            window.electronAPI.runTrimal({
                markers,
                params,
                alignmentMode: window._alignmentMode || (isCodonMode ? "codon" : "nt"),
                codonMode: isCodonMode,
                cdsFastas: isCodonMode ? window._codonOriginalDna : undefined,
                trimTool,
                clipkitMode,
            });
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
        document.getElementById("runAlignmentBtn").addEventListener("click", async () => {
            if (state.selectedGenes.size === 0) { alert("No markers selected."); return; }

            const isCodonMode = document.getElementById("alignSeqCodon")?.checked === true;
            const isAaMode = !isCodonMode && document.getElementById("alignSeqAa")?.checked === true;
            const codonOptions = isCodonMode ? getCodonRunOptions() : null;
            const geneticCode = isCodonMode
                ? codonOptions.fallbackGeneticCode
                : (document.getElementById("alignGeneticCode")?.value || getDefaultCodonFallbackGeneticCode());

            // Guard: AA mode requires translateFasta IPC (codon mode uses JS translation)
            if (isAaMode && !window.electronAPI?.translateFasta) {
                alert("Protein translation is not available in this build. Cannot run amino acid alignment.");
                return;
            }

            // Track alignment mode for downstream steps (set early so modal can read it)
            window._codonMode = isCodonMode;
            window._alignmentMode = isCodonMode ? "codon" : isAaMode ? "aa" : "nt";

            // Get human-readable genetic code name from the select element
            const geneticCodeName = (() => {
                const sel = document.getElementById("alignGeneticCode");
                const opt = sel?.querySelector(`option[value="${geneticCode}"]`);
                return opt ? opt.textContent.trim() : `NCBI ${geneticCode}`;
            })();

            syncPythonRuntimeCard();
            syncCodonTrimToolUi();

            // Save alignment metadata globally for pipeline tracking
            window._alignmentGeneticCode = geneticCode;
            window._alignmentGeneticCodeName = geneticCodeName;
            window._codonInternalStopPolicy = codonOptions?.internalStopPolicy || null;
            window._codonAllowPseudo = codonOptions?.allowPseudo === true;

            // Build MAFFT params
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
            if (isAaMode || isCodonMode) params.push("--amino");

            // ── Codon-aware: open modal first, validate genes async with progress bar ──
            if (isCodonMode) {
                openAnalysisModal("mafft", [...state.selectedGenes], {
                    alignMode: "codon",
                    isTranslated: true,
                    geneticCode,
                    geneticCodeName,
                    mafftMode: "--amino",
                    postProcess: "trimAl -backtrans",
                });
                const logEl = document.getElementById("analysisModalLog");
                const appendLog = (msg) => { logEl.textContent += msg + "\n"; logEl.scrollTop = logEl.scrollHeight; };
                const bar = document.getElementById("analysisModalProgressBar");
                const lbl = document.getElementById("analysisModalProgressLabel");
                const pct = document.getElementById("analysisModalProgressPct");

                setTimeout(() => {
                    const validationLogs = [];
                    const logFn = (msg) => { validationLogs.push(msg); appendLog(msg); console.log("[codon]", msg); };

                    appendLog("=== CDS extraction ===");
                    const candidateContext = buildCodonCandidates(codonOptions, logFn);
                    const geneEntries = [...candidateContext.perGeneCandidates.entries()];
                    const totalGenes = geneEntries.length;
                    const cdsSelectedMarkers = candidateContext.summary.cdsSelectedMarkers || [];
                    const nonCdsMarkers = candidateContext.summary.nonCdsMarkers || [];

                    if (nonCdsMarkers.length) {
                        logFn(`[CDS extraction] Non-CDS markers skipped: ${nonCdsMarkers.join(", ")}`);
                    }

                    if (totalGenes === 0) {
                        appendLog(cdsSelectedMarkers.length === 0
                            ? i18nText("step5.modal.codon.noMarkers", "No CDS markers are currently selected for Codon-aware alignment.")
                            : i18nText("step5.modal.codon.noCandidates", "No CDS candidates remained after Codon-aware extraction."));
                        document.getElementById("analysisModalClose").classList.remove("hidden");
                        document.getElementById("analysisModalIcon").className = "fa-solid fa-xmark text-red-500 text-lg";
                        document.getElementById("analysisModalTitle").textContent = i18nText("step5.modal.codon.noSequencesTitle", "No CDS sequences found");
                        window._codonMode = false;
                        window._alignmentMode = "nt";
                        return;
                    }

                    appendLog(i18nText("step5.modal.codon.validationHeading", "=== CDS validation + translation ==="));
                    appendLog(i18nText("step5.modal.codon.validationProgress", "Validating {count} CDS gene(s)...", { count: totalGenes }));

                    const codonPairs = new Map();
                    let gIdx = 0;

                    function validateNext() {
                        if (gIdx >= totalGenes) {
                            const processedMarkers = [...codonPairs.keys()].sort();
                            const insufficientGenes = cdsSelectedMarkers.filter((geneName) => !codonPairs.has(geneName)).sort();

                            logFn(i18nText("step5.modal.codon.summaryHeading", "=== Codon-aware summary ==="));
                            logFn(i18nText("step5.modal.codon.summary.processed", "CDS markers to process: {markers}", { markers: processedMarkers.length ? processedMarkers.join(", ") : "none" }));
                            logFn(i18nText("step5.modal.codon.summary.nonCds", "Non-CDS markers skipped: {markers}", { markers: nonCdsMarkers.length ? nonCdsMarkers.join(", ") : "none" }));
                            logFn(i18nText("step5.modal.codon.summary.insufficient", "Genes with fewer than two valid CDS sequences: {markers}", { markers: insufficientGenes.length ? insufficientGenes.join(", ") : "none" }));

                            if (codonPairs.size === 0) {
                                appendLog(i18nText("step5.modal.codon.noValidLog", "No valid CDS sequences found for Codon-aware alignment."));
                                document.getElementById("analysisModalClose").classList.remove("hidden");
                                document.getElementById("analysisModalIcon").className = "fa-solid fa-xmark text-red-500 text-lg";
                                document.getElementById("analysisModalTitle").textContent = i18nText("step5.modal.codon.noValidTitle", "No valid CDS sequences");
                                window._codonMode = false;
                                window._alignmentMode = "nt";
                                return;
                            }

                            window._codonOriginalDna = {};
                            window._codonPipelineStats = {};
                            window._codonSequenceMeta = {};
                            window._codonValidationLogs = validationLogs;
                            window._codonAllowBacktranslationWarnings = codonOptions?.allowBacktranslationWarnings === true;
                            const files = {};
                            for (const [geneName, pair] of codonPairs) {
                                window._codonOriginalDna[geneName] = pair.cdsDnaFasta;
                                window._codonPipelineStats[geneName] = pair.stats || {};
                                window._codonSequenceMeta[geneName] = pair.sequenceMeta || {};
                                files[geneName] = pair.proteinFasta;
                            }
                            logFn("=== MAFFT protein alignment ===");
                            logFn(`Sending ${codonPairs.size} gene(s) to MAFFT (--amino)`);

                            if (window._phaseTracker) {
                                window._phaseTracker.complete("validation", `${totalGenes} / ${totalGenes} genes validated`);
                                window._phaseTracker.activate("mafft");
                            }
                            document.getElementById("analysisModalTitle").textContent = "Running MAFFT\u2026";

                            window.electronAPI.runAnalysis({ files, params, threads });
                            return;
                        }

                        const [geneName, seqs] = geneEntries[gIdx++];
                        const pair = processCodonGeneEntry(geneName, seqs, codonOptions, logFn);
                        if (pair?.eligible) codonPairs.set(geneName, pair);

                        const pctVal = Math.round(gIdx / totalGenes * 100);
                        if (window._phaseTracker) window._phaseTracker.update("validation", pctVal, `${gIdx} / ${totalGenes} genes validated`);

                        setTimeout(validateNext, 0);
                    }

                    setTimeout(validateNext, 0);
                }, 50);
                return;
            }

            // ── Nucleotide / Amino-acid: open modal first, then generate FASTAs ──────────
            window._codonOriginalDna = null;
            window._codonSequenceMeta = null;
            window._codonValidationLogs = null;
            window._codonAllowBacktranslationWarnings = false;

            openAnalysisModal("mafft", [...state.selectedGenes], {
                alignMode: isAaMode ? "aa" : "nt",
                isTranslated: isAaMode,
                geneticCode: isAaMode ? geneticCode : null,
                geneticCodeName: isAaMode ? geneticCodeName : null,
                mafftMode: isAaMode ? "--amino" : "nucleotide mode",
                postProcess: null,
            });

            const logEl = document.getElementById("analysisModalLog");
            const appendLog = (msg) => { logEl.textContent += msg + "\n"; logEl.scrollTop = logEl.scrollHeight; };
            appendLog("Preparing FASTA files\u2026");
            if (isAaMode) appendLog(`Using genetic code: ${geneticCodeName}`);

            setTimeout(async () => {
                const rawFastaFiles = generateFastaFiles();
                if (rawFastaFiles.size === 0) {
                    appendLog("No sequences to align. Check records and gene selection.");
                    document.getElementById("analysisModalClose").classList.remove("hidden");
                    document.getElementById("analysisModalIcon").className = "fa-solid fa-xmark text-red-500 text-lg";
                    document.getElementById("analysisModalTitle").textContent = "No sequences to align";
                    return;
                }

                const validRaw = {};
                const skippedRaw = [];
                for (const [name, content] of rawFastaFiles) {
                    const rawCount = countFastaSequences(content);
                    if (rawCount < 2) skippedRaw.push({ name, count: rawCount });
                    else validRaw[name] = { content, rawCount };
                }

                if (Object.keys(validRaw).length === 0) {
                    const detail = skippedRaw.map(s => `\u2022 ${s.name}: ${s.count} sequence(s)`).join("\n");
                    appendLog(`No markers have enough sequences (minimum 2 required):\n${detail}`);
                    document.getElementById("analysisModalClose").classList.remove("hidden");
                    document.getElementById("analysisModalIcon").className = "fa-solid fa-xmark text-red-500 text-lg";
                    document.getElementById("analysisModalTitle").textContent = "Not enough sequences";
                    return;
                }

                for (const { name, count } of skippedRaw) {
                    appendLog(`${name}: skipped, only ${count} sequence${count !== 1 ? "s" : ""}`);
                }

                if (isAaMode) {
                    appendLog(`Preparing FASTA files for ${Object.keys(validRaw).length} marker(s)\u2026`);
                    appendLog("Translating sequences with seqkit\u2026");
                    try {
                        const files = {};
                        for (const [name, { content, rawCount }] of Object.entries(validRaw)) {
                            appendLog(`  Translating ${name} (${rawCount} sequences)\u2026`);
                            const result = await window.electronAPI.translateFasta({ content, code: geneticCode });
                            const translatedContent = result || content;
                            const translatedCount = countFastaSequences(translatedContent);
                            if (translatedCount !== rawCount) {
                                throw new Error(`${name}: translation changed sequence count from ${rawCount} to ${translatedCount}`);
                            }
                            files[name] = translatedContent;
                            appendLog(`${name}: ${rawCount} raw sequences \u2192 ${translatedCount} translated \u2192 MAFFT`);
                        }
                        appendLog("Running MAFFT with --amino\u2026");
                        window.electronAPI.runAnalysis({ files, params, threads });
                    } catch (e) {
                        appendLog(`\nTranslation failed: ${e.message || e}`);
                        document.getElementById("analysisModalClose").classList.remove("hidden");
                        document.getElementById("analysisModalIcon").className = "fa-solid fa-xmark text-red-500 text-lg";
                        document.getElementById("analysisModalTitle").textContent = "Translation Failed";
                    }
                } else {
                    const files = {};
                    for (const [name, { content, rawCount }] of Object.entries(validRaw)) {
                        files[name] = content;
                        appendLog(`${name}: ${rawCount} sequences \u2192 MAFFT`);
                    }
                    window.electronAPI.runAnalysis({ files, params, threads });
                }
            }, 50);
        });
    }
});

// ========================================================================
// Electron: Analysis Modal Helpers
// ========================================================================

let _analysisMarkers = [];

// Add a completed full-size phase bar to the stack and reset the active bar.
function createPhaseTracker(phases) {
    const stack = document.getElementById("analysisProgressStack");
    if (!stack) return null;
    stack.innerHTML = "";
    const tracker = { phases: {} };
    for (const { id, label } of phases) {
        const row = document.createElement("div");
        row.className = "mb-2.5";
        row.innerHTML = `<div class="flex justify-between text-sm mb-1">`
            + `<span class="flex items-center gap-1.5 text-gray-400">`
            + `<i class="fa-regular fa-circle text-gray-300 text-xs" data-phase-icon></i>`
            + `<span data-phase-text>${label}</span></span>`
            + `<span class="text-gray-400 tabular-nums" data-phase-pct>\u2014</span></div>`
            + `<div class="w-full bg-gray-100 rounded-full h-2">`
            + `<div class="bg-gray-200 h-2 rounded-full" data-phase-bar style="width:0%"></div></div>`;
        stack.appendChild(row);
        const spanEl = row.querySelector("[data-phase-icon]").parentElement;
        tracker.phases[id] = {
            iconEl: row.querySelector("[data-phase-icon]"),
            spanEl,
            textEl: row.querySelector("[data-phase-text]"),
            pctEl: row.querySelector("[data-phase-pct]"),
            barEl: row.querySelector("[data-phase-bar]"),
        };
    }
    tracker.activate = function (id) {
        const p = this.phases[id];
        if (!p) return;
        p.iconEl.className = "fa-solid fa-spinner fa-spin text-splace-blue-500 text-xs";
        p.spanEl.className = "flex items-center gap-1.5 text-splace-blue-700";
        p.pctEl.className = "text-splace-blue-600 tabular-nums";
        p.pctEl.textContent = "0%";
        p.barEl.className = "bg-splace-blue-600 h-2 rounded-full";
        p.barEl.style.width = "0%";
    };
    tracker.update = function (id, pctVal, text) {
        const p = this.phases[id];
        if (!p) return;
        p.barEl.style.width = pctVal + "%";
        p.pctEl.textContent = pctVal + "%";
        if (text !== undefined) p.textEl.textContent = text;
    };
    tracker.complete = function (id, text) {
        const p = this.phases[id];
        if (!p) return;
        p.iconEl.className = "fa-solid fa-check text-green-500 text-xs";
        p.spanEl.className = "flex items-center gap-1.5 text-green-700 font-medium";
        p.pctEl.className = "text-green-600 font-medium tabular-nums";
        p.pctEl.textContent = "100%";
        p.barEl.className = "bg-green-500 h-2 rounded-full";
        p.barEl.style.width = "100%";
        if (text !== undefined) p.textEl.textContent = text;
    };
    return tracker;
}

function openAnalysisModal(phase, markers, options) {
    const opts = options || {};
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
    const badges = document.getElementById("analysisModalBadges");
    const sectionCopy = document.getElementById("analysisModalSectionCopy");

    icon.className = "fa-solid fa-spinner fa-spin text-splace-blue-600 text-lg";

    // Title
    if (phase === "mafft") {
        title.textContent = opts.alignMode === "codon"
            ? i18nText("step5.modal.codon.title.running", "Running Codon-aware Alignment\u2026")
            : "Running MAFFT Alignment\u2026";
    } else if (phase === "trimAl") {
        title.textContent = opts.alignMode === "codon"
            ? i18nText("step5.modal.codon.title.trim", "Running trimAl (protein-guided)\u2026")
            : "Running trimAl\u2026";
    } else {
        title.textContent = "Running\u2026";
    }

    // Subtitle
    if (phase === "mafft") {
        if (opts.alignMode === "codon") {
            sub.textContent = i18nText("step5.modal.codon.subtitle.mafft", "CDS extraction · CDS validation · MAFFT --amino · back-translation to codon DNA");
        } else if (opts.alignMode === "aa") {
            sub.textContent = "Amino Acid mode \u00b7 Translation ON \u00b7 MAFFT --amino";
        } else {
            sub.textContent = `Nucleotide mode · Translation OFF · ${markers.length} marker${markers.length !== 1 ? 's' : ''} queued`;
        }
    } else {
        if (phase === "trimAl" && opts.alignMode === "codon") {
            sub.textContent = i18nText("step5.modal.codon.subtitle.trim", "Protein alignment · trimAl -backtrans · codon-preserving DNA");
        } else {
            sub.textContent = `${markers.length} marker${markers.length !== 1 ? "s" : ""} queued`;
        }
    }

    // Section copy
    if (sectionCopy) {
        if (phase === "trimAl" && opts.alignMode === "codon") {
            sectionCopy.textContent = i18nText("step5.modal.codon.copy.trim", "trimAl receives the protein alignment from MAFFT and uses it to guide codon-block trimming. Output is back-translated to codon-preserving DNA.");
        } else if (opts.alignMode === "codon") {
            sectionCopy.textContent = i18nText("step5.modal.codon.copy.mafft", "CDS sequences are extracted, translated, aligned as proteins, then back-translated to codon-preserving DNA.");
        } else if (opts.alignMode === "aa") {
            sectionCopy.textContent = opts.geneticCodeName
                ? `Sequences are being translated before alignment using genetic code ${opts.geneticCodeName}.`
                : "Sequences are being translated before alignment.";
        } else {
            sectionCopy.textContent = "DNA/RNA sequences are aligned directly.";
        }
    }

    // Mode badges
    if (badges) {
        if (opts.alignMode) {
            const modeLabel = { codon: "CDS", aa: "Amino Acid", nt: "Nucleotide" }[opts.alignMode] || opts.alignMode;
            const modeColor = { codon: "violet", aa: "blue", nt: "green" }[opts.alignMode] || "slate";
            const tranLabel = opts.isTranslated ? (opts.alignMode === "codon" ? "ON · CDS to protein" : "ON") : "OFF";
            const tranColor = opts.isTranslated ? "amber" : "slate";
            const mafftVal = opts.mafftMode || "auto";
            const mafftColor = (opts.mafftMode === "--amino") ? "indigo" : "teal";

            const items = [
                { label: "Sequence type", value: modeLabel, color: modeColor },
                { label: "Translation", value: tranLabel, color: tranColor },
                ...(opts.isTranslated && opts.geneticCodeName
                    ? [{ label: "Genetic code", value: opts.geneticCodeName, color: "teal" }]
                    : []),
                { label: "MAFFT mode", value: mafftVal, color: mafftColor },
                ...(opts.postProcess ? [{ label: "Post-processing", value: opts.postProcess, color: "orange" }] : []),
            ];

            const colorMap = {
                violet: "bg-violet-100 text-violet-700 border-violet-200",
                blue: "bg-blue-100 text-blue-700 border-blue-200",
                green: "bg-green-100 text-green-700 border-green-200",
                amber: "bg-amber-100 text-amber-700 border-amber-200",
                slate: "bg-slate-100 text-slate-600 border-slate-200",
                indigo: "bg-indigo-100 text-indigo-700 border-indigo-200",
                teal: "bg-teal-100 text-teal-700 border-teal-200",
                orange: "bg-orange-100 text-orange-700 border-orange-200",
            };

            badges.innerHTML = items.map(b => {
                const cls = colorMap[b.color] || colorMap.slate;
                return `<span class="inline-flex items-center gap-1 px-2 py-0.5 rounded-full text-xs font-medium border ${cls}"><span class="font-semibold opacity-75">${b.label}:</span><span>${b.value}</span></span>`;
            }).join("");
            badges.classList.remove("hidden");
        } else {
            badges.classList.add("hidden");
        }
    }

    bar.style.width = "0%";
    lbl.textContent = `0 / ${markers.length} markers`;
    pct.textContent = "0%";
    log.textContent = "";
    close.classList.add("hidden");

    // Build all phase bars upfront so they all appear at once
    const phaseList = [];
    if (phase === "mafft") {
        if (opts.alignMode === "codon") {
            phaseList.push(
                { id: "validation", label: i18nText("step5.modal.phase.validation", "CDS extraction + validation") },
                { id: "mafft", label: i18nText("step5.modal.phase.mafftProtein", "MAFFT protein alignment (--amino)") },
                { id: "reading", label: i18nText("step5.modal.phase.reading", "Reading alignments") },
                { id: "backtrans", label: i18nText("step5.modal.phase.backtrans", "Back-translating to codon DNA") }
            );
        } else {
            phaseList.push(
                { id: "mafft", label: "MAFFT alignment" },
                { id: "reading", label: i18nText("step5.modal.phase.reading", "Reading alignments") }
            );
        }
    } else if (phase === "trimAl") {
        phaseList.push(
            { id: "trimal", label: opts.alignMode === "codon" ? i18nText("step5.modal.phase.trimal", "trimAl -backtrans (codon-preserving)") : "trimAl trimming" }
        );
    }
    window._phaseTracker = createPhaseTracker(phaseList);
    if (window._phaseTracker && phaseList.length > 0) {
        window._phaseTracker.activate(phaseList[0].id);
    }
    window._activeIpcPhase = phase === "trimAl" ? "trimal" : "mafft";

    list.innerHTML = markers.map(m =>
        `<div id="modal-marker-${CSS.escape(m)}" class="flex items-center gap-2 px-2 py-1 rounded text-xs text-gray-600">
            <i class="fa-solid fa-clock text-gray-300 w-3"></i>
            <span class="font-mono">${m}</span>
        </div>`
    ).join("");

    setModalEscapeEnabled("analysisModal", false);
    modal.classList.remove("hidden");
    refreshManagedModalState("analysisModal");

    // Codon mode: show CDS validation logs collected before modal opened
    if (window._codonMode && window._codonValidationLogs?.length) {
        log.textContent = "=== CDS Validation ===\n" +
            window._codonValidationLogs.join("\n") +
            "\n=== MAFFT ===\n";
        log.scrollTop = log.scrollHeight;
        refreshManagedModalState("analysisModal");
    }
}

function registerManagedModal(config) {
    const modal = document.getElementById(config.modalId);
    const scrollRegion = document.getElementById(config.scrollRegionId);
    const footer = document.getElementById(config.footerId);
    const cue = document.getElementById(config.scrollCueId);
    const closeButton = document.getElementById(config.closeButtonId);
    if (!modal || !scrollRegion || !footer || !cue || !closeButton) return;

    const refresh = () => {
        const scrollable = scrollRegion.scrollHeight > scrollRegion.clientHeight + 6;
        const atBottom = scrollRegion.scrollTop + scrollRegion.clientHeight >= scrollRegion.scrollHeight - 6;
        cue.classList.toggle("hidden", !scrollable || atBottom);
        footer.classList.toggle("has-scroll-cue", scrollable && !atBottom);
    };

    scrollRegion.addEventListener("scroll", refresh, { passive: true });
    if (typeof ResizeObserver !== "undefined") {
        const resizeObserver = new ResizeObserver(() => requestAnimationFrame(refresh));
        resizeObserver.observe(scrollRegion);
        if (scrollRegion.firstElementChild) resizeObserver.observe(scrollRegion.firstElementChild);
    }
    if (typeof MutationObserver !== "undefined") {
        const mutationObserver = new MutationObserver(() => requestAnimationFrame(refresh));
        mutationObserver.observe(scrollRegion, { childList: true, subtree: true, characterData: true });
    }

    modal._managedModal = { ...config, refresh };
    requestAnimationFrame(refresh);
}

function refreshManagedModalState(modalId) {
    const modal = document.getElementById(modalId);
    modal?._managedModal?.refresh?.();
}

function setModalEscapeEnabled(modalId, enabled) {
    const modal = document.getElementById(modalId);
    if (modal) modal.dataset.escapeEnabled = enabled ? "true" : "false";
}

function closeManagedModal(modalId) {
    const modal = document.getElementById(modalId);
    if (!modal || modal.classList.contains("hidden") || modal.dataset.escapeEnabled !== "true") {
        return false;
    }
    const closeButton = document.getElementById(modal._managedModal?.closeButtonId || "");
    if (closeButton && !closeButton.classList.contains("hidden")) {
        closeButton.click();
        return true;
    }
    modal.classList.add("hidden");
    return true;
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
        refreshManagedModalState("analysisModal");
    }

    if (data.total > 0) {
        const p = Math.round((data.done / data.total) * 100);
        bar.style.width = p + "%";
        lbl.textContent = `${data.done} / ${data.total} markers`;
        pct.textContent = p + "%";
        if (window._phaseTracker) {
            const ipcPhaseId = window._activeIpcPhase || "mafft";
            window._phaseTracker.update(ipcPhaseId, p, `${data.done} / ${data.total} markers`);
        }
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

function attachPipelineToResults(markerResults, mode) {
    const isAa = mode === "aa";
    for (const mr of markerResults) {
        mr.pipeline = {
            mode,
            rawDataType: "nucleotide",
            mafftInputType: isAa ? "protein" : "nucleotide",
            mafftMode: isAa ? "--amino" : "nucleotide",
            alignedDataType: isAa ? "protein" : "nucleotide",
            finalMatrixType: isAa ? "protein" : "nucleotide",
            geneticCode: isAa ? (window._alignmentGeneticCode || null) : null,
            geneticCodeName: isAa ? (window._alignmentGeneticCodeName || null) : null,
            stopCodonPolicy: null,
            backtranslation: { status: "not applicable", method: "none" },
            trimming: { status: "not run", inputType: null, outputType: null, strategy: null, codonPreserving: false },
            codonIntegrity: { checked: false },
        };
    }
}

function updatePipelineForTrimming(markerResults) {
    const isCodon = window._alignmentMode === "codon";
    for (const mr of markerResults) {
        // Restore pipeline from the MAFFT phase (main.js re-creates result objects as shallow copies)
        const savedPipeline = window._markerPipelines?.[mr.marker];
        if (savedPipeline) {
            mr.pipeline = { ...savedPipeline };
        } else if (!mr.pipeline) {
            // Fallback: build minimal pipeline from current mode
            const mode = window._alignmentMode || "nt";
            const isAa = mode === "aa";
            mr.pipeline = {
                mode,
                rawDataType: "nucleotide",
                mafftInputType: isAa ? "protein" : "nucleotide",
                mafftMode: isAa ? "--amino" : "nucleotide",
                alignedDataType: isAa ? "protein" : "nucleotide",
                finalMatrixType: isAa ? "protein" : "nucleotide",
                geneticCode: window._alignmentGeneticCode || null,
                geneticCodeName: window._alignmentGeneticCodeName || null,
                stopCodonPolicy: null,
                backtranslation: { status: "not applicable", method: "none" },
                codonIntegrity: { checked: false },
            };
        }
        const trimMeta = mr.trimMetadata || {};
        const trimValidation = mr.trimValidation || null;
        const inputType = mr.pipeline.mafftInputType === "protein" ? "protein" : "nucleotide";
        const outputType = isCodon ? (trimMeta.outputType || "codon-preserving DNA") : (mr.pipeline.alignedDataType || "nucleotide");
        mr.pipeline.trimming = {
            status: mr.trimmedContent ? "done" : (mr.trimError ? "failed" : "not run"),
            inputType,
            outputType,
            strategy: isCodon
                ? ((trimMeta.trimmingTool || (mr.trimMethod === "trimal-backtrans" ? "trimal" : "splace-internal")) === "trimal"
                    ? "protein-guided backtranslation"
                    : "protein-guided codon fallback")
                : "column trimming",
            method: isCodon ? (mr.trimMethod === "trimal-backtrans" ? "trimAl -backtrans" : "approximate JavaScript fallback") : "trimAl",
            codonPreserving: isCodon,
            approximate: isCodon && mr.trimMethod !== "trimal-backtrans",
            trimmingTool: trimMeta.trimmingTool || (isCodon ? (mr.trimMethod === "trimal-backtrans" ? "trimal" : "splace-internal") : "trimal"),
            trimmingMode: trimMeta.trimmingMode || (isCodon ? (mr.trimMethod === "trimal-backtrans" ? "protein-guided-backtrans" : "protein-guided-codon-block-fallback") : "standard"),
            inputAlignment: trimMeta.inputAlignment || (isCodon ? `${mr.marker}_protein_aligned.fasta` : `${mr.marker}_aligned.fasta`),
            backtransInput: trimMeta.backtransInput || (isCodon ? `${mr.marker}_cds_validated.fasta` : null),
            outputAlignment: trimMeta.outputAlignment || (isCodon ? `${mr.marker}_codon_trimmed.fasta` : `${mr.marker}_trimmed.fasta`),
            codonIntegrityValidated: trimMeta.codonIntegrityValidated === true,
            validationIssues: trimValidation?.issues || [],
        };
        mr.pipeline.trimmingMethod = mr.pipeline.trimming.method;
        if (isCodon && mr.trimMethod) {
            mr.pipeline.backtranslation = {
                status: "done",
                method: mr.trimMethod === "trimal-backtrans" ? "trimAl -backtrans" : "approximate JavaScript fallback",
                approximate: mr.trimMethod !== "trimal-backtrans",
                inputAlignment: mr.pipeline.trimming.inputAlignment,
                backtransInput: mr.pipeline.trimming.backtransInput,
                outputAlignment: mr.pipeline.trimming.outputAlignment,
            };
            mr.pipeline.backtranslationMethod = mr.pipeline.backtranslation.method;
        }
        if (isCodon && (mr.trimmedContent || trimValidation)) {
            const trimSeqs = mr.trimmedContent ? parseFasta(mr.trimmedContent) : [];
            const integrity = trimSeqs.length
                ? getCodonIntegrityChecks(trimSeqs)
                : {
                    equalLengths: trimValidation?.equalLengths ?? false,
                    allLengthsMultipleOf3: trimValidation?.allLengthsMultipleOf3 ?? false,
                    gapBlocksAreTriplets: trimValidation?.gapBlocksAreTriplets ?? false,
                };
            mr.pipeline.codonIntegrity = {
                equalLengths: trimValidation?.equalLengths ?? integrity.equalLengths,
                allLengthsMultipleOf3: trimValidation?.allLengthsMultipleOf3 ?? integrity.allLengthsMultipleOf3,
                gapBlocksAreTriplets: trimValidation?.gapBlocksAreTriplets ?? integrity.gapBlocksAreTriplets,
                headersMatch: trimValidation?.headersMatch,
                dnaCompatible: trimValidation?.dnaCompatible,
                frameDisrupted: trimValidation?.frameDisrupted,
                codonIntegrityValidated: trimMeta.codonIntegrityValidated === true,
                checked: true,
            };
        }
    }
}

function finalizeAnalysisModal(result) {
    const icon = document.getElementById("analysisModalIcon");
    const title = document.getElementById("analysisModalTitle");
    const close = document.getElementById("analysisModalClose");
    const bar = document.getElementById("analysisModalProgressBar");
    const lbl = document.getElementById("analysisModalProgressLabel");
    const pct = document.getElementById("analysisModalProgressPct");
    const logEl = document.getElementById("analysisModalLog");
    const log = (msg) => { if (logEl) { logEl.textContent += msg + "\n"; logEl.scrollTop = logEl.scrollHeight; } console.log("[codon]", msg); };

    bar.style.width = "100%";
    lbl.textContent = `${result.aligned} / ${result.total} markers`;
    pct.textContent = "100%";

    // Async per-marker back-translation with progress bar phase.
    function runBacktranslation(markers, onDone) {
        icon.className = "fa-solid fa-spinner fa-spin text-splace-blue-500 text-lg";
        title.textContent = i18nText("step5.modal.codon.backtransTitle", "Back-translating\u2026");
        log("=== Back-translation ===");
        bar.style.width = "0%";
        lbl.textContent = `0 / ${markers.length} markers`;
        pct.textContent = "0%";
        let bIdx = 0;
        const fatalReports = [];

        function next() {
            if (bIdx >= markers.length) {
                if (fatalReports.length) {
                    icon.className = "fa-solid fa-triangle-exclamation text-red-600 text-lg";
                    title.textContent = "Codon-aware validation failed";
                    close.classList.remove("hidden");
                    setModalEscapeEnabled("analysisModal", true);
                    refreshManagedModalState("analysisModal");
                    log("=== Codon-aware validation failed ===");
                    fatalReports.forEach((report) => {
                        report.issues.forEach((issue) => log(`[${report.marker}] ERROR: ${issue}`));
                    });
                    close.onclick = () => { document.getElementById("analysisModal").classList.add("hidden"); };
                    return;
                }
                if (window._phaseTracker) {
                    window._phaseTracker.complete("backtrans", `${markers.length} / ${markers.length} markers back-translated`);
                }
                onDone();
                return;
            }

            const mr = markers[bIdx++];
            const backtranslationResult = applyBacktranslationToResults([mr], {
                log,
                allowCompatibilityWarnings: window._codonAllowBacktranslationWarnings === true,
            });
            const report = backtranslationResult.reports[0];
            if (report && !report.success) {
                fatalReports.push(report);
            }

            const pctVal = Math.round((bIdx / markers.length) * 100);
            bar.style.width = pctVal + "%";
            lbl.textContent = `${bIdx} / ${markers.length} markers`;
            pct.textContent = pctVal + "%";
            if (window._phaseTracker) window._phaseTracker.update("backtrans", pctVal, `${bIdx} / ${markers.length} markers`);
            setTimeout(next, 0);
        }

        setTimeout(next, 0);
    }

    async function syncCodonOutputs(markerResults) {
        if (!window._codonMode || !window.electronAPI?.saveCodonOutputs || !window._alignmentOutputDir) {
            return markerResults;
        }

        title.textContent = i18nText("step5.modal.codon.saveOutputsTitle", "Saving codon-aware outputs\u2026");
        log(i18nText("step5.modal.codon.saveOutputsStart", "Syncing codon-aware DNA/protein files to disk\u2026"));

        const syncResult = await window.electronAPI.saveCodonOutputs({
            outputDir: window._alignmentOutputDir,
            markerResults,
        });

        if (!syncResult?.success) {
            log(i18nText("step5.modal.codon.saveOutputsError", "ERROR: could not save codon-aware outputs to disk ({error})", { error: syncResult?.error || "unknown error" }));
            return markerResults;
        }

        log(i18nText("step5.modal.codon.saveOutputsSuccess", "Saved codon-aware outputs to {path}", { path: syncResult.outputDir }));
        return syncResult.markerResults || markerResults;
    }

    if (result.phase === "mafft") {
        if (result.markerResults) {
            window._alignmentOutputDir = result.outputDir;
            const markers = result.markerResults;

            // Complete MAFFT phase, activate reading phase
            if (window._phaseTracker) {
                window._phaseTracker.complete("mafft", `${result.aligned} / ${result.total} markers aligned`);
                window._phaseTracker.activate("reading");
            }
            icon.className = "fa-solid fa-spinner fa-spin text-splace-blue-500 text-lg";
            title.textContent = "Reading alignments\u2026";
            lbl.textContent = `0 / ${markers.length} markers`;

            let mIdx = 0;
            function readNext() {
                if (mIdx >= markers.length) {
                    if (window._phaseTracker) {
                        window._phaseTracker.complete("reading", `${markers.length} / ${markers.length} alignments read`);
                    }
                    if (!window._codonMode) {
                        attachPipelineToResults(markers, window._alignmentMode || "nt");
                        window._markerPipelines = {};
                        for (const mr of markers) {
                            if (mr.pipeline) window._markerPipelines[mr.marker] = { ...mr.pipeline };
                        }
                        finishMafft(markers);
                    } else {
                        if (window._phaseTracker) window._phaseTracker.activate("backtrans");
                        runBacktranslation(markers, async () => {
                            const syncedMarkers = await syncCodonOutputs(markers);
                            window._markerPipelines = {};
                            for (const mr of syncedMarkers) {
                                if (mr.pipeline) window._markerPipelines[mr.marker] = { ...mr.pipeline };
                            }
                            finishMafft(syncedMarkers);
                        });
                    }
                    return;
                }
                mIdx++;
                const pctVal = Math.round(mIdx / markers.length * 100);
                bar.style.width = pctVal + "%";
                lbl.textContent = `${mIdx} / ${markers.length} markers`;
                pct.textContent = pctVal + "%";
                if (window._phaseTracker) window._phaseTracker.update("reading", pctVal, `${mIdx} / ${markers.length} markers`);
                setTimeout(readNext, 0);
            }

            function finishMafft(finalMarkers) {
                renderAlignmentResults(finalMarkers, false);
                renderConcatSection(finalMarkers, false);
                const trimalCard = document.getElementById("trimalPromptCard");
                trimalCard.classList.remove("hidden");
                if (window._codonMode) {
                    const body = trimalCard.querySelector("[data-i18n='step5.trimal.prompt.body']");
                    if (body) body.textContent = i18nText("step5.trimal.prompt.codon.body", "In Codon-aware mode, trimAl runs on the protein alignment, then back-translates to codon-preserving DNA. Gaps are always codon-aligned.");
                }
                setBadge("step6c-badge", "done");
                updateStepBadges();
                bar.style.width = "100%";
                lbl.textContent = `${result.aligned} / ${result.total} markers`;
                pct.textContent = "100%";
                icon.className = "fa-solid fa-check text-green-600 text-lg";
                title.textContent = `Alignment complete \u2014 ${result.aligned}/${result.total} markers aligned`;
                close.classList.remove("hidden");
                setModalEscapeEnabled("analysisModal", true);
                refreshManagedModalState("analysisModal");
                close.onclick = () => {
                    document.getElementById("analysisModal").classList.add("hidden");
                    document.getElementById("trimalPromptCard").scrollIntoView({ behavior: "smooth", block: "start" });
                };
            }

            setTimeout(readNext, 50);
        } else {
            icon.className = "fa-solid fa-check text-green-600 text-lg";
            title.textContent = `Alignment complete \u2014 ${result.aligned}/${result.total} markers aligned`;
            close.classList.remove("hidden");
            setModalEscapeEnabled("analysisModal", true);
            refreshManagedModalState("analysisModal");
            close.onclick = () => { document.getElementById("analysisModal").classList.add("hidden"); };
        }
    } else if (result.phase === "trimal") {
        if (result.markerResults) {
            setModalEscapeEnabled("analysisModal", false);
            const markers = result.markerResults;
            refreshManagedModalState("analysisModal");
            icon.className = "fa-solid fa-spinner fa-spin text-splace-blue-500 text-lg";
            // Complete trimAl phase
            if (window._phaseTracker) window._phaseTracker.complete("trimal", `${result.aligned} / ${result.total} markers trimmed`);

            if ((window._alignmentMode === "codon" || window._codonMode) && (result.success !== true || (result.failedMarkers?.length || 0) > 0)) {
                icon.className = "fa-solid fa-triangle-exclamation text-red-600 text-lg";
                title.textContent = "Codon-aware trimming failed";
                close.classList.remove("hidden");
                setModalEscapeEnabled("analysisModal", true);
                refreshManagedModalState("analysisModal");
                if (result.failedMarkers?.length) {
                    log("=== Codon-aware trim validation failed ===");
                    result.failedMarkers.forEach((failure) => {
                        failure.issues.forEach((issue) => log(`[${failure.marker}] ERROR: ${issue}`));
                    });
                }
                close.onclick = () => { document.getElementById("analysisModal").classList.add("hidden"); };
                return;
            }

            function finishTrimal() {
                updatePipelineForTrimming(markers);
                renderAlignmentResults(markers, true);
                refreshManagedModalState("analysisModal");
                renderConcatSection(markers, true);
                document.getElementById("alignmentResults").classList.remove("hidden");
                setBadge("step6d-badge", "done");
                updateStepBadges();
                bar.style.width = "100%";
                lbl.textContent = `${result.aligned} / ${result.total} markers`;
                pct.textContent = "100%";
                icon.className = "fa-solid fa-check text-green-600 text-lg";
                title.textContent = `Trimming complete \u2014 ${result.aligned}/${result.total} markers trimmed`;
                close.classList.remove("hidden");
                setModalEscapeEnabled("analysisModal", true);
                refreshManagedModalState("analysisModal");
                close.onclick = () => {
                    document.getElementById("analysisModal").classList.add("hidden");
                    document.getElementById("alignmentResults").scrollIntoView({ behavior: "smooth", block: "start" });
                };
            }

            // trimAl -backtrans already outputs codon-preserving DNA — no back-translation needed
            title.textContent = "Processing results\u2026";
            setTimeout(finishTrimal, 50);
        } else {
            icon.className = "fa-solid fa-check text-green-600 text-lg";
            title.textContent = `Trimming complete \u2014 ${result.aligned}/${result.total} markers trimmed`;
            close.classList.remove("hidden");
            setModalEscapeEnabled("analysisModal", true);
            refreshManagedModalState("analysisModal");
            close.onclick = () => { document.getElementById("analysisModal").classList.add("hidden"); };
        }
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

// Build a { header -> seq } lookup from a FASTA string.
function parseFastaToMap(fastaText) {
    const map = {};
    for (const s of parseFasta(fastaText)) map[s.name] = s.seq;
    return map;
}

// Count FASTA sequences by counting '>' header lines.
function countFastaSequences(fasta) {
    return (fasta.match(/^>/gm) || []).length;
}

function getCodonIntegrityChecks(seqs) {
    const lengths = seqs.map((seq) => seq.seq.length);
    const firstLength = lengths[0] || 0;
    return {
        equalLengths: lengths.every((length) => length === firstLength),
        allLengthsMultipleOf3: lengths.every((length) => length % 3 === 0),
        gapBlocksAreTriplets: seqs.every((seq) =>
            (seq.seq.match(/-+/g) || []).every((run) => run.length % 3 === 0)
        ),
        length: firstLength,
    };
}

// Back-translate one aligned protein sequence to DNA using the original
// (unaligned) CDS codons.  Each '-' in the protein → '---' in the DNA.
function translateCodonWithGeneticCode(codon, geneticCode = "1") {
    const codeKey = String(geneticCode || "1");
    const table = GENETIC_CODE_TABLES[codeKey] || CODON_TABLE;
    return table[codon] || "X";
}

function backTranslateProteinAlignmentDetailed(alignedProteinSeq, originalDna, geneticCode = "1") {
    const dna = originalDna.toUpperCase().replace(/U/g, "T").replace(/[^ACGTNRYSWKMBDHV]/g, "");
    const codons = [];
    for (let i = 0; i + 2 < dna.length; i += 3) codons.push(dna.slice(i, i + 3));
    let idx = 0, output = "";
    let mismatchCount = 0;
    let missingCodonCount = 0;
    const issues = [];

    for (let pos = 0; pos < alignedProteinSeq.length; pos++) {
        const aa = String(alignedProteinSeq[pos] || "").toUpperCase();
        if (aa === "-" || aa === ".") output += "---";
        else {
            const codon = codons[idx++];
            if (!codon) {
                output += "NNN";
                missingCodonCount++;
                issues.push(`AA position ${pos + 1}: no CDS codon available for ${aa}`);
                continue;
            }
            const translatedAa = translateCodonWithGeneticCode(codon, geneticCode);
            if (aa !== "X" && aa !== "?" && aa !== translatedAa) {
                mismatchCount++;
                issues.push(`AA position ${pos + 1}: protein ${aa} but ${codon} translates to ${translatedAa} with ${getGeneticCodeName(geneticCode)}`);
            }
            output += codon;
        }
    }

    const consumedCodons = Math.min(idx, codons.length);
    const leftoverCodons = Math.max(0, codons.length - consumedCodons);
    if (leftoverCodons > 0) {
        issues.push(`CDS retained ${leftoverCodons} unused codon(s) after back-translation`);
    }

    return {
        dna: output,
        mismatchCount,
        missingCodonCount,
        consumedCodons,
        totalCodons: codons.length,
        leftoverCodons,
        frameDisrupted: missingCodonCount > 0 || leftoverCodons > 0,
        issues,
        geneticCode: String(geneticCode || "1"),
    };
}

function backTranslateProteinAlignment(alignedProteinSeq, originalDna) {
    return backTranslateProteinAlignmentDetailed(alignedProteinSeq, originalDna).dna;
}

function summarizeCodonIntegrityIssues(integrity, frameDisrupted = false) {
    const issues = [];
    if (!integrity.equalLengths) issues.push("sequences have unequal lengths");
    if (!integrity.allLengthsMultipleOf3) issues.push(`alignment length ${integrity.length} is not divisible by 3`);
    if (!integrity.gapBlocksAreTriplets) issues.push("gap runs are not multiples of 3");
    if (frameDisrupted) issues.push("codon frame disruption detected during back-translation");
    return issues;
}

// After MAFFT returns protein alignments, back-translate every marker's
// alignedContent to DNA in-place using the stored original CDS sequences.
function applyBacktranslationToResults(markerResults, options = {}) {
    const logEl = options.logEl || document.getElementById("analysisModalLog");
    const log = options.log || ((msg) => {
        if (logEl) { logEl.textContent += msg + "\n"; logEl.scrollTop = logEl.scrollHeight; }
        console.log("[codon]", msg);
    });
    const allowCompatibilityWarnings = options.allowCompatibilityWarnings === true;
    const reports = [];

    for (const mr of markerResults) {
        const dnaCdsContent = options.rawFastas?.[mr.marker] || window._codonOriginalDna?.[mr.marker];
        const sequenceMeta = options.sequenceMeta?.[mr.marker] || window._codonSequenceMeta?.[mr.marker] || {};
        if (!dnaCdsContent || !mr.alignedContent) {
            const issues = ["missing validated CDS DNA or protein alignment"];
            log(`[${mr.marker}] ERROR: ${issues[0]}`);
            reports.push({ marker: mr.marker, success: false, fatal: true, issues, warnings: [] });
            continue;
        }

        const dnaMap = parseFastaToMap(dnaCdsContent);
        const protSeqs = parseFasta(mr.alignedContent);
        let dnaFasta = "", warnings = 0;
        let mismatchCount = 0;
        let frameDisrupted = false;
        const compatibilityIssues = [];
        const frameIssues = [];

        for (const ps of protSeqs) {
            const dna = dnaMap[ps.name];
            if (!dna) {
                log(`[${mr.marker}] WARNING: no DNA for header "${ps.name}"`);
                warnings++;
                frameIssues.push(`${ps.name}: missing validated CDS DNA`);
                continue;
            }

            const geneticCode = String(
                sequenceMeta[ps.name]?.geneticCode ||
                window._codonPipelineStats?.[mr.marker]?.geneticCode ||
                window._alignmentGeneticCode ||
                getDefaultCodonFallbackGeneticCode()
            );
            const detailed = backTranslateProteinAlignmentDetailed(ps.seq, dna.replace(/-/g, ""), geneticCode);
            mismatchCount += detailed.mismatchCount;
            if (detailed.frameDisrupted) {
                frameDisrupted = true;
                frameIssues.push(`${ps.name}: consumed ${detailed.consumedCodons}/${detailed.totalCodons} codons`);
            }
            if (detailed.issues.length) {
                detailed.issues.forEach((issue) => {
                    const formatted = `${ps.name}: ${issue}`;
                    if (issue.includes("protein ")) compatibilityIssues.push(formatted);
                    else frameIssues.push(formatted);
                });
            }
            dnaFasta += `>${ps.name}\n${detailed.dna}\n`;
        }
        if (!dnaFasta) {
            const issues = ["back-translation produced no sequences"];
            log(`[${mr.marker}] ERROR: ${issues[0]}`);
            reports.push({ marker: mr.marker, success: false, fatal: true, issues, warnings: [] });
            continue;
        }

        const seqs = parseFasta(dnaFasta);
        const dnaLen = seqs[0]?.seq.length || 0;
        const gapCount = seqs.reduce((a, s) => a + (s.seq.match(/-/g) || []).length, 0);
        const totalChars = seqs.reduce((a, s) => a + s.seq.length, 0);
        const integrity = getCodonIntegrityChecks(seqs);
        const integrityIssues = summarizeCodonIntegrityIssues(integrity, frameDisrupted);
        const fatalIssues = [
            ...frameIssues,
            ...integrityIssues,
            ...(allowCompatibilityWarnings ? [] : compatibilityIssues),
        ];
        const warningIssues = [
            ...(allowCompatibilityWarnings ? compatibilityIssues : []),
            ...(warnings > 0 ? [`${warnings} sequences skipped because DNA headers were missing`] : []),
        ];
        const success = fatalIssues.length === 0;
        const pipelineStats = window._codonPipelineStats?.[mr.marker] || {};

        mr.pipeline = {
            ...(mr.pipeline || {}),
            mode: "codon",
            rawDataType: "cds-dna",
            mafftInputType: "protein",
            mafftMode: "--amino",
            alignedDataType: "codon-dna",
            finalMatrixType: "codon-preserving DNA",
            geneticCode: pipelineStats.geneticCode || window._alignmentGeneticCode || getDefaultCodonFallbackGeneticCode(),
            geneticCodeName: pipelineStats.geneticCodeName || window._alignmentGeneticCodeName || getGeneticCodeName(window._alignmentGeneticCode || getDefaultCodonFallbackGeneticCode()),
            stopCodonPolicy: {
                terminalStopsRemoved: pipelineStats.terminalStopsRemoved || 0,
                incompleteStopFragmentsRemoved: pipelineStats.incompleteStopFragmentsRemoved || 0,
                internalStopsDetected: pipelineStats.internalStopsDetected || 0,
                internalStopPolicy: pipelineStats.internalStopPolicy || window._codonInternalStopPolicy || "exclude",
            },
            backtranslation: {
                status: success || allowCompatibilityWarnings ? "done" : "failed",
                method: "JavaScript codon-aware back-translation",
                approximate: true,
                compatibilityWarnings: compatibilityIssues.length,
                compatibilityPolicy: allowCompatibilityWarnings ? "continue-with-warning" : "block",
                frameDisrupted,
            },
            trimming: mr.pipeline?.trimming || {
                status: "not run",
                inputType: null,
                outputType: null,
                strategy: null,
                method: null,
                codonPreserving: null,
                approximate: false,
            },
            codonIntegrity: {
                equalLengths: integrity.equalLengths,
                allLengthsMultipleOf3: integrity.allLengthsMultipleOf3,
                gapBlocksAreTriplets: integrity.gapBlocksAreTriplets,
                frameDisrupted,
                checked: true,
            },
        };

        if (success || allowCompatibilityWarnings) {
            mr.proteinAlignedContent = mr.alignedContent;
            mr.alignedContent = dnaFasta;
            mr.rawContent = dnaCdsContent;
            mr.alignedStats = {
                numSeqs: seqs.length,
                length: dnaLen,
                avgLen: dnaLen,
                gapPct: totalChars > 0 ? ((gapCount / totalChars) * 100).toFixed(1) : "0.0",
            };
        }

        if (compatibilityIssues.length) {
            compatibilityIssues.forEach((issue) => {
                log(`[${mr.marker}] ${allowCompatibilityWarnings ? "WARN" : "ERROR"}: ${issue}`);
            });
        }
        if (frameIssues.length) {
            frameIssues.forEach((issue) => log(`[${mr.marker}] ERROR: ${issue}`));
        }
        if (integrityIssues.length) {
            integrityIssues.forEach((issue) => log(`[${mr.marker}] ERROR: ${issue}`));
        }
        if (warningIssues.length && allowCompatibilityWarnings) {
            log(`[${mr.marker}] WARN: continuing with ${warningIssues.length} back-translation warning(s) by user request`);
        }
        if (success || allowCompatibilityWarnings) {
            log(`[${mr.marker}] ✓ Back-translation: ${seqs.length} seq × ${dnaLen} bp${mismatchCount > 0 ? ` | ${mismatchCount} compatibility warning(s)` : ""}`);
        }

        reports.push({
            marker: mr.marker,
            success,
            fatal: !success,
            issues: fatalIssues,
            warnings: warningIssues,
            mismatchCount,
            integrity,
            frameDisrupted,
        });
    }

    return {
        markerResults,
        reports,
        success: reports.every((report) => report.success),
        errors: reports.filter((report) => !report.success),
    };
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
    "3": { ...CODON_TABLE, ATA: "M", TGA: "W", CTT: "T", CTC: "T", CTA: "T", CTG: "T" },
    "4": { ...CODON_TABLE, TGA: "W" },
    "5": { ...CODON_TABLE, ATA: "M", TGA: "W", AGA: "S", AGG: "S" },
    "6": { ...CODON_TABLE, TAA: "Q", TAG: "Q" },
    "9": { ...CODON_TABLE, AAA: "N", AGA: "S", TGA: "W" },
    "11": CODON_TABLE,
    "13": { ...CODON_TABLE, AGA: "G", AGG: "G", ATA: "M", TGA: "W" },
    "14": { ...CODON_TABLE, AAA: "N", AGA: "S", TAA: "Y", TGA: "W" },
    "16": { ...CODON_TABLE, TAG: "L" },
    "21": { ...CODON_TABLE, ATA: "M", TGA: "W", AGA: "S", AGG: "S", AAA: "N" },
    "24": { ...CODON_TABLE, AGA: "S", AGG: "K", TGA: "W" },
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

function renderAlignmentPreview(container, sequences, mode, rawSeqs, geneName) {
    if (!container) return;
    const rows = sequences || [];
    if (!rows.length) {
        container.innerHTML = `<div class="alignment-preview"><div class="alignment-preview-empty">${i18nText("step5.results.preview.empty", "No sequences available for this view.")}</div></div>`;
        return;
    }

    const totalStops = mode === "aa" ? rows.reduce((acc, seq) => acc + (seq.stopCount || 0), 0) : 0;
    const summaryKey = mode === "aa" ? "step5.results.preview.summaryStops" : "step5.results.preview.summary";
    const summaryFallback = mode === "aa"
        ? "{seqs} sequences rendered; {stops} possible stop codon(s)"
        : "{seqs} sequences rendered";
    const summary = i18nText(summaryKey, summaryFallback, { seqs: rows.length, stops: totalStops });

    // Helper: extract start or stop codon (first/last 3 non-gap bases) from a raw DNA sequence
    function getDnaCodon(rawSeq, fromEnd) {
        if (!rawSeq || !rawSeq.seq) return { codon: "---", incomplete: false };
        const dna = rawSeq.seq.replace(/-/g, "").toUpperCase();
        if (!dna) return { codon: "---", incomplete: false };
        const raw = fromEnd ? dna.slice(-3) : dna.slice(0, 3);
        const incomplete = raw.length < 3;
        const codon = fromEnd ? raw.padStart(3, "-") : raw.padEnd(3, "-");
        return { codon, incomplete };
    }

    // Helper: render a 3-base codon as colored spans using NT_COLORS
    function renderCodon(codon, incomplete) {
        const title = incomplete ? ' title="Incomplete codon"' : "";
        return Array.from(codon).map(ch => {
            const col = incomplete && ch === "-" ? "#fef3c7" : (NT_COLORS[ch] || "#e2e8f0");
            const fg = (ch === "-") ? (incomplete ? "#d97706" : "#94a3b8") : "rgba(0,0,0,0.72)";
            return `<span class="alignment-residue" style="background:${col};color:${fg};"${title}>${ch}</span>`;
        }).join("");
    }

    const codonStats = window._codonPipelineStats?.[geneName];

    const previewRows = rows.slice(0, 12).map((seq, i) => {
        const length = seq.seq.length;
        const stopDisplay = mode === "aa"
            ? String(seq.stopCount || 0)
            : i18nText("step5.results.preview.none", "none");
        const rawSeq = rawSeqs ? rawSeqs[i] : seq;
        const { codon: startCodon, incomplete: startIncomplete } = getDnaCodon(rawSeq, false);
        const { codon: stopCodon, incomplete: stopIncomplete } = getDnaCodon(rawSeq, true);

        const stopCodonHtml = codonStats
            ? `<span style="font-size:0.68rem;padding:1px 5px;background:#f1f5f9;border-radius:3px;color:#94a3b8;white-space:nowrap;" title="Stop codon removed during preprocessing">removed</span>`
            : renderCodon(stopCodon, stopIncomplete);

        return `
            <div class="alignment-preview-row">
                <div class="alignment-preview-name" title="${seq.name}">${seq.name}</div>
                <div class="alignment-preview-metric">${length}</div>
                <div class="alignment-preview-metric">${stopDisplay}</div>
                <div class="alignment-preview-codon">${renderCodon(startCodon, startIncomplete)}</div>
                <div class="alignment-preview-codon">${stopCodonHtml}</div>
            </div>`;
    }).join("");

    const alignedNames = new Set(rows.map(s => s.name));
    const missing = [];
    if (geneName && state.records && state.records.length > 0) {
        for (const record of state.records) {
            const tax = extractTaxonomy(record);
            const header = buildHeader(record, tax);
            if (!alignedNames.has(header)) {
                missing.push({
                    name: header,
                    organism: record.editedOrganism || record.organism || record.accession || ""
                });
            }
        }
    }

    const missingHtml = missing.length > 0 ? `
        <div class="alignment-preview-missing">
            <div class="alignment-preview-missing-title">Missing (${missing.length})</div>
            <div class="alignment-preview-missing-list">
                ${missing.map(m => `<div class="alignment-preview-missing-row" title="${m.name}">${m.organism || m.name}</div>`).join("")}
            </div>
        </div>` : "";

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
                    <div>Start Codon</div>
                    <div>Stop Codon</div>
                </div>
                ${previewRows}
            </div>
            ${missingHtml}
        </div>`;
}

function getAlignmentViewOptions(mode) {
    if (mode === "aa") {
        return {
            colors: AA_COLORS,
            isProtein: true,
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

async function renderAlignmentTab(markerData, tab, mode) {
    const contentMap = {
        raw: markerData.rawContent,
        aligned: markerData.alignedContent,
        trimmed: markerData.trimmedContent,
    };
    // For codon mode "protein-guide" view: show protein alignment in the aligned tab
    let content;
    if (mode === "protein-guide" && tab === "aligned" && markerData.proteinAlignedContent) {
        content = markerData.proteinAlignedContent;
    } else {
        content = contentMap[tab];
    }
    const wrap = document.getElementById(`vizwrap-${markerData.id}-${tab}`);
    const preview = document.getElementById(`preview-${markerData.id}-${tab}`);
    if (!wrap || !preview || !content) return;
    const sequences = parseFasta(content);
    const geneticCode = document.getElementById(`genetic-code-${markerData.id}-${tab}`)?.value || "2";

    let displaySequences;
    if (mode === "protein-guide") {
        // Content is already a protein alignment — render with AA colors directly
        renderAlignmentPlotly(wrap, sequences, getAlignmentViewOptions("aa"));
        renderAlignmentPreview(preview, sequences, "aa", sequences, markerData.id);
        wrap.dataset.rendered = "1";
        wrap.dataset.view = mode;
        return;
    }
    if (mode === "aa") {
        if (window.isElectron && window.electronAPI?.translateFasta) {
            // Strip alignment gaps; seqkit translate handles the rest
            const strippedFasta = sequences
                .filter(s => s.seq.replace(/-/g, "").length > 0)
                .map(s => `>${s.name}\n${s.seq.replace(/-/g, "")}`).
                join("\n");
            try {
                const proteinFasta = await window.electronAPI.translateFasta({ content: strippedFasta, code: geneticCode });
                const proteinSeqs = parseFasta(proteinFasta || "");
                displaySequences = proteinSeqs.map(s => ({
                    name: s.name,
                    seq: s.seq,
                    stopCount: (s.seq.match(/\*/g) || []).length,
                }));
            } catch (e) {
                displaySequences = translateAlignmentSequences(sequences, geneticCode);
            }
        } else {
            displaySequences = translateAlignmentSequences(sequences, geneticCode);
        }
    } else {
        displaySequences = sequences;
    }

    renderAlignmentPlotly(wrap, displaySequences, getAlignmentViewOptions(mode));
    renderAlignmentPreview(preview, displaySequences, mode, sequences, markerData.id);
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
    if (!ctx) return;

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
    const BASE_COMP_H = 44;
    let nameW = 190;
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
    let cellW = 20;

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
        const avail = seqOuter.clientWidth || container.clientWidth - nameW;
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
    rulerSpacer.style.cssText = `width:${nameW}px;flex-shrink:0;height:${RULER_H}px;background:#f3f4f6;border-right:2px solid #d1d5db;`;
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
    bodyRow.style.cssText = "display:flex;flex:1;overflow-y:auto;min-height:0;align-items:flex-start;";
    container.appendChild(bodyRow);

    const nameDiv = document.createElement("div");
    nameDiv.style.cssText = `width:${nameW}px;flex-shrink:0;overflow-x:hidden;overflow-y:visible;background:#f3f4f6;border-right:2px solid #d1d5db;`;
    bodyRow.appendChild(nameDiv);

    const nameCanvas = document.createElement("canvas");
    nameCanvas.width = nameW;
    nameCanvas.style.cssText = "display:block;";
    nameDiv.appendChild(nameCanvas);

    // ── resize handle ──────────────────────────────────────────────────────
    const resizeHandle = document.createElement("div");
    resizeHandle.style.cssText = "width:5px;flex-shrink:0;background:transparent;cursor:col-resize;position:relative;z-index:2;";
    resizeHandle.title = "Drag to resize sequence name column";
    bodyRow.appendChild(resizeHandle);

    const seqOuter = document.createElement("div");
    seqOuter.style.cssText = "overflow-x:auto;overflow-y:visible;flex:1;cursor:grab;";
    bodyRow.appendChild(seqOuter);

    const mainCanvas = document.createElement("canvas");
    mainCanvas.style.cssText = "display:block;";
    seqOuter.appendChild(mainCanvas);

    // ── consensus row ──────────────────────────────────────────────────────
    const consSeqRow = document.createElement("div");
    consSeqRow.style.cssText = "display:flex;flex-shrink:0;border-top:2px solid #cbd5e1;";
    container.appendChild(consSeqRow);

    const consSeqNameDiv = document.createElement("div");
    consSeqNameDiv.style.cssText = `width:${nameW}px;flex-shrink:0;height:${CONS_SEQ_H}px;background:#dbeafe;border-right:2px solid #d1d5db;display:flex;align-items:center;justify-content:flex-end;padding-right:8px;`;
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
    consBarNameDiv.style.cssText = `width:${nameW}px;flex-shrink:0;height:${CONS_BAR_H}px;background:#f3f4f6;border-right:2px solid #d1d5db;display:flex;align-items:center;justify-content:flex-end;padding-right:8px;`;
    consBarNameDiv.innerHTML = `<span style="font-size:10px;color:#6b7280;font-weight:600;letter-spacing:.04em;">${labels.conservation}</span>`;
    consBarRow.appendChild(consBarNameDiv);

    const consBarOuter = document.createElement("div");
    consBarOuter.style.cssText = "overflow:hidden;flex:1;";
    consBarRow.appendChild(consBarOuter);

    const consBarCanvas = document.createElement("canvas");
    consBarCanvas.height = CONS_BAR_H;
    consBarCanvas.style.cssText = "display:block;";
    consBarOuter.appendChild(consBarCanvas);

    // ── base composition row ───────────────────────────────────────────────
    const baseCompRow = document.createElement("div");
    baseCompRow.style.cssText = "display:flex;flex-shrink:0;border-top:1px solid #e5e7eb;";
    container.appendChild(baseCompRow);

    const baseCompNameDiv = document.createElement("div");
    baseCompNameDiv.style.cssText = `width:${nameW}px;flex-shrink:0;height:${BASE_COMP_H}px;background:#f3f4f6;border-right:2px solid #d1d5db;display:flex;align-items:center;justify-content:flex-end;padding-right:8px;`;
    baseCompNameDiv.innerHTML = `<span style="font-size:10px;color:#6b7280;font-weight:600;letter-spacing:.04em;">COMPOSITION</span>`;
    baseCompRow.appendChild(baseCompNameDiv);

    const baseCompOuter = document.createElement("div");
    baseCompOuter.style.cssText = "overflow:hidden;flex:1;";
    baseCompRow.appendChild(baseCompOuter);

    const baseCompCanvas = document.createElement("canvas");
    baseCompCanvas.height = BASE_COMP_H;
    baseCompCanvas.style.cssText = "display:block;";
    baseCompOuter.appendChild(baseCompCanvas);

    // ── helper: sync all name-column widths after resize ──────────────────
    function updateNameWidth() {
        nameDiv.style.width = nameW + "px";
        nameCanvas.width = nameW;
        rulerSpacer.style.width = nameW + "px";
        consSeqNameDiv.style.width = nameW + "px";
        consBarNameDiv.style.width = nameW + "px";
        baseCompNameDiv.style.width = nameW + "px";
    }

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
        rCtx.textAlign = "center";
        // tick interval: minimum every 2 positions, then nice steps
        const minTickPx = 35;
        const rawStep = Math.ceil(minTickPx / cellW);
        const niceSteps = [2, 5, 10, 25, 50, 100, 200, 500, 1000, 2000, 5000];
        const tickStep = niceSteps.find(s => s >= rawStep) || rawStep;
        for (let p = 0; p < seqLen; p++) {
            if (p === 0 || (p + 1) % tickStep === 0) {
                const x = p * cellW + cellW / 2;
                rCtx.fillStyle = "#9ca3af";
                rCtx.fillRect(p * cellW + cellW / 2 - 0.5, RULER_H - 6, 1, 6);
                rCtx.fillStyle = "#374151";
                rCtx.fillText(p + 1, x, RULER_H - 9);
            }
        }

        // --- Names ---
        nameCanvas.height = H;
        const nCtx = nameCanvas.getContext("2d");
        nCtx.fillStyle = "#f3f4f6";
        nCtx.fillRect(0, 0, nameW, H);
        const nameFontSize = Math.min(12, cellH - 2);
        nCtx.font = `${nameFontSize}px system-ui,sans-serif`;
        nCtx.textBaseline = "middle";
        nCtx.textAlign = "right";
        for (let si = 0; si < nSeq; si++) {
            const y = si * cellH;
            if (si % 2 === 0) { nCtx.fillStyle = "#eef0f5"; nCtx.fillRect(0, y, nameW, cellH); }
            nCtx.fillStyle = "#374151";
            const name = sequences[si].name;
            const maxChars = Math.floor((nameW - 10) / (nameFontSize * 0.58));
            nCtx.fillText(name.length > maxChars ? name.slice(0, maxChars - 1) + "…" : name, nameW - 6, y + cellH / 2);
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

        // --- Base composition: sequence logo (letters scaled by frequency) ---
        baseCompCanvas.width = W;
        const bcCtx = baseCompCanvas.getContext("2d");
        bcCtx.fillStyle = "#f9fafb";
        bcCtx.fillRect(0, 0, W, BASE_COMP_H);
        const refSize = BASE_COMP_H;
        const isProtein = !!options.isProtein;
        // Nucleotide-specific vivid colors; proteins use palette
        const ntCompColors = { A: "#1AC253", T: "#E11218", G: "#ECA918", C: "#4272EE" };
        const compOrder = isProtein
            ? ["L", "A", "G", "V", "E", "S", "I", "K", "R", "D", "T", "P", "N", "F", "Q", "Y", "M", "H", "C", "W", "*"]
            : ["A", "T", "G", "C"];
        const compColors = isProtein ? palette : ntCompColors;
        for (let pi = 0; pi < seqLen; pi++) {
            const counts = {};
            for (const k of compOrder) counts[k] = 0;
            for (let si = 0; si < nSeq; si++) {
                const ch = (sequences[si].seq[pi] || "-").toUpperCase();
                const mapped = (!isProtein && ch === "U") ? "T" : ch;
                if (mapped in counts) counts[mapped]++;
            }
            const total = compOrder.reduce((s, k) => s + counts[k], 0);
            if (total === 0) continue;
            // Sort: most frequent drawn first (bottom)
            const sorted = compOrder.map(k => [k, counts[k]]).filter(([, v]) => v > 0)
                .sort((a, b) => b[1] - a[1]);
            let yPos = BASE_COMP_H;
            for (const [res, count] of sorted) {
                const frac = count / total;
                const bh = Math.max(1, Math.round(frac * BASE_COMP_H));
                // Only draw letter if tall enough to be legible (≥6px)
                if (cellW >= 4 && bh >= 6) {
                    bcCtx.save();
                    const cx = pi * cellW + cellW / 2;
                    const cy = yPos - bh / 2;
                    const yScale = bh / refSize;
                    const xScale = Math.min(1.2, (cellW - gapX) / (refSize * 0.65));
                    bcCtx.translate(cx, cy);
                    bcCtx.scale(xScale, yScale);
                    bcCtx.fillStyle = compColors[res] || "#cbd5e1";
                    bcCtx.font = `bold ${refSize}px monospace`;
                    bcCtx.textAlign = "center";
                    bcCtx.textBaseline = "middle";
                    bcCtx.fillText(res, 0, 0);
                    bcCtx.restore();
                } else {
                    bcCtx.fillStyle = compColors[res] || "#cbd5e1";
                    bcCtx.fillRect(pi * cellW, yPos - bh, Math.max(1, cellW - gapX), bh);
                }
                yPos -= bh;
            }
        }
    }

    // Initial draw
    redraw();

    // ── scroll sync ────────────────────────────────────────────────────────
    seqOuter.addEventListener("scroll", () => {
        rulerOuter.scrollLeft = seqOuter.scrollLeft;
        consSeqOuter.scrollLeft = seqOuter.scrollLeft;
        consBarOuter.scrollLeft = seqOuter.scrollLeft;
        baseCompOuter.scrollLeft = seqOuter.scrollLeft;
    });
    bodyRow.addEventListener("scroll", () => {
        nameDiv.scrollTop = bodyRow.scrollTop;
    });

    // ── resize handle drag ─────────────────────────────────────────────────
    let isResizing = false, resizeStartX = 0, resizeStartW = nameW;
    resizeHandle.addEventListener("mousedown", (e) => {
        isResizing = true; resizeStartX = e.clientX; resizeStartW = nameW;
        document.body.style.cursor = "col-resize";
        document.body.style.userSelect = "none";
        e.preventDefault();
    });
    document.addEventListener("mousemove", (e) => {
        if (!isResizing) return;
        nameW = Math.max(80, Math.min(420, resizeStartW + (e.clientX - resizeStartX)));
        updateNameWidth();
        redraw();
    });
    document.addEventListener("mouseup", () => {
        if (!isResizing) return;
        isResizing = false;
        document.body.style.cursor = "";
        document.body.style.userSelect = "";
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

// ── Pipeline display helpers ──────────────────────────────────────────────

function buildPipelineFlowHtml(pipeline, hasTrimmed) {
    if (!pipeline) return "";
    const mode = pipeline.mode || "nt";
    const isCodon = mode === "codon";
    const isAa = mode === "aa";

    const node = (icon, label, color) => {
        const cls = {
            gray: "bg-gray-50 border-gray-200 text-gray-600",
            blue: "bg-blue-50 border-blue-200 text-blue-700",
            amber: "bg-amber-50 border-amber-200 text-amber-700",
            violet: "bg-violet-50 border-violet-200 text-violet-700",
            green: "bg-green-50 border-green-200 text-green-700",
            teal: "bg-teal-50 border-teal-200 text-teal-700",
        }[color] || "bg-gray-50 border-gray-200 text-gray-600";
        return `<div class="flex-shrink-0 border rounded-lg px-2.5 py-1.5 text-xs font-medium ${cls}"><i class="fa-solid ${icon} mr-1"></i>${label}</div>`;
    };
    const arrow = `<div class="flex-shrink-0 text-gray-300 text-sm mx-1 self-center">→</div>`;

    let steps = [];

    if (isCodon) {
        steps = [
            node("fa-dna", "CDS (DNA)", "gray"), arrow,
            node("fa-right-left", "Translate", "amber"), arrow,
            node("fa-circle-dot", "Protein", "gray"), arrow,
            node("fa-align-left", "MAFFT --amino", "blue"), arrow,
            node("fa-align-justify", "Aligned protein", "gray"), arrow,
            node("fa-rotate-left", "Back-translate", "violet"), arrow,
            node("fa-dna", "Codon DNA", "green"),
        ];
        if (hasTrimmed) {
            steps.push(arrow, node("fa-scissors", "trimAl -backtrans", "teal"), arrow, node("fa-dna", "Trimmed codon DNA", "green"));
        }
    } else if (isAa) {
        steps = [
            node("fa-dna", "DNA", "gray"), arrow,
            node("fa-right-left", "Translate", "amber"), arrow,
            node("fa-align-left", "MAFFT --amino", "blue"), arrow,
            node("fa-align-justify", "Aligned protein", "green"),
        ];
        if (hasTrimmed) {
            steps.push(arrow, node("fa-scissors", "trimAl", "teal"), arrow, node("fa-align-justify", "Trimmed protein", "green"));
        }
    } else {
        steps = [
            node("fa-dna", "DNA", "gray"), arrow,
            node("fa-align-left", "MAFFT", "blue"), arrow,
            node("fa-align-justify", "Aligned DNA", "green"),
        ];
        if (hasTrimmed) {
            steps.push(arrow, node("fa-scissors", "trimAl", "teal"), arrow, node("fa-align-justify", "Trimmed DNA", "green"));
        }
    }

    return `<div class="flex items-center flex-wrap gap-y-1 py-2 mb-3 overflow-x-auto">${steps.join("")}</div>`;
}

function buildPipelineBadgesHtml(pipeline) {
    if (!pipeline) return "";
    const badges = [];
    // Mode
    const modeColors = { codon: ["Codon-aware", "violet"], aa: ["Amino acid", "blue"], nt: ["Nucleotide", "teal"] };
    const [modeLabel, modeColor] = modeColors[pipeline.mode] || ["Unknown", "gray"];
    badges.push(`<span class="inline-flex items-center gap-1 px-2 py-0.5 rounded-full text-xs font-semibold bg-${modeColor}-100 text-${modeColor}-800"><i class="fa-solid fa-dna text-[0.6rem]"></i>${modeLabel}</span>`);
    // Input type
    const inputLabel = pipeline.rawDataType === "cds-dna" ? "CDS DNA" : "DNA";
    badges.push(`<span class="inline-flex items-center px-2 py-0.5 rounded-full text-xs bg-gray-100 text-gray-600">Input: ${inputLabel}</span>`);
    // MAFFT mode
    const mafftLabel = pipeline.mafftMode === "--amino" ? "MAFFT --amino" : "MAFFT DNA";
    badges.push(`<span class="inline-flex items-center gap-1 px-2 py-0.5 rounded-full text-xs bg-indigo-100 text-indigo-700"><i class="fa-solid fa-align-left text-[0.6rem]"></i>${mafftLabel}</span>`);
    // Genetic code
    if (pipeline.geneticCodeName) {
        badges.push(`<span class="inline-flex items-center px-2 py-0.5 rounded-full text-xs bg-teal-50 text-teal-700 border border-teal-200">${pipeline.geneticCodeName}</span>`);
    }
    // Back-translation
    if (pipeline.backtranslation?.status === "done") {
        badges.push(`<span class="inline-flex items-center gap-1 px-2 py-0.5 rounded-full text-xs bg-amber-100 text-amber-700"><i class="fa-solid fa-right-left text-[0.6rem]"></i>${pipeline.backtranslation.method}</span>`);
    }
    // Stop codons
    const pol = pipeline.stopCodonPolicy;
    if (pol) {
        if (pol.terminalStopsRemoved > 0) {
            badges.push(`<span class="inline-flex items-center gap-1 px-2 py-0.5 rounded-full text-xs bg-red-100 text-red-700"><i class="fa-solid fa-stop text-[0.6rem]"></i>${pol.terminalStopsRemoved} terminal stop${pol.terminalStopsRemoved !== 1 ? "s" : ""} removed</span>`);
        }
        if (pol.incompleteStopFragmentsRemoved > 0) {
            badges.push(`<span class="inline-flex items-center gap-1 px-2 py-0.5 rounded-full text-xs bg-orange-100 text-orange-700"><i class="fa-solid fa-scissors text-[0.6rem]"></i>${pol.incompleteStopFragmentsRemoved} T/TA fragment${pol.incompleteStopFragmentsRemoved !== 1 ? "s" : ""} removed</span>`);
        }
    }
    // Final matrix type
    const finalColors = { "codon-preserving DNA": "green", protein: "blue", nucleotide: "cyan" };
    const fColor = finalColors[pipeline.finalMatrixType] || "gray";
    badges.push(`<span class="inline-flex items-center px-2 py-0.5 rounded-full text-xs bg-${fColor}-100 text-${fColor}-800">Matrix: ${pipeline.finalMatrixType}</span>`);
    // Trimming
    if (pipeline.trimming?.status === "done") {
        badges.push(`<span class="inline-flex items-center gap-1 px-2 py-0.5 rounded-full text-xs bg-orange-100 text-orange-700"><i class="fa-solid fa-scissors text-[0.6rem]"></i>${pipeline.trimming.strategy || "trimAl"}${pipeline.trimming.codonPreserving ? " · frame-safe" : ""}</span>`);
    }
    // Codon integrity
    if (pipeline.codonIntegrity?.checked) {
        const ok =
            pipeline.codonIntegrity.equalLengths !== false &&
            pipeline.codonIntegrity.allLengthsMultipleOf3 !== false &&
            pipeline.codonIntegrity.gapBlocksAreTriplets !== false &&
            pipeline.codonIntegrity.headersMatch !== false &&
            pipeline.codonIntegrity.dnaCompatible !== false &&
            pipeline.codonIntegrity.frameDisrupted !== true;
        if (ok) {
            badges.push(`<span class="inline-flex items-center gap-1 px-2 py-0.5 rounded-full text-xs bg-green-100 text-green-700"><i class="fa-solid fa-check text-[0.6rem]"></i>Frame OK</span>`);
        } else {
            const issues = [];
            if (!pipeline.codonIntegrity.equalLengths) issues.push("unequal lengths");
            if (!pipeline.codonIntegrity.allLengthsMultipleOf3) issues.push("length ÷ 3 ≠ 0");
            if (!pipeline.codonIntegrity.gapBlocksAreTriplets) issues.push("gaps not in triplets");
            if (pipeline.codonIntegrity.headersMatch === false) issues.push("header mismatch");
            if (pipeline.codonIntegrity.dnaCompatible === false) issues.push("non-DNA output");
            if (pipeline.codonIntegrity.frameDisrupted === true) issues.push("frame disrupted");
            badges.push(`<span class="inline-flex items-center gap-1 px-2 py-0.5 rounded-full text-xs bg-red-100 text-red-700"><i class="fa-solid fa-triangle-exclamation text-[0.6rem]"></i>Frame: ${issues.join(", ")}</span>`);
        }
    }
    return badges.join("");
}

function renderAlignmentResults(markerResults, hasTrimal) {
    state.alignmentResultsData = { markerResults, hasTrimal };
    const accordion = document.getElementById("alignmentResultsAccordion");
    if (!accordion) return;
    accordion.innerHTML = "";
    document.getElementById("alignmentResults").classList.remove("hidden");

    for (const r of markerResults) {
        const id = "marker-result-" + r.marker.replace(/[^a-zA-Z0-9]/g, "_");
        // ── Per-mode metadata for display ──────────────────────────────────
        const pMode = r.pipeline?.mode || window._alignmentMode || "nt";
        const isCodon = pMode === "codon";
        const isAa = pMode === "aa";
        const alignUnit = isAa ? "aa" : "bp";
        const rawUnit = isCodon ? "bp CDS" : isAa ? "aa" : "bp";

        // Smart aligned header stat
        const aLen = r.alignedStats?.length ?? "?";
        const aSeqs = r.alignedStats?.numSeqs ?? "?";
        let headerStat;
        if (isCodon && typeof aLen === "number") {
            headerStat = `${aSeqs} seq \u00b7 ${aLen} bp${aLen % 3 === 0 ? ` \u00b7 ${aLen / 3} codons` : " \u26a0"}`;
        } else {
            headerStat = `${aSeqs} seq \u00b7 ${aLen} ${alignUnit}`;
        }

        // Smart trimStat
        const trimStat = r.trimmedStats
            ? (() => {
                const tLen = r.trimmedStats.length;
                let statStr;
                if (isCodon && typeof tLen === "number") {
                    statStr = `${r.trimmedStats.numSeqs} seq \u00b7 ${tLen} bp${tLen % 3 === 0 ? ` \u00b7 ${tLen / 3} codons` : " \u26a0"}`;
                } else {
                    statStr = `${r.trimmedStats.numSeqs} seq \u00b7 ${tLen} ${alignUnit}`;
                }
                return `<div class="text-center"><div class="text-xs text-gray-400 mb-0.5">${i18nText("step5.results.stats.trimmed", "Trimmed")}</div><div class="text-xs font-semibold text-gray-700">${statStr}</div><div class="text-xs text-gray-400">${r.trimmedStats.gapPct}% ${i18nText("step5.results.stats.gaps", "gaps")}</div></div>`;
            })()
            : `<div class="text-center text-xs text-gray-300 italic">${i18nText("step5.results.stats.noTrimming", "No trimming")}</div>`;

        // Codon integrity alert
        let integrityAlert = "";
        if (r.pipeline?.codonIntegrity?.checked) {
            const ci = r.pipeline.codonIntegrity;
            const issues = [];
            if (!ci.equalLengths) issues.push("alignment contains unequal sequence lengths");
            if (!ci.allLengthsMultipleOf3) issues.push("alignment length not multiple of 3");
            if (!ci.gapBlocksAreTriplets) issues.push("gaps not in triplets");
            if (ci.headersMatch === false) issues.push("sequence names do not match the validated CDS records");
            if (ci.dnaCompatible === false) issues.push("trimmed output is not DNA-compatible");
            if (ci.frameDisrupted === true) issues.push("codon frame disruption detected");
            if (issues.length > 0) {
                integrityAlert = `<div class="mt-2 flex items-start gap-2 rounded-lg border border-red-200 bg-red-50 px-3 py-2"><i class="fa-solid fa-triangle-exclamation text-red-500 mt-0.5 text-xs flex-shrink-0"></i><span class="text-xs text-red-700 font-medium">Frame warning: ${issues.join("; ")}</span></div>`;
            }
        }

        // Pipeline badges and flow
        const pipelineBadgesHtml = buildPipelineBadgesHtml(r.pipeline);
        const pipelineFlowHtml = buildPipelineFlowHtml(r.pipeline, !!r.trimmedStats);

        // Mode-aware view buttons per tab
        const buildViewBtns = (tabKey) => {
            if (isCodon && tabKey === "aligned") {
                return `<div class="alignment-view-toggle">
                    <button class="alignment-view-btn active" type="button" data-marker="${r.marker}" data-tab="${tabKey}" data-view="nt">Codon DNA</button>
                    ${r.proteinAlignedContent ? `<button class="alignment-view-btn" type="button" data-marker="${r.marker}" data-tab="${tabKey}" data-view="protein-guide">Protein guide</button>` : ""}
                </div>`;
            }
            if (isCodon) {
                return `<div class="alignment-view-toggle"><button class="alignment-view-btn active" type="button" data-marker="${r.marker}" data-tab="${tabKey}" data-view="nt">Codon DNA</button></div>`;
            }
            if (isAa) {
                return `<div class="alignment-view-toggle"><button class="alignment-view-btn active" type="button" data-marker="${r.marker}" data-tab="${tabKey}" data-view="aa">${i18nText("step5.results.view.aa", "Protein view")}</button></div>`;
            }
            return `<div class="alignment-view-toggle">
                <button class="alignment-view-btn active" type="button" data-marker="${r.marker}" data-tab="${tabKey}" data-view="nt">${i18nText("step5.results.view.nt", "DNA view")}</button>
                <button class="alignment-view-btn" type="button" data-marker="${r.marker}" data-tab="${tabKey}" data-view="aa">${i18nText("step5.results.view.aa", "Translate to amino acids")}</button>
            </div>`;
        };

        // Default initial view per mode
        const defaultView = isAa ? "aa" : "nt";

        // Aligned stat for grid cell
        let alignedStatStr;
        if (isCodon && typeof aLen === "number") {
            alignedStatStr = `${aSeqs} seq \u00b7 ${aLen} bp${aLen % 3 === 0 ? ` \u00b7 ${aLen / 3} codons` : " \u26a0"}`;
        } else {
            alignedStatStr = `${aSeqs} seq \u00b7 ${aLen} ${alignUnit}`;
        }

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
                        ${buildViewBtns(t.key)}
                        <div class="alignment-toolbar-meta">
                            ${!isCodon ? `<span class="alignment-view-hint flex items-center gap-1 text-amber-600"><i class="fa-solid fa-lightbulb text-amber-500"></i>${i18nText("step5.results.view.hint", "Translate to amino acids to check for premature stop codons before phylogenetic inference.")}</span>` : ""}
                            <div class="alignment-genetic-code">
                                <label for="genetic-code-${id}-${t.key}">${i18nText("step5.results.view.geneticCode", "Genetic code")}</label>
                                <select id="genetic-code-${id}-${t.key}" class="alignment-genetic-code-select" data-marker="${r.marker}" data-tab="${t.key}">
                                    <option value="1">${i18nText("step5.results.view.code.standard", "Standard (NCBI 1)")}</option>
                                    <option value="2" selected>${i18nText("step5.results.view.code.vertMito", "Vertebrate mitochondrial (NCBI 2)")}</option>
                                    <option value="3">Yeast mitochondrial (NCBI 3)</option>
                                    <option value="4">Mold mitochondrial (NCBI 4)</option>
                                    <option value="5">${i18nText("step5.results.view.code.invMito", "Invertebrate mitochondrial (NCBI 5)")}</option>
                                    <option value="6">Ciliate nuclear (NCBI 6)</option>
                                    <option value="9">Echinoderm mitochondrial (NCBI 9)</option>
                                    <option value="11">${i18nText("step5.results.view.code.bacterial", "Bacterial, archaeal and plant plastid (NCBI 11)")}</option>
                                    <option value="13">Ascidian mitochondrial (NCBI 13)</option>
                                    <option value="14">Alternative flatworm mitochondrial (NCBI 14)</option>
                                    <option value="16">Chlorophycean mitochondrial (NCBI 16)</option>
                                    <option value="21">Trematode mitochondrial (NCBI 21)</option>
                                    <option value="24">Pterobranchia mitochondrial (NCBI 24)</option>
                                </select>
                            </div>
                        </div>
                    </div>
                    <div id="vizwrap-${id}-${t.key}" data-view="${defaultView}"></div>
                    <div id="preview-${id}-${t.key}"></div>
                </div>
            </div>`;
        }).join("");

        const panel = document.createElement("div");
        panel.className = "bg-white rounded-xl border border-gray-200 overflow-hidden";
        panel.innerHTML = `
            <button class="w-full px-5 py-3 flex items-center justify-between text-left hover:bg-gray-50 transition-colors" onclick="toggleMarkerResult('${id}')">
                <div class="flex items-center gap-3 min-w-0">
                    <span class="step-badge badge-done flex-shrink-0" style="font-size:0.55rem;width:1rem;height:1rem;">\u2713</span>
                    <span class="font-semibold text-sm text-gray-800">${r.marker}</span>
                    <span class="text-xs text-gray-400">${headerStat} ${i18nText("step5.results.stats.alignedSuffix", "aligned")}</span>
                </div>
                <i id="${id}-chevron" class="fa-solid fa-chevron-down text-gray-400 text-xs transition-transform flex-shrink-0"></i>
            </button>
            <div id="${id}-body" class="hidden px-5 pb-5">
                ${pipelineBadgesHtml ? `<div class="flex flex-wrap gap-1.5 pt-3 pb-1">${pipelineBadgesHtml}</div>` : ""}
                ${pipelineFlowHtml}
                ${integrityAlert}
                <div class="grid grid-cols-3 gap-3 mb-4 pt-2">
                    <div class="text-center"><div class="text-xs text-gray-400 mb-0.5">${i18nText("step5.results.stats.raw", "Raw")}</div><div class="text-xs font-semibold text-gray-700">${r.rawStats?.numSeqs ?? "?"} seq \u00b7 ${r.rawStats?.avgLen ?? "?"} ${rawUnit} ${i18nText("step5.results.stats.avg", "avg")}</div></div>
                    <div class="text-center"><div class="text-xs text-gray-400 mb-0.5">${i18nText("step5.results.stats.aligned", "Aligned")}</div><div class="text-xs font-semibold text-gray-700">${alignedStatStr}</div><div class="text-xs text-gray-400">${r.alignedStats?.gapPct ?? "?"}% ${i18nText("step5.results.stats.gaps", "gaps")}</div></div>
                    ${trimStat}
                </div>
                <div class="flex gap-2 mb-3">${tabBtns}</div>
                <div class="bg-gray-50 rounded-lg p-2">${panels}</div>
            </div>`;
        accordion.appendChild(panel);

        // Render raw tab immediately (it's visible by default)
        requestAnimationFrame(() => {
            if (r.rawContent) renderAlignmentTab({ id, rawContent: r.rawContent, alignedContent: r.alignedContent, trimmedContent: r.trimmedContent, proteinAlignedContent: r.proteinAlignedContent, pipeline: r.pipeline }, "raw", defaultView);
        });

        // Store content references for lazy rendering on tab switch / accordion open
        panel._markerData = {
            id,
            rawContent: r.rawContent,
            alignedContent: r.alignedContent,
            trimmedContent: r.trimmedContent,
            proteinAlignedContent: r.proteinAlignedContent,
            pipeline: r.pipeline,
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
                    const defView = wrap.dataset.view || (data.pipeline?.mode === "aa" ? "aa" : "nt");
                    requestAnimationFrame(() => renderAlignmentTab(data, tab, defView));
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

window.toggleThirdCodonWarning = function () {
    const warning = document.getElementById("excludeThirdCodonWarning");
    const checked = document.getElementById("excludeThirdCodonPosition")?.checked === true;
    if (warning) warning.classList.toggle("hidden", !checked);
};

function buildConcatenatedFasta(alignedSeqs) {
    return (alignedSeqs || []).map((seq) => `>${seq.name}\n${seq.seq}`).join("\n") + "\n";
}

function removeThirdCodonPositionsFromSequence(sequence) {
    if (!sequence) return "";
    if (sequence.length % 3 !== 0) {
        throw new Error(`Cannot remove third codon positions from sequence length ${sequence.length}: length is not divisible by 3`);
    }
    let trimmed = "";
    for (let idx = 0; idx < sequence.length; idx += 3) {
        const codon = sequence.slice(idx, idx + 3);
        if (codon.length !== 3) {
            throw new Error(`Incomplete codon chunk encountered at offset ${idx + 1}`);
        }
        trimmed += codon[0] + codon[1];
    }
    return trimmed;
}

function buildNoThirdCodonMatrix(concatResult) {
    const geneSeqs = {};
    const geneLens = {};
    const partitions = [];
    const logs = [];
    let pos = 1;
    let removedThirdPositions = 0;

    for (const gene of concatResult.geneOrder) {
        const originalLength = concatResult.geneLens[gene];
        if (!originalLength) continue;
        if (originalLength % 3 !== 0) {
            throw new Error(`Cannot exclude third codon positions for ${gene}: alignment length ${originalLength} is not divisible by 3`);
        }
        const noThirdLength = (originalLength / 3) * 2;
        const removedForGene = originalLength / 3;
        removedThirdPositions += removedForGene;
        logs.push(`[No-3rd] [${gene}] removed ${removedForGene} third positions (${originalLength} bp → ${noThirdLength} bp)`);

        geneSeqs[gene] = {};
        for (const [species, seq] of Object.entries(concatResult.geneSeqs[gene] || {})) {
            geneSeqs[gene][species] = removeThirdCodonPositionsFromSequence(seq);
        }
        geneLens[gene] = noThirdLength;
        partitions.push({
            gene,
            start: pos,
            end: pos + noThirdLength - 1,
            len: noThirdLength,
            originalLength,
            noThirdLength,
            removedThirdPositions: removedForGene,
            thirdCodonPositionExcluded: true,
            codonPositionMode: "positions 1 and 2 only",
        });
        pos += noThirdLength;
    }

    return {
        ...concatResult,
        alignedSeqs: concatResult.alignedSeqs.map((seq) => ({
            name: seq.name,
            seq: removeThirdCodonPositionsFromSequence(seq.seq),
        })),
        partitions,
        totalLen: pos - 1,
        geneSeqs,
        geneLens,
        thirdCodonPositionExcluded: true,
        originalLength: concatResult.totalLen,
        noThirdLength: pos - 1,
        removedThirdPositions,
        codonPositionMode: "positions 1 and 2 only",
        logs,
    };
}

function resolveCodonPartitionMode(codonPositions, excludeThirdPositions) {
    if (!codonPositions) return "genes";
    return excludeThirdPositions ? "positions12" : "all";
}

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

    return { alignedSeqs, partitions, totalLen: pos - 1, species: finalSpecies, missingReport, removedCount: allowMissing ? 0 : missingReport.length, geneSeqs, geneLens, geneOrder, allSpecies };
}

function validateCodonPartitions(markerResults, useTrimmed) {
    const reports = markerResults.map((markerResult) => {
        const content = useTrimmed && markerResult.trimmedContent ? markerResult.trimmedContent : markerResult.alignedContent;
        const seqs = parseFasta(content || "");
        if (!seqs.length) {
            return {
                marker: markerResult.marker,
                valid: false,
                issues: ["no sequences available"],
                length: 0,
            };
        }

        const integrity = getCodonIntegrityChecks(seqs);
        const frameDisrupted = markerResult.pipeline?.codonIntegrity?.frameDisrupted === true || markerResult.pipeline?.backtranslation?.frameDisrupted === true;
        const issues = summarizeCodonIntegrityIssues(integrity, frameDisrupted);

        return {
            marker: markerResult.marker,
            valid: issues.length === 0,
            issues,
            length: integrity.length,
        };
    });

    return {
        valid: reports.every((report) => report.valid),
        reports,
        errors: reports.filter((report) => !report.valid),
    };
}

function showCodonPartitionValidationError(errors) {
    const warn = document.getElementById("concatMissingWarning");
    if (!warn) return;

    warn.classList.remove("hidden");
    warn.className = "text-xs rounded-lg px-3 py-2 border text-red-700 bg-red-50 border-red-200";
    warn.innerHTML = `<strong>Codon partition validation failed.</strong> ${errors.map((error) => `${escHtml(error.marker)}: ${escHtml(error.issues.join("; "))}`).join(" | ")}`;
}

function buildNexusString(alignedSeqs, partitions, outgroups, datatype, codonPartitionMode) {
    const dtStr = ((datatype || "dna")).toLowerCase();
    const ntax = alignedSeqs.length;
    const nchar = alignedSeqs[0]?.seq.length || 0;
    // Pad names for alignment
    const maxNameLen = Math.max(...alignedSeqs.map(s => s.name.length));
    const matrix = alignedSeqs.map(s =>
        `    ${s.name.padEnd(maxNameLen + 2)}${s.seq}`
    ).join("\n");

    const partitionMode = codonPartitionMode === true ? "all" : (codonPartitionMode || "genes");
    let charsets;
    if (partitionMode === "all") {
        charsets = partitions.flatMap(p => {
            const geneLen = p.len || (p.end - p.start + 1);
            if (geneLen % 3 !== 0) {
                throw new Error(`Cannot create codon charsets for ${p.gene}: alignment length ${geneLen} is not divisible by 3`);
            }
            return [
                `    charset ${p.gene}_pos1 = ${p.start}-${p.end}\\3;`,
                `    charset ${p.gene}_pos2 = ${p.start + 1}-${p.end}\\3;`,
                `    charset ${p.gene}_pos3 = ${p.start + 2}-${p.end}\\3;`,
            ];
        }).join("\n");
    } else if (partitionMode === "positions12") {
        charsets = partitions.flatMap(p => {
            const geneLen = p.len || (p.end - p.start + 1);
            if (geneLen % 2 !== 0) {
                throw new Error(`Cannot create 1+2 codon charsets for ${p.gene}: alignment length ${geneLen} is not divisible by 2`);
            }
            return [
                `    charset ${p.gene}_pos1 = ${p.start}-${p.end}\\2;`,
                `    charset ${p.gene}_pos2 = ${p.start + 1}-${p.end}\\2;`,
            ];
        }).join("\n");
    } else {
        charsets = partitions.map(p =>
            `    charset ${p.gene} = ${p.start}-${p.end};`
        ).join("\n");
    }

    const outgroupBlock = outgroups && outgroups.length
        ? `\nbegin assumptions;\n    outgroup ${outgroups.join(" ")};\nend;\n`
        : "";

    return `#NEXUS

begin data;
    dimensions ntax=${ntax} nchar=${nchar};
    format datatype=${dtStr} missing=? gap=-;
    matrix
${matrix}
    ;
end;

begin sets;
${charsets}
end;
${outgroupBlock}`;
}

function buildPartitionString(partitions, codonPartitionMode, datatype) {
    const prefix = (datatype === "protein") ? "AA" : "DNA";
    const partitionMode = codonPartitionMode === true ? "all" : (codonPartitionMode || "genes");
    if (partitionMode === "all") {
        return partitions.flatMap(p => {
            const geneLen = p.len || (p.end - p.start + 1);
            if (geneLen % 3 !== 0) {
                throw new Error(`Cannot create codon partitions for ${p.gene}: alignment length ${geneLen} is not divisible by 3`);
            }
            return [
                `${prefix}, ${p.gene}_pos1 = ${p.start}-${p.end}\\3`,
                `${prefix}, ${p.gene}_pos2 = ${p.start + 1}-${p.end}\\3`,
                `${prefix}, ${p.gene}_pos3 = ${p.start + 2}-${p.end}\\3`,
            ];
        }).join("\n");
    }
    if (partitionMode === "positions12") {
        return partitions.flatMap(p => {
            const geneLen = p.len || (p.end - p.start + 1);
            if (geneLen % 2 !== 0) {
                throw new Error(`Cannot create 1+2 codon partitions for ${p.gene}: alignment length ${geneLen} is not divisible by 2`);
            }
            return [
                `${prefix}, ${p.gene}_pos1 = ${p.start}-${p.end}\\2`,
                `${prefix}, ${p.gene}_pos2 = ${p.start + 1}-${p.end}\\2`,
            ];
        }).join("\n");
    }
    return partitions.map(p => `${prefix}, ${p.gene} = ${p.start}-${p.end}`).join("\n");
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

    // Show codon partition option only when in codon mode
    const codonPartitionOpt = document.getElementById("codonPartitionOption");
    if (codonPartitionOpt) codonPartitionOpt.classList.toggle("hidden", !window._codonMode);
    const excludeThirdOpt = document.getElementById("excludeThirdCodonOption");
    const excludeThirdCheckbox = document.getElementById("excludeThirdCodonPosition");
    if (excludeThirdOpt) excludeThirdOpt.classList.toggle("hidden", !window._codonMode);
    if (!window._codonMode && excludeThirdCheckbox) excludeThirdCheckbox.checked = false;
    toggleThirdCodonWarning();

    const generationLog = document.getElementById("concatGenerationLog");
    if (generationLog) {
        generationLog.classList.add("hidden");
        generationLog.textContent = "";
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
    const activeFiles = window._nexusData?.activeFiles || {
        nexus: "concatenated.nex",
        partition: "partitions.txt",
        prefix: "concat_tree",
    };

    const parts = [
        "iqtree3",
        "-s", `${dir}/${activeFiles.nexus}`,
        "-p", `${dir}/${activeFiles.partition}`,
        "--prefix", `${dir}/${activeFiles.prefix}`,
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

    // Detect alignment mode for NEXUS datatype and codon partitions
    const isCodonMode = window._codonMode === true;
    const isAaMode = window._alignmentMode === "aa";
    const datatype = isAaMode ? "protein" : "dna";
    const codonPositions = isCodonMode && (document.getElementById("codonPartitionByPosition")?.checked ?? false);
    const excludeThirdPositions = isCodonMode && (document.getElementById("excludeThirdCodonPosition")?.checked ?? false);

    const result = buildConcat(_lastMarkerResultsUI, useTrimmed, allowMissing);

    if (isCodonMode) {
        const codonValidation = validateCodonPartitions(_lastMarkerResultsUI, useTrimmed);
        if (!codonValidation.valid) {
            showCodonPartitionValidationError(codonValidation.errors);
            return;
        }
    }

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

    let nexus, partition;
    let fasta;
    let activeResult = result;
    let activeFiles = { nexus: "concatenated.nex", partition: "partitions.txt", fasta: "concatenated.fasta", prefix: "concat_tree" };
    let fullMatrix = null;
    let noThirdMatrix = null;
    let additionalFiles = {};
    try {
        const fullPartitionMode = resolveCodonPartitionMode(codonPositions, false);
        fullMatrix = {
            ...result,
            nexus: buildNexusString(result.alignedSeqs, result.partitions, outgroups, datatype, fullPartitionMode),
            partition: buildPartitionString(result.partitions, fullPartitionMode, datatype),
            fasta: buildConcatenatedFasta(result.alignedSeqs),
            files: { nexus: "concatenated.nex", partition: "partitions.txt", fasta: "concatenated.fasta", prefix: "concat_tree" },
            thirdCodonPositionExcluded: false,
            codonPositionMode: codonPositions ? "all positions" : "gene partitions",
        };

        if (excludeThirdPositions) {
            const noThirdPartitionMode = resolveCodonPartitionMode(codonPositions, true);
            const noThirdResult = buildNoThirdCodonMatrix(result);
            noThirdMatrix = {
                ...noThirdResult,
                nexus: buildNexusString(noThirdResult.alignedSeqs, noThirdResult.partitions, outgroups, "dna", noThirdPartitionMode),
                partition: buildPartitionString(noThirdResult.partitions, noThirdPartitionMode, "dna"),
                fasta: buildConcatenatedFasta(noThirdResult.alignedSeqs),
                files: { nexus: "concatenated_no3rd.nex", partition: "partitions_no3rd.txt", fasta: "concatenated_no3rd.fasta", prefix: "concat_tree_no3rd" },
            };
            activeResult = noThirdMatrix;
            activeFiles = noThirdMatrix.files;
            additionalFiles = {
                [fullMatrix.files.nexus]: fullMatrix.nexus,
                [fullMatrix.files.partition]: fullMatrix.partition,
                [fullMatrix.files.fasta]: fullMatrix.fasta,
            };
        } else {
            activeResult = fullMatrix;
        }

        nexus = activeResult.nexus;
        partition = activeResult.partition;
        fasta = activeResult.fasta;
    } catch (error) {
        showCodonPartitionValidationError([{ marker: "matrix", issues: [error.message] }]);
        return;
    }

    // Store for IQ-TREE section
    window._nexusData = {
        nexus,
        partition,
        fasta,
        result: activeResult,
        outgroups,
        activeFiles,
        fullMatrix,
        noThirdMatrix,
        additionalFiles,
        thirdCodonPositionExcluded: excludeThirdPositions,
    };

    // Show sequence × marker heatmap table
    const tbl = document.getElementById("partitionTable");
    if (tbl) {
        const markers = activeResult.partitions.map(p => p.gene);
        const allSp = activeResult.allSpecies || [];
        const heatColor = (cov) => {
            const r = Math.round(254 - 34 * cov);
            const g = Math.round(226 + 26 * cov);
            const b = Math.round(226 + 5 * cov);
            return `rgb(${r},${g},${b})`;
        };
        const heatRows = allSp.map(sp => {
            const nameTd = `<td style="position:sticky;left:0;z-index:1;background:#f8fafc;font-weight:600;padding:4px 10px;border:1px solid #e5e7eb;white-space:nowrap;">${escHtml(sp)}</td>`;
            const dataTds = markers.map(gene => {
                const seq = activeResult.geneSeqs[gene]?.[sp];
                if (!seq) return `<td style="padding:4px 8px;border:1px solid #e5e7eb;background:#f3f4f6;color:#9ca3af;text-align:center;">-</td>`;
                const nonGap = seq.replace(/-/g, "").length;
                const gaps = seq.length - nonGap;
                const cov = activeResult.geneLens[gene] > 0 ? nonGap / activeResult.geneLens[gene] : 0;
                return `<td style="padding:4px 8px;border:1px solid #e5e7eb;background:${heatColor(cov)};text-align:center;color:#1f2937;">${nonGap} (${gaps})</td>`;
            }).join("");
            return `<tr>${nameTd}${dataTds}</tr>`;
        }).join("");
        tbl.innerHTML = `
            <div style="margin-bottom:6px;font-size:0.72rem;color:#6b7280;display:flex;align-items:center;gap:6px;">
                Cells: <strong>bp (gaps)</strong>
                <span style="display:inline-flex;align-items:center;gap:4px;">
                    <span style="display:inline-block;width:60px;height:8px;background:linear-gradient(to right,#fee2e2,#dcfce7);border-radius:2px;"></span>
                    <span>low → high coverage</span>
                    <span style="background:#f3f4f6;color:#9ca3af;padding:1px 6px;border-radius:3px;">-</span><span>missing</span>
                </span>
            </div>
            <div style="overflow:auto;max-height:420px;">
                <table style="border-collapse:collapse;font-size:0.72rem;white-space:nowrap;">
                    <thead>
                        <tr>
                            <th style="position:sticky;left:0;top:0;z-index:3;background:#f3f4f6;padding:6px 10px;border:1px solid #e5e7eb;font-weight:700;text-align:left;">Sequence</th>
                            ${markers.map(g => `<th style="position:sticky;top:0;z-index:2;background:#f3f4f6;padding:6px 8px;border:1px solid #e5e7eb;font-weight:700;text-align:center;" title="${escHtml(g)}">${escHtml(g)}</th>`).join("")}
                        </tr>
                    </thead>
                    <tbody>${heatRows}</tbody>
                </table>
            </div>`;
        tbl.classList.remove("hidden");
    }

    // Show summary
    const summary = document.getElementById("concatSummary");
    if (summary) {
        summary.textContent = `${activeResult.alignedSeqs.length} species · ${activeResult.partitions.length} genes · ${activeResult.totalLen} bp${excludeThirdPositions ? ` · positions 1 and 2 only (removed ${activeResult.removedThirdPositions} third positions)` : ""}`;
        summary.classList.remove("hidden");
    }

    const generationLog = document.getElementById("concatGenerationLog");
    if (generationLog) {
        if (excludeThirdPositions && activeResult.logs?.length) {
            generationLog.textContent = `${activeResult.logs.join("\n")}\n[No-3rd] total removed third positions: ${activeResult.removedThirdPositions}`;
            generationLog.classList.remove("hidden");
        } else {
            generationLog.classList.add("hidden");
            generationLog.textContent = "";
        }
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
    a.download = window._nexusData.activeFiles?.nexus || "concatenated.nex";
    a.click();
};

window.downloadPartition = function () {
    if (!window._nexusData) return;
    const blob = new Blob([window._nexusData.partition], { type: "text/plain" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = window._nexusData.activeFiles?.partition || "partitions.txt";
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
        fasta: window._nexusData.fasta,
        params,
        threads,
        outgroup,
        perGene,
        geneFiles,
        fileNames: window._nexusData.activeFiles,
        extraFiles: window._nexusData.additionalFiles,
    });
};

function openIqtreeModal() {
    const modal = document.getElementById("iqtreeModal");
    if (!modal) return;
    modal.classList.remove("hidden");
    setModalEscapeEnabled("iqtreeModal", false);
    document.getElementById("iqtreeModalLog").textContent = "";
    document.getElementById("iqtreeModalTitle").textContent = "Running IQ-TREE…";
    document.getElementById("iqtreeModalClose").classList.add("hidden");
    document.getElementById("iqtreeModalIcon").className = "fa-solid fa-spinner fa-spin text-splace-blue-600 text-lg";
    refreshManagedModalState("iqtreeModal");
}

if (window.electronAPI?.onIqtreeProgress) {
    window.electronAPI.onIqtreeProgress((data) => {
        const log = document.getElementById("iqtreeModalLog");
        if (log && data.message) {
            log.textContent += data.message + "\n";
            // scroll the overflow container (parent of <pre>), not the <pre> itself
            const logContainer = log.parentElement;
            if (logContainer) logContainer.scrollTop = logContainer.scrollHeight;
            refreshManagedModalState("iqtreeModal");
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
        setModalEscapeEnabled("iqtreeModal", true);
        refreshManagedModalState("iqtreeModal");
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

