// ── SPLACE Citation Registry ──────────────────────────────────────────────────
// Centralised citation data for all tools used in SPLACE.
// Provides window.CITATIONS (data) and window.renderCitationsCard() (renderer).

window.CITATIONS = {

    splace: {
        id: "splace",
        name: "SPLACE",
        authors: "Oliveira RR, Vasconcelos S, Oliveira G",
        title: "SPLACE: A tool to automatically SPLit, Align, and ConcatenatE genes for phylogenomic inference of several organisms",
        journal: "Frontiers in Bioinformatics",
        year: "2022",
        volume: "2",
        doi: "10.3389/fbinf.2022.1074802",
        get cite() {
            return `${this.authors} (${this.year}) ${this.title}. ${this.journal}, ${this.volume}. https://doi.org/${this.doi}`;
        },
    },

    syngenes: {
        id: "syngenes",
        name: "SynGenes",
        authors: "Rabelo et al.",
        title: "SynGenes: a Python class for standardizing nomenclatures of mitochondrial and chloroplast genes and a web form for enhancing searches for evolutionary analyses",
        journal: "BMC Bioinformatics",
        year: "2024",
        volume: "25",
        pages: "160",
        doi: "10.1186/s12859-024-05781-y",
        get cite() {
            return `${this.authors} (${this.year}) ${this.title}. ${this.journal} ${this.volume}:${this.pages}. https://doi.org/${this.doi}`;
        },
    },

    datafishing: {
        id: "datafishing",
        name: "dataFishing",
        authors: "Rabelo et al.",
        title: "dataFishing: An efficient Python tool and user-friendly web-form for mining mitochondrial and chloroplast sequences, taxonomic, and biodiversity data",
        journal: "Ecological Informatics",
        year: "2025",
        volume: "85",
        pages: "102970",
        doi: "10.1016/j.ecoinf.2024.102970",
        get cite() {
            return `${this.authors} (${this.year}) ${this.title}. ${this.journal} ${this.volume}:${this.pages}. https://doi.org/${this.doi}`;
        },
    },

    mafft: {
        id: "mafft",
        name: "MAFFT",
        authors: "Katoh K, Standley DM",
        title: "MAFFT Multiple Sequence Alignment Software Version 7: Improvements in Performance and Usability",
        journal: "Molecular Biology and Evolution",
        year: "2013",
        volume: "30",
        issue: "4",
        pages: "772–780",
        doi: "10.1093/molbev/mst010",
        get cite() {
            return `${this.authors} (${this.year}) ${this.title}. ${this.journal}, ${this.volume}(${this.issue}):${this.pages}. https://doi.org/${this.doi}`;
        },
    },

    trimal: {
        id: "trimal",
        name: "trimAl",
        authors: "Capella-Gutiérrez S, Silla-Martínez JM, Gabaldón T",
        title: "trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses",
        journal: "Bioinformatics",
        year: "2009",
        volume: "25",
        issue: "15",
        pages: "1972–1973",
        doi: "10.1093/bioinformatics/btp348",
        get cite() {
            return `${this.authors} (${this.year}) ${this.title}. ${this.journal}, ${this.volume}(${this.issue}):${this.pages}. https://doi.org/${this.doi}`;
        },
    },

    iqtree: {
        id: "iqtree",
        name: "IQ-TREE 2",
        authors: "Minh BQ, Schmidt HA, Chernomor O, Schrempf D, Woodhams MD, von Haeseler A, Lanfear R",
        title: "IQ-TREE 2: New Models and Methods for Phylogenetic Inference",
        journal: "Molecular Biology and Evolution",
        year: "2020",
        volume: "37",
        issue: "5",
        pages: "1530–1534",
        doi: "10.1093/molbev/msaa015",
        get cite() {
            return `${this.authors} (${this.year}) ${this.title}. ${this.journal}, ${this.volume}(${this.issue}):${this.pages}. https://doi.org/${this.doi}`;
        },
    },
};

window.getCitationTooltip = function (key) {
    const citation = window.CITATIONS && window.CITATIONS[key];
    if (!citation) return "";
    const lang = (window.SPLACE_I18N && window.SPLACE_I18N.lang) || document.documentElement.lang || "en";
    const labels = {
        en: { reference: "Reference" },
        pt: { reference: "Referência" },
        es: { reference: "Referencia" },
    };
    const label = (labels[lang] || labels.en).reference;
    const doiUrl = `https://doi.org/${citation.doi}`;
    const volInfo = [
        citation.volume,
        citation.issue ? `(${citation.issue})` : "",
        citation.pages ? `:${citation.pages}` : "",
    ].join("");
    return `
        <div class="font-bold text-white">${label}</div>
        <div class="text-white"><strong>${citation.name}</strong></div>
        <div class="text-white">${citation.authors} (${citation.year}). <em>${citation.title}</em>. ${citation.journal}${volInfo ? ` ${volInfo}` : ""}.</div>
        <div class="mt-1 text-white"><a href="${doiUrl}" target="_blank" rel="noopener" class="text-white">https://doi.org/${citation.doi}</a></div>`;
};

// ── Copy-citation helper (called from data-attribute buttons) ─────────────────
window.copyCitation = function (btn) {
    const text = btn.dataset.cite;
    if (!text) return;
    navigator.clipboard?.writeText(text).then(() => {
        const orig = btn.innerHTML;
        btn.innerHTML = '<i class="fa-solid fa-check mr-1"></i>Copied';
        btn.classList.replace('bg-splace-blue-50', 'bg-green-50');
        btn.classList.replace('border-splace-blue-200', 'border-green-200');
        btn.classList.replace('text-splace-blue-700', 'text-green-700');
        setTimeout(() => {
            btn.innerHTML = orig;
            btn.classList.replace('bg-green-50', 'bg-splace-blue-50');
            btn.classList.replace('border-green-200', 'border-splace-blue-200');
            btn.classList.replace('text-green-700', 'text-splace-blue-700');
        }, 2000);
    });
};

// ── Render a citation card into a container element ───────────────────────────
// keys: array of CITATIONS keys to include, e.g. ["splace","syngenes","mafft"]
window.renderCitationsCard = function (container, keys) {
    if (!container) return;

    const items = keys.map(k => window.CITATIONS[k]).filter(Boolean);
    if (!items.length) return;

    // Build rows — citation text stored in data-cite to avoid attribute quoting issues
    const rows = items.map(c => {
        const doiUrl = `https://doi.org/${c.doi}`;
        const volInfo = [
            c.volume,
            c.issue ? `(${c.issue})` : "",
            c.pages ? `:${c.pages}` : "",
        ].join("");
        return `
        <div class="bg-white rounded-lg border border-splace-blue-100 p-3">
            <div class="flex items-center justify-between gap-3">
                <div class="min-w-0">
                    <span class="inline-block text-sm font-bold text-splace-blue-700 mb-1">${c.name}</span>
                    <p class="text-sm text-gray-700 leading-relaxed">
                        ${c.authors} (${c.year})
                        <em>${c.title}</em>.
                        ${c.journal}${volInfo ? ` ${volInfo}` : ""}.
                        DOI: <a href="${doiUrl}" target="_blank" rel="noopener"
                                class="text-splace-blue-600 hover:text-splace-blue-800 hover:underline break-all">${c.doi}</a>
                    </p>
                </div>
                <button
                    class="flex-shrink-0 inline-flex items-center gap-1.5 px-2.5 py-1.5 rounded-lg bg-splace-blue-50 hover:bg-splace-blue-100 border border-splace-blue-200 text-splace-blue-700 text-xs font-semibold transition-colors self-center"
                    title="Copy citation"
                    onclick="copyCitation(this)">
                    <i class="fa-regular fa-copy"></i>Copy
                </button>
            </div>
        </div>`;
    }).join("\n");

    // Inject HTML, then set data-cite on each button (avoids any escaping in attributes)
    container.innerHTML = `
    <div class="mt-4 border border-splace-blue-200 rounded-xl bg-splace-blue-50 p-4">
        <div class="flex items-center gap-2 mb-3">
            <i class="fa-solid fa-graduation-cap text-splace-blue-600 text-lg"></i>
            <h3 class="text-base font-semibold text-splace-blue-800">How to Cite</h3>
            <span class="ml-auto text-xs text-splace-blue-500">If you use these results, please cite the tools below</span>
        </div>
        <div class="space-y-2">
            ${rows}
        </div>
    </div>`;

    // Attach citation text via JS property — no HTML-escaping issues
    container.querySelectorAll("button[onclick='copyCitation(this)']").forEach((btn, i) => {
        btn.dataset.cite = items[i].cite;
    });
};

// ── Auto-populate containers on page load ─────────────────────────────────────
(function () {
    function populate() {
        // Web citation block (inside Export FASTA section)
        const webEl = document.getElementById("citationsWeb");
        if (webEl) window.renderCitationsCard(webEl, ["splace", "syngenes", "datafishing"]);

        // Desktop citation block (inside IQ-TREE results section)
        const desktopEl = document.getElementById("citationsDesktop");
        if (desktopEl) window.renderCitationsCard(desktopEl, ["splace", "syngenes", "datafishing", "mafft", "trimal", "iqtree"]);
    }
    if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", populate);
    else populate();
})();
