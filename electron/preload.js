const { contextBridge, ipcRenderer } = require("electron");

let spawnSync = null;
try {
    ({ spawnSync } = require("child_process"));
} catch {
    spawnSync = null;
}

let fs = null;
let path = null;
let pathToFileURL = null;
try {
    fs = require("fs");
    path = require("path");
    ({ pathToFileURL } = require("url"));
} catch {
    fs = null;
    path = null;
    pathToFileURL = null;
}

function runCommandCaptureLocal(cmd, args, options = {}) {
    if (!spawnSync) {
        return {
            status: 1,
            error: "child_process is unavailable in the preload context.",
            stdout: "",
            stderr: "",
        };
    }
    const result = spawnSync(cmd, args, {
        encoding: "utf8",
        timeout: options.timeout || 120000,
        shell: false,
        windowsHide: true,
    });
    return {
        status: result.status,
        error: result.error ? result.error.message : null,
        stdout: result.stdout || "",
        stderr: result.stderr || "",
    };
}

function pythonCandidatesLocal() {
    if (process.platform === "win32") {
        return [
            { cmd: "py", args: ["-3"] },
            { cmd: "python", args: [] },
            { cmd: "python3", args: [] },
        ];
    }
    return [
        { cmd: "python3", args: [] },
        { cmd: "python", args: [] },
    ];
}

function findPythonLocal() {
    for (const candidate of pythonCandidatesLocal()) {
        const probe = runCommandCaptureLocal(candidate.cmd, [...candidate.args, "--version"], { timeout: 10000 });
        const version = `${probe.stdout || ""}${probe.stderr || ""}`.trim();
        if (!probe.error && probe.status === 0 && /Python\s+\d+/i.test(version)) {
            return { cmd: candidate.cmd, prefixArgs: candidate.args, version };
        }
    }
    return null;
}

function runPythonLocal(py, args, options = {}) {
    return runCommandCaptureLocal(py.cmd, [...py.prefixArgs, ...args], options);
}

function checkPythonEnvironmentLocal(ipcErrorMessage = "") {
    const py = findPythonLocal();
    if (!py) {
        return {
            pythonFound: false,
            clipkitFound: false,
            ipcFallback: true,
            ipcErrorMessage,
            error: "Python was not found in PATH. On Windows, check that py -3 or python works in a new terminal.",
        };
    }
    const clipkitProbe = runPythonLocal(py, ["-c", "import clipkit, sys; print(getattr(clipkit, '__version__', 'installed'))"], { timeout: 10000 });
    return {
        pythonFound: true,
        clipkitFound: clipkitProbe.status === 0,
        ipcFallback: true,
        ipcErrorMessage,
        executable: [py.cmd, ...py.prefixArgs].join(" "),
        version: py.version.replace(/^Python\s*/i, ""),
        clipkitVersion: clipkitProbe.status === 0 ? (clipkitProbe.stdout || "installed").trim() : null,
        clipkitError: clipkitProbe.status === 0 ? null : (clipkitProbe.stderr || clipkitProbe.error || "ClipKIT import failed"),
    };
}

function installPythonPackagesLocal({ packages } = {}, ipcErrorMessage = "") {
    const py = findPythonLocal();
    if (!py) {
        return {
            success: false,
            ipcFallback: true,
            ipcErrorMessage,
            error: "Python was not found in PATH. On Windows, check that py -3 or python works in a new terminal.",
        };
    }
    const safePackages = (Array.isArray(packages) && packages.length ? packages : ["clipkit"])
        .map(String)
        .filter((pkg) => /^[A-Za-z0-9_.-]+$/.test(pkg));
    if (!safePackages.length) {
        return { success: false, ipcFallback: true, ipcErrorMessage, error: "No valid packages were requested." };
    }
    const result = runPythonLocal(py, ["-m", "pip", "install", "--user", "--upgrade", ...safePackages], { timeout: 600000 });
    return {
        success: result.status === 0,
        ipcFallback: true,
        ipcErrorMessage,
        output: `${result.stdout || ""}${result.stderr || ""}`.trim(),
        error: result.status === 0 ? null : (result.stderr || result.error || `pip exited with status ${result.status}`),
    };
}

async function invokeWithLocalFallback(primaryChannel, fallbackChannel, params, localFallback) {
    try {
        return await ipcRenderer.invoke(primaryChannel, params);
    } catch (primaryError) {
        try {
            return await ipcRenderer.invoke(fallbackChannel, params);
        } catch (fallbackError) {
            return localFallback(primaryError?.message || fallbackError?.message || String(primaryError || fallbackError || "IPC unavailable"));
        }
    }
}


function highchartsAssetCandidatesLocal(relativePath) {
    if (!fs || !path) return [];
    return [
        path.join(__dirname, "node_modules", "highcharts", relativePath),
        path.join(__dirname, "..", "electron", "node_modules", "highcharts", relativePath),
        path.join(process.cwd ? process.cwd() : __dirname, "node_modules", "highcharts", relativePath),
    ];
}

function findHighchartsAssetLocal(relativePath) {
    for (const candidate of highchartsAssetCandidatesLocal(relativePath)) {
        try {
            if (candidate && fs.existsSync(candidate)) return candidate;
        } catch {
            // Try next location.
        }
    }
    throw new Error(`Highcharts asset not found locally: ${relativePath}`);
}

function getHighchartsAssetUrlsLocal(ipcErrorMessage = "") {
    if (!pathToFileURL || !path) throw new Error(ipcErrorMessage || "Highcharts IPC and local file URL fallback are unavailable.");
    return {
        ipcFallback: true,
        ipcErrorMessage,
        highcharts: pathToFileURL(findHighchartsAssetLocal("highcharts.js")).toString(),
        heatmap: pathToFileURL(findHighchartsAssetLocal(path.join("modules", "heatmap.js"))).toString(),
        accessibility: pathToFileURL(findHighchartsAssetLocal(path.join("modules", "accessibility.js"))).toString(),
    };
}

function loadHighchartsAssetsLocal(ipcErrorMessage = "") {
    if (!fs || !path) throw new Error(ipcErrorMessage || "Highcharts IPC and local read fallback are unavailable.");
    return {
        ipcFallback: true,
        ipcErrorMessage,
        highcharts: fs.readFileSync(findHighchartsAssetLocal("highcharts.js"), "utf8"),
        heatmap: fs.readFileSync(findHighchartsAssetLocal(path.join("modules", "heatmap.js")), "utf8"),
        accessibility: fs.readFileSync(findHighchartsAssetLocal(path.join("modules", "accessibility.js")), "utf8"),
    };
}

async function invokeHighchartsUrls() {
    try {
        return await ipcRenderer.invoke("highcharts-asset-urls");
    } catch (primaryError) {
        try {
            return await ipcRenderer.invoke("splace-highcharts-asset-urls");
        } catch (fallbackError) {
            return getHighchartsAssetUrlsLocal(primaryError?.message || fallbackError?.message || String(primaryError || fallbackError || "IPC unavailable"));
        }
    }
}

async function invokeHighchartsAssets() {
    try {
        return await ipcRenderer.invoke("load-highcharts-assets");
    } catch (primaryError) {
        try {
            return await ipcRenderer.invoke("splace-load-highcharts-assets");
        } catch (fallbackError) {
            return loadHighchartsAssetsLocal(primaryError?.message || fallbackError?.message || String(primaryError || fallbackError || "IPC unavailable"));
        }
    }
}

contextBridge.exposeInMainWorld("isElectron", true);

contextBridge.exposeInMainWorld("electronAPI", {
    getCpuCount: () => ipcRenderer.invoke("get-cpu-count"),
    selectGenbankInputs: () => ipcRenderer.invoke("select-genbank-inputs"),
    getHighchartsAssetUrls: () => invokeHighchartsUrls(),
    loadHighchartsAssets: () => invokeHighchartsAssets(),

    checkPythonEnvironment: async () => invokeWithLocalFallback(
        "python-status",
        "splace-python-status",
        undefined,
        (message) => checkPythonEnvironmentLocal(message)
    ),

    installPythonPackages: async (params) => invokeWithLocalFallback(
        "python-install-packages",
        "splace-python-install-packages",
        params,
        (message) => installPythonPackagesLocal(params, message)
    ),

    runAnalysis: (params) => ipcRenderer.send("run-analysis", params),
    runTrimal: (params) => ipcRenderer.send("run-trimal", params),
    saveCodonOutputs: (params) => ipcRenderer.invoke("save-codon-outputs", params),

    onAnalysisProgress: (callback) => {
        ipcRenderer.on("analysis-progress", (_event, data) => callback(data));
    },
    onAnalysisDone: (callback) => {
        ipcRenderer.on("analysis-done", (_event, result) => callback(result));
    },

    runIqtree: (params) => ipcRenderer.send("run-iqtree", params),
    onIqtreeProgress: (callback) => {
        ipcRenderer.on("iqtree-progress", (_event, data) => callback(data));
    },
    onIqtreeDone: (callback) => {
        ipcRenderer.on("iqtree-done", (_event, result) => callback(result));
    },

    saveZip: (params) => ipcRenderer.invoke("save-zip", params),
    translateFasta: (params) => ipcRenderer.invoke("translate-fasta", params),
});
