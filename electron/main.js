const { app, BrowserWindow, ipcMain, dialog, shell } = require("electron");
const path = require("path");
const fs = require("fs");
const os = require("os");
const { spawn, spawnSync } = require("child_process");
const { pathToFileURL } = require("url");


// -----------------------------------------------------------------------
// Desktop-only Highcharts assets
// -----------------------------------------------------------------------
function highchartsAssetCandidates(relativePath) {
    return [
        path.join(__dirname, "node_modules", "highcharts", relativePath),
        path.join(__dirname, "..", "electron", "node_modules", "highcharts", relativePath),
        path.join(process.resourcesPath || "", "node_modules", "highcharts", relativePath),
        path.join(process.resourcesPath || "", "app.asar.unpacked", "node_modules", "highcharts", relativePath),
        path.join(process.resourcesPath || "", "app.asar", "node_modules", "highcharts", relativePath),
    ];
}

function findHighchartsAsset(relativePath) {
    for (const candidate of highchartsAssetCandidates(relativePath)) {
        try {
            if (candidate && fs.existsSync(candidate)) return candidate;
        } catch {
            // Try the next location.
        }
    }
    throw new Error(`Highcharts asset not found: ${relativePath}`);
}

function readHighchartsAsset(relativePath) {
    return fs.readFileSync(findHighchartsAsset(relativePath), "utf8");
}

function getHighchartsAssetUrls() {
    return {
        highcharts: pathToFileURL(findHighchartsAsset("highcharts.js")).toString(),
        heatmap: pathToFileURL(findHighchartsAsset(path.join("modules", "heatmap.js"))).toString(),
        accessibility: pathToFileURL(findHighchartsAsset(path.join("modules", "accessibility.js"))).toString(),
    };
}

function getHighchartsAssetCode() {
    return {
        highcharts: readHighchartsAsset("highcharts.js"),
        heatmap: readHighchartsAsset(path.join("modules", "heatmap.js")),
        accessibility: readHighchartsAsset(path.join("modules", "accessibility.js")),
    };
}

ipcMain.handle("highcharts-asset-urls", async () => getHighchartsAssetUrls());
ipcMain.handle("splace-highcharts-asset-urls", async () => getHighchartsAssetUrls());
ipcMain.handle("load-highcharts-assets", async () => getHighchartsAssetCode());
ipcMain.handle("splace-load-highcharts-assets", async () => getHighchartsAssetCode());

// -----------------------------------------------------------------------
// Desktop-only Python tools for optional ClipKIT trimming
// -----------------------------------------------------------------------
function runCommandCapture(cmd, args, options = {}) {
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

function pythonCandidates() {
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

function findPython() {
    for (const candidate of pythonCandidates()) {
        const probe = runCommandCapture(candidate.cmd, [...candidate.args, "--version"], { timeout: 10000 });
        const version = `${probe.stdout || ""}${probe.stderr || ""}`.trim();
        if (!probe.error && probe.status === 0 && /Python\s+\d+/i.test(version)) {
            return { cmd: candidate.cmd, prefixArgs: candidate.args, version };
        }
    }
    return null;
}

function runPython(py, args, options = {}) {
    return runCommandCapture(py.cmd, [...py.prefixArgs, ...args], options);
}

function checkPythonEnvironment() {
    const py = findPython();
    if (!py) {
        return { pythonFound: false, clipkitFound: false, error: "Python was not found in PATH." };
    }
    const clipkitProbe = runPython(py, ["-c", "import clipkit, sys; print(getattr(clipkit, '__version__', 'installed'))"], { timeout: 10000 });
    return {
        pythonFound: true,
        executable: [py.cmd, ...py.prefixArgs].join(" "),
        version: py.version.replace(/^Python\s*/i, ""),
        clipkitFound: clipkitProbe.status === 0,
        clipkitVersion: clipkitProbe.status === 0 ? (clipkitProbe.stdout || "installed").trim() : null,
        clipkitError: clipkitProbe.status === 0 ? null : (clipkitProbe.stderr || clipkitProbe.error || "ClipKIT import failed"),
    };
}

ipcMain.handle("python-status", async () => checkPythonEnvironment());
ipcMain.handle("splace-python-status", async () => checkPythonEnvironment());

ipcMain.handle("splace-python-install-packages", async (_event, params = {}) => installPythonPackages(params));

async function installPythonPackages({ packages } = {}) {
    const py = findPython();
    if (!py) return { success: false, error: "Python was not found in PATH." };
    const safePackages = (Array.isArray(packages) && packages.length ? packages : ["clipkit"])
        .map(String)
        .filter((pkg) => /^[A-Za-z0-9_.-]+$/.test(pkg));
    if (!safePackages.length) return { success: false, error: "No valid packages were requested." };

    const result = runPython(py, ["-m", "pip", "install", "--user", "--upgrade", ...safePackages], { timeout: 600000 });
    return {
        success: result.status === 0,
        output: `${result.stdout || ""}${result.stderr || ""}`.trim(),
        error: result.status === 0 ? null : (result.stderr || result.error || `pip exited with status ${result.status}`),
    };
}

ipcMain.handle("python-install-packages", async (_event, params = {}) => installPythonPackages(params));


// -----------------------------------------------------------------------
// Binary resolution: WSL on Windows, system PATH otherwise
// -----------------------------------------------------------------------

const IS_WIN = process.platform === "win32";
const GENBANK_FILE_RE = /\.(gb|gbk|genbank)$/i;

// Resolve a bundled binary or fall back to system PATH
function resolveBin(name) {
    const platform = IS_WIN ? "win" : process.platform === "darwin" ? "mac" : "linux";
    const binDir = app.isPackaged
        ? path.join(process.resourcesPath, "bin", platform)
        : path.join(__dirname, "bin", platform);

    if (IS_WIN) {
        // Try .bat first (e.g. mafft.bat), then .exe
        for (const ext of [".bat", ".exe"]) {
            const local = path.join(binDir, name + ext);
            if (fs.existsSync(local)) return local;
        }
    } else {
        const local = path.join(binDir, name);
        if (fs.existsSync(local)) {
            try { fs.chmodSync(local, 0o755); } catch { }
            return local;
        }
    }
    return name;
}

function collectGenbankFilesFromPath(targetPath, files = []) {
    let stat;
    try {
        stat = fs.statSync(targetPath);
    } catch {
        return files;
    }

    if (stat.isDirectory()) {
        let entries = [];
        try {
            entries = fs.readdirSync(targetPath);
        } catch {
            return files;
        }
        for (const entry of entries) {
            collectGenbankFilesFromPath(path.join(targetPath, entry), files);
        }
        return files;
    }

    if (!stat.isFile() || !GENBANK_FILE_RE.test(targetPath)) {
        return files;
    }

    try {
        files.push({
            name: path.basename(targetPath),
            path: targetPath,
            text: fs.readFileSync(targetPath, "utf8"),
            size: stat.size,
            lastModified: stat.mtimeMs,
        });
    } catch {
        return files;
    }

    return files;
}

// -----------------------------------------------------------------------
// App window
// -----------------------------------------------------------------------
function createWindow() {
    const win = new BrowserWindow({
        width: 1280,
        height: 900,
        minWidth: 1280,
        minHeight: 720,
        title: "SPLACE",
        webPreferences: {
            preload: path.join(__dirname, "preload.js"),
            contextIsolation: true,
            nodeIntegration: false,
            sandbox: false,
        },
    });
    const indexHtml = app.isPackaged
        ? path.join(process.resourcesPath, "docs", "index.html")
        : path.join(__dirname, "..", "docs", "index.html");

    win.webContents.setWindowOpenHandler(({ url }) => {
        if (/^https?:\/\//i.test(url)) {
            shell.openExternal(url);
            return { action: "deny" };
        }
        return { action: "allow" };
    });

    win.webContents.on("will-navigate", (event, url) => {
        const currentUrl = win.webContents.getURL();
        if (url === currentUrl) return;

        event.preventDefault();
        if (/^https?:\/\//i.test(url)) {
            shell.openExternal(url);
        }
    });

    win.loadFile(indexHtml);
}

app.whenReady().then(() => {
    createWindow();
    app.on("activate", () => {
        if (BrowserWindow.getAllWindows().length === 0) createWindow();
    });
});

app.on("window-all-closed", () => {
    if (process.platform !== "darwin") app.quit();
});

// -----------------------------------------------------------------------
// IPC: CPU count
// -----------------------------------------------------------------------
ipcMain.handle("get-cpu-count", () => os.cpus().length);

ipcMain.handle("select-genbank-inputs", async () => {
    const focusedWindow = BrowserWindow.getFocusedWindow();
    const result = await dialog.showOpenDialog(focusedWindow, {
        title: "Select GenBank files or folders",
        properties: ["openFile", "openDirectory", "multiSelections"],
        filters: [
            { name: "GenBank files", extensions: ["gb", "gbk", "genbank"] },
            { name: "All files", extensions: ["*"] },
        ],
    });

    if (result.canceled || !result.filePaths?.length) {
        return [];
    }

    const collected = [];
    for (const selectedPath of result.filePaths) {
        collectGenbankFilesFromPath(selectedPath, collected);
    }

    return collected;
});

// -----------------------------------------------------------------------
// IPC: run MAFFT alignment
// { files: {geneName: fastaContent}, params: string[], threads: number }
// -----------------------------------------------------------------------
ipcMain.on("run-analysis", async (event, { files, params, threads }) => {
    const sender = event.sender;
    const send = (channel, data) => { if (!sender.isDestroyed()) sender.send(channel, data); };

    const outputDir = path.join(os.homedir(), "Documents", "SPLACE", "sequences");
    try { fs.rmSync(outputDir, { recursive: true, force: true }); } catch { }
    fs.mkdirSync(outputDir, { recursive: true });

    const tmpDir = fs.mkdtempSync(path.join(os.tmpdir(), "splace-"));
    const markers = Object.keys(files);
    let done = 0, aligned = 0;
    const markerResults = new Array(markers.length).fill(null);

    // Parallel pool — jobs limited to leave 2 threads free for the OS
    const concurrency = Math.max(1, Math.floor((os.cpus().length - 2) / Math.max(1, threads)));
    let nextIdx = 0;

    async function worker() {
        while (true) {
            const idx = nextIdx++;
            if (idx >= markers.length) break;
            const name = markers[idx];
            const inputFile = path.join(tmpDir, `${name}.fasta`);
            const outputFile = path.join(outputDir, `${name}_aligned.fasta`);
            const rawContent = files[name];

            fs.writeFileSync(inputFile, rawContent, "utf8");

            send("analysis-progress", {
                marker: name, status: "running",
                message: `[${name}] Running MAFFT…`,
                done, total: markers.length,
            });

            const rawSeqs = parseFastaStats(rawContent);
            const { success, alignedContent, error } = await runMafft(
                inputFile, outputFile, params, threads,
                (msg) => send("analysis-progress", { message: `[${name}] ${msg}`, done, total: markers.length })
            );

            done++;
            if (success) aligned++;
            const alignedStats = success ? parseFastaStats(alignedContent) : null;

            send("analysis-progress", {
                marker: name,
                status: success ? "done" : "error",
                info: success ? `${alignedStats.numSeqs} seq · ${alignedStats.length} bp` : error,
                message: success ? `[${name}] ✓ Saved → ${outputFile}` : `[${name}] ✗ ${error}`,
                done, total: markers.length,
            });

            markerResults[idx] = {
                marker: name, rawContent,
                alignedContent: success ? alignedContent : null,
                rawStats: rawSeqs, alignedStats,
                trimmedContent: null, trimmedStats: null,
                outputFile: success ? outputFile : null,
            };
        }
    }

    await Promise.all(Array.from({ length: concurrency }, () => worker()));
    try { fs.rmSync(tmpDir, { recursive: true, force: true }); } catch { }

    const finalResults = markerResults.filter(Boolean);
    send("analysis-done", {
        phase: "mafft",
        success: aligned > 0,
        aligned, total: markers.length,
        outputDir,
        markerResults: finalResults,
    });
    _lastMarkerResults = finalResults;
});

// -----------------------------------------------------------------------
// IPC: run trimAl
// { markers: string[], params: string[] }
// (uses already-aligned files from previous run stored in markerResults)
// -----------------------------------------------------------------------
let _lastMarkerResults = [];

ipcMain.on("run-trimal", async (event, { markers, params, codonMode, cdsFastas, alignmentMode, trimTool = "auto", clipkitMode = "gappy", codonTrimAction = "trim-codons" }) => {
    const sender = event.sender;
    const send = (channel, data) => { if (!sender.isDestroyed()) sender.send(channel, data); };

    let done = 0, trimmed = 0;
    const failedMarkers = [];
    const updatedResults = _lastMarkerResults.map(r => ({ ...r }));
    const concurrency = Math.max(1, os.cpus().length - 2);
    let nextIdx = 0;

    async function worker() {
        while (true) {
            const idx = nextIdx++;
            if (idx >= markers.length) break;
            const name = markers[idx];
            const mr = updatedResults.find(r => r.marker === name);
            const proteinInputFile = mr?.proteinAlignedFile || mr?.outputFile;
            if (!mr || !mr.alignedContent || !(codonMode ? proteinInputFile : mr.outputFile)) {
                done++;
                send("analysis-progress", { marker: name, status: "error", message: `[${name}] No aligned file found`, done, total: markers.length });
                continue;
            }

            const inputFile = codonMode ? proteinInputFile : mr.outputFile;
            const trimmedFile = codonMode
                ? path.join(path.dirname(mr.outputFile || inputFile), `${name}_codon_trimmed.fasta`)
                : inputFile.replace(/_aligned\.fasta$/, "_trimmed.fasta");

            // ── Codon-aware trimming (back-translation) ─────────────────────────
            if (codonMode && cdsFastas && cdsFastas[name]) {
                send("analysis-progress", {
                    marker: name, status: "running",
                    message: `[${name}] Running codon-aware trimming…`,
                    done, total: markers.length,
                });

                const dnaCdsContent = cdsFastas[name];
                let dnaResult = null;
                const logicalTrimMetadata = {
                    inputAlignment: `${name}_protein_aligned.fasta`,
                    backtransInput: `${name}_cds_validated.fasta`,
                    outputAlignment: `${name}_codon_trimmed.fasta`,
                    outputType: "codon-preserving DNA",
                };
                const trimalMetadata = {
                    trimmingTool: "trimal",
                    trimmingMode: "protein-guided-backtrans",
                    ...logicalTrimMetadata,
                };
                const fallbackMetadata = {
                    trimmingTool: "splace-internal",
                    trimmingMode: "protein-guided-codon-block-fallback",
                    ...logicalTrimMetadata,
                };
                const failCodonTrim = (issues, metadata, methodLabel) => {
                    const issueList = (Array.isArray(issues) ? issues : [issues]).filter(Boolean);
                    const summary = issueList.join("; ");
                    mr.trimmedContent = null;
                    mr.trimmedStats = null;
                    mr.trimmedFile = null;
                    mr.trimError = summary;
                    mr.trimValidation = {
                        valid: false,
                        issues: issueList,
                        equalLengths: false,
                        allLengthsMultipleOf3: false,
                        gapBlocksAreTriplets: false,
                        headersMatch: false,
                        dnaCompatible: false,
                        frameDisrupted: true,
                    };
                    mr.trimMetadata = {
                        ...metadata,
                        codonIntegrityValidated: false,
                    };
                    mr.trimMethod = metadata?.trimmingTool === "trimal" ? "trimal-backtrans" : "js-fallback";
                    try { fs.rmSync(trimmedFile, { force: true }); } catch { }
                    issueList.forEach((issue) => send("analysis-progress", { message: `[${name}] VALIDATION ERROR: ${issue}` }));
                    done++;
                    failedMarkers.push({ marker: name, issues: issueList, method: methodLabel || metadata?.trimmingMode || "codon-aware trimming" });
                    send("analysis-progress", {
                        marker: name,
                        status: "error",
                        info: summary,
                        message: `[${name}] ✗ Codon-aware trimming failed`,
                        done,
                        total: markers.length,
                    });
                };

                const requestedTrimTool = String(trimTool || "auto").toLowerCase();
                const requestedCodonAction = String(codonTrimAction || "trim-codons").toLowerCase();
                const positionMatch = requestedCodonAction.match(/^remove-pos-([123])$/);
                if (positionMatch) {
                    const removedPosition = Number(positionMatch[1]);
                    try {
                        send("analysis-progress", { message: `[${name}] Back-translating protein alignment before removing codon position ${removedPosition}…` });
                        const proteinAlignmentContent = mr.alignedContent || fs.readFileSync(inputFile, "utf8");
                        const built = buildCodonAlignedDnaFromProteinAlignment(proteinAlignmentContent, dnaCdsContent);
                        if (!built.dnaFasta) throw new Error("Back-translation produced no DNA output");
                        const validation = validateCodonPreservingTrimmedAlignment(built.dnaFasta, dnaCdsContent);
                        if (!validation.valid) throw new Error(validation.issues.join("; "));
                        const filteredFasta = removeCodonPositionFromFasta(built.dnaFasta, removedPosition);
                        fs.writeFileSync(trimmedFile, filteredFasta);
                        const filteredStats = parseFastaStats(filteredFasta);
                        mr.trimmedContent = filteredFasta;
                        mr.trimmedFile = trimmedFile;
                        mr.trimmedStats = filteredStats;
                        mr.trimMethod = `remove-codon-position-${removedPosition}`;
                        mr.trimValidation = {
                            valid: true,
                            issues: [],
                            equalLengths: true,
                            allLengthsMultipleOf3: false,
                            gapBlocksAreTriplets: true,
                            headersMatch: true,
                            dnaCompatible: true,
                            frameDisrupted: false,
                            positionFiltered: true,
                            removedCodonPosition: removedPosition,
                        };
                        mr.trimMetadata = {
                            trimmingTool: "splace-internal",
                            trimmingMode: `remove-codon-position-${removedPosition}`,
                            removedCodonPosition: removedPosition,
                            positionFilter: true,
                            inputAlignment: `${name}_protein_aligned.fasta`,
                            backtransInput: `${name}_cds_validated.fasta`,
                            outputAlignment: `${name}_codon_position_${removedPosition}_removed.fasta`,
                            outputType: "codon-position-filtered DNA",
                            codonIntegrityValidated: true,
                        };
                        done++;
                        trimmed++;
                        send("analysis-progress", {
                            marker: name,
                            status: "done",
                            info: `${filteredStats.numSeqs} seq · ${filteredStats.length} bp (position ${removedPosition} removed)`,
                            message: `[${name}] ✓ Removed codon position ${removedPosition} → ${filteredStats.numSeqs} seq · ${filteredStats.length} bp`,
                            done,
                            total: markers.length,
                        });
                    } catch (error) {
                        failCodonTrim(error?.message || String(error), {
                            trimmingTool: "splace-internal",
                            trimmingMode: `remove-codon-position-${removedPosition}`,
                            removedCodonPosition,
                            positionFilter: true,
                            ...logicalTrimMetadata,
                            outputType: "codon-position-filtered DNA",
                            codonIntegrityValidated: false,
                        }, `remove codon position ${removedPosition}`);
                    }
                    continue;
                }
                const preferClipkit = requestedTrimTool === "clipkit";
                const allowTrimal = requestedTrimTool !== "clipkit";
                const allowClipkitFallback = requestedTrimTool === "auto" || requestedTrimTool === "clipkit";

                async function runClipkitBacktranslation() {
                    const tmpProtOut = path.join(os.tmpdir(), `splace-clipkit-${Date.now()}-${idx}.fasta`);
                    try {
                        const clipRes = await runClipkitProteinTrim(inputFile, tmpProtOut, clipkitMode, (msg) => send("analysis-progress", { message: `[${name}] ${msg}` }));
                        if (!clipRes.success) return { success: false, error: clipRes.error };
                        const protSeqs = parseFastaArr(clipRes.trimmedContent);
                        const dnaSeqs = parseFastaArr(dnaCdsContent);
                        const dnaMap = {};
                        for (const s of dnaSeqs) dnaMap[s.name] = s.seq.replace(/-/g, "");

                        let dnaFasta = "";
                        let skipped = 0;
                        for (const ps of protSeqs) {
                            const dna = dnaMap[ps.name];
                            if (!dna) { skipped++; continue; }
                            dnaFasta += `>${ps.name}\n${backTranslateJS(ps.seq, dna)}\n`;
                        }
                        if (!dnaFasta) return { success: false, error: "ClipKIT back-translation produced no DNA output" };
                        const validation = validateCodonPreservingTrimmedAlignment(dnaFasta, dnaCdsContent);
                        if (!validation.valid) return { success: false, error: validation.issues };
                        if (skipped > 0) send("analysis-progress", { message: `[${name}] WARNING: ${skipped} sequences were unmatched after ClipKIT` });
                        return {
                            success: true,
                            content: dnaFasta,
                            validation,
                            method: "clipkit-backtrans",
                            metadata: {
                                trimmingTool: "clipkit",
                                trimmingMode: `clipkit-${clipkitMode || "gappy"}-protein-guided-backtrans`,
                                ...logicalTrimMetadata,
                                codonIntegrityValidated: true,
                            },
                        };
                    } finally {
                        try { fs.unlinkSync(tmpProtOut); } catch { }
                    }
                }

                if (preferClipkit) {
                    send("analysis-progress", { message: `[${name}] Using ClipKIT for protein-guided codon trimming…` });
                    const clipResult = await runClipkitBacktranslation();
                    if (clipResult.success) {
                        dnaResult = clipResult;
                        send("analysis-progress", { message: `[${name}] ✓ ClipKIT produced validated codon-preserving DNA` });
                    } else {
                        failCodonTrim(clipResult.error || "ClipKIT failed", { trimmingTool: "clipkit", trimmingMode: `clipkit-${clipkitMode || "gappy"}`, ...logicalTrimMetadata, codonIntegrityValidated: false }, "ClipKIT");
                        continue;
                    }
                }

                // Try trimAl -backtrans
                const trimalBin = resolveBin("trimal");
                if (allowTrimal && !dnaResult) {
                    const tmpCdsFile = (!mr.rawFile || !fs.existsSync(mr.rawFile))
                        ? path.join(os.tmpdir(), `splace-cds-${Date.now()}-${idx}.fasta`)
                        : null;
                    const tmpOutFile = path.join(os.tmpdir(), `splace-bt-${Date.now()}-${idx}.fasta`);
                    const backtransInputFile = mr.rawFile && fs.existsSync(mr.rawFile) ? mr.rawFile : tmpCdsFile;
                    try {
                        if (tmpCdsFile) {
                            fs.writeFileSync(tmpCdsFile, dnaCdsContent, "utf8");
                        }
                        const btArgs = ["-in", inputFile, "-backtrans", backtransInputFile, "-out", tmpOutFile, "-fasta", ...params];
                        const btRes = await new Promise((resolve) => {
                            const proc = spawn(trimalBin, btArgs, { shell: false });
                            let stderr = "";
                            proc.stderr.on("data", d => {
                                stderr += d.toString();
                                const m = d.toString().trim();
                                if (m) send("analysis-progress", { message: `[${name}] ${m}` });
                            });
                            proc.on("close", code => {
                                if (code === 0 && fs.existsSync(tmpOutFile)) {
                                    resolve({ success: true, content: fs.readFileSync(tmpOutFile, "utf8") });
                                } else {
                                    resolve({ success: false, error: stderr.trim() || `Exit code ${code}`, unavailable: false });
                                }
                            });
                            proc.on("error", err => resolve({ success: false, error: err.message, unavailable: isTrimAlUnavailableError(err.message) }));
                        });
                        if (btRes.success) {
                            const validation = validateCodonPreservingTrimmedAlignment(btRes.content, dnaCdsContent);
                            if (!validation.valid) {
                                if (requestedTrimTool === "auto") {
                                    send("analysis-progress", { message: `[${name}] trimAl -backtrans validation failed, trying ClipKIT fallback` });
                                } else {
                                    failCodonTrim(validation.issues, trimalMetadata, "trimAl -backtrans");
                                    continue;
                                }
                            } else {
                                dnaResult = {
                                    content: btRes.content,
                                    method: "trimal-backtrans",
                                    validation,
                                    metadata: {
                                        ...trimalMetadata,
                                        codonIntegrityValidated: true,
                                    },
                                };
                                send("analysis-progress", { message: `[${name}] ✓ trimAl -backtrans produced validated codon-preserving DNA` });
                            }
                        } else if (btRes.unavailable) {
                            send("analysis-progress", { message: `[${name}] trimAl is unavailable (${btRes.error}), using approximate JavaScript fallback` });
                        } else {
                            if (requestedTrimTool === "auto") {
                                send("analysis-progress", { message: `[${name}] trimAl -backtrans failed (${btRes.error}), trying ClipKIT fallback` });
                            } else {
                                failCodonTrim(`trimAl -backtrans failed: ${btRes.error}`, trimalMetadata, "trimAl -backtrans");
                                continue;
                            }
                        }
                    } catch (e) {
                        if (isTrimAlUnavailableError(e.message)) {
                            send("analysis-progress", { message: `[${name}] trimAl is unavailable (${e.message}), using approximate JavaScript fallback` });
                        } else {
                            if (requestedTrimTool === "auto") {
                                send("analysis-progress", { message: `[${name}] trimAl -backtrans error (${e.message}), trying ClipKIT fallback` });
                            } else {
                                failCodonTrim(`trimAl -backtrans error: ${e.message}`, trimalMetadata, "trimAl -backtrans");
                                continue;
                            }
                        }
                    } finally {
                        try { if (tmpCdsFile) fs.unlinkSync(tmpCdsFile); } catch { }
                        try { fs.unlinkSync(tmpOutFile); } catch { }
                    }
                }

                if (!dnaResult && allowClipkitFallback) {
                    send("analysis-progress", { message: `[${name}] Trying ClipKIT as codon-aware trimming fallback…` });
                    const clipResult = await runClipkitBacktranslation();
                    if (clipResult.success) {
                        dnaResult = clipResult;
                        send("analysis-progress", { message: `[${name}] ✓ ClipKIT fallback produced validated codon-preserving DNA` });
                    } else {
                        send("analysis-progress", { message: `[${name}] ClipKIT fallback unavailable or failed: ${Array.isArray(clipResult.error) ? clipResult.error.join("; ") : clipResult.error}` });
                    }
                }

                if (!dnaResult) {
                    // JS fallback: trim protein alignment, then back-translate to DNA
                    send("analysis-progress", { message: `[${name}] Using approximate JavaScript codon-aware fallback` });
                    const protContent = fs.existsSync(inputFile)
                        ? fs.readFileSync(inputFile, "utf8")
                        : (mr.proteinAlignedContent || "");
                    const trimResult = trimAlignmentJS(protContent, params);
                    if (!trimResult.success) {
                        failCodonTrim(`SPLACE internal fallback failed: ${trimResult.error}`, fallbackMetadata, "SPLACE internal fallback");
                        continue;
                    }
                    const protSeqs = parseFastaArr(trimResult.trimmedContent);
                    const dnaSeqs = parseFastaArr(dnaCdsContent);
                    const dnaMap = {};
                    for (const s of dnaSeqs) dnaMap[s.name] = s.seq.replace(/-/g, "");

                    let dnaFasta = "";
                    let skipped = 0;
                    for (const ps of protSeqs) {
                        const dna = dnaMap[ps.name];
                        if (!dna) {
                            send("analysis-progress", { message: `[${name}] WARNING: no DNA found for "${ps.name}"` });
                            skipped++;
                            continue;
                        }
                        dnaFasta += `>${ps.name}\n${backTranslateJS(ps.seq, dna)}\n`;
                    }
                    if (!dnaFasta) {
                        failCodonTrim("back-translation produced no DNA output", fallbackMetadata, "SPLACE internal fallback");
                        continue;
                    }
                    const validation = validateCodonPreservingTrimmedAlignment(dnaFasta, dnaCdsContent);
                    if (!validation.valid) {
                        failCodonTrim(validation.issues, fallbackMetadata, "SPLACE internal fallback");
                        continue;
                    }
                    if (skipped > 0) send("analysis-progress", { message: `[${name}] ⚠ ${skipped} sequences were unmatched before validation` });
                    dnaResult = {
                        content: dnaFasta,
                        method: "js-fallback",
                        validation,
                        metadata: {
                            ...fallbackMetadata,
                            codonIntegrityValidated: true,
                        },
                    };
                    send("analysis-progress", { message: `[${name}] ✓ Approximate JavaScript back-translation produced validated codon-preserving DNA` });
                }

                fs.writeFileSync(trimmedFile, dnaResult.content, "utf8");

                done++;
                trimmed++;
                mr.trimmedContent = dnaResult.content;
                mr.trimMethod = dnaResult.method;
                mr.trimmedStats = parseFastaStats(dnaResult.content);
                mr.trimmedFile = trimmedFile;
                mr.trimValidation = dnaResult.validation || null;
                mr.trimMetadata = dnaResult.metadata || null;
                mr.trimError = null;
                const dnaLen = mr.trimmedStats.length;
                const isDiv3 = mr.trimValidation?.allLengthsMultipleOf3 !== false && dnaLen % 3 === 0;
                send("analysis-progress", {
                    marker: name, status: "done",
                    info: `${mr.trimmedStats.numSeqs} seq · ${dnaLen} bp (DNA codon${isDiv3 ? "" : " ⚠"})`,
                    message: `[${name}] ✓ Codon-preserving DNA (${dnaResult.metadata?.trimmingTool === "trimal" ? "trimAl -backtrans" : dnaResult.metadata?.trimmingTool === "clipkit" ? "ClipKIT" : "SPLACE fallback"}) → ${mr.trimmedStats.numSeqs} seq · ${dnaLen} bp${!isDiv3 ? " ⚠ length not multiple of 3" : ""}`,
                    done, total: markers.length,
                });
                continue;
            }

            // ── Standard trimAl flow ────────────────────────────────────────────
            send("analysis-progress", {
                marker: name, status: "running",
                message: `[${name}] Running trimAl…`,
                done, total: markers.length,
            });

            const { success, trimmedContent, error } = await runTrimal(
                inputFile, trimmedFile, params,
                (msg) => send("analysis-progress", { message: `[${name}] ${msg}`, done, total: markers.length })
            );

            done++;
            if (success) {
                trimmed++;
                mr.trimmedContent = trimmedContent;
                mr.trimMethod = "standard";
                mr.trimmedStats = parseFastaStats(trimmedContent);
            }

            send("analysis-progress", {
                marker: name,
                status: success ? "done" : "error",
                info: success ? `${mr.trimmedStats.numSeqs} seq · ${mr.trimmedStats.length} bp` : error,
                message: success ? `[${name}] ✓ Trimmed → ${trimmedFile}` : `[${name}] ✗ ${error}`,
                done, total: markers.length,
            });
        }
    }

    await Promise.all(Array.from({ length: concurrency }, () => worker()));

    send("analysis-done", {
        phase: "trimal",
        success: codonMode ? (failedMarkers.length === 0 && trimmed === markers.length && trimmed > 0) : trimmed > 0,
        aligned: trimmed, total: markers.length,
        markerResults: updatedResults,
        failedMarkers,
    });
});

// -----------------------------------------------------------------------
// Helper: run a single MAFFT alignment
// -----------------------------------------------------------------------
function runMafft(inputFile, outputFile, params, threads, onLog) {
    return new Promise((resolve) => {
        const bin = resolveBin("mafft");
        const args = [...params, "--thread", String(threads), "--quiet", inputFile];

        // .bat files on Windows must be invoked via cmd /c
        const [cmd, cmdArgs] = (bin.toLowerCase().endsWith(".bat"))
            ? ["cmd", ["/c", bin, ...args]]
            : [bin, args];

        const proc = spawn(cmd, cmdArgs, { shell: false });
        let stdout = "";

        proc.stdout.on("data", (d) => { stdout += d.toString(); });
        proc.stderr.on("data", (d) => {
            const msg = d.toString("utf8").trim();
            if (!msg) return;
            // Filter known Windows setup noise from mafft.bat (chcp, env messages)
            if (/65001|Preparing environment|anti.virus|could not find \/tmp/i.test(msg)) return;
            onLog(msg);
        });

        proc.on("close", (code) => {
            if (code === 0 && stdout.trim()) {
                try {
                    fs.writeFileSync(outputFile, stdout, "utf8");
                    resolve({ success: true, alignedContent: stdout });
                } catch (e) {
                    resolve({ success: false, error: `Write failed: ${e.message}` });
                }
            } else {
                resolve({ success: false, error: `Exit code ${code}` });
            }
        });

        proc.on("error", (err) => {
            resolve({ success: false, error: err.message });
        });
    });
}

// -----------------------------------------------------------------------
// Helper: trim a multiple alignment (pure JS — no external binary required)
//
// Supported params (mirrors trimAl flags):
//   --automated1        heuristic: gap threshold derived from gap distribution
//   -gt <0..1>          max gap fraction per column (default 0.6)
//   -st <0..1>          min similarity score per column (default 0, off)
//   -cons <0..100>      min % sequences that must be kept (default 0, off)
//
// Algorithm:
//   1. Compute gap fraction for each alignment column.
//   2. Compute similarity (conservation) for each column using Shannon entropy
//      normalised to [0,1] so that a perfectly conserved column = 1.
//   3. Drop columns whose gap fraction > gapThreshold OR similarity < simThreshold.
//   4. --automated1: if gap distribution is bimodal (many near-0 and near-1
//      columns) use Jenks-style midpoint split; otherwise use gt=0.6.
// -----------------------------------------------------------------------
function trimAlignmentJS(alignedContent, params) {
    // ── Parse params ───────────────────────────────────────────────────────
    let gapThreshold = 0.6;
    let simThreshold = 0.0;
    let consThreshold = 0;
    let automated1 = false;

    for (let i = 0; i < params.length; i++) {
        const p = params[i];
        if (p === "--automated1" || p === "-automated1") { automated1 = true; }
        else if ((p === "-gt") && params[i + 1]) { gapThreshold = parseFloat(params[++i]); }
        else if ((p === "-st") && params[i + 1]) { simThreshold = parseFloat(params[++i]); }
        else if ((p === "-cons") && params[i + 1]) { consThreshold = parseFloat(params[++i]); }
    }

    // ── Parse FASTA ────────────────────────────────────────────────────────
    const seqs = [];
    let cur = null;
    for (const line of alignedContent.split(/\r?\n/)) {
        if (line.startsWith(">")) {
            if (cur) seqs.push(cur);
            cur = { name: line.slice(1).trim(), seq: "" };
        } else if (cur) {
            cur.seq += line.trim();
        }
    }
    if (cur) seqs.push(cur);
    if (!seqs.length) return { success: false, error: "Empty alignment" };

    const nSeqs = seqs.length;
    const seqLen = seqs[0].seq.length;

    // ── Per-column stats ────────────────────────────────────────────────────
    const gapFrac = new Float32Array(seqLen);
    const simScore = new Float32Array(seqLen);
    const log2 = Math.log(2);

    for (let col = 0; col < seqLen; col++) {
        let gaps = 0;
        const freq = {};
        for (const s of seqs) {
            const c = s.seq[col].toUpperCase();
            if (c === "-" || c === "?" || c === "N") { gaps++; }
            else { freq[c] = (freq[c] || 0) + 1; }
        }
        gapFrac[col] = gaps / nSeqs;

        // Shannon entropy → normalised similarity
        const nonGap = nSeqs - gaps;
        if (nonGap <= 1) {
            simScore[col] = nonGap === 1 ? 1.0 : 0.0;
        } else {
            let H = 0;
            for (const n of Object.values(freq)) {
                const p = n / nonGap;
                H -= p * Math.log(p) / log2;
            }
            const maxH = Math.log(Math.min(nonGap, 4)) / log2; // max 2 bits for DNA
            simScore[col] = maxH > 0 ? Math.max(0, 1 - H / maxH) : 1.0;
        }
    }

    // ── Automated1: derive gap threshold from distribution ──────────────────
    if (automated1) {
        // Sort gap fractions; use the midpoint of the largest gap between
        // the two clusters (low-gap "good" columns vs high-gap "bad" columns)
        const sorted = Array.from(gapFrac).sort((a, b) => a - b);
        let maxGap = 0, splitAt = gapThreshold;
        for (let i = 1; i < sorted.length; i++) {
            const d = sorted[i] - sorted[i - 1];
            if (d > maxGap) { maxGap = d; splitAt = (sorted[i] + sorted[i - 1]) / 2; }
        }
        // Only adopt the split threshold if there is a clear bimodal gap
        // (jump > 0.05) and it falls in [0.1, 0.9]
        if (maxGap > 0.05 && splitAt >= 0.1 && splitAt <= 0.9) {
            gapThreshold = splitAt;
        }
    }

    // ── Select columns to keep ──────────────────────────────────────────────
    let keepCols = [];
    for (let col = 0; col < seqLen; col++) {
        if (gapFrac[col] <= gapThreshold && simScore[col] >= simThreshold) {
            keepCols.push(col);
        }
    }

    // ── Enforce -cons: ensure at least consThreshold % of sequences are
    //    represented in the kept columns (relax gap threshold if needed) ─────
    if (consThreshold > 0 && keepCols.length === 0) {
        // Fallback: keep columns sorted by ascending gap fraction until
        // the minimum representation is satisfied
        const byGap = Array.from({ length: seqLen }, (_, i) => i)
            .sort((a, b) => gapFrac[a] - gapFrac[b]);
        keepCols = byGap.slice(0, Math.max(1, Math.ceil(seqLen * consThreshold / 100)));
    }

    if (!keepCols.length) {
        return { success: false, error: "All columns removed by trimming thresholds" };
    }

    // ── Build trimmed sequences ─────────────────────────────────────────────
    const trimmedFasta = seqs
        .map(s => `>${s.name}\n${keepCols.map(i => s.seq[i]).join("")}`)
        .join("\n") + "\n";

    return { success: true, trimmedContent: trimmedFasta, keptCols: keepCols.length, totalCols: seqLen };
}

// -----------------------------------------------------------------------
// Helper: parse FASTA text into an array of { name, seq } objects
// -----------------------------------------------------------------------
function parseFastaArr(text) {
    const seqs = [];
    let cur = null;
    for (const line of (text || "").split(/\r?\n/)) {
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

function normalizeValidationDna(sequence) {
    return (sequence || "").toUpperCase().replace(/\s+/g, "").replace(/U/g, "T");
}

function splitIntoCodons(sequence) {
    const dna = normalizeValidationDna(sequence);
    const codons = [];
    for (let i = 0; i + 2 < dna.length; i += 3) {
        codons.push(dna.slice(i, i + 3));
    }
    return codons;
}

function isOrderedCodonSubsequence(candidateSequence, originalSequence) {
    const candidateCodons = splitIntoCodons(candidateSequence);
    const originalCodons = splitIntoCodons(originalSequence);
    let searchIndex = 0;

    for (const codon of candidateCodons) {
        while (searchIndex < originalCodons.length && originalCodons[searchIndex] !== codon) {
            searchIndex++;
        }
        if (searchIndex >= originalCodons.length) {
            return false;
        }
        searchIndex++;
    }

    return true;
}

function formatValidationNameList(names) {
    if (!names.length) return "";
    const preview = names.slice(0, 3).join(", ");
    return names.length > 3 ? `${preview} (+${names.length - 3} more)` : preview;
}

function validateCodonPreservingTrimmedAlignment(trimmedContent, originalDnaContent) {
    const trimmedSeqs = parseFastaArr(trimmedContent);
    const originalSeqs = parseFastaArr(originalDnaContent);
    const issues = [];

    if (!trimmedSeqs.length) {
        return {
            valid: false,
            issues: ["trimAl -backtrans produced no DNA sequences"],
            equalLengths: false,
            allLengthsMultipleOf3: false,
            gapBlocksAreTriplets: false,
            headersMatch: false,
            dnaCompatible: false,
            frameDisrupted: true,
            length: 0,
            missingHeaders: originalSeqs.map((seq) => seq.name),
            extraHeaders: [],
            subsequenceIssues: [],
        };
    }

    const originalMap = new Map(originalSeqs.map((seq) => [seq.name, normalizeValidationDna(seq.seq)]));
    const originalNames = originalSeqs.map((seq) => seq.name);
    const trimmedNames = trimmedSeqs.map((seq) => seq.name);
    const originalNameSet = new Set(originalNames);
    const trimmedNameSet = new Set(trimmedNames);
    const missingHeaders = originalNames.filter((name) => !trimmedNameSet.has(name));
    const extraHeaders = trimmedNames.filter((name) => !originalNameSet.has(name));
    const lengths = trimmedSeqs.map((seq) => seq.seq.length);
    const length = lengths[0] || 0;
    const equalLengths = lengths.every((seqLen) => seqLen === length);
    const allLengthsMultipleOf3 = trimmedSeqs.every((seq) => seq.seq.length % 3 === 0);
    const gapBlocksAreTriplets = trimmedSeqs.every((seq) => (seq.seq.match(/-+/g) || []).every((run) => run.length % 3 === 0));
    const dnaCompatible = trimmedSeqs.every((seq) => /^[ACGTNRYSWKMBDHV\-\?]*$/i.test(normalizeValidationDna(seq.seq)));
    const subsequenceIssues = [];

    for (const seq of trimmedSeqs) {
        const original = originalMap.get(seq.name);
        if (!original) continue;
        const trimmedUngapped = normalizeValidationDna(seq.seq).replace(/-/g, "");
        if (!trimmedUngapped.length || trimmedUngapped.length % 3 !== 0) continue;
        if (!isOrderedCodonSubsequence(trimmedUngapped, original)) {
            subsequenceIssues.push(`${seq.name}: trimmed DNA is not a codon-ordered subsequence of the validated CDS`);
        }
    }

    if (!equalLengths) issues.push("trimmed sequences do not all have equal length");
    if (!allLengthsMultipleOf3) issues.push(`final trimmed length ${length} is not divisible by 3`);
    if (!gapBlocksAreTriplets) issues.push("trimmed alignment has gap runs that are not multiples of 3");
    if (!dnaCompatible) issues.push("trimAl -backtrans output contains non-DNA symbols");
    if (missingHeaders.length) issues.push(`missing sequence names from validated CDS: ${formatValidationNameList(missingHeaders)}`);
    if (extraHeaders.length) issues.push(`unexpected sequence names in trimmed alignment: ${formatValidationNameList(extraHeaders)}`);
    if (subsequenceIssues.length) {
        issues.push(...subsequenceIssues.slice(0, 3));
        if (subsequenceIssues.length > 3) {
            issues.push(`... and ${subsequenceIssues.length - 3} more sequence(s) with codon frame disruption`);
        }
    }

    const frameDisrupted =
        !equalLengths ||
        !allLengthsMultipleOf3 ||
        !gapBlocksAreTriplets ||
        !dnaCompatible ||
        missingHeaders.length > 0 ||
        extraHeaders.length > 0 ||
        subsequenceIssues.length > 0;

    return {
        valid: issues.length === 0,
        issues,
        equalLengths,
        allLengthsMultipleOf3,
        gapBlocksAreTriplets,
        headersMatch: missingHeaders.length === 0 && extraHeaders.length === 0,
        dnaCompatible,
        frameDisrupted,
        length,
        missingHeaders,
        extraHeaders,
        subsequenceIssues,
    };
}

function removeCodonPositionFromFasta(codonDnaContent, positionToRemove) {
    const position = Number(positionToRemove);
    if (![1, 2, 3].includes(position)) {
        throw new Error(`Invalid codon position to remove: ${positionToRemove}`);
    }
    const seqs = parseFastaArr(codonDnaContent);
    if (!seqs.length) throw new Error("No sequences available for codon-position filtering");
    const filtered = [];
    const lengths = [];
    for (const seq of seqs) {
        const raw = String(seq.seq || "");
        if (raw.length % 3 !== 0) {
            throw new Error(`${seq.name}: codon-aligned DNA length ${raw.length} is not divisible by 3`);
        }
        let output = "";
        for (let i = 0; i < raw.length; i += 3) {
            const codon = raw.slice(i, i + 3);
            if (codon.length !== 3) throw new Error(`${seq.name}: incomplete codon block at position ${i + 1}`);
            output += codon.split("").filter((_, idx) => idx !== position - 1).join("");
        }
        filtered.push({ name: seq.name, seq: output });
        lengths.push(output.length);
    }
    const equalLengths = lengths.every((len) => len === lengths[0]);
    if (!equalLengths) throw new Error("Position-filtered sequences do not all have equal length");
    return filtered.map((seq) => `>${seq.name}\n${seq.seq}`).join("\n") + "\n";
}

function buildCodonAlignedDnaFromProteinAlignment(proteinAlignmentContent, originalDnaContent) {
    const protSeqs = parseFastaArr(proteinAlignmentContent);
    const dnaSeqs = parseFastaArr(originalDnaContent);
    const dnaMap = {};
    for (const s of dnaSeqs) dnaMap[s.name] = s.seq.replace(/-/g, "");
    let dnaFasta = "";
    let skipped = 0;
    for (const ps of protSeqs) {
        const dna = dnaMap[ps.name];
        if (!dna) { skipped++; continue; }
        dnaFasta += `>${ps.name}
${backTranslateJS(ps.seq, dna)}
`;
    }
    return { dnaFasta, skipped };
}

function isTrimAlUnavailableError(message) {
    return /ENOENT|not found|not recognized|No such file or directory/i.test(message || "");
}

// -----------------------------------------------------------------------
// Helper: back-translate an aligned protein sequence to DNA using the
// original CDS codons.  Each '-' gap in the protein → '---' in DNA.
// -----------------------------------------------------------------------
function backTranslateJS(alignedProteinSeq, originalDna) {
    const dna = originalDna.toUpperCase().replace(/U/g, "T").replace(/[^ACGTNRYSWKMBDHV]/g, "");
    const codons = [];
    for (let i = 0; i + 2 < dna.length; i += 3) {
        codons.push(dna.slice(i, i + 3));
    }
    let idx = 0;
    let output = "";
    for (const aa of alignedProteinSeq) {
        if (aa === "-" || aa === ".") {
            output += "---";
        } else {
            output += codons[idx] || "NNN";
            idx++;
        }
    }
    return output;
}

// -----------------------------------------------------------------------
// Helper: run ClipKIT through Python or the console script
// -----------------------------------------------------------------------
function runClipkitProteinTrim(inputFile, outputFile, mode, onLog) {
    const safeMode = String(mode || "gappy").replace(/[^A-Za-z0-9_-]/g, "") || "gappy";
    const attempts = [];
    const py = findPython();
    if (py) attempts.push({ cmd: py.cmd, args: [...py.prefixArgs, "-m", "clipkit", inputFile, "-o", outputFile, "-m", safeMode], label: "python -m clipkit" });
    attempts.push({ cmd: resolveBin("clipkit"), args: [inputFile, "-o", outputFile, "-m", safeMode], label: "clipkit" });

    return new Promise((resolve) => {
        let idx = 0;
        function nextAttempt(lastError = null) {
            if (idx >= attempts.length) {
                resolve({ success: false, error: lastError || "ClipKIT is not available. Install it with the Python tools button." });
                return;
            }
            const attempt = attempts[idx++];
            onLog(`Running ${attempt.label} -m ${safeMode}`);
            let stderr = "";
            let stdout = "";
            const proc = spawn(attempt.cmd, attempt.args, { shell: false, windowsHide: true });
            proc.stdout.on("data", d => { stdout += d.toString(); });
            proc.stderr.on("data", d => {
                const msg = d.toString().trim();
                stderr += d.toString();
                if (msg) onLog(msg);
            });
            proc.on("close", code => {
                if (code === 0 && fs.existsSync(outputFile)) {
                    resolve({ success: true, trimmedContent: fs.readFileSync(outputFile, "utf8"), stdout, stderr, method: "clipkit" });
                } else {
                    nextAttempt(stderr.trim() || `ClipKIT exit code ${code}`);
                }
            });
            proc.on("error", err => nextAttempt(err.message));
        }
        nextAttempt();
    });
}

// -----------------------------------------------------------------------
// Helper: run trimAl — native binary on Linux/macOS, JS fallback on Windows
// -----------------------------------------------------------------------
function runTrimal(inputFile, outputFile, params, onLog) {
    // Use the bundled native binary only on Linux/macOS — on Windows the exe
    // has missing DLL dependencies, so always fall through to the JS fallback.
    const bin = resolveBin("trimal");
    if (!IS_WIN && bin !== "trimal") {
        return new Promise((resolve) => {
            const proc = spawn(bin, ["-in", inputFile, "-out", outputFile, ...params], { shell: false });
            proc.stderr.on("data", (d) => { const m = d.toString().trim(); if (m) onLog(m); });
            proc.on("close", (code) => {
                if (code === 0 && fs.existsSync(outputFile)) {
                    resolve({ success: true, trimmedContent: fs.readFileSync(outputFile, "utf8") });
                } else {
                    resolve({ success: false, error: `Exit code ${code}` });
                }
            });
            proc.on("error", (err) => resolve({ success: false, error: err.message }));
        });
    }

    // JS fallback — used on Windows where no native trimAl binary is available
    return new Promise((resolve) => {
        try {
            const alignedContent = fs.readFileSync(inputFile, "utf8");
            const result = trimAlignmentJS(alignedContent, params);
            if (!result.success) { resolve({ success: false, error: result.error }); return; }
            fs.writeFileSync(outputFile, result.trimmedContent, "utf8");
            onLog(`Kept ${result.keptCols}/${result.totalCols} columns`);
            resolve({ success: true, trimmedContent: result.trimmedContent });
        } catch (e) {
            resolve({ success: false, error: e.message });
        }
    });
}



// -----------------------------------------------------------------------
// Helper: parse basic FASTA stats
// -----------------------------------------------------------------------
function parseFastaStats(text) {
    const seqs = [];
    let cur = null;
    for (const line of (text || "").split(/\r?\n/)) {
        if (line.startsWith(">")) {
            if (cur !== null) seqs.push(cur);
            cur = 0;
        } else if (cur !== null) {
            cur += line.trim().length;
        }
    }
    if (cur !== null) seqs.push(cur);

    const numSeqs = seqs.length;
    const length = seqs[0] || 0;
    const avgLen = numSeqs > 0 ? Math.round(seqs.reduce((a, b) => a + b, 0) / numSeqs) : 0;

    // Gap percentage of first (aligned) sequence set
    const gapCount = (text || "").replace(/^>.*$/gm, "").replace(/\n/g, "").split("").filter(c => c === "-").length;
    const totalChars = (text || "").replace(/^>.*$/gm, "").replace(/\n/g, "").length;
    const gapPct = totalChars > 0 ? ((gapCount / totalChars) * 100).toFixed(1) : "0.0";

    return { numSeqs, length, avgLen, gapPct };
}

function getCodonOutputFiles(outputDir, marker) {
    return {
        rawFile: path.join(outputDir, `${marker}_cds_validated.fasta`),
        proteinAlignedFile: path.join(outputDir, `${marker}_protein_aligned.fasta`),
        codonAlignedFile: path.join(outputDir, `${marker}_codon_aligned.fasta`),
        codonTrimmedFile: path.join(outputDir, `${marker}_codon_trimmed.fasta`),
    };
}

function persistCodonMarkerOutputs(outputDir, markerResult) {
    const files = getCodonOutputFiles(outputDir, markerResult.marker);
    const legacyAlignedFile = markerResult.outputFile;
    const proteinContent = markerResult.proteinAlignedContent || markerResult.alignedContent || "";

    fs.mkdirSync(outputDir, { recursive: true });

    if (markerResult.rawContent) {
        fs.writeFileSync(files.rawFile, markerResult.rawContent, "utf8");
    }

    if (legacyAlignedFile && fs.existsSync(legacyAlignedFile) && path.resolve(legacyAlignedFile) !== path.resolve(files.proteinAlignedFile)) {
        try { fs.rmSync(files.proteinAlignedFile, { force: true }); } catch { }
        try {
            fs.renameSync(legacyAlignedFile, files.proteinAlignedFile);
        } catch {
            fs.copyFileSync(legacyAlignedFile, files.proteinAlignedFile);
            try { fs.rmSync(legacyAlignedFile, { force: true }); } catch { }
        }
    } else if (proteinContent) {
        fs.writeFileSync(files.proteinAlignedFile, proteinContent, "utf8");
    }

    if (markerResult.alignedContent) {
        fs.writeFileSync(files.codonAlignedFile, markerResult.alignedContent, "utf8");
    }

    if (markerResult.trimmedContent) {
        fs.writeFileSync(files.codonTrimmedFile, markerResult.trimmedContent, "utf8");
    }

    return {
        ...markerResult,
        rawFile: markerResult.rawContent ? files.rawFile : null,
        proteinAlignedFile: proteinContent ? files.proteinAlignedFile : null,
        outputFile: markerResult.alignedContent ? files.codonAlignedFile : markerResult.outputFile,
        trimmedFile: markerResult.trimmedContent ? files.codonTrimmedFile : (markerResult.trimmedFile || null),
    };
}

ipcMain.handle("save-codon-outputs", async (_event, { outputDir, markerResults }) => {
    try {
        if (!outputDir) throw new Error("Missing output directory");
        const updatedResults = (markerResults || []).map((mr) => persistCodonMarkerOutputs(outputDir, mr));
        _lastMarkerResults = updatedResults;
        return { success: true, outputDir, markerResults: updatedResults };
    } catch (error) {
        return { success: false, error: error.message };
    }
});

// -----------------------------------------------------------------------
// IPC: run IQ-TREE3
// { nexus, partition, fasta, params, threads, outgroup, fileNames, extraFiles }
// -----------------------------------------------------------------------
ipcMain.on("run-iqtree", async (event, { nexus, partition, fasta, params, threads, outgroup, perGene, geneFiles, fileNames, extraFiles }) => {
    const sender = event.sender;
    const send = (channel, data) => { if (!sender.isDestroyed()) sender.send(channel, data); };

    const outputDir = path.join(os.homedir(), "Documents", "SPLACE", "iqtree");
    try { fs.rmSync(outputDir, { recursive: true, force: true }); } catch { }
    fs.mkdirSync(outputDir, { recursive: true });

    const activeFileNames = {
        nexus: fileNames?.nexus || "concatenated.nex",
        partition: fileNames?.partition || "partitions.txt",
        fasta: fileNames?.fasta || "concatenated.fasta",
        prefix: fileNames?.prefix || "concat_tree",
    };
    const nexusFile = path.join(outputDir, activeFileNames.nexus);
    const partitionFile = path.join(outputDir, activeFileNames.partition);
    const fastaFile = path.join(outputDir, activeFileNames.fasta);
    const concatPrefix = path.join(outputDir, activeFileNames.prefix);

    const writeOutputFile = (targetPath, content) => {
        if (typeof content !== "string") return;
        fs.mkdirSync(path.dirname(targetPath), { recursive: true });
        fs.writeFileSync(targetPath, content, "utf8");
    };

    writeOutputFile(nexusFile, nexus);
    writeOutputFile(partitionFile, partition);
    writeOutputFile(fastaFile, fasta || "");
    for (const [relativeName, content] of Object.entries(extraFiles || {})) {
        if (!relativeName) continue;
        writeOutputFile(path.join(outputDir, relativeName), content);
    }

    // Resolve IQ-TREE binary: bundled first, then iqtree3 > iqtree2 > iqtree
    function findIqtree() {
        for (const name of ["iqtree3", "iqtree2", "iqtree"]) {
            const resolved = resolveBin(name);
            if (resolved !== name) return resolved;
        }
        return "iqtree3";
    }
    const iqtreeBin = findIqtree();

    // Helper: run one IQ-TREE process and stream output
    function runOne(args, label) {
        return new Promise((resolve) => {
            send("iqtree-progress", { message: `\n[${label}] ${iqtreeBin} ${args.join(" ")}` });
            const proc = spawn(iqtreeBin, args, { shell: true, stdio: ["ignore", "pipe", "pipe"] });
            proc.stdout.on("data", (d) => {
                const msg = d.toString().trim();
                if (msg) send("iqtree-progress", { message: `[${label}] ${msg}` });
            });
            proc.stderr.on("data", (d) => {
                const msg = d.toString().trim();
                if (msg) send("iqtree-progress", { message: `[${label}] ${msg}` });
            });
            proc.on("close", (code) => resolve({ code, label }));
            proc.on("error", (err) => resolve({ code: -1, label, error: err.message }));
        });
    }

    const outgroupArgs = outgroup && outgroup.length ? ["-o", outgroup.join(",")] : [];

    // === Concatenated tree ===
    const concatArgs = [
        "-s", nexusFile,
        ...(partition ? ["-p", partitionFile] : []),
        "--prefix", concatPrefix,
        "-T", String(threads),
        "--redo",
        ...outgroupArgs,
        ...params,
    ];
    const concatResult = await runOne(concatArgs, "Concatenated");

    // === Per-gene trees (optional) ===
    const geneResults = [];
    if (perGene && geneFiles) {
        const genesDir = path.join(outputDir, "gene_trees");
        fs.mkdirSync(genesDir, { recursive: true });
        const geneNames = Object.keys(geneFiles);
        for (const gene of geneNames) {
            const fastaFile = path.join(genesDir, `${gene}.fasta`);
            const genePrefix = path.join(genesDir, `${gene}_tree`);
            fs.writeFileSync(fastaFile, geneFiles[gene], "utf8");
            const geneArgs = [
                "-s", fastaFile,
                "--prefix", genePrefix,
                "-T", String(threads),
                "--redo",
                ...outgroupArgs,
                ...params,
            ];
            const res = await runOne(geneArgs, gene);
            geneResults.push(res);
        }
    }

    // Collect output files
    let files = [];
    try {
        const visibleConcatFastas = new Set(["concatenated.fasta", "concatenated_no3rd.fasta"]);
        files = fs.readdirSync(outputDir)
            .filter(f => !f.endsWith(".fasta") || visibleConcatFastas.has(f))
            .map(f => f);
        if (perGene) {
            const genesDir = path.join(outputDir, "gene_trees");
            try {
                fs.readdirSync(genesDir).forEach(f => files.push(`gene_trees/${f}`));
            } catch { }
        }
    } catch { }

    const success = concatResult.code === 0;
    send("iqtree-done", {
        success,
        outputDir,
        files,
        error: !success ? `Concatenated tree exit code ${concatResult.code}` : null,
        geneResults: geneResults.map(r => ({ gene: r.label, success: r.code === 0 })),
    });
});

// -----------------------------------------------------------------------
// IPC: save ZIP with organised folder structure
// { alignmentDir, iqtreeDir, rawFastas, suggestedName }
//   raw/          — raw per-gene FASTAs (from memory, not disk)
// -----------------------------------------------------------------------
// seqkit translate: nucleotide FASTA content → protein FASTA content
// -----------------------------------------------------------------------
ipcMain.handle("translate-fasta", async (event, { content, code }) => {
    const seqkitBin = resolveBin("seqkit");
    const ts = Date.now();
    const tmpIn = path.join(os.tmpdir(), `splace-tx-in-${ts}.fasta`);
    const tmpOut = path.join(os.tmpdir(), `splace-tx-out-${ts}.fasta`);
    try {
        fs.writeFileSync(tmpIn, content, "utf8");
        const result = spawnSync(seqkitBin, ["translate", "-T", String(code || 1), "-o", tmpOut, tmpIn], {
            encoding: "utf8",
            timeout: 30000,
        });
        if (result.status !== 0) {
            throw new Error(result.stderr || "seqkit translate failed");
        }
        return fs.readFileSync(tmpOut, "utf8");
    } finally {
        try { fs.unlinkSync(tmpIn); } catch { }
        try { fs.unlinkSync(tmpOut); } catch { }
    }
});

//   aligned/      — *_aligned.fasta files from alignmentDir
//   trimmed/      — *_aligned_trimmed.fasta files from alignmentDir (if any)
//   concatenated/ — concatenated matrices and partition files from iqtreeDir
//   trees/        — all other iqtree output files (treefile, log, etc.)
// -----------------------------------------------------------------------
ipcMain.handle("save-zip", async (event, { alignmentDir, iqtreeDir, rawFastas, suggestedName }) => {
    const { filePath, canceled } = await dialog.showSaveDialog({
        defaultPath: path.join(os.homedir(), "Documents", suggestedName || "splace_results.zip"),
        filters: [{ name: "ZIP archive", extensions: ["zip"] }],
    });
    if (canceled || !filePath) return { cancelled: true };

    const stagingDir = fs.mkdtempSync(path.join(os.tmpdir(), "splace-zip-"));
    try {
        // raw/ — write raw FASTA content passed from renderer
        if (rawFastas && Object.keys(rawFastas).length) {
            const rawDir = path.join(stagingDir, "raw");
            fs.mkdirSync(rawDir, { recursive: true });
            for (const [gene, content] of Object.entries(rawFastas)) {
                fs.writeFileSync(path.join(rawDir, `${gene}_cds_validated.fasta`), content, "utf8");
            }
        }

        // aligned/ and trimmed/ — sort files from alignmentDir by suffix
        if (alignmentDir && fs.existsSync(alignmentDir)) {
            const alignedDir = path.join(stagingDir, "aligned");
            fs.mkdirSync(alignedDir, { recursive: true });
            let trimmedDir = null;
            let proteinAlignedDir = null;

            for (const f of fs.readdirSync(alignmentDir)) {
                const src = path.join(alignmentDir, f);
                try { if (!fs.statSync(src).isFile()) continue; } catch { continue; }

                if (f.endsWith("_codon_trimmed.fasta") || f.endsWith("_trimmed.fasta") || f.includes("_aligned_trimmed")) {
                    if (!trimmedDir) {
                        trimmedDir = path.join(stagingDir, "trimmed");
                        fs.mkdirSync(trimmedDir, { recursive: true });
                    }
                    fs.copyFileSync(src, path.join(trimmedDir, f));
                } else if (f.endsWith("_protein_aligned.fasta")) {
                    if (!proteinAlignedDir) {
                        proteinAlignedDir = path.join(stagingDir, "protein_aligned");
                        fs.mkdirSync(proteinAlignedDir, { recursive: true });
                    }
                    fs.copyFileSync(src, path.join(proteinAlignedDir, f));
                } else if (f.endsWith("_codon_aligned.fasta")) {
                    fs.copyFileSync(src, path.join(alignedDir, f));
                } else if (f.endsWith("_cds_validated.fasta")) {
                    // raw/ is already populated from renderer memory to match the UI exactly.
                    continue;
                } else if (f.endsWith(".fasta") || f.endsWith(".fa")) {
                    fs.copyFileSync(src, path.join(alignedDir, f));
                }
            }
        }

        // concatenated/ and trees/ — split iqtreeDir contents
        if (iqtreeDir && fs.existsSync(iqtreeDir)) {
            const concatDir = path.join(stagingDir, "concatenated");
            const treesDir = path.join(stagingDir, "trees");
            fs.mkdirSync(concatDir, { recursive: true });
            fs.mkdirSync(treesDir, { recursive: true });

            const CONCAT_FILES = new Set([
                "concatenated.nex",
                "partitions.txt",
                "concatenated.fasta",
                "concatenated_no3rd.nex",
                "partitions_no3rd.txt",
                "concatenated_no3rd.fasta",
            ]);

            for (const f of fs.readdirSync(iqtreeDir)) {
                const src = path.join(iqtreeDir, f);
                try {
                    const stat = fs.statSync(src);
                    if (stat.isDirectory()) {
                        // gene_trees/ subdirectory → trees/gene_trees/
                        if (f === "gene_trees") {
                            const dest = path.join(treesDir, "gene_trees");
                            fs.mkdirSync(dest, { recursive: true });
                            for (const gf of fs.readdirSync(src)) {
                                const gsrc = path.join(src, gf);
                                if (fs.statSync(gsrc).isFile())
                                    fs.copyFileSync(gsrc, path.join(dest, gf));
                            }
                        }
                        continue;
                    }
                    if (CONCAT_FILES.has(f)) {
                        fs.copyFileSync(src, path.join(concatDir, f));
                    } else {
                        fs.copyFileSync(src, path.join(treesDir, f));
                    }
                } catch { }
            }
        }
    } catch (e) {
        try { fs.rmSync(stagingDir, { recursive: true, force: true }); } catch { }
        return { success: false, error: e.message };
    }

    return new Promise((resolve) => {
        let proc;
        if (IS_WIN) {
            // PowerShell Compress-Archive — available on all modern Windows
            proc = spawn("powershell", [
                "-NoProfile", "-NonInteractive", "-Command",
                `Compress-Archive -Path "${stagingDir}\\*" -DestinationPath "${filePath}" -Force`,
            ], { shell: false, stdio: ["ignore", "pipe", "pipe"] });
        } else {
            proc = spawn("zip", ["-r", filePath, "."], {
                cwd: stagingDir, shell: true, stdio: ["ignore", "pipe", "pipe"],
            });
        }
        proc.on("close", (code) => {
            try { fs.rmSync(stagingDir, { recursive: true, force: true }); } catch { }
            if (code === 0) resolve({ success: true, filePath });
            else resolve({ success: false, error: `zip exit code ${code}` });
        });
        proc.on("error", (err) => {
            try { fs.rmSync(stagingDir, { recursive: true, force: true }); } catch { }
            resolve({ success: false, error: err.message });
        });
    });
});
