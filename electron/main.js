const { app, BrowserWindow, ipcMain, dialog, shell } = require("electron");
const path = require("path");
const fs = require("fs");
const os = require("os");
const { spawn, spawnSync } = require("child_process");

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

ipcMain.on("run-trimal", async (event, { markers, params, codonMode, cdsFastas, alignmentMode }) => {
    const sender = event.sender;
    const send = (channel, data) => { if (!sender.isDestroyed()) sender.send(channel, data); };

    let done = 0, trimmed = 0;
    const updatedResults = _lastMarkerResults.map(r => ({ ...r }));
    const concurrency = Math.max(1, os.cpus().length - 2);
    let nextIdx = 0;

    async function worker() {
        while (true) {
            const idx = nextIdx++;
            if (idx >= markers.length) break;
            const name = markers[idx];
            const mr = updatedResults.find(r => r.marker === name);
            if (!mr || !mr.alignedContent || !mr.outputFile) {
                done++;
                send("analysis-progress", { marker: name, status: "error", message: `[${name}] No aligned file found`, done, total: markers.length });
                continue;
            }

            const inputFile = mr.outputFile;
            const trimmedFile = inputFile.replace(/_aligned\.fasta$/, "_trimmed.fasta");

            // ── Codon-aware trimming (back-translation) ─────────────────────────
            if (codonMode && cdsFastas && cdsFastas[name]) {
                send("analysis-progress", {
                    marker: name, status: "running",
                    message: `[${name}] Running codon-aware trimming…`,
                    done, total: markers.length,
                });

                const dnaCdsContent = cdsFastas[name];
                let dnaResult = null;

                // Try trimAl -backtrans
                const trimalBin = resolveBin("trimal");
                if (trimalBin !== "trimal") {
                    const tmpCdsFile = path.join(os.tmpdir(), `splace-cds-${Date.now()}.fasta`);
                    const tmpOutFile = path.join(os.tmpdir(), `splace-bt-${Date.now()}.fasta`);
                    try {
                        fs.writeFileSync(tmpCdsFile, dnaCdsContent, "utf8");
                        const btArgs = ["-in", inputFile, "-backtrans", tmpCdsFile, "-out", tmpOutFile, "-fasta", ...params];
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
                                    resolve({ success: false, error: stderr.trim() || `Exit code ${code}` });
                                }
                            });
                            proc.on("error", err => resolve({ success: false, error: err.message }));
                        });
                        if (btRes.success) {
                            dnaResult = { content: btRes.content, method: "trimal-backtrans" };
                            send("analysis-progress", { message: `[${name}] ✓ trimAl -backtrans succeeded` });
                        } else {
                            send("analysis-progress", { message: `[${name}] trimAl -backtrans failed (${btRes.error}), using JS fallback` });
                        }
                    } catch (e) {
                        send("analysis-progress", { message: `[${name}] trimAl -backtrans error: ${e.message}, using JS fallback` });
                    } finally {
                        try { fs.unlinkSync(tmpCdsFile); } catch { }
                        try { fs.unlinkSync(tmpOutFile); } catch { }
                    }
                }

                if (!dnaResult) {
                    // JS fallback: trim protein alignment, then back-translate to DNA
                    send("analysis-progress", { message: `[${name}] Using JS codon-aware fallback` });
                    const protContent = fs.existsSync(inputFile) ? fs.readFileSync(inputFile, "utf8") : (mr.alignedContent || "");
                    const trimResult = trimAlignmentJS(protContent, params);
                    if (!trimResult.success) {
                        done++;
                        send("analysis-progress", {
                            marker: name, status: "error",
                            message: `[${name}] Trimming failed: ${trimResult.error}`,
                            done, total: markers.length,
                        });
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
                        done++;
                        send("analysis-progress", {
                            marker: name, status: "error",
                            message: `[${name}] Back-translation produced no output`,
                            done, total: markers.length,
                        });
                        continue;
                    }
                    if (skipped > 0) send("analysis-progress", { message: `[${name}] ⚠ ${skipped} sequences skipped (header mismatch)` });
                    dnaResult = { content: dnaFasta, method: "js-fallback" };
                    send("analysis-progress", { message: `[${name}] ✓ JS back-translation complete` });
                }

                done++;
                trimmed++;
                mr.trimmedContent = dnaResult.content;
                mr.trimMethod = dnaResult.method;
                mr.trimmedStats = parseFastaStats(dnaResult.content);
                const dnaLen = mr.trimmedStats.length;
                const isDiv3 = dnaLen % 3 === 0;
                send("analysis-progress", {
                    marker: name, status: "done",
                    info: `${mr.trimmedStats.numSeqs} seq · ${dnaLen} bp (DNA codon${isDiv3 ? "" : " ⚠"})`,
                    message: `[${name}] ✓ Codon alignment (${dnaResult.method}) → ${mr.trimmedStats.numSeqs} seq · ${dnaLen} bp${!isDiv3 ? " ⚠ length not multiple of 3" : ""}`,
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
        success: trimmed > 0,
        aligned: trimmed, total: markers.length,
        markerResults: updatedResults,
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

// -----------------------------------------------------------------------
// IPC: run IQ-TREE3
// { nexus, partition, params, threads, outgroup }
// -----------------------------------------------------------------------
ipcMain.on("run-iqtree", async (event, { nexus, partition, params, threads, outgroup, perGene, geneFiles }) => {
    const sender = event.sender;
    const send = (channel, data) => { if (!sender.isDestroyed()) sender.send(channel, data); };

    const outputDir = path.join(os.homedir(), "Documents", "SPLACE", "iqtree");
    try { fs.rmSync(outputDir, { recursive: true, force: true }); } catch { }
    fs.mkdirSync(outputDir, { recursive: true });

    const nexusFile = path.join(outputDir, "concatenated.nex");
    const partitionFile = path.join(outputDir, "partitions.txt");
    const concatPrefix = path.join(outputDir, "concat_tree");

    fs.writeFileSync(nexusFile, nexus, "utf8");
    fs.writeFileSync(partitionFile, partition, "utf8");

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
        "-p", partitionFile,
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
        files = fs.readdirSync(outputDir)
            .filter(f => !f.endsWith(".fasta"))
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
//   concatenated/ — concatenated.nex + partitions.txt from iqtreeDir
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
                fs.writeFileSync(path.join(rawDir, `${gene}.fasta`), content, "utf8");
            }
        }

        // aligned/ and trimmed/ — sort files from alignmentDir by suffix
        if (alignmentDir && fs.existsSync(alignmentDir)) {
            const alignedDir = path.join(stagingDir, "aligned");
            fs.mkdirSync(alignedDir, { recursive: true });
            let trimmedDir = null;

            for (const f of fs.readdirSync(alignmentDir)) {
                const src = path.join(alignmentDir, f);
                try { if (!fs.statSync(src).isFile()) continue; } catch { continue; }

                if (f.endsWith("_trimmed.fasta") || f.includes("_aligned_trimmed")) {
                    if (!trimmedDir) {
                        trimmedDir = path.join(stagingDir, "trimmed");
                        fs.mkdirSync(trimmedDir, { recursive: true });
                    }
                    fs.copyFileSync(src, path.join(trimmedDir, f));
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

            const CONCAT_FILES = new Set(["concatenated.nex", "partitions.txt"]);

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
