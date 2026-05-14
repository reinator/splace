const { contextBridge, ipcRenderer } = require("electron");

contextBridge.exposeInMainWorld("isElectron", true);

contextBridge.exposeInMainWorld("electronAPI", {
    getCpuCount: () => ipcRenderer.invoke("get-cpu-count"),
    selectGenbankInputs: () => ipcRenderer.invoke("select-genbank-inputs"),

    runAnalysis: (params) => ipcRenderer.send("run-analysis", params),
    runTrimal: (params) => ipcRenderer.send("run-trimal", params),

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
});

