/* =========================================================================
   SPLACE — UI Chrome
   Splash modal, mini-map, logs drawer, settings popover,
   bilingual i18n, theme/density/font-scale, step locking.
   Loaded BEFORE splace.js so it patches console + DOM logging hooks.
   ========================================================================= */

(function () {
    "use strict";

    // ------------------------------------------------------------------
    // Storage helpers
    // ------------------------------------------------------------------
    const LS = {
        get(k, fallback) {
            try { const v = localStorage.getItem(k); return v == null ? fallback : v; }
            catch { return fallback; }
        },
        set(k, v) { try { localStorage.setItem(k, v); } catch { } },
    };

    // ------------------------------------------------------------------
    // i18n — dictionaries are registered from js/i18n/*.js.
    // English is the default fallback.
    // ------------------------------------------------------------------
    const I18N = window.SPLACE_TRANSLATIONS || { en: {}, pt: {}, es: {} };

    let currentLang = LS.get("splace.lang", "en");

    function t(key, params) {
        let value = (I18N[currentLang] && I18N[currentLang][key]) || I18N.en[key] || key;
        if (params) {
            value = value.replace(/\{(\w+)\}/g, (_, name) => {
                return params[name] == null ? `{${name}}` : String(params[name]);
            });
        }
        return value;
    }
    function applyI18n(root) {
        (root || document).querySelectorAll("[data-i18n]").forEach(el => {
            const key = el.getAttribute("data-i18n");
            const val = t(key);
            if (el.dataset.i18nAttr) {
                // Allow safe HTML tags like <b>, <strong>, <em>, <i>, <br>
                el.setAttribute(el.dataset.i18nAttr, val);
            } else {
                el.innerHTML = val;
            }
        });
        (root || document).querySelectorAll("[data-i18n-tooltip]").forEach(el => {
            const key = el.getAttribute("data-i18n-tooltip");
            el.setAttribute("data-tippy-content", t(key));
        });
        document.documentElement.lang = currentLang;
    }
    window.SPLACE_I18N = {
        t,
        get lang() { return currentLang; },
        get translations() { return I18N; },
        applyI18n,
    };

    function setLanguage(lang) {
        currentLang = I18N[lang] ? lang : "en";
        LS.set("splace.lang", currentLang);
        document.querySelectorAll(".lang-toggle button").forEach(b => {
            b.classList.toggle("active", b.dataset.lang === currentLang);
        });
        applyI18n();
        if (typeof window.initTooltips === "function") {
            window.initTooltips();
        }
        if (typeof window.refreshStep1UiTranslations === "function") {
            window.refreshStep1UiTranslations();
        }
        if (typeof window.refreshRecordsUiTranslations === "function") {
            window.refreshRecordsUiTranslations();
        }
        if (typeof window.refreshAlignmentUiTranslations === "function") {
            window.refreshAlignmentUiTranslations();
        }

        const overlay = document.getElementById("splashOverlay");
        const startButtons = Array.from(document.querySelectorAll("[data-splash-start]"));
        if (overlay && startButtons.length && overlay.style.display !== "none" && !splashTimer) {
            startButtons.forEach((startBtn) => {
                const labelKey = startBtn.dataset.splashStart === "desktop" ? "splash.startDesktop" : "splash.start";
                startBtn.innerHTML = `<span>${t(labelKey)}</span> <i class="fa-solid fa-arrow-right"></i>`;
            });
        }
    }
    window.setLanguage = setLanguage;

    // ------------------------------------------------------------------
    // Theme / density / font scale
    // ------------------------------------------------------------------
    function applyTheme(mode) {
        // mode: "light" | "dark" | "auto"
        let resolved = mode;
        if (mode === "auto") {
            resolved = matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light";
        }
        document.documentElement.setAttribute("data-theme", resolved);
        LS.set("splace.theme", mode);
        document.querySelectorAll('.seg[data-seg="theme"] button').forEach(b => {
            b.classList.toggle("active", b.dataset.val === mode);
        });
    }
    function applyDensity(d) {
        document.documentElement.setAttribute("data-density", d);
        LS.set("splace.density", d);
        document.querySelectorAll('.seg[data-seg="density"] button').forEach(b => {
            b.classList.toggle("active", b.dataset.val === d);
        });
    }
    function applyFontSize(px) {
        const v = Math.max(13, Math.min(20, +px));
        document.documentElement.style.setProperty("--font-base", v + "px");
        LS.set("splace.fontsize", v);
        const display = document.querySelector("#fontSizeValue");
        if (display) display.textContent = v + "px";
        const slider = document.querySelector("#fontSizeSlider");
        if (slider && +slider.value !== v) slider.value = v;
    }

    // Restore saved prefs immediately to avoid FOUC
    applyTheme(LS.get("splace.theme", "light"));
    applyDensity(LS.get("splace.density", "comfortable"));
    applyFontSize(LS.get("splace.fontsize", 16));

    // React to OS theme changes when in auto
    matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => {
        if (LS.get("splace.theme") === "auto") applyTheme("auto");
    });

    // ------------------------------------------------------------------
    // Logs panel — captures console + a public window.splaceLog API
    // ------------------------------------------------------------------
    const logBuffer = [];
    let autoScroll = true;

    function pad(n) { return n < 10 ? "0" + n : "" + n; }
    function nowStamp() {
        const d = new Date();
        return pad(d.getHours()) + ":" + pad(d.getMinutes()) + ":" + pad(d.getSeconds());
    }

    function addLog(level, msg, opts) {
        opts = opts || {};
        const entry = {
            time: nowStamp(),
            iso: new Date().toISOString(),
            level: level || "info",
            msg: typeof msg === "string" ? msg : safeStringify(msg),
        };
        logBuffer.push(entry);
        // cap buffer to last 5000
        if (logBuffer.length > 5000) logBuffer.shift();
        renderLogLine(entry);
        updateLogBadge();
    }
    function safeStringify(v) {
        try {
            if (v instanceof Error) return v.message + (v.stack ? "\n" + v.stack : "");
            return JSON.stringify(v, null, 2);
        } catch { return String(v); }
    }
    function renderLogLine(entry) {
        const body = document.getElementById("logsBody");
        if (!body) return;
        const empty = body.querySelector(".logs-empty");
        if (empty) empty.remove();
        const icons = { info: "fa-circle-info", warn: "fa-triangle-exclamation", error: "fa-circle-xmark", success: "fa-circle-check" };
        const div = document.createElement("div");
        div.className = "log-line " + entry.level;
        div.innerHTML =
            `<span class="log-time">${entry.time}</span>` +
            `<span class="log-icon"><i class="fa-solid ${icons[entry.level] || icons.info}"></i></span>` +
            `<span class="log-msg"></span>`;
        div.querySelector(".log-msg").textContent = entry.msg;
        body.appendChild(div);
        if (autoScroll) body.scrollTop = body.scrollHeight;
    }
    function clearLogs() {
        logBuffer.length = 0;
        const body = document.getElementById("logsBody");
        if (!body) return;
        body.innerHTML = `<div class="logs-empty"><i class="fa-solid fa-bars-staggered"></i><p data-i18n="logs.empty"></p></div>`;
        applyI18n(body);
        updateLogBadge();
    }
    function copyLogs() {
        const text = logBuffer.map(e => `[${e.iso}] ${e.level.toUpperCase()} ${e.msg}`).join("\n");
        if (!text) return;
        navigator.clipboard?.writeText(text).then(() => flashTooltip("logsCopyBtn", "✓"));
    }
    function downloadLogs() {
        const text = logBuffer.map(e => `[${e.iso}] ${e.level.toUpperCase()} ${e.msg}`).join("\n");
        const blob = new Blob([text], { type: "text/plain" });
        const url = URL.createObjectURL(blob);
        const a = document.createElement("a");
        a.href = url;
        a.download = `splace_log_${new Date().toISOString().slice(0, 19).replace(/[:T]/g, "-")}.txt`;
        document.body.appendChild(a); a.click(); a.remove();
        setTimeout(() => URL.revokeObjectURL(url), 1000);
    }
    function flashTooltip(id, text) {
        const el = document.getElementById(id);
        if (!el) return;
        const orig = el.innerHTML;
        el.innerHTML = text;
        setTimeout(() => { el.innerHTML = orig; }, 900);
    }
    function updateLogBadge() {
        const badge = document.getElementById("logsBadge");
        if (!badge) return;
        const errors = logBuffer.filter(e => e.level === "error" || e.level === "warn").length;
        if (errors > 0 && !isLogsOpen()) {
            badge.textContent = errors > 99 ? "99+" : errors;
            badge.classList.remove("hidden");
        } else {
            badge.classList.add("hidden");
        }
    }
    function isLogsOpen() {
        return document.querySelector(".app-shell")?.dataset.logs === "open";
    }

    // Expose
    window.splaceLog = {
        info: (m) => addLog("info", m),
        warn: (m) => addLog("warn", m),
        error: (m) => addLog("error", m),
        success: (m) => addLog("success", m),
        clear: clearLogs,
    };

    // Patch console — preserve original output, mirror to drawer
    const origConsole = {
        log: console.log.bind(console),
        info: console.info.bind(console),
        warn: console.warn.bind(console),
        error: console.error.bind(console),
    };
    console.log = function (...args) {
        origConsole.log(...args);
        if (args.length) addLog("info", args.map(a => typeof a === "string" ? a : safeStringify(a)).join(" "));
    };
    console.info = function (...args) {
        origConsole.info(...args);
        if (args.length) addLog("info", args.map(a => typeof a === "string" ? a : safeStringify(a)).join(" "));
    };
    console.warn = function (...args) {
        origConsole.warn(...args);
        if (args.length) addLog("warn", args.map(a => typeof a === "string" ? a : safeStringify(a)).join(" "));
    };
    console.error = function (...args) {
        origConsole.error(...args);
        if (args.length) addLog("error", args.map(a => typeof a === "string" ? a : safeStringify(a)).join(" "));
    };
    window.addEventListener("error", (e) => addLog("error", `${e.message} (${e.filename}:${e.lineno})`));
    window.addEventListener("unhandledrejection", (e) => addLog("error", "Unhandled: " + (e.reason?.message || e.reason)));

    // ------------------------------------------------------------------
    // Mini-map / step navigator
    // ------------------------------------------------------------------
    // Step config — sectionId is what we scroll to and observe.
    // The unlock predicate is what determines whether the step is reachable.
    const STEPS = [
        { n: 1, key: "step.1", section: "step1Section", unlock: () => true },
        { n: 2, key: "step.2", section: "recordsSection", unlock: () => isVisible("recordsSection") },
        { n: 3, key: "step.3", section: "featureTypesSection", unlock: () => isVisible("featureTypesSection") },
        { n: 4, key: "step.4", section: "genesSection", unlock: () => isVisible("genesSection") },
        {
            n: 5, keyWeb: "step.5.web", keyDesktop: "step.5.desktop",
            section: () => window.isElectron ? "nextStepSection" : "downloadSection",
            unlock: () => isVisible(window.isElectron ? "nextStepSection" : "downloadSection")
        },
        {
            n: 6, key: "step.6", section: "concatSection",
            unlock: () => window.isElectron && !document.getElementById("concatSection")?.classList.contains("hidden"),
            desktopOnly: true
        },
        {
            n: 7, key: "step.7", section: "iqtreeSection",
            unlock: () => window.isElectron && !document.getElementById("iqtreeSection")?.classList.contains("hidden"),
            desktopOnly: true
        },
    ];

    function state() { return window.state || {}; }
    function isVisible(id) {
        const el = document.getElementById(id);
        return !!(el && !el.classList.contains("hidden"));
    }
    function hasProceeded() {
        const ft = document.getElementById("featureTypesSection");
        return ft && !ft.classList.contains("hidden");
    }

    function resolveSection(step) {
        const id = typeof step.section === "function" ? step.section() : step.section;
        return document.getElementById(id);
    }

    function buildMiniMap() {
        const list = document.getElementById("mmList");
        if (!list) return;
        list.innerHTML = "";
        STEPS.forEach((s) => {
            if (s.desktopOnly && !window.isElectron) return;
            const a = document.createElement("a");
            a.className = "mm-step";
            a.dataset.step = s.n;
            const labelKey = s.key || (window.isElectron ? s.keyDesktop : s.keyWeb);
            a.innerHTML = `
                <span class="mm-num"><span>${s.n}</span><i class="fa-solid fa-check"></i></span>
                <span class="mm-label" data-i18n="${labelKey}">${t(labelKey)}</span>
                <span class="mm-meta" data-mm-meta></span>
            `;
            a.addEventListener("click", (e) => {
                e.preventDefault();
                navigateToStep(s.n);
            });
            list.appendChild(a);
        });
        applyI18n(list);
    }

    function navigateToStep(n) {
        const s = STEPS.find(x => x.n === n);
        if (!s) return;
        const el = resolveSection(s);
        const unlocked = s.unlock();
        if (!unlocked || !el || el.classList.contains("hidden")) {
            // brief shake on locked
            const item = document.querySelector(`.mm-step[data-step="${n}"]`);
            if (item) {
                item.animate([
                    { transform: "translateX(0)" }, { transform: "translateX(-4px)" },
                    { transform: "translateX(4px)" }, { transform: "translateX(0)" },
                ], { duration: 220 });
            }
            return;
        }
        el.scrollIntoView({ behavior: "smooth", block: "start" });
        // close minimap drawer on mobile
        if (window.matchMedia("(max-width: 900px)").matches) closeMiniMap();
    }
    window.splaceNavigate = navigateToStep;

    function refreshMiniMap() {
        let activeFound = false;
        let doneCount = 0;
        const visibleSteps = STEPS.filter(s => !(s.desktopOnly && !window.isElectron));
        const total = visibleSteps.length;

        // Find topmost visible section
        const scrollTop = window.scrollY + 120;

        visibleSteps.forEach((s, idx) => {
            const item = document.querySelector(`.mm-step[data-step="${s.n}"]`);
            if (!item) return;
            const el = resolveSection(s);
            const unlocked = s.unlock();
            const visible = el && !el.classList.contains("hidden");

            item.classList.remove("active", "locked", "done");

            if (!unlocked || !visible) {
                item.classList.add("locked");
                item.setAttribute("aria-disabled", "true");
                item.setAttribute("title", t("minimap.locked"));
            } else {
                item.removeAttribute("aria-disabled");
                item.removeAttribute("title");
                // "done" if next step's section is also visible (we've moved past)
                const next = visibleSteps[idx + 1];
                const nextVisible = next && !resolveSection(next)?.classList.contains("hidden") && next.unlock();
                if (nextVisible) {
                    item.classList.add("done");
                    doneCount++;
                }
            }

            // active: section in viewport
            if (visible && el) {
                const rect = el.getBoundingClientRect();
                const offset = rect.top + window.scrollY;
                if (offset <= scrollTop && !activeFound) {
                    // mark previous one as active by checking next
                    const nextEl = visibleSteps[idx + 1] ? resolveSection(visibleSteps[idx + 1]) : null;
                    const nextOffset = nextEl ? nextEl.getBoundingClientRect().top + window.scrollY : Infinity;
                    if (scrollTop < nextOffset) {
                        item.classList.add("active");
                        activeFound = true;
                    }
                }
            }
        });

        // Progress bar
        const finalPct = total > 0 ? Math.round((doneCount / total) * 100) : 0;
        const bar = document.getElementById("mmProgress");
        if (bar) bar.style.width = finalPct + "%";
        const pct = document.getElementById("mmProgressPct");
        if (pct) pct.textContent = finalPct + "%";

        // Apply locking visually to step sections
        applyStepLocking(visibleSteps);
    }

    function applyStepLocking(visibleSteps) {
        visibleSteps.forEach((s, idx) => {
            const el = resolveSection(s);
            if (!el || el.classList.contains("hidden")) return;
            el.classList.add("step-section");
            // Active highlight: closest to top
            const item = document.querySelector(`.mm-step[data-step="${s.n}"]`);
            el.classList.toggle("is-active", item?.classList.contains("active"));
            // We don't visually grey-out a step; locked steps are simply hidden by the original JS.
            // But if a section is rendered AND its unlock is false (rare race), grey it.
            const unlocked = s.unlock();
            el.classList.toggle("locked", !unlocked);
        });
    }

    // ------------------------------------------------------------------
    // Drawer toggles
    // ------------------------------------------------------------------
    function shell() { return document.querySelector(".app-shell"); }
    function setLogs(open) {
        const sh = shell(); if (!sh) return;
        sh.dataset.logs = open ? "open" : "closed";
        document.getElementById("logsToggleBtn")?.classList.toggle("active", open);
        LS.set("splace.logs", open ? "open" : "closed");
        if (open) updateLogBadge();
    }
    function setMiniMap(open) {
        const sh = shell(); if (!sh) return;
        sh.dataset.minimap = open ? "open" : "closed";
        document.getElementById("menuToggleBtn")?.classList.toggle("active", open);
        LS.set("splace.minimap", open ? "open" : "closed");
    }
    function closeMiniMap() { setMiniMap(false); }

    function setSplashScrollLock(locked) {
        document.documentElement.classList.toggle("splash-open", locked);
        document.body.classList.toggle("splash-open", locked);
    }

    // ------------------------------------------------------------------
    // Splash modal
    // ------------------------------------------------------------------
    const SPLASH_FIRST_DELAY = 15; // seconds before button appears on first load
    let splashTimer = null;

    function syncSplashPlatformLayout() {
        const isDesktop = !!window.isElectron;
        const overlay = document.getElementById("splashOverlay");
        const webCard = document.querySelector(".splash-platform-web");
        const webStart = document.querySelector("[data-splash-start='web']");
        const desktopStart = document.querySelector("[data-splash-start='desktop']");
        const desktopDownload = document.querySelector("[data-splash-download='desktop']");

        if (overlay) overlay.classList.toggle("splash-desktop-mode", isDesktop);
        if (webCard) webCard.classList.toggle("hidden", isDesktop);
        if (webStart) webStart.classList.toggle("hidden", isDesktop);
        if (desktopStart) desktopStart.classList.toggle("hidden", !isDesktop);
        if (desktopDownload) desktopDownload.classList.toggle("hidden", isDesktop);
    }

    function showSplash(opts) {
        const isFirst = !!(opts && opts.firstLoad);
        const isDesktop = !!window.isElectron;
        const overlay = document.getElementById("splashOverlay");
        if (!overlay) return;
        setSplashScrollLock(true);
        overlay.style.display = "flex";
        syncSplashPlatformLayout();
        applyI18n(overlay);
        if (typeof window.initTooltips === "function") {
            window.initTooltips();
        }

        const startButtons = Array.from(document.querySelectorAll("[data-splash-start]"));
        if (splashTimer) clearInterval(splashTimer);

        // Se NÃO for a primeira vez (ex: clicou no botão "Help"), libera o botão imediatamente
        if (!isFirst || isDesktop) {
            startButtons.forEach((startBtn) => {
                const labelKey = startBtn.dataset.splashStart === "desktop" ? "splash.startDesktop" : "splash.start";
                startBtn.disabled = false;
                startBtn.classList.remove("splash-action-btn-waiting");
                startBtn.innerHTML = `<span>${t(labelKey)}</span> <i class="fa-solid fa-arrow-right"></i>`;
            });
            return;
        }

        // Se FOR a primeira vez, bloqueia e conta
        startButtons.forEach((startBtn) => { startBtn.disabled = true; });

        let remaining = SPLASH_FIRST_DELAY; // Constante que você já tem (15 segundos)

        const render = () => {
            startButtons.forEach((startBtn) => {
                const timerKey = startBtn.dataset.splashStart === "desktop" ? "splash.timerDesktop" : "splash.timer";
                const labelKey = startBtn.dataset.splashStart === "desktop" ? "splash.startDesktop" : "splash.start";
                const isWebStart = startBtn.dataset.splashStart === "web";
                if (remaining > 0) {
                    startBtn.classList.toggle("splash-action-btn-waiting", isWebStart);
                    startBtn.innerHTML = `<span>${t(timerKey, { n: remaining })}</span> <i class="fa-solid fa-hourglass-half"></i>`;
                } else {
                    startBtn.classList.remove("splash-action-btn-waiting");
                    startBtn.innerHTML = `<span>${t(labelKey)}</span> <i class="fa-solid fa-arrow-right"></i>`;
                }
            });
        };

        render(); // Renderiza estado inicial (15s)

        splashTimer = setInterval(() => {
            remaining--;
            render();
            if (remaining <= 0) {
                clearInterval(splashTimer);
                splashTimer = null;
                startButtons.forEach((startBtn) => { startBtn.disabled = false; });
            }
        }, 1000);
    }

    function hideSplash() {
        if (splashTimer) clearInterval(splashTimer);
        splashTimer = null;
        const overlay = document.getElementById("splashOverlay");
        if (!overlay) return;
        overlay.style.opacity = "0";
        setTimeout(() => {
            overlay.style.display = "none";
            overlay.style.opacity = "";
            setSplashScrollLock(false);
        }, 200);
    }
    window.SPLACE_showSplash = showSplash;

    // ------------------------------------------------------------------
    // Init on DOMContentLoaded
    // ------------------------------------------------------------------
    function init() {
        document.documentElement.dataset.platform = window.isElectron ? "desktop" : "web";
        syncSplashPlatformLayout();

        // Apply translations to chrome
        applyI18n();

        // Lang toggle
        document.querySelectorAll(".lang-toggle button").forEach(b => {
            b.classList.toggle("active", b.dataset.lang === currentLang);
            b.addEventListener("click", () => setLanguage(b.dataset.lang));
        });

        // Settings popover
        const settingsBtn = document.getElementById("settingsBtn");
        const popover = document.getElementById("settingsPopover");
        if (settingsBtn && popover) {
            settingsBtn.addEventListener("click", (e) => {
                e.stopPropagation();
                popover.classList.toggle("open");
                settingsBtn.classList.toggle("active", popover.classList.contains("open"));
            });
            document.addEventListener("click", (e) => {
                if (popover.classList.contains("open") && !popover.contains(e.target) && e.target !== settingsBtn) {
                    popover.classList.remove("open");
                    settingsBtn.classList.remove("active");
                }
            });
        }
        // Theme seg
        document.querySelectorAll('.seg[data-seg="theme"] button').forEach(b => {
            b.classList.toggle("active", b.dataset.val === LS.get("splace.theme", "light"));
            b.addEventListener("click", () => applyTheme(b.dataset.val));
        });
        document.querySelectorAll('.seg[data-seg="density"] button').forEach(b => {
            b.classList.toggle("active", b.dataset.val === LS.get("splace.density", "comfortable"));
            b.addEventListener("click", () => applyDensity(b.dataset.val));
        });
        const slider = document.getElementById("fontSizeSlider");
        if (slider) {
            slider.value = LS.get("splace.fontsize", 16);
            slider.addEventListener("input", () => applyFontSize(slider.value));
        }
        const resetBtn = document.getElementById("fontSizeResetBtn");
        if (resetBtn) {
            resetBtn.addEventListener("click", () => applyFontSize(16));
        }
        applyFontSize(LS.get("splace.fontsize", 16));

        // Drawer toggles
        document.getElementById("menuToggleBtn")?.addEventListener("click", () => {
            const open = shell().dataset.minimap !== "open";
            setMiniMap(open);
        });
        document.getElementById("logsToggleBtn")?.addEventListener("click", () => {
            const open = shell().dataset.logs !== "open";
            setLogs(open);
        });
        document.getElementById("logsCollapseBtn")?.addEventListener("click", () => setLogs(false));
        document.getElementById("logsCopyBtn")?.addEventListener("click", copyLogs);
        document.getElementById("logsDownloadBtn")?.addEventListener("click", downloadLogs);
        document.getElementById("logsClearBtn")?.addEventListener("click", clearLogs);
        document.getElementById("logsBody")?.addEventListener("scroll", (e) => {
            const el = e.target;
            autoScroll = (el.scrollTop + el.clientHeight >= el.scrollHeight - 8);
        });

        // Help button — re-open splash (button visible immediately)
        document.getElementById("helpBtn")?.addEventListener("click", () => showSplash({ firstLoad: false }));

        // Splash actions
        document.querySelectorAll("[data-splash-start]").forEach((btn) => {
            btn.addEventListener("click", hideSplash);
        });
        document.getElementById("splashOverlay")?.addEventListener("click", (e) => {
            if (e.target.id === "splashOverlay") hideSplash();
        });
        document.querySelectorAll(".splash-lang button").forEach(b => {
            b.addEventListener("click", () => {
                setLanguage(b.dataset.lang);
                document.querySelectorAll(".splash-lang button").forEach(x => x.classList.toggle("active", x.dataset.lang === currentLang));
            });
        });
        document.querySelectorAll(".splash-lang button").forEach(b => {
            b.classList.toggle("active", b.dataset.lang === currentLang);
        });

        // Mobile tabs (only relevant <900px)
        document.querySelectorAll(".mobile-tabs button").forEach(b => {
            b.addEventListener("click", () => {
                const target = b.dataset.tab;
                document.querySelectorAll(".mobile-tabs button").forEach(x => x.classList.toggle("active", x === b));
                if (target === "logs") { setLogs(true); }
                else if (target === "records") {
                    const el = document.getElementById("recordsSection");
                    if (el && !el.classList.contains("hidden")) el.scrollIntoView({ behavior: "smooth" });
                } else if (target === "heatmap") {
                    const el = document.getElementById("heatmapSection");
                    if (el && !el.classList.contains("hidden")) el.scrollIntoView({ behavior: "smooth" });
                } else {
                    window.scrollTo({ top: 0, behavior: "smooth" });
                }
            });
        });

        // Build mini-map
        buildMiniMap();

        // Restore drawer states
        const isMobile = matchMedia("(max-width: 900px)").matches;
        const isNarrow = matchMedia("(max-width: 1280px)").matches;
        setMiniMap(LS.get("splace.minimap", isMobile ? "closed" : "open") === "open");
        setLogs(LS.get("splace.logs", isNarrow ? "closed" : "open") === "open");

        // Welcome log
        addLog("info", t("logs.welcome"));

        // Refresh mini-map on scroll & on DOM mutations
        let rafId = null;
        const scheduleRefresh = () => {
            if (rafId) cancelAnimationFrame(rafId);
            rafId = requestAnimationFrame(refreshMiniMap);
        };
        window.addEventListener("scroll", scheduleRefresh, { passive: true });
        window.addEventListener("resize", scheduleRefresh);
        const mo = new MutationObserver(scheduleRefresh);
        mo.observe(document.body, { attributes: true, attributeFilter: ["class"], subtree: true });
        // Also poll occasionally for state changes
        setInterval(refreshMiniMap, 800);
        refreshMiniMap();

        // Always show splash on load (with 15s delay before close button)
        setTimeout(() => showSplash({ firstLoad: true }), 80);

        // Hook into proceedToAnalysis if it exists later — re-render minimap
        const oldProceed = window.proceedToAnalysis;
        window.proceedToAnalysis = function () {
            if (typeof oldProceed === "function") oldProceed.apply(this, arguments);
            scheduleRefresh();
        };
    }

    if (document.readyState === "loading") {
        document.addEventListener("DOMContentLoaded", init);
    } else {
        init();
    }
})();
