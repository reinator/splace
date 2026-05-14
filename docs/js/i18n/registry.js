(function () {
    "use strict";

    const languages = ["en", "pt", "es"];
    const store = window.SPLACE_TRANSLATIONS || Object.create(null);

    languages.forEach((lang) => {
        if (!store[lang]) {
            store[lang] = Object.create(null);
        }
    });

    window.SPLACE_TRANSLATIONS = store;
    window.registerSplaceTranslations = function registerSplaceTranslations(bundle) {
        if (!bundle) return store;

        languages.forEach((lang) => {
            if (bundle[lang]) {
                Object.assign(store[lang], bundle[lang]);
            }
        });

        return store;
    };
})();