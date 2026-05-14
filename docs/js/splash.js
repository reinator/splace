(function () {
    "use strict";

    window.SPLACE_T = (key, params) => {
        const api = window.SPLACE_I18N;
        if (api && typeof api.t === "function") {
            return api.t(key, params);
        }

        let value = (window.SPLACE_TRANSLATIONS?.en && window.SPLACE_TRANSLATIONS.en[key]) || key;
        if (params) {
            value = value.replace(/\{(\w+)\}/g, (_, name) => {
                return params[name] == null ? `{${name}}` : String(params[name]);
            });
        }
        return value;
    };
})();