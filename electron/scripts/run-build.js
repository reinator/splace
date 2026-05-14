const fs = require("node:fs");
const { spawnSync } = require("node:child_process");

const platformFlag = process.argv[2];

if (!platformFlag || !["--win", "--mac", "--linux"].includes(platformFlag)) {
    console.error("Usage: node scripts/run-build.js --win|--mac|--linux");
    process.exit(1);
}

function isWsl() {
    if (process.platform !== "linux") {
        return false;
    }

    try {
        return /microsoft/i.test(fs.readFileSync("/proc/version", "utf8"));
    } catch {
        return false;
    }
}

const env = {
    ...process.env,
    CSC_IDENTITY_AUTO_DISCOVERY: "false",
};

const command = isWsl()
    ? { file: "bash", args: ["build-wsl.sh", platformFlag] }
    : { file: "npx", args: ["electron-builder", platformFlag] };

const result = spawnSync(command.file, command.args, {
    cwd: process.cwd(),
    env,
    stdio: "inherit",
    shell: false,
});

if (result.error) {
    console.error(result.error.message);
    process.exit(1);
}

process.exit(result.status ?? 1);