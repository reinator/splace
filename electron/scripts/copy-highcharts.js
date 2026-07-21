const fs = require("fs");
const path = require("path");

const electronRoot = path.resolve(__dirname, "..");
const repoRoot = path.resolve(electronRoot, "..");
const highchartsRoot = path.join(electronRoot, "node_modules", "highcharts");
const docsVendorRoot = path.join(repoRoot, "docs", "vendor", "highcharts");

const requiredFiles = [
  "highcharts.js",
  path.join("modules", "heatmap.js"),
  path.join("modules", "exporting.js"),
  path.join("modules", "export-data.js"),
  path.join("modules", "accessibility.js"),
];

function copyFile(relativePath) {
  const source = path.join(highchartsRoot, relativePath);
  const target = path.join(docsVendorRoot, relativePath);

  if (!fs.existsSync(source)) {
    throw new Error(`Missing Highcharts asset: ${source}`);
  }

  fs.mkdirSync(path.dirname(target), { recursive: true });
  fs.copyFileSync(source, target);
  console.log(`Copied: ${relativePath.replace(/\\/g, "/")}`);
}

function main() {
  if (!fs.existsSync(highchartsRoot)) {
    console.error(`Highcharts was not found at: ${highchartsRoot}`);
    console.error("Run npm ci inside the electron folder before building.");
    process.exit(1);
  }

  for (const file of requiredFiles) {
    copyFile(file);
  }

  console.log("Highcharts files copied successfully.");
}

main();
