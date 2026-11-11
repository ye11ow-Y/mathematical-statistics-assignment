// Local fallback when native Rscript is unavailable. All statistics run in R.
import fs from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { WebR } from '../tmp/r-runtime/node_modules/webr/dist/webr.mjs';

const root = fileURLToPath(new URL('../', import.meta.url));
const group = process.argv[2] ?? 'q1-q3';
const scripts = { 'q1-q3': 'q1_q3.R', 'q4-q6': 'q4_q6.R', 'q7': 'q7.R' };
if (!Object.hasOwn(scripts, group)) {
  throw new Error('Choose q1-q3, q4-q6 or q7.');
}
const scriptName = scripts[group];
const runtime = path.join(root, 'tmp', 'r-runtime', 'node_modules', 'webr', 'dist') + path.sep;
// webR 0.6.0 resolves a native path before dynamic import. Node on Windows
// requires a file URL here; adapt only this loader expression in the local copy.
const workerPath = path.join(runtime, 'webr-worker.js');
const worker = await fs.readFile(workerPath, 'utf8');
const originalImport = 'import((await import("path")).default.resolve(e))';
const windowsImport = 'import((await import("url")).pathToFileURL((await import("path")).default.resolve(e)).href)';
if (process.platform === 'win32' && worker.includes(originalImport)) {
  await fs.writeFile(workerPath, worker.replace(originalImport, windowsImport));
}
const webR = new WebR({ baseUrl: runtime });
try {
  await webR.init();
  console.log(await webR.evalRString('R.version.string'));
  await webR.FS.mkdir('/work');
  await webR.FS.writeFile(`/work/${scriptName}`, await fs.readFile(path.join(root, 'analysis', scriptName)));
  await webR.evalRVoid(`setwd("/work"); source("${scriptName}")`);
  const files = await webR.evalR(`list.files("output/${group}")`);
  const names = await files.toArray();
  const destination = path.join(root, 'output', group);
  await fs.mkdir(destination, { recursive: true });
  for (const name of names) {
    await fs.writeFile(path.join(destination, name), await webR.FS.readFile(`/work/output/${group}/${name}`));
  }
  console.log(`Exported ${names.length} result files.`);
  console.log(await fs.readFile(path.join(destination, 'calculation_summary.txt'), 'utf8'));
} catch (error) {
  console.error(error.message);
  process.exitCode = 1;
} finally {
  await webR.close();
}
