#!/usr/bin/env node
/**
 * CI helper — stamp package.json with a commit-derived pre-release version.
 *
 * Takes the version already declared in package.json (e.g. "0.9.6") and appends
 * the first 7 characters of the commit being built, producing a semver
 * pre-release version such as "0.9.6-abc1234".
 *
 * vsce names the resulting package "<name>-<version>.vsix", so the artifact
 * automatically becomes "Molecular-Viewer-0.9.6-abc1234.vsix".
 *
 * Derived metadata is written as KEY=VALUE lines to $GITHUB_OUTPUT so later
 * workflow steps can consume it, and echoed to the log for readability.
 */

const fs = require('fs');
const path = require('path');

const sha = process.env.GITHUB_SHA;
if (!sha) {
  console.error('ERROR: GITHUB_SHA is not set.');
  process.exit(1);
}

// First 7 characters of the commit hash (short SHA).
const sha7 = sha.slice(0, 7);

const pkgPath = path.join(__dirname, '..', '..', 'package.json');
const pkg = JSON.parse(fs.readFileSync(pkgPath, 'utf8'));

// Strip any pre-release suffix that may already be present so repeated runs
// always start from the clean base version.
const baseVersion = String(pkg.version).split('-')[0];
const version = `${baseVersion}-${sha7}`;

pkg.version = version;
fs.writeFileSync(pkgPath, JSON.stringify(pkg, null, 2) + '\n');

const outputs = {
  sha7,
  version,
  tag: `v${version}`,
  vsix: `${pkg.name}-${version}.vsix`,
};

const lines = Object.entries(outputs)
  .map(([key, value]) => `${key}=${value}`)
  .join('\n');

if (process.env.GITHUB_OUTPUT) {
  fs.appendFileSync(process.env.GITHUB_OUTPUT, `${lines}\n`);
}

console.log(lines);
