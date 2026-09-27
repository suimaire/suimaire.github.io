// Release inventory and failure modes that complement the existing functional suites.
import assert from 'node:assert/strict';
import { readFile, readdir, writeFile, mkdir } from 'node:fs/promises';
import { resolve, relative } from 'node:path';
import { createHash } from 'node:crypto';
import { chromium, webkit } from './.pcr-tools/node_modules/playwright/index.mjs';
import { emptyRecord, STORAGE_KEY } from '../assets/js/pcr-records.mjs';

const output = resolve('verification.local/pcr-primer-design/phase9');
await mkdir(output, { recursive: true });
const url = (process.env.PCR_SITE_URL || 'http://127.0.0.1:4174') + '/bioinformatics/pcr-primer-design/';
const html = await readFile('bioinformatics/pcr-primer-design.html', 'utf8');
const css = await readFile('assets/css/pcr-worksheet.css', 'utf8');
const modules = (await readdir('assets/js')).filter(n => /^pcr-.*\.mjs$/.test(n));
const sources = Object.fromEntries(await Promise.all(modules.map(async n => [n, await readFile(`assets/js/${n}`, 'utf8')])));
const graph = Object.fromEntries(Object.entries(sources).map(([name, text]) => [name, [...text.matchAll(/from\s+['"]\.\/([^'"]+)['"]/g)].map(m => m[1])]));
const reachable = new Set();
function visit(name) { if (reachable.has(name)) return; reachable.add(name); for (const child of graph[name]) visit(child); }
visit('pcr-worksheet.mjs');
assert.equal(reachable.size, modules.length);
assert.doesNotMatch(Object.values(sources).join('\n'), /console\.log|debugger;|sourceMappingURL|Playwright|\.pcr-tools|test-only/);
const allText = html + await readFile('_layouts/pcr-worksheet.html', 'utf8') + Object.values(sources).join('\n');
const cssOnlyClasses = [...new Set([...css.matchAll(/\.((?:pcr-|final-|rna-|ext-)[a-z0-9-]+)/g)].map(m => m[1]))].filter(name => !allText.includes(name));
const externalLinks = [...new Set([...html.matchAll(/href="(https:\/\/[^"\s]+)"/g)].map(m => m[1]))];
const build = resolve('tests/.pcr-output/site');
async function tree(directory) { return (await Promise.all((await readdir(directory, { withFileTypes: true })).map(async e => e.isDirectory() ? tree(resolve(directory, e.name)) : [resolve(directory, e.name)]))).flat(); }
const builtPaths = (await tree(build)).map(p => relative(build, p).replaceAll('\\', '/'));
assert.ok(!builtPaths.some(p => /(?:^|\/)(verification\.local|_codex|tests|node_modules)\//.test(p)));
const digests = [];
for (const path of ['assets/css/pcr-worksheet.css', 'assets/data/pcr-primer-fixture.json', ...modules.map(n => 'assets/js/' + n)]) {
  const source = await readFile(path), built = await readFile(resolve(build, path));
  assert.ok(source.equals(built), `built asset equals reviewed source: ${path}`);
  digests.push({ path, sha256: createHash('sha256').update(built).digest('hex'), bytes: built.length });
}
const builtHtml = await readFile(resolve(build, 'bioinformatics/pcr-primer-design/index.html'), 'utf8');
assert.doesNotMatch(builtHtml, /\{\{|\{%/);
const failures = [], headings = [];
for (const [engineName, engine] of Object.entries({ chromium, webkit })) {
  const browser = await engine.launch({ headless: true });
  try {
    for (const mode of ['no-js', 'module-failure', 'quota', 'malformed', 'unsupported-answer', 'invalid-coordinates', 'partial-rna']) {
      const context = await browser.newContext({ javaScriptEnabled: mode !== 'no-js' });
      const page = await context.newPage();
      if (mode === 'module-failure') await page.route('**/pcr-worksheet.mjs', r => r.fulfill({ status: 503, body: 'QA injected module failure' }));
      if (mode === 'quota') await context.addInitScript(() => { Storage.prototype.setItem = function () { throw new DOMException('QA quota', 'QuotaExceededError'); }; });
      const record = emptyRecord();
      if (mode === 'unsupported-answer') record.answers['unsupported-field'] = 'unsupported';
      if (mode === 'invalid-coordinates') record.draft = { forward: 'ACGT', reverse: 'ACGT', bindings: { F: { start: '-9', end: '9999', direction: 'right' }, R: null } };
      if (mode === 'partial-rna') record.rnaExtension = { activeSection: 'strategy', answers: { strategy: '부분 기록 유지' } };
      if (['malformed', 'unsupported-answer', 'invalid-coordinates', 'partial-rna'].includes(mode)) await context.addInitScript(({ key, raw }) => localStorage.setItem(key, raw), { key: STORAGE_KEY, raw: mode === 'malformed' ? '{broken' : JSON.stringify(record) });
      await page.goto(url);
      if (['no-js', 'module-failure'].includes(mode)) {
        assert.match(await page.locator('#save-status').innerText(), /자동 저장과 계산을 사용할 수 없습니다/);
        assert.equal(await page.locator('#pcr-main > section > h2, #ext-introduction h2').count(), 8);
      } else {
        await page.waitForSelector('#pcr-worksheet[data-ready=true]');
        if (['malformed', 'unsupported-answer'].includes(mode)) {
          const before = await page.evaluate(key => localStorage.getItem(key), STORAGE_KEY);
          assert.match(await page.locator('#save-status').innerText(), /이전 기록을 읽을 수 없습니다/);
          await page.locator('#first-placement').fill('현재 입력은 계속 작성 가능');
          assert.equal(await page.evaluate(key => localStorage.getItem(key), STORAGE_KEY), before);
        }
        if (mode === 'quota') { await page.locator('#first-placement').fill('저장 공간 부족 시 내보낼 기록'); assert.match(await page.locator('#save-status').innerText(), /공간이 부족/); }
        if (mode === 'invalid-coordinates') { assert.match(await page.locator('#analysis-status').innerText(), /좌표/); assert.equal(await page.locator('#primer-f').inputValue(), 'ACGT'); }
        if (mode === 'partial-rna') {
          await page.locator('#rna-extension > summary').click(); assert.equal(await page.locator('#rna-strategy-answer').inputValue(), '부분 기록 유지');
          const list = await page.locator('h1,h2,h3,h4,h5,h6').evaluateAll(els => els.filter(e => e.getClientRects().length).map(e => ({ level: Number(e.tagName.slice(1)), text: e.textContent.trim() })));
          assert.equal(list.filter(h => h.level === 1).length, 1);
          for (let i = 1; i < list.length; i++) assert.ok(list[i].level <= list[i - 1].level + 1, `heading level: ${list[i].text}`);
          headings.push({ engineName, list });
        }
        // Export remains usable during failure, without writing any production state.
        await page.locator('#record-menu > summary').click(); const download = page.waitForEvent('download'); await page.locator('#export-record').click(); await download;
      }
      failures.push({ engineName, mode, pass: true }); await context.close();
    }
  } finally { await browser.close(); }
}
await writeFile(resolve(output, 'source-and-failure-audit.json'), JSON.stringify({ url, modules: graph, unreachableModules: modules.filter(n => !reachable.has(n)), cssOnlyClasses, externalLinks, builtPaths, digests, headings, failures }, null, 2));
console.log(`PASS release audit: ${modules.length} reachable modules; build exclusions/assets; ${failures.length} failure/restore checks; heading hierarchy`);
