// Run against real Jekyll output. All production-origin requests below are intercepted.
import assert from 'node:assert/strict';
import { mkdir, readFile, writeFile } from 'node:fs/promises';
import { resolve } from 'node:path';
import { chromium, webkit } from './.pcr-tools/node_modules/playwright/index.mjs';
import { PAGE_VIEWS_CONFIG } from '../assets/js/page-views.js';
import { STORAGE_KEY } from '../assets/js/pcr-records.mjs';

const base = process.env.PCR_SITE_URL || 'http://127.0.0.1:4174';
const path = '/bioinformatics/pcr-primer-design/';
const output = resolve('verification.local/pcr-primer-design/site-integration');
const previewMode = process.env.PCR_COUNTER_PREVIEW === 'live' ? 'live' : 'mock';
await mkdir(output, { recursive: true });
const results = [], consoleErrors = [], warnings = [], pageErrors = [], rejections = [];
const screenshots = [];
let checks = 0;
function check(value, message) { assert.ok(value, message); checks++; }
async function instrument(context) {
  await context.exposeBinding('siteRejection', (_, reason) => rejections.push(reason));
  await context.addInitScript(() => window.addEventListener('unhandledrejection', e => window.siteRejection(String(e.reason))));
  context.on('page', page => {
    page.on('pageerror', e => pageErrors.push(e.message));
    page.on('console', m => { if (m.type() === 'error') consoleErrors.push(m.text()); if (m.type() === 'warning') warnings.push(m.text()); });
  });
}
const ready = page => page.waitForSelector('#pcr-worksheet[data-ready=true]');
const counterReady = page => page.waitForSelector('[data-page-views] .page-views__total');
async function shot(page, name) {
  await page.evaluate(() => document.fonts.ready);
  const filename = resolve(output, name);
  await page.screenshot({ path: filename }); screenshots.push(filename);
}
async function smoke(page, label) {
  await ready(page);
  await page.locator('#first-placement').fill(label);
  const saved = JSON.parse(await page.evaluate(key => localStorage.getItem(key), STORAGE_KEY));
  check(saved.answers['first-placement'] === label, 'counter independent from worksheet autosave');
  await page.locator('#coordinate-details > summary').click();
  await page.locator('#range-start').fill('41'); await page.locator('#range-end').fill('60'); await page.locator('#apply-range').click();
  await page.locator('#range-primer').selectOption('R'); await page.locator('#range-direction').selectOption('left');
  await page.locator('#range-start').fill('281'); await page.locator('#range-end').fill('300'); await page.locator('#apply-range').click();
  await page.locator('#manual-sequences > summary').click();
  await page.locator('#analyze-design').click();
  check((await page.locator('#analysis-results').textContent()).includes('260 bp'), 'counter independent from PCR calculation');
  await page.locator('#rna-extension > summary').click();
  await page.locator('[data-rna-section=reflection]').click();
  await page.locator('#rna-reflection-answer').fill(label + ' RNA');
  await page.locator('#record-menu > summary').click();
  const pending = page.waitForEvent('download'); await page.locator('#export-record').click();
  const download = await pending, file = resolve(output, 'failure-roundtrip.json'); await download.saveAs(file);
  const record = JSON.parse(await readFile(file, 'utf8'));
  check(record.rnaExtension.answers.reflection === label + ' RNA', 'RNA and JSON export independent from counter');
}

for (const [engineName, engine] of Object.entries({ chromium, webkit })) {
  const browser = await engine.launch({ headless: true });
  try {
    for (const width of [1440, 768, 390]) {
      const context = await browser.newContext({ viewport: { width, height: 1000 } });
      await instrument(context);
      const page = await context.newPage();
      await page.goto(base + '/');
      const nav = page.locator('#site-nav');
      if (!await nav.isVisible()) await page.locator('#menu-button').click();
      check(await nav.getByRole('link', { name: 'PCR과 프라이머 디자인', exact: true }).count() === 0, 'PCR leaf excluded');
      check(await nav.locator('a').count() === 4, 'portal and three category links only');
      check(await nav.getByRole('link', { name: 'C 생물정보학 · 데이터', exact: true }).count() === 1, 'bioinformatics category preserved');
      const heading = page.locator('#bioinformatics h3').filter({ hasText: 'PCR과 프라이머 디자인' });
      check((await heading.locator('.res__code').innerText()).trim() === 'C2', 'portal code C2');
      if (engineName === 'chromium' && width === 1440) {
        await shot(page, 'portal-sidebar.png');
        await page.locator('#bioinformatics').scrollIntoViewIfNeeded(); await shot(page, 'portal-bioinformatics-section.png');
        await page.locator('#search-input').pressSequentially('PCR', { delay: 80 });
        await page.waitForSelector('#search-results a[href*="/bioinformatics/pcr-primer-design/"]');
        check(await page.locator('#search-results').innerText().then(t => t.includes('PCR과 프라이머 디자인')), 'search UI finds PCR');
        await page.locator('#search-results a[href*="/bioinformatics/pcr-primer-design/"]').first().click(); await ready(page);
        check(new URL(page.url()).pathname === path, 'search result opens PCR');
        await page.goto(base + '/');
      }
      await heading.locator(`a[href="${path}"]`).click(); await ready(page);
      check(new URL(page.url()).pathname === path, 'portal link opens PCR');
      check((await page.reload()).status() === 200, 'reload'); await ready(page);
      check(await page.locator('.pcr-back').getAttribute('href') === '/#bioinformatics', 'back link preserved');
      check(await page.locator('.pcr-eyebrow').innerText() === 'C2 · 생물정보학', 'portal code shown on worksheet');
      check(await page.title() === 'PCR과 프라이머 디자인 | HAFS Biology Lab', 'title preserved');
      check(await page.locator('link[rel=canonical]').getAttribute('href') === 'https://suimaire.github.io' + path, 'canonical preserved');
      check(await page.locator('script[src$="/page-views.js"]').count() === 1, 'one module tag');
      check(await page.locator('[data-page-views]').count() === 1, 'one counter mount');
      check(await page.locator('.site-credit__name').innerText() === 'HAFS Biology Lab', 'shared title');
      check(await page.locator('.site-credit__line').innerText() === 'Teacher-built interactive science tools · CH Park', 'shared author credit');
      check(await page.locator('.pcr-site-footer').evaluate(e => !e.closest('#pcr-worksheet') && !!(document.querySelector('#calculation-scope').compareDocumentPosition(e) & Node.DOCUMENT_POSITION_FOLLOWING)), 'footer after all worksheet content and sources');
      // Existing preview modes label read-only live data or mock data. Neither records views.
      await page.goto(base + path + '?page-views=' + previewMode); await ready(page); await counterReady(page);
      check((await page.locator('[data-page-views]').innerText()).includes(previewMode === 'live' ? '(읽기 전용)' : '(mock)'), 'preview counts visibly labelled');
      check(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth), 'no horizontal overflow');
      const identity = await page.locator('.site-credit').boundingBox(), counter = await page.locator('[data-page-views]').boundingBox();
      check(width > 600 ? counter.x > identity.x + identity.width : counter.y >= identity.y + identity.height, 'responsive footer layout');
      await page.locator('.pcr-site-footer').scrollIntoViewIfNeeded();
      if (engineName === 'chromium') {
        await shot(page, width === 1440 ? 'pcr-page-footer.png' : width === 390 ? 'pcr-page-mobile-footer.png' : 'pcr-page-tablet-footer.png');
        if (width === 1440) {
          await page.locator('#rna-extension > summary').click();
          await page.locator('[data-rna-section=reflection]').click();
          await page.setViewportSize({ width, height: 1550 });
          await page.locator('.pcr-site-footer').scrollIntoViewIfNeeded();
          await shot(page, 'pcr-page-full-bottom.png');
        }
      }
      results.push({ engine: engineName, width, navigation: 'pass', footer: 'pass', overflow: false });
      console.log(`PASS ${engineName} ${width}px navigation/footer`);
      await context.close();
    }

    // Same origin and production mode, but local built files and fake RPC responses only.
    // This exercises actual autostart/dedupe without incrementing the public counter.
    for (const scenario of ['success', 'http-error', 'network-error', 'malformed', 'timeout', 'module-unavailable', 'storage-unavailable']) {
      const context = await browser.newContext(); await instrument(context);
      const calls = [];
      await context.route('https://suimaire.github.io/**', async route => {
        const requestUrl = new URL(route.request().url());
        if (scenario === 'module-unavailable' && requestUrl.pathname.endsWith('/page-views.js')) return route.abort();
        const response = await context.request.get(base + requestUrl.pathname + requestUrl.search);
        await route.fulfill({ response });
      });
      await context.route(PAGE_VIEWS_CONFIG.supabaseUrl + '/**', async route => {
        const name = new URL(route.request().url()).pathname.split('/').pop();
        calls.push({ name, ...route.request().postDataJSON() });
        if (scenario === 'network-error') return route.abort();
        if (scenario === 'timeout') return; // unresolved until the module's own 6s AbortController
        const headers = { 'access-control-allow-origin': '*' };
        await route.fulfill({ status: scenario === 'http-error' ? 503 : 200, headers, contentType: 'application/json', body: JSON.stringify(scenario === 'malformed' ? [] : [{ today_views: 7, total_views: 31 }]) });
      });
      if (scenario === 'storage-unavailable') await context.addInitScript(() => {
        for (const name of ['localStorage', 'sessionStorage']) Object.defineProperty(window, name, { get() { throw new Error('storage denied'); } });
      });
      const page = await context.newPage();
      await page.goto('https://suimaire.github.io' + path); await ready(page);
      if (['success', 'storage-unavailable'].includes(scenario)) {
        await counterReady(page);
        check(calls.length === 1 && calls[0].name === 'record_page_view', 'exactly one initial record even without storage');
        check(calls[0].p_page_key === path, 'PCR-specific normalized key');
        if (scenario === 'success') {
          await page.reload(); await ready(page); await counterReady(page);
          check(calls.length === 2 && calls[1].name === 'get_page_view_counts', 'reload within 30 minutes only reads');
          await page.goto('https://suimaire.github.io' + path + 'index.html?check=1#activity-07'); await counterReady(page);
          check(calls.length === 3 && calls[2].name === 'get_page_view_counts' && calls[2].p_page_key === path, 'index/query/hash share the same page key');
          await page.goto('https://suimaire.github.io/'); await counterReady(page);
          check(calls.length === 4 && calls[3].name === 'record_page_view' && calls[3].p_page_key === '/', 'portal counted separately');
        } else {
          await page.locator('#first-placement').fill('storage unavailable');
          check((await page.locator('#save-status').innerText()).includes('읽을 수 없습니다'), 'worksheet handles unavailable storage');
        }
      } else {
        if (scenario !== 'module-unavailable') await page.waitForFunction(() => document.querySelector('[data-page-views]')?.hidden && document.querySelector('[data-page-views]').dataset.pageViewsMounted === '1', null, { timeout: 12000 });
        check(await page.locator('[data-page-views]').isHidden(), 'existing quiet failure fallback');
        await smoke(page, `${engineName} ${scenario}`);
        check(calls.length === (scenario === 'module-unavailable' ? 0 : 1), 'no retry or duplicate request on failure');
      }
      results.push({ engine: engineName, scenario, calls, pass: true });
      console.log(`PASS ${engineName} counter ${scenario}`); await context.close();
    }
  } finally { await browser.close(); }
}
check(pageErrors.length === 0, 'no pageerror');
check(rejections.length === 0, 'no unhandled rejection');
await writeFile(resolve(output, 'site-integration-results.json'), JSON.stringify({ checks, results, screenshots, pageErrors, rejections, consoleErrors: [...new Set(consoleErrors)], warnings: [...new Set(warnings)], screenshotCounterMode: previewMode, liveCountIncremented: false }, null, 2));
console.log(`PASS ${checks} site integration assertions; ${results.length} cases; ${screenshots.length} PNGs; pageerror/unhandledrejection 0`);
