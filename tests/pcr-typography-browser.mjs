// Activity 00 v2 copy, Korean word boundaries and technical-content overflow.
import assert from 'node:assert/strict';
import { mkdir, readFile, writeFile } from 'node:fs/promises';
import { resolve } from 'node:path';
import { chromium, webkit } from './.pcr-tools/node_modules/playwright/index.mjs';
import { startPreview } from './pcr-preview.mjs';
import { parseRecord, STORAGE_KEY } from '../assets/js/pcr-records.mjs';

const server = process.env.PCR_TEST_URL ? null : await startPreview(0);
const url = process.env.PCR_TEST_URL || `http://127.0.0.1:${server.address().port}/bioinformatics/pcr-primer-design/`;
const output = resolve('verification.local/pcr-primer-design/activity00-v2');
await mkdir(output, { recursive: true });
const seed = parseRecord(await readFile(resolve('tests/fixtures/pcr/phase6.json'), 'utf8'));
seed.answers['first-placement'] = '';
seed.finalReview.notebookExpanded = true;
const results = [], screenshots = [];
let checks = 0;
const check = (value, label) => { assert.ok(value, label); checks++; };
try {
  for (const [engineName, engine] of Object.entries({ chromium, webkit })) {
    const browser = await engine.launch({ headless: true });
    try {
      for (const width of [1440, 768, 390]) {
        const baseline = checks;
        const context = await browser.newContext({ viewport: { width, height: 1000 } });
        const page = await context.newPage(), errors = [];
        page.on('pageerror', error => errors.push(error.message));
        await page.goto(url);
        await page.evaluate(({ key, seed }) => localStorage.setItem(key, JSON.stringify(seed)), { key: STORAGE_KEY, seed });
        await page.reload(); await page.waitForSelector('#pcr-worksheet[data-ready=true]');
        await page.evaluate(() => document.fonts.ready);
        const support = await page.evaluate(() => ({ pretty: CSS.supports('text-wrap', 'pretty'), balance: CSS.supports('text-wrap', 'balance') }));
        const css = selector => page.locator(selector).first().evaluate(el => {
          const s = getComputedStyle(el);
          return { wordBreak: s.wordBreak, lineBreak: s.lineBreak, textWrap: s.getPropertyValue('text-wrap'), overflowWrap: s.overflowWrap, whiteSpace: s.whiteSpace, overflowX: s.overflowX };
        });
        const overflow = async label => check(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth), label);
        const bounded = async selector => check(await page.locator(selector).first().evaluate(el => {
          const b = el.getBoundingClientRect(); return b.left >= -1 && b.right <= innerWidth + 1;
        }), `${selector} stays inside the viewport`);
        const open = async selector => { if (!await page.locator(selector).evaluate(el => el.open)) await page.locator(`${selector} > summary`).click(); };
        const shot = async (name, selector, last, maxHeight) => {
          if (engineName !== 'chromium') return;
          const path = resolve(output, name);
          if (!last && !maxHeight) await page.locator(selector).screenshot({ path, animations: 'disabled' });
          else {
            await page.locator(selector).scrollIntoViewIfNeeded();
            const first = await page.locator(selector).boundingBox(), end = last ? await page.locator(last).boundingBox() : first;
            const offset = await page.evaluate(() => ({ x: scrollX, y: scrollY }));
            await page.screenshot({ path, fullPage: true, animations: 'disabled', clip: { x: first.x + offset.x, y: first.y + offset.y, width: first.width, height: Math.min(end.y + end.height - first.y, maxHeight || Infinity) } });
          }
          screenshots.push({ name, path, width });
        };
        const wordBoundaries = () => page.locator('#activity-00 h2 + p').evaluate(el => {
          const walker = document.createTreeWalker(el, NodeFilter.SHOW_TEXT), broken = [];
          let words = 0;
          while (walker.nextNode()) {
            const node = walker.currentNode;
            for (const match of node.textContent.matchAll(/[A-Za-z가-힣]+/g)) {
              if (!/[가-힣]/.test(match[0])) continue;
              words++;
              const range = document.createRange(); range.setStart(node, match.index); range.setEnd(node, match.index + match[0].length);
              const lines = new Set([...range.getClientRects()].map(rect => Math.round(rect.top)));
              if (lines.size > 1) broken.push(match[0]);
            }
          }
          return { words, broken };
        });
        for (const selector of ['#activity-00 h2 + p', '#activity-00 .pcr-prompt', '#activity-00 label', '#activity-00 figcaption', '#activity-00 .pcr-transition p', '.pcr-direction-principles li']) {
          const s = await css(selector);
          check(s.wordBreak === 'keep-all' && s.lineBreak === 'strict', `${selector} Korean word boundaries`);
          if (support.pretty) check(s.textWrap === 'pretty', `${selector} progressive pretty wrapping`);
        }
        const heading = await css('#activity-00 h2');
        check(heading.wordBreak === 'keep-all', 'heading keeps Korean words intact');
        if (support.balance) check(heading.textWrap === 'balance', 'heading balances lines');
        const words = await wordBoundaries();
        check(words.words > 8 && words.broken.length === 0, `no syllable split in opening paragraph: ${words.broken}`);
        check((await page.locator('#prediction-b-caption').innerText()).includes('결실 바깥의 보존된 결합 부위'), 'B binding explanation is concrete');
        check((await page.locator('#prediction-ab-concept').innerText()).includes('B에는 이 80 bp가 없으므로'), 'deletion is connected to product length difference');
        check(!(await page.locator('#activity-00').innerText()).match(/(?:260|180)\s*bp/), 'no absolute amplicon sizes in Activity 00');
        check(!(await page.locator('.pcr-map-legend').innerText()).includes('/'), 'legend uses sentences');
        await overflow('Activity 00 has no page overflow');
        const size = width === 1440 ? 'desktop' : width === 390 ? 'mobile' : 'tablet';
        await shot(`00-v2-${size}.png`, '#activity-00');
        await shot(`typography-${size}.png`, '#activity-00', '#prediction-instructions');
        if (width === 1440) await shot('00-v2-ab-concept.png', '.pcr-sample-map', '#prediction-ab-concept');

        // Existing interactive views, including intrinsically wide scientific content.
        const sequence = await css('.pcr-ordered-sequence');
        check(sequence.wordBreak === 'normal' && sequence.overflowWrap === 'anywhere' && sequence.textWrap !== 'pretty', 'sequence keeps independent wrapping');
        await open('#manual-sequences');
        await page.locator('#primer-f').fill('ACGT'.repeat(25));
        check(await page.locator('#primer-f').inputValue() === 'ACGT'.repeat(25), 'long primer value retained');
        await bounded('#primer-f'); await overflow('long workbench sequence is contained');
        await page.getByRole('button', { name: '설계 1 불러오기', exact: true }).click();
        if (width !== 768) await shot(`regression-03-${size}.png`, '#activity-03', null, 1500);
        await page.locator('#review-tab-complementarity').click();
        const alignment = await css('.pcr-review-alignment');
        check(alignment.wordBreak === 'normal' && alignment.whiteSpace === 'pre' && alignment.overflowWrap === 'normal', 'alignment columns never reflow');
        check((await css('.pcr-alignment-scroll')).overflowX === 'auto', 'alignment has its own horizontal scroll');
        await overflow('alignment does not overflow the page');
        if (width !== 768) await shot(`regression-04-${size}.png`, '#review-panel-complementarity');
        await bounded('.pcr-gel'); await overflow('gel fits the page');
        if (width !== 768) await shot(`regression-05-${size}.png`, '#evidence-gel-figure');
        await page.locator('#ext-route-paper').click();
        const longUrl = 'https://example.invalid/' + 'LongExternalIdentifier'.repeat(20);
        await page.locator('#ext-paper-source').fill(longUrl);
        await page.locator('#ext-paper-forward').fill('ACGT'.repeat(25));
        check(await page.locator('#ext-paper-source').inputValue() === longUrl, 'long URL remains intact in field');
        await bounded('#ext-paper-source'); await bounded('#ext-paper-forward'); await overflow('external URL/primer fields fit the page');
        await page.locator('#ext-route-new').click();
        const accession = 'NM_01234567890123456789.123456789';
        await page.locator('#ext-design-target').fill(accession);
        check((await css('#ext-design-target')).wordBreak === 'normal', 'accession field does not inherit Korean wrapping');
        await open('#ext-plan');
        await page.locator('#ext-plan-target').fill(`Gene ID ${'1234567890'.repeat(20)}`);
        await bounded('#ext-design-target'); await bounded('#ext-plan-target'); await overflow('accession/Gene ID fields fit the page');
        if (width !== 768) await shot(`regression-06-${size}.png`, '#ext-route-input-new');
        await page.locator('#final-unknown').fill(`아직 확인하지 않은 내용입니다. ${longUrl} ${'TechnicalIdentifier'.repeat(20)}`);
        if (await page.locator('#final-notebook').isHidden()) await page.locator('#final-notebook-toggle').click();
        await bounded('#final-notebook'); await overflow('notebook wraps long English tokens and URL');
        check((await css('#final-notebook p')).wordBreak === 'keep-all', 'notebook prose uses Korean wrapping');
        if (width !== 768) await shot(`regression-07-${size}.png`, '#final-notebook', null, 1600);
        await open('#rna-extension'); await page.locator('[data-rna-section=transcripts]').click();
        check((await css('#rna-transcripts p')).wordBreak === 'keep-all', 'RNA prose uses Korean wrapping');
        await bounded('#rna-transcript-comparison'); await overflow('RNA transcript diagrams fit the page');
        if (width !== 768) await shot(`regression-rna-${size}.png`, '#rna-transcripts', null, 1500);

        // Force the progressive-enhancement blocks off, approximating older CSS support.
        await page.evaluate(() => {
          for (const sheet of document.styleSheets) {
            for (let i = sheet.cssRules.length - 1; i >= 0; i--) {
              const rule = sheet.cssRules[i];
              if (rule.conditionText?.includes('text-wrap:')) sheet.deleteRule(i);
            }
          }
        });
        check((await css('#activity-00 h2 + p')).wordBreak === 'keep-all', 'fallback retains keep-all');
        check((await css('#activity-00 h2 + p')).textWrap !== 'pretty', 'pretty enhancement disabled for fallback check');
        check((await css('#activity-00 h2')).textWrap !== 'balance', 'balance enhancement disabled for fallback check');
        check((await wordBoundaries()).broken.length === 0, 'fallback still keeps Korean words intact');
        await overflow('fallback has no page overflow');
        check(errors.length === 0, `runtime errors: ${errors.join(', ')}`);
        results.push({ engine: engineName, browserVersion: browser.version(), width, support, words, checks: checks - baseline, passed: true });
        console.log(`PASS typography ${engineName} ${width}px (${checks - baseline} assertions)`);
        await context.close();
      }
    } finally { await browser.close(); }
  }
  await writeFile(resolve(output, 'typography-verification.json'), JSON.stringify({ url, checks, results, screenshots }, null, 2));
  console.log(`PASS ${checks} typography assertions; ${screenshots.length} PNG files`);
} finally { if (server) await new Promise(resolve => server.close(resolve)); }
