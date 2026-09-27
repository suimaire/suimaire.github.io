// Focused Activity 00 behavior, compatibility and automatic PNG evidence.
// PLAYWRIGHT_BROWSERS_PATH=tests/.pcr-tools/browsers node tests/pcr-prediction-browser.mjs
import assert from 'node:assert/strict';
import { mkdir, readFile, writeFile } from 'node:fs/promises';
import { resolve } from 'node:path';
import { chromium, webkit } from './.pcr-tools/node_modules/playwright/index.mjs';
import { startPreview } from './pcr-preview.mjs';
import { STORAGE_KEY, parseRecord } from '../assets/js/pcr-records.mjs';

const server = process.env.PCR_TEST_URL ? null : await startPreview(0);
const url = process.env.PCR_TEST_URL || `http://127.0.0.1:${server.address().port}/bioinformatics/pcr-primer-design/`;
const output = resolve(process.env.PCR_PREDICTION_OUTPUT || 'verification.local/pcr-primer-design/activity00');
await mkdir(output, { recursive: true });
const screenshots = [], results = [];
let checks = 0;
const check = (value, label) => { assert.ok(value, label); checks++; };
const near = (a, b, label) => check(Math.abs(a - b) < 1, `${label}: ${a} vs ${b}`);
const stable = record => { const copy = structuredClone(record); delete copy.updatedAt; return copy; };
try {
  for (const [engineName, engine] of Object.entries({ chromium, webkit })) {
    const browser = await engine.launch({ headless: true });
    try {
      for (const width of [1440, 768, 390]) {
        const baseline = checks;
        const context = await browser.newContext({ viewport: { width, height: 1000 }, hasTouch: width === 390 });
        const page = await context.newPage(), errors = [];
        page.on('pageerror', e => errors.push(e.message));
        const ready = () => page.waitForSelector('#pcr-worksheet[data-ready=true]');
        const saved = () => page.evaluate(key => JSON.parse(localStorage.getItem(key)), STORAGE_KEY);
        const range = async (name, value) => {
          const slider = page.locator(`#prediction-${name}`);
          await slider.focus(); await page.keyboard.press('Home');
          for (let n = 5; n < value; n += 5) await page.keyboard.press('ArrowRight');
        };
        const pair = async (f, r) => { await range('forward', f); await range('reverse', r); };
        const text = () => page.locator('#prediction-feedback').innerText();
        const shot = async (name, selector) => {
          if (engineName !== 'chromium') return;
          const path = resolve(output, name);
          await page.locator(selector).screenshot({ path, animations: 'disabled' });
          screenshots.push({ name, path, width });
        };
        const geometry = async expectedDifference => {
          const a = await page.locator('#prediction-track').boundingBox();
          const b = await page.locator('#prediction-b-track').boundingBox();
          const af = await page.locator('#prediction-forward-binding').boundingBox();
          const bf = await page.locator('#prediction-b-forward-binding').boundingBox();
          const ar = await page.locator('#prediction-reverse-binding').boundingBox();
          const br = await page.locator('#prediction-b-reverse-binding').boundingBox();
          near(b.width / a.width * 420, 340, 'B shares physical bp scale');
          near(af.width, bf.width, 'same F footprint in A/B');
          near(ar.width, br.width, 'same R footprint in A/B');
          const spanA = await page.locator('#prediction-a-span').boundingBox();
          const spanB = await page.locator('#prediction-b-span').boundingBox();
          near((spanA.width - spanB.width) / a.width * 420, expectedDifference, 'deletion changes the interval length');
        };
        await page.goto(url); await ready(); await page.evaluate(() => document.fonts.ready);
        for (const id of ['forward', 'reverse', 'b-forward', 'b-reverse']) check(await page.locator(`#prediction-${id}-marker`).isHidden(), 'no inferred placement before input');
        if (width === 1440) await shot('00-before-placement.png', '#activity-00');
        check((await page.locator('#prediction-c-row').innerText()).includes('결합 여부 판단 보류'), 'C undecided before input');
        const track = page.locator('#prediction-track');
        const box = await track.boundingBox();
        await track[width === 390 ? 'tap' : 'click']({ position: { x: box.width * .15, y: 20 } });
        check(await page.locator('#prediction-forward').inputValue() === '15', 'mouse/touch selects A reference position');
        check(await page.locator('#prediction-b-forward-marker').isVisible(), 'single F appears on B immediately');
        check(await page.locator('#prediction-b-reverse-marker').isHidden(), 'R remains unselected');
        await page.locator('[data-prediction-primer=reverse]').click();
        await track.click({ position: { x: box.width * .75, y: 20 } });
        await page.locator('#prediction-reverse').focus(); await page.keyboard.press('ArrowRight');
        check((await text()).includes('개념 조건 충족'), 'first valid shared pair');
        await geometry(80);
        const state = await saved();
        await page.locator('#prediction-b-track').click();
        assert.deepEqual(await saved(), state); checks++;
        if (width === 1440) {
          await shot('00-desktop-full.png', '#activity-00');
          await shot('00-shared-ab-primers.png', '.pcr-sample-map');
          await shot('00-c-undecided.png', '#prediction-c-row');
        }
        if (width === 390) await shot('00-mobile-full.png', '#activity-00');
        await pair(20, 75);
        check((await text()).includes('개념 조건 충족'), 'another valid pair accepted');
        check(await page.locator('#prediction-forward-marker').evaluate(el => el.style.left) === '20%', 'original A percent preserved');
        await geometry(80);
        await pair(35, 80);
        check((await text()).includes('겹침 주의'), 'deletion interior warning');
        check(await page.locator('#prediction-b-forward-marker').isHidden() && await page.locator('#prediction-b-reverse-marker').isVisible(), 'no fictitious B F binding');
        check((await page.locator('#prediction-b-caption').innerText()).includes('F: 결실과 겹쳐'), 'missing B marker explicitly explained');
        check(await page.locator('#prediction-b-span').isHidden(), 'no B interval when one binding is absent');
        if (width === 1440) await shot('00-overlap-warning.png', '#activity-00');
        if (width === 390) await shot('00-mobile-warning.png', '#activity-00');
        await pair(15, 50);
        check((await text()).includes('겹침 주의'), 'outside center with footprint crossing the edge warns');
        check(await page.locator('#prediction-b-reverse-marker').isHidden(), 'R overlap handled symmetrically');
        await pair(5, 20);
        check((await text()).includes('결실 포함 확인'), 'same-side warning');
        await geometry(0);
        if (width === 1440) await shot('00-outside-warning.png', '#activity-00');
        await pair(80, 15);
        check((await text()).includes('배치 확인'), 'reversed pair warning');
        check(await page.locator('#prediction-a-span').isHidden() && await page.locator('#prediction-b-span').isHidden(), 'no span for outward pair');
        await pair(15, 15);
        check((await text()).includes('배치 확인'), 'coincident pair warning');
        await pair(25, 55);
        check((await text()).includes('개념 조건 충족'), 'valid locations close to both deletion boundaries');
        await page.locator('#first-placement').fill('같은 pair의 두 결합 부위와 결실 위치를 비교했다. C는 뒤 활동에서 확인한다.');
        const beforeReload = await saved();
        await page.reload(); await ready();
        assert.deepEqual(await saved(), beforeReload); checks++;
        check(await page.locator('#prediction-b-forward-marker').isVisible() && (await text()).includes('개념 조건 충족'), 'derived B view restored without saved B fields');
        await page.locator('#record-menu > summary').click();
        const download = page.waitForEvent('download'); await page.locator('#export-record').click();
        const exported = await download, file = resolve(output, `${engineName}-${width}-record.json`); await exported.saveAs(file);
        const exportedRecord = JSON.parse(await readFile(file, 'utf8'));
        assert.deepEqual(stable(exportedRecord), stable(beforeReload)); checks++;
        await pair(35, 80);
        page.once('dialog', d => d.accept()); await page.locator('#import-record').setInputFiles(file);
        await page.waitForFunction(() => document.querySelector('#prediction-forward').value === '25');
        assert.deepEqual(stable(await saved()), stable(beforeReload)); checks++;
        check((await text()).includes('개념 조건 충족'), 'JSON import redraws the shared pair');
        await page.locator('#activity-00 .pcr-explanation > summary').click();
        check((await page.locator('#activity-00 .pcr-explanation').innerText()).includes('같은 연속 결합 부위'), 'explanation includes binding failure');
        await page.addScriptTag({ path: resolve('tests/.pcr-tools/node_modules/axe-core/axe.min.js') });
        const audit = await page.evaluate(async () => (await axe.run('#activity-00', { runOnly: { type: 'tag', values: ['wcag2a', 'wcag2aa', 'wcag21aa', 'wcag22aa', 'best-practice'] } })).violations.map(v => ({ id: v.id, nodes: v.nodes.map(n => n.target) })));
        assert.deepEqual(audit, []); checks++;
        check(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth), 'no horizontal overflow');
        if (engineName === 'chromium' && width === 1440) {
          for (let phase = 1; phase <= 6; phase++) {
            const historicalFile = resolve(`tests/fixtures/pcr/phase${phase}.json`);
            const raw = await readFile(historicalFile, 'utf8'), expected = parseRecord(raw);
            await page.evaluate(({ key, raw }) => localStorage.setItem(key, raw), { key: STORAGE_KEY, raw });
            await page.reload(); await ready();
            assert.deepEqual(parseRecord(JSON.stringify(await saved())), expected); checks++;
            await page.locator('#record-menu > summary').click();
            page.once('dialog', d => d.accept()); await page.locator('#import-record').setInputFiles(historicalFile);
            await page.waitForFunction(() => document.querySelector('#import-record').value === '' && document.querySelector('#save-status').textContent.includes('자동 저장'));
            assert.deepEqual(stable(await saved()), stable(expected)); checks++;
            const before = stable(await saved()); delete before.initialPrimerPrediction;
            await pair(20, 75);
            const after = stable(await saved()); delete after.initialPrimerPrediction;
            assert.deepEqual(after, before); checks++;
          }
        }
        check(errors.length === 0, `no runtime errors: ${errors.join(', ')}`);
        results.push({ engine: engineName, width, checks: checks - baseline, accessibilityViolations: audit.length, passed: true });
        console.log(`PASS Activity 00 ${engineName} ${width}px (${checks - baseline} assertions)`);
        await context.close();
      }
    } finally { await browser.close(); }
  }
  await writeFile(resolve(output, 'verification.json'), JSON.stringify({ url, generatedAt: new Date().toISOString(), checks, results, screenshots }, null, 2));
  console.log(`PASS ${checks} Activity 00 assertions; ${screenshots.length} PNG files in ${output}`);
} finally { if (server) await new Promise(resolve => server.close(resolve)); }
