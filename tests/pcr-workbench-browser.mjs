// Real interactions plus PNG evidence. Uses the same project-local Playwright as the regression suite.
import assert from 'node:assert/strict';
import { mkdir, readFile, writeFile } from 'node:fs/promises';
import { resolve } from 'node:path';
import { chromium, webkit } from './.pcr-tools/node_modules/playwright/index.mjs';
import { STORAGE_KEY } from '../assets/js/pcr-records.mjs';
const url = process.env.PCR_TEST_URL || 'http://127.0.0.1:4173/bioinformatics/pcr-primer-design/';
const output = resolve(process.env.PCR_VERIFICATION_ROOT || 'verification.local/pcr-primer-design', 'phase2');
await mkdir(output, { recursive: true });
const fixture = JSON.parse(await readFile(new URL('../assets/data/pcr-primer-fixture.json', import.meta.url), 'utf8'));
const evidence = [], results = [];
let checks = 0;
const check = (condition, message) => { assert.ok(condition, message); checks++; };
for (const [engineName, engine] of Object.entries({ chromium, webkit })) {
  const browser = await engine.launch({ headless: true });
  for (const width of [1440, 768, 390]) {
    const context = await browser.newContext({ viewport: { width, height: 1150 }, hasTouch: width === 390 });
    const page = await context.newPage(), errors = [], baseline = checks;
    page.on('pageerror', e => errors.push(e.message));
    page.on('requestfailed', request => errors.push(`${request.url()}: ${request.failure()?.errorText}`));
    await page.goto(url);
    try { await page.waitForSelector('#pcr-worksheet[data-ready=true]'); }
    catch (error) { throw new Error(`${error.message}\n${await page.locator('#save-status').textContent()}\n${await page.locator('#analysis-status').textContent()}\n${errors.join('\n')}`); }
    const screenshot = async (name, selector = '#activity-03') => {
      if (engineName !== 'chromium' || (width !== 1440 && !['03-mobile.png', '03-tablet.png'].includes(name))) return;
      await page.locator(selector).screenshot({ path: resolve(output, name), animations: 'disabled' });
      evidence.push({ name, path: resolve(output, name), viewport: { width, height: 1150 } });
    };
    const open = async id => { if (!await page.locator(id).evaluate(el => el.open)) await page.locator(`${id} > summary`).click(); };
    const close = async id => { if (await page.locator(id).evaluate(el => el.open)) await page.locator(`${id} > summary`).click(); };
    const fine = async (name, start, end, direction) => {
      await page.locator(`[data-select-primer=${name}]`).click(); await open('#coordinate-details');
      await page.locator('#range-start').fill(String(start)); await page.locator('#range-end').fill(String(end));
      if (direction) await page.locator('#range-direction').selectOption(direction);
      await page.locator('#apply-range').click(); await close('#coordinate-details');
    };
    const saved = async () => JSON.parse(await page.evaluate(key => localStorage.getItem(key), STORAGE_KEY));
    const resultText = () => page.locator('#analysis-results').innerText();
    const primerText = name => page.locator(`#primer-summary-${name.toLowerCase()}`).innerText();

    // 00: real range-key interactions and exact marker synchronization.
    await page.locator('#prediction-forward').focus(); await page.keyboard.press('Home');
    await page.keyboard.press('ArrowRight'); await page.keyboard.press('ArrowRight');
    await page.locator('#prediction-reverse').focus(); await page.keyboard.press('End');
    await page.keyboard.press('ArrowLeft'); await page.keyboard.press('ArrowLeft'); await page.keyboard.press('ArrowLeft');
    for (const [name, value] of [['forward', 15], ['reverse', 80]]) {
      check(await page.locator(`#prediction-${name}`).inputValue() === String(value), `00 ${name} slider`);
      check(await page.locator(`#prediction-${name}-marker`).evaluate(el => el.style.left) === `${value}%`, `00 ${name} marker`);
      check(await page.locator(`#prediction-${name}-marker`).isVisible(), `00 ${name} visible`);
    }
    await screenshot('00-overview.png', '#activity-00');
    for (const cycle of [1, 2, 3]) {
      await page.locator(`[data-comparison-cycle="${cycle}"]`).click();
      const description = await page.locator('#cycle-products-description').innerText();
      check(description.startsWith(`${cycle}주기`), `cycle ${cycle} explanation matches button`);
      check(await page.locator(`[data-comparison-cycle="${cycle}"]`).getAttribute('aria-pressed') === 'true', `cycle ${cycle} active`);
      check((await page.locator('#cycle-products-svg').getAttribute('aria-label')) === description, `cycle ${cycle} accessible diagram`);
      if (cycle !== 2) await screenshot(`01-cycle-${cycle}.png`, '.pcr-cycle-comparison');
    }
    check((await page.locator('.pcr-direction-principles').innerText()).includes('두 결합 부위 사이의 표적 영역'), 'correct 02 synthesis explanation');
    await page.locator('#direction-complement').fill('TCAGGCAT'); await page.locator('#check-complement').click();
    await page.locator('#reverse-complement').click(); await page.locator('#direction-reverse').fill('TACGGACT'); await page.locator('#check-direction').click();
    check((await page.locator('#direction-feedback').innerText()).includes('일치합니다'), 'complement exercise intact');
    await screenshot('02-direction.png', '#activity-02');

    // Empty -> F pointer/touch endpoints -> R keyboard. Nothing is seeded into page state.
    await page.locator('#activity-03').scrollIntoViewIfNeeded();
    check((await resultText()).includes('기다리고'), 'empty calculation');
    check(await page.locator('.pcr-base').count() === (width === 390 ? 30 : 60), 'bounded zoom window');
    await screenshot('03-workbench-empty.png');
    const zoom = page.locator('#sequence-window-slider'); await zoom.focus(); await page.keyboard.press('Home');
    for (let i = 0; i < 40; i++) await page.keyboard.press('ArrowRight');
    for (const position of [41, 60]) {
      const cell = page.locator(`[data-base="${position}"]`);
      if (width === 390) await cell.tap(); else await cell.click();
    }
    check((await saved()).draft.forward === fixture.candidatePairs.P1.forward, 'F two endpoint selection');
    check((await primerText('F')).includes('20 nt / GC 50.0%'), 'F live stats');
    check((await saved()).draft.reverse === '', 'F selection does not invent R');
    await screenshot('03-forward-selected.png');
    await page.locator('[data-select-primer=R]').click();
    check(await page.locator('[data-select-primer=R]').getAttribute('aria-pressed') === 'true', 'R explicit mode state');
    await zoom.focus(); await page.keyboard.press('End');
    for (let i = 0; i < (width === 390 ? 110 : 80); i++) await page.keyboard.press('ArrowLeft');
    const rStart = page.locator('[data-base="281"]'); await rStart.focus();
    check(await rStart.evaluate(el => getComputedStyle(el).outlineStyle !== 'none'), 'visible keyboard focus');
    await page.keyboard.press('Enter'); for (let i = 0; i < 19; i++) await page.keyboard.press('ArrowRight'); await page.keyboard.press('Enter');
    check((await saved()).draft.reverse === fixture.candidatePairs.P1.reverse, 'R keyboard selection and reverse complement');
    check((await primerText('R')).includes('TCGTACGATCGTAGCCTGAA') && (await primerText('R')).includes('reverse complement'), 'R shows reference and order transformation');
    check((await resultText()).includes('260 bp') && (await resultText()).includes('180 bp') && (await resultText()).includes('예상 산물 없음'), 'live ABC fixture results');
    check((await page.locator('#map-description').textContent()).includes('41–300, 260 bp'), 'map accessible product range');
    check((await page.locator('#deletion-message').innerText()).includes('포함되어'), 'deletion between primers');
    check(await page.locator('#activity-03 #gel-results').count() === 0, 'no gel in workbench');
    check(await page.locator('#activity-03 .pcr-map-primer').count() > 2, 'dynamic structural comparison');
    await page.locator('#show-initial-prediction').click();
    check((await page.locator('#map-prediction').textContent()).includes('00 예상 F'), '00 prediction overlay');
    check((await saved()).draft.forward === fixture.candidatePairs.P1.forward, 'overlay does not change current design');
    await screenshot(width === 390 ? '03-mobile.png' : width === 768 ? '03-tablet.png' : '03-complete-design.png');
    check(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth), 'no viewport overflow');
    const left = await page.locator('.pcr-workspace').boundingBox(), right = await page.locator('#current-design').boundingBox();
    check(width <= 1150 ? right.y > left.y + left.height - 1 : right.x > left.x + left.width, 'responsive workbench layout');
    check(await page.locator('.pcr-base').evaluateAll(elements => elements.every(el => el.getBoundingClientRect().width >= 43.5 && el.getBoundingClientRect().height >= 44)), 'touch sized sequence targets');
    await page.locator('#design-reason').fill('결실 양옆에서 시작한 첫 설계'); await page.locator('#save-design').click();
    const originalDesign = (await saved()).designs[0];

    // Drag actual bases. Include all engines; touch users have the independently tested two-tap path.
    await page.locator('[data-select-primer=F]').click();
    const first = page.locator('[data-base="41"]'), last = page.locator('[data-base="60"]');
    await first.scrollIntoViewIfNeeded(); const from = await first.boundingBox(), to = await last.boundingBox();
    await page.mouse.move(from.x + from.width / 2, from.y + from.height / 2); await page.mouse.down();
    await page.mouse.move(to.x + to.width / 2, to.y + to.height / 2, { steps: 15 }); await page.mouse.up();
    check((await saved()).draft.bindings.F.start === '41' && (await saved()).draft.bindings.F.end === '60', 'drag commits actual range');
    // Keyboard across zoom edge and Shift extension.
    await page.locator('[data-base="41"]').focus(); await page.keyboard.press('End');
    const edge = await page.evaluate(() => Number(document.activeElement.dataset.base));
    await page.keyboard.press('Enter'); await page.keyboard.press('ArrowRight'); await page.keyboard.press('Enter');
    check((await saved()).draft.bindings.F.end === String(edge + 1), 'keyboard selects across zoom boundary');
    await page.keyboard.press('Shift+ArrowRight');
    check((await saved()).draft.bindings.F.end === String(edge + 2), 'Shift keyboard extends range');
    await page.keyboard.press('Escape');

    await fine('F', 111, 130, 'right');
    check((await page.locator('#binding-message').innerText()).includes('Forward') && (await page.locator('#binding-message').innerText()).includes('존재하지 않습니다'), 'B binding deletion explained');
    check((await saved()).draft.bindings.F.start === '111', 'overlap retained');
    check((await resultText()).includes('190 bp') && !(await resultText()).includes('110 bp'), 'overlap is not A minus 80');
    check((await page.locator('.pcr-product-values dd').nth(1).innerText()) === '예상 산물 없음', 'B has no overlap product');
    await screenshot('03-deletion-overlap.png');
    await page.locator('#design-reason').fill('Forward를 결실 경계 위로 옮김'); await page.locator('#save-design').click();
    check(JSON.stringify((await saved()).designs[0]) === JSON.stringify(originalDesign), 'saved design immutable after edits');
    await screenshot('03-saved-designs.png', '.pcr-design-records');

    await fine('F', 321, 340, 'right'); await fine('R', 41, 60, 'left');
    check((await page.locator('#placement-message').innerText()).includes('inward-facing'), 'reversed F/R order explained');
    check((await saved()).draft.bindings.F.start === '321', 'invalid layout is not silently fixed');
    await screenshot('03-invalid-design.png');
    await fine('F', 41, 60, 'right'); await fine('R', 281, 300, 'right');
    check((await page.locator('#placement-message').innerText()).includes('inward-facing'), 'parallel orientation explained');
    await fine('R', 51, 70, 'left');
    check((await page.locator('#placement-message').innerText()).includes('겹칩니다'), 'overlapping primer regions explained');
    await fine('F', 0, 430, 'right');
    check((await page.locator('#analysis-status').innerText()).includes('좌표'), 'out of range explanation');
    check((await saved()).draft.bindings.F.start === '0' && (await saved()).draft.bindings.F.end === '430', 'out of range values retained');
    check(!(await resultText()).includes('bp'), 'invalid edit clears previous product');
    await page.reload(); await page.waitForSelector('#pcr-worksheet[data-ready=true]');
    check((await saved()).draft.bindings.F.start === '0', 'invalid draft restored without correction');
    check(await page.locator('#prediction-forward').inputValue() === '15', 'initial prediction survives reload');
    check((await saved()).workbench.showPrediction, 'prediction overlay preference restored');
    await fine('F', 1, 101, 'right'); check((await primerText('F')).includes('101 nt'), 'long primer preserved with limit');
    await fine('F', 41, 45, 'right'); check((await primerText('F')).includes('5 nt'), 'short primer allowed');
    await fine('F', 231, 250, 'right'); await fine('R', 281, 300, 'left');
    check((await page.locator('#deletion-message').innerText()).includes('포함되지'), 'same side of deletion explained');
    check((await page.locator('.pcr-product-values dd').nth(0).innerText()) === (await page.locator('.pcr-product-values dd').nth(1).innerText()), 'same side A/B equal length');
    await page.locator('#design-reason').fill('두 primer를 결실 오른쪽으로 옮김'); await page.locator('#save-design').click();
    check((await saved()).designs.length === 3 && await page.locator('#save-design').isDisabled(), 'three design limit');
    await page.getByRole('button', { name: '설계 1 불러오기', exact: true }).click();
    check((await resultText()).includes('260 bp'), 'saved design reload recalculates');
    check((await saved()).draft.bindings.F.start === '41', 'saved coordinate snapshot restored');

    await page.locator('#record-menu > summary').click(); const downloadEvent = page.waitForEvent('download'); await page.locator('#export-record').click();
    const download = await downloadEvent, jsonPath = resolve(output, `${engineName}-${width}-roundtrip.json`); await download.saveAs(jsonPath);
    const exported = JSON.parse(await readFile(jsonPath, 'utf8'));
    check(exported.designs.length === 3 && exported.draft.bindings.R.direction === 'left', 'new state exported');
    await fine('F', 111, 130, 'right');
    await open('#record-menu'); page.once('dialog', dialog => dialog.accept()); await page.locator('#import-record').setInputFiles(jsonPath);
    await page.waitForFunction(() => document.querySelector('#analysis-results').textContent.includes('260 bp'));
    check(JSON.stringify((await saved()).designs) === JSON.stringify(exported.designs), 'JSON restores independent snapshots');
    check((await saved()).draft.bindings.F.start === '41', 'JSON restores current coordinates');
    await page.reload(); await page.waitForSelector('#pcr-worksheet[data-ready=true]');
    check((await resultText()).includes('180 bp') && (await saved()).designs.length === 3, 'localStorage restores live workbench');
    await page.emulateMedia({ reducedMotion: 'reduce' });
    check(await page.locator('#sequence-grid').evaluate(el => getComputedStyle(el).animationName) === 'none', 'reduced motion');
    check(errors.length === 0, `no runtime errors: ${errors.join(', ')}`);
    results.push({ engine: engineName, width, checks: checks - baseline, passed: true });
    console.log(`PASS Phase 2 ${engineName} ${width}px (${checks - baseline} assertions)`);
    await context.close();
  }
  await browser.close();
}
await writeFile(resolve(output, 'verification.json'), JSON.stringify({ url, generatedAt: new Date().toISOString(), checks, results, evidence }, null, 2));
console.log(`PASS Phase 2 ${checks} browser assertions. ${evidence.length} PNG files in ${output}`);
