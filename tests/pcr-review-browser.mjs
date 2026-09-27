// PLAYWRIGHT_BROWSERS_PATH=tests/.pcr-tools/browsers node tests/pcr-review-browser.mjs
import assert from 'node:assert/strict';
import { mkdir, readFile, writeFile } from 'node:fs/promises';
import { resolve } from 'node:path';
import { chromium, webkit } from './.pcr-tools/node_modules/playwright/index.mjs';
import { STORAGE_KEY, emptyRecord } from '../assets/js/pcr-records.mjs';
import { primerStats, reverseComplement } from '../assets/js/pcr-core.mjs';
import { complementarity, REVIEW_LENSES } from '../assets/js/pcr-review.mjs';
const fixture = JSON.parse(await readFile(new URL('../assets/data/pcr-primer-fixture.json', import.meta.url)));
const url = process.env.PCR_TEST_URL || 'http://127.0.0.1:4173/bioinformatics/pcr-primer-design/';
const output = resolve(process.env.PCR_VERIFICATION_ROOT || 'verification.local/pcr-primer-design', 'phase3'); await mkdir(output, { recursive: true });
let checks = 0; const results = [], evidence = [];
const check = (value, label) => { assert.ok(value, label); checks++; };
for (const [engineName, engine] of Object.entries({ chromium, webkit })) {
  const browser = await engine.launch({ headless: true });
  for (const width of [1440, 768, 390]) {
    const baseline = checks, context = await browser.newContext({ viewport: { width, height: 1150 }, hasTouch: width === 390 });
    const page = await context.newPage(), errors = [];
    page.on('pageerror', e => errors.push(e.message));
    await page.goto(url); await page.waitForSelector('#pcr-worksheet[data-ready=true]');
    const open = async id => { if (!await page.locator(id).evaluate(el => el.open)) await page.locator(`${id} > summary`).click(); };
    const saved = async () => JSON.parse(await page.evaluate(key => localStorage.getItem(key), STORAGE_KEY));
    const text = id => page.locator(id).innerText();
    const lens = name => page.locator(`#review-tab-${name}`).click();
    const screenshot = async (name, selector = '#activity-04') => {
      if (engineName !== 'chromium' || (width !== 1440 && !['04-mobile.png', '04-tablet.png'].includes(name))) return;
      await page.locator(selector).screenshot({ path: resolve(output, name), animations: 'disabled' });
      evidence.push({ name, path: resolve(output, name), width, selector });
    };
    const fine = async (name, start, end) => {
      await page.locator(`[data-select-primer=${name}]`).click(); await open('#coordinate-details');
      await page.locator('#range-start').fill(String(start)); await page.locator('#range-end').fill(String(end)); await page.locator('#apply-range').click();
    };
    check((await text('#review-status')).includes('03에서 먼저'), 'no design guidance');
    check(await page.locator('#review-my-design').textContent() === '', 'no invented student sequence');
    check(await page.locator('#review-panel-off-target').isHidden(), 'fixed candidates are not first screen');
    await fine('F', 41, 60); await fine('R', 281, 300);
    const pair = fixture.candidatePairs.P1;
    check((await text('#review-my-design')).includes(pair.forward) && (await text('#review-my-design')).includes(pair.reverse), 'real draft pair passed from 03');
    check((await text('#review-my-design')).includes('260 bp') && (await text('#review-my-design')).includes('180 bp'), 'draft products passed');
    for (const [i, seq] of [pair.forward, pair.reverse].entries()) {
      const stats = primerStats(seq);
      check((await page.locator('#review-length-values tbody tr').nth(0).locator('td').nth(i).innerText()) === `${stats.length} nt`, 'length uses sequence');
      check((await page.locator('#review-length-values tbody tr').nth(1).locator('td').nth(i).innerText()) === `${Math.round(stats.gcPercent)}%`, 'GC uses sequence');
    }
    await page.locator('#review-length-choice').selectOption('no'); await page.locator('#review-length-reason').fill('조성이 같아도 결합 위치와 상보 구간은 다를 수 있다.');
    await screenshot('04-my-design.png'); await screenshot('04-length-gc.png', '#review-panel-length-gc');
    await page.locator('#save-design').click();
    const snapshot = (await saved()).designs, originalDraft = (await saved()).draft;
    await page.locator('#review-design').selectOption('1');
    await fine('F', 121, 140);
    check((await text('#review-status')).includes('저장 설계 1'), 'explicit snapshot takes priority');
    check((await text('#review-my-design')).includes(pair.forward), 'saved sequence is independent from edited draft');
    await page.locator('#review-design').selectOption('draft');
    check((await text('#review-my-design')).includes(fixture.candidatePairs.P3.forward), 'draft selector returns current draft');
    await page.locator('#review-design').selectOption('1');
    const protectedDraft = (await saved()).draft;

    await page.locator('#review-tab-length-gc').focus();
    for (const [key, name] of [['ArrowRight', 'tm'], ['End', 'off-target'], ['ArrowRight', 'length-gc'], ['ArrowLeft', 'off-target'], ['Home', 'length-gc']]) {
      await page.keyboard.press(key);
      check(await page.locator(`#review-tab-${name}`).evaluate(el => el === document.activeElement), `${key} moves focus`);
      check(await page.locator(`#review-tab-${name}`).getAttribute('aria-selected') === 'true', 'tab aria selected');
      check(await page.locator('#review-tabs [tabindex="0"]').count() === 1, 'roving tab stop');
    }
    check(await page.locator('#review-tab-length-gc').evaluate(el => getComputedStyle(el).outlineStyle !== 'none'), 'visible keyboard focus');
    await lens('tm');
    const tmText = await text('#review-panel-tm');
    check(tmText.includes('간이 Tm / 교육용 근사') && tmText.includes('2×(A+T)+4×(G+C)'), 'honest Tm label and formula');
    check(tmText.includes('Na+') && tmText.includes('Mg2+') && tmText.includes('dNTP'), 'Tm unmodeled conditions');
    for (const sequence of [pair.forward, pair.reverse]) check((await text('#tm-values')).includes(`${primerStats(sequence).simpleTm} °C`), 'correct Wallace Tm');
    check((await text('#tm-values')).includes(`약 ${Math.abs(primerStats(pair.forward).simpleTm - primerStats(pair.reverse).simpleTm)} °C`), 'absolute delta');
    await screenshot('04-tm.png', '#review-panel-tm');
    await lens('end');
    for (const [i, seq] of [pair.forward, pair.reverse].entries()) {
      check(await page.locator('#review-end-values .pcr-review-end').nth(i).innerText() === seq.slice(-5), 'last five highlighted');
      check((await page.locator('#review-end-values .pcr-end-row').nth(i).innerText()).includes(`마지막 염기: ${seq.slice(-1)}`), 'physical last base');
    }
    await page.locator('#review-end-observation').fill('마지막 5 nt의 염기 배열을 비교했다.');
    await screenshot('04-end-structure.png', '#review-panel-end');
    await lens('complementarity');
    for (const [id, a, b] of [['self-f', pair.forward, pair.forward], ['self-r', pair.reverse, pair.reverse], ['pair', pair.forward, pair.reverse]]) {
      const model = complementarity(a, b), block = `#review-${id}`;
      check((await text(block)).includes(`가장 긴 연속 상보 구간: ${model.longest.length} nt`), `${id} longest match`);
      check((await text(block)).includes(`위: ${id === 'self-r' ? 'R' : 'F'} 5′→3′`), 'top strand label');
      check((await page.locator(`${block} pre`).first().textContent()).includes([...b].reverse().join('').slice(-5)), 'bottom reversed physical sequence');
      check((await page.locator(`${block} pre`).first().textContent()).includes('|'.repeat(model.longest.length)), 'bars depict continuous match');
    }
    check(!/(?:ΔG|dimer Tm|hairpin Tm|발생 확률\s*\d|PASS|FAIL|Excellent|87\/100)/i.test(await text('#review-panel-complementarity')), 'no invented thermodynamics or score');
    await screenshot('04-complementarity.png', '#review-panel-complementarity');
    if (width !== 1440) await screenshot(width === 390 ? '04-mobile.png' : '04-tablet.png');

    await lens('off-target');
    check(await page.locator('#compare-abc').isDisabled(), 'C cannot be first reveal');
    check(await page.locator('#candidate-results').textContent() === '', 'candidate results initially absent');
    await open('#candidate-sequences');
    for (const p of Object.values(fixture.candidatePairs)) check((await text('#candidate-sequences')).includes(p.forward) && (await text('#candidate-sequences')).includes(p.reverse), 'fixed sequence provenance');
    await page.locator('#candidate-prediction').fill('P1과 P2의 차이를 서열과 산물에서 찾아본다.');
    await page.locator('#compare-ab').click();
    check(await page.locator('#candidate-results [data-source=C]').count() === 0, 'no C results in A/B stage');
    check(await page.locator('#review-c-question').isHidden(), 'C explanation deferred');
    await page.locator('#review-ab-reason').fill('A와 B의 길이만으로는 P1과 P2를 구별하기 어렵다.');
    await screenshot('04-ab-only.png', '#review-panel-off-target');
    await page.locator('#compare-abc').click();
    for (const [id, sources] of Object.entries(fixture.expectedExactMatchProducts)) for (const [source, products] of Object.entries(sources)) {
      const actual = await page.locator(`[data-candidate=${id}] [data-source=${source}]`).innerText();
      check(products.length ? products.every(([start, end, length]) => actual.includes(`${source} ${start}~${end}: ${length} bp`)) : actual.includes('예상 산물 없음'), `${id}/${source} exact fixture result`);
    }
    check(await page.locator('#review-p3-explanation').evaluate(el => !el.open), 'P3 answer stays collapsed');
    await page.locator('#candidate-judgment').fill('C에서 추가 산물이 예측되어 비교 범위를 넓혀야 함을 알았다.');
    await screenshot('04-with-background-c.png', '#review-panel-off-target');
    await page.locator('#candidate-negative').fill('결합 부위가 결실되었는지 먼저 구분해야 한다.');
    await open('#review-p3-explanation');
    check((await text('#review-p3-binding')).includes('121~140') && (await text('#review-p3-binding')).includes('결합: 0개'), 'P3 actual binding and deletion');
    await screenshot('04-p3-binding-failure.png', '#review-panel-off-target');
    await page.locator('#review-unresolved').fill('실제 반응 조건과 더 넓은 서열 데이터베이스에서의 결합을 확인한다.');
    check(JSON.stringify((await saved()).designs) === JSON.stringify(snapshot), '04 never changes saved snapshot');
    check(JSON.stringify((await saved()).draft) === JSON.stringify(protectedDraft), '04 never changes draft coordinates or sequences');
    check(await page.locator('#activity-04 a[href="#activity-05"]').count() === 1 && await page.locator('#activity-04 a[href*="primer-blast"]').count() === 0, 'existing learning order preserved');
    for (const name of REVIEW_LENSES) {
      await lens(name);
      check(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth), `${name} no viewport overflow`);
      check(await page.locator(`#review-panel-${name}`).isVisible(), 'selected panel visible');
      check(await page.locator(`#review-tab-${name}`).getAttribute('aria-controls') === `review-panel-${name}`, 'tab controls real panel');
    }
    await page.reload(); await page.waitForSelector('#pcr-worksheet[data-ready=true]');
    check(await page.locator('#review-design').inputValue() === '1', 'selected snapshot restored');
    check(await page.locator('#candidate-results [data-source=C]').count() === 3, 'C stage restored');
    check(await page.locator('#review-length-choice').inputValue() === 'no', 'choice restored');
    check((await page.locator('#review-unresolved').inputValue()).includes('데이터베이스'), 'answer restored');
    check(await page.locator('#review-tab-off-target').getAttribute('aria-selected') === 'true', 'lens restored');
    await open('#record-menu'); const downloadEvent = page.waitForEvent('download'); await page.locator('#export-record').click();
    const download = await downloadEvent, recordPath = resolve(output, `${engineName}-${width}-roundtrip.json`); await download.saveAs(recordPath);
    const exported = JSON.parse(await readFile(recordPath, 'utf8'));
    check(exported.schemaVersion === 1 && exported.review.stage === 2 && exported.review.designId === 1, 'optional v1 review export');
    await page.locator('#review-design').selectOption('draft');
    await open('#record-menu'); page.once('dialog', d => d.accept()); await page.locator('#import-record').setInputFiles(recordPath);
    await page.waitForFunction(() => document.querySelector('#review-design').value === '1');
    check(JSON.stringify((await saved()).review) === JSON.stringify(exported.review), 'review JSON restored');
    check(JSON.stringify((await saved()).designs) === JSON.stringify(snapshot), 'JSON saved snapshot integrity');

    // Legacy records omit the new optional object; old 04 answer ids remain usable.
    const legacy = { ...exported, draft: originalDraft }; delete legacy.review;
    await open('#record-menu'); page.once('dialog', d => d.accept());
    await page.locator('#import-record').setInputFiles({ name: 'legacy.json', mimeType: 'application/json', buffer: Buffer.from(JSON.stringify(legacy)) });
    await page.waitForFunction(() => document.querySelector('#review-design').value === 'draft');
    check(await page.locator('#review-tab-length-gc').getAttribute('aria-selected') === 'true', 'legacy starts with length/GC');
    check(await page.locator('#candidate-results').textContent() === '', 'legacy starts unrevealed');
    check((await page.locator('#candidate-judgment').inputValue()).includes('C에서'), 'legacy old 04 answer preserved');
    await page.reload(); await page.waitForSelector('#pcr-worksheet[data-ready=true]');
    check((await text('#review-my-design')).includes(pair.forward), 'legacy localStorage default valid');

    // Print all observations/answers, keep hidden results hidden on the blank sheet.
    await page.evaluate(() => { window.print = () => window.dispatchEvent(new Event('beforeprint')); });
    await open('#print-menu'); await page.locator('#print-filled').click(); await page.emulateMedia({ media: 'print' });
    check(await page.locator('#review-panel-tm').isVisible() && await page.locator('#review-panel-end').isVisible(), 'filled print includes inactive lenses');
    check((await text('#review-unresolved + .pcr-print-value')).includes('데이터베이스'), 'filled print answer');
    check((await text('#review-length-choice + .pcr-print-value')) === '결론 내릴 수 없다', 'print uses readable choice label');
    if (engineName === 'chromium' && width === 1440) await page.pdf({ path: resolve(output, 'phase3-filled.pdf'), format: 'A4', printBackground: true });
    await page.emulateMedia({ media: 'screen' }); await page.evaluate(() => window.dispatchEvent(new Event('afterprint')));
    await page.locator('#print-blank').click(); await page.emulateMedia({ media: 'print' });
    check(await page.locator('#review-my-design').isHidden() && await page.locator('#review-complementarity-values').isHidden(), 'blank print excludes computed design');
    check((await text('#review-unresolved + .pcr-print-value')) === '', 'blank print excludes answer');
    if (engineName === 'chromium' && width === 1440) await page.pdf({ path: resolve(output, 'phase3-blank.pdf'), format: 'A4', printBackground: true });
    await page.emulateMedia({ media: 'screen' }); await page.evaluate(() => window.dispatchEvent(new Event('afterprint')));
    check((await page.locator('#review-unresolved').inputValue()).includes('데이터베이스'), 'printing retains record');

    // Editing the draft invalidates review values; choosing a snapshot never repairs it.
    await open('#manual-sequences'); await page.locator('#primer-f').fill('N');
    check(await page.locator('#review-my-design').textContent() === '', 'invalid draft clears review');
    check((await text('#review-status')).includes('03에서 먼저'), 'invalid draft guidance');
    await page.locator('#review-design').selectOption('1');
    check((await text('#review-my-design')).includes(pair.forward), 'snapshot still reviewable');
    check((await saved()).draft.forward === 'N', 'snapshot review never changes invalid draft');
    await page.locator('#review-design').selectOption('draft');
    const longF = fixture.templates.A.sequence.slice(0, 100), longR = reverseComplement(fixture.templates.A.sequence.slice(220, 320));
    await page.locator('#primer-f').fill(longF); await page.locator('#primer-r').fill(longR); await lens('complementarity');
    check(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth), '100 nt alignments bounded on all viewports');
    await page.emulateMedia({ media: 'print' });
    check(await page.locator('.pcr-long-alignment > pre').first().isHidden(), 'long print alignment replaced with shared-column windows');
    check(await page.locator('.pcr-print-alignment pre').evaluateAll(elements => elements.every(el => el.textContent.split('\n').every(line => line.length <= 60))), 'long alignment print fragments preserve bounded columns');
    await page.emulateMedia({ media: 'screen' });
    await page.emulateMedia({ reducedMotion: 'reduce' });
    check(await page.locator('#review-panel-complementarity').evaluate(el => getComputedStyle(el).animationName) === 'none', 'reduced motion');
    check(await page.locator('#activity-04 input,#activity-04 textarea,#activity-04 select').evaluateAll(elements => elements.every(el => el.labels.length)), 'all review fields have labels');
    check(!/[\u00b7\u2022\u2027\u2219\u22c5\u30fb\u318d]/u.test(await text('#activity-04')), 'no forbidden dots');
    check(errors.length === 0, `no errors: ${errors.join(', ')}`);
    results.push({ engine: engineName, width, checks: checks - baseline, passed: true });
    console.log(`PASS Phase 3 ${engineName} ${width}px (${checks - baseline} assertions)`); await context.close();
  }
  const nojs = await browser.newContext({ javaScriptEnabled: false }), nojsPage = await nojs.newPage(); await nojsPage.goto(url);
  for (const lens of REVIEW_LENSES) check(await nojsPage.locator(`#review-panel-${lens}`).isVisible(), 'no-JS all questions readable');
  check((await nojsPage.locator('#activity-04').innerText()).includes('종이에'), 'no-JS paper fallback'); await nojs.close();
  const failure = await browser.newContext(), failurePage = await failure.newPage();
  await failurePage.route('**/pcr-primer-fixture.json', route => route.abort()); await failurePage.goto(url);
  await failurePage.waitForFunction(() => document.querySelector('#analysis-status').textContent.includes('로딩 실패'));
  check((await failurePage.locator('#review-status').innerText()).includes('데이터'), 'fixture failure notice');
  await failurePage.locator('#review-unresolved').fill('로딩 실패에도 기록'); await failurePage.locator('#review-tab-off-target').click();
  check(await failurePage.locator('#compare-ab').isDisabled() && await failurePage.locator('#compare-abc').isDisabled(), 'failed fixture controls disabled');
  await failurePage.locator('#record-menu > summary').click(); const d = failurePage.waitForEvent('download'); await failurePage.locator('#export-record').click(); await d; check(true, 'fixture failure still exports');
  await failure.close(); await browser.close();
}
await writeFile(resolve(output, 'verification.json'), JSON.stringify({ url, generatedAt: new Date().toISOString(), checks, results, evidence }, null, 2));

// Contact sheet uses the actual captured PNGs, preserving aspect ratios and readable scale.
const browser = await chromium.launch({ headless: true });
for (const [filename, list, cellWidth] of [
  ['phase3-contact-sheet.png', evidence.filter(e => e.width === 1440 && !['04-mobile.png', '04-tablet.png'].includes(e.name)), 940],
  ['phase3-contact-sheet-mobile.png', evidence.filter(e => e.name === '04-mobile.png' || e.name === '04-tablet.png'), 720]
]) {
  const page = await browser.newPage({ viewport: { width: cellWidth * 2 + 72, height: 1150 }, deviceScaleFactor: 1 });
  const cards = await Promise.all(list.map(async item => `<figure><figcaption>${item.name}</figcaption><img src="data:image/png;base64,${(await readFile(item.path)).toString('base64')}"></figure>`));
  await page.setContent(`<html><head><style>body{margin:24px;background:#fff;color:#2f3634;font:20px sans-serif}main{display:grid;grid-template-columns:repeat(2,${cellWidth}px);gap:24px}figure{margin:0;border-top:1px solid #dedfdd;padding-top:12px}figcaption{margin-bottom:14px}img{display:block;max-width:100%;height:auto}</style></head><body><main>${cards.join('')}</main></body></html>`);
  await page.locator('img').evaluateAll(images => Promise.all(images.map(img => img.decode())));
  await page.screenshot({ path: resolve(output, filename), fullPage: true }); await page.close();
}
await browser.close(); console.log(`PASS Phase 3 ${checks} assertions. PNGs and contact sheets: ${output}`);
