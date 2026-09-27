// Install Playwright and axe-core in excluded tests/.pcr-tools before running.
import assert from 'node:assert/strict';
import { mkdir, readFile, writeFile } from 'node:fs/promises';
import { resolve } from 'node:path';
import { chromium, webkit } from './.pcr-tools/node_modules/playwright/index.mjs';
import { emptyRecord, STORAGE_KEY } from '../assets/js/pcr-records.mjs';
import { EVIDENCE_CASES, caseSummary } from '../assets/js/pcr-evidence.mjs';
const fixture = JSON.parse(await readFile(new URL('../assets/data/pcr-primer-fixture.json', import.meta.url)));
const url = process.env.PCR_TEST_URL || 'http://127.0.0.1:4173/bioinformatics/pcr-primer-design/';
const output = resolve(process.env.PCR_VERIFICATION_ROOT || 'verification.local/pcr-primer-design', 'phase4'); await mkdir(output, { recursive: true });
let checks = 0; const results = [], evidence = [], accessibility = [];
const check = (value, label) => { assert.ok(value, label); checks++; };
for (const [engineName, engine] of Object.entries({ chromium, webkit })) {
  const browser = await engine.launch({ headless: true });
  for (const width of [1440, 768, 390]) {
    const baseline = checks, context = await browser.newContext({ viewport: { width, height: 1150 }, hasTouch: width === 390 });
    const page = await context.newPage(), errors = []; page.on('pageerror', e => errors.push(e.message));
    await page.goto(url); await page.waitForSelector('#pcr-worksheet[data-ready=true]');
    const open = async id => { if (!await page.locator(id).evaluate(el => el.open)) await page.locator(`${id} > summary`).click(); };
    const saved = async () => JSON.parse(await page.evaluate(key => localStorage.getItem(key), STORAGE_KEY));
    const text = id => page.locator(id).innerText();
    const shot = async (name, selector) => {
      if (engineName !== 'chromium' || (width !== 1440 && name !== '05-mobile.png' && name !== '05-tablet.png')) return;
      await page.locator(selector).screenshot({ path: resolve(output, name), animations: 'disabled' });
      evidence.push({ name, path: resolve(output, name), width, selector });
    };
    const fine = async (name, start, end) => {
      await page.locator(`[data-select-primer=${name}]`).click(); await open('#coordinate-details');
      await page.locator('#range-start').fill(String(start)); await page.locator('#range-end').fill(String(end)); await page.locator('#apply-range').click();
    };
    check((await text('#gel-results')).includes('유효한 primer pair'), 'empty student design guidance');
    check(await page.locator('#evidence-lane-sample .pcr-gel-band').count() === 1, 'fixed case available without student design');
    check(!await page.locator('#evidence-additional').evaluate(el => el.open), 'additional case initially collapsed');
    await page.locator('#evidence-case-2').click(); await page.locator('#evidence-lane-ntc').click();
    check((await text('#evidence-selected-lane')).includes('짧은 band'), 'NTC can be observed without student design');
    await page.locator('#evidence-case-1').click();
    await fine('F', 41, 60); await fine('R', 281, 300);
    check((await text('#gel-results')).includes('260 bp') && (await text('#gel-results')).includes('180 bp') && (await text('#gel-results')).includes('예상 산물 없음'), '03 computed A/B/C prediction');
    await shot('05-overview.png', '#evidence-workspace');
    await page.locator('#save-design').click(); await page.locator('#review-design').selectOption('1');
    await fine('F', 121, 140);
    check((await text('#evidence-design-source')).includes('저장 설계 1'), '05 follows 04 explicit snapshot');
    check((await text('#gel-results')).includes('260 bp'), 'snapshot persists after editing draft');
    await page.locator('#review-design').selectOption('draft');
    check(!(await text('#gel-results')).includes('260 bp') && (await text('#gel-results')).includes('180 bp'), '05 follows actual alternate draft');
    await page.locator('#review-design').selectOption('1');
    await page.locator('#review-length-reason').fill('04 기록을 보존한다.');
    const before = await saved();
    const protectedState = JSON.stringify([before.draft, before.designs, before.review, before.answers['review-length-reason']]);

    for (const c of EVIDENCE_CASES) {
      if (c.additional) await open('#evidence-additional');
      await page.locator(`#evidence-${c.id}`).click();
      check(await page.locator(`#evidence-${c.id}`).getAttribute('aria-pressed') === 'true', 'selected case aria');
      check((await text('#evidence-case-expected')).includes('260 bp') && (await text('#evidence-case-expected')).includes('독립'), 'fixed expected size clearly independent');
      check(await page.locator('#evidence-observation-summary').textContent() === caseSummary(c), 'complete equivalent text for gel');
      for (const [id, lane] of Object.entries(c.lanes)) {
        const button = page.locator(`#evidence-lane-${id}`);
        assert.deepEqual(await button.locator('.pcr-gel-band').evaluateAll(els => els.map(el => Number(el.dataset.bp))), lane.bands); checks++;
        check((await button.getAttribute('aria-label')).includes(lane.label) && (await button.getAttribute('aria-label')).includes(lane.observation), 'lane accessible name and observation');
        const box = await button.boundingBox(); check(box.width >= 44 && box.height >= 44, 'lane touch target');
      }
      check(await page.locator('#evidence-gel-figure').evaluate(el => el.textContent.includes('정량값을 뜻하지')), 'brightness warning');
      check(!/(?:오염 확정|primer-dimer 확정|target 확인 완료|PASS|FAIL|confidence score)/i.test(await text('#evidence-case-view')), 'no automatic cause or identity verdict');
      check(await page.locator(`#evidence-${c.id}-answers .pcr-explanation`).evaluate(el => !el.open), 'explanation not shown automatically');
      const lane = c.id === 'case-2' ? 'ntc' : c.id === 'case-3' ? 'positive' : 'sample';
      if (c.id === 'case-1') await shot('05-case1-expected-band.png', '#evidence-case-view');
      if (c.id === 'case-2') await shot('05-case2-ntc-band.png', '#evidence-case-view');
      if (c.id === 'case-3') await shot('05-case3-positive-control-fail.png', '#evidence-case-view');
      const button = page.locator(`#evidence-lane-${lane}`);
      await button.focus(); await page.keyboard.press('Enter');
      check(await button.getAttribute('aria-pressed') === 'true', 'keyboard selects lane');
      check(await button.evaluate(el => el === document.activeElement && getComputedStyle(el).outlineStyle !== 'none'), 'focus retained and visible');
      check(await page.locator('#evidence-gel [aria-pressed=true]').count() === 1, 'one selected lane');
      check((await text('#evidence-selected-lane')).includes(c.lanes[lane].observation), 'selected observation');
      check((await button.innerText()).includes('선택됨'), 'selected state not color alone');
      if (width === 390) { await page.locator('#evidence-lane-marker').tap(); await button.tap(); }
      else { await page.locator('#evidence-lane-marker').click(); await button.click(); }
      if (c.id === 'case-1') await shot('05-case1-lane-selected.png', '#evidence-case-view');
      if (c.id === 'case-3') await shot('05-case3-control-selected.png', '#evidence-case-view');
      const answers = {
        observation: `${c.lanes[lane].label}: ${c.lanes[lane].observation}`,
        interpretation: c.id === 'case-2' ? 'NTC에서 band가 보인다. 오염이나 primer 유래 산물 등 여러 설명을 비교해야 한다.' : c.id === 'case-3' ? 'Positive control에서도 band가 보이지 않아 Sample 음성을 해석하는 데 제한이 있다.' : '크기와 대조군 관찰을 함께 고려하여 가능한 설명을 비교한다.',
        uncertainty: c.id === 'case-2' ? '작은 크기만으로 primer-dimer라고 단정할 수 없다. 시약과 작업 과정을 확인하고 대조 반응을 반복한다.' : 'Gel만으로 산물의 sequence identity나 원인을 하나로 확정할 수 없다.'
      };
      for (const [key, value] of Object.entries(answers)) {
        await page.locator(`#evidence-${c.id}-${key}`).fill(value);
        check((await saved()).answers[`evidence-${c.id}-${key}`] === value, 'case answer saved');
      }
      if (c.id === 'case-2') await shot('05-case2-interpretation.png', '#evidence-case-2-answers');
      if (c.id === 'case-4') await shot('05-case4-multiple-bands.png', '#evidence-case-view');
      check(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth), 'case no horizontal overflow');
      await page.addScriptTag({ path: resolve('tests/.pcr-tools/node_modules/axe-core/axe.min.js') });
      const audit = await page.evaluate(async () => axe.run(document.querySelector('#activity-05'), { runOnly: { type: 'tag', values: ['wcag2a','wcag2aa','wcag21aa','wcag22aa'] } }));
      accessibility.push({ engineName, width, caseId: c.id, violations: audit.violations, incomplete: audit.incomplete.map(r => r.id) });
      check(audit.violations.length === 0, `axe ${engineName}/${width}/${c.id}: ${audit.violations.map(v => v.id + ': ' + v.nodes.map(n => n.target).join(', ')).join('; ')}`);
    }
    const after = await saved();
    check(JSON.stringify([after.draft, after.designs, after.review, after.answers['review-length-reason']]) === protectedState, '05 preserves 03/04 state and answers');
    await page.locator('#evidence-case-3').click();
    await page.locator('#evidence-identity').fill('Band의 크기는 길이에 관한 증거다. 길이가 같아도 sequence는 다를 수 있으므로 정체를 확정할 수 없다.');
    await page.locator('#evidence-controls').fill('대조군은 시료 결과를 해석할 조건을 확인하게 한다. Positive control과 NTC의 관찰을 함께 보아야 한다.');
    await shot('05-final-reflection.png', '#evidence-final');
    if (width === 390) await shot('05-mobile.png', '#evidence-workspace');
    if (width === 768) await shot('05-tablet.png', '#evidence-workspace');
    await page.reload(); await page.waitForSelector('#pcr-worksheet[data-ready=true]');
    check(await page.locator('#evidence-case-3').getAttribute('aria-pressed') === 'true', 'reload active case');
    check(await page.locator('#evidence-lane-positive').getAttribute('aria-pressed') === 'true', 'reload selected lane');
    check((await page.locator('#evidence-case-3-interpretation').inputValue()).includes('제한'), 'reload case answer');
    check((await page.locator('#evidence-case-2-uncertainty').inputValue()).includes('시약'), 'reload inactive case uncertainty');
    check((await page.locator('#evidence-identity').inputValue()).includes('sequence'), 'reload final reflection');
    await open('#record-menu'); const downloadEvent = page.waitForEvent('download'); await page.locator('#export-record').click();
    const file = resolve(output, `${engineName}-${width}-roundtrip.json`); await (await downloadEvent).saveAs(file);
    const exported = JSON.parse(await readFile(file, 'utf8'));
    check(exported.schemaVersion === 1 && exported.evidence.activeCase === 'case-3', 'v1 export');
    await page.locator('#evidence-case-1').click(); await page.locator('#evidence-case-1-interpretation').fill('임시 수정');
    await open('#record-menu'); page.once('dialog', d => d.accept()); await page.locator('#import-record').setInputFiles(file);
    await page.waitForFunction(() => document.querySelector('#evidence-case-3').getAttribute('aria-pressed') === 'true');
    check(JSON.stringify((await saved()).evidence) === JSON.stringify(exported.evidence), 'JSON view round trip');
    check(JSON.stringify((await saved()).answers) === JSON.stringify(exported.answers), 'JSON every answer round trip');

    await page.evaluate(() => { window.print = () => window.dispatchEvent(new Event('beforeprint')); });
    await open('#print-menu'); await page.locator('#print-filled').click(); await page.emulateMedia({ media: 'print' });
    for (const c of EVIDENCE_CASES) {
      check(await page.locator(`#evidence-${c.id}-answers`).isVisible(), 'print every case including inactive');
      check((await text(`#evidence-${c.id}-interpretation + .pcr-print-value`)) === exported.answers[`evidence-${c.id}-interpretation`], 'print case answers');
    }
    check(await page.locator('#evidence-lane-positive').isVisible(), 'printed lane despite shared button hide');
    check((await text('#evidence-case-heading')).includes('상황 3'), 'print current case');
    check((await text('#gel-results')).includes('260 bp'), 'print student prediction');
    check(await page.locator('#evidence-gel-figure').evaluate(el => getComputedStyle(el).breakInside) === 'avoid', 'gel and labels avoid page split');
    if (engineName === 'chromium' && width === 1440) await page.pdf({ path: resolve(output,'phase4-filled.pdf'), format:'A4', printBackground:false });
    await page.emulateMedia({ media:'screen' }); await page.evaluate(() => window.dispatchEvent(new Event('afterprint')));
    await page.locator('#print-blank').click(); await page.emulateMedia({ media:'print' });
    check(await page.locator('#evidence-prediction').isHidden(), 'blank hides student calculation');
    check(await page.locator('#activity-05 .pcr-print-value').evaluateAll(els => els.every(el => !el.textContent)), 'blank contains no student answers');
    check(await page.locator('#evidence-lane-positive .pcr-gel-selection').isHidden(), 'blank hides lane selection');
    check(await page.locator('#evidence-case-1-answers .pcr-evidence-print-summary').isVisible(), 'blank retains question observations');
    if (engineName === 'chromium' && width === 1440) await page.pdf({ path: resolve(output,'phase4-blank.pdf'), format:'A4', printBackground:false });
    await page.emulateMedia({ media:'screen' }); await page.evaluate(() => window.dispatchEvent(new Event('afterprint')));
    check(JSON.stringify((await saved()).answers) === JSON.stringify(exported.answers), 'printing does not alter answers');
    await page.locator('#review-design').selectOption('draft'); await open('#manual-sequences'); await page.locator('#primer-f').fill('N');
    check(!(await text('#gel-results')).includes('180 bp') && (await text('#gel-results')).includes('유효한'), 'invalid student design clears old prediction');
    check(await page.locator('#evidence-lane-marker .pcr-gel-band').count() === 5, 'invalid design leaves fixed observations intact');
    await page.emulateMedia({ reducedMotion:'reduce' });
    check(await page.locator('#evidence-gel').evaluate(el => getComputedStyle(el).animationName) === 'none', 'reduced motion');
    check(errors.length === 0, `no runtime errors ${errors.join(',')}`);
    results.push({ engineName, width, checks:checks-baseline, passed:true });
    console.log(`PASS Phase 4 ${engineName} ${width}px (${checks-baseline} assertions)`); await context.close();
  }
  // Old v1 on both input paths, all three original activity 05 answers kept.
  const legacy = emptyRecord(); delete legacy.evidence;
  legacy.draft = { ...fixture.candidatePairs.P1 };
  legacy.answers = { 'evidence-cases':'이전 통합 기록 <img src=x onerror="window.injected=true">', 'evidence-controls':'이전 대조군', 'evidence-identity':'이전 정체', 'candidate-judgment':'04 보존' };
  const old = await browser.newContext();
  await old.addInitScript(({key, record}) => { if (!localStorage.getItem(key)) localStorage.setItem(key, JSON.stringify(record)); }, {key:STORAGE_KEY,record:legacy});
  const page = await old.newPage(); await page.goto(url); await page.waitForSelector('#pcr-worksheet[data-ready=true]');
  check(await page.locator('#evidence-case-1').getAttribute('aria-pressed') === 'true', 'old localStorage safe case');
  await page.locator('#evidence-legacy > summary').click();
  check(await page.locator('#evidence-cases').inputValue() === legacy.answers['evidence-cases'], 'old 05 text preserved and accessible');
  for (const id of ['evidence-controls','evidence-identity','candidate-judgment']) check(await page.locator(`#${id}`).inputValue() === legacy.answers[id], 'old other answers preserved');
  await page.locator('#record-menu > summary').click(); page.once('dialog',d=>d.accept());
  await page.locator('#import-record').setInputFiles({name:'old-v1.json',mimeType:'application/json',buffer:Buffer.from(JSON.stringify(legacy))});
  await page.waitForFunction(() => document.querySelector('#save-status').textContent.includes('자동 저장'));
  check(await page.locator('#evidence-cases').inputValue() === legacy.answers['evidence-cases'], 'old JSON legacy answer retained');
  check(!await page.evaluate(() => window.injected), 'legacy text stays inert');
  await page.locator('#evidence-legacy > summary').click();
  await page.evaluate(() => { window.print = () => window.dispatchEvent(new Event('beforeprint')); });
  await page.locator('#print-menu > summary').click(); await page.locator('#print-filled').click(); await page.emulateMedia({media:'print'});
  check(await page.locator('#evidence-cases + .pcr-print-value').isVisible(), 'closed legacy included in filled print');
  check((await page.locator('#evidence-cases + .pcr-print-value').innerText()).includes('이전 통합'), 'legacy print text');
  await old.close();
  const nojs = await browser.newContext({javaScriptEnabled:false}), staticPage = await nojs.newPage(); await staticPage.goto(url);
  check((await staticPage.locator('#activity-05 noscript').innerText()).includes('상황 3'), 'no-JS observation fallback');
  check(await staticPage.locator('#evidence-identity').isVisible(), 'no-JS synthesis question'); await nojs.close();
  const failure = await browser.newContext(), offline = await failure.newPage();
  await offline.route('**/pcr-primer-fixture.json', route => route.abort()); await offline.goto(url);
  await offline.waitForFunction(() => document.querySelector('#analysis-status').textContent.includes('로딩 실패'));
  await offline.locator('#evidence-case-2').click(); await offline.locator('#evidence-lane-ntc').click();
  check(await offline.locator('#evidence-lane-ntc .pcr-gel-band').count() === 1, 'fixture failure leaves fixed cases usable');
  await offline.locator('#evidence-case-2-uncertainty').fill('자료 로딩 실패에도 관찰을 기록');
  check(JSON.parse(await offline.evaluate(key => localStorage.getItem(key), STORAGE_KEY)).answers['evidence-case-2-uncertainty'].includes('기록'), 'fixture failure saves evidence');
  await failure.close(); await browser.close();
}
await writeFile(resolve(output,'verification.json'), JSON.stringify({url,generatedAt:new Date().toISOString(),checks,results,evidence,accessibility},null,2));
const browser = await chromium.launch({headless:true});
for (const [filename, list, cellWidth] of [
  ['phase4-contact-sheet.png',evidence.filter(e=>e.width===1440),940],
  ['phase4-contact-sheet-mobile.png',evidence.filter(e=>e.width!==1440),720]
]) {
  const page = await browser.newPage({viewport:{width:cellWidth*2+72,height:1150},deviceScaleFactor:1});
  const panels = await Promise.all(list.map(async e=>`<figure><figcaption>${e.name}</figcaption><img src="data:image/png;base64,${(await readFile(e.path)).toString('base64')}"></figure>`));
  await page.setContent(`<html><head><style>body{margin:24px;background:#fff;color:#2f3634;font:20px sans-serif}main{display:grid;grid-template-columns:repeat(2,${cellWidth}px);gap:24px}figure{margin:0;border-top:1px solid #dedfdd;padding-top:12px}figcaption{margin-bottom:14px}img{display:block;max-width:100%;height:auto}</style></head><body><main>${panels.join('')}</main></body></html>`);
  await page.locator('img').evaluateAll(els=>Promise.all(els.map(el=>el.decode())));
  await page.screenshot({path:resolve(output,filename),fullPage:true}); await page.close();
}
await browser.close(); console.log(`PASS Phase 4 ${checks} assertions. PNGs and contact sheets: ${output}`);
