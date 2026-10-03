// Phase 9: one student's complete journey against the actual Jekyll output.
import assert from 'node:assert/strict';
import { mkdir, readFile, writeFile } from 'node:fs/promises';
import { resolve } from 'node:path';
import { chromium, webkit } from './.pcr-tools/node_modules/playwright/index.mjs';
import { STORAGE_KEY } from '../assets/js/pcr-records.mjs';

const base = process.env.PCR_SITE_URL || 'http://127.0.0.1:4174';
const url = `${base}/bioinformatics/pcr-primer-design/`;
const output = resolve('verification.local/pcr-primer-design/phase9');
await mkdir(output, { recursive: true });
const fixture = JSON.parse(await readFile('assets/data/pcr-primer-fixture.json', 'utf8'));
const results = [], screenshots = [], audits = [];
let representative;
const stable = record => { const copy = structuredClone(record); delete copy.updatedAt; return copy; };

for (const [engineName, engine] of Object.entries({ chromium, webkit })) {
  const browser = await engine.launch({ headless: true });
  try {
    for (const width of [1440, 768, 390]) {
      const context = await browser.newContext({ viewport: { width, height: 1100 }, hasTouch: width === 390 });
      const page = await context.newPage();
      const consoleErrors = [], warnings = [], pageErrors = [], failedRequests = [], badResponses = [], unhandledRejections = [];
      page.on('console', m => { if (m.type() === 'error') consoleErrors.push(m.text()); if (m.type() === 'warning') warnings.push(m.text()); });
      page.on('pageerror', e => pageErrors.push(e.message));
      page.on('requestfailed', r => failedRequests.push({ url: r.url(), error: r.failure()?.errorText }));
      page.on('response', r => { if (r.status() >= 400) badResponses.push({ url: r.url(), status: r.status() }); });
      await context.exposeBinding('releaseRejection', (_, reason) => unhandledRejections.push(reason));
      await context.addInitScript(() => { window.releaseUnhandled = []; window.addEventListener('unhandledrejection', e => { window.releaseUnhandled.push(String(e.reason)); window.releaseRejection(String(e.reason)); }); });
      const ready = () => page.waitForSelector('#pcr-worksheet[data-ready=true]');
      const saved = () => page.evaluate(key => JSON.parse(localStorage.getItem(key)), STORAGE_KEY);
      const press = async (selector, key = 'Enter') => { assert.ok(await page.locator(selector).isVisible(), `keyboard target visible: ${selector}`); await page.locator(selector).focus(); await page.keyboard.press(key); };
      const open = async selector => { if (!await page.locator(selector).evaluate(e => e.open)) await press(`${selector} > summary`); };
      const importRecord = async (file, expected) => {
        await open('#record-menu'); const dialog = page.waitForEvent('dialog');
        await page.locator('#import-record').setInputFiles(file); await (await dialog).accept();
        await page.waitForFunction(({ key, expected }) => {
          const canonical = v => v && typeof v === 'object' ? Array.isArray(v) ? v.map(canonical) : Object.fromEntries(Object.keys(v).filter(k => k !== 'updatedAt').sort().map(k => [k, canonical(v[k])])) : v;
          return JSON.stringify(canonical(JSON.parse(localStorage.getItem(key)))) === JSON.stringify(canonical(expected));
        }, { key: STORAGE_KEY, expected });
      };
      const fill = (selector, value) => page.locator(selector).fill(value);
      const checkOverflow = async () => assert.ok(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth), `${engineName}/${width}: page overflow`);
      const coordinate = async (name, start, end) => {
        await press(`[data-select-primer=${name}]`); await open('#coordinate-details');
        await fill('#range-start', String(start)); await fill('#range-end', String(end)); await press('#apply-range');
      };
      const audit = async state => {
        await page.addScriptTag({ path: resolve('tests/.pcr-tools/node_modules/axe-core/axe.min.js') });
        const a = await page.evaluate(async () => {
          const r = await axe.run('#pcr-worksheet', { runOnly: { type: 'tag', values: ['wcag2a', 'wcag2aa', 'wcag21aa', 'wcag22aa', 'best-practice'] } });
          return { violations: r.violations.map(v => ({ id: v.id, nodes: v.nodes.map(n => n.target) })), incomplete: r.incomplete.map(v => ({ id: v.id, nodes: v.nodes.map(n => n.target) })) };
        });
        audits.push({ engineName, width, state, ...a }); assert.deepEqual(a.violations, []); await checkOverflow();
      };

      // A fresh browser context has empty localStorage. Enter by the real portal link.
      assert.equal((await page.goto(base + '/')).status(), 200);
      assert.equal(await saved(), null);
      const portalHeading = page.locator('#bioinformatics h3').filter({ hasText: 'PCR과 프라이머 디자인' });
      assert.match((await portalHeading.innerText()).replace(/\s+/g, ' '), /^C2\s*PCR과 프라이머 디자인$/);
      assert.doesNotMatch(await portalHeading.locator('..').innerText(), /\d+분/);
      await portalHeading.locator('a[href="/bioinformatics/pcr-primer-design/"]').click(); await ready();
      assert.equal(page.url(), url);
      assert.doesNotMatch(await page.locator('body').innerText(), /이전 (?:기록 확인|작성 기록)|legacy|old v1|migration/i);
      assert.equal(await page.locator('[data-legacy-answers]:visible,#ext-legacy:visible,#final-legacy:visible,#evidence-legacy:visible').count(), 0);
      await audit('fresh');

      // 00: use the actual keyboard range interaction, retaining relative positions.
      await press('#prediction-forward', 'Home'); await page.keyboard.press('ArrowRight'); await page.keyboard.press('ArrowRight');
      await press('#prediction-reverse', 'End'); for (let i = 0; i < 3; i++) await page.keyboard.press('ArrowLeft');
      await fill('#first-placement', '결실 양옆에 primer를 두면 A와 B의 산물 길이를 비교할 수 있을 것으로 예상했다.');
      const initial = (await saved()).initialPrimerPrediction;
      assert.deepEqual(initial, { reference: 'A', units: 'relative-percent', forward: 15, reverse: 80 });
      // 01 and 02: visit every stage and early cycle, then complement -> order direction.
      for (const stage of ['mixture', 'denaturation', 'annealing', 'extension']) await press(`[data-cycle-stage=${stage}]`);
      for (const cycle of [1, 2, 3]) await press(`[data-comparison-cycle="${cycle}"]`);
      await press('[name=cycle-boundary-choice][value=primer-pair]', 'Space');
      await fill('#cycle-template', '2주기에 양 끝이 정해진 단일가닥, 3주기에 정확한 길이의 이중가닥이 생기며 긴 산물도 남는다.');
      await fill('#direction-complement', 'TCAGGCAT'); await press('#check-complement'); await press('#reverse-complement');
      await fill('#direction-reverse', 'TACGGACT'); await press('#check-direction');
      await press('#flip-strands'); await press('[data-arrangement=parallel]'); await press('[data-arrangement=inward]');
      await press('[name=direction-end-choice][value="3"]', 'Space');
      assert.match(await page.locator('#direction-feedback').innerText(), /일치/);

      // 03: two saved snapshots, including a deletion-overlapping binding site.
      await coordinate('F', 41, 60); await coordinate('R', 281, 300);
      assert.equal((await saved()).draft.forward, fixture.candidatePairs.P1.forward);
      assert.equal((await saved()).draft.reverse, fixture.candidatePairs.P1.reverse);
      await fill('#design-reason', '결실 양옆의 결합 부위로 A 260 bp와 B 180 bp를 비교한다.'); await press('#save-design');
      const design1 = structuredClone((await saved()).designs[0]);
      await coordinate('F', 111, 130); await fill('#design-reason', 'F가 결실 경계와 겹치면 B의 완전 일치 결합 부위가 사라지는 대안을 비교했다.'); await press('#save-design');
      assert.match(await page.locator('#analysis-results').innerText(), /190/);
      assert.match(await page.locator('#binding-message').innerText(), /완전 일치 결합 부위가 존재하지/);
      assert.deepEqual((await saved()).designs[0], design1);
      await page.getByRole('button', { name: '설계 1 불러오기', exact: true }).click();
      const designs = structuredClone((await saved()).designs);
      assert.equal(designs.length, 2);

      // 04: all five lenses and A/B then C; keep the reviewed source explicit.
      await page.locator('#review-design').selectOption('1');
      for (const lens of ['length-gc', 'tm', 'end', 'complementarity', 'off-target']) {
        await press(`#review-tab-${lens}`); await checkOverflow();
        if (lens === 'length-gc') await page.locator('#review-length-choice').selectOption('no');
        if (lens === 'end') await fill('#review-end-observation', '3′ 말단과 상보성을 함께 확인하며 GC clamp 하나로 판정하지 않는다.');
      }
      await press('#compare-ab'); await fill('#review-ab-reason', 'A/B만으로 P1과 P2의 배경 결합을 구별할 수 없다.'); await press('#compare-abc');
      await fill('#candidate-judgment', 'C에서는 P2의 220 bp 산물이 예측된다. P3의 B 음성은 target 부재의 증명이 아니다.');
      await fill('#review-unresolved', '더 넓은 서열 범위와 실제 반응 조건을 별도로 확인해야 한다.');

      // 05: observations are fixed classroom evidence, independent of the student's pair.
      for (const n of [1, 2, 3]) {
        await press(`[data-evidence-case="case-${n}"]`);
        for (const lane of ['sample', 'positive', 'ntc']) await press(`#evidence-lane-${lane}`);
        const fields = page.locator('#evidence-case-answers textarea:visible, #evidence-case-answers input:visible');
        for (let i = 0; i < await fields.count(); i++) await fields.nth(i).fill(['크기는 정체를 확정하지 않는다.', 'NTC band 원인은 추가 확인이 필요하다.', 'Positive control 실패로 Sample 음성을 단정할 수 없다.'][n - 1]);
      }
      await fill('#evidence-identity', '같은 크기의 다른 서열도 같은 위치에 band를 만들 수 있다.');
      await fill('#evidence-controls', 'Positive control과 NTC를 함께 확인해 음성과 비특이적 산물의 해석 범위를 판단한다.');

      // 06: synthetic QA transcription only; never submit a real NCBI search.
      await press('#ext-route-mine'); await press('#ext-to-plan');
      for (const [key, value] of Object.entries({ purpose: '외부 기록 흐름 연습 / 실제 검색 미실시', organism: '교육용 인공 서열', target: 'A/B/C 수업 자료', database: 'Custom / 교육용 입력 서열', notes: 'Release QA용 모의 기록이며 실제 NCBI 결과가 아니다.' })) await fill(`[data-external="plan.${key}"]`, value);
      await press('#ext-performed');
      for (const [key, value] of Object.entries({ date: '2026-09-27', organism: '제한 없음 / Custom 입력 범위', database: 'Custom / 교육용 입력 서열', target: '교육용 A/B/C / 모의 결과', forward: design1.forward, reverse: design1.reverse, specificity: '완전 일치 / 모의 검색 조건', other: '실제 NCBI 검색은 실행하지 않은 QA용 입력이다.' })) await fill(`[data-external="conditions.${key}"]`, value);
      await open('#ext-results');
      for (const [key, value] of Object.entries({ forward: design1.forward, reverse: design1.reverse, product: '260', tmF: '59.2', tmR: '59.8', observations: '실제 검색 미실시 / QA용 모의 보고', other: '입력과 복원 검증용' })) await fill(`[data-external="candidates.0.${key}"]`, value);
      await page.locator('[data-external="candidates.0.unintended"]').selectOption('none');
      assert.equal((await saved()).externalSearch.selectedCandidate, '');
      await page.locator('#ext-selectedCandidate').selectOption('A');
      await fill('#ext-selectionReason', '입력한 조건과 결과를 함께 보존하기 위한 연습 후보이다.');
      await fill('#ext-claimReflection', '기록된 입력 범위와 조건에서만 해석하며 다른 genome 전체로 일반화하지 않는다.');
      await fill('#ext-wetLabReflection', '실제 PCR 성공과 산물 정체는 아직 알 수 없다.'); await press('#ext-complete');
      assert.equal((await saved()).externalSearch.status, 'recorded');
      assert.match(await page.locator('#ext-claim-text').innerText(), /Custom.*Candidate A/);

      // 07: selecting the current draft freezes it; later draft edits cannot alter it.
      await press('#final-select-draft', 'Space'); const draftSnapshot = structuredClone((await saved()).finalReview.primerSnapshot);
      await coordinate('F', 61, 80); assert.deepEqual((await saved()).finalReview.primerSnapshot, draftSnapshot);
      await press('#final-select-design-1', 'Space');
      await page.getByRole('button', { name: '설계 1 불러오기', exact: true }).click();
      for (const [id, value] of Object.entries({
        'final-question': '80 bp 결실을 가진 두 시료를 산물 길이로 구별할 수 있는가?',
        'final-evidence': '설계 1은 결실 양옆에 결합하며 내부 모형에서 A 260 bp, B 180 bp가 예측된다.',
        'final-revision': '대략적 위치 예측을 실제 서열과 좌표로 바꾸고 결실 경계에 겹치는 대안을 비교했다.',
        'final-control-positive': '알려진 양성 template와 같은 primer/조건을 사용한다.',
        'final-control-negative': 'Template 없는 NTC를 함께 둔다.',
        'final-control-additional': '필요하면 추출과 반응 억제 대조를 검토한다.',
        'final-control-limit': 'Positive control 실패 시 sample의 target 부재를 결론 내릴 수 없다.',
        'final-unknown': '실제 증폭과 산물의 서열 정체는 미확인이다.',
        'final-assessment': '교육용 내부 모형의 후보 설계이며 실제 PCR로 검증한 pair는 아니다.'
      })) await fill(`#${id}`, value);
      await press('#final-notebook-toggle');
      assert.match(await page.locator('#final-notebook').innerText(), /미수행/);
      assert.deepEqual((await saved()).designs, designs);
      await audit('core-notebook');

      // RNA: every strategy, both controls and all five candidate binding structures.
      await open('#rna-extension');
      for (const step of ['rna', 'rt', 'cdna', 'pcr']) await press(`[data-rna-option=flowStep][data-value=${step}]`);
      await page.locator('#rna-template').selectOption('cdna');
      await press('[data-rna-section=genomic]'); await press('[data-rna-section=strategy]');
      for (const strategy of ['same-exon', 'intron-spanning', 'junction']) {
        await press(`[data-rna-option=primerStrategy][data-value=${strategy}]`);
        if (strategy === 'intron-spanning') for (const value of ['short', 'long']) await press(`[data-rna-option=intronExample][data-value=${value}]`);
        if (strategy === 'junction') await press('[data-rna-option=junctionPosition][data-value=junction]');
        await checkOverflow();
      }
      await fill('#rna-strategy-answer', 'Same-exon은 gDNA도 증폭 가능하다. Intron과 junction 전략도 조건과 전체 서열 맥락을 확인한다.');
      await press('[data-rna-section=control]');
      for (const value of ['case1', 'case2']) await press(`[data-rna-option=rtControlCase][data-value=${value}]`);
      await fill('#rna-control-answer', 'RT(-) signal은 DNA 유래 가능성을 검토할 근거이며 gDNA 오염 확정은 아니다.');
      await press('[data-rna-section=transcripts]');
      for (const value of ['exon2', 'exon3', 'exon4', 'junction23', 'junction34']) await press(`[data-rna-option=transcriptTarget][data-value=${value}]`);
      await fill('#rna-transcripts-answer', '선택 exon 또는 인접 junction과 공통 R 부위를 확인하며 transcript별 검출 범위를 구별한다.');
      await fill('#rna-specific-answer', '후보 junction이 다른 transcript에도 있는지와 전체 amplicon 구조를 함께 검토한다.');
      await audit('rna-transcripts');
      await press('[data-rna-section=reflection]');
      await fill('#rna-reflection-answer', 'RNA에서 cDNA를 만들고, gDNA 가능성과 RT(-)를 고려한다. 두 primer의 결합 부위 및 transcriptome/genome 맥락을 확인해야 한다.');
      const before = await saved();
      assert.deepEqual(before.initialPrimerPrediction, initial); assert.deepEqual(before.designs, designs);
      assert.equal(before.finalReview.selectedDesignSource, 'design-1');
      assert.equal(before.finalReview.primerSnapshot.forward, design1.forward);
      assert.equal(before.review.stage, 2);
      for (const n of [1, 2, 3]) assert.equal(before.evidence.selectedLanes[`case-${n}`], 'ntc');
      await page.reload(); await ready(); assert.deepEqual(stable(await saved()), stable(before));
      await open('#record-menu');
      const downloaded = page.waitForEvent('download'); await press('#export-record');
      const file = resolve(output, `${engineName}-${width}-record.json`); await (await downloaded).saveAs(file);
      const exported = JSON.parse(await readFile(file, 'utf8'));
      assert.deepEqual(stable(exported), stable(before));
      await open('#record-menu'); page.once('dialog', d => d.accept()); await press('#reset-record');
      await page.waitForFunction(key => JSON.parse(localStorage.getItem(key)).designs.length === 0, STORAGE_KEY);
      await importRecord(file, before);
      assert.deepEqual(stable(await saved()), stable(before));
      await page.reload(); await ready(); assert.deepEqual(stable(await saved()), stable(before));
      await writeFile(resolve(output, `${engineName}-${width}-state-comparison.json`), JSON.stringify({ before: stable(before), after: stable(await saved()), equal: true }, null, 2));

      // Direct hash, reload, history, scroll-driven TOC and portal return preserve state.
      await page.goto(url + '#activity-03'); await ready();
      await page.locator('#activity-03 .pcr-transition a').click(); await page.reload(); await ready();
      assert.equal(await page.locator('nav a[aria-current]').getAttribute('href'), '#activity-04');
      await page.goBack(); await page.waitForFunction(() => document.querySelector('nav a[aria-current]')?.hash === '#activity-03');
      await page.goForward(); await page.waitForFunction(() => document.querySelector('nav a[aria-current]')?.hash === '#activity-04');
      for (const id of ['activity-00', 'activity-05', 'activity-07']) {
        await page.locator(`#${id}`).evaluate(e => window.scrollTo(0, e.getBoundingClientRect().top + scrollY - 20));
        await page.waitForFunction(id => document.querySelector('nav a[aria-current]')?.hash === '#' + id, id);
      }
      await page.locator('.pcr-back').click(); assert.equal(new URL(page.url()).pathname, '/');
      await page.locator('#bioinformatics h3 a').filter({ hasText: 'PCR과 프라이머 디자인' }).click(); await ready();
      assert.deepEqual(stable(await saved()), stable(before));
      await open('#worksheet-toc'); await page.locator('nav a[href="#rna-extension"]').click();
      await page.waitForFunction(() => document.querySelector('nav a[aria-current]')?.hash === '#rna-extension');
      await page.reload(); await ready(); assert.ok(await page.locator('#rna-extension').evaluate(e => e.open));
      await checkOverflow();

      // Real keyboard events on each representative interaction, with visible focus.
      await press('[data-rna-section=strategy]'); await press('[data-rna-option=primerStrategy][data-value=same-exon]');
      await page.keyboard.press('Tab'); await page.keyboard.press('Space');
      assert.equal((await saved()).rnaExtension.primerStrategy, 'intron-spanning');
      await press('[data-rna-section=transcripts]'); await press('[data-rna-option=transcriptTarget][data-value=exon2]');
      await page.keyboard.press('Tab'); await page.keyboard.press('Space');
      assert.equal((await saved()).rnaExtension.transcriptTarget, 'exon3');
      await press('#ext-route-mine'); await page.keyboard.press('ArrowRight'); assert.equal((await saved()).externalSearch.route, 'paper');
      await page.keyboard.press('Home'); assert.equal((await saved()).externalSearch.route, 'mine');
      await press('#review-tab-length-gc'); await page.keyboard.press('ArrowRight'); assert.equal((await saved()).review.lens, 'tm');
      await press('#final-select-design-1', 'Space'); await page.keyboard.press('ArrowRight'); assert.ok(await page.locator('#final-select-design-2').isChecked());
      assert.notEqual(await page.locator('#final-select-design-2').evaluate(e => getComputedStyle(e).outlineStyle), 'none');
      await press('#evidence-lane-sample'); await page.keyboard.press('Tab'); await page.keyboard.press('Space');
      assert.equal((await saved()).evidence.selectedLanes['case-3'], 'positive');
      await press('[data-select-primer=F]');
      await page.locator('#sequence-grid button').first().focus(); await page.keyboard.press('Home'); await page.keyboard.press('Enter'); await page.keyboard.press('ArrowRight'); await page.keyboard.press('Enter');
      assert.equal(Number((await saved()).draft.bindings.F.end) - Number((await saved()).draft.bindings.F.start), 1);
      // Restoring via the real import UI returns the scientific record after keyboard checks.
      await importRecord(file, before);
      assert.deepEqual(stable(await saved()), stable(before));

      await page.emulateMedia({ reducedMotion: 'reduce' });
      await press('[data-cycle-stage=extension]');
      const animations = await page.locator('#pcr-worksheet *').evaluateAll(elements => elements.filter(e => { const s = getComputedStyle(e); return s.animationName !== 'none' && parseFloat(s.animationDuration) > 0.01; }).map(e => e.id || e.className));
      assert.deepEqual(animations, []);
      const nodeCount = await page.locator('#pcr-worksheet *').count();
      for (let n = 0; n < 12; n++) { await press('[data-cycle-stage=mixture]'); await press('[data-cycle-stage=extension]'); }
      assert.equal(await page.locator('#pcr-worksheet *').count(), nodeCount);
      assert.deepEqual(await page.evaluate(() => window.releaseUnhandled), []);
      const result = { engineName, width, roundtrip: true, reload: true, keyboard: true, consoleErrors, warnings, pageErrors, unhandledRejections, failedRequests, badResponses, nodeCount };
      results.push(result); await writeFile(resolve(output, 'release-verification.json'), JSON.stringify({ url, results, audits }, null, 2));
      assert.deepEqual(consoleErrors, []); assert.deepEqual(pageErrors, []); assert.deepEqual(unhandledRejections, []); assert.deepEqual(failedRequests, []); assert.deepEqual(badResponses, []);
      if (engineName === 'chromium' && width === 1440) representative = before;
      await context.close(); console.log(`PASS release journey ${engineName} ${width}px`);
    }
  } finally { await browser.close(); }
}

// Only after all journeys pass, capture the final build's public-facing representative views.
const browser = await chromium.launch({ headless: true });
try {
  for (const width of [1440, 768, 390]) {
    const context = await browser.newContext({ viewport: { width, height: 1100 } });
    await context.addInitScript(({ key, record }) => localStorage.setItem(key, JSON.stringify(record)), { key: STORAGE_KEY, record: representative });
    const page = await context.newPage(); await page.goto(url); await page.waitForSelector('#pcr-worksheet[data-ready=true]'); await page.evaluate(() => document.fonts.ready);
    const shot = async (name, selector, height = 1800) => {
      const path = resolve(output, name);
      if (selector) {
        await page.locator(selector).scrollIntoViewIfNeeded(); const b = await page.locator(selector).boundingBox(), offset = await page.evaluate(() => ({ x: scrollX, y: scrollY }));
        await page.screenshot({ path, fullPage: true, clip: { x: b.x + offset.x, y: b.y + offset.y, width: b.width, height: Math.min(b.height, height) }, animations: 'disabled' });
      } else await page.screenshot({ path, fullPage: name.includes('fullpage'), animations: 'disabled' });
      screenshots.push({ name, path, width });
    };
    if (width === 1440) {
      await page.evaluate(() => scrollTo(0, 0)); await shot('release-top.png');
      await shot('release-00.png', '#activity-00'); await shot('release-01-cycle.png', '#activity-01'); await shot('release-02-direction.png', '#activity-02');
      await shot('release-03-workbench.png', '#activity-03', 2400);
      await page.locator('#review-tab-complementarity').click(); await shot('release-04-review.png', '#activity-04', 2400);
      await page.locator('[data-evidence-case="case-2"]').click(); await page.locator('#evidence-lane-ntc').click(); await shot('release-05-gel.png', '#activity-05');
      await page.locator('#ext-plan').evaluate(e => e.open = true); await shot('release-06-primer-blast.png', '#activity-06'); await shot('release-06-result.png', '#ext-after-search', 2400);
      await shot('release-07-summary.png', '#activity-07'); await shot('release-07-final-notebook.png', '#final-notebook', 30000);
    } else {
      for (const n of [3, 4, 5, 6, 7]) await shot(`release-${width}-0${n}.png`, `#activity-0${n}`, 2400);
      if (width === 390) { await shot('release-mobile-core.png', '#activity-03', 2200); await shot('release-mobile-notebook.png', '#final-notebook', 30000); }
    }
    await page.locator('#rna-extension').evaluate(e => e.open = true);
    await page.locator('[data-rna-section=strategy]').click();
    await shot(width === 1440 ? 'release-rna-extension.png' : width === 390 ? 'release-mobile-rna.png' : 'release-tablet-rna.png', '#rna-extension', 3000);
    await page.locator('[data-rna-section=transcripts]').click();
    await shot(width === 1440 ? 'release-rna-transcripts.png' : `release-${width}-rna-transcripts.png`, '#rna-extension', 3000);
    if (width === 1440) await shot('course-fullpage-desktop.png');
    await context.close();
  }
  const chosen = ['release-top.png', 'release-03-workbench.png', 'release-04-review.png', 'release-05-gel.png', 'release-06-primer-blast.png', 'release-07-summary.png', 'release-rna-extension.png', 'release-rna-transcripts.png', 'release-mobile-core.png'];
  const page = await browser.newPage({ viewport: { width: 1592, height: 1100 } });
  const panels = await Promise.all(chosen.map(async name => `<figure><figcaption>${name}</figcaption><img src="data:image/png;base64,${(await readFile(resolve(output, name))).toString('base64')}"></figure>`));
  await page.setContent(`<html lang="en"><style>body{margin:24px;font:18px sans-serif;color:#2f3634}main{display:grid;grid-template-columns:repeat(3,500px);gap:22px;align-items:start}figure{margin:0;border-top:1px solid #ddd;padding-top:12px}figcaption{margin-bottom:12px}img{width:500px;height:690px;object-fit:contain;object-position:top;background:#f8f9f8}</style><main>${panels.join('')}</main></html>`);
  await page.locator('img').evaluateAll(els => Promise.all(els.map(e => e.decode())));
  await page.screenshot({ path: resolve(output, 'release-contact-sheet.png'), fullPage: true });
  await writeFile(resolve(output, 'release-verification.json'), JSON.stringify({ url, results, audits, screenshots, contactSheet: resolve(output, 'release-contact-sheet.png') }, null, 2));
} finally { await browser.close(); }
console.log('PASS final release screenshots and contact sheet');
