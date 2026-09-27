import assert from 'node:assert/strict';
import { mkdir, readFile, writeFile } from 'node:fs/promises';
import { resolve } from 'node:path';
import { chromium, webkit } from './.pcr-tools/node_modules/playwright/index.mjs';
import { parseRecord, STORAGE_KEY } from '../assets/js/pcr-records.mjs';
const url = process.env.PCR_TEST_URL || 'http://127.0.0.1:4173/bioinformatics/pcr-primer-design/';
const output = resolve('verification.local/pcr-primer-design/phase8');
await mkdir(output, { recursive: true });
const histories = await Promise.all([1, 2, 3, 4, 5, 6].map(async n => JSON.parse(await readFile(`tests/fixtures/pcr/phase${n}.json`, 'utf8'))));
// Phase 7 preserved Phase 6's record format. Include a Phase 7 style legacy RNA answer.
const phase7 = structuredClone(histories[5]); phase7.answers['rna-plan'] = 'RNA에서 cDNA로\n이전 RT(-) 계획 <원문>';
histories.push(phase7);
let checks = 0; const screenshots = [], audits = [];
const check = (v, m) => { assert.ok(v, m); checks++; };
const core = record => { const result = structuredClone(record); delete result.rnaExtension; delete result.updatedAt; return result; };
for (const [engineName, engine] of Object.entries({ chromium, webkit })) {
  const browser = await engine.launch({ headless: true });
  for (const width of [1440, 768, 390]) {
    const context = await browser.newContext({ viewport: { width, height: 1100 }, hasTouch: width === 390 });
    const page = await context.newPage(), errors = []; page.on('pageerror', e => errors.push(e.message));
    const ready = () => page.waitForSelector('#pcr-worksheet[data-ready=true]');
    const open = async id => { if (!await page.locator(id).evaluate(e => e.open)) await page.locator(`${id} > summary`).click(); };
    const section = async name => { await open('#rna-extension'); await page.locator(`[data-rna-section=${name}]`).click(); };
    const option = (key, value) => page.locator(`[data-rna-option=${key}][data-value="${value}"]`);
    const saved = () => page.evaluate(key => JSON.parse(localStorage.getItem(key)), STORAGE_KEY);
    const shot = async (name, selector, maxHeight) => {
      if (engineName !== 'chromium' || width === 768) return;
      const filename = width === 390 ? name.replace('.png', '-mobile.png') : name;
      await page.evaluate(() => document.fonts.ready);
      const target = page.locator(selector); await target.scrollIntoViewIfNeeded();
      const box = await target.boundingBox(), scroll = await page.evaluate(() => ({ x: scrollX, y: scrollY }));
      await page.screenshot({ path: resolve(output, filename), fullPage: true, clip: { x: box.x + scroll.x, y: box.y + scroll.y, width: box.width, height: Math.min(box.height, maxHeight || box.height) }, animations: 'disabled' });
      screenshots.push({ name: filename, width, path: resolve(output, filename) });
    };
    const audit = async name => {
      check(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth), `${name} no page overflow`);
      const result = await page.evaluate(async () => {
        const r = await axe.run('#rna-extension', { runOnly: { type: 'tag', values: ['wcag2a', 'wcag2aa', 'wcag21aa', 'wcag22aa', 'best-practice'] } });
        return { violations: r.violations.map(v => ({ id: v.id, nodes: v.nodes.map(n => n.target) })), incomplete: r.incomplete.map(v => ({ id: v.id, nodes: v.nodes.map(n => n.target) })) };
      });
      audits.push({ engineName, width, name, ...result }); check(!result.violations.length, JSON.stringify(result));
    };
    await page.goto(url); await ready();
    if (width === 1440) {
      await open('#record-menu'); const download = page.waitForEvent('download'); await page.locator('#export-record').click();
      const path = resolve(output, `${engineName}-unused.json`); await (await download).saveAs(path);
      const unused = parseRecord(await readFile(path, 'utf8')); check(Object.keys(unused.rnaExtension.answers).length === 0, 'unused RNA record exports');
    }
    await page.addScriptTag({ path: resolve('tests/.pcr-tools/node_modules/axe-core/axe.min.js') });
    // Exercise the 06 link and full navigation, not a programmatically opened screenshot.
    await open('#ext-plan'); await page.locator('#ext-rna-help a').click(); check(await page.locator('#rna-extension').evaluate(e => e.open), '06 opens extension');
    await shot('RNA-overview.png', '#rna-extension', 1150);
    await section('flow');
    for (const step of ['rna', 'rt', 'cdna', 'pcr']) {
      await option('flowStep', step).focus(); await page.keyboard.press('Enter');
      check(await option('flowStep', step).getAttribute('aria-pressed') === 'true', `flow ${step} selected`);
      check(await option('flowStep', step).evaluate(e => getComputedStyle(e).outlineStyle !== 'none'), 'visible keyboard focus');
      check((await saved()).rnaExtension.flowStep === step, 'step saved');
    }
    await page.locator('#rna-template').selectOption('rna'); check((await page.locator('#rna-template-feedback').innerText()).includes('출발 물질'), 'template corrective explanation');
    await page.locator('#rna-template').selectOption('cdna');
    check((await saved()).rnaExtension.answers.template === 'cdna', 'cDNA answer saved');
    check((await page.locator('#rna-template-feedback').innerText()).includes('cDNA가 PCR template'), 'template concept');
    await option('flowStep', 'cdna').click(); await audit('flow'); await shot('RNA-cdna-flow.png', '#rna-flow');
    if (width === 1440) {
      await open('#record-menu'); const download = page.waitForEvent('download'); await page.locator('#export-record').click();
      const path = resolve(output, `${engineName}-partial.json`); await (await download).saveAs(path);
      page.once('dialog', d => d.accept()); await page.locator('#import-record').setInputFiles(path);
      await page.waitForFunction(() => document.querySelector('#save-status').textContent.includes('자동 저장'));
      check((await saved()).rnaExtension.answers.template === 'cdna', 'partial RNA JSON re-import');
    }
    await section('genomic');
    const structureText = await page.locator('#rna-structure-comparison').innerText();
    check(structureText.includes('Genomic DNA') && structureText.includes('Mature mRNA / cDNA'), 'both structures');
    check(await page.locator('#rna-structure-comparison svg[aria-label]').count() === 2, 'structure text alternatives');
    if (width === 390) {
      const scroller = page.locator('#rna-structure-comparison .rna-diagram-scroll').first();
      check(await scroller.evaluate(e => e.scrollWidth > e.clientWidth), 'internal diagram scroll without page scaling');
      await scroller.focus(); await page.keyboard.press('ArrowRight');
      await page.waitForFunction(() => document.querySelector('#rna-structure-comparison .rna-diagram-scroll').scrollLeft > 0);
      check(true, 'keyboard scroll');
    }
    await audit('genomic'); await shot('RNA-genomic-vs-cdna.png', '#rna-genomic');
    await section('strategy'); const before = core(await saved());
    await option('primerStrategy', 'same-exon').click();
    check((await page.locator('#rna-strategy-comparison').innerText()).includes('동일한 크기'), 'same exon binding on both');
    await shot('RNA-strategy-same-exon.png', '#rna-strategy');
    await option('primerStrategy', 'intron-spanning').focus(); await page.keyboard.press('Space');
    for (const intron of ['short', 'long']) {
      await option('intronExample', intron).click();
      const text = await page.locator('#rna-strategy-comparison').innerText();
      check(text.includes('짧은 product') && text.includes('intron') && text.includes('product가 가능'), 'conditional intron products');
      check(!/불가능|증폭되지/.test(text), 'no absolute exclusion');
    }
    await shot('RNA-strategy-intron-spanning.png', '#rna-strategy');
    await option('primerStrategy', 'junction').click();
    await option('junctionPosition', 'exon').click(); check((await page.locator('#rna-placement').innerText()).includes('F는 Exon 1'), 'F inside exon');
    await option('junctionPosition', 'junction').focus(); await page.keyboard.press('Enter');
    check((await page.locator('#rna-strategy-comparison').innerText()).includes('연속된 F 결합 서열'), 'junction continuous cDNA site');
    check((await page.locator('#rna-strategy-comparison').innerText()).includes('동일한 연속 F 결합 부위가 없'), 'gDNA separated site');
    check(await page.locator('#rna-strategy-comparison .rna-separated').count() === 1, 'gDNA split binding visual');
    await page.locator('#rna-strategy-answer').fill('긴 intron의 길이 차이와 spliced junction의 연속 결합 부위를 활용한다.');
    await audit('strategy'); await shot('RNA-strategy-junction.png', '#rna-strategy');
    await section('control'); await option('rtControlCase', 'case1').click();
    check(await page.locator('#rna-minus-signal').innerText() === 'band 없음', 'case1 RT minus'); await shot('RNA-rt-control.png', '#rna-control');
    await option('rtControlCase', 'case2').focus(); await page.keyboard.press('Space');
    check(await page.locator('#rna-minus-signal').innerText() === 'band 있음', 'case2 RT minus');
    check((await page.locator('#rna-control-meaning').innerText()).includes('확정은 아니며'), 'not contamination diagnosis');
    await page.locator('#rna-control-answer').fill('RNA-derived cDNA가 아닌 DNA template나 다른 실험 문제를 함께 검토한다.');
    await audit('control'); await shot('RNA-rt-minus-signal.png', '#rna-control');
    await section('transcripts');
    const expected = { exon2: [true, true, false], exon3: [true, false, true], exon4: [true, true, true], junction23: [true, false, false], junction34: [true, false, true] };
    for (const [target, values] of Object.entries(expected)) {
      await option('transcriptTarget', target).focus(); await page.keyboard.press('Enter');
      check(await option('transcriptTarget', target).getAttribute('aria-pressed') === 'true', 'transcript selected aria');
      const texts = await page.locator('.rna-transcript-row p').allTextContents();
      texts.forEach((text, i) => check(text.includes(values[i] ? '결합 가능한 구조' : '결합 구조가 없'), `${target} transcript ${i + 1}`));
      if (target === 'exon4') { check((await page.locator('#rna-transcript-meaning').innerText()).includes('합쳐질 수'), 'pooled signal caveat'); await shot('RNA-transcript-variants.png', '#rna-transcripts'); }
    }
    await option('transcriptTarget', 'junction23').click(); await shot('RNA-transcript-selection.png', '#rna-transcripts');
    await page.locator('#rna-transcripts-answer').fill('결합 가능한 transcript subset에 따라 관찰 신호의 범위가 달라진다.');
    await page.locator('#rna-specific-answer').fill('목표 transcript에만 있는 exon 또는 연속 junction을 찾아 전체 transcriptome과 비교한다.');
    await audit('transcripts');
    await section('reflection');
    const reflection = 'RNA → cDNA와 exon 연결을 확인한다.\nRT(-)로 DNA 기여를 검토하고 transcript 범위를 기록한다. <img src=x onerror="window.rnaInjected=true">';
    await page.locator('#rna-reflection-answer').fill(reflection.split(' <img')[0]);
    check(JSON.stringify(core(await saved())) === JSON.stringify(before), 'RNA does not change Core state');
    await audit('reflection'); await shot('RNA-final-reflection.png', '#rna-reflection');
    await page.locator('#rna-reflection-answer').fill(reflection);
    await page.reload(); await ready(); await open('#rna-extension');
    check(await page.locator('#rna-reflection').isVisible(), 'active section reload');
    check(await page.locator('#rna-reflection-answer').inputValue() === reflection, 'reflection reload');
    const complete = await saved();
    check(complete.rnaExtension.primerStrategy === 'junction' && complete.rnaExtension.rtControlCase === 'case2' && complete.rnaExtension.transcriptTarget === 'junction23', 'all choices reload');
    await open('#record-menu'); const downloaded = page.waitForEvent('download'); await page.locator('#export-record').click();
    const filename = resolve(output, `${engineName}-${width}-roundtrip.json`); await (await downloaded).saveAs(filename);
    const exported = JSON.parse(await readFile(filename, 'utf8')); check(JSON.stringify(exported.rnaExtension) === JSON.stringify(complete.rnaExtension), 'new export complete');
    await page.locator('#rna-reflection-answer').fill('임시'); await open('#record-menu'); page.once('dialog', d => d.accept()); await page.locator('#import-record').setInputFiles(filename);
    await page.waitForFunction(() => document.querySelector('#rna-reflection-answer').value.startsWith('RNA →'));
    check(!await page.evaluate(() => window.rnaInjected), 'RNA imported text inert');
    await open('#print-menu'); await page.evaluate(() => { window.print = () => window.dispatchEvent(new Event('beforeprint')); });
    await page.locator('#print-filled').click(); await page.emulateMedia({ media: 'print' });
    for (const id of ['flow', 'genomic', 'strategy', 'control', 'transcripts', 'reflection']) check(await page.locator(`#rna-${id}`).isVisible(), 'all RNA sections printable');
    check(await page.locator('#rna-reflection-answer + .pcr-print-value').innerText() === reflection, 'written print answer');
    await page.emulateMedia({ media: 'screen' }); await page.evaluate(() => window.dispatchEvent(new Event('afterprint')));
    check(await page.locator('[data-rna-panel]:visible').count() === 1, 'print restores selected panel');
    await page.locator('#print-blank').click(); await page.emulateMedia({ media: 'print' });
    for (const id of ['rna-template', 'rna-reflection-answer', 'rna-control-answer']) check(await page.locator(`#${id} + .pcr-print-value`).innerText() === '', 'blank print excludes answer');
    await page.emulateMedia({ media: 'screen' }); await page.evaluate(() => window.dispatchEvent(new Event('afterprint')));
    check(JSON.stringify((await saved()).rnaExtension) === JSON.stringify(complete.rnaExtension), 'printing never alters saved RNA');
    if (width === 1440) for (const [index, raw] of histories.entries()) {
      await page.evaluate(({ key, raw }) => localStorage.setItem(key, JSON.stringify(raw)), { key: STORAGE_KEY, raw }); await page.reload(); await ready();
      check(await page.locator('#rna-reflection-answer').inputValue() === (raw.answers['rna-plan'] || ''), `phase${index + 1} storage migration`);
      await open('#record-menu'); page.once('dialog', d => d.accept());
      await page.locator('#import-record').setInputFiles({ name: `phase${index + 1}.json`, mimeType: 'application/json', buffer: Buffer.from(JSON.stringify(raw)) });
      await page.waitForFunction(() => document.querySelector('#save-status').textContent.includes('자동 저장'));
      check(JSON.stringify(core(await saved())) === JSON.stringify(core(parseRecord(JSON.stringify(raw)))), `phase${index + 1} core unchanged`);
      if (raw.answers['rna-plan']) {
        await section('reflection'); await open('#rna-legacy');
        check(await page.locator('#rna-plan').inputValue() === raw.answers['rna-plan'], 'legacy verbatim visible');
        await page.locator('#rna-reflection-answer').fill('새 reflection');
        check((await saved()).answers['rna-plan'] === raw.answers['rna-plan'], 'editing never overwrites original');
        if (index === 6) {
          await page.locator('#rna-legacy > summary').click(); await page.locator('#rna-extension > summary').click();
          await page.evaluate(() => { window.print = () => window.dispatchEvent(new Event('beforeprint')); });
          await open('#print-menu'); await page.locator('#print-filled').click(); await page.emulateMedia({ media: 'print' });
          check(await page.locator('#rna-plan + .pcr-print-value').innerText() === raw.answers['rna-plan'], 'legacy original prints alongside new reflection');
          await page.emulateMedia({ media: 'screen' }); await page.evaluate(() => window.dispatchEvent(new Event('afterprint')));
          check(!await page.locator('#rna-extension').evaluate(e => e.open) && !await page.locator('#rna-legacy').evaluate(e => e.open), 'print restores closed legacy and extension');
        }
      }
    }
    // Reject malformed optional state without losing the current record.
    const previous = await saved(); const invalid = structuredClone(previous); invalid.rnaExtension.transcriptTarget = 'exon99';
    await open('#record-menu'); await page.locator('#import-record').setInputFiles({ name: 'invalid-rna.json', mimeType: 'application/json', buffer: Buffer.from(JSON.stringify(invalid)) });
    await page.waitForFunction(() => document.querySelector('#save-status').textContent.includes('불러오기 실패'));
    check(JSON.stringify(await saved()) === JSON.stringify(previous), 'invalid RNA import keeps current record');
    await page.emulateMedia({ reducedMotion: 'reduce' }); await section('strategy');
    check(await page.locator('#rna-strategy svg').first().evaluate(e => getComputedStyle(e).animationName) === 'none', 'reduced motion');
    const buttons = await page.locator('#rna-extension button:visible').evaluateAll(els => els.map(e => e.getBoundingClientRect().height));
    check(buttons.every(height => height >= 44), 'touch targets at least 44px');
    if (width === 390) { await option('primerStrategy', 'junction').click(); await option('junctionPosition', 'junction').click(); await page.locator('#rna-strategy').scrollIntoViewIfNeeded(); if (engineName === 'chromium') await page.screenshot({ path: resolve(output, 'RNA-mobile.png'), animations: 'disabled' }); }
    check(!errors.length, JSON.stringify(errors));
    console.log(`PASS Phase 8 ${engineName} ${width}px`); await context.close();
  }
  await browser.close();
}
const browser = await chromium.launch({ headless: true });
for (const width of [1440, 390]) {
  const chosen = screenshots.filter(s => s.width === width && /RNA-(cdna-flow|strategy-|rt-minus|transcript-)/.test(s.name));
  const cell = width === 390 ? 390 : 760;
  const page = await browser.newPage({ viewport: { width: cell * 2 + 72, height: 1100 } });
  const panels = await Promise.all(chosen.map(async s => `<figure><figcaption>${s.name}</figcaption><img src="data:image/png;base64,${(await readFile(s.path)).toString('base64')}"></figure>`));
  await page.setContent(`<html lang="en"><style>body{margin:24px;font:18px sans-serif;color:#2f3634}main{display:grid;grid-template-columns:repeat(2,${cell}px);gap:24px;align-items:start}figure{margin:0;border-top:1px solid #ddd;padding-top:12px}figcaption{margin-bottom:12px}img{width:100%;height:auto}</style><main>${panels.join('')}</main></html>`);
  await page.locator('img').evaluateAll(els => Promise.all(els.map(e => e.decode())));
  await page.screenshot({ path: resolve(output, `phase8-contact-sheet${width === 390 ? '-mobile' : ''}.png`), fullPage: true }); await page.close();
}
await browser.close();
await writeFile(resolve(output, 'verification.json'), JSON.stringify({ url, checks, audits, screenshots }, null, 2));
console.log(`PASS Phase 8 ${checks} assertions, ${audits.length} accessibility audits`);
