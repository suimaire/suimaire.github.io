import assert from 'node:assert/strict';
import {readFile,writeFile,mkdir} from 'node:fs/promises';
import {resolve} from 'node:path';
import {chromium,webkit} from './.pcr-tools/node_modules/playwright/index.mjs';
import {parseRecord,STORAGE_KEY} from '../assets/js/pcr-records.mjs';
const url=process.env.PCR_TEST_URL || 'http://127.0.0.1:4173/bioinformatics/pcr-primer-design/';
const output=resolve(process.env.PCR_VERIFICATION_ROOT || 'verification.local/pcr-primer-design','phase7'); await mkdir(output,{recursive:true});
const records=await Promise.all([1,2,3,4,5,6].map(async n=>JSON.parse(await readFile(`tests/fixtures/pcr/phase${n}.json`,'utf8'))));
const titles=['두 시료를 어떻게 구별할 것인가?','PCR 한 주기에서는 무엇이 달라질까?','두 primer의 3′ 말단은 어디를 향할까?','직접 primer를 배치하기','숫자가 적절하면 좋은 primer인가?','예상과 실험 증거는 같은 것인가?','실제 데이터베이스에서는 어떻게 검토할까?','최종 설계 기록'];
let checks=0; const screenshots=[],audits=[],metrics=[];
const check=(v,m)=>{assert.ok(v,m);checks++;};
const cKeys=['cycle-selector','direction-reason','design-deletion','design-change','review-length-reason','candidate-negative'];
for(const [engineName,engine] of Object.entries({chromium,webkit})) {
 const browser=await engine.launch({headless:true});
 for(const width of [1440,1024,768,390]) {
  const context=await browser.newContext({viewport:{width,height:1100},hasTouch:width===390});
  const page=await context.newPage(),errors=[];page.on('pageerror',e=>errors.push(e.message));
  const ready=()=>page.waitForSelector('#pcr-worksheet[data-ready=true]');
  const open=async id=>{if(!await page.locator(id).evaluate(e=>e.open))await page.locator(`${id} > summary`).click();};
  const audit=async state=>{
   await page.addScriptTag({path:resolve('tests/.pcr-tools/node_modules/axe-core/axe.min.js')});
   const result=await page.evaluate(async()=>{const r=await axe.run('#pcr-worksheet',{runOnly:{type:'tag',values:['wcag2a','wcag2aa','wcag21aa','wcag22aa','best-practice']}});return {violations:r.violations.map(v=>({id:v.id,nodes:v.nodes.map(n=>n.target)})),incomplete:r.incomplete.map(v=>({id:v.id,nodes:v.nodes.map(n=>n.target)}))};});
   audits.push({engineName,width,state,...result});check(!result.violations.length,JSON.stringify(result));
   check(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth),'page overflow');
  };
  const shoot=async(name,selector,maxHeight=1450)=>{
   if(engineName!=='chromium')return;
   await page.evaluate(()=>document.fonts.ready);
   if(selector){
    await page.locator(selector).scrollIntoViewIfNeeded();const b=await page.locator(selector).boundingBox(),scroll=await page.evaluate(()=>({x:scrollX,y:scrollY}));
    await page.screenshot({path:resolve(output,name),fullPage:true,clip:{x:b.x+scroll.x,y:b.y+scroll.y,width:b.width,height:Math.min(b.height,maxHeight)},animations:'disabled'});
   }else await page.screenshot({path:resolve(output,name),fullPage:name.includes('fullpage'),animations:'disabled'});
   screenshots.push({name,path:resolve(output,name),width});
  };
  await page.goto(url);await ready();await page.evaluate(()=>document.fonts.ready);
  const visible=await page.locator('body').innerText();check(!/(?:약\s*)?\d+분/.test(visible),'no time limits');check(!/[\u00b7\u2022\u2027\u2219\u22c5\u30fb\u318d]/u.test(visible),'no decorative dots');
  check(await page.locator('h1').count()===1,'single h1');
  for(let n=0;n<8;n++) {
   const id=`activity-${String(n).padStart(2,'0')}`;
   check(await page.locator(`#${id} h2`).first().innerText()===titles[n],`${id} heading`);
   check((await page.locator(`nav a[href="#${id}"]`).textContent()).includes(titles[n]),`${id} TOC title`);
  }
  for(const id of cKeys)check(await page.locator(`#${id}`).isHidden(),`${id} no redundant blank prompt`);
  for(const id of ['candidate-prediction','review-ab-reason'])check(await page.locator(`#${id}`).getAttribute('rows')==='1','short answer size');
  for(const id of ['final-review','final-experimental','final-external'])check(!await page.locator(`#${id} > details`).evaluate(e=>e.open),'07 details closed');
  check(await page.locator('#final-external-summary').isVisible(),'external status exposed');
  await audit('empty');
  // Each transition is a real hash link; reload and history restore the exact heading.
  for(let n=2;n<7;n++) {
   await page.locator(`#activity-0${n} .pcr-transition a`).click();
   check(new URL(page.url()).hash===`#activity-0${n+1}`,'transition hash');
   await page.waitForFunction(id=>document.querySelector('nav a[aria-current]')?.hash===id,`#activity-0${n+1}`);
   check(await page.locator(`#activity-0${n+1} h2`).first().evaluate(e=>e.getBoundingClientRect().top>=0),'anchor heading visible');
  }
  await page.reload();await ready();check(await page.locator('nav a[aria-current]').getAttribute('href')==='#activity-07','hash reload current');
  await page.goBack();await page.waitForFunction(()=>document.querySelector('nav a[aria-current]')?.hash==='#activity-06');check(true,'history anchor');
  for(const n of [1,5,3,0,7]) {
   await page.evaluate(n=>window.scrollTo(0,document.querySelector(`#activity-0${n}`).getBoundingClientRect().top+scrollY-20),n);
   await page.waitForFunction(n=>document.querySelector('nav a[aria-current]')?.hash===`#activity-0${n}`,n);
   check(await page.locator('nav a[aria-current]').count()===1,'one scroll current');
  }
  if(width===1440)for(const [index,record] of records.entries()) {
   // Load both historical browser storage and exported JSON through the real UI.
   await page.evaluate(({record,key})=>localStorage.setItem(key,JSON.stringify(record)),{record,key:STORAGE_KEY});await page.reload();await ready();
   check(!(await page.locator('#save-status').innerText()).includes('읽을 수 없습니다'),`phase ${index+1} storage accepted`);
   for(const [key,value] of Object.entries(record.answers)) {
    const field=page.locator(`[data-answer][id="${key}"]`);
    if(await field.count())check(await field.inputValue()===value,`phase ${index+1} answer ${key}`);
   }
   await open('#record-menu');page.once('dialog',d=>d.accept());
   await page.locator('#import-record').setInputFiles({name:`phase${index+1}.json`,mimeType:'application/json',buffer:Buffer.from(JSON.stringify(record))});
   await page.waitForFunction(()=>document.querySelector('#save-status').textContent.includes('자동 저장'));
   const saved=await page.evaluate(key=>JSON.parse(localStorage.getItem(key)),STORAGE_KEY),expected=parseRecord(JSON.stringify(record));
   for(const key of ['answers','draft','designs','initialPrimerPrediction','review','evidence','externalSearch','finalReview'])check(JSON.stringify(saved[key])===JSON.stringify(expected[key]),`phase ${index+1} JSON ${key}`);
   if(index===0)for(const id of cKeys.filter(k=>record.answers[k])) {
    if(id.startsWith('design-'))await open('#design-reflection');
    if(id==='candidate-negative'){await page.locator('#review-tab-off-target').click();await page.locator('#compare-ab').click();await page.locator('#compare-abc').click();}
    if(id==='review-length-reason')await page.locator('#review-tab-length-gc').click();
    const parent=page.locator(`#${id}`).locator('..');check(await parent.isVisible(),`${id} legacy wrapper exposed`);await parent.locator('summary').click();check(await page.locator(`#${id}`).isVisible(),'old answer readable');
   }
  }
  // A populated historical Phase 6 record provides the same representative state at every width.
  const representative=structuredClone(records[5]);
  for(const key of cKeys)delete representative.answers[key];
  await page.evaluate(({record,key})=>localStorage.setItem(key,JSON.stringify(record)),{record:representative,key:STORAGE_KEY});await page.goto(url);await ready();
  check(await page.locator('#final-question').inputValue()===representative.finalReview.researchQuestion,'essential 07 reflection retained');
  await page.locator('#review-tab-off-target').click();await page.locator('#compare-ab').click();await page.locator('#compare-abc').click();await page.locator('#review-tab-length-gc').click();
  await page.locator('#evidence-lane-sample').click();
  await audit('populated');
  if(width===1440) {
   await page.evaluate(()=>scrollTo(0,0));await shoot('course-top-desktop.png');await shoot('course-toc-desktop.png','.pcr-sidebar nav');
   for(let n=0;n<8;n++)await shoot(`${String(n).padStart(2,'0')}-polished.png`,`#activity-${String(n).padStart(2,'0')}`);
   await shoot('course-fullpage-desktop.png');
   const inventory=await page.locator('textarea,input:not([type=file]),select').evaluateAll(els=>els.map(e=>({id:e.id,name:e.name,type:e.type,rows:e.rows||null,label:[...(e.labels||[])].map(l=>l.textContent.trim()).join(' / '),answer:e.hasAttribute('data-answer'),external:e.dataset.external,final:e.dataset.final,legacy:Boolean(e.closest('[data-legacy-answers],#ext-legacy,#final-legacy,#evidence-legacy,#legacy-first-negative'))})));
   await writeFile(resolve(output,'input-inventory.json'),JSON.stringify(inventory,null,2));
  }
  if(width===390) {
   await page.evaluate(()=>scrollTo(0,0));await shoot('course-mobile-top.png');
   for(const n of [3,5,6,7])await shoot(`0${n}-mobile-polished.png`,`#activity-0${n}`);
   await shoot('03-mobile-current-design.png','#current-design',30000);
  }
  if(width===768)for(const n of [3,4,7])await shoot(`0${n}-tablet-polished.png`,`#activity-0${n}`);
  metrics.push({engineName,width,editorHeight:await page.locator('#final-editor').evaluate(e=>e.getBoundingClientRect().height),pageHeight:await page.evaluate(()=>document.documentElement.scrollHeight)});
  // Read-only report has ten sections, complete reflections, and one detailed final F/R.
  await page.locator('#final-notebook-toggle').focus();await page.keyboard.press('Enter');
  check(await page.locator('#final-notebook').isVisible(),'keyboard opens notebook');check(await page.locator('#final-notebook input,#final-notebook textarea').count()===0,'report read only');
  check(await page.locator('#final-notebook > section > h4').count()===10,'ten notebook sections');
  for(const value of [representative.finalReview.researchQuestion,representative.finalReview.revisionReflection,representative.answers['candidate-judgment'],representative.answers['evidence-identity']])check((await page.locator('#final-notebook').innerText()).includes(value),'complete report record');
  for(const sequence of [representative.draft.forward,representative.draft.reverse])check((await page.locator('#final-notebook').innerText()).split(sequence).length===2,'final sequence shown once');
  await audit('notebook');if(width===1440)await shoot('07-final-notebook.png','#final-notebook',30000);
  if(width===390)await shoot('07-final-notebook-mobile.png','#final-notebook',30000);
  if(width===1440) {
   const matched=structuredClone(representative);matched.externalSearch=structuredClone(records[4].externalSearch);
   Object.assign(matched.externalSearch.candidates[0],matched.draft);matched.finalReview.notebookExpanded=true;
   await open('#record-menu');page.once('dialog',d=>d.accept());
   await page.locator('#import-record').setInputFiles({name:'same-pair.json',mimeType:'application/json',buffer:Buffer.from(JSON.stringify(matched))});
   await page.waitForFunction(()=>document.querySelector('#final-external-summary').textContent.includes('최종 pair와 같습니다'));
   for(const seq of [matched.draft.forward,matched.draft.reverse])check((await page.locator('#final-notebook').innerText()).split(seq).length===2,'matched external pair references final sequence once');
   check((await page.locator('#final-notebook').innerText()).includes(matched.externalSearch.conditions.database),'external conditions kept');
   check((await page.locator('#final-external-summary').innerText()).includes('기록된 조건'),'claim exposed while details collapsed');
   matched.externalSearch.candidates[0].forward='ACGTACGTACGTACGTACGT';
   page.once('dialog',d=>d.accept());await page.locator('#import-record').setInputFiles({name:'different-pair.json',mimeType:'application/json',buffer:Buffer.from(JSON.stringify(matched))});
   await page.waitForFunction(()=>document.querySelector('#final-notebook').textContent.includes('ACGTACGTACGTACGTACGT'));
   check((await page.locator('#final-notebook').innerText()).includes('검색 결과로 간주하지 마세요'),'different external pair retained and distinguished');
   // Print folding is presentation only; imported old prose remains printable.
   page.once('dialog',d=>d.accept());await page.locator('#import-record').setInputFiles({name:'old-answers.json',mimeType:'application/json',buffer:Buffer.from(JSON.stringify(records[0]))});
   await page.waitForFunction(()=>!document.querySelector('#cycle-selector').closest('details').hidden);
   await page.evaluate(()=>{window.print=()=>window.dispatchEvent(new Event('beforeprint'));});
   await open('#print-menu');await page.locator('#print-filled').click();await page.emulateMedia({media:'print'});
   for(const id of cKeys.filter(k=>records[0].answers[k])){check(await page.locator(`#${id} + .pcr-print-value`).isVisible(),'old answer included in written print');check(await page.locator(`#${id} + .pcr-print-value`).innerText()===records[0].answers[id],'old print text intact');}
   await page.emulateMedia({media:'screen'});await page.evaluate(()=>window.dispatchEvent(new Event('afterprint')));
   check(!await page.locator('#cycle-selector').locator('..').evaluate(e=>e.open),'print restores folded state');
  }
  await page.emulateMedia({reducedMotion:'reduce'});check(await page.evaluate(()=>matchMedia('(prefers-reduced-motion: reduce)').matches),'reduced motion');
  check(!errors.length,JSON.stringify(errors));await context.close();console.log(`PASS Phase 7 ${engineName} ${width}px`);
 }
 await browser.close();
}
const browser=await chromium.launch({headless:true});
for(const [name,items,cellWidth] of [
 ['phase7-contact-sheet.png',screenshots.filter(s=>/^\d\d-polished/.test(s.name)),600],
 ['phase7-contact-sheet-mobile.png',screenshots.filter(s=>/^(course-mobile-top|\d\d-mobile-polished)/.test(s.name)),390]
]) {
 const page=await browser.newPage({viewport:{width:cellWidth*2+72,height:1100}});
 const panels=await Promise.all(items.map(async s=>`<figure><figcaption>${s.name}</figcaption><img src="data:image/png;base64,${(await readFile(s.path)).toString('base64')}"></figure>`));
 await page.setContent(`<html lang="en"><style>body{margin:24px;font:18px sans-serif;color:#2f3634}main{display:grid;grid-template-columns:repeat(2,${cellWidth}px);gap:24px}figure{margin:0;border-top:1px solid #ddd;padding-top:12px}figcaption{margin-bottom:12px}img{width:100%;height:auto}</style><main>${panels.join('')}</main></html>`);
 await page.locator('img').evaluateAll(els=>Promise.all(els.map(e=>e.decode())));await page.screenshot({path:resolve(output,name),fullPage:true});await page.close();
}
await browser.close();
await writeFile(resolve(output,'verification.json'),JSON.stringify({url,checks,audits,metrics,screenshots},null,2));
console.log(`PASS Phase 7 ${checks} assertions, ${audits.length} full-page accessibility audits`);
