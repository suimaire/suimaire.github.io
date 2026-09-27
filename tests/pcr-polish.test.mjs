import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {parseRecord,STORAGE_KEY} from '../assets/js/pcr-records.mjs';
const read=p=>readFileSync(new URL('../'+p,import.meta.url),'utf8');
test('course and portal remove learning durations without changing activity anchors',()=>{
 const html=read('bioinformatics/pcr-primer-design.html'),portal=read('index.md').match(/<li class="portal-resource">\s*<h3><a[^>]+>PCR과 프라이머 디자인[\s\S]*?<\/li>/)[0];
 assert.doesNotMatch(html+portal,/(?:약\s*)?\d+분|80분 수업용/);
 for(let i=0;i<8;i++)assert.equal((html.match(new RegExp(`id="activity-0${i}"`,'g'))||[]).length,1);
 assert.match(html,/id="rna-extension"[^>]*><summary>확장 활동/);
});
test('Phase 1 through 6 fixture records preserve every answer and saved design in schema v1',()=>{
 for(let i=1;i<=6;i++){
  const raw=JSON.parse(read(`tests/fixtures/pcr/phase${i}.json`)),restored=parseRecord(JSON.stringify(raw));
  assert.equal(restored.schemaVersion,1);assert.deepEqual(restored.answers,raw.answers);assert.deepEqual(restored.designs,raw.designs);
  if(raw.finalReview)assert.deepEqual(restored.finalReview,raw.finalReview);
  assert.deepEqual(parseRecord(JSON.stringify(restored)),restored);
 }
 assert.equal(STORAGE_KEY,'hafs:pcr-primer:v1');
});
