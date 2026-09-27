import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { emptyRecord, parseRecord, STORAGE_KEY } from '../assets/js/pcr-records.mjs';
import { emptyRnaExtension, validateRnaExtension, RNA_OPTIONS, RNA_TRANSCRIPTS, strategyModel, transcriptBinding } from '../assets/js/pcr-rna.mjs';
const read = p => readFileSync(new URL('../' + p, import.meta.url), 'utf8');
test('RNA optional v1 state accepts unused and partial activity without changing core records', () => {
  for (let n = 1; n <= 6; n++) {
    const raw = JSON.parse(read(`tests/fixtures/pcr/phase${n}.json`));
    const parsed = parseRecord(JSON.stringify(raw));
    assert.deepEqual(parsed.answers, raw.answers); assert.deepEqual(parsed.designs, raw.designs);
    assert.equal(parsed.rnaExtension.activeSection, 'flow');
    if (raw.answers['rna-plan'] !== undefined) assert.equal(parsed.rnaExtension.answers.reflection, raw.answers['rna-plan']);
  }
  const raw = emptyRecord(); raw.rnaExtension = { activeSection: 'strategy', answers: { strategy: 'intron 길이도 검토' } };
  assert.equal(parseRecord(JSON.stringify(raw)).rnaExtension.primerStrategy, 'same-exon');
  assert.equal(raw.schemaVersion, 1); assert.equal(STORAGE_KEY, 'hafs:pcr-primer:v1');
});
test('RNA legacy migration retains verbatim original and never overwrites edited or intentionally empty reflection', () => {
  const raw = emptyRecord(); delete raw.rnaExtension;
  raw.answers['rna-plan'] = '이전 RNA 답안\n<특수 문자와 개행>';
  let parsed = parseRecord(JSON.stringify(raw));
  assert.equal(parsed.rnaExtension.answers.reflection, raw.answers['rna-plan']);
  parsed.rnaExtension.answers.reflection = '';
  parsed = parseRecord(JSON.stringify(parsed));
  assert.equal(parsed.rnaExtension.answers.reflection, ''); assert.deepEqual(parsed.answers, raw.answers);
});
test('RNA completed state and all selection combinations round trip independently', () => {
  const raw = emptyRecord(), core = structuredClone(raw); delete core.rnaExtension;
  raw.rnaExtension.answers = { template: 'cdna', strategy: 'intron', control: 'DNA 가능성', transcripts: '범위', specific: 'junction', reflection: '연결\n'.repeat(500) };
  for (const [key, options] of Object.entries(RNA_OPTIONS)) for (const option of options) {
    raw.rnaExtension[key] = option;
    const parsed = parseRecord(JSON.stringify(raw)); assert.deepEqual(parsed, raw);
    delete parsed.rnaExtension; assert.deepEqual(parsed, core);
  }
});
test('malformed RNA imports reject invalid structures choices answer keys and oversized text', () => {
  for (const value of [null, [], 'rna', { flowStep: 'dna' }, { transcriptTarget: 'exon5' }, { answers: [] }, { answers: { template: 'dna' } }, { answers: { reflection: 4 } }, { answers: { reflection: 'a'.repeat(12001) } }, { answers: { injected: 'text' } }, JSON.parse('{"__proto__":{}}')]) {
    assert.throws(() => validateRnaExtension(value));
  }
});
test('same-exon intron-spanning and junction models explain conditional template binding', () => {
  assert.match(strategyModel('same-exon').genomic, /동일한 크기/);
  for (const length of ['short', 'long']) {
    const m = strategyModel('intron-spanning', length);
    assert.match(m.cdna, /짧은 product/); assert.match(m.genomic, /intron.*product가 가능/);
    assert.match(m.meaning, /extension time.*polymerase.*PCR 조건/);
    assert.doesNotMatch(m.genomic, /불가능|증폭되지/);
  }
  assert.match(strategyModel('junction', 'long', 'exon').meaning, /아직 Exon 1 내부/);
  assert.match(strategyModel('junction').cdna, /연속된 F 결합 서열/);
  assert.match(strategyModel('junction').genomic, /동일한 연속 F 결합 부위가 없/);
  assert.match(strategyModel('junction').meaning, /보장하지/);
});
test('transcript targets depend on exon presence and adjacency, with shared R on exon 4', () => {
  const expected = { exon2: [true, true, false], exon3: [true, false, true], exon4: [true, true, true], junction23: [true, false, false], junction34: [true, false, true] };
  for (const [target, result] of Object.entries(expected)) assert.deepEqual(RNA_TRANSCRIPTS.map(exons => transcriptBinding(exons, target).present), result);
  assert.equal(transcriptBinding([1, 2, 4, 3], 'junction23').present, false);
});
test('RNA module keeps qPCR optional and avoids unsafe DOM or runtime network requests', () => {
  const html = read('bioinformatics/pcr-primer-design.html').split('<details id="rna-extension"')[1].split('<footer')[0];
  assert.match(html, /RT-PCR 자체가 정량 PCR을 뜻하지는/);
  assert.match(html, /더 알아보기 \/ qPCR/);
  assert.doesNotMatch(html, /ΔΔCt|ΔCt|standard curve|probe chemistry|validated|PASS|FAIL|best design/);
  const sources = ['assets/js/pcr-rna.mjs', 'assets/js/pcr-rna-view.mjs'].map(read).join('');
  assert.doesNotMatch(sources, /innerHTML|insertAdjacentHTML|eval\(|fetch\(|XMLHttpRequest|[\u00b7\u2022\u2027\u2219\u22c5\u30fb\u318d]/u);
  assert.match(html, /id="rna-plan" data-answer rows="4" readonly/);
  assert.deepEqual(emptyRnaExtension().answers, {});
});
