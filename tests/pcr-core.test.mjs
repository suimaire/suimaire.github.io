import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import * as core from '../assets/js/pcr-core.mjs';
import { emptyRecord, parseRecord, validateRecord } from '../assets/js/pcr-records.mjs';
const fixture = JSON.parse(readFileSync(new URL('../assets/data/pcr-primer-fixture.json', import.meta.url)));
const { templates, candidatePairs, expectedExactMatchProducts } = fixture;
test('reverse complement and involution', () => {
  assert.equal(core.reverseComplement('AGTCCGTA'), 'TACGGACT');
  for (const s of ['ACGT', 'AGTCCGTA', templates.A.sequence]) assert.equal(core.reverseComplement(core.reverseComplement(s)), s);
  assert.equal(core.complement('AGTCCGTA'), 'TCAGGCAT');
});
test('fixture lengths and exact deletion', () => {
  assert.deepEqual(Object.values(templates).map(t => t.sequence.length), [420, 340, 420]);
  assert.equal(templates.B.sequence, templates.A.sequence.slice(0, 120) + templates.A.sequence.slice(200));
});
test('P1 length and GC', () => {
  for (const p of Object.values(candidatePairs.P1)) {
    assert.equal(core.primerStats(p).length, 20); assert.equal(core.primerStats(p).gcPercent, 50);
  }
});
for (const [id, pair] of Object.entries(candidatePairs)) for (const [source, template] of Object.entries(templates)) {
  test(`${id}/${source} exact coordinates and inclusive lengths`, () => {
    const actual = core.analyzeTemplate(template.sequence, pair.forward, pair.reverse, source);
    assert.deepEqual(actual.products.map(p => [p.start + 1, p.end, p.length]), expectedExactMatchProducts[id][source]);
    for (const p of actual.products) assert.equal(p.sequence.length, p.length);
  });
}
test('wrong R order sequence does not produce intended product', () => {
  assert.equal(core.analyzeTemplate(templates.A.sequence, candidatePairs.P1.forward, 'TCGTACGATCGTAGCCTGAA').products.length, 0);
});
test('repeated binding sites including overlapping matches', () => {
  assert.deepEqual(core.bindingSites('AAAAA', 'AAA').map(h => h.start), [0, 1, 2]);
  const result = core.analyzeTemplate('AAAGGAAACCCTTT', 'AAA', 'GGG');
  assert.ok(result.products.length > 1);
  assert.ok(result.products.some(p => p.bindings.some(b => b.right.primer === b.left.primer)));
});
test('same primer in both directions, FF/RR and deduplicated products preserve binding alternatives', () => {
  const result = core.analyzeTemplate('AGTCGGGACT', 'AGTC', 'AGTC');
  assert.equal(result.products.length, 1); assert.equal(result.products[0].bindings.length, 4);
  assert.deepEqual(new Set(result.products[0].bindings.map(b => b.right.primer + b.left.primer)), new Set(['FF', 'FR', 'RF', 'RR']));
  assert.deepEqual(core.bindingSites('ATAT', 'AT').map(h => h.direction), ['right', 'right', 'left', 'left']);
});
test('outward, overlap, adjacency and range boundaries', () => {
  assert.equal(core.analyzeTemplate('GACTGGAGTC', 'AGTC', 'AGTC').products.length, 0);
  const hit = (start, end, direction) => ({ start, end, direction });
  assert.throws(() => core.productFromHits('ACGT', hit(0, 3, 'right'), hit(2, 4, 'left')), /지원 범위/);
  assert.equal(core.productFromHits('ACGT', hit(0, 2, 'right'), hit(2, 4, 'left')).length, 4);
  for (const args of [[0, 20, 420], [10, 9, 420], [1, 421, 420], [1.5, 20, 420]]) assert.throws(() => core.toInterval(...args));
  assert.equal(core.primerFromRange(templates.A.sequence, 281, 300, 'left'), candidatePairs.P1.reverse);
});
test('normalization and explicit input errors', () => {
  assert.equal(core.normalizeSequence(' ac g\nT '), 'ACGT');
  for (const p of ['', ' ', 'ACNT', '>id\nACGT', 'A'.repeat(101)]) assert.throws(() => core.primerStats(p));
  assert.throws(() => core.normalizeSequence('AC N'), /3번 N/);
  assert.throws(() => core.analyzeTemplate('A'.repeat(1000) + 'T'.repeat(1000), 'A', 'T'), /100,000/);
});
test('mixtures keep sources even with same lengths', () => {
  const p = candidatePairs.P2, result = core.analyzeAll(templates, p.forward, p.reverse);
  for (const source of ['A', 'B']) {
    const mixed = core.combineProducts(result[source].products, result.C.products, result.C.products);
    assert.deepEqual(mixed.map(x => x.source), [source, 'C']);
  }
  const same = core.combineProducts([{ ...result.A.products[0], source: 'A' }], [{ ...result.A.products[0], source: 'C' }]);
  assert.equal(same.length, 2);
  assert.ok(core.gelPosition(180) > core.gelPosition(260));
});
test('record round trip, size, versions, malicious keys and immutable snapshots', () => {
  const state = emptyRecord(); state.answers.first = '<img src=x onerror=alert(1)>';
  state.designs.push({ id: 1, createdAt: state.updatedAt, ...candidatePairs.P1, prediction: '예측', reason: '이유', unresolved: '미확인' });
  const restored = parseRecord(JSON.stringify(state), ['first']);
  assert.deepEqual(restored, state); state.draft.forward = 'AAAA'; assert.equal(restored.designs[0].forward, candidatePairs.P1.forward);
  assert.throws(() => parseRecord(' '.repeat(1000001)));
  assert.throws(() => validateRecord({ ...state, schemaVersion: 2 }));
  assert.throws(() => parseRecord(JSON.stringify({ ...state, answers: JSON.parse('{"__proto__":"x"}') })));
  assert.throws(() => validateRecord({ ...state, designs: [{ ...state.designs[0], forward: 'ACN' }] }));
});
