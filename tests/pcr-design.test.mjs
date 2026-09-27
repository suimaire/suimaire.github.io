import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { inspectDesign, sequenceFromBinding, cloneBindings } from '../assets/js/pcr-design.mjs';
import { emptyRecord, parseRecord, validateRecord } from '../assets/js/pcr-records.mjs';
const fixture = JSON.parse(readFileSync(new URL('../assets/data/pcr-primer-fixture.json', import.meta.url)));
const binding = (start, end, direction) => ({ start: String(start), end: String(end), direction });
function draft(f, r) {
  return { forward: f ? sequenceFromBinding(fixture.templates.A.sequence, f) : '', reverse: r ? sequenceFromBinding(fixture.templates.A.sequence, r) : '', bindings: { F: f, R: r } };
}
const normal = () => draft(binding(41, 60, 'right'), binding(281, 300, 'left'));
test('selected ranges reuse exact fixture calculations, reverse complement and GC', () => {
  const d = normal(), model = inspectDesign(d, fixture);
  assert.equal(d.reverse, fixture.candidatePairs.P1.reverse);
  assert.equal(model.primers.R.reference, 'TCGTACGATCGTAGCCTGAA');
  for (const p of Object.values(model.primers)) { assert.equal(p.stats.length, 20); assert.equal(p.stats.gcPercent, 50); }
  assert.equal(model.selectedProduct.length, 260);
  assert.deepEqual(Object.values(model.result).map(r => r.products.map(p => p.length)), [[260], [180], []]);
  assert.match(model.deletionMessage, /포함되어/);
});
test('deletion inside, partial overlap, full overlap and adjacent binding boundaries', () => {
  for (const [start, end] of [[111, 130], [121, 140], [191, 210]]) {
    const m = inspectDesign(draft(binding(start, end, 'right'), binding(281, 300, 'left')), fixture);
    assert.equal(m.result.B.products.length, 0); assert.ok(m.result.A.products.length);
    assert.match(m.bindingMessage, /Forward.*결실.*존재하지 않습니다/);
    assert.match(m.deletionMessage, /겹칩니다/);
  }
  const reverseOverlap = inspectDesign(draft(binding(41, 60, 'right'), binding(191, 210, 'left')), fixture);
  assert.equal(reverseOverlap.result.B.products.length, 0); assert.match(reverseOverlap.bindingMessage, /Reverse/);
  const edges = inspectDesign(draft(binding(101, 120, 'right'), binding(201, 220, 'left')), fixture);
  assert.equal(edges.result.A.products[0].length - edges.result.B.products[0].length, 80);
  assert.match(edges.deletionMessage, /포함되어/);
  const outside = inspectDesign(draft(binding(231, 250, 'right'), binding(281, 300, 'left')), fixture);
  assert.equal(outside.result.A.products[0].length, outside.result.B.products[0].length);
  assert.match(outside.deletionMessage, /포함되지 않습니다/);
});
test('invalid orientation, reversed order and overlapping intervals stay unchanged', () => {
  for (const [f, r, text] of [
    [binding(321, 340, 'right'), binding(41, 60, 'left'), /inward-facing/],
    [binding(41, 60, 'right'), binding(281, 300, 'right'), /inward-facing/],
    [binding(41, 60, 'left'), binding(281, 300, 'right'), /inward-facing/],
    [binding(41, 60, 'right'), binding(51, 70, 'left'), /겹칩니다/]
  ]) {
    const d = draft(f, r), original = structuredClone(d), m = inspectDesign(d, fixture);
    assert.equal(m.selectedProduct, null); assert.match(m.placement, text); assert.deepEqual(d, original);
  }
});
test('incomplete, out of range, long and short choices are retained with explicit limits', () => {
  const partial = inspectDesign(draft(binding(41, 60, 'right'), null), fixture);
  assert.equal(partial.result, null); assert.equal(partial.primers.F.stats.length, 20);
  for (const [start, end] of [['', 60], [0, 60], [401, 430], [70, 41], [1.5, 20]]) {
    const d = normal(); d.bindings.F = binding(start, end, 'right'); d.forward = '';
    const m = inspectDesign(d, fixture); assert.equal(m.result, null); assert.match(m.errors[0], /좌표/);
    const record = emptyRecord(); record.draft = d;
    assert.deepEqual(parseRecord(JSON.stringify(record)).draft, d);
  }
  const long = inspectDesign(draft(binding(1, 101, 'right'), binding(281, 300, 'left')), fixture);
  assert.equal(long.primers.F.stats.length, 101); assert.match(long.errors[0], /100 nt/); assert.equal(long.result, null);
  const short = inspectDesign(draft(binding(41, 45, 'right'), binding(281, 300, 'left')), fixture);
  assert.equal(short.primers.F.stats.length, 5); assert.ok(short.result);
});
test('legacy sequences infer a unique binding without altering draft or suppressing C products', () => {
  for (const id of ['P1', 'P2', 'P3']) {
    const original = { ...fixture.candidatePairs[id] }, m = inspectDesign(original, fixture);
    assert.ok(m.primers.F.inferred && m.primers.R.inferred);
    for (const source of ['A', 'B', 'C']) assert.deepEqual(m.result[source].products.map(p => [p.start + 1, p.end, p.length]), fixture.expectedExactMatchProducts[id][source]);
    assert.deepEqual(original, fixture.candidatePairs[id]);
  }
});
test('repeated binding sequences are not assigned invented coordinates', () => {
  const m = inspectDesign({ forward: 'A', reverse: fixture.candidatePairs.P1.reverse }, fixture);
  assert.equal(m.primers.F.binding, null); assert.ok(m.primers.F.hits.length > 1); assert.ok(m.result);
});
test('workbench and independent saved bindings round trip with legacy fields', () => {
  const record = emptyRecord(); record.draft = normal(); record.workbench = { mode: 'R', windowStart: 276, showPrediction: true };
  record.designs.push({ id: 1, createdAt: record.updatedAt, ...record.draft, bindings: cloneBindings(record.draft.bindings), prediction: '이전 예측', reason: '위치 변경', unresolved: '미확인' });
  assert.deepEqual(parseRecord(JSON.stringify(record)), record);
  record.draft.bindings.F.start = '111'; assert.equal(record.designs[0].bindings.F.start, '41');
  const old = emptyRecord(); delete old.workbench;
  assert.equal(validateRecord(old).workbench.mode, 'F');
  for (const patch of [{ mode: 'X' }, { windowStart: 0 }, { windowStart: 421 }, { showPrediction: 'yes' }]) assert.throws(() => validateRecord({ ...record, workbench: { ...record.workbench, ...patch } }));
  for (const b of [{ start: '<script>', end: '20', direction: 'left' }, { start: '1', end: '20', direction: 'up' }, { start: 1, end: 20, direction: 'left' }]) assert.throws(() => validateRecord({ ...record, draft: { ...record.draft, bindings: { F: b, R: null } } }));
});
test('mismatched imported coordinates cannot display stale products', () => {
  const d = normal(); d.bindings.F.start = '51';
  const m = inspectDesign(d, fixture); assert.equal(m.result, null); assert.match(m.errors[0], /일치하지 않습니다/);
});
