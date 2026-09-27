import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { emptyRecord, parseRecord, validateRecord } from '../assets/js/pcr-records.mjs';
import { finalCandidates, selectFinalDesign, inspectFinal, roughPosition, primerReviewRows, productText, externalSummary, limitationRows, LEGACY_FINAL } from '../assets/js/pcr-final.mjs';
import { sequenceFromBinding } from '../assets/js/pcr-design.mjs';
const fixture = JSON.parse(readFileSync(new URL('../assets/data/pcr-primer-fixture.json', import.meta.url)));
const pair = (f = 41, r = 281) => {
  const bindings = { F: { start: String(f), end: String(f + 19), direction: 'right' }, R: { start: String(r), end: String(r + 19), direction: 'left' } };
  return { forward: sequenceFromBinding(fixture.templates.A.sequence, bindings.F), reverse: sequenceFromBinding(fixture.templates.A.sequence, bindings.R), bindings };
};
const withDesigns = n => { const s = emptyRecord(); s.draft = pair(); s.designs = Array.from({ length: n }, (_, i) => ({ id: i + 1, createdAt: s.updatedAt, ...pair(41 + i * 35), reason: `수정 ${i}`, prediction: '', unresolved: '' })); return s; };
test('07 initial positions remain qualitative percentages, including absent values', () => {
  assert.equal(roughPosition(null), '기록 없음');
  for (let v = 5; v <= 95; v += 5) { assert.match(roughPosition(v), /대략적 위치/); assert.doesNotMatch(roughPosition(v), /bp|\d/); }
  assert.match(roughPosition(15), /왼쪽/); assert.match(roughPosition(35), /안쪽/); assert.match(roughPosition(80), /오른쪽/);
});
test('07 offers exactly existing valid saved designs and draft in saved order without choosing', () => {
  for (let n = 0; n <= 3; n++) {
    const s = withDesigns(n), before = structuredClone(s);
    assert.deepEqual(finalCandidates(s, fixture).map(c => c.source), [...Array.from({ length: n }, (_, i) => `design-${i + 1}`), 'draft']);
    assert.deepEqual(s, before); assert.equal(s.finalReview.selectedDesignSource, '');
  }
  const s = emptyRecord(); assert.deepEqual(finalCandidates(s, fixture), []);
  s.draft = pair(); s.draft.bindings.R.direction = 'right'; assert.deepEqual(finalCandidates(s, fixture), []);
  s.draft = pair(41, 51); assert.deepEqual(finalCandidates(s, fixture), []);
  s.draft = pair(); s.draft.forward = 'AC?'; assert.deepEqual(finalCandidates(s, fixture), []);
  assert.deepEqual(finalCandidates(withDesigns(3), null), []);
});
test('07 explicit selection of each saved design is independent and preserves all prior state', () => {
  for (let n = 1; n <= 3; n++) {
    const s = withDesigns(3), before = structuredClone(s);
    selectFinalDesign(s, `design-${n}`, fixture);
    assert.equal(s.finalReview.primerSnapshot.forward, s.designs[n - 1].forward);
    for (const key of Object.keys(before).filter(k => k !== 'finalReview')) assert.deepEqual(s[key], before[key]);
    s.draft = pair(130); assert.equal(s.finalReview.primerSnapshot.forward, before.designs[n - 1].forward);
    assert.deepEqual(parseRecord(JSON.stringify(s)).finalReview, s.finalReview);
  }
});
test('07 final draft snapshot survives edits, including nested bindings, and explicit reselection updates it', () => {
  const s = withDesigns(1); selectFinalDesign(s, 'draft', fixture);
  const snapshot = structuredClone(s.finalReview.primerSnapshot);
  s.draft.bindings.F.start = '50'; s.draft.forward = 'AAA';
  assert.deepEqual(s.finalReview.primerSnapshot, snapshot);
  assert.throws(() => selectFinalDesign(s, 'draft', fixture));
  assert.deepEqual(s.finalReview.primerSnapshot, snapshot);
  s.draft = pair(111); selectFinalDesign(s, 'draft', fixture); assert.notEqual(s.finalReview.primerSnapshot.forward, snapshot.forward);
});
test('07 F/R coordinates, stats, products and 04 calculations reuse the existing model', () => {
  const s = withDesigns(1); selectFinalDesign(s, 'design-1', fixture);
  const m = inspectFinal(s.finalReview.primerSnapshot, fixture);
  assert.equal(m.primers.F.stats.sequence, fixture.candidatePairs.P1.forward);
  assert.equal(m.primers.R.stats.sequence, fixture.candidatePairs.P1.reverse);
  assert.equal(m.primers.F.binding.start, '41'); assert.equal(m.primers.R.binding.end, '300');
  assert.equal(m.primers.F.stats.length, 20); assert.equal(m.primers.F.stats.gcPercent, 50);
  assert.deepEqual(productText(m.result), [['A', '260 bp'], ['B', '180 bp'], ['C', '예상 산물 없음']]);
  const rows = primerReviewRows(m); assert.match(rows[1][1], /Forward 60 °C \/ Reverse 60 °C \/ 차이 0 °C/);
  assert.match(rows[2][1], /3|GACCA/); assert.match(rows[3][1], /연속 상보/);
  assert.doesNotMatch(JSON.stringify(rows), /NaN|undefined|ΔG/);
});
test('07 external records retain status, route, candidate, off-target, date and bounded claim', () => {
  const s = withDesigns(1), e = s.externalSearch;
  assert.equal(externalSummary(e).rows[0][1], '미실시'); assert.equal(externalSummary(e).claim, '');
  e.status = 'performed'; e.searchRoute = 'paper'; e.selectedCandidate = 'A';
  Object.assign(e.conditions, { date: '2026-09-27', organism: '교육용 organism', database: '교육용 DB', target: 'target 예시' });
  Object.assign(e.candidates[0], { ...fixture.candidatePairs.P1, product: '260', unintended: 'reported', observations: 'accession 예시' });
  const summary = externalSummary(e, pair());
  assert.match(summary.note, /F\/R이 최종 pair와 같습니다/); assert.match(summary.claim, /2026-09-27.*교육용 DB.*보고되었다고/);
  assert.match(JSON.stringify(summary.rows), /accession 예시/); assert.match(JSON.stringify(summary.rows), /논문의 primer/);
  assert.match(externalSummary(e, pair(111)).note, /일치가 확인되지/);
  e.status = 'unperformed'; assert.doesNotMatch(JSON.stringify(externalSummary(e)), /accession 예시|교육용 DB/);
});
test('07 limitations never convert simulated gel, simple complementarity or external status to wet-lab validation', () => {
  const s = emptyRecord();
  for (const status of ['unperformed', 'performed', 'recorded']) {
    s.externalSearch.status = status; const rows = limitationRows(s);
    assert.match(rows[0][1], /미수행/); assert.match(rows[1][1], /미수행.*가상/); assert.equal(rows[2][1], '미확인'); assert.match(rows[3][1], /미검증.*실제 검증이 아닙니다/);
  }
});
test('07 old v1 and Phase 1–5 migrate clear meanings while preserving all legacy originals', () => {
  for (let phase = 0; phase <= 5; phase++) {
    const s = withDesigns(phase ? 1 : 0); delete s.finalReview;
    if (phase < 5) delete s.externalSearch; if (phase < 4) delete s.evidence; if (phase < 3) delete s.review;
    if (phase < 2) delete s.workbench; if (phase < 1) { delete s.introView; delete s.initialPrimerPrediction; }
    s.answers = Object.fromEntries(Object.keys(LEGACY_FINAL).map(k => [k, `<b>${k} 원문</b>`]));
    const restored = parseRecord(JSON.stringify(s), Object.keys(LEGACY_FINAL));
    assert.deepEqual(restored.answers, s.answers); assert.equal(restored.finalReview.researchQuestion, s.answers['final-question']);
    assert.equal(restored.finalReview.finalRationale, s.answers['final-evidence']); assert.equal(restored.finalReview.revisionReflection, s.answers['final-revision']);
    assert.equal(restored.finalReview.primerSnapshot, null); assert.equal(restored.finalReview.controls.positive, ''); assert.equal(restored.finalReview.finalAssessment, '');
    restored.finalReview.researchQuestion = ''; assert.equal(parseRecord(JSON.stringify(restored)).finalReview.researchQuestion, '');
  }
});
test('07 all controls, long text, expansion and selected pair round trip with schema v1', () => {
  const s = withDesigns(2); selectFinalDesign(s, 'draft', fixture);
  Object.assign(s.finalReview, { researchQuestion: '연구 질문', finalRationale: '근거\n'.repeat(1000), revisionReflection: '수정', finalAssessment: '판단', otherLimitations: '기타', notebookExpanded: true });
  for (const k of Object.keys(s.finalReview.controls)) s.finalReview.controls[k] = `${k} 대조군`;
  assert.deepEqual(parseRecord(JSON.stringify(s)), s); assert.equal(s.schemaVersion, 1);
});
test('07 malformed imports fail before replacing state or silently inventing snapshots', () => {
  const s = withDesigns(1); selectFinalDesign(s, 'design-1', fixture);
  for (const patch of [{ selectedDesignSource: 'best' }, { notebookExpanded: 'yes' }, { controls: {} }, { primerSnapshot: null }, { selectedAt: 'bad' }, { unknown: 'field' }, { finalAssessment: 'a'.repeat(12001) }, { selectedDesignSource: '' }]) assert.throws(() => validateRecord({ ...s, finalReview: { ...s.finalReview, ...patch } }));
  for (const patch of [{ forward: 'AC?' }, { bindings: { F: null } }, { reverse: 7 }, { extra: 'x' }]) assert.throws(() => validateRecord({ ...s, finalReview: { ...s.finalReview, primerSnapshot: { ...s.finalReview.primerSnapshot, ...patch } } }));
});
