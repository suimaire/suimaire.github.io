import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { EVIDENCE_CASES, EVIDENCE_LANES, emptyEvidence, evidenceBandPosition, caseSummary } from '../assets/js/pcr-evidence.mjs';
import { emptyRecord, validateRecord, parseRecord } from '../assets/js/pcr-records.mjs';
import { reviewDesign } from '../assets/js/pcr-review.mjs';
const fixture = JSON.parse(readFileSync(new URL('../assets/data/pcr-primer-fixture.json', import.meta.url)));

test('fixed observations distinguish size evidence, NTC signal, failed positive control and multiple products', () => {
  assert.deepEqual(EVIDENCE_CASES.map(c => Object.values(c.lanes).map(l => l.bands)), [
    [[500,400,300,200,100],[260],[260],[]], [[500,400,300,200,100],[260],[260],[70]],
    [[500,400,300,200,100],[],[],[]], [[500,400,300,200,100],[420,260,140],[260],[]]
  ]);
  for (const c of EVIDENCE_CASES) {
    assert.ok(c.purpose && c.expectedControlBehavior.positive && c.explanation && c.questions.length === 2);
    assert.equal(c.expectedSize, 260);
    for (const id of Object.keys(EVIDENCE_LANES)) assert.ok(c.lanes[id].observation && caseSummary(c).includes(c.lanes[id].label));
  }
  assert.equal(EVIDENCE_CASES.filter(c => c.additional).length, 1);
});
test('log positioning has equal log intervals and shorter fragments move farther within the lane', () => {
  const sizes = [500,420,400,300,260,200,140,100,70,50], positions = sizes.map(evidenceBandPosition);
  assert.ok(positions.every((p, i) => p > 0 && p < 100 && (i === 0 || p > positions[i-1])));
  assert.ok(Math.abs((evidenceBandPosition(100)-evidenceBandPosition(200))-(evidenceBandPosition(200)-evidenceBandPosition(400))) < 1e-10);
  for (const value of [0,49,501,NaN,Infinity]) assert.throws(() => evidenceBandPosition(value));
});
test('student selection and fixed observations are independent and never mutate draft, snapshots or 04', () => {
  const record = emptyRecord(), fixedBefore = JSON.stringify(EVIDENCE_CASES);
  record.draft = { ...fixture.candidatePairs.P3 };
  record.designs = [{ id: 1, createdAt: record.updatedAt, ...fixture.candidatePairs.P2, prediction: '', reason: '', unresolved: '' }];
  record.review.designId = 1;
  const before = JSON.stringify(record), model = reviewDesign(record, fixture);
  assert.equal(model.result.C.products[0].length, 220);
  assert.equal(JSON.stringify(record), before);
  record.review.designId = null;
  assert.equal(reviewDesign(record, fixture).result.A.products[0].length, 180);
  record.draft.forward = 'N'; assert.equal(reviewDesign(record, fixture), null);
  assert.equal(JSON.stringify(EVIDENCE_CASES), fixedBefore);
});
test('old v1 retains all old 05 answers and gets safe new defaults', () => {
  const record = emptyRecord(); delete record.evidence;
  record.answers = { 'evidence-cases': '과거 해석', 'evidence-controls': '과거 대조군', 'evidence-identity': '과거 정체' };
  const result = parseRecord(JSON.stringify(record), Object.keys(record.answers));
  assert.deepEqual(result.evidence, emptyEvidence()); assert.deepEqual(result.answers, record.answers);
  assert.equal(result.schemaVersion, 1);
});
test('every case selection and answer survives a JSON round trip with independent defaults', () => {
  const record = emptyRecord(); record.evidence.activeCase = 'case-3';
  for (const c of EVIDENCE_CASES) {
    record.evidence.selectedLanes[c.id] = 'positive';
    for (const key of ['observation','interpretation','uncertainty']) record.answers[`evidence-${c.id}-${key}`] = `${c.label} ${key}`;
  }
  assert.deepEqual(parseRecord(JSON.stringify(record), Object.keys(record.answers)), record);
  assert.equal(emptyRecord().evidence.selectedLanes['case-1'], null);
});
test('invalid optional evidence values and unexpected answer keys are rejected', () => {
  for (const value of [null, [], { activeCase:'case-7', selectedLanes:{} }, { activeCase:'case-1', selectedLanes:{} }, { ...emptyEvidence(), selectedLanes:{ ...emptyEvidence().selectedLanes, 'case-1':'constructor' } }]) {
    assert.throws(() => validateRecord({ ...emptyRecord(), evidence: value }), /가상 관찰/);
  }
  const record = emptyRecord(); record.answers['evidence-case-1-interpretation'] = 'x';
  assert.throws(() => validateRecord(record, []), /알 수 없는/);
});
