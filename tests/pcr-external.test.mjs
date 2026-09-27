import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { emptyRecord, parseRecord } from '../assets/js/pcr-records.mjs';
import { emptyExternalSearch, emptyCandidate, currentPrimer, claimScope, completionIssues, comparisonRecorded, tmDifference, validateExternalSearch } from '../assets/js/pcr-external.mjs';
const fixture = JSON.parse(readFileSync(new URL('../assets/data/pcr-primer-fixture.json', import.meta.url)));
const complete = () => {
  const s = emptyExternalSearch(); s.status = 'performed'; s.searchRoute = 'paper';
  s.conditions = { ...s.conditions, date: '2026-09-27', organism: 'Homo sapiens', database: 'Recorded database', target: 'recorded target', specificity: 'recorded settings', forward: 'ACGTACGT', reverse: 'TGCATGCA' };
  s.candidates[0] = { ...emptyCandidate(), forward: 'ACGTACGT', reverse: 'TGCATGCA', product: '150', tmF: '60.2', tmR: '61.7', unintended: 'none' }; s.selectedCandidate = 'A'; return s;
};
test('06 reads the same explicit 04 snapshot or valid draft without mutation', () => {
  const record = emptyRecord(); assert.equal(currentPrimer(record, fixture), null);
  record.draft = { ...fixture.candidatePairs.P1 }; const before = structuredClone(record);
  assert.equal(currentPrimer(record, fixture).forward, record.draft.forward); assert.deepEqual(record, before);
  record.designs = [{ id: 1, ...fixture.candidatePairs.P2 }]; record.review.designId = 1;
  assert.equal(currentPrimer(record, fixture).forward, fixture.candidatePairs.P2.forward);
  record.draft.forward = 'invalid'; assert.equal(currentPrimer(record, fixture).forward, fixture.candidatePairs.P2.forward);
  record.review.designId = null; assert.equal(currentPrimer(record, fixture), null); assert.equal(currentPrimer(before, null), null);
});
test('old v1 and Phase 1 through 4 records preserve prior states and legacy 06 answers', () => {
  const optional = ['initialPrimerPrediction', 'introView', 'workbench', 'review', 'evidence'];
  for (let phase = 0; phase <= 4; phase++) {
    const old = emptyRecord(); delete old.externalSearch;
    const keep = phase === 0 ? 0 : phase + 1;
    optional.slice(keep).forEach(key => delete old[key]);
    old.answers = { 'external-mode': '새 후보 설계', 'external-status': '실시 후 학생이 기록', 'external-plan': 'old plan', 'evidence-cases': 'old evidence' };
    const restored = parseRecord(JSON.stringify(old)); assert.deepEqual(restored.externalSearch, emptyExternalSearch());
    assert.deepEqual(restored.answers, old.answers); for (const key of optional) if (old[key]) assert.deepEqual(restored[key], old[key]);
  }
});
test('all routes, plans, conditions, candidates and reflections round trip independently', () => {
  for (const route of ['mine','paper','new']) {
    const record = emptyRecord(); record.externalSearch = complete(); const s = record.externalSearch;
    s.route = route; s.paper.source = 'DOI or note'; s.paper.target = ''; s.design.type = 'fasta'; s.design.target = '>target\nACGT';
    s.plan = { purpose: 'purpose', organism: 'species', database: 'db', target: 'accession', notes: 'scope' };
    s.candidates.push({ ...s.candidates[0], unintended: 'reported', observations: 'other accession' });
    s.selectedCandidate = 'B'; s.selectionReason = 'reason'; s.claimReflection = 'limits'; s.wetLabReflection = 'experiment'; s.status = 'recorded';
    assert.deepEqual(parseRecord(JSON.stringify(record)), record);
  }
});
test('unfinished numeric, sequence and date text persists safely without completing', () => {
  const s = complete(); s.conditions.date = '2026-02-30'; s.candidates[0].product = '150?'; s.candidates[0].tmF = '6.'; s.paper.forward = 'AC?';
  assert.deepEqual(validateExternalSearch(s), s); assert.ok(completionIssues(s).length >= 3); assert.equal(claimScope(s), '');
  s.status = 'recorded'; assert.throws(() => validateExternalSearch(s));
});
test('invalid optional structure, unknown keys and enums cannot enter state', () => {
  for (const patch of [null, [], { route: 'bad' }, { status: 'verified' }, { candidates: [] }, { candidates: [emptyCandidate(), emptyCandidate(), emptyCandidate()] }, { selectedCandidate: 'B' }, { paper: { source: 10 } }, { extra: 'unexpected' }, { status: 'performed' }]) {
    assert.throws(() => validateExternalSearch(patch === null || Array.isArray(patch) ? patch : { ...emptyExternalSearch(), ...patch }));
  }
  const s = complete(); s.candidates[0].unintended = 'pass'; assert.throws(() => validateExternalSearch(s));
});
test('claim always binds student report to date organism database and selected candidate', () => {
  const s = complete(); const statement = claimScope(s);
  for (const v of ['2026-09-27', 'Homo sapiens', 'Recorded database', 'Candidate A', '기록했습니다']) assert.ok(statement.includes(v));
  assert.doesNotMatch(statement, /완전 특이적|PASS|인증|off-target이 없습니다/);
  s.candidates.push({ ...s.candidates[0], unintended: 'reported' }); s.selectedCandidate = 'B'; assert.match(claimScope(s), /보고되었다고/);
  s.candidates[1].unintended = 'unclear'; assert.match(claimScope(s), /해석하지 못했다고/);
  s.status = 'unperformed'; assert.equal(claimScope(s), ''); s.status = 'performed'; s.selectedCandidate = ''; assert.equal(claimScope(s), '');
});
test('comparison needs two recorded candidates and only subtracts reported temperatures', () => {
  const s = complete(); assert.equal(comparisonRecorded(s), false); assert.equal(tmDifference(s.candidates[0]), '1.5 °C');
  s.candidates.push(emptyCandidate()); assert.equal(comparisonRecorded(s), false);
  s.candidates[1] = { ...s.candidates[0] }; assert.equal(comparisonRecorded(s), true); assert.equal(s.selectedCandidate, 'A');
  s.candidates[1].tmF = ''; assert.equal(tmDifference(s.candidates[1]), '기록 전');
});
