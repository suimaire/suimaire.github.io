import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { primerStats } from '../assets/js/pcr-core.mjs';
import { emptyRecord, parseRecord } from '../assets/js/pcr-records.mjs';
import { reviewDesign, endFeatures, complementarity, alignmentLines, emptyReview } from '../assets/js/pcr-review.mjs';
const fixture = JSON.parse(readFileSync(new URL('../assets/data/pcr-primer-fixture.json', import.meta.url)));
const design = (id, pair) => ({ id, createdAt: new Date().toISOString(), ...pair, prediction: '', reason: '', unresolved: '' });

test('04 reads valid draft, explicit saved design takes priority, never invents a fallback', () => {
  const state = emptyRecord(); assert.equal(reviewDesign(state, fixture), null);
  state.draft = { ...fixture.candidatePairs.P1 };
  let model = reviewDesign(state, fixture);
  assert.equal(model.designId, null); assert.equal(model.result.A.products[0].length, 260);
  state.designs = [design(1, fixture.candidatePairs.P3)]; state.review.designId = 1;
  const before = JSON.stringify(state);
  model = reviewDesign(state, fixture);
  assert.equal(model.designId, 1); assert.equal(model.primers.F.stats.sequence, fixture.candidatePairs.P3.forward);
  assert.deepEqual(model.result.B.products, []); assert.equal(JSON.stringify(state), before);
  state.review.designId = null; state.draft.forward = 'N';
  assert.equal(reviewDesign(state, fixture), null); assert.equal(reviewDesign(state, null), null);
});

test('invalid coordinate metadata cannot leak stale review values', () => {
  const state = emptyRecord(); state.draft = { ...fixture.candidatePairs.P1, bindings: { F: { start: '1', end: '20', direction: 'right' }, R: null } };
  assert.equal(reviewDesign(state, fixture), null);
});

test('composition, Wallace Tm, delta and last five bases come from actual sequences', () => {
  const f = primerStats('AATGCC'), r = primerStats('TTTGGCC');
  assert.equal(f.length, 6); assert.equal(f.gcPercent, 50); assert.equal(f.simpleTm, 18);
  assert.equal(r.simpleTm, 22); assert.equal(Math.abs(f.simpleTm - r.simpleTm), 4);
  assert.deepEqual(endFeatures('AATGCC'), { sequence: 'ATGCC', lastBase: 'C', gcCount: 3, runs: [{ sequence: 'CC', start: 4, length: 2 }] });
  assert.deepEqual(endFeatures(' g '), { sequence: 'G', lastBase: 'G', gcCount: 1, runs: [] });
  assert.equal(endFeatures('AAAAAA').runs[0].length, 5);
});

test('antiparallel matching uses reverse complement and handles zero and full matches', () => {
  assert.equal(complementarity('AAAA').longest, null);
  const pair = complementarity('AAGC', 'GCTT');
  assert.deepEqual(pair.longest, { offset: 0, aStart: 0, bRCStart: 0, length: 4, aThreePrime: true, bThreePrime: true });
  assert.deepEqual(alignmentLines(pair), { top: 'AAGC', bars: '||||', bottom: 'TTCG', topPad: 0, bottomPad: 0, matchStart: 0 });
  assert.equal(complementarity('ATGCAT').longest.length, 6);
  assert.equal(alignmentLines(complementarity('AAAA')), null);
});

test('terminal involvement distinguishes a 3′ from b 3′ and shorter terminal alternatives', () => {
  const fEnd = complementarity('AAAC', 'AGGA');
  assert.equal(fEnd.threePrime.aThreePrime, true); assert.equal(fEnd.threePrime.bThreePrime, false);
  const rEnd = complementarity('AGGA', 'AAAC');
  assert.equal(rEnd.threePrime.aThreePrime, false); assert.equal(rEnd.threePrime.bThreePrime, true);
  const internal = complementarity('CCCCAG', 'GGGGTTTTC');
  assert.equal(internal.longest.length, 4);
  assert.equal(internal.longest.aThreePrime, false); assert.equal(internal.longest.bThreePrime, false);
  assert.equal(internal.bothThreePrime.length, 1);
  assert.equal(internal.bothThreePrime.aThreePrime, true); assert.equal(internal.bothThreePrime.bThreePrime, true);
});

test('all short sequences agree with an independent antiparallel substring oracle', () => {
  const bases = 'ACGT', comp = { A: 'T', C: 'G', G: 'C', T: 'A' };
  const sequences = [...bases].flatMap(a => [...bases].flatMap(b => [...bases].map(c => a + b + c)));
  for (const a of sequences) for (const b of sequences) {
    let max = 0, maxEnd = 0, maxBoth = 0;
    for (let i = 0; i < a.length; i++) for (let j = 0; j < b.length; j++) {
      let n = 0;
      while (i + n < a.length && j - n >= 0 && comp[a[i + n]] === b[j - n]) {
        n++; max = Math.max(max, n);
        if (i + n === a.length || j === b.length - 1) maxEnd = Math.max(maxEnd, n);
        if (i + n === a.length && j === b.length - 1) maxBoth = Math.max(maxBoth, n);
      }
    }
    const result = complementarity(a, b);
    assert.equal(result.longest?.length ?? 0, max, `${a}/${b}`);
    assert.equal(result.threePrime?.length ?? 0, maxEnd);
    assert.equal(result.bothThreePrime?.length ?? 0, maxBoth);
    if (result.longest) {
      const lines = alignmentLines(result), run = result.longest;
      for (let i = lines.matchStart; i < lines.matchStart + run.length; i++) assert.equal(comp[lines.top[i]], lines.bottom[i]);
    }
  }
});

test('both fixture self observations and pair observations are deterministic and contain no thermodynamic values', () => {
  for (const pair of Object.values(fixture.candidatePairs)) for (const [a, b] of [[pair.forward, pair.forward], [pair.reverse, pair.reverse], [pair.forward, pair.reverse]]) {
    const result = complementarity(a, b);
    assert.deepEqual(result, complementarity(a, b));
    assert.deepEqual(Object.keys(result), ['a', 'b', 'longest', 'threePrime', 'bothThreePrime']);
    assert.ok(result.longest.length > 0);
  }
  assert.throws(() => complementarity('N')); assert.throws(() => complementarity('A'.repeat(101)));
});

test('optional review settings and answers round trip without changing snapshots or schema v1', () => {
  const state = emptyRecord(); state.designs = [design(1, fixture.candidatePairs.P1)];
  state.review = { lens: 'off-target', designId: 1, stage: 2 };
  state.answers = { 'review-length-choice': 'no', 'review-unresolved': '더 넓은 검색', 'candidate-prediction': '이전 답안' };
  assert.deepEqual(parseRecord(JSON.stringify(state)), state);
  const legacy = { ...state }; delete legacy.review;
  const restored = parseRecord(JSON.stringify(legacy));
  assert.deepEqual(restored.review, emptyReview()); assert.deepEqual(restored.designs, state.designs);
  assert.deepEqual(restored.answers, state.answers);
  for (const patch of [{ lens: 'unknown' }, { stage: -1 }, { stage: '2' }, { designId: 2 }, { designId: '1' }]) assert.throws(() => parseRecord(JSON.stringify({ ...state, review: { ...state.review, ...patch } })));
  assert.throws(() => parseRecord(JSON.stringify({ ...state, review: null })));
  assert.throws(() => parseRecord(JSON.stringify({ ...state, answers: { 'review-length-choice': 'wrong' } })));
});
