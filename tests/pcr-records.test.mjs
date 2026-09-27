import test from 'node:test';
import assert from 'node:assert/strict';
import { emptyRecord, parseRecord, validateRecord } from '../assets/js/pcr-records.mjs';

test('legacy v1 records keep answers, draft, designs and three-step meaning', () => {
  for (const cycle of [0, 1, 2]) {
    const old = emptyRecord(); delete old.initialPrimerPrediction; delete old.introView;
    old.cycle = cycle; old.answers = { 'first-placement': '예측', 'first-negative': '대조군 필요' };
    old.draft.forward = 'ACGT';
    const restored = parseRecord(JSON.stringify(old), Object.keys(old.answers));
    assert.deepEqual(restored.answers, old.answers); assert.deepEqual(restored.draft, old.draft);
    assert.equal(restored.introView.stage, ['denaturation', 'annealing', 'extension'][cycle]);
    assert.equal(restored.initialPrimerPrediction.forward, null);
  }
});
test('initial prediction and intro choices survive JSON independently of final designs', () => {
  const state = emptyRecord(); state.initialPrimerPrediction.forward = 15; state.initialPrimerPrediction.reverse = 80;
  state.answers = { 'first-placement': '결실 양옆', 'cycle-boundary-choice': 'primer-pair', 'direction-end-choice': '3' };
  state.introView = { stage: 'extension', comparisonCycle: 3, flipped: true, arrangement: 'parallel', complementConfirmed: true, reverseShown: true };
  const restored = parseRecord(JSON.stringify(state), Object.keys(state.answers));
  assert.deepEqual(restored, state);
  state.initialPrimerPrediction.forward = 30; state.draft.forward = 'AAAA';
  assert.equal(restored.initialPrimerPrediction.forward, 15);
});
test('invalid initial prediction coordinates and reference cannot enter saved state', () => {
  const state = emptyRecord();
  for (const forward of [-5, 0, 100, 11, 12.5, '15', NaN, Infinity, {}, undefined]) {
    assert.throws(() => validateRecord({ ...state, initialPrimerPrediction: { ...state.initialPrimerPrediction, forward } }));
  }
  for (const patch of [{ reference: 'B' }, { units: 'bp' }, { reverse: 421 }]) assert.throws(() => validateRecord({ ...state, initialPrimerPrediction: { ...state.initialPrimerPrediction, ...patch } }));
});
test('intro state and multiple choice inputs are strictly validated', () => {
  const state = emptyRecord();
  for (const patch of [{ stage: 'other' }, { comparisonCycle: 4 }, { flipped: 'true' }, { arrangement: 'other' }, { reverseShown: 1 }, { complementConfirmed: null }]) assert.throws(() => validateRecord({ ...state, introView: { ...state.introView, ...patch } }));
  for (const key of ['cycle-boundary-choice', 'direction-end-choice']) assert.throws(() => validateRecord({ ...state, answers: { [key]: '<script>' } }));
});
