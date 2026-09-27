import test from 'node:test';
import assert from 'node:assert/strict';
import { predictSharedPair } from '../assets/js/pcr-intro.mjs';
import { emptyRecord, parseRecord } from '../assets/js/pcr-records.mjs';

const model = (forward, reverse) => predictSharedPair({ forward, reverse });
const near = (actual, expected) => assert.ok(Math.abs(actual - expected) < 1e-9, `${actual} != ${expected}`);

test('00 has no default answer and diagnoses incomplete, overlap and reversed pairs', () => {
  assert.equal(model(null, null).state, 'incomplete');
  assert.equal(model(15, null).reverse.b, null);
  assert.equal(model(null, 80).state, 'incomplete');
  for (const value of [30, 35, 40, 45, 50]) {
    assert.equal(model(value, null).state, 'overlap');
    assert.equal(model(value, 80).forward.b, null);
    assert.equal(model(15, value).reverse.b, null);
  }
  for (const pair of [[80, 15], [15, 15], [15, 20]]) assert.equal(model(...pair).state, 'order');
});

test('00 accepts multiple pairs and maps the deletion on a common physical scale', () => {
  let valid = 0, outside = 0;
  for (let f = 5; f <= 95; f += 5) for (let r = 5; r <= 95; r += 5) {
    const m = model(f, r);
    if (m.state === 'valid') {
      valid++;
      assert.ok(f <= 25 && r >= 55);
      // Physical positions in B: unchanged on the left, 80 bp shorter on the right.
      near(m.forward.b / 100 * 340, f / 100 * 420);
      near(m.reverse.b / 100 * 340, r / 100 * 420 - 80);
      near((r - f) / 100 * 420 - (m.reverse.b - m.forward.b) / 100 * 340, 80);
    }
    if (m.state === 'outside') {
      outside++;
      near((r - f) / 100 * 420, (m.reverse.b - m.forward.b) / 100 * 340);
    }
  }
  assert.equal(valid, 45);
  assert.ok(outside > 1);
  // Center is outside the deletion, but the illustrative binding band crosses its edge.
  assert.equal(model(15, 50).state, 'overlap');
  assert.equal(model(25, 55).state, 'valid');
});

test('00 derives the B display without changing the historical record format or later state', () => {
  const record = emptyRecord();
  record.initialPrimerPrediction.forward = 15; record.initialPrimerPrediction.reverse = 80;
  record.answers['first-placement'] = '같은 pair를 비교한다. C는 아직 판단하지 않는다.';
  const before = JSON.stringify(record);
  predictSharedPair(record.initialPrimerPrediction);
  assert.equal(JSON.stringify(record), before);
  assert.deepEqual(parseRecord(before), record);
  delete record.initialPrimerPrediction;
  assert.equal(predictSharedPair(parseRecord(JSON.stringify(record)).initialPrimerPrediction).state, 'incomplete');
});
