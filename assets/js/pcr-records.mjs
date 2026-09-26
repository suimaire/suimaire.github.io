import { normalizeSequence } from './pcr-core.mjs';
export const STORAGE_KEY = 'hafs:pcr-primer:v1';
export const DATA_VERSION = 'hafs-pcr-synthetic-20260926-v1';
export const MAX_IMPORT_BYTES = 1000000;
export const INTRO_CHOICE_KEYS = ['cycle-boundary-choice', 'direction-end-choice'];
export const PCR_STAGES = ['mixture', 'denaturation', 'annealing', 'extension'];
export function emptyIntroView() {
  return { stage: 'mixture', comparisonCycle: 1, flipped: false, arrangement: 'inward', complementConfirmed: false, reverseShown: false };
}
export function emptyRecord() {
  return { schemaVersion: 1, dataVersion: DATA_VERSION, updatedAt: new Date().toISOString(), answers: {}, draft: { forward: '', reverse: '' }, designs: [], cycle: 0,
    initialPrimerPrediction: { reference: 'A', units: 'relative-percent', forward: null, reverse: null }, introView: emptyIntroView() };
}
const object = x => x !== null && typeof x === 'object' && !Array.isArray(x);
function string(x, max = 12000) { if (typeof x !== 'string' || x.length > max) throw new Error('기록 문자열의 형식 또는 길이가 올바르지 않습니다.'); return x; }
function date(x) { string(x, 40); if (!Number.isFinite(Date.parse(x))) throw new Error('기록 시각이 올바르지 않습니다.'); return x; }
export function validateRecord(value, allowedKeys) {
  if (!object(value) || value.schemaVersion !== 1 || value.dataVersion !== DATA_VERSION) throw new Error('지원하지 않는 학습지 또는 데이터 버전입니다.');
  if (!object(value.answers) || !object(value.draft) || !Array.isArray(value.designs) || value.designs.length > 3) throw new Error('기록 구조가 올바르지 않습니다.');
  const answers = {};
  for (const [key, answer] of Object.entries(value.answers)) {
    if (!/^[a-z][a-z0-9-]{0,60}$/.test(key) || ['constructor', 'prototype'].includes(key) || (allowedKeys && !allowedKeys.includes(key))) throw new Error('알 수 없는 기록 항목입니다.');
    answers[key] = string(answer);
  }
  if (![0, 1, 2].includes(value.cycle)) throw new Error('PCR 단계 값이 올바르지 않습니다.');
  // Optional additions to v1: old records retain their original three-step position.
  const initialPrimerPrediction = value.initialPrimerPrediction ?? emptyRecord().initialPrimerPrediction;
  if (!object(initialPrimerPrediction) || initialPrimerPrediction.reference !== 'A' || initialPrimerPrediction.units !== 'relative-percent' ||
      !['forward', 'reverse'].every(key => initialPrimerPrediction[key] === null || (Number.isInteger(initialPrimerPrediction[key]) && initialPrimerPrediction[key] >= 5 && initialPrimerPrediction[key] <= 95 && initialPrimerPrediction[key] % 5 === 0))) throw new Error('초기 프라이머 위치 예측이 올바르지 않습니다.');
  const introView = value.introView ?? { ...emptyIntroView(), stage: PCR_STAGES[value.cycle + 1] };
  if (!object(introView) || !PCR_STAGES.includes(introView.stage) || ![1, 2, 3].includes(introView.comparisonCycle) ||
      !['inward', 'parallel'].includes(introView.arrangement) || !['flipped', 'complementConfirmed', 'reverseShown'].every(key => typeof introView[key] === 'boolean')) throw new Error('도입 활동 보기 설정이 올바르지 않습니다.');
  for (const [key, options] of [['cycle-boundary-choice', ['polymerase', 'primer-pair', 'dntp', 'buffer']], ['direction-end-choice', ['5', '3']]]) {
    if (answers[key] !== undefined && !options.includes(answers[key])) throw new Error('확인 질문의 선택값이 올바르지 않습니다.');
  }
  const designs = value.designs.map((d, i) => {
    if (!object(d) || d.id !== i + 1) throw new Error('설계 순서가 올바르지 않습니다.');
    return { id: d.id, createdAt: date(d.createdAt), forward: normalizeSequence(string(d.forward, 100)), reverse: normalizeSequence(string(d.reverse, 100)), prediction: string(d.prediction), reason: string(d.reason), unresolved: string(d.unresolved) };
  });
  return { schemaVersion: 1, dataVersion: DATA_VERSION, updatedAt: date(value.updatedAt), answers,
    draft: { forward: string(value.draft.forward, 1000), reverse: string(value.draft.reverse, 1000) }, designs, cycle: value.cycle,
    initialPrimerPrediction: { reference: 'A', units: 'relative-percent', forward: initialPrimerPrediction.forward, reverse: initialPrimerPrediction.reverse },
    introView: { stage: introView.stage, comparisonCycle: introView.comparisonCycle, flipped: introView.flipped, arrangement: introView.arrangement, complementConfirmed: introView.complementConfirmed, reverseShown: introView.reverseShown } };
}
export function parseRecord(raw, allowedKeys) {
  if (typeof raw !== 'string' || new TextEncoder().encode(raw).length > MAX_IMPORT_BYTES) throw new Error('가져올 JSON은 1 MB 이하여야 합니다.');
  let value;
  try { value = JSON.parse(raw); } catch { throw new Error('올바른 JSON 파일이 아닙니다.'); }
  return validateRecord(value, allowedKeys);
}
