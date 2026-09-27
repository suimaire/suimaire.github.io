import { inspectDesign, cloneBindings } from './pcr-design.mjs';
import { endFeatures, complementarity } from './pcr-review.mjs';
import { ROUTES, UNINTENDED, claimScope, CLAIM_LIMIT } from './pcr-external.mjs';

export const LEGACY_FINAL = {
  'final-question': '연구 질문', 'final-f': '이전 Forward / 5′→3′', 'final-r': '이전 Reverse / 5′→3′',
  'final-products': '예상 산물과 출처', 'final-evidence': '관찰 근거와 최종 선택 이유',
  'final-controls': '대조군 계획', 'final-unknown': '확인하지 못한 사항',
  'final-revision': '처음 설명에서 수정한 점', 'final-limits': '계산과 실험의 한계'
};
export const emptyFinalReview = () => ({
  selectedDesignSource: '', selectedAt: '', primerSnapshot: null,
  researchQuestion: '', finalRationale: '', revisionReflection: '',
  controls: { positive: '', negative: '', additional: '', interpretationLimit: '' },
  otherLimitations: '', finalAssessment: '', notebookExpanded: false
});
export function validateFinalReview(value, answers, bindingFields) {
  const defaults = emptyFinalReview();
  if (value === undefined) {
    // Only unambiguous meanings migrate. Keep every original answer unchanged.
    for (const [key, old] of Object.entries({ researchQuestion: 'final-question', finalRationale: 'final-evidence', revisionReflection: 'final-revision', otherLimitations: 'final-unknown' })) defaults[key] = answers[old] || '';
    return defaults;
  }
  const fail = () => { throw new Error('최종 연구 노트의 기록 구조 또는 선택값이 올바르지 않습니다.'); };
  const object = x => x && typeof x === 'object' && !Array.isArray(x);
  const shape = (x, keys) => object(x) && Object.keys(x).every(k => keys.includes(k)) && keys.every(k => Object.hasOwn(x, k));
  const str = (x, max = 12000) => { if (typeof x !== 'string' || x.length > max) fail(); return x; };
  if (!shape(value, Object.keys(defaults)) || !shape(value.controls, Object.keys(defaults.controls)) || typeof value.notebookExpanded !== 'boolean') fail();
  const out = { ...defaults };
  for (const key of Object.keys(defaults).filter(k => typeof defaults[k] === 'string')) out[key] = str(value[key]);
  if (!['', 'draft', 'design-1', 'design-2', 'design-3'].includes(out.selectedDesignSource)) fail();
  out.controls = Object.fromEntries(Object.keys(defaults.controls).map(k => [k, str(value.controls[k])]));
  out.notebookExpanded = value.notebookExpanded;
  if (out.selectedDesignSource) {
    if (!Number.isFinite(Date.parse(str(value.selectedAt, 40))) || !shape(value.primerSnapshot, ['forward', 'reverse', 'bindings'])) fail();
    const p = value.primerSnapshot;
    if (![p.forward, p.reverse].every(x => typeof x === 'string' && /^[ACGT]{1,100}$/.test(x))) fail();
    if (!object(p.bindings) || Object.keys(p.bindings).some(k => !['F', 'R'].includes(k))) fail();
    out.primerSnapshot = { forward: p.forward, reverse: p.reverse, ...bindingFields(p.bindings) };
  } else if (value.primerSnapshot !== null || value.selectedAt !== '') fail();
  return out;
}
export const sourceLabel = source => source === 'draft' ? '선택 당시의 현재 초안' : source ? `설계 ${source.slice(-1)}` : '선택 없음';
export const roughPosition = value => value === null || value === undefined ? '기록 없음' : value < 29 ? '결실 구간 왼쪽 / 대략적 위치' : value <= 48 ? '결실 구간 안쪽 / 대략적 위치' : '결실 구간 오른쪽 / 대략적 위치';
export function inspectFinal(design, fixture) {
  if (!fixture || !design) return null;
  const model = inspectDesign(design, fixture);
  return model.result && !model.errors.length && model.selectedProduct ? model : null;
}
export function finalCandidates(state, fixture) {
  return [...state.designs.map(d => ({ source: `design-${d.id}`, label: `설계 ${d.id}`, design: d })), { source: 'draft', label: '현재 초안', design: state.draft }]
    .filter(c => inspectFinal(c.design, fixture));
}
export function selectFinalDesign(state, source, fixture) {
  const selected = finalCandidates(state, fixture).find(c => c.source === source);
  if (!selected) throw new Error('03에서 유효한 결합 위치와 primer pair를 먼저 설계하세요.');
  const model = inspectFinal(selected.design, fixture);
  // Store only the pair and its coordinates, not derived calculations or reflections.
  state.finalReview.selectedDesignSource = source;
  state.finalReview.selectedAt = new Date().toISOString();
  state.finalReview.primerSnapshot = { forward: model.primers.F.stats.sequence, reverse: model.primers.R.stats.sequence,
    bindings: cloneBindings(Object.fromEntries(['F', 'R'].map(k => [k, model.primers[k].binding]))) };
}
export const productText = result => result ? Object.entries(result).map(([name, r]) => [name, r.products.length ? [...new Set(r.products.map(p => p.length))].join(', ') + ' bp' : '예상 산물 없음']) : [['A / B / C', '계산 자료 없음']];
export function primerReviewRows(model) {
  if (!model) return [['최종 pair 계산', '최종 설계 선택 또는 교육 데이터 확인이 필요합니다.']];
  const f = model.primers.F.stats, r = model.primers.R.stats;
  const ends = s => { const e = endFeatures(s); return `${e.sequence} / 마지막 염기 ${e.lastBase} / 말단 5 nt의 G+C ${e.gcCount}개`; };
  const comp = (a, b = a) => { const c = complementarity(a, b); return `가장 긴 연속 상보 구간 ${c.longest?.length || 0} nt / 3′ 말단 포함 ${c.threePrime?.length || 0} nt`; };
  return [
    ['Length / GC', `Forward ${f.length} nt / ${Math.round(f.gcPercent)}%\nReverse ${r.length} nt / ${Math.round(r.gcPercent)}%`],
    ['간이 Tm / 2(A+T)+4(G+C)', `Forward ${f.simpleTm} °C / Reverse ${r.simpleTm} °C / 차이 ${Math.abs(f.simpleTm - r.simpleTm)} °C`],
    ['3′ 말단', `Forward ${ends(f.sequence)}\nReverse ${ends(r.sequence)}`],
    ['자기 상보성', `Forward ${comp(f.sequence)}\nReverse ${comp(r.sequence)}`],
    ['F/R 상보성', comp(f.sequence, r.sequence)]
  ];
}
export const externalStatus = s => ({ unperformed: '미실시', performed: '외부 검색 실행 / 결과 기록 중', recorded: '외부 검색 결과 기록 완료' })[s.status];
export function externalSummary(s, pair) {
  const rows = [['상태', externalStatus(s)]];
  if (s.status === 'unperformed') return { rows, claim: '', note: '외부 검토 미실시. 검색 기록은 현재 설계의 근거에 포함되지 않습니다.' };
  const candidate = s.candidates[s.selectedCandidate === 'A' ? 0 : s.selectedCandidate === 'B' ? 1 : -1];
  rows.push(['Route', ROUTES[s.searchRoute]], ...[['date', '검색 날짜'], ['organism', 'Organism'], ['database', 'Database'], ['target', 'Target'], ['specificity', 'Specificity 설정']].map(([k, label]) => [label, s.conditions[k]]),
    ['선택한 candidate', s.selectedCandidate ? `Candidate ${s.selectedCandidate}` : '기록 없음']);
  if (candidate) rows.push(['Reported unintended target', UNINTENDED[candidate.unintended]], ['보고된 accession / 관찰', candidate.observations], ['선택한 candidate의 F / 5′→3′', candidate.forward], ['선택한 candidate의 R / 5′→3′', candidate.reverse], ['외부 도구에 보고된 product', candidate.product ? `${candidate.product} bp` : '기록 없음']);
  rows.push(['후보 선택 근거', s.selectionReason], ['학생의 주장 범위 reflection', s.claimReflection], ['학생의 실제 실험 한계 reflection', s.wetLabReflection]);
  const norm = x => x.replace(/\s/g, '').toUpperCase();
  const same = pair && candidate && ['forward', 'reverse'].every(k => norm(pair[k]) === norm(candidate[k]));
  return { rows, claim: claimScope(s), note: `${same ? '선택한 외부 후보의 F/R이 최종 pair와 같습니다.' : '선택한 외부 후보와 최종 pair의 일치가 확인되지 않았습니다. 이 기록을 최종 pair의 검색 결과로 간주하지 마세요.'} 이 학습지는 외부 결과를 직접 가져오거나 검증하지 않습니다. ${CLAIM_LIMIT}` };
}
export const limitationRows = state => [
  ['실제 PCR에서 증폭', '미수행 / 이 학습지에서 실제 PCR을 수행하지 않았습니다.'],
  ['실제 gel의 단일 product', '미수행 / 05는 수업용 가상 자료입니다.'],
  ['Amplicon의 실제 서열 정체', '미확인'],
  ['실제 반응 조건의 primer-dimer / hairpin', '미검증 / 04의 단순 서열 상보성 계산은 실제 검증이 아닙니다.'],
  ['외부 데이터베이스 특이성 검토', `${externalStatus(state.externalSearch)} / 학생의 기록 상태이며 최종 pair의 특이성 판정이 아닙니다.`]
];
