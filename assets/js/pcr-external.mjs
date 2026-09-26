import { reviewDesign } from './pcr-review.mjs';

export const ROUTES = { mine: '내가 설계한 primer 검토', paper: '논문의 primer 다시 검토', new: '새로운 primer 후보 설계' };
export const SEARCH_STATUS = { unperformed: '미실시', performed: '외부 검색 실행', recorded: '결과 기록 완료' };
export const UNINTENDED = { '': '기록 전', none: '보고되지 않음', reported: '보고됨', unclear: '결과를 해석하지 못함' };
export const CLAIM_LIMIT = '이 기록만으로 검색하지 않은 organism, 다른 database, 다른 annotation, 다른 search setting, 실제 wet-lab PCR 조건에 대해서까지 비표적 증폭이 없다고 결론 내릴 수는 없습니다.';
export const emptyCandidate = () => ({ forward: '', reverse: '', product: '', tmF: '', tmR: '', unintended: '', observations: '', other: '' });
export const emptyExternalSearch = () => ({
  route: 'mine', status: 'unperformed', searchRoute: '',
  paper: { source: '', species: '', purpose: '', forward: '', reverse: '', target: '' },
  design: { type: 'accession', target: '', organism: '', purpose: '일반 PCR', min: '', max: '' },
  plan: { purpose: '', organism: '', target: '', database: '', notes: '' },
  conditions: { date: '', organism: '', database: '', target: '', forward: '', reverse: '', product: '', specificity: '', other: '' },
  sourcePrimer: { origin: '', forward: '', reverse: '' },
  candidates: [emptyCandidate()], selectedCandidate: '', selectionReason: '', claimReflection: '', wetLabReflection: ''
});

// Keep unfinished text intact. Structural validation is separate from scientific/input feedback.
export function validateExternalSearch(value) {
  if (value === undefined) return emptyExternalSearch();
  const fail = () => { throw new Error('외부 검색 기록의 구조 또는 선택값이 올바르지 않습니다.'); };
  const shape = (raw, defaults) => {
    if (!raw || typeof raw !== 'object' || Array.isArray(raw) || Object.keys(raw).some(k => !Object.hasOwn(defaults, k))) return fail();
    const out = {};
    for (const [key, defaultValue] of Object.entries(defaults)) {
      const v = raw[key];
      if (Array.isArray(defaultValue)) {
        if (!Array.isArray(v) || v.length < 1 || v.length > 2) return fail();
        out[key] = v.map(c => shape(c, emptyCandidate()));
      } else if (typeof defaultValue === 'object') out[key] = shape(v, defaultValue);
      else { if (typeof v !== 'string' || v.length > (key === 'target' ? 60000 : 12000)) return fail(); out[key] = v; }
    }
    return out;
  };
  const out = shape(value, emptyExternalSearch());
  if (!Object.hasOwn(ROUTES, out.route) || !Object.hasOwn(SEARCH_STATUS, out.status) || !['', ...Object.keys(ROUTES)].includes(out.searchRoute) ||
      !['accession', 'fasta'].includes(out.design.type) || !['일반 PCR', '발현 분석', '기타'].includes(out.design.purpose) ||
      !['', 'A', ...(out.candidates.length > 1 ? ['B'] : [])].includes(out.selectedCandidate) || out.candidates.some(c => !Object.hasOwn(UNINTENDED, c.unintended))) return fail();
  if (out.status !== 'unperformed' && !out.searchRoute) return fail();
  if (out.status === 'recorded' && completionIssues(out).length) return fail();
  return out;
}
export function currentPrimer(state, fixture) {
  const model = reviewDesign(state, fixture);
  if (!model) return null;
  const design = model.designId === null ? state.draft : state.designs.find(d => d.id === model.designId);
  return { origin: model.designId === null ? '03 현재 초안' : `04에서 선택한 저장 설계 ${model.designId}`, forward: design.forward.replace(/\s/g, '').toUpperCase(), reverse: design.reverse.replace(/\s/g, '').toUpperCase() };
}
export const validDate = value => /^\d{4}-\d{2}-\d{2}$/.test(value) && Number.isFinite(Date.parse(value)) && new Date(value).toISOString().slice(0, 10) === value;
export const validPrimer = value => /^[ACGTRYSWKMBDHVN]+$/i.test(value.replace(/\s/g, ''));
export const validProduct = value => /^\d+$/.test(value) && Number(value) > 0 && Number.isSafeInteger(Number(value));
export const validTm = value => /^-?\d+(\.\d+)?$/.test(value) && Number.isFinite(Number(value));
export function tmDifference(candidate) {
  return validTm(candidate.tmF) && validTm(candidate.tmR) ? `${Number(Math.abs(Number(candidate.tmF) - Number(candidate.tmR)).toFixed(2))} °C` : '기록 전';
}
export const candidateRecorded = c => validPrimer(c.forward) && validPrimer(c.reverse) && validProduct(c.product) && c.unintended !== '';
export function completionIssues(s) {
  const issues = [];
  if (!validDate(s.conditions.date)) issues.push('검색 날짜');
  for (const key of ['organism', 'database', 'target', 'specificity']) if (!s.conditions[key].trim()) issues.push({ organism: 'Organism / 제한 없음 여부', database: 'Database', target: 'Target / template 또는 미제공', specificity: 'Specificity 주요 설정' }[key]);
  if (!s.selectedCandidate) issues.push('주장에 사용할 후보 선택');
  s.candidates.forEach((c, i) => {
    if (!candidateRecorded(c)) issues.push(`Candidate ${'AB'[i]}의 F/R, 양의 정수 product size, unintended target 상태`);
    if ([c.tmF, c.tmR].some(v => v && !validTm(v))) issues.push(`Candidate ${'AB'[i]}의 Reported Tm 숫자`);
  });
  for (const key of ['forward', 'reverse']) if (s.conditions[key] && !validPrimer(s.conditions[key])) issues.push('검색에 사용한 primer 서열');
  if (s.searchRoute !== 'new' && (!validPrimer(s.conditions.forward) || !validPrimer(s.conditions.reverse))) issues.push('기존 pair 검토에 사용한 Forward와 Reverse');
  return issues;
}
export function claimScope(s) {
  if (s.status === 'unperformed') return '';
  const c = s.candidates[s.selectedCandidate === 'A' ? 0 : s.selectedCandidate === 'B' ? 1 : -1];
  if (!validDate(s.conditions.date) || !s.conditions.organism.trim() || !s.conditions.database.trim() || !c?.unintended) return '';
  const conclusion = { none: 'unintended PCR target이 보고되지 않았다고', reported: 'unintended PCR target이 보고되었다고', unclear: 'unintended PCR target 결과를 해석하지 못했다고' }[c.unintended];
  return `${s.conditions.date}에 기록한 ${s.conditions.organism} / ${s.conditions.database} 검색 조건에서는 Candidate ${s.selectedCandidate}의 ${conclusion} 기록했습니다.`;
}
export const comparisonRecorded = s => s.candidates.length === 2 && s.candidates.every(candidateRecorded);
export const planText = s => `검색 목적: ${s.plan.purpose}\nTarget organism: ${s.plan.organism}\nIntended target: ${s.plan.target}\nDatabase 계획: ${s.plan.database}\n범위 선택 이유: ${s.plan.notes}`;
