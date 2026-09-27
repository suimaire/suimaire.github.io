import { ROUTES, SEARCH_STATUS, UNINTENDED, CLAIM_LIMIT, currentPrimer, emptyCandidate, validDate, validPrimer, validProduct, validTm, tmDifference, completionIssues, claimScope, comparisonRecorded, planText } from './pcr-external.mjs';

const node = (tag, text = '', className = '') => { const el = document.createElement(tag); el.textContent = text; if (className) el.className = className; return el; };
const get = (value, path) => path.split('.').reduce((v, key) => v?.[key], value);
const set = (value, path, text) => { const keys = path.split('.'); const key = keys.pop(); get(value, keys.join('.'))[key] = text; };

export function initializeExternal(root, getState, getFixture, save, onChange) {
  const $ = id => root.querySelector(`#${id}`), state = () => getState().externalSearch;
  const fields = [...root.querySelectorAll('[data-external]')];
  fields.forEach(f => { if (f.tagName !== 'SELECT') f.maxLength = f.dataset.external === 'design.target' ? 60000 : 12000; });
  function assign(path, value) { if (path.includes('.')) set(state(), path, value); else state()[path] = value; }
  function changed() { save(); renderDerived(); onChange(); }
  function source() {
    const s = state();
    return s.route === 'mine' ? currentPrimer(getState(), getFixture()) : s.route === 'paper' ? { origin: s.paper.source, forward: s.paper.forward, reverse: s.paper.reverse } : null;
  }
  function renderSource() {
    const primer = currentPrimer(getState(), getFixture());
    $('ext-primer-origin').textContent = primer?.origin || '';
    $('ext-my-pair').hidden = !primer; $('ext-no-primer').hidden = !!primer;
    $('ext-my-forward').textContent = primer?.forward || ''; $('ext-my-reverse').textContent = primer?.reverse || '';
  }
  function feedback(f) {
    const path = f.dataset.external, value = f.value;
    let message = '';
    if (value && /\.(forward|reverse)$/.test(path) && !validPrimer(value)) message = '염기 서열만 기록하세요. 입력한 초안은 보존됩니다.';
    if (value && (/^candidates\.\d\.product$/.test(path) || /^design\.(min|max)$/.test(path)) && !validProduct(value)) message = 'bp는 양의 정수로 기록하세요. 초안은 보존됩니다.';
    if (value && /\.tm[FR]$/.test(path) && !validTm(value)) message = 'Reported Tm은 숫자로 기록하세요. 초안은 보존됩니다.';
    if (value && path === 'conditions.date' && !validDate(value)) message = '실제 날짜를 YYYY-MM-DD 형식으로 기록하세요.';
    if (path === 'design.max' && validProduct(value) && validProduct(state().design.min) && Number(value) < Number(state().design.min)) message = 'max는 min 이상이어야 합니다. 초안은 보존됩니다.';
    f.setAttribute('aria-invalid', String(!!message)); $(`${f.id}-help`).textContent = message;
  }
  function comparison() {
    const s = state(), target = $('ext-comparison-table'); target.replaceChildren();
    $('ext-comparison-note').textContent = comparisonRecorded(s) ? '두 후보의 외부 기록을 나란히 비교합니다. 선택은 학생이 근거와 함께 기록합니다.' : '후보 비교를 기록하지 않음';
    if (!comparisonRecorded(s)) return;
    target.append(comparisonTable(s));
  }
  function comparisonTable(s) {
    const table = node('table'); table.append(node('caption', '학생이 기록한 Candidate A / Candidate B'));
    const head = node('thead'), tr = node('tr');
    for (const title of ['항목', 'Candidate A', 'Candidate B']) { const th = node('th', title); th.scope = 'col'; tr.append(th); } head.append(tr); table.append(head);
    const body = node('tbody');
    for (const [label, read] of [ ['Product size', c => `${c.product} bp`], ['Reported F Tm', c => c.tmF ? `${c.tmF} °C` : '미기록'], ['Reported R Tm', c => c.tmR ? `${c.tmR} °C` : '미기록'], ['F/R Tm difference', tmDifference], ['Unintended target', c => `${UNINTENDED[c.unintended]}\n${c.observations}`], ['Position / 기타', c => c.other || '미기록'] ]) {
      const row = node('tr'), th = node('th', label); th.scope = 'row'; row.append(th, ...s.candidates.map(c => node('td', read(c)))); body.append(row);
    }
    table.append(body); return table;
  }
  function renderDerived() {
    const s = state(); renderSource();
    for (const route of Object.keys(ROUTES)) {
      const selected = route === s.route, button = $(`ext-route-${route}`);
      button.setAttribute('aria-pressed', String(selected)); button.querySelector('[data-route-selected]').textContent = selected ? '선택됨' : '';
      $(`ext-route-input-${route}`).hidden = !selected;
    }
    $('ext-status-text').textContent = `외부 검색 / ${SEARCH_STATUS[s.status]}`;
    $('ext-after-search').hidden = s.status === 'unperformed';
    $('ext-search-context').textContent = `아래 결과를 기록한 경로: ${ROUTES[s.searchRoute] || '미기록'}${s.sourcePrimer.origin ? ` / ${s.sourcePrimer.origin}` : ''}. 준비 화면을 바꿔도 이 검색 기록은 유지됩니다.`;
    $('ext-candidate-1').hidden = s.candidates.length < 2; $('ext-add-candidate').hidden = s.candidates.length === 2;
    $('ext-selectedCandidate').querySelector('option[value=B]').disabled = s.candidates.length < 2;
    $('ext-claim-text').textContent = claimScope(s) || (s.status === 'unperformed' ? '외부 검색 미실시' : '검색 날짜, organism, 데이터베이스, 후보와 unintended target 상태를 기록하면 검색 범위를 담은 문장이 나타납니다.');
    $('ext-legacy').hidden = !Object.entries(getState().answers).some(([key, value]) => key.startsWith('external-') && value);
    comparison();
    fields.forEach(feedback);
  }
  function render() {
    fields.forEach(f => { f.value = get(state(), f.dataset.external) ?? ''; });
    $('ext-copy-status').textContent = ''; $('ext-copy-fallback').hidden = true; $('ext-completion-help').textContent = '';
    renderDerived();
  }
  root.addEventListener('input', event => {
    const f = event.target;
    if (!f.matches('[data-external]')) return;
    assign(f.dataset.external, f.value);
    // Edited evidence requires another explicit completion, but never erases the draft.
    if (state().status === 'recorded' && /^(conditions|candidates|selectedCandidate)/.test(f.dataset.external)) state().status = 'performed';
    $('ext-completion-help').textContent = ''; changed();
  });
  for (const button of root.querySelectorAll('[data-route]')) button.addEventListener('click', () => { state().route = button.dataset.route; changed(); });
  // Native buttons support Tab/Enter/Space; arrows additionally cycle this compact selector.
  $('ext-routes').addEventListener('keydown', event => {
    const buttons = [...$('ext-routes').querySelectorAll('button')], index = buttons.indexOf(event.target);
    if (index < 0 || !['ArrowLeft', 'ArrowRight', 'ArrowUp', 'ArrowDown', 'Home', 'End'].includes(event.key)) return;
    event.preventDefault(); const next = event.key === 'Home' ? 0 : event.key === 'End' ? 2 : (index + (['ArrowLeft', 'ArrowUp'].includes(event.key) ? 2 : 1)) % 3;
    buttons[next].click(); buttons[next].focus();
  });
  $('ext-to-plan').addEventListener('click', () => { $('ext-plan').open = true; $('ext-plan-purpose').focus(); });
  $('ext-performed').addEventListener('click', () => {
    const s = state(); s.status = 'performed'; s.searchRoute = s.route;
    s.sourcePrimer = source() || { origin: ROUTES[s.route], forward: '', reverse: '' };
    // Actual conditions are manually transcribed, never inferred from a plan or link click.
    $('ext-conditions').open = true; changed(); $('ext-conditions-date').focus();
  });
  $('ext-unperformed').addEventListener('click', () => { state().status = 'unperformed'; changed(); });
  $('ext-add-candidate').addEventListener('click', () => { if (state().candidates.length < 2) state().candidates.push(emptyCandidate()); if (state().status === 'recorded') state().status = 'performed'; render(); changed(); $('ext-candidates-1-forward').focus(); });
  $('ext-complete').addEventListener('click', () => {
    const issues = completionIssues(state());
    if (issues.length) { $('ext-completion-help').textContent = `아직 기록할 항목: ${issues.join(', ')}. 입력한 초안은 저장되어 있습니다.`; return; }
    state().status = 'recorded'; changed(); $('ext-completion-help').textContent = '결과 기록 완료. 학생의 기록 상태이며 외부 결과를 자동 확인한 표시는 아닙니다.';
  });
  for (const button of root.querySelectorAll('[data-copy]')) button.addEventListener('click', async () => {
    const s = state(), primer = source(), kind = button.dataset.copy;
    const value = kind === 'plan' ? planText(s) : kind === 'target' ? (s.route === 'paper' ? s.paper.target : s.design.target) : kind === 'pair' ? (primer ? `Forward / 5′→3′\n${primer.forward}\nReverse / 5′→3′\n${primer.reverse}` : '') : primer?.[kind] || '';
    $('ext-copy-fallback').hidden = true;
    if (!value) { $('ext-copy-status').textContent = '복사할 내용을 먼저 입력하세요.'; return; }
    try { await navigator.clipboard.writeText(value); $('ext-copy-status').textContent = '복사됨'; }
    catch { $('ext-copy-status').textContent = '아래 텍스트를 선택해 직접 복사하세요.'; $('ext-copy-fallback').textContent = value; $('ext-copy-fallback').hidden = false; }
    // Place feedback next to the triggering tool even when the plan is far below the route.
    button.after($('ext-copy-status'), $('ext-copy-fallback'));
  });
  function preparePrint(blank) {
    const s = state(), target = $('ext-print'); target.replaceChildren();
    const add = (label, value = '') => {
      const wrap = node('div', '', 'ext-print-entry');
      const sequence = /Forward|Reverse|^Target$|내 primer/.test(label) ? ' pcr-sequence' : '';
      wrap.append(node('h4', label), node('p', blank ? '\u00a0' : value || '미기록', `pcr-record-text${sequence}`)); target.append(wrap);
    };
    add('선택 경로', ROUTES[s.route]); add('외부 검색', SEARCH_STATUS[s.status]);
    if (!blank && s.searchRoute) add('결과 기록의 경로', ROUTES[s.searchRoute]);
    const sourceLabels = { source: '출처', species: 'Species', purpose: '연구 목적', forward: 'Forward / 5′→3′', reverse: 'Reverse / 5′→3′', target: 'Target' };
    if (!blank && s.route === 'paper') for (const [key, label] of Object.entries(sourceLabels)) add(label, s.paper[key]);
    if (!blank && s.route === 'new') for (const [key, label] of Object.entries({ type: '입력 종류', target: 'Target', organism: 'Organism', purpose: 'PCR 목적', min: 'Product min bp', max: 'Product max bp' })) add(label, s.design[key]);
    if (!blank && s.route === 'mine') { const p = source(); add('내 primer / 현재 준비', p ? `${p.origin}\nF ${p.forward}\nR ${p.reverse}` : '유효한 primer pair 없음'); }
    for (const [key, label] of Object.entries({ purpose: '검색 목적', organism: 'Target organism 계획', target: 'Intended target 계획', database: '검색 데이터베이스 계획', notes: '검색 범위 선택 이유' })) add(label, s.plan[key]);
    if (blank || s.status !== 'unperformed') {
      target.append(node('h3', '실제 검색 조건'));
      for (const f of fields.filter(f => f.dataset.external.startsWith('conditions.'))) add(root.querySelector(`label[for="${f.id}"]`).textContent, f.value);
      target.append(node('h3', 'Primer-BLAST에서 보고된 결과 / 학생이 기록한 값'));
      const candidates = blank ? [emptyCandidate(), emptyCandidate()] : s.candidates;
      candidates.forEach((c, i) => { target.append(node('h4', `Candidate ${'AB'[i]}`)); for (const [key, label] of Object.entries({ forward: 'Forward / 5′→3′', reverse: 'Reverse / 5′→3′', product: 'Reported product size / bp', tmF: 'Reported F Tm / °C', tmR: 'Reported R Tm / °C', unintended: 'Unintended PCR target', observations: '주요 off-target 또는 관찰', other: 'Position / 기타' })) add(label, key === 'unintended' ? UNINTENDED[c[key]] : c[key]); });
      if (!blank && comparisonRecorded(s)) target.append(comparisonTable(s)); else add('후보 비교', '후보 비교를 기록하지 않음');
      add('선택한 후보', s.selectedCandidate ? `Candidate ${s.selectedCandidate}` : '미선택'); add('선택 근거', s.selectionReason);
      add('주장 범위', claimScope(s) || '주장 범위 기록에 필요한 정보 미완성');
    } else target.append(node('p', '외부 검색 미실시. 보관된 결과 초안은 현재 검색의 증거로 출력하지 않습니다.'));
    target.append(node('p', CLAIM_LIMIT, 'pcr-small'));
    add('Unintended target이 보고되지 않았다는 결과를 어떤 범위까지 주장할 수 있는가?', s.claimReflection);
    add('실제 PCR에서는 데이터베이스 검색 결과만으로 무엇을 아직 알 수 없는가?', s.wetLabReflection);
    if (!blank) {
      const old = Object.entries(getState().answers).filter(([key, value]) => key.startsWith('external-') && value);
      if (old.length) { target.append(node('h3', '이전 활동 06 기록')); for (const [key, value] of old) add(root.querySelector(`label[for="${key}"]`).textContent, value); }
    }
  }
  return { render, renderSource, preparePrint };
}
