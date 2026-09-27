import { inspectDesign } from './pcr-design.mjs';
import { LEGACY_FINAL, sourceLabel, roughPosition, inspectFinal, finalCandidates, selectFinalDesign, productText, primerReviewRows, externalSummary, limitationRows } from './pcr-final.mjs';

const GEL_NOTE = '05의 gel 자료는 수업용 가상 자료이며 이 primer pair의 실제 실험 결과가 아닙니다.';
const ASSESSMENT_NOTE = '실제 PCR 성공, 실제 amplicon identity, 검색하지 않은 sequence에 대한 specificity는 별도의 검증이 필요합니다.';
const REVIEW_ANSWERS = [['길이와 GC에 관한 학생의 이유', 'review-length-reason'], ['3′ 말단에 관한 학생의 관찰', 'review-end-observation'], ['제한된 A/B/C 비교에 관한 학생의 판단', 'candidate-judgment'], ['실제 실험 전에 아직 확인할 것 / 04 마지막 답안', 'review-unresolved']];
const CONTROL_LABELS = { positive: 'Positive control', negative: 'NTC 또는 negative control', additional: '필요하다면 추가 control', interpretationLimit: '대조군이 예상대로 나오지 않을 때 내릴 수 없는 결론' };
const node = (tag, text = '', className = '') => { const n = document.createElement(tag); n.textContent = text; if (className) n.className = className; return n; };
const paragraph = text => node('p', text || '기록 없음', 'pcr-record-text');
function entries(rows) {
  const dl = node('dl', '', 'final-entries');
  for (const [label, value] of rows) { const row = node('div', '', 'final-entry'); row.append(node('dt', label), node('dd', value || '기록 없음', label.includes('5′→3′') ? 'pcr-sequence final-sequence' : '')); dl.append(row); }
  return dl;
}
const section = (title, ...children) => { const s = node('section', '', 'final-note-section'); s.append(node('h4', title), ...children); return s; };
function initial(state, heading = '처음 예측 / 활동 00') {
  const s = section(heading), p = state.initialPrimerPrediction;
  // The track retains percentages. The only bp label is the known deletion interval.
  const diagram = node('div', '', 'final-map'); diagram.setAttribute('aria-hidden', 'true');
  diagram.append(node('span', 'A DNA / 축약 구조', 'final-map-label'), node('span', '', 'final-map-deletion'), node('span', '121~200 결실 구간', 'final-map-deletion-label'));
  for (const [name, key] of [['F', 'forward'], ['R', 'reverse']]) if (p[key] !== null) {
    const mark = node('span', `${name} 예상`, `final-map-marker final-map-${name}`); mark.style.left = `${p[key]}%`; diagram.append(mark);
  }
  s.append(diagram, entries([['Forward', roughPosition(p.forward)], ['Reverse', roughPosition(p.reverse)], ['처음 설명', state.answers['first-placement']]]), node('p', '00에 저장된 대략적 예측입니다. 정밀 bp 좌표를 뜻하지 않습니다.', 'pcr-small'));
  return s;
}
function history(state, fixture) {
  if (!state.designs.length) return paragraph('저장된 설계 없음');
  const list = node('ol', '', 'final-history'); list.setAttribute('aria-label', '저장된 설계 순서');
  for (const design of state.designs) {
    const li = node('li'), model = fixture ? inspectDesign(design, fixture) : null;
    li.append(node('h4', `설계 ${design.id}`));
    const coords = ['F', 'R'].map(k => { const b = model?.primers[k].binding || design.bindings?.[k]; return `${k} ${b ? `${b.start}–${b.end}` : '좌표 기록 없음'}`; }).join(' / ');
    li.append(paragraph(coords), entries(productText(model?.result)), entries([['선택 또는 수정 이유', design.reason], ['아직 확인하지 못한 점', design.unresolved]]));
    if (model?.errors.length) li.append(paragraph('저장된 원문은 보존되어 있습니다. 현재 계산 범위에서 유효한 pair인지 03에서 확인하세요.'));
    list.append(li);
  }
  return list;
}
function pairContent(state, fixture, compact = false) {
  const f = state.finalReview, pair = f.primerSnapshot, model = inspectFinal(pair, fixture);
  const body = node('div');
  if (!pair) { body.append(paragraph('최종 설계 선택 전 / 03에서 설계한 후 사용할 pair를 직접 선택하세요.')); return body; }
  body.append(node('p', `${sourceLabel(f.selectedDesignSource)} / 선택 시점의 pair 보존`, 'pcr-small'));
  if (!model) body.append(paragraph('교육 데이터 또는 pair의 결합 위치를 확인할 수 없어 계산값을 표시하지 않습니다. 선택한 주문 서열은 보존합니다.'));
  for (const [k, name, key] of [['F', 'Forward', 'forward'], ['R', 'Reverse', 'reverse']]) {
    const p = model?.primers[k], b = pair.bindings[k], s = section(name);
    s.append(paragraph(`A 좌표 ${b ? `${b.start}–${b.end}` : '기록 없음'}${p ? ` / ${p.stats.length} nt / GC ${Math.round(p.stats.gcPercent)}%` : ''}`));
    if (!compact) { const sequence = node('div', `5′–${pair[key]}–3′`, 'pcr-sequence final-sequence'); sequence.setAttribute('role', 'group'); sequence.setAttribute('aria-label', `${name} 주문 서열 5′에서 3′: ${pair[key]}`); s.append(sequence); }
    body.append(s);
  }
  if (compact) body.append(node('p', '실제 binding sequence를 확인한 내부 완전 일치 모형', 'pcr-small'), entries(productText(model?.result)));
  return body;
}
function reviewContent(state, fixture) {
  const body = node('div'), source = state.review.designId === null ? '03의 현재 draft' : `저장 설계 ${state.review.designId}`;
  body.append(paragraph('계산 출처: 07에서 선택해 보존한 최종 pair에 04와 동일한 계산 함수를 적용했습니다.'), entries(primerReviewRows(inspectFinal(state.finalReview.primerSnapshot, fixture))),
    node('p', '간이 Tm과 단순 서열 상보성 관찰입니다. 정밀 열역학 분석이나 실제 primer-dimer / hairpin 검증이 아닙니다.', 'pcr-small'),
    section('04에서 학생이 남긴 기록', paragraph(`현재 04 검토 대상: ${source}. 아래 답안은 04의 현재 기록이며, 답안을 쓸 당시의 pair는 별도로 저장되지 않았습니다.`), entries(REVIEW_ANSWERS.map(([label, key]) => [label, state.answers[key]]))));
  return body;
}
function evidenceContent(state) {
  return section('실험 증거를 해석할 때 주의할 점', paragraph(GEL_NOTE), entries([
    ['예상 크기의 band 하나만으로 sequence identity를 확정할 수 없는 이유', state.answers['evidence-identity']],
    ['대조군이 필요한 이유', state.answers['evidence-controls']]
  ]));
}
function externalContent(state) {
  const summary = externalSummary(state.externalSearch, state.finalReview.primerSnapshot);
  return section('학생이 기록한 외부 검색 결과', entries(summary.rows), section('Claim scope / 기록된 조건에서의 주장 범위', paragraph(summary.claim || '기록 없음')), node('p', summary.note, 'pcr-small'));
}
function legacyContent(state) {
  const rows = Object.entries(LEGACY_FINAL).filter(([key]) => state.answers[key]).map(([key, label]) => [label, state.answers[key]]);
  return rows.length ? entries(rows) : null;
}

export function initializeFinalReview(root, getState, getFixture, save) {
  const $ = id => root.querySelector(`#${id}`), fields = [...root.querySelectorAll('[data-final]')];
  let choiceKey = '', renderedKey = '', isPrinting = false;
  const getValue = (review, path) => path.split('.').reduce((v, k) => v[k], review);
  for (const field of fields) {
    field.maxLength = 12000;
    field.addEventListener('input', () => {
      const path = field.dataset.final.split('.'), record = getState().finalReview;
      if (path.length === 1) record[path[0]] = field.value; else record[path[0]][path[1]] = field.value;
      save();
    });
  }
  function choose(source) { selectFinalDesign(getState(), source, getFixture()); save(); render(); $('final-selection-status').textContent = `${sourceLabel(source)}의 primer pair를 선택 시점의 값으로 보존했습니다.`; }
  $('final-design-options').addEventListener('change', e => { if (e.target.matches('[name=final-design]')) choose(e.target.value); });
  $('final-refresh-draft').addEventListener('click', () => choose('draft'));
  $('final-notebook-toggle').addEventListener('click', () => { getState().finalReview.notebookExpanded = !getState().finalReview.notebookExpanded; save(); render(); });
  function renderChoices(state, fixture) {
    const options = finalCandidates(state, fixture), f = state.finalReview;
    const key = JSON.stringify(options.map(c => c.source));
    if (choiceKey !== key) {
      choiceKey = key;
      $('final-design-options').replaceChildren(...options.map(c => {
        const label = node('label', '', 'final-design-option'), radio = node('input');
        radio.type = 'radio'; radio.name = 'final-design'; radio.value = c.source; radio.id = `final-select-${c.source}`;
        radio.checked = c.source === f.selectedDesignSource; label.htmlFor = radio.id;
        label.append(radio, node('span', c.label)); return label;
      }));
    }
    for (const radio of $('final-design-options').querySelectorAll('input')) radio.checked = radio.value === f.selectedDesignSource;
    const current = options.find(c => c.source === f.selectedDesignSource);
    const changed = f.primerSnapshot && (!current || ['forward', 'reverse'].some(k => current.design[k].replace(/\s/g, '').toUpperCase() !== f.primerSnapshot[k]) || (current.design.bindings && JSON.stringify(current.design.bindings) !== JSON.stringify(f.primerSnapshot.bindings)));
    $('final-selection-note').textContent = !options.length ? '선택 가능한 설계 없음. 03에서 유효한 primer pair를 설계하세요.' : '사용할 설계를 직접 선택하세요. 선택한 pair는 이후 03을 수정해도 유지됩니다.';
    $('final-snapshot-note').textContent = f.primerSnapshot ? `${sourceLabel(f.selectedDesignSource)} / ${f.selectedAt} 선택. ${changed ? '원본이 변경되었거나 현재 후보에서 사용할 수 없습니다. 최종 pair는 선택 당시의 값으로 유지합니다.' : '최종 pair는 선택 시점의 서열과 좌표입니다.'}` : '아직 최종 설계를 선택하지 않았습니다.';
    $('final-refresh-draft').hidden = !(f.selectedDesignSource === 'draft' && changed && options.some(c => c.source === 'draft'));
  }
  function notebook(blank = false) {
    const state = getState(), f = state.finalReview, fixture = getFixture(), report = $('final-notebook');
    report.replaceChildren(node('p', 'PCR과 프라이머 디자인', 'pcr-section-label'), node('h3', '최종 연구 노트'));
    const past = () => paragraph('이전 활동의 기록이 여기에 요약됩니다.');
    const answer = (label, value, rows = 3) => { const s = section(label, paragraph(blank ? '' : value)); s.classList.add('final-response'); if (blank) { s.lastChild.textContent = ''; s.lastChild.classList.add('final-writing-space'); s.lastChild.style.minHeight = `${rows * 1.7}em`; } return s; };
    report.append(
      section('1. 연구 질문', answer('이 primer pair를 이용해 무엇을 구별하거나 확인하려는가?', f.researchQuestion)),
      section('2. 처음 예측', blank ? past() : initial(state)),
      section('3. 설계 변화', blank ? past() : history(state, fixture)),
      section('4. 최종 primer pair', blank ? past() : pairContent(state, fixture), answer('최종 선택 근거', f.finalRationale, 4)),
      section('5. 예상 PCR products', blank ? past() : entries(productText(inspectFinal(f.primerSnapshot, fixture)?.result)), node('p', '제공된 A/B/C 서열 안의 내부 완전 일치 모형입니다.', 'pcr-small')),
      section('6. Primer 특성 검토', blank ? past() : reviewContent(state, fixture)),
      section('7. 실험 증거 해석 원칙', blank ? past() : evidenceContent(state)),
      section('8. 외부 database 검토', blank ? past() : externalContent(state)),
      section('9. 대조군 계획', ...Object.entries(CONTROL_LABELS).map(([k, label]) => answer(label, f.controls[k], k === 'interpretationLimit' ? 3 : 2))),
      section('10. 아직 확인하지 못한 것', blank ? past() : entries(limitationRows(state)), answer('기타 확인할 것', f.otherLimitations, 2)),
      section('11. 처음 예측에서 수정된 점', answer('무엇을 바꾸었으며 왜 바꾸었는가?', f.revisionReflection, 4)),
      section('12. 현재 근거의 범위에서의 판단', answer('현재 확보한 계산, 관찰, 외부 검색 기록의 범위에서 이 설계를 어떻게 평가할 수 있는가?', f.finalAssessment, 4), node('p', ASSESSMENT_NOTE, 'pcr-small'))
    );
    if (!blank) { const legacy = legacyContent(state); if (legacy) report.append(section('이전 07 기록 / 보존한 원문', legacy)); }
  }
  function render() {
    if (isPrinting) return;
    const state = getState(), fixture = getFixture(), f = state.finalReview;
    for (const field of fields) if (field !== document.activeElement && field.value !== getValue(f, field.dataset.final)) field.value = getValue(f, field.dataset.final);
    renderChoices(state, fixture);
    const key = JSON.stringify([state.initialPrimerPrediction, state.answers, state.designs, state.review, state.externalSearch, f, Boolean(fixture)]);
    if (key !== renderedKey) {
      renderedKey = key;
      const comparison = $('final-comparison'); comparison.replaceChildren(initial(state), section('최종 설계 / 활동 07', pairContent(state, fixture, true)));
      $('final-history-content').replaceChildren(history(state, fixture));
      $('final-primer-content').replaceChildren(pairContent(state, fixture));
      $('final-products-content').replaceChildren(entries(productText(inspectFinal(f.primerSnapshot, fixture)?.result)));
      $('final-review-content').replaceChildren(reviewContent(state, fixture));
      $('final-evidence-content').replaceChildren(evidenceContent(state));
      $('final-external-content').replaceChildren(externalContent(state));
      $('final-limitations-content').replaceChildren(entries(limitationRows(state)));
      const legacy = legacyContent(state); $('final-legacy').hidden = !legacy; $('final-legacy-content').replaceChildren(...(legacy ? [legacy] : []));
      notebook();
    }
    $('final-notebook').hidden = !f.notebookExpanded;
    $('final-notebook-toggle').setAttribute('aria-expanded', String(f.notebookExpanded));
    $('final-notebook-toggle').textContent = f.notebookExpanded ? '최종 연구 노트 접기' : '최종 연구 노트 보기';
  }
  root.dataset.finalReady = 'true';
  return { render, restore() { $('final-selection-status').textContent = ''; for (const field of fields) field.value = getValue(getState().finalReview, field.dataset.final); renderedKey = ''; render(); }, preparePrint(blank) { isPrinting = true; notebook(blank); }, finishPrint() { isPrinting = false; renderedKey = ''; render(); } };
}
