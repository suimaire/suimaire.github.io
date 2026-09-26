import { EVIDENCE_CASES, EVIDENCE_LANES, evidenceBandPosition, caseSummary } from './pcr-evidence.mjs';
import { reviewDesign } from './pcr-review.mjs';

const node = (tag, text = '', className = '') => {
  const el = document.createElement(tag); el.textContent = text; el.className = className; return el;
};

// Mount before the worksheet collects answer fields, so all answers share its v1 path.
export function mountEvidence(root) {
  const $ = id => root.querySelector(`#${id}`);
  for (const c of EVIDENCE_CASES) {
    const button = node('button', c.label); button.type = 'button'; button.id = `evidence-${c.id}`;
    button.dataset.evidenceCase = c.id; button.setAttribute('aria-controls', 'evidence-case-view');
    if (c.additional) $('evidence-additional').append(button);
    else $('evidence-navigation').insertBefore(button, $('evidence-additional'));
    const panel = node('section', '', 'pcr-evidence-answers'); panel.id = `evidence-${c.id}-answers`;
    const heading = node('h3', `${c.label} / 관찰에서 해석으로`); heading.id = `${panel.id}-title`;
    panel.setAttribute('aria-labelledby', heading.id); panel.append(heading);
    panel.append(node('p', `수업용 가상 관찰 / 사례의 예상 산물 ${c.expectedSize} bp. ${caseSummary(c)}`, 'pcr-evidence-print-summary'));
    const selected = node('p', '', 'pcr-evidence-print-selection'); selected.dataset.caseSelection = c.id; panel.append(selected);
    for (const question of c.questions) panel.append(node('p', question, 'pcr-evidence-question'));
    for (const [key, labelText, hint] of [
      ['observation', '관찰', '선택한 lane에서 보이는 사실을 짧게 적으세요.'],
      ['interpretation', '해석', '이 관찰과 일치하는 설명은 무엇인가요?'],
      ['uncertainty', '아직 확정할 수 없는 점', '남은 불확실성과 다음 확인을 기록하세요.']
    ]) {
      const id = `evidence-${c.id}-${key}`, label = node('label', labelText); label.htmlFor = id;
      const input = node(key === 'observation' ? 'input' : 'textarea'); input.id = id; input.dataset.answer = '';
      if (key !== 'observation') input.rows = 2;
      const help = node('p', hint, 'pcr-small'); help.id = `${id}-help`; input.setAttribute('aria-describedby', help.id);
      panel.append(label, help, input);
    }
    const explanation = node('details', '', 'pcr-explanation');
    explanation.append(node('summary', `${c.label} 해석의 범위 살펴보기`), node('p', c.explanation)); panel.append(explanation);
    $('evidence-case-answers').append(panel);
  }
  const grid = $('evidence-gel');
  const axis = node('div', '', 'pcr-gel-axis'); axis.setAttribute('aria-hidden', 'true');
  axis.append(node('span', 'bp', 'pcr-gel-lane-name'));
  const scale = node('span', '', 'pcr-gel-track');
  for (const bp of EVIDENCE_CASES[0].lanes.marker.bands) {
    const label = node('span', String(bp), 'pcr-gel-tick'); label.style.top = `${evidenceBandPosition(bp)}%`; scale.append(label);
  }
  axis.append(scale); grid.append(axis);
  for (const [id, lane] of Object.entries(EVIDENCE_LANES)) {
    const button = node('button', '', 'pcr-gel-lane'); button.type = 'button'; button.id = `evidence-lane-${id}`; button.dataset.lane = id;
    button.setAttribute('aria-controls', 'evidence-selected-lane');
    const name = node('span', lane.label, 'pcr-gel-lane-name'), track = node('span', '', 'pcr-gel-track'); track.setAttribute('aria-hidden', 'true');
    track.append(node('span', '', 'pcr-gel-well'), node('span', '', 'pcr-gel-bands'));
    button.append(name, track, node('span', '', 'pcr-gel-selection')); grid.append(button);
  }
  $('evidence-navigation').hidden = false; $('evidence-case-view').hidden = false;
}

export function initializeEvidence(root, getState, getFixture) {
  const $ = id => root.querySelector(`#${id}`);
  let legacyOpenBeforePrint = null;
  function renderPrediction() {
    const model = reviewDesign(getState(), getFixture()), target = $('gel-results'); target.replaceChildren();
    $('evidence-design-source').textContent = model ? `${model.designId === null ? '03 현재 초안' : `04에서 선택한 저장 설계 ${model.designId}`} / 완전 일치 모형` : '';
    if (!model) {
      target.append(node('p', '03에서 유효한 primer pair를 먼저 설계하면 여기에 계산된 예상 산물이 표시됩니다.'));
      if (!getFixture()) target.append(node('p', '서열 데이터가 준비되지 않았습니다. 아래 고정 가상 관찰 자료는 사용할 수 있습니다.', 'pcr-small'));
      return;
    }
    const list = node('dl', '', 'pcr-evidence-products');
    for (const [source, result] of Object.entries(model.result)) {
      const item = node('div'); item.append(node('dt', source), node('dd', result.products.length ? [...new Set(result.products.map(p => p.length))].map(bp => `${bp} bp`).join(', ') : '예상 산물 없음')); list.append(item);
    }
    target.append(list);
  }
  function renderSelection() {
    const view = getState().evidence, c = EVIDENCE_CASES.find(c => c.id === view.activeCase), selected = view.selectedLanes[c.id];
    for (const [id, lane] of Object.entries(c.lanes)) {
      const button = $(`evidence-lane-${id}`), active = id === selected;
      button.setAttribute('aria-pressed', String(active)); button.setAttribute('aria-label', `${lane.label} lane. ${lane.observation}. ${lane.description}`);
      button.querySelector('.pcr-gel-selection').textContent = active ? '선택됨' : '';
    }
    const target = $('evidence-selected-lane'); target.replaceChildren();
    if (selected) {
      const lane = c.lanes[selected];
      target.append(node('strong', `선택한 lane / ${lane.label}`), node('p', lane.description, 'pcr-small'), node('p', `관찰 / ${lane.observation}`));
    } else target.append(node('p', 'Lane을 선택해 관찰하세요. 마우스, 터치 또는 Tab과 Enter/Space를 사용할 수 있습니다.'));
    for (const c of EVIDENCE_CASES) {
      const lane = view.selectedLanes[c.id];
      root.querySelector(`[data-case-selection="${c.id}"]`).textContent = lane ? `선택한 lane: ${c.lanes[lane].label}` : '선택한 lane: 기록 없음';
    }
  }
  function render() {
    renderPrediction();
    const view = getState().evidence, c = EVIDENCE_CASES.find(c => c.id === view.activeCase);
    for (const item of EVIDENCE_CASES) {
      $(`evidence-${item.id}`).setAttribute('aria-pressed', String(item.id === c.id));
      $(`evidence-${item.id}-answers`).hidden = item.id !== c.id;
    }
    $('evidence-additional').open = c.additional;
    $('evidence-case-heading').textContent = `${c.label} / 수업용 가상 관찰`;
    $('evidence-case-expected').textContent = `이 고정 사례의 계산된 예상: ${c.expectedSize} bp 산물 하나. 학생 설계와 독립된 조건입니다.`;
    $('evidence-observation-summary').textContent = caseSummary(c);
    for (const [id, lane] of Object.entries(c.lanes)) {
      const target = $(`evidence-lane-${id}`).querySelector('.pcr-gel-bands'); target.replaceChildren();
      for (const bp of lane.bands) {
        const band = node('span', '', 'pcr-gel-band'); band.dataset.bp = bp; band.style.top = `${evidenceBandPosition(bp)}%`; target.append(band);
      }
    }
    renderSelection();
    $('evidence-legacy').hidden = !getState().answers['evidence-cases'];
  }
  function preparePrint() {
    if (legacyOpenBeforePrint === null) legacyOpenBeforePrint = $('evidence-legacy').open;
    if (!$('evidence-legacy').hidden) $('evidence-legacy').open = true;
  }
  function finishPrint() {
    if (legacyOpenBeforePrint !== null) $('evidence-legacy').open = legacyOpenBeforePrint;
    legacyOpenBeforePrint = null;
  }
  return { render, renderPrediction, renderSelection, preparePrint, finishPrint };
}
