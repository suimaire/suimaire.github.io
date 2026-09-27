import { analyzeAll } from './pcr-core.mjs';
import { STORAGE_KEY, DATA_VERSION, MAX_IMPORT_BYTES, emptyRecord, parseRecord, INTRO_CHOICE_KEYS } from './pcr-records.mjs';

import { initializeIntro } from './pcr-intro.mjs';
import { initializeWorkbench } from './pcr-workbench.mjs';
import { initializeReview } from './pcr-review-view.mjs';
import { mountEvidence, initializeEvidence } from './pcr-evidence-view.mjs';
import { cloneBindings, inspectDesign } from './pcr-design.mjs';
import { initializeExternal } from './pcr-external-view.mjs';
import { LEGACY_FINAL } from './pcr-final.mjs';
import { initializeFinalReview } from './pcr-final-view.mjs';

const root = document.querySelector('#pcr-worksheet');
if (root) initialize().catch(error => {
  root.querySelector('#save-status').textContent = `도구를 시작하지 못했습니다: ${error.message} 기본 문항을 읽고 브라우저 인쇄를 이용하세요.`;
});
async function initialize() {
  const $ = id => root.querySelector(`#${id}`);
  mountEvidence(root);
  const fields = [...root.querySelectorAll('[data-answer]')];
  const keys = [...fields.map(f => f.id), ...INTRO_CHOICE_KEYS, ...Object.keys(LEGACY_FINAL)];
  for (const field of fields) if (field.tagName !== 'SELECT') field.maxLength = 12000;
  let state = emptyRecord(), fixture, currentResult = null, storageCorrupt = false, printBlank = false;
  let finalReviewView;
  const defaults = Object.fromEntries(fields.map(f => [f.id, f.value]));
  const textNode = (tag, text, className) => {
    const el = document.createElement(tag); el.textContent = text;
    if (className) el.className = className;
    return el;
  };
  function status(id, message, error = false) { $(id).textContent = message; $(id).dataset.error = String(error); }
  function on(id, action) {
    $(id).addEventListener('click', async () => {
      try { await action(); } catch (error) { status('analysis-status', error.message, true); }
    });
  }
  function ready() { if (!fixture) throw new Error('교육 데이터를 불러오지 못했습니다. 연결을 확인한 뒤 새로고침하세요. 작성 기록은 내보낼 수 있습니다.'); }
  function save() {
    state.updatedAt = new Date().toISOString();
    renderFinal();
    if (storageCorrupt) return;
    try { localStorage.setItem(STORAGE_KEY, JSON.stringify(state)); status('save-status', '이 브라우저에 자동 저장했습니다. 중요한 기록은 JSON으로 내보내세요.'); }
    catch { status('save-status', '브라우저 저장이 차단되었거나 공간이 부족합니다. 현재 기록은 자동 저장되지 않습니다. 기록 내보내기로 보관하세요.', true); }
  }
  try {
    const raw = localStorage.getItem(STORAGE_KEY);
    if (raw) state = parseRecord(raw, keys);
    status('save-status', raw ? '이 브라우저의 이전 기록을 복원했습니다.' : '입력하면 이 브라우저에 자동 저장합니다.');
  } catch (error) {
    storageCorrupt = true;
    status('save-status', `이전 기록을 읽을 수 없습니다: ${error.message} 기존 저장값을 덮어쓰지 않습니다. 현재 작업은 기록 내보내기로 보관하세요.`, true);
  }
  const intro = initializeIntro(root, () => state, save);
  finalReviewView = initializeFinalReview(root, () => state, () => fixture, save);
  const review = initializeReview(root, () => state, () => fixture, save);
  const evidence = initializeEvidence(root, () => state, () => fixture);
  const external = initializeExternal(root, () => state, () => fixture, save, renderFinal);
  const workbench = initializeWorkbench(root, () => state, () => fixture, save, model => {
    currentResult = model.result ? { forward: state.draft.forward.replace(/\s/g, '').toUpperCase(), reverse: state.draft.reverse.replace(/\s/g, '').toUpperCase(), result: model.result } : null;
    evidence.renderPrediction();
    external.renderSource();
    review.render();
    renderFinal();
  });
  function reflectState() {
    for (const f of fields) f.value = state.answers[f.id] ?? defaults[f.id];
    $('primer-f').value = state.draft.forward; $('primer-r').value = state.draft.reverse;
    intro.render(); workbench.resetSelection(); invalidate(); renderDesigns(); renderFinal(); evidence.render(); external.render();
    finalReviewView.restore();
  }
  root.addEventListener('input', event => {
    const f = event.target;
    if (f.matches('[data-answer]')) { state.answers[f.id] = f.value; save(); renderFinal(); }
    if (f.id === 'primer-f' || f.id === 'primer-r') {
      state.draft = { ...state.draft, forward: $('primer-f').value, reverse: $('primer-r').value, bindings: cloneBindings(state.draft.bindings) };
      state.draft.bindings[f.id === 'primer-f' ? 'F' : 'R'] = null; workbench.cancelSelection(); invalidate(); save();
    }
  });
  root.addEventListener('change', event => {
    if (event.target.id === 'review-design') { evidence.renderPrediction(); external.renderSource(); }
  });
  root.addEventListener('click', event => {
    const caseButton = event.target.closest('[data-evidence-case]');
    const laneButton = event.target.closest('[data-lane]');
    if (caseButton) {
      state.evidence.activeCase = caseButton.dataset.evidenceCase;
      evidence.render(); save();
    } else if (laneButton) {
      state.evidence.selectedLanes[state.evidence.activeCase] = laneButton.dataset.lane;
      evidence.renderSelection(); save();
    }
  });
  function setDraft(forward, reverse, bindings) {
    state.draft = { forward, reverse, ...(bindings ? { bindings: cloneBindings(bindings) } : {}) };
    $('primer-f').value = forward; $('primer-r').value = reverse; workbench.resetSelection(); invalidate(); save();
  }
  function invalidate() {
    currentResult = null;
    if (fixture) workbench.refresh();
    else evidence.renderPrediction();
    review.render();
    external.renderSource();
  }
  const productLabel = p => `${p.source} ${p.start + 1}~${p.end}: ${p.length} bp`;
  const summarize = products => products.length ? products.map(productLabel).join('\n') : '이 완전 일치 모형에서 산물 미검출';
  function calculate() {
    ready(); workbench.refresh();
    if (!currentResult) throw new Error($('analysis-status').textContent || '먼저 F와 R을 선택하세요.');
    return currentResult;
  }
  on('analyze-design', calculate);
  on('save-design', () => {
    if (state.designs.length >= 3) return;
    const computed = currentResult || calculate();
    state.designs.push({ id: state.designs.length + 1, createdAt: new Date().toISOString(), forward: computed.forward, reverse: computed.reverse,
      prediction: $('design-prediction').value, reason: $('design-reason').value, unresolved: $('design-unresolved').value, ...(state.draft.bindings ? { bindings: cloneBindings(state.draft.bindings) } : {}) });
    save(); renderDesigns(); renderFinal(); review.render(); evidence.renderPrediction(); status('design-status', `설계 ${state.designs.length}을 별도로 보존했습니다. 이후 입력 변경은 이 설계를 덮어쓰지 않습니다.`);
  });
  function designSummary(design, withButton = false) {
    const block = textNode('div', '', 'pcr-saved');
    block.append(textNode('h3', `설계 ${design.id}`), textNode('p', `F 5′-${design.forward}-3′\nR 5′-${design.reverse}-3′`, 'pcr-sequence'));
    if (fixture) {
      try {
        const result = analyzeAll(fixture.templates, design.forward, design.reverse);
        const description = withButton ? Object.entries(result).map(([source, r]) => `${source} ${r.products.length ? `${[...new Set(r.products.map(p => p.length))].join(', ')} bp` : '예상 산물 없음'}`).join('\n') : Object.values(result).map(r => summarize(r.products)).join('\n');
        block.append(textNode('p', `완전 일치 모형의 예상 산물\n${description}`));
      }
      catch (error) { block.append(textNode('p', `계산 범위 오류: ${error.message}`)); }
    }
    if (withButton && fixture) {
      const info = inspectDesign(design, fixture).primers;
      block.append(textNode('p', ['F', 'R'].map(name => `${name} ${info[name].binding ? `${info[name].binding.start}–${info[name].binding.end}` : '서열 기반'}`).join(' / '), 'pcr-small'));
    }
    for (const [label, key] of [['계산 전 예상', 'prediction'], ['선택 또는 수정 이유', 'reason'], ['미확인', 'unresolved']]) if (!withButton || design[key]) block.append(textNode('p', `${label}: ${design[key] || '기록 없음'}`));
    if (withButton) { const b = textNode('button', `설계 ${design.id} 불러오기`); b.type = 'button'; b.addEventListener('click', () => { setDraft(design.forward, design.reverse, design.bindings); workbench.focus(); }); block.append(b); }
    return block;
  }
  function renderDesigns() {
    $('saved-designs').replaceChildren(...state.designs.map(d => designSummary(d, true)));
    $('save-design').textContent = state.designs.length === 3 ? '설계 1~3 보존됨' : `현재 설계 저장 / ${state.designs.length + 1}번`;
    $('save-design').disabled = state.designs.length === 3;
  }
  function renderFinal() {
    finalReviewView?.render();
  }
  on('export-record', () => {
    state.updatedAt = new Date().toISOString();
    const url = URL.createObjectURL(new Blob([JSON.stringify(state, null, 2)], { type: 'application/json' }));
    const a = textNode('a', ''); a.href = url; a.download = `pcr-notebook-${state.updatedAt.slice(0, 10)}.json`; root.append(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  });
  $('import-record').addEventListener('change', async event => {
    const file = event.target.files[0]; if (!file) return;
    try {
      if (file.size > MAX_IMPORT_BYTES) throw new Error('가져올 JSON은 1 MB 이하여야 합니다.');
      const imported = parseRecord(await file.text(), keys);
      for (const f of fields.filter(f => f.tagName === 'SELECT')) if (imported.answers[f.id] !== undefined && ![...f.options].some(o => o.value === imported.answers[f.id])) throw new Error(`${f.id} 선택값이 올바르지 않습니다.`);
      if (!window.confirm('현재 학습지 기록을 불러온 기록으로 바꿀까요? 보관이 필요하면 먼저 내보내세요.')) return;
      state = imported; storageCorrupt = false; reflectState(); save();
    } catch (error) { status('save-status', `불러오기 실패: ${error.message} 현재 기록은 유지했습니다.`, true); }
    finally { event.target.value = ''; }
  });
  on('reset-record', () => {
    if (!window.confirm('이 학습지의 답안과 설계 1~3을 모두 초기화할까요? 다른 학습지의 기록은 유지합니다.')) return;
    try { localStorage.removeItem(STORAGE_KEY); } catch { /* In-memory reset still works with blocked storage. */ }
    state = emptyRecord(); storageCorrupt = false; reflectState();
    $('direction-feedback').textContent = ''; workbench.resetSelection(); save();
  });
  function preparePrint() {
    root.dataset.printBlank = String(printBlank);
    root.querySelectorAll('.pcr-print-value').forEach(el => el.remove());
    for (const field of root.querySelectorAll('textarea, input:not([type=file]):not([type=radio]):not([type=range]), select')) {
      const value = field.id === 'review-length-choice' ? (field.value ? field.selectedOptions[0].textContent : '') : field.value;
      const mirror = textNode('div', printBlank ? '' : value, 'pcr-print-value'); mirror.dataset.multiline = String(field.tagName === 'TEXTAREA'); mirror.dataset.rows = field.rows || 1; field.after(mirror);
    }
    intro.preparePrint(printBlank, textNode);
    evidence.preparePrint();
    external.preparePrint(printBlank);
    finalReviewView.preparePrint(printBlank);
  }
  window.addEventListener('beforeprint', preparePrint);
  window.addEventListener('afterprint', () => { printBlank = false; delete root.dataset.printBlank; root.querySelectorAll('.pcr-print-value').forEach(el => el.remove()); evidence.finishPrint(); finalReviewView.finishPrint(); });
  for (const [id, blank] of [['print-filled', false], ['print-blank', true]]) on(id, () => { printBlank = blank; preparePrint(); window.print(); });
  const menus = [...root.querySelectorAll('.pcr-menu')];
  for (const menu of menus) {
    menu.addEventListener('toggle', () => { if (menu.open) menus.filter(other => other !== menu).forEach(other => { other.open = false; }); });
    menu.addEventListener('keydown', event => { if (event.key === 'Escape') { menu.open = false; menu.querySelector('summary').focus(); } });
  }
  document.addEventListener('click', event => menus.forEach(menu => { if (!menu.contains(event.target)) menu.open = false; }));
  const tocLinks = [...root.querySelectorAll('.pcr-sidebar nav a')];
  function markCurrentActivity() {
    const current = [...tocLinks].reverse().find(link => root.querySelector(link.hash).getBoundingClientRect().top <= 180) || tocLinks[0];
    for (const link of tocLinks) { if (link === current) link.setAttribute('aria-current', 'location'); else link.removeAttribute('aria-current'); }
  }
  let scrollScheduled = false;
  window.addEventListener('scroll', () => {
    if (scrollScheduled) return;
    scrollScheduled = true;
    requestAnimationFrame(() => { markCurrentActivity(); scrollScheduled = false; });
  }, { passive: true });
  markCurrentActivity();
  const wide = matchMedia('(min-width: 901px)'); $('worksheet-toc').open = wide.matches;
  wide.addEventListener('change', () => { $('worksheet-toc').open = wide.matches; });
  reflectState();
  try {
    const response = await fetch(root.dataset.fixtureUrl);
    if (!response.ok) throw new Error(`자료 응답 ${response.status}`);
    fixture = await response.json();
    if (fixture.id !== DATA_VERSION) throw new Error('교육 데이터 버전 불일치');
    const templateContainer = $('template-sequences'); templateContainer.replaceChildren();
    for (const [id, t] of Object.entries(fixture.templates)) templateContainer.append(textNode('h3', `${id} / ${t.length} nt / 상단 5′→3′`), textNode('p', t.sequence, 'pcr-sequence'));
    workbench.refresh(); renderDesigns(); renderFinal(); root.dataset.ready = 'true';
  } catch (error) { fixture = null; review.render(); evidence.renderPrediction(); status('analysis-status', `교육 데이터 로딩 실패: ${error.message} 답안 기록, JSON과 인쇄는 계속 사용할 수 있습니다.`, true); }
}
