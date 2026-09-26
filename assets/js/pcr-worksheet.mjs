import { analyzeAll, complement, reverseComplement, normalizeSequence, primerStats, primerFromRange, combineProducts, gelPosition } from './pcr-core.mjs';
import { STORAGE_KEY, DATA_VERSION, MAX_IMPORT_BYTES, emptyRecord, parseRecord } from './pcr-records.mjs';

const root = document.querySelector('#pcr-worksheet');
if (root) initialize().catch(error => {
  root.querySelector('#save-status').textContent = `도구를 시작하지 못했습니다: ${error.message} 기본 문항을 읽고 브라우저 인쇄를 이용하세요.`;
});
async function initialize() {
  const $ = id => root.querySelector(`#${id}`);
  const fields = [...root.querySelectorAll('[data-answer]')];
  const keys = fields.map(f => f.id);
  for (const field of fields) if (field.tagName !== 'SELECT') field.maxLength = 12000;
  let state = emptyRecord(), fixture, currentResult = null, storageCorrupt = false, anchor = null, columns = 40, printBlank = false;
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
  function reflectState() {
    for (const f of fields) f.value = state.answers[f.id] ?? defaults[f.id];
    $('primer-f').value = state.draft.forward; $('primer-r').value = state.draft.reverse;
    renderCycle(); invalidate(); renderDesigns(); renderFinal();
  }
  root.addEventListener('input', event => {
    const f = event.target;
    if (f.matches('[data-answer]')) { state.answers[f.id] = f.value; save(); renderFinal(); }
    if (f.id === 'primer-f' || f.id === 'primer-r') {
      state.draft = { forward: $('primer-f').value, reverse: $('primer-r').value }; save(); invalidate();
    }
  });
  const cycles = [
    ['1. 변성 / 예시 95 °C', '가열로 상보적인 두 가닥이 분리됩니다. DNA의 당-인산 골격을 잘라내는 단계가 아닙니다.'],
    ['2. 결합 / 예시 55 °C', '온도를 낮추면 프라이머가 상보적인 주형 구간에 결합할 수 있습니다. F와 R의 3′ 끝이 서로 안쪽을 향하는 기본 예입니다.'],
    ['3. 신장 / 예시 72 °C', 'DNA 중합효소가 dNTP를 사용해 새 가닥의 3′ 끝에만 염기를 추가합니다. 새 가닥은 5′→3′로 자랍니다. 첫 주기에는 반대쪽 프라이머 위치에서 자동으로 잘리지 않고 더 길게 이어집니다.']
  ];
  function renderCycle() {
    $('cycle-title').textContent = cycles[state.cycle][0]; $('cycle-description').textContent = cycles[state.cycle][1];
    $('cycle-primers').setAttribute('visibility', state.cycle >= 1 ? 'visible' : 'hidden');
    $('cycle-primer-ends').setAttribute('visibility', state.cycle === 1 ? 'visible' : 'hidden');
    $('cycle-extension').setAttribute('visibility', state.cycle === 2 ? 'visible' : 'hidden');
  }
  on('cycle-next', () => { state.cycle = (state.cycle + 1) % 3; renderCycle(); save(); });
  on('cycle-reset', () => { state.cycle = 0; renderCycle(); save(); });
  on('check-direction', () => {
    try {
      const c = normalizeSequence($('direction-complement').value), r = normalizeSequence($('direction-reverse').value);
      status('direction-feedback', `상보 서열: ${c === complement('AGTCCGTA') ? '일치합니다' : '다시 확인하세요'}. 주문용 역상보 서열: ${r === reverseComplement('AGTCCGTA') ? '일치합니다' : '다시 확인하세요'}. 3′와 5′ 방향을 함께 읽으세요.`);
    } catch (error) { status('direction-feedback', error.message, true); }
  });
  function setDraft(forward, reverse) {
    state.draft = { forward, reverse }; $('primer-f').value = forward; $('primer-r').value = reverse; save(); invalidate();
  }
  function invalidate() {
    currentResult = null;
    $('analysis-results').replaceChildren(textNode('p', '현재 입력의 계산 결과가 없습니다. 산물 계산을 눌러 확인하세요.'));
    $('gel-results').replaceChildren(textNode('p', '활동 03에서 현재 서열로 산물을 계산하세요.'));
    $('map-bindings').replaceChildren();
    status('analysis-status', '서열을 수정한 경우 다시 계산하세요. 저장된 설계는 그대로 보존됩니다.');
    renderTm();
  }
  function selectedRange(start, end) {
    $('range-start').value = Math.min(start, end); $('range-end').value = Math.max(start, end);
    for (const cell of root.querySelectorAll('.pcr-base')) cell.setAttribute('aria-pressed', String(+cell.dataset.base >= Math.min(start, end) && +cell.dataset.base <= Math.max(start, end)));
    status('range-status', `${Math.min(start, end)}~${Math.max(start, end)}번 선택. 신장 방향을 확인한 뒤 주문 서열로 넣으세요.`);
  }
  function selectEndpoint(index) {
    if (anchor === null) { anchor = index; selectedRange(index, index); status('range-status', `${index}번에서 선택 시작. 끝 좌표를 누르세요.`); }
    else { selectedRange(anchor, index); anchor = null; }
  }
  function renderSequence() {
    if (!fixture) return;
    const grid = $('sequence-grid'), width = grid.clientWidth;
    if (!width) return;
    columns = Math.max(8, Math.min(40, Math.floor(width / 24)));
    const sequence = fixture.templates.A.sequence, lower = complement(sequence);
    const focused = document.activeElement?.dataset.base;
    grid.replaceChildren();
    for (let offset = 0; offset < sequence.length; offset += columns) {
      const row = textNode('div', '', 'pcr-sequence-row'), coords = textNode('div', '', 'pcr-row-coordinates');
      coords.append(textNode('span', `${offset + 1}`), textNode('span', `${Math.min(offset + columns, sequence.length)}`));
      const bases = textNode('div', '', 'pcr-bases'); bases.style.setProperty('--pcr-columns', columns);
      for (let i = offset; i < Math.min(offset + columns, sequence.length); i++) {
        const cell = textNode('button', '', 'pcr-base'); cell.type = 'button'; cell.dataset.base = i + 1;
        cell.tabIndex = i + 1 === +(focused || 1) ? 0 : -1;
        cell.setAttribute('aria-label', `A ${i + 1}번 상단 ${sequence[i]} 하단 ${lower[i]}`);
        cell.setAttribute('aria-pressed', String(i + 1 >= +$('range-start').value && i + 1 <= +$('range-end').value));
        cell.append(textNode('span', sequence[i]), textNode('span', lower[i])); bases.append(cell);
      }
      row.append(coords, bases); grid.append(row);
    }
    if (focused) grid.querySelector(`[data-base="${focused}"]`)?.focus({ preventScroll: true });
  }
  let dragStart = null, suppressClick = false;
  $('sequence-grid').addEventListener('pointerdown', event => {
    const cell = event.target.closest('[data-base]');
    if (cell && event.pointerType === 'mouse') dragStart = +cell.dataset.base;
  });
  root.addEventListener('pointerup', event => {
    const cell = event.target.closest('[data-base]');
    if (dragStart !== null && cell && dragStart !== +cell.dataset.base) { selectedRange(dragStart, +cell.dataset.base); anchor = null; suppressClick = true; }
    dragStart = null;
  });
  root.addEventListener('pointercancel', () => { dragStart = null; });
  $('sequence-grid').addEventListener('click', event => {
    if (suppressClick) { suppressClick = false; return; }
    const cell = event.target.closest('[data-base]'); if (cell) selectEndpoint(+cell.dataset.base);
  });
  $('sequence-grid').addEventListener('keydown', event => {
    const cell = event.target.closest('[data-base]'); if (!cell) return;
    const move = { ArrowLeft: -1, ArrowRight: 1, ArrowUp: -columns, ArrowDown: columns }[event.key];
    if (move !== undefined || event.key === 'Home' || event.key === 'End') {
      event.preventDefault();
      const index = event.key === 'Home' ? 1 : event.key === 'End' ? 420 : Math.max(1, Math.min(420, +cell.dataset.base + move));
      const next = $('sequence-grid').querySelector(`[data-base="${index}"]`);
      cell.tabIndex = -1; next.tabIndex = 0; next.focus();
    }
  });
  $('sequence-details').addEventListener('toggle', renderSequence);
  let sequenceWidth = 0;
  new ResizeObserver(entries => {
    const width = Math.floor(entries[0].contentRect.width);
    if (width !== sequenceWidth) { sequenceWidth = width; renderSequence(); }
  }).observe($('sequence-grid'));
  on('apply-range', () => {
    try {
      ready(); const start = Number($('range-start').value), end = Number($('range-end').value);
      const p = primerFromRange(fixture.templates.A.sequence, start, end, $('range-direction').value);
      const label = $('range-primer').value;
      setDraft(label === 'F' ? p : state.draft.forward, label === 'R' ? p : state.draft.reverse);
      selectedRange(start, end); status('range-status', `${label} 주문 서열 5′-${p}-3′를 넣었습니다. 산물 계산을 눌러 실제 결합 방향을 확인하세요.`);
    } catch (error) { status('range-status', error.message, true); }
  });
  function table(headers, rows, caption) {
    const el = textNode('table', ''), head = document.createElement('thead'), body = document.createElement('tbody'), hr = document.createElement('tr');
    if (caption) el.append(textNode('caption', caption));
    for (const h of headers) { const th = textNode('th', h); th.scope = 'col'; hr.append(th); }
    head.append(hr);
    for (const row of rows) { const tr = document.createElement('tr'); for (const value of row) tr.append(textNode('td', String(value))); body.append(tr); }
    el.append(head, body); return el;
  }
  const productLabel = p => `${p.source} ${p.start + 1}~${p.end}: ${p.length} bp`;
  const summarize = products => products.length ? products.map(productLabel).join('\n') : '이 완전 일치 모형에서 산물 미검출';
  function statsTable(forward, reverse) {
    return table(['프라이머', '길이', 'GC 비율', '3′ 쪽 5 nt', 'Tm'], [forward, reverse].map((p, i) => {
      const s = primerStats(p); return [i ? 'R' : 'F', `${s.length} nt`, `${s.gcPercent.toFixed(1)}%`, s.threePrime, '정밀 계산 전'];
    }), '현재 입력의 기본 조성');
  }
  function calculate() {
    ready(); const forward = normalizeSequence(state.draft.forward), reverse = normalizeSequence(state.draft.reverse);
    const result = analyzeAll(fixture.templates, forward, reverse);
    currentResult = { forward, reverse, result };
    const container = $('analysis-results'); container.replaceChildren(statsTable(forward, reverse));
    for (const [source, data] of Object.entries(result)) {
      container.append(textNode('h3', `인공 DNA ${source}`), textNode('p', summarize(data.products), 'pcr-record-text'));
      if (data.overlappingPairs) container.append(textNode('p', `겹치는 결합 쌍 ${data.overlappingPairs}개는 지원 범위 밖이므로 산물에서 제외했습니다.`));
      const details = document.createElement('details'); details.append(textNode('summary', `${source}의 모든 결합 위치와 산물 서열`));
      if (data.hits.length) details.append(table(['이름', '결합 구간', '주형 가닥', '신장'], data.hits.map(h => [h.primer, `${h.start + 1}~${h.end}`, h.strand === 'lower' ? '하단' : '상단', h.direction === 'right' ? '→ 오른쪽' : '← 왼쪽'])));
      else details.append(textNode('p', '완전 일치 결합 위치 미검출'));
      for (const p of data.products) {
        details.append(textNode('p', `${productLabel(p)} / 가능한 쌍: ${p.bindings.map(b => `${b.right.primer}→ ←${b.left.primer}`).join(', ')}`), textNode('p', `5′-${p.sequence}-3′`, 'pcr-sequence'));
      }
      container.append(details);
    }
    drawMap(result.A.hits); renderGel(result); renderTm();
    status('analysis-status', '현재 F와 R의 모든 완전 일치 결합을 계산했습니다. 입력 이름으로 방향을 강제하지 않습니다.');
    return currentResult;
  }
  on('analyze-design', calculate);
  function svgNode(tag, attrs, text) {
    const node = document.createElementNS('http://www.w3.org/2000/svg', tag);
    for (const [key, val] of Object.entries(attrs)) node.setAttribute(key, val);
    if (text !== undefined) node.textContent = text; return node;
  }
  function drawMap(hits) {
    const layer = $('map-bindings'); layer.replaceChildren();
    for (const h of hits) {
      const x = 30 + h.start / 420 * 660, width = (h.end - h.start) / 420 * 660, y = h.primer === 'F' ? 108 : 126;
      layer.append(svgNode('rect', { x, y, width, height: 5, fill: h.primer === 'F' ? '#0d5751' : 'white', stroke: '#0d5751' }));
      layer.append(svgNode('text', { x, y: y - 3, 'font-size': 12 }, `${h.primer}${h.direction === 'right' ? '→' : '←'}`));
    }
  }
  function renderGel(result) {
    const lanes = { A: result.A.products, B: result.B.products, C: result.C.products, 'A+C': combineProducts(result.A.products, result.C.products), 'B+C': combineProducts(result.B.products, result.C.products) };
    const svg = svgNode('svg', { viewBox: '0 0 700 310', role: 'img', 'aria-label': '완전 일치 모형의 예상 산물 전기영동 개념도. 실제 실험 결과가 아닙니다.' });
    svg.append(svgNode('title', {}, '완전 일치 모형의 예상 산물'));
    for (const size of [5000, 1000, 500, 200, 100, 20]) { const y = gelPosition(size); svg.append(svgNode('text', { x: 0, y: y + 5 }, `${size}`), svgNode('path', { d: `M55 ${y}H680`, stroke: '#ddd' })); }
    svg.append(svgNode('text', { x: 5, y: 286 }, 'bp'));
    Object.entries(lanes).forEach(([name, products], i) => {
      const x = 115 + i * 115; svg.append(svgNode('text', { x: x - 17, y: 304 }, name), svgNode('rect', { x: x - 23, y: 7, width: 46, height: 6, fill: 'none', stroke: '#555' }));
      for (const size of new Set(products.map(p => p.length))) if (size >= 20) svg.append(svgNode('path', { d: `M${x - 25} ${gelPosition(size)}h50`, stroke: '#222', 'stroke-width': 4 }));
    });
    $('gel-results').replaceChildren(svg, table(['시료', '출처를 보존한 산물 목록'], Object.entries(lanes).map(([name, list]) => [name, summarize(list)]), '완전 일치 모형의 예상 산물'));
    if (Object.values(lanes).flat().some(p => p.length < 20)) $('gel-results').append(textNode('p', '20 bp 미만 산물은 그림 범위 밖입니다. 위 목록에는 보존되어 있습니다.'));
  }
  function renderTm() {
    try { $('tm-values').replaceChildren(table(['현재 입력', '조성 기반 간이 Tm'], ['forward', 'reverse'].map((key, i) => [i ? 'R' : 'F', `${primerStats(state.draft[key]).simpleTm} °C`]))); }
    catch { $('tm-values').textContent = '유효한 현재 F와 R 서열을 모두 입력하면 간이값이 표시됩니다.'; }
  }
  function compareCandidates(sources) {
    ready(); const rows = Object.entries(fixture.candidatePairs).map(([id, pair]) => {
      const result = analyzeAll(fixture.templates, pair.forward, pair.reverse);
      return [id, ...sources.map(source => summarize(result[source].products))];
    });
    $('candidate-results').replaceChildren(table(['후보', ...sources], rows, '완전 일치 모형의 예상 산물'));
    for (const [id, pair] of Object.entries(fixture.candidatePairs)) {
      const b = textNode('button', `${id}을 편집 입력으로 가져오기`); b.type = 'button';
      b.addEventListener('click', () => { setDraft(pair.forward, pair.reverse); status('design-status', `${id}을 편집 입력에 넣었습니다. 저장 설계는 변경하지 않았습니다.`); $('primer-f').focus(); });
      $('candidate-results').append(b);
    }
  }
  on('compare-ab', () => compareCandidates(['A', 'B'])); on('compare-abc', () => compareCandidates(['A', 'B', 'C']));
  on('save-design', () => {
    if (state.designs.length >= 3) return;
    const computed = currentResult || calculate();
    state.designs.push({ id: state.designs.length + 1, createdAt: new Date().toISOString(), forward: computed.forward, reverse: computed.reverse,
      prediction: $('design-prediction').value, reason: $('design-reason').value, unresolved: $('design-unresolved').value });
    save(); renderDesigns(); renderFinal(); status('design-status', `설계 ${state.designs.length}을 별도로 보존했습니다. 이후 입력 변경은 이 설계를 덮어쓰지 않습니다.`);
  });
  function designSummary(design, withButton = false) {
    const block = textNode('div', '', 'pcr-saved');
    block.append(textNode('h3', `설계 ${design.id}`), textNode('p', `F 5′-${design.forward}-3′\nR 5′-${design.reverse}-3′`, 'pcr-sequence'));
    if (fixture) {
      try { const result = analyzeAll(fixture.templates, design.forward, design.reverse); block.append(textNode('p', `완전 일치 모형의 예상 산물\n${Object.values(result).map(r => summarize(r.products)).join('\n')}`)); }
      catch (error) { block.append(textNode('p', `계산 범위 오류: ${error.message}`)); }
    }
    for (const [label, key] of [['계산 전 예상', 'prediction'], ['선택 또는 수정 이유', 'reason'], ['미확인', 'unresolved']]) block.append(textNode('p', `${label}: ${design[key] || '기록 없음'}`));
    if (withButton) { const b = textNode('button', `설계 ${design.id}을 편집 입력으로 복사`); b.type = 'button'; b.addEventListener('click', () => { setDraft(design.forward, design.reverse); $('primer-f').focus(); }); block.append(b); }
    return block;
  }
  function renderDesigns() {
    $('saved-designs').replaceChildren(...state.designs.map(d => designSummary(d, true)));
    $('save-design').textContent = state.designs.length === 3 ? '설계 1~3 보존됨' : `설계 ${state.designs.length + 1} 저장`;
    $('save-design').disabled = state.designs.length === 3;
  }
  function renderFinal() {
    const target = $('final-comparison');
    target.replaceChildren(textNode('h3', '처음 예측'), textNode('p', `배치: ${state.answers['first-placement'] || '기록 없음'}\n음성 해석: ${state.answers['first-negative'] || '기록 없음'}`, 'pcr-record-text'));
    target.append(...state.designs.map(d => designSummary(d)));
    target.append(textNode('h3', '관찰과 판단의 출처'), textNode('p', `내부 계산을 본 뒤의 판단: ${state.answers['candidate-judgment'] || '기록 없음'}\n수업용 가상 자료 해석: ${state.answers['evidence-cases'] || '기록 없음'}\n학생이 기록한 외부 검색: ${state.answers['external-status'] || '미실시'}\n외부 후보 비교: ${state.answers['external-comparison'] || '기록 없음'}`, 'pcr-record-text'));
  }
  async function copy(text) {
    try { await navigator.clipboard.writeText(text); status('copy-status', '출처와 함께 클립보드에 복사했습니다.'); }
    catch { $('copy-status').replaceChildren(textNode('span', '클립보드를 사용할 수 없습니다. 아래 텍스트를 선택해 복사하세요.'), textNode('pre', text, 'pcr-sequence')); }
  }
  on('copy-synthetic', () => { const f = normalizeSequence(state.draft.forward), r = normalizeSequence(state.draft.reverse); return copy(`학습용 인공 서열 설계. 자연 유전자 또는 NCBI 검증 결과가 아님.\nF 5′→3′: ${f}\nR 5′→3′: ${r}`); });
  on('copy-external', () => copy(`학생이 외부 결과에서 기록한 값. 학습지 자체 검증 아님.\n실시 여부: ${$('external-status').value}\naccession.version: ${$('external-accession').value}\nF 5′→3′: ${$('external-f').value}\nR 5′→3′: ${$('external-r').value}`));
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
    $('candidate-results').replaceChildren(); $('direction-feedback').textContent = ''; $('range-start').value = ''; $('range-end').value = ''; anchor = null;
    renderSequence(); save();
  });
  function preparePrint() {
    root.dataset.printBlank = String(printBlank);
    root.querySelectorAll('.pcr-print-value').forEach(el => el.remove());
    for (const field of root.querySelectorAll('textarea, input:not([type=file]), select')) {
      const mirror = textNode('div', printBlank ? '' : field.value, 'pcr-print-value'); mirror.dataset.multiline = String(field.tagName === 'TEXTAREA'); field.after(mirror);
    }
  }
  window.addEventListener('beforeprint', preparePrint);
  window.addEventListener('afterprint', () => { printBlank = false; delete root.dataset.printBlank; root.querySelectorAll('.pcr-print-value').forEach(el => el.remove()); });
  for (const [id, blank] of [['print-filled', false], ['print-blank', true]]) on(id, () => { printBlank = blank; preparePrint(); window.print(); });
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
    $('candidate-sequences').append(table(['후보', 'F 주문 서열 5′→3′', 'R 주문 서열 5′→3′'], Object.entries(fixture.candidatePairs).map(([id, p]) => [id, p.forward, p.reverse]), '고정 교육 후보'));
    renderSequence(); renderDesigns(); renderFinal(); root.dataset.ready = 'true';
  } catch (error) { fixture = null; status('analysis-status', `교육 데이터 로딩 실패: ${error.message} 답안 기록, JSON과 인쇄는 계속 사용할 수 있습니다.`, true); }
}
