import { inspectDesign, sequenceFromBinding, PRIMER_KEYS, cloneBindings } from './pcr-design.mjs';

export function initializeWorkbench(root, getState, getFixture, save, onComputed) {
  const $ = id => root.querySelector(`#${id}`);
  const node = (tag, text = '', className = '') => {
    const el = document.createElement(tag); el.textContent = text; el.className = className; return el;
  };
  const svgNode = (tag, attrs = {}, text = '') => {
    const el = document.createElementNS('http://www.w3.org/2000/svg', tag);
    for (const [key, value] of Object.entries(attrs)) el.setAttribute(key, value);
    el.textContent = text; return el;
  };
  const label = (parent, x, y, text, anchor = 'start', className = '') => parent.appendChild(svgNode('text', { x, y, 'text-anchor': anchor, class: className }, text));
  const path = (parent, d, className) => parent.appendChild(svgNode('path', { d, class: className }));
  const direction = (parent, from, to, y, name) => {
    const sign = to > from ? 1 : -1;
    path(parent, `M${from} ${y}H${to}M${to - sign * 6} ${y - 5}L${to} ${y}L${to - sign * 6} ${y + 5}`, `pcr-map-primer pcr-primer-${name}`);
  };
  let model = null, anchor = null, drag = null, suppressClick = false, columns = 12, focusPosition = 1;
  const view = () => getState().workbench;
  const size = () => matchMedia('(max-width: 520px)').matches ? 30 : 60;
  const selected = () => model?.primers[view().mode]?.binding;
  const boundsText = binding => binding ? `${binding.start || '?'}–${binding.end || '?'}` : '선택 전';
  function syncControls() {
    const binding = selected();
    $('range-primer').value = view().mode;
    $('range-start').value = binding?.start ?? ''; $('range-end').value = binding?.end ?? '';
    $('range-direction').value = binding?.direction ?? (view().mode === 'F' ? 'right' : 'left');
  }
  function setWindow(start, focus = false) {
    const length = getFixture()?.templates.A.sequence.length || 420;
    view().windowStart = Math.max(1, Math.min(length - size() + 1, start));
    save(); renderSequence();
    if (focus) $('sequence-grid').querySelector(`[data-base="${focusPosition}"]`)?.focus();
  }
  function chooseMode(name) {
    view().mode = name; anchor = null; drag = null;
    syncControls();
    const binding = selected();
    const start = Number(binding?.start);
    if (Number.isInteger(start) && start >= 1 && start <= (getFixture()?.templates.A.sequence.length || 420)) setWindow(Math.max(1, start - 5));
    refresh(); save();
  }
  root.querySelectorAll('[data-select-primer]').forEach(button => button.addEventListener('click', () => chooseMode(button.dataset.selectPrimer)));
  $('range-primer').addEventListener('change', event => chooseMode(event.target.value));
  function applyBinding(binding, keepControls = false) {
    if (!getFixture()) return;
    const state = getState(), key = PRIMER_KEYS[view().mode];
    // Preserve the other primer's explicit selection or unambiguous legacy binding.
    state.draft.bindings = cloneBindings(state.draft.bindings);
    state.draft.bindings[view().mode] = binding;
    try { state.draft[key] = sequenceFromBinding(getFixture().templates.A.sequence, binding); }
    catch { state.draft[key] = ''; } // Invalid coordinates must never reuse a previous sequence/result.
    anchor = null;
    refresh(!keepControls); save();
  }
  const applyCoordinates = () => applyBinding({ start: $('range-start').value, end: $('range-end').value, direction: $('range-direction').value }, true);
  for (const id of ['range-start', 'range-end']) $(id).addEventListener('input', applyCoordinates);
  $('range-direction').addEventListener('change', applyCoordinates);
  $('apply-range').addEventListener('click', () => {
    applyCoordinates();
    const start = Number($('range-start').value);
    if (start >= 1 && start <= (getFixture()?.templates.A.sequence.length || 420)) setWindow(Math.max(1, start - 5));
  });
  $('show-initial-prediction').addEventListener('click', () => { view().showPrediction = !view().showPrediction; renderMap(); save(); });
  root.addEventListener('input', event => { if (event.target.id.startsWith('prediction-') && view().showPrediction) renderMap(); });
  $('prediction-track').addEventListener('click', () => { if (view().showPrediction) renderMap(); });
  $('sequence-prev').addEventListener('click', () => setWindow(view().windowStart - size()));
  $('sequence-next').addEventListener('click', () => setWindow(view().windowStart + size()));
  $('sequence-window-slider').addEventListener('input', event => setWindow(Number(event.target.value)));

  function paintSelection(start, end) {
    for (const cell of $('sequence-grid').querySelectorAll('[data-base]')) cell.setAttribute('aria-pressed', String(Number(cell.dataset.base) >= Math.min(start, end) && Number(cell.dataset.base) <= Math.max(start, end)));
  }
  function commitSelection(start, end) {
    applyBinding({ start: String(Math.min(start, end)), end: String(Math.max(start, end)), direction: $('range-direction').value });
  }
  function endpoint(index) {
    focusPosition = index;
    if (anchor === null) { anchor = index; paintSelection(index, index); $('range-status').textContent = `${index}번에서 선택 시작. 끝 염기를 선택하세요. Escape로 취소합니다.`; }
    else commitSelection(anchor, index);
  }
  const grid = $('sequence-grid');
  grid.addEventListener('pointerdown', event => {
    suppressClick = false;
    const cell = event.target.closest('[data-base]');
    if (cell && event.isPrimary && event.button === 0 && event.pointerType !== 'touch') {
      drag = { start: Number(cell.dataset.base), end: Number(cell.dataset.base), id: event.pointerId, moved: false };
      suppressClick = false;
      grid.setPointerCapture(event.pointerId);
    }
  });
  grid.addEventListener('pointermove', event => {
    if (!drag) return;
    const cell = document.elementFromPoint(event.clientX, event.clientY)?.closest('[data-base]');
    if (!cell || !grid.contains(cell)) return;
    drag.end = Number(cell.dataset.base); drag.moved ||= drag.end !== drag.start;
    if (drag.moved) { paintSelection(drag.start, drag.end); $('range-status').textContent = `${Math.min(drag.start, drag.end)}–${Math.max(drag.start, drag.end)} 선택 중`; }
  });
  grid.addEventListener('pointerup', event => {
    if (!drag) return;
    const current = drag; drag = null;
    if (grid.hasPointerCapture(event.pointerId)) grid.releasePointerCapture(event.pointerId);
    suppressClick = true; // The captured click targets the grid, so handle both click and drag here.
    focusPosition = current.end;
    if (current.moved) commitSelection(current.start, current.end); else endpoint(current.start);
    grid.querySelector(`[data-base="${focusPosition}"]`)?.focus({ preventScroll: true });
  });
  grid.addEventListener('pointercancel', () => { drag = null; suppressClick = false; refresh(); });
  grid.addEventListener('click', event => {
    if (suppressClick) { suppressClick = false; return; }
    const cell = event.target.closest('[data-base]'); if (cell) endpoint(Number(cell.dataset.base));
  });
  grid.addEventListener('keydown', event => {
    const cell = event.target.closest('[data-base]'); if (!cell) return;
    const current = Number(cell.dataset.base), start = view().windowStart, end = Math.min(getFixture().templates.A.sequence.length, start + size() - 1);
    if (event.key === 'Escape') { event.preventDefault(); anchor = null; refresh(); return; }
    if (event.key === 'Enter' || event.key === ' ') { event.preventDefault(); endpoint(current); return; }
    const delta = { ArrowLeft: -1, ArrowRight: 1, ArrowUp: -columns, ArrowDown: columns }[event.key];
    if (delta === undefined && !['Home', 'End'].includes(event.key)) return;
    event.preventDefault();
    focusPosition = event.key === 'Home' ? start : event.key === 'End' ? end : Math.max(1, Math.min(getFixture().templates.A.sequence.length, current + delta));
    if (event.shiftKey && anchor === null) anchor = current;
    if (focusPosition < start || focusPosition > end) setWindow(focusPosition < start ? focusPosition : focusPosition - size() + 1);
    for (const button of grid.querySelectorAll('[data-base]')) button.tabIndex = Number(button.dataset.base) === focusPosition ? 0 : -1;
    grid.querySelector(`[data-base="${focusPosition}"]`)?.focus();
    if (event.shiftKey) {
      const original = anchor; commitSelection(original, focusPosition); anchor = original;
    }
  });

  function renderSequence() {
    const fixture = getFixture(); if (!fixture || !grid.clientWidth) return;
    columns = Math.max(1, Math.min(14, Math.floor(grid.clientWidth / 44)));
    const start = Math.min(view().windowStart, fixture.templates.A.sequence.length - size() + 1), end = start + size() - 1;
    view().windowStart = start;
    const focused = document.activeElement?.dataset.base;
    if (focusPosition < start || focusPosition > end) focusPosition = start;
    grid.replaceChildren();
    for (let offset = start; offset <= end; offset += columns) {
      const row = node('div', '', 'pcr-sequence-row'), bases = node('div', '', 'pcr-bases');
      bases.style.setProperty('--pcr-columns', columns);
      for (let position = offset; position <= Math.min(end, offset + columns - 1); position++) {
        const base = fixture.templates.A.sequence[position - 1], cell = node('button', '', 'pcr-base');
        cell.type = 'button'; cell.dataset.base = position;
        cell.tabIndex = position === focusPosition ? 0 : -1;
        cell.setAttribute('aria-label', `${view().mode} 선택, A ${position}번 ${base}${position >= fixture.templates.B.deletion.start && position <= fixture.templates.B.deletion.end ? ', B 결실 영역' : ''}`);
        cell.append(node('span', String(position), 'pcr-base-coordinate'), node('span', base));
        cell.dataset.deletion = String(position >= fixture.templates.B.deletion.start && position <= fixture.templates.B.deletion.end);
        bases.append(cell);
      }
      row.append(bases); grid.append(row);
    }
    const binding = selected();
    paintSelection(anchor ?? Number(binding?.start), anchor ?? Number(binding?.end));
    $('sequence-window').textContent = `A ${start}–${end} / ${fixture.templates.A.sequence.length} nt`;
    $('sequence-window-slider').max = fixture.templates.A.sequence.length - size() + 1;
    $('sequence-window-slider').value = start;
    $('sequence-window-slider').setAttribute('aria-valuetext', `${start}–${end}번 확대`);
    $('sequence-prev').disabled = start === 1; $('sequence-next').disabled = end === fixture.templates.A.sequence.length;
    if (focused) grid.querySelector(`[data-base="${focused}"]`)?.focus({ preventScroll: true });
  }
  function renderMap() {
    if (!getFixture() || !model) return;
    const fixture = getFixture(), length = fixture.templates.A.sequence.length, deletion = fixture.templates.B.deletion;
    const x = position => 35 + position / length * 530;
    const reference = $('map-reference'), bindings = $('map-bindings'), prediction = $('map-prediction');
    reference.replaceChildren(); bindings.replaceChildren(); prediction.replaceChildren();
    label(reference, 35, 24, '1'); label(reference, 565, 24, `${length} bp`, 'end');
    path(reference, 'M35 66H565', 'pcr-strand');
    reference.append(svgNode('rect', { x: x(deletion.start - 1), y: 53, width: x(deletion.end) - x(deletion.start - 1), height: 26, class: 'pcr-map-deletion' }));
    label(reference, (x(deletion.start - 1) + x(deletion.end)) / 2, 100, `${deletion.start}–${deletion.end} 결실`, 'middle');
    for (const [name, info] of Object.entries(model.primers)) {
      const binding = info.binding;
      if (!binding || !Number.isInteger(Number(binding.start)) || !Number.isInteger(Number(binding.end)) || Number(binding.start) < 1 || Number(binding.end) < Number(binding.start) || Number(binding.end) > length) continue;
      const start = x(Number(binding.start) - 1), end = x(Number(binding.end)), y = name === 'F' ? 134 : 172;
      direction(bindings, binding.direction === 'right' ? start : end, binding.direction === 'right' ? end : start, y, name);
      const center = Math.max(90, Math.min(510, (start + end) / 2));
      label(bindings, center, y - 10, binding.direction === 'right' ? `5′ ${name} → 3′` : `3′ ← ${name} 5′`, 'middle', `pcr-map-primer-label pcr-primer-${name}`);
    }
    if (model.selectedProduct) {
      const product = model.selectedProduct;
      path(bindings, `M${x(product.start)} 195V203H${x(product.end)}V195`, 'pcr-map-amplicon');
      label(bindings, Math.max(130, Math.min(470, (x(product.start) + x(product.end)) / 2)), 223, `${product.length} bp amplicon`, 'middle');
    } else label(bindings, 300, 217, '선택 위치의 amplicon 미정', 'middle');
    if (view().showPrediction) for (const [name, key] of Object.entries(PRIMER_KEYS)) {
      const value = getState().initialPrimerPrediction[key]; if (value === null) continue;
      const px = 35 + value / 100 * 530;
      path(prediction, `M${px} 29V78`, 'pcr-initial-guide');
      label(prediction, Math.max(70, Math.min(530, px)), name === 'F' ? 42 : 57, `00 예상 ${name}`, 'middle', 'pcr-initial-label');
    }
    $('show-initial-prediction').setAttribute('aria-pressed', String(view().showPrediction));
    const predictions = Object.entries(PRIMER_KEYS).filter(([, key]) => getState().initialPrimerPrediction[key] !== null);
    $('map-prediction-legend').hidden = !view().showPrediction || !predictions.length;
    $('map-prediction-legend').textContent = '점선 / 00의 대략적 예측: ' + predictions.map(([name, key]) => `${name} ${getState().initialPrimerPrediction[key]}% 위치`).join(' / ');
    const description = `A ${length} bp. 결실 ${deletion.start}–${deletion.end}. F ${boundsText(model.primers.F.binding)}, R ${boundsText(model.primers.R.binding)}. ${model.selectedProduct ? `예상 amplicon ${model.selectedProduct.start + 1}–${model.selectedProduct.end}, ${model.selectedProduct.length} bp.` : model.placement || model.errors.join(' ')}`;
    $('map-description').textContent = description;
    $('map-caption').textContent = view().showPrediction ? (Object.values(PRIMER_KEYS).some(key => getState().initialPrimerPrediction[key] !== null) ? '점선: 00의 대략적 예측 / 파란 F, 빨간 R 화살표: 현재 선택. 좌표는 양 끝을 포함합니다.' : '00에서 저장한 초기 위치 예측이 없습니다. 현재 선택은 그대로 유지됩니다.') : '윤곽 구간: B에서 결실된 영역 / 화살표: 합성 방향. 좌표는 양 끝을 포함합니다.';
  }
  function renderPrimer(name) {
    const info = model.primers[name], container = $(name === 'F' ? 'primer-summary-f' : 'primer-summary-r');
    container.replaceChildren(node('h4', name === 'F' ? 'Forward primer' : 'Reverse primer'));
    const meta = info.binding ? `A ${boundsText(info.binding)}${info.inferred ? ' / 서열에서 확인' : ''}` : info.hits.length ? `A 결합 위치 ${info.hits.length}개` : info.stats ? 'A의 완전 일치 결합 부위 없음' : '결합 위치 선택 전';
    container.append(node('p', meta, 'pcr-primer-coordinate'));
    if (info.stats) {
      container.append(node('p', `${info.stats.length} nt / GC ${info.stats.gcPercent.toFixed(1)}%`, 'pcr-primer-metrics'));
      if (name === 'R' && info.reference) {
        container.append(node('p', 'Reference 결합 부위', 'pcr-small'), node('p', `5′ ${info.reference} 3′`, 'pcr-sequence pcr-reference-sequence'));
        container.append(node('p', info.binding.direction === 'left' ? '↓ reverse complement' : '↓ 같은 서열 / 오른쪽 합성', 'pcr-small'));
      }
      container.append(node('p', '주문 서열', 'pcr-small'), node('p', `5′ ${info.stats.sequence} 3′`, 'pcr-sequence pcr-ordered-sequence'));
      if (info.stats.length < 18 || info.stats.length > 25) container.append(node('p', '일반적인 출발 길이 18–25 nt와 다릅니다. 선택은 유지되며 활동 04에서 검토합니다.', 'pcr-small'));
    } else if (!info.error) container.append(node('p', '확대 서열에서 시작과 끝을 선택하세요.', 'pcr-small'));
    if (info.error) container.append(node('p', info.error, 'pcr-small'));
  }
  function renderProducts() {
    const results = $('analysis-results'), comparison = $('product-comparison');
    results.replaceChildren(); comparison.replaceChildren();
    if (!model.result) { results.append(node('p', model.errors.length ? '현재 입력으로 계산할 수 없습니다.' : 'F와 R 선택을 기다리고 있습니다.')); return; }
    const list = node('dl', '', 'pcr-product-values'); results.append(list);
    for (const [source, data] of Object.entries(model.result)) {
      const sizes = [...new Set(data.products.map(p => p.length))];
      list.append(node('dt', source), node('dd', sizes.length ? `${sizes.slice(0, 8).join(', ')} bp${sizes.length > 8 ? ` 외 ${sizes.length - 8}개 길이` : ''}` : '예상 산물 없음'));
      const block = node('div', '', 'pcr-product-source'); block.append(node('h4', source));
      if (!data.products.length) block.append(node('p', '현재 완전 일치 모형에서 대응되는 primer pair 없음', 'pcr-small'));
      for (const product of data.products.slice(0, 3)) {
        const svg = svgNode('svg', { viewBox: '0 0 300 82', role: 'img', 'aria-label': `${source} ${product.start + 1}–${product.end}, ${product.length} bp` });
        const width = 244 * product.length / getFixture().templates.A.sequence.length, start = 28, end = start + width;
        const binding = product.bindings[0];
        path(svg, `M${start} 26H${end}`, 'pcr-strand');
        direction(svg, start, start + (binding.right.end - binding.right.start) / product.length * width, 26);
        direction(svg, end, end - (binding.left.end - binding.left.start) / product.length * width, 26);
        label(svg, start, 15, `${binding.right.primer} →`); label(svg, Math.max(start + 65, end), 15, `← ${binding.left.primer}`, 'end');
        path(svg, `M${start} 40V48H${end}V40`, 'pcr-map-amplicon');
        label(svg, start, 74, `${product.length} bp / ${product.start + 1}–${product.end}`);
        block.append(svg);
      }
      if (data.products.length > 3) block.append(node('p', `산물 ${data.products.length}개 중 3개 구조 표시`, 'pcr-small'));
      comparison.append(block);
    }
    const details = node('details'), summary = node('summary', '결합과 산물 계산 내역'); details.append(summary);
    for (const [source, data] of Object.entries(model.result)) {
      details.append(node('p', `${source}: 완전 일치 결합 ${data.hits.length}개 / 산물 ${data.products.length}개`, 'pcr-small'));
      for (const p of data.products.slice(0, 20)) details.append(node('p', `${source} ${p.start + 1}~${p.end}: ${p.length} bp / ${[...new Set(p.bindings.map(b => `${b.right.primer}→ ←${b.left.primer}`))].join(', ')}`, 'pcr-small'));
      if (data.products.length > 20) details.append(node('p', `전체 ${data.products.length}개를 계산했으며 좌표 목록은 처음 20개만 표시합니다. 더 긴 primer로 결합 부위를 좁혀 보세요.`, 'pcr-small'));
    }
    results.append(details);
  }
  function refresh(controls = true) {
    if (!getFixture()) return null;
    model = inspectDesign(getState().draft, getFixture());
    $('primer-f').value = getState().draft.forward; $('primer-r').value = getState().draft.reverse;
    if (controls) syncControls();
    for (const button of root.querySelectorAll('[data-select-primer]')) {
      const name = button.dataset.selectPrimer;
      button.setAttribute('aria-pressed', String(name === view().mode));
      button.textContent = `${name} 선택 ${(model.primers[name].binding?.direction ?? (name === 'F' ? 'right' : 'left')) === 'right' ? '→' : '←'}`;
    }
    $('selection-title').textContent = view().mode === 'F' ? 'Forward primer 선택' : 'Reverse primer 결합 부위 선택';
    $('selection-title').dataset.mode = view().mode; $('sequence-grid').dataset.mode = view().mode;
    const binding = selected();
    $('range-status').textContent = binding ? `${view().mode} ${boundsText(binding)} 선택${model.primers[view().mode].stats ? ` / ${model.primers[view().mode].stats.length} nt` : ''}. 시작과 끝을 다시 선택해 옮길 수 있습니다.` : '시작 염기와 끝 염기를 차례로 선택하세요.';
    $('selection-direction').textContent = binding ? (binding.direction === 'left' ? '3′ ←──────── 5′ / 왼쪽으로 합성' : '5′ ────────→ 3′ / 오른쪽으로 합성') : '';
    $('placement-message').textContent = model.placement;
    $('deletion-message').textContent = model.deletionMessage; $('binding-message').textContent = model.bindingMessage;
    $('placement-explanation').hidden = !model.placement && !model.deletionMessage && !model.bindingMessage;
    $('analysis-status').textContent = model.errors.join(' ') || (model.result ? '현재 주문 서열의 완전 일치 결과입니다.' : '두 primer를 선택하면 즉시 계산합니다.');
    renderMap(); renderSequence(); renderPrimer('F'); renderPrimer('R'); renderProducts();
    onComputed(model); return model;
  }
  let lastWidth = 0;
  new ResizeObserver(entries => {
    const width = Math.floor(entries[0].contentRect.width);
    if (width !== lastWidth) { lastWidth = width; renderSequence(); }
  }).observe(grid);
  return { refresh, resetSelection() { anchor = null; drag = null; focusPosition = 1; refresh(); }, focus() { $('sequence-grid [tabindex="0"]')?.focus(); }, cancelSelection() { anchor = null; } };
}
