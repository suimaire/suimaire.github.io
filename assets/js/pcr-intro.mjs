import { complement, reverseComplement, normalizeSequence } from './pcr-core.mjs';
import { PCR_STAGES } from './pcr-records.mjs';

// Display-only geometry on the existing relative-percent scale. The short binding
// bands span one slider step; they are NOT a selected nucleotide/primer length.
const DELETION_LEFT = 120 / 420 * 100;
const DELETION_RIGHT = 200 / 420 * 100;
const B_SCALE = 340 / 420;
const ROUGH_HALF_WIDTH = 2.5;
export function predictSharedPair(prediction) {
  const binding = value => {
    if (value === null) return { value, overlap: false, b: null };
    const start = value - ROUGH_HALF_WIDTH, end = value + ROUGH_HALF_WIDTH;
    const overlap = end > DELETION_LEFT && start < DELETION_RIGHT;
    return { value, start, end, overlap, b: overlap ? null : (value - (start >= DELETION_RIGHT ? DELETION_RIGHT - DELETION_LEFT : 0)) / B_SCALE };
  };
  const forward = binding(prediction.forward), reverse = binding(prediction.reverse);
  const complete = forward.value !== null && reverse.value !== null;
  const ordered = complete && forward.end < reverse.start;
  const enclosesDeletion = ordered && forward.end <= DELETION_LEFT && reverse.start >= DELETION_RIGHT;
  const shared = complete && !forward.overlap && !reverse.overlap;
  const state = forward.overlap || reverse.overlap ? 'overlap' : !complete ? 'incomplete' : !ordered ? 'order' : !enclosesDeletion ? 'outside' : 'valid';
  return { forward, reverse, ordered, shared, enclosesDeletion, state };
}

// These are conceptual drawings only. They do not run the sequence calculator.
export function initializeIntro(root, getState, save) {
  const $ = id => root.querySelector(`#${id}`);
  const all = selector => [...root.querySelectorAll(selector)];
  const view = () => getState().introView;
  const ns = 'http://www.w3.org/2000/svg';
  const svgNode = (tag, attributes = {}, text) => {
    const node = document.createElementNS(ns, tag);
    for (const [key, value] of Object.entries(attributes)) node.setAttribute(key, value);
    if (text !== undefined) node.textContent = text;
    return node;
  };
  const line = (svg, x1, y1, x2, y2, className = 'pcr-strand') => svg.appendChild(svgNode('path', { d: `M${x1} ${y1}L${x2} ${y2}`, class: className }));
  const label = (svg, x, y, text, anchor = 'start', className = '') => svg.appendChild(svgNode('text', { x, y, 'text-anchor': anchor, class: className }, text));
  const arrow = (svg, x1, y, x2, className) => {
    const sign = x2 > x1 ? 1 : -1;
    line(svg, x1, y, x2, y, className);
    svg.append(svgNode('path', { d: `M${x2 - sign * 9} ${y - 6}L${x2} ${y}L${x2 - sign * 9} ${y + 6}`, class: className }));
  };
  function strand(svg, y, reference) {
    line(svg, 45, y, 555, y);
    label(svg, 7, y + 6, reference ? '5′' : '3′'); label(svg, 568, y + 6, reference ? '3′' : '5′');
  }
  function primer(svg, five, three, y, name, labelAbove, extending = false) {
    arrow(svg, five, y, three, 'pcr-intro-primer');
    const ly = y + (labelAbove ? -16 : 29);
    label(svg, five, ly, '5′', 'middle');
    if (!extending) label(svg, three, ly, '3′', 'middle', 'pcr-three-prime');
    label(svg, (five + three) / 2, ly, name, 'middle');
  }
  const choose = (selector, selected, key) => all(selector).forEach(button => button.setAttribute('aria-pressed', String(button.dataset[key] === String(selected))));
  const roughRegion = value => value === null ? '아직 표시하지 않음' : value < DELETION_LEFT ? '결실 구간의 왼쪽' : value <= DELETION_RIGHT ? '결실 구간 안쪽' : '결실 구간의 오른쪽';
  let activePrimer = 'forward';
  function renderPrediction() {
    const model = predictSharedPair(getState().initialPrimerPrediction);
    const bDescriptions = [];
    for (const name of ['forward', 'reverse']) {
      const { value, start, overlap, b } = model[name];
      const shortName = name === 'forward' ? 'F' : 'R';
      const marker = $(`prediction-${name}-marker`), slider = $(`prediction-${name}`);
      const bMarker = $(`prediction-b-${name}-marker`);
      marker.hidden = value === null;
      marker.style.left = `${value ?? 50}%`;
      marker.dataset.overlap = String(overlap);
      marker.textContent = name === 'forward' ? `F${overlap ? ' !' : ''} →` : `← R${overlap ? ' !' : ''}`;
      const description = `${roughRegion(value)}${overlap ? ' / 결실과 겹침 주의' : value !== null ? ' / B에도 대응 위치 표시' : ''}`;
      marker.setAttribute('aria-label', `A의 ${shortName}: ${description}`);
      bMarker.hidden = b === null;
      bMarker.style.left = `${b ?? 50}%`;
      bMarker.setAttribute('aria-label', `B의 ${shortName}: A에서 고른 같은 primer가 인식하는 보존된 결합 부위`);
      for (const [id, left, width, hidden] of [
        [`prediction-${name}-binding`, start, ROUGH_HALF_WIDTH * 2, value === null],
        [`prediction-b-${name}-binding`, b - ROUGH_HALF_WIDTH / B_SCALE, ROUGH_HALF_WIDTH * 2 / B_SCALE, b === null]
      ]) {
        const band = $(id);
        band.hidden = hidden; band.style.left = `${left ?? 0}%`; band.style.width = `${width}%`;
        band.dataset.overlap = String(overlap);
      }
      bDescriptions.push(`${shortName}: ${value === null ? '아직 위치를 정하지 않았습니다' : overlap ? '결실과 겹쳐 같은 연속 결합 부위를 표시할 수 없습니다' : 'B에 남아 있는 결합 부위에 같은 primer가 결합합니다'}`);
      slider.value = value ?? 50;
      slider.setAttribute('aria-valuetext', description);
      $(`prediction-${name}-value`).textContent = description;
    }
    for (const [sample, f, r, visible] of [
      ['a', model.forward.value, model.reverse.value, model.ordered],
      ['b', model.forward.b, model.reverse.b, model.ordered && model.shared]
    ]) {
      const span = $(`prediction-${sample}-span`);
      span.hidden = !visible; span.style.left = `${f ?? 0}%`; span.style.width = `${visible ? r - f : 0}%`;
    }
    const messages = {
      incomplete: 'F와 R의 A 기준 위치를 모두 정해 보세요. 보존된 위치를 고르면 같은 primer의 B 대응 위치도 함께 표시됩니다.',
      overlap: '겹침 주의: 한 primer의 개략적 결합 범위가 121~200 결실 구간과 겹칩니다. B에서 같은 연속 부위에 결합하기 어려워 해당 표식을 표시하지 않았습니다. 두 primer 모두 A와 B에 결합할 수 있는 위치인지 확인하세요.',
      order: '배치 확인: F가 왼쪽, R이 오른쪽에서 서로를 향하고 두 결합 부위가 떨어져 있어야 합니다. 결실 구간 전체가 두 primer 사이에 포함되는지도 확인하세요.',
      outside: '결실 포함 확인: 두 primer는 A와 B의 대응 위치에 있지만, 121~200 결실 구간 전체가 두 결합 부위 사이에 들어가지 않습니다. 이 배치의 primer 사이 구간은 A와 B에서 같으므로 결실로 인한 80 bp 산물 차이를 볼 수 없습니다.',
      valid: '개념 조건 충족: 같은 F/R이 A와 B의 보존된 대응 부위에 있고, 121~200 결실 구간 전체가 두 결합 부위 사이에 들어갑니다. 이 가정에서 B의 PCR 산물은 A보다 80 bp 짧습니다. 가능한 배치는 여러 곳이며, 실제 서열과 품질은 뒤 활동에서 확인합니다.'
    };
    $('prediction-feedback').textContent = messages[model.state];
    $('prediction-feedback').dataset.state = ['incomplete', 'valid'].includes(model.state) ? model.state : 'warning';
    $('prediction-b-caption').textContent = model.shared
      ? 'A와 B에 동일한 F/R primer를 사용합니다. B에서도 결실 바깥의 보존된 결합 부위에 같은 primer가 결합합니다.'
      : `${bDescriptions.join('. ')}.`;
    $('prediction-track').setAttribute('aria-label', `A 기준의 대략적 예측 위치. Forward: ${roughRegion(model.forward.value)}. Reverse: ${roughRegion(model.reverse.value)}. 아래 조절 막대로 같은 primer pair의 위치를 바꿀 수 있습니다.`);
    $('prediction-b-track').setAttribute('aria-label', `B에 자동 대응한 같은 primer pair. ${bDescriptions.join('. ')}.`);
    $('legacy-first-negative').hidden = !getState().answers['first-negative'];
  }
  function place(name, value) {
    getState().initialPrimerPrediction[name] = Math.max(5, Math.min(95, Math.round(value / 5) * 5));
    renderPrediction(); save();
  }
  for (const name of ['forward', 'reverse']) $(`prediction-${name}`).addEventListener('input', event => place(name, Number(event.target.value)));
  all('[data-prediction-primer]').forEach(button => button.addEventListener('click', () => {
    activePrimer = button.dataset.predictionPrimer;
    choose('[data-prediction-primer]', activePrimer, 'predictionPrimer');
  }));
  $('prediction-track').addEventListener('click', event => {
    const bounds = event.currentTarget.getBoundingClientRect();
    place(activePrimer, (event.clientX - bounds.left) / bounds.width * 100);
  });

  const stages = {
    mixture: ['반응 혼합', '상보적인 두 DNA 가닥과 프라이머, 효소, 합성 재료를 함께 준비합니다. 아래의 F와 R은 아직 주형에 결합하지 않은 상태입니다.', []],
    denaturation: ['변성 / 예시 95°C', '상보적인 두 가닥의 간격이 벌어집니다. DNA의 phosphodiester backbone, 즉 당-인산 골격을 끊는 과정이 아닙니다.', ['template']],
    annealing: ['결합 / 예시 60°C', 'Forward와 Reverse primer가 각각 상보적인 위치에 결합합니다. 두 위치가 증폭될 구간의 경계를 정합니다. 실제 annealing temperature는 프라이머와 반응 조건에 따라 달라집니다.', ['forward', 'reverse']],
    extension: ['신장 / 예시 72°C', 'Pol은 primer의 3′ OH에서 연장이 시작되는 위치를 표시합니다. 새 가닥은 5′→3′로 자랍니다. 첫 주기에는 반대편 primer 위치를 지나 더 긴 산물이 생길 수 있습니다.', ['polymerase', 'dntp']]
  };
  function renderCycle() {
    const stage = view().stage, [title, description, mixture] = stages[stage], svg = $('cycle-svg');
    $('cycle-title').textContent = title; $('cycle-description').textContent = description;
    svg.replaceChildren(); svg.dataset.stage = stage; svg.setAttribute('aria-label', `${title}. ${description}`);
    const mixed = stage === 'mixture';
    strand(svg, mixed ? 94 : 42, true); strand(svg, mixed ? 126 : 231, false);
    if (mixed) {
      for (let x = 65; x < 550; x += 25) line(svg, x, 98, x, 122, 'pcr-pair-guide');
      primer(svg, 100, 180, 205, 'F', true); primer(svg, 500, 420, 205, 'R', true);
      label(svg, 300, 221, 'Pol', 'middle');
    }
    if (stage === 'annealing' || stage === 'extension') {
      const extending = stage === 'extension';
      primer(svg, 120, 210, 197, 'F', true, extending);
      primer(svg, 480, 390, 77, 'R', false, extending);
      svg.append(svgNode('path', { d: 'M120 130V142H480V130', class: 'pcr-boundary-guide' }));
      if (!extending) label(svg, 300, 169, '두 primer가 정하는 경계', 'middle');
      if (extending) {
        arrow(svg, 210, 197, 550, 'pcr-new-dna pcr-line-growth');
        arrow(svg, 390, 77, 50, 'pcr-new-dna pcr-line-growth');
        label(svg, 553, 180, '3′', 'middle'); label(svg, 48, 110, '3′', 'middle');
        for (const [x, y] of [[212, 167], [352, 89]]) {
          svg.append(svgNode('rect', { x, y, width: 40, height: 25, rx: 2, class: 'pcr-pol' }));
          label(svg, x + 20, y + 19, 'Pol', 'middle', 'pcr-pol-label');
        }
        label(svg, 300, 165, '새 가닥 5′→3′', 'middle');
      }
    }
    choose('[data-cycle-stage]', stage, 'cycleStage');
    all('[data-mixture]').forEach(item => { item.dataset.active = String(mixture.includes(item.dataset.mixture)); });
    $('cycle-prev').disabled = stage === 'mixture'; $('cycle-next').disabled = stage === 'extension';
  }
  function setStage(stage) {
    view().stage = stage;
    getState().cycle = Math.max(0, PCR_STAGES.indexOf(stage) - 1); // Preserve the legacy field's meaning.
    renderCycle(); save();
  }
  all('[data-cycle-stage]').forEach(button => button.addEventListener('click', () => setStage(button.dataset.cycleStage)));
  $('cycle-prev').addEventListener('click', () => setStage(PCR_STAGES[Math.max(0, PCR_STAGES.indexOf(view().stage) - 1)]));
  $('cycle-next').addEventListener('click', () => setStage(PCR_STAGES[Math.min(3, PCR_STAGES.indexOf(view().stage) + 1)]));
  $('cycle-reset').addEventListener('click', () => setStage('mixture'));

  const products = {
    1: '1주기: 한쪽 끝만 primer로 정해진 긴 새 가닥이 만들어집니다. 반대쪽 끝은 아직 다른 primer 위치로 정해지지 않았습니다.',
    2: '2주기: 한쪽 끝이 이미 primer로 정해진 긴 새 가닥이 다음 주기의 주형이 됩니다. 반대쪽 primer에서 합성하면 이 주형의 끝까지 복사하므로 정확한 길이의 단일가닥이 처음 생깁니다. 긴 산물도 함께 존재합니다.',
    3: '3주기: 정확한 길이의 이중가닥 산물이 처음 생깁니다. 이후 주기를 거치며 축적되고, 긴 산물도 함께 존재합니다.'
  };
  function renderProducts() {
    const cycle = view().comparisonCycle, svg = $('cycle-products-svg');
    svg.replaceChildren(); svg.setAttribute('aria-label', products[cycle]);
    if (cycle === 1) {
      label(svg, 30, 24, '원래 주형 / 긴 가닥');
      line(svg, 45, 58, 555, 58); label(svg, 12, 64, '3′'); label(svg, 565, 64, '5′');
      label(svg, 30, 100, '새 가닥 / 한쪽 끝만 정해짐');
      line(svg, 140, 136, 210, 136, 'pcr-intro-primer'); arrow(svg, 210, 136, 550, 'pcr-new-dna');
      label(svg, 128, 142, '5′', 'end'); label(svg, 565, 142, '3′');
      line(svg, 140, 147, 140, 185, 'pcr-pair-guide'); label(svg, 140, 210, 'primer에서 시작', 'middle');
      line(svg, 440, 119, 440, 155, 'pcr-boundary-guide');
      label(svg, 425, 184, '반대 primer 위치를 넘어 합성', 'middle');
      label(svg, 300, 254, '반대쪽 끝은 아직 다른 primer로 정해지지 않음', 'middle');
    } else if (cycle === 2) {
      label(svg, 30, 24, '1주기의 새 가닥 → 이번 주기의 주형');
      line(svg, 140, 62, 550, 62); label(svg, 125, 68, '5′', 'end'); label(svg, 565, 68, '3′');
      line(svg, 440, 140, 370, 140, 'pcr-intro-primer'); arrow(svg, 370, 140, 140, 'pcr-new-dna');
      label(svg, 125, 146, '3′', 'end'); label(svg, 455, 146, '5′');
      for (const x of [140, 440]) line(svg, x, 46, x, 170, 'pcr-boundary-guide');
      label(svg, 440, 109, '반대 primer에서 시작', 'middle');
      label(svg, 140, 199, '주형의 끝에서 종료', 'middle');
      label(svg, 300, 254, '양 끝이 정해진 단일가닥이 생김', 'middle');
    } else {
      svg.append(svgNode('path', { d: 'M140 50V36H440V50', class: 'pcr-pair-guide' }));
      label(svg, 290, 24, '정확한 목표 길이', 'middle');
      for (const x of [140, 440]) line(svg, x, 56, x, 166, 'pcr-boundary-guide');
      line(svg, 140, 91, 440, 91, 'pcr-new-dna'); line(svg, 140, 129, 440, 129, 'pcr-new-dna');
      label(svg, 125, 97, '5′', 'end'); label(svg, 455, 97, '3′'); label(svg, 125, 135, '3′', 'end'); label(svg, 455, 135, '5′');
      for (let x = 155; x < 440; x += 25) line(svg, x, 97, x, 123, 'pcr-pair-guide');
      line(svg, 140, 91, 210, 91, 'pcr-intro-primer'); line(svg, 370, 129, 440, 129, 'pcr-intro-primer');
      label(svg, 290, 196, '정확한 길이 / 이중가닥', 'middle');
      label(svg, 290, 254, '이후 축적됨 / 긴 산물도 함께 존재', 'middle');
    }
    $('cycle-products-description').textContent = products[cycle];
    choose('[data-comparison-cycle]', cycle, 'comparisonCycle');
  }
  all('[data-comparison-cycle]').forEach(button => button.addEventListener('click', () => {
    view().comparisonCycle = Number(button.dataset.comparisonCycle); renderProducts(); save();
  }));

  function renderDirection() {
    const { flipped, arrangement } = view(), svg = $('direction-svg'), inward = arrangement === 'inward';
    svg.replaceChildren();
    const referenceY = flipped ? 229 : 48, complementY = flipped ? 48 : 229;
    strand(svg, referenceY, true); strand(svg, complementY, false);
    label(svg, 45, flipped ? 258 : 22, 'Reference strand');
    const fy = flipped ? 88 : 188, ry = inward ? (flipped ? 188 : 88) : fy;
    primer(svg, 110, 190, fy, 'F', !flipped);
    primer(svg, inward ? 490 : 395, inward ? 410 : 475, ry, 'R', inward ? flipped : !flipped);
    arrow(svg, 195, fy, 290, 'pcr-new-dna');
    arrow(svg, inward ? 405 : 480, ry, inward ? 310 : 552, 'pcr-new-dna');
    svg.append(svgNode('path', { d: 'M110 130V139H490V130', class: 'pcr-boundary-guide' }));
    label(svg, 300, 165, '표적 내부', 'middle');
    const description = inward ? '두 primer의 3′ 말단에서 서로를 향해 합성이 진행됩니다. 검정 화살표는 각 3′ 말단에서의 합성 방향입니다.' : '잘못된 예에서는 두 primer가 같은 가닥에 결합해 오른쪽으로 합성합니다. R의 3′ 말단도 바깥을 향하므로 표적 구간의 양쪽 경계를 정의할 수 없습니다.';
    $('arrangement-description').textContent = description;
    svg.setAttribute('aria-label', `${inward ? '정상 배치' : '잘못된 배치'}. ${description}`);
    $('direction-caption').textContent = `이 그림에서는 ${flipped ? '아래쪽' : '위쪽'} 가닥을 reference strand로 5′→3′ 방향으로 그렸습니다. ${flipped ? '두 가닥의 위아래 위치만 바꾸었으며 각 가닥의 5′, 3′ 방향은 유지했습니다.' : '위아래 위치보다 5′, 3′ 방향과 합성 방향이 중요합니다.'}`;
    choose('[data-arrangement]', arrangement, 'arrangement'); $('flip-strands').setAttribute('aria-pressed', String(flipped));
  }
  all('[data-arrangement]').forEach(button => button.addEventListener('click', () => { view().arrangement = button.dataset.arrangement; renderDirection(); save(); }));
  $('flip-strands').addEventListener('click', () => { view().flipped = !view().flipped; renderDirection(); save(); });

  function correctComplement() {
    try { return normalizeSequence($('direction-complement').value) === complement('AGTCCGTA'); } catch { return false; }
  }
  function renderComplement() {
    const confirmed = view().complementConfirmed && correctComplement();
    $('reverse-step').hidden = !confirmed;
    $('reverse-demonstration').hidden = !confirmed || !view().reverseShown;
    $('complement-feedback').textContent = confirmed ? '상보 서열이 일치합니다. 이제 같은 가닥을 주문 방향으로 읽어 보세요.' : '';
    $('direction-feedback').textContent = '';
  }
  $('check-complement').addEventListener('click', () => {
    view().complementConfirmed = correctComplement(); view().reverseShown = false; renderComplement();
    if (!view().complementConfirmed) $('complement-feedback').textContent = '다시 확인하세요. A–T, G–C를 짝짓고 왼쪽 3′와 오른쪽 5′ 방향을 유지하세요.';
    save();
  });
  $('direction-complement').addEventListener('input', () => { view().complementConfirmed = false; view().reverseShown = false; renderComplement(); save(); });
  $('direction-reverse').addEventListener('input', () => { $('direction-feedback').textContent = ''; });
  $('reverse-complement').addEventListener('click', () => { view().reverseShown = true; renderComplement(); save(); });
  $('check-direction').addEventListener('click', () => {
    try {
      $('direction-feedback').textContent = normalizeSequence($('direction-reverse').value) === reverseComplement('AGTCCGTA') ? '주문용 역상보 서열이 일치합니다. F와 R 모두 5′→3′로 기록합니다.' : '다시 확인하세요. 아래 상보 가닥을 오른쪽 5′에서 왼쪽 3′로 읽어 보세요.';
    } catch (error) { $('direction-feedback').textContent = error.message; }
  });
  function renderChoices() {
    for (const fieldset of all('[data-choice]')) {
      const selected = getState().answers[fieldset.dataset.choice];
      fieldset.querySelectorAll('input').forEach(input => { input.checked = input.value === selected; });
    }
    const a = getState().answers;
    $('cycle-choice-feedback').textContent = !a['cycle-boundary-choice'] ? '' : a['cycle-boundary-choice'] === 'primer-pair' ? '맞습니다. primer pair의 결합 위치가 양 끝을 정합니다.' : '각 요소의 역할을 다시 살펴보세요. 주형에 결합해 합성 시작점을 정하는 요소는 무엇일까요?';
    $('direction-choice-feedback').textContent = !a['direction-end-choice'] ? '' : a['direction-end-choice'] === '3' ? '맞습니다. 3′ 말단이 내부를 향하며 이 말단에서 새 DNA 합성이 진행됩니다.' : 'DNA polymerase가 어느 말단의 OH에서 연장하는지 다시 살펴보세요.';
  }
  all('[data-choice] input').forEach(input => input.addEventListener('change', () => { getState().answers[input.name] = input.value; renderChoices(); save(); }));
  return {
    render() { renderPrediction(); renderCycle(); renderProducts(); renderDirection(); renderComplement(); renderChoices(); },
    preparePrint(blank, textNode) {
      for (const fieldset of all('[data-choice]')) {
        const selected = fieldset.querySelector('input:checked');
        fieldset.append(textNode('div', blank || !selected ? '' : selected.parentElement.textContent.trim(), 'pcr-print-value'));
      }
      for (const name of ['forward', 'reverse']) {
        $(`prediction-${name}`).after(textNode('div', blank ? '' : roughRegion(getState().initialPrimerPrediction[name]), 'pcr-print-value'));
      }
    }
  };
}
