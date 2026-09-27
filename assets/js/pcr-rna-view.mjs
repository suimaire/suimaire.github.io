import { RNA_TRANSCRIPTS, strategyModel, transcriptBinding } from './pcr-rna.mjs';

const node = (tag, text = '', className = '') => {
  const el = document.createElement(tag); el.textContent = text;
  if (className) el.className = className;
  return el;
};
const svgNode = (tag, attributes, text = '') => {
  const el = document.createElementNS('http://www.w3.org/2000/svg', tag);
  for (const [key, value] of Object.entries(attributes)) el.setAttribute(key, value);
  el.textContent = text; return el;
};
// A fixed readable canvas scrolls within its figure on narrow screens.
function structure(label, genomic, { strategy, intron = 'long', position = 'junction', summary = '' } = {}) {
  const figure = node('figure', '', 'rna-figure');
  figure.append(node('h4', label));
  const scroll = node('div', '', 'rna-diagram-scroll');
  scroll.tabIndex = 0; scroll.setAttribute('role', 'region'); scroll.setAttribute('aria-label', `${label} exon 구조 / 좁은 화면에서 좌우로 스크롤`);
  const svg = svgNode('svg', { viewBox: '0 0 620 130', role: 'img', 'aria-label': `${label}: ${genomic ? 'Exon 1, intron, Exon 2, intron, Exon 3' : 'Exon 1, Exon 2, Exon 3이 intron 없이 연결됨'}. ${summary}` });
  const starts = genomic ? [20, 250, 480] : [20, 140, 260];
  starts.forEach((x, i) => {
    svg.append(svgNode('rect', { x, y: 36, width: 120, height: 36, class: 'pcr-outline' }), svgNode('text', { x: x + 60, y: 60, 'text-anchor': 'middle' }, `Exon ${i + 1}`));
    if (genomic && i < 2) svg.append(svgNode('path', { d: `M${x + 120} 54H${starts[i + 1]}`, class: 'pcr-strand' }), svgNode('text', { x: (x + 120 + starts[i + 1]) / 2, y: 26, 'text-anchor': 'middle', class: 'rna-intron-label' }, strategy === 'intron-spanning' && i === 0 ? `${intron === 'long' ? '긴' : '짧은'} intron` : 'intron'));
  });
  function primer(x1, x2, name, right) {
    const tip = right ? x2 : x1, back = right ? tip - 7 : tip + 7;
    svg.append(svgNode('path', { d: `M${x1} 99H${x2} M${back} 94L${tip} 99L${back} 104`, class: 'rna-primer' }), svgNode('text', { x: (x1 + x2) / 2, y: 125, 'text-anchor': 'middle', class: 'rna-primer-label' }, name));
  }
  if (strategy) {
    if (strategy === 'junction' && position === 'junction') {
      if (!genomic) primer(105, 175, 'F', true);
      else {
        svg.append(svgNode('path', { d: 'M105 99H140 M250 99H285', class: 'rna-separated' }), svgNode('text', { x: 195, y: 126, 'text-anchor': 'middle' }, 'F 결합 서열 분리'));
      }
    } else primer(30, 70, 'F', true);
    const rx = strategy === 'same-exon' ? 90 : starts[1] + 70;
    primer(rx, rx + 40, 'R', false);
  }
  scroll.append(svg); figure.append(scroll);
  if (summary) figure.append(node('figcaption', summary));
  return figure;
}

export function initializeRna(root, getState, save) {
  const $ = id => root.querySelector(`#${id}`), activity = $('rna-extension');
  const rna = () => getState().rnaExtension;
  const panels = [...activity.querySelectorAll('[data-rna-panel]')];
  const fields = [...activity.querySelectorAll('[data-rna-answer]')];
  fields.forEach(f => { if (f.tagName === 'TEXTAREA') f.maxLength = 12000; });
  let printOpen;
  function sections() {
    for (const panel of panels) panel.hidden = panel.dataset.rnaPanel !== rna().activeSection;
    activity.querySelectorAll('[data-rna-section]').forEach(b => b.setAttribute('aria-pressed', String(b.dataset.rnaSection === rna().activeSection)));
  }
  function render() {
    const state = rna(); sections();
    activity.querySelectorAll('[data-rna-option]').forEach(b => b.setAttribute('aria-pressed', String(state[b.dataset.rnaOption] === b.dataset.value)));
    for (const f of fields) f.value = state.answers[f.dataset.rnaAnswer] ?? '';
    const flowText = {
      rna: '1 / RNA: 발현을 살펴볼 출발 물질입니다. 일반적인 RT-PCR에서는 먼저 cDNA로 전환합니다.',
      rt: '2 / reverse transcription: reverse transcriptase가 RNA를 바탕으로 cDNA를 만듭니다.',
      cdna: '3 / cDNA: RNA에서 만든 DNA이며, 다음 PCR 단계에 들어가는 template입니다.',
      pcr: '4 / PCR: primer가 cDNA에 결합하고 DNA amplicon이 증폭됩니다. qPCR도 RNA 자체를 직접 증폭하는 단계는 아닙니다.'
    };
    $('rna-flow-description').textContent = flowText[state.flowStep];
    const activeIndex = { rna: 0, rt: 1, cdna: 2, pcr: 3 }[state.flowStep];
    activity.querySelectorAll('.rna-flow-comparison .rna-flow-list')[1].querySelectorAll('li').forEach((li, i) => { li.dataset.active = String(i === activeIndex); });
    templateFeedback();
    const model = strategyModel(state.primerStrategy, state.intronExample, state.junctionPosition);
    $('rna-intron-options').hidden = state.primerStrategy !== 'intron-spanning';
    $('rna-junction-options').hidden = state.primerStrategy !== 'junction';
    $('rna-placement').textContent = model.placement;
    $('rna-strategy-meaning').textContent = model.meaning;
    $('rna-strategy-comparison').replaceChildren(...[false, true].map(genomic => structure(genomic ? 'Genomic DNA (gDNA)' : 'cDNA', genomic, {
      strategy: state.primerStrategy, intron: state.intronExample, position: state.junctionPosition,
      summary: `${model.placement} ${genomic ? model.genomic : model.cdna}`
    })));
    $('rna-minus-signal').textContent = state.rtControlCase === 'case1' ? 'band 없음' : 'band 있음';
    $('rna-control-meaning').textContent = state.rtControlCase === 'case1'
      ? 'Case 1: 이 조건에서 RT(-)의 band는 관찰되지 않았습니다. RNA-derived cDNA 신호라는 해석을 뒷받침하지만, DNA가 전혀 없거나 산물의 정체가 확인되었다는 뜻은 아닙니다.'
      : 'Case 2: RNA-derived cDNA가 아닌 DNA template의 존재, 예를 들어 gDNA carryover 등을 검토할 근거가 됩니다. gDNA contamination 확정은 아니며 기타 오염이나 실험 문제도 검토해야 합니다.';
    renderTranscripts();
  }
  function templateFeedback() {
    const answer = rna().answers.template;
    $('rna-template-feedback').textContent = !answer ? '' : answer === 'cdna' ? 'cDNA가 PCR template입니다. RNA → reverse transcription → cDNA의 연결을 확인했습니다.' : 'RNA는 출발 물질, 단백질은 다른 분자, dNTP는 합성 재료입니다. Reverse transcription 뒤에 만들어진 cDNA가 PCR template가 됩니다.';
  }
  function renderTranscripts() {
    const target = rna().transcriptTarget;
    const models = RNA_TRANSCRIPTS.map(exons => transcriptBinding(exons, target));
    $('rna-transcript-summary').textContent = `선택한 F 영역: ${models[0].label}. 아래에서 해당 exon 또는 연속 junction을 확인하세요.`;
    $('rna-transcript-comparison').replaceChildren(...RNA_TRANSCRIPTS.map((exons, index) => {
      const model = models[index], row = node('div', '', 'rna-transcript-row');
      row.append(node('h4', `Transcript ${index + 1}`));
      const track = node('div', '', 'rna-transcript-track'); track.setAttribute('role', 'img');
      track.setAttribute('aria-label', `Transcript ${index + 1}: Exon ${exons.join(' → ')}. ${model.reason}`);
      for (const exon of exons) {
        const block = node('span', `Exon ${exon}`, 'rna-exon');
        if (model.present && model.region.includes(exon)) block.dataset.target = 'true';
        track.append(block);
      }
      row.append(track, node('p', model.reason, 'pcr-small')); return row;
    }));
    $('rna-transcript-meaning').textContent = target === 'exon4'
      ? '공통 Exon 4를 이용하면 이 예의 세 transcript 모두에 결합 가능한 구조가 있습니다. 해당 primer pair가 결합 가능한 여러 transcript의 신호가 합쳐질 수 있으며, gene 전체의 절대 발현량을 뜻하지는 않습니다.'
      : `${models.filter(m => m.present).length}개 transcript에 선택한 결합 구조가 있습니다. 일부 transcript에만 존재하는 exon 또는 junction을 이용하면 특정 isoform 또는 transcript subset을 더 선택적으로 검출할 수 있습니다. 구조의 존재만으로 실제 검출이나 isoform 특이성이 보장되지는 않습니다.`;
  }
  activity.addEventListener('click', event => {
    const section = event.target.closest('[data-rna-section]');
    const option = event.target.closest('[data-rna-option]');
    if (section) { rna().activeSection = section.dataset.rnaSection; sections(); save(); }
    if (option) { rna()[option.dataset.rnaOption] = option.dataset.value; render(); save(); }
  });
  activity.addEventListener('input', event => {
    const field = event.target.closest('[data-rna-answer]'); if (!field) return;
    rna().answers[field.dataset.rnaAnswer] = field.value; templateFeedback(); save();
  });
  // Opening a hash link into the extension must expose its content, including from 06.
  const revealHash = () => { if (location.hash === '#rna-extension') activity.open = true; };
  root.addEventListener('click', event => { if (event.target.closest('a[href="#rna-extension"]')) activity.open = true; });
  window.addEventListener('hashchange', revealHash); revealHash();
  $('rna-structure-comparison').replaceChildren(
    structure('Genomic DNA (gDNA)', true, { summary: 'Exon 1과 2, Exon 2와 3 사이에 intron이 있습니다.' }),
    structure('Mature mRNA / cDNA', false, { summary: 'Mature mRNA의 Exon 1–2–3 연결 구조가 cDNA에 반영됩니다. 이 그림의 exon 사이에는 intron이 없습니다.' })
  );
  return {
    render,
    preparePrint(blank) {
      if (printOpen === undefined) printOpen = activity.open;
      activity.open = true; panels.forEach(p => { p.hidden = false; });
      // Existing print mirrors handle prose; show the selected option's label, not its key.
      const field = $('rna-template'), mirror = field.nextElementSibling;
      if (mirror?.classList.contains('pcr-print-value')) mirror.textContent = blank || !field.value ? '' : field.selectedOptions[0].textContent;
      if (blank) $('rna-template-feedback').hidden = true;
    },
    finishPrint() { if (printOpen !== undefined) activity.open = printOpen; printOpen = undefined; $('rna-template-feedback').hidden = false; sections(); }
  };
}
