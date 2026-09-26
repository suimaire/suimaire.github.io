import { analyzeAll } from './pcr-core.mjs';
import { REVIEW_LENSES, reviewDesign, endFeatures, complementarity, alignmentLines } from './pcr-review.mjs';

export function initializeReview(root, getState, getFixture, save) {
  const $ = id => root.querySelector(`#${id}`);
  const node = (tag, text = '', className = '') => {
    const el = document.createElement(tag); el.textContent = text; el.className = className; return el;
  };
  const view = () => getState().review;
  const labels = { F: 'Forward', R: 'Reverse' };
  const productsText = products => products.length ? products.map(p => `${p.source} ${p.start + 1}~${p.end}: ${p.length} bp`).join('\n') : '예상 산물 없음';
  function metricTable(headers, rows, caption) {
    const table = node('table'), head = node('thead'), body = node('tbody'), tr = node('tr');
    table.append(node('caption', caption));
    for (const h of headers) { const th = node('th', h); th.scope = 'col'; tr.append(th); }
    head.append(tr);
    for (const row of rows) {
      const tr = node('tr');
      row.forEach((text, i) => { const cell = node(i ? 'td' : 'th', text); if (!i) cell.scope = 'row'; tr.append(cell); });
      body.append(tr);
    }
    table.append(head, body); return table;
  }
  function chooseLens(lens, focus = false) {
    view().lens = lens;
    for (const name of REVIEW_LENSES) {
      const selected = lens === name, button = $(`review-tab-${name}`);
      button.setAttribute('aria-selected', String(selected)); button.tabIndex = selected ? 0 : -1;
      $(`review-panel-${name}`).hidden = !selected;
    }
    if (focus) $(`review-tab-${lens}`).focus();
  }
  $('review-tabs').hidden = false;
  for (const lens of REVIEW_LENSES) {
    const button = $(`review-tab-${lens}`);
    button.addEventListener('click', () => { chooseLens(lens); save(); });
    button.addEventListener('keydown', event => {
      const i = REVIEW_LENSES.indexOf(lens);
      const next = { ArrowRight: (i + 1) % 5, ArrowLeft: (i + 4) % 5, Home: 0, End: 4 }[event.key];
      if (next === undefined) return;
      event.preventDefault(); chooseLens(REVIEW_LENSES[next], true); save();
    });
  }
  $('review-design').addEventListener('change', () => {
    view().designId = $('review-design').value === 'draft' ? null : Number($('review-design').value);
    render(); save();
  });
  function sequenceLine(sequence, markEnd = false) {
    const p = node('p', '', 'pcr-sequence pcr-review-sequence');
    p.append('5′ ');
    if (markEnd) {
      const count = Math.min(5, sequence.length), end = node('span', sequence.slice(-count), 'pcr-review-end');
      end.title = `3′ 말단의 마지막 ${count} nt`; p.append(sequence.slice(0, -count), end);
    } else p.append(sequence);
    p.append(' 3′'); return p;
  }
  function alignmentBlock(analysis, run, first, second) {
    if (!run) return node('p', '이 단순 정렬에서 연속 상보 구간이 관찰되지 않습니다.');
    const lines = alignmentLines(analysis, run), block = node('div');
    const region = node('div', '', 'pcr-alignment-scroll'); region.tabIndex = 0;
    region.setAttribute('role', 'region'); region.setAttribute('aria-label', `${first}와 ${second}의 역평행 정렬. 가로로 스크롤하여 전체 서열 확인`);
    const pre = node('pre', '', 'pcr-structure pcr-review-alignment');
    function strand(prefix, text, suffix) {
      const start = lines.matchStart, end = start + run.length;
      pre.append(prefix, text.slice(0, start), node('span', text.slice(start, end), 'pcr-alignment-match'), text.slice(end), suffix);
    }
    // Both sequence lines have the same four-character prefix. Bottom is b reversed,
    // not reverse-complemented: the visible paired bases must be complementary.
    strand('5′  ', lines.top, '  3′\n');
    pre.append('    ', lines.bars, '\n');
    strand('3′  ', lines.bottom, '  5′'); region.append(pre);
    // Long alignments remain scrollable on screen and split at shared columns on paper.
    // Window fragments do not get invented 5′/3′ end labels.
    const columns = Math.max(lines.top.length, lines.bottom.length);
    if (columns > 60) {
      region.classList.add('pcr-long-alignment');
      const print = node('div', '', 'pcr-print-alignment');
      for (let start = 0; start < columns; start += 60) {
        const end = Math.min(start + 60, columns);
        print.append(node('p', `정렬 열 ${start + 1}~${end} / 위 ${first} 5′→3′ / 아래 ${second} 3′→5′`, 'pcr-small'),
          node('pre', `${lines.top.slice(start, end)}\n${lines.bars.slice(start, end)}\n${lines.bottom.slice(start, end)}`, 'pcr-structure'));
      }
      region.append(print);
    }
    block.append(node('p', `위: ${first} 5′→3′ / 아래: ${second} 3′→5′`, 'pcr-small'), region,
      node('p', `${run.length} nt 연속 상보 / ${first} 3′ ${run.aThreePrime ? '포함' : '미포함'} / ${second} 3′ ${run.bThreePrime ? '포함' : '미포함'} / 양쪽 3′ ${run.aThreePrime && run.bThreePrime ? '모두 포함' : '동시 포함 아님'}`, 'pcr-small'));
    return block;
  }
  function renderComplementarity(f, r) {
    const target = $('review-complementarity-values'); target.replaceChildren();
    for (const [id, title, a, b, first, second] of [
      ['self-f', 'Forward 자기 상보성', f, f, 'F', 'F 사본'],
      ['self-r', 'Reverse 자기 상보성', r, r, 'R', 'R 사본'],
      ['pair', 'F/R 사이 상보성', f, r, 'F', 'R']
    ]) {
      const analysis = complementarity(a, b), section = node('section', '', 'pcr-complementarity-row'); section.id = `review-${id}`;
      section.append(node('h4', title), node('p', `가장 긴 연속 상보 구간: ${analysis.longest?.length ?? 0} nt`), alignmentBlock(analysis, analysis.longest, first, second));
      const terminal = analysis.bothThreePrime || analysis.threePrime;
      if (terminal && terminal !== analysis.longest) {
        const details = node('details'); details.append(node('summary', `3′ 말단이 관여하는 별도 구간 / ${terminal.length} nt`), alignmentBlock(analysis, terminal, first, second)); section.append(details);
      }
      if (!analysis.threePrime) section.append(node('p', '3′ 말단을 포함하는 연속 상보 구간은 이 모형에서 관찰되지 않습니다.', 'pcr-small'));
      target.append(section);
    }
  }
  function renderMyDesign() {
    const state = getState(), fixture = getFixture(), model = reviewDesign(state, fixture);
    const select = $('review-design'), draft = node('option', '03의 현재 초안'); draft.value = 'draft';
    select.replaceChildren(draft);
    for (const design of state.designs) { const option = node('option', `저장 설계 ${design.id}`); option.value = String(design.id); select.append(option); }
    select.value = view().designId === null ? 'draft' : String(view().designId);
    for (const id of ['review-my-design', 'review-length-values', 'tm-values', 'review-end-values', 'review-complementarity-values']) $(id).replaceChildren();
    if (!model) {
      $('review-status').textContent = fixture ? '03에서 먼저 primer를 설계하세요. 유효한 F/R 초안 또는 저장 설계를 선택하면 여기에 표시됩니다.' : '교육 데이터가 준비되지 않아 설계를 검토할 수 없습니다. 답안 기록, JSON과 인쇄는 사용할 수 있습니다.';
      return;
    }
    $('review-status').textContent = model.designId === null ? '03의 현재 초안을 검토 중입니다. 03에서 수정하면 함께 갱신됩니다.' : `저장 설계 ${model.designId}을 검토 중입니다. 03의 초안과 저장된 원본은 변경하지 않습니다.`;
    const pair = node('div', '', 'pcr-review-pair');
    for (const name of ['F', 'R']) {
      const stats = model.primers[name].stats, section = node('div');
      section.append(node('h4', `${labels[name]} primer / ${name}`), sequenceLine(stats.sequence), node('p', `${stats.length} nt / GC ${Math.round(stats.gcPercent)}%`, 'pcr-small')); pair.append(section);
      const end = endFeatures(stats.sequence), endRow = node('div', '', 'pcr-end-row'); endRow.dataset.primer = name;
      endRow.append(node('h4', `${labels[name]} / ${name}`), sequenceLine(stats.sequence, true), node('p', `3′ end: ${end.sequence} / 마지막 염기: ${end.lastBase}`, 'pcr-sequence'), node('p', `표시한 ${end.sequence.length} nt의 GC ${end.gcCount}개 / 연속 동일 염기: ${end.runs.length ? end.runs.map(run => `${run.sequence} (${run.length} nt)`).join(', ') : '없음'}`, 'pcr-small'));
      $('review-end-values').append(endRow);
    }
    const products = node('div', '', 'pcr-review-products');
    for (const [source, result] of Object.entries(model.result)) products.append(node('p', `${source}  ${result.products.length ? [...new Set(result.products.map(p => p.length))].map(length => `${length} bp`).join(', ') : '예상 산물 없음'}`));
    $('review-my-design').append(pair, node('p', '완전 일치 모형의 예상 산물', 'pcr-small'), products);
    const f = model.primers.F.stats, r = model.primers.R.stats;
    $('review-length-values').append(metricTable(['관찰값', 'Forward / F', 'Reverse / R'], [['Length', `${f.length} nt`, `${r.length} nt`], ['GC', `${Math.round(f.gcPercent)}%`, `${Math.round(r.gcPercent)}%`]], '같은 primer pair의 길이와 GC'));
    const tm = node('div', '', 'pcr-tm-values');
    for (const [label, value] of [['Forward', f.simpleTm], ['Reverse', r.simpleTm], ['차이', Math.abs(f.simpleTm - r.simpleTm)]]) {
      const block = node('div'); block.append(node('h4', label), node('p', `${label === '차이' ? '약' : '간이 추정'} ${value} °C`)); tm.append(block);
    }
    $('tm-values').append(tm); renderComplementarity(f.sequence, r.sequence);
  }
  function renderCandidates() {
    const fixture = getFixture(), stage = view().stage;
    $('compare-ab').disabled = !fixture;
    $('compare-abc').disabled = !fixture || stage < 1;
    $('review-ab-question').hidden = stage < 1;
    $('review-c-question').hidden = stage < 2;
    $('review-candidate-status').textContent = !fixture ? '교육 데이터가 준비되지 않아 비교할 수 없습니다.' : stage === 0 ? '먼저 A와 B에서 비교하세요. 배경 C의 결과는 아직 공개하지 않았습니다.' : stage === 1 ? 'A와 B만 공개했습니다. 판단 근거를 기록한 뒤 배경 C를 확인하세요.' : 'A, B와 배경 C를 함께 공개했습니다.';
    $('candidate-results').replaceChildren(); $('candidate-sequence-values').replaceChildren(); $('review-p3-binding').replaceChildren();
    if (!fixture) return;
    for (const [id, pair] of Object.entries(fixture.candidatePairs)) {
      const seq = node('div', '', 'pcr-candidate-sequence');
      seq.append(node('h4', id), node('p', `F 5′ ${pair.forward} 3′\nR 5′ ${pair.reverse} 3′`, 'pcr-sequence'));
      $('candidate-sequence-values').append(seq);
    }
    if (!stage) return;
    const sources = stage === 1 ? ['A', 'B'] : ['A', 'B', 'C'];
    const scale = Math.max(...Object.values(fixture.templates).map(template => template.sequence.length));
    for (const [id, pair] of Object.entries(fixture.candidatePairs)) {
      const result = analyzeAll(fixture.templates, pair.forward, pair.reverse);
      const row = node('section', '', 'pcr-candidate-row'); row.dataset.candidate = id;
      row.append(node('h4', id));
      const lanes = node('div', '', 'pcr-candidate-lanes');
      for (const source of sources) {
        const lane = node('div'); lane.dataset.source = source;
        lane.append(node('h5', source), node('p', productsText(result[source].products), 'pcr-small'));
        for (const product of result[source].products) {
          const bar = node('div', '', 'pcr-product-bracket'); bar.style.width = `${100 * product.length / scale}%`;
          bar.setAttribute('aria-hidden', 'true'); lane.append(bar);
        }
        lanes.append(lane);
      }
      row.append(lanes); $('candidate-results').append(row);
      if (id === 'P3' && stage === 2) {
        const deletion = fixture.templates.B.deletion;
        const hits = result.A.hits.filter(hit => hit.primer === 'F');
        $('review-p3-binding').append(node('p', `A에서 P3 F 결합 위치: ${hits.map(hit => `${hit.start + 1}~${hit.end}`).join(', ')} / B의 결실: A 좌표 ${deletion.start}~${deletion.end}`, 'pcr-sequence'), node('p', `B에서 P3 F 완전 일치 결합: ${result.B.hits.filter(hit => hit.primer === 'F').length}개`, 'pcr-small'));
      }
    }
    $('candidate-results').append(node('p', '완전 일치 모형의 예상 산물입니다. 선은 같은 bp 축척으로 산물 길이를 나타냅니다. 실제 실험 결과가 아닙니다.', 'pcr-small'));
  }
  for (const [id, stage] of [['compare-ab', 1], ['compare-abc', 2]]) $(id).addEventListener('click', () => {
    if (!getFixture() || (stage === 2 && view().stage === 0)) return;
    view().stage = stage; $('review-p3-explanation').open = false; renderCandidates(); save();
  });
  function render() { chooseLens(view().lens); renderMyDesign(); renderCandidates(); }
  return { render };
}
