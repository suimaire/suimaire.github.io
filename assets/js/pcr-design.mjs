import { analyzeAll, bindingSites, normalizeSequence, primerStats, productFromHits, reverseComplement, toInterval, LIMITS } from './pcr-core.mjs';

export const PRIMER_KEYS = { F: 'forward', R: 'reverse' };
export const emptyWorkbench = () => ({ mode: 'F', windowStart: 1, showPrediction: false });
export const cloneBindings = bindings => Object.fromEntries(['F', 'R'].map(name => [name, bindings?.[name] ? { ...bindings[name] } : null]));

// Keep coordinates as entered, including incomplete and invalid edits. Never clamp a design.
export function sequenceFromBinding(sequence, binding) {
  if (!binding) return '';
  const interval = toInterval(Number(binding.start), Number(binding.end), sequence.length);
  const selected = sequence.slice(interval.start, interval.end);
  return binding.direction === 'left' ? reverseComplement(selected) : selected;
}

export function inspectPrimer(draft, name, fixture) {
  const raw = draft[PRIMER_KEYS[name]], explicit = draft.bindings?.[name];
  const info = { name, raw, binding: explicit || null, inferred: false, hits: [], stats: null, reference: '', error: '' };
  try {
    if (explicit) {
      const interval = toInterval(Number(explicit.start), Number(explicit.end), fixture.templates.A.sequence.length);
      info.reference = fixture.templates.A.sequence.slice(interval.start, interval.end);
      if (sequenceFromBinding(fixture.templates.A.sequence, explicit) !== raw.replace(/\s/g, '').toUpperCase()) throw new Error('좌표와 주문 서열이 일치하지 않습니다. 좌표를 다시 적용하거나 서열을 직접 수정하세요.');
    }
    if (!raw && !explicit) return info;
    const sequence = normalizeSequence(raw, LIMITS.template);
    const gc = [...sequence].filter(base => base === 'G' || base === 'C').length;
    info.stats = { sequence, length: sequence.length, gcPercent: gc * 100 / sequence.length };
    if (sequence.length > LIMITS.primer) throw new Error(`선택은 유지했습니다. ${sequence.length} nt는 계산 한도 ${LIMITS.primer} nt를 넘습니다.`);
    info.stats = primerStats(sequence);
    info.hits = bindingSites(fixture.templates.A.sequence, sequence, name);
    if (!explicit && info.hits.length === 1) {
      const hit = info.hits[0];
      info.binding = { start: String(hit.start + 1), end: String(hit.end), direction: hit.direction };
      info.inferred = true;
      info.reference = fixture.templates.A.sequence.slice(hit.start, hit.end);
    }
  } catch (error) { info.error = error.message; }
  return info;
}

export function inspectDesign(draft, fixture) {
  const primers = Object.fromEntries(['F', 'R'].map(name => [name, inspectPrimer(draft, name, fixture)]));
  const errors = Object.values(primers).filter(p => p.error).map(p => `${p.name}: ${p.error}`);
  const deletion = fixture.templates.B.deletion;
  const model = { primers, result: null, selectedProduct: null, placement: '', deletionMessage: '', bindingMessage: '', errors };
  if (!errors.length && primers.F.stats && primers.R.stats) {
    try { model.result = analyzeAll(fixture.templates, draft.forward, draft.reverse); }
    catch (error) { errors.push(error.message); }
  }
  const pair = Object.values(primers).map(p => p.binding && !p.error ? { start: Number(p.binding.start) - 1, end: Number(p.binding.end), direction: p.binding.direction } : null);
  if (pair.every(Boolean)) {
    const [f, r] = pair;
    if (f.start < r.end && r.start < f.end) model.placement = 'F와 R 결합 영역이 겹칩니다. 이 모형에서는 겹친 영역으로 하나의 amplicon을 계산하지 않습니다.';
    else {
      const right = pair.find(p => p.direction === 'right'), left = pair.find(p => p.direction === 'left');
      try {
        if (!right || !left) throw new Error('orientation');
        model.selectedProduct = productFromHits(fixture.templates.A.sequence, right, left);
      } catch { model.placement = '현재 두 primer는 하나의 inward-facing amplicon을 정의하지 않습니다. 두 3′ 말단과 합성 방향을 확인하세요.'; }
    }
    const [first, second] = [...pair].sort((a, b) => a.start - b.start);
    const between = first.end < deletion.start && second.start >= deletion.end;
    const overlapping = ['F', 'R'].filter((name, i) => pair[i].start < deletion.end && pair[i].end >= deletion.start);
    if (overlapping.length) {
      model.deletionMessage = `${overlapping.join('/') } 결합 부위가 ${deletion.start}~${deletion.end} 결실 영역과 겹칩니다. 결실을 두 primer 사이에 두는 배치와 다릅니다.`;
      model.bindingMessage = overlapping.map(name => {
        const full = name === 'F' ? 'Forward' : 'Reverse';
        const alternatives = model.result?.B.hits.filter(h => h.primer === name).length || 0;
        return `B에서는 ${full} primer의 선택한 결합 부위 일부 또는 전부가 결실됩니다. ${alternatives ? `다른 위치의 완전 일치 결합 ${alternatives}개는 별도로 계산합니다.` : '해당 primer의 완전 일치 결합 부위가 존재하지 않습니다.'}`;
      }).join(' ');
    } else if (between) {
      model.deletionMessage = `두 primer 사이에 ${deletion.start}~${deletion.end} 결실 영역이 포함되어 있습니다.${model.selectedProduct ? '' : ' 다만 현재 합성 방향으로는 이 구간의 산물을 정의할 수 없습니다.'}`;
    } else {
      model.deletionMessage = `${deletion.start}~${deletion.end} 결실 영역이 두 primer 사이에 포함되지 않습니다. 현재 선택한 두 위치에서는 A와 B의 ${deletion.end - deletion.start + 1} bp 차이가 산물 길이 차이로 나타나지 않습니다.`;
    }
  } else if (!errors.length) {
    model.placement = primers.F.stats && primers.R.stats ? '주문 서열의 결합 위치가 없거나 여러 개입니다. 아래 실제 완전 일치 결과를 확인하세요.' : 'F와 R 결합 부위를 모두 선택하면 예상 amplicon을 표시합니다.';
  }
  return model;
}
