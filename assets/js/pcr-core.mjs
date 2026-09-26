// Display coordinates: 1-based inclusive. Internal intervals: 0-based half-open.
export const LIMITS = Object.freeze({ template: 5000, primer: 100, raw: 30000, pairs: 100000 });
export function normalizeSequence(raw, max = LIMITS.primer) {
  if (typeof raw !== 'string' || raw.length > LIMITS.raw) throw new Error(`입력은 ${LIMITS.raw}자 이하여야 합니다.`);
  const sequence = raw.replace(/\s/g, '').toUpperCase();
  if (!sequence.length) throw new Error('서열을 입력하세요. 빈 입력은 계산하지 않습니다.');
  if (sequence.length > max) throw new Error(`서열은 ${max} nt까지 지원합니다.`);
  const invalid = [...sequence].flatMap((base, i) => /[ACGT]/.test(base) ? [] : [`${i + 1}번 ${base}`]);
  if (invalid.length) throw new Error(`공백을 제외한 서열 위치: ${invalid.slice(0, 8).join(', ')}. A/C/G/T만 지원하며 FASTA 헤더와 N은 지원하지 않습니다.`);
  return sequence;
}
export function complement(raw) {
  return [...normalizeSequence(raw, LIMITS.template)].map(b => ({ A: 'T', T: 'A', C: 'G', G: 'C' })[b]).join('');
}
export function reverseComplement(raw) { return [...complement(raw)].reverse().join(''); }
export function primerStats(raw) {
  const sequence = normalizeSequence(raw);
  const gc = [...sequence].filter(b => b === 'G' || b === 'C').length;
  return { sequence, length: sequence.length, gcPercent: 100 * gc / sequence.length, threePrime: sequence.slice(-5), simpleTm: 2 * (sequence.length - gc) + 4 * gc };
}
export function toInterval(start, end, length) {
  if (!Number.isInteger(start) || !Number.isInteger(end) || start < 1 || end < start || end > length) {
    throw new Error(`좌표는 1~${length} 범위의 정수이고 시작 ≤ 끝이어야 합니다.`);
  }
  return { start: start - 1, end };
}
export function primerFromRange(raw, start, end, direction) {
  const sequence = normalizeSequence(raw, LIMITS.template);
  const range = toInterval(start, end, sequence.length);
  if (!['right', 'left'].includes(direction)) throw new Error('신장 방향을 선택하세요.');
  return normalizeSequence(direction === 'right' ? sequence.slice(range.start, range.end) : reverseComplement(sequence.slice(range.start, range.end)));
}
function positions(sequence, needle) {
  const found = [];
  for (let from = 0, at; (at = sequence.indexOf(needle, from)) !== -1; from = at + 1) found.push(at);
  return found;
}
export function bindingSites(raw, primer, label = 'F') {
  const sequence = normalizeSequence(raw, LIMITS.template), p = normalizeSequence(primer);
  return ['right', 'left'].flatMap(direction => positions(sequence, direction === 'right' ? p : reverseComplement(p))
    .map(start => ({ primer: label, start, end: start + p.length, direction, strand: direction === 'right' ? 'lower' : 'upper' })));
}
export function productFromHits(raw, right, left) {
  const sequence = normalizeSequence(raw, LIMITS.template);
  for (const hit of [right, left]) toInterval(hit.start + 1, hit.end, sequence.length);
  if (right.direction !== 'right' || left.direction !== 'left') throw new Error('안쪽을 향하는 반대 가닥의 결합 쌍이 필요합니다.');
  if (right.end > left.start) throw new Error('겹치거나 바깥쪽을 향한 결합 구간은 이 모형의 지원 범위 밖입니다.');
  return { start: right.start, end: left.end, length: left.end - right.start, sequence: sequence.slice(right.start, left.end) };
}
export function analyzeTemplate(raw, forward, reverse, source = '') {
  const sequence = normalizeSequence(raw, LIMITS.template);
  const hits = [...bindingSites(sequence, forward, 'F'), ...bindingSites(sequence, reverse, 'R')];
  const rights = hits.filter(h => h.direction === 'right'), lefts = hits.filter(h => h.direction === 'left');
  if (rights.length * lefts.length > LIMITS.pairs) throw new Error('반복 결합 조합이 100,000개를 넘습니다. 더 긴 프라이머로 범위를 좁히세요. 일부 결과만 표시하지 않습니다.');
  const unique = new Map();
  let overlappingPairs = 0;
  for (const right of rights) for (const left of lefts) {
    if (right.end > left.start) { if (right.start < left.end && left.start < right.end) overlappingPairs++; continue; }
    const key = `${right.start}:${left.end}`;
    if (!unique.has(key)) unique.set(key, { ...productFromHits(sequence, right, left), source, bindings: [] });
    unique.get(key).bindings.push({ right: { ...right }, left: { ...left } });
  }
  return { hits, overlappingPairs, products: [...unique.values()].sort((a, b) => a.start - b.start || a.end - b.end) };
}
export function analyzeAll(templates, forward, reverse) {
  return Object.fromEntries(Object.entries(templates).map(([id, data]) => [id, analyzeTemplate(data.sequence, forward, reverse, id)]));
}
// Equal-length bands may overlap visually; each distinct source and interval stays in the list.
export function combineProducts(...lists) {
  const unique = new Map();
  for (const p of lists.flat()) {
    const key = `${p.source}:${p.start}:${p.end}`;
    if (!unique.has(key)) unique.set(key, { ...p, bindings: [...p.bindings] });
  }
  return [...unique.values()];
}
export function gelPosition(bp, min = 20, max = 5000) {
  if (!Number.isFinite(bp) || bp < min || bp > max) throw new Error('전기영동 표시 범위는 20~5000 bp입니다.');
  return 30 + 220 * (Math.log10(max) - Math.log10(bp)) / (Math.log10(max) - Math.log10(min));
}
