import { normalizeSequence, primerStats, reverseComplement } from './pcr-core.mjs';
import { inspectDesign } from './pcr-design.mjs';

export const REVIEW_LENSES = ['length-gc', 'tm', 'end', 'complementarity', 'off-target'];
export const emptyReview = () => ({ lens: 'length-gc', designId: null, stage: 0 });

// Read the explicitly selected snapshot, otherwise the current draft. Never change either.
export function reviewDesign(state, fixture) {
  const designId = state.review?.designId ?? null;
  const design = designId === null ? state.draft : state.designs.find(d => d.id === designId);
  if (!fixture || !design) return null;
  const model = inspectDesign(design, fixture);
  if (!model.result || model.errors.length) return null;
  return { designId, primers: model.primers, result: model.result };
}

export function endFeatures(raw) {
  const { threePrime } = primerStats(raw);
  return {
    sequence: threePrime, lastBase: threePrime.slice(-1),
    gcCount: [...threePrime].filter(base => /[GC]/.test(base)).length,
    runs: [...threePrime.matchAll(/(A{2,}|C{2,}|G{2,}|T{2,})/g)].map(m => ({ sequence: m[0], start: m.index + 1, length: m[0].length }))
  };
}

// Ungapped, antiparallel sequence observation only. No energy model or hairpin geometry.
// offset is the start of reverseComplement(b) in a's 0-based coordinate system.
// bRCStart=0 is b's physical 3′ end, not b's 5′ end.
export function complementarity(rawA, rawB = rawA) {
  const a = normalizeSequence(rawA), b = normalizeSequence(rawB), rc = reverseComplement(b);
  const runs = [];
  for (let offset = 1 - b.length; offset < a.length; offset++) {
    let start = null;
    const left = Math.max(0, offset), right = Math.min(a.length, offset + b.length);
    const finish = end => {
      if (start === null) return;
      runs.push({ offset, aStart: start, bRCStart: start - offset, length: end - start,
        aThreePrime: end === a.length, bThreePrime: start === offset });
      start = null;
    };
    for (let i = left; i < right; i++) {
      if (a[i] === rc[i - offset]) { if (start === null) start = i; }
      else finish(i);
    }
    finish(right);
  }
  const ends = run => Number(run.aThreePrime) + Number(run.bThreePrime);
  runs.sort((x, y) => y.length - x.length || ends(y) - ends(x) || Math.abs(x.offset) - Math.abs(y.offset) || x.offset - y.offset || x.aStart - y.aStart);
  return { a, b, longest: runs[0] ?? null,
    threePrime: runs.find(r => r.aThreePrime || r.bThreePrime) ?? null,
    bothThreePrime: runs.find(r => r.aThreePrime && r.bThreePrime) ?? null };
}

export function alignmentLines(analysis, run = analysis.longest) {
  if (!run) return null;
  const origin = Math.min(0, run.offset), topPad = Math.max(0, -origin), bottomPad = run.offset - origin;
  return {
    top: ' '.repeat(topPad) + analysis.a,
    bars: ' '.repeat(run.aStart - origin) + '|'.repeat(run.length),
    bottom: ' '.repeat(bottomPad) + [...analysis.b].reverse().join(''),
    topPad, bottomPad, matchStart: run.aStart - origin
  };
}
