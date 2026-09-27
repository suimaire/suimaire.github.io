// Educational structures, not sequence analysis or amplification predictions.
export const RNA_OPTIONS = {
  activeSection: ['flow', 'genomic', 'strategy', 'control', 'transcripts', 'reflection'],
  flowStep: ['rna', 'rt', 'cdna', 'pcr'],
  primerStrategy: ['same-exon', 'intron-spanning', 'junction'],
  intronExample: ['short', 'long'],
  junctionPosition: ['exon', 'junction'],
  rtControlCase: ['case1', 'case2'],
  transcriptTarget: ['exon2', 'exon3', 'exon4', 'junction23', 'junction34']
};
export const RNA_ANSWERS = ['template', 'strategy', 'control', 'transcripts', 'specific', 'reflection'];
export const RNA_TRANSCRIPTS = [[1, 2, 3, 4], [1, 2, 4], [1, 3, 4]];
export function emptyRnaExtension() {
  return { activeSection: 'flow', flowStep: 'rna', primerStrategy: 'same-exon', intronExample: 'long', junctionPosition: 'exon', rtControlCase: 'case1', transcriptTarget: 'exon2', answers: {} };
}
export function validateRnaExtension(value, legacyAnswers = {}) {
  const result = emptyRnaExtension();
  if (value !== undefined) {
    if (!value || typeof value !== 'object' || Array.isArray(value)) throw new Error('RNA 확장 기록이 올바르지 않습니다.');
    if (Object.keys(value).some(key => key !== 'answers' && !Object.hasOwn(RNA_OPTIONS, key))) throw new Error('알 수 없는 RNA 확장 항목입니다.');
    for (const [key, options] of Object.entries(RNA_OPTIONS)) {
      if (value[key] !== undefined && !options.includes(value[key])) throw new Error('RNA 확장 선택값이 올바르지 않습니다.');
      if (value[key] !== undefined) result[key] = value[key];
    }
    if (value.answers !== undefined) {
      if (!value.answers || typeof value.answers !== 'object' || Array.isArray(value.answers)) throw new Error('RNA 답안 구조가 올바르지 않습니다.');
      for (const [key, answer] of Object.entries(value.answers)) {
        if (!RNA_ANSWERS.includes(key) || typeof answer !== 'string' || answer.length > 12000) throw new Error('RNA 답안 항목 또는 길이가 올바르지 않습니다.');
        if (key === 'template' && !['', 'rna', 'cdna', 'protein', 'dntp'].includes(answer)) throw new Error('RNA template 답안이 올바르지 않습니다.');
        result.answers[key] = answer;
      }
    }
  }
  // Copy only on the first migration. The old answer remains verbatim in answers.
  if (value === undefined && legacyAnswers['rna-plan'] !== undefined) result.answers.reflection = legacyAnswers['rna-plan'];
  return result;
}
export function strategyModel(strategy, intron = 'long', position = 'junction') {
  if (strategy === 'same-exon') return {
    placement: 'F와 R 모두 Exon 1 내부에서 서로 마주봅니다.',
    cdna: '두 primer 모두 결합 가능한 구조입니다. Exon 1 내부 product가 가능합니다.',
    genomic: '같은 Exon 1 결합 부위가 있어 동일한 크기의 product가 가능합니다.',
    meaning: '이 배치만으로는 cDNA 유래 증폭과 gDNA 유래 증폭을 구별하기 어려울 수 있습니다. gDNA 분석이 목적이라면 이 배치 자체가 문제가 되지는 않습니다.'
  };
  if (strategy === 'intron-spanning' || (strategy === 'junction' && position === 'exon')) return {
    placement: 'F는 Exon 1, R은 Exon 2에 놓여 서로 마주봅니다.',
    cdna: 'Exon 1과 Exon 2가 연결되어 상대적으로 짧은 product가 가능합니다.',
    genomic: intron === 'long' ? '긴 intron을 포함하는 훨씬 긴 product가 가능합니다. 현재 PCR 조건에서 증폭이 불리할 수 있습니다.' : '짧은 intron을 포함하는 더 긴 product가 가능합니다. gDNA도 함께 증폭될 수 있습니다.',
    meaning: strategy === 'junction' ? 'F는 아직 Exon 1 내부에 있어 두 template 모두에 결합 가능한 구조입니다. R은 Exon 2에 그대로 두고 F를 junction으로 옮겨 연속 결합 부위를 비교하세요.' : '충분히 긴 intron을 사이에 두면 현재 PCR 조건에서 gDNA product의 증폭이 불리해지거나 cDNA product와 크기로 구별하는 데 도움이 될 수 있습니다. 실제 결과는 intron 길이, extension time, polymerase와 PCR 조건에 따라 달라집니다.'
  };
  return {
    placement: 'F 하나가 Exon 1–2 junction을 가로지르고 R은 Exon 2에서 마주봅니다.',
    cdna: 'Splicing된 Exon 1–2에 연속된 F 결합 서열이 존재하며 R도 결합 가능한 구조입니다.',
    genomic: '이 intron을 포함하는 gDNA에는 동일한 연속 F 결합 부위가 없습니다. 두 exon 부분은 intron으로 떨어져 있습니다.',
    meaning: 'Junction을 가로지르는 배치는 cDNA를 구별하는 데 도움이 됩니다. 실제 특이성은 primer 서열과 genome context를 함께 검토해야 합니다. 유사 유전자나 가공 위유전자 등에서는 다른 결합 가능성이 있어 gDNA 증폭 배제를 보장하지 않습니다.'
  };
}
export function transcriptBinding(exons, target) {
  const region = { exon2: [2], exon3: [3], exon4: [4], junction23: [2, 3], junction34: [3, 4] }[target];
  if (!region) throw new Error('알 수 없는 transcript 영역입니다.');
  const present = region.length === 1 ? exons.includes(region[0]) : exons.some((e, i) => e === region[0] && exons[i + 1] === region[1]);
  const label = region.length === 1 ? `Exon ${region[0]}` : `Exon ${region.join('–')} junction`;
  return { present, region, label, reason: present ? `${label}이 있어 선택한 F와 공통 Exon 4의 R이 결합 가능한 구조입니다.` : `${label}이 없어 선택한 F의 결합 구조가 없습니다.` };
}
