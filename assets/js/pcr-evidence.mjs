import { gelPosition } from './pcr-core.mjs';

// Fixed teaching observations, independent of every student sequence and calculation.
export const EVIDENCE_LANES = {
  marker: { label: 'M', description: 'DNA size marker / 길이를 비교하는 DNA ladder' },
  sample: { label: 'Sample', description: '해석하려는 시료' },
  positive: { label: 'Positive control', description: '해당 primer pair와 PCR 조건이 정상적으로 작동할 때 예상 band가 나와야 하는 대조군' },
  ntc: { label: 'NTC', description: 'No-template control / template DNA를 넣지 않은 대조군' }
};
const marker = [500, 400, 300, 200, 100];
const defineCase = (number, purpose, bands, observations, questions, explanation) => ({
  id: `case-${number}`, label: `상황 ${number}`, additional: number === 4, purpose,
  expectedSize: 260,
  expectedControlBehavior: { marker: '기준 band 5개', positive: '약 260 bp band 하나', ntc: 'band 없음' },
  lanes: Object.fromEntries(Object.entries(EVIDENCE_LANES).map(([id, lane]) => [id, {
    ...lane, bands: id === 'marker' ? [...marker] : bands[id],
    observation: id === 'marker' ? '500, 400, 300, 200, 100 bp 기준 band가 보임' : observations[id]
  }])), questions, explanation
});
export const EVIDENCE_CASES = [
  defineCase(1, '산물 길이에 대한 증거와 sequence identity를 구별한다.',
    { sample: [260], positive: [260], ntc: [] },
    { sample: '예상 크기 부근에 band 하나가 보임 (약 260 bp)', positive: '예상 크기 부근에 band 하나가 보임 (약 260 bp)', ntc: 'band가 보이지 않음' },
    ['이 결과는 이 사례의 계산된 예상과 일치하는가?', 'Sample에 band가 하나이고 예상 크기와 비슷하다는 사실만으로 그 DNA가 목표 sequence라고 확정할 수 있는가?'],
    '예상 크기와의 일치는 target product라는 해석을 지지하는 증거가 될 수 있습니다. 그러나 같은 크기의 다른 산물이나 비특이적 산물 가능성을 gel의 크기 정보만으로 완전히 배제할 수는 없습니다. Band size와 molecular identity는 같은 정보가 아닙니다.'),
  defineCase(2, 'NTC 관찰과 가능한 원인, 추가 확인을 구분한다.',
    { sample: [260], positive: [260], ntc: [70] },
    { sample: '예상 크기 부근에 band 하나가 보임 (약 260 bp)', positive: '예상 크기 부근에 band 하나가 보임 (약 260 bp)', ntc: '100 bp 기준보다 아래쪽에 짧은 band 하나가 보임' },
    ['NTC에 band가 나타났습니다. 이 관찰 하나만으로 contamination과 primer-dimer 중 하나를 확정할 수 있는가?', '현재 자료에서 가능한 설명을 적고, 추가로 무엇을 확인해야 하는지 기록하세요.'],
    'NTC band는 정상적인 음성 대조 결과가 아닙니다. Contamination, primer-derived product, 비특이적 amplification 또는 기타 실험 문제를 고려할 수 있습니다. 작은 band는 primer-dimer 가능성을 생각하게 하지만, 위치 하나만으로 원인이나 molecular identity를 확정하지는 않습니다. 시약과 작업 과정, 반복 대조 반응 및 산물의 성격을 추가로 확인할 수 있습니다.'),
  defineCase(3, 'Positive control이 증폭되지 않을 때 Sample 음성 해석의 제한을 찾는다.',
    { sample: [], positive: [], ntc: [] },
    { sample: 'band가 보이지 않음', positive: 'band가 보이지 않음', ntc: 'band가 보이지 않음' },
    ['Sample에 band가 없습니다. 이 시료가 target-negative라고 결론 내릴 수 있는가?', '어떤 control의 결과가 정상이어야 Sample의 음성 결과를 해석할 수 있는가? Positive control lane을 직접 선택해 확인하세요.'],
    'Positive control이 정상적으로 증폭되지 않았다면 PCR reaction, reagent, thermal cycling, primer performance 등 여러 문제를 고려해야 합니다. Sample의 band absence를 곧바로 target absence로 해석할 수 없습니다. Positive control의 예상 band와 NTC의 band 부재를 함께 확인하고, 필요하면 시료별 추출 상태와 반응 억제를 살피는 내부 대조도 검토합니다. 이 자료로 원인을 하나로 특정할 수는 없습니다.'),
  defineCase(4, '하나의 desired product 예상과 여러 관찰 산물의 차이를 해석한다.',
    { sample: [420, 260, 140], positive: [260], ntc: [] },
    { sample: '예상 위치 부근 외에도 band가 보여 총 세 개임 (약 420, 260, 140 bp)', positive: '예상 크기 부근에 band 하나가 보임 (약 260 bp)', ntc: 'band가 보이지 않음' },
    ['하나의 desired product를 예상했는데 여러 band가 보입니다. 어떤 종류의 문제를 시사할 수 있는가?', '각 band의 sequence identity에 관해 현재 gel만으로 무엇을 말할 수 있는가?'],
    '여러 band는 비특이적 amplification 가능성을 고려하게 합니다. 계산에서 하나의 desired product를 예상했더라도 관찰에서는 여러 DNA product가 나타날 수 있습니다. 각 band의 sequence identity는 gel만 보고 확정할 수 없습니다.')
];
export const emptyEvidence = () => ({ activeCase: 'case-1', selectedLanes: Object.fromEntries(EVIDENCE_CASES.map(c => [c.id, null])) });
// Reuse the core's log scale with an explicit 50–500 bp teaching window.
// Output is a percentage of the 280-unit lane, not a measured migration distance.
export const evidenceBandPosition = bp => gelPosition(bp, 50, 500) / 280 * 100;
export const caseSummary = c => Object.values(c.lanes).map(lane => `${lane.label}: ${lane.observation}`).join(' / ');
