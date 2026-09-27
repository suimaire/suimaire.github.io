# PCR 프라이머 디자인 Phase 4 검증 보고서

검증일: 2026-09-27 (KST). 기준 commit: `deb572c0897a5155928045978036054206c9fb9e`. 작업 branch: `feature/pcr-primer-design`.

활동 05 ‘예상 산물과 실험 증거 구별하기’를 선택 가능한 gel viewer 중심으로 재설계했다. 학생의 계산된 예상과 고정 수업용 가상 관찰을 분리하고, 관찰 / 해석 / 아직 확정할 수 없는 점을 각각 기록한다. 자동 채점이나 원인 확정 표시는 없다.

## 1. 시작 상태와 보호 범위

- 시작 HEAD는 요청한 Phase 3 commit `deb572c`였다. 작업 트리에는 미추적 `_codex/`만 있었다.
- 최초 원격 조회는 sandbox 네트워크 제한으로 연결되지 않았다. 허용된 네트워크 실행으로 재확인한 `origin/feature/pcr-primer-design`의 실제 원격 HEAD도 `deb572c0897a5155928045978036054206c9fb9e`였다. Phase 3은 이미 push된 상태였다.
- 저장소 및 적용되는 상위 디렉터리에서 `AGENTS.md`를 찾지 못했다. Phase 3 보고서, draft/snapshot, review 선택, core, fixture, 저장/인쇄, 기존 테스트를 확인했다.
- `_codex/`와 사용자 개인 파일은 수정하거나 staging하지 않았다. `main` 병합 또는 push는 수행하지 않는다.
- 기준 commit과 비교하여 활동 05 바깥의 HTML 전체가 동일함을 검사했다. 00~04, 06~07, RNA의 문항과 핵심 상호작용을 재설계하지 않았다.
- `pcr-core.mjs`, `pcr-design.mjs`, `pcr-workbench.mjs`, `pcr-review.mjs`, `pcr-review-view.mjs`, `pcr-intro.mjs`, 기존 JSON fixture는 기준 commit과 동일하다.

## 2. 수정 파일

| 파일 | 변경 |
| --- | --- |
| `bioinformatics/pcr-primer-design.html` | 활동 05 예상/관찰 분리, viewer, 종합 질문, legacy/no-JS 안내 |
| `assets/css/pcr-worksheet.css` | 05 전용 gel/lane, 반응형, 인쇄 스타일 |
| `assets/js/pcr-evidence.mjs` | 고정 사례, lane 정의, log 위치 계산 연결, 보기 기본값 |
| `assets/js/pcr-evidence-view.mjs` | 05 DOM 생성, 예상 읽기, 사례/lane 표시, 답안과 인쇄 연결 |
| `assets/js/pcr-records.mjs` | 선택적 v1 `evidence` 기본값 및 검증 |
| `assets/js/pcr-worksheet.mjs` | 옛 gel 렌더러 제거, 05 갱신/저장/복원/인쇄 연결, 07 해석 요약 연결 |
| `tests/pcr-evidence.test.mjs` | 새 단위 테스트 6개 |
| `tests/pcr-evidence-browser.mjs` | 새 브라우저/접근성/저장/인쇄 검사, PNG/contact sheet 생성 |
| `tests/pcr-browser.mjs` | 옛 혼합 gel 검사를 계산된 예상 검사로 변경, no-JS 선택자 범위 명시 |
| `tests/pcr-static.test.mjs` | 새 모듈의 금지 문자와 안전한 DOM 검사 포함 |
| `tests/PCR_PRIMER_DESIGN_PHASE4.md` | 이 보고서 |

## 3. 기존 gel 조사

기존 `pcr-worksheet.mjs`의 `renderGel(result)`는 학생 초안의 계산 결과로 A, B, C, A+C, B+C lane과 산물 표를 만들었다. `combineProducts`는 출처와 구간을 보존했고 같은 크기의 그림 band는 겹쳐 표시했다. 그 아래 세 가상 관찰은 별도 HTML 표와 통합 textarea에 있었다.

기존 `pcr-core.mjs`의 실제 함수는 다음 식이다.

```js
gelPosition(bp, min = 20, max = 5000)
// 30 + 220 * (log10(max) - log10(bp)) / (log10(max) - log10(min))
```

이번에는 옛 예상 gel과 표를 제거하고 A/B/C의 계산된 예상 목록을 위쪽에 둔다. 핵심 core 함수와 `combineProducts` 자체는 변경하지 않았다. 고정 가상 관찰 gel을 학생의 product 계산에서 만들지 않는다.

## 4. 학생 설계 연결

05는 04의 기존 읽기 함수 `reviewDesign(state, fixture)`를 그대로 사용한다. `review.designId`가 명시되면 그 저장 snapshot을, `null`이면 03의 현재 draft를 읽는다. `inspectDesign`의 실제 유효성 확인과 `analyzeAll` 결과에서 A/B/C 산물 길이를 가져온다. 학생 결과를 hard-code하지 않았다.

- 03 편집/불러오기, 04 설계 선택, JSON import와 reset에 따라 예상 목록을 갱신한다.
- snapshot을 선택한 상태에서 초안을 수정해도 해당 snapshot의 예상이 유지된다.
- 유효한 설계가 없으면 안내를 표시한다. 좌표/서열 불일치나 잘못된 염기가 있으면 이전 계산값을 남기지 않는다.
- 설계 부재 또는 기존 JSON fixture fetch 실패 때도 고정 가상 관찰 활동을 사용할 수 있다.
- 활동 05는 draft, bindings, designs, review 및 04 answers를 변경하지 않는다.

P1/P2/P3 고정 결과도 기존 core 테스트로 보호했다: P1 A/B/C = 260/180/없음, P2 = 260/180/220, P3 = 180/없음 (B)/없음 (C). 모든 단위는 bp다.

## 5. 고정 가상 관찰 데이터

`pcr-evidence.mjs`의 `EVIDENCE_CASES`가 id, label, additional, purpose, expectedSize, expectedControlBehavior, lanes, questions, explanation을 관리한다. 각 lane에는 label, description, bands, observation이 있다. 사용자 설계나 저장 파일에서 band를 주입하지 않는다.

모든 사례의 독립된 계산 예상은 260 bp 산물 하나다. M은 500/400/300/200/100 bp 기준 band다. 아래 수치는 그림 배치를 위한 고정 교육 값이며 실제 측정값이 아니다. 화면의 관찰 표현은 ‘약’ 또는 기준 band와의 상대 위치를 사용한다.

| 상황 | Sample | Positive control | NTC | 목적 |
| --- | --- | --- | --- | --- |
| 1 | 260 | 260 | 없음 | 길이 일치와 sequence identity 구별 |
| 2 | 260 | 260 | 70 | NTC signal의 복수 설명과 추가 확인 |
| 3 | 없음 | 없음 | 없음 | positive control 부재에 따른 음성 해석 제한 |
| 4 / 추가 사례 | 420, 260, 140 | 260 | 없음 | 여러 산물의 가능성과 identity 불확실성 |

상황별 expectedControlBehavior에는 marker 기준 band, positive 약 260 bp 하나, NTC band 없음이 명시된다. 이 기대와 다른 관찰을 오류 판정이나 자동 점수로 변환하지 않는다.

## 6. band 위치 계산과 한계

`evidenceBandPosition(bp)`는 기존 `gelPosition(bp, 50, 500) / 280 * 100`을 호출한다. 반환값은 CSS lane 높이에 대한 백분율이다. 작은 fragment일수록 아래로 이동하며 같은 길이의 band는 lane 간 같은 높이에 그려진다. 두 배 길이 차이는 log 축에서 같은 위치 차이다.

이는 50~500 bp의 개념적 표시 범위다. Agarose percentage, voltage, running time, buffer, DNA conformation을 모델링하지 않는다. 일정한 5 px 실선 두께와 같은 색을 사용하며 수율/초기 DNA 양의 정량값을 저장하거나 표현하지 않는다. 화면과 인쇄 설명에 위치 모형과 밝기 해석의 한계를 명시한다.

## 7. 상황별 추론

상황 1은 예상 크기의 단일 band와 정상 기대에 맞는 controls를 제시한다. 크기의 일치는 해석을 지지할 수 있지만 같은 길이의 다른 산물을 배제하지 못한다는 설명은 접힌 해설에 있다.

상황 2는 NTC의 짧은 band를 직접 선택하게 한다. Contamination, primer-derived product, 비특이적 amplification, 기타 실험 문제를 가능한 설명으로 제시한다. 작은 위치를 primer-dimer의 확정 근거로 표시하지 않는다.

상황 3은 Sample과 Positive control, NTC 모두 band가 없는 관찰이다. Positive control lane을 선택하는 질문을 제공한다. PCR reaction, reagent, thermal cycling, primer performance 등의 가능성과 Sample 음성 해석의 제한을 해설하되 특정 원인을 단정하지 않는다.

상황 4는 기본적으로 접힌 ‘추가 사례’에서 연다. Sample에 세 band가 있고 controls는 기대에 맞는다. 비특이적 증폭 가능성과 각 산물의 정체가 미확인이라는 점을 구분한다.

## 8. Lane interaction과 답안

Gel은 HTML/CSS로 만든다. lane 전체가 native button이어서 마우스, 터치, Tab과 Enter/Space로 선택한다. 얇은 청록 경계, ‘선택됨’ 텍스트, `aria-pressed`를 함께 갱신한다. 아래 live region에는 lane 이름, 대조군 역할, 관찰을 표시한다. 자동 정답/오답은 없다.

각 사례에는 관찰 input 한 줄, 해석 textarea 2줄, 아직 확정할 수 없는 점 textarea 2줄이 있다. 사례 전환 시 답안 DOM을 교체하지 않아 입력과 focus를 불필요하게 잃지 않는다. 해설은 수동으로 펼친다. 마지막 두 종합 문항은 각각 3줄이며 기존 `evidence-identity`, `evidence-controls` 키를 유지한다.

마지막 질문은 band의 sequence identity와 다른 genomic/transcript sequence에 대한 specificity를 더 살필 방법을 묻는다. 일반적인 활동 06 이동 링크를 사용한다. 07은 새 상황별 해석을 기존 출처 요약에 읽기 전용으로 포함하며, 새 해석이 없으면 이전 통합 기록을 사용한다.

## 9. State, JSON, localStorage

`schemaVersion: 1`, `dataVersion`, localStorage 키 `hafs:pcr-primer:v1`을 유지했다. 선택적 최상위 객체를 추가한다.

```json
"evidence": {
  "activeCase": "case-1",
  "selectedLanes": {
    "case-1": null, "case-2": null,
    "case-3": null, "case-4": null
  }
}
```

허용 lane은 marker/sample/positive/ntc 또는 null이다. 알려지지 않은 case/lane과 잘못된 객체는 import에서 거부한다. 객체가 없는 old v1은 안전한 기본값을 받는다.

12개 새 답안은 기존 `answers`에 `evidence-case-{1~4}-{observation|interpretation|uncertainty}`로 저장한다. 새 답안 필드를 기존 수집 절차 전에 생성하여 동일한 문자열 길이 제한, allowed-key 검사, 자동 저장, 내보내기/가져오기, 인쇄 경로를 사용한다.

기존 `evidence-cases`는 ‘이전 활동 05의 통합 해석 기록’에서 확인하고 수정할 수 있다. 값이 있는 경우만 보이며 조용히 삭제하지 않는다. 기존 두 종합 답안도 그대로 유지한다. 두 브라우저 엔진에서 old v1 JSON과 localStorage, HTML 모양의 텍스트 안전성까지 확인했다.

## 10. 인쇄

작성본에는 계산된 예상과 선택 출처, 현재 case gel, 모든 case의 관찰 자료/선택 lane/답안, 최종 종합 답안을 포함한다. 화면에서 접힌 legacy 답안도 인쇄 준비 때 잠시 펼치고 종료 후 원래 상태로 돌린다.

빈 학습지는 학생 계산, 선택 lane 표시, legacy 기록과 답안을 제외한다. 고정 사례의 관찰 자료와 질문은 유지한다. 인쇄용 mirror만 비우며 메모리/localStorage 답안을 지우지 않는다.

Gel band는 배경색이 아닌 검은 경계선으로 그려져 배경 인쇄 옵션 없이도 보인다. Figure 단위 `break-inside: avoid`로 lane label과 band를 같은 페이지에 둔다. 활동 05는 새 페이지에서 시작하며 각 사례는 새 페이지, 짧은 종합 문항 묶음은 한 페이지에 유지한다. 긴 답안은 잘라내지 않고 이어서 출력한다.

실제 Chromium PDF를 생성하고 Poppler로 관련 페이지를 PNG 렌더링했다. 작성본 PDF의 14개 활동 05 답안이 모두 포함되는지 텍스트 추출로 확인했고, 빈 PDF에는 작성한 해석/불확실성/종합 답안이 없음을 확인했다. PDF/인쇄용 PNG도 `verification.local`에만 둔다.

## 11. 접근성과 작은 화면

- lane 이름과 관찰/역할을 accessible label로 제공한다. 전체 관찰 요약은 접힌 텍스트 영역에서도 읽을 수 있고 figure의 설명으로 연결된다.
- 기본 button keyboard 동작, 유지되는 focus, visible outline, `aria-pressed`, 비색상 선택 문구를 검사했다.
- 모든 lane은 폭과 높이 44 px 이상이며 모바일에서 실제 tap을 실행했다. 작은 화면의 Positive control은 줄바꿈하며 페이지 가로 넘침이 없다.
- 대조군 설명은 선택 영역 또는 열 수 있는 안내에서 보이므로 hover에 의존하지 않는다.
- 새 animation을 추가하지 않았고 기존 reduced-motion 규칙을 유지했다.
- JavaScript가 없으면 세 필수 상황의 관찰 텍스트, 대조군 의미, 종이 기록 안내와 종합 질문을 읽을 수 있다.
- axe-core 4.11.1로 WCAG 2 A/AA, 2.1 AA, 2.2 AA 태그를 실행했다. 2개 엔진 × 3개 화면 폭 × 4개 사례 = 24회에서 violations 0, incomplete 0이었다.

## 12. 테스트 추가와 전체 결과

| 검증 | 결과 |
| --- | --- |
| `node --test tests/*.test.mjs` | 64개 통과 (기존 58 + 새 6) |
| PCR 모듈 10개 구문 검사 | 통과 |
| 기존 browser / 00~03 / 저장 / 06~07 / RNA / 인쇄 / no-JS | 536 assertions 통과 |
| Phase 2 workbench browser | 402 assertions 통과 |
| Phase 3 review browser | 708 assertions 통과 |
| Phase 4 evidence browser | 926 assertions 통과 |
| axe 활동 05 검사 | 24회, 위반 0, 미완료 0 |
| 실제 Jekyll 3.10.0 | 기존 원격 테마로 16개 페이지 빌드 |
| 포털/직접 URL/새로고침/drag/Day 1~5/fixture 실패 통합 | 통과 |
| 보호 범위, `git diff --check`, Git/Jekyll 산출물 제외 | 통과 |

브라우저는 Chromium과 WebKit, 화면은 1440×1150 / 768×1150 / 390×1150이다. 기존 browser suite 높이는 1000이다. 기존 no-JS 테스트의 전역 `noscript` 선택자는 새 05 fallback 추가로 두 요소를 찾았으므로 기존 header 요소로 범위를 명시했다. 기존 검사를 삭제하지 않았다.

추가 테스트는 사례별 band/대조군 구성, log 순서, 설계와 관찰의 독립성, old v1, 잘못된 state, 모든 답안 왕복을 검사한다. 브라우저에서는 실제 좌표 설계, 저장 snapshot 선택, 다른 draft로 변경, 사례 전환, lane 선택, 답안 입력, reload, JSON, legacy, print, invalid draft, fixture 실패, no-JS를 실행한다.

재실행 예시:

```powershell
# 기존 tests/.pcr-tools Playwright 환경. 추가 패키지: axe-core@4.11.1
$env:PLAYWRIGHT_BROWSERS_PATH = Join-Path (Get-Location) 'tests/.pcr-tools/browsers'
$env:PCR_TEST_URL = 'http://127.0.0.1:4174/bioinformatics/pcr-primer-design/'
node --test tests/*.test.mjs
node tests/pcr-browser.mjs
node tests/pcr-workbench-browser.mjs
node tests/pcr-review-browser.mjs
node tests/pcr-evidence-browser.mjs
node tests/pcr-integration.mjs
```

4174는 실제 Jekyll 출력 `tests/.pcr-output/site/`를 제공하는 로컬 서버다. 별도 미리보기는 기존 `tests/pcr-preview.mjs`의 4173을 사용할 수 있다.

## 13. 자동 PNG와 절대 경로

아래 파일은 Playwright가 실제 case 전환, lane 선택과 입력을 수행한 후 생성했다. overview도 03의 실제 좌표 선택으로 예상값을 만든 상태다. 단순 최초 페이지 캡처만 사용하지 않았다.

| 파일 | 절대 경로 |
| --- | --- |
| overview | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase4/05-overview.png` |
| 상황 1 예상 크기 band | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase4/05-case1-expected-band.png` |
| 상황 1 Sample 선택 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase4/05-case1-lane-selected.png` |
| 상황 2 NTC band | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase4/05-case2-ntc-band.png` |
| 상황 2 기록 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase4/05-case2-interpretation.png` |
| 상황 3 band 부재 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase4/05-case3-positive-control-fail.png` |
| 상황 3 Positive control 선택 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase4/05-case3-control-selected.png` |
| 상황 4 multiple bands | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase4/05-case4-multiple-bands.png` |
| 최종 종합 기록 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase4/05-final-reflection.png` |
| mobile | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase4/05-mobile.png` |
| tablet | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase4/05-tablet.png` |

`verification.json`에 URL, 생성 시각, engine/viewport별 결과, PNG 경로, 접근성 결과를 기록한다. JSON 왕복 파일과 인쇄 확인용 PDF/PNG도 같은 제외 폴더에 있다.

## 14. Contact sheet와 제외 확인

- Desktop: `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase4/phase4-contact-sheet.png`
- Tablet/mobile: `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase4/phase4-contact-sheet-mobile.png`

실제 PNG 원본을 파일명과 함께 2열로 배치한다. Desktop은 원본 약 941 px에 대해 셀 940 px이므로 텍스트를 과도하게 축소하지 않는다. 작은 화면은 셀 720 px 안에서 원본 비율을 유지한다. 원본을 자르지 않는다.

기존 `.gitignore`와 `_config.yml`이 `verification.local/`을 제외한다. `git check-ignore` 및 실제 Jekyll 출력에서 제외를 확인했다. PNG, PDF, verification JSON, 테스트 의존성, `_codex`는 커밋/push하지 않는다.

## 15. 과학 문구 확인과 수동 검증 한계

계산식과 사례 수치는 실제 구현을 읽어 기록했다. 과학적 설명은 [Addgene gel electrophoresis protocol](https://www.addgene.org/protocols/gel-electrophoresis/), [NEB NTC amplification 설명](https://www.neb.com/en-gb/faqs/2016/11/15/why-do-i-see-amplification-curves-in-my-ntc-samples), [NEB PCR troubleshooting guide](https://www.neb.com/en/tools-and-resources/troubleshooting-guides/taq-pcr-kit-troubleshooting-guide)를 대조했다. NTC 자료는 qPCR 맥락의 설명이며 이 활동에 qPCR 정량 모델을 도입한 것은 아니다.

실제 학생 실험, 물리 프린터, 실제 iPhone/iPad Safari, 화면 읽기 프로그램의 실제 발화는 검증하지 않았다. 브라우저 엔진 자동화, native keyboard/touch, 접근성 검사, DOM/매체 검사, 실제 PDF 렌더링과 텍스트 검사를 수행했다. 이 gel은 정밀한 전기영동 simulator나 실험 검증 결과가 아니다.

모든 구현/검증 후 위 수정 파일만 현재 feature branch에 커밋하고 같은 원격 branch에 push한다. 최종 commit hash와 push 결과는 작업 완료 답변에 제시한다.
