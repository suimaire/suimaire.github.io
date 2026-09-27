# PCR 프라이머 디자인 Phase 6 검증 보고서

검증일: 2026-09-27 (KST). 작업 branch: `feature/pcr-primer-design`.

활동 07을 이전 기록을 자동으로 모으고 학생이 최종 실험 계획을 서술하는 연구 노트로 재설계했다. 최종 pair는 학생이 직접 선택하며 선택 당시의 서열과 결합 좌표를 별도로 보존한다.

## 1. 기준 commit과 작업 전 조사

- 시작 HEAD: `25073abebb0009811425e063347cf06ccdde35f7`. 요청한 Phase 5 기준과 정확히 일치했다.
- 시작 branch: `feature/pcr-primer-design`, 추적 branch: `origin/feature/pcr-primer-design`.
- 저장소와 적용 가능한 상위 경로에서 `AGENTS.md`가 발견되지 않았다.
- 최초 작업 트리에는 개인 미추적 `_codex/`만 있었다. 개인 파일과 기존 verification.local 산출물은 수정하거나 staging하지 않았다.
- Phase 5 보고서, 전체 state와 import/export, localStorage, 00 예측과 설명, 03 draft와 bindings 및 저장 설계, 04 선택과 답안, 05 종합 reflection, 06 검색 기록과 claim scope, 기존 07, 인쇄, 기존 테스트를 조사했다.
- `main`에는 merge 또는 push하지 않았다. 완료 commit과 기능 branch push 결과는 최종 응답에 별도로 보고한다.

## 2. 수정 파일

| 파일 | 변경 |
| --- | --- |
| `bioinformatics/pcr-primer-design.html` | 07 자동 요약, 선택, 핵심 서술, 구조화된 대조군, 읽기 전용 노트 |
| `assets/css/pcr-worksheet.css` | 07 비교, 이력, 서열, 모바일, 인쇄 스타일 |
| `assets/js/pcr-final.mjs` | optional state, legacy migration, 선택 snapshot, 이전 데이터 읽기 adapter |
| `assets/js/pcr-final-view.mjs` | 안전한 DOM 렌더링, 최종 선택과 입력, report 및 인쇄 |
| `assets/js/pcr-records.mjs` | optional finalReview 기본값과 import 검증 |
| `assets/js/pcr-worksheet.mjs` | 저장/복원/변경/인쇄 시 07 연결과 legacy key 허용 |
| `tests/pcr-final.test.mjs` | Phase 6 단위 검증 |
| `tests/pcr-final-browser.mjs` | Phase 6 브라우저, 접근성, JSON, 인쇄, PNG/contact sheet |
| `tests/pcr-browser.mjs` | 07 인쇄 검사를 새 연구 노트 출력에 연결 |
| `tests/pcr-external-browser.mjs` | 07 외부 요약의 새 위치 검사 |
| `tests/pcr-workbench-browser.mjs`, `tests/pcr-review-browser.mjs`, `tests/pcr-evidence-browser.mjs`, `tests/pcr-external-browser.mjs` | 기존 산출물을 보존하기 위한 선택적 `PCR_VERIFICATION_ROOT` |
| `tests/pcr-static.test.mjs` | 신규 모듈 금지 문자와 안전한 DOM 검사 |
| `tests/PCR_PRIMER_DESIGN_PHASE6.md` | 이 보고서 |

## 3. 기존 07 구조

기존 `final-comparison`은 00의 서술, 저장 설계, 일부 판단과 06 상태를 보여주었다. 그러나 F/R과 예상 산물도 학생이 다시 입력해야 했고, 최종 설계 source 선택이나 draft 변경에서 독립된 최종 pair가 없었다. `final-*` 답안 9개는 answers에 저장되었다. 인쇄는 textarea 옆에 전체 텍스트를 복제하는 기존 공통 방식이었다.

## 4. 00 초기 예측 연결

`initialPrimerPrediction`의 reference A, units relative-percent와 forward/reverse를 그대로 읽는다. `first-placement`는 최초 설명으로 읽는다. 위치는 00과 같은 구간 기준으로 왼쪽/안쪽/오른쪽 및 '대략적 위치'라고 표시한다. 상대값을 정밀 bp 좌표로 바꾸지 않는다. 축약 구조도도 원래 백분율 위치를 사용하며 121~200은 알려진 결실 구간의 label일 뿐 예측 좌표가 아니다.

00 값은 기존 역사 기록 필드를 읽는 live reference다. 07에서는 편집하거나 대체하지 않는다. 학생이 00 자체를 수정하면 07도 현재 00 기록을 읽는다. 값이 없으면 '기록 없음'이며 과거 생각을 새로 작성하는 입력란은 없다.

## 5. 03 설계 history

`designs`의 순서를 유지하는 의미 있는 ordered list를 사용한다. 실제 저장한 1~3개만 출력한다. 각 항목에 F/R 좌표, 기존 계산으로 구한 A/B/C 산물, 저장 당시의 reason과 unresolved를 표시한다. 설계가 없으면 '저장된 설계 없음'이다. 계산할 수 없는 예전 설계도 원문을 삭제하지 않는다. 03의 snapshot에는 쓰기를 하지 않는다.

## 6. 최종 설계 선택 정책

기존 `inspectDesign`에서 오류 없이 결과와 inward-facing 선택 산물을 정의할 수 있는 저장 설계와 현재 draft만 native radio로 제공한다. 유효하지 않은 draft, 겹치는 위치, 잘못된 방향이나 불완전 pair는 후보에서 제외된다. 아무 후보도 자동 선택하지 않는다. 'best', winner, 점수나 검증 완료 판정은 없다.

`selectedDesignSource`에 `design-1`, `design-2`, `design-3`, `draft` 중 학생이 고른 source를 저장한다. source가 달라지는 경우에도 해당 선택은 학생의 명시적인 조작에서만 발생한다.

## 7. current draft와 최종 pair snapshot

모든 최종 선택 시 `primerSnapshot`에 forward, reverse, F/R bindings만 깊은 복사로 저장한다. `selectedAt`도 기록한다. 계산값, 기존 설계 이유, 04~06 답안은 복제하지 않는다.

**03 draft를 이후 수정해도 07의 최종 서열과 좌표는 바뀌지 않는다.** 현재 원본과 보존한 pair가 달라졌거나 원본이 더 이상 유효하지 않으면 그 사실을 문장으로 알린다. 유효한 변경 draft를 반영하려면 '변경된 현재 draft로 최종 선택 갱신'을 눌러야 한다. 무효 draft는 갱신 후보로 제공하지 않는다. saved design도 source reference와 작은 pair snapshot을 함께 보존하여 import/추후 원본 상태에 관계없이 최종 pair가 조용히 바뀌지 않도록 했다.

선택한 pair의 length, GC, 주문 서열과 A/B/C 산물은 이 snapshot을 기존 계산 함수로 읽어 표시한다. fixture가 없으면 저장한 서열은 보존하고 계산값은 표시하지 않는다.

## 8. 04 설계 검토 요약

최종 pair에 `inspectDesign`/primerStats, `endFeatures`, `complementarity`를 그대로 적용한다. Length/GC, 간이 Tm과 차이, 3′ 말단, 자기 상보성, F/R 상보성을 표시한다. GC 반올림과 간이 Tm도 04와 일치한다. 새로운 열역학 공식을 만들지 않았다.

그 아래에 **04 학생 기록**을 별도 출처로 둔다. 현재 04 검토 대상(저장 설계 또는 draft)을 표시하며, 과거 답안을 쓸 당시의 pair는 기존 schema에 별도로 묶여 있지 않음을 설명한다. 길이/GC 이유, 3′ 관찰, 제한된 A/B/C 판단, 마지막 `review-unresolved`를 그대로 읽는다. 답안이 없으면 '기록 없음'이다. 최종 pair와 04 답안이 같은 설계를 가리킨다고 단정하지 않는다.

## 9. 05 실험 증거 해석

`evidence-identity`, `evidence-controls` 종합 답안만 중심으로 가져온다. 모든 case와 gel을 다시 출력하지 않는다. '05의 gel 자료는 수업용 가상 자료이며 이 primer pair의 실제 실험 결과가 아닙니다.'를 화면과 report에 명시한다. 현재 05에 별도의 최종 추가 검증 아이디어 field가 없어 없는 답을 만들지 않는다.

## 10. 06 외부 database 검토

현재 externalSearch의 상태, 실제 searchRoute, 날짜, organism, database, target, specificity 설정, 학생이 선택한 candidate, unintended 상태와 관찰/accession, 외부 후보 F/R 및 product, 선택 이유와 두 reflection을 읽는다. 기존 `claimScope`를 재사용한다. 미실시에는 미실시를 명확하게 표시하고 이전 후보 초안을 현재 증거로 출력하지 않는다. 실행만 한 상태도 '결과 기록 중'으로 구분한다.

표제는 '학생이 기록한 외부 검색 결과'다. 학습지가 NCBI를 직접 fetch하거나 결과를 검증했다고 표현하지 않는다. 선택한 외부 candidate와 최종 pair의 F/R이 일치하는지도 구분한다. 다르거나 비교할 자료가 없으면 그 검색을 최종 pair의 결과로 간주하지 않도록 안내한다. 일치할 때도 외부 기록의 진위를 검증한 것은 아니며 검색 범위 제한을 유지한다.

## 11. 대조군 계획

Positive, NTC 또는 negative, optional additional control을 각각 2줄로 둔다. 각 대조군이 예상대로 나오지 않았을 때 내릴 수 없는 결론은 3줄이다. 다른 활동의 답을 자동으로 대조군 계획에 채우지 않는다.

## 12. 아직 확인하지 못한 것

실제 PCR/gel 미수행, amplicon identity 미확인, 실제 조건의 primer-dimer/hairpin 미검증은 이 학습지 범위의 자동 상태다. 06 상태는 현재 학생 기록에 따라 표시한다. 04의 상보성 계산과 05의 가상 gel을 실제 수행 완료로 바꾸지 않는다. 외부 기록 완료도 최종 pair의 specificity 판정이 아니라고 표시한다. 학생은 기타 확인할 것과 현재 근거 범위에서의 판단을 별도로 서술한다.

## 13. 읽기 전용 최종 연구 노트

페이지 안의 button으로 펼치고 접는 article이며 12개 요청 항목을 순서대로 담는다. 최종 선택 근거는 최종 pair 항목에 포함한다. 이전 활동 내용과 학생 서술을 안전한 DOM textContent로 구성한다. report에는 편집 입력란이나 모달이 없다. 흰 배경, 청록 heading, 가는 구분선, monospace 서열을 사용한다. Desktop 비교는 좌우, mobile 비교는 세로 배치다.

## 14. Legacy 07 migration

`finalReview`가 없는 기록을 읽을 때에만 다음을 복사한다.

| 이전 answers key | 새 의미 |
| --- | --- |
| final-question | researchQuestion |
| final-evidence | finalRationale |
| final-revision | revisionReflection |
| final-unknown | otherLimitations |

기존 `answers`는 하나도 삭제하거나 덮어쓰지 않는다. F/R, 예상 산물, 통합 대조군과 이전 한계 문항은 의미를 임의로 새 field에 나누지 않고 접힌 '이전 07 기록'에 보존한다. 명확히 mapping한 항목의 원문도 함께 보관한다. 이전 F/R을 최종 pair로 자동 선택하지 않는다. 이미 finalReview가 존재하면 빈 답안도 학생의 현재 기록으로 존중하여 legacy에서 재주입하지 않는다.

## 15. State 변경

`schemaVersion: 1`, `dataVersion`, localStorage key를 유지하고 optional `finalReview`만 추가했다.

```text
finalReview
  selectedDesignSource, selectedAt
  primerSnapshot: null | { forward, reverse, bindings: { F, R } }
  researchQuestion, finalRationale, revisionReflection
  controls: { positive, negative, additional, interpretationLimit }
  otherLimitations, finalAssessment, notebookExpanded
```

최종 pair는 선택 시 snapshot, 03 이력은 기존 immutable saved snapshots, 00은 기존 역사 기록 필드, 04~06 서술과 상태는 현재 기록을 읽는 live reference다. 화면과 report도 같은 adapter를 이용한다.

## 16. JSON/localStorage 호환

`hafs:pcr-primer:v1`을 유지했다. Old v1 및 Phase 1~5 optional field 조합을 복원하고 기존 값을 보존한다. 새로운 객체의 허용 key, 필수 구조, 문자열 길이, 선택 enum, boolean, 시각, snapshot 서열과 bindings 구조를 검사한다. 최대 JSON 1 MB 및 기존 answers 검증도 유지한다. 오류 import는 현재 기록을 바꾸지 않는다.

최종 pair, 모든 서술과 대조군, expansion 상태의 JSON round trip과 reload를 검증했다. 가져오기 직전 작성란에 초점이 남아 있어도 새 기록 값이 확실하게 복원되도록 연결했다. 기존 localStorage 접근 차단/손상 대응은 유지한다.

## 17. 인쇄

작성본은 접힘 여부와 무관하게 07의 읽기 전용 연구 노트를 출력하고 편집 화면은 중복 출력하지 않는다. 학생 답, source가 표시된 최종 pair, 06 상태, controls, limitations를 포함한다. 긴 서술은 고정 높이로 자르지 않고 여러 페이지로 자연스럽게 이어진다. 제목과 다음 내용의 분리, 개별 서열과 요약 행의 분리를 피하며, 설계 이력은 한 설계 단위로 유지한다.

빈 학습지는 별도로 조립한 report에서 과거 기록을 placeholder로 바꾸고 서술란을 비운다. 현재 학생 값과 legacy를 넣지 않는다. 인쇄 후 메모리/화면/localStorage는 유지된다. no-JS에서도 07의 종이 작성란을 인쇄할 수 있다.

Chromium의 배경 인쇄 없는 A4 PDF를 생성하고 pypdf 텍스트 검사와 Poppler PNG 시각 검수를 병행했다. 최종 페이지 수와 확인 범위는 아래 검증 결과에 기록한다.

## 18. 접근성과 반응형

- Native radio의 Space, 방향키, Tab 선택과 visible focus를 검사했다. 선택할 때 radio DOM을 불필요하게 다시 만들지 않아 초점이 유지된다.
- 44 px 이상 선택 행, 터치로 펼치기/접기, label 연결, aria-expanded/controls, 서열 group의 accessible name을 제공한다.
- 초기 그림에는 바로 아래 텍스트 대안이 있고 history는 ordered list다. 상태는 색만으로 구별하지 않는다.
- Desktop 1440, tablet 768, mobile 390 px에서 page horizontal overflow가 없다. 서열은 내부에서 안전하게 줄바꿈한다.
- 새 animation이 없으며 reduced motion 상태도 검사한다.
- axe WCAG 2 A/AA, 2.1 AA, 2.2 AA를 두 엔진과 세 폭에서 빈 기록, 완성/펼침, legacy/부분 기록에 적용했다.

## 19. 테스트 실행

```powershell
$env:PLAYWRIGHT_BROWSERS_PATH = Join-Path (Get-Location) 'tests/.pcr-tools/browsers'
$env:PCR_VERIFICATION_ROOT = Join-Path (Get-Location) 'verification.local/pcr-primer-design/phase6/regression'
node --test tests/*.test.mjs
node tests/pcr-browser.mjs
node tests/pcr-workbench-browser.mjs
node tests/pcr-review-browser.mjs
node tests/pcr-evidence-browser.mjs
node tests/pcr-external-browser.mjs
$env:PCR_TEST_URL = 'http://127.0.0.1:4174/bioinformatics/pcr-primer-design/'
node tests/pcr-final-browser.mjs
node tests/pcr-integration.mjs
```

기존 보조 Ruby 환경은 optional server의 eventmachine이 없어 Phase 5와 같은 실제 Jekyll `Site#process` 경로로 빌드했다. 이번 작업용 제외 파일에서는 cache 출력 경로만 tests 아래로 지정했으며 기존 개인 보조 파일은 편집하지 않았다. 기존 remote theme 다운로드가 허용된 실행에서 Jekyll 3.10.0이 16 pages를 생성했다. 4174는 이 실제 출력을 제공한다. 기존 회귀는 최신 소스를 제공하는 4173에서, Phase 6 최종 검증과 통합 검사는 4174에서 수행했다.

## 20. 전체 검증 결과

| 검증 | 결과 |
| --- | --- |
| 전체 `node --test tests/*.test.mjs` | 81개 통과 / 기존 71 + 신규 10 |
| PCR 모듈 구문 | 14개 통과 |
| 기본 browser / 00~03 / RNA / 저장 / 인쇄 / no-JS | 536 assertions 통과 |
| Phase 2 workbench browser | 402 assertions 통과 |
| Phase 3 review browser | 708 assertions 통과 |
| Phase 4 evidence browser | 926 assertions 통과 |
| Phase 5 external browser | 512 assertions 통과 |
| Phase 6 final browser | 676 assertions 통과 |
| 브라우저 합계 | 3,760 assertions 통과 |
| Phase 6 axe | 2 engines × 3 widths × 3 states = 18회 / violations 0, incomplete 0 |
| 실제 Jekyll 3.10.0 | 기존 remote theme로 16 pages 빌드 |
| 포털 1.3.2 / 직접 URL / 새로고침 / drag / Day 1~5 / fixture 실패 | 통합 검사 통과 |
| no-JS 07 읽기와 인쇄 작성란 | 두 엔진 통과 |
| JSON / reload / legacy / 긴 답안 / invalid import | 통과 |
| 금지 문자 / 안전한 DOM / `git diff --check` | 통과 |

Phase 5 기준과 비교해 07 밖의 HTML 전체(00~06 및 RNA)는 동일하다. Core, design, workbench, review, review-view, intro, evidence, evidence-view, external, external-view의 기존 10개 모듈도 동일하다. Fixture 및 기존 계산은 수정하지 않았다. 07 조작 전후에 draft, saved designs, 04/05 state, 00 예측 및 06 외부 기록의 불변성을 검사했다.

작성본에는 수천 자의 선택 근거를 넣어 여러 페이지로 이어지는 마지막 문장까지 확인했다. 최종 F/R은 잘리지 않고 출력되고, 빈 인쇄에는 연구 질문, 대조군, 외부 DB 이름 등 테스트 학생 값이 없다. PDF의 07 모든 관련 페이지를 Poppler 이미지로 확인했다. 이 테스트 상태에서 전체 작성본은 32쪽이며 07은 23~32쪽, 전체 빈 학습지는 22쪽이며 07은 19~21쪽이다. 작성본 마지막 페이지에는 기존 참고문헌도 이어진다. 페이지 수는 앞 활동의 답안과 긴 답안 분량에 따라 달라진다.

## 21. 생성 PNG 절대 경로

모든 요청 PNG는 Playwright가 실제 UI로 예측 위치, 세 설계, 선택, 04/05 답안, 06 검색 상태와 후보, 07 서술을 구성한 후 생성한다. 외부 결과 예시는 실제 NCBI 결과를 주장하지 않는 UI 테스트 데이터다.

| 화면 | 절대 경로 |
| --- | --- |
| 07 도입과 흐름 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase6/07-overview.png` |
| 설계 1~3 이력 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase6/07-design-history.png` |
| 처음과 최종 비교 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase6/07-initial-vs-final.png` |
| 최종 F/R과 산물 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase6/07-final-primer.png` |
| 04 요약 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase6/07-primer-review.png` |
| 05 해석 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase6/07-experimental-interpretation.png` |
| 06 완료 기록 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase6/07-external-review.png` |
| 06 미실시 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase6/07-no-external-search.png` |
| 대조군 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase6/07-controls.png` |
| 미확인 사항 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase6/07-limitations.png` |
| 최종 연구 노트 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase6/07-final-notebook.png` |
| Mobile | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase6/07-mobile.png` |
| Tablet | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase6/07-tablet.png` |

## 22. Contact sheet와 제외 정책

- Desktop: `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase6/phase6-contact-sheet.png`
- Tablet/mobile: `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase6/phase6-contact-sheet-mobile.png`
- 검증 상세: `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase6/verification.json`

2열 contact sheet에서 desktop cell은 940 px로 설정하여 원본 약 941 px의 주요 화면을 읽을 수 있게 유지했다. 긴 연구 노트도 끝까지 포함한다. 기존 Phase 2~5 캡처를 덮어쓰지 않고 이번 회귀 결과는 phase6/regression 아래에 생성했다.

기존 `.gitignore`와 `_config.yml`의 verification.local 제외를 유지한다. 실제 Jekyll 출력에서 verification.local 및 _codex가 없는 것도 검사했다. PNG/PDF/contact sheet/검증 데이터와 개인 파일은 commit하지 않는다.

## 23. 남은 수동 검증 한계

자동 검증은 Chromium/WebKit, 세 viewport, keyboard/touch, axe, 실제 생성 PDF의 텍스트와 PNG 렌더링을 대상으로 했다. 물리 프린터, 실제 iPhone/iPad와 화면 읽기 프로그램의 발화, 수업 현장에서의 학습 효과는 별도 확인이 필요하다.

학습지는 실제 PCR, 실제 gel, amplicon identity, 열역학적 dimer/hairpin 및 외부 결과의 진위를 검증하지 않는다. 06의 실제 검색을 학생 대신 실행하지 않았다. 이 Phase는 기존 과학 계산 및 검색 워크플로를 재설계하지 않고 그 기록을 정확한 출처와 한계로 모으는 데 한정했다.
