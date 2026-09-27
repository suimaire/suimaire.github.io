# PCR 프라이머 디자인 Phase 3 검증 보고서

검증일: 2026-09-27 (KST). 기준 커밋: `09fd27e7ee55fb37f3f90cac72dac2fa2df9205d`. 작업 브랜치: `feature/pcr-primer-design`.

활동 04를 03의 실제 학생 설계에서 시작하는 다섯 검토 관점으로 재설계했다. 길이/GC, 간이 Tm, 3′ 말단, 자기/F-R 상보성, 비표적 결합을 구분한다. 점수나 합격 판정은 만들지 않았다. 고정 P1/P2/P3는 비표적 결합 관점에서 A/B → 배경 C 순서로 공개한다.

## 1. 사전 확인과 보호 범위

- 시작 시 현재 브랜치와 HEAD가 요청과 일치했다. 기존 미추적 항목은 `_codex/`뿐이었다. 사용자 미추적 파일과 `_codex/`는 수정하거나 커밋하지 않았다.
- 저장소, 상위 경로와 하위 소스 경로에서 적용되는 `AGENTS.md`를 찾지 못했다. Phase 2 보고서, workbench state/snapshot, 04 HTML/JS/CSS, core, fixture, JSON/localStorage, 04~07 연결과 기존 테스트를 읽었다.
- 기준 커밋과 비교하여 00~03 HTML, 05~07/RNA/출처 HTML이 동일함을 확인했다. `pcr-core.mjs`, `pcr-design.mjs`, `pcr-workbench.mjs`, `pcr-intro.mjs`, fixture도 동일하다.
- 기존 `.gitignore`와 `_config.yml`이 `verification.local/`을 이미 제외한다. 새 PNG/JSON/PDF/contact sheet를 커밋하거나 배포하지 않는다. 실제 Jekyll 출력에도 해당 디렉터리가 없다.

## 2. 수정 파일

| 파일 | 변경 |
| --- | --- |
| `bioinformatics/pcr-primer-design.html` | 활동 04의 내 설계, 다섯 관점, 단계별 비교, 짧은 답안과 05 연결 |
| `assets/css/pcr-worksheet.css` | 04 전용 얇은 선, 서열/말단/정렬, 반응형 탭, 인쇄 스타일 |
| `assets/js/pcr-review.mjs` | 새 순수 함수: 설계 읽기, 말단 특징, 연속 상보 구간과 정렬 |
| `assets/js/pcr-review-view.mjs` | 새 04 UI 모듈: 관점 전환, 계산 표시, A/B/C 공개 |
| `assets/js/pcr-records.mjs` | v1의 선택적 `review` 검증과 기본값 |
| `assets/js/pcr-worksheet.mjs` | 03 변경/복원/저장과 review 갱신 연결, 기존 Tm/후보 UI 통합, 선택 답안 인쇄 |
| `tests/pcr-review.test.mjs` | 새 단위 테스트 8개 |
| `tests/pcr-review-browser.mjs` | 새 브라우저 검증, PNG/contact sheet 자동 생성 |
| `tests/pcr-browser.mjs` | 새 탭 진입에 맞춘 기존 회귀 절차와 RNA 저장/인쇄 검사 추가 |
| `tests/pcr-static.test.mjs` | 새 모듈의 금지 문자와 안전한 DOM 표시 검사 |
| `tests/PCR_PRIMER_DESIGN_PHASE3.md` | 이 보고서 |

## 3. primerStats와 Tm 조사 결과

실제로 확인한 `assets/js/pcr-core.mjs:16`의 `primerStats`는 먼저 `normalizeSequence`로 서열을 검증하고 G/C 개수를 센다. 19행의 반환식은 다음과 같다.

```js
{ sequence, length: sequence.length,
  gcPercent: 100 * gc / sequence.length,
  threePrime: sequence.slice(-5),
  simpleTm: 2 * (sequence.length - gc) + 4 * gc }
```

따라서 Tm은 `2×(A+T)+4×(G+C)`라는 염기 개수 기반 Wallace 근사다. nearest-neighbor parameter set이나 반응 조건을 받는 thermodynamic engine은 없다. Na+, Mg2+, primer 농도, dNTP, 인접 염기의 효과를 계산하지 않는다. 이 core는 변경하지 않았다.

기존 자기/F-R 상보성 기능은 실제 입력 분석이 아니라 고정 서열과 그림이었다. 기존 `reverseComplement`는 정규화된 서열의 상보 염기를 만든 뒤 순서를 뒤집는다. 새 상보성 분석은 이 함수를 재사용한다.

## 4. 03 → 04 데이터 연결

`reviewDesign(state, fixture)`가 다음 우선순위로 읽는다.

1. 04의 ‘검토할 설계’에서 명시적으로 고른 `review.designId`의 저장 snapshot.
2. 명시적인 저장 설계 선택이 없다면 `state.draft`.
3. 유효한 F/R이 없거나 좌표와 서열이 불일치하면 ‘03에서 먼저 primer를 설계하세요’ 안내.

저장 설계가 존재한다는 이유만으로 자동 선택하지 않는다. `inspectDesign`을 사용하므로 03과 같은 서열/좌표 검증과 `analyzeAll` 결과를 읽는다. 내 설계의 F/R, 길이/GC와 A/B/C 산물은 현재 선택한 데이터에서 계산한다. 유효한 서열이지만 예상 산물이 없는 설계도 관찰할 수 있다.

04의 선택/탭/공개/답안 조작은 `draft`, `bindings`, `designs[]`를 변경하지 않는다. 03 편집 중에도 04에서 명시적으로 선택한 snapshot은 그대로 유지된다. ‘현재 초안’을 선택하면 다시 초안을 따른다. 과거 후보 가져오기 버튼은 제거해 고정 사례가 04 조작으로 학생 초안을 덮어쓰지 않게 했다. 기존 저장 1~3 회귀 검증에서는 03의 원래 서열 편집 입력을 사용한다.

## 5. 활동 04 구조

제목과 흰 배경/짙은 회색/청록 한 가지 강조색을 유지했다. 내 설계를 먼저 표시하고 관점 목록, 줄 형태의 탭, 현재 관점 하나, 마지막 미확인 질문 순서로 구성한다. 큰 카드, gradient, 그림자, 게임 요소, 별점, 합격 배지와 가운데 점을 추가하지 않았다.

탭은 `tablist`/`tab`/`tabpanel`, `aria-controls`, `aria-labelledby`, `aria-selected`와 하나의 `tabindex=0`을 사용한다. 좌우 화살표, Home, End로 선택과 focus가 함께 이동한다. 좁은 화면에서는 탭이 줄바꿈된다. 계산 내용이 길면 서열 정렬 영역만 가로 스크롤할 수 있으며 영역 이름과 keyboard focus를 제공한다. JavaScript가 없을 때는 다섯 관점의 질문이 모두 정적으로 표시된다.

## 6. 길이와 GC

현재 같은 F/R의 `primerStats.length`와 `gcPercent`를 나란히 표시한다. GC는 정수 %로 표시한다. 약 18~25 nt와 40~60%는 일반적인 출발점으로만 안내한다. 자동 판정 없이 선택 질문과 2줄 이유 입력을 둔다. 인쇄에는 내부 선택값 `no` 대신 ‘결론 내릴 수 없다’ 같은 읽을 수 있는 문구가 출력된다.

## 7. Tm 구현과 한계

기존 간이 Tm을 새 관점 안으로 통합했다. F와 R의 `simpleTm`, 절대 차이 `Math.abs(F.simpleTm - R.simpleTm)`를 정수 °C로 표시한다. ‘간이 Tm / 교육용 근사’, 실제 식, 반영하지 않는 조건을 같은 화면에 명시한다. 두 값의 차이가 하나의 annealing 조건을 맞추기 어렵게 할 수 있음을 설명하며 적정 annealing 온도나 반응 성공을 자동 결정하지 않는다.

## 8. 3′ 말단

각 주문 서열은 5′→3′로 표시하며 마지막 최대 5 nt를 굵기와 청록 밑줄로 강조한다. 실제 마지막 염기, 표시한 말단의 GC 개수, 2개 이상 연속한 동일 염기 구간을 계산한다. 5 nt 미만이면 전체 서열만 표시한다. GC clamp의 유무나 반복 길이를 점수로 합산하지 않는다. 관찰 답안은 2줄이다.

## 9. 자기 상보성

각 primer와 자신의 reverse complement를 모든 가능한 겹침 위치에서 비교한다. gap과 mismatch를 건너뛰지 않고 일치가 끊기는 지점마다 연속 구간을 끝낸다. 가장 긴 구간 하나를 대표로 표시한다. 길이가 같으면 3′ 말단을 더 많이 포함하는 구간, 이동량이 작은 정렬, 좌표 순서로 결정한다.

이는 같은 서열을 가진 두 사본의 역평행 비교다. 같은 한 분자의 hairpin 구조나 고리 크기를 예측하지 않는다. 실제 형성 확률, 자유에너지, 열역학적 안정성을 계산하지 않음을 표시한다.

## 10. F/R 상보성

F와 reverse complement(R)을 같은 방법으로 비교한다. F의 마지막 좌표가 구간 끝에 있으면 F 3′ 포함, reverse complement(R)의 첫 좌표가 구간 시작이면 R 3′ 포함이다. 양쪽 3′를 동시에 포함하는지도 기록한다.

화면 아래 가닥에는 R의 **역순 서열**을 3′→5′로 그린다. reverse complement 자체를 아래 가닥에 그려 동일 염기끼리 결합하는 것처럼 보이지 않도록 했다. 밑줄과 `|`는 선택한 연속 구간에만 표시한다.

전체 최장 구간이 내부에 있어도 더 짧은 말단 상보성을 놓치지 않도록, 3′가 관여하는 구간과 양쪽 3′가 관여하는 구간을 별도로 찾는다. 대표와 다르면 접힌 상세 정렬로 제공하며 양쪽 3′ 구간을 우선 보여준다. 계산 결과에는 에너지, 확률, 구조 Tm과 품질 점수 필드가 없다. 긴 정렬은 인쇄 시 같은 열 기준으로 60열씩 나누며 조각 경계에 가짜 5′/3′ 말단을 붙이지 않는다.

## 11. P1/P2/P3 비교와 출처

고정 후보와 기준 결과는 `assets/data/pcr-primer-fixture.json`에 있다. 기존 정적 테스트가 `_codex/pcr_primer_teaching_fixture.json`과 JSON 내용의 동일성을 확인한다. fixture의 SHA-256 검증도 유지했다. 결과는 `candidatePairs`와 실제 A/B/C 서열을 `analyzeAll`에 넣어 계산한다. `expectedExactMatchProducts`는 검증 기준으로만 사용한다.

| 고정 사례 | A | B | C |
| --- | --- | --- | --- |
| P1 | 41~300, 260 bp | 41~220, 180 bp | 예상 산물 없음 |
| P2 | 71~330, 260 bp | 71~250, 180 bp | 101~320, 220 bp |
| P3 | 121~300, 180 bp | 예상 산물 없음 | 예상 산물 없음 |

처음에는 결과를 표시하지 않는다. ‘A와 B에서 비교’ 후에는 A/B만 표시하고 판단 근거를 기록한다. 이후 ‘배경 C까지 확인’으로 C를 추가한다. 산물의 짧은 bracket은 모든 주형에 같은 bp 축척을 사용한다. 후보 서열은 별도의 접힌 영역에서 볼 수 있다.

P3 질문 뒤 해설은 학생이 열어야 보인다. F의 A 결합 위치, B 결실 좌표, B의 F 완전 일치 결합 개수도 fixture에서 계산해 표시한다. ‘target 자체의 부재’와 ‘결합 부위 소실에 따른 증폭 실패’를 구별한다. 실제 데이터베이스 검색의 필요성을 설명하되 04에서 NCBI 검색은 실행하지 않는다. 마지막 이동 링크는 05로 향한다.

## 12. State 변경

`schemaVersion: 1`, `dataVersion`, localStorage 키 `hafs:pcr-primer:v1`을 유지했다. 선택적 최상위 필드만 추가한다.

```json
"review": { "lens": "length-gc", "designId": null, "stage": 0 }
```

`lens`는 다섯 고정 관점 중 하나, `designId`는 `null` 또는 존재하는 저장 설계 ID, `stage`는 0(미공개)/1(A/B)/2(A/B/C)이다. 잘못된 값과 존재하지 않는 ID는 import에서 거부한다. 정렬, 통계와 산물은 저장하지 않고 선택한 실제 서열로 다시 계산한다.

새 답안 키는 기존 `answers` 안의 `review-length-choice`, `review-length-reason`, `review-end-observation`, `review-ab-reason`, `review-unresolved`다. 기존 `candidate-prediction`, `candidate-judgment`, `candidate-negative`는 유지해 과거 답안과 07의 판단 요약을 보존한다.

## 13. JSON/localStorage 호환

`review`가 없는 이전 v1에는 기본 관점/초안/미공개 상태를 제공한다. 기존 답안과 저장 snapshot은 그대로 둔다. 새 필드는 기존 import 허용 키/값 검증에 포함되며, 새 JSON과 localStorage는 선택한 저장 설계, 관점, 공개 단계와 답안을 복원한다. 04 필드 입력은 기존 자동 저장 경로를 사용한다. 저장 차단과 fixture 로딩 실패 시에도 현재 답안의 JSON 내보내기와 인쇄가 가능하다.

## 14. 테스트 범위

- 기존 단위 테스트를 모두 유지했다. 새 8개는 초안/저장 선택, 좌표 불일치, 길이/GC/Tm/말단, 0개/완전/내부/3′ 상보성, 결정적 정렬과 v1 호환을 검사한다.
- 3 nt 서열 64개를 서로 조합한 4,096쌍에 대해 독립적인 역평행 substring 탐색과 최장/말단/양쪽 말단 결과를 대조했다. 표시한 상보 염기도 직접 대조했다.
- 기존 browser suite의 00~03, 05 gel, 06 검색 기록, 07 긴 답안, 3개 저장 설계, JSON, 저장 차단, 초기화 취소, no-JS와 인쇄 검사를 유지했다. RNA 기록의 reload/JSON/작성본/빈 인쇄 검사를 추가했다.
- Phase 2 workbench suite는 변경하지 않고 실행했다.
- 새 Phase 3 suite는 두 엔진 × 세 viewport에서 실제 03 좌표 선택, 저장 설계 우선순위, 계산값, 다섯 탭, keyboard/focus/aria, 공개 단계, 이전 v1, JSON, snapshot/초안 불변, 인쇄, 100 nt 정렬과 reduced motion을 검사한다. 엔진별 no-JS 및 fixture 실패도 확인한다.
- PNG는 초기 화면 캡처만 하지 않고 테스트가 해당 조작을 수행한 후 생성한다. PNG 원본을 이용해 두 열 contact sheet를 자동 생성한다.

## 15. 실행 결과

| 검증 | 결과 |
| --- | --- |
| `node --test tests/*.test.mjs` | 58개 통과 |
| PCR JS 모듈 8개 구문 검사 | 통과 |
| 기존 browser + RNA 회귀 | 536개 assertion 통과 |
| Phase 2 workbench 회귀 | 402개 assertion 통과 |
| Phase 3 browser | 708개 assertion 통과, PNG 10개와 contact sheet 2개 생성 |
| 실제 Jekyll 3.10.0 빌드 | 기존 원격 테마로 16개 페이지 생성 |
| 포털/직접 주소/새로고침/드래그/fixture 실패 통합 | 통과, Day 1~5 링크 포함 |
| 코드 보호 범위 | 00~03, 05~07/RNA HTML 및 core/fixture 동일 |
| Git diff / 산출물 제외 | `git diff --check` 통과, Git 및 Jekyll 제외 확인 |

viewport는 1440×1150, 768×1150, 390×1150이며 기존 browser suite의 높이는 1000이다. Chromium과 WebKit을 사용했다. 최초 Jekyll 빌드는 sandbox 네트워크 차단으로 기존 테마를 받지 못했고, 필요한 네트워크 권한으로 다시 실행해 실제 빌드를 완료했다. 코드/테마를 모킹하지 않았다. 정렬 출력 테스트에서 발견한 JavaScript의 `-0`은 0으로 정규화했다.

재실행 예시:

```powershell
$env:PLAYWRIGHT_BROWSERS_PATH = (Join-Path (Get-Location) 'tests/.pcr-tools/browsers')
$env:PCR_TEST_URL = 'http://127.0.0.1:4174/bioinformatics/pcr-primer-design/'
node --test tests/*.test.mjs
node tests/pcr-browser.mjs
node tests/pcr-workbench-browser.mjs
node tests/pcr-review-browser.mjs
node tests/pcr-integration.mjs
```

4174는 실제 Jekyll 출력 `tests/.pcr-output/site/`의 로컬 서버다. 학습지 전용 미리보기는 `node tests/pcr-preview.mjs`와 4173을 사용한다. 이 환경은 Phase 1/2의 로컬 Ruby 도구로 실제 Jekyll `Site#process`를 실행했다.

## 16. 자동 생성 PNG 절대 경로

공통 경로: `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase3/`.

| 상태 | 절대 경로 |
| --- | --- |
| 03 초안으로 시작한 내 설계 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase3/04-my-design.png` |
| Length/GC와 짧은 답안 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase3/04-length-gc.png` |
| 간이 Tm와 한계 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase3/04-tm.png` |
| F/R 3′ 말단 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase3/04-end-structure.png` |
| self 및 F/R 정렬 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase3/04-complementarity.png` |
| A/B만 공개 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase3/04-ab-only.png` |
| 배경 C까지 공개 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase3/04-with-background-c.png` |
| P3 답안과 결합 부위 해설 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase3/04-p3-binding-failure.png` |
| mobile / 선택된 저장 설계와 상보성 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase3/04-mobile.png` |
| tablet / 선택된 저장 설계와 상보성 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase3/04-tablet.png` |

`verification.json`에는 실행 URL, 생성 시각, 엔진/폭별 assertion 수와 PNG 경로를 기록한다. 엔진/폭별 JSON 왕복 파일, `phase3-filled.pdf`, `phase3-blank.pdf`도 같은 제외 폴더에 있다.

## 17. Contact sheet

- Desktop 8개 상태: `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase3/phase3-contact-sheet.png`
- Tablet/mobile: `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase3/phase3-contact-sheet-mobile.png`

실제 PNG를 두 열로 배치하고 파일명을 위에 표시했다. Desktop 셀은 940 px로 원본 941 px와 거의 같은 크기이며, 작은 화면 원본은 확대하지 않는다. 비율을 유지하며 일부를 잘라내지 않는다. 주요 개별 PNG를 직접 열어 내 설계, 말단, 정렬, 후보 비교와 작은 화면의 배치/잘림을 검수했다.

## 18. 수동 검증 한계와 과학 문구 확인

실제 iPhone/iPad/macOS Safari, 실제 화면 읽기 프로그램의 발화, 물리 프린터는 수동 검증하지 않았다. 브라우저 엔진 자동화, touch/keyboard, DOM/인쇄 매체와 PDF 생성으로 확인했다. 열역학, 실제 PCR 성공, hairpin 구조, 실제 유전체 특이성은 계산 범위가 아니다.

과학 문구는 기존 코드/fixture와 함께 제조사 및 NCBI의 1차 안내를 확인했다. 3′ 상호 상보성 검토는 [NEB primer design 안내](https://www.neb.com/en-gb/nebinspired-blog/proven-tips-for-pcr-primer-design), 검색 범위에 따른 specificity 검토는 [NCBI Primer-BLAST 안내](https://www.ncbi.nlm.nih.gov/guide/howto/design-pcr-primers/)를 참조했다. Tm 계산법에 대한 보고는 외부 문서에서 추측한 것이 아니라 위 core의 실제 식에 근거한다.

모든 검증이 통과한 뒤 이번 수정 파일만 `feature/pcr-primer-design`에 커밋하고 같은 원격 브랜치로 push한다. `main` 병합/push는 수행하지 않는다.
