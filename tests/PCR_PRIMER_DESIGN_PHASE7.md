# PCR과 프라이머 디자인 Phase 7 검증 보고서

검증일: 2026-09-27 (KST). 작업 branch: `feature/pcr-primer-design`.

Phase 7은 기존 학습 기능과 과학 계산을 보존하면서 00~07의 문구, 위계, 연결과 반응형 화면을 정리했다. 수업 시간을 맞추기 위한 내용 축소가 아니며 새 활동, 계산 엔진, API, 외부 서비스 또는 의존성을 추가하지 않았다.

## 1. 기준과 작업 전 조사

- 시작 HEAD: `fcfa16f6206c62307f611af0ba66665ea7bd2266`.
- 시작 branch: `feature/pcr-primer-design`. Remote: `https://github.com/suimaire/suimaire.github.io.git`.
- 시작 `git status --short`: `?? _codex/`만 존재했다. 해당 폴더와 개인 파일을 수정하거나 staging하지 않았다.
- 저장소 및 적용 가능한 상위 경로에서 AGENTS.md가 발견되지 않았다.
- Phase 1~6 보고서를 읽고 기존 단계별 역할을 확인했다: 공통 UI/00~02, 03 workbench, 04 검토 lens, 05 gel, 06 외부 검색, 07 최종 연구 노트.
- 전체 HTML 00~07/RNA/footer, layout metadata, 전용 JS, 공통 CSS, state schema, JSON/localStorage, 포털 카드, 버튼/details/label, 인쇄 코드와 기존 테스트를 조사했다.
- `_codex/`, `verification.local/`, PNG/PDF와 개인 보조 파일은 commit 대상에서 제외했다. Main merge/push는 수행하지 않는다. 최종 commit 및 feature branch push 결과는 작업 완료 응답에 기록한다.

## 2. 수정 파일과 보호 범위

| 파일 | 변경 |
| --- | --- |
| `bioinformatics/pcr-primer-design.html` | 제목/목차/시간, 반복 문항의 이전 기록 처리, 연결, 제한 범위, 07 상세 접기, 약어와 접근성 이름 |
| `index.md` | 1.3.2 카드의 ‘약 80분’만 제거 |
| `assets/css/pcr-worksheet.css` | 화면용 목차, details, 위계, tablet stack, 모바일 지도/서열/gel, 07/노트 스타일 |
| `assets/js/pcr-worksheet.mjs` | 이전 답안 표시, 짧은 저장 상태, 기존 scroll 계산 재사용, hash 복원, 인쇄 시 이전 답안의 접힘 복원 |
| `assets/js/pcr-workbench.mjs` | 모바일 초기 예측 legend |
| `assets/js/pcr-review-view.mjs` | 첫 overview의 중복 Length/GC 제거 |
| `assets/js/pcr-evidence.mjs` | 고정 자료의 일반 문장 용어 정리 |
| `assets/js/pcr-external.mjs`, `assets/js/pcr-external-view.mjs` | 일반 문장의 데이터베이스 용어 정리 |
| `assets/js/pcr-final.mjs` | 초안/서열/데이터베이스 표시 문구 |
| `assets/js/pcr-final-view.mjs` | 편집 요약, 10항목 보고서, 같은 외부 F/R 중복 생략, heading 위계 |
| 기존 browser 테스트 7개 | 삭제한 새 답안 요구, 접힌 상세, 10항목 노트, tablet stack에 맞는 검사; 나머지 회귀 유지 |
| `tests/pcr-polish.test.mjs`, `tests/pcr-polish-browser.mjs` | 시간/anchor/legacy/전체 UX/시각 검증 |
| `tests/fixtures/pcr/phase1.json`~`phase6.json`, `README.md` | 각 역사 버전으로 직렬화한 실제 import fixture와 출처 |
| 이 보고서 | 검토 판단, 검증 및 산출물 |

Core, design, review 계산, records schema/validator, intro interaction, 교육용 fixture와 전용 layout metadata는 유지했다. Final/external의 계산·검증 함수 수정은 없으며 표시 문자열만 바꿨다. 기존 인쇄 CSS 블록은 변경하지 않았다.

## 3. 시간 표기 제거

활동 label의 `/ 5분`, `/ 10분`, `/ 20분`, 목차 시간과 RNA의 `약 15분`, 포털 카드의 `약 80분`을 제거했다. 페이지와 해당 포털 카드에 학습 분량을 시간으로 제한하는 문자열이 없다. 상단 metadata에도 시간을 넣지 않았다.

## 4. 제목 체계

| 번호 | 최종 제목 |
| --- | --- |
| 00 | 두 시료를 어떻게 구별할 것인가? |
| 01 | PCR 한 주기에서는 무엇이 달라질까? |
| 02 | 두 primer의 3′ 말단은 어디를 향할까? |
| 03 | 직접 primer를 배치하기 |
| 04 | 숫자가 적절하면 좋은 primer인가? |
| 05 | 예상과 실험 증거는 같은 것인가? |
| 06 | 실제 데이터베이스에서는 어떻게 검토할까? |
| 07 | 최종 설계 기록 |

06 설명 첫 문장에 NCBI Primer-BLAST를 유지했다. 06은 본 활동이며 RNA만 ‘확장 활동 / RNA 발현을 보려면?’으로 구분한다.

## 5. 목차와 현재 위치

왼쪽 목차를 유지하고 번호와 제목을 분리해 긴 제목의 두 번째 줄이 제목 위치에 맞춰지도록 했다. 현재 활동은 얇은 청록선, 굵은 글자와 `aria-current="location"`로 표시한다. 카드 배경이나 progress UI는 없다.

기존 수동 스크롤 계산과 requestAnimationFrame 제한을 재사용한다. Section 참조를 한 번 수집해 재검색을 줄였다. 별도 scroll listener를 더하지 않고 hashchange/pageshow에서도 같은 갱신 함수를 사용한다. 폰트와 복원 기록이 배치된 후 초기 URL hash 위치를 맞춘다. 기존 `scroll-margin-top: 20px`에서 제목이 잘리지 않아 값을 유지했다.

## 6. 학생 입력 A/B/C 검토

A는 사고와 근거의 서술, B는 짧은 답/선택/조작값, C는 별도 새 답안을 요구할 필요가 없는 반복 또는 이전 버전 보존 항목으로 분류했다. 모든 DOM 입력의 ID, label, 유형과 분류는 부록에 기록한다. 여러 radio 선택지는 한 문항으로 묶는다. 외부 route/candidate처럼 선택에 따라 나타나는 필드도 포함한다.

00 초기 이유, 03 설계 수정 이유, 04 A/B→C 판단 변화, 05 관찰/해석/불확실성, 06 목적/범위/claim/선택 근거, 07의 모든 핵심 reflection과 대조군 계획은 유지했다. 단순히 입력 개수를 일정 비율 줄이려는 기준은 사용하지 않았다.

## 7. 실제 제거·축소한 반복 질문

| 기존 key | 결정과 이유 |
| --- | --- |
| `cycle-selector` | C. primer pair 선택 직후 같은 이유를 다시 쓰는 요구를 제거. 선택 feedback과 주기 해설로 확인 |
| `direction-reason` | C. 3′ 말단 선택 뒤 원리를 재서술하는 요구 제거. 방향 원리와 reverse complement 실습 유지 |
| `design-deletion` | C. 00 초기 이유 및 실시간 결실 설명과 반복. 03의 짧은 ‘확인한 것’으로 통합 |
| `design-change` | C. 설계마다 남기는 `design-reason`과 반복. 수정 이유의 주 기록은 저장 설계에 유지 |
| `review-length-reason` | C. Length/GC 선택 후 동일한 한계를 다시 입력하는 요구 제거. lens 해설과 마지막 판단 유지 |
| `candidate-negative` | C. 03 결실 겹침과 05 음성 결과 해석의 반복. P3 인접 설명과 해설에 통합 |
| `candidate-prediction`, `review-ab-reason` | B. 2줄에서 1줄 시작 높이로 축소. textarea와 key는 유지해 기존 줄바꿈 손실 방지 |

6개 C 항목은 DOM과 allowed answer key를 삭제하지 않았다. 값이 있는 과거 기록을 불러오면 접힌 ‘이전 설명 기록’에 나타난다. 새 빈 기록에는 추가 작성란이 나타나지 않는다. 기존 값은 import/export, localStorage, 작성본 인쇄에서 그대로 보존된다. 각 field를 독립적으로 보존해 일부 답안만 있던 기록도 읽을 수 있다.

## 8. 용어와 약어

일반 문장의 sequence/database/specificity는 서열/데이터베이스/특이성으로, current draft/draft는 현재 초안/초안으로 정리했다. Primer, Forward/Reverse primer, primer pair, amplicon, off-target, Primer-BLAST, Positive control, NTC, Tm은 문맥에 맞게 유지했다. 실제 도구 label인 Database, Specificity, RefSeq RNA, accession.version, Gene ID와 공식 참고문헌 제목은 보존한다. 학생이 쓴 과거 답안의 문구는 자동 치환하지 않는다.

PCR, F/R, bp, nt, GC, Tm, dNTP, NTC, cDNA/gDNA의 의미를 도입부 또는 첫 관련 조작 근처에 짧게 설명했다. 금지된 가운데 점 검사를 전체 소스와 실제 표시 문장에 적용했다.

## 9. 버튼과 해설

Details는 제목만 쓰고 공통 UI가 ‘펼치기 + / 접기 −’를 표시한다. 해설은 ‘해설 / 질문 제목 / 펼치기’의 얇은 선 구조를 유지한다. 큰 회색 배경 카드는 추가하지 않았다.

정답 판정은 ‘상보 서열 채점’, ‘주문 서열 채점’, 단순 변환은 ‘5′→3′로 뒤집기’, 좌표 조작은 ‘좌표 적용’으로 구별했다. 외부 열기 ↗, 복사, 선택, 저장 설계 불러오기, 계산 동작은 기존 의미를 유지한다. ‘최종 연구 노트 보기/접기’는 읽기 전용 결과를 여는 동작이다.

## 10. 과학적 한계 A/B/C 재배치

| 등급 | 위치와 유지 내용 |
| --- | --- |
| A: 조작 바로 옆 | 01 온도 예시와 초기 cycle 모형, 02 방향 연습용 8 nt, 03 미지원 입력/겹침과 실제 결실 상태, 04 Wallace 근사 및 농도·열역학 미반영, 05 가상 gel·정량 아님, 06 학생 외부 기록·검색 scope |
| B: 활동 끝 | 03 primer 위치/결실 겹침 정리, 04 계산과 실험의 구별, 05 크기와 서열 특이성의 구별, 06 근거를 모으는 연결 |
| C: 전체 공통 | 교육용 A/B/C, 완전 일치 선형 모형, 실제 PCR 미수행 및 성공 보장 불가를 하단 ‘자료와 계산 범위’에 통합 |

00 시작/해설과 03 현재 설계의 반복 범위 문장을 줄였다. 03 산물 설명은 F/R, F/F, R/R 탐색과 선택 위치/다른 위치 산물의 구별에 집중한다. 05의 중복 도입과 종합 문항 앞의 동일 질문을 제거했다. 07의 최종 판단과 읽기 전용 노트에는 독립적으로 읽는 데 필요한 근거 한계를 유지했다.

## 11. 활동 간 연결과 공통 흐름

00→01→02→03→04→05→06→07에 문장과 text link를 추가했다. 특히 02 방향 규칙→실제 서열, 03 위치→primer 특성, 04 계산→실험 증거, 05 크기→특이성, 06 내부/외부 기록→최종 판단을 연결한다. 기존 04/05의 연결 문장은 교체해 이중 표시하지 않았다.

03 끝에만 ‘이번 활동에서 확인한 것’을 두 문장 분량으로 도입했다. 모든 활동에 동일한 bullet/카드나 추가 답안을 강요하지 않았다. Workbench, gel, 외부 검색 기록의 고유 구조는 보존했다.

## 12. 03 polish

중복 subtitle을 제거하고 현재 설계의 예상 산물을 읽기 쉽게 했다. 좁은 desktop과 tablet(1150px 이하)에서 workbench와 현재 설계를 세로로 배치한다. 390px에서 주요 지도 label은 확대하고 초기 예측의 작은 F/R 글자는 지도 아래 범례로 옮겼다. 좌표의 작은 글자는 12px로 조정했다. 과학 정보와 00 예측 점선은 유지한다.

선택, 드래그/키보드/터치, F/R 변환, A/B/C 산물, 결실 겹침, invalid design, 저장 설계는 이전 모형을 그대로 쓴다.

## 13. 04 polish

Overview의 F/R 서열과 산물은 남기고 길이/GC 숫자 반복은 첫 lens에 집중했다. 3′ 마지막 5 nt는 밑줄과 글자 간격으로 강조했다. 단순 상보성 정렬은 16px monospace와 얇은 왼쪽 선으로 주변 설명에서 구별한다. 긴 정렬은 기존 내부 가로 스크롤 및 인쇄 분할을 유지한다. 모든 lens와 P1/P2/P3의 단계적 비교, C 공개 후 판단은 보존한다.

## 14. 05 polish

04에서 이어지는 목적을 한 문장으로 만들고 같은 질문을 종합 문항 앞에서 다시 반복하지 않게 했다. 상황 navigation, 선택됨 문구, aria-pressed와 schematic gel은 그대로다. 관찰/해석/불확실성의 답안 구역을 얇은 구분선으로 나눴고 모바일 lane 제목은 13px로 유지한다. 개별 case의 기록은 서로 다른 관찰/해석이므로 삭제하지 않았다.

## 15. 06 polish

세 route와 검색 전 계획/실제 조건/학생 외부 결과의 차이를 보존했다. 일반 문장을 한국어로 정리하고 후보 선택 label의 반복 질문을 하나로 줄였다. 필드 열 간격과 section 여백을 조정했다. Claim scope는 청록 구분선과 본문 크기로 드러낸다. 결과의 진위나 최종 pair의 특이성을 자동 인증하지 않으며 새 NCBI 기능을 추가하지 않았다.

## 16. 07 편집 화면 압축

항상 보이는 것은 처음 예측, 좌표 중심의 최종 비교, 저장 설계 순서/수정 이유, 최종 pair 선택, F/R 주문 서열/내부 예상 산물, 연구 질문, 최종 근거, 수정 reflection, 04 핵심 수치, 05 증거 구분, 외부 상태/claim/최종 pair 일치 여부, 대조군, 미확인 사항과 최종 판단이다.

04 세부 수치와 학생 기록, 05 세부 해석, 06 조건과 후보 상세는 기본적으로 닫힌 details로 옮겼다. 기록은 자동으로 가져오며 다시 쓰도록 요구하지 않는다. 상단의 반복 flow 목록과 최종 비교 안의 상세 수치/산물 중복도 줄였다. 설계 이력은 편집기에서 수정 이유 중심의 짧은 행, 읽기 노트에서는 당시 예상/미확인 사항까지 보여준다.

동일한 합성 Phase 6 fixture를 기준 HEAD와 Phase 7에서 렌더링해 07 기본 편집 높이를 비교했다. 1440px: 약 7031→5235px, 390px: 약 8428→6368px. 이는 해당 기록/폰트의 결과이며 줄일 비율을 목표로 작업한 것이 아니다. 상세 기록은 삭제하지 않았으며 입력 내용에 따라 높이는 달라진다.

## 17. 최종 연구 노트

읽기 전용 article을 ‘최종 연구 노트 보기’로 연다. 기본 편집기에 전체 보고서를 펼쳐 중복 표시하지 않는다. 보고서는 10개 항목이다: 연구 질문, 처음 예측과 설계 변화, 최종 primer pair, 내부 예상 산물, primer 특성, 증거 해석, 외부 데이터베이스, 대조군, 미확인 사항, 현재 판단. 처음 생각에서의 수정은 설계 변화에 통합했다.

최종 F/R은 pair 항목에서 한 번 상세히 보여준다. 외부 후보가 같은 F/R이면 조건과 결과를 표시하고 같은 서열 문자열은 다시 쓰지 않는다. 다르면 외부 후보 서열을 보존하고 최종 pair의 검색 결과로 간주하지 말라는 구별을 유지한다. Report의 중첩 heading을 h3→h4→h5/h6로 맞췄다. 인증서 표현은 없다.

## 18. Desktop/tablet/mobile 및 시각 검증

전체 max-width는 바꾸지 않았다. 설명 문장에 읽기 폭을 적용하고 03/04/05/06/07별 내부 배치를 조정했다. 07 처음/최종 비교는 800px 이하에서 stack한다. 390px에서 F/R, 긴 서열, gel의 네 lane, 표/후보 비교, 노트를 실제 렌더링했다. 1024px의 좁은 desktop도 추가했다.

자동 생성된 활동별 대표 화면, mobile/tablet, 최종 노트, full-page 및 두 contact sheet를 직접 열어 활동 간격, 정보 밀도, 제목/선의 일관성과 작은 화면의 글자·겹침을 확인했다. 긴 전체 PNG는 문서 리듬 확인용이며 작은 텍스트 판독용이 아니다. 별도 모바일 현재 설계 PNG로 지도 아래 패널까지 확인한다.

## 19. 접근성

두 엔진/네 폭에서 keyboard focus, 단일 aria-current, hash/back/reload, native details, 기존 tab/sequence/lane/form 조작과 touch target을 검사했다. 표시 상태는 밑줄/선택됨/텍스트와 aria 속성을 함께 제공한다. 새 animation은 없다. 비의미적 div의 aria-label에는 group 역할을 추가했다.

Phase 7 전체 화면의 빈 기록/작성 기록/노트 상태 24회 axe 검사에서 violations 0이다. 기존 04~07 전용 axe 회귀도 유지한다. 전체 화면에서는 SVG text와 결실 기호의 `color-contrast`를 axe가 자동 확정하지 못해 incomplete로 남긴다. 해당 글자는 흰 배경의 진한 회색/청록색이며 PNG와 지정 색상을 별도로 확인했다. 사용된 네 텍스트 색의 흰 배경 대비는 약 5.05:1~12.37:1이다. 이를 자동 검사 완료로 숨기지 않는다. 실제 화면 읽기 프로그램의 발화는 수동 검증 범위에 남는다.

## 20. Legacy/state/성능

SchemaVersion 1, dataVersion 및 `hafs:pcr-primer:v1`을 유지했다. 새 저장 state는 없다. 새 details 접힘은 저장하지 않으며 기존 notebookExpanded를 재사용한다. 새로운 질문을 숨겨도 필드와 allowed key를 보존한다. JSON 가져오기 시 문자열, 개행, saved snapshot과 모든 이전 optional state를 검사했다.

Phase 1/6은 각 역사 commit의 실제 `emptyRecord`/선택 함수를 사용해 만든 재구성 fixture, Phase 2~5는 기존 브라우저 export를 역사 commit의 실제 parser로 재직렬화한 fixture다. 실제 학생 개인정보나 실제 NCBI 결과가 아니다. 각 버전의 원래 optional 구조를 보존한 파일을 두 엔진에서 localStorage와 JSON으로 각각 불러왔다. 출처는 fixture README에 명시했다.

기존 scroll listener 하나와 requestAnimationFrame 제한을 유지한다. 새 큰 dependency, framework, 외부 요청 또는 계산은 없다. 웹 학습지의 유일한 자동 데이터 fetch는 기존 교육용 fixture다.

## 21. 테스트와 인쇄 회귀

최종 테스트는 최신 실제 Jekyll 3.10.0의 16-page 빌드(4174)에서 수행한다. 기존 bundled Ruby의 optional server 문제를 피하는 기존 실제 `Site#process` 보조 경로를 사용했다. 공개 remote theme를 받는 데 sandbox 네트워크 승인이 필요했으며 승인된 빌드는 성공했다. Production 설정/테마는 변경하지 않았다.

| 검증 | 결과 |
| --- | --- |
| `node --test tests/*.test.mjs` | 83개 통과 (기존 81 + Phase 7 2) |
| 기본 browser / 00~03 / RNA / 저장 / 인쇄 / no-JS | 548 assertions |
| Phase 2 workbench | 402 assertions |
| Phase 3 review | 720 assertions |
| Phase 4 evidence | 926 assertions |
| Phase 5 external | 512 assertions |
| Phase 6 final | 670 assertions |
| Phase 7 전체 UX | 768 assertions, 전체 화면 axe 24회 |
| 브라우저 합계 | 4546 assertions |
| Jekyll / portal / direct URL / refresh / drag / Day 1~5 / fixture failure | 통과 |
| syntax / forbidden dots / safe DOM / diff whitespace | 통과 |

작성본/빈 학습지, 긴 답안, no-JS 작성란, 인쇄 후 값과 접힘 복원을 회귀했다. 새로 접은 이전 답안은 beforeprint에서 일시적으로 열고 afterprint에서 복원하여 작성본에서 사라지지 않게 했다. 빈 인쇄에는 학생 값을 넣지 않는다. 기존 print CSS redesign, page break 조정, print PNG polish와 PDF 페이지 수 최적화는 수행하지 않았다.

실행 방법:

```powershell
$env:PLAYWRIGHT_BROWSERS_PATH = Join-Path (Get-Location) 'tests/.pcr-tools/browsers'
$env:PCR_TEST_URL = 'http://127.0.0.1:4174/bioinformatics/pcr-primer-design/'
$env:PCR_VERIFICATION_ROOT = Join-Path (Get-Location) 'verification.local/pcr-primer-design/phase7/regression'
node --test tests/*.test.mjs
node tests/pcr-browser.mjs
node tests/pcr-workbench-browser.mjs
node tests/pcr-review-browser.mjs
node tests/pcr-evidence-browser.mjs
node tests/pcr-external-browser.mjs
node tests/pcr-final-browser.mjs
node tests/pcr-polish-browser.mjs
node tests/pcr-integration.mjs
```

## 22. 과학 문구 검수와 공식 자료

5′→3′ 합성, primer의 3′ OH, reverse complement, 2주기의 exact-length 단일가닥과 3주기의 이중가닥, 결실 겹침, 일반적인 길이/GC 출발점, Wallace 근사, 단순 정렬과 실제 열역학, gel 크기와 서열 정체, Positive control/NTC, 검색 scope 및 실제 PCR 미수행을 재검토했다. 올바른 기존 과학 설명과 계산은 스타일 통일을 위해 다시 만들지 않았다.

2026-09-27에 [NCBI Primer-BLAST](https://www.ncbi.nlm.nih.gov/tools/primer-blast/), [NCBI 공식 사용 안내](https://www.ncbi.nlm.nih.gov/guide/howto/design-pcr-primers/), [Primer3 Manual](https://primer3.org/manual.html), [Addgene primer 설계 안내](https://www.addgene.org/protocols/primer-design/)를 대조했다. 기존 pair/template 입력, 선택 database의 검색 범위, Custom database 예외, 농도/열역학 계산과 일반적 출발점의 의미를 확인했다. 사이트의 기존 공식 참고문헌 목록은 유지했다. 실제 검색을 제출하거나 새로운 외부 결과를 생성하지 않았다.

## 23. 남은 수동 검증 한계

실제 모바일 기기, Safari 앱 자체, 물리 프린터, 보조공학 도구의 실제 발화 및 수업 현장의 학습 효과는 이번 자동 검증에 포함되지 않는다. SVG 대비에 대한 axe incomplete는 위에 별도로 기록했다. 과학적 성공 여부, 외부 검색 기록의 진위와 실제 wet-lab PCR 결과는 이 polish가 검증하는 범위가 아니다.

## 24. PNG와 검증 자료 절대 경로

아래 PNG는 실제 UI를 자동 렌더링한 결과다. 테스트용 합성 기록이며 실제 외부 검색 결과로 표시하지 않는다. 모든 PNG 및 contact sheet는 Git에서 제외한다.

| 화면 | 절대 경로 |
| --- | --- |
| 00-polished.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase7/00-polished.png` |
| 01-polished.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase7/01-polished.png` |
| 02-polished.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase7/02-polished.png` |
| 03-mobile-current-design.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase7/03-mobile-current-design.png` |
| 03-mobile-polished.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase7/03-mobile-polished.png` |
| 03-polished.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase7/03-polished.png` |
| 03-tablet-polished.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase7/03-tablet-polished.png` |
| 04-polished.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase7/04-polished.png` |
| 04-tablet-polished.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase7/04-tablet-polished.png` |
| 05-mobile-polished.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase7/05-mobile-polished.png` |
| 05-polished.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase7/05-polished.png` |
| 06-mobile-polished.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase7/06-mobile-polished.png` |
| 06-polished.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase7/06-polished.png` |
| 07-final-notebook-mobile.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase7/07-final-notebook-mobile.png` |
| 07-final-notebook.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase7/07-final-notebook.png` |
| 07-mobile-polished.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase7/07-mobile-polished.png` |
| 07-polished.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase7/07-polished.png` |
| 07-tablet-polished.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase7/07-tablet-polished.png` |
| course-fullpage-desktop.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase7/course-fullpage-desktop.png` |
| course-mobile-top.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase7/course-mobile-top.png` |
| course-toc-desktop.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase7/course-toc-desktop.png` |
| course-top-desktop.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase7/course-top-desktop.png` |
| phase7-contact-sheet-mobile.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase7/phase7-contact-sheet-mobile.png` |
| phase7-contact-sheet.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase7/phase7-contact-sheet.png` |

Contact sheet: `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase7/phase7-contact-sheet.png` 및 `phase7-contact-sheet-mobile.png`. Full-page: `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase7/course-fullpage-desktop.png`.

구조화된 결과와 화면 높이는 같은 폴더의 `verification.json`, 입력 원본은 `input-inventory.json`, 기준 화면 높이는 `baseline-metrics.json`에 보관했다. `regression/phase2`~`phase6`에는 기존 회귀 산출물을 분리했다.

## 25. 입력 항목 전체 검토 부록

아래 표는 모든 route/candidate와 legacy를 포함한 입력 inventory다. C의 기존 답안은 삭제하지 않는다. A 중 검색 목적처럼 한 문장으로 충분한 항목은 기존 단문 입력을 유지하며, B 중 기존 개행 보존이 필요한 항목은 작은 textarea를 유지한다. Radio는 문항 단위로 묶었으며 최종 pair 선택지는 저장 설계 수에 따라 늘어난다.

대표 기록에서 DOM 입력 125개, radio를 문항으로 묶은 검토 단위 121개: A 34 / B 63 / C 24. 파일 선택은 기록 관리 도구이므로 사고 기록 분류에서 제외했다. 숨겨진 외부 경로와 선택적 추가 gel case도 포함한다. 기존 07 legacy 9개 key는 별도 입력을 새로 만들지 않고 읽기 전용 원문으로 보존한다.

| ID 또는 radio name | 문항/값 | 형태 | 분류 | 처리 |
| --- | --- | --- | --- | --- |
| `prediction-forward` | Forward 대략적 위치 | range | B | 짧은 답/선택/조작 유지 |
| `prediction-reverse` | Reverse 대략적 위치 | range | B | 짧은 답/선택/조작 유지 |
| `first-placement` | 왜 이 위치에 놓았나요? | textarea / 3줄 | A | 핵심 서술 유지 |
| `first-negative` | 한 시료에서 밴드가 보이지 않으면 해당 DNA가 없다고 결론 내릴 수 있는가? | textarea / 2줄 | C | 이전 값이 있을 때만 보존/열람 |
| `cycle-boundary-choice` | DNA polymerase / primer pair / dNTP / buffer | radio | B | 짧은 답/선택/조작 유지 |
| `cycle-selector` | 그렇게 생각한 이유를 1~2줄로 적어 보세요. | textarea / 2줄 | C | 이전 값이 있을 때만 보존/열람 |
| `cycle-template` | 새로 합성된 DNA가 다음 주기의 주형이 되면서 산물의 길이와 구성은 어떻게 달라지는가? | textarea / 3줄 | A | 핵심 서술 유지 |
| `direction-complement` | 1. 아래 가닥의 complement를 왼쪽 3′에서 오른쪽 5′로 쓰세요. | text | B | 짧은 답/선택/조작 유지 |
| `direction-reverse` | Reverse complement / 주문 서열을 기록하세요. | text | B | 짧은 답/선택/조작 유지 |
| `direction-end-choice` | 5′ / 3′ | radio | B | 짧은 답/선택/조작 유지 |
| `direction-reason` | 왜 그런가? | textarea / 3줄 | C | 이전 값이 있을 때만 보존/열람 |
| `sequence-window-slider` | 확대 위치 이동 | range | B | 짧은 답/선택/조작 유지 |
| `range-primer` | 조정할 primer | select-one | B | 짧은 답/선택/조작 유지 |
| `range-direction` | 합성 방향 | select-one | B | 짧은 답/선택/조작 유지 |
| `range-start` | 시작 좌표 | number | B | 짧은 답/선택/조작 유지 |
| `range-end` | 끝 좌표 | number | B | 짧은 답/선택/조작 유지 |
| `primer-f` | F 주문 서열 / 5′→3′ | textarea / 2줄 | B | 짧은 답/선택/조작 유지 |
| `primer-r` | R 주문 서열 / 5′→3′ | textarea / 2줄 | B | 짧은 답/선택/조작 유지 |
| `design-reason` | 이전 설계에서 무엇을 바꾸었는가? | textarea / 2줄 | A | 핵심 서술 유지 |
| `design-prediction` | 계산 전 예상 결과 | textarea / 2줄 | A | 핵심 서술 유지 |
| `design-unresolved` | 아직 확인하지 못한 점 | textarea / 2줄 | A | 핵심 서술 유지 |
| `design-deletion` | 두 시료의 길이 차이가 산물에 나타나도록 하려면 두 프라이머 사이에 어떤 구간이 포함되어야 하는가? | textarea / 3줄 | C | 이전 값이 있을 때만 보존/열람 |
| `design-change` | 결합 위치를 바꾸었을 때 예상 산물이 달라진 이유를 설명하라. | textarea / 3줄 | C | 이전 값이 있을 때만 보존/열람 |
| `review-design` | 검토할 설계 | select-one | B | 짧은 답/선택/조작 유지 |
| `review-length-choice` | Length와 GC만으로 이 primer pair가 좋은 설계라고 결론 내릴 수 있는가? | select-one | B | 짧은 답/선택/조작 유지 |
| `review-length-reason` | 그렇게 생각한 이유 | textarea / 2줄 | C | 이전 값이 있을 때만 보존/열람 |
| `review-end-observation` | 두 primer의 3′ 말단에서 눈에 띄는 특징이 있는가? | textarea / 2줄 | B | 짧은 답/선택/조작 유지 |
| `candidate-prediction` | P1, P2, P3 중 어떤 차이가 중요할 것으로 예상하는가? | textarea / 1줄 | B | 1줄 시작 textarea / 개행 보존 |
| `review-ab-reason` | 현재 정보만으로 P1과 P2를 구별할 근거가 충분한가? | textarea / 1줄 | B | 1줄 시작 textarea / 개행 보존 |
| `candidate-judgment` | A와 B만 보았을 때와 C까지 확인했을 때 판단이 어떻게 달라졌는가? | textarea / 3줄 | A | 핵심 서술 유지 |
| `candidate-negative` | P3에서 B의 산물이 예측되지 않는다는 사실만으로 B에 해당 DNA가 없다고 결론 내릴 수 있는가? | textarea / 2줄 | C | 이전 값이 있을 때만 보존/열람 |
| `review-unresolved` | 현재 primer를 실제 실험에 사용하기 전에 아직 무엇을 확인해야 하는가? | textarea / 3줄 | A | 핵심 서술 유지 |
| `evidence-case-1-observation` | 관찰 | text | B | 짧은 답/선택/조작 유지 |
| `evidence-case-1-interpretation` | 해석 | textarea / 2줄 | A | 핵심 서술 유지 |
| `evidence-case-1-uncertainty` | 아직 확정할 수 없는 점 | textarea / 2줄 | A | 핵심 서술 유지 |
| `evidence-case-2-observation` | 관찰 | text | B | 짧은 답/선택/조작 유지 |
| `evidence-case-2-interpretation` | 해석 | textarea / 2줄 | A | 핵심 서술 유지 |
| `evidence-case-2-uncertainty` | 아직 확정할 수 없는 점 | textarea / 2줄 | A | 핵심 서술 유지 |
| `evidence-case-3-observation` | 관찰 | text | B | 짧은 답/선택/조작 유지 |
| `evidence-case-3-interpretation` | 해석 | textarea / 2줄 | A | 핵심 서술 유지 |
| `evidence-case-3-uncertainty` | 아직 확정할 수 없는 점 | textarea / 2줄 | A | 핵심 서술 유지 |
| `evidence-case-4-observation` | 관찰 | text | B | 짧은 답/선택/조작 유지 |
| `evidence-case-4-interpretation` | 해석 | textarea / 2줄 | A | 핵심 서술 유지 |
| `evidence-case-4-uncertainty` | 아직 확정할 수 없는 점 | textarea / 2줄 | A | 핵심 서술 유지 |
| `evidence-identity` | 예상 크기의 band 하나가 보였다는 사실만으로 해당 band가 목표 서열이라고 확정할 수 없는 이유는 무엇인가? | textarea / 3줄 | A | 핵심 서술 유지 |
| `evidence-controls` | PCR 결과를 해석할 때 대조군이 필요한 이유를 설명하세요. | textarea / 3줄 | A | 핵심 서술 유지 |
| `evidence-cases` | 이전 상황별 해석 | textarea / 3줄 | C | 이전 값이 있을 때만 보존/열람 |
| `ext-paper-source` | 출처 / 논문 citation, DOI 또는 메모 | text | B | 짧은 답/선택/조작 유지 |
| `ext-paper-species` | Species | text | B | 짧은 답/선택/조작 유지 |
| `ext-paper-purpose` | 연구 목적 | text | A | 핵심 서술 유지 |
| `ext-paper-forward` | Forward / 5′→3′ | text | B | 짧은 답/선택/조작 유지 |
| `ext-paper-reverse` | Reverse / 5′→3′ | text | B | 짧은 답/선택/조작 유지 |
| `ext-paper-target` | 알고 있다면 target accession / 선택 입력 | text | B | 짧은 답/선택/조작 유지 |
| `ext-design-type` | 입력 종류 | select-one | B | 짧은 답/선택/조작 유지 |
| `ext-design-organism` | Organism | text | B | 짧은 답/선택/조작 유지 |
| `ext-design-target` | Target | textarea / 2줄 | B | 짧은 답/선택/조작 유지 |
| `ext-design-purpose` | PCR 목적 | select-one | A | 핵심 서술 유지 |
| `ext-design-min` | 원하는 product size / min bp | text | B | 짧은 답/선택/조작 유지 |
| `ext-design-max` | 원하는 product size / max bp | text | B | 짧은 답/선택/조작 유지 |
| `ext-plan-purpose` | 검색 목적 / 무엇을 확인하려는가? | text | A | 핵심 서술 유지 |
| `ext-plan-organism` | Target organism | text | B | 짧은 답/선택/조작 유지 |
| `ext-plan-target` | Intended target / accession, Gene ID 또는 설명 | text | B | 짧은 답/선택/조작 유지 |
| `ext-plan-database` | 검색 데이터베이스 계획 | text | B | 짧은 답/선택/조작 유지 |
| `ext-plan-notes` | 왜 이 범위를 검색하려는가? | textarea / 2줄 | A | 핵심 서술 유지 |
| `ext-conditions-date` | 검색 날짜 / YYYY-MM-DD | text | B | 짧은 답/선택/조작 유지 |
| `ext-conditions-organism` | Organism / 제한 없음이면 그 사실 기록 | text | B | 짧은 답/선택/조작 유지 |
| `ext-conditions-database` | Database / 실제 선택한 이름 | text | B | 짧은 답/선택/조작 유지 |
| `ext-conditions-target` | Target / template / 미제공이면 그 사실 기록 | text | B | 짧은 답/선택/조작 유지 |
| `ext-conditions-forward` | 검색에 사용한 Forward / 5′→3′ | text | B | 짧은 답/선택/조작 유지 |
| `ext-conditions-reverse` | 검색에 사용한 Reverse / 5′→3′ | text | B | 짧은 답/선택/조작 유지 |
| `ext-conditions-product` | Product size 조건 | text | B | 짧은 답/선택/조작 유지 |
| `ext-conditions-specificity` | Specificity 관련 주요 설정 | text | B | 짧은 답/선택/조작 유지 |
| `ext-conditions-other` | 기타 변경한 조건 | text | B | 짧은 답/선택/조작 유지 |
| `ext-candidates-0-forward` | Forward / 5′→3′ | text | B | 짧은 답/선택/조작 유지 |
| `ext-candidates-0-reverse` | Reverse / 5′→3′ | text | B | 짧은 답/선택/조작 유지 |
| `ext-candidates-0-product` | Reported product size / bp | text | B | 짧은 답/선택/조작 유지 |
| `ext-candidates-0-tmF` | Reported Tm / F °C | text | B | 짧은 답/선택/조작 유지 |
| `ext-candidates-0-tmR` | Reported Tm / R °C | text | B | 짧은 답/선택/조작 유지 |
| `ext-candidates-0-unintended` | 보고된 unintended PCR target | select-one | B | 짧은 답/선택/조작 유지 |
| `ext-candidates-0-observations` | 주요 off-target accession 또는 관찰 | textarea / 2줄 | B | 짧은 답/선택/조작 유지 |
| `ext-candidates-0-other` | Primer position 또는 기타 특징 | text | B | 짧은 답/선택/조작 유지 |
| `ext-candidates-1-forward` | Forward / 5′→3′ | text | B | 짧은 답/선택/조작 유지 |
| `ext-candidates-1-reverse` | Reverse / 5′→3′ | text | B | 짧은 답/선택/조작 유지 |
| `ext-candidates-1-product` | Reported product size / bp | text | B | 짧은 답/선택/조작 유지 |
| `ext-candidates-1-tmF` | Reported Tm / F °C | text | B | 짧은 답/선택/조작 유지 |
| `ext-candidates-1-tmR` | Reported Tm / R °C | text | B | 짧은 답/선택/조작 유지 |
| `ext-candidates-1-unintended` | 보고된 unintended PCR target | select-one | B | 짧은 답/선택/조작 유지 |
| `ext-candidates-1-observations` | 주요 off-target accession 또는 관찰 | textarea / 2줄 | B | 짧은 답/선택/조작 유지 |
| `ext-candidates-1-other` | Primer position 또는 기타 특징 | text | B | 짧은 답/선택/조작 유지 |
| `ext-selectedCandidate` | 주장과 선택 근거를 기록할 pair | select-one | B | 짧은 답/선택/조작 유지 |
| `ext-selectionReason` | 검색 조건과 결과를 근거로 이 후보를 선택한 이유는 무엇인가? | textarea / 3줄 | A | 핵심 서술 유지 |
| `ext-claimReflection` | Primer-BLAST에서 unintended target이 보고되지 않았다는 결과를 어떤 범위까지 주장할 수 있는가? | textarea / 3줄 | A | 핵심 서술 유지 |
| `ext-wetLabReflection` | 실제 PCR에서는 데이터베이스 검색 결과만으로 무엇을 아직 알 수 없는가? | textarea / 3줄 | A | 핵심 서술 유지 |
| `external-mode` | 이번 검토 경로 | select-one | C | 이전 값이 있을 때만 보존/열람 |
| `external-status` | 외부 검색 실시 여부 | select-one | C | 이전 값이 있을 때만 보존/열람 |
| `external-purpose` | 분석 목적 | text | C | 이전 값이 있을 때만 보존/열람 |
| `external-organism` | 생물종 | text | C | 이전 값이 있을 때만 보존/열람 |
| `external-template` | 주형 종류 | select-one | C | 이전 값이 있을 때만 보존/열람 |
| `external-accession` | accession.version | text | C | 이전 값이 있을 때만 보존/열람 |
| `external-database` | 검색 데이터베이스 | text | C | 이전 값이 있을 때만 보존/열람 |
| `external-date` | 검색 날짜 | text | C | 이전 값이 있을 때만 보존/열람 |
| `external-f` | 외부 후보 F / 5′→3′ | textarea / 2줄 | C | 이전 값이 있을 때만 보존/열람 |
| `external-r` | 외부 후보 R / 5′→3′ | textarea / 2줄 | C | 이전 값이 있을 때만 보존/열람 |
| `external-product` | 도구가 보고한 예상 산물 길이 / bp | text | C | 이전 값이 있을 때만 보존/열람 |
| `external-tm` | 도구가 보고한 F/R Tm / °C | text | C | 이전 값이 있을 때만 보존/열람 |
| `external-settings` | 검색 설정, 농도 조건, 검토한 전사체 범위와 프라이머 출처 | textarea / 3줄 | C | 이전 값이 있을 때만 보존/열람 |
| `external-offtargets` | 보고된 비표적 산물 | textarea / 3줄 | C | 이전 값이 있을 때만 보존/열람 |
| `external-comparison` | 이전 후보 비교와 선택 이유 | textarea / 3줄 | C | 이전 값이 있을 때만 보존/열람 |
| `external-plan` | 추가 검증 계획 / 접속하지 못한 경우 미실시 이유와 검토 계획 | textarea / 3줄 | C | 이전 값이 있을 때만 보존/열람 |
| `final-select-design-1` | 설계 1 | radio | B | 짧은 답/선택/조작 유지 |
| `final-select-draft` | 현재 초안 | radio | B | 짧은 답/선택/조작 유지 |
| `final-question` | 이 primer pair를 이용해 무엇을 구별하거나 확인하려는가? | textarea / 3줄 | A | 핵심 서술 유지 |
| `final-evidence` | 어떤 관찰과 계산을 근거로 이 primer pair를 최종 선택했는가? | textarea / 4줄 | A | 핵심 서술 유지 |
| `final-revision` | 00에서의 처음 예측과 비교해 무엇을 바꾸었으며 왜 바꾸었는가? | textarea / 4줄 | A | 핵심 서술 유지 |
| `final-control-positive` | Positive control | textarea / 2줄 | A | 핵심 서술 유지 |
| `final-control-negative` | NTC 또는 negative control | textarea / 2줄 | A | 핵심 서술 유지 |
| `final-control-additional` | 필요하다면 추가 control | textarea / 2줄 | A | 핵심 서술 유지 |
| `final-control-limit` | 각 대조군이 예상대로 나오지 않았을 때 어떤 결론을 내릴 수 없게 되는가? | textarea / 3줄 | A | 핵심 서술 유지 |
| `final-unknown` | 기타 확인할 것 | textarea / 2줄 | A | 핵심 서술 유지 |
| `final-assessment` | 현재 확보한 계산, 관찰, 외부 검색 기록의 범위에서 이 설계를 어떻게 평가할 수 있는가? | textarea / 4줄 | A | 핵심 서술 유지 |
| `rna-plan` | RNA 발현을 보기 위해 설계와 대조군 계획을 어떻게 바꾸겠는가? | textarea / 4줄 | A | 핵심 서술 유지 |
