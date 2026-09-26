# PCR 프라이머 디자인 Phase 5 검증 보고서

검증일: 2026-09-27 (KST). 작업 branch: `feature/pcr-primer-design`.

활동 06을 검색 목적과 범위를 준비하고, 학생이 외부 결과와 주장 한계를 기록하는 활동으로 재설계했다. 첫 화면에는 세 경로와 해당 준비 정보가 나오며, 계획과 외부 도구 안내는 펼쳐서 사용한다. 결과 기록은 학생의 검색 실행 표시 후 사용할 수 있다.

## 1. 기준 commit과 시작 상태

- 시작 HEAD: `5fa41c7dba4c9567c06812225e984baf7043093b` (요청한 Phase 4 commit과 일치).
- 현재 branch와 원격 `origin/feature/pcr-primer-design`도 동일한 기준 commit임을 확인했다.
- 저장소와 적용되는 상위 디렉터리에 `AGENTS.md`가 없었다.
- 최초 작업 트리에는 개인 미추적 `_codex/`만 있었다. 해당 파일과 개인 파일은 수정하거나 staging하지 않았다.
- Phase 4 보고서, 기존 06, draft/saved designs, 04 선택 규칙, 05 state, JSON/localStorage, 인쇄와 전체 테스트를 조사했다.
- `main` 병합 또는 push는 수행하지 않는다. 최종 commit hash와 기능 브랜치 push 결과는 완료 응답에 별도로 제시한다.

## 2. 수정 파일

| 파일 | 변경 |
| --- | --- |
| `bioinformatics/pcr-primer-design.html` | 활동 06의 경로, 준비, 공통 기록, 개념도, legacy |
| `assets/css/pcr-worksheet.css` | 06 선택 행, 서열, 모바일, 인쇄 전용 스타일 |
| `assets/js/pcr-external.mjs` | optional state, 검증, 학생 primer 읽기, 주장과 비교 |
| `assets/js/pcr-external-view.mjs` | 경로와 입력, copy, 수동 상태, 출력, 인쇄 |
| `assets/js/pcr-records.mjs` | 기존 v1에 optional externalSearch 추가 |
| `assets/js/pcr-worksheet.mjs` | 갱신/복원/저장/인쇄 연결, 07 상태 요약 |
| `tests/pcr-external.test.mjs` | 7개 새 단위 테스트 |
| `tests/pcr-external-browser.mjs` | 새 브라우저, 접근성, JSON, 인쇄, PNG/contact sheet 검증 |
| `tests/pcr-browser.mjs` | 기존 06 조작을 새 경로/계획 선택자로 갱신 |
| `tests/pcr-static.test.mjs` | 신규 모듈의 금지 문자와 안전한 DOM 검사 포함 |
| `tests/PCR_PRIMER_DESIGN_PHASE5.md` | 이 보고서 |

## 3. 직접 확인한 현재 NCBI 동작

2026-09-27에 공식 제출 페이지와 공식 안내를 직접 열어 확인했다. 실제 검색 제출은 하지 않았다.

| 확인 항목 | 구현 판단과 근거 |
| --- | --- |
| 두 primer 입력 | Template 없이도 기존 pair 검색을 시작할 수 있다. 공식 제출 페이지의 PCR Template 도움말에서 확인했다. |
| Template와 F/R 모두 입력 | 기존 pair의 specificity check를 수행한다. 새로운 pair를 설계하는 경우와 구분했다. 공식 How-to에 명시되어 있다. |
| Target 입력 | 현재 제출 페이지는 accession, GI, FASTA를 안내한다. UI는 accession.version과 FASTA를 제공한다. Gene ID 직접 입력 지원은 확인되지 않아 관련 nucleotide accession을 찾는 출발점으로만 안내한다. |
| RefSeq | 가능한 경우 RefSeq accession을 사용하라는 현재 공식 안내를 반영했다. |
| 검색 범위 | Organism과 database를 목적에 맞게 정한다. 제출 페이지의 Custom database에서는 organism 필드가 적용되지 않는 예외도 안내했다. |
| Pair의 가능성 | 단일 primer match와 pair product를 구별한다. F/R 외 F/F, R/R도 고려한다는 공식 설명을 확인했다. |
| RNA 확장 | RefSeq mRNA template를 요구하는 exon junction/intron inclusion 설정이 존재한다. 짧은 선택 안내만 추가했다. |

Pair만 제공한 경우에는 의도한 표적을 해석할 맥락을 학생이 살펴야 하고, 알려진 template를 제공하면 그 정체에 관한 정보를 줄 수 있다는 설명은 공식 How-to와 검색 팁에 근거한 교육적 정리다. NCBI 버튼의 위치나 모든 기본값을 복제하지 않았다.

## 4. 사용한 공식 자료와 자료 간 차이

- [NCBI Primer-BLAST 현재 제출 페이지](https://www.ncbi.nlm.nih.gov/tools/primer-blast/): target, F/R, database, organism, specificity, RNA 설정과 Tm 관련 도움말을 확인했다.
- [NCBI How-to: Design PCR primers and check them for specificity](https://www.ncbi.nlm.nih.gov/guide/howto/design-pcr-primers/): 기존 primer와 새 설계 경로, template와 두 primer의 specificity-only 동작, 범위 선택 원칙을 확인했다.
- [NCBI Primer-BLAST 설명](https://www.ncbi.nlm.nih.gov/tools/primer-blast/primerinfo.html): Primer3, BLAST와 정렬을 통한 pair 검토 및 F/R, F/F, R/R 조합을 확인했다.
- [NCBI Tips for finding specific primers](https://www.ncbi.nlm.nih.gov/tools/primer-blast/search_tips.html): RefSeq 정체 정보, organism 제한, 목적에 맞는 범위의 중요성을 대조했다. 페이지 자체의 마지막 수정 표기는 2020-11-19다.

How-to와 오래된 팁의 database 명칭 일부가 현재 제출 페이지와 다르다. 따라서 학습지에 고정 DB 메뉴나 오래된 명칭의 추천 목록을 만들지 않았다. 실제 선택한 이름과 범위를 학생이 기록한다. 테스트 PNG의 `RefSeq RNA`는 현재 제출 페이지에서 확인되는 이름을 쓴 입력 예시이며, 모든 PCR 목적에 대한 추천이 아니다.

## 5. Route A / 내 primer

04가 사용하는 `reviewDesign(state, fixture)`를 재사용한다. 명시적으로 선택한 saved design을 읽고, 선택이 없으면 현재 유효한 draft를 읽는다. 유효성은 기존 `inspectDesign`을 거친다. 유효하지 않거나 fixture가 없으면 안내만 표시한다.

F/R을 5′→3′로 표시하고 개별/함께 복사를 제공한다. P1/P2/P3를 자동으로 학생 설계에 넣지 않는다. 교육용 인공 서열임을 명시하며 자연 유전자의 표적 기록으로 자동 변환하지 않는다. 03 편집, 저장본 선택, import/reset에 따라 갱신한다. 원본 draft, saved designs와 04 선택을 변경하지 않는다.

## 6. Route B / 논문의 primer

출처 한 줄, species, 연구 목적, F/R, 선택적 target accession을 기록한다. 문헌에 보고된 사실과 현재 검색 조건의 해석을 구별하며, 기존 논문 primer에 문제가 있다는 판단을 만들지 않는다. 원문 입력과 빈 target도 저장한다.

## 7. Route C / 새 후보 설계 준비

Accession.version 또는 FASTA, organism, PCR 목적, product min/max를 준비한다. FASTA는 2줄 입력 영역에서 시작하며 여러 줄을 보존한다. Min/max의 숫자 오류나 역전은 설명하고 초안을 유지한다. Primer 생성 기능은 없으며 외부 사이트에 전달할 자료만 준비한다.

## 8. 외부 링크와 copy UX

기존 공식 URL을 활동 06 HTML의 링크에서 사용한다. JavaScript에는 외부 URL을 추가로 hard-code하지 않았다. 새 탭, `noopener noreferrer`, 새 탭임을 밝히는 accessible name을 적용했다. Query string으로 primer, organism, database를 보내지 않는다. iframe, 자동 submit, NCBI 결과 fetch 또는 결과 삽입이 없다.

내 primer F/R/pair, 문헌 primer F/R/pair/target, 새 설계 target, 검색 계획을 복사한다. 성공은 `복사됨`, 실패는 선택 가능한 텍스트로 표시한다. 단일 inline feedback은 누른 도구 옆에 표시된다.

## 9. 검색 전 계획

Purpose, organism, intended target, database 계획, 범위 선택 이유의 다섯 항목이다. 먼저 계획을 열어 작성한 뒤 같은 영역의 공식 링크를 사용한다. 실제 검색 조건은 계획에서 자동 복사하지 않는다. NCBI에서 확인한 조건을 별도로 적도록 하여 계획과 실행을 구별한다.

## 10. 상태와 결과 기록

`unperformed`, `performed`, `recorded`를 각각 미실시, 외부 검색 실행, 결과 기록 완료로 표시한다. 링크 클릭으로 상태가 바뀌지 않는다. 학생이 실행 버튼을 누르면 조건과 결과 영역을 사용할 수 있다. 검색 날짜도 자동 생성하지 않는다.

조건은 날짜, organism, database, target/template, 사용한 F/R, product size 조건, specificity 설정, 기타 변경 사항이다. 새 후보 설계에서 F/R을 제공하지 않았다면 비워 둘 수 있다. 기존 pair 검토의 기록 완료에는 사용한 F/R이 필요하다.

결과 완료 버튼은 기록에 필요한 항목과 기본 입력 형식만 확인한다. 후보의 진위나 과학적 타당성을 자동 확인하지 않는다. 조건/후보/선택을 고치면 완료 상태를 실행 상태로 돌리고 초안을 보존한다. Tm이 보고되지 않았다면 빈칸을 허용한다.

## 11. Off-target 기록

각 후보에 보고되지 않음, 보고됨, 결과를 해석하지 못함의 상태와 accession/관찰을 둔다. 초기값은 기록 전이다. 이를 합격, 인증, 완전한 특이성 같은 자동 판정으로 변환하지 않는다. 후보별로 저장하여 A와 B의 결과를 섞지 않는다.

## 12. Specificity 개념도

Intended target, 하나의 primer만 결합하는 unintended sequence, 방향/거리 조건에 따라 pair product가 가능한 unintended sequence를 HTML/CSS 선과 화살표로 표현했다. 04 A/B/C 비교와 실제 database 검색의 연결을 설명한다. 실제 PCR 발생이나 Primer-BLAST 알고리즘 전체를 재현하는 그림이 아님을 명시한다.

## 13. Claim scope

학생이 선택한 후보, 검색 날짜, organism, database, unintended target 상태를 문장에 포함한다. 결과가 보고되지 않았다는 학생 기록으로 표현하며 검색 범위 밖의 부재를 단정하지 않는다. 결과를 해석하지 못한 경우도 그 사실을 표현한다. 미실시이거나 필수 범위 정보가 없으면 주장 문장을 만들지 않는다.

바로 아래에는 검색하지 않은 organism, 다른 database/annotation/settings, wet-lab 조건으로 일반화할 수 없다는 한계를 둔다. 최종 문항은 주장 범위와 실제 실험에서 아직 모르는 점을 각각 3줄로 기록한다.

## 14. 후보 비교

최대 두 후보 A/B를 기록한다. 두 후보에 F/R, product size, unintended 상태가 기록되어야 비교 표를 표시한다. 그 전에는 `후보 비교를 기록하지 않음`이다.

표는 product size, Reported F/R Tm, 해당 두 수치의 절대 차이, unintended record, 기타 특징을 나란히 보여준다. 차이는 입력된 Reported Tm의 단순 차이며 내부 간이 Tm을 재계산하지 않는다. 초기 selectedCandidate는 빈 값이고 학생이 직접 고른다. 자동 winner, score 또는 후보 1번 추천이 없다.

## 15. State 변경과 독립성

`schemaVersion: 1`, dataVersion, localStorage 키 `hafs:pcr-primer:v1`을 유지했다. 최상위 optional `externalSearch`만 추가했다.

```text
externalSearch
  route, status, searchRoute
  paper: source, species, purpose, forward, reverse, target
  design: type, target, organism, purpose, min, max
  plan: purpose, organism, target, database, notes
  conditions: date, organism, database, target, forward, reverse,
              product, specificity, other
  sourcePrimer: origin, forward, reverse
  candidates[1..2]: forward, reverse, product, tmF, tmR,
                    unintended, observations, other
  selectedCandidate, selectionReason, claimReflection, wetLabReflection
```

준비 화면의 현재 route와 실행 확인 시 기록한 searchRoute/sourcePrimer를 구분한다. 준비 화면을 바꿔도 이미 기록한 결과를 자동으로 다른 경로의 결과라고 바꾸지 않는다. 이 학습지는 검색 기록 하나를 위한 구조이며 별도의 다중 검색 이력 관리자는 아니다.

07은 미실시를 `외부 검토 미실시`로 읽으며 완료/실행 상태와 후보 비교 기록을 요약한다. 07 자체의 문항은 그대로다.

## 16. Legacy 보존

기존 `external-*` 답안 키와 값은 answers에 그대로 유지한다. 값이 있으면 접힌 이전 활동 06 기록을 표시하고 수정/인쇄/내보내기할 수 있다. 과거의 실시 여부를 새 실행 확인으로 자동 간주하지 않는다. 기존 후보 1번 선택 문항의 편향된 label은 중립적인 이전 후보 비교 label로 바꾸고 답안 값은 보존했다.

## 17. JSON/localStorage 호환

이전 v1, Phase 1/2/3/4에 해당하는 optional field 조합을 모두 import해 기존 상태와 답안을 보존하는지 검사했다. externalSearch가 없으면 안전한 미실시 기본값을 받는다.

새 객체는 허용 키, 객체/배열 형태, 최대 후보 수 2, 문자열 길이, route/status/type/selection enum을 검사한다. 완성 상태인데 필요한 기록이 누락된 import는 거부한다. 입력 중인 잘못된 숫자나 날짜는 문자열 초안으로 저장하고 오류를 설명한다. import 실패 시 현재 기록을 유지한다. 기존 전체 JSON 1 MB 제한도 유지한다.

브라우저 두 엔진에서 export→import, reload, 미완성 입력, invalid import, legacy localStorage/JSON, HTML 형태의 legacy 텍스트를 검사했다. Imported text는 DOM의 textContent/value로 처리한다.

## 18. 인쇄

06 전용 인쇄 출력을 만들어 선택 경로, 준비 자료, 목적/범위, 실제 조건, 날짜, 사용 F/R, 외부 후보, off-target, 비교, 주장과 reflection을 포함한다. 접힌 화면 상태와 무관하게 작성 기록을 출력한다. Copy/link/button과 접기 UI는 출력하지 않는다.

미실시 출력에는 미실시를 명시하고 보관된 후보 초안을 현재 검색의 증거로 출력하지 않는다. 빈 학습지는 공통 질문과 두 후보 기록 공간만 출력하며 학생 값과 legacy 값을 넣지 않는다. 인쇄 후 메모리/localStorage는 유지된다.

Chromium에서 배경 인쇄 없이 작성본/빈 학습지/미실시 PDF를 생성했다. Poppler PNG 렌더링과 pypdf 텍스트 검사를 병행했다. 최종 비교 표는 한 페이지에 유지하고 07은 다음 페이지에서 시작한다. 작성본의 06은 19~22쪽, 빈 학습지는 16~18쪽, 미실시는 19~20쪽에서 확인했다. 번호는 이 테스트의 앞 활동 작성 상태에 따른 전체 학습지 페이지 번호다.

## 19. 접근성과 반응형

- Native button의 Tab/Enter/Space와 추가 화살표/Home/End 조작을 지원한다.
- 선택됨 텍스트와 aria-pressed를 함께 갱신하고 visible focus를 유지한다.
- Form label과 id, 오류/도움말의 aria-describedby, aria-invalid를 연결한다.
- Copy와 외부 링크에 명확한 이름을 제공한다.
- 터치 가능한 선택 행은 44 px 이상이며 모바일에서 실제 tap을 검사했다.
- Desktop 1440, tablet 768, mobile 390 px에서 가로 넘침이 없다.
- Hover에 의존하지 않으며 새 animation이 없다. Reduced motion 검사도 수행했다.
- axe-core 4.11.1의 WCAG 2 A/AA, 2.1 AA, 2.2 AA로 2 engines × 3 widths × 4 states = 24회 검사에서 violations 0, incomplete 0이다.
- Fixture 실패 때도 B/C와 기록 기능을 사용하며 no-JS에는 종이 기록 안내를 제공한다.

## 20. 전체 테스트

| 검증 | 결과 |
| --- | --- |
| 전체 `node --test tests/*.test.mjs` | 71개 통과 / 기존 64 + 새 7 |
| PCR 모듈 12개 구문 검사 | 통과 |
| 기본 browser / 00~03 / 06~07 / RNA / 저장 / 인쇄 / no-JS | 536 assertions 통과 |
| Phase 2 workbench browser | 402 assertions 통과 |
| Phase 3 review browser | 708 assertions 통과 |
| Phase 4 evidence browser | 926 assertions 통과 |
| Phase 5 external browser | 512 assertions 통과 |
| 합계 | 브라우저 3,084 assertions 통과 |
| Phase 5 axe | 24회, violations 0, incomplete 0 |
| 실제 Jekyll 3.10.0 | 기존 remote theme로 16 pages 빌드 |
| 포털 번호/직접 URL/새로고침/drag/Day 1~5/fixture 실패 | 통과 |
| 보호 범위/금지 문자/DOM/`git diff --check` | 통과 |

기존 browser suite의 06 조작만 새 경로/계획 필드로 바꿨으며 검사 수를 줄이지 않았다. Intro browser helper는 기본 browser suite에서 실행된다.

현재 Windows의 별도 Ruby 환경에서는 optional server용 eventmachine이 없어 일반 `tests/pcr-build.rb`의 gem activation이 실패한다. 이전 Phase에서 마련한 제외 파일 `tests/.pcr-tools/build.rb`의 실제 Jekyll `Site#process` 경로로 빌드했다. Jekyll 소스 수정이나 가짜 빌드를 사용하지 않았다. 원격 theme 다운로드는 네트워크 권한을 받은 실행에서 성공했다.

```powershell
$env:PLAYWRIGHT_BROWSERS_PATH = Join-Path (Get-Location) 'tests/.pcr-tools/browsers'
$env:PCR_TEST_URL = 'http://127.0.0.1:4174/bioinformatics/pcr-primer-design/'
node --test tests/*.test.mjs
node tests/pcr-browser.mjs
node tests/pcr-workbench-browser.mjs
node tests/pcr-review-browser.mjs
node tests/pcr-evidence-browser.mjs
node tests/pcr-external-browser.mjs
node tests/pcr-integration.mjs
```

4174는 실제 Jekyll 출력 `tests/.pcr-output/site/`를 제공한다. 기존 Phase 2~4 suite는 최신 소스를 즉시 읽는 4173 미리보기에서 실행했고, Phase 5와 기본 회귀/통합 검사는 최종 실제 Jekyll 출력에서도 실행했다.

## 21. 회귀 보호와 검증 결론

기준 commit과 비교해 06 바깥의 HTML 전체가 동일하다. Core, design, workbench, review, review-view, intro, evidence, evidence-view의 기존 8개 모듈도 동일하다. fixture와 기존 계산은 변경하지 않았다.

활동 06의 입력을 실제로 조작한 뒤 draft, saved designs, review, evidence, introView, 초기 예측 값이 그대로임을 브라우저에서 확인했다. 07 연결은 상태를 읽는 요약에 한정된다. 기존 RNA 확장을 재설계하지 않았다.

과학적 합격 판정, 자동 winner, 외부 결과 자동 확인 기능을 넣지 않았다. 내부 간이 Tm과 Reported Tm의 출처를 구분했다. 검색 범위와 실제 실험의 불확실성은 마지막 기록에 남긴다.

## 22. 자동 생성 PNG 절대 경로

각 PNG는 Playwright가 실제 route 선택, primer 좌표 설계, 입력, copy, 수동 상태 확인과 후보 추가를 수행한 상태에서 생성했다. 테스트의 입력은 UI 확인용 예시이며 실제 NCBI 결과를 재현한 데이터가 아니다.

| 상태 | 절대 경로 |
| --- | --- |
| 도입과 세 경로 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase5/06-overview.png` |
| 내 primer 연결 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase5/06-my-primer-route.png` |
| 문헌 primer | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase5/06-paper-primer-route.png` |
| 새 후보 준비 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase5/06-new-design-route.png` |
| 검색 전 계획 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase5/06-search-plan.png` |
| Specificity 개념 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase5/06-specificity-concept.png` |
| 외부 결과 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase5/06-result-record.png` |
| Unintended target | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase5/06-unintended-target.png` |
| 주장 범위 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase5/06-claim-scope.png` |
| 후보 비교 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase5/06-candidate-comparison.png` |
| 검색 미실시 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase5/06-no-search.png` |
| Mobile | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase5/06-mobile.png` |
| Tablet | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase5/06-tablet.png` |

## 23. Contact sheet와 산출물 제외

- Desktop: `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase5/phase5-contact-sheet.png`
- Tablet/mobile: `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase5/phase5-contact-sheet-mobile.png`

원본 PNG를 파일명과 함께 2열로 배치한다. Desktop 원본 약 941 px에 대해 cell 940 px을 사용해 과도하게 축소하지 않는다. Mobile의 긴 입력 화면도 끝까지 캡처한다. `verification.json`에는 URL, 시간, 엔진/폭별 결과, 접근성 검사와 PNG 경로가 있다.

기존 `.gitignore`와 `_config.yml`의 `verification.local/` 제외를 유지한다. PNG, contact sheet, 테스트 JSON/PDF, 캡처 보조 파일, 개인 `_codex/`는 커밋하지 않는다. 실제 Jekyll 출력에도 verification.local이 없다.

## 24. 수동 검증 한계

NCBI의 현재 공식 페이지와 안내를 확인했지만 실제 검색을 학생 대신 제출하지 않았다. 이 구현은 외부 검색 결과의 진위, 실제 PCR 성공, primer의 전역적 유일성이나 wet-lab 성능을 검증하지 않는다.

Chromium/WebKit 자동화, keyboard/touch, axe, 실제 인쇄 PDF의 텍스트 추출과 PNG 시각 검수를 수행했다. 물리 프린터, 실제 iPhone/iPad Safari와 화면 읽기 프로그램의 실제 발화, 수업 현장의 학습 효과는 검증하지 않았다. NCBI 화면과 DB는 이후 바뀔 수 있으므로 위치 중심 매뉴얼 대신 입력과 범위 기록을 남겼다.
