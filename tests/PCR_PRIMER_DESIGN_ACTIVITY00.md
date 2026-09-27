# 활동 00: A와 B에 같은 primer pair 적용

검증일: 2026-09-27. 작업 브랜치: `codex/activity00-shared-primer`.
기준 커밋: `deed9900f8f00d5dd026df5957f7c88dd8d3fdf1`.
Push 및 공개 사이트 배포는 수행하지 않았다. 최종 커밋 해시는 작업 응답에 기재한다.

## 변경 파일

| 파일 | 변경 내용 |
| --- | --- |
| `bioinformatics/pcr-primer-design.html` | 활동 00 질문, A/B 도식, C 판단 보류, 조작 안내, 서술형 질문, 저장 안내, 접이식 해설 |
| `assets/js/pcr-intro.mjs` | A 기준 위치에서 B 표시를 도출하는 개념 모형, 표시 갱신, 겹침/결실 미포함/방향 피드백 |
| `assets/css/pcr-worksheet.css` | A/B 동일 축척, 같은 F/R 결합 범위와 사이 구간, 경고 표시 |
| `tests/pcr-prediction.test.mjs` | 361개 위치 조합을 포함한 개념 모형 검사와 저장 상태 불변성 검사 |
| `tests/pcr-prediction-browser.mjs` | 키보드/마우스/터치, A/B 도형 길이, 경고, C, 저장 호환, 접근성, PNG 자동 생성 |
| `tests/PCR_PRIMER_DESIGN_ACTIVITY00.md` | 이 검증 보고서 |

기존 미추적 `_codex/` 및 `tests/PCR_PRIMER_DESIGN_DEPLOYMENT.md`는 수정하거나 커밋에 포함하지 않았다. PNG와 실행 결과 JSON은 기존 ignored `verification.local/`, 실행용 파일은 ignored `tests/.pcr-output/`에 보관한다.

## 개념과 화면 동작

- 같은 Forward/Reverse primer pair를 A와 B에 사용한다는 전제를 질문과 도식 범례에 명시했다. 정답은 하나의 좌표가 아니며 두 결합 부위 사이에 121~200 결실 구간 전체가 포함되어야 한다고 안내한다.
- 학생은 A 기준 위치만 선택한다. 보존된 부위의 F/R은 B에도 자동 표시한다. B 도식을 클릭해도 선택이나 저장 상태는 바뀌지 않는다.
- A는 420 bp, B는 340 bp를 같은 축척으로 그린다. 결실 왼쪽의 대응 위치는 그대로, 오른쪽은 80 bp만큼 당겨진다. 같은 primer의 개략적 결합 띠는 A/B에서 같은 화면 폭이다. 두 primer 사이의 옅은 구간도 B에서 짧아진다.
- 처음에는 F/R이 모두 숨겨져 있다. F만 고르면 A/B의 F만 나타난다. 위치를 이동하면 해당 표식과 구간, 안내가 즉시 바뀐다.
- 결실 내부 또는 경계에 결합 띠가 겹치면 A 표식에 `!`가 붙고 B의 해당 표식과 사이 구간은 숨겨진다. B 캡션과 경고가 같은 연속 결합 부위를 표시할 수 없는 이유를 설명한다.
- 둘 다 결실의 같은 쪽이면 A/B 표식은 유지하되 사이 구간의 길이가 같음을 보여 주고 결실 미포함을 안내한다. F/R 순서가 뒤바뀌거나 결합 띠가 겹치거나 맞닿으면 배치를 다시 확인하도록 안내한다.
- C는 “배경 DNA / A와 다른 서열”, “결합 여부 판단 보류”로 표시한다. 표식 부재가 결합 불가를 뜻하지 않으며 specificity는 뒤 활동에서 검토한다고 명시한다.
- 서술형 질문은 A/B 결합 가능성, 결실 위치, 산물 길이, C에 대해 지금 판단할 수 있는 범위를 설명하도록 유도한다. 예시 답안이나 유일한 정답 좌표를 제시하지 않는다.

개략적 결합 띠의 폭은 기존 슬라이더 한 단계인 5%이고 중심 양쪽으로 2.5%씩이다. 실제 primer 길이나 nt 좌표를 선택한 것으로 취급하지 않는다. 이 한계와 경계 겹침 주의를 화면에 명시했다. 모델은 염기서열 일치 검색이나 실제 증폭을 계산하지 않는다. 이 표시 모형에서는 361개 F/R 조합 중 서로 다른 45개 배치가 결실 전체를 사이에 두는 조건을 충족한다.

## 저장 호환 및 변경 범위

- `hafs:pcr-primer:v1`, schemaVersion 1, dataVersion을 유지했다.
- `initialPrimerPrediction`의 `reference: 'A'`, `units: 'relative-percent'`, `forward`, `reverse`를 그대로 사용한다. 새 B 좌표나 판정 결과를 저장하지 않고 기존 값에서 다시 그린다.
- `answers['first-placement']`와 과거 `first-negative` 기록을 보존했다. 기존 기록을 새 정답으로 바꾸거나 보정하지 않는다.
- Phase 1~6 JSON fixture 각각의 localStorage 복원과 실제 파일 입력을 통한 가져오기를 검사했다. 활동 00 위치를 바꿔도 후속 활동 상태, 설계, 답안, 최종 snapshot은 불변이었다.
- HTML의 활동 01 이후와 intro 모듈의 활동 01/02 부분은 기준 커밋과 동일하다. 기존 core/design/workbench/review/evidence/external/final/RNA 계산, 저장 validator, 교육용 fixture는 변경하지 않았다.
- 빈 학습지 인쇄 시 새로 추가한 예측 띠와 동적 피드백도 숨겨 기존 기록 노출 규칙을 유지한다. 인쇄 레이아웃을 별도로 검증하거나 재설계하지 않았다.

## 검증 결과

| 검사 | 결과 |
| --- | --- |
| 기존 전체 단위 검사와 신규 개념 검사 | 93 PASS, 실패 0 |
| 활동 00 집중 브라우저 검사 | 282 assertions PASS |
| 브라우저와 화면 폭 | Chromium / WebKit × 1440 / 768 / 390 px, 6/6 PASS |
| F/R 표시 전후, 단독 F, 마우스/터치/키보드 조작 | PASS |
| 서로 다른 유효 배치와 A/B 구간의 80 bp 축척 차이 | PASS |
| 결실 내부, R의 경계 겹침, 같은 쪽, 역순, 같은 위치 | PASS |
| B에서 독립 조작 불가, C 판단 보류 안내 | PASS |
| 현재 기록 내보내기/가져오기/새로고침 | 6/6 PASS, updatedAt 외 상태 동일 |
| 과거 Phase 1~6 localStorage 복원 및 JSON 가져오기 | 6종 모두 PASS, 후속 상태 보존 |
| 활동 00 접근성 및 가로 넘침 | 접근성 위반 0, 넘침 0 |
| 실제 Jekyll 3.10.0 빌드 | 기존 테마와 전용 layout으로 16 pages 성공 |
| 기존 전체 학생 흐름 회귀 검사 | Chromium/WebKit × 3개 폭, 6/6 PASS |
| 전체 흐름 접근성 | 18 audits, 위반 0 |
| 전체 흐름 runtime/console error, 실패 요청, HTTP 오류 | 모두 0 |
| 최종 빌드의 HTML 활동 00, CSS, intro JS와 source 대조 | 일치 |
| PNG 시각 검토 | A/B 같은 표식, B 길이 단축, 모바일, 경고, C 안내 확인 |

기존 학생 흐름의 실제 assertions를 그대로 실행했으며, 별도 대량 screenshot 단계만 제외했다. 활동 03의 정밀 배치와 설계 저장, 04의 품질 검토, 05의 virtual gel, 06의 external validation 기록, 07의 최종 snapshot, RNA 확장, 키보드, 새로고침과 JSON round trip을 포함한다. 이번 요청의 PNG는 새 전용 검사에서 생성한다. 마지막 경고 띠의 글자 가림 보정 후 실제 Jekyll 출력으로 활동 00의 6개 환경 검사를 다시 통과하고 PNG를 다시 생성했다.

재실행:

```powershell
node --test tests/*.test.mjs
$env:PLAYWRIGHT_BROWSERS_PATH = (Resolve-Path 'tests/.pcr-tools/browsers').Path
node tests/pcr-prediction-browser.mjs
```

`PCR_TEST_URL`을 지정하면 실제 빌드 서버를 대상으로 검사하며, 생략하면 기존 worksheet preview를 임시 포트에 띄워 종료 시 닫는다. 최종 캡처는 실제 Jekyll 출력 `http://127.0.0.1:4184/bioinformatics/pcr-primer-design/`에서 생성했다.

검증 데이터:

- `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/activity00/verification.json`
- `D:/Codex/260926 bioinformatics/tests/.pcr-output/activity00/regression-verification.json`

## 자동 생성 PNG 절대 경로

| 화면 | PNG 절대 경로 |
| --- | --- |
| 전체 활동 00 데스크톱 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/activity00/00-desktop-full.png` |
| 전체 활동 00 모바일 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/activity00/00-mobile-full.png` |
| A/B 같은 F/R 동시 표시 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/activity00/00-shared-ab-primers.png` |
| 결실 겹침 경고 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/activity00/00-overlap-warning.png` |
| C 판단 보류 설명 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/activity00/00-c-undecided.png` |
| F/R 표시 전 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/activity00/00-before-placement.png` |
| 결실 미포함 경고 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/activity00/00-outside-warning.png` |
| 모바일 겹침 경고 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/activity00/00-mobile-warning.png` |

## 남은 우려 / 미검증 사항

- 공개 사이트에 push/deploy하지 않았으므로 이번 수정의 공개 URL 반영은 검증하지 않았다.
- 활동 00은 결실 구조와 개략적 결합 범위에 대한 교육용 예측이다. 실제 primer 서열 특이성, 품질 또는 실험 성공을 확인한 것으로 표현하지 않는다.
- 기존 `NanumSquareNeoR.ttf` / `NanumSquareNeoB.ttf` 파일의 decode 경고가 Chromium에서 발생하며 대체 글꼴로 표시된다. 이번 변경에서 해당 파일은 수정하지 않았다. 그 외 browser console error나 runtime 오류는 없었다.
- 요청에 따라 인쇄물 레이아웃과 수업 시간 구성의 추가 검증은 하지 않았다.
