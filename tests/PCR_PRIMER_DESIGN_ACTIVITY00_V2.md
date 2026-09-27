# 활동 00 문구 및 한국어 줄바꿈 보정

검증일: 2026-09-27. 브랜치: `codex/activity00-shared-primer`.
기준 커밋: `3b5ad2c9a4eafc28eb25bdab4307f33af046289a`.
현재 브랜치에만 로컬 커밋한다. Main 변경, push, 공개 배포는 수행하지 않았다.

## 수정 범위

| 파일 | 변경 내용 |
| --- | --- |
| `bioinformatics/pcr-primer-design.html` | 활동 00의 legend, A/B 캡션, 대응 결합 부위와 좌표 설명, 80 bp 차이의 원인, C 안내 |
| `assets/js/pcr-intro.mjs` | B의 상태별 캡션과 접근성 설명 문구 |
| `assets/css/pcr-worksheet.css` | 학습지 전반의 한국어 본문/제목 줄바꿈 및 서열/코드 예외 |
| `tests/pcr-prediction-browser.mjs` | 자연스럽게 수정된 C 안내에 맞춘 기존 assertion |
| `tests/pcr-typography-browser.mjs` | 실제 단어 줄바꿈, computed CSS, 긴 값/후속 활동 overflow, fallback, PNG 생성 검사 |
| `tests/PCR_PRIMER_DESIGN_ACTIVITY00_V2.md` | 이 보고서 |

새 조작 기능, 계산 또는 저장 필드는 추가하지 않았다. 활동 01 이후 HTML은 기준 커밋과 동일하다. Primer 대응 위치 및 판정 계산, 후속 활동 계산, localStorage와 JSON schema는 그대로다. 기존 미추적 `_codex/` 및 `tests/PCR_PRIMER_DESIGN_DEPLOYMENT.md`는 수정하거나 커밋에 포함하지 않았다.

## 학생용 문구

- 상단 범례: “A와 B에는 같은 primer pair를 사용합니다. 옅은 영역은 두 primer 사이의 구간입니다.”
- A/B 모두에 결합 부위를 표시할 수 있을 때: “A와 B에 동일한 F/R primer를 사용합니다. B에서도 결실 바깥의 보존된 결합 부위에 같은 primer가 결합합니다.”
- B의 결실 오른쪽 결합 부위는 bp 좌표가 80만큼 앞당겨짐을 설명한다. 비교하는 것은 같은 숫자 좌표가 아니라 같은 primer가 인식하는 결합 부위임을 명시한다.
- 도식 아래에, 두 primer가 결실 양옆의 보존된 부위에 결합하면 A에는 121~200의 80 bp가 포함되고 B에는 없으므로 같은 pair의 B 산물이 80 bp 짧아진다고 설명한다. 잘못된 배치에서도 무조건 80 bp 차이가 난다고 읽히지 않도록 조건을 먼저 제시했다.
- C는 다른 서열의 배경 DNA이며 아직 결합 여부를 판단하지 않는다고 안내한다. 표식 부재가 결합 불가를 뜻하지 않고, 특이성은 뒤 활동에서 검토한다는 의미를 유지했다.
- 활동 00에 260 bp/180 bp와 같은 절대 PCR 산물 길이는 추가하지 않았다. “같은 축척”은 하단 보조 설명으로 옮겼다.

## Typography 규칙

모든 규칙은 `#pcr-worksheet` 내부에 한정했다. 다른 포털/강좌의 CSS에는 영향을 주지 않는다. 글자 크기나 본문 폭을 바꾸지 않았다.

일반 문단, 질문, 해설, label, list item, caption, transition, RNA 본문, 작은 안내와 기록 텍스트에는 `word-break: keep-all`과 `line-break: strict`를 적용한다. 공백 없는 토큰이 컨테이너보다 길 때만 안전하게 나뉘도록 `overflow-wrap: anywhere`를 함께 둔다. 지원하는 브라우저에서만 `@supports (text-wrap: pretty)`로 본문의 마지막 줄을 개선한다. 제목은 `keep-all`을 유지하고 `@supports (text-wrap: balance)`로 줄 길이를 균형 있게 조정한다.

입력 필드와 `.pcr-sequence`, `.final-sequence`, code, monospace 서열/좌표 등에는 `word-break: normal`, `line-break: auto`, `overflow-wrap: anywhere`, `text-wrap: wrap`을 명시한다. 입력 필드는 기존 폭 안에 머무르게 한다. 공백 정렬이 중요한 `pre`, `.pcr-structure`, `pre code`는 `white-space: pre`, `text-wrap: nowrap`, `overflow-wrap: normal`로 별도 처리해 열을 유지한다. 기존 `.pcr-alignment-scroll`의 가로 스크롤을 보존했다.

레이아웃용 수동 `<br>`, `&nbsp;`, 특정 문장 전용 non-breaking span을 추가하지 않았다. 활동 00의 기존 `<br>` 두 개는 자연스러운 문장으로 정리하면서 제거했다.

브라우저 동작을 검토할 때 참고한 공식 문서:

- [MDN word-break](https://developer.mozilla.org/en-US/docs/Web/CSS/Reference/Properties/word-break): keep-all의 CJK 단어 경계 처리.
- [MDN text-wrap](https://developer.mozilla.org/en-US/docs/Web/CSS/Reference/Properties/text-wrap): pretty/balance의 목적과 구현별 차이.
- [MDN overflow-wrap](https://developer.mozilla.org/en-US/docs/Web/CSS/Reference/Properties/overflow-wrap): 긴 토큰의 안전한 줄바꿈.
- [MDN line-break](https://developer.mozilla.org/en-US/docs/Web/CSS/Reference/Properties/line-break): CJK 문장의 줄바꿈 규칙.

## 검증 결과

| 검사 | 결과 |
| --- | --- |
| 기존 단위 검사 | 93 PASS |
| 활동 00 표시/조작/경고/저장 호환 검사 | 282 assertions PASS |
| Typography 및 긴 값 검사 | 294 assertions PASS |
| 화면 폭 | Desktop 1440px / tablet 768px / mobile 390px |
| 브라우저 | Chromium 153.0.8010.12 / WebKit 26.6, 각 3개 화면 폭 |
| 첫 문단의 실제 한국어 단어 경계 | 각 환경의 15개 한국어 포함 단어에서 내부 줄바꿈 0 |
| pretty/balance computed CSS 및 지원 확인 | 두 시험 엔진에서 모두 지원, 의도한 값 적용 |
| progressive enhancement fallback | 두 @supports 블록을 제거한 상태에서도 keep-all 유지, 한글 단어 분리 0, 페이지 가로 넘침 0 |
| 긴 기술 값 | 100 nt primer, 긴 URL, accession.version, 200자리 Gene ID, 긴 영문 식별자 보존 및 가로 넘침 0 |
| 후속 표시 | 03 workbench / 04 alignment / 05 gel / 06 외부 값 / 07 notebook / RNA transcript 검사 PASS |
| 실제 Jekyll 3.10.0 빌드 | 기존 테마와 layout으로 16 pages 성공 |
| 기존 전체 학생 흐름 | 두 브라우저 × 3개 폭, 6/6 PASS |
| 전체 흐름 접근성 | 18 audits, 위반 0 |
| 전체 흐름 console/runtime 오류 및 실패 요청/HTTP 오류 | 모두 0 |
| 소스 및 실제 빌드 대조 | 활동 00 HTML, CSS, intro JS 일치 |
| PNG | 아래 19개 자동 생성 후 주요 화면 시각 검수 |

기존 전체 흐름 assertions는 변경 없이 실행했다. 출력 폴더만 `tests/.pcr-output/activity00-v2`로 바꾸고, 별도 대량 PNG 단계 대신 이번 전용 캡처를 사용했다. 활동 03~07과 RNA의 입력/조작, JSON 내보내기·가져오기, 새로고침, 최종 snapshot 보존을 포함한다. 최종 typography 검사와 아래 PNG는 실제 Jekyll 빌드 서버에서 생성했다.

재실행:

```powershell
node --test tests/*.test.mjs
$env:PLAYWRIGHT_BROWSERS_PATH = (Resolve-Path 'tests/.pcr-tools/browsers').Path
$env:PCR_PREDICTION_OUTPUT = 'verification.local/pcr-primer-design/activity00-v2/prediction-regression'
node tests/pcr-prediction-browser.mjs
node tests/pcr-typography-browser.mjs
```

`PCR_TEST_URL`을 지정하면 해당 빌드 서버를 검사하며, 생략하면 임시 worksheet preview를 실행하고 종료 시 닫는다.

검증 원본:

- `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/activity00-v2/typography-verification.json`
- `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/activity00-v2/prediction-regression/verification.json`
- `D:/Codex/260926 bioinformatics/tests/.pcr-output/activity00-v2/regression-verification.json`

## PNG 절대 경로

| 화면 | 경로 |
| --- | --- |
| 활동 00 전체 데스크톱 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/activity00-v2/00-v2-desktop.png` |
| A/B 개념 도식과 설명 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/activity00-v2/00-v2-ab-concept.png` |
| 활동 00 전체 모바일 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/activity00-v2/00-v2-mobile.png` |
| 데스크톱 상단 줄바꿈 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/activity00-v2/typography-desktop.png` |
| 모바일 상단 줄바꿈 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/activity00-v2/typography-mobile.png` |
| 활동 00 전체 태블릿 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/activity00-v2/00-v2-tablet.png` |
| 태블릿 상단 줄바꿈 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/activity00-v2/typography-tablet.png` |
| 03 데스크톱 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/activity00-v2/regression-03-desktop.png` |
| 03 모바일 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/activity00-v2/regression-03-mobile.png` |
| 04 데스크톱 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/activity00-v2/regression-04-desktop.png` |
| 04 모바일 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/activity00-v2/regression-04-mobile.png` |
| 05 데스크톱 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/activity00-v2/regression-05-desktop.png` |
| 05 모바일 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/activity00-v2/regression-05-mobile.png` |
| 06 데스크톱 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/activity00-v2/regression-06-desktop.png` |
| 06 모바일 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/activity00-v2/regression-06-mobile.png` |
| 07 데스크톱 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/activity00-v2/regression-07-desktop.png` |
| 07 모바일 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/activity00-v2/regression-07-mobile.png` |
| RNA 데스크톱 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/activity00-v2/regression-rna-desktop.png` |
| RNA 모바일 | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/activity00-v2/regression-rna-mobile.png` |

별도 `prediction-regression` 하위 폴더에는 기존 활동 00 검사가 생성하는 표시 전, 겹침/결실 미포함 경고, C 설명 등의 PNG 8개도 보존했다.

## 브라우저 호환성과 남은 사항

- pretty/balance는 보조 개선이다. 미지원 시 기본 줄바꿈으로 동작하며 keep-all과 긴 값 overflow 처리는 유지된다. 시험에서는 두 엔진 모두 지원했고, 미지원 경로는 enhancement 규칙을 제거해 확인했다. 실제 구형 브라우저 및 Firefox 실행은 하지 않았다.
- 글꼴과 브라우저의 조판 방식에 따라 정확한 줄 위치가 달라질 수 있다. 매우 긴 토큰은 넘침 방지를 위해 예외적으로 단어 안에서 나뉜다. 모든 문단의 마지막 줄을 특정 형태로 강제하지 않는다.
- 기존 NanumSquareNeo R/B 글꼴 파일의 decode 경고가 계속 있어 시스템 대체 글꼴로 검수했다. 글꼴 파일은 이번 범위에서 수정하지 않았다.
- Push/deploy하지 않았으므로 공개 사이트 반영은 아직 없다. 로컬 PNG를 먼저 검수할 수 있다.
