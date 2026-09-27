# PCR 학습지 1차 개선 결과

검증일: 2026-09-27. 작업 대상: `feature/pcr-primer-design`.

범위는 공통 UI와 활동 00~02이다. 활동 03~07, RNA 확장 및 자료 출처의 HTML은 작업 전과 동일하다. 계산 모듈과 교육용 fixture도 변경하지 않았다. 기존 전역 스타일, 포털과 강좌 페이지는 변경하지 않았다.

## 1. 수정한 파일

| 파일 | 변경 내용 |
| --- | --- |
| `bioinformatics/pcr-primer-design.html` | 제목, 기록/인쇄 메뉴, 목차, 활동 00~02 |
| `assets/css/pcr-worksheet.css` | 공통 UI, 답안 높이, 도입 활동 반응형 도식, 접근성, 인쇄 |
| `assets/js/pcr-intro.mjs` | 새 파일. 초기 위치 예측, 단계/주기 보기, 방향 비교, 상보/역상보 학습 |
| `assets/js/pcr-records.mjs` | 선택적인 초기 예측 및 보기 상태 추가, 이전 기록 호환, 입력 검증 |
| `assets/js/pcr-worksheet.mjs` | 도입 활동과 저장 연결, 메뉴 동작, 현재 목차 표시, 새 답안 인쇄 |
| `tests/pcr-intro-browser.mjs` | 새 파일. 00~02 상호작용과 접근성 검증 |
| `tests/pcr-records.test.mjs` | 새 파일. 기존 기록 호환과 추가 데이터 검증 |
| `tests/pcr-browser.mjs` | 새 학습 흐름, 메뉴, 기존 기록 불러오기, 인쇄 회귀 검사 |
| `tests/pcr-integration.mjs` | 기록 메뉴를 통한 내보내기 검사 반영 |
| `tests/pcr-static.test.mjs` | 새 모듈의 금지 문자와 안전한 DOM 표시 검사 |
| `tests/PCR_PRIMER_DESIGN_PHASE1.md` | 이 보고서 |

## 2. 공통 UI 변경

- 제목 계층을 `1.3.2 생물정보학`, 학습지 제목, 부제와 학습 목적 순서로 정리했다.
- 상단에는 닫힌 `기록`, `인쇄` 메뉴만 보인다. JSON 내보내기/불러오기, 초기화, 두 인쇄 기능을 유지했다. 메뉴는 키보드로 열 수 있고 Escape로 닫으면 제목에 초점을 돌려준다.
- 현재 목차는 얇은 청록색 선, 굵은 글자, `aria-current`로 표시한다. RNA 링크는 `확장 활동`으로 구분했다.
- 짧은 답 2줄, 근거 3줄, 넓은 기록 4줄의 높이를 구별했다. 기존 문항의 답안과 구조는 유지한다.
- 해설은 얇은 위아래 구분선과 `해설`, 질문 제목, `펼치기/접기`로 표현한다.
- 흰 배경, 기존 글꼴, 짙은 회색 본문과 청록색 하나를 유지했다. 점수, 진행률, 배지와 그림자는 추가하지 않았다.

## 3. 활동 00 변경

A 420 bp, B 80 bp 결실 후 340 bp, C 배경 DNA를 수평 구조도로 표시한다. A 선을 클릭/터치하거나 두 개의 기본 range 입력을 조작해 F/R 위치를 예측한다. 키보드 화살표로도 움직일 수 있다. 초기값에는 학생이 정하지 않은 프라이머를 표시하지 않는다.

이 단계에서는 정밀 좌표나 Tm을 계산하지 않는다. 화면에는 결실 구간의 왼쪽/안쪽/오른쪽만 안내하며 정답 배치를 강제하지 않는다. 설명 답안은 3줄이다.

기록의 선택적인 `initialPrimerPrediction`은 다음 의미를 갖는다.

```json
{
  "reference": "A",
  "units": "relative-percent",
  "forward": 15,
  "reverse": 80
}
```

두 값은 A 그림 전체에 대한 상대 위치로, 5부터 95까지 5 간격이며 미지정은 `null`이다. bp 좌표 또는 실제 결합 서열이 아니다. 설명은 기존 `answers["first-placement"]`에 보존한다. 향후 07에서 이 필드들을 읽어 최종 설계와 비교할 수 있도록 준비했으며 이번에는 07 UI를 변경하지 않았다.

기존 `first-negative` 답안은 삭제하지 않는다. 이전 내용이 있으면 00의 접힌 문항에서 확인/수정할 수 있고, 기존 최종 비교와 JSON에도 남는다.

## 4. 활동 01 변경

- 반응 혼합, 변성, 결합, 신장 버튼 및 이전/다음 이동을 제공한다.
- 변성에서는 두 가닥의 실제 도식 간격이 벌어진다. 결합에서는 상보적인 위치의 F/R과 경계 가이드가 나타난다. 신장에서는 primer의 3′ OH에서 출발하는 Pol 표식과 검정 새 가닥이 나타난다.
- 모든 가닥에 필요한 5′/3′ 방향을 표시했다. 신장 시 primer의 원래 3′ 위치를 새 가닥의 현재 말단으로 오인하지 않도록 Pol을 연장 시작 위치로 설명한다.
- 혼합물 표를 간결한 목록으로 바꾸고 단계의 주요 요소를 굵기와 밑줄로 강조했다.
- 1/2/3주기를 비교하며 긴 산물, 정확한 길이 단일가닥, 정확한 길이 이중가닥이 나타나는 과정을 구별한다. 분자 수나 비율을 정량화하지 않는다.
- 영역의 양 끝을 정하는 요소를 고르는 질문, 짧은 이유, 다음 주기 산물의 변화 설명을 제공한다.

95°C, 60°C, 72°C는 예시 조건이다. backbone 절단과 가닥 분리를 구별하고, annealing 조건의 의존성, Taq의 작동 온도에 대한 한계, 교육용 도식의 한계를 명시했다. 온도 설명은 [NEB PCR 최적화 안내](https://www.neb.com/en/tools-and-resources/usage-guidelines/guidelines-for-pcr-optimization-with-thermophilic-dna-polymerases)를, 초기 산물의 개념은 [NCBI Bookshelf, Studying DNA](https://www.ncbi.nlm.nih.gov/books/NBK21129/)의 긴 산물/짧은 산물 설명을 참고했다.

## 5. 활동 02 변경

primer의 3′ OH에서 시작하는 5′→3′ 합성, 내부를 향하는 두 3′ 말단, 두 주문 서열 모두 5′→3′라는 세 개념을 중심에 두었다.

정상 배치와 두 primer가 같은 방향을 향하는 잘못된 배치를 비교한다. 각 primer의 5′/3′와 합성 화살표를 표시한다. 가닥 뒤집기는 위아래 위치만 바꾸며 reference strand의 방향과 현재 위치를 설명한다.

8 nt 예제는 먼저 `TCAGGCAT` 상보 서열을 확인한 뒤 두 번째 단계를 여는 방식이다. 이어 같은 가닥을 뒤집어 주문 방향의 `TACGGACT`를 확인하고 기록한다. 첫 답을 수정하면 두 번째 단계의 확인 상태를 해제하지만 기존 작성 문자열은 보존한다. 오른쪽 reference 결합 부위라는 맥락을 설명하고, 말단 선택과 3줄 이유 문항을 제공한다.

## 6. 추가 또는 수정한 테스트

- 이전 v1 JSON에서 답안, draft, designs와 기존 3단계 위치가 복원되는지 검사한다.
- 초기 예측의 독립적인 JSON 왕복, 상대 위치 범위/간격, reference/units, 보기 상태의 열거값과 불리언, 선택 답안의 허용값을 검증한다.
- Chromium/WebKit의 1440, 768, 390px에서 위치 클릭/터치/키보드, 단계 도식 변화, Pol, 주기 비교, 가닥 뒤집기, 잘못된 배치, 상보→역상보 조작과 복원을 검사한다.
- 메뉴 키보드 동작, 보이는 focus, `aria-pressed`와 현재 목차, 작은 화면 overflow, `prefers-reduced-motion`을 검사한다.
- 03의 위치 선택/계산, 04 후보 비교, 05 gel, 06 기록, 07 설계 보존과 긴 답안을 기존 흐름으로 계속 검사한다.
- 자동 저장, JSON 왕복과 잘못된 파일 거부, 이전 작성본의 브라우저 복원/불러오기, 저장 차단, 자료 로딩 실패, 초기화 취소 및 다른 저장 키 보존을 검사한다.
- 작성본의 예측/선택 답 포함, 빈 학습지의 예측/선택 답/확인한 서열 제외, 인쇄 후 원래 기록 보존을 검사한다.

## 7. 실행한 검증과 결과

| 검사 | 결과 |
| --- | --- |
| `node --test tests/*.test.mjs` | 42개 통과. 기존 계산 테스트 포함 |
| 새 모듈과 학습지 모듈의 구문 검사 | 통과 |
| 실제 Jekyll 출력 대상 브라우저 검사 | 500개 assertion 통과. 2개 엔진 x 3개 화면 크기 |
| Jekyll 사이트 빌드 | 3.10.0, 기존 원격 테마로 16개 페이지 생성 |
| 포털 통합 | 1.3.2 번호, 직접 URL/새로고침, Day 1~5, 드래그 선택, fixture 실패 시 기록/내보내기 통과 |
| 범위 보존 | 03 이후 HTML 동일, 계산 모듈과 fixture 미변경 확인 |
| 인쇄 PDF | 테스트 긴 답안 포함 작성본 19쪽, 빈 학습지 13쪽. 00 한 장 배치, 긴 답안 마지막 문장, 빈 답안과 금지 문자 검사 통과 |
| 시각 검수 | 데스크톱/모바일 도식, 도입 활동의 인쇄 페이지를 PNG로 확인 |
| `git diff --check` | 통과 |

재실행 시 브라우저 경로와 실제 빌드 URL을 지정한다.

```powershell
$env:PLAYWRIGHT_BROWSERS_PATH = (Join-Path (Get-Location) 'tests/.pcr-tools/browsers')
$env:PCR_TEST_URL = 'http://127.0.0.1:4174/bioinformatics/pcr-primer-design/'
node tests/pcr-browser.mjs
node tests/pcr-integration.mjs
```

이 Windows 환경은 일반 Ruby gem 활성화 시 optional server 의존성인 eventmachine이 누락되어 있다. 이전에 준비된 `tests/.pcr-tools/build.rb`로 설치된 gem의 경로를 로드하여 실제 Jekyll `Site#process`를 실행했다. 소스나 테마를 모킹하지 않았다. 원격 테마 다운로드에는 네트워크 허용이 필요했다. 도구, 빌드, 스크린샷과 PDF는 제외된 `tests/.pcr-tools/`, `tests/.pcr-output/`에 두며 커밋하지 않는다.

## 8. 남아 있는 문제와 검증 한계

검사 범위에서 남은 기능 오류는 없다. 실제 iPhone/iPad, macOS Safari, 실제 프린터와 화면 읽기 프로그램을 이용한 수동 검사는 수행하지 않았다. WebKit, 터치 입력 에뮬레이션과 PDF 출력으로 확인했다. 인쇄 페이지 수는 학생의 답안 길이에 따라 달라진다.

초기 예측 좌표를 07 UI에 표시하는 연결은 다음 단계의 작업이다. 실제 PCR 성공률, 정밀 thermodynamic 값이나 유전체 특이성을 새로 계산하지 않았다. 외부 NCBI 검색도 실행하지 않았다.

이번 변경의 커밋/push 대상은 `feature/pcr-primer-design`뿐이다. main 병합 및 main push, 공개 Pages 배포는 작업 범위가 아니다.
