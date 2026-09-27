# PCR 학습지 구현 및 검증 기록

검증 날짜: 2026-09-26

## 범위

기존 Jekyll + Just the Docs 저장소 안에 `/bioinformatics/pcr-primer-design/`를 추가했다. 포털의 자동 번호 생성 방식을 유지하여 생물정보학 기초 다음에 1.3.2로 표시한다. 작업 브랜치는 `feature/pcr-primer-design`이다. 기존 강좌, Day 1~5, 전역 스타일과 전역 메뉴 구현은 변경하지 않았다.

`_codex/`의 제공 파일 두 개는 수정하지 않았다. 배포용 fixture는 원본 JSON과 동일하다. 새 clone의 테스트는 배포용 fixture를 사용하고, 정규화한 JSON의 SHA-256으로 제공 데이터와의 동일성을 확인한다. 원본 `_codex/`가 있으면 두 JSON도 직접 비교한다.

## 생성하거나 수정한 파일

| 파일 | 역할 |
| --- | --- |
| `index.md` | 기존 목록 형식으로 새 학습지 항목 추가 |
| `_layouts/pcr-worksheet.html` | 한국어 문서, 전용 자산, 전역 메뉴와 분리한 Jekyll 레이아웃 |
| `bioinformatics/pcr-primer-design.html` | 활동 00~07, RNA 확장, 자료 출처, 문항과 기록란 |
| `assets/css/pcr-worksheet.css` | 학습지 안으로 격리한 반응형 및 인쇄 스타일 |
| `assets/data/pcr-primer-fixture.json` | 제공된 고정 인공 서열과 기준 결과 |
| `assets/js/pcr-core.mjs` | 정규화, 역상보, GC, 좌표, 모든 방향의 결합과 산물 탐색 |
| `assets/js/pcr-records.mjs` | 기록 버전과 JSON 구조 및 크기 검증 |
| `assets/js/pcr-worksheet.mjs` | 조작, 그림, 설계 보존, 자동 저장, JSON, 인쇄 |
| `tests/pcr-core.test.mjs` | 계산 및 기록 회귀 테스트 |
| `tests/pcr-static.test.mjs` | fixture, 금지 문자, 안전한 표시, 포털 연결 검사 |
| `tests/pcr-browser.mjs` | Chromium/WebKit의 주요 사용 흐름과 인쇄 검사 |
| `tests/pcr-integration.mjs` | 실제 Jekyll 출력의 번호, URL, 기존 강좌 링크 검사 |
| `tests/pcr-preview.mjs` | Ruby 없이 가능한 학습지 전용 로컬 미리보기 |
| `tests/pcr-build.rb` | 기존 사이트 전체를 생성하는 Jekyll 빌드 검사 |
| `tests/.gitignore` | 임시 도구, 브라우저, 빌드와 검증 산출물 제외 |
| `tests/PCR_PRIMER_DESIGN.md` | 이 구현 및 검증 기록 |

## 구현한 기능

1. 80분 활동 00~07과 접힌 15분 RNA 발현 확장. 학생의 첫 예측, 수정된 설명, 외부 결과와 가상 자료를 별도로 기록한다.
2. PCR의 변성, 결합, 신장을 직접 이동한다. 5′/3′ 표기와 새 가닥의 3′ 말단 첨가를 표시하고 초기 주기의 긴 산물을 설명한다.
3. 양 가닥 확대 보기에서 마우스 드래그, 두 지점 클릭, 터치, 화살표와 Enter로 구간을 선택한다. 숫자 좌표 및 주문 서열 직접 입력도 지원한다. 좁은 화면은 서열 열 수를 줄인다.
4. F/R 이름에 방향을 강제하지 않고 양 방향의 모든 완전 일치 위치와 F/R, F/F, R/R 산물을 탐색한다. 비중첩 조건, 양쪽 결합 구간을 포함한 산물 길이, 반복 결합의 대안 정보를 보존한다.
5. 설계 1~3은 저장 시점의 서열과 설명을 보존하며 이후 편집으로 덮어쓰지 않는다. 최종 기록에서 초기 설명과 나란히 검토할 수 있다.
6. P1~P3의 A/B 비교와 C 추가 비교, 기본 조성 분석, 접힌 간이 Tm, 명시적으로 구분된 고정 이차구조 예시를 제공한다.
7. 흑백 전기영동 개념도에 A, B, C, A+C, B+C를 표시한다. 같은 길이 밴드가 겹쳐도 출처가 있는 산물 목록을 유지한다.
8. 공식 Primer-BLAST 링크, 두 검토 경로, 출처를 포함하는 복사, 검색 조건과 후보 비교 기록을 제공한다. 외부 검색은 자동 실행하지 않는다.
9. 전용 localStorage 키, 저장 실패 안내, JSON 내보내기와 검증된 불러오기, 확인을 받는 현재 학습지 초기화를 지원한다.
10. 작성본과 빈 학습지를 구별해 인쇄한다. 긴 답안 전체를 별도 인쇄 요소로 출력하고, 해설은 제외하며 원래 기록은 유지한다.

## 실행한 검증

### 계산 및 정적 검사

```powershell
node --test tests/*.test.mjs
node --check assets/js/pcr-worksheet.mjs
git diff --check
```

기존 15개와 새 계산/기록 19개, 새 정적 검사 4개를 합해 **38개 통과**. 제공된 A/B/C 길이, B 결실, P1~P3의 9개 기준 산물 결과, 잘못된 R, 반복 결합, 동일 프라이머의 양 방향, F/F와 R/R, 경계 오류와 겹침, 미지원 문자, 혼합 시료의 출처 및 JSON 왕복을 확인했다.

### 브라우저

테스트 도구는 사이트 런타임에 추가하지 않는다. 재실행 준비:

```powershell
npm install --prefix tests/.pcr-tools --no-audit --no-fund playwright@1.63.0
$env:PLAYWRIGHT_BROWSERS_PATH = (Join-Path (Get-Location) 'tests/.pcr-tools/browsers')
node tests/.pcr-tools/node_modules/playwright/cli.js install chromium webkit
node tests/pcr-browser.mjs
```

학습지 전용 미리보기와 실제 Jekyll 출력에 각각 실행했다. Chromium과 WebKit에서 1440, 768, 390px의 **236개 검증 통과**. 실제 Jekyll 출력 대상으로 실행할 때는 다음 주소를 설정한다.

```powershell
$env:PCR_TEST_URL = 'http://127.0.0.1:4174/bioinformatics/pcr-primer-design/'
node tests/pcr-browser.mjs
node tests/pcr-integration.mjs
```

포털에서 1.3.2가 정확히 한 번 표시됨, 기존 1.3.1 유지, 링크 이동, 직접 주소와 새로고침 HTTP 200, Day 1~5의 기존 링크 HTTP 200, 마우스 드래그, fixture 로딩 실패 시 기록과 내보내기가 유지됨을 확인했다.

그 밖에 터치 및 키보드 선택, 좌표 입력, 서열 방향, 오래된 계산 결과 제거, 저장 설계 보존, 새로고침 복원, JSON 왕복, 악성 HTML 문자열의 비실행, 저장 차단, JavaScript 비활성화 안내, 빈 학습지 인쇄 시 기록 유지, 초기화 취소와 다른 저장 키 보존을 확인했다.

### Jekyll 빌드

Windows에 Ruby가 없어 테스트 폴더에 임시 Ruby 3.3.12를 준비했다. 서버용 EventMachine 네이티브 빌드는 개발 도구 부재로 설치되지 않았으나 정적 사이트 생성에는 필요하지 않아, 실제 Jekyll `Site#process` API를 사용했다. Jekyll 또는 테마 소스를 수정하거나 모킹하지 않았다.

Jekyll 3.10.0, jekyll-remote-theme 0.4.3, jekyll-seo-tag 2.9.0, jekyll-include-cache 0.2.2와 기존 `pmarsceill/just-the-docs` 원격 테마로 **16개 페이지 빌드 성공**. 출력은 `tests/.pcr-output/site/`다. 새 파일의 출력 경로는 `bioinformatics/pcr-primer-design/index.html`이다.

일반 Ruby/Jekyll 환경에서는 필요한 gem이 설치된 상태로 다음 스크립트를 실행한다.

```powershell
ruby tests/pcr-build.rb
python -m http.server 4174 --bind 127.0.0.1 --directory tests/.pcr-output/site
```

### 인쇄 및 시각 검사

스크린샷과 PDF는 `tests/.pcr-output/`에 남겼으며 Git에 포함하지 않는다. 데스크톱과 모바일 서열 정렬을 실제 이미지로 확인했다. Chromium 인쇄 PDF는 긴 테스트 답안 포함 작성본 17쪽, 빈 학습지 12쪽이었다. 마지막 답안 문장 보존, 금지 문자 부재, 빈 학습지의 학생 답안 제외를 텍스트로 검사하고 대표 인쇄 페이지를 Poppler로 렌더링해 확인했다. 페이지 수는 학생의 답안 길이에 따라 달라진다.

## 가장 간단한 로컬 미리보기

현재 저장소 폴더에서 다음 명령을 실행한다. Node.js 외의 설치는 필요 없다.

```powershell
node tests/pcr-preview.mjs
```

주소: http://127.0.0.1:4173/bioinformatics/pcr-primer-design/

이 서버는 학습지 전용이다. 전체 포털 미리보기와 배포 확인을 대신하지 않는다. 전체 Jekyll 출력은 위 4174 서버로 확인한다.

## 아직 검증하지 않은 범위

1. 실제 iPhone/iPad 및 macOS Safari 장치, 실제 프린터, 화면 읽기 프로그램을 이용한 수동 검토는 하지 않았다. WebKit과 인쇄 PDF로 확인했다.
2. NCBI에 프라이머를 제출하거나 검색을 실행하지 않았다. 공식 링크와 사용 설명 자료만 확인했다.
3. 정밀 Tm, 열역학적 안정성, 실제 PCR 성공률과 유전체 특이성은 구현 및 검증 범위가 아니다. 간이 Tm을 합격 판정이나 결합 온도 추천에 사용하지 않는다.
4. 공개 Pages 배포는 수행하지 않았다. 원격 기본 브랜치는 main이며 저장소에 별도 Actions workflow는 없다. 인증 없는 Pages 설정 API가 404를 반환해 실제 배포 대상 브랜치 설정은 확정하지 못했다. main 병합 및 main push는 사용자 요청에 따라 수행하지 않는다.

공개 반영 시의 경로는 `https://suimaire.github.io/bioinformatics/pcr-primer-design/`이며, 작업 브랜치 push를 공개 반영 완료로 취급하지 않는다.
