# PCR navigation / site footer integration

검증일: 2026-09-27 (KST). 작업 브랜치: `fix/pcr-site-integration`.

## 1. 기준 branch / commit

시작 시 로컬 `main`은 `1bca1ad`였고 PCR 페이지가 없었다. `git fetch origin`으로 최신 상태를 확인한 뒤, `origin/main`의 `a344f52f8bb997d7319c9ef308e00c3e95eb7ccb`에서 작업 브랜치를 만들었다. 로컬 main은 이동하거나 수정하지 않았다. Main merge 및 공개 배포는 수행하지 않았다.

## 2. Sidebar에 PCR 링크가 생긴 원인

Jekyll은 front matter가 있는 `bioinformatics/pcr-primer-design.html`을 `site.html_pages`에 포함한다. 이 페이지에는 `title`, `nav_order: 3`이 있고 `parent`와 `nav_exclude`가 없었다. 전용 layout 사용 여부와 무관하게 전역 navigation의 최상위 페이지 후보가 된다.

로컬 override인 `_includes/components/sidebar.html`은 홈페이지에서만 `components/site_nav.html`에 `hide_home_nav_branches=true`를 전달한다. 해당 include는 `_config.yml`의 `home_nav_hidden_branches`에 등록된 `생물정보학 기초` branch를 숨긴 후 테마의 `components/nav/pages.html`로 렌더링한다. PCR은 이 branch와 무관한 독립 페이지여서 그대로 노출되었다. 세 상위 분류는 `nav_external_links`에서 생성된다.

## 3. 제거 방식과 hierarchy 보존

PCR front matter에 `nav_exclude: true`만 추가했다. 기존 navigation override는 수정하지 않았다. 실제 Jekyll 결과의 홈페이지 메뉴는 과학 수업 포털과 세 상위 분류, 총 4개 링크다. PCR leaf는 없다.

포털 `#bioinformatics` 본문의 두 번째 자료 링크 및 자동 번호 `1.3.2`는 그대로다. 전용 permalink, title, description, canonical, 페이지의 `1.3.2 생물정보학` 표기 및 `← HAFS Biology Lab 포털` 링크를 유지했다.

## 4. 검색 포함 여부

`search_exclude`를 추가하지 않았다. 실제 빌드의 `assets/js/search-data.json`에 PCR 문서와 절이 포함된다. 홈페이지 검색창에 `PCR`을 키보드로 입력해 검색 결과를 확인하고 결과 클릭으로 해당 페이지에 진입했다. Sidebar 제외와 검색 제외는 독립적으로 동작한다.

## 5. 기존 author identity source

기존 site name은 `_config.yml`의 `title`에 있었고, 동일한 두 줄이 `index.md`의 masthead와 `_includes/footer_custom.html`에도 직접 적혀 있었다.

이번 변경에서는 site name에 `site.title`, 제작자 문구에 새 공통 설정 `site.author_credit`를 사용한다. 두 문구는 각각 **HAFS Biology Lab**, **Teacher-built interactive science tools · CH Park** 그대로다. 포털 masthead와 공통 footer가 동일 설정을 참조한다. PCR 전용 HTML에는 제작자 문구를 중복 작성하지 않았다.

## 6. 기존 counter 구현 / 데이터 소스

- 공통 include: `_includes/footer_custom.html`.
- JS module: `assets/js/page-views.js` (변경 없음).
- DOM mount: `p.page-views[data-page-views]`, 초기 `hidden`.
- 표시 span: `.page-views__today`, `.page-views__sep`, `.page-views__total`.
- 데이터 소스: `https://szbmpsvxzrnewyzqiokr.supabase.co/rest/v1/rpc/`.
- RPC: `record_page_view`, `get_page_view_counts`; body는 `{ "p_page_key": "…" }`.
- 응답 필드: `kst_date`, `today_views`, `total_views`.
- 기존 공개용 publishable key와 REST 요청을 그대로 사용한다. 신규 서비스, DB, SDK, 키를 만들지 않았다.
- 모듈이 참조하는 DB 정의도 기존 인접 저장소의 `supabase/migrations/20260913_page_views.sql`에서 읽어 확인했다. daily 테이블의 key는 `(page_key, view_date)`, total 테이블의 key는 `page_key`다. 날짜 계산은 DB의 `Asia/Seoul` 기준이다. 해당 저장소는 수정하지 않았다.

## 7. Counter의 실제 의미

**기존 counter는 page-specific 조회수다. Site-wide 값이나 고유 방문자 수가 아니다.** PCR은 `/bioinformatics/pcr-primer-design/`, 포털은 `/`로 분리된다. Today는 KST 당일, total은 해당 페이지의 누적 조회수다.

pathname을 소문자로 바꾸고 query/hash, `index.html`, `.html` 등을 정규화해 마지막 `/`를 통일한다. 따라서 PCR의 직접 URL, `index.html` URL, query/hash는 같은 key를 사용한다. 같은 브라우저에서 같은 key의 30분 이내 재방문은 증가 없이 현재 값만 읽는다.

## 8. PCR layout과 footer 연결

기존 `_layouts/pcr-worksheet.html`은 기본 layout을 상속하지 않고 전용 HTML과 학습지 모듈만 출력했다. 공통 footer를 호출하지 않았고 front matter에 `page_views: false`도 있었다.

전용 layout을 유지하면서 `#pcr-worksheet`의 닫는 태그 뒤에 `.pcr-site-footer`를 추가하고 기존 `footer_custom.html`을 한 번 include했다. `page_views: false`를 제거해 기존 기본값을 사용한다. 모든 활동, RNA 확장과 `자료와 계산 범위` 이후에 footer가 나온다. 학습지 모듈의 DOM 범위 밖에 있어 답안/출력 처리와 분리된다.

기존 `.site-credit__name` / `.site-credit__line` 스타일을 `override.css`에서 `site-footer.css`로 옮겨 일반 페이지와 PCR이 함께 로드한다. 기존 13px/12px 크기, 중립색 `#687470`, teal `#0d5751`, 얇은 선 `#e5e3dd` 및 포털 font stack을 사용한다. PCR에만 desktop/tablet 좌우 배치, mobile 세로 배치를 적용했다. 카드, 그림자 또는 배경 박스는 추가하지 않았다.

## 9. Double counting 방지

빌드된 PCR HTML에서 counter module script 1개와 mount 1개를 확인했다. 기본 layout은 호출하지 않으며 별도 수동 init이나 loader를 추가하지 않았다. 기존 모듈은 `window.__hafsPageViewsStarted`와 `DOMContentLoaded`의 `once`로 자동 초기화한다.

두 브라우저 엔진에서 production origin을 모사하되 모든 site/RPC 요청을 가로채 로컬 빌드와 모의 응답으로 검증했다. 첫 load는 `record_page_view` 1회, 30분 이내 reload는 `get_page_view_counts` 1회였다. `index.html?check=1#activity-07`도 같은 PCR key로 읽기만 하고, 포털은 별도 `/` key로 기록했다. 저장소 사용 불가 상태에서도 첫 load 요청은 1회다. 실제 공개 RPC의 증가 요청은 보내지 않았다.

## 10. Loading failure / fallback

기존 방식 그대로, 로딩 중에는 `조회수 불러오는 중…`, 실패 시 counter 줄을 비우고 조용히 숨긴다. Identity는 계속 표시된다. 요청 timeout은 6초이고 오류는 catch된다. HTTP 503, 네트워크 중단, 잘못된 응답, timeout 및 module 자체 로딩 실패를 시험했다.

각 실패 상태에서도 답안 입력/자동 저장, 260 bp primer 계산, RNA 답안과 JSON 내보내기가 동작했다. `pageerror` 및 `unhandledrejection`은 0개다. 의도적으로 실패시킨 요청의 network error와 기존 `console.warn`은 남지만 interaction을 중단하지 않는다.

## 11. 수정 파일

| 파일 | 내용 |
| --- | --- |
| `_config.yml` | 공통 `author_credit` 설정 |
| `index.md` | masthead에서 공통 identity 참조 |
| `bioinformatics/pcr-primer-design.html` | nav_exclude 추가, page_views 비활성화 제거; 본문 불변 |
| `_includes/footer_custom.html` | 공통 설정 참조와 identity wrapper |
| `_includes/head_custom.html` | 공통 footer CSS 로드 |
| `_layouts/pcr-worksheet.html` | 공통 footer include 및 CSS 연결 |
| `assets/css/override.css` | 기존 identity 스타일을 공통 CSS로 이동 |
| `assets/css/site-footer.css` | 기존 identity 스타일과 PCR footer 배치 |
| `tests/pcr-browser.mjs`, `tests/pcr-polish-browser.mjs` | 학습지의 장식 점 금지 검사를 `#pcr-worksheet`로 한정; 요청된 공통 저작자 문구의 `·` 허용 |
| `tests/pcr-preview.mjs` | Ruby 없는 간이 학습지 preview에서도 실제 공통 footer/config 반영 |
| `tests/pcr-site-integration.mjs` | 실제 Jekyll navigation/search/footer 및 counter 실패·중복 검증과 PNG 생성 |
| `tests/PCR_PRIMER_DESIGN_SITE_INTEGRATION.md` | 본 보고서 |

## 12. 테스트 결과

| 검사 | 결과 |
| --- | --- |
| 실제 Jekyll 3.10.0 / 기존 remote theme 빌드 | 16 pages PASS |
| 기존 Node 단위 검사 | 93/93 PASS |
| 기존 PCR 브라우저 회귀 | 548 assertions PASS; Chromium/WebKit × 1440/768/390px |
| 기존 Jekyll 통합 검사 | 번호, 링크, 직접 URL/index.html, reload, 좌표 드래그, Day 1–5, fixture failure PASS |
| 신규 site integration | 178 assertions / 20 cases PASS |
| 반응형 footer / 가로 overflow | 두 엔진 × 3 크기 PASS, overflow 없음 |
| Counter failure / 중복 집계 | 두 엔진 × 7 상태 PASS |
| 런타임 오류 | pageerror 0, unhandledrejection 0 |
| 변경 범위 비교 | front matter 이후 PCR HTML 본문 불변; 19개 보호 자산 불변 |
| 실제 counter 읽기 | 기존 `get_page_view_counts`가 PCR key에 today 0 / total 0 응답 |
| PNG 검수 | desktop/mobile 및 contact sheet를 직접 확인 |

19개 보호 자산은 기존 PCR module 16개, worksheet CSS, teaching fixture, 공통 counter JS다. 활동 00–07, gel, Primer-BLAST 기록, final notebook, RNA, localStorage, JSON schema 및 계산 로직을 수정하지 않았다. 기존 548개 검사는 계산, 저장/복원, JSON round trip, RNA, 인쇄, 기존 기록 및 저장 차단 상태를 포함한다.

재실행 예 (필요한 Ruby gems 및 기존 테스트용 Playwright 설치 후):

```powershell
ruby tests/pcr-build.rb
python -m http.server 4174 --bind 127.0.0.1 --directory tests/.pcr-output/site
# 별도 터미널
node --test tests/*.test.mjs
$env:PCR_TEST_URL = 'http://127.0.0.1:4174/bioinformatics/pcr-primer-design/'
node tests/pcr-browser.mjs
node tests/pcr-integration.mjs
$env:PCR_COUNTER_PREVIEW = 'live' # 기존 RPC 읽기 전용; 생략하면 명시된 mock 모드
node tests/pcr-site-integration.mjs
```

이번 환경은 기존 인접 작업의 Ruby 3.3.12 및 Chromium/WebKit runtime을 재사용했다. Ruby의 선택적 server native dependency 없이 실제 `Jekyll::Site#process`를 사용했으며 Jekyll/테마 소스는 수정하지 않았다. 실행 로그와 요청별 증거는 PNG와 같은 로컬 검증 폴더에 있다.

## 13. PNG 절대 경로

캡처의 `today 0 / total 0`은 기존 counter의 **실제 읽기 전용 응답**이며 hard-code가 아니다. `(읽기 전용)`은 기존 개발 모드가 자동 표시한다. Production에서는 이 suffix 없이 같은 today/total span을 표시한다.

- `D:/Codex/260901 science/suimaire.github.io/verification.local/pcr-primer-design/site-integration/portal-sidebar.png`
- `D:/Codex/260901 science/suimaire.github.io/verification.local/pcr-primer-design/site-integration/portal-bioinformatics-section.png`
- `D:/Codex/260901 science/suimaire.github.io/verification.local/pcr-primer-design/site-integration/pcr-page-footer.png`
- `D:/Codex/260901 science/suimaire.github.io/verification.local/pcr-primer-design/site-integration/pcr-page-full-bottom.png`
- `D:/Codex/260901 science/suimaire.github.io/verification.local/pcr-primer-design/site-integration/pcr-page-mobile-footer.png`
- `D:/Codex/260901 science/suimaire.github.io/verification.local/pcr-primer-design/site-integration/pcr-page-tablet-footer.png`
- `D:/Codex/260901 science/suimaire.github.io/verification.local/pcr-primer-design/site-integration/site-integration-contact-sheet.png`

PNG, `verification.local`, 테스트용 runtime/output, 개인 파일은 commit에 포함하지 않는다.

## 14. 남은 검증 한계

- 이 변경을 main에 merge하거나 GitHub Pages에 배포하지 않았다. 결과는 작업 branch의 실제 로컬 Jekyll 빌드를 기준으로 한다.
- 공개 데이터는 읽기만 했다. 실제 DB 증가/실시간 배포는 이번 검사 범위 밖이며 production 분기와 중복 집계는 요청을 가로챈 브라우저 시험으로 검증했다.
- 30분 제한은 기존 브라우저 storage 기반이다. 다른 브라우저, 저장소 삭제/차단, 동시 탭 사이의 원자적 중복 제거를 보장하는 고유 방문자 추적 방식이 아니다. 기존 동작을 유지했다.
- 실제 장치의 Safari/Firefox는 실행하지 않았다. Chromium/WebKit에서 desktop/tablet/mobile viewport로 검사했다.
- 기존 NanumSquareNeo font 파일이 비어 있어 decode 경고가 있으며 시스템 대체 글꼴로 표시된다. 이번 범위에서 font 자산을 변경하지 않았다.
