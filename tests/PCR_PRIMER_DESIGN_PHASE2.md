# PCR 프라이머 디자인 Phase 2 검증 보고서

검증일: 2026-09-27 (KST). 기준 커밋: `41a02c9`. 작업 브랜치: `feature/pcr-primer-design`.

활동 00~02의 소규모 개선과 활동 03의 sequence workbench 재설계를 완료했다. 기존 계산 모듈, 교육용 fixture, 활동 04~07 및 RNA 확장 HTML은 기준 커밋과 동일하다. `main` 병합이나 push는 수행하지 않는다.

## 1. 사전 확인과 수정 파일

작업 시작 시 현재 브랜치와 기준 커밋이 일치했으며, 기존 미추적 `_codex/`만 있었다. 해당 사용자 파일은 수정하거나 커밋에 포함하지 않았다. 저장소와 상위 경로에서 적용되는 `AGENTS.md`는 발견되지 않았다. HTML/Liquid 레이아웃, CSS, JS, v1 저장 구조, JSON 가져오기/내보내기, fixture, 기존 테스트와 Phase 1 보고서를 확인했다.

| 파일 | 수정 내용 |
| --- | --- |
| `.gitignore` | `verification.local/`의 Git 추적 제외 |
| `_config.yml` | 로컬 검증 PNG/JSON의 Jekyll 배포 제외 |
| `bioinformatics/pcr-primer-design.html` | 활동 02 과학 문구, 주기 비교 SVG 크기, 활동 03 workbench 구조 |
| `assets/css/pcr-worksheet.css` | 00 marker 연결선, 03 화면/터치/인쇄 스타일 |
| `assets/js/pcr-intro.mjs` | 00 위치의 접근성 설명, 1~3주기 비교 그림 |
| `assets/js/pcr-design.mjs` | 새 파일. 선택 좌표 해석, 기존 서열의 위치 확인, 결실/배치 설명 adapter |
| `assets/js/pcr-workbench.mjs` | 새 파일. 지도, 확대 서열, 선택 조작, 실시간 설계/산물 표시 |
| `assets/js/pcr-records.mjs` | 선택적인 보기 상태와 결합 위치 기록의 검증/호환 |
| `assets/js/pcr-worksheet.mjs` | workbench와 자동 저장, 기존 계산/gel, 설계 1~3 연결 |
| `tests/pcr-design.test.mjs` | 새 파일. 선택/결실/잘못된 설계/이전 기록의 단위 테스트 8개 |
| `tests/pcr-workbench-browser.mjs` | 새 파일. 402개 검증, 실제 조작 후 PNG 생성과 검증 manifest |
| `tests/pcr-browser.mjs` | 기존 500개 회귀 검증을 새 선택/접힌 보조 입력에 연결 |
| `tests/pcr-integration.mjs` | 기본 공개 상태인 확대 서열에서 드래그 검증 |
| `tests/pcr-static.test.mjs` | 새 모듈에도 금지 문자 및 안전한 DOM 표시 검사 적용 |
| `tests/PCR_PRIMER_DESIGN_PHASE2.md` | 이 보고서 |

## 2. 활동 00 변경

기존 slider와 F/R marker의 `input` 즉시 연동을 유지하고, marker와 DNA 선 사이에 가는 연결선을 추가했다. 구조도의 접근성 설명도 현재 F/R의 대략적 위치와 함께 갱신한다. bp 좌표, 주문 서열, GC, Tm을 계산하지 않는다. 초기 예측은 기존 `initialPrimerPrediction`에 별도로 보존한다. 03의 `00의 초기 예측 보기` 버튼은 이 값만 읽어 옅은 회색 점선 marker를 표시한다.

## 3. 활동 01 변경

PCR cycle viewer는 유지했다. 아래 비교 그림에서 1주기는 원래 긴 주형과 primer에서 시작해 반대 primer 위치를 넘어가는 새 가닥을 구별한다. 2주기는 그 긴 새 가닥이 주형이 되고, 반대 primer에서 합성한 가닥이 이미 정해진 주형 끝에 도달하는 모습을 그린다. 3주기는 양 끝이 정해진 이중가닥과 목표 길이 bracket을 보여 준다. 실제 분자 비율을 나타내지 않는 교육용 모형이라는 설명과 효율 한계 문구는 유지했다.

## 4. 활동 02 변경

3′ 말단 자체가 최종 산물의 바깥 경계를 정한다는 오해를 피하도록 다음 의미로 고쳤다: 두 3′ 말단이 서로를 향해야 각 primer에서 시작하는 합성이 두 결합 부위 사이의 표적 영역을 향해 진행된다. 5′/3′ 그림, reference strand, 가닥 뒤집기와 complement → reverse complement 학습을 유지했다.

## 5. 활동 03 구조와 primer 선택 UX

- 넓은 화면은 약 2:1의 DNA workbench / 현재 설계 패널이다. 작은 화면에서는 세로로 쌓인다. 기존 목차와 함께 폭이 좁아지는 중간 desktop 구간도 세로 배치한다.
- 전체 지도는 fixture에서 A의 길이와 B의 결실 구간을 읽고, F/R 위치, 합성 방향, 선택 위치에서 정의되는 amplicon bracket을 표시한다.
- 확대 서열은 desktop/tablet 60 nt, mobile 30 nt만 표시한다. 염기별 좌표, monospace 서열, 최소 44px 선택 칸, 이전/다음 버튼과 위치 slider를 제공한다. 결실 영역의 염기는 중립색 점선으로 구별한다.
- F/R 선택 모드, 두 지점 클릭/터치, 마우스/펜 드래그, 화살표 이동 + Enter/Space, Shift+화살표 확장, 확대 경계를 넘는 키보드 선택을 지원한다. Escape는 진행 중인 선택을 취소한다.
- 좌표와 합성 방향은 `좌표로 미세 조정`에 접었다. 입력 즉시 갱신하며, 적용 버튼은 선택 위치로 확대 화면을 옮긴다. 시작/끝 역전, 범위 밖 값, 잘못된 F/R 배치나 합성 방향을 자동 교정하지 않는다.
- F와 R 주문 서열은 모두 5′→3′다. R에는 reference 결합 부위와 역상보 변환 과정을 함께 표시한다. 기본 R 선택은 자동으로 역상보를 계산한다. 보조 도구에서 합성 방향을 바꾼 경우 그 방향을 그대로 계산/표시한다.
- 현재 설계는 좌표, 길이, GC, 주문 서열, A/B/C 산물과 결실 관계를 즉시 보여 준다. 03에는 gel을 그리지 않는다. Tm 및 상보성 검토는 04에 둔다.
- 저장 설계는 최대 3개이며, 좌표/서열/A/B/C 결과와 짧은 변경 이유로 비교한다. 불러오기는 편집본만 바꾸고 저장된 snapshot은 보존한다. 이전 큰 답안 필드들은 접힌 추가 설명에 남겨 과거 기록을 지우지 않는다.

## 6. 계산 로직과 과학적 확인

`pcr-core.mjs`와 fixture는 변경하지 않았다. 새 adapter가 `bindingSites`, `primerStats`, `reverseComplement`, `productFromHits`, `analyzeAll`을 재사용한다. 핵심 탐색은 모든 F/R, F/F, R/R의 완전 일치 조합을 계속 계산한다. 지도는 학생이 선택한 한 쌍, 예상 산물은 실제 주문 서열의 모든 완전 일치 결과로 구분한다. 여러 결합 위치를 가진 이전 서열에는 임의의 대표 좌표를 부여하지 않는다.

| 실제 검사 배치 (A 좌표) | A | B | 관찰 |
| --- | --- | --- | --- |
| F 41–60 / R 281–300 | 260 bp | 180 bp | 결실이 primer 사이에 포함됨 |
| F 111–130 / R 281–300 | 190 bp | 예상 산물 없음 | F 결합 부위가 결실과 부분적으로 겹침 |
| F 121–140 / R 281–300 | 180 bp | 예상 산물 없음 | F 결합 부위가 결실에 완전히 포함됨 |
| F 231–250 / R 281–300 | 70 bp | 70 bp | 두 primer 모두 결실 오른쪽 |

B는 A 산물에서 단순히 80을 빼지 않고 B 서열에서 직접 탐색한다. 결실이 선택한 결합 부위를 손상시키지만 다른 위치의 완전 일치 결합이 있는 경우에는 대안 결합을 따로 설명한다. 고정 P1/P2/P3의 A/B/C 9개 기준 결과도 모두 유지된다. 범위 밖/뒤집힌 좌표와 100 nt 초과 선택은 현재 입력을 유지하면서 계산 제한을 안내한다. 이런 초안은 JSON/localStorage에 보존되며, 기존 제한을 넘어서는 주문 서열은 저장 설계로 확정하지 않는다. 짧은 primer도 허용하고 길이 검토 안내를 표시한다.

## 7. 저장 state와 이전 JSON/localStorage 호환

`schemaVersion: 1`, `dataVersion`, localStorage 키 `hafs:pcr-primer:v1`, 기존 답안과 `forward`/`reverse` 문자열을 유지했다. 추가 필드는 선택적이다.

```json
{
  "workbench": { "mode": "R", "windowStart": 281, "showPrediction": true },
  "draft": {
    "forward": "AGTCGATGCTACGTTGACCA",
    "reverse": "TTCAGGCTACGATCGTACGA",
    "bindings": {
      "F": { "start": "41", "end": "60", "direction": "right" },
      "R": { "start": "281", "end": "300", "direction": "left" }
    }
  }
}
```

`designs[]`에도 저장 시점의 선택적인 `bindings`를 복사한다. 좌표는 문자열로 남겨 미완성/잘못된 편집 값을 보존한다. 정상 범위로 몰래 바꾸지 않는다. 구조, 방향 열거값, 길이와 보기 설정은 import 시 검증한다. 좌표와 주문 서열이 맞지 않는 기록은 이전 결과를 표시하지 않고 불일치를 안내한다.

이전 v1에는 기본 보기 설정을 제공하고, 서열이 A에서 단 하나의 결합 위치를 갖는 경우에만 UI에 위치를 확인해 표시한다. 기존 저장 설계의 서열/설명/시각은 보존된다. 직접 서열을 편집하면 해당 primer의 선택 메타데이터만 해제해 새 서열과 낡은 좌표가 혼용되지 않게 한다. 기존 04 후보 비교, 05 gel, 06 기록, 07 최종 비교가 같은 서열 필드를 계속 사용한다.

## 8. 테스트와 실행 결과

| 검사 | 결과 |
| --- | --- |
| `node --test tests/*.test.mjs` | 50개 통과: 기존 42개 + 새 8개 |
| PCR JS 모듈 5개 구문 검사 | 통과 |
| 기존 Chromium/WebKit 브라우저 회귀 | 500개 assertion 통과 |
| 새 workbench 브라우저 검증 | 402개 assertion 통과 |
| 브라우저/viewport | 두 엔진 각각 1440×1150, 768×1150, 390×1150 (기존 suite 높이는 1000) |
| Jekyll 3.10.0 실제 빌드 | 기존 원격 테마로 16개 페이지 생성 |
| 포털 통합 | 번호 1.3.2, 직접 주소/새로고침, 드래그, Day 1~5 링크, fixture 실패 시 기록/JSON 유지 통과 |
| 보호 범위 | 04 이후 HTML, 계산 core, fixture가 기준 커밋과 동일 |
| 저장 호환 | 이전 v1, 3개 snapshot, 좌표/방향, 미완성/잘못된 초안, JSON 왕복과 reload 통과 |
| 접근성 | 보이는 focus, 모드 `aria-pressed`, 명확한 염기 label, 키보드/터치, 44px target, reduced motion, 가로 overflow 검사 통과 |
| 인쇄 회귀 | 작성본/빈 학습지 구분, 긴 답안, 저장 기록 및 marker 숨김, 인쇄 후 원래 기록 유지 통과 |
| Git diff/출력 제외 | `git diff --check` 통과. PNG는 Git 및 Jekyll 배포에서 제외 |

새 검증은 F/R 선택, 두 지점 조작, drag, 키보드와 확대 경계 이동, 좌표 미세 조정, 역상보/GC, 결실 포함/부분 겹침/완전 겹침, 합성 방향 오류/서로 겹친 primer/범위 오류/길이 한도, A/B/C, 저장 1~3과 복원, JSON/localStorage, 00 overlay를 확인한다. 기존 suite는 00~02 학습 조작과 04~07, 이전 기록, 저장 차단, 잘못된 import, 초기화 취소, no-JS 및 인쇄를 계속 검증한다.

브라우저 suite를 동시에 실행하던 중 새 suite의 초기 페이지 로딩이 한 번 시간 초과했다. 실제 Jekyll 출력을 대상으로 직렬 재실행하여 전 항목 통과를 확인했다. 초기 로딩 실패 시 저장/계산 안내와 네트워크 오류를 출력하도록 검증 도구를 보강했다.

재실행 예시:

```powershell
$env:PLAYWRIGHT_BROWSERS_PATH = (Join-Path (Get-Location) 'tests/.pcr-tools/browsers')
$env:PCR_TEST_URL = 'http://127.0.0.1:4174/bioinformatics/pcr-primer-design/'
node --test tests/*.test.mjs
node tests/pcr-browser.mjs
node tests/pcr-workbench-browser.mjs
node tests/pcr-integration.mjs
```

4174는 실제 Jekyll 출력인 `tests/.pcr-output/site/`를 제공하는 서버다. 학습지 전용 미리보기에서는 `node tests/pcr-preview.mjs` 실행 후 `PCR_TEST_URL`을 4173으로 바꾼다. 이 Windows 환경의 빌드는 Phase 1에서 준비된 Ruby와 실제 Jekyll `Site#process` 도구를 사용했다. 원격 테마 다운로드에 필요한 네트워크 권한을 사용했다.

## 9. 자동 생성 PNG와 절대 경로

아래 12개 PNG는 Playwright가 실제 slider/버튼/염기/키보드/좌표 입력을 조작한 후 생성했다. 테스트용 overlay를 넣지 않았다. 활동 또는 해당 도식/비교 영역 전체를 캡처하여 화면 아래 정보도 확인할 수 있다. PNG는 커밋하지 않는다.

| PNG | 절대 경로 |
| --- | --- |
| 00-overview.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase2/00-overview.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase2/00-overview.png>) |
| 01-cycle-1.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase2/01-cycle-1.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase2/01-cycle-1.png>) |
| 01-cycle-3.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase2/01-cycle-3.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase2/01-cycle-3.png>) |
| 02-direction.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase2/02-direction.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase2/02-direction.png>) |
| 03-workbench-empty.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase2/03-workbench-empty.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase2/03-workbench-empty.png>) |
| 03-forward-selected.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase2/03-forward-selected.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase2/03-forward-selected.png>) |
| 03-complete-design.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase2/03-complete-design.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase2/03-complete-design.png>) |
| 03-deletion-overlap.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase2/03-deletion-overlap.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase2/03-deletion-overlap.png>) |
| 03-saved-designs.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase2/03-saved-designs.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase2/03-saved-designs.png>) |
| 03-invalid-design.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase2/03-invalid-design.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase2/03-invalid-design.png>) |
| 03-tablet.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase2/03-tablet.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase2/03-tablet.png>) |
| 03-mobile.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase2/03-mobile.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase2/03-mobile.png>) |

화면별 의미: 00은 F/R 예측 marker, 01은 첫 주기 긴 산물/세 번째 주기 정확한 이중가닥, 02는 수정 문구와 상보→역상보 활동, 03은 빈 설계/F만 선택/정상 설계/결실 경계 겹침/잘못된 F/R 배치/두 저장 설계/모바일/태블릿이다. 정상 설계 PNG에는 오른쪽 계산값과 00 초기 예측 overlay도 포함한다.

같은 출력 폴더의 `verification.json`은 최종 실행 URL, 시각, 엔진/viewport별 성공 결과와 PNG 절대 경로를 기록한다. JSON 왕복 검증 산출물도 해당 제외 폴더에 있다.

## 10. 수동 검증 한계와 브랜치 전달

실제 iPhone/iPad와 macOS Safari, 화면 읽기 프로그램의 실제 발화, 펜 장치, 실제 프린터는 수동 검증하지 않았다. Chromium/WebKit 및 touch/keyboard 자동화, 인쇄 PDF 생성/DOM 검사로 확인했다. 대표 desktop/mobile PNG, 주기 도식, 결실 겹침/잘못된 배치와 저장 비교를 시각 검수했다. 실제 PCR 성공률, 열역학, 유전체 특이성이나 외부 Primer-BLAST 검색은 이번 구현 범위가 아니다.

모든 검증 결과를 확인한 후 이 변경과 보고서를 `feature/pcr-primer-design`에 커밋하고 해당 브랜치만 원격에 push한다. `main` 병합/push나 공개 사이트 배포는 수행하지 않는다.
