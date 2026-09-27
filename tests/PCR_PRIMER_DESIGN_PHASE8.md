# PCR과 프라이머 디자인 Phase 8 검증 보고서

검증일: 2026-09-27 (KST). RNA 확장 활동만 재설계했다. Core Course 00~07의 기존 계산, 데이터 및 핵심 DOM은 유지했다.

## 1. 기준 commit / branch / 원격 상태

- 작업 branch: `feature/pcr-primer-design`.
- 시작 HEAD: `b15b2e23572f44ffff7191092cba56bc0c8073d0` (Phase 7).
- 시작 상태: `?? _codex/`만 존재. 해당 사용자 폴더를 수정하거나 staging하지 않았다.
- 저장소와 상위 `D:/`, `D:/Codex/`에 적용할 AGENTS.md가 없었다.
- 시작 remote-tracking 비교: 0 behind / 0 ahead. 이어서 실제 `git ls-remote`로 원격 feature도 위 Phase 7 commit임을 확인했다.
- 읽기 시점 원격 main: `1bca1ad14247c2a1973df3981e3e57c8b0190b2a`.
- main merge/push는 하지 않는다. 최종 Phase 8 commit hash 및 feature push 성공 여부는 완료 응답에 기록한다. 이 보고서가 자신의 commit hash를 포함하게 만들지는 않는다.

## 2. 수정 파일과 보호 범위

| 파일 | 변경 |
| --- | --- |
| `bioinformatics/pcr-primer-design.html` | A~E/reflection 확장 구조, 짧은 질문, legacy 원문, 06 연결 링크 |
| `assets/js/pcr-rna.mjs` | 독립 RNA 상태/검증, primer 구조 설명, transcript exon 인접성 판별 |
| `assets/js/pcr-rna-view.mjs` | 단계/설계/배치/RT case/transcript 조작, 그림과 텍스트 대안, 인쇄 복원 |
| `assets/js/pcr-records.mjs` | optional v1 RNA 상태와 legacy migration 연결 |
| `assets/js/pcr-worksheet.mjs` | RNA 초기화, restore, 기존 인쇄 수명 주기에 연결 |
| `assets/css/pcr-worksheet.css` | RNA 범위의 화면/도식/모바일 규칙을 뒤에 추가 |
| `tests/pcr-rna.test.mjs` | 상태, migration, 잘못된 입력, 구조적 검출 범위 등 7개 단위 테스트 |
| `tests/pcr-rna-browser.mjs` | 실제 UI 조작, 저장/JSON/인쇄/접근성/viewport, PNG/contact sheet 생성 |
| `tests/pcr-browser.mjs` | 기존 RNA 입력/내보내기/인쇄 검사를 새 reflection 위치로 연결 |
| `tests/PCR_PRIMER_DESIGN_PHASE8.md` | 이 보고서 |

Phase 7과 소스를 비교해 Core 00~07 HTML이 06의 `ext-rna-help` 연결 부분을 제외하면 동일함을 확인했다. 기존 CSS 전체가 새 CSS의 prefix이며, 기존 print CSS도 그대로다. Core/design/workbench/intro/review/evidence/external/final 계산 모듈, 교육용 fixture, layout, portal은 수정하지 않았다. 새 외부 API, 서비스, 라이브러리, 프레임워크를 추가하지 않았다.

## 3. 기존 RNA 활동 조사

기존 `#rna-extension`은 접힌 확장 활동으로, 두 exon/gDNA/cDNA SVG, junction/intron 설명, RT(-)/qPCR 짧은 문장, `answers['rna-plan']` textarea 하나와 해설로 이루어져 있었다. 독립 RNA JS/state는 없었고 공통 `[data-answer]` 자동 저장과 인쇄를 사용했다.

Phase 7 보고서, 전체 state validator, localStorage/JSON 흐름, 01 cycle, 02 방향, 03 workbench, 05 controls, 06 RNA 안내, 07 notebook, 공통 UI와 모든 기존 테스트를 조사했다. 기존 목차의 ‘확장 활동 / RNA 발현을 보려면?’ 및 `#rna-extension` anchor를 유지했다. Core 사이에 필수 번호를 추가하지 않았다.

## 4. RNA → cDNA 구현

DNA template → PCR → amplicon과 RNA → reverse transcription → cDNA → PCR 또는 qPCR → amplicon을 비교한다. 네 단계 버튼을 눌러 현재 단계와 짧은 설명을 확인한다. 긴 textarea 대신 RNA/cDNA/단백질/dNTP 선택 질문을 두었다.

일반적인 RT-PCR에서 reverse transcription으로 만든 cDNA를 PCR template로 사용한다고 설명한다. 모든 DNA polymerase에 대한 보편적 금지 표현은 쓰지 않았다. RT-PCR의 RT와 one-step workflow는 짧은 접힌 보충 설명으로 둔다. 별도 animation은 없다.

## 5. Genomic DNA 비교

Genomic DNA의 Exon 1–intron–Exon 2–intron–Exon 3과 mature mRNA/cDNA의 연결된 Exon 1–2–3을 비교한다. Splicing된 mature RNA의 exon 연결이 cDNA에 반영된다는 점을 primer 설계와 연결한다.

gDNA가 시료에 남을 가능성을 다루며 항상 오염되어 있다고 단정하지 않는다. 같은 primer pair가 두 template에 결합하면 신호를 RNA 유래만으로 해석하기 어려울 수 있음을 설명한다.

## 6. Same-exon strategy

F/R을 같은 Exon 1 내부에 마주보게 표시한다. cDNA/gDNA에 같은 결합 부위가 있어 동일 크기 product가 가능하다는 모형이다. 목적에 따라 의미가 달라지므로 나쁜 설계나 자동 부적합 판정을 내리지 않는다.

## 7. Intron-spanning strategy

F를 Exon 1, R을 Exon 2에 놓는다. cDNA에서는 상대적으로 짧고, gDNA에서는 intron을 포함하는 더 긴 product가 가능하다. 짧은/긴 intron 두 예를 선택한다. 도식은 길이 축척이 아니며 label과 결과 설명이 바뀐다.

긴 intron이 현재 조건에서 증폭을 불리하게 하거나 크기 구별을 도울 수 있다고 설명한다. 실제 결과는 intron 길이, extension time, polymerase와 PCR 조건에 따라 달라진다. 성공 확률이나 ‘gDNA 증폭 불가능’을 계산하지 않는다.

## 8. Junction strategy

R을 Exon 2에 고정하고 F를 Exon 1 내부 또는 Exon 1–2 junction에 직접 배치하는 작은 버튼 조작이다. Junction 상태에서 하나의 F 화살표가 cDNA의 연결 경계를 실제로 가로지른다. gDNA 쪽에는 떨어진 결합 서열 부분을 점선으로 표시한다.

해당 intron을 포함하는 gDNA에 동일한 연속 F 결합 부위가 없다는 교육용 예다. 실제 특이성에는 primer 서열, genome context, 유사 유전자/가공 위유전자 등의 검토가 필요하며 gDNA 증폭 배제를 보장하지 않는다고 표시한다. 03 전체 workbench를 복제하지 않았다.

## 9. RT(+)/RT(-)

Reverse transcriptase를 넣거나 생략한 같은 RNA sample을 비교한다. Case 1은 RT(+)에만 band, Case 2는 둘 다 band가 있는 고정 교육 자료다. 작은 signal table을 사용해 별도 gel component를 만들지 않았다.

Case 2는 RNA-derived cDNA가 아닌 DNA template의 존재, 예를 들어 gDNA carryover 등을 검토할 근거다. gDNA contamination으로 확정하지 않고 기타 오염/실험 문제도 검토한다. Case 1 역시 DNA가 전혀 없다는 증명이나 산물 정체 확인으로 해석하지 않는다. 05 control 해석과 연결했다.

## 10. Transcript variant

Gene X의 exon 구조를 T1=1–2–3–4, T2=1–2–4, T3=1–3–4로 표시한다. R은 공통 Exon 4 하류에 고정하고 F의 영역을 선택한다. Exon 4 선택 시 F는 같은 exon 내 R 상류라는 전제를 명시했다.

| F 선택 | 결합 구조가 있는 transcript |
| --- | --- |
| Exon 2 | T1, T2 |
| Exon 3 | T1, T3 |
| Exon 4 | T1, T2, T3 |
| Exon 2–3 junction | T1 |
| Exon 3–4 junction | T1, T3 |

각 행에 실제 exon 구성과 해당 exon/junction 존재 또는 부재의 이유를 같이 표시한다. Junction은 두 exon의 단순 존재가 아니라 인접성을 확인한다. 공통 exon에서는 여러 transcript 신호가 합쳐질 수 있지만 gene 전체의 절대 발현량을 뜻하지 않는다. 특정 exon/junction도 실제 isoform specificity를 보장하지 않으며 전체 transcriptome/서열 검토가 필요하다.

Gene 측정 범위에 관한 3줄 질문과 특정 transcript를 위한 2줄 질문을 두고, 마지막 4줄 reflection은 ‘내 primer가 실제로 무엇을 측정하는가?’로 마무리한다.

## 11. qPCR 범위 제한

하단 접힌 ‘더 알아보기 / qPCR에서 발현량을 비교하려면?’에는 정량을 위한 추가 고려사항 한 문단만 있다. Ct/Cq 계산, ΔCt/ΔΔCt, reference gene activity, efficiency 계산기, probe chemistry, standard curve를 구현하지 않았다. RT-PCR을 정량 PCR과 동일시하거나 qPCR이 RNA 자체를 직접 증폭한다고 설명하지 않는다.

## 12. 06과 중복 정리

기존 `ext-rna-help`의 긴 설명을 ‘RNA 발현 분석에서 exon/intron을 고려하려면 확장 활동을 참고하세요’ 링크로 교체했다. 클릭하면 접힌 확장 활동도 열린다. 기존 검색 route, 조건/후보 기록, 실제 Primer-BLAST 링크/기능은 변경하지 않았다.

## 13. State 변경

`schemaVersion: 1`, `dataVersion`, `hafs:pcr-primer:v1`을 유지했다. Optional `rnaExtension`을 추가한다.

```text
rnaExtension
  activeSection: flow | genomic | strategy | control | transcripts | reflection
  flowStep: rna | rt | cdna | pcr
  primerStrategy: same-exon | intron-spanning | junction
  intronExample: short | long
  junctionPosition: exon | junction
  rtControlCase: case1 | case2
  transcriptTarget: exon2 | exon3 | exon4 | junction23 | junction34
  answers: template, strategy, control, transcripts, specific, reflection
```

누락된 RNA 상태/필드는 기본값을 적용한다. 문자열 답안은 기존과 같은 최대 12,000자다. 잘못된 enum, 타입, 답안 key, 알 수 없는 RNA state key는 import 전에 거절한다. Core 답안/설계/계산 state를 변경하지 않는다. Hover나 임시 animation 상태는 저장하지 않는다.

## 14. Legacy migration

이전 `answers['rna-plan']`을 삭제하지 않는다. RNA state가 없는 기록의 첫 migration에서 의미가 대응되는 새 `rnaExtension.answers.reflection`에 복사한다. 원문은 읽기 전용 ‘이전 RNA 확장 기록’ 접힌 영역에 남긴다.

새 state가 이미 있으면 migration을 반복하지 않아 수정한 답안 또는 의도적으로 비운 답안을 덮어쓰지 않는다. 원문의 개행/특수문자도 보존한다. 새 reflection을 수정해도 이전 원문과 Core 답안은 바뀌지 않는다.

## 15. JSON / localStorage

단계/primer strategy/배치/RT case/transcript 영역/reflection은 기존 저장 함수를 통해 자동 저장한다. 기존 export/import에 같은 JSON으로 포함한다.

미사용, 일부 진행, 완료 RNA 기록의 export/import를 검사했다. 완료 기록의 reload, export 후 임시 수정, re-import로 선택과 답안을 복원한다. 잘못된 RNA import는 현재 기록을 유지한다. 입력은 DOM textContent/value로만 처리해 HTML처럼 실행하지 않는다.

Phase 1~6은 기존 역사 fixture를 localStorage와 실제 파일 가져오기 UI 양쪽으로 검사했다. Phase 7은 Phase 6와 record schema/validator를 바꾸지 않았으므로 Phase 6 fixture에 이전 RNA 답안을 넣은 재구성 Phase 7 record를 사용했다. 실제 학생 기록이라고 주장하지 않는다. Core state는 각 역사 parser 호환 결과와 동일하고 이전 RNA 원문도 보존된다.

## 16. Responsive

1440px desktop, 768px tablet, 390px mobile을 두 엔진에서 검사했다. 기존 Phase 7 전체 회귀는 1024px도 포함한다. 모바일에서 flow/strategy 선택을 세로 배치한다.

Exon 도식은 최소 540px의 읽을 수 있는 내부 canvas를 유지하고 해당 그림 영역만 가로 스크롤한다. 페이지 전체를 축소하지 않는다. Transcript 행은 각 exon block과 label이 보이는 flex 구조다. 모든 RNA 절에서 page horizontal overflow가 없고 그림 영역을 키보드로 스크롤할 수 있다.

## 17. 접근성

Native button/select/textarea, 연결 label, aria-pressed/aria-controls, polite feedback, visible focus, 44px 이상 touch target을 검사했다. 현재 절은 하나만 노출하며 활성 절을 저장한다. Exon/intron/junction/결합 유무는 label, 선, 구조 설명과 함께 제공하므로 색만으로 구별하지 않는다.

SVG는 각 template와 primer 위치를 설명하는 접근성 이름을 갖고, transcript별 exon 구성/선택 이유도 text alternative에 포함한다. Reduced-motion에서 새 animation은 없다.

RNA 36회 axe 검사에서 violations 0. 이 중 12회는 SVG text의 color-contrast를 자동 확정하지 못해 incomplete이며 나머지 오류로 숨기지 않았다. 실제 PNG와 지정 색을 확인했다. 진한 텍스트/중립 윤곽, 흰 배경, 기존 청록 primer만 사용한다. 기존 전체 화면 Phase 7 axe 24회도 violations 0. 실제 screen reader 발화와 물리 터치는 수동 검증 한계다.

## 18. 테스트

| 검증 | 결과 |
| --- | --- |
| 전체 `node --test tests/*.test.mjs` | 90개 통과: 기존 83 + RNA 7 |
| 기본 browser / intro / 저장 / 인쇄 / no-JS | 548 assertions |
| Phase 2 workbench | 402 assertions |
| Phase 3 review | 720 assertions |
| Phase 4 evidence | 926 assertions |
| Phase 5 external | 512 assertions |
| Phase 6 final notebook | 670 assertions |
| Phase 7 polish | 768 assertions, 전체 화면 axe 24회 |
| Phase 8 RNA | 554 assertions, RNA axe 36회 |
| 브라우저 assertions 합계 | 5,100 |
| 실제 Jekyll 3.10.0 | 16-page build 성공 |
| 통합 | 포털 번호, 직접 URL/refresh, 드래그, Day 1~5, fixture failure 통과 |
| 정적 확인 | syntax, 금지 표현/안전한 DOM, Core 보호 범위, diff whitespace 통과 |

Phase 8 신규 검사에는 단계/cDNA 답안, 두 구조, 세 strategy, 두 intron 예, F junction 배치, RT case, 5개 transcript 영역의 15개 구조 판별, 모든 답안/reload, 미사용/부분/완료 JSON, 역사 기록, XSS inert text, invalid import 보존, print와 키보드/overflow가 포함된다.

초기 병렬 회귀 중 06 legacy 페이지 준비에 timeout이 한 차례 발생했다. 같은 기록을 별도로 재현했을 때 정상 초기화됐고 이후 06 전체 512 assertions가 통과했다. 테스트 timeout 값을 늘리거나 검사를 생략하지 않았다.

## 19. 전체 검증 결과와 과학적 한계

완료 범위는 RNA 확장, Core 회귀, 저장 호환, 반응형/접근성, 기존 인쇄 회귀, 자동 PNG/contact sheet다. 실제 Jekyll 산출물 `http://127.0.0.1:4174/bioinformatics/pcr-primer-design/`에서 검증했다. 기존 공개 remote theme를 읽기 위한 네트워크 권한으로 빌드했으며 production 설정은 변경하지 않았다.

과학적으로 단순화한 부분은 고정 Gene X exon 구성, 길이/실서열 없는 도식, exon block 기반 결합 구조 판별, 고정 RT band 예다. 실제 조건 의존적인 부분은 gDNA product 증폭 여부, intron 크기 구분, junction primer의 실제 특이성, isoform specificity 및 band의 정체다. 후자는 확인 완료나 실험 성공으로 자동 판정하지 않는다. PCR/qPCR의 실제 반응이나 정량을 실행하지 않았다.

참고 확인: [NCBI Primer-BLAST exon/intron 및 variant 설정](https://www.ncbi.nlm.nih.gov/tools/primer-blast/index.cgi?GROUP_TARGET=on), [NCBI Primer-BLAST 원 논문](https://pmc.ncbi.nlm.nih.gov/articles/3412702/), [Thermo Fisher RT-PCR control 안내](https://assets.thermofisher.com/TFS-Assets/LSG/manuals/MAN0012720_20x20uL_rxns_Maxima_H_Minus_FirstStrand_cDNA_Syn_UG.pdf). 이 자료는 개념과 한계를 확인하기 위한 것으로 새 runtime API가 아니다.

인쇄는 기존 스타일을 재설계하지 않았다. 기존 textarea/select mirror에 RNA 답안도 포함하며, 인쇄 중만 전체 RNA 절을 노출하고 이후 현재 절/접힘을 복원한다. 작성본은 새 답안과 이전 원문을 보존하고 빈 학습지는 답안 값을 비운다. 기본 PDF 생성 회귀와 RNA 전용 인쇄 미디어 검사를 통과했다. 페이지 분할/인쇄 미감 polish는 하지 않았다.

## 20. PNG 절대 경로

모든 PNG는 Playwright의 실제 버튼/선택/입력 조작 이후 자동 생성했다. 보안 검증 문자열과 별도로 마지막 reflection 대표 PNG에는 학습 답안 예시만 표시했다.

| PNG | 절대 경로 |
| --- | --- |
| RNA-overview.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase8/RNA-overview.png` |
| RNA-cdna-flow.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase8/RNA-cdna-flow.png` |
| RNA-genomic-vs-cdna.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase8/RNA-genomic-vs-cdna.png` |
| RNA-strategy-same-exon.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase8/RNA-strategy-same-exon.png` |
| RNA-strategy-intron-spanning.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase8/RNA-strategy-intron-spanning.png` |
| RNA-strategy-junction.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase8/RNA-strategy-junction.png` |
| RNA-rt-control.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase8/RNA-rt-control.png` |
| RNA-rt-minus-signal.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase8/RNA-rt-minus-signal.png` |
| RNA-transcript-variants.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase8/RNA-transcript-variants.png` |
| RNA-transcript-selection.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase8/RNA-transcript-selection.png` |
| RNA-final-reflection.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase8/RNA-final-reflection.png` |
| RNA-mobile.png | `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase8/RNA-mobile.png` |

같은 폴더에 개별 절의 `-mobile.png`도 생성했다. `RNA-mobile.png`는 모바일 실제 viewport, 개별 절 PNG는 그 절 전체 캡처다.

## 21. Contact sheet 절대 경로

- `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase8/phase8-contact-sheet.png`
- `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase8/phase8-contact-sheet-mobile.png`

두 contact sheet는 2열로 RNA→cDNA, 세 설계, RT(-) signal, 공통 exon 및 특정 junction 선택을 비교한다. 원본 크기는 desktop 1592×4471px, mobile 852×7756px다. 확대해 읽을 수 있는 PNG이며 개별 원본도 함께 제공한다. 직접 열어 도식/글자/겹침과 선택 상태를 확인했다.

구조화된 결과는 같은 폴더의 `verification.json`, 실행 요약은 `unit.log`, `base.log`, `workbench.log`, `review.log`, `evidence.log`, `external.log`, `final.log`, `polish.log`, `integration.log`, `rna.log`에 있다. Phase 2~6 기존 산출물은 `phase8/regression/`에 분리했다. 기존 Phase 7 script의 contact sheet는 기존 phase7 출력 경로를 사용한다.

## 22. 남은 수동 검증 한계 / 실행과 commit 범위

- NVDA/VoiceOver의 실제 발화, 실기기 Safari/터치, 교실 투사와 학생 사용성은 자동 검사로 대체하지 않았다.
- SVG contrast incomplete는 PNG/명시된 색을 별도로 확인했으나 axe의 완전 자동 통과라고 보고하지 않는다.
- 실제 PCR, DNA 잔존 검출한계, full transcriptome specificity 및 정량 해석은 구현 범위 밖이다.
- 인쇄 기능의 값/접힘 회귀만 검증했으며 새 pagination/인쇄 디자인 polish는 수행하지 않았다.

재현 명령:

```powershell
$env:PLAYWRIGHT_BROWSERS_PATH = Join-Path (Get-Location) 'tests/.pcr-tools/browsers'
$env:PCR_TEST_URL = 'http://127.0.0.1:4174/bioinformatics/pcr-primer-design/'
$env:PCR_VERIFICATION_ROOT = Join-Path (Get-Location) 'verification.local/pcr-primer-design/phase8/regression'
node --test tests/*.test.mjs
node tests/pcr-browser.mjs
node tests/pcr-workbench-browser.mjs
node tests/pcr-review-browser.mjs
node tests/pcr-evidence-browser.mjs
node tests/pcr-external-browser.mjs
node tests/pcr-final-browser.mjs
node tests/pcr-polish-browser.mjs
node tests/pcr-rna-browser.mjs
node tests/pcr-integration.mjs
```

Source/test/report만 명시적으로 staging한다. `verification.local/`, PNG/PDF/로그, `_codex/`, 도구 캐시 및 개인 파일을 commit하지 않는다. 최종 commit은 현재 feature branch에만 만들고 main에는 merge/push하지 않는다.
