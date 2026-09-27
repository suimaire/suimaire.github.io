# PCR과 프라이머 디자인 Phase 9 Release Candidate QA

검증일: 2026-09-27 (KST). 대상: Core Course 00–07 및 RNA 확장. 판정: **MERGE READY**. 검토한 범위에서 blocking 과학/기능 오류를 발견하지 않았고, 아래 실제 결함 수정 후 실제 Jekyll 출력에 대한 전체 회귀 검사와 학생 E2E를 다시 통과했다. main merge/push 및 production 배포는 수행하지 않았다.

## 1. 기준 commit / branch / remote

| 항목 | QA 시작 시 확인값 |
| --- | --- |
| Branch | `feature/pcr-primer-design` |
| HEAD | `e3ee4ebf0d4dd095e8309a0fa0cdeae83a1c6144` |
| 실제 remote feature HEAD | `e3ee4ebf0d4dd095e8309a0fa0cdeae83a1c6144` |
| origin/main 및 실제 remote main HEAD | `1bca1ad14247c2a1973df3981e3e57c8b0190b2a` |
| Remote | `https://github.com/suimaire/suimaire.github.io.git` |
| Working tree | `?? _codex/`만 존재 |

실제 원격은 `git ls-remote`로 확인했다. 저장소와 상위 `D:/`, `D:/Codex/`에 적용할 AGENTS.md가 없었다. 이 보고서가 자신의 commit hash를 포함하게 만들지 않으며, 최종 RC hash와 feature push 결과는 완료 응답에 기록한다.

## 2. 조사와 보호 범위

Phase 1–8 보고서, 기존 종합 보고서, portal index, 전용 layout, 학습지 전체 HTML/CSS, 16개 PCR JS 모듈, synthetic fixture, state validation/migration, localStorage, import/export, 기존 단위/브라우저 검사를 조사했다. 학생의 조작에 따라 바뀌는 설명도 실제 UI의 단계별 흐름과 렌더링 소스에서 읽었다.

`/verification.local/`은 Git ignore 및 Jekyll exclude 대상이다. 기존 Phase 1–8 검증 자료는 보존하고, 사용자가 이번 작업에 지정한 `verification.local/pcr-primer-design/phase9/`에 새 증거를 생성했다. Phase 7/8 테스트도 출력 root를 선택할 수 있도록 두 줄만 조정하여 새 회귀 결과를 `phase9/regression/`에 격리했다.

`_codex/`는 Git ignore 대상이 아니므로 무차별 staging하지 않았다. 사용자 파일은 수정하거나 commit하지 않았다. Jekyll 출력에도 `_codex/`, `tests/`, `verification.local/`, `node_modules/`가 없음을 확인했다. 이번 실행이 만든 루트 Sass cache는 제거했다. 로컬 test runtime/output은 기존 `tests/.gitignore` 적용 범위다.

## 3. 변경 범위 요약

Phase 9의 production 수정은 CSS 한 규칙, 02 feedback 두 문장, HTML의 다섯 label/문구에 한정된다. 계산, state schema, fixture, RNA 동작, layout 구조 및 portal은 변경하지 않았다.

| 분류 | Phase 9 변경 |
| --- | --- |
| CSS | WebKit에서 방향키로 이동한 native radio에도 기존 focus outline 적용 |
| JS | 숨겨진 legacy 이유 입력란에 답을 요구하던 02 feedback 수정 |
| HTML | 새 기록의 일반 details에서 불필요한 이전 기록 표기 제거, 조사/용어/공식 문서 제목 수정 |
| 기존 tests | Phase 7/8 검증 출력 root 설정 가능 |
| Release tests | 전체 학생 E2E, state snapshot, console, 접근성, 최종 PNG 및 source/failure audit |
| 보고서 | 이 문서 |

새 activity, 계산 엔진, 외부 API, animation, dependency, 디자인 체계 또는 대규모 refactor는 없다. main 대비 feature 전체 파일 목록과 diff stat은 34절에 포함한다.

## 4. 과학 문구와 용어 audit

00의 초기 가설 → 01 기작 → 02 방향 → 03 설계 → 04 특성 검토 → 05 증거 해석 → 06 검색 기록 → 07 판단의 연결을 확인했다. 계산된 예상, 수업용 가상 관찰, 학생이 옮긴 외부 검색 기록, 실제 PCR 수행 상태가 서로 구별된다.

primer, Forward/Reverse primer, amplicon, off-target, Tm, Positive control, NTC, 서열, 데이터베이스, 특이성, genomic DNA (gDNA), cDNA, reverse transcription, exon/intron, exon-exon junction, transcript variant를 문맥에 맞게 유지했다. 자연스러운 전문 영어를 일괄 번역하지 않았다. 최종 갱신 버튼의 `현재 draft`만 기존 UI의 `현재 초안`과 맞췄다. NCBI 문서의 고유 제목에 잘못 섞인 한국어는 원래 영어 제목으로 복원했다.

## 5. PCR 기작 및 초기 product audit

- DNA 합성은 5′→3′이고 primer의 3′ 말단에서 연장한다. 3′ OH와 새 nucleotide 연결 설명이 그림의 화살표와 일치한다.
- 변성은 가닥 분리로 표현하며 DNA backbone 절단으로 설명하지 않는다.
- 95/60/72°C 표시는 교육용 조건 예다. Annealing 온도를 보편적 고정 법칙으로 만들거나 Taq가 72°C에서만 작동한다고 설명하지 않는다.
- 1주기에는 한쪽 끝만 primer로 정의된 긴 가닥이 생길 수 있고, 2주기에 양 끝이 정의된 단일가닥, 3주기에 정확한 길이의 이중가닥이 나타나는 설명을 확인했다. 처음부터 모든 산물이 목표 길이라고 하지 않는다.
- 1–3주기 그림은 분자 수 비율의 정확한 재현이나 항상 100%인 PCR 효율을 주장하지 않는다. 정확한 기존 문구는 유지했다.

## 6. Primer direction / reverse complement audit

Reference의 `5′-AGTCCGTA-3′`에 대한 정렬된 complement `3′-TCAGGCAT-5′`와 주문 방향의 reverse complement `5′-TACGGACT-3′`를 구별한다. Complement 문자열을 그대로 5′→3′ 주문 서열이라고 하지 않는다.

가닥 뒤집기, inward/parallel 배치를 실제로 전환했다. F/R의 위/아래 위치는 절대 규칙이 아니며, reference 방향과 결합 가닥을 기준으로 설명한다. 두 primer의 3′가 안쪽을 향해야 두 primer 위치로 경계가 정의되는 산물을 만들 수 있다는 맥락이 유지된다. 모든 주문 서열 표시는 5′→3′이다.

새 세션에서 보이지 않는 legacy 이유란을 가리키던 feedback의 추가 작성 요구를 제거했다. 정답 판단과 서열 계산은 바꾸지 않았다.

## 7. 활동 03 계산 audit

실제 A/B/C 서열에서 완전 일치 위치를 다시 찾는 구현과 기존 core/design 회귀 검사를 확인했다. B는 A에서 80 bp를 무조건 빼지 않는다. F/R뿐 아니라 가능한 방향 조합과 여러 결합 위치를 검사하는 기존 제한된 교육용 모형이다. 불일치 허용 BLAST나 thermodynamic 모델로 표현하지 않는다.

| 입력/사례 | 확인 결과 |
| --- | --- |
| F 41–60, R 281–300 | F `AGTCGATGCTACGTTGACCA`, R `TTCAGGCTACGATCGTACGA`; 각 20 nt, GC 50%; A 260 bp, B 180 bp, C 예상 산물 없음 |
| F 111–130, R 281–300 | A 190 bp; 결실 경계와 겹친 F의 B 결합 부위 소실 → B 예상 산물 없음 |
| F 121–140, R 281–300 | A 180 bp; B 예상 산물 없음 |
| F 231–250, R 281–300 | 결실 뒤에 있는 pair → A/B 모두 70 bp |
| Invalid direction / overlap / 범위 오류 | 조용히 보정하여 정상 pair처럼 표시하지 않고 이유를 안내 |
| 직접 입력 / multiple hits | 유일한 위치와 다중 위치 구별, R 주문 서열의 reverse complement 및 가능한 산물 확인 |

정상 설계와 결실 경계 대안을 저장한 뒤 정상 설계를 불러오는 흐름도 최종 E2E에서 실행했다. Fixture 정규화 JSON의 SHA-256은 기존 테스트의 `3f404112df85687ccb27a813069129cfb60b6907ea824a102033f6d594490589`와 일치한다.

## 8. 활동 04 primer review audit

18–25 nt와 GC 40–60%는 일반적 출발점이며 합격선이 아니다. 간이 Tm label과 구현은 Wallace `2(A+T)+4(G+C)` 근사로 일치한다. P1은 F/R 각 60°C지만 정밀 Tm, 반응 성공 또는 특이성 판정으로 확장하지 않는다.

3′ GC clamp를 절대 규칙으로 만들지 않는다. Self/F–R complementarity는 gap 없는 antiparallel 상보 run을 살피는 교육용 지표이며, 실제 hairpin 자유에너지, ΔG, dimer 발생 확률로 표시하지 않는다. 서로 다른 Tm을 비교할 때도 실제 조건 검토가 필요하다는 한계가 있다.

P1: A 260/B 180/C 없음, P2: A 260/B 180/C 220, P3: A 180/B 없음/C 없음을 fixture 및 실제 UI와 대조했다. P3의 B 음성을 target 자체의 부재로 단정하지 않는다. 다섯 lens와 A/B → C 공개를 실행했고 저장된 원본 설계를 바꾸지 않는다.

## 9. 활동 05 gel audit

학생 설계의 계산된 예상과 세 가지 고정 수업용 가상 관찰을 구별한다. 각 case에서 Sample, Positive control, NTC lane을 선택하고 관찰/해석/불확실성 답안을 작성했다.

예상 크기의 single band만으로 sequence identity를 확정하지 않는다. NTC의 작은 band도 contamination 또는 primer-dimer 중 하나로 단정하지 않는다. Positive control이 실패한 case 3에서는 sample 음성을 target-negative로 해석할 수 없다고 안내한다. Multiple band의 정체는 gel만으로 정하지 않는다. Band 밝기도 초기 template 양의 직접 정량값이 아니다. 가상 관찰이 학생의 실제 primer 실험 결과라는 표현은 없다.

## 10. 활동 06 Primer-BLAST audit

NCBI 공식 [Primer-BLAST](https://www.ncbi.nlm.nih.gov/tools/primer-blast/), [입력 설명](https://www.ncbi.nlm.nih.gov/tools/primer-blast/primerinfo.html), [공식 사용 안내](https://www.ncbi.nlm.nih.gov/guide/howto/design-pcr-primers/)를 확인했다. 기존 F/R을 모두 제공하면 template 없이 검토할 수 있고, template/accession도 함께 제공하면 intended target 맥락을 지정할 수 있다. Template와 두 primer를 제공한 경우 기존 pair의 specificity 검토라는 설명이 현재 안내와 일치한다. 새 pair 설계에는 target/template가 필요하다. Custom database 등 선택한 검색 범위의 의미를 확인했다.

학습지는 NCBI 검색을 실행하지 않는다. 미실시/진행/기록 완료 상태와 학생이 입력한 조건/결과를 보존한다. Candidate A(첫 후보)를 자동 best로 선택하지 않으며, E2E에서도 명시적 선택 전에는 선택값이 비어 있음을 확인했다. Reported Tm과 내부 Wallace 근사는 별개다. Claim scope에는 기록된 organism/database와 조건 맥락이 있으며, 다른 organism/database/annotation/search setting/wet-lab 조건까지 off-target 부재를 일반화하지 않는다. NCBI verified/validated 자동 판정은 없다.

NCBI 화면의 픽셀 위치를 전제로 하지 않는다. URL과 입력 의미를 안내한다. `특이성 검토을` 오타를 `특이성 검토 설정을`로 수정했다. E2E 데이터의 검색 조건/관찰에는 **실제 검색 미실시 / QA용 모의 기록**을 명시했고 실제 BLAST 요청은 전송하지 않았다.

## 11. 활동 07 final notebook audit

00의 위치는 relative-percent의 대략적 예측으로 보존된다. 03의 bp 좌표를 소급하여 채우지 않는다. 저장한 설계는 저장 당시 sequence, 좌표, 이유를 유지하는 snapshot이다. Current draft를 final로 선택한 뒤 F를 61–80으로 바꾸어도 final snapshot이 바뀌지 않음을 확인하고, 다시 저장 설계 1을 최종 선택했다.

04의 현재 답안은 pair별 실험 검증으로 둔갑하지 않으며, 05는 수업용 가상 관찰, 06은 학생의 외부 기록으로 명시된다. Notebook은 실제 PCR 미수행 상태를 유지하고 자동 validated 결론을 만들지 않는다. Positive/NTC/추가 control, 해석 제한, 미확인 사항, 최종 판단을 작성하고 읽기 전용 notebook에서 확인했다. Editor와 펼친 전체 notebook의 의도된 재표시는 서로 다른 기능이며 자동으로 중복 legacy 편집란을 노출하지 않는다.

## 12. RNA workflow와 세 전략 audit

일반적인 RT-PCR의 RNA → reverse transcription → cDNA → PCR 흐름과 splicing된 mature RNA/cDNA의 exon 연결을 확인했다. 모든 DNA polymerase가 어떤 조건에서도 RNA를 읽을 수 없다는 보편 명제로 서술하지 않는다. 확장 활동은 qPCR 정량 분석 기능으로 확대되지 않았다.

Same-exon pair는 cDNA와 gDNA 양쪽에 결합 부위가 있을 수 있다. Intron-spanning은 intron 길이/extension time/polymerase/반응 조건의 영향을 받으므로 gDNA 증폭 절대 불가능으로 표현하지 않는다. 짧은/긴 intron 예를 모두 전환했다. Exon-exon junction도 절대 cDNA specificity를 보장하지 않는다.

RT(-) 두 case를 확인했다. RT(-) signal은 DNA 유래 가능성 등을 검토할 근거이며 gDNA contamination 확정 진단이 아니다. RT(-) 음성도 모든 오염 가능성을 없애는 증명으로 쓰지 않는다. 모든 RNA 답안/선택/단계가 새로고침과 JSON 복원 후 유지되었다.

## 13. Transcript variant / junction 특별 audit

Exon 2/3/4, junction 2–3/3–4의 다섯 선택을 비교했다. 도식에는 공통 R 부위를 Exon 4로 고정하고 F 후보와 inward pair 구조를 함께 보여 준다. 선택 exon 또는 인접 junction의 존재를 보는 교육 모형이며 실제 증폭/검출, isoform specificity 판정과 구별한다. 실제 pair 두 결합 위치와 전체 amplicon 구조, transcriptome 검토가 필요하다는 기존 설명이 충분하여 중복 경고를 추가하지 않았다.

Junction 설명은 유사 유전자와 processed pseudogene 등 genome 맥락 때문에 gDNA 신호 배제를 보장하지 않는다고 명시한다. Variant 비교에서 같은 junction을 공유하는 transcript가 있을 수 있고 전체 transcriptome에서 실제 특이성을 확인해야 한다는 한계도 확인했다. 새 강의를 추가하거나 RNA 엔진을 바꾸지 않았다.

## 14. Dead visible UI / 이전 기록

빈 localStorage로 시작했을 때 `[data-legacy-answers]`, `#ext-legacy`, `#final-legacy`, `#evidence-legacy`가 보이지 않음을 여섯 E2E에서 확인했다. 03의 일반 details 제목에 항상 붙던 ‘이전 기록 확인’/‘이전 작성 기록’을 제거했다. 실제 오래된 답이 있을 때만 나타나는 기존 wrapper와 migration은 유지했다.

중복 좌표 form, 04의 옛 고정 편집 표, 05의 과거 답안 표, 06의 옛 긴 form, 07 legacy F/R editing, RNA의 과거 단일 답안은 새 기록의 활성 UI로 나타나지 않는다. Phase 1–6 JSON 및 Phase 7 방식의 legacy RNA 기록을 storage와 import로 복원하는 기존 검사를 재실행했다. 원문 답안은 보존되며 사용자 입력을 HTML로 실행하지 않는다. 07의 ‘이전 기록’은 앞선 00–06 학습 기록을 뜻하는 자연스러운 설명이므로 유지했다.

## 15. Dead code / CSS / debug artifact / dependency audit

`pcr-worksheet.mjs`에서 시작하는 import graph로 PCR JS 16개 모두에 도달한다. 미사용 module은 없다. Render/restore/print/event 연결을 읽고 활동 반복 조작을 수행했다. 과거 data migration 분기를 추측으로 제거하지 않았다.

문자열 대조가 지목한 CSS-only 후보는 `.pcr-forward`, `.pcr-row-coordinates`, `.pcr-workbench-subtitle`, `.final-flow`다. 동적 selector/이전 markup까지 고려해야 하므로 이 결과만으로 삭제 안전성을 단정하지 않았다. 충돌하는 학생 UI나 실패를 일으키지 않아 제한된 QA 범위에서 유지했다. Phase별 media/print 규칙은 각 활동의 반응형 및 인쇄 목적을 갖고 있으며, `[hidden]`, reduced-motion, print의 `!important`를 일괄 삭제하지 않았다. 화면별 충돌이나 실패를 발견하지 않아 breakpoint도 유지했다.

Production 소스에서 console.log, debugger, sourceMappingURL, Playwright/test-tool marker 및 노출되는 TEMP/TODO/FIXME 개발 흔적을 점검했다. QA 모의 학생 기록은 테스트와 ignored evidence에만 있으며 배포 asset에 삽입하지 않았다. Built module/CSS/fixture는 소스와 byte 단위로 일치하고 SHA-256을 `source-and-failure-audit.json`에 기록했다. Production은 기존 vanilla ES modules/CSS/JSON이며 새 package/dependency를 추가하지 않았다.

## 16. Production-like Jekyll build

실제 `_config.yml`, `pmarsceill/just-the-docs` remote theme, 전용 `pcr-worksheet` layout, 비어 있는 baseurl, permalink, local asset 경로로 **Jekyll 3.10.0 / Ruby 3.3.12** build를 수행했다. 최종 production 수정 후 재build 결과는 16 pages 성공이다. Remote theme도 실제로 읽었다. Mock HTML을 release 대상으로 삼지 않았다.

실행 환경에 이미 있던 `tests/.pcr-tools/build.rb`는 설치된 gem lib 경로를 적재하여 Jekyll의 `Site#process`를 호출한다. 선택적 서버 native dependency가 필요 없는 실제 정적 빌드이며 Jekyll 소스 자체를 바꾸지 않는다. 출력은 `tests/.pcr-output/site`이다. 기록된 명령:

```powershell
& 'tests/.pcr-tools/ruby/rubyinstaller-3.3.12-1-x64/bin/ruby.exe' --disable-gems tests/.pcr-tools/build.rb
& 'C:/Users/BIO/AppData/Local/Programs/Python/Python312/python.exe' -m http.server 4174 --bind 127.0.0.1 --directory tests/.pcr-output/site
```

Build 로그: `verification.local/pcr-primer-design/phase9/build.log`. Built HTML에 미해결 Liquid가 없고 전용 CSS, 16 modules, fixture가 모두 제공된다. 실제 Pages CI/deploy를 실행한 결과로 표현하지 않는다.

## 17. Portal / public path / direct URL / reload

Portal의 1.3.1 뒤에 **1.3.2 PCR과 프라이머 디자인**이 한 번만 자동 번호로 표시된다. 웹 학습지 설명, 링크, 같은 resource 구조와 시간 제한 표현 부재를 확인했다. 제목 링크 직접 클릭 → 학습지 → `← HAFS Biology Lab 포털` → 재진입이 정상이며 기록은 보존된다. 기존 Day 1–5 링크도 HTTP 200이다.

로컬 실제 build의 `/bioinformatics/pcr-primer-design/`와 `/bioinformatics/pcr-primer-design/index.html`을 확인했다. Activity hash 직접 진입, reload, browser back/forward, transition link와 스크롤에 따른 현재 목차 갱신, 390px 모바일 목차, RNA hash 재진입/펼침을 검사했다. 404나 asset path 실패가 없다. Canonical 공개 경로는 `https://suimaire.github.io/bioinformatics/pcr-primer-design/`이지만, 이 주소가 이번 feature commit을 서비스한다고 확인하지 않았다. 공개 배포 완료를 주장하지 않는다.

## 18. 외부 링크 전체 목록

2026-09-27에 학생 HTML의 9개 고유 HTTPS 링크를 목록화하고 GET 응답을 확인했다. 아래 8개는 HTTP 200이다. NCBI 내용 대조에는 공식 원문을 사용했다.

| 링크 | 확인 |
| --- | --- |
| [NCBI Primer-BLAST](https://www.ncbi.nlm.nih.gov/tools/primer-blast/) | 200 |
| [NCBI 사용 안내](https://www.ncbi.nlm.nih.gov/guide/howto/design-pcr-primers/) | 200 |
| [NCBI 입력 정보](https://www.ncbi.nlm.nih.gov/tools/primer-blast/primerinfo.html) | 200 |
| [Primer3 Manual](https://primer3.org/manual.html) | 200 |
| [CSHL PCR animation](https://dnalc.cshl.edu/resources/animations/pcr.html) | 200 |
| [Addgene primer design](https://www.addgene.org/protocols/primer-design/) | 200 |
| [NC Anchor PCR lab sheet](https://www.ncanchor.org/sites/default/files/documents/pcr-lab-sheet.pdf) | 200 / PDF |
| [ASBMB student workbook](https://www.asbmb.org/getmedia/988b145d-926b-4368-9536-bb09912d15dd/hopes-6-8-student-workbook.pdf) | 200 / PDF |
| [Clinical Chemistry 71(6):634–651](https://academic.oup.com/clinchem/article/71/6/634/8119148) | 자동 GET 403. 공식 검색 결과와 PubMed에서 같은 논문의 현존 확인 |

OUP 직접 자동 접근은 거부되었지만, 출판사 공식 검색 결과와 [PubMed 40272429](https://pubmed.ncbi.nlm.nih.gov/40272429/)에서 같은 DOI `10.1093/clinchem/hvaf043`, 권/페이지/URL을 확인했다. 404나 삭제로 단정하지 않고 공식 링크를 유지했다. 일반 브라우저에서 본문이 열리는지는 남은 확인 한계이며, 9개 모두 자동 HTTP 200이었다고 보고하지 않는다. 기계 판독 가능한 응답 목록은 `external-links.json`이다.

## 19. Release student E2E

`tests/pcr-release-browser.mjs`로 Chromium/WebKit × 1440/768/390px, 각각 새 context의 빈 localStorage에서 다음 흐름을 실행했다. 최종 production 수정 후 여섯 조합 모두 처음부터 재실행하여 PASS했다.

1. 실제 portal의 1.3.2 링크로 진입한다.
2. 00 F/R slider를 keyboard로 15%/80%로 옮기고 초기 이유를 쓴다.
3. 01의 네 단계와 1–3주기, 02의 complement/reverse complement/가닥 뒤집기/inward·parallel을 조작한다.
4. 03 정상 pair를 설계 1로 저장하고, F111–130의 결실 경계 pair를 설계 2로 저장한 뒤 설계 1을 복원한다.
5. 04의 다섯 lens, A/B 비교와 C 공개, 선택 이유/한계를 기록한다.
6. 05의 세 case와 전체 lane을 조작하고 관찰/해석/불확실성 및 reflection을 쓴다.
7. 06 Route A 계획, 명시된 모의 외부 결과, 후보의 직접 선택, claim scope를 기록한다.
8. 07에서 draft 선택 후 draft를 바꾸어 snapshot 불변을 확인하고 설계 1을 final로 선택한다. Controls/reflection/notebook을 완성한다.
9. RNA flow 네 단계, gDNA 비교, 세 전략, 짧은/긴 intron, RT(-) 두 case, 다섯 transcript 후보, 최종 reflection을 기록한다.
10. Reload, JSON export, UI의 reset 확인, 같은 JSON을 UI에서 import, 다시 reload 후 전체 state를 비교한다.
11. Direct hash, back/forward, scroll TOC, portal 왕복, 대표 keyboard 조작, reduced motion 및 반복 조작을 확인한다.

상태를 직접 덮어써서 ‘조작 완료’로 만드는 시나리오가 아니다. 입력/선택/저장/복원은 실제 UI를 사용한다. 완성 기록의 직접 재주입은 최종 스크린샷 생성에만 사용했다. 외부 검색은 실행하지 않고 교육용 fixture와 QA 모의 기록을 사용했다.

## 20. State integrity / export-import roundtrip

초기 예측, immutable saved designs, draft/saved 분리, final snapshot, 04 답안, 05 각 case/lane/답안, 06 route/계획/조건/결과/claim, 07 선택, RNA 전체 state를 비교했다. Export 직전, download한 JSON, import 직후, reload 후의 정규화된 객체가 deep equality로 일치한다.

비교에서 제외한 것은 저장 시각인 `updatedAt`뿐이다. 학습 데이터, 설계 시각, 선택 당시 snapshot을 제외하지 않았다. 여섯 조합의 `*-record.json`과 `*-state-comparison.json`에 before/after/equal을 저장했다. Import 확인 창과 비동기 복원 완료를 기다린 뒤 비교하여 중간 state를 통과로 처리하지 않는다. Reset은 이 학습지의 key만 대상으로 하며, 무관한 localStorage 보호는 기존 검사에서도 확인했다.

## 21. No-JS / load failure / malformed state

Release audit는 두 engine × 일곱 조건의 14개 검사를 통과했다. No-JS와 module 503에서는 자동 저장/계산을 사용할 수 없는 이유와 정적 본문이 남는다. Quota 초과는 설명을 표시하며 현재 입력을 export할 수 있다. Malformed JSON과 unsupported answer field는 기존 raw 기록을 덮어쓰지 않고 현재 입력을 계속 작성/export할 수 있다. Invalid 좌표는 계산란에서 이유를 보여 준다. Partial RNA는 작성한 답을 유지하고 부족한 state에 기본값을 제공한다.

기존 검사에서도 fixture 로딩 실패 시 설명/답안 저장/export, storage 접근 거부, clipboard 거부 시 선택 가능한 전사용 텍스트, 오래되거나 불완전한 기록 복원을 재확인했다. 장애 주입 사례의 의도된 오류는 정상 학생 E2E의 console 집계와 구분한다.

## 22. 최종 테스트 결과와 재현

| 실행 대상 | 최종 결과 |
| --- | --- |
| `node --test tests/pcr-*.test.mjs` | 90 tests PASS, fail/skip/cancel 0 |
| base + intro browser | 548 assertions PASS |
| workbench | 402 PASS |
| review | 720 PASS |
| evidence | 926 PASS |
| external | 512 PASS |
| final | 670 PASS |
| polish | 768 PASS, 24 accessibility audits |
| RNA | 554 PASS, 36 accessibility audits |
| Jekyll integration | portal 번호/direct URL/reload/drag/Day 1–5/fixture failure PASS |
| Release E2E | 6 full journeys PASS, 18 accessibility audits |
| Release source/failure audit | 16 reachable modules, 14 failure checks, build 자산/제외/heading PASS |

기존 browser assertion은 합계 5,100이다. 숫자 증가를 목표로 삼지 않고 신규 검증은 전체 E2E와 release 증거에 집중했다. 최종 전체 회귀 실행도 끝까지 exit 0이었다. 기존 Node 24.20.0, Playwright 1.63.0, axe 설치를 사용했다.

```powershell
$env:PLAYWRIGHT_BROWSERS_PATH = (Resolve-Path 'tests/.pcr-tools/browsers').Path
$env:PCR_TEST_URL = 'http://127.0.0.1:4174/bioinformatics/pcr-primer-design/'
$env:PCR_SITE_URL = 'http://127.0.0.1:4174'
$env:PCR_VERIFICATION_ROOT = 'verification.local/pcr-primer-design/phase9/regression'
& 'D:/nodejs/node.exe' --test tests/pcr-*.test.mjs
& 'D:/nodejs/node.exe' tests/pcr-browser.mjs
& 'D:/nodejs/node.exe' tests/pcr-workbench-browser.mjs
& 'D:/nodejs/node.exe' tests/pcr-review-browser.mjs
& 'D:/nodejs/node.exe' tests/pcr-evidence-browser.mjs
& 'D:/nodejs/node.exe' tests/pcr-external-browser.mjs
& 'D:/nodejs/node.exe' tests/pcr-final-browser.mjs
& 'D:/nodejs/node.exe' tests/pcr-polish-browser.mjs
& 'D:/nodejs/node.exe' tests/pcr-rna-browser.mjs
& 'D:/nodejs/node.exe' tests/pcr-integration.mjs
& 'D:/nodejs/node.exe' tests/pcr-release-browser.mjs
& 'D:/nodejs/node.exe' tests/pcr-release-audit.mjs
```

각 suite 로그와 JSON 증거는 `verification.local/pcr-primer-design/phase9/`에 있다. 기존 Phase를 재실행한 이미지는 그 아래 `regression/`에 저장했다. 로컬 서버는 loopback에만 바인딩했다.

## 23. Desktop visual QA

1440px에서 00 DNA 지도/slider, 01 cycle/초기 product, 02 방향, 03 workbench, 04 complementarity/A·B·C, 05 gel, 06 plan/results, 07 비교/notebook, RNA 전략/transcript를 실제 이미지로 확인했다. 제목, 입력란, 도식 label, 테두리, 빈 panel, 의도하지 않은 잘림/겹침을 점검했다.

Fullpage는 1440×35511px이다. 전체 section 간격/반복과 마지막 RNA 위치를 확인했다. Notebook을 펼친 완성 기록이어서 07이 길지만, 편집 시 접힌 detail과 의도적인 전체 기록 보기를 구별한다. Notebook의 글자는 전체 이미지 축소만으로 판단하지 않고 pair/외부 기록/최종 판단의 확대 section 이미지에서도 읽었다. 대표 section PNG의 높이 상한에 따른 하단 clip은 스크린샷의 부분 캡처이며 페이지 내용의 유실이 아니다.

## 24. Tablet visual QA

768px에서 03의 세로 배치 workbench, 04 candidate comparison, 05 네 lane, 06 forms/후보, 07 초기 예측과 최종 pair, RNA exon/transcript를 확인했다. Form 열 너비, 줄바꿈, 선택 조작이 유지된다. 기존 breakpoint 수정이 필요한 과밀 배치나 페이지 가로 overflow는 발견하지 않았다.

## 25. Mobile visual QA

390px에서 전체 E2E를 실행하고 단계별 document scrollWidth가 viewport를 넘지 않는지 확인했다. Sequence는 제한된 창과 이동 조작으로 읽을 수 있고, gel lane 선택, external form, candidate comparison, final notebook의 주문 서열/본문도 사용할 수 있다. 글씨를 극단적으로 줄이는 수정을 하지 않았다.

RNA 그림에는 label을 뭉개지 않기 위한 의도적인 내부 가로 스크롤 영역이 있다. 페이지 전체 가로 overflow와 구별하며, focus와 방향키로 내부 이동하는 동작은 기존 RNA 검사로 확인했다. Exon block과 transcript label의 끝도 실제 이미지로 점검했다.

## 26. Accessibility / heading / keyboard

기존 접근성 검사를 포함한 전체 suite가 PASS했다. Release에서는 fresh/core-notebook/rna-transcripts의 세 상태 × 여섯 조합, 총 18회의 axe (WCAG 2 A/AA, 2.1 AA, 2.2 AA, best-practice)를 실행하여 violations 0이었다. Phase 7의 24회, Phase 8의 36회도 위반 0이다.

**자동 판정의 한계는 남긴다.** Release의 각 audit에 color-contrast incomplete가 있으며, 대부분 SVG text와 결실 경계 label로 37–46 nodes/audit이다. Violation 0을 모든 색의 자동 보장으로 해석하지 않는다. 그림을 실제 이미지로 읽고, 최종 기록에서 SVG text의 computed fill이 본문색/강조색인 점도 확인했다. 이를 screen reader 검사나 모든 상태의 수동 WCAG 적합성 인증이라고 부르지 않는다.

h1은 하나이고 Core는 h2, 그 내부는 h3 이후이며 RNA도 h2/h3을 가진다. Audit에서 heading level의 건너뛰기를 검사했고 notebook의 중첩 제목도 유지했다. `lang=ko`, main, nav 이름, skip link, input label, status/live 영역, 그림의 텍스트 대안을 확인했다.

실제 keyboard event로 slider, details, cycle, sequence, review tab, gel lane, route, final radio, RNA 전략/transcript를 조작했다. Tab/Space, 방향키, Home/Enter를 사용했다. 긴 E2E의 다수 조작은 대상에 focus한 뒤 키를 누르는 방식이며, 처음부터 끝까지 Tab만으로 이동했다고 주장하지 않는다. 보충으로 Chromium 1440/390px에서 처음부터 각각 38회 연속 Tab을 실행해 skip link→portal→기록/인쇄/목차→활동의 논리적 순서와 모든 대상의 outline을 확인했다 (`visual-review-details.json`).

WebKit의 final radio에서 ArrowRight 뒤 `:focus-visible`이 false가 되며 outline이 사라지는 실제 버그를 재현했다. Native radio의 `:focus`에도 기존 3px outline을 적용해 수정하고 두 engine/모든 너비에서 재확인했다. 실제 screen reader 발화는 미검증이다.

## 27. Reduced motion / performance sanity / console

Reduced-motion 환경에서 기존 CSS가 transition/animation/scroll behavior를 억제하는지 확인하고, PCR 단계를 조작하여 긴 animation이 남지 않는지 검사했다. 새 motion은 추가하지 않았다.

PCR 단계를 12회 왕복해 DOM node 수가 늘지 않았다. 복원한 기록 상태의 node 수는 desktop/tablet 2316, mobile 2224로 반복 전후 동일하다. 반복 scroll/활동 이동/재진입에서 handler 중복 실행 증상은 없고 sequence window도 제한된다. 이는 성능 sanity check이며 heap profiler로 memory leak 부재를 증명하거나 대규모 benchmark를 수행한 것은 아니다.

여섯 E2E 모두 console.error, pageerror, unhandledrejection, requestfailed, HTTP 4xx/5xx는 **0**이다. 외부 링크로 이동하지 않는 정상 흐름의 결과다. Console.warning은 Chromium 각 journey 40개, WebKit 0개다. 40개는 reload/왕복 시 반복된 아래 기존 font 경고뿐이었다.

```text
Failed to decode downloaded font: .../assets/fonts/NanumSquareNeoR.ttf (또는 B.ttf)
OTS parsing error: file less than 4 bytes
```

두 font는 main/시작 HEAD 모두 동일 blob `d3f5a12faa99758192ecc4ed3fc22c9249232e86`, 각각 2바이트의 기존 공용 파일이다. HTTP는 200이지만 font로 유효하지 않아 CSS의 Malgun Gothic/system-ui fallback이 적용된다. 이번 실제 fallback 이미지의 가독성과 조작/저장을 확인했다. 공용 font 교체는 다른 페이지까지 바꾸므로 이번 QA에서 변경하지 않고 비차단 기존 제약으로 명시한다. Warning 0이라고 보고하지 않는다.

## 28. 실제 수정한 문제 — 중요도 순

1. **Keyboard focus 누락**: WebKit에서 final 설계 radio를 방향키로 옮기면 outline 소실. CSS native radio focus fallback으로 수정하고 여섯 E2E 및 전체 회귀를 재실행했다.
2. **없는 입력란으로 안내**: 02 feedback이 새 세션에서는 숨겨진 legacy 이유란에 추가 작성을 요구했다. 잘못된 작성 요구를 없애고 3′ 말단의 합성 설명을 유지했다.
3. **불필요한 이전 기록 제목 노출**: 03 일반 details에 항상 붙던 이전 기록 안내를 일반 제목으로 바꾸고 실제 과거 기록 wrapper는 유지했다.
4. **표기 불일치**: `특이성 검토을` 문법, `현재 draft` 용어, NCBI 공식 제목의 혼합 표기를 수정했다.

계산/데이터 손상이나 blocking 과학 오류는 발견하지 않아 계산 엔진과 정확한 기존 과학 문구를 바꾸지 않았다.

## 29. 남은 수동 검증 한계

- Playwright WebKit은 실제 iPhone/iPad/macOS Safari 수동 확인과 다르다.
- 실기기 touch/OS 입력 지원/screen reader 발화는 미검증이다. Axe incomplete의 SVG 색 판정에는 자동 검사 한계가 있다.
- Wet-lab PCR, 실제 amplicon sequence, 실제 NCBI 검색 결과의 진위를 검증하지 않았다. 외부 결과는 명시된 QA 모의 기록이다.
- 공개 GitHub Pages에 이번 commit을 배포하거나 실제 서비스 여부를 확인하지 않았다. 검사 대상은 실제 로컬 Jekyll 출력이다.
- OUP 논문 직접 자동 GET은 403이며 일반 브라우저에서 본문이 열리는지는 미확인이다. 공식 metadata의 현존은 확인했다.
- 기존 2바이트 공용 font의 Chromium 경고를 유지했다. Fallback으로 검증했으며 NanumSquareNeo 실제 font가 표시된다고 주장하지 않는다.
- Remote theme은 저장소의 기존 미고정 참조다. 미래 theme 업데이트까지 포함한 무기한 재현성을 보장하지 않는다.

## 30. Merge readiness 판단

**MERGE READY**. 이번 검사 범위에서 blocking bug/과학적 오류가 없으며, 최종 production-like build, portal integration, 여섯 전체 E2E, state roundtrip, final snapshot, responsive, 기존 전체 회귀, 접근성 violation 0, keyboard 수정 후 확인, console blocking error 0, diff 범위 점검, release 증거 생성을 충족한다.

위 수동 한계와 비차단 기존 font 경고를 명시한 판단이다. 기술적 QA 상태이며 자동 merge/deploy 승인이 아니다. 현재 feature만 commit/push 대상으로 하고, main checkout/merge/push 및 production deploy는 수행하지 않는다.

## 31. Release PNG 절대 경로

다음은 여섯 journey 통과 후 같은 최종 build에서 생성했다. 대표 PNG 30개와 contact sheet, 추가 notebook 판독용 detail PNG 6개를 저장했다. QA 입력은 공개 source에 포함하지 않는다.

| PNG | 절대 경로 |
| --- | --- |
| release-top.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-top.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-top.png>) |
| release-00.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-00.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-00.png>) |
| release-01-cycle.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-01-cycle.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-01-cycle.png>) |
| release-02-direction.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-02-direction.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-02-direction.png>) |
| release-03-workbench.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-03-workbench.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-03-workbench.png>) |
| release-04-review.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-04-review.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-04-review.png>) |
| release-05-gel.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-05-gel.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-05-gel.png>) |
| release-06-primer-blast.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-06-primer-blast.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-06-primer-blast.png>) |
| release-06-result.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-06-result.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-06-result.png>) |
| release-07-summary.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-07-summary.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-07-summary.png>) |
| release-07-final-notebook.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-07-final-notebook.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-07-final-notebook.png>) |
| release-rna-extension.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-rna-extension.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-rna-extension.png>) |
| release-rna-transcripts.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-rna-transcripts.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-rna-transcripts.png>) |
| course-fullpage-desktop.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/course-fullpage-desktop.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/course-fullpage-desktop.png>) |
| release-768-03.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-768-03.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-768-03.png>) |
| release-768-04.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-768-04.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-768-04.png>) |
| release-768-05.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-768-05.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-768-05.png>) |
| release-768-06.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-768-06.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-768-06.png>) |
| release-768-07.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-768-07.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-768-07.png>) |
| release-tablet-rna.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-tablet-rna.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-tablet-rna.png>) |
| release-768-rna-transcripts.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-768-rna-transcripts.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-768-rna-transcripts.png>) |
| release-390-03.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-390-03.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-390-03.png>) |
| release-390-04.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-390-04.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-390-04.png>) |
| release-390-05.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-390-05.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-390-05.png>) |
| release-390-06.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-390-06.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-390-06.png>) |
| release-390-07.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-390-07.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-390-07.png>) |
| release-mobile-core.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-mobile-core.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-mobile-core.png>) |
| release-mobile-notebook.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-mobile-notebook.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-mobile-notebook.png>) |
| release-mobile-rna.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-mobile-rna.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-mobile-rna.png>) |
| release-390-rna-transcripts.png | [D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-390-rna-transcripts.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-390-rna-transcripts.png>) |

## 32. Release contact sheet / fullpage

Contact sheet는 3열로 top, 03, 04, 05, 06, 07, RNA 전략, RNA transcript, mobile을 배치했다. 실제 공개 후보의 대표 상태를 보여 준다.

- [release-contact-sheet.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-contact-sheet.png>)
- [course-fullpage-desktop.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/course-fullpage-desktop.png>)
- [release-07-final-notebook.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-07-final-notebook.png>)
- [release-rna-extension.png](<D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/release-rna-extension.png>)

## 33. Machine-readable 증거

Evidence root: `D:/Codex/260926 bioinformatics/verification.local/pcr-primer-design/phase9/`.

| 파일 | 내용 |
| --- | --- |
| `release-verification.json` | 여섯 engine/width 결과, console/network 목록, 18 axe 결과/incomplete, PNG 목록 |
| `chromium-{1440,768,390}-record.json` / `webkit-{1440,768,390}-record.json` | 실제 UI에서 download한 완성 기록 |
| 같은 prefix의 `-state-comparison.json` | before/after/equal |
| `source-and-failure-audit.json` | module graph, CSS 후보, asset hash, build 목록, heading, 14 failure 결과 |
| `external-links.json` | 전체 외부 URL의 HTTP 응답 |
| `visual-review-details.json` | 실제 연속 Tab 순서, focus outline, SVG text 색 |
| `build.log`, `unit.log`, `*-browser.log`, `integration.log`, `release.log`, `audit.log` | 최종 실행 로그 |
| `regression/phase2`–`phase8` | 기존 suite의 이번 재실행 증거 |

## 34. origin/main 대비 final diff / commit 제외

시작 시 main 대비 62개 파일이었으며 Phase 9 report와 두 release script를 더해 최종 65개다. 전체 목록:

| 분류 | 파일 수 |
| --- | --- |
| Portal integration | 1 |
| Worksheet source / fixture | 3 |
| CSS | 1 |
| JS | 16 |
| Tests | 30 |
| Reports | 11 |
| Config | 3 |

### Portal integration

```text
index.md
```

### Worksheet source / fixture

```text
_layouts/pcr-worksheet.html
assets/data/pcr-primer-fixture.json
bioinformatics/pcr-primer-design.html
```

### CSS

```text
assets/css/pcr-worksheet.css
```

### JS

```text
assets/js/pcr-core.mjs
assets/js/pcr-design.mjs
assets/js/pcr-evidence-view.mjs
assets/js/pcr-evidence.mjs
assets/js/pcr-external-view.mjs
assets/js/pcr-external.mjs
assets/js/pcr-final-view.mjs
assets/js/pcr-final.mjs
assets/js/pcr-intro.mjs
assets/js/pcr-records.mjs
assets/js/pcr-review-view.mjs
assets/js/pcr-review.mjs
assets/js/pcr-rna-view.mjs
assets/js/pcr-rna.mjs
assets/js/pcr-workbench.mjs
assets/js/pcr-worksheet.mjs
```

### Tests

```text
tests/fixtures/pcr/phase1.json
tests/fixtures/pcr/phase2.json
tests/fixtures/pcr/phase3.json
tests/fixtures/pcr/phase4.json
tests/fixtures/pcr/phase5.json
tests/fixtures/pcr/phase6.json
tests/pcr-browser.mjs
tests/pcr-build.rb
tests/pcr-core.test.mjs
tests/pcr-design.test.mjs
tests/pcr-evidence-browser.mjs
tests/pcr-evidence.test.mjs
tests/pcr-external-browser.mjs
tests/pcr-external.test.mjs
tests/pcr-final-browser.mjs
tests/pcr-final.test.mjs
tests/pcr-integration.mjs
tests/pcr-intro-browser.mjs
tests/pcr-polish-browser.mjs
tests/pcr-polish.test.mjs
tests/pcr-preview.mjs
tests/pcr-records.test.mjs
tests/pcr-release-audit.mjs
tests/pcr-release-browser.mjs
tests/pcr-review-browser.mjs
tests/pcr-review.test.mjs
tests/pcr-rna-browser.mjs
tests/pcr-rna.test.mjs
tests/pcr-static.test.mjs
tests/pcr-workbench-browser.mjs
```

### Reports

```text
tests/fixtures/pcr/README.md
tests/PCR_PRIMER_DESIGN_PHASE1.md
tests/PCR_PRIMER_DESIGN_PHASE2.md
tests/PCR_PRIMER_DESIGN_PHASE3.md
tests/PCR_PRIMER_DESIGN_PHASE4.md
tests/PCR_PRIMER_DESIGN_PHASE5.md
tests/PCR_PRIMER_DESIGN_PHASE6.md
tests/PCR_PRIMER_DESIGN_PHASE7.md
tests/PCR_PRIMER_DESIGN_PHASE8.md
tests/PCR_PRIMER_DESIGN_PHASE9_RELEASE_QA.md
tests/PCR_PRIMER_DESIGN.md
```

### Config

```text
_config.yml
.gitignore
tests/.gitignore
```

<!-- DIFF_STAT_START -->
```text
 .gitignore                                   |   1 +
 _config.yml                                  |   1 +
 _layouts/pcr-worksheet.html                  |  16 +
 assets/css/pcr-worksheet.css                 | 635 +++++++++++++++++++
 assets/data/pcr-primer-fixture.json          |  98 +++
 assets/js/pcr-core.mjs                       |  80 +++
 assets/js/pcr-design.mjs                     |  80 +++
 assets/js/pcr-evidence-view.mjs              | 118 ++++
 assets/js/pcr-evidence.mjs                   |  46 ++
 assets/js/pcr-external-view.mjs              | 144 +++++
 assets/js/pcr-external.mjs                   |  78 +++
 assets/js/pcr-final-view.mjs                 | 178 ++++++
 assets/js/pcr-final.mjs                      |  97 +++
 assets/js/pcr-intro.mjs                      | 224 +++++++
 assets/js/pcr-records.mjs                    |  85 +++
 assets/js/pcr-review-view.mjs                | 183 ++++++
 assets/js/pcr-review.mjs                     |  63 ++
 assets/js/pcr-rna-view.mjs                   | 138 ++++
 assets/js/pcr-rna.mjs                        |  64 ++
 assets/js/pcr-workbench.mjs                  | 276 ++++++++
 assets/js/pcr-worksheet.mjs                  | 246 ++++++++
 bioinformatics/pcr-primer-design.html        | 911 +++++++++++++++++++++++++++
 index.md                                     |  12 +
 tests/.gitignore                             |   2 +
 tests/PCR_PRIMER_DESIGN.md                   | 116 ++++
 tests/PCR_PRIMER_DESIGN_PHASE1.md            | 113 ++++
 tests/PCR_PRIMER_DESIGN_PHASE2.md            | 148 +++++
 tests/PCR_PRIMER_DESIGN_PHASE3.md            | 190 ++++++
 tests/PCR_PRIMER_DESIGN_PHASE4.md            | 206 ++++++
 tests/PCR_PRIMER_DESIGN_PHASE5.md            | 239 +++++++
 tests/PCR_PRIMER_DESIGN_PHASE6.md            | 217 +++++++
 tests/PCR_PRIMER_DESIGN_PHASE7.md            | 373 +++++++++++
 tests/PCR_PRIMER_DESIGN_PHASE8.md            | 229 +++++++
 tests/PCR_PRIMER_DESIGN_PHASE9_RELEASE_QA.md | 555 ++++++++++++++++
 tests/fixtures/pcr/README.md                 |  14 +
 tests/fixtures/pcr/phase1.json               |  42 ++
 tests/fixtures/pcr/phase2.json               | 111 ++++
 tests/fixtures/pcr/phase3.json               |  79 +++
 tests/fixtures/pcr/phase4.json               |  95 +++
 tests/fixtures/pcr/phase5.json               | 145 +++++
 tests/fixtures/pcr/phase6.json               | 161 +++++
 tests/pcr-browser.mjs                        | 159 +++++
 tests/pcr-build.rb                           |  15 +
 tests/pcr-core.test.mjs                      |  77 +++
 tests/pcr-design.test.mjs                    |  87 +++
 tests/pcr-evidence-browser.mjs               | 197 ++++++
 tests/pcr-evidence.test.mjs                  |  62 ++
 tests/pcr-external-browser.mjs               | 182 ++++++
 tests/pcr-external.test.mjs                  |  66 ++
 tests/pcr-final-browser.mjs                  | 186 ++++++
 tests/pcr-final.test.mjs                     | 103 +++
 tests/pcr-integration.mjs                    |  43 ++
 tests/pcr-intro-browser.mjs                  |  81 +++
 tests/pcr-polish-browser.mjs                 | 151 +++++
 tests/pcr-polish.test.mjs                    |  20 +
 tests/pcr-preview.mjs                        |  30 +
 tests/pcr-records.test.mjs                   |  36 ++
 tests/pcr-release-audit.mjs                  |  80 +++
 tests/pcr-release-browser.mjs                | 298 +++++++++
 tests/pcr-review-browser.mjs                 | 217 +++++++
 tests/pcr-review.test.mjs                    | 105 +++
 tests/pcr-rna-browser.mjs                    | 200 ++++++
 tests/pcr-rna.test.mjs                       |  69 ++
 tests/pcr-static.test.mjs                    |  25 +
 tests/pcr-workbench-browser.mjs              | 180 ++++++
 65 files changed, 9478 insertions(+)
```
<!-- DIFF_STAT_END -->

`git diff --check`와 staged 차이를 점검했다. Portal/worksheet/CSS/JS/fixture/tests/reports/config 이외의 무관한 변경은 없다. PNG, `verification.local/`, `_codex/`, runtime/cache, 개인 파일은 commit에 포함하지 않는다. 최종 commit 뒤에도 `origin/main...HEAD`의 stat/file 목록을 재대조한다.
