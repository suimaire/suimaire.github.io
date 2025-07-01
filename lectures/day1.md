---
layout: default
title: "Day 1"
nav_order: 2
---

# Day 1: 기본 개념 소개 및 환경 구축


##  1부: 생물학 기본 개념과 Python 소개 (50분)

###  (15분) 생물학 기본 용어와 개념 소개
  ---

## 1. 분자생물학 & 중심원리

| 용어 | 한국어 해석 | 설명 |
|-----------|-------------|-----------|
| gene | 유전자 | **DNA 기능 단위** |
| genome | 유전체 | 한 생물의 모든 유전정보 |
| chromosome | 염색체 | DNA와 히스톤 단백질이 응축된 형태  |
| DNA / RNA | 디엔에이 / 알엔에이 | 핵산의 두 종류 |
| replication | 복제 | **DNA → DNA** |
| transcription | 전사 | **DNA → mRNA** |
| translation | 번역 | **mRNA → 단백질** |
| codon | 코돈 | mRNA 3 염기 묶음 |
| start codon | 개시코돈 | 보통 **AUG** |
| stop codon | 종결코돈 | **UAA / UAG / UGA** |
| reading frame | 리딩 프레임 | 3염기씩 자르는 기준 |
| ORF (open reading frame) | 개방된 리딩 프레임 | 시작 ∼ 종결코돈 사이 |
| mutation | 돌연변이 | 염기서열 변화 |
| point mutation | 점 돌연변이 | 염기 1개 변화 |
| insertion / deletion | 삽입 / 결실 | 염기 추가 · 삭제 |
| frameshift | 프레임 이동 | (변이 수≠3의 배수) 에 의한 프레임 틀어짐 |
| missense / nonsense / silent | 미스센스 / 넌센스 / 침묵 | 아미노산 변경 · 조기 종결 · 무변화 |

  ---

## 2. NCBI & 변이 데이터베이스 (Natural Center for Biotechnology Information)

| 영어 용어 | 한국어 대응 | 핵심 의미 |
|-----------|-------------|-----------|
| accession number | 등록번호 | 서열 · 논문 고유 ID |
| FASTA / GenBank / RefSeq | 포맷(형식) · DB 이름 | 서열 · 주석 저장 방식 |
| locus | 로커스 | 유전자의 물리적 위치 |
| annotation | 주석 | 기능 · 구조 정보 |
| BLAST | 블라스트 | 서열 유사성 탐색 툴(구글 이미지 검색의 생물정보학 버전) |
| variant / variation | 변이 | 표준서열 대비 차이 |
| consequence (type) | 결과 유형 | 변이가 끼치는 분자적 영향 |
| condition | 임상 조건 | 질환·표현형 명칭 |
| classification | 분류 | Pathogenic / Benign 등 |
| ClinVar / dbSNP | 클린바 / dbSNP | 임상(의학판 위키백과) / 단일염기 변이 DB(어디에서 A→G가 흔한가?) |
| allele / genotype / phenotype | 대립유전자 / 유전자형 / 표현형 | 하디-바인베르크 연결 |

  ---

## 3. Google Colab & Python 셸

| 영어 용어 | 한국어 대응 | 용도 |
|-----------|-------------|------|
| `!` (bang) | 셸 명령 접두사 | `!ls`, `!pip install` |
| `pip` | 패키지 관리자 | 라이브러리 설치 |
| `cd` / `ls` / `pwd` | 디렉터리 이동 · 목록 확인 · 경로 확인 |
| install / import | 설치 / 가져오기 | 패키지 준비 → 사용 |
| notebook / cell | 노트북 / 셀 | 코드·텍스트 실행 단위 |
| interpreter | 인터프리터 | 코드 실행기 |
| shell / terminal | 셸 / 터미널 | 명령줄 창 |
| variable / function / module | 변수 / 함수 / 모듈 | Python 기본 단위 |
| list / dict | 리스트 / 딕셔너리 | 주요 자료형 |
| `for` / `if` | 반복 / 조건 | 제어 구문 |
| indentation | 들여쓰기 | 블록 구분 규칙 |
| stdout / stderr | 표준출력 / 표준오류 | 결과·오류 창구 |

---

## 4. Biopython 핵심 객체 (2차시 예고)

| 모듈·함수 | 한국어 대응 | 간단 설명 |
|-----------|-------------|-----------|
| `Bio` | 최상위 패키지 | 서브모듈 루트 |
| `Seq`, `SeqRecord` | 서열 / 서열 기록 | 서열 + 메타데이터 저장 |
| `SeqIO` | 서열 입출력 | FASTA·GenBank 파싱 |
| `translate()` | 번역 함수 | DNA/RNA → 단백질 |
| `Entrez` | NCBI API | 원격 서열 다운로드 |
| `AlignIO` | 정렬 입출력 | 서열 정렬 파일 파싱 |

  ---

###  (10분) 영문 타자 연습과 영단어 소개
- https://tadaktadak.co.kr/taja/en_sentence.html (타닥타닥)
- https://support.hancomtaja.com/2fa61bff-472b-4a91-bbe7-f60f1d3f78ab　(한컴 타자)
  
###  (15분) 단백질 기초 개념과 아미노산 약자 설명
### 표준 20 아미노산 — 구조식 · 약어 · 코돈 (로컬 이미지 버전)

| 구조식<sup>†</sup> | Amino acid | 3-letter | 1-letter | Codon(s) (RNA) |
|:---:|:---|:---:|:---:|:---|
| <img src="/assets/img/amino/alanine.png"          width="110"/> | Alanine 알라닌       | Ala | **A** | GCU · GCC · GCA · GCG |
| <img src="/assets/img/amino/arginine.png"         width="110"/> | Arginine 아르지닌       | Arg | **R** | CGU · CGC · CGA · CGG · AGA · AGG |
| <img src="/assets/img/amino/asparagine.png"       width="110"/> | Asparagine 아스파라진    | Asn | **N** | AAU · AAC |
| <img src="/assets/img/amino/asparatate.png"    width="110"/> | Aspartic acid 아스파르트산 | Asp | **D** | GAU · GAC |
| <img src="/assets/img/amino/cysteine.png"         width="110"/> | Cysteine 시스테인      | Cys | **C** | UGU · UGC |
| <img src="/assets/img/amino/glutamate.png"    width="110"/> | Glutamic acid 글루탐산 | Glu | **E** | GAA · GAG |
| <img src="/assets/img/amino/glutamine.png"        width="110"/> | Glutamine 글루타민      | Gln | **Q** | CAA · CAG |
| <img src="/assets/img/amino/glycine.png"          width="110"/> | Glycine 글라이신(글리신)        | Gly | **G** | GGU · GGC · GGA · GGG |
| <img src="/assets/img/amino/histidine.png"        width="110"/> | Histidine 히스티딘     | His | **H** | CAU · CAC |
| <img src="/assets/img/amino/isoleucine.png"       width="110"/> | Isoleucine 아이소류신(이소류신)    | Ile | **I** | AUU · AUC · AUA |
| <img src="/assets/img/amino/leucine.png"          width="110"/> | Leucine 류신       | Leu | **L** | UUA · UUG · CUU · CUC · CUA · CUG |
| <img src="/assets/img/amino/lysine.png"           width="110"/> | Lysine 라이신        | Lys | **K** | AAA · AAG |
| <img src="/assets/img/amino/methionine.png"       width="110"/> | Methionine 메싸이오닌(메티오닌)    | Met | **M** | AUG |
| <img src="/assets/img/amino/phenylalanine.png"    width="110"/> | Phenylalanine 페닐알라닌 | Phe | **F** | UUU · UUC |
| <img src="/assets/img/amino/proline.png"          width="110"/> | Proline 프롤린       | Pro | **P** | CCU · CCC · CCA · CCG |
| <img src="/assets/img/amino/serine.png"           width="110"/> | Serine 세린        | Ser | **S** | UCU · UCC · UCA · UCG · AGU · AGC |
| <img src="/assets/img/amino/threonine.png"        width="110"/> | Threonine 트레오닌     | Thr | **T** | ACU · ACC · ACA · ACG |
| <img src="/assets/img/amino/tryptophane.png"       width="110"/> | Tryptophane 트립토판    | Trp | **W** | UGG |
| <img src="/assets/img/amino/tyrosine.png"         width="110"/> | Tyrosine 타이로신(티로신)      | Tyr | **Y** | UAU · UAC |
| <img src="/assets/img/amino/valine.png"           width="110"/> | Valine 발린        | Val | **V** | GUU · GUC · GUA · GUG |



###  (10분) 전사(transcription), 번역(translation) 설명
- 유튜브 영상 활용

---

##  2부: Python과 Google Colab 소개 (50분)

###  (15분) Google Colab 환경 구축
- [colab.research.google.com](https://colab.research.google.com)
- 간단한 코드 실행: `print("Hello, Bioinformatics!")`

###  (15분) Python 기본 자료형 소개
- 변수, 문자열, 숫자 자료형 및 연산 예시

###  (15분) Python 코드 실습 (아미노산 출력 예시)

```python
amino_acids = ["Ala", "Gly", "Val", "Leu"]
for aa in amino_acids:
    print(aa)
```
- 최적화 개요  
- Gradient Descent  
- Learning Rate 스케줄링  

**실습 자료**  
- [Colab 노트북](lectures/day1_notebook.ipynb)  
