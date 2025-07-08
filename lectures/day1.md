---
layout: default
title: "Day 1"
nav_order: 2
---

# Day 1 기본 개념 소개 및 환경 구축
1부: 생물학 기본 개념과 Python 소개 (50분)

  ---
## 영문 타자 연습과 영단어 소개

-   https://support.hancomtaja.com/2fa61bff-472b-4a91-bbe7-f60f1d3f78ab   (한컴 타자)
-   https://tadaktadak.co.kr/taja/en_sentence.html   (타닥타닥)
　 
  ---

## 중심 원리 및 아미노산(Amino acids) 소개

### 전사와 번역

- [from DNA to protein - 3D](https://www.youtube.com/watch?v=gG7uCskUOrA)

### 20 종류의 아미노산 : 구조식 · 약어 · 코돈

| 구조식<sup>†</sup> | Amino acid | 3-letter | 1-letter | Codon(s) (RNA) |
|:---:|:---|:---:|:---:|:---|
| <img src="/assets/img/amino/alanine.png"          width="110"/> | Alanine<br>알라닌                        | Ala | **A** | GCU · GCC · GCA · GCG |
| <img src="/assets/img/amino/arginine.png"         width="110"/> | Arginine<br>아르지닌                     | Arg | **R** | CGU · CGC · CGA · CGG · <br>AGA · AGG |
| <img src="/assets/img/amino/asparagine.png"       width="110"/> | Asparagine<br>아스파라진                 | Asn | **N** | AAU · AAC |
| <img src="/assets/img/amino/asparatate.png"       width="110"/> | Aspartic acid<br>아스파르트산            | Asp | **D** | GAU · GAC |
| <img src="/assets/img/amino/cysteine.png"         width="110"/> | Cysteine<br>시스테인                     | Cys | **C** | UGU · UGC |
| <img src="/assets/img/amino/glutamate.png"        width="110"/> | Glutamic acid<br>글루탐산                | Glu | **E** | GAA · GAG |
| <img src="/assets/img/amino/glutamine.png"        width="110"/> | Glutamine<br>글루타민                    | Gln | **Q** | CAA · CAG |
| <img src="/assets/img/amino/glycine.png"          width="110"/> | Glycine<br>글라이신<br>(글리신)           | Gly | **G** | GGU · GGC · GGA · GGG |
| <img src="/assets/img/amino/histidine.png"        width="110"/> | Histidine<br>히스티딘                    | His | **H** | CAU · CAC |
| <img src="/assets/img/amino/isoleucine.png"       width="110"/> | Isoleucine<br>아이소류신<br>(이소류신)    | Ile | **I** | AUU · AUC · AUA |
| <img src="/assets/img/amino/leucine.png"          width="110"/> | Leucine<br>류신                          | Leu | **L** | UUA · UUG · <br>CUU · CUC · CUA · CUG |
| <img src="/assets/img/amino/lysine.png"           width="110"/> | Lysine<br>라이신                         | Lys | **K** | AAA · AAG |
| <img src="/assets/img/amino/methionine.png"       width="110"/> | Methionine<br>메싸이오닌<br>(메티오닌)    | Met | **M** | AUG |
| <img src="/assets/img/amino/phenylalanine.png"    width="110"/> | Phenylalanine<br>페닐알라닌              | Phe | **F** | UUU · UUC |
| <img src="/assets/img/amino/proline.png"          width="110"/> | Proline<br>프롤린                        | Pro | **P** | CCU · CCC · CCA · CCG |
| <img src="/assets/img/amino/serine.png"           width="110"/> | Serine<br>세린                           | Ser | **S** | UCU · UCC · UCA · UCG · <br>AGU · AGC |
| <img src="/assets/img/amino/threonine.png"        width="110"/> | Threonine<br>트레오닌                    | Thr | **T** | ACU · ACC · ACA · ACG |
| <img src="/assets/img/amino/tryptophane.png"      width="110"/> | Tryptophane<br>트립토판                  | Trp | **W** | UGG |
| <img src="/assets/img/amino/tyrosine.png"         width="110"/> | Tyrosine<br>타이로신<br>(티로신)          | Tyr | **Y** | UAU · UAC |
| <img src="/assets/img/amino/valine.png"           width="110"/> | Valine<br>발린                           | Val | **V** | GUU · GUC · GUA · GUG |



 ---

## 생물학 기본 용어 및 Bioinformatics Camp를 위한 용어 정리
### 분자생물학 & 중심원리

| 용어 | 한국어 해석 | 설명 · 용도 |
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

### NCBI & 변이 데이터베이스 (Natural Center for Biotechnology Information)

| 영어 용어 | 한국어 대응 | 핵심 의미 |
|-----------|-------------|-----------|
| accession number | 등록번호 | 서열 · 논문 고유 ID |
| FASTA / GenBank / RefSeq | 포맷(형식) · DB 이름 | 서열 · 주석 저장 방식 |
| locus | 로커스 | 유전자의 물리적 위치 |
| annotation | 주석 | 기능 · 구조 정보<br>염기 서열에 담긴 생물학적 의미와<br>기능적 특징을 '주석'처럼 표시 |
| BLAST | 블라스트 | 서열 유사성 탐색 툴<br>(구글 이미지 검색의 생물정보학 버전) |
| variant / variation | 변이 | 표준서열 대비 차이 |
| consequence (type) | 결과 유형 | 변이가 끼치는 분자적 영향 |
| condition | 임상 조건 | 질환·표현형 명칭 |
| classification | 분류 | Pathogenic / Benign 등 |
| ClinVar / dbSNP | 클린바 / dbSNP | 임상(의학판 위키백과) /<br>단일염기 변이 DB(어디에서 A→G가 흔한가?) |
| allele / genotype / phenotype | 대립유전자 / 유전자형 /<br>표현형 | 하디-바인베르크 연결 |

  ---

### Google Colab & Python 셸 (참고용)

| 용어 | 한국어 해석 | 설명 · 용도 |
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





 

---

## Google Colab 환경 구축

### Access Google Colab
- [colab.research.google.com](https://colab.research.google.com) 링크 접속
- 학교 구글 계정 로그인
  
### Create a New Notebook
- `File` 클릭 → `New Notebook` 클릭해 새 노트북 생성
- `+코드` 버튼을 누르면 코드를 입력할 수 있는 노트북이 나옵니다. <br>(기본으로 1개의 코드 셀이 생성되어 있습니다.)
- print 명령어를 사용해 봅시다.
 `print("Hello, Bioinformatics!")
  
### 필요한 library를 설치하자!
Google Colab은 많은 library들이 사전에 적용되어 있으나, 추가로 필요한 몇 가지의 library들을 
설치할 <br>필요가 있습니다. 예를 들어 `biopython`과 `scikit-bio`와 같은 library들이 있습니다.
<br>설치는 다음과 같은 명령어로 수행할 수 있습니다. `!pip install` 커맨드를 직접 셀에 입력하면 됩니다.
```python
!pip install biopython
```
```python
!pip install scikit-bio
```

### library를 import하고 설치를 확인하자!
새로운 코드 셀에 방금 설치한 library를 import하여 정확히 설치가 완료됐는지 확인합시다.
```python
# Importing necessary libraries
import Bio
import skbio

print("Biopython version:", Bio.__version__)
print("scikit-bio version:", skbio.__version__)
```
1.84 / 0.6.2 

## Python 기본 자료형 소개
- 변수, 문자열, 숫자 자료형 및 연산 예시

### Python 코드 실습 (아미노산 출력 예시)
- 코드 셀에 다음을 입력하세요.

```python
for i in [1, 2, 3, 4, 5]:
    print(i)
```

- for는 특정 명령을 반복하라고 할 때 사용합니다. 예를 들어 [1, 2, 3, 4, 5]라는 리스트에서
- 첫 번째엔 i=1, 두 번째엔 i=2, ... 가 대입되며, 그 때마다 print가 실행됩니다.
- for 변수 in (반복 가능한 개체) 의 형태로 사용하며, 아래 들여쓰기 된 부분을 반복합니다.

```python
amino_acids = ["Ala", "Gly", "Val", "Leu"]
for aa in amino_acids:
    print(aa)
```

## Google Colab 접속 및 환경 세팅
- 새 코드 셀 생성, 아래 내용 작성
  
  ```python
  from Bio.Seq import Seq
  print("Biopython 설치 성공!")
  ```

## Biopython 다루기 (기초)
  - 진행 순서
    ```python

    # DNA 서열 입력
    dna_seq = Seq("___") # 밑줄을 지우고 무작위로 A, T, G, C를 수십 개 넣어서 염기 서열을 작성해 보세요.

    # DNA → RNA 전사
    rna_seq = dna_seq.transcribe()
    print("RNA 서열:", rna_seq)

    # RNA → 아미노산 번역
    protein_seq = rna_seq.translate()
    print("단백질 서열:", protein_seq)
    ```

### Reading Frame 변화에 따른 아미노산 변화 확인
  - 목표: Reading Frame의 중요성을 이해하기
  - 진행 순서
    ```python

    # 원래 Reading Frame (0부터 시작)
    print("원래:", dna_seq.translate())

    # 1부터 시작
    print("Frame +1:", dna_seq[1:].translate())

    # 2부터 시작
    print("Frame +2:", dna_seq[2:].translate())
    ```
  - 세 종류의 Frame 결과를 비교하고, 왜 정확한 Frame이 중요한지 이유를 생각하고 발표해 봅시다.

### DNA 서열 변형 및 상보적 서열 확인
  - 목표: complement() 함수를 통해 DNA 서열의 구조 이해
  ```python
  # DNA 상보적 서열 염기
  complement_seq = data_seq.complement()
  print("원본 DNA 서열:", dna_seq)
  print("상보적 DNA 서열:", complement_seq)
  ```

---
