---
layout: default
title: "Day 2"
nav_order: 3
---

# Day 2: Biopython 환경 준비 및 기본 실습

## 미맹 키트를 이용한 미각 감지 여부 확인
  - 목표: 유전자형과 표현형의 관계의 직관적 이해
  - 진행 순서<br>
    a. 농도 별 PTC 용액을 면봉에 찍어 혀에 살짝 닿게 하여 쓴맛 감지 여부를 확인<br>
    b. 결과를 기록하여 학생 수 분포 확인

  - NCBI에서 미맹 관련 유전자 서열 불러오기
  - (PTC taste receptor gene, TAS2R38)


### Access Google Colab
- [colab.research.google.com](https://colab.research.google.com) 링크 접속
- 학교(혹은 개인) 구글 계정 로그인
  
### Create a New Notebook
- `File` 클릭 → `New Notebook` 클릭해 새 노트북 생성
- 기본 코드 셀에 !pip install biopython을 작성하여 biopython을 설치한 뒤 새로운 코드 셀을 생성
- 아래와 같이 작성
  
  ```python
  from Bio import Entrez, SeqIO
  from Bio.Seq import Seq

  # Bio는 python에서 생물정보학 작업을 도와주는 Biopython 패키지이다.
  # Entrez 모듈은 NCBI 데이터베이스와 통신(데이터 요청)을, SeqIO 모듈은 염기서열(FASTA, GenBank 등) 파일을 읽고 쓰는 기능을 담당한다.


  # 1) 쓴 맛을 느끼는 사람인가요? 아닌가요?
  ptc_group = "____"   # "밑줄"에 쓴맛을 느끼면 "taster", 못 느끼면 "nontaster" 입력

  # 2) 그룹에 따라 자동으로 accession 선택
  if ptc_group == "taster":
      accession = "AY258598"   # PAV (쓴맛 잘 느낌)
  else:
      accession = "AY258597"   # AVI (쓴맛 못 느낌)


  Entrez.email = "your.email@example.com"  # NCBI 서버는 요청자의 신원파악을 위해 이메일을 요구한다.(자신의 실제 이메일주소 입력)
  handle = Entrez.efetch(
          db="nucleotide",     # db="nucleotide" : 핵산(nucleotide) 데이터베이스에서 찾겠다는 의미
          id=accession,        # id=accession    : 요청할 접근번호이며, TAS2R38 유전자의 대표 서열 번호인 AY258597 또는 AY258598을 입력된다.
          rettype="fasta",     # rettype="fasta" : FASTA 형식으로 파일을 요청
          retmode="text"       # retmode="text"  : text 형태로 결과를 받겠다는 의미
)                              # 이렇게 가독성을 위해 인수별로 줄바꿈을 하기도 한다.
                        
                        
                        
  tas2r38_seq = SeqIO.read(handle, "fasta") # SeqIO.read()는 handle에서 한 개의 FASTA 레코드를 읽어서 python의 객체(미지수)로 설정(정의)한다는 의미
  handle.close()        # 서버와 통신하던 통로(handle)을 닫아서, 컴퓨터의 소모 자원을 절약

  print(f"TAS2R38 유전자 서열 ({ptc_group}): {tas2r38_seq.seq[:200]} …")  
     # tas2r38_seq.seq는 전체 염기서열(ATGC...)이다.
     # [:200] 은 처음 200개만 잘라서 print하라는 의미이다.
     # 문자열 앞에 f를 붙이면 f-string이라고 부른다.
     # f-string을 사용해야 {} 내부의 변수, 연산 등을 실제 값으로 평가해서 문자열에 삽입해준다.
     # f-string을 사용하지 않고 큰 따옴표만을 사용한다면, python을 그 내부를 그대로 텍스트로 인식한다.
     # 따라서 {snp_nucleotide}조차 그대로 단순한 문자의 나열 {, s, n, ..., }으로 인식하여 위에서 정의한 값을 불러오지 않는다.
  ```

## TAS2R38 서열 분석 및 SNP 위치 확인
  - 목표: 미각 유전자에서 실제로 발견되는 SNP를 확인
  - SNP : **S**ingle **N**ucleotide **p**olymorphism의 약자로, 유전 정보에서 한 자리 염기만 다른 유전 형질
  - 마치 백과사전에서 하나의 오타가 난 것과 같은 아주 작은 차이
    <br>
  - 진행 순서<br>
    a. TAS2R38의 알려진 SNP 위치(예: 145번째 뉴클레오타이드) 확인
    ```python

    # TAS2R38 유전자 SNP 위치(100번째 염기) 확인
    snp_position = 100 - 1                           # python에서 염기서열(sequence)은 0부터 시작하는 위치(index)로 관리된다.
                                                     # 즉, 첫 번째 염기는 seq[0], 두 번째는 seq[1] ... 이다.
                                                     # 우리가 '100번째' 염기를 찾고 싶으면 '99' 위치를 사용해야한다.
    
    snp_nucleotide = tas2r38_seq.seq[snp_position]   # tas2r38_seq.seq는 전체 염기서열(ATGC...)이므로, 
                                                     # 뒤에 [snp_position]을 붙이면 문자열의 99번째 위치, 즉 100번째 문자를 꺼내온다.
                                                     # 꺼낸 문자는 A, T, G, C 중 하나이므로 이를 snp_nucleotide라는 변수에 저장하라는 의미이다.
    
    print(f"SNP 위치의 뉴클레오타이드(100번째): {snp_nucleotide}")  
    ```

## 서로 다른 두 FASTA 염기 서열에서 다른 부분을 빠르게 비교해보자
  - 목표: 우리는 200개의 염기를 들여다보며 어디가 다른지 찾고 있을 시간이 없다.
  - Biopython 명령어를 이용하여 다른 부분을 빠르게 찾아보자.
  - 코드 셀에 다음과 같이 작성한다.

```python
from Bio import Entrez, SeqIO

def fetch_sequence(accession, email="your.email@example.com"):

    # 다음 줄과 같이 큰따옴표 3개를 양옆에 붙여 주석을 작성하기도 한다.
    """주어진 accession의 FASTA를 NCBI에서 가져와 SeqRecord로 반환하라는 명령이다"""
    

    Entrez.email = email
    handle = Entrez.efetch(
        db="nucleotide",
        id=accession,
        rettype="fasta",
        retmode="text"
    )
    record = SeqIO.read(handle, "fasta")   
    handle.close()
    return record                          
    """
    SeqIO.read()는 handle에서 한 개의 FASTA 레코드를 읽어서 python의 객체(미지수)로 설정(정의)한다는 의미
    호출한 쪽으로 SeqRecord 객체를 돌려준다는 의미이며, SeqRecord 객체는 .seq 형식으로 염기 서열 정보가 들어있다.
    """

def find_differences(rec1, rec2):
    """
    두 SeqRecord를 비교해서, 같은 인덱스(1-based)에서
    염기가 다른 위치만 [(위치, rec1_염기, rec2_염기), ...] 형태로 반환.
    """
    seq1 = rec1.seq
    seq2 = rec2.seq
    diffs = []
    # 두 서열 길이가 같다고 가정한다. 서열 길이가 다르면 짧은 쪽까지 비교한다.
    for i, (nt1, nt2) in enumerate(zip(seq1, seq2), start=1):    # enumerate(..., start=1) : 1부터 시작하는 위치 번호를 매긴다.
                                           # zip(seq1, seq2)로 두 서열을 같은 길이만큼 한 글자씩 묶어서 처리하게 만든다.
        if nt1 != nt2:                     # 두 서열의 같은 위치 염기가 다르면, diffs 리스트에 (위치, rec_염기, rec_염기)를 추가한다는 의미
            diffs.append((i, nt1, nt2))
    return diffs                           # 염기가 달랐던 위치와 각 서열의 염기를 담은 리스트를 반환한다는 의미

# 다음과 같은 코드를 작성하여 비교 출력
""" ───────────────────── 비교 출력 ──────────────────────────── """
# 1) 두 accession 불러오기
rec_avi = fetch_sequence("AY258597")   # nontaster (AVI)
rec_pav = fetch_sequence("AY258598")   # taster    (PAV)

# 2) 차이나는 위치 찾기 : 위에서 정의한 find_differences를 호출하여 염기가 다른 모든 위치를 diff_positions에 저장한다.
diff_positions = find_differences(rec_avi, rec_pav)

# 3) 결과 출력
print("서열 비교 결과 (위치: AVI → PAV):")
for pos, nt_avi, nt_pav in diff_positions:   # for : 아래 들여쓰기한 부분을 반복해서 실행하라는 명령
    print(f"  {pos}번째: {nt_avi} → {nt_pav}")    # for 변수 in (반복 가능한 개체) 의 형태로 사용한다.

```

---
### 질문에 대해 답해봅시다.

  - 다른 친구들에 비해 나는 PTC 용액에 얼마나 민감하였나요?
  - 오늘 수행한 활동 결과 그 차이는 어디에서 기인한 것으로 생각되나요?
  - 어느 SNP에서 차이가 발생하나요?
  - 그렇다면, **몇 번째 아미노산**이 바뀌어서 쓴 맛 민감성에 차이가 발생하는 걸까요?

---

## 마무리(심화)
  - 락타아제 지속성(Lactase persistence)과 음주 후 피부가 붉어지는 것(Alcohol flushing)에 대해서도 위의 활동을 참고하여 조사해 봅시다.
  - db를 "snp"로 변경하고, id의 변수 값을 하나로 고정한 뒤, rettype을 "flt"로 변경해주면 SNP 주변의 염기서열을 두 개의 FASTA 파일로 돌려줍니다.
  - 이후 SeqIO.read 대신 SeqIO.parse를 사용하면 반환된 텍스트에서 FASTA 형식을 여러 개 읽어 옵니다.
  - 추가로 기존의 record에 대해 정의한 코드를 records = list(SeqIO.parse(handle, "fasta")) 과 같이 변경하고 records를 return 합니다.
  - 이후 사용 예시의 함수를 참고하여 작성하면 락타아제 지속성 및 Alcohol flushing의 염기서열을 직접 비교할 수 있습니다.
  - 이 설명을 읽고 직접 코드를 작성해 보세요!

