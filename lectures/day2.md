---
layout: default
title: "Day 2"
nav_order: 3
---

# Day 2: Biopython 환경 준비 및 기본 실습

## Google Colab 접속 및 환경 세팅
- 진행 순서<br>
  a. https://colab.research.google.com 접속<br>
  b. 새 노트북 만들기<br>
  c. 첫 번째 코드 셀에 아래 명령어 입력 및 실행:
  
  ```python
  !pip install biopython
  ```
     
  d. 설치 완료 후 아래 코드 실행하여 정상 설치 여부 확인:
  
  ```python
  from Bio.Seq import Seq
  print("Biopython 설치 성공!")
  ```

## Biopython 다루기 (기초)
  - 진행 순서
    ```python

    # DNA 서열 입력
    dna_seq = Seq("___")

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

## 미맹 키트를 이용한 미각 감지 여부 확인
  - 목표: 유전자형과 표현형의 관계의 직관적 이해
  - 진행 순서<br>
    a. 농도 별 PTC 용액을 면봉에 찍어 혀에 살짝 닿게 하여 쓴맛 감지 여부를 확인<br>
    b. 결과를 기록하여 학생 수 분포 확인

  - NCBI에서 미맹 관련 유전자 서열 불러오기
  - (PTC taste receptor gene, TAS2R38)
    
  ```python
  from Bio import Entrez, SeqIO   # Bio는 python에서 생물정보학 작업을 도와주는 Biopython 패키지이다.
                                  # Entrez 모듈은 NCBI 데이터베이스와 통신(데이터 요청)을, SeqIO 모듈은 염기서열(FASTA, GenBank 등) 파일을 읽고 쓰는 기능을 담당한다.

  Entrez.email = "your.email@example.com"  # NCBI 서버는 요청자의 신원파악을 위해 이메일을 요구한다.(자신의 실제 이메일주소 입력)
  handle = Entrez.efetch(db="nucleotide", id="AY258597", rettype="fasta", retmode="text")
                        # db="nucleotide" : 핵산(nucleotide) 데이터베이스에서 찾겠다는 의미
                        # id="AY258597"   : 요청할 접근번호이며, TAS2R38 유전자의 대표 서열 번호인 AY258597을 입력한다.
                        # rettype="fasta" : FASTA 형식으로 파일을 요청
                        # retmode="text"  : text 형태로 결과를 받겠다는 의미

  tas2r38_seq = SeqIO.read(handle, "fasta")  # SeqIO.read()는 handle에서 한 개의 FASTA 레코드를 읽어서 python의 객체(미지수)로 설정(정의)한다는 의미
  handle.close()        # 서버와 통신하던 통로(handle)을 닫아서, 컴퓨터의 소모 자원을 절약

  print("TAS2R38 유전자 서열:", tas2r38_seq.seq[:100], "...")  # tas2r38_seq.seq는 전체 염기서열(ATGC...)이다.
                                                              # [:100] 은 처음 100개만 잘라서 print하라는 의미이다.
  ```

## TAS2R38 서열 분석 및 SNP 위치 확인
  - 목표: 미각 유전자에서 실제로 발견되는 SNP를 확인
  - 진행 순서<br>
    a. TAS2R38의 알려진 SNP 위치(예: 145번째 뉴클레오타이드) 확인
    ```python

    # TAS2R38 유전자 SNP 위치(145번째 염기) 확인
    snp_position = 145 - 1                           # python에서 염기서열(sequence)은 0부터 시작하는 위치(index)로 관리된다.
                                                     # 즉, 첫 번째 염기는 seq[0], 두 번째는 seq[1] ... 이다.
                                                     # 우리가 '145번째' 염기를 찾고 싶으면 '144' 위치를 사용해야한다.
    
    snp_nucleotide = tas2r38_seq.seq[snp_position]   # tas2r38_seq.seq는 전체 염기서열(ATGC...)이므로, 
                                                     # 뒤에 [snp_position]을 붙이면 문자열의 144번째 위치, 즉 145번째 문자를 꺼내온다.
                                                     # 꺼낸 문자는 A, T, G, C 중 하나이므로 이를 snp_nucleotide라는 변수에 저장하라는 의미이다.
    
    print(f"SNP 위치의 뉴클레오타이드(145번째): {snp_nucleotide}")  # 문자열 앞에 f를 붙이면 f-string이라고 부른다.
                                                                  # f-string을 사용해야 {} 내부의 변수, 연산 등을 실제 값으로 평가해서 문자열에 삽입해준다.
                                                                  # f-string을 사용하지 않고 큰 따옴표만을 사용한다면, python을 그 내부를 그대로 텍스트로 인식한다.
                                                                  # 따라서 {snp_nucleotide}조차 그대로 단순한 문자의 나열 {, s, n, ..., }으로 인식하여 위에서 정의한 값을 불러오지 않는다.
    ```

## PTC 


  - 다른 친구들에 비해 나는 PTC 용액에 얼마나 민감하였나요?
  - 오늘 수행한 활동 결과 그 차이는 어디에서 기인한 것으로 생각되나요?



