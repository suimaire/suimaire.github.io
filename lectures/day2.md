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
  from Bio import Entrez, SeqIO
  Entrez.email = "your.email@example.com"
  handle = Entrez.efetch(db="nucleotide", id="AY258597", rettype="fasta", retmode="text")
  tas2r38_seq = SeqIO.read(handle, "fasta")
  handle.close()

  print("TAS2R38 유전자 서열:", tas2r38_seq.seq[:100], "...")
  ```

## TAS2R38 서열 분석 및 SNP 위치 확인
  - 목표: 미각 유전자에서 실제로 발견되는 SNP를 확인
  - 진행 순서<br>
    a. TAS2R38의 알려진 SNP 위치(예: 145번째 뉴클레오타이드) 확인
    ```python

    # TAS2R38 유전자 SNP 위치(145번째 염기) 확인
    snp_position = 145 - 1
    snp_nucleotide = tas2r38_seq.seq[snp_position]
    print(f"SNP 위치의 뉴클레오타이드(145번째): {snp_nucleotide}")
    ```

  - 다른 친구들에 비해 나는 PTC 용액에 얼마나 민감하였나요?
  - 오늘 수행한 활동 결과 그 차이는 어디에서 기인한 것으로 생각되나요?



