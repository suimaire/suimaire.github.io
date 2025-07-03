---
layout: default
title: "Day 2"
nav_order: 3
---

# Day 2: Biopython 환경 준비 및 기본 실습

## 1. Google Colab 접속 및 환경 세팅
- 진행 순서
  1. https://colab.research.google.com 접속
  2. 새 노트북 만들기
  3. 첫 번째 코드 셀에 아래 명령어 입력 및 실행:
     ```python
     !pip install biopython
     ```
  4. 설치 완료 후 아래 코드 실행하여 정상 설치 여부 확인:
     ```python
     from Bio.Seq import Seq
     print("Biopython 설치 성공!")
     ```

## 2. Biopython 다루기 (기초)
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
## 3. Reading Frame 변화에 따른 아미노산 변화 확인
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

## 4. DNA 서열 변형 및 상보적 서열 확인
  - 목표: complement() 함수를 통해 DNA 서열의 구조 이해
  ```python
  # DNA 상보적 서열 염기
  complement_seq = data_seq.complement()
  print("원본 DNA 서열:", dna_seq)
  print("상보적 DNA 서열:", complemetn_seq)
  ```

**실습 자료**  
- [Colab 노트북](lectures/day1_notebook.ipynb)  
