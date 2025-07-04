---
layout: default
title: "Day 3"
nav_order: 4
---

# Day 3 낫 모양 적혈구 빈혈증 탐색

---
## 낫 모양 적혈구 빈혈증

  - 헤모글로빈: 적혈구에 존재하는 산소 운반 단백질
  - 낫 모양의 적혈구를 특징으로 하는 헤모글로빈의 유전성 유전자 이상으로 발생
  - 이상 적혈구의 과다 파괴가 원인이 되는 만성 빈혈
  - 낫 모양 적혈구는 생존 기간이 짧으며, 서로 달라붙어 딱딱한 형태로 작은 모세혈관을 막기 쉬움
  - 모세혈관이 막힌 부위에서는 조직으로의 산소 전달이 어려움

---
### Google Colab Notebook 열기
  - [GOOGLE COLAB](https://colab.research.google.com)
  - 왼쪽 위 [파일] → Drive의 새 노트북

### biopython 분석 준비
  - 첫 번째 코드 셀에 다음과 같이 입력

    ```python
    !pip install biopython py3Dmol ipywidgets
    from Bio.Seq import Seq
    from Bio import Align
    import matplotlib.pyplot as plt
    from ipywidgets import interact
    import py3Dmol
    ```

  - 코드 셀 왼쪽의 [셀 실행] (ctrl + enter) 클릭 후 대기
  - Successfully installed biopython-1.85 jedi-0.19.2 py3Dmol-2.5.1 결과값이 나오면<br>정상 설치가 된 것

### 노트 목차 나누기 및 새로운 코드 셀 생성

  - 왼쪽 메뉴의 [목차] → [+섹션] 클릭
  - 추가로 생성할 코드 셀의 제목 입력 (예시: 정상 서열 vs 변이 서열)
  - 새롭게 섹션을 만들고 제목을 입력했다면, 그 부근에 마우스를 올려 [+코드] 클릭

### NCBI에서 정상 서열과 낫 모양 적혈구 빈혈증 환자 서열 찾기
  - NCBI에서 정상인과 환자의 염기 서열을 찾아보자!
  - [NCBI](https://www.ncbi.nlm.nih.gov) 클릭하여 염기서열 데이터베이스 접속
  - NCBI 메인 페이지 상단 검색창에 정상-베타-글로빈 유전자 검색
    - 다음과 같은 검색어 입력

    ```css
    hemoglobin beta Homo sapiens normal
    ```

  - 여러 결과들 중 [nucleotide]를 선택하여 DNA 서열만 표시되도록 설정
  - 여러 서열들이 나오므로, 그 중 하나를 선택하기 <br>(정상 서열 추천 ID: NM_000518.5)
  - 서열 ID를 클릭하여 상세 정보 페이지로 이동 # NM_으로 시작하는 ID는 공식적인 표준 유전자 서열을 의미

### FASTA 형식으로 서열 보기
  - 화면 상단 메뉴에서 "FASTA" 클릭
  - 아래와 같은 형식의 서열이 나오면 해당 페이지를 유지한 채 <br>다시 Google Colab Notebook으로 이동

  - 두 번째 코드 셀에 아래 표를 참고하여 다음과 같이 입력
  - | 용어 | 한국어 해석 |
  - | PairwiseAligner | 두 서열을 비교하는 도구 |
  - | match | 서열이 같은 글자일 때 점수 부여 |
  - | mismatch | 서열이 다른 글자일 때 벌점 혹은 0점 부여 |
  - | gap | 빈칸, 서열을 맞추기 위한 공백<br>gap이 많으면 서열이 많이 다르다는 의미 |
  - | gap penalty | 서열을 맞추기 위해 gap을 넣을 때마다 주는 벌점 |
  - | ㅣ(세로선) | 두 서열에서 같은 글자가 위치할 때 표시 |
  - | sorted() | 여러 결과를 정리해 특정 기준에 따라 정렬 |

    ```python
    # 필요한 모듈 임포트
    from Bio.Seq import Seq
    from Bio import Align 

    # DNA 서열 정의 
    normal_dna = Seq(" ")       # 큰 따옴표를 반드시 입력하고, 따옴표 사이에 NCBI에서 찾은 서열을 입력
    mutated_dna = Seq(" ")       # 돌연변이가 일어난 낫 모양 적혈구 빈혈증 환자의 서열을 찾아서 입력

    # 단백질 서열로 변환 (기존 코드와 동일)
    normal_prot = normal_dna.transcribe().translate()
    mutated_prot = mutated_dna.transcribe().translate()

    print("\n--- 단백질 서열 정렬 결과 ---")

    aligner = Align.PairwiseAligner
    
    
---
