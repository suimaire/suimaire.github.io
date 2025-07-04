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
  - 왼쪽 메뉴의 [목차] → [+섹션] 클릭
  - 추가로 생성할 코드 셀의 제목 입력 (예시: 정상 서열 vs 변이 서열)
  - 새롭게 섹션을 만들고 제목을 입력했다면, 그 부근에 마우스를 올려 [+코드] 클릭
  - 두 번째 코드 셀에 다음과 같이 입력

    ```python
    # 필요한 모듈 임포트
    from Bio.Seq import Seq
    from Bio import Align 

    
---
