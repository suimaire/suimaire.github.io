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
  - 여러 서열들이 나오므로, 그 중 하나를 선택하기 <br>(정상 서열 추천 ID: NM_000518.5 - 51번째 염기부터 89번째 염기)
  - 서열 ID를 클릭하여 상세 정보 페이지로 이동 <br>(NM_으로 시작하는 ID는 공식적인 표준 유전자 서열을 의미)

  - 새 NCBI 창을 열어 검색창에 변이 서열(낫 모양 적혈구 빈혈증 환자) 유전자 검색
    - 다음과 같은 검색어 입력

    ```css
    sickle cell anemia hemoglobin beta Homo sapiens
    ```

  - 마찬가지로 [nucleotide]를 선택하여 DNA 서열을 찾기
  - 여러 서열들 중 하나를 선택<br>(변이 서열 추천 ID: AY356351.1 - 159번째 염기부터 197번째 염기 )
  
### FASTA 형식으로 서열 보기 & biopython 분석
  - 화면 상단 메뉴에서 "FASTA" 클릭
  - 아래와 같은 형식의 서열이 나오면 해당 페이지를 유지한 채 <br>다시 Google Colab Notebook으로 이동

  - 두 번째 코드 셀에 아래 표를 참고하여 다음과 같이 입력

| 용어 | 한국어 해석 |
| Align | Biopython에서 제공하는 기능<br>두 개의 서열을 서로 맞추고 비교할 때 사용하는 도구 |
| PairwiseAligner | 두 서열을 비교하는 도구 |
| match | 서열이 같은 글자일 때 점수 부여 |
| mismatch | 서열이 다른 글자일 때 벌점 혹은 0점 부여 |
| gap | 빈칸, 서열을 맞추기 위한 공백<br>gap이 많으면 서열이 많이 다르다는 의미 |
| gap penalty | 서열을 맞추기 위해 gap을 넣을 때마다 주는 벌점 |
| ㅣ(세로선) | 두 서열에서 같은 글자가 위치할 때 표시 |
| sorted() | 여러 결과를 정리해 특정 기준에 따라 정렬 |

  - 코드 셀
    ```python
    # 필요한 모듈 임포트
    from Bio.Seq import Seq
    from Bio import Align 

    # DNA 서열 정의 
    normal_dna = Seq(" ")       # 큰 따옴표를 반드시 입력하고, 따옴표 사이에 아까 NCBI에서 찾은 서열을 입력
    mutated_dna = Seq(" ")       # 돌연변이가 일어난 낫 모양 적혈구 빈혈증 환자의 서열을 찾아서 입력
                                 # 어떻게 하면 필요한 서열만을 사용할 수 있을까? [__:__]을 사용해보자. 
      
    # 단백질 서열로 변환 (기존 코드와 동일)
    normal_prot = normal_dna.transcribe().translate()
    mutated_prot = mutated_dna.transcribe().translate()
    
    print("\n--- 단백질 서열 정렬 결과 ---")

    aligner = Align.PairwiseAligner() # A = B 꼴로 긴 문자열 명령어를 짧은 명령어로 '재정의'했다는 의미
                                      # PairwiseAligner: 두 서열을 서로 잘 맞추어 비교해주는 기능
    
    alignments = aligner.align(str(normal_prot), str(mutated_prot))
    # 정상 단백질과 변이 단백질을 aligner로 비교한 다음, 모든 결과를 alignments로 '재정의'했다는 의미

    for alignment in sorted(alignments, key=lambda x: x.score, reverse=True):
        print(alignment)
    # 정렬 결과 중 가장 잘 맞는(점수가 높은) 순서로 하나씩 꺼내어 화면에 출력하라는 의미
    ```

  - 정렬 결과가 다음과 같이 나왔다면 올바른 결과가 나온 것입니다.
  
  ```
  - --- 단백질 서열 정렬 결과 ---
    target            0 MVHLTP-EEKSAVT 13
                      0 ||||||-|-||||| 14
    query             0 MVHLTPVE-KSAVT 13
  ```

### [Q1] 정상/변이 단백질 정렬 결과에서 점선(ㅣ)이 사라지는 이유는 무엇일까요?
    
---

## 염기 하나의 변화가 단백질의 3차원 구조 형성에 영향을 미친다?

### 생각해보자!
  
####  하나의 염기 서열 변화가 가져올 수 있는 모든 변화의 경우에 대해 말해보자
      
      - 가장 위험한 하나의 염기 서열 변화는 무엇일까?
  
####  아미노산이 모여서 만들어지는 거대한 분자가 단백질이니까...

      - 아미노산이 바뀌면 단백질의 거대한 형태가 어떻게 바뀔까?

### 직접 단백질의 3차원 구조를 눈으로 확인해보자!
  
  - 왼쪽 메뉴의 [목차] → [+섹션] 클릭
  - 추가로 생성할 코드 셀의 제목 입력 (예시: 시뮬레이션: 세대별 변화 관찰)
  - 새롭게 섹션을 만들고 제목을 입력했다면, 그 부근에 마우스를 올려 [+코드] 클릭

  - 코드 셀에 다음과 같이 입력 

```python
import py3Dmol # 단백질의 3D 구조를 노트북 안에 표시해 주는 도구(py3Dmol)을 불러오라는 명령

def view_polymerization_with_labels():     # view_polymerization_with_labels(): 라는 이름의 함수(프로그램 묶음)을
    """                                    # 만드는 문장으로, 한 번 정의해 두면 나중에 저 이름만 쳐도 밑의 코드를 한 번에 불러올 수 있다.
    위쪽: 정상 헤모글로빈 단량체 (PDB:1A3N)
    아래쪽: 낫형 헤모글로빈 섬유    (PDB:2HBS)             
    β사슬 6번 위치에 Glu6 / Val6 레이블을 확실히 표시       
    """
    examples = [
        ('1A3N', False, 'Glu6', 'red', 'white'),      # examples라는 리스트를 만들어 여러가지 특징을 설정한다.
        ('2HBS', True,  'Val6', 'yellow', 'black'),   # 각각 (PDB 코드, 기능 사용 여부, 레이블(메모), 스틱 색, 글자 색)을 의미한다.
    ]
    
    for pdb, is_fiber, label_text, stick_color, font_color in examples:
    # examples 안에 두 개씩 들어있는 정보를 꺼내 pdb, is_fiber, label_text, stick_color, font_color 라는
    # 변수에 각각 대입해서 사용한다는 명령이다.

        # 1) 뷰어 만들기
        view = py3Dmol.view(query=f'pdb:{pdb}', width=350, height=350)
        view.setBackgroundColor('white')
        # PDB 데이터베이스에서 해당 pdb 코드(ex: 1A3N or 2HBS)의 3D 구조를 불러와 캔버스를 만든다.
        # width/height: 보이는 크기를 350x350 픽셀로 지정한다.
        # 뷰어의 배경색을 white로 설정하여 단백질을 잘 보이게 한다.

        # 2) 리본 스타일
        if is_fiber:
            # 섬유: 연한 회색, 반투명
            view.setStyle({'cartoon': {}}, {'cartoon': {'color':'lightgrey','opacity':0.4}})
            # 섬유 표면
            view.addSurface({'opacity':0.3,'color':'lightgreen'})
        else:
            # 단량체: 스펙트럼 컬러, 불투명
            view.setStyle({'cartoon': {}}, {'cartoon': {'color':'spectrum','opacity':1.0}})

        # is_fiber에 True, False가 대입될 때 어떠한 형태로 그릴지 정의한다.
        # True일 때: 연한 회색, 반투명(opacity=0.4)로 그리라는 의미
        # False일 때: 무지개 색(spectrum), 불투명(opacity=1.0)으로 그리라는 의미


        # 3) 6번 잔기 stick 강조
        view.addStyle(
            {'chain':'B','resi':6},                          # addStyle(A, B): A부분을 B모양으로 표현하라는 함수
            {'stick':{'color': stick_color, 'radius':0.2}}   # 여기서는 B chain의 6번 자리를 stick 모양으로 강조하라는 의미
        )
        
        # 4) 레이블 추가 — 세 개의 인자를 정확히 분리
        view.addLabel(
            label_text,                                       # addLabel(A, B, C):
            {                                                 # A: 그림 안에 표시할 실제 메모(Label) 내용
                'fontColor': font_color,                      # B: 글자 색(fontColor), 배경 색, 글자 크기
                'backgroundColor': stick_color,               
                'fontSize': 12
            },
            {'chain':'B','resi':6,'atom':'CA'}                # C: 메모(Label)을 붙일 위치 
        )                                                     # 여기서는 B chain의 6번 자리의 CA 원자(알파 탄소)
        
        # 5) 6번 잔기로 줌
        view.zoomTo({'chain':'B','resi':6})                   # 선택된 부위(6번 자리) 중심으로 확대해서 보여주라는 의미
        
        # 6) 제목 출력 후 렌더링                               
        title = '정상 단량체 (PDB:1A3N)' if not is_fiber else '낫형 섬유 (PDB:2HBS)'
        print(title)                                          # is_fiber 여부에 따라 화면 위에 title을 다르게 표시하라는 명령
        view.show()

# 함수 호출
view_polymerization_with_labels()                             # 마지막으로 함수 이름 뒤에 ()를 붙여서 지금까지 정의한 모든 과정을 한꺼번에 실행한다.
```

### 마무리
  - 1. 정상/변이 단백질 정렬 결과에서 점선(ㅣ)이 사라지는 이유는 무엇일까요?
    2. 20종류의 아미노산은 각각 친수성, 소수성, 산성, 염기성과 같은 성격이 모두 다릅니다. 이것은 단백질의 구조에 어떤 영향을 미칠 수 있을까요?
    3. 3D 모델을 관찰한 결과 염기 하나의 변화로 해당 아미노산은 단백질의 3차원 구조 중 어디서 어디로 이동하였나요? <br>그리고, 그 이유는 무엇일까요?
