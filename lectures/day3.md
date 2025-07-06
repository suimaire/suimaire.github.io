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
  - 여러 서열들 중 하나를 선택<br>(변이 서열 추천 ID: AY356351.1 - 155번째 염기부터 193번째 염기 )
  
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

    # 단백질 서열로 변환 (기존 코드와 동일)
    normal_prot = normal_dna.transcribe().translate()
    mutated_prot = mutated_dna.transcribe().translate()

    print("\n--- 단백질 서열 정렬 결과 ---")

    aligner = Align.PairwiseAligner # A = B 꼴로 긴 문자열 명령어를 짧은 명령어로 '재정의'했다는 의미
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

## 말라리아와 빈혈의 치명률 변화가 후대로의 형질 전달에 영향을 미친다?

### 만일, 낫 모양 적혈구 빈혈에 걸린 사람이 100% 확률로 결혼 전에 죽는다면?
  - 예상되는 해당 염기 서열이 자손에게 전달 될 확률은?

### 만일, 낫 모양 적혈구 빈혈에 걸린 사람이 50% 확률로 10살, 50% 확률로 40살에 죽는다면?
  - 예상되는 해당 염기 서열이 자손에게 전달 될 확률은?

### 1.2.2의 상황에서, 해당 지역에 말라리아 병원체를 가지는 모기가 900만마리가 들어왔다!
  - 말라리아에 걸렸을 때 100% 확률로 3일 뒤에 죽는다면, 예상되는 빈혈 환자의 염기 서열이 후손에게 전달 될 확률은?

### 눈으로 확인해보자!
  - 새 코드 셀을 만들고, 아래 코드를 입력하자.
    ```python
    # simulate라는 함수를 정의하자! [def] 이후, 괄호 안에 4가지 파라미터를 넣는다. 
    # selection_sickle: SS(낫모양 동형접합)의 생존율 (ex:0.9 → 정상보다 10% 덜 산다.)
    # malaria_mortality: AA(정상 동형접합)의 사망률 (ex:0.2 → 말라리아에 걸린 환자 중 20%가 죽는다.)
    # generations: 시뮬레이션을 진행 할 세대 수
    # initial_S_allele: 최초 S 대립유전자의 빈도(ex:0.1 → 10%)
    
    def simulate(selection_sickle=0.9, malaria_mortality=0.2, generations=30, initial_S_allele=0.1):
    p = 1 - initial_S_allele; q = initial_S_allele
    AA=AS=SS=None
    AA_lst, AS_lst, SS_lst = [], [], []
    for _ in range(generations):
        f_AA, f_AS, f_SS = p*p, 2*p*q, q*q
        w_AA, w_AS, w_SS = 1-malaria_mortality, 1, 1-selection_sickle
        mean_w = f_AA*w_AA + f_AS*w_AS + f_SS*w_SS
        f_AA, f_AS, f_SS = f_AA*w_AA/mean_w, f_AS*w_AS/mean_w, f_SS*w_SS/mean_w
        AA_lst.append(f_AA); AS_lst.append(f_AS); SS_lst.append(f_SS)
        p = f_AA + 0.5*f_AS; q = 1 - p
    plt.figure(figsize=(5,3))
    plt.plot(AA_lst,label="AA"); plt.plot(AS_lst,label="AS"); plt.plot(SS_lst,label="SS")
    plt.xlabel("generation"); plt.ylabel("frequency"); plt.legend(); plt.title("malaria selection pressure")
    plt.show()
    interact(simulate, selection_sickle=(0.5,1.0,0.05), malaria_mortality=(0,0.5,0.05),
        generations=(10,100,10), initial_S_allele=(0.01,0.3,0.01));

    ```
    
---
### 다시
  ```python
  import py3Dmol

  def view_polymerization():
    # --- 1) 정상 단량체 (PDB:1A3N) ---
    v1 = py3Dmol.view(width=400, height=400, query='pdb:1A3N')
    # 리본(cartoon) 스타일로 전체 구조
    v1.setStyle({'cartoon': {}})
    # β체인 B의 6번 잔기(=Glu6)를 빨간 스틱으로 강조
    v1.addStyle({'chain':'B','resi':6}, {'stick':{'color':'red','radius':0.4}})
    v1.setBackgroundColor('white')
    
    # --- 2) 낫형 헤모글로빈이 만든 섬유 (PDB:2HBS) ---
    v2 = py3Dmol.view(width=400, height=400, query='pdb:2HBS')
    # 리본 스타일
    v2.setStyle({'cartoon': {}})
    # β체인 B의 6번 잔기(=Val6)를 노란 스틱으로 강조
    v2.addStyle({'chain':'B','resi':6}, {'stick':{'color':'yellow','radius':0.4}})
    # 전체 섬유 모양이 드러나도록 살짝 투명한 스틱으로 주변 결합부도 표시
    v2.addStyle({'cartoon':{ }}, {'cartoon':{'opacity':0.3}})
    v2.setBackgroundColor('white')
    
    # 두 뷰어를 나란히 출력
    v1.show()
    v2.show()
  ```
