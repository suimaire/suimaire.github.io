---
layout: default
title: "Day 4"
nav_order: 5
---

# Day 4 하디-바인베르크 평형<br>(Hardy-Weinberg Equilibrium) 실습

---
## 하디-바인베르크 유전 평형

  - 이상적인 조건을 가정한 집단(멘델 집단)에서, 시간이 흘러도 대립유전자 빈도와 유전자형 빈도가 유지되는 것
  - 우성인 유전자는 열성인 유전자에 비해 우선적으로 발현되는데, 왜 세대를 거듭하며 우성 유전자를 가지는<br>개체의 수가 열성 유전자를 가지는 개체보다 많아지지 않는가?
  - 열성 유전자를 가지는 개체가 세대를 거듭해도 사라지지 않는 이유는 무엇인가?에 대한 해답을 간접적으로 제시
  - 이상적인 조건을 가정한 집단에서만 유지되는 것으로 실제 생태계에서 이론이 성립하는 경우는 거의 없다.
  
### 하디-바인베르크 유전 평형의 원리

  - 부모의 표현형이 우성이라고 해서 자손에게 우성 유전자만 전달되는 것이 아니다.
  - 어머니와 아버지가 자손에게 제공하는 생식세포는 무작위적으로 부모 본인의 유전 정보의 절반만 제공하는 것이기 때문이다.
  - 어머니에게 우성 유전자 A를 받고, 아버지에게 열성 유전자 a를 받은 자손 개체(Aa)는 표현형이 우성(A)일 뿐,<br>다음 세대에게 A를 물려줄 수도 있고 a를 물려줄 수도 있다.
  - 따라서 개체 수가 일정 수준 이상인 집단에서는, 유전에 별도로 영향을 끼치는 요인 없이는 대립 유전자의 빈도가 유지된다.

### 하디-바인베르크 유전 평형 성립 조건

  - 무작위적 교배: 이성에게 보다 특별하게 매력적인 개체는 존재하지 않으며, 오직 무작위적으로 짝을 이루어 세대를 거듭한다.
  - 돌연변이 없음: 돌연변이는 아주 드물게 발생하지만, 특정 경우에는 기존 형질에 비해 적응 능력이 높거나 낮은 개체가 등장할 수 있다.
  - 자연선택 없음: 우성 및 열성 형질 개체를 비교할 시 생존에 보다 더 유리한 형질은 존재하지 않는다.
  - 충분히 큰 개체군: 집단의 수가 충분히 커야 확률이 의미하는 숫자에 수렴할 수 있다.
  - 유전자 흐름 없음: 주변 다른 집단과의 개체 이동, 교류, 교배는 없다고 가정한다.

### 하디-바인베르크 유전 평형

  - 어떤 생물적 특성은 2가지 형질로 나뉜다. (ex: 발가락이 5개인 사람과 6개인 사람)
  - 발가락이 5개인 형질이 우성이며 A이고, 발가락이 6개인 형질이 열성이며 a라고 정의한다.
  - 발가락이 5개인 암컷과 수컷(모두 Aa 개체)이 자손을 낳는다면, 그 자손에서 예상되는 대립유전자의 비율은 다음과 같다.
  
|   | A  | a  |
| A | AA | Aa |
| a | Aa | aa |

  - 이제 멘델 집단 내부에서 A 유전자의 빈도가 p, a 유전자의 빈도가 q라고 하자.<br>형질의 종류가 2가지이므로 유전자 빈도의 총합은 p+q=1이다.
  - 예를 들어 p=0.6, q=0.4 일 때 집단 내부에 개체가 100마리 존재한다면, 개체들은 각각 2개의 <br>유전자를 가지므로 총 200개의 A 또는 a가 집단 내부에 존재한다는 의미이다. <br>이 때 p가 0.6이므로, 200개중 120개의 유전자는 A이며 q가 0.4이므로 80개의 유전자는 a이다.
  - 해당 빈도가 멘델 집단 내에서는 세대를 거듭하여도 유지된다는 것이 하디-바인베르크 법칙이다.<br>증명은 해당 아카데미와 관련이 없으므로 생략한다.
  - 위 표를 멘델 집단에서 부모가 자손에게 특정 유전자를 전달할 확률이라고 생각할 수 있다.<br>**p=0.6, q=0.4 일때 확률은 다음과 같다.**

|   | p    | q    |
| p | 0.36 | 0.24 |
| q | 0.24 | 0.16 |

  - 이를 biopython을 이용해 시각적으로 관찰하며 이해해보자.

## Biopython을 이용한 Hardy-Weinberg Equilibrium의 이해
  - 학습 목표
    - 대립유전자 빈도 p, q의 이해
    - 세대 반복 시 유전자형의 빈도가 p² : 2pq : q² 로 수렴함을 확인
    - FASTA/VCF 등 실측 서열에서 대립유전자 빈도를 계산해 가설 검정

### Google Colab Notebook 열기
  - [GOOGLE COLAB](https://colab.research.google.com)
  - 왼쪽 위 [파일] → Drive의 새 노트북

### biopython 분석 준비
  - 첫 번째 코드 셀에 다음과 같이 입력

    ```python
    !pip install biopython ipywidgets
    import random, matplotlib.pyplot as plt
    from ipywidgets import interact
    ```

  - 코드 셀 왼쪽의 [셀 실행] (ctrl + enter) 클릭 후 대기
  - Successfully installed biopython-1.85 jedi-0.19.2 결과값이 나오면 정상 설치가 된 것

### 노트 목차 나누기 및 새로운 코드 셀 생성

  - 왼쪽 메뉴의 [목차] → [+섹션] 클릭
  - 추가로 생성할 코드 셀의 제목 입력 (예시: 기대 유전자형 빈도 함수)
  - 새롭게 섹션을 만들고 제목을 입력했다면, 그 부근에 마우스를 올려 [+코드] 클릭

### 관찰할 시뮬레이션의 조건을 입력

  - 코드 셀에 다음과 같이 입력

    ```python
    def hwe_expected(p: float): # def는 함수를 정의하라는 명령, hwe_expected는 "Hardy-Weinberg expected"의 약자
                                # (p: float) : p라는 값을 넣을 건데, float는 "소수점 숫자"라는 의미 (ex: p=0.6같은 값을 넣을거야)
        q = 1 - p               # p + q = 1 이 되어야 한다!
        return p*p, 2*p*q, q*q  # 해당 함수가 세 가지 값을 한꺼번에 돌려준다는 의미 (p², 2pq, q²)
    print("예시 p=0.6 :", hwe_expected(0.6)) # print는 화면에 내용을 보여주라는 명령, "예시 p=0.6:" 라는 글자 뒤에 결과를 출력
    ```

### 세대별 변화 관찰 시뮬레이션

  - 왼쪽 메뉴의 [목차] → [+섹션] 클릭
  - 추가로 생성할 코드 셀의 제목 입력 (예시: 시뮬레이션: 세대별 변화 관찰)
  - 새롭게 섹션을 만들고 제목을 입력했다면, 그 부근에 마우스를 올려 [+코드] 클릭

  - 코드 셀에 다음과 같이 입력

    ```python
    def simulate_hwe(p_init=0.7, pop_size=1000, generations=5000):  # simulate_hwe라는 이름의 함수 생성
                                                                    # p_init: 첫 세대의 A 대립유전자 빈도(0.7 = 70%)
    p = p_init                                                      # pop_size: 각 세대의 개체 수 / generations: 반복할 세대 수
    AA= []; Aa=[]; aa=[]
    for _ in range(generations):           # _ 는 "반복 횟수 자체는 중요하지 않을 때" 사용 , range(generations)만큼 세대별 분석
        AA_exp, Aa_exp, aa_exp = hwe_expected(p)  # p값을 설정함에 따라 그에 맞는 값을 반환하라는 명령
        AA.append(AA_exp); Aa.append(Aa_exp); aa.append(aa_exp)
        # 무작위 짝짓기 개체 생성
        import random                             # python 내장 모듈 random을 사용하라는 명령
                                                  # 이론적 비율에 따라 실제로 pop_size명의 유전형을 뽑는 효과 발생
        pop = random.choices(['AA','Aa','aa'], weights=[AA_exp, Aa_exp, aa_exp], k=pop_size)
        allele_A = pop.count('AA')*2 + pop.count('Aa')  # allele = 대립유전자
                                                        # 'AA' 개체는 A를 2개, 'Aa' 개체는 A를 1개를 가지므로 이렇게 계산하라는 명령
        p = allele_A/(2*pop_size)                       # p를 구할 때, A 대립유전자 수를 전체 인구수*2로 나눈 값으로 구하라는 명령
    plt.figure()
    plt.plot(AA, label='AA'); plt.plot(Aa, label='Aa'); plt.plot(aa, label='aa')
    plt.xlabel('Generation'); plt.ylabel('Frequency'); plt.title('Hardy–Weinberg simulation')
    plt.legend(); plt.show()                            # AA, Aa, aa 리스트에 담긴 세대별 비율은 선 그래프로 그리라는 명령
                                                        # x축 = 세대 번호, y축 = 유전형 빈도
                                                        # legend = 범례, title = 제목 (마음대로 바꿔도 됩니다)

    interact(simulate_hwe, p_init=(0.1,0.9,0.05), pop_size=(50,1000,50), generations=(5,5000,5));
                                                        # Colab 환경에서 사용되는 interact 명령
                                                        # p_init, pop_size, generations 값을 슬라이더로 바꿔가며 시뮬레이션 갱신
    ```

  - 이론적 예상치를 구한 뒤, 그 예상대로 유한한 개체수에서 실제 뽑기
  - 그 결과로 얻은 실제 빈도로 p값(A 대립유전자 빈도) 계산

---

### 생각해보자

  - 직접 값을 변경해 보며 하디-바인베르크 유전 평형을 직접 체험해 봅시다.
  - 집단 크기가 감소하면 유전자 빈도의 변동성은 어떻게 변화하나요?
  - 이론과 실제 시뮬레이션 수치가 맞지 않다면, 그 이유는 무엇일까요?

---

## 빈혈과 말라리아의 상관 관계 - 선택 압력 시뮬레이션<br>말라리아와 빈혈의 치명률 변화가 후대로의 형질 전달에 영향을 미친다?

  - 지난 시간에 학습한 낫 모양 적혈구 빈혈증과 말라리아의 관계(영향력)에 따라 낫 모양 적혈구 빈혈증 유전자가<br>자손에게 전달되는 경향성을 알아보자.
  - 빈혈에 걸린 사람이 생존 확률, 말라리아에 걸린 사람의 사망 확률을 직접 설정해보며 변화를 관찰해보자.

### 다음 상황에 대해 생각해보자!
   
#### 만일, 낫 모양 적혈구 빈혈에 걸린 사람이 100% 확률로 결혼 전에 죽는다면?
  - 예상되는 해당 염기 서열이 자손에게 전달 될 확률은?

#### 만일, 낫 모양 적혈구 빈혈에 걸린 사람이 50% 확률로 10살, 50% 확률로 40살에 죽는다면?
  - 예상되는 해당 염기 서열이 자손에게 전달 될 확률은?

#### 위의 상황에서, 해당 지역에 말라리아 병원체를 가지는 모기가 900만마리가 들어왔다!
  - 말라리아에 걸렸을 때 100% 확률로 3일 뒤에 죽는다면, 예상되는 빈혈 환자의 염기 서열이 후손에게 전달 될 확률은?

### Biopython을 이용한 Malaria & Sickle cell Selection Pressure 실습

  - 학습 목표
    - 대립유전자 빈도 p, q의 이해
    - 세대 반복 시 유전자형의 빈도가 p² : 2pq : q² 로 수렴함을 확인
    - FASTA/VCF 등 실측 서열에서 대립유전자 빈도를 계산해 가설 검정
   
  - 왼쪽 메뉴의 [목차] → [+섹션] 클릭
  - 추가로 생성할 코드 셀의 제목 입력 (예시: 시뮬레이션: 세대별 변화 관찰)
  - 새롭게 섹션을 만들고 제목을 입력했다면, 그 부근에 마우스를 올려 [+코드] 클릭

  - 코드 셀에 다음과 같이 입력

    ```python

    # simulate라는 함수를 정의하자! [def] 이후, 괄호 안에 4가지 파라미터를 넣는다. 
    # selection_sickle: SS(낫모양 동형접합)의 생존율 (ex:0.9 → 정상보다 10% 덜 산다.)
    # malaria_mortality: AA(정상 동형접합)의 사망률 (ex:0.2 → 말라리아에 걸린 환자 중 20%가 죽는다.)
    # generations: 시뮬레이션을 진행 할 세대 수
    # initial_S_allele: 최초 S 대립유전자의 빈도(ex:0.1 → 10%)
    
    def simulate(selection_sickle=생존율, malaria_mortality=사망률, generations=세대, initial_S_allele=최초값):
    p = 1 - initial_S_allele; q = initial_S_allele    # def를 이용해 simulate라는 이름의 함수를 생성, 4가지 입력값을 받는다고 설정
    AA=AS=SS=None
    AA_lst, AS_lst, SS_lst = [], [], []               # AA, AS, SS는 비율(숫자 값)이 들어올 '자리'이며, 이를 설정
    for _ in range(generations):                      # 위와 마찬가지로 세대 수만큼 반복, '_'는 반복 횟수 자체가 필요 없을 때 사용
        f_AA, f_AS, f_SS = p*p, 2*p*q, q*q            # p*p는 AA 개체의 빈도, 2*p*q는 AS 개체 빈도, q*q는 SS 개체 빈도
        w_AA, w_AS, w_SS = 1-malaria_mortality, 1, 1-selection_sickle  # 각 유전형의 생존율 지정
        # AA는 말라리아에 취약 (정상 적혈구를 가진 사람)
        # AS는 정상과 변이 유전자가 섞인 사람이지만, 말라리아 저항성을 얻음과 동시에 빈혈 증상은 거의 없다고 가정
        # SS는 낫 모양 적혈구 빈혈증 환자로, 심한 빈혈과 조직 손상 발생 및 말라리아 면역
    
        mean_w = f_AA*w_AA + f_AS*w_AS + f_SS*w_SS  # 평균(mean) 생존율 계산
                                                    # 한 세대 전체에서 "살아남을 확률의 합"을 이용
    
        f_AA, f_AS, f_SS = f_AA*w_AA/mean_w, f_AS*w_AS/mean_w, f_SS*w_SS/mean_w
        # f_(유전자형)을 각각의 유전형 빈도 × 그 유전형의 생존율 ÷ 전체 평균 생존율로 정의
        # 다음 세대의 유전형 빈도 합을 1로 맞추기 위한 작업  
    
        AA_lst.append(f_AA); AS_lst.append(f_AS); SS_lst.append(f_SS) 
        # 계산하여 정의한 다음 세대의 유전형 빈도를 리스트에 순서대로 담는 명령어
    
        p = f_AA + 0.5*f_AS; q = 1 - p
        # 다음 세대의 대립유전자 빈도 결정
    
    plt.figure(figsize=(5,3))
    plt.plot(AA_lst,label="AA"); plt.plot(AS_lst,label="AS"); plt.plot(SS_lst,label="SS")
    plt.xlabel("generation"); plt.ylabel("frequency"); plt.legend(); plt.title("malaria selection pressure")
    plt.show()
    interact(simulate, selection_sickle=(0.5,1.0,0.05), malaria_mortality=(0,0.5,0.05),
        generations=(10,100,10), initial_S_allele=(0.01,0.3,0.01));

    # 그래프 그리기: plt.figure(figsize=(5,3)) → 그림의 크기 지정
    # plt.plot(AA_lst,...): AA 빈도, AS 빈도, SS 빈도를 세대별로 연결
    # plt.xlabel, plt.ylabel, plt.title : 축 이름과 제목 설정
    # plt.legend(): 그래프 안에 'AA', 'AS', 'SS' 라벨 표시
    # plt.show(): 그래프를 화면에 출력하라는 명령
    # interact(...): Google Colab에서 사용할 수 있는 기능으로, 슬라이더를 움직여 parameter를 변경할 수 있게 하는 도구
    # (0.5,1.0,0.05) 등은 슬라이더의 최솟값, 최댓값, 단계를 의미

    ```

### 마무리

  - 1. 시뮬레이션 그래프에서 AS 곡선이 안정적으로 유지되는 생물학적 의미를 해석하고, 친구들과 이야기 해봅시다.
