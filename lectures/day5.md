---
layout: default
title: "Day 5"
nav_order: 6
---

# Day 5 COVID-19 변이 실습 및 응용

---
## 기본 개념
  - 렉틴(lectin): 탄수화물과 강한 친화력과 높은 특이성을 가지고 결합할 수 있는 단백질들의 총칭, 모든 생물체에 존재
  - 렉틴은 세포가 주변 세포들을 인식할 수 있도록 하는 기능 및 신호 전달과 같은 기능을 수행한다.

### 바이러스의 감염에 관여하는 헤마글루티닌(hemagglutinin)
  - 독감 바이러스(influenza virus)를 포함하는 여러 가지 동물 바이러스들은 숙주세포의 세포막에 존재하는 탄수화물과 상호작용하여 숙주세포에 부착된다.
  - 해당 탄수화물의 끝 부분(말단)에는 **시알산(cialic acid)** 라고 불리는 물질이 존재한다.
  - 시알산을 특이적으로 분해하는 효소의 이름은 **뉴라미니데이스(neuraminidase)** 이다.
  - 이 때 숙주세포의 탄수화물(올리고당; 단당류가 여러 개 모인 것)과 결합하는 독감바이러스의 렉틴(lectin)은 **헤마글루티닌**이라고 한다.
  - 인간 인플루엔자 바이러스의 외피 단백질에는 H_N_ (_은 숫자)로 표시하는 외피 단백질들이 존재한다.
  - H는 헤마글루티닌, N은 뉴라미니데이스를 의미한다. 이들은 인플루엔자 바이러스의 가장 바깥쪽 표면에 돌출되어 있기 때문에 스파이크(spike)라고 부르기도 하며 바이러스 감염과 증식에 주요 역할을 한다.
  - 헤마글루티닌은 숙주 세포(감염시킬 목표 개체의 세포) 표면의 시알산과 특이적 결합을 함으로써 감염을 일으킨다.
  - 이러한 결합은 바이러스가 세포 내에서 증식한 후, 합체되어 다시 빠져나올 때에도 일어난다.
  - 만일, 바이러스의 헤마글루티닌이 숙주 세포막의 시알산과 결합한 채로 있으면 바이러스가 숙주 세포로부터 분리되지 못한다.
  - 이 때 인플루엔자 바이러스의 **뉴라미니데이스** 가 헤마글루티닌과 결합한 시알산을 절단함으로써 바이러스가 최종적으로 세포로부터 분리될 수 있게 해 준다.
  - 항바이러스제인 **타미플루(Oseltamivir)** 는 뉴라미니데이스 수용체에 결합하여 해당 효소의 작용을 억제한다.
  - 이는 감염된 세포로부터 바이러스 방출을 못하게 막음으로써 감염이 전파 및 확산되지 못하게 하는 기전을 통해 바이러스에 대한 치료효과를 보인다.

### 계속 변이되는 바이러스
  - [델타 우세종 된 뒤 백신 효과 66%로 감소](https://www.youtube.com/watch?v=Fk-ebDSKOns)
  - [오미크론 재감염 위험, 델타보다 2.4배 높아](https://www.youtube.com/watch?v=yXr7o2RjNb8)
  - [오미크론, 델타보다 감염력 월등... 3차 접종 필요](https://www.youtube.com/watch?v=_jHxR5PkGH0)

### Google Colab Notebook 열기
  - [GOOGLE COLAB](https://colab.research.google.com)
  - 왼쪽 위 [파일] → Drive의 새 노트북

### 스파이크 단백질의 실제 서열 변이를 Biopython으로 확인해 보자.
  - 코드 셀에 다음과 같이 입력

```python
!pip install biopython		
from Bio import Entrez, SeqIO		

Entrez.email = "your.email@example.com"

def fetch_spike(acc):
    handle = Entrez.efetch(db="nucleotide", id=acc,
    rettype="fasta", retmode="text")
    return SeqIO.read(handle, "fasta").seq


wuhan = fetch_spike("NC_045512.2") # Wuhan 전체 게놈에서 S 유전자 서열만 추출
delta = fetch_spike("OL955326.1") # 델타 변이 스파이크
omicron= fetch_spike("OM283822.1") # 오미크론 변이 스파이크
```

```python
""" 서열 정렬 & 변이 위치 찾기 """

from Bio import pairwise2
from Bio.pairwise2 import format_alignment

# Wuhan vs Delta
aln1 = pairwise2.align.globalms(
    wuhan, delta, 2, -1, -0.5, -0.1)[0]
print(format_alignment(*aln1))

# Wuhan vs Omicron 실습
# 학생들이 ACC 넣어봤다가 출력해 보기
```

### RBD(319-541aa) 영역에서 변이 확인

```python
import matplotlib.pyplot as plt

variants = ["Wuhan","Delta","Omicron"]
folds    = [1.0, 3.3, 20.0]  # 예시: 중화능 저하 배수

plt.figure()
plt.bar(variants, folds)
plt.ylabel("Neutralization titer fold reduction")
plt.title("Variant vs Neutralization Drop")
plt.show()
```

  - bar chart로 얼마나 중화능이 떨어졌는지 비교
  - RBD 변이가 많을수록 백신 효능이 저하되는 이유에 대해서 생각해봅시다.
