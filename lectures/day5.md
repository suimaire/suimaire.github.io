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

## 스파이크 단백질의 실제 서열 변이를 Biopython으로 확인해 보자.
  - 코드 셀에 다음과 같이 입력

```python
!pip install biopython py3Dmol		
from Bio import Entrez, SeqIO

!apt-get -y -qq install clustalo
```

```python
Entrez.email = "너네_이메일"

def fetch_spike(acc):
    handle = Entrez.efetch(db="nucleotide", id=acc,
    rettype="fasta", retmode="text")
    return SeqIO.read(handle, "fasta").seq


wuhan = fetch_spike("NC_045512.2") # Wuhan 전체 게놈에서 S 유전자 서열만 추출
delta = fetch_spike("OL955326.1") # 델타 변이 스파이크
omicron = fetch_spike("OM283822.1") # 오미크론 변이 스파이크
```

```python
""" 서열 정렬 & 변이 위치 찾기 """

from Bio import pairwise2
from Bio.pairwise2 import format_alignment

start_nt, end_nt = 21563-1, 25384
wuhan_spike = wuhan [start_nt:end_nt]
delta_spike = delta [start_nt:end_nt]
omicron_spike = omicron[start_nt:end_nt]

# Wuhan vs Delta
aln1 = pairwise2.align.globalms(
    wuhan_spike, delta_spike, 2, -1, -0.5, -0.1)[0]
print(format_alignment(*aln1))

```
### 생각해보기
- 추가로 변수명을 바꾸어 Wuhan vs Omicron 을 비교해 봅시다.
- 어떤 변이가 가장 치명적인 변이였나요? 무엇을 통해 알 수 있었나요?


### RBD(319-541aa) 영역에서 변이 확인

```python
import matplotlib.pyplot as plt

 # bar chart 만들어 보기
 # RBD : Spike 단백질에서 세포 수용체(ACE2)에 결합하는 부위로, 바이러스가 우리 세포에 침투할 때 중요한 부위
 # 따라서 백신 및 항체는 주로 RBD 부위를 겨냥하여 우리 몸을 보호한다.

variants = ["Wuhan","_____","_____"] # 원하는 변이가 있다면 추가로 입력해도 됩니다.
folds    = [1.0, ___ , ___ ]         # variants의 예시 - Delta , Omicron 등
"""
# 예시: 중화능 저하 배수
# ___에 원하는 변이 명 및 중화능 저하 배수를 스스로 찾아서 넣어봅시다.
# 포털 사이트에 'Pubmed' 검색 후 접속
# (비교하고 싶은 변이명) neutralization titer fold reduction 검색
# 논문을 읽으며 중화능 차이를 찾고 코드 내부에 입력
# 참고: 논문 본문 및 도표의 WA-1(또는 Wuhan) 대비 omicron의 배수 값을 읽어내면 됩니다.
"""
plt.figure()
plt.bar(variants, folds)
plt.ylabel("Neutralization titer fold reduction")
plt.title("Variant vs Neutralization Drop")
plt.show()
"""
중화능 1 : 같은 항원에 대한 면역반응을 일으키기 위해 같은 농의 항체가 필요하다
중화능 20 : 같은 항원에 대한 면역반응을 일으키기 위해 20배 농도의 항체가 필요하다

"""
```
####
- bar chart로 얼마나 중화능이 떨어졌는지 비교해 봅시다.
- RBD 변이가 많을수록 백신 효능이 저하되는 이유에 대해서 생각해봅시다.
- 실제 논문 및 보고서 등에서 다른 변이(BQ.1.1 등)의 fold 값을 찾아 리스트에 추가하고 그래프를 업데이트 해봅시다.

---

## 바이러스 표면 Spike py3Dmol 시각화

```python
import py3Dmol

# Omicron RBD 변이(예시)
omicron_resi = [339,371,373,375,417,440,446,484,493,498,501,505]

# RBD·항체 체인 매핑 (PDB ID ➜ 체인명)
# 오미크론 변이가 생긴 RBD 부분의 잔기 번호(residue, 아미노산 위치)를 리스트에 담아 둡니다.
# 이 위치들에 빨간 구(sphere)를 그려 강조할 거예요.

info = {
    "6W41": {  # Wuhan RBD–CR3022
        "rbd":  "A",
        "ab":   ["H","L"]
    },
    "7T9L": {  # Omicron RBD–CR3022
        "rbd":  "A",
        "ab":   ["H","L"]
    }
}

    # PDB ID(단백질 데이터베이스 식별자)별로
    # RBD가 속한 체인 이름(chain ID)과
    # 항체의 Heavy/light 체인 정보를 저장해 둡니다.

def draw_complex(pdb_id, title, highlight=False):
    """pdb_id (str) : 6W41 or 7T9L
       highlight=True : Omicron 변이 구 표시"""
    """
    [draw_complex]라는 이름의 함수를 정의합니다.

    호출할 때 pdb_id로 어떤 구조를 불러올지 지정하고,

    title으로 표시할 글자를 주고,

    highlight=True로 설정하면(오미크론일 때) 변이 부위를 강조합니다.
    """

    # 1) 뷰어 만들기
    view = py3Dmol.view(query=f"pdb:{pdb_id}", width=520, height=420)
    view.setBackgroundColor("white")
    """
    py3Dmol.view로 3D 시각화 창을 만들고,

    query="pdb:6W41"처럼 PDB 서버에서 구조 파일을 가져옵니다.

    창 크기(width×height)와 배경색을 설정합니다.
    """

    # 2) 모든 단백질 cartoon 기본 (연회색)
    view.setStyle({"protein":"true"}, {"cartoon":{"color":"gainsboro","opacity":1.0}})

    # 3) RBD 체인만 파란색으로 재도색 (굵기 증가)
    rbd_chain = info[pdb_id]["rbd"]
    view.addStyle(
        {"chain": rbd_chain},
        {"cartoon":{"color":"deepskyblue","thickness":1.0,"opacity":1.0}}
    )

    # 4) 항체 체인은 주황색으로 재도색
    for ch in info[pdb_id]["ab"]:
        view.addStyle(
            {"chain": ch},
            {"cartoon":{"color":"sandybrown","opacity":1.0}}
        )

    # 5) Omicron 변이 구 강조
    if highlight:
        for res in omicron_resi:
            view.addStyle(
                {"chain": rbd_chain, "resi": res},
                {"sphere":{"color":"red","radius":0.5}}
            )
    """
    highlight가 True이면(오미크론 그릴 때)

    RBD 체인의 지정된 잔기(resi) 위치마다

    **빨간 구(sphere)**를 반지름 0.5 크기로 그려 변이 부위를 표시합니다.  
    """

    # 6) 구조 전체가 들어오도록 줌
    view.zoomTo()

    # 7) 제목 레이블(좌상단)
    view.addLabel(
        title,
        {
          "fontColor":"white",
          "backgroundColor":"black",
          "fontSize":14,
          "inFront":True
        },
        {"screenOffset":{"x":10,"y":10}}
    )

    # 8) 결과값 그림 표출
    view.show()


# 실제로 그리라는 명령
draw_complex("6W41", "Wuhan RBD–CR3022")          # Wuhan
draw_complex("7T9L", "Omicron RBD–CR3022", True)  # Omicron (변이 구 표시)
"""
이렇게 draw_complex를 두 번 부르면

첫 번째는 우한 구조만 기본 색으로,

두 번째는 오미크론 변이 부위를 빨간 구로 강조해 그려 줍니다.
"""
```

### Spike의 모양이 변하는 것은 면역 반응에 어떤 영향을 미칠까요?

## Biopython으로 COVID-19 변이 계통수(tree) 제작

- 목표: 실제 Wuhan → Alpha → Delta → Omicron 유전자를 불러와 계통수를 통해 비교
- 탐구 과정
  1) 변이별 Genbank acc.ID를 확인 후, NCBI에서 이를 불러옵시다.
  2) NCBI Nucleotide 검색 - 각 변이의 GenBank ID 복사
- 코드 셀 추가, 다음과 같은 코드 작성

```python
# 1) Clustal Omega 실행 ★최소 안정 옵션만 남김
cmd = """
clustalo -i spike.fasta -o aligned.fasta
         --seqtype=DNA
         --guidetree-out tree.dnd
         --force
""".split()
import subprocess, textwrap, shlex
subprocess.run(cmd, check=True)


# 2) 변이 GenBank ID
seq_ids = {
    "Wuhan"  : "NC_045512.2",
    "Alpha"  : "OK091006",
    "Delta"  : "OM061695",
    "Omicron": "OL672836"
}

# 3) Spike 유전자 좌표 (nt 21563~25384)만 다운
Entrez.email = "너네_이메일"
records = []
for name, acc in seq_ids.items():
    handle = Entrez.efetch(db="nucleotide", id=acc,
                           rettype="fasta", retmode="text",
                           seq_start=21563, seq_stop=25384)   # ← Spike
    rec = SeqIO.read(handle, "fasta")
    rec.id = name; rec.description = ""
    records.append(rec)
SeqIO.write(records, "spike.fasta", "fasta")

# 4) Clustal Omega 실행 (threads & DNA 모드 지정)
cmd = [
    "clustalo",
    "-i", "spike.fasta",
    "-o", "aligned.fasta",
    "--seqtype=DNA",          # DNA 모드
    "--guidetree-out", "tree.dnd",
    "--force"                 # 기존 파일 덮어쓰기
]
subprocess.run(cmd, check=True)

# 5) 트리 시각화
tree = Phylo.read("tree.dnd", "newick")
Phylo.draw(tree, label_colors={
    "Wuhan":"black", "Alpha":"blue",
    "Delta":"orange", "Omicron":"red"
})
```

### 생각해보기
- 방금 작성한 코드 셀을 이용하여, BA.5 / XBB.1.5 FASTA를 추가한 뒤 다시 셀을 실행해봅시다.
- 트리를 캡쳐한 뒤 최종 활동지에 업로드 합시다. (구글 클래스룸)
- 5번의 수업으로 느낀 점을 활동지에 작성합시다. (구글 클래스룸)

```python
"""
** 고생 많았습니다!! **
"""
```
