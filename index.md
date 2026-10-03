---
layout: default
title: 과학 수업 포털
nav_order: 1
description: 질문하고, 관찰하고, 데이터로 설명하는 과학 수업 공간
og_title: HAFS BIOLOGY LAB
og_description: EXPLORE · MODEL · ANALYZE · EXPLAIN
og_image: /assets/images/hafs-biology-lab-og.png
og_image_alt: HAFS BIOLOGY LAB
---

{%- comment -%}
  포털 홈. 학생이 "내 수업"을 먼저 고르고, 그 수업의 자료를 여는 구조다.
  스타일은 assets/css/portal.css 의 "메인 포털" 부분 참고.

  자료 번호 — 수업 시간에 "B4 열어 보세요"처럼 부르기 위한 번호
    수업마다 글자 하나, 자료마다 그 뒤에 1, 2, 3 … 을 붙인다.
      A  #ecosystem      통합과학2      .course--eco  (초록)
      B  #molecular      기초 생화학    .course--mol  (파랑)
      C  #bioinformatics 생물정보학     .course--bio  (주황)
    · 번호는 화면의 배지(.res__code)와 칸의 id 가 같다.
      suimaire.github.io/#B4 로 열면 그 자료로 바로 내려가 테두리가 강조된다.
    · 이미 수업에서 쓴 번호가 바뀌지 않도록, 새 자료는 각 수업의 맨 뒤에 추가한다.
    · 수업 글자는 사이드바(_config.yml nav_external_links 의 code)와 맞춘다.
    · PCR 학습지(bioinformatics/pcr-primer-design.html) 상단에도 "C2" 가 적혀 있다.
  #ecosystem 등 영역 id 는 사이드바와 PCR 학습지의 "← 포털" 링크가 쓴다.

  자료 추가 위치
    B 분자·생화학 탐구 → _data/molecular_explorers.yml 맨 아래에만 항목 추가 (번호 자동)
    A, C               → 아래 해당 <ul class="res-grid"> 맨 뒤에 <li class="res"> 추가하고
                         번호(id, .res__code)와 맨 위 "수업 고르기"의 번호 범위·개수를 고친다.

  자료 한 칸(.res)은 제목 링크 하나가 칸 전체를 덮는다. 따로 버튼 링크를 두지 않고
  .res__go 는 눈에 보이는 안내 문구일 뿐이다(스크린리더에는 번호와 제목이 한 번 읽힌다).

  heading-numbers.js 가 붙이는 1.1 / 1.2.1 번호는 이 페이지에서는 CSS 로 숨긴다.
{%- endcomment -%}

{%- assign explorer_count = site.data.molecular_explorers | size -%}
{%- assign prelearning_count = 0 -%}
{%- for item in site.data.molecular_explorers -%}
  {%- if item.prelearning -%}{%- assign prelearning_count = prelearning_count | plus: 1 -%}{%- endif -%}
{%- endfor -%}
{%- assign molecular_last = explorer_count | plus: prelearning_count -%}

<div class="science-portal">
  <header class="home-intro">
    <div class="home-intro__text">
      <p class="home-intro__brand">{{ site.title | escape }}</p>
      <p class="home-intro__credit">{{ site.author_credit | escape }}</p>
      <h1 id="portal-title">생명과학 수업 자료실</h1>
      <p class="home-intro__lead">
        수업 시간에 쓰는 시뮬레이션, 3D 분자 모델, 학습지를 수업별로 모아 두었습니다.
        지금 듣고 있는 수업을 고르면 바로 시작할 수 있어요.
      </p>
      <p class="home-intro__hint">
        모두 웹 브라우저에서 열리며 따로 설치할 프로그램은 없습니다.
        생물정보학 실습만 Google Colab을 사용합니다.
      </p>
    </div>
    {% include home/population-figure.html %}
  </header>

  <nav class="home-picker" aria-labelledby="home-picker-title">
    <div class="home-picker__head">
      <p class="home-picker__title" id="home-picker-title">어느 수업 자료를 찾나요?</p>
      <p class="home-picker__help">수업 시간에 들은 번호(예: <strong>B4</strong>)를 찾아 열어도 됩니다.</p>
    </div>
    <ul class="home-picker__list">
      <li>
        <a class="home-picker__item course--eco" href="#ecosystem">
          <span class="home-picker__letter" aria-hidden="true">A</span>
          <span class="home-picker__course">통합과학2</span>
          <span class="home-picker__topic">생태계 · 개체군 모델링</span>
          <span class="home-picker__count"><b>A1–A3</b> 학습지&nbsp;1 · 시뮬레이션&nbsp;2</span>
        </a>
      </li>
      <li>
        <a class="home-picker__item course--mol" href="#molecular">
          <span class="home-picker__letter" aria-hidden="true">B</span>
          <span class="home-picker__course">기초 생화학</span>
          <span class="home-picker__topic">분자 · 생화학 탐구</span>
          <span class="home-picker__count"><b>B1–B{{ molecular_last }}</b> 탐색기&nbsp;{{ explorer_count }}{% if prelearning_count > 0 %} · 사전&nbsp;학습&nbsp;{{ prelearning_count }}{% endif %}</span>
        </a>
      </li>
      <li>
        <a class="home-picker__item course--bio" href="#bioinformatics">
          <span class="home-picker__letter" aria-hidden="true">C</span>
          <span class="home-picker__course">생물정보학 특강</span>
          <span class="home-picker__topic">생물정보학 · 데이터</span>
          <span class="home-picker__count"><b>C1–C2</b> 5일&nbsp;강좌 · 웹&nbsp;학습지&nbsp;1</span>
        </a>
      </li>
    </ul>
  </nav>

  <div id="courses">

    <section class="course course--eco" id="ecosystem" aria-labelledby="ecosystem-module-title">
      <div class="course__head">
        <p class="course__tab"><span class="course__letter">A</span>통합과학2</p>
        <h2 id="ecosystem-module-title">생태계 · 개체군 모델링</h2>
        <p class="course__note">
          서로 다른 두 모형으로 생태계의 변화를 관찰하고, 그래프에 나타난 관계를 학습지에서 해석합니다.
        </p>
      </div>

      <div class="res res--lead" id="A1">
        <p class="res__kind">{% include home/icon.html name="worksheet" %}학습지 · 약 20분 <span class="res__flag">여기서 시작</span></p>
        <h3 class="res__title"><span class="res__code">A1</span> <a href="{{ '/predator-prey-worksheet/' | relative_url }}">포식자와 피식자, 누가 먼저 변할까?</a></h3>
        <p class="res__desc">
          두 시뮬레이션에서 관찰한 개체군 변화의 순서와 시간 지연을 그래프로 해석하고 하나의 설명으로 연결합니다.
        </p>
        <span class="res__go" aria-hidden="true">학습지 시작하기 <span class="arrow">→</span></span>
      </div>

      <p class="course__sub">학습지에서 함께 쓰는 시뮬레이션</p>
      <ul class="res-grid">
        <li class="res" id="A2">
          <p class="res__kind">{% include home/icon.html name="graph" %}시뮬레이션 · 개체군 동역학 모형</p>
          <h3 class="res__title"><span class="res__code">A2</span> <a href="{{ '/predator-prey-simulation/' | relative_url }}">포식자-피식자 동역학 실험실</a></h3>
          <p class="res__desc">
            개별 개체가 아니라 포식자와 피식자 개체군 전체의 시간에 따른 변화와 시간 지연을 그래프로 탐구합니다.
          </p>
          <span class="res__go" aria-hidden="true">시뮬레이션 열기 <span class="arrow">→</span></span>
        </li>
        <li class="res" id="A3">
          <p class="res__kind">{% include home/icon.html name="graph" %}시뮬레이션 · 개체 기반 모형</p>
          <h3 class="res__title"><span class="res__code">A3</span> <a href="{{ '/predator-prey-simulation-2/' | relative_url }}">토끼와 늑대 숲 생태계</a></h3>
          <p class="res__desc">
            개별 토끼와 늑대의 이동 · 먹이 · 번식 조건을 조절하고, 그 행동이 모여 전체 생태계 변화를 만드는 과정을 관찰합니다.
          </p>
          <span class="res__go" aria-hidden="true">시뮬레이션 열기 <span class="arrow">→</span></span>
        </li>
      </ul>
    </section>

    <section class="course course--mol" id="molecular" aria-labelledby="molecular-module-title">
      <div class="course__head">
        <p class="course__tab"><span class="course__letter">B</span>2026 HAFS Elective Track · Basic Biochemistry</p>
        <h2 id="molecular-module-title">분자 · 생화학 탐구</h2>
        <p class="course__note">
          생체분자의 구조와 화학적 성질을 3D 모델, 그래프, 시뮬레이션으로 직접 조작하며 분자 간 상호작용과 생화학적 기전을 탐구하는 수업자료입니다.
        </p>
      </div>

      <ol class="res-list">
        {%- assign n = 0 %}
        {%- for item in site.data.molecular_explorers %}
        {%- comment -%} 사전 학습은 연결된 탐색기 바로 앞 칸. 번호도 그 순서대로 붙는다 {%- endcomment -%}
        {%- if item.prelearning %}
        {%- assign n = n | plus: 1 %}
        {%- assign next = n | plus: 1 %}
        <li class="res res--prep" id="B{{ n }}">
          <p class="res__kind">{% include home/icon.html name="cards" %}{{ item.prelearning.label }} · B{{ next }} 전에 해 보세요</p>
          <h3 class="res__title"><span class="res__code">B{{ n }}</span> <a href="{{ item.prelearning.url }}">{{ item.prelearning.title }}</a></h3>
          <p class="res__desc">{{ item.prelearning.description }}</p>
          <p class="res__tags">{{ item.prelearning.meta }}</p>
          <span class="res__go" aria-hidden="true">{{ item.prelearning.cta }} <span class="arrow">→</span></span>
        </li>
        {%- endif %}
        {%- assign n = n | plus: 1 %}
        <li class="res" id="B{{ n }}">
          <p class="res__kind">{% include home/icon.html name="molecule" %}{{ item.label }}{% if item.status and item.status != "운영 중" %} <span class="res__flag">{{ item.status }}</span>{% endif %}</p>
          <h3 class="res__title"><span class="res__code">B{{ n }}</span> <a href="{{ item.url }}">{{ item.title }}</a></h3>
          <p class="res__desc">{{ item.description }}</p>
          {%- if item.meta %}
          <p class="res__tags">{{ item.meta | join: " · " }}</p>
          {%- endif %}
          <span class="res__go" aria-hidden="true">{{ item.cta }} <span class="arrow">→</span></span>
        </li>
        {%- endfor %}
      </ol>
    </section>

    <section class="course course--bio" id="bioinformatics" aria-labelledby="bioinformatics-module-title">
      <div class="course__head">
        <p class="course__tab"><span class="course__letter">C</span>2025 HSHS 생물정보학 특강</p>
        <h2 id="bioinformatics-module-title">생물정보학 · 데이터</h2>
        <p class="course__note">
          공개 생명과학 데이터와 분석 도구를 활용해 생명 현상을 탐구하는 수업자료입니다.
        </p>
      </div>

      <ul class="res-grid">
        <li class="res" id="C1">
          <p class="res__kind">{% include home/icon.html name="book" %}강좌 · 5일 과정</p>
          <h3 class="res__title"><span class="res__code">C1</span> <a href="{{ '/bioinformatics/' | relative_url }}">생물정보학 기초</a></h3>
          <p class="res__desc">
            Biopython과 공개 생명과학 데이터를 활용해 서열, 단백질 구조, 유전 평형과 변이를 탐구합니다.
          </p>
          <p class="res__tags">Google Colab · 탐구 활동 중심</p>
          <span class="res__go" aria-hidden="true">강좌 안내 보기 <span class="arrow">→</span></span>
        </li>
        <li class="res" id="C2">
          <p class="res__kind">{% include home/icon.html name="worksheet" %}웹 학습지</p>
          <h3 class="res__title"><span class="res__code">C2</span> <a href="{{ '/bioinformatics/pcr-primer-design/' | relative_url }}">PCR과 프라이머 디자인</a></h3>
          <p class="res__desc">
            프라이머의 결합 위치와 방향을 바꾸며 PCR 산물을 예측하고, 설계한 프라이머의 적합성을 검토합니다.
          </p>
          <span class="res__go" aria-hidden="true">학습지 열기 <span class="arrow">→</span></span>
        </li>
      </ul>
    </section>

  </div>
</div>
