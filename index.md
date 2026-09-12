---
layout: default
title: 과학 수업 포털
nav_order: 1
description: 질문하고, 관찰하고, 데이터로 설명하는 과학 수업 공간
---

{%- comment -%}
  포털 홈. 자료 목록은 "카드"가 아니라 편집된 목록으로 조판한다.
  스타일은 assets/css/portal.css 참고.

  번호(1.1, 1.2.1 …)는 문서에 적지 않는다. assets/js/heading-numbers.js 가
  h1 → 1, h2 → 1.1, h3 → 1.1.1 순서대로 자동으로 붙인다. 따라서
    · 절은 반드시 h2, 개별 자료는 반드시 h3
    · 자료를 넣고 빼면 번호는 알아서 다시 매겨진다.

  자료 추가 위치
    1.2 분자·생화학 3D 탐구 → _data/molecular_explorers.yml 에만 항목 추가
    1.1 / 1.3             → 아래 해당 <ol class="portal-resources"> 안에 <li> 추가
{%- endcomment -%}

<div class="science-portal">
  <header class="portal-masthead">
    <p class="portal-masthead__org">HAFS Science Learning Lab</p>
    <h1 id="portal-title">데이터로 확인하는 과학</h1>
    <p class="portal-masthead__lead">
      생명과학 수업에서 사용하는 시뮬레이션, 3D 분자 모델, 데이터 탐구 자료를 모았습니다.
      관찰하고, 비교하고, 데이터로 설명하는 과학 수업 공간입니다.
    </p>
  </header>

  <nav class="portal-contents" aria-label="자료 영역 바로가기">
    <a href="#ecosystem"><span class="portal-contents__num">1.1</span>생태계·개체군 모델링</a>
    <a href="#molecular"><span class="portal-contents__num">1.2</span>분자·생화학 3D 탐구</a>
    <a href="#bioinformatics"><span class="portal-contents__num">1.3</span>생물정보학·데이터</a>
  </nav>

  <div id="courses">

    <section class="portal-part" id="ecosystem" aria-labelledby="ecosystem-module-title">
      <h2 id="ecosystem-module-title">생태계·개체군 모델링</h2>
      <p class="portal-part__note">
        서로 다른 두 모형으로 생태계의 변화를 관찰하고, 그래프에 나타난 관계를 학습지에서 해석합니다.
      </p>

      <ol class="portal-resources">
        <li class="portal-resource portal-resource--compact">
          <h3><a href="{{ '/predator-prey-worksheet/' | relative_url }}">포식자와 피식자, 누가 먼저 변할까?</a></h3>
          <p class="portal-resource__meta">
            <span class="portal-resource__kind">자료 해석 학습지</span>
            <span class="portal-resource__status">약 20분</span>
          </p>
          <p class="portal-resource__desc">
            두 시뮬레이션에서 관찰한 개체군 변화의 순서와 시간 지연을 그래프로 해석하고 하나의 설명으로 연결합니다.
          </p>
          <p class="portal-resource__action">
            <a href="{{ '/predator-prey-worksheet/' | relative_url }}">학습지 시작하기 <span class="arrow" aria-hidden="true">→</span></a>
          </p>
        </li>

        <li class="portal-resource">
          <h3><a href="{{ '/predator-prey-simulation/' | relative_url }}">포식자-피식자 동역학 실험실</a></h3>
          <p class="portal-resource__meta">
            <span class="portal-resource__kind">개체군 동역학 모형</span>
            <span class="portal-resource__status">운영 중</span>
          </p>
          <p class="portal-resource__desc">
            개별 개체가 아니라 포식자와 피식자 개체군 전체의 시간에 따른 변화와 시간 지연을 그래프로 탐구합니다.
          </p>
          <p class="portal-resource__action">
            <a href="{{ '/predator-prey-simulation/' | relative_url }}">시뮬레이션 열기 <span class="arrow" aria-hidden="true">→</span></a>
          </p>
        </li>

        <li class="portal-resource">
          <h3><a href="{{ '/predator-prey-simulation-2/' | relative_url }}">토끼와 늑대 숲 생태계</a></h3>
          <p class="portal-resource__meta">
            <span class="portal-resource__kind">개체 기반 모형</span>
            <span class="portal-resource__status">운영 중</span>
          </p>
          <p class="portal-resource__desc">
            개별 토끼와 늑대의 이동·먹이·번식 조건을 조절하고, 그 행동이 모여 전체 생태계 변화를 만드는 과정을 관찰합니다.
          </p>
          <p class="portal-resource__action">
            <a href="{{ '/predator-prey-simulation-2/' | relative_url }}">시뮬레이션 열기 <span class="arrow" aria-hidden="true">→</span></a>
          </p>
        </li>
      </ol>
    </section>

    <section class="portal-part" id="molecular" aria-labelledby="molecular-module-title">
      <h2 id="molecular-module-title">분자·생화학 3D 탐구</h2>
      <p class="portal-part__note">
        생체분자의 구조를 직접 회전하고, 분자 간 상호작용과 생화학적 기전을 3D 모델과 애니메이션으로 탐구하는 수업자료입니다.
      </p>

      <ol class="portal-resources">
        {%- for item in site.data.molecular_explorers %}
        <li class="portal-resource">
          <h3><a href="{{ item.url }}">{{ item.title }}</a></h3>
          <p class="portal-resource__meta">
            <span class="portal-resource__kind">{{ item.label }}</span>
            <span class="portal-resource__status">{{ item.status }}</span>
          </p>
          <p class="portal-resource__desc">{{ item.description }}</p>
          {%- if item.meta %}
          <p class="portal-resource__tags">{{ item.meta | join: " · " }}</p>
          {%- endif %}
          <p class="portal-resource__action">
            <a href="{{ item.url }}">{{ item.cta }} <span class="arrow" aria-hidden="true">→</span></a>
          </p>
        </li>
        {%- endfor %}
      </ol>
    </section>

    <section class="portal-part" id="bioinformatics" aria-labelledby="bioinformatics-module-title">
      <h2 id="bioinformatics-module-title">생물정보학·데이터</h2>
      <p class="portal-part__note">
        공개 생명과학 데이터와 분석 도구를 활용해 생명 현상을 탐구하는 수업자료입니다.
      </p>

      <ol class="portal-resources">
        <li class="portal-resource">
          <h3><a href="{{ '/bioinformatics/' | relative_url }}">생물정보학 기초</a></h3>
          <p class="portal-resource__meta">
            <span class="portal-resource__kind">5일 강좌</span>
            <span class="portal-resource__status">운영 중</span>
          </p>
          <p class="portal-resource__desc">
            Biopython과 공개 생명과학 데이터를 활용해 서열, 단백질 구조, 유전 평형과 변이를 탐구합니다.
          </p>
          <p class="portal-resource__tags">Google Colab · 탐구 활동 중심</p>
          <p class="portal-resource__action">
            <a href="{{ '/bioinformatics/' | relative_url }}">강좌 안내 보기 <span class="arrow" aria-hidden="true">→</span></a>
          </p>
        </li>
      </ol>
    </section>

  </div>
</div>
