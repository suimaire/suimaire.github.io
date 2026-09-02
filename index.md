---
layout: default
title: 과학 수업 포털
nav_order: 1
description: 질문하고, 관찰하고, 데이터로 설명하는 과학 수업 공간
---

<div class="science-portal">
  <section class="portal-hero portal-hero--compact" aria-labelledby="portal-title">
    <div class="portal-hero__copy">
      <p class="portal-kicker">HAFS SCIENCE LEARNING LAB</p>
      <h1 id="portal-title">데이터로 확인하는 과학</h1>
      <p class="portal-hero__lead">
        관찰하고, 비교하고, 데이터로 설명하는 과학 수업 공간
      </p>
    </div>
    <a class="portal-button portal-button--primary" href="#courses">수업자료 보기 <span aria-hidden="true">↓</span></a>
  </section>

  <section class="portal-section" id="courses" aria-labelledby="courses-title">
    <div class="portal-section__heading">
      <p class="portal-eyebrow">COURSE LIBRARY</p>
      <h2 id="courses-title">수업자료(공개 중)</h2>
    </div>

    <section class="ecosystem-module" aria-labelledby="ecosystem-module-title">
      <div class="ecosystem-module__heading">
        <div>
          <p class="portal-eyebrow">01 · ECOSYSTEM &amp; POPULATION MODELING</p>
          <h3 id="ecosystem-module-title">생태계와 개체군 모델링</h3>
          <p>서로 다른 두 모형으로 생태계의 변화를 관찰하고, 그래프에 나타난 관계를 학습지에서 해석합니다.</p>
        </div>
        <span class="ecosystem-module__sequence">모형 탐구 <b aria-hidden="true">→</b> 자료 해석</span>
      </div>

      <div class="portal-course-grid ecosystem-model-grid" aria-label="생태계와 개체군 모델링 시뮬레이션">
        <article class="course-card course-card--coming course-card--available ecosystem-card">
          <div class="course-card__topline">
            <span class="ecosystem-card__number">01 · SIMULATION</span>
            <span class="course-card__status">운영 중</span>
          </div>
          <span class="ecosystem-card__type">개체 기반 모형</span>
          <h3>토끼와 늑대 숲 생태계</h3>
          <p>개별 토끼와 늑대의 이동·먹이·번식 조건을 조절하고, 그 행동이 모여 전체 생태계 변화를 만드는 과정을 관찰합니다.</p>
          <a class="course-card__link course-card__link--dark" href="{{ '/predator-prey-simulation-2/' | relative_url }}">시뮬레이션 바로가기 <span aria-hidden="true">→</span></a>
        </article>

        <article class="course-card course-card--coming course-card--available ecosystem-card">
          <div class="course-card__topline">
            <span class="ecosystem-card__number">02 · SIMULATION</span>
            <span class="course-card__status">운영 중</span>
          </div>
          <span class="ecosystem-card__type">개체군 동역학 모형</span>
          <h3>포식자-피식자 동역학 실험실</h3>
          <p>개별 개체가 아니라 포식자와 피식자 개체군 전체의 시간에 따른 변화와 시간 지연을 그래프로 탐구합니다.</p>
          <a class="course-card__link course-card__link--dark" href="{{ '/predator-prey-simulation/' | relative_url }}">시뮬레이션 바로가기 <span aria-hidden="true">→</span></a>
        </article>
      </div>

      <div class="ecosystem-connector" aria-hidden="true">
        <span></span><b>두 모형에서 관찰한 패턴을 설명으로 연결</b><span></span>
      </div>

      <article class="course-card ecosystem-followup">
        <div class="ecosystem-followup__copy">
          <div class="course-card__topline">
            <span class="ecosystem-card__number">03 · FOLLOW-UP ACTIVITY</span>
            <span class="ecosystem-followup__time">약 20분</span>
          </div>
          <span class="ecosystem-card__type ecosystem-card__type--light">자료 해석 활동</span>
          <h3>포식자와 피식자, 누가 먼저 변할까?</h3>
          <p>두 시뮬레이션에서 관찰한 개체군 변화의 순서와 시간 지연을 그래프로 해석하고 하나의 설명으로 연결합니다.</p>
        </div>
        <a class="ecosystem-followup__link" href="{{ '/predator-prey-worksheet/' | relative_url }}">학습지 시작하기 <span aria-hidden="true">→</span></a>
      </article>
    </section>

    <section class="bioinformatics-module" aria-labelledby="bioinformatics-module-title">
      <div class="bioinformatics-module__heading">
        <p class="portal-eyebrow">02 · BIOINFORMATICS &amp; DATA</p>
        <h3 id="bioinformatics-module-title">생물정보학과 데이터</h3>
        <p>공개 생명과학 데이터와 분석 도구를 활용해 생명 현상을 탐구하는 수업자료입니다.</p>
      </div>

      <article class="course-card course-card--featured">
        <div class="course-card__topline">
          <span class="course-card__subject">생명과학 · 데이터</span>
          <span class="course-card__status">운영 중</span>
        </div>
        <h3>생물정보학 기초</h3>
        <p>Biopython과 공개 생명과학 데이터를 활용해 서열, 단백질 구조, 유전 평형과 변이를 탐구합니다.</p>
        <div class="course-card__meta" aria-label="강좌 정보">
          <span>5일 과정</span><span>Google Colab</span><span>탐구 활동 중심</span>
        </div>
        <a class="course-card__link" href="{{ '/bioinformatics/' | relative_url }}">강좌 안내 보기 <span aria-hidden="true">→</span></a>
      </article>
    </section>
  </section>
</div>
