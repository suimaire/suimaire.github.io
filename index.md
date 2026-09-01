---
layout: default
title: 과학 수업 포털
nav_order: 1
description: 질문하고, 관찰하고, 데이터로 설명하는 과학 수업 공간
---

<div class="science-portal">
  <section class="portal-hero" aria-labelledby="portal-title">
    <div class="portal-hero__copy">
      <p class="portal-kicker">HAFS SCIENCE LEARNING LAB</p>
      <h1 id="portal-title">데이터로 확인하는 과학</h1>
      <p class="portal-hero__lead">
        시각 자료를 관찰하고 분석하며 설명을 만들어 가는
        과학 수업 공간
      </p>
      <div class="portal-actions">
        <a class="portal-button portal-button--primary" href="#courses">수업 둘러보기</a>
        <a class="portal-button portal-button--quiet" href="{{ '/bioinformatics/' | relative_url }}">생물정보학 살펴보기</a>
      </div>
    </div>

    <div class="portal-hero__visual" aria-label="질문, 탐구, 설명으로 이어지는 수업 흐름">
      <div class="orbit orbit--one"></div>
      <div class="orbit orbit--two"></div>
      <div class="science-node science-node--question"><strong>01</strong><span>질문</span></div>
      <div class="science-node science-node--explore"><strong>02</strong><span>탐구</span></div>
      <div class="science-node science-node--explain"><strong>03</strong><span>설명</span></div>
      <div class="science-core"><span>SCIENCE</span><strong>LAB</strong></div>
    </div>
  </section>

  <section class="portal-section" id="courses" aria-labelledby="courses-title">
    <div class="portal-section__heading">
      <p class="portal-eyebrow">COURSE LIBRARY</p>
      <h2 id="courses-title">수업자료(공개 중)</h2>
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

    <div class="portal-course-grid" aria-label="공개 중인 시뮬레이션 수업">
      <article class="course-card course-card--coming course-card--available">
        <div class="course-card__topline">
          <span class="course-card__subject">생태 · 모델링</span>
          <span class="course-card__status">운영 중</span>
        </div>
        <h3>토끼와 늑대 숲 생태계</h3>
        <p>토끼와 늑대의 개체군 조건을 조절하고 시간에 따른 생태계 변화를 직접 관찰합니다.</p>
        <a class="course-card__link course-card__link--dark" href="{{ '/predator-prey-simulation-2/' | relative_url }}">시뮬레이션 바로가기 <span aria-hidden="true">→</span></a>
      </article>

      <article class="course-card course-card--coming course-card--available">
        <div class="course-card__topline">
          <span class="course-card__subject">생태 · 개체군</span>
          <span class="course-card__status">운영 중</span>
        </div>
        <h3>포식자-피식자 동역학 실험실</h3>
        <p>포식자와 피식자 개체군이 시간에 따라 어떻게 달라지는지 시뮬레이션으로 관찰합니다.</p>
        <a class="course-card__link course-card__link--dark" href="{{ '/predator-prey-simulation/' | relative_url }}">시뮬레이션 바로가기 <span aria-hidden="true">→</span></a>
      </article>
    </div>
  </section>
</div>
