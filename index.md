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
      <h1 id="portal-title">데이터로 확인하는 과학</span></h1>
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
      <h2 id="courses-title">지금 만날 수 있는 수업</h2>
      <p>기존 생물정보학 자료는 원래 주소 그대로 보존되어 있습니다.</p>
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

    <div class="portal-course-grid" aria-label="준비 중인 수업 영역">
      <article class="course-card course-card--coming">
        <div class="course-card__topline">
          <span class="course-card__subject">생태 · 모델링</span>
          <span class="course-card__status course-card__status--muted">준비 중</span>
        </div>
        <h3>생태계 시뮬레이션</h3>
        <p>개체군의 상호작용을 조절하고 시간에 따른 변화를 해석하는 탐구 수업을 준비하고 있습니다.</p>
      </article>

      <article class="course-card course-card--coming">
        <div class="course-card__topline">
          <span class="course-card__subject">과학 · 탐구</span>
          <span class="course-card__status course-card__status--muted">확장 예정</span>
        </div>
        <h3>다음 질문을 위한 자리</h3>
        <p>새로운 수업을 같은 방식으로 모을 수 있도록 포털의 구조를 열어 두었습니다.</p>
      </article>
    </div>
  </section>

  <section class="portal-section portal-method" aria-labelledby="method-title">
    <div class="portal-section__heading">
      <p class="portal-eyebrow">HOW WE LEARN</p>
      <h2 id="method-title">답보다 먼저, 좋은 질문을 만듭니다</h2>
    </div>

    <ol class="method-list">
      <li><strong>01</strong><div><h3>질문하기</h3><p>현상에서 궁금한 점을 찾고 검증 가능한 질문으로 다듬습니다.</p></div></li>
      <li><strong>02</strong><div><h3>탐구하기</h3><p>관찰, 실험, 공개 데이터와 모델을 이용해 근거를 모읍니다.</p></div></li>
      <li><strong>03</strong><div><h3>설명하기</h3><p>결과를 비교하고 한계까지 살펴보며 나만의 설명을 완성합니다.</p></div></li>
    </ol>
  </section>

  <aside class="portal-preservation" aria-labelledby="preservation-title">
    <div>
      <p class="portal-eyebrow">CONTENT PRESERVATION</p>
      <h2 id="preservation-title">기존 수업 자료를 안전하게 보존하고 있습니다</h2>
      <p>현재 생물정보학 Day 1–5 문서와 이미지 주소는 바꾸지 않았습니다. 포털 초안은 그 위에 새로운 입구만 더한 상태입니다.</p>
    </div>
    <a href="{{ '/lectures/day1.html' | relative_url }}">기존 Day 1 열기 <span aria-hidden="true">→</span></a>
  </aside>
</div>
