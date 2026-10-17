---
layout: default
title: 학생 의견 게시판
description: HAFS Biology Lab 학생 의견과 개선 제안을 공유하는 공개 게시판
nav_exclude: true
---

<link rel="stylesheet" href="{{ '/assets/css/feedback.css' | relative_url }}">

<section class="feedback-board" data-feedback-board aria-labelledby="feedback-title">
  <div class="feedback-board__intro">
    <a class="feedback-back" href="{{ '/' | relative_url }}">← 수업 자료실로 돌아가기</a>
    <h1 id="feedback-title">학생 의견 게시판</h1>
    <p class="feedback-board__help">수업 자료를 더 좋은 방향으로 바꾸는 아이디어를 함께 나누는 공간입니다. 다른 학생의 의견도 읽을 수 있습니다.</p>
  </div>
  <p class="feedback-board__rules">
    작성한 의견은 승인 없이 즉시 공개됩니다. 익명 또는 이름·닉네임을 선택할 수 있지만, 이름은 신원 인증을 의미하지 않습니다.
    타인의 이름·학번·연락처·사진 등 개인정보를 게시하지 마세요.
  </p>
  <div class="feedback-board__controls">
    <strong>최근 의견</strong>
    <button type="button" class="feedback-button feedback-button--primary" data-feedback-compose>＋ 의견 남기기</button>
  </div>
  <p data-feedback-notice class="feedback-notice" role="status" aria-live="polite" hidden></p>
  <ul class="feedback-list feedback-board__list" data-feedback-list aria-label="공개 의견 목록"></ul>
  <button type="button" class="feedback-button feedback-board__more" data-feedback-more hidden>이전 의견 더 보기</button>
  <p class="feedback-board__admin"><a href="{{ '/feedback/admin/' | relative_url }}">교사 관리</a></p>
</section>

<script type="module" src="{{ '/assets/js/feedback.js' | relative_url }}"></script>
