---
layout: default
title: 의견 게시판 관리
description: HAFS Biology Lab 의견 게시판 교사 관리 페이지
nav_exclude: true
---

<link rel="stylesheet" href="{{ '/assets/css/feedback.css' | relative_url }}">

<section class="feedback-admin" data-feedback-admin aria-labelledby="feedback-admin-title">
  <a class="feedback-back" href="{{ '/feedback/' | relative_url }}">← 의견 게시판으로 돌아가기</a>
  <h1 id="feedback-admin-title">의견 게시판 관리</h1>
  <p class="feedback-board__help">교사 계정으로 로그인하면 게시글 숨김, 답변, 상태 변경과 삭제를 관리할 수 있습니다.</p>

  <form data-admin-login>
    <label class="feedback-field">교사 이메일<input name="email" type="email" autocomplete="username" required></label>
    <label class="feedback-field">비밀번호<input name="password" type="password" autocomplete="current-password" required></label>
    <button type="submit" class="feedback-button feedback-button--primary">교사 로그인</button>
    <p class="feedback-notice" role="status" aria-live="polite" data-admin-auth-notice hidden></p>
  </form>
  <div data-admin-panel hidden>
    <div class="feedback-admin__head">
      <strong>게시글 관리 · 최근 50건</strong>
      <div class="feedback-actions">
        <button type="button" class="feedback-button" data-admin-refresh>새로고침</button>
        <button type="button" class="feedback-button feedback-button--secondary" data-admin-logout>로그아웃</button>
      </div>
    </div>
    <p class="feedback-notice" role="status" aria-live="polite" data-admin-notice hidden></p>
    <div class="feedback-admin__list" data-admin-list></div>
    <button type="button" class="feedback-button feedback-board__more" data-admin-more hidden>이전 의견 더 보기</button>
  </div>
</section>

<script type="module" src="{{ '/assets/js/feedback-admin.js' | relative_url }}"></script>
