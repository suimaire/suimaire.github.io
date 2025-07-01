document.addEventListener('DOMContentLoaded', () => {
  // 페이지 메인 콘텐츠 영역 선택 (just-the-docs 기준)
  const container = document.querySelector('.page__content') || document.body;
  if (!container) return;

  const counters = [];  // 각 heading 레벨별 카운터
  const headings = container.querySelectorAll('h1, h2, h3, h4, h5');

  headings.forEach(h => {
    const lvl = parseInt(h.tagName.substr(1));  // h1→1, h2→2...
    counters[lvl] = (counters[lvl] || 0) + 1;
    // 하위 레벨 카운터 초기화
    for (let i = lvl+1; i < counters.length; i++) counters[i] = 0;
    // "1.2.3" 형태 문자열 생성
    const num = counters.slice(0, lvl).join('.');
    // 번호 span 삽입
    const span = document.createElement('span');
    span.className = 'heading-number';
    span.textContent = num + ' ';
    h.prepend(span);
  });
});
