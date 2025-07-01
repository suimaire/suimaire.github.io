document.addEventListener('DOMContentLoaded', () => {
  // 페이지 메인 콘텐츠 영역 선택
  const container = document.querySelector('.page__content') || document.body;
  if (!container) return;

  // 0–5 레벨용 카운터(0번은 안 씁니다)
  const counters = new Array(6).fill(0);

  // h1~h5 모두 선택
  const headings = container.querySelectorAll('h1, h2, h3, h4, h5');

  headings.forEach(h => {
    // h1→1, h2→2...
    const lvl = parseInt(h.tagName.substr(1), 10);

    // 해당 레벨 카운터 1 증가
    counters[lvl] += 1;
    // 하위 레벨 카운터 초기화
    for (let i = lvl + 1; i < counters.length; i++) {
      counters[i] = 0;
    }

    // 1번 인덱스부터 현재 레벨까지 번호 문자열 생성
    // h1 → slice(1,2) → [1] → "1"
    // h2 → slice(1,3) → [1,1] → "1.1"
    const num = counters.slice(1, lvl + 1).join('.');

    // 번호 span 삽입
    const span = document.createElement('span');
    span.className = 'heading-number';
    span.textContent = num + ' ';
    h.prepend(span);
  });
});
