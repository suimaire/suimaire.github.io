/*
 * HAFS 과학 수업 포털 공통 조회수 표시 ("today 37 · total 1,284")
 *
 * 이 파일 하나를 포털(Jekyll)과 모든 학습 앱(별도 GitHub Pages 저장소)이 함께 씁니다.
 *   포털  : _includes/footer_custom.html 이 로드
 *   각 앱 : index.html 의 짧은 loader 가 https://suimaire.github.io/assets/js/page-views.js 를 로드
 * 표시 위치는 페이지의 [data-page-views] 요소입니다. 없으면 body 맨 끝에 작은 줄을 만듭니다.
 *
 * DB 스키마: suimaire/predator-prey-simulation-2 저장소 supabase/migrations/20260913_page_views.sql
 *
 * 집계 규칙
 *   - page_key = URL pathname 정규화 (query/hash 제외, index.html·.html 제거, 끝 '/' 통일, 소문자)
 *   - https://suimaire.github.io 에서만 증가. localhost·127.0.0.1·file:// 등은 절대 증가시키지 않음
 *   - 같은 브라우저에서 같은 page_key 를 30분 이내 다시 열면 증가 없이 현재 값만 읽음
 *   - 오늘/전체 계산과 KST 날짜 판정은 전부 DB 함수가 수행
 *
 * 개발 환경 표시 (production 이 아닌 곳)
 *   기본       : "조회수 · 개발 환경에서는 집계하지 않음"
 *   ?page-views=mock : 가짜 숫자로 UI 확인
 *   ?page-views=live : 실제 숫자를 읽기 전용으로 표시(증가 없음)
 *
 * 이 모듈의 어떤 실패도 페이지로 전파되지 않습니다. 모든 비동기 작업은 catch 되고,
 * 네트워크 요청에는 timeout 이 있으며, 실패하면 표시 줄을 조용히 숨깁니다.
 */

export const PAGE_VIEWS_CONFIG = Object.freeze({
  supabaseUrl: 'https://szbmpsvxzrnewyzqiokr.supabase.co',
  // publishable(anon) key — 브라우저 공개용. 권한은 DB 의 GRANT/RLS 가 제한합니다.
  // service_role / secret key 는 절대 여기에 넣지 마세요.
  publishableKey: 'sb_publishable_bOM16W0j9Dh6U72M_fNp8Q_HwRZm9tF',
  productionHost: 'suimaire.github.io',
  dedupeWindowMs: 30 * 60 * 1000,
  timeoutMs: 6000,
  storagePrefix: 'hafs-page-views:v1:',
});

const MAX_KEY_LENGTH = 200;
const KEY_PATTERN = /^\/([a-z0-9_~%-][a-z0-9._~%-]*\/)*$/;

/**
 * URL(또는 pathname) → canonical page key. 형식에 맞지 않으면 null.
 *   https://suimaire.github.io/                         → /
 *   https://suimaire.github.io/bioinformatics           → /bioinformatics/
 *   https://suimaire.github.io/bioinformatics/?x=1#foo  → /bioinformatics/
 *   /lectures/day1.html                                 → /lectures/day1/
 *   /carbohydrate-3d-explorer/index.html                → /carbohydrate-3d-explorer/
 */
export function normalizePageKey(input) {
  let path;
  try {
    const raw = String(input);
    const origin = 'https://' + PAGE_VIEWS_CONFIG.productionHost;
    // '/...' 는 항상 경로로 취급합니다('//x' 가 호스트로 해석되지 않도록).
    path = new URL(raw.startsWith('/') ? origin + raw : raw, origin).pathname;
  } catch {
    return null;
  }
  path = path.toLowerCase().replace(/\/{2,}/g, '/');
  path = path.replace(/\/index\.html?$/, '/').replace(/\.html?$/, '');
  if (!path.endsWith('/')) path += '/';
  if (path.length > MAX_KEY_LENGTH || !KEY_PATTERN.test(path)) return null;
  return path;
}

/** 'record' | 'read' | 'mock' | 'disabled' */
export function resolveMode(location) {
  const param = safeSearchParam(location, 'page-views');
  if (location.protocol === 'https:' && location.hostname === PAGE_VIEWS_CONFIG.productionHost) {
    return param === 'off' ? 'disabled' : 'record';
  }
  if (param === 'mock') return 'mock';
  if (param === 'live' && /^https?:$/.test(location.protocol)) return 'read';
  return 'disabled';
}

function safeSearchParam(location, name) {
  try {
    return new URLSearchParams(location.search || '').get(name);
  } catch {
    return null;
  }
}

/** 마지막 집계 시각이 dedupe 창 안이면 true. 시계가 거꾸로 간 경우(미래 시각)는 창 밖으로 봅니다. */
export function isWithinDedupeWindow(lastCountedAt, now, windowMs = PAGE_VIEWS_CONFIG.dedupeWindowMs) {
  if (!Number.isFinite(lastCountedAt)) return false;
  const age = now - lastCountedAt;
  return age >= 0 && age < windowMs;
}

export function formatCounts(today, total) {
  const f = (value) => Number(value).toLocaleString('ko-KR');
  return `today ${f(today)} · total ${f(total)}`;
}

function openStorage() {
  for (const name of ['localStorage', 'sessionStorage']) {
    try {
      const store = globalThis[name];
      const probe = PAGE_VIEWS_CONFIG.storagePrefix + 'probe';
      store.setItem(probe, '1');
      store.removeItem(probe);
      return store;
    } catch {
      /* 사용 불가 — 다음 후보 */
    }
  }
  return null;
}

/**
 * Supabase RPC 호출. REST fetch 만 사용합니다(SDK 불필요).
 * @returns {Promise<{today:number,total:number}>}
 */
export async function callCounterRpc(name, pageKey, { fetchImpl = globalThis.fetch, timeoutMs = PAGE_VIEWS_CONFIG.timeoutMs } = {}) {
  const controller = typeof AbortController === 'function' ? new AbortController() : null;
  const timer = controller ? setTimeout(() => controller.abort(), timeoutMs) : null;
  try {
    const response = await fetchImpl(`${PAGE_VIEWS_CONFIG.supabaseUrl}/rest/v1/rpc/${name}`, {
      method: 'POST',
      headers: {
        apikey: PAGE_VIEWS_CONFIG.publishableKey,
        'Content-Type': 'application/json',
        Accept: 'application/json',
      },
      body: JSON.stringify({ p_page_key: pageKey }),
      signal: controller ? controller.signal : undefined,
      keepalive: name === 'record_page_view',
      credentials: 'omit',
    });
    if (!response.ok) {
      const error = new Error(`page views ${name} HTTP ${response.status}`);
      error.httpStatus = response.status;
      throw error;
    }
    const rows = await response.json();
    const row = Array.isArray(rows) ? rows[0] : rows;
    const today = Number(row && row.today_views);
    const total = Number(row && row.total_views);
    if (!Number.isFinite(today) || !Number.isFinite(total)) throw new Error('page views: unexpected response');
    return { today, total };
  } finally {
    if (timer) clearTimeout(timer);
  }
}

/**
 * 30분 dedupe 를 적용해 조회수를 기록하거나 읽습니다.
 * 기록 요청을 보내기 직전에 시각을 저장하므로 요청 중 새로고침해도 두 번 세지 않습니다.
 * 서버가 명시적으로 거절(HTTP 오류)한 경우에만 이전 값으로 되돌려 다음 방문에 다시 시도합니다.
 */
export async function loadCounts({ pageKey, mode, storage = openStorage(), now = Date.now(), rpc = callCounterRpc }) {
  if (mode === 'mock') return { today: 24, total: 1392, counted: false };
  if (mode === 'read') return { ...(await rpc('get_page_view_counts', pageKey)), counted: false };
  if (mode !== 'record') throw new Error('page views disabled');

  const storageKey = PAGE_VIEWS_CONFIG.storagePrefix + pageKey;
  let previous = null;
  try {
    previous = storage ? storage.getItem(storageKey) : null;
  } catch {
    previous = null;
  }
  if (previous !== null && isWithinDedupeWindow(Number.parseInt(previous, 10), now)) {
    return { ...(await rpc('get_page_view_counts', pageKey)), counted: false };
  }

  try {
    if (storage) storage.setItem(storageKey, String(now));
  } catch {
    /* 저장 실패해도 집계는 진행 */
  }
  try {
    return { ...(await rpc('record_page_view', pageKey)), counted: true };
  } catch (error) {
    if (error && error.httpStatus && storage) {
      try {
        if (previous === null) storage.removeItem(storageKey);
        else storage.setItem(storageKey, previous);
      } catch {
        /* 무시 */
      }
    }
    throw error;
  }
}

const STYLE_ID = 'hafs-page-views-style';
const STYLE = `
.page-views{font-size:12px;line-height:1.5;color:var(--page-views-color,#56626b);font-variant-numeric:tabular-nums;letter-spacing:0;overflow-wrap:anywhere}
.page-views--standalone{display:block;margin:0;padding:12px 16px 20px;text-align:center}
`;

function ensureStyle(doc) {
  if (doc.getElementById(STYLE_ID)) return;
  const style = doc.createElement('style');
  style.id = STYLE_ID;
  style.textContent = STYLE;
  (doc.head || doc.documentElement).appendChild(style);
}

function waitForMount(doc, selector, timeoutMs) {
  const found = doc.querySelector(selector);
  if (found) return Promise.resolve(found);
  return new Promise((resolve) => {
    let observer = null;
    const finish = (element) => {
      if (observer) observer.disconnect();
      clearTimeout(timer);
      resolve(element);
    };
    const timer = setTimeout(() => finish(doc.querySelector(selector)), timeoutMs);
    if (typeof MutationObserver === 'function') {
      observer = new MutationObserver(() => {
        const element = doc.querySelector(selector);
        if (element) finish(element);
      });
      observer.observe(doc.documentElement, { childList: true, subtree: true });
    }
  });
}

function createStandaloneMount(doc) {
  const element = doc.createElement('p');
  element.className = 'page-views--standalone';
  element.setAttribute('data-page-views', '');
  doc.body.appendChild(element);
  return element;
}

/**
 * 공통 진입점. 반환 Promise 는 절대 reject 되지 않습니다.
 * @param {{pageKey?: string, mount?: Element|string, location?: Location, document?: Document}} [options]
 */
export async function initPageViewCounter(options = {}) {
  try {
    const doc = options.document || document;
    const loc = options.location || window.location;
    const mode = options.mode || resolveMode(loc);
    const pageKey = normalizePageKey(options.pageKey || loc.pathname);

    ensureStyle(doc);
    // 네트워크 요청은 mount 를 기다리지 않고 바로 시작합니다.
    const countsPromise = pageKey && mode !== 'disabled'
      ? loadCounts({ pageKey, mode })
      : Promise.reject(new Error(mode === 'disabled' ? 'page views disabled' : 'invalid page key'));
    countsPromise.catch(() => {}); // mount 를 기다리는 동안 unhandled rejection 방지

    let element = typeof options.mount === 'object' && options.mount ? options.mount : null;
    if (!element) {
      const selector = typeof options.mount === 'string' ? options.mount : '[data-page-views]';
      element = await waitForMount(doc, selector, 10000);
    }
    if (!element) element = createStandaloneMount(doc);
    if (element.getAttribute('data-page-views-mounted') === '1') return;
    element.setAttribute('data-page-views-mounted', '1');
    element.classList.add('page-views');

    if (mode === 'disabled') {
      element.textContent = '조회수 · 개발 환경에서는 집계하지 않음';
      element.hidden = false;
      return;
    }

    element.textContent = '조회수 불러오는 중…';
    element.hidden = false;
    try {
      const { today, total } = await countsPromise;
      element.textContent = formatCounts(today, total) + (mode === 'mock' ? ' (mock)' : mode === 'read' ? ' (읽기 전용)' : '');
    } catch (error) {
      element.hidden = true;
      element.textContent = '';
      console.warn('[page-views] 조회수를 불러오지 못했습니다.', error && error.message ? error.message : error);
    }
  } catch (error) {
    try {
      console.warn('[page-views] 초기화 실패', error);
    } catch {
      /* 무시 */
    }
  }
}

// <script type="module" src=".../page-views.js"> 로 로드되면 자동 실행합니다(한 페이지에 한 번만).
if (typeof window !== 'undefined' && typeof document !== 'undefined' && !window.__hafsPageViewsStarted) {
  window.__hafsPageViewsStarted = true;
  const start = () => { initPageViewCounter(); };
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', start, { once: true });
  else start();
}
