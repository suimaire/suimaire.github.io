// node --test tests/
import test from 'node:test';
import assert from 'node:assert/strict';
import {
  normalizePageKey,
  resolveMode,
  isWithinDedupeWindow,
  formatCounts,
  loadCounts,
  callCounterRpc,
  PAGE_VIEWS_CONFIG,
} from '../assets/js/page-views.js';

test('normalizePageKey: 요구사항 예시', () => {
  assert.equal(normalizePageKey('https://suimaire.github.io/'), '/');
  assert.equal(normalizePageKey('https://suimaire.github.io'), '/');
  assert.equal(normalizePageKey('https://suimaire.github.io/bioinformatics'), '/bioinformatics/');
  assert.equal(normalizePageKey('https://suimaire.github.io/bioinformatics/'), '/bioinformatics/');
  assert.equal(normalizePageKey('https://suimaire.github.io/bioinformatics/?x=123#foo'), '/bioinformatics/');
  assert.equal(normalizePageKey('https://suimaire.github.io/lipid-3d-explorer/?foo=1#section2'), '/lipid-3d-explorer/');
});

test('normalizePageKey: Jekyll .html, index.html, 중복 슬래시, 대소문자', () => {
  assert.equal(normalizePageKey('/lectures/day1.html'), '/lectures/day1/');
  assert.equal(normalizePageKey('/lectures/day1'), '/lectures/day1/');
  assert.equal(normalizePageKey('/carbohydrate-3d-explorer/index.html'), '/carbohydrate-3d-explorer/');
  assert.equal(normalizePageKey('/index.html'), '/');
  assert.equal(normalizePageKey('//predator-prey-simulation-2//'), '/predator-prey-simulation-2/');
  assert.equal(normalizePageKey('/Protein-3D-Explorer/'), '/protein-3d-explorer/');
  assert.equal(normalizePageKey('/한글/'), '/%ed%95%9c%ea%b8%80/');
});

test('normalizePageKey: 비정상 입력은 null (서버 형식과 동일)', () => {
  assert.equal(normalizePageKey('/.hidden/'), null);
  assert.equal(normalizePageKey('/' + 'a'.repeat(250)), null);
  assert.equal(normalizePageKey('/a"b/'), '/a%22b/'); // 인코딩되어 안전한 문자만 남음
});

test('normalizePageKey: 따옴표/태그는 인코딩되어 안전한 문자만 남음', () => {
  const key = normalizePageKey('/<script>alert(1)</script>/');
  assert.ok(key === null || /^[a-z0-9._~%\/-]+$/.test(key));
});

test('resolveMode: production 에서만 record', () => {
  const loc = (href) => new URL(href);
  assert.equal(resolveMode(loc('https://suimaire.github.io/')), 'record');
  assert.equal(resolveMode(loc('https://suimaire.github.io/x/?page-views=off')), 'disabled');
  assert.equal(resolveMode(loc('http://suimaire.github.io/')), 'disabled');
  assert.equal(resolveMode(loc('http://localhost:5173/predator-prey-simulation/')), 'disabled');
  assert.equal(resolveMode(loc('http://127.0.0.1:4000/')), 'disabled');
  assert.equal(resolveMode(loc('file:///D:/site/index.html')), 'disabled');
  assert.equal(resolveMode(loc('file:///D:/site/index.html?page-views=live')), 'disabled');
  assert.equal(resolveMode(loc('http://localhost:5173/?page-views=mock')), 'mock');
  assert.equal(resolveMode(loc('http://localhost:5173/?page-views=live')), 'read');
  assert.equal(resolveMode(loc('https://suimaire-github-io.example.com/')), 'disabled');
});

test('isWithinDedupeWindow: 30분', () => {
  const now = 10_000_000;
  const win = PAGE_VIEWS_CONFIG.dedupeWindowMs;
  assert.equal(win, 30 * 60 * 1000);
  assert.equal(isWithinDedupeWindow(NaN, now), false);
  assert.equal(isWithinDedupeWindow(now, now), true);
  assert.equal(isWithinDedupeWindow(now - win + 1, now), true);
  assert.equal(isWithinDedupeWindow(now - win, now), false);
  assert.equal(isWithinDedupeWindow(now + 60_000, now), false);
});

test('formatCounts: ko-KR 천 단위', () => {
  assert.equal(formatCounts(37, 1284), '오늘 조회 37 · 전체 조회 1,284');
  assert.equal(formatCounts(0, 1234567), '오늘 조회 0 · 전체 조회 1,234,567');
});

function memoryStorage() {
  const map = new Map();
  return {
    getItem: (k) => (map.has(k) ? map.get(k) : null),
    setItem: (k, v) => map.set(k, String(v)),
    removeItem: (k) => map.delete(k),
    map,
  };
}

function fakeServer() {
  const counts = new Map();
  const calls = [];
  const rpc = async (name, key) => {
    calls.push(name);
    const c = counts.get(key) || { today: 0, total: 0 };
    if (name === 'record_page_view') { c.today += 1; c.total += 1; counts.set(key, c); }
    return { ...c };
  };
  return { rpc, calls };
}

test('loadCounts: 첫 방문 +1, 30분 이내 새로고침은 읽기만, 30분 후 다시 +1, 다른 페이지는 별도', async () => {
  const storage = memoryStorage();
  const server = fakeServer();
  const t0 = 1_700_000_000_000;
  const run = (pageKey, now) => loadCounts({ pageKey, mode: 'record', storage, now, rpc: server.rpc });

  assert.deepEqual(await run('/a/', t0), { today: 1, total: 1, counted: true });
  assert.deepEqual(await run('/a/', t0 + 1000), { today: 1, total: 1, counted: false });
  assert.deepEqual(await run('/a/', t0 + 29 * 60_000), { today: 1, total: 1, counted: false });
  assert.deepEqual(await run('/b/', t0 + 2000), { today: 1, total: 1, counted: true });
  assert.deepEqual(await run('/a/', t0 + 30 * 60_000), { today: 2, total: 2, counted: true });
  assert.deepEqual(server.calls, ['record_page_view', 'get_page_view_counts', 'get_page_view_counts', 'record_page_view', 'record_page_view']);
});

test('loadCounts: 서버 HTTP 오류면 dedupe 표시를 되돌림, 네트워크 timeout 이면 유지', async () => {
  const storage = memoryStorage();
  const httpFail = async () => { const e = new Error('500'); e.httpStatus = 500; throw e; };
  await assert.rejects(loadCounts({ pageKey: '/a/', mode: 'record', storage, now: 1, rpc: httpFail }));
  assert.equal(storage.getItem(PAGE_VIEWS_CONFIG.storagePrefix + '/a/'), null);

  const netFail = async () => { throw new Error('aborted'); };
  await assert.rejects(loadCounts({ pageKey: '/a/', mode: 'record', storage, now: 5, rpc: netFail }));
  assert.equal(storage.getItem(PAGE_VIEWS_CONFIG.storagePrefix + '/a/'), '5');
});

test('loadCounts: disabled/read/mock 은 절대 record 하지 않음', async () => {
  const server = fakeServer();
  await assert.rejects(loadCounts({ pageKey: '/a/', mode: 'disabled', storage: memoryStorage(), rpc: server.rpc }));
  await loadCounts({ pageKey: '/a/', mode: 'read', storage: memoryStorage(), rpc: server.rpc });
  await loadCounts({ pageKey: '/a/', mode: 'mock', storage: memoryStorage(), rpc: server.rpc });
  assert.deepEqual(server.calls, ['get_page_view_counts']);
});

test('loadCounts: storage 사용 불가여도 동작', async () => {
  const server = fakeServer();
  const broken = { getItem() { throw new Error('denied'); }, setItem() { throw new Error('denied'); }, removeItem() {} };
  assert.deepEqual(await loadCounts({ pageKey: '/a/', mode: 'record', storage: broken, now: 1, rpc: server.rpc }), { today: 1, total: 1, counted: true });
});

test('callCounterRpc: 요청 형식, publishable key 만 사용, 응답 파싱', async () => {
  let seen;
  const fetchImpl = async (url, init) => {
    seen = { url, init };
    return new Response(JSON.stringify([{ kst_date: '2026-09-13', today_views: 37, total_views: 1284 }]), { status: 200 });
  };
  assert.deepEqual(await callCounterRpc('record_page_view', '/x/', { fetchImpl }), { today: 37, total: 1284 });
  assert.equal(seen.url, 'https://szbmpsvxzrnewyzqiokr.supabase.co/rest/v1/rpc/record_page_view');
  assert.equal(seen.init.method, 'POST');
  assert.equal(seen.init.body, '{"p_page_key":"/x/"}');
  assert.match(seen.init.headers.apikey, /^sb_publishable_/);
  assert.equal(seen.init.headers.Authorization, undefined);
});

test('callCounterRpc: HTTP 오류와 timeout 은 reject', async () => {
  await assert.rejects(
    callCounterRpc('get_page_view_counts', '/x/', { fetchImpl: async () => new Response('{}', { status: 404 }) }),
    (e) => e.httpStatus === 404,
  );
  const hang = (url, init) => new Promise((_, reject) => init.signal.addEventListener('abort', () => reject(new Error('aborted'))));
  await assert.rejects(callCounterRpc('get_page_view_counts', '/x/', { fetchImpl: hang, timeoutMs: 20 }), /aborted/);
});

test('소스에 service_role/secret key 가 없음', async () => {
  const { readFile } = await import('node:fs/promises');
  const source = await readFile(new URL('../assets/js/page-views.js', import.meta.url), 'utf8');
  assert.doesNotMatch(source, /sb_secret_|service_role"|eyJhbGciOi/);
});
