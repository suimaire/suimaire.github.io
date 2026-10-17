import { FEEDBACK_CONFIG } from "./feedback-config.js";

export const feedbackReady = Boolean(
  /^https:\/\/[a-z0-9-]+\.supabase\.co\/functions\/v1\/feedback-api$/.test(FEEDBACK_CONFIG.functionUrl) &&
  /^https:\/\/[a-z0-9-]+\.supabase\.co$/.test(FEEDBACK_CONFIG.supabaseUrl) &&
  FEEDBACK_CONFIG.publishableKey.startsWith("sb_publishable_") &&
  FEEDBACK_CONFIG.turnstileSiteKey
);

export function feedbackClientId() {
  const key = "hafs-feedback-client-id";
  try {
    let value = localStorage.getItem(key);
    if (!value) {
      value = crypto.randomUUID();
      localStorage.setItem(key, value);
    }
    return value;
  } catch {
    if (!window.__feedbackTemporaryClientId) window.__feedbackTemporaryClientId = crypto.randomUUID();
    return window.__feedbackTemporaryClientId;
  }
}

async function responseJson(res) {
  let data = {};
  try { data = await res.json(); } catch { /* generic failure below */ }
  if (!res.ok) throw new Error(data.error || "서버 요청에 실패했습니다.");
  return data;
}

export async function feedbackApi(action, { params = {}, body = null, token = null } = {}) {
  if (!feedbackReady) throw new Error("의견 게시판 연결을 준비하고 있습니다.");
  const url = new URL(FEEDBACK_CONFIG.functionUrl);
  url.searchParams.set("action", action);
  for (const [key, value] of Object.entries(params)) url.searchParams.set(key, String(value));

  const headers = { apikey: FEEDBACK_CONFIG.publishableKey };
  if (body !== null) headers["Content-Type"] = "application/json";
  if (token) headers.Authorization = "Bearer " + token;
  const res = await fetch(url, {
    method: body === null ? "GET" : "POST",
    headers, body: body === null ? undefined : JSON.stringify(body),
    cache: "no-store",
    credentials: "omit",
  });
  return responseJson(res);
}

let turnstilePromise = null;
function ensureTurnstileScript() {
  if (window.turnstile) return Promise.resolve(window.turnstile);
  if (!turnstilePromise) {
    turnstilePromise = new Promise((resolve, reject) => {
      const script = document.createElement("script");
      script.src = "https://challenges.cloudflare.com/turnstile/v0/api.js?render=explicit";
      script.async = true;
      script.onload = () => window.turnstile ? resolve(window.turnstile) : reject(new Error("보안 인증을 불러오지 못했습니다."));
      script.onerror = () => reject(new Error("보안 인증 서비스에 연결하지 못했습니다."));
      document.head.append(script);
    });
  }
  return turnstilePromise;
}

export async function initFeedbackCaptcha(mount, onToken) {
  if (!feedbackReady) throw new Error("보안 인증 설정을 준비 중입니다.");
  const turnstile = await ensureTurnstileScript();
  const id = turnstile.render(mount, {
    sitekey: FEEDBACK_CONFIG.turnstileSiteKey,
    theme: "auto",
    callback: token => onToken(token),
    "expired-callback": () => onToken(""),
    "error-callback": () => onToken(""),
  });
  return () => {
    onToken("");
    try { turnstile.reset(id); } catch { /* widget may already have been removed */ }
  };
}

export async function teacherLogin(email, password) {
  if (!feedbackReady) throw new Error("의견 게시판 연결을 준비하고 있습니다.");
  const response = await fetch(FEEDBACK_CONFIG.supabaseUrl + "/auth/v1/token?grant_type=password", {
    method: "POST",
    headers: {
      apikey: FEEDBACK_CONFIG.publishableKey,
      "Content-Type": "application/json",
    },
    body: JSON.stringify({ email, password }),
    cache: "no-store",
    credentials: "omit",
  });
  return responseJson(response);
}
