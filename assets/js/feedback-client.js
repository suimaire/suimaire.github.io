import { FEEDBACK_CONFIG } from "./feedback-config.js";

export const feedbackReady = Boolean(
  /^https:\/\/[a-z0-9-]+\.supabase\.co\/functions\/v1\/feedback-api$/.test(FEEDBACK_CONFIG.functionUrl) &&
  /^https:\/\/[a-z0-9-]+\.supabase\.co$/.test(FEEDBACK_CONFIG.supabaseUrl) &&
  FEEDBACK_CONFIG.publishableKey.startsWith("sb_publishable_")
);

// The server retains its own independent write-enable secret. Never rely
// on this UI flag for actual authorization or write protection.
export const feedbackPostingEnabled = feedbackReady && FEEDBACK_CONFIG.postingEnabled === true;

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
