import "jsr:@supabase/functions-js/edge-runtime.d.ts";
import { createClient } from "npm:@supabase/supabase-js@2.95.0";

// Public anonymous API. Deploy with verify_jwt=false: admin JWTs are verified
// explicitly below; all DB operations use a server-only secret key.
const ALLOWED_ORIGINS = new Set([
  "https://suimaire.github.io",
  "http://localhost:4000",
  "http://127.0.0.1:4000",
]);
const CATEGORIES = new Set(["bug", "content", "idea", "question", "other"]);
const REASONS = new Set(["spam", "abuse", "privacy", "other"]);
const STATUSES = new Set(["received", "reviewing", "done"]);
const UUID = /^[0-9a-f]{8}-[0-9a-f]{4}-[1-8][0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}$/i;
const CODE = /^(portal|[ABC][1-9][0-9]?)$/;
const POST_FIELDS = "id,author_name,is_anonymous,category,resource_code,title,body,status,created_at,updated_at";
const encoder = new TextEncoder();

const supabaseUrl = Deno.env.get("SUPABASE_URL") || "";
const secretKeys = JSON.parse(Deno.env.get("SUPABASE_SECRET_KEYS") || "{}");
const serviceKey = secretKeys.default || Deno.env.get("SUPABASE_SERVICE_ROLE_KEY") || "";
if (!supabaseUrl || !serviceKey) throw new Error("Supabase server credentials missing");
const db = createClient(supabaseUrl, serviceKey, {
  auth: { autoRefreshToken: false, persistSession: false },
});

class HttpError extends Error {
  constructor(public status: number, message: string) { super(message); }
}
function json(req: Request, data: unknown, status = 200): Response {
  const origin = req.headers.get("origin") || "";
  const headers = new Headers({
    "Content-Type": "application/json; charset=utf-8",
    "Cache-Control": "no-store",
    "Vary": "Origin",
    "X-Content-Type-Options": "nosniff",
    "Access-Control-Allow-Methods": "GET,POST,OPTIONS",
    "Access-Control-Allow-Headers": "apikey,authorization,content-type",
  });
  if (ALLOWED_ORIGINS.has(origin)) headers.set("Access-Control-Allow-Origin", origin);
  return new Response(JSON.stringify(data), { status, headers });
}
const clean = (value: unknown) => typeof value === "string" ? value.trim() : "";
function sized(value: unknown, min: number, max: number, label: string): string {
  const v = clean(value);
  if (v.length < min || v.length > max) throw new HttpError(400, label + "의 길이를 확인해 주세요.");
  return v;
}
function uuid(value: unknown): string {
  const id = clean(value);
  if (!UUID.test(id)) throw new HttpError(400, "올바르지 않은 게시글 번호입니다.");
  return id;
}
async function sha256(value: string): Promise<string> {
  const bytes = await crypto.subtle.digest("SHA-256", encoder.encode(value));
  return Array.from(new Uint8Array(bytes), b => b.toString(16).padStart(2, "0")).join("");
}
function randomHex(length: number): string {
  const bytes = crypto.getRandomValues(new Uint8Array(length));
  return Array.from(bytes, b => b.toString(16).padStart(2, "0")).join("");
}
async function passwordHash(password: string, salt: string): Promise<string> {
  const key = await crypto.subtle.importKey("raw", encoder.encode(password), "PBKDF2", false, ["deriveBits"]);
  const digest = await crypto.subtle.deriveBits({
    name: "PBKDF2", salt: encoder.encode(salt), iterations: 180000, hash: "SHA-256",
  }, key, 256);
  return Array.from(new Uint8Array(digest), b => b.toString(16).padStart(2, "0")).join("");
}
function constantTimeEqual(a: string, b: string): boolean {
  if (a.length !== b.length) return false;
  let diff = 0;
  for (let i = 0; i < a.length; i++) diff |= a.charCodeAt(i) ^ b.charCodeAt(i);
  return diff === 0;
}
function requester(req: Request): string {
  // Read the rightmost hop where the gateway supplies x-forwarded-for.
  // This is best-effort abuse control, not identity verification.
  const forwarded = req.headers.get("x-forwarded-for") || "";
  const hop = forwarded.split(",").map(s => s.trim()).filter(Boolean).pop();
  return (hop || req.headers.get("cf-connecting-ip") || "unavailable").slice(0, 100);
}
async function consume(req: Request, action: string, clientId: string,
                       browserMax: number, windowSeconds: number): Promise<void> {
  const network = requester(req);
  const browser = /^[a-f0-9-]{32,40}$/i.test(clientId) ? clientId : "missing";
  // Per-browser quota prevents accidental repeated submissions. A much looser
  // network quota avoids blocking a whole classroom behind one public IP.
  const rules: Array<[string, string, number, number]> = [
    [action + "-browser", "browser:" + browser + ":" + network, browserMax, windowSeconds],
    [action + "-network-hour", "network:" + network, 120, 3600],
    [action + "-network-day", "network:" + network, 320, 86400],
  ];
  for (const [scope, key, max, seconds] of rules) {
    const digest = await sha256(serviceKey + ":" + key);
    const { data, error } = await db.rpc("feedback_consume_quota", {
      p_key_hash: digest, p_action: scope, p_window_seconds: seconds, p_max_hits: max,
    });
    if (error) throw error;
    if (data !== true) throw new HttpError(429, "요청이 너무 많습니다. 잠시 후 다시 시도해 주세요.");
  }
}
function rejectBotTrap(value: unknown): void {
  // Simple honeypot. This deters naive automated form submissions but not
  // determined scripts; the database-enforced quotas remain necessary.
  if (clean(value)) throw new HttpError(400, "올바르지 않은 입력입니다.");
}
async function readBody(req: Request): Promise<Record<string, unknown>> {
  if (!(req.headers.get("content-type") || "").startsWith("application/json"))
    throw new HttpError(415, "JSON 요청만 허용합니다.");
  const body = await req.text();
  if (body.length > 12000) throw new HttpError(413, "입력 크기가 너무 큽니다.");
  let parsed: unknown;
  try { parsed = JSON.parse(body); } catch { throw new HttpError(400, "입력 내용을 읽을 수 없습니다."); }
  if (!parsed || typeof parsed !== "object" || Array.isArray(parsed))
    throw new HttpError(400, "잘못된 입력 형식입니다.");
  return parsed as Record<string, unknown>;
}
async function moderator(req: Request): Promise<void> {
  const auth = req.headers.get("authorization") || "";
  const jwt = auth.match(/^Bearer (.+)$/i)?.[1];
  if (!jwt) throw new HttpError(401, "교사 로그인이 필요합니다.");
  const { data: { user }, error } = await db.auth.getUser(jwt);
  if (error || !user) throw new HttpError(401, "로그인이 만료되었습니다.");
  const { data, error: adminError } = await db.from("feedback_moderators")
    .select("user_id").eq("user_id", user.id).maybeSingle();
  if (adminError) throw adminError;
  if (!data) throw new HttpError(403, "관리자 권한이 없습니다.");
}
async function postsList(offset: number, admin: boolean): Promise<unknown> {
  const limit = admin ? 50 : 15;
  let query = db.from("feedback_posts").select(POST_FIELDS + (admin ? ",hidden" : ""))
    .order("created_at", { ascending: false })
    .order("id", { ascending: false })
    .range(offset, offset + limit);
  if (!admin) query = query.eq("hidden", false);
  const { data, error } = await query;
  if (error) throw error;
  const items = (data || []).slice(0, limit);
  const ids = items.map(p => p.id);
  const replies = ids.length ? await db.from("feedback_replies")
    .select("post_id,body,created_at,updated_at").in("post_id", ids) : { data: [], error: null };
  if (replies.error) throw replies.error;
  const replyMap = new Map((replies.data || []).map(r => [r.post_id, r]));
  const reports = admin && ids.length ? await db.from("feedback_reports").select("post_id,reason").in("post_id", ids) : { data: [], error: null };
  if (reports.error) throw reports.error;
  const reportCounts = new Map<string,number>();
  const reasons = new Map<string,string[]>();
  for (const item of reports.data || []) {
    reportCounts.set(item.post_id, (reportCounts.get(item.post_id) || 0) + 1);
    reasons.set(item.post_id, [...(reasons.get(item.post_id) || []), item.reason]);
  }
  return {
    posts: items.map(p => ({ ...p, reply: replyMap.get(p.id) || null,
      ...(admin ? { report_count: reportCounts.get(p.id) || 0, report_reasons: reasons.get(p.id) || [] } : {}) })),
    has_more: (data || []).length > limit,
  };
}
async function visiblePost(id: string): Promise<unknown> {
  const { data, error } = await db.from("feedback_posts").select(POST_FIELDS)
    .eq("id", id).eq("hidden", false).maybeSingle();
  if (error) throw error;
  if (!data) throw new HttpError(404, "게시글을 찾을 수 없습니다.");
  const { data: reply, error: replyError } = await db.from("feedback_replies")
    .select("body,created_at,updated_at").eq("post_id", id).maybeSingle();
  if (replyError) throw replyError;
  return { post: { ...data, reply: reply || null } };
}
async function handle(req: Request): Promise<Response> {
  const url = new URL(req.url);
  const action = url.searchParams.get("action") || "";
  if (req.method === "GET" && action === "list") {
    const offset = Math.min(10000, Math.max(0, Number.parseInt(url.searchParams.get("offset") || "0") || 0));
    return json(req, await postsList(offset, false));
  }
  if (req.method === "GET" && action === "detail")
    return json(req, await visiblePost(uuid(url.searchParams.get("id"))));
  if (req.method !== "POST") throw new HttpError(405, "허용되지 않은 요청입니다.");

  const b = await readBody(req);
  const clientId = clean(b.client_id);
  if (action === "create") {
    if (Deno.env.get("FEEDBACK_PUBLIC_POSTING_ENABLED") === "false")
      throw new HttpError(503, "새 의견 작성이 일시 중지되었습니다.");
    rejectBotTrap(b.website);
    await consume(req, "create", clientId, 3, 600);
    const anonymous = b.is_anonymous !== false;
    const author = anonymous ? "익명" : sized(b.author_name, 1, 24, "작성자");
    const category = clean(b.category);
    const resource = clean(b.resource_code || "portal");
    if (!CATEGORIES.has(category) || !CODE.test(resource))
      throw new HttpError(400, "의견 분류 또는 자료 번호를 확인해 주세요.");
    const title = sized(b.title, 2, 80, "제목");
    const body = sized(b.body, 5, 1500, "내용");
    const password = sized(b.password, 6, 64, "게시글 비밀번호");
    const salt = randomHex(16);
    const hash = await passwordHash(password, salt);
    const { data: post, error } = await db.from("feedback_posts").insert({
      author_name: author, is_anonymous: anonymous, category, resource_code: resource,
      title, body,
    }).select(POST_FIELDS).single();
    if (error) throw error;
    const saved = await db.from("feedback_post_credentials")
      .insert({ post_id: post.id, password_salt: salt, password_hash: hash });
    if (saved.error) {
      await db.from("feedback_posts").delete().eq("id", post.id);
      throw saved.error;
    }
    return json(req, { post }, 201);
  }
  if (action === "report") {
    rejectBotTrap(b.website);
    await consume(req, "report", clientId, 5, 3600);
    const id = uuid(b.id);
    const reason = clean(b.reason);
    if (!REASONS.has(reason)) throw new HttpError(400, "신고 사유를 선택해 주세요.");
    await visiblePost(id);
    const reporterHash = await sha256(serviceKey + ":report:" + requester(req) + ":" + clientId);
    const { error } = await db.from("feedback_reports").insert({
      post_id: id, reason, reporter_hash: reporterHash,
    });
    if (error && error.code !== "23505") throw error;
    return json(req, { ok: true });
  }
  if (action === "delete") {
    await consume(req, "delete", clientId, 8, 3600);
    const id = uuid(b.id);
    const password = sized(b.password, 6, 64, "게시글 비밀번호");
    const { data, error } = await db.from("feedback_post_credentials")
      .select("password_salt,password_hash").eq("post_id", id).maybeSingle();
    if (error) throw error;
    if (!data) throw new HttpError(404, "게시글을 찾을 수 없습니다.");
    const computed = await passwordHash(password, data.password_salt);
    if (!constantTimeEqual(computed, data.password_hash))
      throw new HttpError(403, "게시글 비밀번호가 맞지 않습니다.");
    const result = await db.from("feedback_posts").delete().eq("id", id);
    if (result.error) throw result.error;
    return json(req, { ok: true });
  }
  if (action.startsWith("admin-")) {
    await moderator(req);
    if (action === "admin-list") {
      const offset = Math.min(10000, Math.max(0, Number(b.offset) || 0));
      return json(req, await postsList(offset, true));
    }
    const id = uuid(b.id);
    if (action === "admin-reply") {
      const body = sized(b.body, 2, 1500, "답변");
      const { error } = await db.from("feedback_replies").upsert({
        post_id: id, body, updated_at: new Date().toISOString(),
      }, { onConflict: "post_id" });
      if (error) throw error;
      return json(req, { ok: true });
    }
    if (action === "admin-status") {
      const status = clean(b.status);
      if (!STATUSES.has(status)) throw new HttpError(400, "잘못된 상태입니다.");
      const { error } = await db.from("feedback_posts")
        .update({ status, updated_at: new Date().toISOString() }).eq("id", id);
      if (error) throw error;
      return json(req, { ok: true });
    }
    if (action === "admin-hide") {
      if (typeof b.hidden !== "boolean") throw new HttpError(400, "잘못된 숨김 값입니다.");
      const { error } = await db.from("feedback_posts")
        .update({ hidden: b.hidden, updated_at: new Date().toISOString() }).eq("id", id);
      if (error) throw error;
      return json(req, { ok: true });
    }
    if (action === "admin-delete") {
      const { error } = await db.from("feedback_posts").delete().eq("id", id);
      if (error) throw error;
      return json(req, { ok: true });
    }
  }
  throw new HttpError(404, "지원하지 않는 기능입니다.");
}

Deno.serve(async (req: Request) => {
  const origin = req.headers.get("origin");
  if (origin && !ALLOWED_ORIGINS.has(origin))
    return json(req, { error: "허용되지 않은 출처입니다." }, 403);
  if (req.method === "OPTIONS") return json(req, { ok: true });
  try {
    return await handle(req);
  } catch (error) {
    const known = error instanceof HttpError;
    if (!known) console.error("[feedback-api] " + (error instanceof Error ? error.message : "unknown error"));
    return json(req, { error: known ? error.message : "요청을 처리하지 못했습니다." },
      known ? error.status : 500);
  }
});
