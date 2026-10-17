import { feedbackReady, feedbackApi, teacherLogin } from "./feedback-client.js";

const root = document.querySelector("[data-feedback-admin]");
if (root) {
  const login = root.querySelector("[data-admin-login]");
  const panel = root.querySelector("[data-admin-panel]");
  const feed = root.querySelector("[data-admin-list]");
  const notice = root.querySelector("[data-admin-notice]");
  const authNotice = root.querySelector("[data-admin-auth-notice]");
  let token = "";

  const say = (target, message) => { target.hidden = !message; target.textContent = message; };
  const dateText = value => new Date(value).toLocaleString("ko-KR", { dateStyle: "short", timeStyle: "short" });
  function el(tag, cls, text) {
    const e = document.createElement(tag);
    e.className = cls || "";
    if (text !== undefined) e.textContent = text;
    return e;
  }
  function logout(message = "") {
    token = "";
    login.hidden = false;
    panel.hidden = true;
    feed.replaceChildren();
    say(authNotice, message);
  }
  async function mutate(action, payload) {
    try {
      await feedbackApi(action, { body: payload, token });
      await load();
    } catch (error) {
      say(notice, error.message);
    }
  }
  function adminCard(post) {
    const card = el("article", "feedback-admin-card");
    const h = el("h3", "", post.title);
    const meta = el("p", "feedback-post__meta",
      [post.author_name, post.resource_code, post.category, dateText(post.created_at),
        post.report_count + "건 신고", post.hidden ? "숨김" : "공개"].join(" · "));
    const body = el("p", "feedback-detail__body", post.body);
    const row = el("div", "feedback-actions");
    const status = el("select", "feedback-admin-select");
    [["received","접수"],["reviewing","검토 중"],["done","반영 완료"]].forEach(([key,label]) => {
      const option = el("option", "", label);
      option.value = key;
      status.append(option);
    });
    status.value = post.status;
    status.setAttribute("aria-label", "진행 상태");
    status.addEventListener("change", () => mutate("admin-status", { id: post.id, status: status.value }));

    const hidden = el("button", "feedback-button feedback-button--secondary", post.hidden ? "공개 복원" : "글 숨기기");
    hidden.type = "button";
    hidden.addEventListener("click", () => mutate("admin-hide", { id: post.id, hidden: !post.hidden }));

    const deleteButton = el("button", "feedback-button feedback-button--danger", "영구 삭제");
    deleteButton.type = "button";
    deleteButton.addEventListener("click", () => {
      if (confirm("게시글과 답변을 영구 삭제할까요? 복구할 수 없습니다."))
        mutate("admin-delete", { id: post.id });
    });
    row.append(status, hidden, deleteButton);

    const replyLabel = el("label", "feedback-field", "교사 답변");
    const reply = el("textarea", "");
    reply.rows = 3;
    reply.maxLength = 1500;
    reply.value = post.reply?.body || "";
    replyLabel.append(reply);
    const save = el("button", "feedback-button feedback-button--primary", "답변 저장");
    save.type = "button";
    save.addEventListener("click", () => mutate("admin-reply", { id: post.id, body: reply.value.trim() }));
    card.append(h, meta, body, row, replyLabel, save);
    return card;
  }
  async function load() {
    if (!token) return;
    say(notice, "의견을 불러오는 중입니다.");
    try {
      const { posts, has_more } = await feedbackApi("admin-list", { body: { offset: 0 }, token });
      feed.replaceChildren();
      if (!posts.length) feed.append(el("p", "feedback-empty", "접수된 의견이 없습니다."));
      else feed.append(...posts.map(adminCard));
      say(notice, has_more ? "최근 50건만 표시됩니다. 이전 의견은 추후 페이지 이동 기능으로 확장할 수 있습니다." : "");
      login.hidden = true;
      panel.hidden = false;
    } catch (error) {
      if (/로그인|권한|만료/.test(error.message)) logout(error.message);
      else say(notice, error.message);
    }
  }
  login.addEventListener("submit", async event => {
    event.preventDefault();
    if (!feedbackReady) return say(authNotice, "게시판 연결을 준비하고 있습니다.");
    const button = login.querySelector('[type="submit"]');
    button.disabled = true;
    say(authNotice, "로그인하고 있습니다.");
    try {
      const email = login.querySelector('[name="email"]').value.trim();
      const password = login.querySelector('[name="password"]').value;
      const session = await teacherLogin(email, password);
      token = session.access_token;
      login.reset();
      say(authNotice, "");
      await load();
    } catch (error) {
      say(authNotice, error.message);
    } finally {
      button.disabled = false;
    }
  });
  root.querySelector("[data-admin-logout]")?.addEventListener("click", () => logout());
  root.querySelector("[data-admin-refresh]")?.addEventListener("click", load);
  if (!feedbackReady) say(authNotice, "독립 Supabase 프로젝트 연결 및 보안 인증 설정 전입니다.");
}
