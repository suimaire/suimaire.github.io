import {
  feedbackReady, feedbackApi, feedbackClientId, initFeedbackCaptcha,
} from "./feedback-client.js";

const CATEGORY = { bug: "작동 오류", content: "내용 오류", idea: "개선 제안", question: "질문", other: "자유 의견" };
const STATUS = { received: "접수", reviewing: "검토 중", done: "반영 완료" };
const home = document.querySelector("[data-feedback-home]");
const board = document.querySelector("[data-feedback-board]");
const dateText = value => new Date(value).toLocaleDateString("ko-KR", {
  year: "numeric", month: "short", day: "numeric",
});

function element(tag, className = "", value = null) {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (value !== null) node.textContent = value;
  return node;
}
function showNotice(target, text, failure = false) {
  if (!target) return;
  target.textContent = text;
  target.hidden = !text;
  target.classList.toggle("feedback-notice--error", failure);
}
function showDialog(dialog) {
  const previousFocus = document.activeElement;
  document.body.append(dialog);
  dialog.addEventListener("close", () => {
    dialog.remove();
    if (previousFocus instanceof HTMLElement) previousFocus.focus();
  }, { once: true });
  if (typeof dialog.showModal === "function") dialog.showModal();
  else dialog.setAttribute("open", "");
  dialog.querySelectorAll("[data-close]").forEach(button => button.addEventListener("click", () => dialog.close()));
}

function postCard(post) {
  const item = element("li", "feedback-post");
  const button = element("button", "feedback-post__open");
  button.type = "button";
  const heading = element("span", "feedback-post__title", post.title);
  const meta = element("span", "feedback-post__meta");
  meta.textContent = [CATEGORY[post.category] || "의견", post.resource_code === "portal" ? "전체 포털" : post.resource_code,
    post.author_name, dateText(post.created_at)].join(" · ");
  const lower = element("span", "feedback-post__status", STATUS[post.status] || "접수");
  if (post.reply) lower.textContent += " · 교사 답변";
  button.append(heading, meta, lower);
  button.addEventListener("click", () => openDetail(post.id));
  item.append(button);
  return item;
}

function initCompose() {
  const dialog = element("dialog", "feedback-dialog");
  dialog.innerHTML = [
    '<form class="feedback-form" novalidate>',
    '<div class="feedback-dialog__top"><h2>의견 남기기</h2><button type="button" data-close class="feedback-close" aria-label="닫기">×</button></div>',
    '<p class="feedback-help">게시글은 제출 즉시 전체 공개됩니다. 개인정보와 타인의 실명은 게시하지 마세요.</p>',
    '<fieldset class="feedback-fieldset"><legend>작성자 표시</legend>',
    '<label><input type="radio" name="identity" value="anonymous" checked> 익명</label>',
    '<label><input type="radio" name="identity" value="named"> 이름·닉네임</label></fieldset>',
    '<label class="feedback-field feedback-nickname" hidden>표시할 이름 또는 닉네임<input name="author" maxlength="24" autocomplete="nickname" placeholder="24자 이내"></label>',
    '<div class="feedback-form__row"><label class="feedback-field">분류<select name="category">',
    '<option value="idea">개선 제안</option><option value="bug">작동 오류</option>',
    '<option value="content">내용 오류</option><option value="question">질문</option><option value="other">자유 의견</option>',
    '</select></label><label class="feedback-field">관련 자료 번호 (선택)<input name="resource" maxlength="3" placeholder="예: A3 또는 B4"></label></div>',
    '<label class="feedback-field">제목<input name="title" required minlength="2" maxlength="80" placeholder="어떤 의견인가요?"></label>',
    '<label class="feedback-field">내용<textarea name="body" required minlength="5" maxlength="1500" rows="6" placeholder="의견을 구체적으로 적어 주세요."></textarea></label>',
    '<label class="feedback-field">게시글 비밀번호<input name="password" type="password" required minlength="6" maxlength="64" autocomplete="new-password" placeholder="나중에 글을 삭제할 때 필요합니다."></label>',
    '<p class="feedback-help">비밀번호는 서버에서 해시하여 저장합니다. 분실하면 본인이 직접 삭제할 수 없습니다.</p>',
    '<div class="feedback-captcha" data-captcha></div>',
    '<p class="feedback-notice" role="status" aria-live="polite" data-notice hidden></p>',
    '<div class="feedback-actions"><button type="button" data-close class="feedback-button feedback-button--secondary">취소</button>',
    '<button type="submit" class="feedback-button feedback-button--primary">작성 즉시 공개</button></div>',
    '</form>',
  ].join("");
  const form = dialog.querySelector("form");
  const notice = dialog.querySelector("[data-notice]");
  const nickname = dialog.querySelector(".feedback-nickname");
  const nameInput = form.elements.namedItem("author");
  form.addEventListener("change", event => {
    if (event.target.name === "identity") {
      const named = form.elements.namedItem("identity").value === "named";
      nickname.hidden = !named;
      nameInput.required = named;
    }
  });

  let captchaToken = "";
  let resetCaptcha = () => {};
  showDialog(dialog);
  if (!feedbackReady) {
    showNotice(notice, "게시판 서버 연결을 준비하고 있습니다.", true);
    form.querySelector('[type="submit"]').disabled = true;
    return;
  }
  initFeedbackCaptcha(dialog.querySelector("[data-captcha]"), token => { captchaToken = token; })
    .then(reset => { resetCaptcha = reset; })
    .catch(error => showNotice(notice, error.message, true));

  form.addEventListener("submit", async event => {
    event.preventDefault();
    const named = form.elements.namedItem("identity").value === "named";
    const code = String(form.elements.namedItem("resource").value || "").trim().toUpperCase();
    const title = String(form.elements.namedItem("title").value || "").trim();
    const body = String(form.elements.namedItem("body").value || "").trim();
    const password = String(form.elements.namedItem("password").value || "");
    const author = String(nameInput.value || "").trim();
    if (named && (!author || author.length > 24)) return showNotice(notice, "작성자 이름을 확인해 주세요.", true);
    if (code && !/^[ABC][1-9][0-9]?$/.test(code))
      return showNotice(notice, "자료 번호는 A3, B4, C2와 같이 입력해 주세요.", true);
    if (title.length < 2 || title.length > 80 || body.length < 5 || body.length > 1500)
      return showNotice(notice, "제목은 2~80자, 내용은 5~1500자로 입력해 주세요.", true);
    if (password.length < 6 || password.length > 64)
      return showNotice(notice, "게시글 비밀번호는 6~64자로 입력해 주세요.", true);
    if (!captchaToken) return showNotice(notice, "보안 인증을 완료해 주세요.", true);

    const submit = form.querySelector('[type="submit"]');
    submit.disabled = true;
    showNotice(notice, "게시글을 등록하고 있습니다.");
    try {
      const result = await feedbackApi("create", { body: {
        category: form.elements.namedItem("category").value,
        is_anonymous: !named, author_name: author, resource_code: code || "portal",
        title, body, password, captcha_token: captchaToken,
        client_id: feedbackClientId(),
      } });
      dialog.close();
      if (home && !board) {
        location.href = "/feedback/?post=" + encodeURIComponent(result.post.id);
      } else {
        await refreshBoard();
        openDetail(result.post.id);
      }
    } catch (error) {
      showNotice(notice, error.message, true);
      resetCaptcha();
    } finally {
      submit.disabled = false;
    }
  });
}

async function openDetail(id) {
  const dialog = element("dialog", "feedback-dialog");
  dialog.innerHTML = [
    '<article class="feedback-detail"><div class="feedback-dialog__top"><h2 data-title>의견 읽기</h2>',
    '<button type="button" data-close class="feedback-close" aria-label="닫기">×</button></div>',
    '<p class="feedback-post__meta" data-meta></p>',
    '<p class="feedback-detail__body" data-body></p>',
    '<section class="feedback-reply" data-reply hidden><h3>선생님 답변</h3><p data-reply-text></p></section>',
    '<p class="feedback-notice" data-notice role="status"></p>',
    '<div class="feedback-actions"><button type="button" data-report class="feedback-button feedback-button--secondary">신고</button>',
    '<button type="button" data-delete class="feedback-button feedback-button--secondary">내 글 삭제</button>',
    '<button type="button" data-close class="feedback-button feedback-button--primary">닫기</button></div>',
    '</article>',
  ].join("");
  showDialog(dialog);
  const notice = dialog.querySelector("[data-notice]");
  showNotice(notice, "게시글을 불러오고 있습니다.");
  try {
    const { post } = await feedbackApi("detail", { params: { id } });
    dialog.querySelector("[data-title]").textContent = post.title;
    dialog.querySelector("[data-meta]").textContent = [CATEGORY[post.category], post.resource_code === "portal" ? "전체 포털" : post.resource_code,
      post.author_name, dateText(post.created_at), STATUS[post.status]].join(" · ");
    dialog.querySelector("[data-body]").textContent = post.body;
    const replyBox = dialog.querySelector("[data-reply]");
    if (post.reply) {
      replyBox.hidden = false;
      replyBox.querySelector("[data-reply-text]").textContent = post.reply.body;
    }
    showNotice(notice, "");
    dialog.querySelector("[data-report]").addEventListener("click", () => openReport(id));
    dialog.querySelector("[data-delete]").addEventListener("click", () => openDelete(id, dialog));
  } catch (error) {
    showNotice(notice, error.message, true);
    dialog.querySelector("[data-report]").disabled = true;
    dialog.querySelector("[data-delete]").disabled = true;
  }
}
function openReport(id) {
  const dialog = element("dialog", "feedback-dialog feedback-dialog--small");
  dialog.innerHTML = [
    '<form class="feedback-form"><div class="feedback-dialog__top"><h2>게시글 신고</h2>',
    '<button type="button" data-close class="feedback-close" aria-label="닫기">×</button></div>',
    '<label class="feedback-field">신고 사유<select name="reason">',
    '<option value="spam">광고·도배</option><option value="abuse">비방·부적절한 표현</option>',
    '<option value="privacy">개인정보 노출</option><option value="other">기타</option></select></label>',
    '<p class="feedback-help">신고가 접수되면 관리자가 확인합니다. 자동으로 삭제되지는 않습니다.</p>',
    '<div class="feedback-captcha" data-captcha></div>',
    '<p data-notice class="feedback-notice" role="status" hidden></p>',
    '<div class="feedback-actions"><button type="button" data-close class="feedback-button feedback-button--secondary">취소</button>',
    '<button type="submit" class="feedback-button feedback-button--primary">신고하기</button></div></form>',
  ].join("");
  showDialog(dialog);
  let captchaToken = "";
  let resetCaptcha = () => {};
  const notice = dialog.querySelector("[data-notice]");
  initFeedbackCaptcha(dialog.querySelector("[data-captcha]"), token => { captchaToken = token; })
    .then(reset => { resetCaptcha = reset; })
    .catch(error => showNotice(notice, error.message, true));
  dialog.querySelector("form").addEventListener("submit", async event => {
    event.preventDefault();
    if (!captchaToken) return showNotice(notice, "보안 인증을 완료해 주세요.", true);
    const button = dialog.querySelector('[type="submit"]');
    button.disabled = true;
    try {
      await feedbackApi("report", { body: {
        id, reason: dialog.querySelector("select").value,
        captcha_token: captchaToken, client_id: feedbackClientId(),
      } });
      showNotice(notice, "신고가 접수되었습니다.");
      button.textContent = "신고 완료";
    } catch (error) {
      showNotice(notice, error.message, true);
      resetCaptcha();
      button.disabled = false;
    }
  });
}
function openDelete(id, parentDialog) {
  const dialog = element("dialog", "feedback-dialog feedback-dialog--small");
  dialog.innerHTML = [
    '<form class="feedback-form"><div class="feedback-dialog__top"><h2>내 글 삭제</h2>',
    '<button type="button" data-close class="feedback-close" aria-label="닫기">×</button></div>',
    '<p class="feedback-help">게시글과 선생님 답변이 함께 삭제되며 복구할 수 없습니다.</p>',
    '<label class="feedback-field">작성할 때 설정한 비밀번호<input type="password" name="password" required minlength="6" maxlength="64"></label>',
    '<p data-notice class="feedback-notice" role="status" hidden></p>',
    '<div class="feedback-actions"><button type="button" data-close class="feedback-button feedback-button--secondary">취소</button>',
    '<button type="submit" class="feedback-button feedback-button--danger">영구 삭제</button></div></form>',
  ].join("");
  showDialog(dialog);
  const notice = dialog.querySelector("[data-notice]");
  dialog.querySelector("form").addEventListener("submit", async event => {
    event.preventDefault();
    const button = dialog.querySelector('[type="submit"]');
    button.disabled = true;
    try {
      await feedbackApi("delete", { body: {
        id, password: dialog.querySelector('[name="password"]').value,
        client_id: feedbackClientId(),
      } });
      dialog.close();
      parentDialog.close();
      await refreshBoard();
      await refreshHome();
    } catch (error) {
      showNotice(notice, error.message, true);
      button.disabled = false;
    }
  });
}
async function refreshHome() {
  if (!home) return;
  const list = home.querySelector("[data-feedback-recent]");
  list.replaceChildren();
  if (!feedbackReady) return list.append(element("li", "feedback-empty", "게시판 연결을 준비 중입니다."));
  try {
    const { posts } = await feedbackApi("list", { params: { offset: 0 } });
    const recent = posts.slice(0, 2);
    if (!recent.length) list.append(element("li", "feedback-empty", "첫 의견을 남겨 주세요."));
    else list.append(...recent.map(postCard));
  } catch {
    list.append(element("li", "feedback-empty", "의견을 불러오지 못했습니다."));
  }
}
let boardOffset = 0;
async function refreshBoard(append = false) {
  if (!board) return;
  const list = board.querySelector("[data-feedback-list]");
  const more = board.querySelector("[data-feedback-more]");
  const notice = board.querySelector("[data-feedback-notice]");
  if (!append) { boardOffset = 0; list.replaceChildren(); }
  if (!feedbackReady) {
    showNotice(notice, "게시판 서버 연결을 준비하고 있습니다.", true);
    more.hidden = true;
    return;
  }
  more.disabled = true;
  showNotice(notice, "");
  try {
    const { posts, has_more } = await feedbackApi("list", { params: { offset: boardOffset } });
    if (!append && !posts.length) list.append(element("li", "feedback-empty", "아직 등록된 의견이 없습니다. 첫 의견을 남겨 주세요."));
    list.append(...posts.map(postCard));
    boardOffset += posts.length;
    more.hidden = !has_more;
  } catch (error) {
    showNotice(notice, error.message, true);
  } finally {
    more.disabled = false;
  }
}
if (home || board) {
  document.querySelectorAll("[data-feedback-compose]").forEach(button =>
    button.addEventListener("click", initCompose));
  if (home) refreshHome();
  if (board) {
    refreshBoard();
    board.querySelector("[data-feedback-more]")?.addEventListener("click", () => refreshBoard(true));
    const params = new URLSearchParams(location.search);
    if (params.has("post") && feedbackReady) openDetail(params.get("post"));
  }
}
