// Public-only configuration. Configure AFTER inspecting the existing page-view project; keep feedback_* tables separate.
// Never place sb_secret, service_role, password or server-only secrets in this file.
export const FEEDBACK_CONFIG = Object.freeze({
  functionUrl: "",
  supabaseUrl: "",
  publishableKey: "",
});
