// Public-only configuration. Configure AFTER inspecting the existing page-view project; keep feedback_* tables separate.
// Never place sb_secret, service_role, password or server-only secrets in this file.
export const FEEDBACK_CONFIG = Object.freeze({
  // Existing HAFS page-view project. Endpoint is not yet deployed.
  functionUrl: "https://szbmpsvxzrnewyzqiokr.supabase.co/functions/v1/feedback-api",
  supabaseUrl: "https://szbmpsvxzrnewyzqiokr.supabase.co",
  publishableKey: "sb_publishable_bOM16W0j9Dh6U72M_fNp8Q_HwRZm9tF",
  // Fail-closed staging. Turn on only after the deployed server's
  // FEEDBACK_PUBLIC_POSTING_ENABLED=true and the moderation checks pass.
  postingEnabled: false,
});
