// Public-only configuration. Configure AFTER provisioning an isolated feedback project.
// Never place sb_secret, service_role, password or Turnstile secret in this file.
export const FEEDBACK_CONFIG = Object.freeze({
  functionUrl: "",
  supabaseUrl: "",
  publishableKey: "",
  turnstileSiteKey: "",
});
