# HAFS Biology Lab — Public feedback board

Status: **feature branch only.** Do not merge into main or publish until the existing backend and safety checks are reviewed.

## Selected architecture

- Reuse the existing **page-view** Supabase project: predator-prey-simulation-2
  (project ref: szbmpsvxzrnewyzqiokr). Do NOT create another project.
- That project already serves assets/js/page-views.js. Do not change or
  delete its existing page-view tables, functions, settings, or keys.
- ChatGPT's current Supabase connection is to a **different account**.
  Project changes require its owning account in the dashboard or reconnecting.
- Cloudflare Turnstile is **not required and not used** in this version.
- Public reads and immediate anonymous/named writes are allowed; only
  a verified teacher Auth user can moderate.

## Check permissions

1. Open https://supabase.com/dashboard/project/szbmpsvxzrnewyzqiokr
2. Click the organization at top left (cavempt Org) and open
   Organization Settings > Team / Members.
3. Confirm logged-in member's role. Owner or Administrator has project
   management permission. Developer can work with database content but has
   limited project settings; Read-Only cannot make changes.
4. Open SQL Editor and Edge Functions to check availability. Do not run SQL
   mutations until the existing schema is understood.

Reference: https://supabase.com/docs/guides/platform/access-control

## User flow

- Home: latest 2 public posts, full board link, quick compose.
- /feedback/: most recent 15/page, anonymous or named submission,
  teacher replies, reporting and password-protected deletion.
- /feedback/admin/: teacher login, 50-item pagination, reply,
  status changes, hiding/unhiding and deletion.
- All comments are displayed as plain text (no visitor-provided HTML).
- Categories: bug, content, idea, question, other. Resources: portal, A1–C99.

## Authorization

- Browser: static GitHub Pages ES modules.
- Server: Supabase Edge Function feedback-api (verify_jwt=false).
- Six new tables named feedback_*: posts, post_credentials, replies,
  reports, moderators, rate_limits.
- RLS is enabled on every new table and anonymous/authenticated direct
  table operations are denied; server operations use its own secret key.
- Teacher operations require a verified Auth JWT and membership in
  feedback_moderators; never accept a name or browser flag as proof.
- Author post deletion uses salted PBKDF2-SHA256 password hashing.

## No-CAPTCHA tradeoffs

The Cloudflare widget and both Turnstile keys are removed. Protection uses:
- Atomic database quotas via feedback_consume_quota: 3 posts/10m per
  browser+network, 5 reports/hour per browser+network, 120 requests/hour and
  320 requests/day per network.
- Honeypot input for simple generic bots.
- Input validation and safe text-only rendering.
- Reporting, moderator hiding/deletion.
- Emergency pause: Supabase Edge Function server setting
  FEEDBACK_PUBLIC_POSTING_ENABLED = false prevents new posts but allows
  reading and reports.

**Limit:** browser identifiers can be replaced, and network/IP signals can be
unreliable. Deliberate attackers may bypass these measures and cause spam
or consume the project's shared resources. Be prepared to disable posting
and add stronger bot verification if abuse occurs. CORS and honeypots are
NOT proof of a genuine human visitor.

## Deployment on the page-view project

1. First run only this **read-only** SQL in the selected project's SQL Editor:

       SELECT schemaname, tablename FROM pg_tables
       WHERE schemaname = 'public' ORDER BY tablename;

   Check existing table/function names before adding any objects.
2. Review and apply supabase/feedback/schema.sql, which creates only
   feedback_* tables and one feedback_* quota RPC.
3. Run database security advisors. Check that anon and authenticated roles
   cannot directly access feedback_* tables or the private quota RPC.
4. Authentication > Users: create a dedicated teacher account. Confirm the
   current project's Auth configuration before disabling any signup setting.
   Insert the teacher's verified Auth UID into feedback_moderators.
5. Deploy supabase/functions/feedback-api/index.ts, with
   supabase/config.toml specifying verify_jwt=false. Its admin actions
   authenticate the token inside the handler.
6. assets/js/feedback-config.js already contains the existing project URL
   and function URL, but its publishableKey is intentionally empty. Add only
   the existing project's public sb_publishable_... key AFTER the function
   is deployed and its permissions checked. No secret belongs in GitHub.
7. Run all tests below before merging the feature branch to main.

## Retention

feedback_rate_limits stores hashes, not plaintext IP addresses. The schema
adds a created_at timestamp to support cleanup. Decide a retention policy;
for example, an operator-approved SQL cleanup after 7 days:

    DELETE FROM public.feedback_rate_limits
    WHERE created_at < now() - interval '7 days';

Also agree retention and removal policies for published posts, reports and
teacher replies, particularly any student personal information.

## Must-pass tests

1. Existing page-view counter and all A1–C2 links still work.
2. Visitors can read public posts; hidden posts are not exposed.
3. Anonymous and named posts appear immediately after submission.
4. Honeypot, invalid payloads and over-quota writes are rejected.
5. Teacher can authenticate, reply, change status, hide/unhide and delete.
6. No-session and non-moderator JWTs cannot call admin actions.
7. Data API denies anon reads/writes to new tables and private quota RPC.
8. Incorrect password cannot delete; correct one can delete its own post.
9. User input never executes HTML or scripts.
10. Report reasons/counts and deduplication work as specified.
11. Emergency pause blocks new posts while preserving reads and reports.
12. Test iPad Safari, school network, keyboard focus and modal cancellation.
13. Monitor quota and ensure sharing the DB does not disrupt page views.

## Rollback

All code is on feature/public-feedback-board, not main. If problems occur
after public launch, revert the frontend or unset publishableKey. Keep
student data and page-view schema intact; never reset the Supabase project
to remove this feature.

## Official references

- https://supabase.com/docs/guides/platform/access-control
- https://supabase.com/docs/guides/functions/secrets
- https://supabase.com/docs/guides/functions/auth
- https://supabase.com/docs/guides/api/securing-your-api
- https://supabase.com/docs/guides/functions/examples/rate-limiting
