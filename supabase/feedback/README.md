# HAFS Biology Lab — public feedback board

Status: **feature branch only**. Do not merge/deploy until the separate backend and safety checklist below are completed.

## User flow

- Home page: latest 2 public posts, "전체 의견 보기", "의견 남기기".
- /feedback/: public timeline (15/page), anonymous or named posts, view teacher replies,
  report posts, and delete one's own post using a per-post password.
- /feedback/admin/: teacher email/password login, reply, update status, hide/unhide,
  permanently delete. No public registration is needed.
- New posts appear immediately after successful database insertion. No approval step.
- Plain-text rendering via textContent; no HTML or Markdown from visitors.
- Categories: bug / content / idea / question / other. Resource: portal or A1–C99.

## Architecture and data protection

- Frontend: static GitHub Pages + ES modules under assets/js.
- Server: one Supabase Edge Function, feedback-api.
- Database: schema.sql creates posts, credentials, replies, reports, moderators, rate limits.
- All six tables have RLS enabled and no public CRUD grants or policies.
- The Edge Function is the only access path, using the server-side service/secret key.
- Teacher authorization checks a verified Auth user against feedback_moderators.
- Post passwords: random per-post salt + 180,000-iteration PBKDF2/SHA-256.
- Abuse controls: Cloudflare Turnstile server-side Siteverify and atomic browser/network quotas.
- Turnstile secret is mandatory for anonymous submissions. No secret = writes fail closed.
- Browser IDs/IP-based quotas are not a guarantee against deliberate attackers; monitor reports.
- Neither IP addresses nor raw passwords are stored in post tables. Hashed abuse identifiers
  are stored and should have a documented retention/cleanup policy.

## Deployment prerequisites

**Do not use a database that contains student screen-time uploads or English participation
records.** Provision/designate a separate Supabase project for the feedback board first.
The current connected Supabase account does not grant access to the legacy page-view project.

1. Create a separate Supabase project in the desired organization and region after reviewing
   any billing implications. Do not assume a third free project is available.
2. Review and run supabase/feedback/schema.sql on the chosen project. Prefer generating a
   regular timestamped migration through the installed Supabase CLI for repository-managed history.
   Run Supabase security advisors; confirm all six feedback tables deny anon/authenticated access.
3. In Supabase Authentication, create the teacher's account manually. Disable public signup
   unless it is needed for other parts of the dedicated project. Add its verified Auth user UUID
   to public.feedback_moderators with the example INSERT at the end of schema.sql.
4. Create a Cloudflare Turnstile widget for the exact hostname suimaire.github.io.
   Configure the private FEEDBACK_TURNSTILE_SECRET in the Edge Function's **server-only** secrets.
   Never commit the secret key.
5. Deploy supabase/functions/feedback-api/index.ts with the setting in supabase/config.toml.
   The function must use verify_jwt=false because anonymous visitors are allowed;
   admin endpoints verify JWTs and authorization inside the handler.
6. Set the four public values in assets/js/feedback-config.js:
   - functionUrl = https://<project-ref>.supabase.co/functions/v1/feedback-api
   - supabaseUrl = https://<project-ref>.supabase.co
   - publishableKey = the project's **sb_publishable_...** public key
   - turnstileSiteKey = the Turnstile public site key
   NEVER put a sb_secret/service_role key in a browser file.
7. Verify the deployment tests below BEFORE merging this branch into main.

Note: recent Supabase Data API configurations may require explicit exposure/grants for
newly created public-schema tables and functions. The only callable RPC needed by
feedback-api is feedback_consume_quota, granted to service_role only.
Confirm PostgREST is configured to permit the service-role requests needed by the function.

## Must-pass tests before public launch

1. Guest reads latest posts and detail without login; hidden posts do not appear.
2. A visitor can submit both as "익명" and as a named author. Post appears instantly.
3. The stored displayed name is "익명" for anonymous submissions.
4. Missing/invalid/expired/reused Turnstile token is rejected server-side.
5. Input length, invalid resource codes and excessive posting are rejected server-side.
6. HTML injection payload displays as literal plain text, not executable markup.
7. Wrong deletion password does not delete; correct password deletes only its own post.
8. Reports appear in moderator counts, and duplicate reports from the same identifier
   do not count twice.
9. Non-moderator Auth users and missing/expired JWTs cannot call admin actions.
10. A moderator can answer, set statuses, hide/reveal and delete posts.
11. Anonymous Data API requests to ALL feedback tables and RPC are denied.
12. Check mobile Safari/iPad layout, keyboard focus, modal cancellation, page scrolling,
    and behavior on the school's actual network.
13. Confirm existing home course links A1–C2 and page-view counter still work.
14. Observe moderation workload and define data retention, handling of student information,
    content takedown response time and abuse-response procedures.

## Rollback

Before production, the entire feature is isolated to feature/public-feedback-board.
If deployment fails later, revert the feature's frontend merge or clear functionUrl
to disable frontend requests. Do not drop feedback tables if legitimate posts exist:
archive/export responsibly and plan deletion with the operator.

## References

- Supabase Edge secrets: https://supabase.com/docs/guides/functions/secrets
- Edge function JWT configuration: https://supabase.com/docs/guides/functions/function-configuration
- Data API grants and RLS: https://supabase.com/docs/guides/api/securing-your-api
- Turnstile server validation: https://developers.cloudflare.com/turnstile/get-started/server-side-validation/
