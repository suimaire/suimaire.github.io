-- HAFS Biology Lab / public feedback board
-- Apply only to a designated feedback Supabase project. Do not apply to student-record projects.
-- The browser receives no service/secret key. All database access goes through feedback-api.
-- Apply the entire migration as one transaction: failure rolls back all new objects.

BEGIN;

create table public.feedback_posts (
  id uuid primary key default gen_random_uuid(),
  author_name text not null check (char_length(author_name) between 1 and 24),
  is_anonymous boolean not null default true,
  category text not null check (category in ('bug','content','idea','question','other')),
  resource_code text not null default 'portal'
    check (resource_code = 'portal' or resource_code ~ '^[ABC][1-9][0-9]?$'),
  title text not null check (char_length(title) between 2 and 80),
  body text not null check (char_length(body) between 5 and 1500),
  status text not null default 'received'
    check (status in ('received','reviewing','done')),
  hidden boolean not null default false,
  created_at timestamptz not null default now(),
  updated_at timestamptz not null default now()
);

create index feedback_posts_visible_recent_idx
  on public.feedback_posts (created_at desc, id desc) where hidden = false;

create table public.feedback_post_credentials (
  post_id uuid primary key references public.feedback_posts(id) on delete cascade,
  password_salt text not null,
  password_hash text not null
);

create table public.feedback_replies (
  post_id uuid primary key references public.feedback_posts(id) on delete cascade,
  body text not null check (char_length(body) between 2 and 1500),
  created_at timestamptz not null default now(),
  updated_at timestamptz not null default now()
);

create table public.feedback_reports (
  id uuid primary key default gen_random_uuid(),
  post_id uuid not null references public.feedback_posts(id) on delete cascade,
  reason text not null check (reason in ('spam','abuse','privacy','other')),
  reporter_hash text not null,
  created_at timestamptz not null default now(),
  unique(post_id, reporter_hash)
);

create index feedback_reports_post_id_idx
  on public.feedback_reports(post_id);

create table public.feedback_moderators (
  user_id uuid primary key references auth.users(id) on delete cascade,
  created_at timestamptz not null default now()
);

create table public.feedback_rate_limits (
  key_hash text not null,
  action text not null,
  window_bucket bigint not null,
  hit_count integer not null default 1 check (hit_count >= 1),
  created_at timestamptz not null default now(),
  primary key (key_hash, action, window_bucket)
);

create index feedback_rate_limits_created_at_idx
  on public.feedback_rate_limits(created_at);

-- Every table in the exposed public schema is locked down, including credential hashes.
alter table public.feedback_posts enable row level security;
alter table public.feedback_post_credentials enable row level security;
alter table public.feedback_replies enable row level security;
alter table public.feedback_reports enable row level security;
alter table public.feedback_moderators enable row level security;
alter table public.feedback_rate_limits enable row level security;

revoke all on public.feedback_posts from public, anon, authenticated;
revoke all on public.feedback_post_credentials from public, anon, authenticated;
revoke all on public.feedback_replies from public, anon, authenticated;
revoke all on public.feedback_reports from public, anon, authenticated;
revoke all on public.feedback_moderators from public, anon, authenticated;
revoke all on public.feedback_rate_limits from public, anon, authenticated;

grant select, insert, update, delete on public.feedback_posts to service_role;
grant select, insert, update, delete on public.feedback_post_credentials to service_role;
grant select, insert, update, delete on public.feedback_replies to service_role;
grant select, insert, update, delete on public.feedback_reports to service_role;
grant select, insert, update, delete on public.feedback_moderators to service_role;
grant select, insert, update, delete on public.feedback_rate_limits to service_role;

-- One atomic counter. Only the server-side service role can call it.
-- UTC epoch buckets keep the quota independent of application server clocks.
create function public.feedback_consume_quota(
  p_key_hash text, p_action text, p_window_seconds integer, p_max_hits integer
) returns boolean
language plpgsql
security invoker
set search_path = ''
as $$
declare
  v_hits integer;
  v_bucket bigint;
begin
  if p_window_seconds < 30 or p_max_hits < 1
    or char_length(p_key_hash) != 64 or char_length(p_action) > 32 then
    return false;
  end if;
  v_bucket := floor(extract(epoch from now()) / p_window_seconds)::bigint;

  insert into public.feedback_rate_limits(key_hash, action, window_bucket, hit_count)
  values (p_key_hash, p_action, v_bucket, 1)
  on conflict (key_hash, action, window_bucket)
  do update set hit_count = public.feedback_rate_limits.hit_count + 1
  where public.feedback_rate_limits.hit_count < p_max_hits
  returning hit_count into v_hits;
  return v_hits is not null;
end
$$;

revoke all on function public.feedback_consume_quota(text,text,integer,integer)
  from public, anon, authenticated;
grant execute on function public.feedback_consume_quota(text,text,integer,integer)
  to service_role;

-- Moderator bootstrap (after creating a teacher Auth user):
-- insert into public.feedback_moderators(user_id) values ('<verified-auth-user-uuid>');
-- Password hashes, rate limits, reports, moderator IDs must never be returned to visitors.

COMMIT;
