import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, existsSync } from 'node:fs';
import { createHash } from 'node:crypto';
const read = p => readFileSync(new URL('../' + p, import.meta.url), 'utf8');
test('runtime data equals supplied teaching fixture exactly', () => {
  const fixture = JSON.parse(read('assets/data/pcr-primer-fixture.json'));
  assert.equal(createHash('sha256').update(JSON.stringify(fixture)).digest('hex'), '3f404112df85687ccb27a813069129cfb60b6907ea824a102033f6d594490589');
  if (existsSync(new URL('../_codex/pcr_primer_teaching_fixture.json', import.meta.url))) assert.deepEqual(fixture, JSON.parse(read('_codex/pcr_primer_teaching_fixture.json')));
});
test('worksheet source, scripts and styling contain no forbidden decorative dots', () => {
  for (const path of ['bioinformatics/pcr-primer-design.html', '_layouts/pcr-worksheet.html', 'assets/css/pcr-worksheet.css', 'assets/js/pcr-worksheet.mjs', 'assets/js/pcr-core.mjs', 'assets/js/pcr-records.mjs', 'assets/js/pcr-intro.mjs', 'assets/js/pcr-design.mjs', 'assets/js/pcr-workbench.mjs', 'assets/js/pcr-review.mjs', 'assets/js/pcr-review-view.mjs', 'assets/js/pcr-evidence.mjs', 'assets/js/pcr-evidence-view.mjs', 'assets/js/pcr-external.mjs', 'assets/js/pcr-external-view.mjs', 'assets/js/pcr-final.mjs', 'assets/js/pcr-final-view.mjs']) assert.doesNotMatch(read(path), /[\u00b7\u2022\u2027\u2219\u22c5\u30fb\u318d]/u, path);
});
test('no whole-storage clear or unsafe imported HTML rendering', () => {
  const source = read('assets/js/pcr-worksheet.mjs') + read('assets/js/pcr-intro.mjs') + read('assets/js/pcr-workbench.mjs') + read('assets/js/pcr-review-view.mjs') + read('assets/js/pcr-evidence-view.mjs') + read('assets/js/pcr-external-view.mjs') + read('assets/js/pcr-final-view.mjs');
  assert.doesNotMatch(source, /localStorage\.clear|innerHTML|insertAdjacentHTML|eval\(/);
  assert.doesNotMatch(read('_layouts/pcr-worksheet.html'), /head_custom|page-views|heading-numbers|supabase/i);
});
test('dedicated permalink and one automatically numbered portal heading', () => {
  assert.match(read('bioinformatics/pcr-primer-design.html'), /permalink: \/bioinformatics\/pcr-primer-design\//);
  const home = read('index.md');
  assert.equal((home.match(/<h3><a[^>]+>PCR과 프라이머 디자인<\/a><\/h3>/g) || []).length, 1);
  assert.ok(home.indexOf('>생물정보학 기초</a></h3>') < home.indexOf('>PCR과 프라이머 디자인</a></h3>'));
  assert.doesNotMatch(home.match(/<li class="portal-resource">\s*<h3><a[^>]+>PCR과 프라이머 디자인[\s\S]*?<\/li>/)[0], /1\.3\.2/);
});
