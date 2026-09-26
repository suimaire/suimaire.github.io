// Worksheet-only local preview when Ruby is unavailable. This is not a Jekyll build.
import { createServer } from 'node:http';
import { readFile, stat } from 'node:fs/promises';
import { resolve, extname, sep } from 'node:path';
import { fileURLToPath } from 'node:url';
const root = fileURLToPath(new URL('..', import.meta.url));
const pagePath = '/bioinformatics/pcr-primer-design/';
const relative = source => source.replace(/\{\{\s*'([^']+)'\s*\|\s*relative_url\s*\}\}/g, '$1');
export async function previewHtml() {
  const page = (await readFile(resolve(root, 'bioinformatics/pcr-primer-design.html'), 'utf8')).replace(/^---[\s\S]*?---\s*/, '');
  const layout = await readFile(resolve(root, '_layouts/pcr-worksheet.html'), 'utf8');
  return relative(layout.replace('{{ content }}', page).replaceAll('{{ page.title | escape }}', 'PCR과 프라이머 디자인').replaceAll('{{ site.title | escape }}', 'HAFS Biology Lab').replaceAll('{{ page.description | escape }}', 'PCR 학습지 로컬 미리보기').replaceAll('{{ page.url | absolute_url }}', `http://127.0.0.1:4173${pagePath}`));
}
export function startPreview(port = 4173) {
  const server = createServer(async (req, res) => {
    try {
      const url = new URL(req.url, `http://${req.headers.host}`);
      if ([pagePath, pagePath + 'index.html', '/'].includes(url.pathname)) { res.setHeader('content-type', 'text/html; charset=utf-8'); res.end(await previewHtml()); return; }
      if (!url.pathname.startsWith('/assets/')) { res.writeHead(404); res.end('Worksheet-only preview'); return; }
      const path = resolve(root, '.' + decodeURIComponent(url.pathname));
      if (!path.startsWith(resolve(root, 'assets') + sep) || !(await stat(path)).isFile()) throw new Error('Not found');
      const types = { '.mjs': 'text/javascript', '.js': 'text/javascript', '.css': 'text/css', '.json': 'application/json', '.ttf': 'font/ttf' };
      res.setHeader('content-type', types[extname(path)] || 'application/octet-stream'); res.end(await readFile(path));
    } catch { res.writeHead(404); res.end('Not found'); }
  });
  return new Promise(resolveServer => server.listen(port, '127.0.0.1', () => resolveServer(server)));
}
if (process.argv[1] && resolve(process.argv[1]) === fileURLToPath(import.meta.url)) {
  await startPreview(); console.log(`Worksheet preview: http://127.0.0.1:4173${pagePath}`);
}
