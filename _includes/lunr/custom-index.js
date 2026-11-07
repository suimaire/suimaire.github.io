{%- comment -%}
  The default lunr.trimmer discards Korean characters at the edges of tokens.
  Preserve Hangul terms in the indexing pipeline so Korean title/keyword
  searches work. The search pipeline already retains Hangul tokens.
  This hook runs once per indexed document, so configure only on the first.
{%- endcomment -%}
if (i === '0') {
  this.pipeline.remove(lunr.trimmer);
}
