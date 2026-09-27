# Phase 1–6 v1 compatibility fixtures

These are synthetic teaching/test records, never real student records or actual NCBI results.

| Fixture | Historical implementation | Origin |
| --- | --- | --- |
| phase1.json | 41a02c9 | Historical `emptyRecord`, populated with introductory and design answers |
| phase2.json | 09fd27e | Existing Phase 2 Chromium 1440 JSON export |
| phase3.json | deb572c | Existing Phase 3 Chromium 1440 JSON export |
| phase4.json | 5fa41c7 | Existing Phase 4 Chromium 1440 JSON export |
| phase5.json | 25073ab | Existing Phase 5 Chromium 1440 JSON export |
| phase6.json | fcfa16f | Historical `emptyRecord` and `selectFinalDesign`, populated with final reflections |

Every file was serialized through that commit's actual `parseRecord` implementation. This removes later optional defaults that earlier regression runs may have added to exported files. Phase 1 and 6 fixtures are reconstructed examples, not preserved original browser exports. Phase 7 loads each fixture from localStorage and imports it through the file picker, then checks every answer, design, optional state and final snapshot. No answer keys or schema versions are changed by the polish.
