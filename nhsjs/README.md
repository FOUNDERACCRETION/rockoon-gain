# NHSJS-format skeleton

This folder is a formatting frame for a later NHSJS submission. It holds the required section order,
12-point single-spaced layout, the figures, tables, equation and code-generated numbers, and the
reference list in the NHSJS style. It holds no body text: the grey NOTES boxes list facts only.

## Why there is no body text

The NHSJS AI usage policy (https://nhsjs.com/ai-usage-policy/, read 18 Nov 2026) allows AI for code,
analysis and literature search, and does not allow AI to draft, paraphrase or rewrite prose, including
passages the authors edit afterwards. The arXiv version's prose was drafted with AI help, so rewording
it sentence by sentence would not meet that policy. Write each section from the notes, the figures and
your own understanding, with the arXiv PDF closed. Figure and table captions are prose too: rewrite them.

## Build

    python3 build.py

writes `nhsjs_standard.pdf` (superscript numbers; the review copy) and `nhsjs_online.pdf` (full citation
in double parentheses at every location). The reference list is built from the citations in the text, in
order of first appearance. Reference texts live in `refs.tex`.

## What NHSJS asks for (LaTeX route)

1. PDF with standard citations and no names, affiliations or acknowledgements, named after the title.
2. A Word file with the online citation format. After the text is final:
   `pandoc nhsjs.tex -o paper.docx` gives a starting point; tables must be rebuilt in Word and figures
   attached as images.
3. The LaTeX source as supplementary information.

Limits: 20 pages including figures and tables; abstract 200-250 words under the four headings.
The skeleton's floats and references take about 9 pages, which leaves about 4,700 words of text.

## Before submitting

- Ask NHSJS whether a paper already on arXiv is eligible. Their pages do not say.
- `scipy` in `refs.tex` needs its full author list (NHSJS: include all author names).
- Only cite what you have read. Jones (1955) is still unread.
