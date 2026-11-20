#!/bin/sh
# Build the arXiv source tarball: the manuscript with its generated inputs at the top level, and the code,
# results and supplement as ancillary files under anc/ (arXiv lists files in anc/ on the abstract page).
set -e
ROOT="$(cd "$(dirname "$0")" && pwd)"
OUT="$ROOT/arxiv_build"
rm -rf "$OUT"; mkdir -p "$OUT/figures" "$OUT/anc/code" "$OUT/anc/results" "$OUT/anc/supplement" "$OUT/anc/paper/figures"
cp "$ROOT"/paper/paper.tex "$ROOT"/paper/results.tex "$ROOT"/paper/tab_*.tex "$OUT/"
cp "$ROOT"/paper/figures/*.pdf "$OUT/figures/"
cp "$ROOT"/code/*.py "$ROOT"/code/thermo_table.json "$OUT/anc/code/"
cp "$ROOT"/results/*.json "$OUT/anc/results/"
cp "$ROOT"/supplement/attribution.pdf "$ROOT"/supplement/attribution.tex "$ROOT"/supplement/audit_2022.pdf "$ROOT"/supplement/audit_2022.tex "$OUT/anc/supplement/"
cp "$ROOT"/paper/paper.tex "$ROOT"/paper/results.tex "$ROOT"/paper/tab_*.tex "$OUT/anc/paper/"
cp "$ROOT"/paper/figures/*.pdf "$OUT/anc/paper/figures/"
cp "$ROOT"/run_all.sh "$ROOT"/requirements.txt "$ROOT"/validate.py "$ROOT"/README.md "$ROOT"/CITATION_LEDGER.md "$OUT/anc/"
# compile once to make sure the tarball builds on its own
cd "$OUT" && pdflatex -interaction=nonstopmode paper.tex >/dev/null && pdflatex -interaction=nonstopmode paper.tex >/dev/null
rm -f paper.aux paper.log paper.out
mv paper.pdf "$ROOT/arxiv_check.pdf"
cd "$OUT" && tar czf "$ROOT/arxiv_submission.tar.gz" .
echo "wrote $ROOT/arxiv_submission.tar.gz"
