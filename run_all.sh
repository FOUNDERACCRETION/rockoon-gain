#!/bin/sh
# Regenerate every number, table and figure, check the paper, and compile it.
set -e
cd "$(dirname "$0")/code"
python rockoon_model.py
python similarity.py
python extra_tests.py
python cea_check.py
python rocketpy_check.py
python revision_tests.py
python realdata_tests.py
python aerobee_check.py
python run_analysis.py
python make_tables.py
python make_figures.py
python -m pytest -q test_model.py
python check_paper.py
cd ../paper
latexmk -pdf -interaction=nonstopmode paper.tex
cd ../supplement
latexmk -pdf -interaction=nonstopmode attribution.tex
latexmk -pdf -interaction=nonstopmode audit_2022.tex
