# Rockoon float altitude and the envelope rule

Code and manuscript for "The mass-optimal float altitude for balloon-launched sounding rockets is robust; the mass saving is not" (Ramnani, Aduru, Reddy, 2026).

## Check validity

    python validate.py           # quick checks in about 20 seconds: physics identities, published reference
                                 # values (NACA Deacon, NASA Super Loki, propellant benchmark), NASA CEA
                                 # and RocketPy cross-checks, headline results against results.json,
                                 # unit tests, paper consistency
    python validate.py --full    # also regenerates every result in a temporary copy and compares
                                 # the generated LaTeX files byte for byte (about two hours on 2 cores)

## Reproduce everything

One command (about an hour on 4 cores):

    pip install -r requirements.txt
    ./run_all.sh

On 18 Nov 2026 `python validate.py --full` rebuilt every result of that version in a fresh copy and reproduced
`paper/results.tex` and every table byte for byte. After the restructuring of 19 Nov 2026 the quick checks pass
and the earlier results are unchanged (see REVISION_NOTES.md), but `--full` has not been repeated: run it before upload.
`code/test_model.py` holds eleven checks
against values that can be worked out by hand (sea-level atmosphere, column mass = p/g, buoyancy
balance, rocket equation in vacuum, energy conservation on the coast, attribution sums to the gain,
the envelope rule against direct minimisation).

Step by step:

Requires Python 3.10+ and a LaTeX installation.

    pip install numpy scipy matplotlib cantera==3.2.0 cea==3.3.4 rocketpy==1.13.0
    cd code
    python rockoon_model.py     # builds thermo_table.json (about 1 minute)
    python similarity.py        # ensemble, optimum altitude, Sobol -> results/similarity.json (about 5 minutes)
    python extra_tests.py       # H3 on four vehicles, H1 at a fixed float altitude -> results/extra.json
    python cea_check.py         # thermochemistry against NASA CEA run directly -> results/cea_check.json
    python rocketpy_check.py    # baseline trajectories against RocketPy -> results/rocketpy_check.json
    python revision_tests.py    # envelope fit, design searches, penalty grid, two-stage and solid comparators
                                # -> results/revision.json (about 20 minutes on 2 cores)
    python realdata_tests.py    # published dry masses and fineness-dependent drag
                                # -> results/realdata.json (about 45 minutes on 2 cores)
    python aerobee_check.py     # Aerobee 150: thermochemistry (Cantera) and trajectories -> results/aerobee_check.json
    python run_analysis.py      # all numbers -> results/results.json, paper/results.tex (about 6 minutes;
                                # calls revision_analysis.py, general_rule.py and realdata_analysis.py)
    python make_tables.py       # table bodies -> paper/tab_*.tex
    python make_figures.py      # figures -> paper/figures/*.pdf
    python check_paper.py       # consistency checks; must print "none" on every line
    cd ../paper
    latexmk -pdf paper.tex

Every computed number in `paper.tex` is a `\R{key}` macro read from `results.tex`, so the text
cannot drift from the code. The random seeds are fixed (20261118, 20261119, 20261120 and 20261121).

## Layout

- `code/rockoon_model.py`  atmosphere, balloon, thermochemistry, nozzle, trajectory, tanks
- `code/run_analysis.py`   every calculation in the paper
- `code/revision_tests.py`, `code/revision_analysis.py`  comparators, penalties, envelope fit, analytic optimum
- `code/general_rule.py`   the envelope rule for any envelope law and price (Sec. 4.3)
- `code/realdata_tests.py`, `code/realdata_analysis.py`  published dry masses, fineness-dependent drag
- `code/aerobee_check.py`  check against the Aerobee 150
- `code/cea_check.py`, `code/rocketpy_check.py`  cross-checks against NASA CEA and RocketPy
- `code/make_tables.py`, `code/make_figures.py`, `code/check_paper.py`
- `paper/`                 LaTeX source, generated tables and figures
- `results/results.json`   all results, including Monte Carlo samples
- `CITATION_LEDGER.md`     verification state of each reference


## Supplements

`supplement/attribution.tex` holds the study of where the apogee gain comes from (attribution, similarity test, robustness, Monte Carlo and Sobol analysis), which the paper summarises in Sec. 6. It reads `paper/results.tex`, the tables and the figures of the paper, and takes its cross-references to the paper from `paper/paper.aux`, so compile the paper first.

`supplement/audit_2022.tex` holds the full audit of the 2022 Project AKS proposal (sixteen errors, their effects and the checks that would have caught them) and the hazards and regulation notes. The paper no longer refers to it. `run_all.sh` compiles both after the paper.

## arXiv submission

`./make_arxiv.sh` builds `arxiv_submission.tar.gz`: the manuscript and its generated inputs at the top level,
and a copy of this repository's code, results, paper source and supplement under `anc/`, which arXiv
publishes as ancillary files.

## arXiv category

Primary physics.app-ph, cross-list physics.class-ph. Not physics.space-ph (space plasma physics).
