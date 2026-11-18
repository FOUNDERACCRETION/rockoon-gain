# Rockoon apogee-gain decomposition

Code and manuscript for "The mass-optimal float altitude for balloon-launched sounding rockets, and where the rockoon advantage comes from" (Ramnani, Aduru, Reddy, 2026).

## Check validity

    python validate.py           # 33 checks in about 20 seconds: physics identities, published reference
                                 # values (NACA Deacon, NASA Super Loki, propellant benchmark), headline
                                 # results against results.json, unit tests, paper consistency
    python validate.py --full    # also regenerates every result in a temporary copy and compares
                                 # the generated LaTeX files byte for byte (about 9 minutes; 46 checks)

## Reproduce everything

One command (about 10 minutes on 4 cores):

    pip install -r requirements.txt
    ./run_all.sh

On 18 Nov 2026 a fresh clone with every generated file deleted was rebuilt this way and gave a
`paper/results.tex` byte-identical to the committed one. `code/test_model.py` holds eight checks
against values that can be worked out by hand (sea-level atmosphere, column mass = p/g, buoyancy
balance, rocket equation in vacuum, energy conservation on the coast, attribution sums to the gain).

Step by step:

Requires Python 3.10+ and a LaTeX installation.

    pip install numpy scipy matplotlib cantera==3.2.0
    cd code
    python rockoon_model.py     # builds thermo_table.json (about 1 minute)
    python similarity.py        # ensemble, optimum altitude, Sobol -> results/similarity.json (about 5 minutes)
    python run_analysis.py      # all numbers -> results/results.json, paper/results.tex (about 2 minutes)
    python make_tables.py       # table bodies -> paper/tab_*.tex
    python make_figures.py      # figures -> paper/figures/*.pdf
    python check_paper.py       # consistency checks; must print "none" on every line
    cd ../paper
    latexmk -pdf paper.tex

Every computed number in `paper.tex` is a `\R{key}` macro read from `results.tex`, so the text
cannot drift from the code. The random seeds are fixed (20261118 and 20261119).

## Layout

- `code/rockoon_model.py`  atmosphere, balloon, thermochemistry, nozzle, trajectory, tanks
- `code/run_analysis.py`   every calculation in the paper
- `code/make_tables.py`, `code/make_figures.py`, `code/check_paper.py`
- `paper/`                 LaTeX source, generated tables and figures
- `results/results.json`   all results, including Monte Carlo samples
- `CITATION_LEDGER.md`     verification state of each reference


## Supplement

`supplement/audit_2022.tex` holds the full audit of the 2022 Project AKS proposal (sixteen errors, their effects and the checks that would have caught them) and the hazards and regulation notes. It reads the same `paper/results.tex`; `run_all.sh` compiles it after the paper.

## arXiv category

Primary physics.app-ph, cross-list physics.class-ph. Not physics.space-ph (space plasma physics).
