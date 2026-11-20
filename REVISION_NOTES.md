# Revision of 19 November 2026: one paper about one result

This pass carries out the eight steps of the journal review. The uploaded zip of 18 November is the "before" state.
`NEXT_STEPS.md` lists what only a person can do. Read its Sec. 1 first: several new source numbers were read by
machine and must be checked by an author before upload.

## What changed in the paper

- **Structure.** The paper is now about the float altitude and the rule that fixes it. Main text 22 pages (was 31),
  29 pages in all (was 35). The attribution, similarity and robustness study (old Secs. 2.7, 4, 5.5, 6 and Appendices
  A-C in part) moved, unchanged apart from cross-references, to `supplement/attribution.tex` (14 pages) and is
  summarised in a new Sec. 6. The hypotheses H1-H3 moved with it.
- **Abstract** rewritten around three findings. Title unchanged.
- **Sec. 4.3, the envelope rule, is new and replaces the old Sec. 5.3.** The optimum condition is derived for any
  envelope law m_e = E(V, m_s) and any price on the envelope (Eq. 7). The old Eq. (6) is its film-law case. New: the
  one-step estimate (Eq. 9) with L_eff = (1/L + alpha_2/H_rho)^-1, Table 5 (four envelope laws) and Fig. 5.
  It agrees with direct minimisation to 6 significant figures, and one step from a 20 km trial is within 0.3 km.
- **A correction.** The old text said the optimum rises with vehicle size and that the fitted envelope law puts it
  "2-3 km higher". With the fitted law the optimum *falls* with size (22.8 to 18.6 km over the aluminium families;
  your own `scenario_fittedNone` shows 22.3 to 19.0 km), and against the film law it moves by -2.0 to +6.5 km
  depending on the family. The rule explains it: the sign of the trend is the sign of 1 - alpha - beta.
- **Latex balloons and price.** Burst-diameter law fitted to nine catalogue balloons (exponent 0.64); clusters
  (alpha = 1) have their optimum at 22-26 km. A price on the balloon system lowers the optimum by at most 4-7 km;
  a cheap envelope raises it to 30-35 km, which is where Okninski et al. found it.
- **Sec. 4.4 is new:** a table of firing altitudes stated by ten programmes and studies.
- **Sec. 3.4 is new: Aerobee 150**, a flown pressure-fed storable-propellant sounding rocket. Nominal peak
  altitudes reproduced to +3.7% and -2.0%; four 1964 flights over-predicted by 7-17%; delivered specific impulse
  0.91 of the ideal value by the paper's method, against the 0.92 assumed. The sentence "The liquid-engine model is
  not validated against flight" is gone from the abstract; the limitation that remains is stated in Sec. 3.5.
- **Sec. 5 (saving), with published dry masses.**
  - Tank mass of the Aerobee 150A sustainer, 0.16 kg per kg of propellant (your aluminium estimate is 1.3 times
    that): saving -6 to +8%, negative for 3 of 7 sizes.
  - Inert mass of the Super Loki booster, 0.32 kg per kg (generic motor: 0.15): saving 6-23%, against -7 to +11%.
  - Skin friction that grows with fineness: every ground vehicle still chooses fineness 25 but is 3-22% heavier;
    saving 23-41% (steel) and 0 to 22% (aluminium), 3-8 points more.
  - The cost paragraph now states the bound instead of leaving r open.
- **Notation.** Envelope area is S (was A, which is the body cross-section); gas ratio is chi (was gamma); gondola
  fraction is f_g (was g); geopotential altitude is h (was H); nozzle half-angle is theta (was alpha). The abstract
  uses H_rho.
- **Validation tables** (Deacon rounds, rockoon flights, closure) moved to Appendix A, with the new Aerobee table.
- **The 2022 audit** is no longer mentioned in the paper. `supplement/audit_2022.tex` is kept, with its title and
  section references updated.

## New and changed files

- `code/general_rule.py` (new): the rule for any envelope law; Table 5; numbers `gr*`, `ltx*`, `clu*`, `fs*`.
- `code/realdata_tests.py` (new, about 45 minutes on 2 cores) and `code/realdata_analysis.py` (new): published
  dry masses, fineness-dependent drag; numbers `rd*` and `ab*`; `paper/tab_aerobee.tex`.
- `code/aerobee_check.py` (new, needs Cantera): thermochemistry and trajectories of the Aerobee 150.
- `code/rockoon_model.py`: `FinenessDrag`, and the options `fin_drag` (in `vehicle_run`) and `cd_fn` (in
  `two_stage_apogee`). Defaults unchanged.
- `code/revision_tests.py`: `_two_stage` accepts an optional tenth argument for fineness-dependent drag.
- `code/run_analysis.py`: calls the two new analysis modules. `code/make_figures.py`: `fig_general`.
- `code/check_paper.py`: also checks the supplement's result keys and the table files.
- `code/test_model.py`: 11 tests (three new). `validate.py`: 52 quick checks (Aerobee and envelope-rule groups new).
- `run_all.sh`, `make_arxiv.sh`: include the new scripts and `supplement/attribution.pdf`.
- `paper/paper.tex` rewritten as above; `supplement/attribution.tex` new; `NEXT_STEPS.md` new.

## What was verified on 19 November

- All 854 result values of the 18 November version and its 14 generated tables are unchanged, byte for byte, after
  `run_analysis.py` and `make_tables.py` were run with the modified model code.
- `python validate.py` (quick): all checks pass. `validate.py --full` was **not** run in this pass; it now takes
  about two hours on 2 cores. Run it before upload.
- With the modified `rockoon_model.py` and `revision_tests.py`, `similarity.py`, `extra_tests.py` and
  `revision_tests.py` were run again in a separate copy: `similarity.json`, `extra.json` and `revision.json` are
  byte-for-byte identical to the uploaded ones. The new code changes no earlier result.

## Not done

- The Van Allen originals (interlibrary loan) and an expert reading: see `NEXT_STEPS.md`.
- Prices: no published source fixed the price ratio; the paper gives a bound instead.
- `nhsjs/` was not updated.

---

# Revision of 18 November 2026

## What changed in the paper

- **New title**: "The mass-optimal float altitude for balloon-launched sounding rockets is robust; the mass saving is not".
- **Sec. 5.3 (new): analytic optimum.** Exact first-order condition: at the optimum the envelope weighs
  (3/2) H_rho |dm_s/dz| (Eq. 6). The closed form (Eq. 7) has length scale L_eff = 4.1-5.1 km. It reproduces
  the full model's optimum to within 0.7 km and predicts the sigma and gondola sensitivities and the mass slope.
  The ground comparator does not enter the condition.
- **Sec. 5.4 (new): whether a rockoon saves mass.**
  - Ground and balloon vehicles are each designed over a wider space (fineness 8-25, T/W 1.5-6, exit pressure down
    to 0.4 of ambient). Every ground optimum sits at the fineness bound; going to 30-35 would cost the rockoon a
    further 2-9%.
  - Central penalty case (10 kg/m2 nozzle wall, 20% gondola, k = 0.1) next to the worst case.
  - Envelope law fitted to 11 published balloons: m_e = a A + b m_s V^(1/3), a = 7.1 g/m2, b = 6.77e-3 /m,
    13% r.m.s. The paper's old law (sigma A + k m) fits them at 44% r.m.s.
  - Two-stage pressure-fed and single- and two-stage solid ground comparators.
  - Cost test rerun with both vehicles designed for their own launch.
  - Monte Carlo of the saving: 4000 samples per family, giving a probability of a saving.
  - Table 11 collects the comparators; Fig. 9 shows the Monte Carlo.
- **Validation.**
  - NASA CEA run directly (Python package `cea` 3.3.4) agrees with the Cantera thermochemistry to 0.14%.
  - RocketPy 1.13 (3- and 6-DOF) agrees with the integrator to 0.03% in apogee.
  - Unfitted rockoon test with the 205 lb gross mass from Encyclopedia Astronautica: the model over-predicts every
    flight. Matching the Iowa flights would need 2.1x the drag curve.
- **Restructuring.**
  - The six-group regression is moved to Appendix B.
  - Fig. 11(a) and the "316.7 s (proposal)" label are removed.
  - The bibliography is renumbered in citation order.
  - Stale sentences in Secs. 6.1, 7.1 and 7.2 are fixed, and "nothing for aluminium" is reworded.
- **Headline numbers that moved.**
  - Saving against the best ground option with the central penalties: 20-34% (steel), -3 to +14% (aluminium),
    -7 to +11% (solid).
  - Probability of a saving: at least 99.9% (steel), 39-99% (aluminium).
  - Optimum: 19-25 km with the fitted envelope law; 15-27 km over the Monte Carlo.

## New and changed files

- `code/revision_tests.py`: design searches, penalty grid, two-stage and solid comparators, envelope fit.
- `code/revision_analysis.py`: numbers, tables and Monte Carlo (called from run_analysis.py).
- `code/cea_check.py` and `code/rocketpy_check.py`: the cross-checks.
- Other changed files:
  - `code/rockoon_model.py`: `balloon_law2`, `two_stage_apogee`.
  - `code/run_analysis.py`, `make_figures.py` and `check_paper.py`.
  - `validate.py`: now 40 checks.
  - `run_all.sh`, `requirements.txt`, `README.md` and `CITATION_LEDGER.md`.
- `make_arxiv.sh` builds `arxiv_submission.tar.gz`, with the code, results and supplement under `anc/`.

## Still for the authors to do

1. **Endorsement** in physics.app-ph from an existing arXiv author.
2. **Make the GitHub repository public** before submitting, or remove the URL. It was private on 18 Nov.
3. **Van Allen & Gottlieb (1954)**, pp. 53-64 of *Rocket Exploration of the Upper Atmosphere* (Boyd and Seaton,
   eds., Pergamon), via interlibrary loan. If it gives the rocket and payload masses, rerun the rockoon comparison
   with them.
4. **Work through `CITATION_LEDGER.md`**: every "Still to do" entry, including the five new references.
5. **The `nhsjs/` version was not updated** in this revision.
