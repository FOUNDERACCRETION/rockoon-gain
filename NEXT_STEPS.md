# What is left for the authors (after the restructuring of 19 Nov 2026)

This file lists what a person has to do before submission. Nothing here can be done by the code.

## 1. Verify the machine-read numbers (do this first)

The Aerobee 150 validation, the latex-balloon law and the table of programme altitudes are new. Their source
numbers were read by a tool from OCR text or web pages, not by a person. arXiv suspends authors for unchecked
machine output, so each number below must be compared with the original by one of you.

| Where it is used | Source | What to check |
|---|---|---|
| `code/aerobee_check.py`, lines `PC ... T_BOOST` | NASA TR R-226 (https://ntrs.nasa.gov/citations/19660005621) | 324.0 psia; 4100 lb; 20.71 lb/s; 50.9 s (Table 1) and 51.5 s (text); area ratio 4.6; O/F 2.56; oxidizer 758.2 lb; fuel 303.3 lb; helium 5.15 lb; "Rocket (empty) 279.1"; booster 18,600 lb, 2.5 s, 104 lb/s; nominal: 150 lb payload, sea level, ogive nose, 158 statute miles, 6580 ft/s at 132,000 ft |
| `code/aerobee_check.py`, `NOMINAL` and `FLIGHTS` | NASA compendium for 1964 (https://ntrs.nasa.gov/citations/19670025553) | 250 lb net: 5870 ft/s at 128,000 ft, 128 statute miles. Flights 4.109 GG, 4.110 GG, 4.118 NA, 4.120 CG: net payload, peak altitude, burnout velocity and altitude. Also read each flight's "effective launch angle" from its figure |
| `code/aerobee_check.py`, `SUSTAINER_ITEMS`, `BOOSTER` | https://open-aerospace.github.io/Aerobee-150 | The thirteen sustainer rows and two booster rows, in kg |
| `code/general_rule.py`, `LATEX` | https://the-rocketman.com/products/high-altitude-weather-balloons-hab | Nine burst diameters in feet. Better: Kaymont's own table |
| `code/general_rule.py`, `FARSIDE`; Table of programmes | StratoCat page and the seven pages cited in the table | One number each; see CITATION_LEDGER.md |

If a number differs, change it in the script, run `./run_all.sh`, and the paper follows.

## 2. Library requests (step 3 of the review)

None of these has a DOI, so ask by title.

1. J. A. Van Allen and M. B. Gottlieb, "The inexpensive attainment of high altitudes with balloon-launched rockets",
   in *Rocket Exploration of the Upper Atmosphere*, ed. R. L. F. Boyd and M. J. Seaton (Pergamon, London, 1954), pp. 53-64.
   **Look for:** gross mass and payload mass of the Deacon rockoons, flight by flight if given; the balloon size and
   mass; the intended firing altitude and why it was chosen; launch angle.
2. J. A. Van Allen, "Balloon-launched rockets for high-altitude research", in *Sounding Rockets*, ed. H. E. Newell, Jr.
   (McGraw-Hill, New York, 1959), Chap. 9, pp. 143-164. **Look for:** the same, plus the Loki and Hawk rockoons.
3. "The Far Side Program", *Navigation* **6**, 34-42 (1958). **Look for:** balloon volume and mass, suspended mass,
   firing altitude of each flight. It would replace the StratoCat page.

With the Deacon masses in hand, set them in `run_analysis.py` (the rockoon block, `rk_ap` and `m_ax`) in place of the
fitted masses. The rockoon check then needs no fitted parameter.

Two papers with DOIs that would strengthen Sec. 4.3 if you can read them:
- Okninski et al., doi:10.1016/j.actaastro.2016.06.010. **Check:** how they size the latex cluster (by burst volume or
  by rated payload), the price they give the balloons, and the exact sentence on the highest altitude being best.
  Sec. 4.3 says their result is the cheap-envelope limit of the same rule; make sure their assumptions agree.
- Shoyama et al., doi:10.2514/1.A35040. **Check:** the launch altitudes they study and any statement of a best one.

## 3. An expert reader (step 8)

Send the PDF to one balloon or sounding-rocket academic before submitting. The same person can be asked for the arXiv
endorsement. A draft:

> Subject: Request for comments on a short paper on rockoon float altitude
>
> Dear Professor ...,
>
> We are three students at CHIREC International School, Hyderabad. We have written a paper on the float altitude that
> minimises the launched mass of a balloon-launched sounding rocket. Its main result is an exact condition for any
> balloon: at the optimum the envelope is the fraction 1/(alpha Lambda/H - beta) of the suspended mass. It gives an
> optimum near 20 km for film balloons and explains why the answer moves so little.
>
> We would be grateful for your view on two points where your experience matters more than our model:
> (1) whether the envelope-mass law we fitted to eleven published zero-pressure balloons is reasonable for a load as
> heavy for its size as a rocket, and (2) whether the balloon data you know would change the exponents we use.
>
> The paper is 20 pages; the code and data are public. We used an AI model heavily in the analysis and writing, as the
> paper declares, and we have checked the derivation and the sources ourselves.
>
> If you find the work sound, we would also be grateful for an arXiv endorsement in physics.app-ph.
>
> Yours sincerely, ...

Before sending, each of you should be able to derive Eq. (general) on paper and explain the three checks of
`code/test_model.py` that test it.

## 4. Prices (step 6): what was found, and why it is not in the paper

No published source gave the price of a zero-pressure envelope or of a sounding rocket per kilogram, so the price
ratio r could not be fixed from the literature. What was found, from search results only (not opened, not cited):
- latex balloons, Kaymont list prices: USD 70 (600 g), 90 (800 g), 110 (1000 g), 130 (1200 g), 135 (1500 g), about
  USD 90-120 per kg;
- helium: USD 0.5-1 per cubic foot (2014), which is about USD 16-32 per kilogram lifted;
- NASA balloon programme: about USD 38 million a year for 10-15 flights (Science, 2017); commercial flights "a few
  hundred thousand dollars";
- UP Aerospace SpaceLoft XL: one mission contracted at USD 1 million (2021).

The paper therefore bounds the effect of price instead of fixing it (Sec. 4.3): no price lowers the optimum by more
than a few kilometres. A quotation from a balloon manufacturer, which the expert reader may be able to give, would
turn the bound into a number.

## 5. Decisions that are yours

- Target journal. Acta Astronautica is the natural first choice: Okninski et al., the study this paper answers most
  directly, appeared there. Advances in Space Research (COSPAR) reaches the scientific-ballooning community, which
  would judge the envelope laws. Format the paper to the chosen journal only after the checks of Sec. 1.

- The attribution study now lives in `supplement/attribution.pdf` and is summarised in Sec. 6. It could become a
  second paper.
- `supplement/audit_2022.pdf` is no longer mentioned in the paper. Keep it as an ancillary file or leave it out.
- The GitHub repository was private when this was written. Make it public before submitting, or remove the URL.
- `nhsjs/` was not updated.
