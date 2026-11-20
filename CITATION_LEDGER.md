# Citation ledger (state on 19 Nov 2026)

Status codes
- OPENED: the source was opened on 18 Nov 2026 and the cited facts were read from it.
- RECORD: bibliographic record confirmed from a publisher/index listing; full text not read.
- CARRIED: carried over from the October 2026 draft; not re-opened in this pass.

Every entry must reach "checked by an author against the original" before upload.

| Key | Status | What it is cited for | Still to do |
|---|---|---|---|
| ludwig2004 | OPENED (machine-read PDF) | 28 Aug 1952 flight: fired 38,000 ft, summit about 200,000 ft (p. 8); cost < USD 1800 vs about USD 25,000 Aerobee (p. 7) | Read pp. 7-8 yourself |
| watson1949 | OPENED (OCR text, garbled in places) | Table I gross weights, impulses, max velocities for rounds 2-9; 70 deg launch; 6.75 in drag diameter; max velocity at 2.9 s (round 6); propellant about 103 lb; rounds 8-9 at lower grain temperature | **Check every Table I number against the scanned page.** Propellant weight and grain temperatures were only partly legible |
| okninski2016 | OPENED by the research pass (full text, author-uploaded copy); abstract page opened separately | 20-40 km launch altitudes, 15 deg planar trajectory, latex balloon clusters, 0.92 Isp factor, nozzle exit limited to stage diameter, total fin span 3d, "highest launch altitude best", negligible relative wind under the balloon | **An author must read the full paper** and confirm each of these statements; the research pass noted an inconsistency between its Tables 3 and 4 |
| li2022 | RECORD (publisher page blocked) | Multi-stage rockoon optimisation for orbit | Confirm authors, volume, pages on the SAGE page |
| jorda2012 | RECORD (IAF listing) | Rockoons for small satellites | Confirm paper code |
| shoyama2021 | OPENED (MDPI open-access page, 18 Nov) | Crane-suspended rail launcher, CMG azimuth control, three model-rocket launches; contains no analysis of apogee against launch altitude | Open access: each author should read it once |
| shoyama2026 | OPENED | Balloon launch experiment; elevation oscillation open | Confirm the "open problem" wording against the paper |
| traveler | OPENED (mirror copy) | +20 % drag tuning; 339,800 ft apogee | Confirm the URL resolves; the copy read was a mirror |
| astronautix | RECORD (read by the research pass, not re-opened) | 285 s, 3415 K, 1720 m/s at O/F 2.61 | Now backed by NASA CEA run directly (cea_check.py); open the page once |
| nioshudmh, nioshno2 | OPENED | Exposure limits, vapour pressure, flash point, flammable range | none |
| nist | OPENED | Liquid enthalpies of formation: UDMH 48.3 kJ/mol, N2O4 -19.56 kJ/mol | none |
| gordon1994, mcbride2002 | RECORD (NTRS listings seen 18 Nov 2026) | CEA method; thermodynamic data | none |
| ussa1976, scipy | CARRIED (standard references) | Atmosphere; integrator | Confirm each record |
| cantera | OPENED | Software version and DOI (Zenodo record is version 3.2.0) | none |
| csbf | OPENED | 0.002 cm film thickness | none |
| asm, astm | CARRIED | Yield strengths 276 MPa and 205 MPa | Confirm values and editions |
| varghese1999 | RECORD (journal page seen) | Vol. 49, No. 1, pp. 71-78; 316.7 s at O/F 1.5 for the 30 % Al gel | none |
| jones2018 | RECORD (NTRS copy seen) | ICES-2018-81; Falcon 9 about USD 2,720 per kg | none |
| iarc | RECORD (IARC and NCBI pages seen) | Vol. 71 (1999); 1,1-dimethylhydrazine Group 2B | none |
| inspace, bva, peso, msihc | CARRIED | Regulatory paragraph | Confirm each, or cut the paragraph |
| bollermann1968 | OPENED (machine-read PDF) | Super Loki Dart masses, thrust, drag tables, computed trajectory (113.4 km), White Sands flights about 125 km | **Check the drag tables and Table 3.1-3.3 against the scan**; the Mach 2.50 entry was illegible and is not used |
| sarigul2005 | OPENED (full text, author-uploaded copy) | Loss terms; Minotaur 36,200 kg; altitudes to 30,490 m; Mach 0-3; >160 POST runs; nozzle kept at 10:1; very little improvement above 15,240 m | Page range now 569-575 (AIAA table of contents, via the research pass); the copy read ends its text on p. 572. Confirm on the AIAA page |
| shapley1953, saltelli2010 | CARRIED (standard references, not opened) | Attribution method; Sobol estimator | Confirm each record. The code uses the Jansen-type estimators described in Saltelli et al. |

okninski2016 additionally supports two statements taken from the research pass, not from my own reading: exit diameter capped at body diameter, and "highest launch altitude best". Confirm both in the full text.
| shoyama2022 | RECORD + abstract | Hybrid rockoon study; abstract attributes the gain to a larger nozzle expansion ratio | **Read the full paper** before upload; contents beyond the abstract are not described in the manuscript |
| corliss1971 | RECORD (NASA and NTRS pages) | Sounding-rocket history | Confirm it mentions the 1952 rockoons where cited |
| vanallen1954, vanallen1959 | RECORD (author's publication list) | Original descriptions; the paper says they were not consulted. 1954 chapter pp. 53-64, 1959 chapter pp. 143-164 | Obtain if possible |
| yajima2009 | RECORD (publisher page) | Zero-pressure balloon principles | none |
| niskanen2013, niskanen2009 | OPENED (technical documentation v13.05, Chapter 3 and Appendix B, machine-read) | Drag build-up equations 3.78-3.94 and the tangent-ogive rule of Appendix B, as implemented in DragBuildUp | **Check each implemented equation against the PDF**; the interpolation between Mach 1 and 1.3 is our cubic through the stated end values and slopes, and the subsonic nose term is fitted at Mach 1 |
| jansen1999 | RECORD (first page seen) | Sobol estimators used in the code | none |
| janhunen2018 | OPENED (publisher page) + research pass (full text) | 18 km chosen qualitatively for a steam-balloon launcher. Now cited as published: Aeronaut. J. 123(1263), 600-616 (2019), doi:10.1017/aer.2019.10 | none |
| nasaldb | OPENED (machine-read PDF) | NASA 29.47 and 59.84 MCF balloons: 3600 lb and 2751 lb, 0.8 and 0.4 mil film | Document carries no number; title taken from its heading. Confirm |
| nist (helium) | OPENED | Helium density 43.571 kg/m3 at 290 K, 30 MPa | none |
| smith1954 | OPENED (machine-read scan of NRL Report 4276 on DTIC; Deacon section and data sheets) | Firing altitude and zenith of nine Deacon rockoon flights (SUI-3, 4, 7, 13, 15, 23; NRL-1, 4, 5); nominal payload 50 lb | **Check all nine pairs against the scanned data sheets.** Several other sheets were garbled or had uncertain designations and were left out; the SUI-5/SUI-6 sheets conflict with the summary table. Confirm author and report number on the title page |
| allen1958 | RECORD (NTRS record and report title page seen; summary read) | That the drag group is the parameter of the Allen-Eggers ballistic-entry solution | **Open the report and confirm the form of the parameter before upload** |
| sutton2017 | CARRIED (not opened) | Over-expanded sea-level nozzles; separation at an exit pressure of roughly 0.4 of ambient | **Confirm the criterion, edition and page in the book, or replace with a source you have read** |
| bloostar2015 | RECORD (IAF directory abstract page) | Zero 2 Infinity's Bloostar balloon-assisted launcher for microsatellites; first author A. Schoenmaker | **Check the full author list and paper number on the IAF page** |
| jaxa2026 | OPENED (JAXA press release, 28 May 2026) | AstroX (rockoon launch service) and JAXA joint demonstration of a suspension-type attitude-control device for stratospheric balloons | none |
| tsien1951 | RECORD (Caltech repository abstract; Scilit metadata) | Classical theory of vertical sounding-rocket performance (optimum thrust programming) | Open the paper before upload |
| traveler | OPENED (white paper, Sec. IV.1 p. 10) | 6-DOF sim tuned by -0.9% thrust and +20% drag, attributed to damage seen on recovered vehicle | none |

| cea2026 | OPENED (documentation site nasa.github.io/cea, version 3.3.4; PyPI metadata: NASA, Apache-2.0) | NASA CEA run directly as a cross-check of the thermochemistry | none |
| rocketpy2021 | RECORD (ASCE page: J. Aerosp. Eng. 34(6), 2021; authors Ceotto, Schmitt, Alves, Pezente, Carmo; DOI from the URL) | RocketPy as an independent trajectory code (version 1.13.0 used) | Confirm the first and last author names and the article number on the ASCE page |
| astronautixdeacon | OPENED (Encyclopedia Astronautica, Deacon Rockoon entry) | Gross mass 93 kg (205 lb), payload 17 kg (37 lb) | Secondary source with no references; the entry text also says 18 kg payload. Find a primary source (vanallen1954) |
| aerostar2005 | OPENED (Aerostar PDF, Rev. 17-Feb-05) | Volume, balloon weight and maximum payload of eight small zero-pressure balloons (Tables 1 and 2) | Check the eight rows against the PDF |
| rainwater2004 | OPENED (Aerostar copy of the paper) + RECORD (ADS: Adv. Space Res. 33, 1648, 2004, doi:10.1016/j.asr.2003.07.037) | 2.1 million m3 design, 774 kg, 60 kg payload (Fig. 4) | Confirm the end page on the publisher page |
Removed 18 Nov 2026: jones1955 (never obtained; cited only to say so). Removed from the paper to the supplement: varghese1999, jones2018.

Moved to the supplement (supplement/audit_2022.tex) with the audit and regulation text, no longer cited in the paper: varghese1999, jones2018, inspace, bva, peso, msihc. iarc and the two NIOSH entries remain in the paper's safety note.

Corrections applied on 18 Nov 2026 from the second research pass: Sarigul-Klijn pages 569-575; 6061-T6 276 MPa relabelled as typical; Van Allen & Gottlieb pages added. **arXiv endorsement: physics subject classes are separate endorsement domains, so a physics.gen-ph endorsement does not cover physics.app-ph.**

Removed since the October draft because they could not be verified: Van Allen & Gottlieb 1954 (page range),
NASA NSS 1740.16 (cancelled standard; hydrogen flammability numbers dropped), McDowell 2018, the USC news
items, and the web propellant table.

## Added on 19 Nov 2026 (restructuring pass)

New status code
- MACHINE: opened on 19 Nov 2026 through a web-fetch tool that returns a machine summary of the page or of the OCR text of a
  scanned PDF. Quotations were requested verbatim, but **no human has compared them with the scan**. Every number below
  marked MACHINE must be checked by an author against the original before upload.

| Key | Status | What it is cited for | Still to do |
|---|---|---|---|
| busse1966 | MACHINE (NTRS PDF of NASA TR R-226, OCR text) | Aerobee 150A propulsion data: chamber pressure 324.0 psia, sea-level thrust 4100 lb, flow 20.71 lb/s, duration 51.5 s (Table 1: 50.9 s, used in the code), area ratio 4.6, O/F 2.56, Isp 198 s; propellant weights 758.2 lb oxidizer, 303.3 lb fuel, 5.15 lb helium; tank working pressure 500 psig; "Rocket (empty) 279.1" lb; booster 18,600 lb for 2.5 s, 104 lb/s; nominal performance at sea level with a 150 lb payload and ogive nose: 158 statute miles, 6580 ft/s at 132,000 ft | **Check every one of these numbers against the scanned pages.** Confirm the year (1966 assumed from the NTRS number) and add page numbers. Table 3 of the report (weights) was garbled in the OCR and is not used |
| busse1967 | MACHINE (NTRS PDF, OCR text; two separate extractions agreed on the numbers used) | Nominal Aerobee 150 performance with a 250 lb net payload: 5870 ft/s at 128,000 ft, peak 128 statute miles. Flights from White Sands in 1964 with stated net payloads: 4.109 GG (231.6 lb, 130.7 mi, 5900 ft/s), 4.110 GG (236.5 lb, 129 mi, 5600 ft/s at 130,000 ft), 4.118 NA (314.8 lb, 97.5 mi), 4.120 CG (356.9 lb, 89 mi, 4650 ft/s at 108,000 ft) | **Check all four flights and the nominal figures against the scan** (the text of 4.109 GG prints "13 statute miles"; its figure gives 130.70). **Find the report number** (the OCR text shows none) and add it to the reference. Launch angles are in the per-flight figures ("effective launch angle"); adding them would let the flights be flown at their real angles |
| openaerospace | MACHINE (web page) | Item-by-item masses of the Aerobee 150A sustainer and booster in kg (tank assembly 74.156, thrust chamber 14.56, fins 12.701, and so on; booster 153.314 + 118.841) | Secondary compilation with no references. Its oxidizer mass (343.914 kg) equals the 758.2 lb of busse1966 and its dry items sum to within 1% of the 279.1 lb empty weight, which is why it is used. **Replace with a primary weight table if TR R-226 Table 3 can be read from the scan** |
| kaymont | MACHINE (vendor page listing the Kaymont/Totex HAB series) | Burst diameters of nine latex balloons, 350-4000 g: 13.5, 19.8, 23.0, 25.8, 28.3, 31.0, 34.6, 42.7, 50.59 ft | Vendor page, not the manufacturer. **Check against Kaymont's own specification table and cite that instead if it agrees** |
| farside | MACHINE (StratoCat page) | Far Side balloon: 3,750,000 ft3, 2-mil polyethylene, 1500 lb; flight table lists a 2400 lb payload; planned 80,000-100,000 ft | Secondary. The page lists its sources (Navigation 6, 34-42, 1958; TIME, 7 Aug 1957). The firing altitude of each flight is not stated on the page; the paper says the vehicle was designed to fire at about 30 km. **Open the page; prefer the 1958 Navigation paper if obtainable** |
| arca2010 | MACHINE (Universe Today article, N. Atkinson, 5 Oct 2010) | Helen 2 fired from a balloon at 14,000 m; 18 km planned for the larger system | Press article. Open once |
| eoportal | MACHINE (eoPortal page) | Bloostar ignited "at an altitude of 20-25 km"; prototype test of 1 March 2017 at 25 km (a figure caption on the page says about 22 km); 4.9 t; methane and oxygen | Open once; note the 22/25 km discrepancy |
| purdue2018 | MACHINE (Purdue news release, 4 Apr 2018) | Leo Aerospace launches the rocket 11 miles (18 km) up | none beyond opening it |
| ansys2021 | MACHINE (Ansys blog, 14 Dec 2021; the URL redirects to ansys.synopsys.com) | SpaceRyde balloon lifts the rocket to 30 km | Company-supplied blog. Open once and update the URL if it has moved |
| b2space | MACHINE (company page) | Balloon lifts the rocket to 35 km; payloads up to 200 kg | Open once |
| parsonson2025 | MACHINE (European Spaceflight, 22 Nov 2025) | B2Space test: rocket launched at 21.5 km from the Canary Islands; 35 km for the operational system | Open once; replace the ?p= URL by the article's permanent link |
| jiji2024 | MACHINE (Jiji Press on Nippon.com, 22 Dec 2024) | AstroX rocket to be launched from a balloon about 20 km up | Open once. A May 2026 Japanese report gives 20-25 km for the FOX2 mission; cite it instead if you can read it |

Not cited, but used in the code and in REVISION_NOTES: the Black Brant V motor (loaded 2803 lb, propellant 2243 lb). These two
numbers come from a search-engine excerpt of the NASA Sounding Rocket Program Handbook vehicle sheet and were **not** read
from the document. They are kept out of the paper for that reason. If an author confirms them in the handbook, the result
(results.json, `realdata/solid_BlackBrantV`) can be cited next to the Super Loki value.

Changes to existing entries on 19 Nov 2026
- sutton2017 is now cited in Sec. 5 (design search) for the separation limit of 0.4 of ambient; still CARRIED, still to confirm.
- bollermann1968 is now also cited for the inert mass of the Super Loki booster, (28.32 - 16.87 - 6.12)/16.87 = 0.32 kg per kg
  of propellant; the three masses are those already in the ledger.
- allen1958, shapley1953, jansen1999, saltelli2010 are cited once each in Sec. 6 of the paper and again in the supplement.
- The supplement (supplement/attribution.tex) has its own reference list, built from the same entries.

