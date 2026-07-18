# Consolidated references — the weaning–survival paradox manuscript

Every citation gathered across the four models' calibration passes plus the
paradox-framing literature. Organised by the role it plays in the paper, with
the evidence grade used throughout the programme (A = directly reported primary
measurement; B = derived/interpolated here; C = secondary/textbook; X = not
found / unverified) and a one-line note on how the manuscript uses it.

Verification status: paradox-framing set (§1) surfaced and abstract-verified via
Consensus, Aug 2013–2024 span; mechanism/calibration set (§2–5) retrieved and
where possible page-verified during the model calibration passes (see
`docs/CALIBRATION.md` for the arithmetic). PMIDs/DOIs given where confirmed.

---

## 1. The paradox itself (Introduction)

| # | Citation | What it establishes | Grade |
|---|---|---|---|
| P1 | **Siempos I et al. Effect of early versus late or no tracheostomy on mortality and pneumonia. Lancet Respir Med 2014;2(2):150–8.** (13 RCTs, 2434 pts) | ICU mortality OR 0.80 (95% CI 0.59–1.09, NS); **1-year mortality RR 0.93 (0.85–1.02, NS)**; VAP lower (OR 0.60). The anchor for "no survival benefit." | A (meta) |
| P2 | **Chorath K et al. Association of Early vs Late Tracheostomy With Pneumonia and Ventilator Days. JAMA Otolaryngol Head Neck Surg 2021;147(5):450–9.** (17 RCTs, 3145 pts) | Early trach: **+1.74 ventilator-free days, −6.25 ICU days, VAP OR 0.59; mortality OR 0.66 (0.38–1.15, NS)**. The anchor for "shortens the process, not the endpoint." | A (meta) |
| P3 | **Han R et al. Effect of tracheotomy timing on patients receiving mechanical ventilation: a meta-analysis. PLOS ONE 2024;19.** (21 RCTs) | MV −2.77 d, ICU −6.36 d, **all-cause mortality RR 0.86 (0.73–1.00, NS)**. Confirms P1/P2 in the most recent large synthesis. | A (meta) |
| P4 | **Merola R et al. Timing of Tracheostomy in ICU Patients: a systematic review and meta-analysis. Life 2024;14.** (19 RCTs, 3586 pts) | The DISSENT: a modest mortality reduction (RR −0.151, p=0.040) alongside shorter ICU stay and MV. Cited for honesty — the mortality signal is not uniformly null, though TSA says inconclusive. | A (meta) |
| P5 | **Jubran A et al. Long-Term Outcome after Prolonged Mechanical Ventilation. Am J Respir Crit Care Med 2019;199(12):1508–16.** (n=315, LTACH) | 1-yr survival 66.9%; **Pimax 41.3 cmH₂O (53% predicted), maintained; handgrip 21.5% predicted, severely impaired.** Third independent MIP anchor AND direct support that respiratory strength ≠ limb strength. | A |
| P6 | **Chelluri L et al. Long-term mortality and quality of life after prolonged mechanical ventilation. Crit Care Med 2004;32(1):61–9.** (n=817) | 44% alive at 1 yr; 57% need caregiver assistance. Framing: mortality is a blunt endpoint; "meaningful survival" is the relevant target. | A |
| P7 | **Ludski J et al. Long-term outcomes of critically ill patients requiring prolonged mechanical ventilation: a systematic review. 2023.** (24 studies) | Mortality 57%/69% at 12/48 mo; only 30.2% discharged home; ≤39% psychiatric sequelae. Reinforces P6. | A (review) |

**Introduction arc:** early tracheostomy reliably shortens ventilation and ICU
stay (P2, P3) and reduces VAP (P1) but does not change mortality (P1–P3; P4 the
qualified dissent) — and mortality is itself a blunt endpoint over a population
whose survivors are largely dependent (P5–P7). The paper asks *why the process
moves and the endpoint does not*, and answers it mechanistically.

---

## 2. Model 1 — respiratory mechanics (load axis)

| # | Citation | Use | Grade |
|---|---|---|---|
| M1.1 | **Flevari AG et al. Anaesth Intensive Care 2011;39(3):410–7. PMID 21675060.** | ETT Rohrer coefficients K1/K2, directly in our form. Sensitivity arm. | A |
| M1.2 | **Guttmann J et al. Anesthesiology 1993;79(3):503–13. PMID 8363076.** | ETT + tracheostomy on ONE rig (primary coefficient set); refit to Rohrer form. Reproduces its own 70% trach/ETT check. ⚠ its published "K2" is a power-law exponent, not a quadratic coefficient. | A (data) / B (refit) |
| M1.3 | **Davis K et al. Arch Surg 1999;134(1):59–62. PMID 9927132.** | Apparatus dead space per device; in-vivo WOB ordering extubated > ETT > trach. | A |
| M1.4 | **Nunn JF, Campbell EJM, Peckett BW. J Appl Physiol 1959;14(2):174–6. PMID 13641137.** | Upper-airway dead space 72±32 mL (n=6 cadavers) — the origin of "trach halves dead space." | A (weak) |
| M1.5 | **Chadda K et al. Intensive Care Med 2002;28(12):1761–7. PMID 12447520.** | Only in-vivo non-intubated vs tracheostomised ΔV_D (74 mL) and WOB ratio (1.32). Validation target. | A |
| M1.6 | **Mohr AM et al. J Trauma 2001;51(5):843–8. PMID 11706329.** | No V_D/V_T change after tracheostomy (0.51→0.51). Null validation of the cancellation. | A |
| M1.7 | **Joseph MJ et al. Anaesth Intensive Care 2013;41(2):216–21. PMID 23530788.** | "The myth of dead space": 41%→40%, p=0.75. Null validation. | A |
| M1.8 | **Radford EP. J Appl Physiol 1955;7(4):451–60. PMID 13233138.** | Anatomic dead space ≈2.2 mL/kg. | A |
| M1.9 | **Hart MC et al. J Appl Physiol 1963;18(3):519–22. PMID 31094493.** | Height-based dead-space scaling (better than weight). | A |
| M1.10 | **Pryor 2016; Carter 2013** | Trach 7.0 leads — Pryor unusable (flow^3.5–4), Carter access-blocked. Why TRACH_7_0 is grade X. | X |

## 3. Model 1b — CO₂ kinetics (gas-exchange axis)

Shares M1.3–M1.7 for dead space. Additional:

| # | Citation | Use | Grade |
|---|---|---|---|
| M1b.1 | Alveolar gas equation (standard); constant 0.863 (STPD/BTPS). | PaCO₂ = 0.863·VCO₂/V_A. | C |
| M1b.2 | Body CO₂-store and dissociation physiology (Nunn's/West standard). | Two-compartment τ_body vs τ_lung. S_body reinterpreted per-kg. | C |
| M1b.3 | Critical-illness VCO₂ 300–400 mL/min (metabolic axis). | Metabolic rate as the third severity axis. | C |

## 4. Model 2 — load–capacity bifurcation (tipping-point axis)

| # | Citation | Use | Grade |
|---|---|---|---|
| M2.1 | **Roussos CS, Macklem PT. J Appl Physiol 1977;43(2):189–97. PMID 893274.** | Pdicrit ≈ 0.4, a PURE pressure ratio → correct `u_crit`. Overturns the spec's TTdi-anchored 0.15. | A |
| M2.2 | **Bellemare F, Grassino A. J Appl Physiol 1982;53(5):1190–5. PMID 7174413.** | Tension-time index and its 0.15 threshold — and why it carries a duty-cycle factor `u=L/C` does not. | A |
| M2.3 | **Vassilakopoulos T et al. Am J Respir Crit Care Med 1998;158(2):378–85. PMID 9700110.** | Weaning discrimination Pi/Pimax 0.31 (success) vs 0.46 (fail); MIP 53.8/42.3. Calibrates `s`; validation target. | A |
| M2.4 | **Laghi F et al. J Appl Physiol 1995;79(2):539–46. PMID 7592215.** | 24-h diaphragm recovery curve → `α`, and the biphasic-recovery limitation. | A |
| M2.5 | **Laghi F et al. Am J Respir Crit Care Med 2003;167(2):120–7. PMID 12411288.** | Weaning failure WITHOUT low-frequency fatigue — the challenge to a fatigue model, confronted in Discussion. | A |
| M2.6 | **Bertoni M, Goligher EC et al. Crit Care 2020;24:106. PMID 32204729.** | Normal P_mus 5–10 cmH₂O; P0.1 windows. Load-scale sanity. | C (review) |

## 5. Model 2b — VIDD capacity dynamics (capacity axis)

| # | Citation | Use | Grade |
|---|---|---|---|
| M2b.1 | **Jaber S et al. Am J Respir Crit Care Med 2011;183(3):364–71. PMID 20813887.** | TwPtr −32% at 6 d (force, cmH₂O) → `k_deg` = 0.064/day. Primary anchor. Also force −32% vs CSA −39% (the C∝CSA data point). | A |
| M2b.2 | **Zambon M et al. Crit Care Med 2016;44(7):1347–52. PMID 26992064.** | Atrophy rate by ventilatory mode → the MONOTONIC g(A). Its own conclusion: linear in support, no U. | A |
| M2b.3 | **Lecronier M et al. Ann Intensive Care 2022;12(1):34. PMID 35403916.** | Septic diaphragm force +19% over 4 d (reversible) → contradicts strong-H3; two groups → solves `k_syn`. | A |
| M2b.4 | **Demoule A et al. Am J Respir Crit Care Med 2013;188(2):213–9. PMID 23641946.** | Sepsis = admission OFFSET (−3.74 cmH₂O), no Day1→3 decline → sepsis as C(0) offset, not a drain. | A |
| M2b.5 | **Levine S et al. N Engl J Med 2008;358(13):1327–35. PMID 18367735.** | Rapid disuse atrophy is real (53–57% CSA, 18–69 h). EXCLUDED from rate calibration (window too wide). | A (phenomenon) / X (rate) |
| M2b.6 | **Goligher EC et al. Am J Respir Crit Care Med 2018;197(2):204–13. PMID 28930478.** | U-shape in OUTCOME (TFdi 15–30% shortest ventilation) — supports the `ushape` variant, not the default. | A |
| M2b.7 | **Goligher EC et al. Am J Respir Crit Care Med 2015;192(9):1080–8. PMID 26167730.** | Thickness change by effort; no post-extubation regrowth; function non-monotonic in thickness. | A |
| M2b.8 | **Schepens T et al. Crit Care 2015;19:422. PMID 26639081.** | 9/20/26% loss at 24/48/72 h — first-order (exponential) form check. | B |
| M2b.9 | **Orozco-Levi M et al. Am J Respir Crit Care Med 2001;164(9):1734–9. PMID 11719318.** | Load-induced sarcomere injury exists in humans — a_injury > 0, magnitude unconstrained. | A (existence) / X (magnitude) |
| M2b.10 | **Hermans G et al. Crit Care 2010;14(4):R127. PMID 20594319.** | Logarithmic (vs exponential) decline dissent (n=10). Form-check caveat. | B |
| M2b.11 | **Yamada T et al. Lung 2024.** | Thickness did NOT independently predict MIP (n=109) — contrary evidence for C∝CSA. | A |

## 6. Model 2b — the unmeasured sedation/activity link (A-gap)

| # | Citation | Use | Grade |
|---|---|---|---|
| S1 | **Young D et al. (TracMan). JAMA 2013;309(20):2121–9. PMID 23695482.** | Early trach 5 vs 8 sedation-days; no mortality benefit. RCT support for Link 1 (trach → less sedation). | A |
| S2 | **Nieszkowska A et al. Crit Care Med 2005;33(11):2527–33. PMID 16276177.** | Same-patient before/after: midazolam −84%, fentanyl −92%, heavy-sedation time 7→1 h/day. | A (confounded) |
| S3 | **Trouillet JL et al. Ann Intern Med 2011;154(6):373–83. PMID 21403073.** | RCT: less IV sedation in trach arm under randomisation. De-confounds Link 1. | A |
| S4 | **Meng L et al. Clin Respir J 2016.** | Meta: sedation duration WMD −5.99 d. | A (meta) |
| S5 | *(gap)* No study measures diaphragm activity ETT vs trach. | Link 2 (sedation → activity) UNMEASURED. The A-gap is a sensitivity variable, not a constant. | X |

---

## 7. How the reference base maps to the claims

- **Paradox is real and mortality-null:** P1–P4 (with P4 the honest dissent).
- **Mortality is a blunt endpoint:** P5–P7.
- **Every mechanism number is graded**, and the load axis (M1/M2) is anchored on
  A-grade primary data; the capacity axis (M2b) has one A-grade rate anchor
  (Jaber) but rests critically on one X-grade unmeasured quantity (the A-gap, S5).
- **Three independent MIP anchors converge:** Vassilakopoulos 42.3/53.8 (M2.3),
  Jubran 41.3 (P5), and the ATS/ERS norms — all in the 40–55 cmH₂O band the
  model's capacity grid uses.
- **Two hypotheses were overturned by the literature**, not by us: H3 (M2b.3,
  M2b.4) and the degradation U-shape (M2b.2 vs M2b.6). Both spec versions are
  retained as switchable modes.
