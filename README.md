# The Tracheostomy Weaning–Survival Paradox: A Computational Account

MATLAB implementation of four coupled models that together explain why
tracheostomy reliably shortens mechanical ventilation yet does not change
survival.

**The paradox.** Meta-analyses agree: tracheostomy shortens ventilation time
and ICU stay (it helps weaning) but does not move long-term survival. Why does
an intervention that improves the *process* leave the *endpoint* untouched?

**The thesis.** The tracheostomy moves one thing — the respiratory *load* —
and it moves it less than assumed and least where disease is worst. It can flip
a weaning outcome only inside a narrow band of neuromuscular capacity whose
position depends on lung mechanics, gas exchange, and metabolic rate at once.
It may also raise muscle *capacity* by enabling lighter sedation — but by how
much has never been measured. What it cannot touch is the disease-driven
elastic load and the competing causes of death. So it optimises the process
and averages to nothing on the endpoint.

Four models make that argument quantitative, each on its own axis, sharing one
repository and one calibration file.

---

## The four models, and how they couple

```
                          target_VA (metabolic rate)
                                   │
   ┌───────────────┐   V_E, V_T    ▼   P_mus, f_device   ┌───────────────┐
   │  Model 1b     │ ───────────▶ Model 1 ─────────────▶ │   Model 2     │
   │  CO₂ kinetics │              respiratory            │  load–capacity│
   │  (+co2)       │              mechanics (+wob)       │  bifurcation  │
   └───────────────┘                                     │   (+lc)       │
      gas-exchange axis          load axis                └───────┬───────┘
                                                        C_max(t)  │  fold
                                                          ┌───────▼───────┐
                                                          │   Model 2b    │
                                                          │  VIDD capacity│
                                                          │   (+vidd)     │
                                                          └───────────────┘
                                                          capacity axis
```

| Model | Axis | Question | Thesis figure |
|---|---|---|---|
| **1** `+wob` | respiratory load | What share of the work of breathing is the device — and how is it diluted by disease? | `fig_device_fraction_heatmap` |
| **1b** `+co2` | gas exchange | How much ventilatory demand does the device's dead space save, and how is *that* diluted? | `fig_co2_dilution_grid` |
| **2** `+lc` | tipping point | Can removing that load flip the weaning outcome, and for whom? | `fig_rescue_window` |
| **2b** `+vidd` | muscle capacity | Does the device also raise capacity by unloading sedation, and does that dominate? | `fig_vidd_dynamic_rescue` |

**The couplings are real, not narrative.** Loads are never hard-coded:

- **1b → 1** — `co2.coupling` computes the ventilation needed to hold PaCO₂,
  then hands it to `wob.simulateModeB` via its `Ventilation=` port. Model 1
  computes the *mechanics* of delivering it; Model 1b decides *how much*.
- **1 → 2** — `lc.coupling` calls `wob.simulateEffort` and nothing else. The
  load `L` is Model 1's `P_mus`, apportioned by `f_device`.
- **2b → 2** — `vidd.couplingToModel2` supplies `C_max(t)` and drives Model 2's
  fold through time via `lc.rescueWindowTrajectory`, which calls the named API
  `lc.rescueWindow` (Model 2 spec §2.5b). The fold is solved in exactly one
  place; 2b never re-derives it.

Everything shares `docs/CALIBRATION.md` and the same evidence grades.

---

## Reproduce

The MATLAB MCP server attaches to a *shared* desktop instance. To avoid
clobbering another session's workspace, run headless — an isolated process:

```bash
matlab -batch "runAll"          # tests + 19 experiments + 4 summaries (~4 min)
matlab -batch "runTests('all')" # 121 tests only (~10 s + startup)
```

```matlab
runAll(Only="m2b")        % one group: "m1" | "m1b" | "m2" | "m2b" | "all"
runAll(Tests=false)       % skip the suite
runTests("m1b")           % one suite
```

Outputs land in `results/`: figures (300 dpi PNG + vector PDF), tables (CSV),
and four `summary*.md` files **regenerated from the models at write time** — no
number in them is transcribed, so a summary can never drift from the code.

Requires MATLAB R2025b with the Statistics and Machine Learning Toolbox
(`sobolset`) and Optimization Toolbox (`fmincon`, `fzero`). No network access.

---

## Layout

```
config/    params_m1.json  params_m2.json  params_co2.json  params_vidd.json
+wob/      Model 1  — devices, waveform, mechanics, metrics, Sobol
+co2/      Model 1b — dead space, alveolar gas eqn, two-compartment CO₂, coupling
+lc/       Model 2  — load-capacity dynamics, folds, rescue API, coupling
+vidd/     Model 2b — VIDD capacity dynamics, U/monotone degradation, coupling
+viz/      shared figure style
experiments/   exp1-5, co2_exp1-4, m2_exp1-5, vidd_exp1-5   (19)
tests/         tModel1 (29)  tModel2 (34)  tModel1b (30)  tModel2b (28)  = 121
docs/          CALIBRATION.md — provenance and evidence grade for every number
results/       figures, tables, summary.md / _1b / _m2 / _2b
runAll.m  runTests.m  writeSummary.m
```

---

## What each model finds

### Model 1 — the device is a small, disease-diluted share of the load
`f_device = WOB_device / WOB_total` for an ETT falls monotonically from **0.36**
in a healthy lung to **0.11** in a severe one, in both directions across the
(C_rs × R_aw) grid. The device is the smallest part of the problem exactly
where the problem is largest. **H3 correction:** the dead-space benefit the
spec expected (~72 mL) *cancels* between ETT and trach — both bypass the same
airway — leaving only ~12 mL of apparatus volume, which two negative clinical
studies confirm.

### Model 1b — the gas-exchange benefit is diluted the same way
The ventilatory saving of ETT→trach is **3.7%** of demand in the mildest case
and **1.2%** in the sickest — the gas-exchange analogue of Model 1's dilution,
on a new axis. The device step is only **19%** of the non-intubated→trach step,
so the 72 mL literature describes the wrong comparison for a clinician holding
an ETT. The two-compartment model earns its keep: the half-time to a new PaCO₂
after tracheostomy is **~55 min**, set by the body CO₂ stores (τ ≈ 35 min), not
the lung (τ ≈ 0.6 min) — a single compartment would be wrong by 100×. And the
dead-space benefit lands in *effort*, not blood gas — which is why studies
measuring PaCO₂ find nothing, and why that null is not evidence of no benefit.

### Model 2 — the rescue window is narrow and slides, it does not shrink
The load–capacity system is bistable with a saddle-node fold. Removing the
device's load flips a weaning outcome only inside a band of capacity **~3 cmH₂O
wide (≈7% of the clinical C_max range)**, whose position depends on compliance,
airway resistance *and* metabolic rate. An unselected trial puts few patients
in it and averages to nothing — the paradox, quantified. The system is
scale-invariant in capacity, which gives the rescue window a closed form and
means bistability never disappears with falling capacity (a *wedge*, not a
cusp). Reproduces Vassilakopoulos 1998's actual wean/fail discrimination.

### Model 2b — the device may also raise capacity, but the size is unmeasured
Capacity is reversible and activity-sensitive (`k_deg` from Jaber 2010, `k_syn`
from Lecronier 2022's two groups). Put in Model 2's units, the tracheostomy's
capacity effect *can* exceed its load effect — the two channels cross at an
activity gap of only **~0.04** — but the magnitude rides on the ETT→trach
increase in diaphragm activity, which **has never been measured**. At a small,
sedation-plausible gap the channels are comparable (~2×); the headline "5×"
needs the assumed doubling and is reported as an upper bound.

**This reframes the paradox.** If capacity is reversible and the device helps
it, "capacity is beyond rescue" cannot be the explanation. The account shifts
to what the models show: a narrow load-axis window, a capacity benefit
contingent on sedation practice actually changing, and — still to come — the
device's costs and competing causes of death.

---

## Deviations from the build specs

Each changed a *result*, not just an implementation. Full reasoning and
citations in [docs/CALIBRATION.md](docs/CALIBRATION.md); the specs carry
revision logs.

**Model 1 / 2**

| Spec said | Model does | Why |
|---|---|---|
| `u_crit = 0.15` (tension-time index) | **0.4** | TTdi carries a duty-cycle factor `u = L/C` does not; 0.15 understates ~3×. Anchor: Roussos & Macklem 1977. With 0.15 the model fails the patients who actually weaned. |
| ΔV_D ≈ 70 mL drives H3 | **it cancels** | Both tubes bypass the same airway; only ~12 mL apparatus differs. Confirmed null by Mohr 2001, Joseph 2013. |
| pseudo-arclength continuation | **closed form** | The branch inverts exactly; folds are analytic. Cannot jump branches. |
| cusp in (L, C_max) | **a wedge** | Scale invariance → fold loci are lines through the origin (residual 1e-14). |
| window shrinks with C_max | **it slides** | width = ΔL/l_high, and ΔL is near-constant across severity. |
| severity = (C_rs, R_aw) | **+ metabolic rate** | At normal V_A the model can't reach observed failure loads; Sobol ranks V_A the top driver. |
| `β`, `s` placeholders | **calibrated** | Laghi 1995 (β) and Vassilakopoulos 1998 (s). |

**Model 1b / 2b**

| Spec said | Model does | Why |
|---|---|---|
| alveolar dead space ∝ V_T | **∝ (V_T − V_D_series)** | As a fraction of V_T the device vanishes from ventilatory demand entirely, making H3 vacuous. |
| `S_body = 15` (no unit) | **15 mL/(kg·mmHg)** | Absolute reading gives a 0.5 min body time constant, contradicting the spec's own "minutes to hours". |
| H3: sepsis is an activity-independent drain | **a reversible offset** | Lecronier 2022: septic force *rose* 19% over 4 d. Demoule 2013: an admission offset, no Day1→Day3 decline. |
| degradation g(A) is U-shaped | **monotonic (default)** | Zambon 2016, the only study of atrophy *rate* by support, is monotone. The U exists in outcome (Goligher 2018), not in dC/dt. |
| `k_deg` from Levine 2008 | **from Jaber 2010** | Levine's 18–69 h window implies a rate 4–17× every other study; it was the obvious anchor and would have been wrong by an order of magnitude. |

Both spec versions of the overturned hypotheses are kept as switchable modes
(`g_mode`, `d_mode`) so the contrast is a result, not a silent choice. And
across the board: `params.json`, not `.yaml` — MATLAB's `readstruct` supports
json/xml only.

---

## External validations

Checks, not fits, except the one marked as a calibration target.

| check | source | observed | model |
|---|---|---|---|
| ETT vs trach resistance, matched ID 8.0, same rig | Guttmann 1993 | trach 70% of ETT | **70.7%** |
| WOB ratio, non-intubated vs tracheostomised | Chadda 2002 (n=9) | 1.32 | **1.36** |
| ΔV_D, non-intubated vs tracheostomised | Chadda 2002 | 74 mL | **66 mL** |
| WOB ordering extubated > ETT > trach | Davis 1999 | that order | **that order** |
| dead-space change ETT→trach | Mohr 2001, Joseph 2013 | null | **null** |
| diaphragm force loss, 6 d controlled MV | Jaber 2010 | −32% | **−31%** |
| sepsis capacity offset | Lecronier 2022 vs Demoule 2013 | 33% vs 36% (two methods) | **agree to 1.3 pp** |
| weaning discrimination, Pi/Pimax 0.31 vs 0.46 | Vassilakopoulos 1998 | wean / fail | **wean / fail** *(calibration target)* |

---

## Honest boundaries

- **Individual `f_device` values are not precise.** Sobol ranks `K2_scale`
  second (ST ≈ 0.29); K2 is a grade-B refit of Guttmann's curves. The *trend*
  and *order of magnitude* hold; a decimal place does not.
- **The A-gap (Model 2b) is unmeasured and load-bearing.** Diaphragm activity
  has never been compared ETT vs trach. Direction is RCT-supported; magnitude
  is not. This is the single most valuable future measurement for the programme.
- **`s` (Model 2) is bounded only from below.** `TRACH_7_0` has no published
  coefficients (grade X, excluded from every verdict).
- **Declared assumptions:** capacity ∝ cross-sectional area (contrary evidence:
  Yamada 2024); twitch-pressure rates transferred to the MIP scale as fractional
  changes; the support→A mapping.
- **Laghi 2003 challenges the fatigue premise** (weaning failure without
  diaphragm fatigue) and is confronted in `summary_m2.md`, not omitted.

---

## Not yet built

The programme spine names two more pillars, both outside what this MATLAB
environment can verify:

- **Model 3 — cuff–trachea contact / FSI** (the device's *cost*: ischaemia,
  stenosis). Needs contact mechanics, which MATLAB's PDE Toolbox does not
  provide; its spec targets FEniCSx or Abaqus.
- **Study C — target-trial emulation + causal ML** (who actually benefits).
  Needs MIMIC-IV credentialed data access.

Together with the device's costs and competing mortality, these are where the
"why no survival benefit" account is completed. The four models here establish
that the *process* benefit is real, small, disease-diluted on three axes, and
contingent on one unmeasured number.
```
