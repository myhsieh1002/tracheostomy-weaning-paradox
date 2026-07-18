# Parameter Calibration & Provenance

Every number in `config/params_m1.json` and `config/params_m2.json` is
recorded here with its source and evidence grade. Nothing in this file is
invented; where the literature does not supply a value, that is stated
plainly and the parameter is carried as a sensitivity variable rather than
a constant.

Evidence grades used throughout:

| Grade | Meaning |
|---|---|
| **A** | Directly reported primary measurement |
| **B** | Derived or interpolated here from a primary source (arithmetic shown) |
| **C** | Secondary assertion / textbook — traceable but not itself a measurement |
| **X** | Not found in the literature; remains `TODO-verify` |

---

## 1. Three findings that changed the model

The calibration pass did not merely fill in numbers — it overturned two
parameter choices in the build specs and confirmed a third concern. These
are recorded first because they change what the model *says*, not just what
it computes.

### 1.1 `u_crit = 0.15` is an invalid anchoring → **0.4**

The Model 2 spec anchors the fatigue threshold to the Bellemare–Grassino
tension-time index, TTdi ≈ 0.15. That anchoring does not hold, because

```
TTdi = (Pdi/Pdimax) × (TI/Ttot)        ← carries a duty-cycle factor
u    =  L/C                             ← the model's utilisation: no duty cycle
```

The two are not the same quantity. Using 0.15 for `u` would understate the
threshold by roughly a factor of three.

The correct anchor is **Roussos & Macklem 1977**, who measured the
sustainable pressure ratio *directly* as a pure Pdi/Pdimax with no timing
term — structurally identical to `u = L/C` — and found **Pdicrit ≈ 0.4**.

Three independent lines converge:

| Route | Value | Grade |
|---|---|---|
| Roussos & Macklem 1977 — Pdicrit measured directly | 0.40 | A |
| Bellemare & Grassino TTdi 0.15 ÷ (TI/Ttot ≈ 0.34) | 0.44 | B |
| Vassilakopoulos 1998 — Pi/Pimax 0.46 (fail) vs 0.31 (success) | ~0.4 | A |

**Consequence for the code:** `u_crit = 0.4` pairs with `L` defined as the
**mean inspiratory** P_mus. The model also exposes a duty-weighted `L`,
which reproduces TTdi and must then be paired with `u_crit = 0.15`. The two
parameterisations are consistent descriptions of the same physiology; they
must not be mixed. `lc.utilisation` enforces the pairing.

### 1.2 The dead-space mechanism behind H3 largely cancels

The Model 1 spec treats ΔV_D ≈ 70 mL as the mechanism behind H3. That
figure is the upper-airway volume bypassed **relative to the
non-intubated state** — not relative to an ETT.

An ETT and a tracheostomy tube are both sited with the tip in the
mid-trachea. **Both bypass the same anatomic region.** In the ETT-vs-trach
contrast the bypass term cancels, and the residual difference is only the
apparatus volume: ~15–24 mL (ETT) against ~5–6 mL (trach), i.e. **10–18 mL**
against a tidal volume of ~450 mL.

**This prediction is externally validated.** Two studies measured dead space
before and after tracheostomy in already-intubated patients and found
nothing:

| Study | n | Result |
|---|---|---|
| Mohr 2001 (PMID 11706329) | 42 | V_D/V_T 0.51 ± 0.10 → 0.51 ± 0.11 — no change |
| Joseph 2013 (PMID 23530788) | 24 | Dead-space fraction 41 ± 12.6% → 40 ± 14.6%, p = 0.75 |

Joseph's title says it outright: *"Tracheostomy in the critically ill: the
myth of dead space."* The model predicts a null and the clinical
measurements report a null. This is Model 1's first external check, and it
strengthens the programme thesis rather than weakening it: the dead-space
axis is one more on which the device moves less than assumed.

**Do not mis-cite the negative studies.** Mohr and Joseph compare ETT→trach,
so they bear on the apparatus term only, *not* on ΔV_D from upper-airway
bypass. Only **Chadda 2002** compares the non-intubated state with the
tracheostomised state, and it is the sole in vivo estimate of ΔV_D itself.

### 1.3 The evidence base for ΔV_D is far thinner than its ubiquity implies

The canonical "tracheostomy halves the dead space" reduces to a single
arithmetic ratio, 72 ÷ (72 + 66) = 52%, from **Nunn, Campbell & Peckett
1959** — **six cadavers**, SD ±32 mL on a mean of 72 mL (CV ≈ 44%). Every
secondary repetition traced during calibration leads back to that one
study. The only corroboration is Chadda 2002 (n = 9, in vivo), which
independently gives 74 mL.

Two further problems:

- Nunn's "extrathoracic" volume **over-counts** what a stoma bypasses: it
  includes the cervical trachea, part of which remains in circuit distal to
  the stoma. 72 mL is an **upper bound**, not an estimate.
- Within-subject lability is of the same order as the quantity itself: jaw
  depression + neck flexion moved it by −31.4 mL, protrusion + extension by
  +39.7 mL (Nunn 1959).

**ΔV_D is therefore carried as a sensitivity parameter (40 / 72 / 90 mL),
never as a constant.** The claim "up to 150 mL" appears in secondary
sources but **no primary source was found**; it is not used.

---

## 2. Model 1 parameters

### 2.1 Dead space

| Parameter | Value | Source | Grade |
|---|---|---|---|
| `Vd_anatomic_mL` | 150 (≈2.2 mL/kg at 68 kg) | Radford 1955, PMID 13233138 | A |
| — preferred if height known | `V_D = 7.585 × Ht(cm)^2.363 × 1e-4` | Hart 1963, PMID 31094493 (n=73, r=0.917) | A |
| `Vd_upper_bypassed_mL` | 72 (range **40 / 72 / 90**) | Nunn 1959, PMID 13641137 (n=6 cadavers); Chadda 2002, PMID 12447520 (n=9, 74 mL) | A, but weak |

Note: "150 mL" is a pedagogic round number; Pierson sources it to *Egan's
Fundamentals* (grade C). Radford or Hart are the citable primaries.

### 2.2 Apparatus dead space

All from **Davis K et al., Arch Surg 1999;134(1):59–62** (PMID 9927132),
tabulated in Pierson 2005 Table 2. Direct internal-volume measurement.

| Device | Volume (mL) | Grade | Note |
|---|---|---|---|
| ETT ID 7.0, 34.5 cm | 15 | A | published |
| ETT ID 7.5 | ~18 | B | interpolated 7.0→8.5 |
| ETT ID 8.0 | ~21 | B | interpolated 7.0→8.5 |
| ETT ID 8.5, 36.5 cm | 24 | A | published |
| Trach ID 7.0, 12 cm | 5 | A | published |
| Trach ID 8.0 | ~5.5 | B | interpolated 7.0→8.5 |
| Trach ID 8.5, 12 cm | 6 | A | published |

Only ID 7.0 and 8.5 are published; 7.5 and 8.0 are interpolated here.
Manufacturer datasheets publish ID/OD/length but **not** internal volume.

⚠️ **Length consistency:** Davis measured *uncut* tubes (ETT 34.5 cm). The
build spec's device table assumed 27 cm. Volume and resistance both scale
with length, so the device geometry must be self-consistent — see
`§2.4 Open issues`.

Davis 1999 also reports, in vivo (n=20): WOB/min 8.9 ± 2.9 → 6.6 ± 1.4 J/min
after tracheostomy (p < 0.04); airway resistance 9.4 ± 4.1 → 6.3 ± 4.5
cmH₂O/(L/s) (p < 0.07). These are direct validation targets for Mode B.

### 2.3 Rohrer coefficients

Two published sources supply K1/K2, and the choice between them is not a
matter of taste.

**The naming trap.** Guttmann's own published "K1/K2" are a **power law**,
ΔP = K1·V̇^K2, in which K2 is a dimensionless **exponent** (1.94–2.03) — not
the quadratic coefficient of our Rohrer form. Same letters, different
objects. Citing his K2 as a quadratic coefficient would be a straight error.
His true Rohrer fits were **never published**: the paper states (p.510) that
"coefficient tables based on quadratic approximation are available from the
authors on request". The values below are therefore a **refit** of his
published curves over 0.05–1.5 L/s (rms < 0.12 cmH₂O, well inside his own
0.5 cmH₂O inter-model spread).

| Source | Form | Grade | Devices |
|---|---|---|---|
| **Flevari 2011**, PMID 21675060 | `R_ETT = K1 + K2·V̇`, **algebraically identical** to ours — directly usable, no conversion | A | ETT only (Portex Blue Line) |
| **Guttmann 1993**, PMID 8363076 | power law; **refitted here** to Rohrer form | B | ETT **and** trach, one rig |

#### Which set is primary, and why it matters

**`guttmann_within_rig` is primary.** The scientific question *is* the
ETT-vs-tracheostomy contrast, so within-rig consistency outranks the grade
of any single number. Guttmann is the only source measuring both device
classes on one rig, so any rig-specific bias — swivel connector, artificial
trachea, flow waveform — is **common to both arms and cancels from the
contrast**. Pairing Flevari's ETT with Guttmann's trach would inject a
systematic between-rig difference (1.10–1.47×, rising with ID, consistent
with Guttmann's 15 mm connector) directly into the quantity of interest.

**The choice validates itself.** Guttmann-within-rig reproduces his own
published check — ETT 8.0 = 6.55 vs trach 8.0 = 4.63 cmH₂O at 1 L/s, ratio
**0.707** against his reported "trach ≈ 70% of ETT". The mixed-source
pairing gives 0.84 and does not reproduce it.

**The bias runs the safe way.** Guttmann gives the *larger* ETT–trach gap,
hence a larger ΔL_device and a *wider* rescue window — i.e. it is the choice
most favourable to tracheostomy. Using it makes the programme's "the device
moves little" conclusion conservative rather than assumed. `flevari_ett` is
retained as the sensitivity arm; swapping sets moves f_device by ~0.01.

#### Values

| Device | Guttmann (primary) K1 / K2 | Flevari (sensitivity) K1 / K2 | Grade |
|---|---|---|---|
| ETT 7.0 | 0.10 / 11.01 | 0.97 / 9.17 | B / A |
| ETT 7.5 | 0.30 / 8.09 | 0.92 / 6.01 | B / A |
| ETT 8.0 | 0.36 / 6.19 | 0.83 / 4.65 | B / A |
| Trach 8.0 (9.0 cm) | 0.0 / 4.63 | — (no trach measured) | B |
| **Trach 7.0** | **0.0 / 7.2** | — | **X — UNVERIFIED** |

**K1 ≈ 0 is not a fitting artefact.** Guttmann states explicitly that "there
is no need for a linear term ... because flow separation and turbulence
prevail", and his Rohrer fit was equally good (rms 0.23 vs 0.22). Flevari's
authors likewise disclaim their K1 as having "only mathematical value" (the
intercept lies outside their tested flow range). The unconstrained trach fit
returns −0.12; it is **clamped to 0** rather than shipping a negative value.
The Sobol analysis independently confirms this: `ST(K1_scale) ≈ 0.001`.

**Trach 7.0 has no publishable values.** Its K2 is extrapolated here from
Guttmann's trach K1 ~ D^−3.52 trend (r² = 0.9999) **beyond his measured
range** — his smallest trach is 8.0. Pryor 2016 (Respir Care 61(5):607–14)
does report a measured 7.0, but its ΔP scales as V̇^3.5–4.0, which turbulence
cannot produce; forcing Rohrer onto it gives K1 = −8.4 (negative ΔP below
0.5 L/s), and it is a *double-lumen* tube (inner cannula ~3× resistance),
not a bare 7.0. **Carter 2013** (*Anaesthesia*, 6–10 mm trach bench) is the
best remaining lead; it needs institutional access. Until then TRACH_7_0 is
grade X, drawn hatched, and excluded from every verdict.

### 2.4 Metabolic rate is a severity axis — the spec omitted it

`target_VA = 4.2 L/min` corresponds to a **normal** VCO₂ of ~200 mL/min.
Critical illness runs 300–400 mL/min (Model 1b's own spec), and V_A scales
with VCO₂ at fixed PaCO₂, giving V_A ≈ 6.3–8.4.

This is not a detail. Held at 4.2, the model's **maximum load anywhere on
the disease grid is 13.6 cmH₂O** — it cannot reach the 19.5 cmH₂O measured
in failing weaning patients (Vassilakopoulos 1998). The rescue window then
falls entirely *below* the observed MIP range, producing a spurious "the
device never matters for real patients" result that is an artefact of
assuming normal metabolism. Sobol subsequently ranks `target_VA` the
**largest single driver of total load** (ST = 0.56).

Base case **6.3** (VCO₂ 300); grid [4.2, 6.3, 8.4]; grade B (derived from
the VCO₂ proportionality, not measured directly).

### 2.5 The non-intubated arm needs upper-airway resistance

exp3's third arm is the intact airway. Modelling it as "no tube, therefore
no resistance" is wrong and **inverts the comparison against Davis 1999**,
who found WOB *higher* extubated (1.2 J/L) than via an ETT (0.81) or trach
(0.77).

The natural upper airway resists. **Chadda 2002** measured it in vivo and
reported it "did not differ from in vitro tracheostomy tube resistance" —
which is both the licence for using the trach's coefficients for the
`NATIVE_UPPER_AIRWAY` pseudo-device, and the reason Chadda's mouth-vs-trach
difference is attributable to dead space rather than resistance. With this
in place the model reproduces Davis's ordering and Chadda's WOB ratio to 3%.

### 2.6 Open issues

- **Device length is heterogeneous across sources.** Guttmann's ETTs are
  30.8–32.3 cm and trachs 9.0–10.5 cm; Davis's dead-space tubes are 34.5 cm
  ETT and 12 cm trach. Resistance and dead space therefore come from
  different tubes. Lengths are recorded per their own source rather than
  averaged into a fiction. This is a real limitation.
- **Connector/HME/T-piece volume** is not in Davis's numbers and is not
  modelled; Pierson notes it further erodes the trach's apparatus advantage.
- **Sub-stomal cervical trachea** — the segment distal to the stoma that
  remains in circuit — is **unquantified in any source found**. It is the
  main reason the ΔV_D base case should sit *below* Nunn's 72 mL. Worth
  stating in the manuscript. (In practice ΔV_D cancels from the ETT-vs-trach
  contrast and Sobol gives it ST ≈ 0.006, so nothing downstream turns on it.)

---

## 3. Model 2 parameters

| Parameter | Value | Range | Source | Grade |
|---|---|---|---|---|
| `u_crit` | **0.4** | 0.35–0.5 | Roussos & Macklem 1977, PMID 893274 — Pdicrit measured as pure Pdi/Pdimax | A |
| `C_max` grid | **[20,30,40,50,60,70]** cmH₂O | — | Vassilakopoulos 1998, PMID 9700110: MIP 42.3 ± 12.7 (fail) / 53.8 ± 15.1 (success) | A |
| `alpha` | 0.15 /hr | 0.1–0.3 | Laghi 1995, PMID 7592215 — 24 h recovery curve | B / semi-constrained |
| `beta` | 1.0 /hr | 0.5–2 | **No published value exists** — inferred from a 35% capacity fall over an assumed 15–45 min Tlim | X |
| `s` | 25 | — | Not addressed in the literature; constrained only indirectly (Tlim > 45 min at u ≈ u_crit) | X |
| `L` scale | 5–25 cmH₂O; 15–20 in weaning | — | Vassilakopoulos 1998 (derived: Pi = Pi/Pimax × Pimax) | B |

**`C_max` grid correction.** The spec's [40, 50, 60, 70, 80] is skewed high:
70–80 exceeds the weaning population, and the grid has **no resolution over
20–40 cmH₂O — precisely where weaning failure lives** (failure-group mean
−1 SD ≈ 29.6; ICU-acquired weakness routinely gives MIP 20–40).

**`alpha` is a model-structure finding, not just a parameter.** Laghi's
recovery is **not single-exponential**: a fast phase (τ ≈ 3–5 h, α ≈ 0.2–0.35
/hr) to ~8 h, then a residual ~16% deficit still unresolved at 24 h. The
model's `α(C_max − C)` term forces complete exponential recovery, so **no
single α reproduces the data** — fitting the early phase overshoots the 24 h
endpoint; fitting 24 h misses the fast phase. Reported as a limitation
rather than hidden.

| Time | Pditw (cmH₂O) | % baseline | % of deficit recovered |
|---|---|---|---|
| baseline | 38.9 ± 1.1 | 100% | — |
| +10 min | 25.1 ± 0.6 | 64.5% | 0% |
| +1 h | 27.6 ± 0.9 | 71.0% | 18% |
| +8 h | 31.6 ± 1.1 | 81.2% | 47% |
| +24 h | 32.7 ± 1.2 | 84.1% | 55% (still < baseline, p < 0.01) |

**`beta` is unconstrained.** No publication reports a fatigue rate constant
or anything dimensionally equivalent. It must not be presented as
literature-derived.

### 3.1 Validation target — Vassilakopoulos 1998

The model reproduces the actual clinical discrimination with `u_crit = 0.4`:

| | Pi (cmH₂O) | Pimax (cmH₂O) | u = L/C | Predicted | Observed |
|---|---|---|---|---|---|
| Failure | 19.5 | 42.3 | **0.46** | > u_crit → fail | fail |
| Success | 16.7 | 53.8 | **0.31** | < u_crit → wean | success |

Implemented as a regression test (`tests/tValidation.m`).

### 3.2 The Laghi 2003 challenge — must be addressed head-on

**Laghi F et al. 2003 (PMID 12411288)** found that weaning failure was **not**
accompanied by low-frequency diaphragm fatigue: twitch Pdi was essentially
unchanged after failed SBTs (8.9 → 9.4 cmH₂O). Distress mandated
reinstitution of ventilation *before* fatigue developed, and **weakness —
not fatigue — dominated**.

This is a direct challenge to a load-vs-capacity *fatigue* model of weaning
failure: the mechanism the model describes may be pre-empted in clinical
reality. It should be confronted in the Discussion, not omitted. One
defensible reading is that the model's low-C attractor represents a state
the patient is removed from *before* reaching — i.e. the bifurcation governs
the decision to stop, not the physiological collapse itself.

Related: `C_max` itself decays on a **days** timescale — an order of
magnitude slower than α — via ventilator-induced diaphragmatic dysfunction
(Levine 2008, PMID 18367735: 18–69 h inactivity → 53–57% myofiber CSA loss;
Goligher 2015, PMID 26167730: 44% lost >10% diaphragm thickness in week 1).
This supports the M2c slow-variable extension.

---

## 4. Outstanding

- [x] Rohrer `K1`/`K2` — done for four of five devices (§2.3)
- [ ] **`TRACH_7_0` K2 — no publishable value.** Blocked on Carter 2013
      (*Anaesthesia*), which needs institutional access. Currently grade X
      and excluded from all conclusions.
- [ ] Device length heterogeneity across sources (§2.6) — documented, not
      resolved; would need a single bench study measuring both resistance
      and internal volume on the same tubes.
- [ ] Connector/HME dead space — not modelled.
- [ ] `beta` — the Laghi-derived anchor is an interpretation (§3); `s` is
      only bounded from below (§10.5 of the M2 spec). Both remain
      sensitivity parameters.
- [ ] `target_VA` — derived from the VCO₂ proportionality (§2.4), not
      measured. Model 1b would supply it properly.

---

## 5. Do-not-cite list

Claims encountered during calibration that could **not** be traced to a
primary source. Recorded so they are not reintroduced later:

| Claim | Status |
|---|---|
| `Tlim = 0.1 × TTdi^-3.6` | Widely repeated; **not verifiable** against any accessible primary source. Do not cite. |
| ΔV_D "up to 150 mL" | **No primary source found.** Appears to conflate total V_D with the bypassed portion. |
| "150 mL anatomic dead space" attributed to a primary | Textbook round number (Egan's). Cite Radford 1955 or Hart 1963 instead. |
