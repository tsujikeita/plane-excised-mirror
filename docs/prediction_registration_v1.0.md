# A registered prediction: the plane-excised mirror anti-symmetry in future low-noise CMB polarization
## plane-excised-mirror repository — prediction registration v1.0 (2026-08-15)

**Purpose.** We register, in advance of any suitable data, a falsifiable polarization
prediction derived from the temperature characterization of the plane-excised mirror
anti-symmetry (companion paper; Tsuji 2026c). This document freezes the axis, the
statistic, the hypothesis benchmarks, and the decision rule. It is committed publicly
so that future low-noise polarization data (e.g. LiteBIRD) can adjudicate the origin
of the anomaly without analyst freedom.

**Provenance.** The temperature findings this prediction rests on: the mirror
anti-symmetry of the CMB about a single axis deepens to p = (1.4–5)×10⁻⁴ when the
Galactic plane is excised (10⁴–10⁵-realization ensembles); the axis agrees to pixel
resolution (≤3.7°) across 10 published maps × 6 processing configurations; the carrier
is ℓ = 2–4, whose mirror-symmetric component is suppressed to 0.09–0.11 of the
Gaussian expectation, frequency-independently (70–143 GHz, ±10%); half-mission split
systematics account for ≤ 0.02% of the missing symmetric power.

## 1. Frozen quantities

- **Axis** (mirror-plane normal, ± identified):
  HEALPix N_side=16 RING pixel **1134**, Galactic (l, b) = **(81.56°, +14.48°)**
  — equivalently (261.56°, −14.48°); consistent with the published full-sky
  temperature axis. Frozen from the unit-vector mean of the PR3 axes at the
  N_side=16 / common-mask configuration (rule recorded in the companion repository).
- **Statistic**: S⁺_E — the mask-aware, fixed-axis mirror-symmetric power of the
  E-mode scalar map at N_side = 16, band-limited to ℓ ≤ 32 with the Planck
  degrade convention (Gaussian 640′ × pixel window), pixel pairs restricted to
  both-valid pairs under the analysis mask, no monopole/dipole removal.
  Implementation frozen in `src/plane_mirror.py::FixedAxisMirror` (this repository).
- **Tail**: lower (suppressed symmetric power = anti-symmetry).
- **Mask rule** (dataset-dependent, rule frozen): the mission's recommended
  polarization analysis mask, degraded to N_side = 16 by thresholding at 0.9;
  the p-value shall additionally be reported for threshold 0.5 as a robustness line.
- **Null**: the mission's end-to-end signal+noise simulations (≥ 400 realizations
  preferred); Gaussian ΛCDM E-mode ensembles as a stated secondary null.

## 2. Hypothesis benchmarks (computed 2026-08-15, cosmic-variance limit)

- **H_corr (Gaussian sky; anomaly shared only via TE correlation):** expected shift
  z ≈ 0.5 — *not* detectable even noise-free. A null polarization result is
  therefore fully consistent with H_corr and does not disfavour the temperature
  anomaly itself.
- **H_geo (geometric origin, e.g. topology/anisotropy: E inherits the ℓ = 2–4
  symmetric-power suppression at the temperature rate 0.10):** expected z ≈ 1.7 in
  the cosmic-variance limit at f_sky ≈ 0.72. With CV-limited full-sky E at
  ℓ = 2–8 (LiteBIRD-class), the test discriminates H_geo at the ≥ 90% detection
  probability level for p < 0.05.

## 3. Decision rule (frozen)

Compute the single p-value of S⁺_E at the frozen axis on the frozen-rule mask
against the mission simulations. No axis scan, no configuration scan, no statistic
selection. Interpretation: p < 0.01 supports a geometric (mirror-anti-symmetric)
component in polarization; p ≥ 0.01 places an upper limit on the H_geo suppression
factor (to be reported with the measured value and its binomial interval, floors as
bounds, per the tail-lottery reporting checklist).

## 4. Why Planck cannot settle this

For Planck PR4, the NPIPE large-angle polarization transfer suppression and the
low-ℓ E noise reduce both benchmarks below testability (the H_corr benchmark is
already below testability in the noise-free limit). This registration therefore
targets the next generation of large-scale polarization data.

## Change log
- v1.0 (2026-08-15): frozen. **Public commit `63aa853`**,
  repository `github.com/tsujikeita/plane-excised-mirror`,
  path `docs/prediction_registration_v1.0.md`. Permanent link:
  https://github.com/tsujikeita/plane-excised-mirror/blob/63aa853/docs/prediction_registration_v1.0.md
  (This change-log entry is the only amendment after freezing.)
- Correction 1 (2026-08-15): two benchmark descriptions in §2 are corrected;
  no frozen quantity (axis, statistic, mask rule, null, decision rule) changes.
  (i) The H_geo benchmark sky fraction reads f_sky ≈ 0.72; the mask actually used
  is the common Galactic mask at N_side=16, f_sky = 0.65 (0.72 was the N_side=64
  value of a companion analysis). (ii) The detection-probability statement is
  verified and sharpened by Monte Carlo: in the cosmic-variance limit the median
  expected p under H_geo is < 10^-3, with P(p < 0.01) = 92% at f_sky = 0.65 and
  97% full-sky (z = 1.7 and 2.0 respectively); the non-Gaussian null makes the
  z-score alone an understatement.
