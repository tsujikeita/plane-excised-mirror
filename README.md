# The plane-excised mirror anti-symmetry of the CMB

Companion repository to K. Tsuji, *"The plane-excised mirror anti-symmetry of the
CMB: characterization and a registered prediction for future polarization data"*
(JCAP, submitted; arXiv:XXXX.XXXXX). Third paper of a series; companions:
`lowl-anomaly-audit` (the preprocessing audit) and `tail-lottery` (sub-percent
p-value stability).

**Registered prediction.** The polarization consequence of a geometric origin is
frozen -- axis, statistic, mask rule, null, decision rule -- in advance of suitable
data, at public commit **63aa853**:
`docs/prediction_registration_v1.0.md`
(permanent link: github.com/tsujikeita/plane-excised-mirror/blob/63aa853/docs/prediction_registration_v1.0.md).
A later commit appends Correction 1 (two benchmark descriptions sharpened; no
frozen quantity changed).

## Headline numbers
- Depth, measured below the Monte-Carlo floor: p = (1.4-5)e-4 under the common
  Galactic mask (1e5/1e4-realization ensembles), 10 maps x 2 resolutions
- One axis: all 60 (map, configuration) best axes within 3.7 deg; equals the
  published full-sky axis (l,b) = (264, -17) [antipode (81.6, +14.5)]
- Carrier l = 2-4: symmetric power suppressed to 0.09-0.11, frequency-independent
  (70-143 GHz, +/-10%)
- Symmetric residual localized near the Galactic anticenter (175, +18), 5x null
- Half-mission split systematics: <= 0.02% of the missing symmetric power
- Polarization: H_corr z = 0.5 (CV limit) -> Planck cannot decide (pre-registered
  NO-GO honored); H_geo median expected p < 1e-3, P(p<0.01) = 92-97% (CV limit)
  -> LiteBIRD-class data can

## Repository map
| Path | Contents |
|---|---|
| `docs/` | The registered prediction (frozen v1.0 + Correction 1) |
| `notebooks/` | A1: characterization (axes, floor-break, bands, exclusion, frequency) / A1b: half-mission split null test |
| `src/` | `plane_mirror.py` (MirrorBatch, FixedAxisMirror, exclusion scan) + shared `phase2_core.py` |
| `scripts/` | `make_figures.py` (figs 1-4 from `data/`) / `forecast_polarization.py` (Sec. 4 benchmarks + fig 5) |
| `data/` | All characterization tables and the exclusion maps (see `data/README.md`) |
| `figures/` | Paper figures (PNG+PDF), regenerable via `scripts/` |

## Quick start
pip install numpy==2.0.2 healpy==1.20.0 scipy pandas matplotlib camb
git clone --depth 1 https://github.com/LauraHerold/CMBanom.git   # masks + maps
python scripts/make_figures.py
python scripts/forecast_polarization.py

## License / citation / acknowledgment
See LICENSE (to be added). Please cite the paper. Developed in extensive
collaboration with Claude (Anthropic).
