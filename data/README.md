# data/

- **a1_axes.csv** -- best anti-symmetry axis per (map, nside, mask): axis pixel
  (RING), Galactic (l, b), S+ minimum, landscape width (fraction of directions
  within 5% of the minimum). 60 rows.
- **a1_consensus_axis.json** -- the frozen consensus axis (pixel 1134,
  l = 81.56, b = +14.48; rule recorded) used by all fixed-axis analyses and by
  the registered prediction.
- **a1_floorbreak.csv** -- masked-sky depths measured with enlarged ensembles
  (N16: 1e5, N32: 1e4 realizations): k, N, p, 68% intervals per map and mask.
- **a1_bands.csv** -- fixed-axis mirror-symmetric power per multipole band
  (2-4, 5-8, 9-16, 17-32, 33-64): data value, Gaussian null median, ratio,
  lower-tail p (1000 realizations). 10 maps.
- **a1_exclusion.npz** -- exclusion-scan maps Delta ln S+ (N_side=8 positions,
  15-deg discs, mirror-symmetric excision) per map, plus the Gaussian null
  |Delta| median scale (`null_absmed`).
- **a1b_absplit.csv** -- half-mission split test: half-sum reproduction
  (l=2-4 ratio, axis offset) and half-difference systematics fraction per method.
