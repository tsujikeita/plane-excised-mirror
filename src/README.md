# src/

- **plane_mirror.py** -- theme-A tools, all verified against the audit reference
  implementation (max relative deviation 1.2e-6):
  `MirrorBatch` (direction scan over map batches; N16 33 ms/map, N32 0.9 s/map,
  enabling 1e5-realization ensembles), `FixedAxisMirror` (O(npix) fixed-axis
  symmetric power; band decomposition, band-sum completeness 5.5e-9; exclusion
  scan), `scan_S` (full S(n) landscapes), `with_mask` (mask swap sharing the
  reflection table), axis helpers (+/- identification).
- **phase2_core.py** -- shared audit core (masks, transfer conventions, CRN
  master realizations, MirrorStat reference implementation); identical to the
  copy in the companion repositories.
