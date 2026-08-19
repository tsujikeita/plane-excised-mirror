# notebooks/ -- run order (Google Colab + Drive)

1. **PlaneMirror_A1_characterization_v0.2** -- staged and checkpointed:
   AXIS (best axes and landscapes, 10 maps x 6 configurations) ->
   FLOORBREAK (null ensembles to 1e5/1e4 at the masked anchors; ~6 h, resumable)
   -> BANDS (consensus-axis freeze + band decomposition vs 1000-realization null)
   -> EXCLUSION (disc-excision scan) -> FREQ (Nofi consistency).
   Note: on Colab the notebook now *requires* a successful Drive mount and stops
   otherwise (v0.2 fix; a silent local fallback previously endangered long runs).
2. **PlaneMirror_A1b_ABsplit_v0.2** -- split-systematics null test: scans PLA for
   split component-separated maps (the 2015 half-mission IDs are the ones that
   exist publicly), verifies half-sum reproduction, and bounds the
   half-difference contribution.

Both embed the verified modules verbatim (self-contained runs).
