# SHAPTURE 0.0.0.9007

- Replaced across-feature min-max scaling in the evidence score with a fixed
  median-anchor strategy for unbounded SHAP importance and absolute expression
  log2FC.
- Preserved detection-rate difference, within-condition expression
  consistency, Value-SHAP correlation, Cross-patient support, and model
  reproducibility on their natural 0-1 scales.
- Applied the same fixed log2FC anchor inside Cross-patient bootstraps and kept
  raw SHAP importance for its bootstrap retention calculation.
- Added `docs/score_anchor_summary.csv` so each task records its anchors and
  fixed candidate-feature population.
