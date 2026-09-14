# SHAPTURE 0.0.0.9009

- Renamed the multi-folder plotting recipe to
  `pipeline/run_shapture_visualization_loop_multiple_folder.R` and made folder
  discovery independent of comparison mode or folder naming.
- Both editable visualization recipes now define `top_feature_show`,
  `top_feature_label`, `top_feature_show_extended`, `width_per_feature`, and
  `height_per_feature` once near the top. Figure dimensions grow from the
  actual number of displayed features within plot-specific safety limits.
- Automatic sizing changes only graphics-device dimensions. Saved analysis
  tables, SHAP values, component scores, final scores, and rankings are
  unchanged.

# SHAPTURE 0.0.0.9008

- Cross-patient bootstraps can now use memory-aware PSOCK parallel execution
  on macOS, Windows, and Linux. The runtime limits workers using available CPU
  and RAM, validates a real worker before the full run, recycles workers in
  small batches, and replays the complete bootstrap serially after any setup,
  worker, or output-validation failure. Bootstrap plans, formulas, and results
  are unchanged.
- Each task now records the applied bootstrap backend, worker limits, memory
  estimates, timing, and any fallback in
  `docs/cross_patient_bootstrap_runtime.csv` and `run_summary.csv`.
- Feature-importance plots now color feature names by the saved model-expression
  pattern, using the same group colors and legend order as the beeswarm.
  Scoring, ranking, and plotted score values are unchanged.

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
