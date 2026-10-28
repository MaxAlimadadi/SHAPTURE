# SHAPTURE 0.0.0.9018

- Condition-level waterfall plots now label every displayed feature, independent
  of `top_feature_label`. This prevents unlabeled contributions from collapsing
  onto the same row when fewer labels than displayed features are requested.
  Both editable visualization recipes use `top_feature_show` for the waterfall
  display count. Force plots still use `top_feature_label` for selective labels.
- This is a visualization-only fix; analysis, SHAP values, scores, and rankings
  are unchanged. Previously saved results can be redrawn without rerunning models.

# SHAPTURE 0.0.0.9017

- Seurat and SingleCellExperiment adapters now accept source split columns
  whose names overlap SHAPTURE's prepared-data `condition` or `sample_id`
  fields. Conflicting split names receive a deterministic internal alias,
  returned as `dataset_res$split_by_col` with a source-to-prepared
  `column_mapping` table.
- The editable data-preparation templates now pass the returned prepared split
  name into analysis configuration. Nonconflicting inputs retain their prior
  names, values, column order, and analytical behavior. The script-config
  resolver also translates the source split name for compatibility with older
  copied templates.

# SHAPTURE 0.0.0.9016

- Dependency plots now add a 3% visual margin around the observed SHAP and
  feature-value ranges. This prevents boundary hexagons from being clipped
  while leaving grouped values, GAM fits, and analysis results unchanged.

# SHAPTURE 0.0.0.9015

- Added the editable `binned_beeswarm_hexagon_size` setting to both
  visualization recipes and a validated `hexagon_size` argument to the plot
  builder. The default remains `0.98`.
- Dependency-plot explanatory text now has a dedicated top header margin so it
  stays outside the SHAP-feature coordinate panel.

# SHAPTURE 0.0.0.9014

- Reduced the filled-hexagon size in the equal-count binned beeswarm so its
  Tukey layout remains readable without heavy symbol overlap. Equal-count
  groups, dependency plots, and all analytical results are unchanged.

# SHAPTURE 0.0.0.9013

- Binned beeswarm and SHAP-value dependency plots now draw filled hexagons as
  geometric polygons instead of Unicode font symbols. This fixes empty-square
  symbols on graphics devices without the required glyph and does not change
  equal-count groups, averages, smooths, scores, or rankings.

# SHAPTURE 0.0.0.9012

- Added a display-only equal-count binned beeswarm and one SHAP-value
  dependency PDF for every feature actually shown in the beeswarm. The
  original unit-level beeswarm is unchanged.
- The two new plots reuse deterministic feature-specific groups formed jointly
  from SHAP and model-seen feature values. They target 500 groups, require at
  least 10 units per group, and show original units below 500 valid units.
- Dependency smooths continue to use original unbinned values. Binning does
  not change correlations, component scores, final scores, or rankings.
- Beeswarm filenames and automatic figure sizing now use the number of saved
  features actually displayed when it is smaller than the requested count.

# SHAPTURE 0.0.0.9011

- Global visualization PDFs now size themselves from the matrices actually
  displayed: the cross-level landscape uses both selected-feature and
  split-level counts, while the direction-aware Jaccard matrix uses its
  split-level count on both axes.
- The editable global visualization recipe keeps visible `width` and `height`
  controls. Leave either value as `NULL` for automatic sizing or replace it
  with an exact size in inches. This changes graphics-device dimensions only;
  feature selection, scores, Jaccard values, and clustering are unchanged.

# SHAPTURE 0.0.0.9010

- Density plots now arrange their actual available features in a near-square
  facet grid and allocate consistent space to each panel. This prevents plots
  with few features from becoming unnecessarily wide and plots with many
  features from becoming disproportionately tall.
- Density sizing remains a post-analysis visualization change. Analysis
  tables, SHAP values, component scores, final scores, and rankings are
  unchanged.

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
