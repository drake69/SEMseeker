# semseeker NEWS

## semseeker 0.99.5 (development)

### Breaking changes

- **The extent a number is valid over is a coordinate of its own, `SCOPE`
  (AI-255).** A burden over the whole sample and a burden per gene are the same
  marker reduced over different extents. They used to live in artefacts of
  different *shape* — a sibling CSV with samples down the rows, and a pivot with
  areas down the rows — and that difference in shape is what hid the difference
  in meaning. Both are now pivots, and the name says which is which:

  ```
  <MARKER>_<FIGURE>_<SCOPE>_<AREA>_<SUBAREA>_<AGGREGATION>_<GENOME_BUILD>

  MUTATIONS_HYPER_SAMPLE_PROBE_WHOLE_SUM_HG19     burden of the whole sample
  MUTATIONS_HYPER_SAMPLE_GENE_WHOLE_SUM_HG19      same, over gene probes only
  DELTAS_HYPO_INSTANCE_GENE_TSS200_MEDIAN_HG19    median per gene, TSS200 window
  SIGNAL_BETA_INSTANCE_PROBE_WHOLE_VALUE_HG19     the beta value per probe
  ```

  `SCOPE = SAMPLE` collapses the region class to one number per sample — such an
  artefact is one row tall — while `SCOPE = INSTANCE` keeps one row per gene,
  island or probe.

- **The aggregation requested at `SCOPE = INSTANCE` is now the one you get.** It
  was validated and then dropped on the way to the read, so asking for the
  median of a region returned its mean: the file existed, its name said nothing
  about which operator had produced it, and nothing complained. Result folders
  from earlier versions do not match the new names and recompute once.

- **`VALUE` names the identity.** A `PROBE` or `POSITION` row already holds a
  single position, so there is nothing to reduce; `VALUE` says so, where before
  the aggregation segment was simply absent. Its absence is again an error
  everywhere else.

- **`SAMPLE_STATS_RESULT.csv` is gone.** Its content is the `SCOPE = SAMPLE`
  artefacts, and the readable per-sample table is composed on read by
  `sem_study_summary_get(regions = ...)`, which joins them onto the sample
  sheet. `sem_study_summary_get()` with no arguments behaves as before.

- **`semseeker(sample_stats_scopes = ...)` was removed.** Which region classes
  you want is no longer a decision to be made before the run: the artefacts are
  built when they are asked for. The error that used to say *"produce it with
  `semseeker(sample_stats_scopes = ...)` and rerun the analysis"* is gone —
  changing your mind now costs one scan of the position pivot instead of a whole
  SEM run. Ask for a region class at analysis time with
  `association_analysis(inference_details$scopes)` or
  `sem_study_summary_get(regions = ...)`.

- **`MODE_LOW` / `MODE_HIGH` are spelled `MODELOW` / `MODEHIGH`**, and they are
  admissible only at `SCOPE = SAMPLE`. Per instance a gene holds about nineteen
  probes and a TSS200 window two to five: a two-peak density estimate on those
  is shaped by the bandwidth rather than by the data, and the old code returned
  that number instead of refusing. A minimum-numerosity guard now returns `NA`
  rather than a plausible-looking value.

- **`depth_analysis` no longer selects a code path.** The granularity of an
  artefact is the `(AREA, SUBAREA)` pair, and those pairs are only *partially*
  ordered — a TSS200 window and an open-sea stretch refine neither each other
  nor anything between them — so an integer scale projected a lattice onto a
  line. `1` still means `SCOPE = SAMPLE` and anything above it `INSTANCE`; the
  value is carried into the results as a label.

- **`N_PROBES` moved to `SAMPLE_SHEET_RESULT.csv`.** It is a property of the
  imputation — how many positions of that sample survived the treatment of
  missing values — not an aggregation of a marker, so it sits with the
  descriptive properties of the sample. Nothing is lost for the density: the
  `MEAN` of a binary marker *is* the density, denominator included.

- **Every computed quantity is now named by four coordinates (AI-248).** A
  number produced by SEMseeker is the reduction of a set of genomic positions,
  and until now the operator that reduced them was implicit — one per marker, so
  the marker and its direction identified the value. A scope can carry several
  reductions, so the operator became an axis of its own and the names changed
  accordingly:

  ```
  <SCOPE>_<MARKER>_<FIGURE>_<AGGREGATION>
  ```

  | before | after |
  |---|---|
  | `SAMPLE_MUTATIONS_HYPER` | `SAMPLE_MUTATIONS_HYPER_SUM` |
  | `SAMPLE_DELTAS_HYPO` | `SAMPLE_DELTAS_HYPO_MEAN` |
  | `SAMPLE_MEDIAN` | `SAMPLE_SIGNAL_BETA_MEDIAN` |
  | `SAMPLE_MODE_LOW` | `SAMPLE_SIGNAL_BETA_MODELOW` |
  | `SAMPLE_N_PROBES` | moved to `SAMPLE_SHEET_RESULT.csv` as `N_PROBES` — see the AI-255 entry above |

- **The figure of `SIGNAL` is the scale of the values, `BETA` or `MVALUE`.** It
  used to be `MEAN`, a placeholder that made the key unique while describing a
  way of aggregating — which now has an axis of its own. Pivot files carry it,
  so a beta run and an M-value run no longer write the same file name and
  overwrite each other in one folder. WGBS, ONT and PacBio stay `BETA`: a
  methylation fraction from reads is the same bounded scale as an array beta
  value, and the package compares the two on purpose.

- **Pivot file names carry the aggregation**, e.g.
  `MUTATIONS_HYPER_CHR_CYTOBAND_SUM_HG19.parquet`. Result folders produced by
  earlier versions do not match the new names and recompute once. This buys a
  guarantee that did not exist: the name used to say nothing about which
  operator had produced the file, so an existing pivot was reused on trust.

- **`inference_details$aggregation` is required.** A request that does not name
  the reduction does not identify what it wants tested. A malformed request (no
  aggregation, or a name outside the taxonomy) stops the run; a request that no
  marker of the run admits — the median of a 0/1 marker, say — drops that row
  with a warning naming it and lets the other rows proceed.

- **Results carry an `AGGREGATION` column**, part of the row identity along with
  marker, figure, area and subarea, in the deduplication, in the resume match
  and in the cross-study overlaps.

### New features

- **Descriptors on any scope.** The signal descriptors are no longer a separate
  path: `SIGNAL` is a marker like the others, so `sample_stats_scopes` now
  produces its median, mean, variance, IQR and — on the beta scale — its two
  modes restricted to a region class, e.g. `GENE_TSS1500_SIGNAL_BETA_MEDIAN`.

- **Density.** The mean of a binary marker over a scope is the fraction of its
  positions classified as epimutations. Unlike the raw burden it is comparable
  between region classes of different size, which the burden is not.

- Which reductions a marker admits is declared in one place. Continuous markers
  (`SIGNAL`, `DELTAS`, `DELTAR`) take the full descriptor set; counts take the
  sum and the mean, because their median and IQR are degenerate on a vector of
  zeros and ones. The two modes assume a bounded bimodal scale and are therefore
  restricted to the signal on the beta scale.


- **Per-sample burden restricted to a region class (AI-223 slice 2a).** The
  statistics sibling can now carry a scope other than the whole sample:
  `semseeker(sample_stats_scopes = c("SAMPLE", "GENE_TSS1500"))` adds
  `GENE_TSS1500_<MARKER>_<FIGURE>`, the burden computed over the probes of that
  region class only. Any registered `(AREA, SUBAREA)` pair of the run is a
  legal scope; `SAMPLE` is always produced. A scope that does not resolve stops
  the run instead of being ignored.

  Each position is counted **once**, even when the annotation maps it to
  several genes. The burden is computed from the POSITION pivot restricted to
  the probes of the region, not by summing the per-gene pivot — summing the
  latter would count a multi-gene probe once per gene.

- **Region scopes as dependent variable at `depth = 1`.** The new
  `inference_details$scopes` column (several separated by `"+"`, default
  `"SAMPLE"`) selects which scopes `association_analysis()` tests at depth 1.
  A region scope is still one value per sample, so its rows keep `DEPTH = 1`
  and `SUBAREA = "SAMPLE"`, with the scope in `AREA` (e.g. `GENE_TSS1500`) and
  `AREA_OF_TEST` unchanged. Requesting a scope that was never produced is an
  error naming the scope, not a silently empty result.

### Documentation

- The vignettes now show the names this version writes: `SIGNAL` pivots carry
  the scale of the run (`SIGNAL_BETA_*`, or `SIGNAL_MVALUE_*`), file names are
  uppercase, and the optional `AGGREGATION` segment is part of the documented
  pattern. `sem_coverage_analysis_report()` describes its `signal_data` path
  the same way.

## semseeker 0.99.4

### Breaking changes

- **The per-sample burden moved out of `SAMPLE_SHEET_RESULT.csv` into a new
  sibling file, `SAMPLE_STATS_RESULT.csv` (AI-223).** The columns that used to
  be appended to the sample sheet (`MUTATIONS_HYPER`, `DELTAS_HYPO`, …, and
  `PROBES_COUNT`) were aggregated over every probe with no genomic filter —
  that is the *sample* scope of the new artefact — and are now written there as
  `SAMPLE_<MARKER>_<FIGURE>` and `SAMPLE_N_PROBES`. The sample sheet stays the
  description of the study; the two files join on `Sample_ID`, and
  `sem_study_summary_get()` performs that join for you, so analyses inside the
  package are unaffected. Code that read the burden straight from
  `SAMPLE_SHEET_RESULT.csv` must read the sibling instead. *(Superseded in
  0.99.5: the sibling file no longer exists, the join is composed from the
  `SCOPE = SAMPLE` artefacts, and `N_PROBES` returned to the sample sheet.)*

### New features

- **Per-sample signal descriptors (AI-223).** Alongside the burden, the new
  sibling carries `SAMPLE_MEDIAN`, `SAMPLE_MEAN`, `SAMPLE_VARIANCE`,
  `SAMPLE_IQR` and `SAMPLE_N_PROBES` for every sample. On the beta scale it
  also carries the two modes of the bimodal distribution, `SAMPLE_MODE_LOW` and
  `SAMPLE_MODE_HIGH`, estimated as the highest density peak on each side of
  0.5. On M-values the distribution is unimodal and the split carries no
  meaning, so the two columns are omitted rather than filled with a
  meaningless number. Column names are composed by a single shared helper, so
  producer and consumer cannot drift apart. Region scopes
  (`<AREA>[_<SUBAREA>]_*`) are the next step; the naming already accommodates
  them.

- **Coverage is now a mandatory pre-step of every SEM analysis (AI-074).**
  Coverage used to be an opt-in report, so a run could go straight to SEM
  detection on input that barely overlaps the reference annotation (long-read
  positions against an Illumina manifest, a batch with a probe mismatch, the
  wrong genome build) and produce thresholds computed on a handful of probes.
  Every run now: writes the coverage charts, logs a coverage banner
  (`input_positions` / `covered_by_reference` / percentage), writes a
  `COVERAGE_GATE.json` sidecar for audit, and **stops** when coverage falls
  below the new `coverage_minimum` parameter (default 80%). The charts are
  produced even when the gate rejects the run — they are what shows which
  regions were missing. Coordinate-based technologies (WGBS, long reads)
  derive their annotation from the positions themselves, so the threshold is
  not applied there. Lower `coverage_minimum` explicitly for deliberate
  cross-technology runs.
- The coverage charts are now computed over the full annotation set instead of
  the areas selected for the run. `areas` defaults to `POSITION` only, which
  the coverage report skips, so on a default run the charts were previously
  never produced.

### Documentation

- `inst/CITATION` now carries the Zenodo DOI `10.5281/zenodo.5095416` instead
  of a Bioconductor DOI that does not resolve (the package is not on
  Bioconductor yet), and the README no longer says a Zenodo DOI "will be
  registered" — the archive has existed since 2021 (AI-119).

### Bug fixes

- **A sample sheet that shares no identifier with the computed pivots no
  longer produces a silently empty result (AI-083).**
  `sem_study_summary_total()` merged the per-sample burden onto the sample
  sheet by name with `all.x = TRUE`: when the two sides had no identifier in
  common, every burden column landed as `NA` while `PROBES_COUNT` stayed
  populated, and each `depth = 1` inference downstream then failed with
  "data are not the same size". The merge now stops with a message naming
  both sides, and reports partially missing samples as a warning. The
  accumulator inside the function is also a proper local binding now: it was
  resolved with `exists()`, which reaches the global environment.
- The same `exists()` pattern was removed from `sem_coverage_analysis()`,
  which could also raise "object 'tot_result' not found" when no area
  produced a usable table.

- **Sample identifiers containing `-`, `.` or spaces are now handled correctly
  (AI-224).** Sample-column names of the signal matrix are normalised
  (uppercase, non-alphanumeric characters replaced by `_`), but the sample
  sheet identifiers were used as-is when selecting those columns. With
  identifiers such as `C3L-00001-06` the two sides never matched:
  `util_exploratory_analysis()` wrote a cleaned signal parquet holding the
  `PROBE` column only, and the SEM pipeline stopped with
  `undefined columns selected`. Datasets whose identifiers are already
  alphanumeric (e.g. `GSM123456`) were unaffected. Both sides are now
  normalised through a single, idempotent internal helper before any
  name-based column selection; distinct identifiers that would collapse onto
  the same normalised name, and identifiers without a matching signal column,
  are reported explicitly instead of failing later.

## semseeker 0.99.3

### Documentation

- Aligned the delta-metric documentation in the vignettes and README with the
  implementation: the `DELTAS`/`DELTAR`/`DELTAP`/`DELTARP`/`DELTAQ`/`DELTARQ`
  descriptions now reflect the relative-delta ratio and its equal-width and
  quantile ranked variants. Corrected the pivot file-name examples to
  `Data/Pivots/<MARKER>/<MARKER>_<FIGURE>_<AREA>_<SUBAREA>_<build>.parquet` and
  updated a stale `enrichment_analysis()` reference.

### Bug fixes

- **`plot_box_plot()`: fixed R CMD check ERROR under ggplot2 >= 4.0.**
  `ggpubr::stat_compare_means(label = "p.format")` builds an internal
  `after_stat(create_p_label(...))` expression that ggplot2 >= 4.0 cannot
  resolve when ggpubr is loaded via `::` but not attached, causing
  "could not find function 'create_p_label'" (test-2-box-plot.R). The
  comparison p-value is now computed directly from `stats`
  (`t.test`/`wilcox.test`/`kruskal.test`/`anova`) and drawn with
  `ggplot2::annotate()`, removing the ggpubr rendering-path dependency.
  Behaviour change: the `kruskal.test` branch now shows the overall test
  p-value only; the per-pair `ggpubr` brackets have been dropped. Also fixed
  the `unit=`/`units=` partial-argument-match warning in `ggsave()`.

### Dependencies

- Bumped pinned GitHub Actions: `actions/cache` 4→6, `actions/upload-artifact`
  4→7, `actions/dependency-review-action` 4→5, `codecov/codecov-action` 4→7.

## semseeker 0.99.2

### Breaking changes

- **LESIONS detection now uses genomic distance, not probe count (AI-092).**
  The legacy `sliding_window_size` parameter (probe-count based, default 11)
  has been REMOVED. It is replaced by `LESIONS_BP` (default 2000), the
  maximum bp distance for two probes to be considered part of the same
  enrichment window. Same registration pattern as `DELTAQ_Q` / `DELTAP_B`:
  registered in `core_init_env()`, surfaces through the `Q` column of
  `keys_markers_default_discrete`, exposed via `semseeker(LESIONS_BP = ...)`.
  Rationale: probe-density varies dramatically across the genome and across
  array platforms (450K ~485k probes, EPIC ~865k, LONGREAD bedmethyl
  variable). A fixed 11-probe window covers ~70kb on average on 450K but
  ~500bp inside a CpG island — biologically incomparable. `LESIONS_BP=2000`
  matches the typical span of a CpG island or DMR and yields a single
  semantics across platforms. Migration: any caller passing
  `sliding_window_size=N` must replace it with `LESIONS_BP=M` (no equivalence
  formula; the metric has changed). Multi-window sensitivity will be tackled
  in AI-091 (vector-valued `LESIONS_BP`).

### New features

- **AI-190: CpG-island subareas aligned to Illumina `Relation_to_Island`.**
  The `ISLAND` area now exposes all six Illumina contexts plus the whole
  neighbourhood: `WHOLE`, `ISLAND`, `N_SHORE`, `S_SHORE`, `N_SHELF`,
  `S_SHELF`, `OPENSEA`. Previously the island core (`Island`) and the
  open-sea compartment (`OpenSea`) were lost. `ISLAND_WHOLE` is **redefined**
  to mean the whole island neighbourhood (core + shores + shelves, ±4 kb),
  mirroring `GENE_WHOLE`; use the new `ISLAND` subarea for the core alone.
  `OPENSEA` groups each open-sea CpG by the inter-neighbourhood genomic gap
  that contains it, labelled `chr:start-end` (never spanning a chromosome).
  Semantics are centralised in `island_opensea.R` and shared by both the
  Illumina (`anno_probe_annotation_build`) and coordinate/AnnotationHub
  (`anno_area_granges_build`) backends. See the getting-started vignette.

- **AI-044: binomial_bulk family + goodness-of-fit metrics extension.**
  New `family_test = "binomial_bulk"` dispatches to `assoc_glm_model_bulk()` for
  bulk per-probe logistic regression via `Rfast::glm_logistic` (parallelised
  with `foreach %dorng%`), ~10-20× faster than the per-probe `stats::glm`
  path. Drop-in replacement for `family_test = "binomial"` at PROBE-level
  inference details (LESIONS, MUTATIONS). Same legacy schema: per-coef
  PVALUE/ESTIMATE columns plus top-level PVALUE/PVALUE_ADJ.

  Standard errors are derived from the Fisher information at the MLE
  (`Var(β̂) = (X' diag(p(1-p)) X)^{-1}`), matching `stats::glm` Wald output.

  `io_data_preparation()` gains a universal degenerate-burden filter: columns
  with `var(Y) == 0` are dropped before reaching the model. Critical for
  LESIONS @ PROBE where ~92% of probes are all-zero across samples
  (manifest-aligned pivot, retained for positional join with annotations).
  The filter applies to every family — binomial GLM no longer diverges on
  constant Y, limma/voom no longer produces NaN t-stats, polynomial no
  longer fits rank-deficient designs.

  Ten new metrics registered in `metrics_properties.rda`:

  | Metric | Engine | Direction | Notes |
  |---|---|---|---|
  | `T_STAT_MODERATED` | limma_2, voom_2 | Higher better | eBayes-moderated t-stat per coef |
  | `B_STATISTIC` | limma_2, voom_2 | Higher better | log-odds of differential expression (lods) |
  | `F_STAT_MODERATED` | limma_2, voom_2 | Higher better | joint F across non-intercept coefs (parabolic test) |
  | `POSTERIOR_RESIDUAL_VAR` | limma_2, voom_2 | Lower better | s2.post diagnostic |
  | `MCFADDEN_R2` | binomial_bulk | Higher better | 1 − devi/null_devi pseudo-R² |
  | `NAGELKERKE_R2` | binomial_bulk | Higher better | scaled Cox-Snell pseudo-R² |
  | `C_STATISTIC_AUC` | binomial_bulk | Higher better | discrimination (= AUC) |
  | `DEVIANCE_RATIO` | binomial_bulk | Lower better | devi / null_devi |
  | `AIC_VALUE` | all GLM | Lower better | uppercase canonical of legacy `aic_value` |
  | `BIC_VALUE` | all GLM | Lower better | new |

  The metrics replace R²/R²_adj for limma/voom (lmFit doesn't return R²
  natively, and the eBayes-moderated diagnostics are more informative)
  and provide the analogous goodness-of-fit signal for binomial logistic
  regression (R² doesn't apply to {0,1} outcomes).

- **Session provenance metadata** (C-06).
  Every `semseeker()` run now writes `session_metadata.json` to the result
  folder at the start of analysis:
  `{"genome_build":"hg19","tech":"K850","semseeker_version":"0.99.0","created":"...","sample_n":120}`
  The file is human-readable and parseable without R, suitable for pipeline
  auditing and automated compatibility checks.
  Pivot files (parquet) additionally receive a sidecar `*_meta.json` with the
  same build/tech stamp; pivot file names now include the genome build as a
  suffix before the extension (e.g. `MUTATIONS_HYPER_GENE_TSS1500_hg19.parquet`)
  as belt-and-suspenders provenance.
  Inference CSV outputs gain two constant columns — `GENOME_BUILD` and `TECH`
  — that survive any downstream merge or stack operation.
  New exported function `core_check_session_compatibility(session_list)`:
  **stops** if `genome_build` differs across sessions (physically incomparable
  coordinates); **warns** if `tech` differs (cross-array meta-analysis is valid
  on the probe intersection but must be explicit).
  `assoc_intra_study_association_replication()` now calls this guard internally: stops
  with a clear message when origin results carry a different `GENOME_BUILD` than
  the current session.

### Breaking changes

- **`semseeker()` no longer auto-converts M-values to beta.**
  The `auto_convert_mvalues` parameter has been removed and the
  `.looks_like_mvalues()` helper deleted. Input values are now passed
  through unchanged; if you need conversion, call `sem_mvalue_to_beta()`
  explicitly before `semseeker()`. The beta-vs-M-value flag is still
  detected and stored in `ssEnv$beta` by `core_get_meth_tech()`, so downstream
  code that consults that flag continues to work.

- **`start_fresh` defaults to `FALSE`.**
  `core_init_env()` no longer deletes the result folder before starting. Existing
  results are preserved unless the caller explicitly passes `start_fresh = TRUE`.
  The Shiny UI exposes this as a checkbox ("Delete result folder before running",
  unchecked by default).

### Bug fixes

- **macOS: tests default to `multisession` instead of `multicore`** (E-14).
  `multicore` (fork) is unsafe on macOS with Polars' C++ thread pool — forked
  children can be killed by Mach exceptions. `setup.R` now selects `multisession`
  on Darwin, `multicore` on Linux. All `%dorng%` foreach bodies now call
  `core_update_session_info(ssEnv)` as their first statement to populate the worker's
  `.pkgglobalenv` — required because `multisession` workers are fresh R processes
  where the session singleton is empty.

- **`chr` prefix mismatch in Polars join** (E-13).
  `io_dump_sample_as_bed_file()` prepends `chr` to chromosome names, but
  `signal_thresholds` retains bare numbers from probe annotations. The inner
  join in `util_join_values_to_thresholds()` returned 0 rows → no mutations detected.
  Fix: `util_join_values_to_thresholds()` now strips the `chr` prefix from both sides
  before joining.

- **`exists("signal_data")` scoped to local environment in `sem_analyze_batch()` and
  `sem_analyze_population()`** (E-01).
  The previous `exists("signal_data")` used `inherits = TRUE` (R default), which walked
  all the way up to `.GlobalEnv`. If the user had a `signal_data` object in their session
  (the normal case — they load data before calling `semseeker()`), the guard fired
  even though the local copy had already been freed, and the subsequent `rm(signal_data)`
  produced a spurious `"object not found"` warning. In the worst case (interactive session
  where local and global scopes coincide) it could silently delete the user's data.
  Fix: `exists("signal_data", envir = environment(), inherits = FALSE)` +
  `rm("signal_data", envir = environment())` in both files.

- **Polars inner join replaces positional zip in `sem_mutations_get()`, `sem_delta_single_sample()`,
  `sem_deltar_single_sample()`; coverage banner in `sem_analyze_population()`** (A-10).
  The previous implementation sorted both `values` and `thresholds` by CHR/START/END and
  then compared them row-by-row (positional zip). When `signal_thresholds` came from a
  different run (cross-run analysis, e.g. Nanopore sample vs Illumina reference batch passed
  via `populationControlRangeBetaValues`), the two position sets could differ in size or
  overlap, producing silently wrong mutation/delta calls or an out-of-bounds crash.
  Changes:
  - Added `util_join_values_to_thresholds()` — a private helper that performs a Polars lazy inner
    join on `(CHR, START, END)`. Shared by all three per-sample functions to ensure consistent
    intersection logic. Required for Nanopore bedmethyl files (28M+ rows where base-R
    `merge()`/`match()` is too slow).
  - `sem_mutations_get()`, `sem_delta_single_sample()`, `sem_deltar_single_sample()` now use the helper;
    zero-overlap guard returns an empty result rather than crashing.
  - Coverage banner moved from per-sample (inside `sem_mutations_get()`) to per-batch (at the
    start of `sem_analyze_population()` before the per-sample loop). Emitted once per batch:
    `input_positions | beta_range_positions | covered_by_inner_join`.

### Bug fixes (A-09: assoc_bayes_analysis rewrite)

- **`assoc_bayes_analysis()`: 9 bugs fixed** (A-09).
  - **Loop off-by-one** (bug 1): `for (a in length(markers))` iterated only once
    (`a = 2`, only "LESIONS"), silently skipping "MUTATIONS". Fixed with `seq_along()`.
  - **Wrong file path** (bug 2): `read_delim(io_pivot_file_name)` passed the function
    object instead of the local variable `pivot_filename` → runtime connection error.
  - **Missing assignment** (bug 3): `tempDataFrame` used before being assigned from
    `pivot` → "object not found" crash. Added `tempDataFrame <- pivot`.
  - **`subset()` with string literal** (bug 4): `subset(df, "Sample_Group" != "Reference")`
    is always `TRUE` → Reference samples never filtered → contaminated Bayes estimates.
    Fixed to bare column reference `subset(df, Sample_Group != "Reference")`.
  - **Column drop off-by-one** (bug 3 cont.): `[, c(-1,-3)]` dropped the first data
    column or the `independent_variable` column depending on session config. Replaced
    with `colnames != "Sample_ID"` (remove only the merge key).
  - **`exists()` without scope** (bug 5): same pattern as E-01; fixed with
    `inherits = FALSE, envir = environment()`.
  - **Duplicate `max()` computation** (bug 6): `max_P_case` was computed twice,
    `max_P_control` never computed → filtered output missed Control criterion.
  - **Output filename typo** (bug 7): `"bayes_analisys"` → `"assoc_bayes_analysis"`.
  - **`c` as loop variable** (bug 8): foreach variable named `c` shadowed the base
    `c()` function. Renamed to `col_idx`.
  - **Hardcoded thresholds** (bug 9): 0.9 / 0.1 are now `bayes_case_threshold` and
    `bayes_control_threshold` parameters with the original values as defaults.

### Statistical model changes

- **`sem_lesions_get()`: replaced hypergeometric with binomial test** (A-01).
  The sliding window advances one probe at a time, so each probe participates
  in up to `sliding_window_size` consecutive windows — sampling with
  replacement. The hypergeometric distribution assumes sampling without
  replacement and produced inflated (too-small) p-values. The new test is:
  `P(X ≥ ENRICHMENT | Binomial(sliding_window_size, p0))` where
  `p0 = MUTATIONS_COUNT / PROBES_COUNT` is the empirical background
  mutation rate for the grouping unit (gene, chromosome, …). Expected
  impact: more conservative lesion calls, better calibration for samples
  with low or high global methylation variation.

## semseeker 0.11.0

### New features

- Added three pkgdown vignettes: Getting started, Association analysis (all 15+
  model families), Pathway and enrichment analysis (all 6 backends).
- Added GitHub Actions CI matrix: macOS, Ubuntu, Windows (`R-CMD-check.yml`).
- Added test coverage workflow (`test-coverage.yml`) with Codecov upload.
- Added `./ci-local.sh` for local Docker-based CI reproduction (`check` and
  `coverage` modes).

### Bug fixes

- Fixed `future::plan(multicore, workers = 0)` crash when `availableCores()`
  returns 1 (e.g. in covr subprocess): added `max(1L, nCore)` guard in
  `parallel_session.R`.

### Documentation

- Corrected SEM citations: replaced Teschendorff with correct attribution to
  Gentilini et al. 2015 (doi:10.18632/aging.100792) and Corsaro et al. 2023
  (doi:10.3390/cancers15164109).
- Added differential signal analysis section (SIGNAL_MEAN) to
  getting-started vignette.

## semseeker 0.10.0 and earlier

See git log for earlier changes.
