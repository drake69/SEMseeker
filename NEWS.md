# semseeker NEWS

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
  New `family_test = "binomial_bulk"` dispatches to `glm_model_bulk()` for
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
  `intra_study_association_replication()` now calls this guard internally: stops
  with a clear message when origin results carry a different `GENOME_BUILD` than
  the current session.

### Breaking changes

- **`semseeker()` no longer auto-converts M-values to beta.**
  The `auto_convert_mvalues` parameter has been removed and the
  `.looks_like_mvalues()` helper deleted. Input values are now passed
  through unchanged; if you need conversion, call `mvalue_to_beta()`
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

- **`exists("signal_data")` scoped to local environment in `analyze_batch()` and
  `analyze_population()`** (E-01).
  The previous `exists("signal_data")` used `inherits = TRUE` (R default), which walked
  all the way up to `.GlobalEnv`. If the user had a `signal_data` object in their session
  (the normal case — they load data before calling `semseeker()`), the guard fired
  even though the local copy had already been freed, and the subsequent `rm(signal_data)`
  produced a spurious `"object not found"` warning. In the worst case (interactive session
  where local and global scopes coincide) it could silently delete the user's data.
  Fix: `exists("signal_data", envir = environment(), inherits = FALSE)` +
  `rm("signal_data", envir = environment())` in both files.

- **Polars inner join replaces positional zip in `mutations_get()`, `delta_single_sample()`,
  `deltar_single_sample()`; coverage banner in `analyze_population()`** (A-10).
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
  - `mutations_get()`, `delta_single_sample()`, `deltar_single_sample()` now use the helper;
    zero-overlap guard returns an empty result rather than crashing.
  - Coverage banner moved from per-sample (inside `mutations_get()`) to per-batch (at the
    start of `analyze_population()` before the per-sample loop). Emitted once per batch:
    `input_positions | beta_range_positions | covered_by_inner_join`.

### Bug fixes (A-09: bayes_analysis rewrite)

- **`bayes_analysis()`: 9 bugs fixed** (A-09).
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
  - **Output filename typo** (bug 7): `"bayes_analisys"` → `"bayes_analysis"`.
  - **`c` as loop variable** (bug 8): foreach variable named `c` shadowed the base
    `c()` function. Renamed to `col_idx`.
  - **Hardcoded thresholds** (bug 9): 0.9 / 0.1 are now `bayes_case_threshold` and
    `bayes_control_threshold` parameters with the original values as defaults.

### Statistical model changes

- **`lesions_get()`: replaced hypergeometric with binomial test** (A-01).
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
- Added differential signal analysis section (SIGNAL_MEAN / SIGNAL_RANGE) to
  getting-started vignette.

## semseeker 0.10.0 and earlier

See git log for earlier changes.
