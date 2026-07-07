# Run depth=1 (sample-level) association for one marker

Extracted from association_analysis() (was inline at lines 218-290).
Reads per-sample counts already present in study_summary columns,
optionally resumes from a partial results CSV, then applies the stat
model per key (one row per genomic region key).

## Usage

``` r
sem_run_depth1_marker(
  prep,
  keys,
  family_test,
  fileNameResults,
  filter_p_value,
  ssEnv,
  ...
)
```

## Arguments

- prep:

  list returned by sem_prepare_study_for_analysis().

- keys:

  data.frame of keys for this marker (subset of
  ssEnv\$keys_markers_figures).

- family_test:

  character.

- fileNameResults:

  character. Path of the output CSV.

- filter_p_value:

  logical.

- ssEnv:

  list. Session environment from core_get_session_info().

- ...:

  forwarded to assoc_apply_stat_model().

## Value

list(results = data.frame, processed_items = integer). Side effect:
writes the CSV via assoc_analysis_save_results().
