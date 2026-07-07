# Final per-job save with TRANSFORMATION_X annotation

Extracted from association_analysis() (was inline at lines 436-444).
Adds the TRANSFORMATION_X column to the last marker's results and
performs the final CSV save; then emits the JOURNAL closing entry.

## Usage

``` r
util_finalize_job_results(
  results,
  inference_detail,
  family_test,
  filter_p_value,
  fileNameResults,
  start_time,
  processed_items
)
```

## Arguments

- results:

  data.frame from the last marker iteration (may be empty).

- inference_detail:

  single-row data.frame.

- family_test:

  character.

- filter_p_value:

  logical.

- fileNameResults:

  character. Path of the output CSV.

- start_time:

  POSIXct. Job start.

- processed_items:

  integer. Total items processed in this job.

## Value

Invisibly NULL.
