# sem_deltar_single_sample

sem_deltar_single_sample

## Usage

``` r
sem_deltar_single_sample(values, thresholds, sample_detail)
```

## Arguments

- values:

  data.frame of methylation values with columns CHR, START, END and
  signal value in column 4

- thresholds:

  data.frame of signal thresholds (from sem_signal_range_values) with
  columns CHR, START, END, signal_superior_thresholds,
  signal_inferior_thresholds, signal_median_values

- sample_detail:

  named list/row with at least Sample_ID and Sample_Group fields

## Value

invisibly NULL; HYPER and HYPO relative-delta results are written as
bedgraph.gz files. Only positions present in both `values` and
`thresholds` are processed (inner join on CHR, START, END via
`util_join_values_to_thresholds`). Returns `invisible(NULL)` early if
there is no positional overlap.
