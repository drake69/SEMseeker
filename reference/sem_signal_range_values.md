# calculate the range of signal values to define the outlier

calculate the range of signal values to define the outlier

## Usage

``` r
sem_signal_range_values(populationMatrix, batch_id, probe_features)
```

## Arguments

- populationMatrix:

  matrix of methylation for the population under calculation (probes ×
  samples)

- batch_id:

  character string identifying the batch; used to name the cached
  parquet output file

- probe_features:

  data.frame of probe annotations (from anno_probe_features_get) with
  columns PROBE, CHR, START, END

## Value

data.frame of per-probe thresholds with columns
signal_inferior_thresholds, signal_superior_thresholds,
signal_median_values, iqr, q1, q3, PROBE, CHR, START, END
