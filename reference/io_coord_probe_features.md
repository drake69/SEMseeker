# Build a minimal probe_features data frame from synthetic probe IDs.

Used inside `sem_analyze_batch()` for WGBS/LONGREAD data in place of the
Bioconductor-annotation-based `anno_probe_features_get("PROBE")` call.

## Usage

``` r
io_coord_probe_features(probe_ids)
```

## Arguments

- probe_ids:

  Character vector of synthetic probe IDs.

## Value

data.frame with columns PROBE, CHR, START, END.
