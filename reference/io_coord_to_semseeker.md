# Convert a coordinate-based methylation data frame to SEMseeker internal format.

The input is a wide data frame where each row is a CpG position: CHR \|
START \| \[END\] \| sample1 \| sample2 \| ...

## Usage

``` r
io_coord_to_semseeker(df)
```

## Arguments

- df:

  Data frame with CHR and START columns (END optional).

## Value

Data frame with rownames = synthetic probe IDs, columns = samples.

## Details

The output is a data frame with synthetic probe IDs as rownames and
sample columns only (CHR/START/END are consumed and encoded in the
rowname).
