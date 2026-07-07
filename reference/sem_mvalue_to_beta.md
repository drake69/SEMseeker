# Convert M-values to beta values

Applies the standard transformation \\\beta = 2^M / (1 + 2^M)\\
element-wise. Handles numeric matrices and numeric columns of data
frames. Non-numeric columns (e.g. annotation columns in
coordinate-format data frames) are preserved unchanged.

## Usage

``` r
sem_mvalue_to_beta(x, coord_cols = NULL)
```

## Arguments

- x:

  A numeric matrix, numeric vector, or data frame containing M-values.

- coord_cols:

  Optional character vector of column names to preserve unchanged (e.g.
  `c("CHR","START","END")`). If `NULL` (default), all numeric columns
  are converted; all non-numeric columns are preserved.

## Value

An object of the same class as `x` with M-values converted to beta
values.
