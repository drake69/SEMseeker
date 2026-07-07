# Compute Spearman permutation statistic (CPU)

Compute Spearman permutation statistic (CPU)

## Usage

``` r
assoc_compute_spearman_permutation(sig.formula, df, shuffle = FALSE)
```

## Arguments

- sig.formula:

  formula to apply

- df:

  dataframe to use

- shuffle:

  logical; if TRUE, permute the independent variable before computing
  the correlation

## Value

A numeric scalar: the Spearman correlation coefficient between the
burden and the independent variable (used as the permutation test
statistic).
