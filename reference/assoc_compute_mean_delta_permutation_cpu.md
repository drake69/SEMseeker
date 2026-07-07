# Compute mean delta permutation (CPU)

Compute mean delta permutation (CPU)

## Usage

``` r
assoc_compute_mean_delta_permutation_cpu(sig.formula, df, shuffle = FALSE)
```

## Arguments

- sig.formula:

  formula to apply

- df:

  dataframe to use

- shuffle:

  logical; if TRUE, permute the independent variable before fitting

## Value

A numeric scalar: the observed mean difference between the two groups
(used as the test statistic in the permutation test).
