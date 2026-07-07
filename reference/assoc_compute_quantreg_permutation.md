# Title

Title

## Usage

``` r
assoc_compute_quantreg_permutation(sig.formula, df, tau, lqm_control)
```

## Arguments

- sig.formula:

  formula to apply

- df:

  dataframe to use

- tau:

  tau at which apply the wuantile regression

- lqm_control:

  specification of the lqmm package

## Value

A numeric scalar: the permuted quantile regression coefficient (used as
a single draw in the permutation distribution).
