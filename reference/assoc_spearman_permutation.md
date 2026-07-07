# Spearman permutation model

Spearman permutation model

## Usage

``` r
assoc_spearman_permutation(
  family_test,
  sig.formula,
  tempDataFrame,
  independent_variable,
  plot,
  samples_sql_condition = samples_sql_condition,
  key
)
```

## Arguments

- family_test:

  family test string encoding model type and permutation parameters

- sig.formula:

  formula of the model

- tempDataFrame:

  data

- independent_variable:

  name of regressor

- plot:

  logical; if TRUE, generate diagnostic plots

- samples_sql_condition:

  SQL condition string used to filter samples (used for plot file
  naming)

- key:

  named list with AREA, SUBAREA, MARKER and FIGURE identifiers for this
  test

## Value

A numeric p-value from the permutation-based Spearman correlation test.
