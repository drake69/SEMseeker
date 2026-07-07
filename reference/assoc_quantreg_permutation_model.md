# Quantile regression permutation model

Quantile regression permutation model

## Usage

``` r
assoc_quantreg_permutation_model(
  family_test,
  sig.formula,
  tempDataFrame,
  independent_variable,
  transformation_y,
  plot,
  samples_sql_condition = samples_sql_condition,
  key
)
```

## Arguments

- family_test:

  family test string encoding model type, quantile and permutation
  counts

- sig.formula:

  formula of the model

- tempDataFrame:

  data

- independent_variable:

  name of regressor

- transformation_y:

  transformation applied to the dependent variable (for labelling plots)

- plot:

  logical; if TRUE, generate diagnostic plots

- samples_sql_condition:

  SQL condition string used to filter samples (used for plot file
  naming)

- key:

  named list with AREA, SUBAREA, MARKER and FIGURE identifiers for this
  test

## Value

A numeric p-value from the permutation-based quantile regression test.
