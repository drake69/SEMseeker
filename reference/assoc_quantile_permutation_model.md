# assoc_quantile_permutation_model calculate differences between the same quantile of two distribution

assoc_quantile_permutation_model calculate differences between the same
quantile of two distribution

## Usage

``` r
assoc_quantile_permutation_model(
  family_test,
  sig.formula,
  tempDataFrame,
  independent_variable,
  transformation_y,
  plot = plot,
  samples_sql_condition = samples_sql_condition,
  key
)
```

## Arguments

- family_test:

  family quantile

- sig.formula:

  formula of the model

- tempDataFrame:

  data

- independent_variable:

  name of regressor

- transformation_y:

  transformation to apply to the dependent variable

- plot:

  logical; if TRUE, generate diagnostic plots

- samples_sql_condition:

  SQL condition string used to filter samples (used for plot file
  naming)

- key:

  named list with AREA, SUBAREA, MARKER and FIGURE identifiers for this
  test

## Value

A numeric p-value from the permutation-based quantile-difference test.
