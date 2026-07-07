# Quantile regression model (lqm)

Quantile regression model (lqm)

## Usage

``` r
assoc_quantreg_model(
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

  family test string encoding model type and quantile level (e.g.
  "quantreg_0.5")

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

A list (as a data.frame row) with quantile regression model results
including tau, p-value, standard error, regression coefficient,
confidence interval bounds, and model performance metrics.
