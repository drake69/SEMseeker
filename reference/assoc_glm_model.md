# GLM association model

GLM association model

## Usage

``` r
assoc_glm_model(
  family_test,
  tempDataFrame,
  sig.formula,
  transformation_y,
  plot,
  samples_sql_condition = samples_sql_condition,
  key
)
```

## Arguments

- family_test:

  regression model family to apply (e.g. "gaussian", "binomial",
  "poisson")

- tempDataFrame:

  data frame to use for the model

- sig.formula:

  formula to apply the model

- transformation_y:

  transformation applied to the dependent variable (for labelling plots)

- plot:

  logical; if TRUE, generate and save a scatter/fit plot

- samples_sql_condition:

  SQL condition string used to filter samples (used for plot file
  naming)

- key:

  named list with AREA, SUBAREA, MARKER and FIGURE identifiers for this
  test

## Value

A list (as a data.frame row) with GLM fit results including AIC,
per-coefficient p-values and estimates, model performance metrics, and
optionally residual diagnostics.
