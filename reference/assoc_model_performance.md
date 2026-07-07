# Compute model performance metrics

Calculates a comprehensive set of regression performance metrics for
both the training set (fitted vs expected) and, optionally, a held-out
test set (prediction vs prediction_expected). Also detects possible
overfitting by comparing train and test metrics.

## Usage

``` r
assoc_model_performance(
  fitted_values,
  expected_values,
  prediction_values,
  prediction_expected_values
)
```

## Arguments

- fitted_values:

  Numeric vector of model-fitted (training) values.

- expected_values:

  Numeric vector of observed (training) values.

- prediction_values:

  Numeric vector of model predictions on the test set. Pass an empty
  vector ([`c()`](https://rdrr.io/r/base/c.html)) to skip test-set
  metrics.

- prediction_expected_values:

  Numeric vector of observed values for the test set. Ignored when
  `prediction_values` is empty.

## Value

A single-row `data.frame` with columns:

- mse, rmse, mape, mpe, sse, mae:

  Training-set error metrics.

- r_squared, r_squared_adj:

  Training-set goodness-of-fit.

- msle:

  Mean squared log error (training).

- mse_test, rmse_test, ...:

  Same metrics on test set (if provided).

- overfitting:

  Logical; `TRUE` if any test metric is worse than the corresponding
  training metric.

## Examples

``` r
fitted   <- c(1.1, 1.9, 3.2, 3.8)
expected <- c(1,   2,   3,   4  )
SEMseeker:::assoc_model_performance(fitted, expected, c(), c())
#>     mse      rmse       mape         mpe sse  mae r_squared r_squared_adj
#> 1 0.025 0.1581139 0.06466806 -0.01203648 0.1 0.15      0.98          0.97
#>          msle
#> 1 0.001894178
```
