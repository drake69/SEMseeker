# Filter metrics by transformation type

Removes scale-sensitive metrics from the requested set when a non-scale
transformation (e.g. `"log"`, `"sqrt"`) is applied to the dependent
variable. Scale-sensitive metrics (e.g. MAE, RMSE) are meaningless after
a non-linear transformation because their units change.

## Usage

``` r
sem_metrics_filter(metrics, transformation_y)
```

## Arguments

- metrics:

  Character vector of metric names to filter (upper-case).

- transformation_y:

  Character scalar describing the transformation applied to the
  dependent variable. Use `"none"` to return all metrics unchanged,
  `"scale"` to keep all metrics (z-score does not change units), or any
  other value (e.g. `"log"`) to drop scale-affected metrics.

## Value

A sorted character vector of metric names that are valid for the given
transformation.

## Examples

``` r
SEMseeker:::sem_metrics_filter(c("MAE", "RMSE", "COUNT_SIGN"), "none")
#> [1] "MAE"        "RMSE"       "COUNT_SIGN"
SEMseeker:::sem_metrics_filter(c("MAE", "RMSE", "COUNT_SIGN"), "log")
#> [1] "COUNT_SIGN"
```
