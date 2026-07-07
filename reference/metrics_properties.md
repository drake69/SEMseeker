# metrics_properties

Metadata table describing the statistical properties of each SEM metric.
Used internally to determine ranking direction and scaling behaviour
during association analysis and pathway enrichment.

## Usage

``` r
metrics_properties
```

## Format

A data frame with one row per metric and columns including `Metric`
(metric name), `Higher_the_Better` (logical: whether higher values
indicate stronger signal), and `Affected_by_Scaling` (logical: whether
the metric is affected by data transformation).
