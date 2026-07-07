# Log a per-inference-detail journal header

Extracted from association_analysis() (was inline at lines 124-134).
Emits the JOURNAL banner with the prettified inference_detail row
rendered as a kable table.

## Usage

``` r
core_log_inference_header(inference_detail)
```

## Arguments

- inference_detail:

  single-row data.frame.

## Value

Invisibly NULL.
