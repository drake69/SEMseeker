# Title

Title

## Usage

``` r
io_data_preparation(
  family_test,
  transformation_y,
  tempDataFrame,
  independent_variable,
  g_start,
  g_end,
  dototal,
  covariates,
  depth_analysis,
  key,
  transformation_x = "none"
)
```

## Arguments

- family_test:

  test or regression to apply

- transformation_y:

  transformation_y to apply to data

- tempDataFrame:

  data frame to use for test/regression

- independent_variable:

  regressor

- g_start:

  index of the first burden column in tempDataFrame

- g_end:

  index of the last burden column in tempDataFrame

- dototal:

  logical; if TRUE, append a column with the total (row-sum) burden

- covariates:

  vector of covariates to be found in the sample sheet

- depth_analysis:

  1 only sample, 2 chr, 3 alle genomic areas

- key:

  named list with AREA, SUBAREA, MARKER and FIGURE identifiers (used to
  name TOTAL columns)

## Value

A named list with two elements: `tempDataFrame` (the prepared and
optionally transformed data.frame) and `independent_variableLevels` (the
factor levels of the independent variable, or `NULL` for continuous
outcomes).
