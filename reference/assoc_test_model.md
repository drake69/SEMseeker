# Statistical test model dispatcher

Statistical test model dispatcher

## Usage

``` r
assoc_test_model(
  family_test,
  tempDataFrame,
  sig.formula,
  burdenValue,
  independent_variable,
  transformation_y,
  plot,
  samples_sql_condition = samples_sql_condition,
  key
)
```

## Arguments

- family_test:

  which family test to apply (e.g. "wilcoxon", "t.test", "kruskal.test",
  "pearson", "spearman")

- tempDataFrame:

  data frame to use with the test

- sig.formula:

  formula to apply

- burdenValue:

  name of the burden (dependent) column in tempDataFrame

- independent_variable:

  name of the independent variable column

- transformation_y:

  transformation applied to the dependent variable (for labelling plots)

- plot:

  logical; if TRUE, generate and save diagnostic box/scatter plots

- samples_sql_condition:

  SQL condition string used to filter samples (used for plot file
  naming)

- key:

  named list with AREA, SUBAREA, MARKER and FIGURE identifiers for this
  test

## Value

A list (as a data.frame row) with test results including the p-value,
test statistic, effect size, power, and model identifier; the exact
fields depend on the chosen `family_test`.
