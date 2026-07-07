# Describe a data frame with summary statistics per column

Pure helper: no I/O, no side effects. Returns a data frame with one row
per column of `df` and the following fields: Variable, Class,
Missing_Values, Missing_Values_Percent, Unique_Values, Mean, Median,
Min, Max.

## Usage

``` r
util_describe_dataframe(df)
```

## Arguments

- df:

  A data frame to describe.

## Value

A data frame with summary statistics.
