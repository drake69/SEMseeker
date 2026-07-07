# Add or update a column in a data frame

If `col_name` already exists in `df`, its values are overwritten. If it
does not exist, a new column is appended. Handles empty data frames
(zero rows) by returning the new column as a single-column data frame.

## Usage

``` r
util_data_frame_add_column(df, col_name, value)
```

## Arguments

- df:

  A `data.frame`.

- col_name:

  Character scalar: name of the column to add or update.

- value:

  Vector of values to assign to the column. Must be compatible with the
  number of rows in `df` (or any length when `df` is empty).

## Value

The modified `data.frame` with the column added or updated.

## Examples

``` r
df <- data.frame(a = 1:3)
SEMseeker:::util_data_frame_add_column(df, "b", c(4, 5, 6))
#>   a b
#> 1 1 4
#> 2 2 5
#> 3 3 6
```
