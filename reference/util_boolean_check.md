# Coerce a value to logical

Converts any representation of TRUE/FALSE (logical, character, or
numeric) to a proper R logical value. Returns `FALSE` for `NA`, empty
string, or `NULL`.

## Usage

``` r
util_boolean_check(x)
```

## Arguments

- x:

  A scalar: logical, character (`"TRUE"`, `"FALSE"`, `"true"`, `"T"`,
  `"1"`) or numeric (`1`, `0`).

## Value

A single logical value (`TRUE` or `FALSE`).

## Examples

``` r
SEMseeker:::util_boolean_check("TRUE")   # TRUE
#> [1] TRUE
SEMseeker:::util_boolean_check("false")  # FALSE
#> [1] FALSE
SEMseeker:::util_boolean_check(1)        # TRUE
#> [1] TRUE
SEMseeker:::util_boolean_check(NA)       # FALSE
#> [1] FALSE
```
