# Load a named sub-table from an Illumina annotation package

Each data object (e.g. "Locations", "Islands.UCSC", "Other") is stored
as an independent lazy dataset in the package and must be loaded via
[`data()`](https://rdrr.io/r/utils/data.html). Accessing the top-level
S4 wrapper and then reading `@data$Locations` only yields a lazy
descriptor, not the actual table.

## Usage

``` r
.anno_pkg_load_table(pkg, table)
```

## Arguments

- pkg:

  Character scalar: annotation package name.

- table:

  Character scalar: name of the dataset to load (e.g. `"Locations"`).

## Value

A `data.frame` with probe IDs as rownames.
