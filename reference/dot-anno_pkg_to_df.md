# Build a data.frame from an Illumina annotation package

Loads the Locations, Islands.UCSC, and Other sub-tables from the
annotation package and combines them into a single data.frame with one
row per probe. Does not require minfi.

## Usage

``` r
.anno_pkg_to_df(pkg)
```

## Arguments

- pkg:

  Character scalar: annotation package name.

## Value

A `data.frame` with probe IDs as rownames and columns: `chr`, `pos`,
`strand`, `Islands_Name`, `Relation_to_Island`, `UCSC_RefGene_Name`,
`UCSC_RefGene_Group`, and further columns from the Other table.
