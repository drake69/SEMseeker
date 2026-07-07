# Find gene sets unique to each group

Given a named list where each element is a character vector of gene set
identifiers, returns only those identifiers that appear exclusively in
one group and not in any other.

## Usage

``` r
enrich_find_unique_gene_sets(split_list)
```

## Arguments

- split_list:

  A named list of character vectors (one per group/category).

## Value

A named list of the same structure, containing only the gene sets that
are unique to each group. Groups with no unique sets are dropped.
