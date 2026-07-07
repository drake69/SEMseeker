# Write a pivot sidecar metadata file

Writes a small JSON file alongside a pivot file (parquet or csv.gz) that
records the `genome_build` and `tech` of the session that created it.
The sidecar path is constructed by appending `_meta.json` to the pivot
base name (before the extension).

## Usage

``` r
core_pivot_sidecar_write(pivot_path)
```

## Arguments

- pivot_path:

  character. Full path of the pivot file (e.g.
  `.../Pivots/MUTATIONS/MUTATIONS_HYPER_GENE_TSS1500_hg19.parquet`).

## Value

Invisibly `NULL`.
