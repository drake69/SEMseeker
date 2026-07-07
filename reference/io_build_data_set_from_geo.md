# io_build_data_set_from_geo

io_build_data_set_from_geo

## Usage

``` r
io_build_data_set_from_geo(GEOgse, downloadFiles = 0, result_folder, ...)
```

## Arguments

- GEOgse:

  geo accession dataset identification

- downloadFiles:

  0 means download all files from Gene Expression Omnibus (GEO),
  different than zero means how many to download

- result_folder:

  where sample sheet and files will be saved

- ...:

  additional arguments passed to internal helpers

## Value

samplesheet, and sample's file saved and samplesheet csv

## Examples

``` r
result_dir <- tempdir()
if (FALSE) { # \dontrun{
io_build_data_set_from_geo(
  GEOgse        = "GSE55763",
  downloadFiles = 1,
  result_folder = tempdir()
)
} # }
```
