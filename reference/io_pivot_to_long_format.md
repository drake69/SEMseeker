# Get the pivot in long format instead of wide format

Get the pivot in long format instead of wide format

## Usage

``` r
io_pivot_to_long_format(
  marker,
  figure,
  area,
  subarea,
  phenotype_column,
  sample_sheet,
  areas_selection = NULL
)
```

## Arguments

- marker:

  marker to filer HYPER, HYPO, BOTH

- figure:

  DELTAS, DELTAQ,DELTAR, MUTATIONS

- area:

  GENE, DMR ...

- subarea:

  TSS1500 ...

- phenotype_column:

  column from the sample sheet to pair to each sample

- sample_sheet:

  sample sheet of samples

- areas_selection:

  genomic area to select, if NULL all areas will be selected

## Value

the pivot in a long format of 3 columnns, the phontype column with name
phenotype, the value of the marker and the area investigated
