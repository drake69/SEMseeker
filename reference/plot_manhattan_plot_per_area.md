# Manhattan plot of association results per genomic area

Manhattan plot of association results per genomic area

## Usage

``` r
plot_manhattan_plot_per_area(
  marker,
  figure,
  area,
  subarea,
  family,
  adjust_method,
  phenotype,
  only_significant_areas = FALSE
)
```

## Arguments

- marker:

  character. SEM metric to plot (e.g. `"MUTATIONS"`, `"DELTARP"`,
  `"DELTAR"`, `"DELTAQ"`).

- figure:

  character. Mutation direction: `"HYPO"` or `"HYPER"`.

- area:

  character. Genomic area level (e.g. `"GENE"`, `"ISLAND"`, `"DMR"`).

- subarea:

  character. Sub-area within the area (e.g. `"TSS1550"`, `"WHOLE"`).

- family:

  character. Statistical family used in the association analysis (e.g.
  `"wilcoxon"`, `"gaussian"`).

- adjust_method:

  character. Column name of the adjusted p-value to use for colouring
  (e.g. `"BH"`).

- phenotype:

  character. Sample sheet column used to colour points (e.g.
  `"Sample_Group"`).

- only_significant_areas:

  logical. If `TRUE`, show only regions with adjusted p-value \< 0.05
  (default `FALSE`).

## Value

Invisibly `NULL`. A Manhattan plot PNG is saved under
`Charts/MARKER_PER_AREA/` in the active result folder.

## Examples

``` r
result_dir <- tempdir()
if (FALSE) { # \dontrun{
plot_manhattan_plot_per_area(
  marker        = "DELTARP",
  figure        = "HYPO",
  area          = "GENE",
  subarea       = "WHOLE",
  family        = "wilcoxon",
  adjust_method = "BH",
  phenotype     = "Sample_Group"
)
} # }
```
