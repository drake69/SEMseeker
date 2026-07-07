# Preliminary exploratory analysys of the data

Preliminary exploratory analysys of the data

## Usage

``` r
util_exploratory_analysis(
  categorical_variables,
  numerical_variables,
  sample_sheet,
  signal_data,
  max_missed_sample_sheet = 0.3,
  max_missed_signal_data = 0.1,
  sample_id_column = "Sample_ID",
  delete_keyword = "REMOVE",
  exploration_phase = "0",
  mapping_folder = NULL,
  sample_sheet_mapping = NULL,
  values_mapping = c(),
  removal_folder = NULL,
  removal_rules = c(),
  result_folder,
  ...
)
```

## Arguments

- categorical_variables:

  vector of variables or variables (with plus sign) to create
  exploratory pivot

- numerical_variables:

  vector of variables to have summary

- sample_sheet:

  path to the sample sheet

- signal_data:

  path to the signal data

- max_missed_sample_sheet:

  max number of missing values in each sample sheet column

- max_missed_signal_data:

  max number of missing values in each signal data column and row

- sample_id_column:

  name of the column with sample ids

- delete_keyword:

  character. Value in `values_mapping` that marks samples for removal
  (default `"REMOVE"`).

- exploration_phase:

  numeric value to preserve history of exploratory analysis

- mapping_folder:

  character. Optional path to folder containing mapping files; overrides
  `sample_sheet_mapping` / `values_mapping` if provided (default
  `NULL`).

- sample_sheet_mapping:

  rules to rename and remove columns

- values_mapping:

  rules to recode values and remove samples (source missed leave blank
  in the mapping file)

- removal_folder:

  character. Optional path to folder containing removal rule files
  (default `NULL`).

- removal_rules:

  character vector. Inline removal rules applied before mapping (default
  [`c()`](https://rdrr.io/r/base/c.html)).

- result_folder:

  path to the result folder

- ...:

  other parameters to define options for semseeker

## Value

Invisibly `NULL`. Cleaned sample sheet and signal data are written to
the result folder together with exploratory summary reports.

## Examples

``` r
result_dir <- tempdir()
if (FALSE) { # \dontrun{
util_exploratory_analysis(
  categorical_variables = c("Sample_Group", "Sex"),
  numerical_variables   = c("Age"),
  sample_sheet          = sample_sheet,
  signal_data           = beta_matrix
)
} # }
```
