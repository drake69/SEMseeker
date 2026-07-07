# Cross-study meta-analysis of association results

Combines inference results from multiple studies using a random-effects
meta-analysis model (`metagen`). For each unique combination of FIGURE,
SUBAREA and AREA_OF_TEST, the function pools effect sizes (BETA) and
standard errors across studies and reports fixed-effect and
random-effect estimates with heterogeneity statistics.

## Usage

``` r
assoc_inter_study_association_meta_analysis(
  inference_details,
  statistic_parameter = "BETA",
  pvalue_column = "PVALUE_ADJ_ALL_BH",
  studies,
  studies_base_folder,
  result_folder
)
```

## Arguments

- inference_details:

  `data.frame` describing the inference configuration (same format as
  used by `association_analysis`).

- statistic_parameter:

  Character scalar: column name of the effect size estimate in the
  inference results (default `"BETA"`).

- pvalue_column:

  Character scalar: column name of the adjusted p-value (default
  `"PVALUE_ADJ_ALL_BH"`).

- studies:

  Character vector: study identifiers to include.

- studies_base_folder:

  Character scalar: base directory containing per-study result folders.

- result_folder:

  Character scalar: output directory for the meta-analysis results.

## Value

Invisibly returns a `data.frame` with one row per stratum containing
pooled effect estimates, confidence intervals, p-values, and
heterogeneity statistics (\\\tau^2\\, Q-test p-value).

## Details

Requires at least two studies per stratum; strata with fewer studies are
silently skipped.

## Examples

``` r
# Stub: see vignette('imprinting-disorders', package = 'SEMseeker') for a
# runnable Beckwith-Wiedemann workflow on the GSE133774 subset (AI-112b).
invisible(NULL)
```
