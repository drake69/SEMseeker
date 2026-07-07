# Read bedmethyl files (modkit/nanopolish) into a SEMseeker coordinate data frame

Parses one or more bedmethyl files (produced by `modkit pileup` or
`nanopolish call-methylation`) and returns a wide data frame in
SEMseeker coordinate format: `CHR, START, END, sample1, sample2, ...`.
Values are methylation fractions in \[0, 1\]. Positions with coverage
below `min_coverage` are set to `NA`.

## Usage

``` r
io_bedmethyl_read(file_paths, sample_ids = NULL, min_coverage = 5L)
```

## Arguments

- file_paths:

  Character vector of paths to bedmethyl files.

- sample_ids:

  Optional character vector of sample IDs, one per file. If `NULL`
  (default), IDs are derived from the file basenames (extension
  stripped).

- min_coverage:

  Integer. Positions with `N_valid_cov < min_coverage` are dropped.
  Default 5 (common threshold for Nanopore methylation calls).

## Value

A data frame with columns `CHR, START, END` followed by one numeric
column per sample (values in `[0, 1]`, `NA` if dropped by coverage
filter or absent in that sample).

## Details

Expected bedmethyl schema (modkit, tab-separated, no header):

1.  `chrom`

2.  `start_position`

3.  `end_position`

4.  `modified_base_code` (e.g. `m` for 5mC)

5.  `score`

6.  `strand`

7.  `start_code`

8.  `end_code`

9.  `color`

10. `N_valid_cov`

11. `fraction_modified` (percent, 0-100)

12. ... remaining modkit columns (ignored)
