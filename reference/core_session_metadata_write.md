# Write session provenance metadata

Creates (or overwrites) `session_metadata.json` in `result_folder`. The
JSON is human-readable and does not require R to parse, making it
suitable for pipeline auditing and automated compatibility checks across
studies.

## Usage

``` r
core_session_metadata_write(result_folder, sample_n = 0L)
```

## Arguments

- result_folder:

  character. Path to the SEMseeker result folder.

- sample_n:

  integer. Total number of samples across all batches (default `0L`).

## Value

Invisibly returns the metadata list that was serialised to JSON.

## Details

Fields written:

- genome_build:

  Reference assembly used (`"hg19"` / `"hg38"` / `"mm10"`).

- tech:

  Methylation technology detected or declared by the user (`"K850"`,
  `"K450"`, `"K27"`, `"WGBS"`, `"LONGREAD"`, or `""` if not yet
  determined at write time).

- semseeker_version:

  Package version string.

- created:

  ISO-8601 timestamp of the run.

- sample_n:

  Total number of samples in the run.
