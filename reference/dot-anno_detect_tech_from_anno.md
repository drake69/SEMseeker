# Detect Illumina array technology by probe-ID overlap

Queries each installed annotation package and counts how many probe IDs
from `probe_ids` are present in each array's probe list. Returns the
technology key (`"K27"`, `"K450"`, or `"K850"`) with the highest overlap
count.

## Usage

``` r
.anno_detect_tech_from_anno(probe_ids)
```

## Arguments

- probe_ids:

  Character vector of probe identifiers from the signal matrix.

## Value

Named integer vector of overlap counts, or `""` if no packages are
available.

## Details

Returns `""` if none of the annotation packages are installed.
