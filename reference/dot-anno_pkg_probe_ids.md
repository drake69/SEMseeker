# Get probe IDs from an Illumina annotation package

Extracts the complete vector of probe identifiers (e.g. `cg00000029`)
from an annotation package.

## Usage

``` r
.anno_pkg_probe_ids(pkg)
```

## Arguments

- pkg:

  Character scalar: annotation package name.

## Value

Character vector of probe IDs (rownames of the Locations table).
