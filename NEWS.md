# semseeker NEWS

## semseeker 0.12.0 (dev)

### Statistical model changes

- **`lesions_get()`: replaced hypergeometric with binomial test** (A-01).
  The sliding window advances one probe at a time, so each probe participates
  in up to `sliding_window_size` consecutive windows — sampling with
  replacement. The hypergeometric distribution assumes sampling without
  replacement and produced inflated (too-small) p-values. The new test is:
  `P(X ≥ ENRICHMENT | Binomial(sliding_window_size, p0))` where
  `p0 = MUTATIONS_COUNT / PROBES_COUNT` is the empirical background
  mutation rate for the grouping unit (gene, chromosome, …). Expected
  impact: more conservative lesion calls, better calibration for samples
  with low or high global methylation variation.

## semseeker 0.11.0

### New features

- Added three pkgdown vignettes: Getting started, Association analysis (all 15+
  model families), Pathway and enrichment analysis (all 6 backends).
- Added GitHub Actions CI matrix: macOS, Ubuntu, Windows (`R-CMD-check.yml`).
- Added test coverage workflow (`test-coverage.yml`) with Codecov upload.
- Added `./ci-local.sh` for local Docker-based CI reproduction (`check` and
  `coverage` modes).

### Bug fixes

- Fixed `future::plan(multicore, workers = 0)` crash when `availableCores()`
  returns 1 (e.g. in covr subprocess): added `max(1L, nCore)` guard in
  `parallel_session.R`.

### Documentation

- Corrected SEM citations: replaced Teschendorff with correct attribution to
  Gentilini et al. 2015 (doi:10.18632/aging.100792) and Corsaro et al. 2023
  (doi:10.3390/cancers15164109).
- Added differential signal analysis section (SIGNAL_MEAN / SIGNAL_RANGE) to
  getting-started vignette.

## semseeker 0.10.0 and earlier

See git log for earlier changes.
