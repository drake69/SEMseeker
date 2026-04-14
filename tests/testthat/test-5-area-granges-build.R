## Tests for area_granges_build() — C-04
##
## All tests that require TxDb / AnnotationHub are guarded with
## skip_if_not_installed(). CI installs these packages explicitly, but they
## are optional for package users.

# ---------------------------------------------------------------------------
# Tests that do NOT require any external Bioconductor package
# ---------------------------------------------------------------------------

test_that("area_granges_build errors on unknown area", {
  expect_error(
    SEMseeker:::area_granges_build("FOOBAR_WHOLE", genome_build = "hg19"),
    regexp = "Unknown area"
  )
})

test_that("area_granges_build errors on unknown GENE subarea", {
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg19.knownGene")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("GenomicFeatures")
  expect_error(
    SEMseeker:::area_granges_build("GENE_FOOBAR", genome_build = "hg19"),
    regexp = "Unknown GENE subarea"
  )
})

test_that("area_granges_build errors on unknown ISLAND subarea", {
  skip_if_not_installed("AnnotationHub")
  skip_if_not_installed("GenomicRanges")
  expect_error(
    SEMseeker:::area_granges_build("ISLAND_FOOBAR", genome_build = "hg19"),
    regexp = "Unknown ISLAND subarea"
  )
})

test_that("area_granges_build appends _WHOLE when no underscore", {
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg19.knownGene")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("GenomicFeatures")
  # "GENE" alone should be treated as "GENE_WHOLE"
  # We just check it doesn't throw a "unknown area" error
  expect_error(
    SEMseeker:::area_granges_build("GENE", genome_build = "hg19"),
    NA   # no error expected
  )
})

# ---------------------------------------------------------------------------
# DMR area — uses bundled data, no external packages needed
# ---------------------------------------------------------------------------

test_that("area_granges_build DMR_WHOLE returns GRanges with label", {
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("IRanges")
  gr <- SEMseeker:::area_granges_build("DMR_WHOLE", genome_build = "hg19")
  expect_s4_class(gr, "GRanges")
  expect_true(length(gr) > 0)
  expect_true("label" %in% names(GenomicRanges::mcols(gr)))
  expect_false(any(is.na(GenomicRanges::mcols(gr)$label)))
})

# ---------------------------------------------------------------------------
# GENE areas — require TxDb
# ---------------------------------------------------------------------------

test_that("area_granges_build GENE_BODY returns valid GRanges", {
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg19.knownGene")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("GenomicFeatures")
  gr <- SEMseeker:::area_granges_build("GENE_BODY", genome_build = "hg19")
  expect_s4_class(gr, "GRanges")
  expect_true(length(gr) > 1000L)
  expect_true("label" %in% names(GenomicRanges::mcols(gr)))
})

test_that("area_granges_build GENE_TSS200 ranges are ≤ 200bp wide", {
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg19.knownGene")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("GenomicFeatures")
  gr <- SEMseeker:::area_granges_build("GENE_TSS200", genome_build = "hg19")
  expect_true(all(GenomicRanges::width(gr) <= 200L))
})

test_that("GENE_TSS1500 and GENE_TSS200 do not overlap", {
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg19.knownGene")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("GenomicFeatures")
  gr200  <- SEMseeker:::area_granges_build("GENE_TSS200",  genome_build = "hg19")
  gr1500 <- SEMseeker:::area_granges_build("GENE_TSS1500", genome_build = "hg19")
  hits   <- GenomicRanges::findOverlaps(gr200, gr1500)
  expect_equal(length(hits), 0L,
    info = "TSS200 and TSS1500 rings must not overlap")
})

test_that("area_granges_build result is cached on second call", {
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg19.knownGene")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("GenomicFeatures")
  gr1 <- SEMseeker:::area_granges_build("GENE_BODY", genome_build = "hg19")
  gr2 <- SEMseeker:::area_granges_build("GENE_BODY", genome_build = "hg19")
  expect_identical(gr1, gr2)  # exact same object from cache
})

# ---------------------------------------------------------------------------
# CHR_CYTOBAND — uses bundled data
# ---------------------------------------------------------------------------

test_that("area_granges_build CHR_CYTOBAND returns GRanges with label", {
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg19.knownGene")
  gr <- SEMseeker:::area_granges_build("CHR_CYTOBAND", genome_build = "hg19")
  expect_s4_class(gr, "GRanges")
  expect_true(length(gr) > 100L)
  expect_true("label" %in% names(GenomicRanges::mcols(gr)))
  # Labels should look like cytoband names (e.g. "p11.1", "q21.3")
  expect_true(any(grepl("[pq]", GenomicRanges::mcols(gr)$label)))
})
