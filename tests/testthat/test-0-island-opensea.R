# Unit tests for the shared CpG-island / OpenSea semantics (island_opensea.R).
# Pure GRanges logic — no Bioconductor annotation packages or downloads.

testthat::skip_if_not_installed("GenomicRanges")
testthat::skip_if_not_installed("IRanges")

test_that(".islands_gr_from_names parses unique coordinate strings, ignores junk", {
  gr <- SEMseeker:::.islands_gr_from_names(
    c("chr1:1000-2000", "chr1:1000-2000", "chr2:5000-6000", NA, "", "OpenSea"))
  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 2L)  # dedup + drop NA/""/non-coordinate
  expect_setequal(as.character(GenomicRanges::seqnames(gr)), c("chr1", "chr2"))
})

test_that(".opensea_gaps returns inter-neighbourhood gaps labelled by coordinate", {
  islands <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(start = c(10000L, 30000L), end = c(11000L, 31000L)))

  gaps_gr <- SEMseeker:::.opensea_gaps(islands, chrom_ends = list(chr1 = 50000L))

  # neighbourhoods: [6000,15000] and [26000,35000] -> gaps:
  #   [1,5999], [15001,25999], [35001,50000]
  expect_true(all(as.character(GenomicRanges::strand(gaps_gr)) == "*"))
  labs <- GenomicRanges::mcols(gaps_gr)$label
  expect_true("chr1:15001-25999" %in% labs)
  expect_true(all(grepl("^chr1:[0-9]+-[0-9]+$", labs)))
  # No gap crosses into a neighbourhood
  expect_false(any(GenomicRanges::start(gaps_gr) >= 6000 &
                   GenomicRanges::end(gaps_gr)   <= 15000))
})

test_that(".assign_opensea_labels maps each open-sea CpG to its gap; NA inside islands", {
  islands <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(start = c(10000L, 30000L), end = c(11000L, 31000L)))

  labels <- SEMseeker:::.assign_opensea_labels(
    probe_chr = rep("chr1", 4L),
    probe_pos = c(1000L, 20000L, 40000L, 10500L),
    island_gr = islands  # no seqlengths: helper extends universe to last CpG
  )

  expect_true(grepl("^chr1:[0-9]+-[0-9]+$", labels[1]))  # before first island
  expect_equal(labels[2], "chr1:15001-25999")            # between islands
  expect_true(grepl("^chr1:35001-", labels[3]))          # after last island
  expect_true(is.na(labels[4]))                          # inside island core
})
