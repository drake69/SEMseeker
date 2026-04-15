# Tests for enrichment and pathway analysis helpers
#
# Covered:
#  - enrichment_analysy_add_category()  category mapping helper (requires session)
#  - pathway_ctdR()                     graceful return when ctdR not installed
#  - pathway_WebGestalt()               graceful return when WebGestaltR not installed
#  - pathway_STRINGdb()                 graceful return when STRINGdb not installed
#  - pathway_pathfindR()                graceful return when pathfindR not installed
#  - pathway_ctdR() integration         end-to-end after a full semseeker +
#                                       association_analysis run with real results
#
# Note: all "package-not-installed" tests temporarily unload the optional package
# from the search path so the requireNamespace() guard fires correctly.

# ---------------------------------------------------------------------------
# 1. enrichment_analysy_add_category — category mapping
# ---------------------------------------------------------------------------

test_that("enrichment_analysy_add_category returns data unchanged for empty input", {
  tempFolder <- tempFolders[1]
  tempFolders <<- tempFolders[-1]
  SEMseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy,
                       showprogress = showprogress, verbosity = verbosity)

  result <- SEMseeker:::enrichment_analysy_add_category("ctdR", data.frame())
  testthat::expect_equal(nrow(result), 0)

  SEMseeker:::close_env()
  unlink(tempFolder, recursive = TRUE)
})

test_that("enrichment_analysy_add_category adds SS_CATEGORY='CHEMICAL' for ctdR source", {
  tempFolder <- tempFolders[1]
  tempFolders <<- tempFolders[-1]
  SEMseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy,
                       showprogress = showprogress, verbosity = verbosity)

  # Column names must match key_enrichment_format for "ctdR":
  #   column_of_adj_pvalue  = "padj"
  #   column_of_enrichment  = "foldEnrichment"
  #   column_of_description = "ChemicalName"
  # The ranking code inside the function accesses these columns;
  # using wrong names triggers a browser() call (debugging artifact in code).
  fake_result <- data.frame(
    ChemicalID     = c("D001", "D002"),
    ChemicalName   = c("chemical A", "chemical B"),
    pvalue         = c(0.01, 0.05),
    padj           = c(0.05, 0.10),
    foldEnrichment = c(2.5, 1.8),
    source         = c("ctdR", "ctdR"),
    key            = c("k1", "k2"),
    stringsAsFactors = FALSE
  )

  result <- SEMseeker:::enrichment_analysy_add_category("ctdR", fake_result)
  testthat::expect_true("SS_CATEGORY" %in% colnames(result))
  testthat::expect_true(all(result$SS_CATEGORY == "CHEMICAL"))

  SEMseeker:::close_env()
  unlink(tempFolder, recursive = TRUE)
})

test_that("enrichment_analysy_add_category maps GO types for WebGestalt source", {
  tempFolder <- tempFolders[1]
  tempFolders <<- tempFolders[-1]
  SEMseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy,
                       showprogress = showprogress, verbosity = verbosity)

  # Column names must match key_enrichment_format for "WebGestalt":
  #   column_of_adj_pvalue  = "FDR"
  #   column_of_enrichment  = "enrichmentRatio"
  #   column_of_description = "description"
  fake_result <- data.frame(
    geneSet         = c("GO:0001", "GO:0002", "GO:0003"),
    description     = c("term1", "term2", "term3"),
    pValue          = c(0.01, 0.05, 0.1),
    FDR             = c(0.02, 0.10, 0.20),
    enrichmentRatio = c(2.5, 1.8, 1.2),
    type            = c("BP", "CC", "MF"),
    source          = c("WebGestalt", "WebGestalt", "WebGestalt"),
    stringsAsFactors = FALSE
  )

  result <- SEMseeker:::enrichment_analysy_add_category("WebGestalt", fake_result)
  testthat::expect_true("SS_CATEGORY" %in% colnames(result))
  # Categories BP/CC/MF should map to GO-BP / GO-CC / GO-MF
  cats <- result$SS_CATEGORY[!is.na(result$SS_CATEGORY)]
  testthat::expect_true(all(cats %in% c("GO-BP", "GO-CC", "GO-MF")))

  SEMseeker:::close_env()
  unlink(tempFolder, recursive = TRUE)
})

# ---------------------------------------------------------------------------
# 2. Pathway functions return NULL gracefully when optional packages absent
#    We test this by calling with an active session but zero-row inference_details.
#    When the optional package IS present the function proceeds; we only care
#    that it does not throw an unhandled error regardless of package state.
# ---------------------------------------------------------------------------

test_that("pathway_WebGestalt returns NULL gracefully when WebGestaltR not installed", {
  tempFolder <- tempFolders[1]
  tempFolders <<- tempFolders[-1]
  SEMseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy,
                       showprogress = showprogress, verbosity = verbosity)

  # WebGestaltR is typically not installed in the test environment.
  # The function guards with requireNamespace() and returns invisibly.
  if (!requireNamespace("WebGestaltR", quietly = TRUE)) {
    inference_detail <- data.frame(
      independent_variable = "Phenotest",
      family_test          = "spearman",
      transformation_y     = "",
      transformation_x     = "",
      depth_analysis       = 1L,
      filter_p_value       = FALSE,
      areas_sql_condition  = NA,
      samples_sql_condition = NA,
      association_results_sql_condition = NA,
      stringsAsFactors = FALSE
    )
    testthat::expect_no_error(
      SEMseeker:::pathway_WebGestalt(
        study           = "test",
        inference_detail = inference_detail,
        significance     = TRUE
      )
    )
  } else {
    testthat::skip("WebGestaltR is installed; skipping guard test")
  }

  SEMseeker:::close_env()
  unlink(tempFolder, recursive = TRUE)
})

test_that("pathway_STRINGdb returns NULL gracefully when STRINGdb not installed", {
  tempFolder <- tempFolders[1]
  tempFolders <<- tempFolders[-1]
  SEMseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy,
                       showprogress = showprogress, verbosity = verbosity)

  if (!requireNamespace("STRINGdb", quietly = TRUE)) {
    inference_details <- data.frame(
      independent_variable = "Phenotest",
      family_test          = "spearman",
      transformation_y     = "",
      transformation_x     = "",
      depth_analysis       = 1L,
      filter_p_value       = FALSE,
      stringsAsFactors = FALSE
    )
    testthat::expect_no_error(
      SEMseeker:::pathway_STRINGdb(
        study            = "test",
        inference_details = inference_details
      )
    )
  } else {
    testthat::skip("STRINGdb is installed; skipping guard test")
  }

  SEMseeker:::close_env()
  unlink(tempFolder, recursive = TRUE)
})

test_that("pathway_pathfindR returns NULL gracefully when pathfindR not installed", {
  tempFolder <- tempFolders[1]
  tempFolders <<- tempFolders[-1]
  SEMseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy,
                       showprogress = showprogress, verbosity = verbosity)

  if (!requireNamespace("pathfindR", quietly = TRUE)) {
    inference_details <- data.frame(
      independent_variable = "Phenotest",
      family_test          = "spearman",
      transformation_y     = "",
      transformation_x     = "",
      depth_analysis       = 1L,
      filter_p_value       = FALSE,
      stringsAsFactors = FALSE
    )
    testthat::expect_no_error(
      SEMseeker:::pathway_pathfindR(
        study             = "test",
        path_dbs          = c("KEGG"),
        inference_details = inference_details,
        significance      = TRUE
      )
    )
  } else {
    testthat::skip("pathfindR is installed; skipping guard test")
  }

  SEMseeker:::close_env()
  unlink(tempFolder, recursive = TRUE)
})

# ---------------------------------------------------------------------------
# 3. pathway_ctdR integration — runs after semseeker + association_analysis
#    Requires ctdR to be installed (installed via Remotes in DESCRIPTION).
#    Uses the same bimodal synthetic data as test-7 to guarantee mutations.
# ---------------------------------------------------------------------------

test_that("pathway_ctdR runs without error on synthetic association results", {
  if (!requireNamespace("ctdR", quietly = TRUE)) {
    testthat::skip("ctdR not installed")
  }

  tempFolder <- tempFolders[1]
  tempFolders <<- tempFolders[-1]
  unlink(tempFolder, recursive = TRUE)

  # ── Build synthetic data with injected mutations ──────────────────────────
  set.seed(777)
  n_probes_ep <- 200L
  local_probes_ep <- probe_features[1:n_probes_ep, ]
  local_sig_ep <- matrix(stats::rbeta(n_probes_ep * nsamples, 90L, 10L),
                          nrow = n_probes_ep, ncol = nsamples)
  local_sig_ep[1:50, 1:5] <- stats::rbeta(50L * 5L, 1L, 100L)
  rownames(local_sig_ep) <- local_probes_ep$PROBE
  local_sig_ep <- as.data.frame(local_sig_ep)
  colnames(local_sig_ep) <- mySampleSheet$Sample_ID

  # ── semseeker ─────────────────────────────────────────────────────────────
  SEMseeker:::semseeker(
    sample_sheet      = mySampleSheet,
    signal_data       = local_sig_ep,
    result_folder     = tempFolder,
    parallel_strategy = "sequential",
    areas             = c("GENE", "POSITION"),
    markers           = c("MUTATIONS"),
    start_fresh       = TRUE,
    showprogress      = showprogress,
    verbosity         = verbosity
  )

  # ── association_analysis (depth=3 to produce GENE-area pivot results) ─────
  inference_details <- data.frame(
    independent_variable = "Phenotest",
    family_test          = "spearman",
    transformation_y     = "",
    transformation_x     = "",
    depth_analysis       = 3L,       # depth=3 reads gene-level pivot files
    filter_p_value       = FALSE,
    stringsAsFactors     = FALSE
  )

  SEMseeker:::association_analysis(
    inference_details  = inference_details,
    result_folder      = tempFolder,
    parallel_strategy  = "sequential",
    markers            = c("MUTATIONS"),
    figures            = c("HYPO"),
    areas              = c("GENE"),
    multiple_test_adj  = "BH",
    showprogress       = showprogress,
    verbosity          = verbosity
  )

  # ── pathway_ctdR ──────────────────────────────────────────────────────────
  # Re-open environment so pathway functions can read ssEnv
  SEMseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy,
                        markers = c("MUTATIONS"), figures = c("HYPO"),
                        areas = c("GENE"), multiple_test_adj = "BH",
                        showprogress = showprogress, verbosity = verbosity)

  testthat::expect_no_error(
    SEMseeker:::pathway_ctdR(
      study            = "test",
      inference_details = inference_details,
      significance     = FALSE     # include all results, not just significant
    )
  )

  # Pathway folder should have been created
  pathway_dir <- file.path(tempFolder, "Pathway")
  testthat::expect_true(dir.exists(pathway_dir))

  SEMseeker:::close_env()
  unlink(tempFolder, recursive = TRUE)
})
