# AI-061: lazy / low-memory batch path for limma_<N> and voom_<N>.
#
# Why a separate function: apply_stat_model_batch.R receives a fully
# materialised, transposed, sample-merged tempDataFrame that the caller
# (run_depth_n_marker) builds with the standard
#   tempDataFrame <- as.data.frame(pivot_lazy$collect())
#   ... t() ... merge(sample_sheet, ...)
# pipeline. Even at GENE/WHOLE scale (~17 k genes × ~500 samples
# ≈ 70 MB raw) that path peaks around 4×–5× the raw size in R memory
# alone, plus the polars Rust DataFrame that stays alive in parallel,
# plus Arrow IPC buffers from the collect step.  Smoke SIGNAL@PROBE
# (366 k × 500 ≈ 1.4 GB raw) ran to 45 GB before lmFit was even
# called and was OOM-killed at ~60 GB.
#
# This function takes the pivot's polars LazyFrame directly, applies
# the AI-043 area_to_remove filter LAZILY, materialises only ONE
# R matrix (genes × samples) without an intervening data.frame, and
# explicitly rm()+gc()'s the polars-side artifacts before lmFit so
# the peak working set drops to roughly:
#   polars collect      ≈ 1×  raw  (in Rust heap, drops after rm+gc)
#   R y_mat             ≈ 1×  raw  (the single matrix lmFit sees)
#   design + lmFit work  small
# i.e. ~3×–4× the raw size at peak vs ~30× on the legacy path.

#' Polars-lazy batch fit for limma_<N> and voom_<N>.
#'
#' @param pivot_lazy polars LazyFrame. The annotated AREA-level pivot
#'   for one (MARKER, FIGURE, AREA, SUBAREA) key, with at least an
#'   AREA column and one column per sample.
#' @param sample_sheet data.frame. Per-sample metadata returned by
#'   prepare_study_for_analysis(); must contain Sample_ID, the
#'   independent variable, and any covariates referenced.
#' @param family_test character, of the form 'limma_<degree>' or
#'   'voom_<degree>' (optional '_<partition>' suffix is ignored).
#' @param covariates character vector of covariate column names.
#' @param key data.frame row carrying MARKER/FIGURE/AREA/SUBAREA.
#' @param transformation_y character label passed through unchanged.
#' @param independent_variable character, single column name.
#' @param area_to_remove character. AREA values to skip (resume cache).
#' @param ... ignored (kept for caller-symmetry).
#'
#' @return data.frame with one row per kept area, same schema
#'   apply_stat_model_batch() returns. NULL on degenerate input.
#'
#' @keywords internal
#' @noRd
apply_stat_model_batch_lazy <- function(pivot_lazy,
                                         sample_sheet,
                                         family_test,
                                         covariates = NULL,
                                         key,
                                         transformation_y,
                                         independent_variable,
                                         area_to_remove = character(0),
                                         ...) {

  # Parser: <engine>_<degree>[_<partition>]
  parts  <- unlist(strsplit(as.character(family_test), "_"))
  engine <- parts[1]
  if (!engine %in% c("limma", "voom") || length(parts) < 2L) {
    log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"),
              " apply_stat_model_batch_lazy: malformed family_test='",
              family_test, "'")
    return(NULL)
  }
  degree <- suppressWarnings(as.integer(parts[2]))
  if (is.na(degree) || degree < 1L) {
    log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"),
              " apply_stat_model_batch_lazy: invalid degree in '",
              family_test, "'")
    return(NULL)
  }
  partition_percentage <- if (length(parts) >= 3L)
    suppressWarnings(as.numeric(parts[3])) else 1
  if (is.na(partition_percentage)) partition_percentage <- 1

  # Drop genomic-coord and probe-metadata columns LAZILY — no
  # materialisation yet.
  schema_names <- names(pivot_lazy$collect_schema())
  drop_cols <- intersect(schema_names,
                          c("CHR", "START", "END", "PROBE",
                            "K27", "K450", "K850"))
  if (length(drop_cols) > 0L) {
    pivot_lazy <- pivot_lazy$drop(drop_cols)
  }

  # AI-043 resume: filter out areas already in the on-disk CSV. The
  # caller has already applied gsub("-","_") to the area_to_remove
  # values, so we apply the same normalisation to the lazy AREA column
  # before the membership check.
  if (length(area_to_remove) > 0L) {
    pivot_lazy <- pivot_lazy$filter(
      !polars::pl$col("AREA")$str$replace_all("-", "_")$is_in(area_to_remove)
    )
  }

  # Materialise once. The polars DF stays in Rust heap until we drop it
  # explicitly below; that prevents both heaps holding the data
  # simultaneously.
  pivot_df <- pivot_lazy$collect()
  n_genes  <- pivot_df$height
  if (n_genes == 0L) {
    log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
              " apply_stat_model_batch_lazy: nothing left after resume",
              " filter for ", key$MARKER, "/", key$FIGURE, "/",
              key$AREA, "/", key$SUBAREA, ".")
    return(NULL)
  }

  # Pull the AREA column out as an R character vector (small), then
  # work on the value-only DataFrame.
  # NOTE: polars R 1.x exposes series -> R via as.character() / as.vector()
  # directly; there is no $to_r() method — earlier draft using it crashed
  # at runtime with rlang::abort.
  gene_names <- as.character(pivot_df$select("AREA")$to_series())
  pivot_vals <- pivot_df$drop("AREA")
  sample_cols <- names(pivot_vals)
  rm(pivot_df)

  # Convert column-by-column to the final numeric matrix. as.vector() on a
  # numeric polars Series returns an R double vector; we never go through
  # an intermediate data.frame for the value side. The only sample-by-gene
  # allocation alive at once is y_mat itself.
  y_mat <- matrix(NA_real_, nrow = n_genes, ncol = length(sample_cols),
                   dimnames = list(gene_names, sample_cols))
  for (i in seq_along(sample_cols)) {
    y_mat[, i] <- as.vector(pivot_vals$select(sample_cols[i])$to_series())
  }
  rm(pivot_vals)
  gc(verbose = FALSE)
  y_mat[is.na(y_mat)] <- 0

  # Align sample_sheet rows to y_mat columns. Drop samples that don't
  # appear in sample_sheet OR have NA in IV/covariates.
  ss <- sample_sheet
  ss <- ss[match(colnames(y_mat), as.character(ss$Sample_ID)), , drop = FALSE]

  use_iv <- c(independent_variable, covariates)
  use_iv <- use_iv[nzchar(use_iv) & use_iv %in% colnames(ss)]
  keep <- !is.na(ss$Sample_ID) &
          stats::complete.cases(ss[, use_iv, drop = FALSE])
  if (sum(keep) < (degree + length(covariates) + 2L)) {
    log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " apply_stat_model_batch_lazy: too few complete samples (",
              sum(keep), " < ", degree + length(covariates) + 2L, ")")
    return(NULL)
  }
  ss <- ss[keep, , drop = FALSE]
  y_mat <- y_mat[, keep, drop = FALSE]

  # Design matrix: poly(IV, degree, raw=TRUE) + covariates.
  # Name the polynomial columns the same way association_model_polynomial's
  # I(...) formula winds up after name_cleaning — 'I_<IV>_<deg>' — so the
  # CSV columns landed by the build_pname/build_ename gsub chain below
  # match the polynomial CSV schema bit-for-bit instead of carrying the
  # ugly long-form 'STATS_POLY_EVAL_PARSE_TEXT_EQ_<IV>_EQ_RAW_EQ_TRUE_<n>'.
  iv_vec <- as.numeric(ss[, independent_variable])
  poly_mat <- stats::poly(iv_vec, degree, raw = TRUE)
  colnames(poly_mat) <- paste0("I_", independent_variable, "_", seq_len(degree))

  cov_used <- character(0)
  if (length(covariates) > 0L && any(nzchar(covariates))) {
    cov_used <- intersect(covariates[nzchar(covariates)], colnames(ss))
  }
  if (length(cov_used) > 0L) {
    cov_mat <- as.matrix(ss[, cov_used, drop = FALSE])
    design  <- cbind(`(Intercept)` = 1, poly_mat, cov_mat)
  } else {
    design <- cbind(`(Intercept)` = 1, poly_mat)
  }

  # Fit. y_mat is the only sample-by-gene matrix we hold; everything
  # else (design, fit$coefficients, fit$p.value) is at most degree+1+|cov|
  # columns wide.
  if (engine == "voom") {
    voom_obj <- tryCatch(limma::voom(y_mat, design),
                          error = function(e) {
                            log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"),
                                      " voom failed: ", conditionMessage(e))
                            NULL
                          })
    if (is.null(voom_obj)) return(NULL)
    fit <- limma::lmFit(voom_obj, design)
    r_model_label <- "limma::voom+lmFit+eBayes"
    rm(voom_obj); gc(verbose = FALSE)
  } else {
    fit <- tryCatch(limma::lmFit(y_mat, design),
                     error = function(e) {
                       log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"),
                                 " lmFit failed: ", conditionMessage(e))
                       NULL
                     })
    if (is.null(fit)) return(NULL)
    r_model_label <- "limma::lmFit+eBayes"
  }
  area_cols <- rownames(y_mat)
  rm(y_mat); gc(verbose = FALSE)

  fit <- limma::eBayes(fit)

  # Build the result data.frame in the same schema apply_stat_model_batch
  # returns — same column names, so the FDR + selector machinery in
  # association_analysis_save_results() picks them up unchanged.
  coef_names <- colnames(fit$coefficients)

  build_pname <- function(cn) {
    pn <- name_cleaning(paste0(cn, "_pvalue"))
    pn <- name_cleaning(gsub("_STATS_POLY_EVAL_PARSE_TEXT_EQ", "", pn))
    pn <- name_cleaning(gsub("_RAW_EQ_TRUE", "", pn))
    pn <- name_cleaning(gsub("INDEPENDENT_VARIABLE", independent_variable, pn))
    pn
  }
  build_ename <- function(cn) {
    en <- name_cleaning(paste0(cn, "_estimate"))
    en <- name_cleaning(gsub("_STATS_POLY_EVAL_PARSE_TEXT_EQ", "", en))
    en <- name_cleaning(gsub("_RAW_EQ_TRUE", "", en))
    en <- name_cleaning(gsub("INDEPENDENT_VARIABLE", independent_variable, en))
    en
  }
  pnames <- vapply(coef_names, build_pname, character(1))
  enames <- vapply(coef_names, build_ename, character(1))

  cov_label <- if (length(cov_used) > 0L) paste(cov_used, collapse = "+") else ""

  result_temp <- data.frame(
    MARKER               = rep(as.character(key$MARKER),  length(area_cols)),
    FIGURE               = rep(as.character(key$FIGURE),  length(area_cols)),
    AREA                 = rep(as.character(key$AREA),    length(area_cols)),
    SUBAREA              = rep(as.character(key$SUBAREA), length(area_cols)),
    AREA_OF_TEST         = area_cols,
    PL_DEGREE            = degree,
    PL_PERC              = partition_percentage,
    R_MODEL              = r_model_label,
    FAMILY_TEST          = as.character(family_test),
    TRANSFORMATION_Y     = as.character(transformation_y),
    INDEPENDENT_VARIABLE = as.character(independent_variable),
    COVARIATES           = cov_label,
    stringsAsFactors     = FALSE
  )
  for (i in seq_along(coef_names)) {
    result_temp[[pnames[i]]] <- fit$p.value[, i]
    result_temp[[enames[i]]] <- fit$coefficients[, i]
  }

  first_poly_pcol <- pnames[2L]  # coef 1 is (Intercept), coef 2 is poly_1
  if (!is.null(first_poly_pcol) && first_poly_pcol %in% colnames(result_temp)) {
    result_temp$PVALUE <- result_temp[[first_poly_pcol]]
  }

  colnames(result_temp) <- name_cleaning(colnames(result_temp))
  result_temp
}
