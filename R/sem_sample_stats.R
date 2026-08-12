#' Per-sample statistics sibling of the sample sheet (internal)
#'
#' AI-223, slice 1 — scope `SAMPLE`. Produces `SAMPLE_STATS_RESULT.csv` next to
#' `SAMPLE_SHEET_RESULT.csv` in the data folder, one row per sample, wide:
#'
#' \itemize{
#'   \item burden per marker — `SAMPLE_<MARKER>_<FIGURE>`, aggregated over the
#'     `AREA == "POSITION"` keys (i.e. every probe, no genomic filter). The
#'     operator is intrinsic to the marker: discrete markers are summed,
#'     continuous ones averaged;
#'   \item signal descriptors — `SAMPLE_MEDIAN`, `SAMPLE_MEAN`,
#'     `SAMPLE_VARIANCE`, `SAMPLE_IQR`, plus `SAMPLE_MODE_LOW` /
#'     `SAMPLE_MODE_HIGH` on the beta scale only;
#'   \item `SAMPLE_N_PROBES`.
#' }
#'
#' Keeping these out of the sample sheet is deliberate: the sheet stays the
#' description of the study, this file carries what the run computed, and the
#' two join on `Sample_ID`.
#'
#' AI-223 slice 2a — region scopes. `sample_stats_scopes` (a `semseeker()`
#' argument, default `"SAMPLE"`) adds the burden restricted to one registered
#' `(AREA, SUBAREA)` pair, e.g. `"GENE_TSS1500"` produces
#' `GENE_TSS1500_<MARKER>_<FIGURE>`. The `SAMPLE` scope is always produced,
#' whatever the argument says: the depth=1 consumer relies on it.
#'
#' @return Invisibly the path of the file written, or `NULL` when there was
#'   nothing to write.
#' @keywords internal
#' @noRd
#' @importFrom doRNG %dorng%
sem_sample_stats_build <- function() {

  ssEnv <- core_get_session_info()

  burden <- NULL
  scopes <- .sem_stats_scopes_get()
  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
            " Per-sample statistics, scopes: ",
            paste(vapply(scopes, function(s) s$scope, character(1)), collapse = ", "), ".")
  for (scope in scopes) {
    scope_burden <- .sem_sample_burden_get(scope$area, scope$subarea)
    if (is.null(scope_burden))
      next
    burden <- if (is.null(burden)) scope_burden else
      merge(burden, scope_burden, by = "Sample_ID", all = TRUE)
  }

  # Signal descriptors stay at sample level: the region-scoped descriptors are
  # a separate decision (they are not marker/figure pairs, so they do not
  # travel through the keys machinery the way the burden does).
  descriptors <- .sem_sample_descriptors_get()

  if (is.null(burden) && is.null(descriptors)) {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " No per-sample statistics could be computed; ",
              "SAMPLE_STATS_RESULT not written.")
    return(invisible(NULL))
  }

  stats <- if (is.null(burden)) descriptors else
           if (is.null(descriptors)) burden else
           merge(descriptors, burden, by = "Sample_ID", all = TRUE)

  # Sample sheet identifiers are the reference: report what the pivots know
  # about samples the sheet does not mention, and vice versa (AI-083).
  sheet <- sem_study_summary_get(with_sample_stats = FALSE)
  if (!is.null(sheet) && "Sample_ID" %in% colnames(sheet)) {
    ids_sheet <- unique(as.character(sheet$Sample_ID))
    ids_stats <- unique(as.character(stats$Sample_ID))
    if (length(intersect(ids_sheet, ids_stats)) == 0) {
      core_log_event("ERROR: ", format(Sys.time(), "%a %b %d %X %Y"),
                " Per-sample statistics match no sample of the sample sheet: [",
                paste(utils::head(ids_sheet, 4), collapse = ", "), "] vs [",
                paste(utils::head(ids_stats, 4), collapse = ", "), "].")
      stop("sem_sample_stats_build(): no Sample_ID in common between the sample ",
           "sheet and the computed pivots. Sample sheet ids: ",
           paste(utils::head(ids_sheet, 4), collapse = ", "),
           " | pivot columns: ", paste(utils::head(ids_stats, 4), collapse = ", "), ".")
    }
    missing_ids <- setdiff(ids_sheet, ids_stats)
    if (length(missing_ids) > 0)
      core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
                " ", length(missing_ids), " sample(s) have no computed statistics: ",
                paste(utils::head(missing_ids, 10), collapse = ", "),
                if (length(missing_ids) > 10) ", ..." else "")
  }

  stats <- stats[order(stats$Sample_ID), , drop = FALSE]
  file_path <- io_file_path_build(ssEnv$result_folderData, "sample_stats_result", "csv")
  utils::write.csv2(stats, file_path, row.names = FALSE)
  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
            " Per-sample statistics written: ", nrow(stats), " samples x ",
            ncol(stats) - 1, " columns.")

  invisible(file_path)
}

#' Scopes the sibling must carry, resolved against the keys of the run
#'
#' AI-223 slice 2a. `ssEnv$sample_stats_scopes` names scopes the way they
#' appear in the column prefixes (`SAMPLE`, `GENE_TSS1500`, …); this resolves
#' each one back to the `(AREA, SUBAREA)` pair that defines it, using
#' `ssEnv$keys_areas_subareas$COMBINED` as the registry.
#'
#' A scope that does not resolve is an error, not a silent drop: a run that
#' quietly ignores the argument would look identical to one that honoured it,
#' and the researcher would only find out when the column is missing at
#' analysis time.
#'
#' @return list of `list(scope, area, subarea)`, `SAMPLE` always first.
#' @keywords internal
#' @noRd
.sem_stats_scopes_get <- function() {

  ssEnv <- core_get_session_info()
  requested <- ssEnv$sample_stats_scopes
  if (is.null(requested) || length(requested) == 0)
    requested <- "SAMPLE"
  # SAMPLE is not optional: the depth=1 consumer and the historical burden
  # columns both live there.
  requested <- unique(c("SAMPLE", core_name_cleaning(as.character(requested))))

  registry <- ssEnv$keys_areas_subareas
  available <- if (is.null(registry) || !("COMBINED" %in% colnames(registry)))
    character(0) else core_name_cleaning(as.character(registry$COMBINED))

  scopes <- list()
  for (name in requested) {
    if (identical(name, "SAMPLE")) {
      scopes[[length(scopes) + 1L]] <- list(scope = "SAMPLE", area = NULL, subarea = NULL)
      next
    }
    idx <- which(available == name)
    if (length(idx) == 0)
      stop("sample_stats_scopes: unknown scope '", name,
           "'. A scope is an (AREA, SUBAREA) pair of this run. Available: ",
           paste(c("SAMPLE", available), collapse = ", "), ".")
    scopes[[length(scopes) + 1L]] <- list(
      scope   = name,
      area    = as.character(registry$AREA[idx[1]]),
      subarea = as.character(registry$SUBAREA[idx[1]]))
  }
  scopes
}

#' Probe mask of a region scope
#'
#' The set of positions that belong to the scope, one row per position. Built
#' from [anno_probe_features_get()] — the same annotation the pivots use — and
#' deliberately NOT from the annotated `<AREA>/<SUBAREA>` pivot: annotation
#' explodes probes that map to several genes, so summing that pivot would count
#' a probe once per gene. Restricting the POSITION pivot to this mask counts
#' every position exactly once, which is what a per-sample burden restricted to
#' a region class means.
#'
#' @return a lazy polars frame with CHR/START/END normalised the way the
#'   POSITION pivot stores them, or `NULL` when the scope selects nothing.
#' @keywords internal
#' @noRd
.sem_scope_probe_mask <- function(area, subarea) {

  subarea <- if (is.null(subarea) || !nzchar(as.character(subarea))) "WHOLE" else
    as.character(subarea)
  area_subarea <- core_name_cleaning(paste0(as.character(area), "_", subarea))

  probe_features <- anno_probe_features_get(area_subarea)
  if (is.null(probe_features) || nrow(probe_features) == 0)
    return(NULL)
  if (!(area_subarea %in% colnames(probe_features)))
    stop(".sem_scope_probe_mask(): anno_probe_features_get('", area_subarea,
         "') returned no '", area_subarea, "' column.")

  membership <- as.character(probe_features[[area_subarea]])
  keep <- !is.na(membership) & nzchar(trimws(membership))
  mask <- unique(probe_features[keep, c("CHR", "START", "END"), drop = FALSE])
  if (nrow(mask) == 0)
    return(NULL)

  # The POSITION pivot stores CHR without the "chr" prefix and the coordinates
  # as Int32 (anno_create_position_pivots / io_stream_merge_bed); the manifest
  # may carry either form.
  mask$CHR   <- as.character(mask$CHR)
  mask$START <- as.integer(mask$START)
  mask$END   <- as.integer(mask$END)
  mask <- mask[!is.na(mask$START) & !is.na(mask$END), , drop = FALSE]

  polars::as_polars_df(mask)$lazy()$with_columns(
    polars::pl$col("CHR")$cast(polars::pl$String)$str$replace("^(?i)chr", ""),
    polars::pl$col("START")$cast(polars::pl$Int32),
    polars::pl$col("END")$cast(polars::pl$Int32)
  )
}

#' Burden per marker for one scope
#'
#' Moved here from sem_study_summary_total() (AI-223): the columns it used to
#' append to the sample sheet were aggregated over the `POSITION` keys, i.e.
#' they were already the `SAMPLE` scope of this artefact. AI-223 slice 2a adds
#' the region scopes, computed from the same POSITION pivots restricted to the
#' scope's probe mask.
#'
#' @param area,subarea `NULL` for the whole sample (scope `SAMPLE`), otherwise
#'   the registered pair that defines the region scope.
#' @keywords internal
#' @noRd
.sem_sample_burden_get <- function(area = NULL, subarea = NULL) {

  ssEnv <- core_get_session_info()
  keys <- ssEnv$keys_areas_subareas_markers_figures
  keys <- subset(keys, keys$AREA == "POSITION")
  if (nrow(keys) == 0) {
    # Every burden, whatever its scope, is aggregated from the POSITION pivots.
    # A `subareas` argument that excludes WHOLE drops them (util_keys_create),
    # and the sibling would come out empty without saying why.
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " No POSITION key in this run: no burden can be aggregated. ",
              "Include \"WHOLE\" in `subareas` if you restricted it.")
    return(NULL)
  }

  scope <- io_scope_name(area = area, subarea = subarea)
  mask  <- if (identical(scope, "SAMPLE")) NULL else .sem_scope_probe_mask(area, subarea)
  if (!identical(scope, "SAMPLE") && is.null(mask)) {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " Scope ", scope, " selects no position; its burden columns are ",
              "not written.")
    return(NULL)
  }

  result <- NULL
  for (k in seq_len(nrow(keys))) {
    key     <- keys[k, ]
    marker  <- key$MARKER
    figure  <- key$FIGURE

    pivot <- io_read_pivot(marker, figure, key$AREA, key$SUBAREA)
    if (is.null(pivot))
      next

    if (!is.null(mask)) {
      pivot <- pivot$with_columns(
        polars::pl$col("CHR")$cast(polars::pl$String)$str$replace("^(?i)chr", ""),
        polars::pl$col("START")$cast(polars::pl$Int32),
        polars::pl$col("END")$cast(polars::pl$Int32)
      )
      pivot <- pivot$join(mask, on = c("CHR", "START", "END"), how = "semi")
    }

    pivot <- pivot$drop(c("CHR", "START", "END"))
    # The operator is intrinsic to the marker: counts are summed, continuous
    # deviations averaged. DISCRETE is the flag that already encodes it
    # (util_keys_create.R).
    pivot <- if (isTRUE(key$DISCRETE)) pivot$sum() else pivot$mean()
    pivot <- pivot$with_columns(polars::pl$col("*"))

    pivot <- as.data.frame(t(as.data.frame(pivot$collect())))
    colname <- io_burden_colname(scope, marker, figure)
    colnames(pivot) <- colname
    pivot$Sample_ID <- rownames(pivot)

    if (is.null(result)) {
      result <- pivot
    } else {
      result <- result[, !(colnames(result) == colname), drop = FALSE]
      result <- merge(result, pivot, by = "Sample_ID", all = TRUE)
    }
  }

  if (is.null(result))
    return(NULL)

  result[is.na(result)] <- 0
  rownames(result) <- NULL
  result
}

#' Signal descriptors at sample level
#'
#' One pass per sample over the SIGNAL pivot, parallelised with the same
#' foreach/%dorng% idiom used elsewhere in the package. Each worker collects a
#' single sample column (one column ~ n_probes doubles, so memory stays bounded
#' by the number of workers, not by the size of the matrix) and reduces it to
#' one row of statistics.
#'
#' @keywords internal
#' @noRd
.sem_sample_descriptors_get <- function() {

  ssEnv <- core_get_session_info()

  pivot_path <- io_pivot_file_name_parquet("SIGNAL", "MEAN", "PROBE", "WHOLE")
  if (!file.exists(pivot_path)) {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " SIGNAL pivot not found, signal descriptors skipped: ", pivot_path)
    return(NULL)
  }

  lazy <- polars::pl$scan_parquet(pivot_path)
  # The probe identifier column is named AREA in the PROBE pivot
  # (io_signal_save()); every other column is a sample.
  samples <- setdiff(names(lazy$collect_schema()), c("AREA", "PROBE", "CHR", "START", "END"))
  if (length(samples) == 0) {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " SIGNAL pivot carries no sample column; descriptors skipped.")
    return(NULL)
  }

  beta  <- isTRUE(ssEnv$beta)
  stats_wanted <- io_signal_stats(beta)
  colnames_out <- vapply(stats_wanted, function(s) io_stat_colname("SAMPLE", s),
                         character(1))

  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
            " Computing signal descriptors for ", length(samples),
            " samples (scale: ", if (beta) "beta" else "M-value", ").")

  s <- NULL   # quiet R CMD check note
  descriptors <- foreach::foreach(
    s         = seq_along(samples),
    .combine  = rbind,
    .packages = c("SEMseeker", "polars", "stats"),
    .export   = c("samples", "pivot_path", "beta", "stats_wanted", "colnames_out")
  ) %dorng% tryCatch({
    sample_id <- samples[s]
    values <- as.data.frame(
      polars::pl$scan_parquet(pivot_path)$select(sample_id)$collect()
    )[[1]]
    values <- as.numeric(values)
    values <- values[is.finite(values)]

    row <- SEMseeker:::util_signal_descriptors(values, beta = beta)
    out <- data.frame(Sample_ID = sample_id, stringsAsFactors = FALSE)
    for (i in seq_along(stats_wanted)) {
      # NULL would DROP the column instead of filling it, and rows of unequal
      # width break the rbind combine.
      value <- row[[stats_wanted[i]]]
      out[[colnames_out[i]]] <- if (is.null(value)) NA_real_ else value
    }
    out
  }, error = function(e) {
    NULL
  })

  if (is.null(descriptors) || nrow(descriptors) == 0)
    return(NULL)

  rownames(descriptors) <- NULL
  descriptors
}
