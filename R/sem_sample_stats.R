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

  stats  <- NULL
  scopes <- .sem_stats_scopes_get()
  core_log_event("INFO: ", format(Sys.time(), "%a %b %d %X %Y"),
            " Per-sample statistics, scopes: ",
            paste(vapply(scopes, function(s) s$scope, character(1)), collapse = ", "), ".")
  for (scope in scopes) {
    # AI-248: one call covers every marker of the scope and, for each, every
    # aggregation it carries — the signal descriptors included, since SIGNAL is
    # a marker like the others and its figure is the scale of the values.
    scope_stats <- .sem_sample_burden_get(scope$area, scope$subarea)
    if (is.null(scope_stats))
      next
    stats <- if (is.null(stats)) scope_stats else
      merge(stats, scope_stats, by = "Sample_ID", all = TRUE)
  }

  if (is.null(stats)) {
    core_log_event("WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
              " No per-sample statistics could be computed; ",
              "SAMPLE_STATS_RESULT not written.")
    return(invisible(NULL))
  }

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

  # SIGNAL is not an analysis you opt into: it is the data the run was given,
  # and its descriptors — plus N_PROBES, which is metadata of the scope and not
  # of any marker — must be there whatever `markers` was restricted to. Before
  # AI-248 they came from a separate path that read the signal pivot directly;
  # now they travel through the keys, so the key has to be there.
  if (!any(keys$MARKER == "SIGNAL")) {
    signal_key <- keys[1, , drop = FALSE]
    signal_key$MARKER   <- "SIGNAL"
    signal_key$FIGURE   <- io_signal_figure()
    signal_key$DISCRETE <- FALSE
    keys <- rbind(keys, signal_key)
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
    marker  <- as.character(key$MARKER)
    figure  <- as.character(key$FIGURE)

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

    # AI-248: the operator is no longer a boolean intrinsic to the marker. The
    # registry says which aggregations this pair carries; DISCRETE now only
    # decides which one is produced without being asked.
    aggregations <- util_aggregations_allowed(marker, figure,
                                              discrete = isTRUE(key$DISCRETE),
                                              default  = TRUE)
    # The number of usable positions is a property of the SCOPE, not an
    # aggregation of a marker, so it carries no marker/figure in its name. It is
    # read off the signal, which is the only marker present at every position.
    if (identical(marker, "SIGNAL"))
      aggregations <- c(aggregations, "N_PROBES")

    for (aggregation in aggregations) {
      values <- .sem_pivot_aggregate(pivot, aggregation)
      if (is.null(values))
        next
      colname <- if (identical(aggregation, "N_PROBES"))
        io_feature_colname(scope, aggregation = aggregation) else
        io_feature_colname(scope, marker, figure, aggregation)
      colnames(values) <- colname
      values$Sample_ID <- rownames(values)
      rownames(values) <- NULL

      if (is.null(result)) {
        result <- values
      } else {
        result <- result[, !(colnames(result) == colname), drop = FALSE]
        result <- merge(result, values, by = "Sample_ID", all = TRUE)
      }
    }
  }

  if (is.null(result))
    return(NULL)

  rownames(result) <- NULL
  result
}

#' Reduce every sample column of a pivot with one aggregation (internal)
#'
#' AI-248. Two paths, and the difference is not cosmetic:
#'
#' \itemize{
#'   \item `SUM` and `MEAN` are streaming reduces: polars computes them over the
#'     whole frame without materialising a single column in R;
#'   \item everything else needs the distribution — a median must sort, the two
#'     modes must estimate a density — so each sample column is collected in
#'     turn and reduced in R. One column at a time keeps memory bounded by the
#'     number of workers, not by the size of the matrix, which is the same idiom
#'     the rest of the package uses.
#' }
#'
#' @param pivot lazy polars frame, sample columns only (keys already dropped).
#' @param aggregation one name from [util_aggregations_allowed()].
#' @return a one-column data.frame, rows named by sample, or `NULL`.
#' @keywords internal
#' @noRd
.sem_pivot_aggregate <- function(pivot, aggregation) {

  if (aggregation %in% c("SUM", "MEAN")) {
    reduced <- if (identical(aggregation, "SUM")) pivot$sum() else pivot$mean()
    reduced <- reduced$with_columns(polars::pl$col("*"))
    out <- as.data.frame(t(as.data.frame(reduced$collect())))
    # A pivot with no row after the mask reduces to NA on MEAN and 0 on SUM;
    # both are honest answers to "no position in this scope".
    return(out)
  }

  samples <- names(pivot$collect_schema())
  if (length(samples) == 0)
    return(NULL)

  collected <- as.data.frame(pivot$collect())
  out <- data.frame(value = vapply(samples, function(s)
    as.numeric(util_aggregate_values(collected[[s]], aggregation)),
    numeric(1)), stringsAsFactors = FALSE)
  rownames(out) <- samples
  out
}
