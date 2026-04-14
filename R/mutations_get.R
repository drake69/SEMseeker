#' mutations_get
#'
#' @param values values of methylation — data.frame with columns CHR, START, END
#'   and a fourth numeric VALUE column.
#' @param figure figure to get Mutations of HYPO or HYPER methylation
#' @param thresholds threshold to use for comparison — data.frame with columns
#'   CHR, START, END, signal_inferior_thresholds, signal_superior_thresholds.
#' @param sampleName name of the sample
#'
#' @return mutations data.frame with columns CHR, START, END, MUTATIONS (0/1).
#'   Only positions present in BOTH values and thresholds are returned
#'   (inner join on CHR, START, END).
#'
mutations_get <- function(values, figure, thresholds, sampleName) {

  # ── Coverage report ─────────────────────────────────────────────────────────
  # Compute the size of the (CHR, START, END) intersection BEFORE joining so
  # that the user can audit how many positions are shared between this sample
  # and the beta-range reference.  This is especially important for cross-run
  # analysis (e.g. a Nanopore sample vs an Illumina reference population) and
  # for long-read files (28M+ rows) where a silent positional mismatch would
  # otherwise produce wrong mutation calls or an out-of-bounds crash.
  n_input  <- nrow(values)
  n_ranges <- nrow(thresholds)

  # ── Polars inner join on (CHR, START, END) ──────────────────────────────────
  # Replaces the previous sort_by_chr_and_start + positional zip, which
  # silently produced wrong results when the two sets differed in size or
  # overlap.  Polars is mandatory here: Nanopore bedmethyl files can have
  # 28M+ rows — base-R merge()/match() is too slow at that scale.
  values_df_for_join <- data.frame(
    CHR   = as.character(values$CHR),
    START = as.integer(values$START),
    END   = as.integer(values$END),
    VALUE = values[, 4L],
    stringsAsFactors = FALSE
  )

  threshold_cols <- c("CHR", "START", "END",
                      "signal_inferior_thresholds",
                      "signal_superior_thresholds")
  thresholds_df_for_join <- thresholds[,
    intersect(threshold_cols, colnames(thresholds)),
    drop = FALSE
  ]
  thresholds_df_for_join$CHR   <- as.character(thresholds_df_for_join$CHR)
  thresholds_df_for_join$START <- as.integer(thresholds_df_for_join$START)
  thresholds_df_for_join$END   <- as.integer(thresholds_df_for_join$END)

  values_lf     <- polars::as_polars_df(values_df_for_join)$lazy()
  thresholds_lf <- polars::as_polars_df(thresholds_df_for_join)$lazy()

  joined <- values_lf$join(
    thresholds_lf,
    on  = c("CHR", "START", "END"),
    how = "inner"
  )$collect()

  n_covered <- joined$height

  log_event(
    "BANNER: ", format(Sys.time(), "%a %b %d %X %Y"),
    " [mutations_get] sample=", sampleName,
    " figure=", figure,
    " | input_positions=", n_input,
    " | beta_range_positions=", n_ranges,
    " | covered_by_inner_join=", n_covered
  )

  # ── Empty-result guard ──────────────────────────────────────────────────────
  if (n_covered == 0L) {
    log_event(
      "WARNING: ", format(Sys.time(), "%a %b %d %X %Y"),
      " [mutations_get] No overlapping positions between input and beta ranges",
      " — returning empty mutation result for sample=", sampleName
    )
    return(data.frame(
      CHR       = character(),
      START     = integer(),
      END       = integer(),
      MUTATIONS = integer(),
      stringsAsFactors = FALSE
    ))
  }

  # ── Mutation call ────────────────────────────────────────────────────────────
  joined_df <- as.data.frame(joined)

  if (figure == "HYPO") {
    mutation <- as.numeric(joined_df$VALUE < joined_df$signal_inferior_thresholds)
  } else {
    # HYPER
    mutation <- as.numeric(joined_df$VALUE > joined_df$signal_superior_thresholds)
  }

  mutationAnnotated <- data.frame(
    CHR       = joined_df$CHR,
    START     = joined_df$START,
    END       = joined_df$END,
    MUTATIONS = mutation,
    stringsAsFactors = FALSE
  )

  mutation_annotated_sorted <- sort_by_chr_and_start(mutationAnnotated)
  return(mutation_annotated_sorted)
}
