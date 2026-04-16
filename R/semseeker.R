#' Run SEMseeker on methylation data from any supported source
#'
#' Public entry point. Accepts a wide range of inputs — bedmethyl files
#' (modkit / nanopolish), coordinate-based data frames (WGBS / long-read),
#' Illumina probe-indexed matrices, or already-loaded data frames — normalises
#' them to the internal format, optionally converts M-values to beta values,
#' validates the tech / genome_build combination, and delegates to the core
#' pipeline \code{\link{semseeker_core}}.
#'
#' Supported \code{input} forms:
#' \itemize{
#'   \item Character vector of bedmethyl file paths (\code{.bed}/\code{.tsv}/
#'     \code{.bedmethyl}) — parsed via \code{\link{bedmethyl_read}}.
#'   \item Data frame with \code{CHR}/\code{START}[\code{/END}] columns
#'     (WGBS / long-read coordinate format) — normalised via
#'     \code{\link{normalize_signal_input}}.
#'   \item Matrix or data frame with probe-ID rownames (Illumina array) —
#'     passed through unchanged.
#'   \item List of any of the above, one element per batch.
#' }
#'
#' @param input Input signal data. See details for supported forms.
#' @param sample_sheet Data frame (or list of data frames) with a
#'   \code{Sample_ID} column identifying samples.
#' @param result_folder Output directory.
#' @param input_type One of \code{"auto"} (default), \code{"bedmethyl"},
#'   \code{"coord_df"}, \code{"matrix"}. When \code{"auto"}, the type is
#'   inferred from the object class / file extension.
#' @param tech Optional technology label (\code{"K850"}, \code{"K450"},
#'   \code{"K27"}, \code{"WGBS"}, \code{"LONGREAD"}). If \code{NULL} (default)
#'   it is auto-detected downstream by \code{get_meth_tech()}.
#' @param genome_build Reference genome build. One of \code{"hg19"} (default),
#'   \code{"hg38"}, \code{"mm10"}, \code{"legacy"}.
#' @param auto_convert_mvalues If \code{TRUE} (default) and the input values
#'   exceed the \code{[0, 1]} range, they are assumed to be M-values and
#'   converted to beta via \eqn{\beta = 2^M / (1 + 2^M)}. Set to \code{FALSE}
#'   to disable the check (the downstream pipeline will then run on raw values).
#' @param strict_build_check If \code{TRUE} (default), impossible tech /
#'   genome_build combinations (e.g. \code{tech="LONGREAD"} with
#'   \code{genome_build="hg19"}) raise an error. If \code{FALSE}, a warning is
#'   emitted and the pipeline proceeds.
#' @param ... Additional arguments forwarded to \code{\link{semseeker_core}}
#'   and \code{init_env()} (e.g. \code{parallel_strategy}, \code{alpha},
#'   \code{sliding_window_size}, \code{marker}, \code{areas}).
#'
#' @return Invisibly \code{NULL}; writes output files to \code{result_folder}.
#'
#' @examples
#' \dontrun{
#' # Bedmethyl (Nanopore / modkit):
#' semseeker(
#'   input = list.files("bedmethyl/", pattern = "\\.bed$", full.names = TRUE),
#'   sample_sheet = ss,
#'   result_folder = tempdir(),
#'   tech = "LONGREAD",
#'   genome_build = "hg38"
#' )
#'
#' # Illumina beta matrix:
#' semseeker(
#'   input = beta_matrix,
#'   sample_sheet = ss,
#'   result_folder = tempdir()
#' )
#'
#' # WGBS coordinate data frame:
#' semseeker(
#'   input = wgbs_df,              # columns: CHR, START, END, sample1, ...
#'   sample_sheet = ss,
#'   result_folder = tempdir(),
#'   tech = "WGBS"
#' )
#' }
#'
#' @export
semseeker <- function(input,
                      sample_sheet,
                      result_folder,
                      input_type = c("auto", "bedmethyl", "coord_df", "matrix"),
                      tech = NULL,
                      genome_build = "hg19",
                      auto_convert_mvalues = TRUE,
                      strict_build_check = TRUE,
                      ...) {

  input_type <- match.arg(input_type)

  # ---- Step 1: validate tech × genome_build ----------------------------
  .validate_tech_build(tech, genome_build, strict_build_check)

  # ---- Step 2: normalise input to SEMseeker signal_data ----------------
  # Recurse for list inputs (multi-batch), unless the character vector is a
  # list of bedmethyl paths (handled as a single batch inside .dispatch_one).
  signal_data <- .dispatch_one(input, input_type, auto_convert_mvalues)

  # ---- Step 3: delegate to core pipeline -------------------------------
  semseeker_core(
    sample_sheet  = sample_sheet,
    signal_data   = signal_data,
    result_folder = result_folder,
    tech          = if (is.null(tech)) "" else tech,
    genome_build  = genome_build,
    ...
  )
}

# --- Internal helpers --------------------------------------------------------

#' @keywords internal
.dispatch_one <- function(input, input_type, auto_convert_mvalues) {

  # Auto-detect when requested
  if (identical(input_type, "auto")) {
    input_type <- .detect_input_type(input)
  }

  signal_data <- switch(
    input_type,
    bedmethyl = bedmethyl_read(input),
    coord_df  = normalize_signal_input(input),
    matrix    = input,
    stop(".dispatch_one(): unknown input_type '", input_type, "'.")
  )

  # coord_df path: normalize_signal_input() returns probe-ID-indexed df already
  # bedmethyl path: returns coord df → also run through normalize_signal_input
  if (identical(input_type, "bedmethyl")) {
    signal_data <- normalize_signal_input(signal_data)
  }

  # ---- M-value detection & conversion --------------------------------
  if (isTRUE(auto_convert_mvalues) && .looks_like_mvalues(signal_data)) {
    message("semseeker(): detected M-values (|x| > 1); converting to beta via 2^M / (1 + 2^M).")
    signal_data <- mvalue_to_beta(signal_data, coord_cols = c("CHR","START","END"))
  }

  signal_data
}

#' @keywords internal
.detect_input_type <- function(input) {
  # List of inputs → apply per element; dispatcher delegates one-by-one
  # (semseeker_core itself iterates lists, so we just need each element
  # normalised before being wrapped again.)
  if (is.character(input)) {
    exts <- tolower(tools::file_ext(input))
    if (all(exts %in% c("bed", "tsv", "bedmethyl", "gz")))
      return("bedmethyl")
    stop(".detect_input_type(): character input with unsupported extensions: ",
         paste(unique(exts), collapse = ", "))
  }
  if (is.data.frame(input) && is_coord_format(input))
    return("coord_df")
  if (is.matrix(input) || is.data.frame(input))
    return("matrix")
  stop(".detect_input_type(): cannot infer input_type for object of class '",
       class(input)[1L], "'. Pass input_type explicitly.")
}

#' @keywords internal
.looks_like_mvalues <- function(signal_data) {
  coord_cols <- c("CHR", "START", "END")
  numeric_cols <- if (is.data.frame(signal_data)) {
    setdiff(colnames(signal_data), coord_cols)
  } else {
    colnames(signal_data)
  }
  if (!length(numeric_cols)) return(FALSE)

  sample_rows <- seq_len(min(10000L, nrow(signal_data)))
  chunk <- if (is.data.frame(signal_data)) {
    as.matrix(signal_data[sample_rows, numeric_cols, drop = FALSE])
  } else {
    signal_data[sample_rows, , drop = FALSE]
  }
  mx <- suppressWarnings(max(abs(chunk), na.rm = TRUE))
  is.finite(mx) && mx > 1
}

#' @keywords internal
.validate_tech_build <- function(tech, genome_build, strict) {
  if (is.null(tech) || !nzchar(tech)) return(invisible(NULL))

  bad <- (identical(tech, "LONGREAD") && identical(genome_build, "hg19"))

  if (bad) {
    msg <- paste0(
      "tech = 'LONGREAD' with genome_build = 'hg19' is almost certainly wrong: ",
      "long-read methylation pipelines (Nanopore / PacBio) align to GRCh38. ",
      "Set genome_build = 'hg38' in init_env() / semseeker()."
    )
    if (isTRUE(strict)) stop(msg) else warning(msg)
  }

  invisible(NULL)
}
