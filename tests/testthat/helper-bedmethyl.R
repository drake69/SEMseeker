## ============================================================================
## helper-bedmethyl.R
##
## Derives bedmethyl files (modkit-compatible 11-column layout) from a beta
## matrix + probe-feature data frame. Used by the cross-format convergence
## test to simulate WGBS and ONT/PacBio input from the same biological signal
## that the Illumina-array test path consumes.
##
## See R/bedmethyl_read.R for the reader and modkit pileup documentation
## (Oxford Nanopore Technologies) for the canonical column layout.
## ============================================================================

#' Write a bedmethyl file (modkit-compatible) from a beta matrix and probe coords
#'
#' @param beta_values Numeric matrix or data.frame of beta values in [0, 1]
#'   with probes in rows and samples in columns. Row names are probe IDs.
#' @param probe_features Data frame with columns PROBE, CHR, START, END,
#'   matching the rows of \code{beta_values}.
#' @param sample_id Character. Which sample column to write.
#' @param file_path Output path (will be overwritten if existing).
#' @param depth Integer. Synthetic coverage depth used to convert beta → Nmod
#'   (Nmod = round(beta * depth), Ncanonical = depth - Nmod). Default 10L.
#' @param mod_code Modification code per modkit (default "m" = 5-methylcytosine).
#' @return Invisibly, the file_path written.
make_bedmethyl_file <- function(beta_values, probe_features, sample_id,
                                file_path, depth = 10L, mod_code = "m") {
  stopifnot(sample_id %in% colnames(beta_values))
  stopifnot(all(c("PROBE", "CHR", "START", "END") %in% colnames(probe_features)))

  ## Align probe_features to beta_values via PROBE
  aligned <- probe_features[match(rownames(beta_values), probe_features$PROBE), ]
  betas   <- as.numeric(beta_values[, sample_id])

  keep <- !is.na(betas) & !is.na(aligned$CHR)
  aligned <- aligned[keep, , drop = FALSE]
  betas   <- betas[keep]

  n_mod        <- as.integer(round(betas * depth))
  n_canonical  <- as.integer(depth - n_mod)
  fraction_pct <- round(betas * 100, 3)

  ## Build the 11-column modkit layout. Columns 7-9 (start_code, end_code,
  ## color) are placeholders that the reader ignores; we put sensible values
  ## so the file is also human-inspectable.
  bedmethyl <- data.frame(
    chrom            = paste0("chr", aligned$CHR),
    start_position   = as.integer(aligned$START - 1L),  # bedmethyl is 0-based
    end_position     = as.integer(aligned$END),
    modified_base    = mod_code,
    score            = 1000L,
    strand           = "+",
    start_code       = as.integer(aligned$START - 1L),
    end_code         = as.integer(aligned$END),
    color            = "255,0,0",
    N_valid_cov      = as.integer(depth),
    fraction_modified = fraction_pct,
    stringsAsFactors = FALSE
  )

  ## Sort by chr/start for predictability
  bedmethyl <- bedmethyl[order(bedmethyl$chrom, bedmethyl$start_position), ]

  utils::write.table(
    bedmethyl, file = file_path, sep = "\t",
    quote = FALSE, row.names = FALSE, col.names = FALSE
  )
  invisible(file_path)
}


#' Convert a beta matrix into a set of bedmethyl files (one per sample)
#'
#' Convenience wrapper that writes one bedmethyl file per sample into
#' \code{out_dir} and returns the named vector of file paths.
#'
#' @param beta_values Beta matrix (probes × samples).
#' @param probe_features As in \code{make_bedmethyl_file}.
#' @param out_dir Output directory (created if missing).
#' @param depth Integer coverage depth. Default 10.
#' @return Named character vector \code{sample_id → file_path}.
make_bedmethyl_per_sample <- function(beta_values, probe_features, out_dir,
                                      depth = 10L) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  samples <- colnames(beta_values)
  paths <- vapply(samples, function(sid) {
    p <- file.path(out_dir, paste0(sid, ".bedmethyl"))
    make_bedmethyl_file(beta_values, probe_features, sample_id = sid,
                        file_path = p, depth = depth)
    p
  }, FUN.VALUE = character(1L))
  names(paths) <- samples
  paths
}
