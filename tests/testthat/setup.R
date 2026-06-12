# num_rows <- 3e^6
# num_cols <- 5200
# populationMatrix <- as.data.frame(matrix(runif(num_rows * num_cols), nrow = num_rows, ncol = num_cols))
# library(future)
# options(future.globals.maxSize = 10 * 1024^3)
Sys.setenv(OBJC_DISABLE_INITIALIZE_FORK_SAFETY='YES')

rm(list = ls())
# DEBUG: trace setup.R progress so we can see (in the macOS run log) at
# which step GEOquery / tcltk loading is triggered. AI-017 tracking.
.trace_step <- function(msg) {
  cat(sprintf("[SETUP-TRACE %s] %s\n",
              format(Sys.time(), "%H:%M:%OS3"), msg))
  flush.console()
}
.trace_step("setup.R BEGIN")
loadNamespace("future")
.trace_step("after loadNamespace(future)")
loadNamespace("stats")
.trace_step("after loadNamespace(stats)")
use_synthetic_data <- TRUE
if (use_synthetic_data)
{
  # 2e4 keeps statistical validity while reducing setup+test time ~3x vs 6e4
  nprobes <<- 2e4
  nsamples <<- 3e1

  # intersample_mean <- 5
  # intersample_sd <- 1
  # perc_epimutation <- 0.05

  Sys.setenv(OBJC_DISABLE_INITIALIZE_FORK_SAFETY='YES')
  # Load the bundled imprinting-aware fixture (20k real EPIC probe IDs around
  # 58 curated human imprinting DMRs). Built once via
  # data-raw/build_test_master_features.R, ships in data/test_master_features.rda.
  # See ?SEMseeker::test_master_features.
  #
  # Using bundled real probe IDs guarantees:
  #  - tech-detection tests pass on every platform (probe-ID overlap → K850),
  #  - no requireNamespace() on IlluminaHumanMethylationEPICanno at test-setup
  #    time, so no minfi → GEOquery → tcltk → libX11 chain → no macOS arm64
  #    segfault path.
  .trace_step("loading bundled test_master_features fixture")
  utils::data("test_master_features", package = "SEMseeker", envir = environment())
  probe_features <- as.data.frame(test_master_features, stringsAsFactors = FALSE)
  nprobes <<- nrow(probe_features)
  .trace_step(sprintf("loaded test_master_features: %d probes", nprobes))

  # test <- stats::rnorm(nsamples,mean = intersample_mean, sd = intersample_sd)
  # mean(test)
  # plot(density(as.numeric(test)))

  set.seed(474693)
  nCpG <- nprobes * 0.9
  nTypeI <- 90
  nTypeII <- 10
  Meth <- as.data.frame(matrix(stats::rbeta(nCpG * nsamples, nTypeI, nTypeII), nrow = nCpG))

  nCpG <- nprobes * 0.1
  nTypeI2 <- 4
  nTypeII2 <- 6
  Unmeth <- as.data.frame(matrix(stats::rbeta(nCpG * nsamples, nTypeI2, nTypeII2), nrow = nCpG))

  signal_data <- rbind(Meth, Unmeth)
  # signal_data <- rep(stats::rnorm(nprobes,mean = intersample_mean, sd = intersample_sd),nsamples)

  # plot(density(as.numeric(signal_data)))

  q1 <-stats::quantile(signal_data, probs = 0.05, na.rm = TRUE)
  q3 <-stats::quantile(signal_data, probs = 0.05, na.rm = TRUE)

  # noise_position <- unique(floor(runif(nsamples*nprobes*perc_epimutation*2, min = 1, max = nsamples*nprobes)))
  # noise_position <- noise_position[1:(nsamples*nprobes*perc_epimutation)]
  # min(noise_position)
  # max(noise_position)
  # median(noise_position)
  # noise <- c(runif(nsamples*nprobes*perc_epimutation/2,min = q3,  max= q3 + 2*(abs(q3)) ),runif(nsamples*nprobes*perc_epimutation/2,min = q1-2*(abs(q1)),max=q1))
  # signal_data[noise_position] <- noise
  # min(noise)
  # max(noise)
  # signal_data <- signal_data + abs(min(noise))
  # min(signal_data)
  # signal_data <- as.data.frame(matrix(data = signal_data,nrow = nprobes,ncol = nsamples, byrow = TRUE))
  # plot(density(as.numeric(signal_data[1,])))

  row.names(signal_data) <- probe_features$PROBE

  Sample_ID <- stringi::stri_rand_strings(nsamples, 15, pattern = "[A-Za-z]")
  colnames(signal_data) <- Sample_ID
  Sample_Group <- c(rep("Control",nsamples/3),rep("Case",nsamples/3),rep("Reference",nsamples/3))
  mySampleSheet <- data.frame(Sample_Group, Sample_ID)

  mySampleSheet$Phenotest <- stats::rnorm(nsamples, mean= 1000, sd= 567)
  mySampleSheet$Group <- c(rep(TRUE,ceiling(nsamples/2)), rep(FALSE,floor(nsamples/2)))
  mySampleSheet$Covariates1 <- stats::rnorm(nsamples, mean= 567, sd= 1000)
  mySampleSheet$Covariates2 <- stats::rnorm(nsamples, mean= 67, sd= 100)

  mySampleSheet_batch <<- list(mySampleSheet, mySampleSheet, mySampleSheet)
  signal_data_batch <<- list(signal_data, signal_data, signal_data)

  q1 <- apply(signal_data, 1, function(x) stats::quantile(x, probs = 0.25, na.rm = TRUE))
  q3 <- apply(signal_data, 1, function(x) stats::quantile(x, probs = 0.75, na.rm = TRUE))
  signal_medians <- apply(signal_data, 1, stats::median)
  iqr <- data.frame(q3 - q1)

  signal_superior_thresholds <- data.frame("HIGH" = q3 + 3 * iqr)
  signal_inferior_thresholds <- data.frame("LOW" = q1 - 3 * iqr)
  colnames(signal_inferior_thresholds) <- "LOW"
  colnames(signal_superior_thresholds) <- "HIGH"

  row.names(signal_superior_thresholds) <- probe_features$PROBE
  row.names(signal_inferior_thresholds) <- probe_features$PROBE

  signal_thresholds <- data.frame("signal_median_values"=signal_medians,
    "signal_inferior_thresholds"=signal_inferior_thresholds,
    "signal_superior_thresholds"=signal_superior_thresholds,
    "iqr" = iqr,
    "q1"=q1,
    "q3"=q3)

  colnames(signal_thresholds) <- c("signal_median_values","signal_inferior_thresholds","signal_superior_thresholds","iqr","q1","q3")

  # Add genomic coordinates so mutations_get / sort_by_chr_and_start can sort thresholds
  signal_thresholds$CHR <- probe_features$CHR
  signal_thresholds$START <- probe_features$START
  signal_thresholds$END <- probe_features$END

  #####
  ### copy as global
  mySampleSheet <<- mySampleSheet
  signal_data <<- (signal_data)
  signal_medians <<- signal_medians
  signal_inferior_thresholds <<- (signal_inferior_thresholds)
  signal_superior_thresholds <<- (signal_superior_thresholds)
  signal_thresholds <<- signal_thresholds
  # perc_epimutation <<- perc_epimutation
  nsamples <<- nsamples
  nprobes <<- nprobes
  # multisession
  # sequential
}



LESIONS_BP <<- 5000L  # AI-092 + AI-044 merged: bp-based window, literature-aligned default 5 kbp (see AI-048).
bonferroni_threshold <<- 0.1
batch_id <<- 1
iqrTimes <<- 3
# "multicore" (fork) is unsafe on macOS with Polars' C++ thread pool — forked
# children can be killed by Mach exceptions.  Use "multisession" when the
# package is installed (R CMD check, CI, devtools::install()).  Fall back to
# "multicore" only when running under devtools::load_all() on non-macOS,
# because multisession workers cannot see load_all()'d internals.
if (Sys.info()[["sysname"]] == "Darwin") {
  parallel_strategy <<- "multisession"
} else {
  parallel_strategy <<- "multicore"
}
markers <<- c("MUTATIONS","DELTAQ","DELTARQ","DELTAP","DELTARP","LESIONS")

if (!use_synthetic_data) {
  signal_data <- readRDS("~/Documents/Dati_Lavoro/beckwith-wiedemann/raw/GSE95486/beta.rds")
  signal_data <- as.data.frame(signal_data)
  signal_data <- signal_data[1:1000,]
  # count rows with all missing values
  nrow_missed <- sum(apply(signal_data, 1, function(x) all(is.na(x))))
  probe_features <<- SEMseeker::probe_annotation_build("K850")[rownames(signal_data), c("CHR","START","PROBE")]
  probe_features$ABSOLUTE <- paste(probe_features$CHR, probe_features$START, sep="_")
  nprobes <<- nrow(signal_data) - nrow_missed
  nsamples <<- ncol(signal_data)
  mySampleSheet <- read.csv2("~/Documents/Dati_Lavoro/beckwith-wiedemann/raw/GSE95486/final_samplesheet.csv")
  mySampleSheet$Covariates1 <- stats::rnorm(nrow(mySampleSheet), mean= 567, sd= 1000)
  mySampleSheet$Covariates2 <- stats::rnorm(nrow(mySampleSheet), mean= 67, sd= 100)
  mySampleSheet$Phenotest <-  mySampleSheet[,3]

  mySampleSheet_batch <<- list(mySampleSheet, mySampleSheet, mySampleSheet)
  signal_data_batch <<- list(signal_data, signal_data, signal_data)
}

# TODO
# recover session stored
#
tmp <- normalizePath(tempdir())
tempFolders <<- paste(tmp,"/semseeker/",stringi::stri_rand_strings(50, 7, pattern = "[A-Za-z0-9]"),sep="")


check_execution_context <- function() {
  calls <- sys.calls()
  if (any(sapply(calls, function(x) "test_file" %in% names(x)))) {
    showprogress <<- FALSE
    verbosity <<- 1
    # message("Called from testthat")
  } else {
    showprogress <<- TRUE
    verbosity <<- 4
    # message("Called from source or directly")
  }
}

check_execution_context()

# TODO
# recover session stored
#
tmp <- tempdir()
tempFolders <<- paste(tmp,"/semseeker/",stringi::stri_rand_strings(50, 7, pattern = "[A-Za-z0-9]"),sep="")
