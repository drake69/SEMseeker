# TODO:
# test_that("quantreg_summary", {
#
#   
#
#   nitem <- 1e4
#   nsamples <- 4
#   tempDataFrame <- rnorm(nitem*nsamples,mean = 0.5, sd = 0.7)
#   tempDataFrame <- as.data.frame(matrix(tempDataFrame,nitem,nsamples))
#
#   Sample_ID <- stri_rand_strings(nsamples, 14, pattern = "[A-Za-z]")
#   colnames(tempDataFrame) <- Sample_ID
#
#   independent_variable <- paste(Sample_ID[2], sep="")
#   sig.formula <- stats::as.formula(paste0(Sample_ID[1] ,"~",Sample_ID[2],"+",Sample_ID[3], "+", Sample_ID[4], sep=""))
#   lqm_control <- list(loop_tol_ll = 1e-5, loop_max_iter = 5000, verbose = F )
#   n_permutations_test <- 5000
#   tau <- 0.5
#
#   model.x <-  suppressMessages(lqmm::lqm(sig.formula, tau=tau,  data=as.data.frame(tempDataFrame) , na.action = stats::na.omit, control = lqm_control))
#   model.x.boot <- suppressMessages(lqmm::boot(model.x, R = n_permutations_test))
#   signal_full <- suppressMessages(summary(model.x.boot)[independent_variable,"Value"])
#   tt <- as.data.frame((as.matrix.data.frame(model.x.boot)))
#   colnames(tt) <- colnames(model.x.boot)
#   boot_vector <- stats::na.omit(tt[,independent_variable])
#   boot.bca <- quantreg_summary(boot_vector, signal_full, as.data.frame(tempDataFrame), sig.formula, tau, independent_variable, lqm_control = lqm_control)
#
#   testthat::expect_true(sum(is.null(boot.bca))==0)
#   SEMseeker:::close_env()
# })
#
