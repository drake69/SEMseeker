# test_that("cluster", {
#
#   # source("./tests/testthat/setup.R")
#   tempFolder <- "~/Documents/tmp_semseeker/test_1"
#   parallel_strategy <- "multicore"
#   # 10.161.5.13
#   # 192.168.37.207
#   ssEnv <- SEMseeker:::init_env(tempFolder, parallel_strategy = parallel_strategy,
#     areas = c("GENE"), subareas= c("BODY"), showprogress=FALSE, start_fresh=TRUE,
#     cluster_workers=c("localhost"))
#
#   library(future)
#
#   funcA <- function()
#   {
#     SEMseeker:::get_session_info()
#   }
#
#   funcB <- function()
#   {
#     funcA()
#   }
#
#   funcC <- function()
#   {
#     funcB()
#   }
#
#   res <- foreach::foreach(i = 1:10, .combine = rbind, .export =c("funcA","funcB","funcC","SEMseeker:::get_session_info")) %dorng% {
#     funcA()
#   }
#
#   res <- foreach::foreach(i = 1:10, .combine = rbind, .export =c("funcA","funcB","funcC")) %dorng% {
#     funcB()
#   }
#
#   res <- foreach::foreach(i = 1:10, .combine = rbind, .export =c("funcA","funcB","funcC")) %dorng% {
#     funcC()
#   }
#
#   testthat::expect_true(length(res)>10)
#
#   gg1 %<-% {
#     # p <- ggplot(data = d, aes(x = carat, y = price)) +
#     #   geom_point(aes(text = paste("Clarity:", clarity)), size = 4) +
#     #   geom_smooth(aes(colour = cut, fill = cut)) + facet_wrap(~ cut)
#     # ggplotly(p)
#     SEMseeker:::get_session_info()
#   }
#   testthat::expect_true(gg1$cluster_workers=="192.168.37.207")
#   testthat::expect_true(gg1$result_folderData==paste(tempFolder,"/Data",sep=""))
# })
