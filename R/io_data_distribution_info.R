# data_distribution_info <- function(family_test, tempDataFrame, burdenValue, independent_variable)
# {
#   family_test_missed <- TRUE
#   if(family_test=="binomial" | family_test=="wilcoxon" | family_test=="jsd" | family_test=="t.test"
#     | grepl("mean-permutation",family_test) | family_test=="chisq.test" | family_test=="fisher.test" | family_test=="kruskal.test")
#   {
#     tempDataFrame[,independent_variable] <- as.factor(tempDataFrame[,independent_variable])
#     bartlett_pvalue <- stats::bartlett.test( stats::as.formula(paste0(burdenValue,"~", independent_variable, sep="")),
#       data= as.data.frame(tempDataFrame) )
#     family_test_missed <- FALSE
#   }
#
#   if(family_test=="gaussian" | family_test=="spearman" | family_test=="kendall" | family_test=="pearson"
#     | grepl("quantreg",family_test) | grepl("poisson",family_test) | grepl("spearman-permutation",family_test)
#     | grepl("quantile-permutation",family_test) | grepl("polynomial",family_test) | grepl("exp",family_test) | grepl("log",family_test))
#   {
#     localDataFrame <- data.frame("depVar"=tempDataFrame[, burdenValue],"indepVar"=1 )
#     localDataFrame <- rbind( localDataFrame,  data.frame("depVar"=tempDataFrame[, independent_variable],"indepVar"=2 ))
#     bartlett_pvalue <- stats::bartlett.test( stats::as.formula("depVar ~ indepVar"), data= localDataFrame )
#     family_test_missed <- FALSE
#   }
#
#   if(family_test_missed)
#     stop("data_distribution_info: familyt test is missed")
#
#   return(bartlett_pvalue$p.value)
# }
