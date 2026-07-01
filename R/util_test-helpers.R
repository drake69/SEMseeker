# if(getRversion() >= "2.15.1")  utils::globalVariables(c( "ssEnv$result_folderData","chartFolder","ssEnv$session_folder","computationCluster",
#                                                          "ssEnv$result_folderChart","ssEnv$result_folderInference","ssEnv$result_folderEuristic",
#                                                          "ssEnv$keys_markers","ssEnv$keys_figures","ssEnv$keys_sample_groups","ssEnv$keys",
#                                                          "probe_features"), add = FALSE)

# my_env <- new.env(parent = emptyenv())
# functionToExport <<- ssEnv$result_folderData <<- chartFolder <<- ssEnv$session_folder <<- computationCluster <<- NULL
# ssEnv$result_folderChart <<- ssEnv$result_folderInference <<- ssEnv$result_folderEuristic <<-NULL
# ssEnv$keys_markers <<- ssEnv$keys_figures <<- ssEnv$keys_sample_groups <<- ssEnv$keys <<- probe_features <<- foreachIndex <<- NULL

# utils::globalVariables( names= c( "ssEnv$result_folderData","chartFolder","ssEnv$session_folder","computationCluster",
#                           "ssEnv$result_folderChart","ssEnv$result_folderInference","ssEnv$result_folderEuristic",
#                           "ssEnv$keys_markers","ssEnv$keys_figures","ssEnv$keys_sample_groups","ssEnv$keys",
#                           "probe_features","foreachIndex"), add = FALSE, package = "semseeker")

 # ssEnv <- new.env(parent = emptyenv())
# ssEnv$result_folderData <- NULL
# chartFolder  <- NULL
# ssEnv$session_folder   <- NULL
# computationCluster <- NULL
# ssEnv$result_folderChart  <- NULL
# ssEnv$result_folderInference  <- NULL
# ssEnv$result_folderEuristic   <- NULL
# ssEnv$keys_markers  <- NULL
# ssEnv$keys_figures  <- NULL
# ssEnv$keys_sample_groups   <- NULL
# ssEnv$keys   <- NULL
# probe_features  <- NULL
# foreachIndex  <- NULL
#


# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("BiocGenerics")
# BiocManager::install("XVector")
# BiocManager::install("S4Vectors")
# BiocManager::install("Biostrings")
# BiocManager::install("KEGGREST")
# BiocManager::install("GEOquery")
# BiocManager::install("limma")
# BiocManager::install("zlibbioc")
# BiocManager::install("Biobase")
# BiocManager::install("igraph")
# BiocManager::install("RBGL")
#
# to develope in windows is needed RTools
# https://cran.rstudio.com/bin/windows/Rtools/rtools42/rtools.html
# choco install rtools
#
# BiocManager::install("GenomeInfoDbData")
# or installed by tar
