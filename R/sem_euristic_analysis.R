sem_euristic_analysis_phenolizer <- function(result_folder, rule, cutoff)
{

  # cutoff mean
  # cases mutated control not mutated
  ssEnv <- core_init_env(result_folder)

  geneAnnotatedFile <- utils::read.csv2(file.path(ssEnv$result_folder_data , "GENE_ANNOTATED.csv"))
  # geneAnnotatedFile <-subset(geneAnnotatedFile,geneAnnotatedFile$POPULATION != "Reference")

  markers <- c("MUTATIONS","LESIONS")
  for ( m in seq_along(markers))
  {
    marker <- markers[m]
    # marker <- "LESIONS"
    geneAnnotated <-subset(geneAnnotatedFile,geneAnnotatedFile$MARKER == marker)
    tempDataFrame <- as.data.frame(reshape2::dcast(data = geneAnnotated, SAMPLEID  ~ GENE, value.var = "freq", sum))
    tempDataFrame <- as.data.frame(t(tempDataFrame))
    tempDataFrame <- apply(tempDataFrame[],1, as.numeric)
    cutOff <- mean(tempDataFrame, na.rm = TRUE)
    geneAnnotated <- geneAnnotated[geneAnnotated$freq>cutOff, ]

    geneAnnotated <-unique(geneAnnotated[, c("GENE", "POPULATION", "MARKER", "SAMPLEID")])
    geneAnnotated$freq <- 1

    tempDataFrame <- as.data.frame(reshape2::dcast(data = geneAnnotated, POPULATION  ~ GENE, value.var = "freq", sum))

    geneMutated <- as.data.frame(t(tempDataFrame))
    names(geneMutated) <- geneMutated[1, ]
    geneMutated <- geneMutated[-1, ]
    geneMutated$gene <- rownames(geneMutated)
    # head(geneMutated[order(geneMutated$Case, decreasing = TRUE),])

    # caseCutOff <- mean(as.numeric(geneMutated$Case))
    # geneMutated <- geneMutated[(as.numeric(geneMutated$Case) > caseCutOff & as.numeric(geneMutated$Control) == 0), ]
    geneMutated <- geneMutated[(as.numeric(geneMutated$Case) > 0 & as.numeric(geneMutated$Control) == 0), ]

    geneMutated$gene <- rownames(geneMutated)
    # utils::write.table(
    #   row.names(geneMutated),
    #   file.path(ssEnv$result_folder_data, paste0("/Euristic/Phenolyzer/", marker,"_genes_only_case.csv",sep="")),
    #   row.names = FALSE,
    #   col.names = FALSE ,
    #   quote = FALSE,
    #   sep = "\t"
    # )
    # utils::write.table(
    #   geneMutated[order(geneMutated$Case, decreasing = TRUE),] ,
    #   file.path(ssEnv$result_folder_data, paste0("/Euristic/Phenolyzer/", marker, "_genes_only_case_all_details.csv", sep="")),
    #   row.names = FALSE,
    #   col.names = TRUE ,
    #   quote = FALSE,
    #   sep = ";"
    # )
    #
  }
}
