enrichment_analysy_add_category <- function(source, data)
{

  if(nrow(data)==0)
    return(data)

  ssEnv <- get_session_info()
  #TO DO: normalize naming of differente data sources
  if(source=="phenolyzer")
    return(data)

  data$SS_CATEGORY <- NA
  if(any(source %in% c("WebGestalt","Phenolyzer_WebGestalt")))
  {
    #
    data$type <- toupper(data$type)
    data[data$type=='BP',"SS_CATEGORY"] <- "GO-BP"
    data[data$type=='CC',"SS_CATEGORY"] <- "GO-CC"
    data[data$type=='MF',"SS_CATEGORY"] <- "GO-MF"
  }
  if(source=="ctdR")
  {
    data$source <- toupper(data$source)
    data[,"SS_CATEGORY"] <- "CHEMICAL"
  }

  if(source=="pathfindR")
  {
    #
    data$source <- toupper(data$source)
    data[data$source=='GO-BP',"SS_CATEGORY"] <- "GO-BP"
    data[data$source=='GO-MF',"SS_CATEGORY"] <- "GO-MF"
    data[data$source=='KEGG',"SS_CATEGORY"] <- "KEGG-NTA"
    data[data$source=='REACTOME',"SS_CATEGORY"] <- "REACTOME-NTA"
    data$PVALUE_ADJ_ALL_FDR <-stats::p.adjust(data[,"highest_p"], method = "fdr")
  }

  if(any(source %in% c("STRINGdb","Phenolyzer_STRINGdb")))
  {
    data$category <- toupper(data$category)
    data$SS_CATEGORY <- data$category
    data[data$category=='PROCESS',"SS_CATEGORY"] <- "GO-BP"
    data[data$category=='FUNCTION',"SS_CATEGORY"] <- "GO-MF"
    data[data$category=='COMPONENT',"SS_CATEGORY"] <- "GO-CC"
    data[data$category=='KEGG',"SS_CATEGORY"] <- "KEGG-ORA"
  }

  with_category <- data[!is.na(data$SS_CATEGORY),]
  categories <- na.omit(unique(with_category$SS_CATEGORY))

  key_enrichment_format <- ssEnv$key_enrichment_format


  # Normalization functions
  # normalize_minimize <- function(x) (max(x) - x) / (max(x) - min(x))
  # normalize_maximize <- function(x) (x - min(x)) / (max(x) - min(x))

  # rank per SS_CATEGORY
  data <- data[is.na(data$SS_CATEGORY),]
  for (c in categories)
  {
    if(!any(c %in% with_category$SS_CATEGORY))
      next
    to_rank <- with_category[with_category$SS_CATEGORY==c,]
    if(nrow(to_rank)==0)
      next
    # fdr
    column_to_rank <- key_enrichment_format[key_enrichment_format$label==source,"column_of_adj_pvalue"]
    if (!any(column_to_rank %in% colnames(to_rank))) next
    if(max(to_rank[,column_to_rank], na.rm=TRUE) == min(to_rank[,column_to_rank], na.rm=TRUE))
      to_rank$SS_RANK_FDR <- 1
    else
    {
      to_rank$SS_RANK_FDR <- normalize_minimize(to_rank[,column_to_rank])
      to_rank$SS_RANK_FDR <- rank(-to_rank$SS_RANK_FDR, ties.method = "min")
    }
    to_rank$SIGNIFICATIVE <- (to_rank[,column_to_rank] < as.numeric(ssEnv$alpha))
    # get max rank of column_of_adj_pvalue when is lower than alpha
    max_rank_pvalue <- max(to_rank$SS_RANK_FDR[to_rank[,column_to_rank] < as.numeric(ssEnv$alpha)], na.rm = TRUE)
    # replace for all values of SS_RANK_FDR that are higher than alpha
    to_rank$SS_RANK_FDR[to_rank[,column_to_rank] >= as.numeric(ssEnv$alpha)] <- max_rank_pvalue + 1

    column_to_rank <- key_enrichment_format[key_enrichment_format$label==source,"column_of_enrichment"]
    if (!any(column_to_rank %in% colnames(to_rank))) next
    if(max(to_rank[,column_to_rank], na.rm=TRUE) == min(to_rank[,column_to_rank], na.rm=TRUE))
      to_rank$SS_RANK_ENRICHMENT <- 1
    else
    {
      to_rank$SS_RANK_ENRICHMENT <- normalize_maximize(to_rank[,column_to_rank])
      to_rank$SS_RANK_ENRICHMENT <- rank(-to_rank$SS_RANK_ENRICHMENT, ties.method = "min")
    }

    # do a total score
    to_rank$SS_RANK <- to_rank$SS_RANK_FDR + to_rank$SS_RANK_ENRICHMENT

    to_rank$SS_RANK <- rank(-to_rank$SS_RANK, ties.method = "min")
    data <- plyr::rbind.fill(data,to_rank)
  }

  data$SS_RANK <- round(data$SS_RANK,2)

  return(data)


}
