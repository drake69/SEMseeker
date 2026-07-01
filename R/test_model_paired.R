test_model_paired <- function (family_test, tempDataFrame, sig.formula,burdenValue,independent_variable , transformation_y, plot , samples_sql_condition="", key)
{
  ssEnv <- core_get_session_info()
  res <- data.frame(pvalue=NA)

  tempDataFrameOriginal <- tempDataFrame
  to_split <- family_test
  family_test <- strsplit(to_split, "@")[[1]][1]
  pairing_variable <- strsplit(to_split, "@")[[1]][2]

  # create a tempDataFrame using the pairing_variable to create the pairs and the independent_variable for the pairs
  pairing_category <- unique(tempDataFrame[,independent_variable])
  if(length(pairing_category)!=2)
  {
    core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " I skip this test because the pairing variable has more than 2 levels." )
    return(res)
  }
  first_category <- as.character(pairing_category[1])
  second_category <- as.character(pairing_category[2])

  # create the pairs
  first_cat_values <- tempDataFrame[tempDataFrame[,independent_variable]==first_category, c(pairing_variable, burdenValue)]
  colnames(first_cat_values) <- c(pairing_variable, first_category)

  second_cat_values <- tempDataFrame[tempDataFrame[,independent_variable]==second_category, c(pairing_variable, burdenValue)]
  colnames(second_cat_values) <- c(pairing_variable, second_category)

  tempDataFrame <- merge(first_cat_values, second_cat_values, by=pairing_variable, all=TRUE)
  # drop not complete pairs
  tempDataFrame <- tempDataFrame[complete.cases(tempDataFrame),]

  # check ifthere are enough pairs
  if(nrow(tempDataFrame)<2)
  {
    core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " I skip this test because there are not enough pairs." )
    return(res)
  }

  tempDataFrameOriginal <- tempDataFrameOriginal[tempDataFrameOriginal[,pairing_variable] %in% tempDataFrame[,pairing_variable],]

  sig.formula <- as.formula(paste(first_category, "~", second_category))

  if(family_test=="wilcoxon.paired")
  {
    if (plot)
      plot_box_plot(tempDataFrameOriginal, independent_variable,burdenValue, transformation_y, family_test, samples_sql_condition, key)

    result_w  <- suppressWarnings(stats::wilcox.test(tempDataFrame[,first_category], tempDataFrame[,second_category], exact=TRUE, paired = TRUE))
    res$pvalue <- result_w$p.value

    #
    res[1,"Wilcox_Value"] <- result_w$statistic

    res$r_model <- "stats_wilcox.test.paired"
    # dep_var <- strsplit(gsub("\ ","",as.character(sig.formula)),"~")
    # SPLIT <- split(tempDataFrame[,dep_var[[2]]], tempDataFrame[,dep_var[[3]]])
    #
    # ## calculate the JSD
    # sample1 <- SPLIT[[1]]
    # sample2 <- SPLIT[[2]]
    #
    # # Combine the unique elements from both samples to create a common event space
    # common_events <- unique(c(sample1, sample2))
    #
    # # Create adjusted frequency tables for both samples
    # frequency_table1_adjusted <- tabulate(match(sample1, common_events), nbins = length(common_events))
    # frequency_table2_adjusted <- tabulate(match(sample2, common_events), nbins = length(common_events))
    #
    # # Convert adjusted frequencies to probabilities
    # probability_distribution1_adjusted <- frequency_table1_adjusted / sum(frequency_table1_adjusted)
    # probability_distribution2_adjusted <- frequency_table2_adjusted / sum(frequency_table2_adjusted)
    #
    # # Calculate the Jensen-Shannon distance
    # res$jsd <- suppressMessages(suppressWarnings(philentropy::JSD(rbind(probability_distribution1_adjusted, probability_distribution2_adjusted))))
    #
    # # Calculate effect size
    # es_res <- effsize::VD.A(SPLIT[[1]], SPLIT[[2]])
    # res$effect_size_estimate <- es_res$estimate
    # res$effect_size_magnitude <- es_res$magnitude
    # # res$cohen_d <- effsize::cohen.d(SPLIT[[1]], SPLIT[[2]], pooled=TRUE, paired=FALSE, na.rm=TRUE)
    #
    # # Calculate rank-biserial correlation as effect size
    # res$RANK_BISERIAL_CORRELATION <- result_w$statistic / (length(SPLIT[[1]]) * length(SPLIT[[2]]))
    # # Calculate power
    # power_result = pwr::pwr.t2n.test(d = res$effect_size_estimate, n1 = length(SPLIT[[1]]), n2=length(SPLIT[[2]]), sig.level = as.numeric(ssEnv$alpha), power = NULL)
    # res$power <- power_result$power
  }


  if(family_test=="t.test.paired")
  {
    if (plot)
      plot_box_plot(tempDataFrameOriginal, independent_variable,burdenValue, transformation_y, family_test, samples_sql_condition, key)

    result_w  <-stats::t.test(formula= sig.formula, data = as.data.frame(tempDataFrame), paired = TRUE)
    res$pvalue <- result_w$p.value
    res$r_model <- "stats_t.test"
    dep_var <- strsplit(gsub("\ ","",as.character(sig.formula)),"~")
    SPLIT <- split(tempDataFrame[,dep_var[[2]]], tempDataFrame[,dep_var[[3]]])
    res$statistic_parameter <- mean(SPLIT[[1]]) - mean(SPLIT[[2]])
    power_result <- pwr::pwr.t2n.test(d = res$statistic_parameter, n1 = length(SPLIT[[1]]), n2=length(SPLIT[[2]]), sig.level = as.numeric(ssEnv$alpha), power = NULL)
    res$power <- power_result$power
  }

  return(res)

}
