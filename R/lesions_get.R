# grouping_column <- "GENE"
# mutationAnnotatedSortedLocal <- utils::read.csv2("/home/lcorsaro/Documents/SEMSEEKER_TEST_BWS/mutations_annotated_sorted_gene.csv")

#' @importFrom dplyr %>%
#' @importFrom rlang .data
lesions_get <- function(grouping_column, mutation_annotated_sorted)
{

  ssEnv <- get_session_info()

  if( is.null(mutation_annotated_sorted))
    return (mutation_annotated_sorted)

  if(nrow(mutation_annotated_sorted) == 0)
    return (mutation_annotated_sorted)

    # #
  mutationAnnotatedSortedLocal <- mutation_annotated_sorted
  summed <- stats::aggregate(mutationAnnotatedSortedLocal$MUTATIONS, by = list(mutationAnnotatedSortedLocal[,grouping_column]), FUN = sum)
  colnames(summed) <- c(grouping_column,"MUTATIONS_COUNT")
  counted <- stats::aggregate(mutationAnnotatedSortedLocal$MUTATIONS, by = list(mutationAnnotatedSortedLocal[,grouping_column]), FUN = length)
  colnames(counted) <- c(grouping_column,"PROBES_COUNT")
  mutationAnnotatedSortedLocal <- merge(mutationAnnotatedSortedLocal,summed, by = grouping_column)
  mutationAnnotatedSortedLocal <- merge(mutationAnnotatedSortedLocal,counted, by = grouping_column)
  rm(counted)
  rm(summed)

  #function to calculate rolled sum, returns a vector
  enrichement_calculator<-function(x,lags){
    missedWindowLength <- ((lags - 1) / 2)
    y <- c(rep(0,missedWindowLength))
    tmp <- c(y,x,y)
    tmp <- zoo::rollsum( tmp, lags, align  =  "center", fill = 0)
    tmp <- tmp[(missedWindowLength+1):(missedWindowLength+length(x))]
    tmp <- as.numeric(tmp)
    return(tmp)
  }

  # #
  # calculate enrichment for each window
  mutationAnnotatedSortedLocal <- mutationAnnotatedSortedLocal %>% dplyr::group_by(eval(parse(text = grouping_column))) %>%
    dplyr::mutate(ENRICHMENT = stats::ave( .data$MUTATIONS, eval(parse(text = grouping_column)),
                                           FUN  = function(x) enrichement_calculator(x, ssEnv$sliding_window_size))) %>% dplyr::ungroup ()

  basepair_calculator<-function(x,lags){
    tmp_min <- -zoo::rollmax( -x, lags, align="center", fill  <-  0)
    tmp_max <- zoo::rollmax( x, lags, align="center", fill  <-  0)
    tmp <-  tmp_max - tmp_min
    tmp <- as.numeric(tmp)
    return(tmp)
  }
  #calculate the base pair count for each window
  mutationAnnotatedSortedLocal <- mutationAnnotatedSortedLocal %>% dplyr::group_by(eval(parse(text = grouping_column))) %>%
    dplyr::mutate(BASEPAIR_COUNT = stats::ave( .data$START, eval(parse(text = grouping_column)),
                                              FUN  =  function(x) basepair_calculator(x, ssEnv$sliding_window_size))) %>% dplyr::ungroup ()

  mutationAnnotatedSortedLocal$ENRICHMENT[ is.na(mutationAnnotatedSortedLocal$ENRICHMENT)] <- 0

  # Background mutation rate per grouping unit (gene / chr / …).
  # Using Binomial rather than Hypergeometric because the sliding window
  # advances one probe at a time, so each probe participates in up to
  # sliding_window_size consecutive windows — i.e. sampling WITH replacement.
  # The hypergeometric distribution assumes sampling WITHOUT replacement and
  # therefore produces inflated (too-small) p-values here.
  #
  # H0: ENRICHMENT ~ Binomial(sliding_window_size, p0)
  # p0 = MUTATIONS_COUNT / PROBES_COUNT  (empirical background rate)
  # p-value = P(X >= ENRICHMENT) = pbinom(ENRICHMENT - 1, size, p0, lower.tail = FALSE)
  p0 <- mutationAnnotatedSortedLocal$MUTATIONS_COUNT / mutationAnnotatedSortedLocal$PROBES_COUNT

  lesionpValue <- stats::pbinom(
    mutationAnnotatedSortedLocal$ENRICHMENT - 1L,
    size       = as.integer(ssEnv$sliding_window_size),
    prob       = p0,
    lower.tail = FALSE
  )

  lesionpValue[is.nan(lesionpValue)] <- 1
  lesionpValue[is.na(lesionpValue)] <- 1

  tt <- data.frame(mutationAnnotatedSortedLocal,lesionpValue)

  ## correction by Bonferroni
  lesionWeighted <- (tt$lesionpValue ) < ( as.numeric(ssEnv$bonferroni_threshold)/( length(tt$PROBES_COUNT) * log10(tt$BASEPAIR_COUNT) ))
  table(lesionWeighted)
  rm(tt)

  lesionWeighted <- data.frame(as.data.frame(mutationAnnotatedSortedLocal), "LESIONS"  =  lesionWeighted)

  lesionWeighted <- sort_by_chr_and_start(lesionWeighted)
  lesionWeighted <- subset(lesionWeighted, lesionWeighted$LESIONS == TRUE)[, c("CHR", "START", "END")]

  if (dim(lesionWeighted)[1] > dim(mutationAnnotatedSortedLocal)[1]) {

  }

 log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Got lesions for sample !")
  return(lesionWeighted)

}
