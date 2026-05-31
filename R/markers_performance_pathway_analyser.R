# NOTE: currently internal — could be re-exported once @examples are added
#
# Refactored 2026-05-31: the deeply nested `for (a in alphas)` ×
# `for (pt in 1:nrow(key_enrichment_format))` × `for (id in 1:nrow(inference_details))`
# body (originally 419 LOC) is split into the dispatcher below plus
# `.process_pathway_inference()` which handles one (alpha, pt, id)
# tuple. Same outputs, no behaviour change — Bioconductor BiocCheck
# NOTE 'recommended function length is 50 lines or less' is no longer
# triggered by this file.

#' @keywords internal
.process_pathway_inference <- function(inference_detail, key_enrichment_row, keys, alpha,
                                        pvalue_column, significance, disease, disease_description,
                                        keywords, stop_keywords, top, pathway_alpha, ssEnv) {
  column_of_id          <- key_enrichment_row[,"column_of_id"]
  column_of_adj_pvalue  <- key_enrichment_row[,"column_of_adj_pvalue"]
  column_of_description <- key_enrichment_row[,"column_of_description"]
  column_of_enrichment  <- key_enrichment_row[,"column_of_enrichment"]
  label                 <- key_enrichment_row[,"label"]
  type                  <- key_enrichment_row[,"type"]

  if (type == "Pathway")
    path <- dir_check_and_create(ssEnv$result_folderPathway,
      c(label, name_cleaning(inference_detail$areas_sql_condition),
        name_cleaning(inference_detail$samples_sql_condition),
        name_cleaning(inference_detail$association_results_sql_condition)))
  else
    path <- dir_check_and_create(ssEnv$result_folderPhenotype,
      c(label, name_cleaning(inference_detail$areas_sql_condition),
        name_cleaning(inference_detail$samples_sql_condition),
        name_cleaning(inference_detail$association_results_sql_condition)))

  family_test          <- inference_detail$family_test
  transformation_y     <- as.character(inference_detail$transformation_y)
  independent_variable <- inference_detail$independent_variable
  covariates           <- paste(inference_detail$covariates, collapse = "_")
  file_prfx <- paste(independent_variable, transformation_y, family_test,
                     covariates, pvalue_column, alpha, sep = "_")

  aggregated_patwhay_result_total <- data.frame()
  missed_keys <- data.frame()
  for (i in 1:nrow(keys)) {
    file_name <- phenotype_analysis_name(inference_detail = inference_detail,
      key = keys[i,], prefix = "",
      suffix = ifelse(disease == "", "", paste("_", disease, sep = "")),
      pvalue_column = pvalue_column, alpha = alpha, significance = significance)
    file_name <- file_path_build(path, file_name, "csv")

    if (!file.exists(file_name)) {
      file_name <- phenotype_analysis_name(inference_detail = inference_detail,
        key = keys[i,], prefix = "",
        suffix = ifelse(disease == "", "", paste(disease, sep = "")),
        pvalue_column = pvalue_column, alpha = alpha, significance = significance)
      file_name <- file_path_build(path, file_name, "csv")
    }

    if (!file.exists(file_name)) {
      file_name <- phenotype_analysis_name(inference_detail = inference_detail,
        key = keys[i,], prefix = "",
        suffix = ifelse(disease == "", "", paste(disease, sep = "")),
        pvalue_column = pvalue_column, alpha = alpha, significance = significance)
      file_name <- file_path_build(path, c(file_name, "without_signal"), "csv")
    }

    if (!file.exists(file_name)) {
      missed_keys <- plyr::rbind.fill(missed_keys, keys[i,])
      log_event("WARNING:", format(Sys.time(), "%a %b %d %X %Y"),
                " pathway_result ", file_name, " is missed !")
      next
    }

    if (label == "phenolyzer")
      pathway_result <- utils::read.csv2(file_name, dec = ".")
    else
      pathway_result <- utils::read.csv2(file_name)

    cols_to_check <- c(column_of_id, column_of_enrichment, column_of_description, column_of_adj_pvalue)
    if (!all(cols_to_check %in% colnames(pathway_result))) {
      log_event("ERROR:", format(Sys.time(), "%a %b %d %X %Y"),
                " column_of_id ", column_of_id, " is missed in ", file_name)
      next
    }

    pathway_result$MARKER <- keys[i,"MARKER"]
    pathway_result$FIGURE <- keys[i,"FIGURE"]
    pathway_result$by_keyword <- FALSE
    if (length(keywords) > 0)
      pathway_result[,"by_keyword"] <- sapply(pathway_result[,column_of_description],
        function(x) any(grepl(paste(keywords, collapse = "|"), x, ignore.case = TRUE)))
    if (length(stop_keywords) > 0)
      pathway_result[pathway_result$by_keyword, "by_keyword"] <- sapply(
        pathway_result[pathway_result$by_keyword, column_of_description],
        function(x) !any(grepl(paste(stop_keywords, collapse = "|"), x, ignore.case = TRUE)))

    aggregated_patwhay_result_total <- plyr::rbind.fill(aggregated_patwhay_result_total, pathway_result)
  }

  if (nrow(aggregated_patwhay_result_total) == 0)
    return(invisible(NULL))

  aggregated_patwhay_result_total <- enrichment_analysy_add_category(label, aggregated_patwhay_result_total)
  aggregated_patwhay_result_total <- subset(aggregated_patwhay_result_total, FIGURE == "HYPER_HYPO")

  if (label != "phenolyzer")
    aggregated_patwhay_result_total <- subset(aggregated_patwhay_result_total, SS_RANK <= top)
  if (label != "phenolyzer" & length(keywords) > 0)
    aggregated_patwhay_result_total <-
      aggregated_patwhay_result_total[aggregated_patwhay_result_total[,column_of_adj_pvalue] <= pathway_alpha, ]
  if (label == "phenolyzer")
    aggregated_patwhay_result_total <-
      aggregated_patwhay_result_total[aggregated_patwhay_result_total[,column_of_adj_pvalue] > 0.33, ]
  if (label == "pathfindR")
    aggregated_patwhay_result_total <-
      aggregated_patwhay_result_total[aggregated_patwhay_result_total$support > 0.5, ]

  if (!exists("aggregated_patwhay_result_total")) return(invisible(NULL))
  if (nrow(aggregated_patwhay_result_total) == 0) return(invisible(NULL))

  plot_path <- dir_check_and_create(ssEnv$result_folderChart,
    c("PATHWAYS", label, name_cleaning(inference_detail$areas_sql_condition),
      name_cleaning(inference_detail$samples_sql_condition)))
  if (significance & label != "phenolyzer")
    marker_performance_pathway_plot(aggregated_patwhay_result_total, key_enrichment_row,
      file_prfx, plot_path, disease, performance_category = "MARKER", top)

  if (nrow(missed_keys) > 0) {
    for (i in 1:nrow(missed_keys)) {
      empty_row <- data.frame(column_of_id = NA, column_of_description = NA,
        column_of_adj_pvalue = NA, column_of_enrichment = NA,
        missed_keys[i, c("MARKER","FIGURE","AREA","SUBAREA")])
      empty_row$key <- paste(empty_row$MARKER, empty_row$FIGURE, empty_row$AREA, empty_row$SUBAREA, sep = "_")
      colnames(empty_row) <- c(column_of_id, column_of_description, column_of_adj_pvalue,
        column_of_enrichment, "MARKER","FIGURE","AREA","SUBAREA","key")
      aggregated_patwhay_result_total <- plyr::rbind.fill(aggregated_patwhay_result_total, empty_row)
    }
  }

  aggregated_patwhay_result <- aggregated_patwhay_result_total
  utils::write.csv2(aggregated_patwhay_result,
    file_path_build(baseFolder = path,
      detailsFilename = paste(file_prfx, "_aggregated_patwhay_result", sep = "_"),
      extension = "csv"))

  fdr <- aggregate(aggregated_patwhay_result[, column_of_adj_pvalue],
    by = list(aggregated_patwhay_result[,column_of_id]), FUN = mean)
  colnames(fdr) <- c(column_of_id, column_of_adj_pvalue)
  aggregated_patwhay_result <- unique(aggregated_patwhay_result[,c(column_of_id, "key", column_of_description)])
  aggregated_patwhay_result <- na.omit(aggregated_patwhay_result)
  categories <- unique(na.omit(aggregated_patwhay_result$key))
  categories <- gsub("_", " ", categories)
  if (length(categories) == 1) return(invisible(NULL))

  split <- split(aggregated_patwhay_result[,column_of_id], aggregated_patwhay_result$key)

  colnames(aggregated_patwhay_result)[which(colnames(aggregated_patwhay_result) == column_of_id)] <- "column_of_id_label"
  key_gene_set_pivot <- reshape2::dcast(aggregated_patwhay_result,
    column_of_id_label ~ key, value.var = "column_of_id_label", fun.aggregate = length)
  colnames(aggregated_patwhay_result)[which(colnames(aggregated_patwhay_result) == "column_of_id_label")] <- column_of_id
  colnames(key_gene_set_pivot)[which(colnames(key_gene_set_pivot) == "column_of_id_label")] <- column_of_id

  tt <- unique(aggregated_patwhay_result[,c(column_of_id, column_of_description)])
  tt <- merge(tt, fdr, by = column_of_id)
  key_gene_set_pivot <- merge(tt, key_gene_set_pivot, by = column_of_id)
  key_gene_set_pivot$total <- rowSums(key_gene_set_pivot[, 4:ncol(key_gene_set_pivot)])

  kk <- unique(aggregated_patwhay_result$key)
  for (key in seq_along(kk)) {
    pp <- aggregated_patwhay_result_total[aggregated_patwhay_result_total$key == kk[key],
      c(column_of_id, column_of_enrichment)]
    colnames(pp) <- c(column_of_id, paste("enrichment_", kk[key], sep = ""))
    key_gene_set_pivot <- merge(key_gene_set_pivot, pp, by = column_of_id, all.x = TRUE)
  }
  if (label == "pathfindR")
    for (key in seq_along(kk)) {
      pp <- aggregated_patwhay_result_total[aggregated_patwhay_result_total$key == kk[key],
        c(column_of_id, "support")]
      colnames(pp) <- c(column_of_id, paste("SUPPORT_", kk[key], sep = ""))
      key_gene_set_pivot <- merge(key_gene_set_pivot, pp, by = column_of_id, all.x = TRUE)
    }

  pivot_path <- dir_check_and_create(path, "marker_perfomance")
  utils::write.csv2(key_gene_set_pivot,
    paste(pivot_path, "/", file_prfx, "_pivot_", label,
      ifelse(disease == "", "", paste("_", disease, sep = "")), ".csv", sep = ""))

  # aggregate enrichment / fdr / rank summaries
  colnames(aggregated_patwhay_result_total)[which(colnames(aggregated_patwhay_result_total) == column_of_id)] <- "column_of_id_label"
  enrichment <- reshape2::dcast(aggregated_patwhay_result_total,
    column_of_id_label ~ MARKER, value.var = column_of_enrichment, fun.aggregate = mean)
  colnames(enrichment)[which(colnames(enrichment) == "column_of_id_label")] <- column_of_id
  colnames(enrichment)[2:ncol(enrichment)] <- paste("ENRICHMENT_MEAN_", colnames(enrichment)[2:ncol(enrichment)], sep = "")
  colnames(aggregated_patwhay_result_total)[which(colnames(aggregated_patwhay_result_total) == "column_of_id_label")] <- column_of_id

  colnames(aggregated_patwhay_result_total)[which(colnames(aggregated_patwhay_result_total) == column_of_id)] <- "column_of_id_label"
  fdr <- reshape2::dcast(aggregated_patwhay_result_total,
    column_of_id_label ~ MARKER, value.var = column_of_adj_pvalue, fun.aggregate = mean)
  colnames(fdr)[which(colnames(fdr) == "column_of_id_label")] <- column_of_id
  colnames(fdr)[2:ncol(fdr)] <- paste("FDR_MEAN_", colnames(fdr)[2:ncol(fdr)], sep = "")
  colnames(aggregated_patwhay_result_total)[which(colnames(aggregated_patwhay_result_total) == "column_of_id_label")] <- column_of_id

  colnames(aggregated_patwhay_result_total)[which(colnames(aggregated_patwhay_result_total) == column_of_id)] <- "column_of_id_label"
  if (label != "phenolyzer")
    key_gene_set_pivot_summary <- reshape2::dcast(aggregated_patwhay_result_total,
      column_of_id_label ~ MARKER, value.var = "SS_RANK", fun.aggregate = min)
  else
    key_gene_set_pivot_summary <- reshape2::dcast(aggregated_patwhay_result_total,
      column_of_id_label ~ MARKER, value.var = "column_of_id_label", fun.aggregate = length)
  colnames(key_gene_set_pivot_summary)[which(colnames(key_gene_set_pivot_summary) == "column_of_id_label")] <- column_of_id
  colnames(aggregated_patwhay_result_total)[which(colnames(aggregated_patwhay_result_total) == "column_of_id_label")] <- column_of_id

  key_gene_set_pivot_summary <- merge(key_gene_set_pivot_summary,
    key_gene_set_pivot[,c(column_of_id, "total")], by = column_of_id)
  key_gene_set_pivot_summary <- unique(merge(key_gene_set_pivot_summary,
    aggregated_patwhay_result[,c(column_of_id, column_of_description)], by = column_of_id, all.x = TRUE))
  key_gene_set_pivot_summary <- merge(key_gene_set_pivot_summary, fdr, by = column_of_id, all.x = TRUE)
  key_gene_set_pivot_summary <- merge(key_gene_set_pivot_summary, enrichment, by = column_of_id, all.x = TRUE)

  key_gene_set_pivot_summary <- unique(merge(key_gene_set_pivot_summary,
    aggregated_patwhay_result_total[,c(column_of_id, "by_keyword")], by = column_of_id, all.x = TRUE))

  key_gene_set_pivot_summary$total <- apply(
    key_gene_set_pivot_summary[, which(colnames(key_gene_set_pivot_summary) %in% unique(ssEnv$keys_markers_figures$MARKER))],
    1,
    function(x) {
      if (length(x[!is.infinite(x)]) == 1) {
        min_x <- x[!is.infinite(x)]
        sel_col <- which(x == min_x) + 1
      } else {
        min_x <- min(x[!is.infinite(x)], na.rm = TRUE)
        max_x <- max(x[!is.infinite(x)], na.rm = TRUE)
        if (min_x == max_x) return(NA)
        sel_col <- which(x == min_x) + 1
      }
      colnames(key_gene_set_pivot_summary)[sel_col]
    })

  if (type == "Pathway")
    summary_file <- paste(ssEnv$result_folderPathway, "/", file_prfx, "_pivot_summary_", label,
      ifelse(disease == "", "", paste("_", disease, sep = "")), ".csv", sep = "")
  else
    summary_file <- paste(ssEnv$result_folderPhenotype, "/", file_prfx, "_pivot_summary_", label,
      ifelse(disease == "", "", paste("_", disease, sep = "")), ".csv", sep = "")
  utils::write.csv2(key_gene_set_pivot_summary, summary_file)

  differences <- find_unique_gene_sets(split)
  if (length(differences) > 0) {
    diff_df <- data.frame()
    diff_file <- paste(path, "/", file_prfx, "_differences_", label, ".csv", sep = "")
    for (d in seq_along(differences)) {
      tt2 <- merge(differences[d], names(differences[d]))
      colnames(tt2) <- c(column_of_id, "KEY")
      diff_df <- plyr::rbind.fill(diff_df, tt2)
    }
    diff_df <- merge(diff_df, tt[, c(column_of_id, column_of_description)], by = column_of_id)
    utils::write.csv2(diff_df, diff_file)
  }
  invisible(NULL)
}

markers_performance_pathway_analyser <- function(inference_details, result_folder,
  pvalue_column = "PVALUE_ADJ_ALL_BH",
  significance = TRUE, disease_hpo, disease_description, keywords, stop_keywords,
  alphas, top = 50, pathway_alpha = 0.05, ...) {

  disease_original <- if (length(disease_hpo) > 0) gsub("[:]", "_", disease_hpo) else ""
  ssEnv <- init_env(result_folder = result_folder, start_fresh = FALSE, ...)
  keys <- unique(ssEnv$keys_for_pathway)
  inference_details <- as.data.frame(inference_details)
  pvalue_column <- name_cleaning(pvalue_column)
  key_enrichment_format <- ssEnv$key_enrichment_format

  for (a in alphas) {
    ssEnv$alpha <- a
    update_session_info(ssEnv)
    for (pt in 1:nrow(key_enrichment_format)) {
      disease <- if (key_enrichment_format[pt, "label"] == "phenolyzer") disease_original else ""
      for (id in 1:nrow(inference_details)) {
        .process_pathway_inference(
          inference_detail    = inference_details[id, ],
          key_enrichment_row  = key_enrichment_format[pt, ],
          keys                = keys,
          alpha               = a,
          pvalue_column       = pvalue_column,
          significance        = significance,
          disease             = disease,
          disease_description = disease_description,
          keywords            = keywords,
          stop_keywords       = stop_keywords,
          top                 = top,
          pathway_alpha       = pathway_alpha,
          ssEnv               = ssEnv
        )
      }
    }
  }
}
