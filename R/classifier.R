
#' Classify Consensus Transcriptomic Subtypes of Sepsis
#'
#' Classify patients into molecular subtypes using gene expression data from two batches (e.g., core dataset and external cohort).
#'
#' @param exp_core_g Optional. Expression matrix of core samples (genes x samples). Row names are Ensembl gene IDs. Defaults to internal training data.
#' @param core_samples Optional. Data frame of sample labels with a column `CTS` for subtype classification. Row names match column names in exp_core_g. Defaults to internal labels.
#' @param new_expr_data A gene expression matrix (genes x samples) for the new dataset to classify. If a data frame, the first column should be gene identifiers with the column name "gene". Row names must be Ensembl gene IDs matching those in `exp_core_g`.
#' @param heatmap_file Optional file path to save the heatmap PDF. If NULL, the heatmap is printed to screen.
#' @param silhouette_file Optional file path to save the silhouette plot PDF. If NULL, the silhouette plot is returned as an object only.
#' @return A list with predicted subtypes, corrected expression matrices, and silhouette metrics.
#' @export
run_subtype_classifier <- function(new_expr_data, exp_core_g = NULL, core_samples = NULL, heatmap_file = NULL, silhouette_file = NULL) {
  library(sva)
  library(mixOmics)
  library(randomForest)
  library(pheatmap)
  library(RColorBrewer)
  library(cluster)
  library(dplyr)

  if (is.null(exp_core_g) || is.null(core_samples)) {
    data("exp_core_g", package = "ConsensusTranscriptomicSubtype", envir = environment())
    data("core_samples", package = "ConsensusTranscriptomicSubtype", envir = environment())
  }

  g <- data.frame(gene = rownames(exp_core_g))
  core_samples <- core_samples[colnames(exp_core_g), , drop = FALSE]
  cts <- as.factor(core_samples$CTS)

  if ("gene" %in% colnames(new_expr_data)) {
    rownames(new_expr_data) <- new_expr_data$gene
    new_expr_data <- new_expr_data[, -1]
  }

  new_expr_data <- as.matrix(new_expr_data)
  new_expr_data_gm <- new_expr_data[rownames(new_expr_data) %in% g$gene, ]

  exp_core_df <- as.data.frame(exp_core_g)
  exp_core_df$gene <- rownames(exp_core_g)
  new_expr_df <- as.data.frame(new_expr_data_gm)
  new_expr_df$gene <- rownames(new_expr_data_gm)

  batch_df <- rbind(
    data.frame(id = colnames(exp_core_g), batch = "A"),
    data.frame(id = colnames(new_expr_data_gm), batch = "B")
  )

  exp_all <- merge(exp_core_df, new_expr_df, by = "gene")
  rownames(exp_all) <- exp_all$gene
  exp_all_m <- data.matrix(exp_all[, -1])

  exp_all_combat <- ComBat(dat = exp_all_m, batch = batch_df$batch, par.prior = TRUE, prior.plots = FALSE)
  exp_core_combat <- exp_all_combat[, colnames(exp_core_g)]
  new_expr_combat <- exp_all_combat[, colnames(new_expr_data_gm)]

  rf <- randomForest(cts ~ ., data = t(exp_core_combat), proximity = TRUE)
  pred_new <- predict(rf, t(new_expr_combat))
  pred_df <- data.frame(Sample = colnames(new_expr_combat), CTS = pred_new)
  pred_df <- pred_df[!pred_df$Sample %in% c("gene"), ]

  ordered_pred_df <- pred_df %>% arrange(CTS)
  new_expr_ordered <- new_expr_combat[, ordered_pred_df$Sample]

  ann_colors <- list(CTS = c("1" = "royalblue", "2" = "#B2DF8A", "3" = "orange"))

  if (!is.null(heatmap_file)) {
    pdf(heatmap_file, height = 5, width = 8)
    on.exit(dev.off(), add = TRUE)
  }

  pheatmap(
    new_expr_ordered,
    color = colorRampPalette(rev(brewer.pal(11, "PuOr")))(50),
    clustering_distance_rows = "correlation",
    cluster_rows = TRUE, cluster_cols = FALSE,
    show_rownames = TRUE, show_colnames = FALSE,
    annotation_col = ordered_pred_df["CTS", drop = FALSE],
    annotation_colors = ann_colors,
    scale = "row", treeheight_row = 0, border_color = NA
  )

  rf_new <- randomForest(pred_new ~ ., data = t(new_expr_combat), proximity = TRUE)
  sil <- silhouette(as.numeric(pred_new), as.dist(1 - rf_new$proximity))

  if (!is.null(silhouette_file)) {
    pdf(silhouette_file, height = 6, width = 6)
    plot(sil, col = c("royalblue", "#B2DF8A", "orange"))
    dev.off()
  }

  return(list(
    predictions = pred_df,
    expression_corrected = list(core = exp_core_combat, new = new_expr_combat),
    silhouette = sil,
    rf_model = rf
  ))
}
