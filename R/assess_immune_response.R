#' Assess predicted immune response performance
#'
#' `assess_immune_response` generates a roc curve plot and a barplot showing the average
#' (across tasks) area under the ROC curve (AUC) values for each quantitative descriptor on
#' the patients' response provided. For comparison, a gold standard is plotted based on the
#' average of the tasks.
#'
#' @importFrom ROCR prediction performance plot
#' @importFrom grDevices gray.colors pdf dev.off
#' @importFrom stats aggregate median sd
#' @import ggplot2
#' @importFrom graphics legend par title
#'
#' @export
#'
#' @param predictions_immune_response A list containing the predictions for each quantitative descriptor and for each task.
#' @param real_patient_response A character vector with two factors (Non-responders = NR, Responders = R).
#' @param RNA_tpm A numeric matrix of patients' gene expression data as tpm values.
#' @param output_file_path A character string pointing to a directory to save the plots returned by the function.
#' @param list_gold_standards A character string of task names to be considered as gold standards for comparison.
#' @param cancer_type A character string indicating which cancer-specific model should be used to compute the predictions.
#' @param TMB_values A numeric vector containing patients' tumor mutational burden (TMB) values.
#' @param easier_with_TMB A logical flag indicating whether to apply refined approach using the combination of easier predictions and tumor mutational burden..
#'
#' @return If easier_with_TMB is set to FALSE, two figures (roc curve and bar plots) are directly saved in the path specified in output_file_path.
#' If easier_with_TMB is set to TRUE, an additional plot is returned displaying an integrated approach that uses both immune response and tumor
#' mutational burden (TMB) to predict patients' response.
#'
#' @examples
#' # Example: Mariathasan cohort (Mariathasan et al., Nature, 2018)
#' if (!requireNamespace("BiocManager", quietly = TRUE))
#'  install.packages("BiocManager")
#'
#' BiocManager::install(c("biomaRt",
#'  "circlize",
#'  "ComplexHeatmap",
#'  "corrplot",
#'  "DESeq2",
#'  "dplyr",
#'  "DT",
#'  "edgeR",
#'  "ggplot2",
#'  "limma",
#'  "lsmeans",
#'  "reshape2",
#'  "spatstat",
#'  "survival",
#'  "plyr"))
#'
#' install.packages("Downloads/IMvigor210CoreBiologies_1.0.0.tar.gz", repos = NULL)
#' library(IMvigor210CoreBiologies)
#'
#' data(cds)
#' mariathasan_data <- preprocess_mariathasan(cds)
#' gene_count <- mariathasan_data$counts
#' gene_tpm <- mariathasan_data$tpm
#' rm(cds)
#'
#' # Computation of cell fractions
#' cell_fractions <- compute_cell_fractions(RNA_tpm = gene_tpm)
#'
#' # Computation of pathway scores
#' pathway_activity <- compute_pathways_scores(RNA_counts = gene_count,
#' remove_genes_ICB_proxies = TRUE)
#'
#' # Computation of TF activity
#' tf_activity <- compute_TF_activity(RNA_tpm = gene_tpm,
#' remove_genes_ICB_proxies = FALSE)
#'
#' # Computation of ligand-receptor pair weights
#' lrpair_weights <- compute_LR_pairs(RNA_tpm = gene_tpm,
#' remove_genes_ICB_proxies = FALSE,
#' cancer_type = "pancan")
#'
#' # Computation of cell-cell interaction scores
#' ccpair_scores <- compute_CC_pairs_grouped(lrpairs = lrpair_weights,
#' cancer_type = "pancan")
#'
#' # Predict patients' immune response
#' predictions_immune_response <- predict_immune_response(pathways = pathway_activity,
#' immunecells = cell_fractions,
#' tfs = tf_activity,
#' lrpairs = lrpair_weights,
#' ccpairs = ccpair_scores,
#' cancer_type = "SKCM",
#' verbose = TRUE)
#'
#' # retrieve clinical response
#' clinical_data <- mariathasan_data$clinical
#' patient_response <- clinical_data[, "Best Confirmed Overall Response"]
#' patient_response <- gsub("CR","R", patient_response)
#' patient_response <- gsub("PD","NR", patient_response)
#' names(patient_response) <- rownames(clinical_data)
#' # retrieve TMB
#' TMB <- clinical_data[, "FMOne mutation burden per MB"]
#' names(TMB) <- rownames(clinical_data)
#'
#' # Assess patient-specific likelihood of response to ICB therapy
#' assess_immune_response(predictions_immune_response = predictions_immune_response, real_patient_response = patient_response,
#' RNA_tpm = gene_tpm,
#' output_file_path = "../figures",
#' cancer_type = "SKCM",
#' TMB_values = TMB)
assess_immune_response <- function(predictions_immune_response = NULL,
                                   real_patient_response,
                                   RNA_tpm,
                                   output_file_path,
                                   list_gold_standards,
                                   cancer_type,
                                   TMB_values,
                                   easier_with_TMB = FALSE) {
  if (missing(cancer_type)) stop("cancer type needs to be specified")
  if (is.null(predictions_immune_response)) stop("none predictions found")
  if (missing(TMB_values)){
    TMB_available <- FALSE
    easier_with_TMB <- FALSE
  }else{
    TMB_available <- TRUE
    if (anyNA(TMB_values)) warning("NA values were found in TMB data, patients with NA values are removed from the analysis")
    message(paste0(", considering ", length(TMB_values[!is.na(TMB_values)]), " patients out of ", length(TMB_values)))
    patients_to_keep <- names(TMB_values[!is.na(TMB_values)])
    TMB_values <- TMB_values[patients_to_keep]
    real_patient_response <- real_patient_response[patients_to_keep]
    RNA_tpm <- RNA_tpm[, patients_to_keep]
  }
  # Check that folder exists, create folder otherwise
  if (dir.exists(output_file_path) == FALSE) {
    dir.create(file.path(output_file_path), showWarnings = FALSE)
    warning(paste0(
      sapply(strsplit(output_file_path, "/", fixed = TRUE), tail, 1),
      " folder does not exist, creating ", sapply(strsplit(output_file_path, "/", fixed = TRUE), tail, 1), " folder"
    ))
  }
  # Initialize variables
  views <- names(predictions_immune_response)
  tasks <- names(predictions_immune_response[[1]])
  # All views
  all_color_views <- vector("character", length = length(views))
  all_color_views <- c(
    "#52ac68", "#6a70d7", "#bbb442",
    "#5b3788", "#72a646"
  )
  names(all_color_views) <- views
  # Patient response labels
  labels <- matrix(real_patient_response,
    nrow = length(real_patient_response), ncol = length(tasks),
    dimnames = list(colnames(RNA_tpm), tasks)
  )
  # AUC predictions (when response available)
  if (missing(real_patient_response) == FALSE) {
    if (all(levels(as.factor(real_patient_response)) %in% c("NR", "R")) == FALSE) {
      stop("real_patient_response factor levels are not NR and R")
    }
    # Compute gold standards
    default_list_gold_standards <- c("CYT", "Roh_IS", "chemokines", "Davoli_IS", "IFNy", "Ayers_expIS", "Tcell_inflamed", "RIR", "TLS")
    if (missing(list_gold_standards)) {
      list_gold_standards <- default_list_gold_standards
    }
    gold_standards <- compute_gold_standards(RNA_tpm, list_gold_standards)
    # Assess correlation between chemokines and the other correlated tasks
    cor_tasks <- list_gold_standards[!list_gold_standards %in% c("IMPRES", "MSI")]
    cor_tasks <- cor_tasks[!cor_tasks %in% c("Ock_IS")] # Unfeasible computation
    tasks_values <- gold_standards[, cor_tasks]
    tasks_cormatrix <- cor(tasks_values)
    cor_sign <- sign(tasks_cormatrix[, "chemokines"])
    cor_sign <- cor_sign[names(cor_sign) != "chemokines"]
    if (all(cor_sign == -1)) {
      tasks_values[, "chemokines"] <- -tasks_values[, "chemokines"]
    }
    tasks_values <- as.data.frame(tasks_values)
    # ---------------------------#
    # Predictions single views #
    # ---------------------------#
    ROC_pred <- lapply(views, function(spec_view) {
        ROC_pred <- sapply(tasks, function(spec_task) {
            df <- predictions_immune_response[[spec_view]][[spec_task]]
            if(TMB_available) df[patients_to_keep, ]
            # check patients match
            df <- df[match(rownames(labels), rownames(df)), ]
            df_runs <- rowMeans(df)
        })
        pred <- ROCR::prediction(ROC_pred, labels, label.ordering = c("NR", "R"))
        perf <- ROCR::performance(pred, "tpr", "fpr")
        AUC <- unlist(ROCR::performance(pred, "auc")@y.values); names(AUC) <- tasks
        return(list(Curve = list(perf), Barplot = list(AUC)))
    })
    names(ROC_pred) <- views
    # Collect derived signatures predictions into data.frame #
    AUC_data <- do.call(rbind, lapply(views, function(spec_view) {
        return(data.frame(
          AUC = as.numeric(unlist(ROC_pred[[spec_view]]$Barplot)),
          View = spec_view,
          Task = names(unlist(ROC_pred[[spec_view]]$Barplot)),
          Run = "average"
        )
        )
    }))
    # Average across tasks
    AUC_mean_sd_run_tasks <- do.call(data.frame,
                                    aggregate(AUC ~ View, data = AUC_data, FUN = function(x) c(mean = mean(x), sd = sd(x))))
    # ---------------------------#
    # Predictions ensemble view #
    # ---------------------------#
    ensemble_df <- lapply(views, function(spec_view){
      ensemble_df <- sapply(tasks, function(spec_task){
        df <- predictions_immune_response[[spec_view]][[spec_task]]
        if(TMB_available) df[patients_to_keep, ]
        df_runs <- rowMeans(df)
      })
      return(ensemble_df)
    })
    names(ensemble_df) <- views
    overall_df <- sapply(tasks, function(spec_task){
      overall_df <- apply(cbind(ensemble_df$pathways[, spec_task],
                                ensemble_df$immunecells[, spec_task],
                                ensemble_df$tfs[, spec_task],
                                ensemble_df$lrpairs[, spec_task],
                                ensemble_df$ccpairs[, spec_task]), 1, mean)
    })
    # check patients match
    overall_df <- overall_df[match(rownames(labels), rownames(overall_df)), ]
    pred <- ROCR::prediction(overall_df, labels, label.ordering = c("NR", "R"))
    perf <- ROCR::performance(pred,"tpr","fpr")
    AUC <- unlist(ROCR::performance(pred, "auc")@y.values); names(AUC) <- tasks
    ensemble_ROC_pred <- list(ensemble=list(Curve = list(perf), Barplot = list(AUC)))
    # Collect predictions into data.frame #
    AUC_data_ensemble <- data.frame(
      AUC = as.numeric(unlist(ensemble_ROC_pred$ensemble$Barplot)),
      View = "ensemble",
      Task = names(unlist(ensemble_ROC_pred$ensemble$Barplot))
    )
    # Average across tasks
    AUC_mean_sd_ensemble_run_tasks <- do.call(data.frame,
                                         aggregate(AUC ~ View, data = AUC_data_ensemble, FUN = function(x) c(mean = mean(x), sd = sd(x))))
    # ---------------------------#
    # Predictions gold standards #
    # ---------------------------#
    # check patients match
    tasks_values <- tasks_values[match(rownames(labels), rownames(tasks_values)), ]
    pred <- ROCR::prediction(tasks_values, labels[, colnames(tasks_values)], label.ordering = c("NR", "R"))
    perf <- ROCR::performance(pred,"tpr","fpr")
    AUC  <- unlist(ROCR::performance(pred, "auc")@y.values); names(AUC) <- colnames(tasks_values)
    goldstandard_ROC_pred <- list(gold_standard=list(Curve = list(perf), Barplot = list(AUC)))

    AUC_data_goldstandard <- data.frame(
      AUC = as.numeric(unlist(goldstandard_ROC_pred$gold_standard$Barplot)),
      View = "gold_standard",
      Task = names(unlist(goldstandard_ROC_pred$gold_standard$Barplot))
    )
    # Average across tasks
    AUC_mean_sd_goldstandard_run_tasks <- do.call(data.frame,
                                                  aggregate(AUC ~ View, data = AUC_data_goldstandard, FUN = function(x) c(mean = mean(x), sd = sd(x))))
    # ------------------------------------------------------#
    # Combine single views, ensemble and gold standard
    # ------------------------------------------------------#
    AUC_mean_sd_all_run_tasks <- rbind(AUC_mean_sd_run_tasks, AUC_mean_sd_ensemble_run_tasks, AUC_mean_sd_goldstandard_run_tasks)
    ROC_all_run_tasks <- c(ROC_pred, ensemble_ROC_pred, goldstandard_ROC_pred); names(ROC_all_run_tasks) <- c(names(all_color_views), "ensemble", "gold_standard")
    AUC_mean_sd_all_run_tasks$View <- factor(AUC_mean_sd_all_run_tasks$View,
                                             levels = c(names(all_color_views), "ensemble", "gold_standard"))

    # Colors gold standard and ensemble
    color_gold_standard <- "gray82"; names(color_gold_standard) <- "gold_standard"
    color_ensemble <- "gold2"; names(color_ensemble) <- "ensemble"
    # ------------------------------------------------------#
    # Add TMB if available
    # ------------------------------------------------------#
    if (TMB_available == TRUE){
      # check patients match
      TMB_values <- TMB_values[match(rownames(labels), names(TMB_values))]
      pred <- ROCR::prediction(TMB_values, labels[,1], label.ordering = c("NR", "R"))
      perf <- ROCR::performance(pred, "tpr", "fpr")
      AUC  <- unlist(ROCR::performance(pred, "auc")@y.values)
      TMB_ROC_pred <- list(TMB=list(Curve = list(perf), Barplot = list(AUC)))

      # Collect predictions into data.frame #
      AUC_data_TMB <- data.frame(
        AUC = as.numeric(unlist(TMB_ROC_pred$TMB$Barplot)),
        View = "TMB",
        Task = "TMB"
      )
      # Average across tasks
      AUC_mean_sd_TMB_run_tasks <- do.call(data.frame,
                                                    aggregate(AUC ~ View, data = AUC_data_TMB, FUN = function(x) c(mean = mean(x), sd = sd(x))))
      AUC_mean_sd_all_run_tasks <- rbind(AUC_mean_sd_all_run_tasks, AUC_mean_sd_TMB_run_tasks)
      ROC_all_run_tasks <- c(ROC_all_run_tasks, TMB_ROC_pred); names(ROC_all_run_tasks) <- c(names(all_color_views), "ensemble", "gold_standard", "TMB")
      AUC_mean_sd_all_run_tasks$View <- factor(AUC_mean_sd_all_run_tasks$View,
                                               levels = c(names(all_color_views), "ensemble", "gold_standard", "TMB"))
      # Colors TMB
      color_TMB <- "salmon"; names(color_TMB) <- "salmon"
    }
    # *******************************************
    # Plot AUC values using barplot
    # *******************************************
    n_R <- table(real_patient_response)[["R"]]
    n_NR <- table(real_patient_response)[["NR"]]
    gg <- ggplot2::ggplot(AUC_mean_sd_all_run_tasks, ggplot2::aes(x = .data$View, y = round(.data$AUC.mean, 2), fill = .data$View)) +
      ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(), color = "white") +
      if (TMB_available){
        ggplot2::scale_fill_manual(values = c(
          as.vector(all_color_views), as.vector(color_ensemble),
          as.vector(color_gold_standard), as.vector(color_TMB)
        ), guide = FALSE)
      }else{
        ggplot2::scale_fill_manual(values = c(
          as.vector(all_color_views), as.vector(color_ensemble),
          as.vector(color_gold_standard)
        ), guide = FALSE)
      }
    gg + ggplot2::scale_x_discrete(labels = c(
        "ensemble" = "Ensemble",
        "immunecells" = "Cell fractions",
        "pathways" = "Pathways",
        "tfs" = "TFs",
        "lrpairs" = "LR pairs",
        "ccpairs" = "CC pairs",
        "gold_standard" = "Tasks (gold standard)"
        )) +
      ggplot2::theme(panel.grid = ggplot2::element_blank(), panel.background = ggplot2::element_rect(fill = NA)) +
      ggplot2::theme_bw() +
      ggplot2::ylim(0, 1) +
      ggplot2::ylab("Area under the curve (AUC)") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = round(.data$AUC.mean, 2) - .data$AUC.sd, ymax = round(.data$AUC.mean, 2) + .data$AUC.sd), width = .3, color="black", position = ggplot2::position_dodge(0.9)) +
      ggplot2::geom_text(ggplot2::aes(label= round(.data$AUC.mean, 2)), stat = "identity", color="black", size = 4, angle = 90, hjust = -0.5, position = ggplot2::position_dodge(0.9)) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(size=12, angle = 45, vjust = 1, hjust = 1, color = "black"),
        axis.text.y = ggplot2::element_text(size = 12, color = "black"),
        axis.title.y = ggplot2::element_text(size = 12), axis.title.x = ggplot2::element_blank(),
        legend.position = "right", legend.text = ggplot2::element_text(size = 10),
        legend.box.background = ggplot2::element_rect(color = "black", size = 0.3),
        legend.box.margin = ggplot2::margin(0.5, 0.5, 0.5, 0.5)
      ) +
      ggplot2::labs(title = paste0(" n=", length(real_patient_response), " (R=", n_R, "; ", "NR=", n_NR, ")"))

    ggplot2::ggsave(file.path(output_file_path, "auc_barplot.pdf"), width = 7, height = 5)

    # *******************************************
    # Plot ROC curve
    # *******************************************
    pdf(file.path(output_file_path, "roc_curve.pdf"), width = 8, height = 8)
    par(cex.axis = 1.6, mar = c(5, 5, 5, 5), col.lab = "black")

    # Single views
    ROCR::plot(ROC_all_run_tasks$pathways$Curve[[1]],
      avg = "threshold", col = all_color_views["pathways"], lwd = 3, type = "S",
      cex.lab = 1.6, ylab = "True Positive Rate", xlab = "False Positive Rate"
    )

    ROCR::plot(ROC_all_run_tasks$immunecells$Curve[[1]],
      avg = "threshold", col = all_color_views["immunecells"], lwd = 3, type = "S",
      cex.lab = 1.6, ylab = "True Positive Rate", xlab = "False Positive Rate", add = TRUE
    )

    ROCR::plot(ROC_all_run_tasks$tfs$Curve[[1]],
      avg = "threshold", col = all_color_views["tfs"], lwd = 3, type = "S",
      cex.lab = 1.6, ylab = "True Positive Rate", xlab = "False Positive Rate", add = TRUE
    )

    ROCR::plot(ROC_all_run_tasks$lrpairs$Curve[[1]],
      avg = "threshold", col = all_color_views["lrpairs"], lwd = 3, type = "S",
      cex.lab = 1.6, ylab = "True Positive Rate", xlab = "False Positive Rate", add = TRUE
    )

    ROCR::plot(ROC_all_run_tasks$ccpairs$Curve[[1]],
      avg = "threshold", col = all_color_views["ccpairs"], lwd = 3, type = "S",
      cex.lab = 1.6, ylab = "True Positive Rate", xlab = "False Positive Rate", add = TRUE
    )

    # ensemble
    ROCR::plot(ROC_all_run_tasks$ensemble$Curve[[1]],
      avg = "threshold", col = color_ensemble, lwd = 3, type = "S",
      lty = 1, add = TRUE)

    # gold standard
    ROCR::plot(ROC_all_run_tasks$gold_standard$Curve[[1]],
      avg = "threshold", col = color_gold_standard, lwd = 3, type = "S",
      lty = 1, add = TRUE)

    # TMB
    if (TMB_available){
      ROCR::plot(ROC_all_run_tasks$TMB$Curve[[1]], col = color_TMB, lwd = 3, type = "S",
        lty = 1, add = TRUE)

      legend(
        x = 0.65, y = 0.27,
        legend = c(
          paste0("Pathways", " (", round(subset(AUC_mean_sd_all_run_tasks, View == "pathways")$AUC.mean, 2), ")"),
          paste0("Cell fractions", " (", round(subset(AUC_mean_sd_all_run_tasks, View == "immunecells")$AUC.mean, 2), ")"),
          paste0("TFs", " (", round(subset(AUC_mean_sd_all_run_tasks, View == "tfs")$AUC.mean, 2), ")"),
          paste0("LR pairs", " (", round(subset(AUC_mean_sd_all_run_tasks, View == "lrpairs")$AUC.mean, 2), ")"),
          paste0("CC pairs", " (", round(subset(AUC_mean_sd_all_run_tasks, View == "ccpairs")$AUC.mean, 2), ")"),
          paste0("Ensemble", " (", round(subset(AUC_mean_sd_all_run_tasks, View == "ensemble")$AUC.mean, 2), ")"),
          paste0("Tasks", " (", round(subset(AUC_mean_sd_all_run_tasks, View == "gold_standard")$AUC.mean, 2), ")"),
          paste0("TMB", " (", round(subset(AUC_mean_sd_all_run_tasks, View == "TMB")$AUC.mean, 2), ")")
        ),
        col = c(
          all_color_views, as.vector(color_ensemble), as.vector(color_gold_standard), color_TMB
        ), lty = 1, lwd = 3, cex = 0.9, bty = "n"
      )
    }else{
      legend(
        x = 0.65, y = 0.27,
        legend = c(
          paste0("Pathways", " (", round(subset(AUC_mean_sd_all_run_tasks, View == "pathways")$AUC.mean, 2), ")"),
          paste0("Cell fractions", " (", round(subset(AUC_mean_sd_all_run_tasks, View == "immunecells")$AUC.mean, 2), ")"),
          paste0("TFs", " (", round(subset(AUC_mean_sd_all_run_tasks, View == "tfs")$AUC.mean, 2), ")"),
          paste0("LR pairs", " (", round(subset(AUC_mean_sd_all_run_tasks, View == "lrpairs")$AUC.mean, 2), ")"),
          paste0("CC pairs", " (", round(subset(AUC_mean_sd_all_run_tasks, View == "ccpairs")$AUC.mean, 2), ")"),
          paste0("Ensemble", " (", round(subset(AUC_mean_sd_all_run_tasks, View == "ensemble")$AUC.mean, 2), ")"),
          paste0("Tasks", " (", round(subset(AUC_mean_sd_all_run_tasks, View == "gold_standard")$AUC.mean, 2), ")")
        ),
        col = c(
          all_color_views, as.vector(color_ensemble), as.vector(color_gold_standard)
        ), lty = 1, lwd = 3, cex = 0.9, bty = "n"
      )
    }
    title(main = paste0("n=", length(real_patient_response), " (R=", n_R, "; ", "NR=", n_NR, ")"))
    dev.off()

    # *******************************************
    # Integrated score
    # *******************************************
    if (combo_easier_and_TMB == TRUE){
      rp_df <- data.frame(response = real_patient_response,
                          prediction_easier = apply(overall_df, 1, mean),
                          TMB = TMB_values)
      # Categorize TMB
      if (length(unique(rp_df$TMB))>3){
        rp_df$TMBcat <- categorize_TMB(rp_df$TMB)
        # rp.df$TMB <- categorize.TMB(rp.df$TMB, thresholds = c(100,400)) # if I want specify the thresholds
      }
      # compute the integrated score for different penalties #
      pred_combined <- rp_df$prediction_easier
      AUC_combined_v <- sapply(seq(from=0, to=1, by=0.1), function(p){

        pred_combined[rp_df$TMBcat==1] <- pred_combined[rp_df$TMBcat==1] - p
        pred_combined[rp_df$TMBcat==3] <- pred_combined[rp_df$TMBcat==3] + p

        pred <- ROCR::prediction(pred_combined, rp_df$response)
        AUC_combined <-  ROCR::performance(pred, measure = "auc")
        AUC_combined_v <- AUC_combined@y.values[[1]]
      })
      # compute the integrated score as weighted average #
      linear_func <- function(x) {
        min_x <- min(x)
        max_x <- max(x)
        x01 = (x - min_x) / (max_x - min_x)
        return(x01)
      }
      pred_lin <- linear_func(rp_df$prediction_easier)
      TMB_lin <- linear_func(rp_df$TMB)

      AUC_averaged_v <- sapply(seq(from=0, to=1, by=0.1), function(p){
        pred_averaged <- apply(cbind((1-p)*pred_lin, (p)*TMB_lin), 1, mean)
        pred <- ROCR::prediction(pred_averaged, rp_df$response)
        AUC_averaged <- ROCR::performance(pred, measure = "auc")
        AUC_averaged_v <- AUC_averaged@y.values[[1]]
      })

      pdf(file.path(output_file_path, "easier_tmb_combo_auc.pdf"), width = 8, height = 8)
      par(cex.axis = 1.6, mar = c(5, 5, 5, 5), col.lab = "black")

      plot(seq(from=0, to=1, by=0.1), AUC_combined_v, xlab="Penalty or Relative weight", ylab="Area under the curve (AUC)", type = "b", col="#c15050", lty = 1, pch = 19, lwd = 3, ylim=c(0.5,1),  cex.lab = 1.6)
      lines(seq(from=0, to=1, by=0.1), AUC_averaged_v, xlab="Penalty or Relative weight", ylab="Area under the curve (AUC)", type = "b", col="#693c72", lty = 1, pch = 19, lwd = 3, ylim=c(0.5,1),  cex.lab = 1.6)
      # TMB
      abline(h = AUC_mean_sd_TMB_run_tasks$AUC.mean, col = color_TMB)
      # easier (ensemble)
      abline(h = AUC_mean_sd_ensemble_run_tasks$AUC.mean, col = color_ensemble)

      legend(
        x = 0.55, y = 1,
        legend = c("Penalized score", "Weighted average", "EaSIeR", "TMB"
        ),
        col = c(
          "#c15050", "#693c72",  as.vector(color_ensemble), color_TMB
        ), lty = 1, lwd = 3, cex = 1.3, bty = "n"
      )
      dev.off()
    }

  } else {
    stop("No patients' response provided")
    # all_scores_df <- data.frame(weighted_average = pred_averaged,
    #                             penalized_score = pred_combined,
    #                             patient = names(TMB_values))
    # all_scores_df <- reshape2::melt(all_scores_df)
    # names(all_scores_df) <- c("patient", "approach", "pred")
    #
    # all_scores_df$approach <- factor(all_scores_df$approach, levels = c("penalized_score", "weighted_average"))
    #
    # ggplot2::ggplot(all_scores_df, ggplot2::aes(x=pred, y=patient, fill = approach, color = approach)) +
    #   ggplot2::geom_point(size = 2) +
    #   ggplot2::scale_fill_manual(name = "Approach",
    #                              labels = as.character(unique(all_scores_df$approach)),
    #                              values = c("#c15050", "#693c72")) +
    #   ggplot2::scale_color_manual(name = "Approach",
    #                               labels = as.character(unique(all_scores_df$approach)),
    #                               values = c("#c15050", "#693c72")) +
    #   ggplot2::theme(panel.grid = ggplot2::element_blank(), panel.background = ggplot2::element_rect(fill = NA)) +
    #   ggplot2::geom_vline(xintercept = 0) +
    #   ggplot2::theme_linedraw() +
    #   #ggplot2::xlim(-3, 3) +
    #   ggplot2::theme(axis.text.x = ggplot2::element_text(size=10,face="bold", angle = 0, vjust = 0.5, hjust=0.5, color = "black"), axis.title.x = ggplot2::element_blank(),
    #                  axis.text.y = ggplot2::element_text(size=10), axis.ticks.y = ggplot2::element_blank(),
    #                  axis.title.y = ggplot2::element_blank(), panel.background = ggplot2::element_rect(fill = "white"),
    #                  legend.position = "right", legend.direction = "vertical", panel.grid.major.y = ggplot2::element_line(linetype = 1, colour = "gray"),
    #                  legend.text = ggplot2::element_text(size=10), legend.title = ggplot2::element_text(size = 10, face="bold", vjust = 0.5),
    #                  strip.background = ggplot2::element_rect(fill = "white", colour = "white"), strip.text = ggplot2::element_text(size=10)) +
    #   ggplot2::labs(x = "prediction", y = "patients", title = "Scoring patients")
    #
    # ggsave(paste0("/Users/Oscar/Desktop/Kim_patient_score_original.pdf"), width = 12, height = 10)
  }
}

