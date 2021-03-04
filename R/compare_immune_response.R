#' Comparison of the actual predictions based on different metrics with real patient response.
#'
#' Plot ROC curve and barplot with the area under the ROC curve.
#' If list of gold standards is not provided, the function used a default one.
#'
#' `compare_response` plots ROC curves and barplots showing the accuracy of the predictions
#' on the real patient response data. It receives as input the predicted immune response as well
#' as the real patient response (they should be provided with the same order of samples). Gold
#' standards are also plotted to consider them as reference. To compute these gold standards,
#' transcriptomics data in different formats is provided. An output file path should be given
#' to save the plots as pdfs. An additional input must be a list with the gold standards the
#' user wish to compare the prediction results.
#'
#' @importFrom ROCR prediction performance plot
#' @importFrom grDevices gray.colors pdf dev.off
#' @importFrom stats aggregate median sd
#' @import ggplot2
#' @importFrom graphics legend par title
#'
#' @export
#'
#' @param predictions_immune_response list that contains the predictions for each model
#' @param real_patient_response vector with two factors (NR,R)
#' @param RNA_tpm numeric matrix with data
#' @param output_file_path string with a file name
#' @param list_gold_standards string with gold standards names
#' @param cancertype specify cancer type
#'
#' @return ROC curves plots and barplots showing AUC values.
#'
#' @examples
#' # Example: Riaz
#' data("Riaz_data")
#'
#' # Computation of cell fractions
#' cell_fractions <- compute_cell_fractions(RNA_tpm = Riaz_data$tpm_RNAseq)
#'
#' # Computation of pathway scores
#' pathway_activity <- compute_pathways_scores(
#'   RNA_counts = Riaz_data$raw_counts_RNAseq,
#'   remove_genes_ICB_proxies = TRUE)
#'
#' # Computation of TF activity
#' tf_activity <- compute_TF_activity(
#'   RNA_tpm = Riaz_data$tpm_RNAseq,
#'   remove_genes_ICB_proxies = FALSE)
#'
#' # Computation of ligand-receptor pair weights
#' lrpairs_weights <- compute_LR_pairs(
#'   RNA_tpm = Riaz_data$tpm_RNAseq,
#'   remove_genes_ICB_proxies = FALSE,
#'   cancertype = "pancan")
#'
#' # Computation of cell-cell interaction scores
#' ccpairs_scores <- compute_CC_pairs_grouped(
#'   lrpairs = lrpairs_weights, # TODO: I removed the $LRpairs as it was not present?
#'   cancertype = "pancan")
#'
#' # Predict patients' immune response
#' # TODO: error raised in predict_immune_response
#' #####predictions_immune_response <- predict_immune_response(
#' #####  pathways = pathway_activity,
#' #####  immunecells = cell_fractions,
#' #####  lrpairs = lrpairs_weights,
#' #####  tfs = tf_activity,
#' #####  ccpairs = ccpairs_scores,
#' #####  include_pairwise_combos = FALSE,
#' #####  cancertype = "SKCM")
#' #####
#' ###### Assess patient-specific likelihood of response to ICB therapy
#' #####compare_immune_response(
#' #####  predictions_immune_response = predictions_immune_response,
#' #####  real_patient_response = Riaz_data$patient_response,
#' #####  RNA_tpm = Riaz_data$tpm_RNAseq,
#' #####  output_file_path = "/Users/Oscar/Desktop/Riaz",
#' #####  cancertype = "SKCM")
compare_immune_response <- function(predictions_immune_response = NULL,
                                    real_patient_response,
                                    RNA_tpm,
                                    output_file_path,
                                    list_gold_standards,
                                    cancertype) {
  if (missing(cancertype)) stop("cancer type needs to be specified")
  if (is.null(predictions_immune_response)) stop("none predictions found")

  # Check that folder exists, create folder otherwise
  if (dir.exists(output_file_path) == FALSE) {
    dir.create(file.path(output_file_path), showWarnings = FALSE)
    warning(paste0(
      sapply(strsplit(output_file_path, "/", fixed = TRUE), tail, 1),
      " folder does not exist, creating ", sapply(strsplit(output_file_path, "/", fixed = TRUE), tail, 1), " folder"
    ))
  }

  # Initialize variables
  AUC_median <- AUC_sd <- Model <- Task <- View <- NULL
  view_combinations <- names(predictions_immune_response)
  algorithms <- names(predictions_immune_response[[1]])
  tasks <- c(names(predictions_immune_response[[1]][[1]]), "common_mean", "common_median")
  models <- names(predictions_immune_response[[1]][[1]][[1]])
  # sel_views <- names(predictions_immune_response)

  # Derived signature colors:
  # Spatial info: "#ff983d"
  # Protein: "darkorange2"
  # Pathways: "#6CD8CB"
  # Immune cells: "#CC6AF2"
  # Pathways + Immune cells: "#0392f2"
  # Cytokine or LR pairs: "salmon"
  # Cell-cell pairs:"#9E1D3F"
  # Transcription factors: "#01AD3F"
  # Transcriptomics: "black"
  # Consensus: grey (just a bit darker than now, can have same line thikness as the other)

  # All views
  all_color_views <- vector("character", length = length(view_combinations))
  all_color_views <- c(
    "#52ac68", "#6a70d7", "#bbb442",
    "#5b3788", "#72a646", "#c972c4",
    "#43c8ac", "#d25543", "#6d8dd7",
    "#cd8232", "#b1457b", "#9b843b",
    "#ba4758", "#a24e2e"
  )
  names(all_color_views) <- view_combinations

  # Selected views
  # sel_color_views <- all_color_views[names(all_color_views) %in% sel_views]

  # Tasks
  color_tasks <- vector("character", length = 14)
  color_tasks <- toupper(c(
    "#7763ce", "#69b444", "#c754bf", "#c7ad3e", "#6484c8", "#d35238", "#4ac0cd",
    "#d74278", "#4fa471", "#9b4d8a", "#7e843b", "#d48dca", "#c07c42", "#bd5f68"
  ))
  names(color_tasks) <- c(
    "CYT", "Ock_IS", "IPS", "IMPRES", "Roh_IS", "chemokines", "Davoli_IS", "IFNy", "Ayers_expIS",
    "Tcell_inflamed", "TIDE", "MSI", "RIR", "TLS"
  )

  # Patient response labels
  labels <- matrix(real_patient_response,
    nrow = length(real_patient_response), ncol = 100,
    dimnames = list(colnames(RNA_tpm), seq(1, 100, 1))
  )

  # AUC predictions (when response available)

  if (missing(real_patient_response) == FALSE) {
    if (all(levels(as.factor(real_patient_response)) %in% c("NR", "R")) == FALSE) {
      stop("real_patient_response factor levels are not NR and R")
    }

    # Compute gold standards
    default_list_gold_standards <- names(color_tasks)
    if (missing(list_gold_standards)) {
      list_gold_standards <- default_list_gold_standards
    }
    gold_standards <- compute_gold_standards(RNA_tpm, list_gold_standards, cancertype, output_file_path)

    # Assess correlation between chemokines and the other correlated tasks
    cor_tasks <- names(color_tasks)[!names(color_tasks) %in% c("IPS", "IMPRES", "TIDE", "MSI")]
    cor_tasks <- cor_tasks[!cor_tasks %in% c("Ock_IS")] # Unfeasible computation
    tasks_values <- do.call(cbind, lapply(cor_tasks, function(X) {
      tmp <- as.numeric(unlist(gold_standards[[X]]))
    }))
    colnames(tasks_values) <- cor_tasks
    rownames(tasks_values) <- colnames(gold_standards$CYT)
    tasks_cormatrix <- cor(tasks_values)
    cor_sign <- sign(tasks_cormatrix[, "chemokines"])
    cor_sign <- cor_sign[names(cor_sign) != "chemokines"]
    if (all(cor_sign == -1)) {
      tasks_values[, "chemokines"] <- -tasks_values[, "chemokines"]
    }
    tasks_values <- as.data.frame(tasks_values)

    # Tasks normalization
    all_tasks_values_norm <- as.data.frame(standardization(tasks_values))
    gold_standards_scaled <- all_tasks_values_norm

    # Predictions
    ROC_info <- lapply(c(view_combinations), function(view) {
      ROC_info <- lapply(algorithms, function(alg) {
        ROC_info <- lapply(c(names(predictions_immune_response[[view]][[alg]]), "common_mean", "common_median"), function(task) {
          ROC_info <- lapply(models, function(model) {
            if (alg != "BEMKL") {
              if (task %in% "common_mean") {
                df <- list(matrix(0, nrow(labels), ncol = 100))

                for (run in 1:ncol(labels)) {
                  df[[1]][, run] <- apply(rbind(
                    predictions_immune_response[[view]][[alg]][["CYT"]][[model]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["IS"]][[model]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["RohIS"]][[model]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["chemokine"]][[model]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["IS_Davoli"]][[model]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["IFny"]][[model]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["ExpandedImmune"]][[model]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["T_cell_inflamed"]][[model]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["resF.down"]][[model]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["TLS"]][[model]][[1]][, run]
                  ), 2, mean)
                }
              } else if (task %in% "common_median") {
                df <- list(matrix(0, nrow(labels), ncol = 100))

                for (run in 1:ncol(labels)) {
                  df[[1]][, run] <- apply(rbind(
                    predictions_immune_response[[view]][[alg]][["CYT"]][[model]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["IS"]][[model]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["RohIS"]][[model]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["chemokine"]][[model]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["IS_Davoli"]][[model]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["IFny"]][[model]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["ExpandedImmune"]][[model]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["T_cell_inflamed"]][[model]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["resF.down"]][[model]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["TLS"]][[model]][[1]][, run]
                  ), 2, median)
                }
              } else {
                df <- predictions_immune_response[[view]][[alg]][[task]][[model]]
              }
            } else {
              if (task %in% "common_mean") {
                df <- list(matrix(0, nrow(labels), ncol = 100))

                for (run in 1:ncol(labels)) {
                  df[[1]][, run] <- apply(rbind(
                    predictions_immune_response[[view]][[alg]][["CYT"]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["IS"]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["RohIS"]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["chemokine"]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["IS_Davoli"]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["IFny"]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["ExpandedImmune"]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["T_cell_inflamed"]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["resF.down"]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["TLS"]][[1]][, run]
                  ), 2, mean)
                }
              } else if (task %in% "common_median") {
                df <- list(matrix(0, nrow(labels), ncol = 100))

                for (run in 1:ncol(labels)) {
                  df[[1]][, run] <- apply(rbind(
                    predictions_immune_response[[view]][[alg]][["CYT"]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["IS"]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["RohIS"]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["chemokine"]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["IS_Davoli"]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["IFny"]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["ExpandedImmune"]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["T_cell_inflamed"]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["resF.down"]][[1]][, run],
                    predictions_immune_response[[view]][[alg]][["TLS"]][[1]][, run]
                  ), 2, median)
                }
              } else {
                df <- predictions_immune_response[[view]][[alg]][[task]]
              }
            }

            pred <- ROCR::prediction(df[[1]], labels, label.ordering = c("NR", "R"))
            perf <- ROCR::performance(pred, "tpr", "fpr")
            AUC <- unlist(ROCR::performance(pred, "auc")@y.values)

            data_ROC <- list(perf)
            Barplot <- list(AUC)
            return(list(Curve = data_ROC, Barplot = Barplot))
          })
          names(ROC_info) <- models
          return(ROC_info)
        })
        names(ROC_info) <- c(names(predictions_immune_response[[view]][[alg]]), "common_mean", "common_median")
        return(ROC_info)
      })
      names(ROC_info) <- algorithms
      return(ROC_info)
    })
    names(ROC_info) <- c(view_combinations)

    # Predictions #
    overall_types <- "overall_mean_single"
    overall_ROC_info <- lapply(overall_types, function(main_view){
      overall_ROC_info <- lapply(algorithms, function(alg){
        overall_ROC_info <- lapply("overall", function(task){
          overall_ROC_info <- lapply(models, function(model){
            df_overall <- lapply(view_combinations, function(view){
              df <- matrix(0, nrow = nrow(labels), ncol = 100)

              if (alg != "BEMKL"){ # Only for now RMTLR implemented
                if (main_view %in% c("overall_mean_single")){
                  df_1 <- rowMeans(predictions_immune_response[[view]][[alg]][[1]][[model]][[1]])
                  df_2 <- rowMeans(predictions_immune_response[[view]][[alg]][[2]][[model]][[1]])
                  df_3 <- rowMeans(predictions_immune_response[[view]][[alg]][[3]][[model]][[1]])
                  df_4 <- rowMeans(predictions_immune_response[[view]][[alg]][[4]][[model]][[1]])
                  df_5 <- rowMeans(predictions_immune_response[[view]][[alg]][[5]][[model]][[1]])
                  df_6 <- rowMeans(predictions_immune_response[[view]][[alg]][[6]][[model]][[1]])
                  df_7 <- rowMeans(predictions_immune_response[[view]][[alg]][[7]][[model]][[1]])
                  df_8 <- rowMeans(predictions_immune_response[[view]][[alg]][[8]][[model]][[1]])
                  df_9 <- rowMeans(predictions_immune_response[[view]][[alg]][[9]][[model]][[1]])
                  df_10 <- rowMeans(predictions_immune_response[[view]][[alg]][[10]][[model]][[1]])

                  df_overall <- cbind(df_1, df_2,df_3,df_4,df_5,df_6, df_7, df_8, df_9, df_10)
                }
              }else{
                if (main_view %in% c("overall_mean_single", "overall_mean_combo", "overall_mean_all")){
                  df_1 <- rowMeans(predictions_immune_response[[view]][[alg]][[1]][[1]])
                  df_2 <- rowMeans(predictions_immune_response[[view]][[alg]][[2]][[1]])
                  df_3 <- rowMeans(predictions_immune_response[[view]][[alg]][[3]][[1]])
                  df_4 <- rowMeans(predictions_immune_response[[view]][[alg]][[4]][[1]])
                  df_5 <- rowMeans(predictions_immune_response[[view]][[alg]][[5]][[1]])
                  df_6 <- rowMeans(predictions_immune_response[[view]][[alg]][[6]][[1]])
                  df_7 <- rowMeans(predictions_immune_response[[view]][[alg]][[7]][[1]])
                  df_8 <- rowMeans(predictions_immune_response[[view]][[alg]][[8]][[1]])
                  df_9 <- rowMeans(predictions_immune_response[[view]][[alg]][[9]][[1]])
                  df_10 <- rowMeans(predictions_immune_response[[view]][[alg]][[10]][[1]])

                  df_overall <- cbind(df_1, df_2,df_3,df_4,df_5,df_6, df_7, df_8, df_9, df_10)
                }
              }
              return(df_overall)
            })

            names(df_overall) <- view_combinations
            overall_df <- matrix(0, nrow = nrow(labels), ncol = 10)
            for (task in 1:10){
              overall_df[,task] <- apply(rbind(
                df_overall$Pathways.cor[, task], df_overall$ImmuneCells[, task], df_overall$LRpairs.spec.pc[, task],
                                               df_overall$TFs[, task], df_overall$CCpairsGroupedScores.spec.pc[, task]
              ), 2, mean)
            }

            pred <- ROCR::prediction(overall_df, labels[, 1:10], label.ordering = c("NR", "R"))
            perf <- ROCR::performance(pred,"tpr","fpr")
            AUC <- unlist(ROCR::performance(pred, "auc")@y.values)

            data_ROC <- list(perf)
            Barplot <- list(AUC)
            return(list(Curve = data_ROC, Barplot = Barplot))
          })
          names(overall_ROC_info) <- models
          return(overall_ROC_info)
        })
        names(overall_ROC_info) <- "overall"
        return(overall_ROC_info)
      })
      names(overall_ROC_info) <- algorithms
      return(overall_ROC_info)
    })
    names(overall_ROC_info) <- overall_types

    # Combine ROC data for views and overall
    ROC_info <- c(ROC_info, overall_ROC_info)

    # Gold standards #
    ROC_info_GS <- lapply(colnames(gold_standards_scaled), function(GS) {
      pred_GS <- ROCR::prediction(gold_standards_scaled[,GS], labels[,1], label.ordering = c("NR", "R"))
      perf_GS <- ROCR::performance(pred_GS, "tpr", "fpr")
      AUC_GS <- unlist(ROCR::performance(pred_GS, "auc")@y.values)

      data_ROC_GS <- list(perf_GS)
      Barplot_GS <- list(AUC_GS)

      return(list(Curve = data_ROC_GS, Barplot = Barplot_GS))
    })
    names(ROC_info_GS) <- paste0(colnames(gold_standards_scaled), ".GS")

    # Collect derived signatures predictions into data.frame #
    AUC_mean_sd <- do.call(rbind, lapply(names(ROC_info), function(view) {
      AUC_mean_sd <- do.call(rbind, lapply(names(ROC_info[[view]]), function(alg) {
        AUC_mean_sd <- do.call(rbind, lapply(names(ROC_info[[view]][[alg]]), function(task) {
          AUC_mean_sd <- do.call(rbind, lapply(names(ROC_info[[view]][[alg]][[task]]), function(model) {
            View <- rep(view, times = 1)
            Alg <- rep(alg, times = 1)
            Model <- rep(model, times = 1)
            Task <- rep(task, times = 1)
            cv_iter <- "average"

            AUC_data <- data.frame(
              Alg = Alg,
              Model = Model,
              AUC = as.numeric(unlist(ROC_info[[view]][[alg]][[task]][[model]]$Barplot)),
              View = View,
              Task = Task,
              Iteration = cv_iter
            )
            AUC_mean_sd <- AUC_data

            return(AUC_mean_sd)
          }))
          return(AUC_mean_sd)
        }))
        return(AUC_mean_sd)
      }))
      return(AUC_mean_sd)
    }))

    # Calculate mean and sd across tasks
    AUC_mean_sd_across_tasks <- do.call(
      data.frame,
      aggregate(AUC ~ Alg + Model + View + Iteration,
        data = AUC_mean_sd, FUN = function(x)
          c(
            median = mean(x),
            sd = sd(x)
            )
        )
      )

    # Collect gold standards predictions into data.frame #
    AUC_mean_sd_GS <- do.call(rbind, lapply(names(ROC_info_GS), function(view) {
      Task <- rep(sapply(strsplit(view, split = ".", fixed = TRUE), head, 1), times = 1)
      Alg <- rep(algorithms, times = 1)
      Model <- rep("1se.mse", times = 1)
      View <- rep("Gold_Standard", times = 1)
      cv_iter <- "average"

      AUC_data <- data.frame(
        Alg = Alg,
        Model = Model,
        AUC = as.numeric(unlist(ROC_info_GS[[view]]$Barplot)),
        View = View,
        Task = Task,
        Iteration = cv_iter
      )
      AUC_mean_sd_GS <- AUC_data

      return(AUC_mean_sd_GS)
    }))

    # Calculate mean and sd across tasks
    AUC_mean_sd_GS_across_tasks <- do.call(
      data.frame,
      aggregate(AUC ~ Alg + Model + View + Iteration,
                data = AUC_mean_sd_GS, FUN = function(x) {
                  c(
                    median = median(x),
                    sd = sd(x)
                  )
                }
      )
    )

    AUC_mean_sd <- rbind(AUC_mean_sd_across_tasks, AUC_mean_sd_GS_across_tasks)
    AUC_mean_sd$View <- factor(AUC_mean_sd$View, levels = unique(AUC_mean_sd$View))

    # Keep more regularized model
    AUC_mean_sd_RMTLR <- subset(AUC_mean_sd, Model == "1se.mse")
    AUC_mean_sd_RMTLR$View <- factor(AUC_mean_sd_RMTLR$View,
                                   levels = c(names(all_color_views), overall_types, "Gold_Standard"))

    # Colors Gold standards and Overall
    color_gold_standard <- "gray82"; names(color_gold_standard) <- "Gold Standard"
    color_overalls <- "gold2"
    names(color_overalls) <- c( "Overall (mean) single views")

    n_R <- table(real_patient_response)[["R"]]
    n_NR <- table(real_patient_response)[["NR"]]

    # *******************************************
    # Barplot AUC values

    ggplot2::ggplot(AUC_mean_sd_RMTLR, ggplot2::aes(x = View, y = round(AUC_median, 2), fill = View,  alpha = Alg)) +
      ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(), color = "white") +
      ggplot2::scale_fill_manual(values = c(
        as.vector(all_color_views), as.vector(color_overalls),
        as.vector(color_gold_standard)
      ), guide = FALSE) +
      ggplot2::scale_alpha_manual(values = c(1),
         labels = c("RMTLR"),
         name = "Algorithm",
         guide = FALSE) +
      ggplot2::scale_x_discrete(labels = c(
        "overall_mean_single" = "Ensemble",
        "ImmuneCells" = "Cell fractions",
        "Pathways.cor" = "Pathways",
        "TFs" = "TFs",
        "LRpairs.spec.pc" = "LR pairs",
        "CCpairsGroupedScores.spec.pc" = "CC pairs",
        "ImmuneCells_CCpairsGroupedScores.spec.pc" = "Cell fractions + CC pairs",
        "ImmuneCells_LRpairs.spec.pc" = "Cell fractions + LR pairs",
        "ImmuneCells_TFs" = "Cell fractions + TFs",
        "Pathways.cor_CCpairsGroupedScores.spec.pc" = "Pathways + CC pairs",
        "Pathways.cor_ImmuneCells" = "Pathways + Cell fractions",
        "Pathways.cor_LRpairs.spec.pc" = "Pathways + LR pairs",
        "Pathways.cor_TFs" = "Pathways + TFs",
        "TFs_CCpairsGroupedScores.spec.pc" = "TFs + CC pairs",
        "TFs_LRpairs.spec.pc" = "TFs + LR pairs",
        "Gold_Standard" = "Tasks (gold standard)",
        "Pathways_ImmuneCells_TFs_LRpairs" = "Combo\nall views"
        )) +
      ggplot2::theme(panel.grid = ggplot2::element_blank(), panel.background = ggplot2::element_rect(fill = NA)) +
      ggplot2::theme_bw() +
      ggplot2::ylim(0, 1) +
      ggplot2::ylab("Area under the curve (AUC)") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = round(AUC_median, 2) - AUC_sd, ymax = round(AUC_median, 2) + AUC_sd), width = .3, color="black", position = ggplot2::position_dodge(0.9)) +
      ggplot2::geom_text(ggplot2::aes(label= round(AUC_median, 2)), stat = "identity", color="black", size = 4, angle = 90, hjust = -0.5, position = ggplot2::position_dodge(0.9)) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(size=12, angle = 45, vjust = 1, hjust = 1, color = "black"),
        axis.text.y = ggplot2::element_text(size = 12, color = "black"),
        axis.title.y = ggplot2::element_text(size = 12), axis.title.x = ggplot2::element_blank(),
        legend.position = "right", legend.text = ggplot2::element_text(size = 10),
        legend.box.background = ggplot2::element_rect(color = "black", size = 0.3),
        legend.box.margin = ggplot2::margin(0.5, 0.5, 0.5, 0.5)
      ) +
      ggplot2::labs(title = paste0("n=", length(real_patient_response), " (R=", n_R, "; ", "NR=", n_NR, ")"))

    ggplot2::ggsave(file.path(output_file_path, "plot_auc_barplot_values.pdf"), width = 10, height = 8)

    # *******************************************
    # Plot ROC curve

    all_color_views <- c("#6CD8CB", "#CC6AF2", "#01AD3F", "salmon", "#9E1D3F")
    names(all_color_views) <- c("Pathways.cor", "ImmuneCells", "TFs", "LRpairs.spec.pc", "CCpairsGroupedScores.spec.pc")

    pdf(file.path(output_file_path, "plot_roc_curves.pdf"), width = 10, height = 10)
    par(cex.axis = 1.6, mar = c(5, 5, 5, 5), col.lab = "black")

    # Single views
    ROCR::plot(ROC_info$Pathways.cor$RMTLR$common_mean[["1se.mse"]]$Curve[[1]],
      avg = "threshold", col = all_color_views["Pathways.cor"], lwd = 3, type = "S",
      cex.lab = 1.6, ylab = "True Positive Rate", xlab = "False Positive Rate"
    )

    ROCR::plot(ROC_info$ImmuneCells$RMTLR$common_mean[["1se.mse"]]$Curve[[1]],
      avg = "threshold", col = all_color_views["ImmuneCells"], lwd = 3, type = "S",
      cex.lab = 1.6, ylab = "True Positive Rate", xlab = "False Positive Rate", add = TRUE
    )

    ROCR::plot(ROC_info$TFs$RMTLR$common_mean[["1se.mse"]]$Curve[[1]],
      avg = "threshold", col = all_color_views["TFs"], lwd = 3, type = "S",
      cex.lab = 1.6, ylab = "True Positive Rate", xlab = "False Positive Rate", add = TRUE
    )

    ROCR::plot(ROC_info$LRpairs.spec.pc$RMTLR$common_mean[["1se.mse"]]$Curve[[1]],
      avg = "threshold", col = all_color_views["LRpairs.spec.pc"], lwd = 3, type = "S",
      cex.lab = 1.6, ylab = "True Positive Rate", xlab = "False Positive Rate", add = TRUE
    )

    ROCR::plot(ROC_info$CCpairsGroupedScores.spec.pc$RMTLR$common_mean[["1se.mse"]]$Curve[[1]],
      avg = "threshold", col = all_color_views["CCpairsGroupedScores.spec.pc"], lwd = 3, type = "S",
      cex.lab = 1.6, ylab = "True Positive Rate", xlab = "False Positive Rate", add = TRUE
    )

    # overall views
    ROCR::plot(ROC_info$overall_mean$RMTLR$overall[["1se.mse"]]$Curve[[1]],
      avg = "threshold", col = color_overalls[1], lwd = 3, type = "S",
      lty = 1, add = TRUE)

    # gold standard
    pred_GS <- ROCR::prediction(all_tasks_values_norm, labels[,1:9], label.ordering = c("NR", "R"))
    perf_GS  <- ROCR::performance(pred_GS,"tpr","fpr")
    AUC_GS <- unlist(ROCR::performance(pred_GS, "auc")@y.values)

    data_ROC_GS <- list(pred_GS)
    Barplot_GS <- list(AUC_GS)
    ROC_info_GS <- list(Curve = data_ROC_GS, Barplot = Barplot_GS)

    ROCR::plot(ROC_info_GS$Curve[[1]],
      avg = "threshold", col = color_gold_standard[1], lwd = 3, type = "S",
      lty = 1, add = TRUE)

    # abline(a=0, b=1, lty = 3, lwd = 2, col = "antiquewhite4")
    legend(
      x = 0.72, y = 0.57,
      legend = c(
        paste0("Pathways", " (", round(subset(AUC_mean_sd_RMTLR, View == "Pathways.cor" & Alg == "RMTLR")$AUC_median, 2), ")"),
        paste0("Cell fractions", " (", round(subset(AUC_mean_sd_RMTLR, View == "ImmuneCells" & Alg == "RMTLR")$AUC_median, 2), ")"),
        paste0("TFs", " (", round(subset(AUC_mean_sd_RMTLR, View == "TFs" & Alg == "RMTLR")$AUC_median, 2), ")"),
        paste0("LR pairs", " (", round(subset(AUC_mean_sd_RMTLR, View == "LRpairs.spec.pc" & Alg == "RMTLR")$AUC_median, 2), ")"),
        paste0("CC pairs", " (", round(subset(AUC_mean_sd_RMTLR, View == "CCpairsGroupedScores.spec.pc" & Alg == "RMTLR")$AUC_median, 2), ")"),
        paste0("\nEnsemble", " (", round(subset(AUC_mean_sd_RMTLR, View == "overall_mean_single")$AUC_median, 2), ")"),
        paste0("\nTasks", " (", round(subset(AUC_mean_sd_RMTLR, View == "Gold_Standard")$AUC_median, 2), ")")

      ),

      col = c(
        all_color_views[c("Pathways.cor", "ImmuneCells", "TFs", "LRpairs.spec.pc", "CCpairsGroupedScores.spec.pc")],
        as.vector(color_overalls)[1], as.vector(color_gold_standard)[1]
      ), lty = 1, lwd = 3, cex = 0.9, bty = "n"
    )
    title(main = paste0("n=", length(real_patient_response), " (R=", n_R, "; ", "NR=", n_NR, ")"))

    dev.off()
  } else {
    stop("No patients' response provided")
  }
}
