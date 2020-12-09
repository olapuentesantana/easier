#' Comparison of the actual predictions based on different metrics with real patient response.
#'
#' Plot ROC curve and barplot with the area under the ROC curve.
#' If list of gold standards is not provided, the function used a default one.
#'
#' \code{compare_response} plots ROC curves and barplots showing the accuracy of the predictions
#'  on the real patient response data. It receives as input the predicted immune response as well
#'  as the real patient response (they should be provided with the same order of samples). Gold
#'  standards are also plotted to consider them as reference. To compute these gold standards,
#'  transcriptomics data in different formats is provided. An output file path should be given
#'  to save the plots as pdfs. An additional input must be a list with the gold standards the
#'  user wish to compare the prediction results.
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
#' @param RNA.tpm numeric matrix with data
#' @param output_file_path string with a file name
#' @param list_gold_standards string with gold standards names
#' @param cancertype specify cancer type
#'
#' @return ROC curves plots and barplots showing AUC values.
#'
compare_immune_response <- function(predictions_immune_response = NULL,
                                    real_patient_response,
                                    RNA.tpm,
                                    output_file_path,
                                    list_gold_standards,
                                    cancertype){

  try(if(missing(cancertype)) stop("cancer type needs to be specified"))
  try(if(is.null(predictions_immune_response)) stop("none predictions found"))

  # Check that folder exists, create folder otherwise
  if(dir.exists(output_file_path) == F) {
    dir.create(file.path(output_file_path), showWarnings = FALSE)
    warning(paste0(sapply(strsplit(output_file_path, "/", fixed = T), tail , 1),
                   " folder does not exist, creating ", sapply(strsplit(output_file_path, "/", fixed = T), tail , 1), " folder"))
  }

  # Initialize variables
  AUC.median <- AUC.sd <- Model <- Task <- View <- NULL
  view_combinations <- names(predictions_immune_response)
  algorithms <- names(predictions_immune_response[[1]])
  tasks <- c(names(predictions_immune_response[[1]][[1]]), "common_mean", "common_median")
  models <- names(predictions_immune_response[[1]][[1]][[1]])
  #sel_views <- names(predictions_immune_response)

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
  all_color_views <- c("#52ac68","#6a70d7","#bbb442",
                      "#5b3788","#72a646","#c972c4",
                      "#43c8ac","#d25543","#6d8dd7",
                      "#cd8232","#b1457b","#9b843b",
                      "#ba4758","#a24e2e")
  names(all_color_views) <- view_combinations

  # all_color_views <- c("#6CD8CB","#CC6AF2","salmon","#01AD3F","#9E1D3F")
  # names(all_color_views) <- c("Pathways","ImmuneCells", "TFs", "LRpairs.spec.pc","CCpairsGroupedScores.spec.pc")


  # Selected views
  # sel_color_views <- all_color_views[names(all_color_views) %in% sel_views]

  # Tasks
  color_tasks <- vector("character", length = 14)
  color_tasks <- toupper(c("#7763ce", "#69b444", "#c754bf", "#c7ad3e", "#6484c8", "#d35238", "#4ac0cd",
                            "#d74278", "#4fa471", "#9b4d8a", "#7e843b", "#d48dca", "#c07c42", "#bd5f68"))
  names(color_tasks) <- c("CYT", "Ock_IS", "IPS", "IMPRES", "Roh_IS", "chemokines", "Davoli_IS", "IFNy", "Ayers_expIS",
                          "Tcell_inflamed", "TIDE", "MSI", "RIR", "TLS")

  # Patient response labels
  labels <- matrix(real_patient_response, nrow = length(real_patient_response), ncol = 100,
                   dimnames = list(colnames(RNA.tpm), seq(1, 100, 1)))

  # AUC predictions (when response available)

  if (missing(real_patient_response) == F){

    try(if(all(levels(as.factor(real_patient_response)) %in% c("NR", "R")) == F) {
      stop("real_patient_response factor levels are not NR and R") })

    # Compute gold standards
    default_list_gold_standards <- c(names(color_tasks), "CTLA4", "PD1", "PDL1")
    if(missing(list_gold_standards)){list_gold_standards = default_list_gold_standards}
    gold.standards <- compute_gold_standards(RNA.tpm, list_gold_standards, cancertype, output_file_path)

    # *******************************
    # Assess correlation between chemokines and the other correlated tasks
    # correlated tasks
    cor_tasks <- names(color_tasks)[!names(color_tasks) %in% c("IPS", "IMPRES", "TIDE", "MSI")]
    cor_tasks <- cor_tasks[!cor_tasks %in% c("Ock_IS")] # Unfeasible computation
    tasks_values <- do.call(cbind, lapply(cor_tasks, function(X){
      tmp <- as.numeric(unlist(gold.standards[[X]]))
    }))
    colnames(tasks_values) <- cor_tasks
    rownames(tasks_values) <- colnames(gold.standards$CYT[[1]])
    tasks_cormatrix <- cor(tasks_values)
    cor_sign <- sign(tasks_cormatrix[,"chemokines"])
    cor_sign <- cor_sign[names(cor_sign) != "chemokines"]
    if (all(cor_sign == -1)){
      tasks_values[,"chemokines"] <- -tasks_values[,"chemokines"]
    }
    tasks_values <- as.data.frame(tasks_values)

    # Overall mean and median tasks measure
    tasks_values$consensus_mean <- tasks_values$consensus_median
    tasks_values$consensus_median <- apply(tasks_values, 1, median)
    tasks_values$consensus_mean <- apply(tasks_values, 1, mean)

    gold.standards_unscaled <- sapply(colnames(tasks_values), function(X){
      tmp <- matrix(as.numeric(tasks_values[[X]]), nrow = nrow(tasks_values), ncol = 100,
                    dimnames = list(colnames(tasks_values[[X]])))
      return(list(tmp))
    })

    # *******************************
    # Scale tasks (according to our analysis, scaling does not affect prediction)

    # tasks_values_scaled <- data.frame(standarization(tasks_values))
    # rownames(tasks_values_scaled) <- colnames(RNA.tpm)
    # colnames(tasks_values_scaled) <- cor_tasks

    # tasks_values_scaled$consensus_mean <- tasks_values_scaled$consensus_median
    # tasks_values_scaled$consensus_median <- apply(tasks_values_scaled, 1, median)
    # tasks_values_scaled$consensus_mean <- apply(tasks_values_scaled, 1, mean)
    #
    # gold.standards_scaled <- sapply(colnames(tasks_values_scaled), function(X){
    #   tmp <- matrix(as.numeric(tasks_values_scaled[[X]]), nrow = nrow(tasks_values_scaled), ncol = 100,
    #                 dimnames = list(colnames(tasks_values_scaled[[X]])))
    #   return(list(tmp))
    # })

    # *******************************************
    # Comparison with transcriptomics data (only us)

    # transcript <- t(RNA.tpm)
    # view_info <- c(transcript = 'gaussian')
    # view_name <- paste(names(view_info), collapse="_")
    # view_data <- lapply(tolower(names(view_info)), function(x) as.data.frame(get(x)))
    # names(view_data) <- names(view_info)
    #
    # # Read tcga model
    # file <- dir(path = paste0("output/PanCancer_draft_v1/", cancertype, "/group_cor_tasks"),
    #             pattern = "all_cv_res_", full.names = T, recursive = F)
    # which_file <- grep(pattern = paste0("_with_cor_tasks_Transcript",".RData"), file, fixed = T)
    # load(file[which_file])
    #
    # # Predict immune response using model parameters
    # pred_alg_EN <- easier::predict_with_multitaskelasticnet(view_name = view_name,
    #                                                         view_info = view_info,
    #                                                         view_data = view_data,
    #                                                         learned_model = all_cv_res[["Multi_Task_EN"]])
    #
    # pred_alg_BE <- easier::predict_with_bemkl(view_name = view_name,
    #                                           view_info = view_info,
    #                                           view_data = view_data,
    #                                           learned_model = all_cv_res[["BEMKL"]])
    #
    # predictions_immune_response$transcript <- list(Multi_Task_EN = pred_alg_EN, BEMKL = pred_alg_BE)

    # *******************************************
    # Predictions #
    ROC_info <- lapply(c(view_combinations), function(view){
      ROC_info <- lapply(algorithms, function(alg){
        ROC_info <- lapply(c(names(predictions_immune_response[[view]][[alg]]), "common_mean", "common_median"), function(task){
          ROC_info <- lapply(models, function(model){

            if (alg != "BEMKL"){
              if (task %in% "common_mean"){

                df <- list(matrix(0,nrow(labels), ncol = 100))

                for (run in 1:ncol(labels)) {

                  df[[1]][,run] <- apply(rbind(predictions_immune_response[[view]][[alg]][["CYT"]][[model]][[1]][,run],
                                               predictions_immune_response[[view]][[alg]][["IS"]][[model]][[1]][,run],
                                               predictions_immune_response[[view]][[alg]][["RohIS"]][[model]][[1]][,run],
                                               predictions_immune_response[[view]][[alg]][["chemokine"]][[model]][[1]][,run],
                                               predictions_immune_response[[view]][[alg]][["IS_Davoli"]][[model]][[1]][,run],
                                               predictions_immune_response[[view]][[alg]][["IFny"]][[model]][[1]][,run],
                                               predictions_immune_response[[view]][[alg]][["ExpandedImmune"]][[model]][[1]][,run],
                                               predictions_immune_response[[view]][[alg]][["T_cell_inflamed"]][[model]][[1]][,run],
                                               predictions_immune_response[[view]][[alg]][["resF.down"]][[model]][[1]][,run],
                                               predictions_immune_response[[view]][[alg]][["TLS"]][[model]][[1]][,run]), 2, mean)
                }

              }else if(task %in% "common_median") {

                  df <- list(matrix(0,nrow(labels), ncol = 100))

                  for (run in 1:ncol(labels)) {

                      df[[1]][,run] <- apply(rbind(predictions_immune_response[[view]][[alg]][["CYT"]][[model]][[1]][,run],
                                                   predictions_immune_response[[view]][[alg]][["IS"]][[model]][[1]][,run],
                                                   predictions_immune_response[[view]][[alg]][["RohIS"]][[model]][[1]][,run],
                                                   predictions_immune_response[[view]][[alg]][["chemokine"]][[model]][[1]][,run],
                                                   predictions_immune_response[[view]][[alg]][["IS_Davoli"]][[model]][[1]][,run],
                                                   predictions_immune_response[[view]][[alg]][["IFny"]][[model]][[1]][,run],
                                                   predictions_immune_response[[view]][[alg]][["ExpandedImmune"]][[model]][[1]][,run],
                                                   predictions_immune_response[[view]][[alg]][["T_cell_inflamed"]][[model]][[1]][,run],
                                                   predictions_immune_response[[view]][[alg]][["resF.down"]][[model]][[1]][,run],
                                                   predictions_immune_response[[view]][[alg]][["TLS"]][[model]][[1]][,run]), 2, median)
                  }

              }else{

                df <- predictions_immune_response[[view]][[alg]][[task]][[model]]

              }

            }else{

              if (task %in% "common_mean"){

                df <- list(matrix(0, nrow(labels), ncol = 100))

                for (run in 1:ncol(labels)) {

                  df[[1]][,run] <- apply(rbind(predictions_immune_response[[view]][[alg]][["CYT"]][[1]][,run],
                                               predictions_immune_response[[view]][[alg]][["IS"]][[1]][,run],
                                               predictions_immune_response[[view]][[alg]][["RohIS"]][[1]][,run],
                                               predictions_immune_response[[view]][[alg]][["chemokine"]][[1]][,run],
                                               predictions_immune_response[[view]][[alg]][["IS_Davoli"]][[1]][,run],
                                               predictions_immune_response[[view]][[alg]][["IFny"]][[1]][,run],
                                               predictions_immune_response[[view]][[alg]][["ExpandedImmune"]][[1]][,run],
                                               predictions_immune_response[[view]][[alg]][["T_cell_inflamed"]][[1]][,run],
                                               predictions_immune_response[[view]][[alg]][["resF.down"]][[1]][,run],
                                               predictions_immune_response[[view]][[alg]][["TLS"]][[1]][,run]), 2, mean)
                }

              }else if(task %in% "common_median") {

                df <- list(matrix(0,nrow(labels), ncol = 100))

                for (run in 1:ncol(labels)) {

                  df[[1]][,run] <- apply(rbind(predictions_immune_response[[view]][[alg]][["CYT"]][[1]][,run],
                                               predictions_immune_response[[view]][[alg]][["IS"]][[1]][,run],
                                               predictions_immune_response[[view]][[alg]][["RohIS"]][[1]][,run],
                                               predictions_immune_response[[view]][[alg]][["chemokine"]][[1]][,run],
                                               predictions_immune_response[[view]][[alg]][["IS_Davoli"]][[1]][,run],
                                               predictions_immune_response[[view]][[alg]][["IFny"]][[1]][,run],
                                               predictions_immune_response[[view]][[alg]][["ExpandedImmune"]][[1]][,run],
                                               predictions_immune_response[[view]][[alg]][["T_cell_inflamed"]][[1]][,run],
                                               predictions_immune_response[[view]][[alg]][["resF.down"]][[1]][,run],
                                               predictions_immune_response[[view]][[alg]][["TLS"]][[1]][,run]), 2, median)
                }

              }else{

                df <- predictions_immune_response[[view]][[alg]][[task]]

             }
            }

            pred <- ROCR::prediction(df[[1]], labels, label.ordering = c("NR", "R"))
            perf <- ROCR::performance(pred,"tpr","fpr")
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

    # save(ROC_info, file = "../../../Desktop/comb.Gide_Auslander_pre_na/comb.Gide_Auslander_ROC_predictions_count_varstab_tcga.RData")

    # Predictions #
    overall_ROC_info <- lapply(c("overall_median", "overall_mean"), function(main_view){
      overall_ROC_info <- lapply(algorithms, function(alg){
        overall_ROC_info <- lapply("overall", function(task){
          overall_ROC_info <- lapply(models, function(model){
            df_overall <- lapply(view_combinations, function(view){

              df <- matrix(0, nrow = nrow(labels), ncol = 100)

              if (alg != "BEMKL"){

                if (main_view == "overall_median"){

                  for (run in 1:ncol(labels)) {

                    df[,run] <- apply(rbind(predictions_immune_response[[view]][[alg]][["CYT"]][[model]][[1]][,run],
                                            predictions_immune_response[[view]][[alg]][["IS"]][[model]][[1]][,run],
                                            predictions_immune_response[[view]][[alg]][["RohIS"]][[model]][[1]][,run],
                                            predictions_immune_response[[view]][[alg]][["chemokine"]][[model]][[1]][,run],
                                            predictions_immune_response[[view]][[alg]][["IS_Davoli"]][[model]][[1]][,run],
                                            predictions_immune_response[[view]][[alg]][["IFny"]][[model]][[1]][,run],
                                            predictions_immune_response[[view]][[alg]][["ExpandedImmune"]][[model]][[1]][,run],
                                            predictions_immune_response[[view]][[alg]][["T_cell_inflamed"]][[model]][[1]][,run],
                                            predictions_immune_response[[view]][[alg]][["resF.down"]][[model]][[1]][,run],
                                            predictions_immune_response[[view]][[alg]][["TLS"]][[model]][[1]][,run]), 2, median)
                  }

                  df_overall <- apply(df, 1, median)

                }else if (main_view == "overall_mean"){

                  for (run in 1:ncol(labels)) {

                    df[,run] <- apply(rbind(predictions_immune_response[[view]][[alg]][["CYT"]][[model]][[1]][,run],
                                            predictions_immune_response[[view]][[alg]][["IS"]][[model]][[1]][,run],
                                            predictions_immune_response[[view]][[alg]][["RohIS"]][[model]][[1]][,run],
                                            predictions_immune_response[[view]][[alg]][["chemokine"]][[model]][[1]][,run],
                                            predictions_immune_response[[view]][[alg]][["IS_Davoli"]][[model]][[1]][,run],
                                            predictions_immune_response[[view]][[alg]][["IFny"]][[model]][[1]][,run],
                                            predictions_immune_response[[view]][[alg]][["ExpandedImmune"]][[model]][[1]][,run],
                                            predictions_immune_response[[view]][[alg]][["T_cell_inflamed"]][[model]][[1]][,run],
                                            predictions_immune_response[[view]][[alg]][["resF.down"]][[model]][[1]][,run],
                                            predictions_immune_response[[view]][[alg]][["TLS"]][[model]][[1]][,run]), 2, mean)
                  }
                  df_overall <- apply(df, 1, mean)
                }
              }else{

                if (main_view == "overall_median"){

                  for (run in 1:ncol(labels)) {

                    df[,run] <- apply(rbind(predictions_immune_response[[view]][[alg]][["CYT"]][[1]][,run],
                                            predictions_immune_response[[view]][[alg]][["IS"]][[1]][,run],
                                            predictions_immune_response[[view]][[alg]][["RohIS"]][[1]][,run],
                                            predictions_immune_response[[view]][[alg]][["chemokine"]][[1]][,run],
                                            predictions_immune_response[[view]][[alg]][["IS_Davoli"]][[1]][,run],
                                            predictions_immune_response[[view]][[alg]][["IFny"]][[1]][,run],
                                            predictions_immune_response[[view]][[alg]][["ExpandedImmune"]][[1]][,run],
                                            predictions_immune_response[[view]][[alg]][["T_cell_inflamed"]][[1]][,run],
                                            predictions_immune_response[[view]][[alg]][["resF.down"]][[1]][,run],
                                            predictions_immune_response[[view]][[alg]][["TLS"]][[1]][,run]), 2, median)
                  }

                  df_overall <- apply(df, 1, median)

                }else if (main_view == "overall_mean"){

                    for (run in 1:ncol(labels)) {

                      df[,run] <- apply(rbind(predictions_immune_response[[view]][[alg]][["CYT"]][[1]][,run],
                                              predictions_immune_response[[view]][[alg]][["IS"]][[1]][,run],
                                              predictions_immune_response[[view]][[alg]][["RohIS"]][[1]][,run],
                                              predictions_immune_response[[view]][[alg]][["chemokine"]][[1]][,run],
                                              predictions_immune_response[[view]][[alg]][["IS_Davoli"]][[1]][,run],
                                              predictions_immune_response[[view]][[alg]][["IFny"]][[1]][,run],
                                              predictions_immune_response[[view]][[alg]][["ExpandedImmune"]][[1]][,run],
                                              predictions_immune_response[[view]][[alg]][["T_cell_inflamed"]][[1]][,run],
                                              predictions_immune_response[[view]][[alg]][["resF.down"]][[1]][,run],
                                              predictions_immune_response[[view]][[alg]][["TLS"]][[1]][,run]), 2, mean)
                  }
                  df_overall <- apply(df, 1, mean)
                }
             }
             return(df_overall)
            })
            names(df_overall) <- view_combinations

            if (main_view == "overall_median"){

              overall_df <- apply(rbind(df_overall$Pathways.cor, df_overall$ImmuneCells, df_overall$LRpairs.spec.pc,
                                        df_overall$TFs, df_overall$CCpairsGroupedScores.spec.pc), 2, median)
            }else if (main_view == "overall_mean"){

              overall_df <- apply(rbind(df_overall$Pathways.cor, df_overall$ImmuneCells, df_overall$LRpairs.spec.pc,
                                        df_overall$TFs, df_overall$CCpairsGroupedScores.spec.pc), 2, mean)
            }

            pred <- ROCR::prediction(overall_df, labels[,1], label.ordering = c("NR", "R"))
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
    names(overall_ROC_info) <- c("overall_median", "overall_mean")

    # Combine views and overall
    ROC_info <- c(ROC_info, overall_ROC_info)

    # Gold standards #
    gold_standards <- gold.standards_unscaled

    ROC_info.GS <- lapply(names(gold_standards), function(GS){

      pred.GS <- ROCR::prediction(gold_standards[[GS]], labels, label.ordering = c("NR", "R"))
      perf.GS  <- ROCR::performance(pred.GS,"tpr","fpr")
      AUC.GS <- unlist(ROCR::performance(pred.GS, "auc")@y.values)

      data_ROC.GS <- list(perf.GS)
      Barplot.GS <- list(AUC.GS)

      return(list(Curve = data_ROC.GS, Barplot = Barplot.GS))
    })
    names(ROC_info.GS) <- paste0(names(gold_standards),".GS")

    # Predictions #
    AUC.mean.sd <-  do.call(rbind, lapply(names(ROC_info), function(view){
      AUC.mean.sd <-  do.call(rbind, lapply(names(ROC_info[[view]]), function(alg){
        AUC.mean.sd <-  do.call(rbind, lapply(names(ROC_info[[view]][[alg]]), function(task){
          AUC.mean.sd <-  do.call(rbind, lapply(names(ROC_info[[view]][[alg]][[task]]), function(model){

            View <- rep(view, times = 100)
            Alg <- rep(alg, times = 100)
            Model <- rep(model, times = 100)
            Task <- rep(task, times = 100)
            cv_iter <- seq(1,100)

            AUC.data <- data.frame(Alg = Alg,
                                   Model = Model,
                                   AUC = as.numeric(unlist(ROC_info[[view]][[alg]][[task]][[model]]$Barplot)),
                                   View = View,
                                   Task = Task,
                                   iteration = cv_iter)

            AUC.mean.sd <- do.call(data.frame,
                                   aggregate(AUC ~ Alg + Model + View + Task,
                                             data = AUC.data,
                                             FUN = function(x) c(median = median(x),
                                                                 sd = sd(x))))

            return(AUC.mean.sd)
          }))
        return(AUC.mean.sd)
        }))
      return(AUC.mean.sd)
      }))
      return(AUC.mean.sd)
    }))

    # Gold standards #
    AUC.mean.sd_GS <-  do.call(rbind, lapply(names(ROC_info.GS), function(view){

      View <- rep(sapply(strsplit(view, split = ".", fixed = T), head, 1), times = 100)
      Alg <- rep("any", times = 100)
      Model <- rep("1se.mse", times = 100)
      Task <- rep("Gold_Standard", times = 100)
      cv_iter <- seq(1,100)

      AUC.data <- data.frame(Alg = Alg,
                             Model = Model,
                             AUC = as.numeric(unlist(ROC_info.GS[[view]]$Barplot)),
                             View = View,
                             Task = Task,
                             iteration = cv_iter)

      AUC.mean.sd_GS <- do.call(data.frame,
                                aggregate(AUC ~ Alg + Model + View + Task,
                                          data = AUC.data, FUN = function(x) c(median = median(x),
                                                                               sd = sd(x))))

      return(AUC.mean.sd_GS)
    }))
    AUC.mean.sd <- rbind(AUC.mean.sd, AUC.mean.sd_GS)

    AUC.mean.sd$Task <- factor(AUC.mean.sd$Task, levels = unique(AUC.mean.sd$Task))
    AUC.mean.sd$View <- factor(AUC.mean.sd$View, levels = unique(AUC.mean.sd$View))

    # Keep more regularized model
    AUC.mean.sd_reg <- subset(AUC.mean.sd, Model == "1se.mse")

    # *******************************************
    # Symplify for us (subset) #
    keep <- c("overall_mean", "overall_median", "consensus_mean", "consensus_median")
    AUC.mean.sd_reg_a <- subset(AUC.mean.sd_reg, View %in% keep)
    AUC.mean.sd_reg_b <- subset(AUC.mean.sd_reg, Task %in% c("common_median", "common_mean"))

    AUC.mean.sd_reg_final <- rbind(AUC.mean.sd_reg_a, AUC.mean.sd_reg_b)

    AUC.mean.sd_reg_final$View <- factor(AUC.mean.sd_reg_final$View,
                                         levels = c("overall_median","overall_mean", "consensus_median", "consensus_mean",
                                                    names(all_color_views)))

    AUC.mean.sd_reg_final$Task <- factor(AUC.mean.sd_reg_final$Task,
                                         levels = c("overall", "Gold_Standard", "common_median", "common_mean"))

    color_consensus <- c("gray82", "gray57"); names(color_consensus) <- c("Overall (median) tasks", "Overall (mean) tasks")
    color_overalls <- c("gold1","gold3"); names(color_overalls) <- c("Overall (median) views", "Overall (mean) views")

    n_R <- table(real_patient_response)[["R"]]
    n_NR <- table(real_patient_response)[["NR"]]

    # *******************************************
    # Barplot AUC values

    # We go for common mean in Gide + Auslander (SD as R) #
    AUC.mean.sd_reg_final_sub <- subset(AUC.mean.sd_reg_final, Task %in% c("overall", "Gold_Standard", "common_mean")
                                        & !View %in% c("overall_median", "consensus_median"))


    ggplot2::ggplot(AUC.mean.sd_reg_final_sub, aes(x=View, y=round(AUC.median, 4), fill=View)) + # alpha = Task
      ggplot2::geom_bar(stat="identity", position = position_dodge(), color = "white") +
      ggplot2::scale_fill_manual(values = c(as.vector(color_overalls)[2], as.vector(color_consensus)[2],
                                            as.vector(all_color_views)), guide = F) +
      #ggplot2::scale_alpha_manual(values = c(1,1,0.6,0.8)) +
      ggplot2::scale_x_discrete(labels = c("overall_mean" = "Overall\n(mean) views", "overall_median" = "Overall\n(median) views",
                                           "consensus_median" = "Overall\n(median) tasks", "consensus_mean" = "Overall\n(mean) tasks",
                                           "ImmuneCells" = "Cell fractions",
                                           "Pathways.cor" = "Pathways",
                                           "TFs" = "TFs",
                                           "LRpairs.spec.pc" = "L-R pairs",
                                           "CCpairsGroupedScores.spec.pc" = "C-C pairs",
                                           "ImmuneCells_CCpairsGroupedScores.spec.pc" = "Cell fractions + C-C pairs",
                                           "ImmuneCells_LRpairs.spec.pc" = "Cell fractions + L-R pairs",
                                           "ImmuneCells_TFs" = "Cell fractions + TFs",
                                           "Pathways.cor_CCpairsGroupedScores.spec.pc" = "Pathways + C-C pairs",
                                           "Pathways.cor_ImmuneCells" = "Pathways + Cell fractions",
                                           "Pathways.cor_LRpairs.spec.pc" = "Pathways + L-R pairs",
                                           "Pathways.cor_TFs" = "Pathways + TFs",
                                           "TFs_CCpairsGroupedScores.spec.pc" = "TFs + C-C pairs",
                                           "TFs_LRpairs.spec.pc" = "TFs + L-R pairs",
                                           "Multi_Task_EN" = "Gold_Standard",
                                           "Pathways_ImmuneCells_TFs_LRpairs" = "Combo\nall views",
                                           "CTLA4" = "\nCTLA4", "PD1" = "\nPD1", "PDL1" = "\nPDL1")) +
      ggplot2::theme(panel.grid = element_blank(), panel.background = element_rect(fill = NA)) +
      ggplot2::theme_bw()+
      ggplot2::ylim(0, 1) +
      ggplot2::ylab("Area under the curve (AUC)") +
      ggplot2::geom_errorbar(aes(ymin = round(AUC.median,4) - AUC.sd, ymax = round(AUC.median,4) + AUC.sd), width = .5, color="black", position = position_dodge(0.9)) +
      ggplot2::geom_text(aes(label= round(AUC.median, 4)), stat = "identity", color="black", size = 3, angle = 90, hjust = -0.8, position = position_dodge(0.9)) +
      ggplot2::theme(axis.text.x = element_text(size=12, angle = 45, vjust = 0.7, hjust = 0.7, color = "black"),
                     axis.text.y = element_text(size=12, color = "black"),
                     axis.title.y = element_text(size=12), axis.title.x = element_blank(),
                     legend.position = "right", legend.text = element_text(size = 10),
                     legend.box.background = element_rect(color="black", size=0.3),
                     legend.box.margin = margin(0.5, 0.5, 0.5, 0.5)) +
    ggplot2::labs(title = paste0("n=", length(real_patient_response)," (R=",n_R,"; ","NR=",n_NR,")"))

    ggplot2::ggsave(paste0(output_file_path,"/AUC_barplot_all_views_and_overall.pdf"), width = 10, height = 10)

    # *******************************************
    # Plot ROC curve

    all_color_views <- c("#6CD8CB","#CC6AF2","#01AD3F","salmon","#9E1D3F")
    names(all_color_views) <- c("Pathways.cor","ImmuneCells", "TFs", "LRpairs.spec.pc","CCpairsGroupedScores.spec.pc")

    pdf(paste0(output_file_path,"/Final/ROC_curve_all_views_and_overall.pdf"), width = 10, height = 10)
    par(cex.axis=1.6, mar = c(5, 5, 5, 5), col.lab="black")

    # Single views
    ROCR::plot(ROC_info$Pathways.cor$Multi_Task_EN$common_mean[["1se.mse"]]$Curve[[1]],
               avg = "threshold", col = all_color_views["Pathways.cor"], lwd = 3, type = "S",
               cex.lab=1.6, ylab="True Positive Rate", xlab="False Positive Rate")

    ROCR::plot(ROC_info$ImmuneCells$Multi_Task_EN$common_mean[["1se.mse"]]$Curve[[1]],
               avg = "threshold", col = all_color_views["ImmuneCells"], lwd = 3, type = "S",
               cex.lab=1.6, ylab="True Positive Rate", xlab="False Positive Rate", add = TRUE)

    ROCR::plot(ROC_info$TFs$Multi_Task_EN$common_mean[["1se.mse"]]$Curve[[1]],
               avg = "threshold", col = all_color_views["TFs"], lwd = 3, type = "S",
               cex.lab=1.6, ylab="True Positive Rate", xlab="False Positive Rate", add = TRUE)

    ROCR::plot(ROC_info$LRpairs.spec.pc$Multi_Task_EN$common_mean[["1se.mse"]]$Curve[[1]],
               avg = "threshold", col = all_color_views["LRpairs.spec.pc"], lwd = 3, type = "S",
               cex.lab=1.6, ylab="True Positive Rate", xlab="False Positive Rate", add = TRUE)

    ROCR::plot(ROC_info$CCpairsGroupedScores.spec.pc$Multi_Task_EN$common_mean[["1se.mse"]]$Curve[[1]],
               avg = "threshold", col = all_color_views["CCpairsGroupedScores.spec.pc"], lwd = 3, type = "S",
               cex.lab=1.6, ylab="True Positive Rate", xlab="False Positive Rate", add = TRUE)


    # ROCR::plot(ROC_info$transcript$Multi_Task_EN$common_mean[["1se.mse"]]$Curve[[1]],
    #            avg = "threshold", col = "black", lwd = 2, type = "S",
    #            cex.lab=2.2, ylab="True Positive Rate", xlab="False Positive Rate", add = TRUE)

    # gold standard - overall mean tasks
    ROCR::plot(ROC_info.GS$consensus_mean.GS$Curve[[1]], col = color_consensus[2], lwd = 3, type = "S", lty = 1, add = TRUE)

    # gold standard - overall mean tasks
    ROCR::plot(ROC_info$overall_mean$Multi_Task_EN$overall[["1se.mse"]]$Curve[[1]], col = color_overalls[2], lwd = 3, type = "S", lty = 1, add = TRUE)

    #abline(a=0, b=1, lty = 3, lwd = 2, col = "antiquewhite4")
    legend(x = 0.72,y = 0.57,
           legend = c(paste0("Pathways","\n(AUC=", round(subset(AUC.mean.sd, Task == "common_mean"
                      & Model == "1se.mse" & View == "Pathways.cor" & Alg == "Multi_Task_EN")$AUC.median, 3),")"),
                      paste0("Cell fractions","\n(AUC=", round(subset(AUC.mean.sd, Task == "common_mean"
                      & Model == "1se.mse" & View == "ImmuneCells" & Alg == "Multi_Task_EN")$AUC.median, 3),")"),
                      paste0("TFs","\n(AUC=", round(subset(AUC.mean.sd, Task == "common_mean"
                      & Model == "1se.mse" & View == "TFs" & Alg == "Multi_Task_EN")$AUC.median, 3),")"),
                      paste0("L-R pairs","\n(AUC=", round(subset(AUC.mean.sd, Task == "common_mean"
                      & Model == "1se.mse" & View == "LRpairs.spec.pc" & Alg == "Multi_Task_EN")$AUC.median, 3),")"),
                      paste0("C-C pairs","\n(AUC=", round(subset(AUC.mean.sd, Task == "common_mean"
                      & Model == "1se.mse" & View == "CCpairsGroupedScores.spec.pc" & Alg == "Multi_Task_EN")$AUC.median, 3),")"),
                      # paste0("Transcriptomics","\n(AUC=", round(subset(AUC.mean.sd, Task == "common_mean"
                      # & Model == "1se.mse" & View == "transcript" & Alg == "Multi_Task_EN")$AUC.median, 3),")"),
                      paste0("\nOverall tasks","\n(AUC=", round(subset(AUC.mean.sd_GS, Model == "1se.mse" & View == "consensus_mean")$AUC.median, 3),")"),
                      paste0("\nOverall views","\n(AUC=", round(subset(AUC.mean.sd, Model == "1se.mse" & View == "overall_mean")$AUC.median, 3),")")),

           col = c(all_color_views[c("Pathways.cor", "ImmuneCells","TFs", "LRpairs.spec.pc", "CCpairsGroupedScores.spec.pc")],
                   as.vector(color_consensus)[2], as.vector(color_overalls)[2]), lty = 1, lwd =4, cex = 1.1, bty = "n")
    title(main = paste0("n=", length(real_patient_response)," (R=",n_R,"; ","NR=",n_NR,")"))

    dev.off()

  }else{
    stop("No patients'response provided")
  }
}


