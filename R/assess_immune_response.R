#' Assess easier score as predictor of patients'
#' immune response
#'
#' Evaluates the predictive performance of easier score as
#' predictor of patients' immune response. Provided that patients'
#' response is available, a roc curve plot and a barplot showing the
#' average (across tasks) area under the ROC curve (AUC) values for each
#' quantitative descriptor. An ensemble model is considered based on
#' the average of the individual descriptors. Also, the average of the
#' gold standard scores is shown. When available information on tumor
#' mutational burden (TMB), this is also displayed.
#'
#' Since both immune response and TMB are essential for effective
#' immunotherapy response, an integrated score is provided given
#' two different approaches based on a applying either a weighted
#' average or penalty to patients' easier score depending on their
#' TMB category.
#'
#' @importFrom ROCR prediction performance plot
#' @importFrom grDevices recordPlot
#' @importFrom stats aggregate median sd
#' @importFrom graphics legend par title segments lines
#' @import ggplot2
#'
#' @export
#'
#' @param predictions_immune_response list containing the predictions
#' for each quantitative descriptor and for each task.
#' @param patient_response character vector with two factors
#' (Non-responders = NR, Responders = R).
#' @param RNA_tpm numeric matrix of patients' gene expression data as
#' tpm values.
#' @param TMB_values numeric vector containing patients' tumor mutational
#' burden (TMB) values.
#' @param easier_with_TMB character string indicating which approach
#' should be used to integrate easier with TMB. One of "weighted_average"
#' (default) or "penalized_score".
#' @param weight_penalty integer value from 0 to 1, which is used to
#' define the weight or penalty for combining easier and TMB scores based
#' on a weighted average or penalized score, in order to derive a score of
#' patient's likelihood of immune response. The default value is 0.5.
#' @param verbose logical flag indicating whether to display messages
#' about the process.
#'
#' @return When patient_response is provided, a roc curve plot and a barplot
#' displaying the average (across tasks) area under the ROC curve (AUC)
#' values is provided. When patient_response is not provided, the easier
#' score is represented as boxplots (10 tasks) for each patient.
#'
#' When patient_response is provided and easier_with_TMB = weighted_average
#' or easier_with_TMB = penalized_score, an scatterplot displays the AUC
#' values of the integrated approach, easier score and TMB. In this case,
#' when patient_response is not provided, the integrated score is represented
#' as a dotplot for each patient.
#'
#' @examples
#' # using a SummarizedExperiment object
#' library(SummarizedExperiment)
#' # Using example exemplary dataset (Mariathasan et al., Nature, 2018)
#' # from easierData. Original processed data is available from
#' # IMvigor210CoreBiologies package.
#' library("easierData")
#'
#' dataset_mariathasan <- easierData::get_Mariathasan2018_PDL1_treatment()
#' RNA_tpm <- assays(dataset_mariathasan)[["tpm"]]
#' RNA_counts <- assays(dataset_mariathasan)[["counts"]]
#' cancer_type <- metadata(dataset_mariathasan)[["cancertype"]]
#'
#' # Select a subset of patients to reduce vignette building time.
#' pat_subset <- c(
#'     "SAM76a431ba6ce1", "SAMd3bd67996035", "SAMd3601288319e",
#'     "SAMba1a34b5a060", "SAM18a4dabbc557"
#' )
#' RNA_counts <- RNA_counts[, colnames(RNA_counts) %in% pat_subset]
#' RNA_tpm <- RNA_tpm[, colnames(RNA_tpm) %in% pat_subset]
#'
#' # Computation of cell fractions  (Finotello et al., Genome Med, 2019)
#' cell_fractions <- compute_cell_fractions(RNA_tpm = RNA_tpm)
#'
#' # Computation of pathway scores (Holland et al., BBAGRM, 2019;
#' # Schubert et al., Nat Commun, 2018)
#' pathway_activities <- compute_pathway_activity(
#'     RNA_counts = RNA_counts,
#'     remove_sig_genes_immune_response = TRUE
#' )
#'
#' # Computation of TF activity (Garcia-Alonso et al., Genome Res, 2019)
#' tf_activities <- compute_TF_activity(
#'     RNA_tpm = RNA_tpm
#' )
#'
#' # Computation of ligand-receptor pair weights
#' lrpair_weights <- compute_LR_pairs(
#'     RNA_tpm = RNA_tpm,
#'     cancer_type = "pancan"
#' )
#'
#' # Computation of cell-cell interaction scores
#' ccpair_scores <- compute_CC_pairs(
#'     lrpairs = lrpair_weights,
#'     cancer_type = "pancan"
#' )
#'
#' # Predict patients' immune response
#' predictions <- predict_immune_response(
#'     pathways = pathway_activities,
#'     immunecells = cell_fractions,
#'     tfs = tf_activities,
#'     lrpairs = lrpair_weights,
#'     ccpairs = ccpair_scores,
#'     cancer_type = cancer_type,
#'     verbose = TRUE
#' )
#'
#' # retrieve clinical response
#' patient_ICBresponse <- colData(dataset_mariathasan)[["BOR"]]
#' names(patient_ICBresponse) <- colData(dataset_mariathasan)[["pat_id"]]
#'
#' # retrieve TMB
#' TMB <- colData(dataset_mariathasan)[["TMB"]]
#' names(TMB) <- colData(dataset_mariathasan)[["pat_id"]]
#'
#' patient_ICBresponse <- patient_ICBresponse[names(patient_ICBresponse) %in% pat_subset]
#' TMB <- TMB[names(TMB) %in% pat_subset]
#'
#' # Assess patient-specific likelihood of response to ICB therapy
#' output_eval_with_resp <- assess_immune_response(
#'     predictions_immune_response = predictions,
#'     patient_response = patient_ICBresponse,
#'     RNA_tpm = RNA_tpm,
#'     TMB_values = TMB,
#'     easier_with_TMB = "weighted_average",
#' )
assess_immune_response <- function(predictions_immune_response = NULL,
                                   patient_response = NULL,
                                   RNA_tpm,
                                   TMB_values,
                                   easier_with_TMB = c("weighted_average", "penalized_score"),
                                   weight_penalty,
                                   verbose = TRUE) {
    if (is.null(predictions_immune_response)) stop("None predictions found")
    if (length(names(table(patient_response))) == 1) {
        stop("A binary class is required to evaluate patient's response")
    }
    if (length(easier_with_TMB) > 1) {
        stop("Please provide only one approach to combine easier with TMB")
    }
    if (missing(TMB_values)) {
        TMB_available <- FALSE
        easier_with_TMB <- "none"
    } else {
        TMB_available <- TRUE
        if (missing(easier_with_TMB)) easier_with_TMB <- "weighted_average"
        if (anyNA(TMB_values)) {
            warning(
                "NA values were found in TMB data, ",
                "patients with NA values are removed from the analysis"
            )
        }
        message(
            "\nConsidering ", length(TMB_values[!is.na(TMB_values)]),
            " patients out of ", length(TMB_values), " with available TMB"
        )
        patients_to_keep <- names(TMB_values[!is.na(TMB_values)])
        TMB_values <- as.numeric(TMB_values[patients_to_keep])
        names(TMB_values) <- patients_to_keep
        if (is.null(patient_response) == FALSE) {
            patient_response <- patient_response[patients_to_keep]
        }
        RNA_tpm <- RNA_tpm[, patients_to_keep]
    }
    # Initialize function
    if (easier_with_TMB == "weighted_average") {
        linear_func <- function(x) {
            min_x <- min(x)
            max_x <- max(x)
            x01 <- (x - min_x) / (max_x - min_x)
            return(x01)
        }
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
    all_color_views <- all_color_views[views]
    # ---------------------------#
    # Obtain Ensemble view predictor #
    # ---------------------------#
    #  Collect views data, and compute ensemble
    ensemble_df <- lapply(views, function(spec_view) {
        ensemble_df <- vapply(tasks, function(spec_task) {
            df <- predictions_immune_response[[spec_view]][[spec_task]]
            if (TMB_available) df <- df[patients_to_keep, ]
            df_runs <- rowMeans(df)
        }, FUN.VALUE = numeric(ncol(RNA_tpm)))
        return(ensemble_df)
    })
    names(ensemble_df) <- views
    overall_df <- vapply(tasks, function(spec_task) {
        overall_df <- apply(cbind(
            ensemble_df$pathways[, spec_task],
            ensemble_df$immunecells[, spec_task],
            ensemble_df$tfs[, spec_task],
            ensemble_df$lrpairs[, spec_task],
            ensemble_df$ccpairs[, spec_task]
        ), 1, mean)
    }, FUN.VALUE = numeric(ncol(RNA_tpm)))
    # AUC predictions (when response available)
    if (is.null(patient_response) == FALSE) {
        if (all(levels(as.factor(patient_response)) %in% c("NR", "R")) == FALSE) {
            stop("patient_response factor levels are not NR and R")
        }
        # Patient response labels
        labels <- matrix(patient_response,
            nrow = length(patient_response), ncol = length(tasks),
            dimnames = list(colnames(RNA_tpm), tasks)
        )
        # Compute scores of immune response and consider them as gold standards
        tasks_values <- compute_scores_immune_response(RNA_tpm)
        if (verbose) message("Scores of immune response computed!")
        # Assess correlation between chemokines and the other correlated tasks
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
            ROC_pred <- vapply(tasks, function(spec_task) {
                df <- predictions_immune_response[[spec_view]][[spec_task]]
                if (TMB_available) df <- df[patients_to_keep, ]
                # check patients match
                df <- df[match(rownames(labels), rownames(df)), ]
                df_runs <- rowMeans(df)
            }, FUN.VALUE = numeric(ncol(RNA_tpm)))
            pred <- ROCR::prediction(ROC_pred, labels, label.ordering = c("NR", "R"))
            perf <- ROCR::performance(pred, "tpr", "fpr")
            AUC <- unlist(ROCR::performance(pred, "auc")@y.values)
            names(AUC) <- tasks
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
            ))
        }))
        # Average across tasks
        AUC_mean_sd_run_tasks <- do.call(
            data.frame,
            stats::aggregate(AUC ~ View,
                data = AUC_data,
                FUN = function(x) c(mean = mean(x), sd = stats::sd(x))
            )
        )
        # ---------------------------#
        # Predictions ensemble view #
        # ---------------------------#
        # check patients match
        overall_df <- overall_df[match(rownames(labels), rownames(overall_df)), ]
        pred <- ROCR::prediction(overall_df, labels, label.ordering = c("NR", "R"))
        perf <- ROCR::performance(pred, "tpr", "fpr")
        AUC <- unlist(ROCR::performance(pred, "auc")@y.values)
        names(AUC) <- tasks
        ensemble_ROC_pred <- list(ensemble = list(Curve = list(perf), Barplot = list(AUC)))
        # Collect predictions into data.frame #
        AUC_data_ensemble <- data.frame(
            AUC = as.numeric(unlist(ensemble_ROC_pred$ensemble$Barplot)),
            View = "ensemble",
            Task = names(unlist(ensemble_ROC_pred$ensemble$Barplot))
        )
        # Average across tasks
        AUC_mean_sd_ensemble_run_tasks <- do.call(
            data.frame,
            stats::aggregate(AUC ~ View,
                data = AUC_data_ensemble,
                FUN = function(x) c(mean = mean(x), sd = stats::sd(x))
            )
        )
        # ---------------------------#
        # Predictions gold standards #
        # ---------------------------#
        # check patients match
        tasks_values <- tasks_values[match(rownames(labels), rownames(tasks_values)), ]
        pred <- ROCR::prediction(tasks_values, labels[, colnames(tasks_values)],
            label.ordering = c("NR", "R")
        )
        perf <- ROCR::performance(pred, "tpr", "fpr")
        AUC <- unlist(ROCR::performance(pred, "auc")@y.values)
        names(AUC) <- colnames(tasks_values)
        goldstandard_ROC_pred <- list(gold_standard = list(
            Curve = list(perf),
            Barplot = list(AUC)
        ))

        AUC_data_goldstandard <- data.frame(
            AUC = as.numeric(unlist(goldstandard_ROC_pred$gold_standard$Barplot)),
            View = "gold_standard",
            Task = names(unlist(goldstandard_ROC_pred$gold_standard$Barplot))
        )
        # Average across tasks
        AUC_mean_sd_goldstandard_run_tasks <- do.call(
            data.frame,
            stats::aggregate(AUC ~ View,
                data = AUC_data_goldstandard,
                FUN = function(x) c(mean = mean(x), sd = stats::sd(x))
            )
        )
        # Combine single views, ensemble and gold standard
        AUC_mean_sd_all_run_tasks <- rbind(
            AUC_mean_sd_run_tasks,
            AUC_mean_sd_ensemble_run_tasks,
            AUC_mean_sd_goldstandard_run_tasks
        )
        ROC_all_run_tasks <- c(ROC_pred, ensemble_ROC_pred, goldstandard_ROC_pred)
        names(ROC_all_run_tasks) <- c(
            names(all_color_views), "ensemble",
            "gold_standard"
        )
        AUC_mean_sd_all_run_tasks$View <- factor(AUC_mean_sd_all_run_tasks$View,
            levels = c(names(all_color_views), "ensemble", "gold_standard")
        )

        # Colors gold standard and ensemble
        color_gold_standard <- "gray82"
        names(color_gold_standard) <- "gold_standard"
        color_ensemble <- "gold2"
        names(color_ensemble) <- "ensemble"
        # ------------------------------------------------------#
        # Predictions TMB (if available)
        # ------------------------------------------------------#
        if (TMB_available == TRUE) {
            # check patients match
            TMB_values <- TMB_values[match(rownames(labels), names(TMB_values))]
            pred <- ROCR::prediction(TMB_values, labels[, 1], label.ordering = c("NR", "R"))
            perf <- ROCR::performance(pred, "tpr", "fpr")
            AUC <- unlist(ROCR::performance(pred, "auc")@y.values)
            TMB_ROC_pred <- list(TMB = list(Curve = list(perf), Barplot = list(AUC)))

            # Collect predictions into data.frame #
            AUC_data_TMB <- data.frame(
                AUC = as.numeric(unlist(TMB_ROC_pred$TMB$Barplot)),
                View = "TMB",
                Task = "TMB"
            )
            # Average across tasks
            AUC_mean_sd_TMB_run_tasks <- do.call(
                data.frame,
                stats::aggregate(AUC ~ View,
                    data = AUC_data_TMB,
                    FUN = function(x) c(mean = mean(x), sd = stats::sd(x))
                )
            )
            AUC_mean_sd_all_run_tasks <- rbind(
                AUC_mean_sd_all_run_tasks,
                AUC_mean_sd_TMB_run_tasks
            )
            ROC_all_run_tasks <- c(ROC_all_run_tasks, TMB_ROC_pred)
            names(ROC_all_run_tasks) <- c(
                names(all_color_views), "ensemble",
                "gold_standard", "TMB"
            )
            AUC_mean_sd_all_run_tasks$View <- factor(AUC_mean_sd_all_run_tasks$View,
                levels = c(names(all_color_views), "ensemble", "gold_standard", "TMB")
            )
            # Colors TMB
            color_TMB <- "salmon"
            names(color_TMB) <- "salmon"
        }
        plot_list <- list()
        # *******************************************
        # Figure: Barplot with AUC values
        # *******************************************
        n_R <- table(patient_response)[["R"]]
        n_NR <- table(patient_response)[["NR"]]
        gg <- ggplot2::ggplot(
            AUC_mean_sd_all_run_tasks,
            ggplot2::aes(
                x = .data$View,
                y = round(.data$AUC.mean, 2),
                fill = .data$View
            )
        ) +
            ggplot2::geom_bar(
                stat = "identity", position = ggplot2::position_dodge(),
                color = "white"
            ) +
            if (TMB_available) {
                ggplot2::scale_fill_manual(values = c(
                    as.vector(all_color_views), as.vector(color_ensemble),
                    as.vector(color_gold_standard), as.vector(color_TMB)
                ), guide = "none")
            } else {
                ggplot2::scale_fill_manual(values = c(
                    as.vector(all_color_views), as.vector(color_ensemble),
                    as.vector(color_gold_standard)
                ), guide = "none")
            }
        ggg <- gg + ggplot2::scale_x_discrete(labels = c(
            "ensemble" = "Ensemble",
            "immunecells" = "Cell fractions",
            "pathways" = "Pathways",
            "tfs" = "TFs",
            "lrpairs" = "LR pairs",
            "ccpairs" = "CC pairs",
            "gold_standard" = "Tasks (gold standard)"
        )) +
            ggplot2::theme(
                panel.grid = ggplot2::element_blank(),
                panel.background = ggplot2::element_rect(fill = NA)
            ) +
            ggplot2::theme_bw() +
            ggplot2::ylim(0, 1) +
            ggplot2::ylab("Area under the curve (AUC)") +
            ggplot2::geom_errorbar(ggplot2::aes(
                ymin = round(.data$AUC.mean, 2) - .data$AUC.sd,
                ymax = round(.data$AUC.mean, 2) + .data$AUC.sd
            ),
            width = .3, color = "black",
            position = ggplot2::position_dodge(0.9)
            ) +
            ggplot2::geom_text(ggplot2::aes(label = round(.data$AUC.mean, 2)),
                stat = "identity",
                color = "black", size = 4, angle = 90, hjust = -0.5,
                position = ggplot2::position_dodge(0.9)
            ) +
            ggplot2::theme(
                axis.text.x = ggplot2::element_text(
                    size = 12, angle = 45, vjust = 1, hjust = 1,
                    color = "black"
                ),
                axis.text.y = ggplot2::element_text(size = 12, color = "black"),
                axis.title.y = ggplot2::element_text(size = 12),
                axis.title.x = ggplot2::element_blank(),
                legend.position = "right", legend.text = ggplot2::element_text(size = 10),
                legend.box.background = ggplot2::element_rect(color = "black", size = 0.3),
                legend.box.margin = ggplot2::margin(0.5, 0.5, 0.5, 0.5)
            ) +
            ggplot2::labs(title = paste0(
                " n=", length(patient_response),
                " (R=", n_R, "; ", "NR=", n_NR, ")"
            ))

        plot_list[[1]] <- print(ggg)
        # *******************************************
        # Figure: ROC curves
        # *******************************************
        all_colors <- c(all_color_views, color_ensemble, color_gold_standard)
        graphics::par(
            cex.axis = 1.3, mar = c(5, 4, 2, 12), col.lab = "black",
            pty = "s", xpd = TRUE
        )
        # Single views & Gold Standard
        ROCR::plot(ROC_all_run_tasks[[1]]$Curve[[1]],
            avg = "threshold", col = all_colors[1], lwd = 2, type = "S",
            cex.lab = 1.3, ylab = "True Positive Rate", xlab = "False Positive Rate",
            bty = "L"
        )
        lapply(setdiff(
            names(ROC_all_run_tasks)[2:length(names(ROC_all_run_tasks))],
            "TMB"
        ), function(descriptor) {
            ROCR::plot(ROC_all_run_tasks[[descriptor]]$Curve[[1]],
                avg = "threshold", col = all_colors[descriptor], lwd = 2, type = "S",
                cex.lab = 1.3, ylab = "True Positive Rate", xlab = "False Positive Rate",
                add = TRUE
            )
        })
        legend_text <- do.call(
            c,
            lapply(setdiff(
                names(ROC_all_run_tasks)[seq_len(length(names(ROC_all_run_tasks)))],
                "TMB"
            ), function(descriptor) {
                paste0(descriptor, " (", round(subset(
                    AUC_mean_sd_all_run_tasks,
                    View == descriptor
                )$AUC.mean, 2), ")")
            })
        )
        # TMB
        if (TMB_available) {
            ROCR::plot(ROC_all_run_tasks$TMB$Curve[[1]],
                col = color_TMB, lwd = 2, type = "S", cex.lab = 1.3,
                lty = 1, add = TRUE
            )
            legend_text <- do.call(
                c,
                lapply(
                    names(ROC_all_run_tasks)[seq_len(length(names(ROC_all_run_tasks)))],
                    function(descriptor) {
                        paste0(descriptor, " (", round(subset(
                            AUC_mean_sd_all_run_tasks,
                            View == descriptor
                        )$AUC.mean, 2), ")")
                    }
                )
            )
            graphics::legend(
                x = "topright", inset = c(-0.5, 0),
                legend = legend_text,
                col = c(all_colors, color_TMB), lty = 1, lwd = 2, cex = 0.8, bty = "n"
            )
        } else {
            graphics::legend(
                x = "topright", inset = c(-0.5, 0),
                legend = legend_text,
                col = all_colors, lty = 1, lwd = 2, cex = 0.8, bty = "n"
            )
        }
        plot_list[[2]] <- grDevices::recordPlot()

        # *******************************************
        # Figure: Scatterplot (easier with TMB)
        # *******************************************
        if (easier_with_TMB != "none") {
            rp_df <- data.frame(
                response = patient_response,
                prediction_easier = apply(overall_df, 1, mean),
                TMB = TMB_values
            )
            # Categorize TMB
            if (length(unique(rp_df$TMB)) > 3) {
                rp_df$TMBcat <- categorize_TMB(rp_df$TMB)
                rp_df$TMB <- rp_df$TMBcat
                # if I want specify the thresholds
                # rp.df$TMB <- categorize.TMB(rp.df$TMB, thresholds = c(100,400))
            }
            # compute the integrated score as weighted average #
            if (easier_with_TMB == "weighted_average") {
                pred_lin <- linear_func(rp_df$prediction_easier)
                TMB_lin <- linear_func(rp_df$TMB)
                AUC_averaged_v <- vapply(seq(from = 0, to = 1, by = 0.1), function(p) {
                    pred_averaged <- apply(cbind((1 - p) * pred_lin, (p) * TMB_lin), 1, mean)
                    pred <- ROCR::prediction(pred_averaged, rp_df$response)
                    AUC_averaged <- ROCR::performance(pred, measure = "auc")
                    AUC_averaged_v <- AUC_averaged@y.values[[1]]
                }, FUN.VALUE = numeric(1))
            }
            # compute the integrated score for different penalties #
            if (easier_with_TMB == "penalized_score") {
                pred_combined <- rp_df$prediction_easier
                AUC_combined_v <- vapply(seq(from = 0, to = 1, by = 0.1), function(p) {
                    pred_combined[rp_df$TMBcat == 1] <- pred_combined[rp_df$TMBcat == 1] - p
                    pred_combined[rp_df$TMBcat == 3] <- pred_combined[rp_df$TMBcat == 3] + p

                    pred <- ROCR::prediction(pred_combined, rp_df$response)
                    AUC_combined <- ROCR::performance(pred, measure = "auc")
                    AUC_combined_v <- AUC_combined@y.values[[1]]
                }, FUN.VALUE = numeric(1))
            }
            pred <- ROCR::prediction(rp_df$prediction_easier, rp_df$response)
            AUC_easier <- ROCR::performance(pred, measure = "auc")
            AUC_easier_v <- AUC_easier@y.values[[1]]
            # scatterplot
            graphics::par(
                cex.axis = 1.3, mar = c(5, 4, 2, 8), col.lab = "black",
                pty = "s", xpd = TRUE
            )
            if (easier_with_TMB == "penalized_score") {
                plot(seq(from = 0, to = 1, by = 0.1), AUC_combined_v,
                    xlab = "Penalty or Relative weight",
                    ylab = "Area under the curve (AUC)",
                    type = "b", col = "#c15050", lty = 1,
                    pch = 19, lwd = 2, ylim = c(0, 1),
                    xlim = c(0, 1), cex.lab = 1.3
                )
                graphics::legend(
                    x = "topright", inset = c(-0.5, 0),
                    legend = c("Penalized score", "EaSIeR", "TMB"),
                    col = c(
                        "#c15050", as.vector(color_ensemble), color_TMB
                    ), lty = 1, lwd = 2, cex = 0.8, bty = "n"
                )
            }
            if (easier_with_TMB == "weighted_average") {
                plot(seq(from = 0, to = 1, by = 0.1), AUC_averaged_v,
                    xlab = "Penalty or Relative weight",
                    ylab = "Area under the curve (AUC)",
                    type = "b", col = "#693c72", lty = 1,
                    pch = 19, lwd = 2, ylim = c(0, 1),
                    xlim = c(0, 1), cex.lab = 1.3
                )
                graphics::legend(
                    x = "topright", inset = c(-0.5, 0),
                    legend = c("Weighted average", "EaSIeR", "TMB"),
                    col = c(
                        "#693c72", as.vector(color_ensemble), color_TMB
                    ), lty = 1, lwd = 2, cex = 0.8, bty = "n"
                )
            }
            # TMB
            graphics::segments(
                x0 = 0, y0 = AUC_mean_sd_TMB_run_tasks$AUC.mean, x1 = 1,
                y1 = AUC_mean_sd_TMB_run_tasks$AUC.mean, col = color_TMB
            )
            # easier (ensemble)
            graphics::segments(
                x0 = 0, y0 = AUC_easier_v, x1 = 1, y1 = AUC_easier_v,
                col = color_ensemble
            )
            plot_list[[3]] <- grDevices::recordPlot()
        }
        return(plot_list)
    } else {
        # what if patients' response is not provided
        # *******************************************
        # Figure: Dotplot (easier prediction)
        # *******************************************
        if (verbose) message("Scoring patients' as real response is not provided\n")
        plot_list <- list()
        rp_df <- data.frame(
            prediction_easier = apply(overall_df, 1, mean),
            patient = rownames(overall_df)
        )
        all_scores_df <- reshape2::melt(overall_df)
        names(all_scores_df) <- c("patient", "approach", "pred")
        all_scores_df$approach <- factor(all_scores_df$approach, levels = tasks)
        # sort patients
        easier <- stats::aggregate(pred ~ patient, data = all_scores_df, FUN = "median")
        easier <- easier[match(sort(easier$pred, decreasing = TRUE), easier$pred), ]
        all_scores_df$patient <- factor(all_scores_df$patient,
            levels = as.character(easier$patient)
        )
        ordering <- all_scores_df$patient
        score_plot <- ggplot2::ggplot(
            all_scores_df,
            ggplot2::aes(x = .data$pred, y = .data$patient)
        ) +
            ggplot2::geom_boxplot() +
            ggplot2::theme(
                panel.grid = ggplot2::element_blank(),
                panel.background = ggplot2::element_rect(fill = NA)
            ) +
            ggplot2::theme(
                axis.text.x = ggplot2::element_text(
                    size = 12, face = "bold", angle = 0,
                    vjust = 0.5, hjust = 0.5, color = "black"
                ),
                axis.title = ggplot2::element_text(size = 14),
                axis.text.y = ggplot2::element_text(size = 12),
                axis.ticks.y = ggplot2::element_blank(),
                panel.background = ggplot2::element_rect(fill = "white"),
                legend.position = "right", legend.direction = "vertical",
                panel.grid.major.y = ggplot2::element_line(linetype = 1, colour = "gray"),
                legend.text = ggplot2::element_text(size = 10),
                legend.title = ggplot2::element_text(size = 10, face = "bold", vjust = 0.5),
                strip.background = ggplot2::element_rect(fill = "white", colour = "white"),
                strip.text = ggplot2::element_text(size = 10)
            ) +
            ggplot2::labs(x = "easier prediction score", y = "patients") +
            ggplot2::guides(fill = "none", color = "none", shape = "none")

        plot_list[[1]] <- print(score_plot)
        # *******************************************
        # Figure: Dotplot (easier with TMB prediction)
        # *******************************************
        if (easier_with_TMB != "none") {
            # default weight_penalty value = 0.5
            if (missing(weight_penalty)) weight_penalty <- 0.5
            # Add TMB #
            rp_df$TMB <- TMB_values
            # Categorize TMB #
            if (length(unique(rp_df$TMB)) > 3) {
                rp_df$TMBcat <- categorize_TMB(rp_df$TMB)
                rp_df$TMB <- rp_df$TMBcat
                # if I want specify the thresholds
                # rp.df$TMB <- categorize_TMB(rp.df$TMB, thresholds = c(100,400))
            }
            # Compute the integrated score as weighted average #
            if (easier_with_TMB == "weighted_average") {
                # Compute the integrated score as weighted average #
                pred_lin <- linear_func(rp_df$prediction_easier)
                names(pred_lin) <- rp_df$patient
                TMB_lin <- linear_func(rp_df$TMB)
                names(TMB_lin) <- rp_df$patient
                pred_averaged_rf <- vapply(seq(from = 0, to = 1, by = 0.1), function(p) {
                    pred_averaged <- apply(cbind((1 - p) * pred_lin, (p) * TMB_lin), 1, mean)
                }, FUN.VALUE = numeric(length(pred_lin)))
                # By default, weight or penalty of 0.5.
                weight_penalty_pos <- match(weight_penalty, (seq_len(11) - 1) / 10)
                rp_df$weighted_average <- pred_averaged_rf[, weight_penalty_pos]
            }
            # Compute the integrated score for different penalties #
            if (easier_with_TMB == "penalized_score") {
                pred_combined <- rp_df$prediction_easier
                names(pred_combined) <- rp_df$patient
                pred_combined_rf <- vapply(seq(from = 0, to = 1, by = 0.1), function(p) {
                    pred_combined[rp_df$TMBcat == 1] <- pred_combined[rp_df$TMBcat == 1] - p
                    pred_combined[rp_df$TMBcat == 3] <- pred_combined[rp_df$TMBcat == 3] + p
                    return(pred_combined)
                }, FUN.VALUE = numeric(numeric(length(pred_combined))))
                weight_penalty_pos <- match(weight_penalty, (seq_len(11) - 1) / 10)
                rp_df$penalized_score <- NA
                rp_df$penalized_score <- pred_combined_rf[, weight_penalty_pos]
            }
            rp_df$patient <- paste0(rownames(rp_df), " (TMBcat=", rp_df$TMBcat, ")")
            rp_df$TMB <- NULL
            rp_df$TMBcat <- NULL
            all_scores_df <- reshape2::melt(rp_df)
            names(all_scores_df) <- c("patient", "approach", "pred")
            # Keep only weighted average
            all_scores_df <- all_scores_df[all_scores_df$approach == "weighted_average", ]
            # sort patients
            all_scores_df$patient <- factor(all_scores_df$patient,
                levels = all_scores_df$patient[match(
                    levels(ordering),
                    vapply(strsplit(all_scores_df$patient, split = " "), head, 1,
                        FUN.VALUE = character(1)
                    )
                )]
            )
            rf_score_plot <- ggplot2::ggplot(
                all_scores_df,
                ggplot2::aes(
                    x = pred,
                    y = .data$patient,
                    shape = .data$approach
                )
            ) +
                ggplot2::geom_point(ggplot2::aes(
                    fill = .data$approach,
                    colour = .data$approach
                ), size = 2) +
                ggplot2::scale_fill_manual(
                    name = "Approach",
                    labels = as.character(unique(all_scores_df$approach)),
                    values = c("#cc8100", "#00bc93")
                ) +
                ggplot2::scale_shape_manual(
                    name = "Approach",
                    labels = as.character(unique(all_scores_df$approach)),
                    values = c(21, 22)
                ) +
                ggplot2::scale_color_manual(
                    name = "Approach",
                    labels = as.character(unique(all_scores_df$approach)),
                    values = c("#cc8100", "#00bc93")
                ) +
                ggplot2::xlim(c(0, 1)) +
                ggplot2::theme(
                    panel.grid = ggplot2::element_blank(),
                    panel.background = ggplot2::element_rect(fill = NA)
                ) +
                ggplot2::theme(
                    axis.text.x = ggplot2::element_text(
                        size = 12, face = "bold", angle = 0,
                        vjust = 0.5, hjust = 0.5, color = "black"
                    ),
                    axis.title = ggplot2::element_text(size = 14),
                    axis.text.y = ggplot2::element_text(size = 12),
                    axis.ticks.y = ggplot2::element_blank(),
                    panel.background = ggplot2::element_rect(fill = "white"),
                    legend.position = "top", legend.direction = "horizontal",
                    panel.grid.major.y = ggplot2::element_line(linetype = 1, colour = "gray"),
                    legend.text = ggplot2::element_text(size = 10),
                    legend.title = ggplot2::element_text(size = 10, face = "bold", vjust = 0.5),
                    strip.background = ggplot2::element_rect(fill = "white", colour = "white"),
                    strip.text = ggplot2::element_text(size = 10),
                    legend.justification = c("right", "top"),
                    legend.key = ggplot2::element_rect(fill = "white")
                ) +
                ggplot2::labs(x = "prediction", y = "patients")

            plot_list[[2]] <- print(rf_score_plot)
        }
        return(plot_list)
    }
}
