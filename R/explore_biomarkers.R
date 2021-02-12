#' Examination of possible relevant biomarkers in cancer dataset
#'
#' \code{explore_biomarkers} provides an overview of the mechanistic signatures features distribution in
#' the dataset (separating responders and non-responders if known) together with information on the
#' feature contribution to the model building.
#'
#' Plot original features distribution together with the
#' feature importance in the model
#'
#' @importFrom grDevices pdf dev.off
#' @importFrom stats aggregate
#' @import ggplot2
#' @importFrom grid unit.pmax grid.newpage grid.draw
#'
#' @export
#'
#' @param pathways numeric matrix with data
#' @param immunecells numeric matrix with data
#' @param lrpairs numeric matrix with data
#' @param cytokinepairs numeric matrix with data
#' @param tfs numeric matrix with data
#' @param ccpairsgroupedscores numeric matrix with data
#' @param ccpairsgroupedpval numeric matrix with data
#' @param cancertype string character
#' @param real_patient_response vector with two factors (NR,R)
#' @param output_file_path string with a file name
#'
#' @return boxplot with features distribution
#'
#' @examples
#' # TODOTODO
explore_biomarkers <- function(pathways,
                               immunecells,
                               lrpairs,
                               cytokinepairs,
                               tfs,
                               ccpairsgroupedscores,
                               ccpairsgroupedpval,
                               cancertype,
                               real_patient_response,
                               output_file_path){

  # Check that folder exists, create folder otherwise
  if(dir.exists(output_file_path) == FALSE) {
    dir.create(file.path(output_file_path), showWarnings = FALSE)
    warning(paste0(sapply(strsplit(output_file_path, "/", fixed = TRUE), tail , 1),
            " folder does not exist, creating ", sapply(strsplit(output_file_path, "/", fixed = TRUE), tail , 1), " folder"))
  }

  try(if(missing(cancertype)) stop("cancer type needs to be specified"))

  # Initialize variables
  views <- c(Pathways = 'gaussian',
             ImmuneCells = 'gaussian',
             LRpairs = 'gaussian',
             CYTOKINEpairs = 'gaussian',
             TFs = 'gaussian',
             CCpairsGroupedScores = 'gaussian',
             CCpairsGroupedPval = 'gaussian')

  view_combinations <- NULL

  algorithm <-  c("Multi_Task_EN") #"BEMKL"

  # Check which views are missing
  miss_views <- c(ifelse(missing(pathways), NA, 1),
                  ifelse(missing(immunecells), NA, 2),
                  ifelse(missing(lrpairs), NA, 3),
                  ifelse(missing(cytokinepairs), NA, 4),
                  ifelse(missing(tfs), NA, 5),
                  ifelse(missing(ccpairsgroupedscores), NA, 6),
                  ifelse(missing(ccpairsgroupedpval), NA, 7))

  # Possible combinations
  possible_combo <- combn(miss_views[1:2], 2)

  # Remove combinations with are not feasible due to missing views
  if(any(is.na(miss_views[1:2]))){
    possible_combo <- possible_combo[!is.na(miss_views[1:2]),]
  }

  # Views single
  view_simples <- lapply(miss_views[!is.na(miss_views)], function(X){
    tmp <- views[X] ; return(tmp)
  })

  # Views combination
  if(is.vector(possible_combo) & length(possible_combo) > 1) {
    view_combinations <- list(do.call(c, lapply(possible_combo, function(X){
      tmp <- views[X] ; return(tmp)
    })))
  }

  # Views combination
  if(is.matrix(possible_combo)) {
    view_combinations <- lapply(1:ncol(possible_combo), function(X){
      tmp <- views[possible_combo[,X]] ; return(tmp)
    })
  }
  view_combinations <- c(view_simples, view_combinations)

  # Remove unavailable combo
  combo_names <- sapply(1:length(view_combinations), function(X) {
    paste(names(view_combinations[[X]]), collapse = "_")
  })

  # Immune cells features curation:
  if (missing(immunecells) == FALSE){
    colnames(immunecells) <-  gsub(".","_", colnames(immunecells), fixed = TRUE)
  }

  # We explore biomarkers for each input separately:
  sapply(1:length(view_combinations), function(X){

      view_info <- view_combinations[[X]]
      view_name <- paste(names(view_info), collapse="_")
      view_data <- do.call(cbind, lapply(tolower(names(view_info)), function(X) as.data.frame(get(X))))

      # To improve visualization, data needs to be normalized (z-score)
      mas_mea_view_data <- apply(view_data, 2, FUN = "mean", na.rm = TRUE)
      mas_std_view_data <- apply(view_data, 2, FUN = "sd", na.rm = TRUE)
      view_data_z <- standarization(view_data, mas_mea_view_data, mas_std_view_data)

      learned_model <- trained_models[[cancertype]][[view_name]]

        # 100 iterations
        summary_iter <- do.call(rbind, lapply(1:length(learned_model), function(iteration){

          tmp_iter_model <- learned_model[[iteration]]$model

          # cv = 1se.mse
          tmp_iter_model_cv <- tmp_iter_model$cv.glmnet.features[["1se.mse"]]

          # # All tasks
          summary_task <- do.call(rbind, lapply(colnames(tmp_iter_model_cv), function(task){

            info <- data.frame(Iteration = iteration,
                               Hyp_model = paste0(tmp_iter_model$cv.glmnet.hyperparameters[["1se.mse"]]$alpha,",",
                                                  round(as.numeric(tmp_iter_model$cv.glmnet.hyperparameters[["1se.mse"]]$lambda),3)),
                               Task = task,
                               Feature = rownames(tmp_iter_model_cv)[2:nrow(tmp_iter_model_cv)],
                               Estimate = tmp_iter_model_cv[2:nrow(tmp_iter_model_cv),task])

            return(info)
          }))

          return(summary_task)
        }))

        len_summary <- nrow(summary_iter)
        summary_iter$DataType <- rep(view_name, len = len_summary)
        summary_iter$CancerType <- rep(cancertype, len =  len_summary)

        # Save .RData
        biomarkers_data <- summary_iter
        save(biomarkers_data, file = paste0(output_file_path,"/biomarkers_data_",view_name,"_",cancertype,".RData"))

        # Boxplots
        #1. Calculate median across iteractions
        median.features <- stats::aggregate(Estimate ~ Feature +  DataType + CancerType + Task,
                                            FUN = "median", na.rm = TRUE, data = summary_iter)
        #2. Calculate median across tasks
        median.features <- stats::aggregate(Estimate ~ Feature + DataType + CancerType,
                                     FUN = "median", na.rm = TRUE, data = median.features)
        #3. Sort median values (not necessary)
        median.features <- median.features[order(abs(median.features$Estimate), decreasing = TRUE),]
        median.features$Feature <- factor(median.features$Feature, levels = rev(unique(median.features$Feature)))

        # Original data (view_data)
        df <- view_data_z
        labels <- real_patient_response ; names(labels) <- rownames(df)

        view_distribution <- do.call(rbind, lapply(colnames(df), function(feature){

          df.feature <- df[,feature]

          view_distribution <- data.frame(Feature = feature,
                                          Patient = names(labels),
                                          Label = labels,
                                          Value = df.feature)
          return(view_distribution)
        }))

        view_distribution$Feature <- factor(view_distribution$Feature, levels = rev(unique(median.features$Feature)))
        view_distribution$Label <- factor(view_distribution$Label, levels = c("NR","R"))

        # Due to interpretability:
        # If a mechanistic signature contains a high number of features, we just keep those that the model is actually using.
        if (length(unique(median.features$Feature)) > 150){
          median.features <- subset(median.features, Estimate != 0)
          view_distribution <- subset(view_distribution, Feature %in% unique(median.features$Feature))
        }

        # barplot
        median.features$Correlation <- sign(median.features$Estimate)
        median.features$Correlation <- gsub("-1", "-", median.features$Correlation, fixed = TRUE)
        median.features$Correlation <- gsub("1", "+", median.features$Correlation, fixed = TRUE)
        median.features$Correlation <- factor(median.features$Correlation, levels = unique(median.features$Correlation))

        barplot <- ggplot2::ggplot(median.features, aes(x =  abs(Estimate), y = Feature, fill = Correlation)) +
          ggplot2::geom_bar(stat="identity", color="white") +
          ggplot2::scale_fill_manual(name = "Correlation",
                                     labels = levels(median.features$Correlation),
                                     values = c("-" = "#BB4444", "+" = "#4477AA", "0" = "gray")) +
          #ggplot2::theme_light() +
          #ggplot2::coord_flip() +
          ggplot2::coord_fixed(ratio = 0.04) +
          ggplot2::theme(panel.grid = element_blank()) +
          ggplot2::scale_y_discrete(labels = c("T_cells_CD8" = "CD8 T cells",
                                               "Macrophages_M2"= "M2",
                                               "B_cells" = "B cells",
                                               "T_cells_regulatory_Tregs" = "Tregs",
                                               "Macrophages_M1"= "M1",
                                               "T_cells_CD4" = "CD4 T cells",
                                               "NK_cells" = "NK cells",
                                               "Dendritic_cells" = "DC cells")) +
          ggplot2::theme(axis.text.x = element_text(size=12, color = "black"), axis.title.y = element_blank(), axis.text.y = element_blank(),
                         axis.ticks.y = element_blank(), axis.ticks.x = element_line(size = 0.5, color = "black"),
                         legend.position = "bottom", legend.direction = "horizontal",
                         legend.box.background = element_rect(color="black", size=0.3),
                         legend.box.margin = margin(0.5, 0.5, 0.5, 0.5),
                         legend.text=element_text(size=12),
                         legend.title = element_text(size=12, face="bold", vjust = 0.5),
                         panel.border = element_blank(), panel.background = element_blank(),
                         plot.margin=unit(c(1,1,1,1), "cm"), axis.line.x = element_line(colour ="black")) + # set negative value for the left margin)
          ggplot2::labs(x = "Biomarker weight")

        boxplot <- ggplot2::ggplot(view_distribution, aes(x = Feature, y = Value, fill = Label, color = Label)) +
          ggplot2::geom_boxplot(alpha = 0.8) +
          ggplot2::geom_point(position = position_jitterdodge()) +
          ggplot2::scale_fill_manual(name = "Label",
                                     labels = levels(view_distribution$Label),
                                     values = c("darkgrey", "black")) +
          ggplot2::scale_color_manual(name = "Label",
                                      labels = levels(view_distribution$Label),
                                      values = c("darkgrey", "black")) +
          # ggplot2::ylim(c(-5, 5)) +
          ggplot2::theme_minimal() +
          ggplot2::coord_flip() +
          ggplot2::theme(panel.grid = element_blank()) +
          ggplot2::scale_x_discrete(labels = c("T_cells_CD8" = "CD8 T cells",
                                               "Macrophages_M2"= "M2",
                                               "B_cells" = "B cells",
                                               "T_cells_regulatory_Tregs" = "Tregs",
                                               "Macrophages_M1"= "M1",
                                               "T_cells_CD4" = "CD4 T cells",
                                               "NK_cells" = "NK cells",
                                               "Dendritic_cells" = "DC cells")) +
          ggplot2::theme(axis.text.x = element_text(size=12, color = "black"), axis.text.y = element_text(size=12),
                         axis.title.y = element_blank(), axis.ticks.y = element_blank(),
                         legend.position = "bottom", legend.direction = "horizontal", axis.ticks.x = element_line(size = 0.5, color = "black"),
                         legend.box.background = element_rect(color="black", size=0.3),
                         legend.box.margin = margin(0.5, 0.5, 0.5, 0.5),
                         legend.text=element_text(size=12),
                         legend.title = element_text(size=12, face="bold", vjust = 0.5),
                         plot.margin=unit(c(1,1,1,1), "cm"), axis.line.x = element_line(colour ="black")) + # set negative value for the right margin
          ggplot2::labs(y = "Z-score")

        if (length(unique(median.features$Feature)) > 30){
          boxplot <- boxplot + ggplot2::theme(axis.text.y = element_text(size = 10))
        } else{
          boxplot <- boxplot + ggplot2::theme(axis.text.y = element_text(size = 12))
        }

        # Combine plots
        g1 <- ggplot2::ggplotGrob(boxplot)
        g2 <- ggplot2::ggplotGrob(barplot)
        g <- cbind(g1, g2, size = "first")
        g$heights <- grid::unit.pmax(g1$heights, g2$heights)
        grid::grid.newpage()
        pdf(paste0(output_file_path,"/biomarkers_plot_",view_name,"_",cancertype,".pdf"), width = 8, height = 12)
        grid::grid.draw(g)
        dev.off()

      })


  }
