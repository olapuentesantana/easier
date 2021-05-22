#' Convert tumor mutational burden (TMB) to categorical.
#'
#' `categorize_TMB` encodes tumor mutational variable (TMB) numerical variable to categorical.
#'
#' @importFrom stats quantile
#'
#' @export
#'
#' @param TMB A numeric vector with tumor mutational burden values.
#' @param thresholds A numeric vector to specify thresholds to be used. Default thresholds are low (<100), moderate (100-400) and high TMB (>400).
#'
#' @return A numeric vector assigning each sample a class from 1 to 3.
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
#'data(cds)
#'mariathasan_data <- preprocess_mariathasan(cds)
#'rm(cds)
#'
#'# retrieve Tumor Mutational Burden
#'TumorMutationalBurden <- clinical_data[, "FMOne mutation burden per MB"]
#'names(TumorMutationalBurden) <- rownames(clinical_data)
#'
#' TMB_cat <- categorize_TMB(TumorMutationalBurden)
#' TMB_cat
categorize_TMB <- function(x, thresholds=NULL) {
  vTert <- stats::quantile(x , c(0:3/3))
  if (is.null(thresholds)){
    x_out <- cut(x, vTert, include.lowest = TRUE, labels = c(1, 2, 3))
  }else if ((is.numeric(thresholds)) & (length(thresholds)==2)){
    # x_out = cut(x, c(min(x, na.rm = T), 100, 400, max(x, na.rm = T)), include.lowest = T, labels = c(1, 2, 3))
    x_out <- cut(x, c(min(x, na.rm = TRUE), 100, 400, max(x, na.rm = TRUE)), include.lowest = FALSE, labels = c(1, 2, 3))
  }
  x_out <- as.numeric(x_out)
  names(x_out) <- names(x)
  return(x_out)
}
