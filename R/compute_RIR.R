#' Compute Immune resistance program
#'
#' `compute_IRP` computes immune resistance program using the code provided by
#' the authors. (Jerby-Arnon et al., 2018).
#' TODOTODOoscar: out of sync name vs doc
#'
#' @importFrom stats na.omit
#'
#' @param RNA_tpm numeric matrix with rows=genes and columns=samples
#' @param verbose A logical value indicating whether to display informative messages
#'
#' @return numeric matrix with rows=samples and columns=IRP
#'
#' @export
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
#' gene_tpm <- mariathasan_data$tpm
#' rm(cds)
#'
#' RIR <- compute_RIR(RNA_tpm = gene_tpm)
#' head(RIR)
compute_RIR <- function(RNA_tpm,
                        verbose = TRUE) {

  # Literature signature
  sig_read <- unique(unlist(res_sig))
  match_sig_read <- match(sig_read, rownames(RNA_tpm))

  if (anyNA(match_sig_read)) {
    warning("differently named or missing signature genes : \n", paste(sig_read[!sig_read %in% rownames(RNA_tpm)], collapse = "\n"), "\n")
    match_sig_read <- stats::na.omit(match_sig_read)
  }

  # Log2 transformation:
  log2_RNA_tpm <- log2(RNA_tpm + 1)

  # Prepare input data
  r <- list()
  r$tpm <- log2_RNA_tpm
  r$genes <- rownames(log2_RNA_tpm)

  # Apply function to calculate OE:
  res_scores <- get_OE_bulk(r, gene_sign = res_sig, verbose = TRUE)

  # Merge as recommend by authors
  res <- cbind.data.frame(
    excF.up = rowMeans(res_scores[, c("exc.up", "exc.seed.up")]),
    excF.down = rowMeans(res_scores[, c("exc.down", "exc.seed.down")]),
    res.up = rowMeans(res_scores[, c("trt.up", "exc.up", "exc.seed.up")]),
    res.down = rowMeans(res_scores[, c("trt.down", "exc.down", "exc.seed.down")]),
    res_scores
  )

  res <- cbind.data.frame(
    resF.up = res[, "res.up"] + res[, "fnc.up"],
    resF.down = res[, "res.down"] + res[, "fnc.down"],
    res
  )

  # Keep that signature considered to be relevant
  keep_sig <- c("resF.down")
  score <- as.matrix(res[, colnames(res) %in% keep_sig])
  rownames(score) <- colnames(log2_RNA_tpm)

  if (verbose) message("RIR score computed")
  return(data.frame(RIR = score, check.names = FALSE))
}
