#' Compute Roh immune score
#'
#' `compute_Roh_IS` computes Roh immune score as the geometric-mean of immune
#' score genes (Roh et al., 2017).
#'
#' @importFrom stats na.omit
#'
#' @param RNA_tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=Roh immune score
#'
#' @export
#'
#' @examples
#' # Example: Riaz
#' data("Riaz_data")
#'
#' Roh_IS <- compute_Roh_IS(RNA_tpm = Riaz_data$tpm_RNAseq)
#' head(Roh_IS)
compute_Roh_IS <- function(RNA_tpm) {

  # Literature genes
  Roh_IS.read <- c(
    "GZMA", "GZMB", "PRF1", "GNLY", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F",
    "HLA-G", "HLA-H", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DPA1",
    "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DRA", "HLA-DRB1",
    "IFNG", "IFNGR1", "IFNGR2", "IRF1", "STAT1", "PSMB9", "CCR5", "CCL3", "CCL4",
    "CCL5", "CXCL9", "CXCL10", "CXCL11", "ICAM1", "ICAM2", "ICAM3", "ICAM4",
    "ICAM5", "VCAM1"
  )

  match_Roh_IS.genes <- match(Roh_IS.read, rownames(RNA_tpm))

  if (anyNA(match_Roh_IS.genes)) {
    warning(c("differenty named or missing signature genes : \n", paste(Roh_IS.read[!Roh_IS.read %in% rownames(RNA_tpm)], collapse = "\n")))
    match_Roh_IS.genes <- stats::na.omit(match_Roh_IS.genes)
  }

  # Subset RNA_tpm
  sub_gene.tpm <- RNA_tpm[match_Roh_IS.genes, ]

  # Pseudocount of 0.01 for all genes
  sub_gene.tpm <- sub_gene.tpm + 0.01

  # Pseudocount of 1 for genes with 0 expr
  if (any(sub_gene.tpm == 0)) sub_gene.tpm[sub_gene.tpm == 0] <- sub_gene.tpm[sub_gene.tpm == 0] + 1

  # Calculation: geometric mean (so-called log-average) [TPM, 0.01 offset]
  score <- apply(sub_gene.tpm, 2, function(X) exp(mean(log(X))))

  message("Roh_IS computed score")
  return(data.frame(Roh_IS = score, check.names = FALSE))
}
