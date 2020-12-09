#' Compute TIDE score
#'
#' \code{compute_TIDE} computes TIDE score using tidepy python module. Input data is processed as recommended
#' by the authors: log2 transformation, quantile normalization and mean centralization (Jiang et al., 2018).
#'
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom utils read.table write.table
#'
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#' @param cancertype string character
#'
#' @return numeric matrix with rows=samples and columns=TIDE score
#'
#' @export
compute.TIDE <- function(RNA.tpm, cancertype, output_file_path){

  # Pre-process #

  # Log2 transformation:
  log2.RNA.tpm <- as.data.frame(log2(RNA.tpm + 1))

  # Quantile normalization
  log2.RNA.tpm.norm <- preprocessCore::normalize.quantiles(as.matrix(log2.RNA.tpm))
  dimnames(log2.RNA.tpm.norm) <- dimnames(log2.RNA.tpm)

  # Mean-centralization: substract gene average across all conditions
  average.gene <- rowMeans(log2.RNA.tpm.norm)
  log2.RNA.tpm.norm <- sweep(log2.RNA.tpm.norm, 1, average.gene, FUN = "-")

  utils::write.table(log2.RNA.tpm.norm, file = paste0(output_file_path, "/log2mas1TPM_norm_", cancertype,".txt"), sep = "\t")

  warning("tidepy need to be installed")

  # Calculation:
  if (cancertype == "SKCM") {
    try(system(paste0("/anaconda3/bin/tidepy ", output_file_path,"/log2mas1TPM_norm_", cancertype,".txt", " -o ",output_file_path,"/output_TIDE_norm_", cancertype,".txt -c Melanoma")))
  }else{
    try(system(paste0("/anaconda3/bin/tidepy ", output_file_path,"/log2mas1TPM_norm_", cancertype,".txt", " -o ",output_file_path,"/output_TIDE_norm_",cancertype,".txt -c Other")))
  }

  TIDE.table.norm <- utils::read.table(file = paste0(output_file_path, "/output_TIDE_norm_", cancertype,".txt"), sep = "\t", header = T, row.names = 1)
  score <- TIDE.table.norm[, "TIDE", drop = FALSE]

  system(paste0("rm ", output_file_path, "/output_TIDE_norm_", cancertype,".txt"))
  system(paste0("rm ", output_file_path, "/log2mas1TPM_norm_", cancertype,".txt"))

  message("TIDE score computed")
  return(data.frame(TIDE = score, check.names = F))
}
