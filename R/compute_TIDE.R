#' Compute TIDE score
#'
#' `compute_TIDE` computes TIDE score using tidepy python module. Input data is
#' processed as recommended by the authors: log2 transformation, quantile
#' normalization and mean centralization (Jiang et al., 2018).
#'
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom utils read.table write.table
#'
#' @param RNA_tpm numeric matrix with rows=genes and columns=samples
#' @param cancertype string character
#' @param output_file_path TODOTODO to define
#' @param verbose A logical value indicating whether to display informative messages
#'
#' @return numeric matrix with rows=samples and columns=TIDE score
#'
#' @export
#'
#' @examples
#' # Example: Riaz
#' data("Riaz_data")
#'
#' ##### TIDE <- compute_TIDE(
#' #####   RNA_tpm= Riaz_data$tpm_RNAseq,
#' #####   cancertype = "SKCM",
#' #####   output_file_path = tempdir())
#' ##### head(TIDE)
#' # TODO: will need to fix the anaconda part...
compute_TIDE <- function(RNA_tpm,
                         cancertype,
                         output_file_path,
                         verbose = TRUE) {

  # Pre-process #

  # Log2 transformation:
  log2_RNA_tpm <- as.data.frame(log2(RNA_tpm + 1))

  # Quantile normalization
  norm_log2_RNA_tpm <- preprocessCore::normalize.quantiles(as.matrix(log2_RNA_tpm))
  dimnames(norm_log2_RNA_tpm) <- dimnames(log2_RNA_tpm)

  # Mean-centralization: substract gene average across all conditions
  average_genes <- rowMeans(norm_log2_RNA_tpm)
  norm_log2_RNA_tpm <- sweep(norm_log2_RNA_tpm, 1, average_genes, FUN = "-")

  utils::write.table(norm_log2_RNA_tpm, file = paste0(output_file_path, "/log2mas1TPM_norm_", cancertype, ".txt"), sep = "\t")

  warning("tidepy need to be installed")

  # Calculation:
  if (cancertype == "SKCM") {
    system(paste0("/anaconda3/bin/tidepy ", output_file_path, "/log2mas1TPM_norm_", cancertype, ".txt", " -o ", output_file_path, "/output_TIDE_norm_", cancertype, ".txt -c Melanoma"))
  } else {
    system(paste0("/anaconda3/bin/tidepy ", output_file_path, "/log2mas1TPM_norm_", cancertype, ".txt", " -o ", output_file_path, "/output_TIDE_norm_", cancertype, ".txt -c Other"))
  }

  TIDE_output <- utils::read.table(file = paste0(output_file_path, "/output_TIDE_norm_", cancertype, ".txt"), sep = "\t", header = TRUE, row.names = 1)
  score <- TIDE_output[, "TIDE", drop = FALSE]

  system(paste0("rm ", output_file_path, "/output_TIDE_norm_", cancertype, ".txt"))
  system(paste0("rm ", output_file_path, "/log2mas1TPM_norm_", cancertype, ".txt"))

  if (verbose) message("TIDE score computed")
  return(data.frame(TIDE = score, check.names = FALSE))
}
