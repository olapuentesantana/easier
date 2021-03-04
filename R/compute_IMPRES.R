#' Compute Immuno-Predictive Score (IMPRES)
#'
#' `compute_IMPRES` computes IMPRES score by applying logical comparison of
#' checkpoint gene pairs (Auslander et al., 2018).
#'
#' @importFrom stats na.omit
#'
#' @param RNA_tpm numeric matrix with rows=genes and columns=samples
#' @param verbose A logical value indicating whether to display informative messages
#'
#' @return numeric matrix with rows=samples and columns=IMPRES score
#'
#' @export
#'
#' @examples
#' # Example: Riaz
#' data("Riaz_data")
#'
#' IMPRES <- compute_IMPRES(RNA_tpm = Riaz_data$tpm_RNAseq)
#' head(IMPRES)
compute_IMPRES <- function(RNA_tpm,
                           verbose = TRUE) {

  # Literature signature
  sig_basis <- data.frame(
    Gene_1 = c(
      "PDCD1", "CD27", "CTLA4", "CD40", "CD86", "CD28", "CD80",
      "CD274", "CD86", "CD40", "CD86", "CD40", "CD28", "CD40", "TNFRSF14"
    ),
    Gene_2 = c(
      "TNFSF4", "PDCD1", "TNFSF4", "CD28", "TNFSF4", "CD86", "TNFSF9",
      "C10orf54", "HAVCR2", "PDCD1", "CD200", "CD80", "CD276", "CD274", "CD86"
    )
  )

  sig_read <- unique(as.vector(as.matrix(sig_basis))) # 15 genes

  # EQUIVALENT : "VISTA" = "C10orf54", "PDL-1" = "CD274", "TIM-3" = "HAVCR2",
  # "PD-1" = "PDCD1", "HVEM" = "TNFRSF14", "OX40L" = "TNFSF4", "CD137L" = "TNFSF9"

  # Some genes might have other name: case for "C10orf54", it's called "VSIR", be carefull
  if (any(rownames(RNA_tpm) == "VSIR")) {
    warning("Gene name changed: C10orf54 instead of VSIR", "\n")
    rownames(RNA_tpm)[which(rownames(RNA_tpm) == "VSIR")] <- "C10orf54"
  }

  # Subset RNA_tpm
  match_F_1 <- match(as.character(sig_basis[, 1]), rownames(RNA_tpm))
  match_F_2 <- match(as.character(sig_basis[, 2]), rownames(RNA_tpm))

  if (anyNA(c(match_F_1, match_F_2))) {
    warning("differenty named or missing signature genes : \n",
      paste(sig_read[!sig_read %in% rownames(RNA_tpm)], collapse = "\n"),
      "\n"
    )
  }

  # Initialize variables
  F_pair_expr_A <- matrix(0, nrow(sig_basis), ncol(RNA_tpm))
  F_pair_expr_B <- matrix(0, nrow(sig_basis), ncol(RNA_tpm))
  F_matrix <- matrix(0, nrow(sig_basis), ncol(RNA_tpm))
  colnames(F_matrix) <- colnames(RNA_tpm)
  score <- vector("numeric", length = ncol(RNA_tpm))
  names(score) <- colnames(RNA_tpm)

  # Log2 transformation:
  log2_RNA_tpm <- as.data.frame(log2(RNA_tpm + 1))

  # Calculation:
  F_pair_expr_A <- log2_RNA_tpm[match_F_1, ]
  F_pair_expr_B <- log2_RNA_tpm[match_F_2, ]

  if (anyNA(F_pair_expr_A + F_pair_expr_B)) {
    remove_pairs <- as.vector(which(is.na(rowSums(F_pair_expr_A + F_pair_expr_B) == TRUE)))
  }

  F_matrix <- F_pair_expr_A > F_pair_expr_B
  if (anyNA(F_matrix)) {
    score <- colSums(F_matrix, na.rm = TRUE)
    score <- (score * nrow(F_matrix)) / (nrow(F_matrix) - length(remove_pairs))
  } else {
    score <- colSums(F_matrix)
  }

  if (verbose) message("IMPRES score computed")
  return(data.frame(IMPRES = score, check.names = FALSE))
}
