#' Compute Immuno-Predictive Score (IMPRES)
#'
#' \code{compute_IMPRES} computes IMPRES score by applying logical comparison of checkpoint gene pairs
#' (Auslander et al., 2018).
#'
#' @importFrom stats na.omit
#'
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=IMPRES score
#'
#' @export
compute.IMPRES <- function(RNA.tpm){

  # Literature genes
  IMPRES.basis <- data.frame(Gene_1 = c("PDCD1","CD27","CTLA4","CD40","CD86", "CD28", "CD80",
                                        "CD274","CD86","CD40","CD86","CD40","CD28","CD40","TNFRSF14"),
                             Gene_2 = c("TNFSF4","PDCD1","TNFSF4","CD28","TNFSF4", "CD86", "TNFSF9",
                                        "C10orf54","HAVCR2","PDCD1","CD200","CD80","CD276","CD274","CD86"))

  IMPRES.read <- unique(as.vector(as.matrix(IMPRES.basis))) # 15 genes

  # EQUIVALENT : "VISTA" = "C10orf54", "PDL-1" = "CD274", "TIM-3" = "HAVCR2",
  # "PD-1" = "PDCD1", "HVEM" = "TNFRSF14", "OX40L" = "TNFSF4", "CD137L" = "TNFSF9"

  # Some genes might have other name: case for "C10orf54", it's called "VSIR", be carefull
  if (any(rownames(RNA.tpm) == "VSIR")){
    cat("Gene name changed: C10orf54 instead of VSIR","\n")
    rownames(RNA.tpm)[which(rownames(RNA.tpm) == "VSIR")] = "C10orf54"
  }

  # Subset RNA.tpm
  match_F_1 <- match(as.character(IMPRES.basis[,1]), rownames(RNA.tpm))
  match_F_2 <- match(as.character(IMPRES.basis[,2]), rownames(RNA.tpm))

  if (anyNA(c(match_F_1, match_F_2))) {
    warning(c("differenty named or missing signature genes : \n", paste(IMPRES.read[!IMPRES.read %in% rownames(RNA.tpm)], collapse = "\n")))
  }

  # Initialize variables
  F_pair_expr_A <- matrix(0, nrow(IMPRES.basis), ncol(RNA.tpm))
  F_pair_expr_B <- matrix(0, nrow(IMPRES.basis), ncol(RNA.tpm))
  IMPRES.matrix <- matrix(0, nrow(IMPRES.basis), ncol(RNA.tpm)) ; colnames(IMPRES.matrix) <- colnames(RNA.tpm)
  score <- vector("numeric", length = ncol(RNA.tpm)) ; names(score) <- colnames(RNA.tpm)

  # Log2 transformation:
  log2.RNA.tpm <- as.data.frame(log2(RNA.tpm + 1))

  # Calculation:
  F_pair_expr_A <- log2.RNA.tpm[match_F_1, ]
  F_pair_expr_B <- log2.RNA.tpm[match_F_2, ]

  if(anyNA(F_pair_expr_A + F_pair_expr_B)) {
    remove_pairs <- as.vector(which(is.na(rowSums(F_pair_expr_A + F_pair_expr_B) == TRUE)))
  }

  IMPRES.matrix <- F_pair_expr_A > F_pair_expr_B
  if(anyNA(IMPRES.matrix)){
    score <- colSums(IMPRES.matrix, na.rm = TRUE)
    score <- (score * nrow(IMPRES.matrix)) / (nrow(IMPRES.matrix) - length(remove_pairs))
  }else{
    score <- colSums(IMPRES.matrix)
  }

  message("IMPRES score computed")
  return(data.frame(IMPRES = score, check.names = FALSE))
}
