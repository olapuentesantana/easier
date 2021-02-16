#' Compute MSI score
#'
#' \code{compute_MSI} computes MSI score by applying logical comparison of MSI-related gene pairs
#' (Fu et al., 2019).
#'
#' @importFrom stats na.omit
#'
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=MSI score
#'
#' @export
#'
#' @examples
#' # TODOTODO
compute.MSI <- function(RNA.tpm){

  # Literature genes: * (CCRN4L in tcga, NOCT approved symbol)
  MSI.basis <- data.frame(Gene_1 = c("HNRNPL","MTA2","CALR","RASL11A","LYG1", "STRN3", "HPSE",
                                     "PRPF39","NOCT","AMFR"),
                          Gene_2 = c("CDC16","VGF","SEC22B","CAB39L","DHRS12", "TMEM192", "BCAS3",
                                     "ATF6","GRM8","DUSP18"))
  MSI.read <- unique(as.vector(as.matrix(MSI.basis))) # 20 genes

  # Some genes might have other name: case for "CCRN4L", it's called "NOCT", be carefull
  if (any(rownames(RNA.tpm) %in% "CCRN4L")){
    cat("Gene name changed: NOCT is approved symbol, not CCRN4L","\n")
    rownames(RNA.tpm)[rownames(RNA.tpm) %in% "CCRN4L"] <- "NOCT"
  }

  # Subset RNA.tpm
  match_F_1 <- match(as.character(MSI.basis[,1]), rownames(RNA.tpm))
  match_F_2 <- match(as.character(MSI.basis[,2]), rownames(RNA.tpm))

  if (anyNA(c(match_F_1,match_F_2))){
    warning(c("differenty named or missing signature genes : \n", paste(MSI.read[!MSI.read %in% rownames(RNA.tpm)], collapse = "\n")))
  }

  # Initialize variables
  F_pair_expr_A <- matrix(0, nrow(MSI.basis), ncol(RNA.tpm))
  F_pair_expr_B <- matrix(0, nrow(MSI.basis), ncol(RNA.tpm))
  MSI.matrix <- matrix(0, nrow(MSI.basis), ncol(RNA.tpm)) ; colnames(MSI.matrix) <- colnames(RNA.tpm)
  remove_pairs <- vector("list", length = ncol(RNA.tpm)) ; names(remove_pairs) <- colnames(RNA.tpm)
  score <- vector("numeric", length = ncol(RNA.tpm)) ; names(score) <- colnames(RNA.tpm)

  # Log2 transformation:
  log2.RNA.tpm <- as.data.frame(log2(RNA.tpm + 1))

  # Calculation:
  F_pair_expr_A <- log2.RNA.tpm[match_F_1, ]
  F_pair_expr_B <- log2.RNA.tpm[match_F_2, ]

  if(anyNA(F_pair_expr_A + F_pair_expr_B)) {
    remove_pairs <- as.vector(which(is.na(rowSums(F_pair_expr_A + F_pair_expr_B) == TRUE)))
  }

  MSI.matrix <- F_pair_expr_A > F_pair_expr_B
  if(anyNA(MSI.matrix)){
    score <- colSums(MSI.matrix, na.rm = TRUE)
    score <- (score * nrow(MSI.matrix)) / (nrow(MSI.matrix) - length(remove_pairs))
  }else{
    score <- colSums(MSI.matrix)
  }

  message("MSI score computed")
  return(data.frame(MSI = score, check.names = FALSE))
}
