#' Compute Immune resistance program
#'
#' \code{compute_IRP} computes immune resistance program using the code provided by the authors.
#' (Jerby-Arnon et al., 2018).
#'
#' @importFrom stats na.omit
#'
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=IRP
#'
#' @export
compute.RIR <- function(RNA.tpm){

  # Literature genes
  RIR.read <- unique(unlist(res.sig))
  match_RIR.read <- match(RIR.read, rownames(RNA.tpm))

  if (anyNA(match_RIR.read)){
    warning(c("differently named or missing signature genes : \n", paste(RIR.read[!RIR.read %in% rownames(RNA.tpm)], collapse = "\n")))
    match_RIR.read <- stats::na.omit(match_RIR.read)
  }

  # Log2 transformation:
  log2.RNA.tpm <- log2(RNA.tpm + 1)

  # Prepare input data
  r <- list()
  r$tpm <- log2.RNA.tpm
  r$genes <- rownames(log2.RNA.tpm)

  # Apply function to calculate OE:
  res.scores <- get_OE_bulk(r, gene.sig = res.sig)

  # Merge as recommend by authors
  res <- cbind.data.frame(excF.up = rowMeans(res.scores[, c("exc.up", "exc.seed.up")]),
                          excF.down = rowMeans(res.scores[, c("exc.down", "exc.seed.down")]),
                          res.up = rowMeans(res.scores[, c("trt.up", "exc.up", "exc.seed.up")]),
                          res.down = rowMeans(res.scores[, c("trt.down", "exc.down", "exc.seed.down")]),
                          res.scores)

  res <- cbind.data.frame(resF.up = res[ , "res.up"] + res[, "fnc.up"],
                          resF.down = res[, "res.down"] + res[, "fnc.down"],
                          res)

  # Keep that signature considered to be relevant
  keep.sig <- c("resF.down")
  score <- as.matrix(res[, colnames(res) %in% keep.sig]); rownames(score) <- colnames(RNA.tpm)

  message("RIR score computed")
  return(data.frame(RIR = score, check.names = FALSE))
}
