#' Compute repressed immune resistance signature (RIR) score
#'
#' This function calculates RIR score defined by combining a set of gene signatures associated with downregulation of T cell exclusion, post-treatment and functional resistance.
#' More info can be found in original work from Jerby-Arnon et al., Cell, 2018.
#'
#' @references Jerby-Arnon, L., Shah, P., Cuoco, M.S., Rodman, C., Su, M.-J., Melms, J.C., Leeson, R., Kanodia, A., Mei, S., Lin, J.-R., et al. (2018).
#' A Cancer Cell Program Promotes T Cell Exclusion and Resistance to Checkpoint Blockade. Cell 175, 984–997.e24. https://doi.org/10.1016/j.cell.2018.09.006.
#'
#' @importFrom stats na.omit
#'
#' @param RNA_tpm data.frame containing TPM values with HGNC symbols in rows and samples in columns.
#' @param verbose logical variable indicating whether to display informative messages.
#'
#' @return A numeric matrix with samples in rows and RIR score in a column.
#'
#' @export
#'
#' @examples
#' # use example dataset from IMvigor210CoreBiologies package (Mariathasan et al., Nature, 2018)
#' data("dataset_mariathasan")
#' gene_tpm <- dataset_mariathasan@tpm
#'
#' # compute RIR signature score
#' RIR <- compute_RIR(RNA_tpm = gene_tpm)
compute_RIR <- function(RNA_tpm,
                        verbose = TRUE) {

  # Literature signature
  sig_read <- unique(unlist(res_sig))
  match_sig_read <- match(sig_read, rownames(RNA_tpm))

  if (anyNA(match_sig_read)) {
    warning("differently named or missing signature genes : \n", paste(sig_read[!sig_read %in% rownames(RNA_tpm)], collapse = "\n"), "\n")
    match_sig_read <- stats::na.omit(match_sig_read)
    # re-annotate genes not found
    out <- reannotate_genes(sig_read[!sig_read %in% rownames(RNA_tpm)])
    sig_read[match(out$old_names[!is.na(out$new_names)], sig_read)] <- out$new_names[!is.na(out$new_names)]
    warning("after gene re-annotation, differently named or missing signature genes : \n", paste(sig_read[!sig_read %in% rownames(RNA_tpm)], collapse = "\n"), "\n")
  }

  new_res_sig <- sapply(names(res_sig), function(X) {
    if (any(is.na(match(out$old_names, res_sig[[X]])) == FALSE)) {
      res_sig[[X]][stats::na.omit(match(out$old_names, res_sig[[X]]))] <- out$new_names[!is.na(match(out$old_names, res_sig[[X]]))]
    }
    return(res_sig[[X]])
  })

  # Log2 transformation:
  log2_RNA_tpm <- log2(RNA_tpm + 1)

  # Prepare input data
  r <- list()
  r$tpm <- log2_RNA_tpm
  r$genes <- rownames(log2_RNA_tpm)

  # Apply function to calculate OE:
  res_scores <- get_OE_bulk(r, gene_sign = new_res_sig, verbose = TRUE)

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
