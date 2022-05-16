#' Compute repressed immune resistance signature (RIR)
#' score
#'
#' Calculates RIR score by combining a set of gene
#' signatures associated with upregulation and
#' downregulation of T cell exclusion, post-treatment
#' and functional resistance.
#' We used the original approach defined in Jerby-Arnon
#' et al., Cell, 2018.
#'
#' The gene signatures were provided by original work:
#' https://github.com/livnatje/ImmuneResistance
#'
#' @references Jerby-Arnon, L., Shah, P., Cuoco, M.S.,
#' Rodman, C., Su, M.-J., Melms, J.C., Leeson, R., Kanodia,
#' A., Mei, S., Lin, J.-R., et al. (2018). A Cancer Cell
#' Program Promotes T Cell Exclusion and Resistance to
#' Checkpoint Blockade. Cell 175, 984–997.e24.
#' https://doi.org/10.1016/j.cell.2018.09.006.
#'
#' @importFrom stats na.omit
#'
#' @param RNA_tpm data.frame containing TPM values with HGNC symbols
#' in rows and samples in columns.
#' @param RIR_program list with gene signatures included in the immune
#' resistance program from Jerby-Arnon et al., 2018.
#'
#' @return A numeric matrix with samples in rows and three RIR scores as
#' columns: "resF_up" (upregulated score), "resF_down" (downregulated score)
#' and "resF" (upregulated score - downregulated score).
#'
compute_RIR <- function(RNA_tpm,
                        RIR_program) {

    # Log2 transformation:
    log2_RNA_tpm <- log2(RNA_tpm + 1)

    # Prepare input data
    r <- list()
    r$tpm <- log2_RNA_tpm
    r$genes <- rownames(log2_RNA_tpm)

    # Apply function to calculate OE:
    res_scores <- get_OE_bulk(r, gene_sign = RIR_program, verbose = TRUE)

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
    keep_sig <- c("resF.down", "resF.up")
    score <- as.data.frame(res[, colnames(res) %in% keep_sig])
    score$resF <- score[,"resF.up"] - score[,"resF.down"]
    rownames(score) <- colnames(log2_RNA_tpm)
    colnames(score) <- gsub(".", "_", colnames(score), fixed = TRUE)

    return(data.frame(RIR = score, check.names = FALSE))
}
