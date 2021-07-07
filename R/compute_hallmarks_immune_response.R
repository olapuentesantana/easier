#' Compute published transcriptomics-based scores of hallmarks of anti-cancer immune response
#' (so-called gold standards)
#'
#' The function computes the gold standard scores required by the user.
#'
#' @export
#'
#' @param RNA_tpm data.frame containing TPM values with HGNC symbols in rows and samples in columns.
#' @param list_hallmarks_of_immune_response character string of task names to be considered as gold standards for comparison.
#'
#' @return A numeric matrix with samples in rows and gold standard scores in columns.
#'
#' @examples
#' # use example dataset from IMvigor210CoreBiologies package (Mariathasan et al., Nature, 2018)
#' data("dataset_mariathasan")
#' gene_tpm <- dataset_mariathasan@tpm
#'
#' # Computation of different hallmarks of the immune response
#' list_hallmarks_of_immune_response <- c(
#'   "CYT", "Roh_IS", "chemokines", "Davoli_IS", "IFNy",
#'   "Ayers_expIS", "Tcell_inflamed", "RIR", "TLS"
#' )
#' tasks_values <- compute_hallmarks_immune_response(
#'   RNA_tpm = gene_tpm,
#'   list_hallmarks_of_immune_response = list_hallmarks_of_immune_response
#' )
compute_hallmarks_immune_response <- function(RNA_tpm,
                                              list_hallmarks_of_immune_response = c("CYT", "Roh_IS", "chemokines", "Davoli_IS", "IFNy", "Ayers_expIS", "Tcell_inflamed", "RIR", "TLS")) {
  hallmarks_immune_response <- sapply(list_hallmarks_of_immune_response, function(X) {
    if ("CYT" == X) {

      # calculate Cytolytic activity #
      CYT <- t(compute_CYT(RNA_tpm))
      return(CYT)
    } else if ("IMPRES" == X) {

      # calculate Impres #
      IMPRES <- t(compute_IMPRES(RNA_tpm))
      return(IMPRES)
    } else if ("Roh_IS" == X) {

      # calculate roh immune signature #
      Roh_IS <- t(compute_Roh_IS(RNA_tpm))
      return(Roh_IS)
    } else if ("chemokines" == X) {

      # calculate chemokine signature #
      chemokines <- t(compute_chemokines(RNA_tpm))
      return(chemokines)
    } else if ("Davoli_IS" == X) {

      # calculate davoli cytotoxic immune signature #
      Davoli_IS <- t(compute_Davoli_IS(RNA_tpm))
      return(Davoli_IS)
    } else if ("IFNy" == X) {

      # calculate ayers IFNy #
      IFNy <- t(compute_IFNy(RNA_tpm))
      return(IFNy)
    } else if ("Ayers_expIS" == X) {

      # calculate ayers expanded immune signature #
      Ayers_expIS <- t(compute_Ayers_expIS(RNA_tpm))
      return(Ayers_expIS)
    } else if ("Tcell_inflamed" == X) {

      # calculate ayers T cell inflamed signature #
      Tcell_inflamed <- t(compute_Tcell_inflamed(RNA_tpm))
      return(Tcell_inflamed)
    } else if ("MSI" == X) {

      # calculate MSI signature #
      MSI <- t(compute_MSI(RNA_tpm))
      return(MSI)
    } else if ("RIR" == X) {

      # calculate MSI signature #
      RIR <- t(compute_RIR(RNA_tpm))
      return(RIR)
    } else if ("TLS" == X) {

      # calculate MSI signature #
      TLS <- t(compute_TLS(RNA_tpm))
      return(TLS)
    }
  })
  rownames(hallmarks_immune_response) <- colnames(RNA_tpm)
  return(hallmarks_immune_response)
}
