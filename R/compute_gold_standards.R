#' Compute gold standards (i.e. tasks)
#'
#' `computation_gold_standards` computes the scores for the gold standards required by the user
#'
#' @export
#'
#' @param RNA_tpm A numeric matrix of patients' gene expression data as tpm values.
#' @param list_gold_standards A character string of task names to be considered as gold standards for comparison.
#'
#' @return A numeric matrix of patients' gold standard values (rows = samples; columns = tasks).
#'
#' @examples
#' # use example dataset from Mariathasan cohort (Mariathasan et al., Nature, 2018)
#' data(cds)
#' mariathasan_data <- preprocess_mariathasan(cds)
#' gene_tpm <- mariathasan_data$tpm
#' rm(cds)
#'
#' # Computation of different hallmarks of the immune response
#' tasks <- c("CYT", "Roh_IS", "chemokines", "Davoli_IS", "IFNy",
#'            "Ayers_expIS", "Tcell_inflamed", "RIR", "TLS")
#' tasks_values <- compute_gold_standards(
#'   RNA_tpm = gene_tpm,
#'   list_gold_standards = tasks)
compute_gold_standards <- function(RNA_tpm,
                                   list_gold_standards = c("CYT", "Roh_IS", "chemokines", "Davoli_IS", "IFNy", "Ayers_expIS", "Tcell_inflamed", "RIR", "TLS")) {

  gold_standards <- sapply(list_gold_standards, function(X) {
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
  rownames(gold_standards) <- colnames(RNA_tpm)
  return(gold_standards)
}
