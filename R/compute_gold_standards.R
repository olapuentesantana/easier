#' Compute gold standards
#'
#' `computation_gold_standards` computes the scores for the gold standards
#' required by the user
#'
#' @export
#'
#' @param RNA_tpm numeric matrix with data, as tpm values
#' @param list_gold_standards string with gold standards names
#' @param cancertype string character
#' @param output_file_path TODOTODO to define - maybe name it more specifically to the function it absolves?
#'
#' @return List with the scores of all the gold standards specified.
#'
#' @examples
#' # Example: Riaz
#' data("Riaz_data")
#'
#' Computation of different hallmarks of the immune response
#' tasks <- c("CYT", "Roh_IS", "chemokines", "Davoli_IS", "IFNy", "Ayers_expIS", "Tcell_inflamed", "RIR", "TLS")
#' immune_response <- compute_gold_standards(RNA_tpm = Riaz_data$tpm_RNAseq,
#' list_gold_standards = tasks,
#' cancertype = "SKCM",
#' output_file_path = "/Users/Oscar/Desktop/Riaz")
compute_gold_standards <- function(RNA_tpm,
                                   list_gold_standards,
                                   cancertype,
                                   output_file_path) {

  # calculate Immune Checkpoint genes expression #
  ICB_genes <- compute_ICB_genes(RNA_tpm)

  gold.standards <- sapply(list_gold_standards, function(X) {
    if ("CYT" == X) {

      # calculate Cytolytic activity #
      CYT <- t(compute_CYT(RNA_tpm))
      return(list(CYT))
    } else if ("IPS" == X) {

      # calculate Immunophenoscore #
      IPS <- t(compute_IPS(RNA_tpm))
      return(list(IPS))
    } else if ("IMPRES" == X) {

      # calculate Impres #
      IMPRES <- t(compute_IMPRES(RNA_tpm))
      return(list(IMPRES))
    } else if ("Roh_IS" == X) {

      # calculate roh immune signature #
      Roh_IS <- t(compute_Roh_IS(RNA_tpm))
      return(list(Roh_IS))
    } else if ("chemokines" == X) {

      # calculate chemokine signature #
      chemokines <- t(compute_chemokines(RNA_tpm))
      return(list(chemokines))
    } else if ("Davoli_IS" == X) {

      # calculate davoli cytotoxic immune signature #
      Davoli_IS <- t(compute_Davoli_IS(RNA_tpm))
      return(list(Davoli_IS))
    } else if ("IFNy" == X) {

      # calculate ayers IFNy #
      IFNy <- t(compute_IFNy(RNA_tpm))
      return(list(IFNy))
    } else if ("Ayers_expIS" == X) {

      # calculate ayers expanded immune signature #
      Ayers_expIS <- t(compute_Ayers_expIS(RNA_tpm))
      return(list(Ayers_expIS))
    } else if ("Tcell_inflamed" == X) {

      # calculate ayers T cell inflamed signature #
      Tcell_inflamed <- t(compute_Tcell_inflamed(RNA_tpm))
      return(list(Tcell_inflamed))
    } else if ("TIDE" == X) {

      # calculate TIDE signature #
      TIDE <- t(compute_TIDE(RNA_tpm, cancertype, output_file_path))
      return(list(TIDE))
    } else if ("MSI" == X) {

      # calculate MSI signature #
      MSI <- t(compute_MSI(RNA_tpm))
      return(list(MSI))
    } else if ("RIR" == X) {

      # calculate MSI signature #
      RIR <- t(compute_RIR(RNA_tpm))
      return(list(RIR))
    } else if ("TLS" == X) {

      # calculate MSI signature #
      TLS <- t(compute_TLS(RNA_tpm))
      return(list(TLS))
    } else if ("CTLA4" == X) {
      CTLA4 <- ICB_genes$CTLA4
      return(list(CTLA4))
    } else if ("PD1" == X) {
      PD1 <- ICB_genes$PD1
      return(list(PD1))
    } else if ("PDL1" == X) {
      PDL1 <- ICB_genes$PDL1
      return(list(PDL1))
    }
  })

  return(gold.standards)
}
