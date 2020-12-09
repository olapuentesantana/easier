#' Compute gold standards
#'
#' \code{computation_gold_standards} computes the scores for the gold standards required by the user
#'
#' @export
#'
#' @param RNA.tpm numeric matrix with data, as tpm values
#' @param list_gold_standards string with gold standards names
#' @param cancertype string character
#'
#' @return List with the scores of all the gold standards specified.
#'
compute_gold_standards <- function(RNA.tpm,
                                   list_gold_standards,
                                   cancertype,
                                   output_file_path){

  # calculate Immune Checkpoint genes expression #
  ICB_genes <- compute_ICB_genes(RNA.tpm)

  gold.standards <- sapply(list_gold_standards, function(X){

    if ("CYT" == X) {

      # calculate Cytolytic activity #
      CYT <- t(compute.CYT(RNA.tpm))
      return(list(CYT))

    }else if("IPS" == X) {

      # calculate Immunophenoscore #
      IPS <- t(compute.IPS(RNA.tpm))
      return(list(IPS))

    }else if("IMPRES" == X) {

      # calculate Impres #
      IMPRES <- t(compute.IMPRES(RNA.tpm))
      return(list(IMPRES))

    }else if("Roh_IS" == X) {

      # calculate roh immune signature #
      Roh_IS <- t(compute.Roh_IS(RNA.tpm))
      return(list(Roh_IS))

    }else if("chemokines" == X) {

      # calculate chemokine signature #
      chemokines <- t(compute.chemokines(RNA.tpm))
      return(list(chemokines))

    }else if("Davoli_IS" == X) {

      # calculate davoli cytotoxic immune signature #
      Davoli_IS <- t(compute.Davoli_IS(RNA.tpm))
      return(list(Davoli_IS))

    }else if("IFNy" == X) {

      # calculate ayers IFNy #
      IFNy <- t(compute.IFNy(RNA.tpm))
      return(list(IFNy))

    }else if("Ayers_expIS" == X) {

      # calculate ayers expanded immune signature #
      Ayers_expIS <- t(compute.Ayers_expIS(RNA.tpm))
      return(list(Ayers_expIS))

    }else if("Tcell_inflamed" == X) {

      # calculate ayers T cell inflamed signature #
      Tcell_inflamed <- t(compute.Tcell_inflamed(RNA.tpm))
      return(list(Tcell_inflamed))

    }else if("TIDE" == X) {

      # calculate TIDE signature #
      TIDE <- t(compute.TIDE(RNA.tpm, cancertype, output_file_path))
      return(list(TIDE))

    }else if("MSI" == X) {

      # calculate MSI signature #
      MSI <- t(compute.MSI(RNA.tpm))
      return(list(MSI))

    }else if("RIR" == X) {

      # calculate MSI signature #
      RIR <- t(compute.RIR(RNA.tpm))
      return(list(RIR))

    }else if("TLS" == X) {

      # calculate MSI signature #
      TLS <- t(compute.TLS(RNA.tpm))
      return(list(TLS))

    }else if("CTLA4" == X) {

      CTLA4 <- ICB_genes$CTLA4
      return(list(CTLA4))

    }else if("PD1" == X) {

      PD1 <- ICB_genes$PD1
      return(list(PD1))

    }else if("PDL1" == X) {

      PDL1 <- ICB_genes$PDL1
      return(list(PDL1))

    }

  })

  return(gold.standards)

}

