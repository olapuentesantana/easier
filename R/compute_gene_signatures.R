#' Compute gold standards
#'
#' \code{computation_gold_standards} computes the scores for the gold standards required by the user
#'
#' @param RNA_tpm numeric matrix with rows=genes and columns=samples
#' @param selected_signatures character string of task names to be considered as gold standards for comparison. (Default scores are computed for: "CYT", "Roh_IS", "chemokines", "Davoli_IS", "IFNy", "Ayers_expIS", "Tcell_inflamed", "RIR", "TLS")
#' @param verbose logical variable indicating whether to display informative messages.
#'
#' @return data.frame containing the gold standards
#' @export
#'
#' @examples
#' # use example dataset from IMvigor210CoreBiologies package (Mariathasan et al., Nature, 2018)
#' data("dataset_mariathasan")
#' gene_tpm <- dataset_mariathasan@tpm
#'
#' # Computation of different hallmarks of the immune response
#' tasks <- c(
#'   "CYT", "Roh_IS", "chemokines", "Davoli_IS", "IFNy",
#'   "Ayers_expIS", "Tcell_inflamed", "RIR", "TLS"
#' )
#' tasks_values <- compute_gold_standards(
#'   RNA_tpm = gene_tpm,
#'   selected_signatures = tasks
#' )
compute_gene_signatures <- function(RNA_tpm, selected_signatures=c("CYT", "Roh_IS", "chemokines", "Davoli_IS", "IFNy", "Ayers_expIS", "Tcell_inflamed", "RIR", "TLS"), verbose = TRUE){
  easier_sigs <- readRDS(file.path(system.file("extdata", "signature_genes.RDS", package = "easier")))

  #check for which selected signatures appropriate functions exist
  sigs <- names(easier_sigs) %in% selected_signatures
  if(verbose) message(c("Following scores can be computed: \n", paste(names(easier_sigs)[sigs], collapse = "\n")))

  result <- lapply(names(easier_sigs)[sigs], function(sig){
    tryCatch(
      {
        if(sig=="Tcell_inflamed"){
          if (any(rownames(RNA_tpm) %in% "C14orf102")){
            message("Gene name changed: NRDE2 is approved symbol, not C14orf102","\n")
            rownames(RNA_tpm)[rownames(RNA_tpm) %in% "C14orf102"] <- "NRDE2"
          }

          # Subset RNA_tpm
          match_genes.housekeeping <- match(easier_sigs$Tcell_inflamed$Housekeeping.read, rownames(RNA_tpm))
          match_genes.predictors <- match(easier_sigs$Tcell_inflamed$Tcell_inflamed.read, rownames(RNA_tpm))

          if (anyNA(c(match_genes.housekeeping, match_genes.predictors))){
            tmp <- c(easier_sigs$Tcell_inflamed$Housekeeping.read, easier_sigs$Tcell_inflamed$Tcell_inflamed.read)
            message(c(paste0("Differenty named or missing signature genes for ",sig,": \n"), paste(tmp[!tmp %in% rownames(RNA_tpm)], collapse = "\n")))
            match_genes.housekeeping <- match_genes.housekeeping[!is.na(match_genes.housekeeping)]
            match_genes.predictors <- match_genes.housekeeping[!is.na(match_genes.housekeeping)]
          }
          do.call(
            paste0("compute_", sig),
            args = list(
              housekeeping = match_genes.housekeeping,
              predictors = match_genes.predictors,
              weights = easier_sigs$Tcell_inflamed$weights,
              RNA_tpm = RNA_tpm,
              verbose = verbose
            )
          )
        } else if (sig=="IMPRES" | sig=="MSI") {
          read <- unique(unlist(easier_sigs[[sig]])) # 15 genes

          # EQUIVALENT : "VISTA" = "C10orf54", "PDL-1" = "CD274", "TIM-3" = "HAVCR2",
          # "PD-1" = "PDCD1", "HVEM" = "TNFRSF14", "OX40L" = "TNFSF4", "CD137L" = "TNFSF9"

          # Some genes might have other name: case for "C10orf54", it's called "VSIR", be carefull
          if (any(rownames(RNA_tpm) %in% "VSIR") & sig == "IMPRES"){
            message("Gene name changed: C10orf54 instead of VSIR","\n")
            rownames(RNA_tpm)[rownames(RNA_tpm) %in% "VSIR"] <- "C10orf54"
          }

          # Some genes might have other name: case for "CCRN4L", it's called "NOCT", be carefull
          if (any(rownames(RNA_tpm) %in% "CCRN4L") & sig == "MSI"){
            message("Gene name changed: NOCT is approved symbol, not CCRN4L","\n")
            rownames(RNA_tpm)[rownames(RNA_tpm) %in% "CCRN4L"] <- "NOCT"
          }

          # Subset RNA_tpm
          match_F_1 <- match(easier_sigs[[sig]]$Gene_1, rownames(RNA_tpm))
          match_F_2 <- match(easier_sigs[[sig]]$Gene_2, rownames(RNA_tpm))

          if (anyNA(c(match_F_1, match_F_2))) {
            message(c("Differenty named or missing signature genes : \n", paste(read[!read %in% rownames(RNA_tpm)], collapse = "\n")))
          }

          do.call(
            "compute_IMPRES_MSI",
            args = list(
              sig = sig,
              len = length(easier_sigs[sig]$Gene_1),
              match_F_1 = match_F_1,
              match_F_2 = match_F_2,
              RNA_tpm = RNA_tpm,
              verbose = verbose
            )
          )

        } else if (sig=="RIR") {
          # TODOTODO add cancertype and file.path as arguments
          do.call(
            paste0("compute_",sig),
            args = list(
              RNA_tpm = RNA_tpm,
              verbose = verbose
            )
          )

        } else {
          # Literature genes
          literature_matches <- match(easier_sigs[[sig]], rownames(RNA_tpm))

          if (anyNA(literature_matches)){
            message(c(paste0("Differenty named or missing signature genes for ",sig,": \n"),
                      paste(easier_sigs[[sig]][!easier_sigs[[sig]] %in% rownames(RNA_tpm)], collapse = "\n")))
            literature_matches <- literature_matches[!is.na(literature_matches)]
          }

          do.call(paste0("compute_", sig), args=list(matches=literature_matches, RNA_tpm=RNA_tpm, verbose=verbose))
        }
      },
      error=function(cond){
        message(paste("The following error occurred while computing sinature of", sig, ":"))
        message(paste(cond, collapse = "/n"))
        df <- data.frame(rep(NA, ncol(RNA_tpm)), row.names = colnames(RNA_tpm))
        names(df)[1] <- sig
        return(df)

      },
      warning=function(cond){
        message(paste("The following warning occurred while computing sinature of", sig, ":"))
        message(paste(cond, collapse ="/n"))
        df <- data.frame(rep(NA, ncol(RNA_tpm)), row.names = colnames(RNA_tpm))
        names(df)[1] <- sig
        return(df)
      }
    )
  })
  return(as.data.frame(result))
}

#' Compute cytolytic activity (CYT) score
#'
#' This function calculates the CYT score as the geometric mean of its signature genes.
#'
#' @references Rooney, M.S., Shukla, S.A., Wu, C.J., Getz, G., and Hacohen, N. (2015). Molecular and genetic properties
#' of tumors associated with local immune cytolytic activity. Cell 160, 48–61. https://doi.org/10.1016/j.cell.2014.12.033.
#'
#' @importFrom stats na.omit
#'
#' @param matches numeric vector indicating the index of signature genes in `RNA_tpm`.
#' @param RNA_tpm data.frame containing TPM values with HGNC symbols in rows and samples in columns.
#' @param verbose logical value indicating whether to display informative messages.
#'
#' @return A numeric matrix with samples in rows and CTY score in a column.
#'
#' @examples
#' # use example dataset from IMvigor210CoreBiologies package (Mariathasan et al., Nature, 2018)
#' data("dataset_mariathasan")
#' gene_tpm <- dataset_mariathasan@tpm
#'
#' # Compute cytolytic activity (Rooney et al, Cell, 2015)
#' CYT <- compute_CYT(RNA_tpm = gene_tpm)
compute_CYT <- function(matches, RNA_tpm, verbose){
  # Subset RNA_tpm
  subset_RNA_tpm <- RNA_tpm[matches, ]

  # Calculation: geometric mean (so-called log-average) [TPM, 0.01 offset]
  score <- as.matrix(apply(subset_RNA_tpm + 0.01, 2, function(X) exp(mean(log(X)))))

  if(verbose) message("CYT score computed")
  return(data.frame(CYT = score, check.names = FALSE))
}

#' Computation of tertiary lymphoid structures signature (TLS) score
#'
#' This function calculates TLS score as the geometric-mean of the expression of its signature genes.
#'
#' @references Cabrita, R., Lauss, M., Sanna, A., Donia, M., Skaarup Larsen, M., Mitra, S., Johansson, I., Phung, B.,
#' Harbst, K., Vallon-Christersson, J., et al. (2020). Tertiary lymphoid structures improve immunotherapy and survival
#' in melanoma. Nature 577, 561–565.
#'
#' @importFrom stats na.omit
#'
#' @param matches numeric vector indicating the index of signature genes in `RNA_tpm`.
#' @param RNA_tpm data.frame containing TPM values with HGNC symbols in rows and samples in columns.
#' @param verbose logical value indicating whether to display informative messages.
#'
#' @return A numeric matrix with samples in rows and TLS score in a column.
#'
#' @examples
#' # use example dataset from IMvigor210CoreBiologies package (Mariathasan et al., Nature, 2018)
#' data("dataset_mariathasan")
#' gene_tpm <- dataset_mariathasan@tpm
#'
#' # compute tertiary lymphoid structures signature (Cabrita et al., Nature, 2020)
#' TLS <- compute_TLS(RNA_tpm = gene_tpm)
compute_TLS <- function(matches, RNA_tpm, verbose){
  # Subset RNA_tpm
  sub_gene.tpm <- RNA_tpm[matches, ]

  # Calculation: geometric mean (so-called log-average) [TPM, 1 offset]
  geom_mean <- apply(sub_gene.tpm, 2, function(X) exp(mean(log2(X + 1))))

  if(verbose) message("TLS score computed")
  return(data.frame(TLS = geom_mean, check.names = FALSE))
}

#' Compute IFNy signature (IFNy) score
#'
#' This function calculates IFNy signature score as the average expression of its signature genes.
#'
#' @references Ayers, M., Lunceford, J., Nebozhyn, M., Murphy, E., Loboda, A., Kaufman, D.R., Albright,
#' A., Cheng, J.D., Kang, S.P., Shankaran, V., et al. (2017). IFN-γ-related mRNA profile predicts clinical
#' response to PD-1 blockade. J. Clin. Invest. 127, 2930–2940. https://doi.org/10.1172/JCI91190.
#'
#' @importFrom stats na.omit
#'
#' @param matches numeric vector indicating the index of signature genes in `RNA_tpm`.
#' @param RNA_tpm data.frame containing TPM values with HGNC symbols in rows and samples in columns.
#' @param verbose logical value indicating whether to display informative messages.
#'
#' @return A numeric matrix with samples in rows and IFNy score in a column.
#'
#' @examples
#' # use example dataset from IMvigor210CoreBiologies package (Mariathasan et al., Nature, 2018)
#' data("dataset_mariathasan")
#' gene_tpm <- dataset_mariathasan@tpm
#'
#' # Compute IFNy signature (Ayers et al., JCI, 2017)
#' IFNy <- compute_IFNy(RNA_tpm = gene_tpm)
compute_IFNy <- function(matches, RNA_tpm, verbose){
  # Log2 transformation:
  log2.RNA_tpm <- log2(RNA_tpm + 1)

  # Subset log2.RNA_tpm
  sub_log2.RNA_tpm <- log2.RNA_tpm[matches, ]

  # Calculation: average of the included genes for the IFN-y signature
  score <- apply(sub_log2.RNA_tpm, 2, mean)

  if(verbose) message("IFNy score computed")
  return(data.frame(IFNy = score, check.names = FALSE))
}

#' Compute Expanded Immune signature (Ayers_expIS) score
#'
#' This function calculates Ayers_expIS score as the average expression of its signature genes.
#'
#' @references Ayers, M., Lunceford, J., Nebozhyn, M., Murphy, E., Loboda, A., Kaufman, D.R., Albright,
#' A., Cheng, J.D., Kang, S.P., Shankaran, V., et al. (2017). IFN-γ-related mRNA profile predicts clinical
#' response to PD-1 blockade. J. Clin. Invest. 127, 2930–2940. https://doi.org/10.1172/JCI91190.
#'
#' @importFrom stats na.omit
#'
#' @param matches numeric vector indicating the index of signature genes in `RNA_tpm`.
#' @param RNA_tpm numeric matrix with rows=genes and columns=samples.
#' @param verbose logical value indicating whether to display informative messages.
#'
#' @return A numeric matrix with rows=samples and columns=Expanded Immune signature score.
#'
#' @examples
#' # use example dataset from IMvigor210CoreBiologies package (Mariathasan et al., Nature, 2018)
#' data("dataset_mariathasan")
#' gene_tpm <- dataset_mariathasan@tpm
#'
#' # Compute expanded immune signature (Ayers et al., JCI, 2017)
#' Ayers_expIS <- compute_Ayers_expIS(gene_tpm)
compute_Ayers_expIS <- function(matches, RNA_tpm, verbose){
  # Log2 transformation:
  log2.RNA_tpm <- log2(RNA_tpm + 1)

  # Subset log2.RNA_tpm
  sub_log2.RNA_tpm  <- log2.RNA_tpm[matches, ]

  # Calculation: average of the included genes for Expanded Immune signature
  score <- apply(sub_log2.RNA_tpm, 2, mean)

  if(verbose) message("Ayers_expIS score computed")
  return(data.frame(Ayers_expIS = score, check.names = FALSE))
}

#' Compute Roh immune score (Roh_IS)
#'
#' This function computes Roh_IS score as the geometric-mean of its signature genes.
#'
#' @references Roh, W., Chen, P.-L., Reuben, A., Spencer, C.N., Prieto, P.A., Miller, J.P., Gopalakrishnan, V.,
#' Wang, F., Cooper, Z.A., Reddy, S.M., et al. (2017). Integrated molecular analysis of tumor biopsies on sequential
#' CTLA-4 and PD-1 blockade reveals markers of response and resistance. Sci. Transl. Med. 9.
#' https://doi.org/10.1126/scitranslmed.aah3560.
#'
#' @importFrom stats na.omit
#'
#' @param matches numeric vector indicating the index of signature genes in `RNA_tpm`.
#' @param RNA_tpm data.frame containing TPM values with HGNC symbols in rows and samples in columns.
#' @param verbose logical value indicating whether to display informative messages.
#'
#' @return A numeric matrix with samples in rows and Roh_IS score in a column.
#'
#' @examples
#' # use example dataset from IMvigor210CoreBiologies package (Mariathasan et al., Nature, 2018)
#' data("dataset_mariathasan")
#' gene_tpm <- dataset_mariathasan@tpm
#'
#' # compute Roh immune score (Roh et al., Sci. Transl. Med., 2017)
#' Roh_IS <- compute_Roh_IS(RNA_tpm = gene_tpm)
compute_Roh_IS <- function(matches, RNA_tpm, verbose){
  # Subset RNA_tpm
  sub_gene.tpm <- RNA_tpm[matches, ]

  # Pseudocount of 0.01 for all genes
  sub_gene.tpm <- sub_gene.tpm + 0.01

  # Pseudocount of 1 for genes with 0 expr
  if(any(sub_gene.tpm == 0)) sub_gene.tpm[sub_gene.tpm == 0] <- sub_gene.tpm[sub_gene.tpm == 0] + 1

  # Calculation: geometric mean (so-called log-average) [TPM, 0.01 offset]
  score <- apply(sub_gene.tpm, 2, function(X) exp(mean(log(X))))

  if(verbose) message("Roh_IS computed score")
  return(data.frame(Roh_IS = score, check.names = FALSE))
}

#' Compute Davoli immune signature (Davoli_IS) score
#'
#' The function calculates Davoli_IS score as the average of the expression of its signature genes after applying rank normalization
#'
#' @references Davoli, T., Uno, H., Wooten, E.C., and Elledge, S.J. (2017). Tumor aneuploidy correlates
#' with markers of immune evasion and with reduced response to immunotherapy. Science 355.
#' https://doi.org/10.1126/science.aaf8399.
#'
#' @importFrom stats na.omit
#'
#' @param matches numeric vector indicating the index of signature genes in `RNA_tpm`.
#' @param RNA_tpm data.frame containing TPM values with HGNC symbols in rows and samples in columns.
#' @param verbose logical value indicating whether to display informative messages
#'
#' @return A numeric matrix with samples in rows and Davoli_IS score in a column.
#'
#' @examples
#' # use example dataset from IMvigor210CoreBiologies package (Mariathasan et al., Nature, 2018)
#' data("dataset_mariathasan")
#' gene_tpm <- dataset_mariathasan@tpm
#'
#' # Compute davoli immune signature (Davoli et al., Science 2017)
#' Davoli_IS <- compute_Davoli_IS(RNA_tpm = gene_tpm)
compute_Davoli_IS <- function(matches, RNA_tpm, verbose){
  # Log2 transformation:
  log2.RNA_tpm <- log2(RNA_tpm + 1)

  # Subset log2.RNA_tpm
  sub_log2.RNA_tpm <- log2.RNA_tpm[matches, ]

  # Calculate rank position for each gene across samples
  ranks_sub_log2.RNA_tpm <- apply(sub_log2.RNA_tpm, 1, rank)

  # Get normalized rank by divided
  ranks_sub_log2.RNA_tpm.norm <- (ranks_sub_log2.RNA_tpm - 1)/(nrow(ranks_sub_log2.RNA_tpm) - 1)

  # Calculation: average of the expression value of all the genes within-sample
  score <- apply(ranks_sub_log2.RNA_tpm.norm, 1, mean)

  if(verbose) message("Davoli_IS score computed")
  return(data.frame(Davoli_IS = score, check.names = FALSE))
}

#' Compute chemokine signature (chemokines) score
#'
#' Computes chemokines score as the PC1 score that results from applying PCA to z-score expression of its signature genes.
#'
#' @references Messina, J.L., Fenstermacher, D.A., Eschrich, S., Qu, X., Berglund, A.E., Lloyd, M.C., Schell, M.J.,
#' Sondak, V.K., Weber, J.S., and Mulé, J.J. (2012). 12-Chemokine gene signature identifies lymph node-like structures
#' in melanoma: potential for patient selection for immunotherapy? Sci. Rep. 2, 765. https://doi.org/10.1038/srep00765.
#'
#' @importFrom stats na.omit prcomp
#'
#' @param matches numeric vector indicating the index of signature genes in `RNA_tpm`.
#' @param RNA_tpm data.frame containing TPM values with HGNC symbols in rows and samples in columns.
#' @param verbose logical value indicating whether to display informative messages.
#'
#' @return A numeric matrix with samples in rows and chemokines score in a column.
#'
#' @examples
#' # use example dataset from IMvigor210CoreBiologies package (Mariathasan et al., Nature, 2018)
#' data("dataset_mariathasan")
#' gene_tpm <- dataset_mariathasan@tpm
#'
#' # Compute chemokine signature (Messina et al., Nat. Sci. Rep., 2012)
#' chemokines <- compute_chemokines(RNA_tpm = gene_tpm)
compute_chemokines <- function(matches, RNA_tpm, verbose){
  # Log2 transformation:
  log2.RNA_tpm <- log2(RNA_tpm + 1)

  # Subset gene_expr
  sub_log2.RNA_tpm <- log2.RNA_tpm[matches, ]

  # calculation: using PCA (Z-score calculated within prcomp)
  chemokine.pca <- stats::prcomp(t(sub_log2.RNA_tpm), center = TRUE, scale = TRUE)
  score <- chemokine.pca$x[, 1]

  if(verbose) message("Chemokines score computed")
  return(data.frame(chemokines = score, check.names = FALSE))
}

#' Compute T cell-inflamed signature (Tcell_inflamed) score
#'
#' This function calculates Tcell_inflamed score as a weighted sum of housekeeping normalized expression of its signature genes.
#' Weightes were available at Table S2B from Cristescu R, et al. Pan-tumor genomic biomarkers for PD-1 checkpoint
#' blockade-based immunotherapy. Science. (2018) 362:eaar3593. doi: 10.1126/science.aar3593.
#'
#' @references Ayers, M., Lunceford, J., Nebozhyn, M., Murphy, E., Loboda, A., Kaufman, D.R., Albright,
#' A., Cheng, J.D., Kang, S.P., Shankaran, V., et al. (2017). IFN-γ-related mRNA profile predicts clinical
#' response to PD-1 blockade. J. Clin. Invest. 127, 2930–2940. https://doi.org/10.1172/JCI91190.
#'
#' @importFrom stats na.omit
#'
#' @param housekeeping numeric vector indicating the index of houskeeping genes in `RNA_tpm`.
#' @param predictors numeric vector indicating the index of predictor genes in `RNA_tpm`.
#' @param weights numeric vector containing the weights.
#' @param RNA_tpm data.frame containing TPM values with HGNC symbols in rows and samples in columns.
#' @param verbose logical value indicating whether to display informative messages
#'
#' @return A numeric matrix with samples in rows and Tcell_inflamed score in a column.
#'
#' @examples
#' # use example dataset from IMvigor210CoreBiologies package (Mariathasan et al., Nature, 2018)
#' data("dataset_mariathasan")
#' gene_tpm <- dataset_mariathasan@tpm
#'
#' # compute T-cell inflamed signature (Ayers et al., JCI, 2017)
#' Tcell_inflamed <- compute_Tcell_inflamed(RNA_tpm = gene_tpm)
compute_Tcell_inflamed <- function(housekeeping, predictors, weights, RNA_tpm, verbose){
  # Log2 transformation:
  log2.RNA_tpm <- log2(RNA_tpm + 1)

  # Subset log2.RNA_tpm
  ## housekeeping
  log2.RNA_tpm.housekeeping <- log2.RNA_tpm[housekeeping, ]
  ## predictors
  log2.RNA_tpm.predictors <- log2.RNA_tpm[predictors, ]
  weights <- weights[rownames(log2.RNA_tpm.predictors)]

  # Housekeeping normalization
  average.log2.RNA_tpm.housekeeping <- apply(log2.RNA_tpm.housekeeping, 2, mean)
  log2.RNA_tpm.predictors.norm <- sweep(log2.RNA_tpm.predictors, 2, average.log2.RNA_tpm.housekeeping, FUN = "-")

  # Calculation: weighted sum of the normalized predictor gene values
  tidy <- match(rownames(log2.RNA_tpm.predictors.norm), names(weights))

  #transform vector to matrix
  weights <- matrix(weights, ncol = 1, dimnames = list(names(weights)))
  score <- t(log2.RNA_tpm.predictors.norm[tidy,]) %*% weights

  if(verbose) message("Tcell_inflamed score computed")
  return(data.frame( Tcell_inflamed = score, check.names = FALSE))
}

#' Compute Immuno-Predictive Score (IMPRES)
#'
#' This function calculates IMPRES score by logical comparison of checkpoint gene pairs expression.
#'
#' @references Auslander,N.,Zhang,G.,Lee,J.S.,Frederick,D.T.,Miao,B.,Moll,T.,Tian,T.,Wei,Z., Madan, S.,
#' Sullivan, R.J., et al. (2018). Robust prediction of response to immune checkpoint blockade therapy in
#' metastatic melanoma. Nat. Med. 24, 1545–1549. https://doi.org/10.1038/s41591-018-0157-9.
#'
#' @importFrom stats na.omit
#'
#' @param sig can be either 'IMPRES' or 'MSI'.
#' @param len the length of gene_1 vector.
#' @param match_F_1 numeric vector indicating the index of signature genes defined in 'gene_1' in `RNA_tpm`.
#' @param match_F_2 numeric vector indicating the index of signature genes defined in 'gene_2' in `RNA_tpm`.
#' @param RNA_tpm data.frame containing TPM values with HGNC symbols in rows and samples in columns.
#' @param verbose logical value indicating whether to display informative messages.
#'
#' @return A numeric matrix with samples in rows and IMPRES score in a column.
#'
#' @examples
#' # use example dataset from IMvigor210CoreBiologies package (Mariathasan et al., Nature, 2018)
#' data("dataset_mariathasan")
#' gene_tpm <- dataset_mariathasan@tpm
#'
#' # Compute IMPRES signature (Auslander et al., Nat.Med., 2018)
#' IMPRES <- compute_IMPRES(RNA_tpm = gene_tpm)
compute_IMPRES_MSI <- function(sig, len, match_F_1, match_F_2, RNA_tpm, verbose){
  # Initialize variables
  F_pair_expr_A <- F_pair_expr_B <- IMPRES.matrix <- matrix(0, len, ncol(RNA_tpm))
  colnames(IMPRES.matrix) <- colnames(RNA_tpm)
  score <- vector("numeric", length = ncol(RNA_tpm))
  names(score) <- colnames(RNA_tpm)

  # Log2 transformation:
  log2.RNA_tpm <- as.data.frame(log2(RNA_tpm + 1))

  # Calculation:
  F_pair_expr_A <- log2.RNA_tpm[match_F_1, ]
  F_pair_expr_B <- log2.RNA_tpm[match_F_2, ]

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

  if(verbose) message(paste(sig,"score computed"))
  df <- data.frame(score, check.names = FALSE)
  names(df)[1] <- sig

  return(df)
}


