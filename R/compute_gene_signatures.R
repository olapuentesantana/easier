selected_signatures <- c("CYT", "IPS", "IMPRES", "Roh_IS", "chemokines", "Davoli_IS", "IFNy",
  "Ayers_expIS", "Tcell_inflamed", "TIDE", "MSI", "TLS")

#'
#'
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#'
#' @return
#' @export
#'
#' @examples
#' # TODOTODO
compute_gene_signatures <- function(RNA.tpm, selected_signatures){

  easier_sigs <- list(CYT=c("GZMA", "PRF1"),
                      TLS=c("CD79B", "CD1D", "CCR6", "LAT", "SKAP1", "CETP", "EIF1AY", "RBP5", "PTGDS"),
                      IFNy=c("IFNG", "STAT1", "CXCL9", "CXCL10", "IDO1", "HLA-DRA"),
                      Ayers_expIS=c("GZMB", "GZMK", "CXCR6", "CXCL10", "CXCL13", "CCL5", "STAT1","CD3D", "CD3E",
                                    "CD2", "IL2RG" , "NKG7", "HLA-E", "CIITA","HLA-DRA", "LAG3", "IDO1", "TAGAP"),
                      Tcell_inflamed=list(
                        Tcell_inflamed.read=c("CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6", "HLA-DQA1",
                                              "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2", "PSMB10", "STAT1", "TIGIT"),
                        Housekeeping.read=c("STK11IP", "ZBTB34", "TBC1D10B", "OAZ1", "POLR2A", "G6PD", "ABCF1", "NRDE2", "UBB", "TBP", "SDHA"),
                        weights=c(CCL5=0.008346, CD27=0.072293, CD274=0.042853, CD276=-0.0239, CD8A=0.031021 ,CMKLR1=0.151253, CXCL9=0.074135,
                                  CXCR6=0.004313, `HLA-DQA1`=0.020091, `HLA-DRB1`=0.058806, `HLA-E`=0.07175, IDO1=0.060679, LAG3=0.123895, NKG7=0.075524, PDCD1LG2=0.003734,
                                  PSMB10=0.032999, STAT1=0.250229, TIGIT=0.084767)),
                      Roh_IS=c("GZMA", "GZMB", "PRF1", "GNLY", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F",
                               "HLA-G", "HLA-H", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DPA1",
                               "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DRA", "HLA-DRB1",
                               "IFNG", "IFNGR1", "IFNGR2", "IRF1", "STAT1", "PSMB9", "CCR5", "CCL3", "CCL4",
                               "CCL5", "CXCL9", "CXCL10", "CXCL11", "ICAM1", "ICAM2", "ICAM3", "ICAM4", "ICAM5", "VCAM1"),
                      Davoli_IS=c("CD247", "CD2", "CD3E", "GZMH", "NKG7", "PRF1", "GZMK"),
                      chemokines=c("CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21",
                                   "CXCL9", "CXCL10", "CXCL11", "CXCL13"))

  sigs <- names(easier_sigs) %in% selected_signatures
  message(c("Following scores can be computed: \n", paste(names(easier_sigs)[sigs], collapse = "\n")))

  result <- lapply(names(easier_sigs)[sigs], function(sig){

    if(sig=="Tcell_inflamed"){
      if (any(rownames(RNA.tpm) %in% "C14orf102")){
        cat("Gene name changed: NRDE2 is approved symbol, not C14orf102","\n")
        rownames(RNA.tpm)[rownames(RNA.tpm) %in% "C14orf102"] <- "NRDE2"
      }

      match_genes.housekeeping <- match(easier_sigs$Tcell_inflamed$Housekeeping.read, rownames(RNA.tpm))
      match_genes.predictors <- match(easier_sigs$Tcell_inflamed$Tcell_inflamed.read, rownames(RNA.tpm))

      if (anyNA(c(match_genes.housekeeping, match_genes.predictors))){
        tmp <- c(Tcell_inflamed.read, Housekeeping.read)
        warning(c(paste0("Differenty named or missing signature genes for ",sig,": \n"), paste(tmp[!tmp %in% rownames(RNA.tpm)], collapse = "\n")))
        match_genes.housekeeping <- match_genes.housekeeping[!is.na(match_genes.housekeeping)]
        match_genes.predictors <- match_genes.housekeeping[!is.na(match_genes.housekeeping)]
      }
      do.call(
        paste0("compute_", sig),
        args = list(
          housekeeping = match_genes.housekeeping,
          predictors = match_genes.predictors,
          weights = easier_sigs$Tcell_inflamed$weights,
          RNA.tpm = RNA.tpm
        )
      )
    } else {
      # Literature genes
      literature_matches <- match(easier_sigs[[sig]], rownames(RNA.tpm))

      if (anyNA(literature_matches)){
        warning(c(paste0("Differenty named or missing signature genes for ",sig,": \n"), paste(easier_sigs[[sig]][!easier_sigs[[sig]] %in% rownames(RNA.tpm)], collapse = "\n")), immediate. = TRUE)
        literature_matches <- literature_matches[!is.na(literature_matches)]
      }

      do.call(paste0("compute_", sig), args=list(matches=literature_matches, RNA.tpm=RNA.tpm))
    }
  })
  as.data.frame(result)
}

compute_gene_signatures(tpm, selected_signatures)


compute_CYT <- function(matches, RNA.tpm){
  # Subset RNA.tpm
  subset_RNA.tpm <- RNA.tpm[matches, ]

  # Calculation: geometric mean (so-called log-average) [TPM, 0.01 offset]
  score <- as.matrix(apply(subset_RNA.tpm + 0.01, 2, function(X) exp(mean(log(X)))))

  message("CYT score computed")
  return(data.frame(CYT = score, check.names = FALSE))
}

compute_TLS <- function(matches, RNA.tpm){
  # Subset RNA.tpm
  sub_gene.tpm <- RNA.tpm[matches, ]

  # Calculation: geometric mean (so-called log-average) [TPM, 1 offset]
  geom_mean <- apply(sub_gene.tpm, 2, function(X) exp(mean(log2(X + 1))))

  message("TLS score computed")
  return(data.frame(TLS = geom_mean, check.names = FALSE))
}

compute_IFNy <- function(matches, RNA.tpm){
  # Log2 transformation:
  log2.RNA.tpm <- log2(RNA.tpm + 1)

  # Subset log2.RNA.tpm
  sub_log2.RNA.tpm <- log2.RNA.tpm[matches, ]

  # Calculation: average of the included genes for the IFN-y signature
  score <- apply(sub_log2.RNA.tpm, 2, mean)

  message("IFNy score computed")
  return(data.frame(IFNy = score, check.names = FALSE))
}

compute_Ayers_expIS <- function(matches, RNA.tpm){
  # Log2 transformation:
  log2.RNA.tpm <- log2(RNA.tpm + 1)

  # Subset log2.RNA.tpm
  sub_log2.RNA.tpm  <- log2.RNA.tpm[matches, ]

  # Calculation: average of the included genes for Expanded Immune signature
  score <- apply(sub_log2.RNA.tpm, 2, mean)

  message("Ayers_expIS score computed")
  return(data.frame(Ayers_expIS = score, check.names = FALSE))
}

compute_Roh_IS <- function(matches, RNA.tpm){
  # Subset RNA.tpm
  sub_gene.tpm <- RNA.tpm[matches, ]

  # Pseudocount of 0.01 for all genes
  sub_gene.tpm <- sub_gene.tpm + 0.01

  # Pseudocount of 1 for genes with 0 expr
  if(any(sub_gene.tpm == 0)) sub_gene.tpm[sub_gene.tpm == 0] <- sub_gene.tpm[sub_gene.tpm == 0] + 1

  # Calculation: geometric mean (so-called log-average) [TPM, 0.01 offset]
  score <- apply(sub_gene.tpm, 2, function(X) exp(mean(log(X))))

  message("Roh_IS computed score")
  return(data.frame(Roh_IS = score, check.names = FALSE))
}

compute_Davoli_IS <- function(matches, RNA.tpm){
  # Log2 transformation:
  log2.RNA.tpm <- log2(RNA.tpm + 1)

  # Subset log2.RNA.tpm
  sub_log2.RNA.tpm <- log2.RNA.tpm[matches, ]

  # Calculate rank position for each gene across samples
  ranks_sub_log2.RNA.tpm <- apply(sub_log2.RNA.tpm, 1, rank)

  # Get normalized rank by divided
  ranks_sub_log2.RNA.tpm.norm <- (ranks_sub_log2.RNA.tpm - 1)/(nrow(ranks_sub_log2.RNA.tpm) - 1)

  # Calculation: average of the expression value of all the genes within-sample
  score <- apply(ranks_sub_log2.RNA.tpm.norm, 1, mean)

  message("Davoli_IS score computed")
  return(data.frame(Davoli_IS = score, check.names = FALSE))
}

compute_chemokines <- function(matches, RNA.tpm){
  # Log2 transformation:
  log2.RNA.tpm <- log2(RNA.tpm + 1)

  # Subset gene_expr
  sub_log2.RNA.tpm <- log2.RNA.tpm[matches, ]

  # calculation: using PCA (Z-score calculated within prcomp)
  chemokine.pca <- stats::prcomp(t(sub_log2.RNA.tpm), center = TRUE, scale = TRUE)
  score <- chemokine.pca$x[, 1]

  message("Chemokines score computed")
  return(data.frame(chemokines = score, check.names = FALSE))
}

compute_Tcell_inflamed <- function(housekeeping, predictors, weights, RNA.tpm){
  # Log2 transformation:
  log2.RNA.tpm <- log2(RNA.tpm + 1)

  # Subset log2.RNA.tpm
  ## housekeeping
  log2.RNA.tpm.housekeeping <- log2.RNA.tpm[housekeeping, ]
  ## predictors
  log2.RNA.tpm.predictors <- log2.RNA.tpm[predictors, ]
  weights <- weights[rownames(log2.RNA.tpm.predictors)]

  # Housekeeping normalization
  average.log2.RNA.tpm.housekeeping <- apply(log2.RNA.tpm.housekeeping, 2, mean)
  log2.RNA.tpm.predictors.norm <- sweep(log2.RNA.tpm.predictors, 2, average.log2.RNA.tpm.housekeeping, FUN = "-")

  # Calculation: weighted sum of the normalized predictor gene values
  tidy <- match(rownames(log2.RNA.tpm.predictors.norm), names(weights))

  #transform vector to matrix
  weights <- matrix(weights, ncol = 1, dimnames = list(names(weights)))
  score <- t(log2.RNA.tpm.predictors.norm[tidy,]) %*% weights

  message("Tcell_inflamed score computed")
  return(data.frame( Tcell_inflamed = score, check.names = FALSE))
}


as.data.frame(c(myvec))



df1 <- compute_CYT(c(8614, 14408), tpm)

df2 <- compute_TLS(c(177, 10077, 4034, 17645, 7978, 1720,16867, 7686, 3381), tpm)

df3 <- compute_IFNy(c(3872,  4389,  7624, 12340,  6340, 17148), tpm)

dfs <- list(df1, df2, df3)

as.data.frame(lapply(dfs, function(df){return(df)}))


