#' Compute gold standards
#'
#' \code{computation_gold_standards} computes the scores for the gold standards required by the user
#'
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#' @param selected_signatures TODOTODO
#'
#' @return TODOTODO
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
                                   "CXCL9", "CXCL10", "CXCL11", "CXCL13"),
                      IMPRES=list(
                        Gene_1 = c("PDCD1", "CD27", "CTLA4", "CD40", "CD86", "CD28", "CD80",
                                   "CD274", "CD86", "CD40", "CD86", "CD40", "CD28", "CD40", "TNFRSF14"),
                        Gene_2 = c("TNFSF4", "PDCD1", "TNFSF4", "CD28", "TNFSF4", "CD86", "TNFSF9",
                                   "C10orf54", "HAVCR2", "PDCD1", "CD200", "CD80", "CD276", "CD274", "CD86")),
                      MSI=list(
                        Gene_1 = c("HNRNPL", "MTA2", "CALR", "RASL11A", "LYG1", "STRN3", "HPSE", "PRPF39", "NOCT","AMFR"),
                        Gene_2 = c("CDC16", "VGF", "SEC22B", "CAB39L", "DHRS12", "TMEM192", "BCAS3", "ATF6","GRM8","DUSP18")),
                      ICB_genes = c("CD274","CTLA4","PDCD1"),
                      TIDE=c(),
                      IPS=c(),
                      RIR=c())

  #check for which selected signatures appropriate functions exist
  sigs <- names(easier_sigs) %in% selected_signatures
  message(c("Following scores can be computed: \n", paste(names(easier_sigs)[sigs], collapse = "\n")))

  result <- lapply(names(easier_sigs)[sigs], function(sig){
    tryCatch(
      {
        if(sig=="Tcell_inflamed"){
          if (any(rownames(RNA.tpm) %in% "C14orf102")){
            message("Gene name changed: NRDE2 is approved symbol, not C14orf102","\n")
            rownames(RNA.tpm)[rownames(RNA.tpm) %in% "C14orf102"] <- "NRDE2"
          }

          # Subset RNA.tpm
          match_genes.housekeeping <- match(easier_sigs$Tcell_inflamed$Housekeeping.read, rownames(RNA.tpm))
          match_genes.predictors <- match(easier_sigs$Tcell_inflamed$Tcell_inflamed.read, rownames(RNA.tpm))

          if (anyNA(c(match_genes.housekeeping, match_genes.predictors))){
            tmp <- c(easier_sigs$Tcell_inflamed$Housekeeping.read, easier_sigs$Tcell_inflamed$Tcell_inflamed.read)
            message(c(paste0("Differenty named or missing signature genes for ",sig,": \n"), paste(tmp[!tmp %in% rownames(RNA.tpm)], collapse = "\n")))
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
        } else if (sig=="IMPRES" | sig=="MSI") {
          read <- unique(unlist(easier_sigs[[sig]])) # 15 genes

          # EQUIVALENT : "VISTA" = "C10orf54", "PDL-1" = "CD274", "TIM-3" = "HAVCR2",
          # "PD-1" = "PDCD1", "HVEM" = "TNFRSF14", "OX40L" = "TNFSF4", "CD137L" = "TNFSF9"

          # Some genes might have other name: case for "C10orf54", it's called "VSIR", be carefull
          if (any(rownames(RNA.tpm) %in% "VSIR") & sig == "IMPRES"){
            message("Gene name changed: C10orf54 instead of VSIR","\n")
            rownames(RNA.tpm)[rownames(RNA.tpm) %in% "VSIR"] <- "C10orf54"
          }

          # Some genes might have other name: case for "CCRN4L", it's called "NOCT", be carefull
          if (any(rownames(RNA.tpm) %in% "CCRN4L") & sig == "MSI"){
            message("Gene name changed: NOCT is approved symbol, not CCRN4L","\n")
            rownames(RNA.tpm)[rownames(RNA.tpm) %in% "CCRN4L"] <- "NOCT"
          }

          # Subset RNA.tpm
          match_F_1 <- match(easier_sigs[[sig]]$Gene_1, rownames(RNA.tpm))
          match_F_2 <- match(easier_sigs[[sig]]$Gene_2, rownames(RNA.tpm))

          if (anyNA(c(match_F_1, match_F_2))) {
            message(c("Differenty named or missing signature genes : \n", paste(read[!read %in% rownames(RNA.tpm)], collapse = "\n")))
          }

          do.call(
            "compute_IMPRES_MSI",
            args = list(
              sig = sig,
              len = length(easier_sigs[sig]$Gene_1),
              match_F_1 = match_F_1,
              match_F_2 = match_F_2,
              RNA.tpm = RNA.tpm
            )
          )

        } else if (sig=="ICB_genes") {
          do.call(
            "compute_ICB_genes",
            args = list(
              genes = easier_sigs[[sig]],
              RNA.tpm = RNA.tpm
            )
          )

        } else if (sig=="TIDE" | sig=="IPS" | sig=="RIR") {
          # TODOTODO add cancertype and file.path as arguments
          do.call(
            paste0("compute.",sig),
            args = list(
              RNA.tpm = RNA.tpm
            )
          )

        } else {
          # Literature genes
          literature_matches <- match(easier_sigs[[sig]], rownames(RNA.tpm))

          if (anyNA(literature_matches)){
            message(c(paste0("Differenty named or missing signature genes for ",sig,": \n"),
                      paste(easier_sigs[[sig]][!easier_sigs[[sig]] %in% rownames(RNA.tpm)], collapse = "\n")))
            literature_matches <- literature_matches[!is.na(literature_matches)]
          }

          do.call(paste0("compute_", sig), args=list(matches=literature_matches, RNA.tpm=RNA.tpm))
        }
      },
      error=function(cond){
        message(paste("The following error occurred while computing sinature of", sig, ":"))
        message(paste(cond, collapse = "/n"))
        df <- data.frame(rep(NA, ncol(RNA.tpm)), row.names = colnames(RNA.tpm))
        names(df)[1] <- sig
        return(df)

      },
      warning=function(cond){
        message(paste("The following warning occurred while computing sinature of", sig, ":"))
        message(paste(cond, collapse ="/n"))
        df <- data.frame(rep(NA, ncol(RNA.tpm)), row.names = colnames(RNA.tpm))
        names(df)[1] <- sig
        return(df)
      }
    )
  })
  return(as.data.frame(result))
}

#' Compute cytolytic activity score
#'
#' \code{compute_CYT} computes cytolytic activity score as the geometric mean of immune cytolytic genes
#' (Rooney et al., 2015).
#'
#' @importFrom stats na.omit
#' @param matches TODOTODO
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=cytolytic activity score
#'
#' @export
#'
#' @examples
#' # TODOTODO
compute_CYT <- function(matches, RNA.tpm){
  # Subset RNA.tpm
  subset_RNA.tpm <- RNA.tpm[matches, ]

  # Calculation: geometric mean (so-called log-average) [TPM, 0.01 offset]
  score <- as.matrix(apply(subset_RNA.tpm + 0.01, 2, function(X) exp(mean(log(X)))))

  message("CYT score computed")
  return(data.frame(CYT = score, check.names = FALSE))
}

#' Compute tertiary lymphoid structures signature
#'
#' \code{compute_TLS} computes TLS signature as the geometric-mean of TLS signature genes
#' (Cabrita et al., 2020).
#'
#' @importFrom stats na.omit
#' @param matches TODOTODO
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#' @return numeric matrix with rows=samples and columns=TLS signature
#' @export
#'
#' @examples
#' # TODOTODO
compute_TLS <- function(matches, RNA.tpm){
  # Subset RNA.tpm
  sub_gene.tpm <- RNA.tpm[matches, ]

  # Calculation: geometric mean (so-called log-average) [TPM, 1 offset]
  geom_mean <- apply(sub_gene.tpm, 2, function(X) exp(mean(log2(X + 1))))

  message("TLS score computed")
  return(data.frame(TLS = geom_mean, check.names = FALSE))
}

#' Compute IFNy signature score
#'
#' \code{compute_IFNy} computes IFNy signature score as the arithmetic mean of genes included
#' in the IFN-γ signature (Ayers et al., JCI, 2017)
#'
#' @importFrom stats na.omit
#' @param matches TODOTODO
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=IFNy signature score
#'
#' @export
#'
#' @examples
#' # TODOTODO
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

#' Compute Expanded Immune signature
#'
#' \code{compute_Ayers_expIS} computes Expanded Immune signature score as the arithmetic mean of genes included
#' in the Expanded Immune signature (Ayers et al., JCI, 2017)
#'
#' @importFrom stats na.omit
#' @param matches TODOTODO
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=Expanded Immune signature score
#'
#' @export
#'
#' @examples
#' # TODOTODO
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

#' Compute Roh immune score
#'
#' \code{compute_roh_IS} computes Roh immune score as the geometric-mean of immune score genes
#' (Roh et al., 2017).
#'
#' @importFrom stats na.omit
#' @param matches TODOTODO
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=Roh immune score
#'
#' @export
#'
#' @examples
#' # TODOTODO
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

#' Compute Davoli immune signature
#'
#' \code{compute_davoli_IS} computes Davoli immune signature as the arithmetic mean of cytotoxic
#' immune infiltrate signature genes, after rank normalization (Davoli et al., 2017).
#'
#' @importFrom stats na.omit
#' @param matches TODOTODO
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=Davoli immune signature
#'
#' @export
#'
#' @examples
#' # TODOTODO
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

#' Compute chemokine score
#'
#' \code{compute_chemokine} computes chemoine score as the PC1 score that results from applying PCA
#' to z-score expression of 12 chemokine genes (Messina et al., 2012).
#'
#' @importFrom stats na.omit prcomp
#' @param matches TODOTODO
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#'
#' @return numeric matrix with rows=samples and columns=chemokine score
#'
#' @export
#'
#' @examples
#' # TODOTODO
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

#' Compute T cell-inflamed signature score
#'
#' \code{compute_ayersTcellInfl} computes T cell-inflamed signature score by taking a weighted sum of
#'  the housekeeping normalized values of the T cell-inflamed signature genes
#'
#' @importFrom stats na.omit
#'
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#' @param housekeeping TODOTODO
#' @param predictors TODOTODO
#' @param weights TODOTODO
#'
#' @return numeric matrix with rows=samples and columns=T cell-inflamed signature score
#'
#' @export
#'
#' @examples
#' # TODOTODO
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

#' Compute Immuno-Predictive Score (IMPRES)
#'
#' \code{compute_IMPRES} computes IMPRES score by applying logical comparison of checkpoint gene pairs
#' (Auslander et al., 2018).
#'
#' @importFrom stats na.omit
#'
#' @param RNA.tpm numeric matrix with rows=genes and columns=samples
#' @param sig TODOTODO
#' @param len TODOTODO
#' @param match_F_1 TODOTODO
#' @param match_F_2 TODOTODO
#'
#' @return numeric matrix with rows=samples and columns=IMPRES score
#'
#' @export
#'
#' @examples
#' # TODOTODO
compute_IMPRES_MSI <- function(sig, len, match_F_1, match_F_2, RNA.tpm){
  # Initialize variables
  F_pair_expr_A <- F_pair_expr_B <- IMPRES.matrix <- matrix(0, len, ncol(RNA.tpm))
  colnames(IMPRES.matrix) <- colnames(RNA.tpm)
  score <- vector("numeric", length = ncol(RNA.tpm))
  names(score) <- colnames(RNA.tpm)

  # Log2 transformation:
  log2.RNA.tpm <- as.data.frame(log2(RNA.tpm + 1))

  # Calculation:
  F_pair_expr_A <- log2.RNA.tpm[match_F_1, ]
  F_pair_expr_B <- log2.RNA.tpm[match_F_2, ]

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

  message(paste(sig,"score computed"))
  df <- data.frame(score, check.names = FALSE)
  names(df)[1] <- sig

  return(df)
}

#' Compute the expression of the immune checkpoints genes
#'
#' \code{computation_ICB_genes} computes the scores for the immune checkpoint genes.
#'
#' @export
#'
#' @param RNA.tpm numeric matrix with data
#' @param genes TODOTODO
#'
#' @return Data.frame with the expression of the immune checkpoint genes
#'
#' @examples
#' # TODOTODO
compute_ICB_genes <- function(genes, RNA.tpm){
  # Extract position genes for GZMA and PRF1
  tmp <- match(genes, rownames(RNA.tpm))

  df <- data.frame(
    PDL1=rep(NA, ncol(RNA.tpm)),
    CTLA4=rep(NA, ncol(RNA.tpm)),
    PD1=rep(NA, ncol(RNA.tpm)))

  # PDL-1 calculation
  if(!is.na(tmp[1])) df$PDL1 = RNA.tpm[tmp[1],]

  # CTLA-4 calculation
  if(!is.na(tmp[2])) df$CTLA4 = RNA.tpm[tmp[2],]

  # PD-1 calculation
  if(!is.na(tmp[3])) df$PD1 = RNA.tpm[tmp[3],]
  message("ICB genes PDL1, CTLA4, PD1 computed")
  return(df)
}

