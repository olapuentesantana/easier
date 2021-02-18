#' Compute transription factors activity
#'
#' `compute_TF_activity` infers transription factor activity from tpm RNA-seq data.
#'
#' @importFrom dorothea run_viper
#' @importFrom stats na.exclude
#' @importFrom dplyr filter
#'
#' @param RNA_tpm numeric matrix of tpm values with rows=genes and columns=samples
#' @param remove_genes_ICB_proxies boolean variable to remove all those genes
#' involved in the computation of ICB proxy's of response
#'
#' @return A list with the following elements:
#'   \describe{
#'     \item{scores}{TF activity matrix with rows=samples and columns=TFs}
#'     \item{transcripts_kept}{vector with available gene names}
#'     \item{transcripts_left}{vector with missing gene names}
#'   }
#' @export
#'
#' @examples
#' # TODOTODO
compute_TF_activity <- function(RNA_tpm,
                                remove_genes_ICB_proxies = FALSE){

  # Gene expression data
  tpm <- RNA_tpm
  genes <- rownames(tpm)

  # HGNC symbols are required
  try(if (any(grep("ENSG00000", genes))) stop("hgnc gene symbols are required", call. = FALSE))

  # Genes to remove according to all ICB proxy's
  if (remove_genes_ICB_proxies) {
    message("Removing signatures genes for proxy's of ICB response  \n")
    idy <- stats::na.exclude(match(cor_genes_to_remove, rownames(tpm)))
    tpm <- tpm[-idy,]
  }

  # Log transformed expression matrix (log2[tpm+1]): expression matrix scaled and recentered.
  gene_expr <- standardization(t(tpm) , mean = TCGA.mean.pancancer, sd = TCGA.sd.pancancer)

  # redefine gene names to match transcripts for viper
  E <- t(gene_expr)
  newNames <- sapply(rownames(E), function(x){
    # strsplit(x, "\\.")[[1]][1]
    zz_tmp <- strsplit(x, "\\.")[[1]]
    paste0(zz_tmp[1:(length(zz_tmp)-1)], collapse = "-")
  })
  rownames(E) <- newNames

  # data extracted from publication
  regulons <- dplyr::filter(dorothea::dorothea_hs, confidence %in% c("A", "B"))
  all_regulated_transcripts <- unique(regulons$target)
  all_tfs <- unique(regulons$tf)

  # check what is the percentage of regulated transcripts and TF that we have in our data
  message(" percentage of regulated transcripts = ", sum(all_regulated_transcripts %in%  rownames(E))*100/length(all_regulated_transcripts), "\n")
  message(" percentage of TF = ", sum(all_tfs %in% rownames(E))*100/length(all_tfs), "\n")

  # TF activity: run viper
  TF_activities <- dorothea::run_viper(input = E, regulons = regulons,
                                       options = list(method = "none", minsize = 4, eset.filter = FALSE, cores = 1, verbose=FALSE))

  # Samples as rows, TFs as columns
  TF_activities <- t(TF_activities)

  # check what is the percentage of genes we have in our data
  genes_kept <- intersect(rownames(E), all_regulated_transcripts)
  genes_left <- setdiff(all_regulated_transcripts, rownames(E))

  # Output list:
  TFs <- list(scores = as.data.frame(TF_activities),
              transcripts_kept = length(genes_kept),
              transcripts_left = length(genes_left))

  message("TF activities computed \n")

  return(TFs)
}
