

## defining signatures here
signature_davoli <- c("CD247", "CD2", "CD3E", "GZMH", "NKG7", "PRF1", "GZMK")

## alternate thought: maybe make a list-like dataset with all signatures?
# usage: easier_sigs$signature_davoli, and that would be returning the content of the current Davoli_IS.read
easier_sigs <- list()
easier_sigs$signature_davoli <- signature_davoli


# if making a dataset, it would also be "exported". Otherwise, it will be "accessed only as internal", i.e. easier:::OBJECT
# Reason: easier (not the package :) )  to maintain, or to extend and refine!
# taking it to the next level: it could be a generic compute_signature function



compute_signature <- function(RNA_tpm,
                              selected_signature) {

  # will do some checks on these objects
  ## RNA_tpm is a matrix-like/df
  ## selected signature must be a vector
  ### -- ideally one could use the generic function right away, therefore the extra checks

  # Log2 transformation:
  log2.RNA.tpm <- log2(RNA_tpm + 1)

  # Subset log2.RNA.tpm
  sub_log2.RNA.tpm <- log2.RNA.tpm[selected_signature, ] # here: it was the specific signature

  # Calculate rank position for each gene across samples
  ranks_sub_log2.RNA.tpm <- apply(sub_log2.RNA.tpm, 1, rank)

  # Get normalized rank by divided
  ranks_sub_log2.RNA.tpm.norm <- (ranks_sub_log2.RNA.tpm - 1)/(nrow(ranks_sub_log2.RNA.tpm) - 1)

  # Calculation: average of the expression value of all the genes within-sample
  score <- apply(ranks_sub_log2.RNA.tpm.norm, 1, mean)

  #### message("Davoli_IS score computed")
  #### return(data.frame(Davoli_IS = score, check.names = FALSE))

  # here still returned as a vector
  ## if needed as a df, could also be handled here (could need an extra param to name the column properly! - here, as well as in the specific signature computing function)
  return(score)

}


compute_Davoli_IS_pimped <- function(RNA.tpm,
                                     sig_davoli = NULL) {

  #### Literature genes
  #### Davoli_IS.read <- c("CD247", "CD2", "CD3E", "GZMH", "NKG7", "PRF1", "GZMK")

  if (is.null(sig_davoli))
    sig_davoli <- easier_sigs$signature_davoli # or just the single one
  # else: keep the selection provided by the user


  match_Davoli_IS.genes <- match(sig_davoli, rownames(RNA.tpm))

  if (anyNA(match_Davoli_IS.genes)){
    warning(c("differenty named or missing signature genes : \n", paste(sig_davoli[!sig_davoli %in% rownames(RNA.tpm)], collapse = "\n")))
    match_Davoli_IS.genes <- stats::na.omit(match_Davoli_IS.genes)
  }

  # then using the generic function
  score <- compute_signature(RNA_tpm = RNA.tpm,
                             selected_signature = match_Davoli_IS.genes)


  # and then again returning the score, just as for the original version


  #### # Log2 transformation:
  #### log2.RNA.tpm <- log2(RNA.tpm + 1)
  ####
  #### # Subset log2.RNA.tpm
  #### sub_log2.RNA.tpm <- log2.RNA.tpm[match_Davoli_IS.genes, ]
  ####
  #### # Calculate rank position for each gene across samples
  #### ranks_sub_log2.RNA.tpm <- apply(sub_log2.RNA.tpm, 1, rank)
  ####
  #### # Get normalized rank by divided
  #### ranks_sub_log2.RNA.tpm.norm <- (ranks_sub_log2.RNA.tpm - 1)/(nrow(ranks_sub_log2.RNA.tpm) - 1)
  ####
  #### # Calculation: average of the expression value of all the genes within-sample
  #### score <- apply(ranks_sub_log2.RNA.tpm.norm, 1, mean)

  message("Davoli_IS score computed")
  return(data.frame(Davoli_IS = score, check.names = FALSE))
}
