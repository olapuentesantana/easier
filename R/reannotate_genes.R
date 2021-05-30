#' Gene re-annotation using HGNC symbols
#'
#' This function performs gene re-annotation using curated data by the HGNC.
#'
#' @param cur_genes A character string containing gene HGNC symbols to be consider for re-annotation.
#'
#' @return A data.frame with the old gene HGNC symbol and the new corresponding gene HGNC symbol.
#'
#' @examples
#' # use example dataset from Mariathasan cohort (Mariathasan et al., Nature, 2018)
#' data(cds)
#' mariathasan_data <- preprocess_mariathasan(cds)
#' gene_count <- mariathasan_data$counts
#' rm(cds)
#'
#' genes_info <- map_genes(cur_genes = rownames(gene_count))
reannotate_genes <- function(cur_genes){

  # cur_genes <- rownames(gene_expr)
  new_genes <- rep(NA,length(cur_genes))
  new_genes2 <- rep(NA,length(cur_genes))
  ind <- match(cur_genes, HGNC$ApprovedSymbol)

  # Current symbols and withdrawn ones
  genes_ind_notNA <- which(!is.na(ind))
  for (i in genes_ind_notNA) {
    genei <- cur_genes[i]
    if (HGNC$Status[ind[i]]=="Approved") {
      new_genes[i] <- cur_genes[i]

    } else if (HGNC$Status[ind[i]]=="EntryWithdrawn") {
      next

    } else {
      W_string <- "symbolwithdrawn,see"
      new_symbol <- gsub(W_string, "", HGNC$ApprovedName[ind[i]])
      new_genes2[i] <- new_symbol
    }
  }

  # Not found as symbols
  genes_ind_NA <- which(is.na(ind))
  for (i in genes_ind_NA) {
    genei <- cur_genes[i]

    # Previos symbol?
    ind1 <- grep(genei,HGNC$PreviousSymbols)
    for (i1 in ind1) {
      array1 <- unlist(strsplit(as.character(HGNC$PreviousSymbols[i1]), "|", fixed = TRUE))
      flag1 <- length(which(array1==genei))>0
      if (flag1) {
        new_symbol <- as.character(HGNC$ApprovedSymbol[i1])
        new_genes2[i] <- new_symbol
      }
    }
    # Synonym?
    ind2 <- grep(genei,HGNC$Synonyms)
    for (i2 in ind2) {
      array2 <- unlist(strsplit(as.character(HGNC$Synonyms[i2]), "|"))
      flag2 <- length(which(array2==genei))>0
      if (flag2) {
        new_symbol <- as.character(HGNC$ApprovedSymbol[i2])
        new_genes2[i] <- new_symbol
      }
    }

  }

  new_genes2[which(new_genes2 %in% setdiff(new_genes, NA))] <- NA
  ind <- intersect(which(is.na(new_genes)),
                   which(!is.na(new_genes2)))
  new_genes[ind] <- new_genes2[ind]

  out <- data.frame(old_names = cur_genes, new_names = new_genes)

  return(out)

}
