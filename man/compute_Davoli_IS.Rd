% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_Davoli_IS.R
\name{compute_Davoli_IS}
\alias{compute_Davoli_IS}
\title{Compute Davoli immune signature (Davoli_IS) score}
\usage{
compute_Davoli_IS(matches, RNA_tpm)
}
\arguments{
\item{matches}{numeric vector indicating the index of signature
genes in \code{RNA_tpm}.}

\item{RNA_tpm}{data.frame containing TPM values with HGNC symbols
in rows and samples in columns.}
}
\value{
A numeric matrix with samples in rows and Davoli_IS score
in a column.
}
\description{
Calculates Davoli_IS score as the average of the expression
of its signature genes after applying rank normalization, as
defined in Davoli et al., Science, 2017.
}
\references{
Davoli, T., Uno, H., Wooten, E.C., and Elledge,
S.J. (2017). Tumor aneuploidy correlates with markers of
immune evasion and with reduced response to immunotherapy.
Science 355. https://doi.org/10.1126/science.aaf8399.
}
