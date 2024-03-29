% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_IMPRES_MSI.R
\name{compute_IMPRES_MSI}
\alias{compute_IMPRES_MSI}
\title{Compute Immuno-Predictive Score (IMPRES) and
Micro Satellite Instability (MSI) status score}
\usage{
compute_IMPRES_MSI(sig, len, match_F_1, match_F_2, RNA_tpm)
}
\arguments{
\item{sig}{can be either 'IMPRES' or 'MSI'.}

\item{len}{the length of gene_1 vector.}

\item{match_F_1}{numeric vector indicating the index of signature
genes defined in 'gene_1' in \code{RNA_tpm}.}

\item{match_F_2}{numeric vector indicating the index of signature
genes defined in 'gene_2' in \code{RNA_tpm}.}

\item{RNA_tpm}{data.frame containing TPM values with HGNC symbols
in rows and samples in columns.}
}
\value{
A numeric matrix with samples in rows and IMPRES score in
a column.
}
\description{
Calculates IMPRES score by logical comparison of
checkpoint gene pairs expression, as defined in
Auslander et al., Nat. Med., 2018.
}
\details{
Calculates MSI status score by logical comparison
of MSI-related gene pairs, as defined in Fu et al.,
BMC Genomics, 2019.
}
\references{
Auslander,N., Zhang,G., Lee,J.S., Frederick, D.T., Miao,
B., Moll,T.,Tian, T., Wei,Z., Madan, S., Sullivan, R.J.,
et al. (2018). Robust prediction of response to immune
checkpoint blockade therapy in metastatic melanoma. Nat.
Med. 24, 1545–1549. https://doi.org/10.1038/s41591-018-0157-9.

Fu, Y., Qi, L., Guo, W., Jin, L., Song, K., You, T.,
Zhang, S., Gu, Y., Zhao, W., and Guo, Z. (2019). A qualitative
transcriptional signature for predicting microsatellite
instability status of right-sided Colon Cancer. BMC Genomics
20, 769.
}
