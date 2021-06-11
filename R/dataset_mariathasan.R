#' Processed data from IMvigor210CoreBiologies R package (Mariathasan et al., Nature, 2018)
#'
#' An exemplary dataset with samples from 192 patients with bladder cancer
#'
#' @references Mariathasan, S., Turley, S., Nickles, D. et al. TGFβ attenuates tumour response
#' to PD-L1 blockade by contributing to exclusion of T cells. Nature 554, 544–548 (2018).
#' https://doi.org/10.1038/nature25501.
#'
#' @details easier ships with an example dataset with samples from 192 patients with bladder cncer.
#' The dataset easier::dataset_mariathasan contains
#' \describe{
#'   \item{name}{cohort name related to the original work from Mariathasan et al., Nature, 2018}
#'   \item{counts}{gene expression matrix, counts derived from bulk RNA-seq}
#'   \item{tpm}{gene expression matrix, tpm values derived from bulk RNA-seq}
#'   \item{response}{patient's response to immunotherapy, defining complete response (CR) patients as responders and progressive disease (PD) patients as non-responders}
#'   \item{TMB}{Tumor mutational burden (TMB) defined by using panel sequencing to estimate the TMB by including synonymous mutations additionally to non-synonymous mutations}
#'   \item{cancertype}{cancer type from which the bulk RNA-seq was obtained}
#' }
#' @source The processed data is made available via IMvigor210CoreBiologies R package
#' (\url{http://research-pub.gene.com/IMvigor210CoreBiologies/}) under the CC-BY license.
"dataset_mariathasan"
