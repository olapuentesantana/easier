# easier run only for cancer types with an associated cancer-specific model...
library(easier)
library(easierData)
library(SummarizedExperiment)

# testing
test_that("easier runs on cancer type with available cancer-specific model",{

  # Using example exemplary dataset (Mariathasan et al., Nature, 2018)
  # from easierData. Original processed data is available from
  # IMvigor210CoreBiologies package.
  dataset_mariathasan <- easierData::get_Mariathasan2018_PDL1_treatment()
  RNA_tpm <- SummarizedExperiment::assays(dataset_mariathasan)[["tpm"]]

  # Select a subset of patients to reduce vignette building time.
  pat_subset <- c(
      "SAM76a431ba6ce1", "SAMd3bd67996035", "SAMd3601288319e",
      "SAMba1a34b5a060", "SAM18a4dabbc557"
  )
  RNA_tpm <- RNA_tpm[, colnames(RNA_tpm) %in% pat_subset]

  # Computation of TF activity (Garcia-Alonso et al., Genome Res, 2019)
  tf_activity <- easier::compute_TF_activity(
      RNA_tpm = RNA_tpm
  )

  # provide right cancer type input
  predictions <- easier::predict_immune_response(
    tfs = tf_activity,
    cancer_type = "BLCA"
  )
  expect_type(predictions, "list")

})

test_that("easier fails to run on cancer type with non-available cancer-specific model",{

  # Using example exemplary dataset (Mariathasan et al., Nature, 2018)
  # from easierData. Original processed data is available from
  # IMvigor210CoreBiologies package.
  dataset_mariathasan <- easierData::get_Mariathasan2018_PDL1_treatment()
  RNA_tpm <- SummarizedExperiment::assays(dataset_mariathasan)[["tpm"]]

  # Select a subset of patients to reduce vignette building time.
  pat_subset <- c(
    "SAM76a431ba6ce1", "SAMd3bd67996035", "SAMd3601288319e",
    "SAMba1a34b5a060", "SAM18a4dabbc557"
  )
  RNA_tpm <- RNA_tpm[, colnames(RNA_tpm) %in% pat_subset]

  # Computation of TF activity (Garcia-Alonso et al., Genome Res, 2019)
  tf_activity <- easier::compute_TF_activity(
    RNA_tpm = RNA_tpm
  )

  # provide wrong cancer type input
  expect_error(
    predictions <- easier::predict_immune_response(
      tfs = tf_activity,
      cancer_type = "LAML"
    )
  )

})
