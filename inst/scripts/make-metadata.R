## validate with `ExperimentHubData::makeExperimentHubMetadata()`
## (do this above pkg directory)

Mariathasan2018_PDL1_treatment <- data.frame(
  stringsAsFactors = FALSE,
  Title = "Mariathasan2018_PDL1_treatment",
  Description = paste(
    "Processed RNA-seq data, quantified as counts and TPM,",
    "for bladder cancer samples pre anti-PDL1 immunotherapy.",
    "Information on patient's tumor mutational burden",
    "and response to anti-PD1 treatment is also provided.",
    "Represented as a SummarizedExperiment;",
    "derived from IMvigor210CoreBiologies R package",
    "(Mariathasan et al., Nature, 2018) under the CC-BY license."),
  BiocVersion = "3.13",
  Genome = NA,
  SourceType = "tar.gz",
  SourceUrl = "http://research-pub.gene.com/IMvigor210CoreBiologies/",
  SourceVersion = "Feb 14 2018",
  Species = "Homo sapiens",
  TaxonomyId = "9606",
  Coordinate_1_based  = NA,
  DataProvider = "Mariathasan S, Turley S, Nickles D et al., “TGF-b attenuates tumor response to PD-L1 blockade by contributing to exclusion of T cells.",
  Maintainer = "Federico Marini <marinif@uni-mainz.de>",
  RDataClass = "SummarizedExperiment",
  DispatchClass = "Rda",
  RDataPath = "easier/Mariathasan2018_PDL1_treatment.Rda")

# write to .csv
write.csv(Mariathasan2018_PDL1_treatment, file = "extdata/metadata.csv", row.names = FALSE)
