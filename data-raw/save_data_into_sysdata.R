# setwd("/Volumes/nonshib-webdav/SystemsImmunoOncology/Mechanistic_signatures_project/")
#
# # Cancer types
# load("analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS.RData")
# PanCancer_names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS)
#
# Views
views <- c(Pathways.cor = 'gaussian',
           ImmuneCells = 'gaussian',
           TFs = 'gaussian',
           LRpairs.spec.pc = 'gaussian',
           CCpairsGroupedScores.spec.pc = 'gaussian')

single_views <- list(views[1], views[2], views[3], views[4], views[5])
combined_views <- lapply(1:10, function(X){
  tmp <- c("gaussian")
  name_tmp <- paste(combn(names(views), m = 2)[,X], collapse = "_")
  names(tmp) <- name_tmp
  return(tmp)
})
combined_views <- combined_views[-10]
#
# # general
# load("data/cor_genes_ICB_proxies.RData")
# load("data/IPS_genes.RData")
# load("data/list_top100pathways_responsive_genes.RData")
# load("data/Ligand_Receptors_Rdata/intercell.network.CC.pairs.grouped.cancer.spec.RData")
# lr.frequency <- LR.frequency
# load("data/resistance.program.RData")
# load("data/grouping_lrpairs_features_info.RData")
# load("data/TCGA_mean_sd.RData")
# #
cor_genes_to_remove
IPSG_read
top_100_per_pathway_responsive_genes
TCGA_mean_pancancer <- TCGA.mean.pancancer
TCGA_sd_pancancer <- TCGA.sd.pancancer
lr_frequency <- lr.frequency
intercell_network_cancer_spec <- intercell.network.cancer.spec
res_sig <- res.sig
grouping_lrpairs_info

#
views_combination <- c(single_views, combined_views)
# views_combination[[15]] <- c(Transcript = 'gaussian')
# Save learned models in a list
alg <- "Multi_Task_EN"
trained_models <- list()
trained_models <- lapply("SKCM", function(CancerType){

  print(CancerType)
  file <- dir(path = paste0("../../../../Desktop/SKCM/"),
                           # output/PanCancer_draft_v1/", CancerType, "/group_cor_normalized_tasks"),
              pattern = "all_cv_res_", full.names = T, recursive = F)

  trained_models <- lapply(1:length(views_combination), function(X){
    DataType <- names(views_combination[[X]])
    print(DataType)
    load(file[grep(pattern = paste0("_with_cor_tasks_", DataType,".RData"), file, fixed = T)])

    model_alg <- all_cv_res[[alg]]
    model <- lapply(1:length(model_alg), function(iter){
      model_alg[[iter]][["performances"]] <- NULL
      model_alg[[iter]][["training_set"]] <- NULL
      model_alg[[iter]][["model"]][["cv.glmnet.mse"]] <- NULL
      return(model_alg[[iter]])
    })

    return(model)
  })
  names(trained_models) <- sapply(1:length(views_combination), function(X) names(views_combination[[X]]))
  return(trained_models)
})
names(trained_models) <- "SKCM"
# save(trained_models, file = "data/learned_model_for_package.RData")
#
# load("data/learned_model_for_package.RData")
#
# setwd("~/ownCloud2/SystemsImmunoOncology/Mechanistic_signatures_project/easier/")
# usethis::use_data(cor_genes_to_remove,
#                   TCGA.mean.pancancer,
#                   TCGA.sd.pancancer,
#                   IPSG_read,
#                   res.sig,
#                   grouping_lrpairs_info,
#                   top_100_per_pathway_responsive_genes,
#                   intercell.network.cancer.spec,
#                   lr.frequency,
#                   trained_models,
#                   internal = TRUE, overwrite = TRUE, compress = "xz")
#
# # Riaz data
# # usethis::use_data(Riaz_data, compress = "xz")
#
#
# # usethis::use_data_raw(Riaz_data, compress = "xz")
#
#
#
# Change immune cell names (to shorter ones) in trained models

old_cellnames <- names(trained_models$SKCM$ImmuneCells[[1]]$mas.mea.learning.X[[1]])
new_cellnames <- c("B", "M1", "M2", "Monocyte", "Neutrophil", "NK", "CD4 T", "CD8+ T", "Treg", "DC", "Other")

# load trained models
load("~/ownCloud2/SystemsImmunoOncology/easier_project/easier_devel/R/sysdata.rda")

views <- names(trained_models$SKCM)
where <- grep("ImmuneCells", views, fixed = TRUE)

for(X in where){
  for(Y in 1:100){

    old_features <- rownames(trained_models[["SKCM"]][[X]][[Y]][["model"]][["cv.glmnet.features"]][["1se.mse"]])

    old_features[match(old_cellnames, old_features)] <- new_cellnames
    rownames(trained_models[["SKCM"]][[X]][[Y]][["model"]][["cv.glmnet.features"]][["1se.mse"]]) <- old_features
    rownames(trained_models[["SKCM"]][[X]][[Y]][["model"]][["cv.glmnet.features"]][["min.mse"]]) <- old_features

    if (sapply(strsplit(views[X], split = "_"), head, 1) == "ImmuneCells") {
      names(trained_models[["SKCM"]][[X]][[Y]][["mas.mea.learning.X"]][[1]]) <- new_cellnames
    } else if (sapply(strsplit(views[X], split = "_"), tail, 1) == "ImmuneCells") {
      names(trained_models[["SKCM"]][[X]][[Y]][["mas.mea.learning.X"]][[2]]) <- new_cellnames
    }
  }
}


setwd("~/ownCloud2/SystemsImmunoOncology/easier_project/easier_devel/")
usethis::use_data(cor_genes_to_remove,
                  TCGA_mean_pancancer,
                  TCGA_sd_pancancer,
                  IPSG_read,
                  res_sig,
                  grouping_lrpairs_info,
                  top_100_per_pathway_responsive_genes,
                  intercell_network_cancer_spec,
                  lr_frequency,
                  trained_models,
                  internal = TRUE, overwrite = TRUE, compress = "xz")


