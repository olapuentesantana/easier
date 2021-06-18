# setwd("/Volumes/nonshib-webdav/SystemsImmunoOncology/Mechanistic_signatures_project/")
#
# # Cancer types
# load("analysis/pre-processing/TCGA_samples_available_screening_with_quanTIseq_IS.RData")
# PanCancer_names <- names(TCGA.samples.pancancer_with_screen_quantiseg_IS)
#

# Old names for Views
views <- c(
  Pathways.cor = "gaussian",
  ImmuneCells = "gaussian",
  TFs = "gaussian",
  LRpairs.spec.pc = "gaussian",
  CCpairsGroupedScores.spec.pc = "gaussian"
)

single_views <- list(views[1], views[2], views[3], views[4], views[5])
combined_views <- lapply(1:10, function(X) {
  tmp <- c("gaussian")
  name_tmp <- paste(combn(names(views), m = 2)[, X], collapse = "_")
  names(tmp) <- name_tmp
  return(tmp)
})
combined_views <- combined_views[-10]

# New names for Views
views_new <- c(
  pathways = "gaussian",
  immunecells = "gaussian",
  tfs = "gaussian",
  lrpairs = "gaussian",
  ccpairs = "gaussian"
)

single_views_new <- list(views_new[1], views_new[2], views_new[3], views_new[4], views_new[5])
combined_views_new <- lapply(1:10, function(X) {
  tmp <- c("gaussian")
  name_tmp <- paste(combn(names(views_new), m = 2)[, X], collapse = "_")
  names(tmp) <- name_tmp
  return(tmp)
})
combined_views_new <- combined_views_new[-10]

# --------------------------------------------------------------- #
# load data needed for easier
setwd("~/ownCloud2/SystemsImmunoOncology/Mechanistic_signatures_project/")

load("data/cor_genes_ICB_proxies.RData")
load("data/Ligand_Receptors_Rdata/intercell.network.CC.pairs.grouped.cancer.spec.RData")
lr_frequency <- LR.frequency
load("data/resistance.program.RData")
load("data/grouping_lrpairs_features_info.RData")
load("data/TCGA_mean_sd.RData")
# #
cor_genes_to_remove
TCGA_mean_pancancer <- TCGA.mean.pancancer
TCGA_sd_pancancer <- TCGA.sd.pancancer
intercell_network_cancer_spec <- intercell.network.cancer.spec
res_sig <- res.sig
grouping_lrpairs_info

#
views_combination <- c(single_views, combined_views)
views_combination_new <- c(single_views_new, combined_views_new)

# views_combination[[15]] <- c(Transcript = 'gaussian')

# --------------------------------------------------------------- #
# Save learned models in a list
alg <- c("Multi_Task_EN")
trained_models <- list()
trained_models <- lapply(c("BLCA"), function(CancerType) {
  print(CancerType)
  file <- dir(
    path = file.path("../../../Desktop/STAD"),
    pattern = "all_cv_res_", full.names = T, recursive = F
  )

  trained_models <- lapply(1:length(views_combination), function(X) {
    DataType <- names(views_combination[[X]])
    print(DataType)
    load(file[grep(pattern = paste0("_with_cor_tasks_", DataType, ".RData"), file, fixed = T)])

    model_alg <- all_cv_res[[alg]]
    model <- lapply(1:length(model_alg), function(iter) {
      model_alg[[iter]][["performances"]] <- NULL
      model_alg[[iter]][["training_set"]] <- NULL
      model_alg[[iter]][["model"]][["cv.glmnet.mse"]] <- NULL
      return(model_alg[[iter]])
    })

    return(model)
  })
  names(trained_models) <- sapply(1:length(views_combination_new), function(X) names(views_combination_new[[X]]))
  return(trained_models)
})

names(trained_models) <- "BLCA"

# --------------------------------------------------------------- #
# Change immune cell names (to shorter ones) in trained models
new_cellnames <- c("B", "M1", "M2", "Monocyte", "Neutrophil", "NK", "CD4 T", "CD8+ T", "Treg", "DC", "Other")

CancerType <- "STAD"
old_cellnames <- names(trained_models[[CancerType]]$immunecells_lrpairs[[1]]$mas.mea.learning.X[[1]])
views <- names(trained_models[[CancerType]])
where <- grep("immunecells", views, fixed = TRUE)

for (X in setdiff(where, 2)) {
  for (Y in 1:100) {
    old_features <- rownames(trained_models[[CancerType]][[X]][[Y]][["model"]][["cv.glmnet.features"]][["1se.mse"]])

    old_features[match(old_cellnames, old_features)] <- new_cellnames
    rownames(trained_models[[CancerType]][[X]][[Y]][["model"]][["cv.glmnet.features"]][["1se.mse"]]) <- old_features
    rownames(trained_models[[CancerType]][[X]][[Y]][["model"]][["cv.glmnet.features"]][["min.mse"]]) <- old_features

    if (sapply(strsplit(views[X], split = "_"), head, 1) == "immunecells") {
      names(trained_models[[CancerType]][[X]][[Y]][["mas.mea.learning.X"]][[1]]) <- new_cellnames
    } else if (sapply(strsplit(views[X], split = "_"), tail, 1) == "immunecells") {
      names(trained_models[[CancerType]][[X]][[Y]][["mas.mea.learning.X"]][[2]]) <- new_cellnames
    }
  }
}

# --------------------------------------------------------------- #
# Restructure sysdata

tasks_new_names <- c("CYT", "Ock_IS", "Roh_IS", "chemokines", "Davoli_IS", "IFNy", "Ayers_expIS", "Tcell_inflamed", "RIR", "TLS")
task_old_names <- c("CYT", "IS", "RohIS", "chemokine", "IS_Davoli", "IFny", "ExpandedImmune", "T_cell_inflamed", "resF.down", "TLS")

# collect optimization models
opt_models <- lapply(names(trained_models), function(cancertype) {
  opt_models <- lapply(names(trained_models[[cancertype]]), function(view) {
    opt_models <- lapply(colnames(trained_models[[cancertype]][[view]][[1]][["model"]][["cv.glmnet.features"]][["1se.mse"]]), function(task) {
      opt_models <- do.call(cbind, lapply(1:100, function(iter) {
        tmp <- trained_models[[cancertype]][[view]][[iter]][["model"]][["cv.glmnet.features"]][["1se.mse"]][, task]

        return(tmp)
      }))
      return(opt_models)
    })
    names(opt_models) <- tasks_new_names
    return(opt_models)
  })
  names(opt_models) <- c("pathways", "immunecells", "tfs", "lrpairs", "ccpairs")
  return(opt_models)
})
names(opt_models) <- c("LUAD", "LUSC", "BLCA", "BRCA", "CESC", "CRC", "GBM", "HNSC", "KIRC", "KIRP", "LIHC", "OV", "PAAD", "PRAD", "SKCM", "STAD", "THCA", "UCEC")

# collect training statistics (mean and sd)
opt_xtrain_stats <- lapply(names(trained_models), function(cancertype) {
  opt_xtrain_stats <- lapply(names(trained_models[[cancertype]]), function(view) {
    opt_xtrain_stats <- lapply(c("mean", "sd"), function(stat) {
      opt_xtrain_stats <- do.call(cbind, lapply(1:100, function(iter) {
        if (stat %in% "mean") {
          tmp <- trained_models[[cancertype]][[view]][[iter]]$mas.mea.learning.X[[1]]
        } else {
          tmp <- trained_models[[cancertype]][[view]][[iter]]$mas.std.learning.X[[1]]
          names(tmp) <- names(trained_models[[cancertype]][[view]][[iter]]$mas.mea.learning.X[[1]])
        }
        return(tmp)
      }))
      return(opt_xtrain_stats)
    })
    names(opt_xtrain_stats) <- c("mean", "sd")
    return(opt_xtrain_stats)
  })
  names(opt_xtrain_stats) <- c("pathways", "immunecells", "tfs", "lrpairs", "ccpairs")
  return(opt_xtrain_stats)
})
names(opt_xtrain_stats) <- c("LUAD", "LUSC", "BLCA", "BRCA", "CESC", "CRC", "GBM", "HNSC", "KIRC", "KIRP", "LIHC", "OV", "PAAD", "PRAD", "SKCM", "STAD", "THCA", "UCEC")


HGNC <- read.csv(file.path("../../Anti_PD1_challenge/Anti-PD1-DREAM-cSysImmunoOnco_system/tmp_data/HGNC_genenames_20170418.txt"),
  header = TRUE, sep = "\t"
)


setwd("~/ownCloud2/SystemsImmunoOncology/easier_project/easier_devel/")
usethis::use_data(cor_genes_to_remove,
  TCGA_mean_pancancer,
  TCGA_sd_pancancer,
  res_sig,
  grouping_lrpairs_info,
  intercell_network_cancer_spec,
  lr_frequency,
  opt_models,
  opt_xtrain_stats,
  HGNC,
  internal = TRUE, overwrite = TRUE, compress = "xz"
)

###################
# add NSCLC model

# Save learned models in a list
alg <- c("Multi_Task_EN")
trained_models <- list()
trained_models <- lapply(c("NSCLC"), function(CancerType) {
  print(CancerType)
  file <- dir(
    path = file.path("../Anti_PD1_challenge/LUAD_LUSC"),
    pattern = "all_cv_res_", full.names = T, recursive = F
  )

  trained_models <- lapply(1:length(single_views), function(X) {
    DataType <- names(single_views[[X]])
    print(DataType)
    load(file[grep(pattern = paste0("_with_cor_tasks_", DataType, ".RData"), file, fixed = T)])

    model_alg <- all_cv_res[[alg]]
    model <- lapply(1:length(model_alg), function(iter) {
      model_alg[[iter]][["performances"]] <- NULL
      model_alg[[iter]][["training_set"]] <- NULL
      model_alg[[iter]][["model"]][["cv.glmnet.mse"]] <- NULL
      return(model_alg[[iter]])
    })

    return(model)
  })
  names(trained_models) <- sapply(1:length(single_views_new), function(X) names(single_views_new[[X]]))
  return(trained_models)
})

names(trained_models) <- "NSCLC"

# --------------------------------------------------------------- #
# Change immune cell names (to shorter ones) in trained models
new_cellnames <- c("B", "M1", "M2", "Monocyte", "Neutrophil", "NK", "CD4 T", "CD8+ T", "Treg", "DC", "Other")

CancerType <- "NSCLC"
old_cellnames <- names(trained_models[[CancerType]]$immunecells[[1]]$mas.mea.learning.X[[1]])
views <- names(trained_models[[CancerType]])
where <- grep("immunecells", views, fixed = TRUE)

for (X in 2) {
  for (Y in 1:100) {
    old_features <- rownames(trained_models[[CancerType]][[X]][[Y]][["model"]][["cv.glmnet.features"]][["1se.mse"]])

    old_features[match(old_cellnames, old_features)] <- new_cellnames
    rownames(trained_models[[CancerType]][[X]][[Y]][["model"]][["cv.glmnet.features"]][["1se.mse"]]) <- old_features
    rownames(trained_models[[CancerType]][[X]][[Y]][["model"]][["cv.glmnet.features"]][["min.mse"]]) <- old_features

    if (sapply(strsplit(views[X], split = "_"), head, 1) == "immunecells") {
      names(trained_models[[CancerType]][[X]][[Y]][["mas.mea.learning.X"]][[1]]) <- new_cellnames
    } else if (sapply(strsplit(views[X], split = "_"), tail, 1) == "immunecells") {
      names(trained_models[[CancerType]][[X]][[Y]][["mas.mea.learning.X"]][[2]]) <- new_cellnames
    }
  }
}

# --------------------------------------------------------------- #
# Restructure sysdata

tasks_new_names <- c("CYT", "Ock_IS", "Roh_IS", "chemokines", "Davoli_IS", "IFNy", "Ayers_expIS", "Tcell_inflamed", "RIR", "TLS")
task_old_names <- c("CYT", "IS", "RohIS", "chemokine", "IS_Davoli", "IFny", "ExpandedImmune", "T_cell_inflamed", "resF.down", "TLS")

# collect optimization models
opt_model_NSCLC <- lapply(names(trained_models), function(cancertype) {
  opt_model_NSCLC <- lapply(names(trained_models[[cancertype]]), function(view) {
    opt_model_NSCLC <- lapply(colnames(trained_models[[cancertype]][[view]][[1]][["model"]][["cv.glmnet.features"]][["1se.mse"]]), function(task) {
      opt_model_NSCLC <- do.call(cbind, lapply(1:100, function(iter) {
        tmp <- trained_models[[cancertype]][[view]][[iter]][["model"]][["cv.glmnet.features"]][["1se.mse"]][, task]
        return(tmp)
      }))
      return(opt_model_NSCLC)
    })
    names(opt_model_NSCLC) <- tasks_new_names
    return(opt_model_NSCLC)
  })
  names(opt_model_NSCLC) <- c("pathways", "immunecells", "tfs", "lrpairs", "ccpairs")
  return(opt_model_NSCLC)
})
names(opt_model_NSCLC) <- c("NSCLC")

# collect training statistics (mean and sd)
opt_xtrain_stats_NSCLC <- lapply(names(trained_models), function(cancertype) {
  opt_xtrain_stats_NSCLC <- lapply(names(trained_models[[cancertype]]), function(view) {
    opt_xtrain_stats_NSCLC <- lapply(c("mean", "sd"), function(stat) {
      opt_xtrain_stats_NSCLC <- do.call(cbind, lapply(1:100, function(iter) {
        if (stat %in% "mean") {
          tmp <- trained_models[[cancertype]][[view]][[iter]]$mas.mea.learning.X[[1]]
        } else {
          tmp <- trained_models[[cancertype]][[view]][[iter]]$mas.std.learning.X[[1]]
          names(tmp) <- names(trained_models[[cancertype]][[view]][[iter]]$mas.mea.learning.X[[1]])
        }
        return(tmp)
      }))
      return(opt_xtrain_stats_NSCLC)
    })
    names(opt_xtrain_stats_NSCLC) <- c("mean", "sd")
    return(opt_xtrain_stats_NSCLC)
  })
  names(opt_xtrain_stats_NSCLC) <- c("pathways", "immunecells", "tfs", "lrpairs", "ccpairs")
  return(opt_xtrain_stats_NSCLC)
})
names(opt_xtrain_stats_NSCLC) <- c("NSCLC")

opt_models$NSCLC <- opt_model_NSCLC$NSCLC
opt_xtrain_stats$NSCLC <- opt_xtrain_stats_NSCLC$NSCLC

setwd("~/ownCloud2/SystemsImmunoOncology/easier_project/easier_devel/")
usethis::use_data(cor_genes_to_remove,
                  TCGA_mean_pancancer,
                  TCGA_sd_pancancer,
                  res_sig,
                  grouping_lrpairs_info,
                  intercell_network_cancer_spec,
                  lr_frequency,
                  opt_models,
                  opt_xtrain_stats,
                  HGNC,
                  internal = TRUE, overwrite = TRUE, compress = "xz"
)



