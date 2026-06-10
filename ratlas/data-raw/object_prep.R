# Initial object prep prior to converting to assay5 and on-disk
# This reduces the dataset to the required needs to for the app
# **This prep is done only once per dataset**

library(Seurat)
dir.create("./lean_datasets")

##### Dataset clean-up #####

### VTA_pain_dataset ###

# load original object from the analysis

VTA_pain_dataset <- readRDS(file = "./analysis_objects/VTA_Pain_final.rds")

## Modify a few metadata variables for consistency and add new ones
VTA_pain_dataset$Sex <- VTA_pain_dataset$sex
VTA_pain_dataset$Veh_CFA <- VTA_pain_dataset$treatment_A
VTA_pain_dataset$Sal_Mor <- VTA_pain_dataset$treatment_B
VTA_pain_dataset$Veh_CFA_Sal_Mor <- VTA_pain_dataset$treatment

VTA_pain_dataset$Veh_CFA_Sex <- paste(VTA_pain_dataset$Veh_CFA, VTA_pain_dataset$Sex, sep = "_")
VTA_pain_dataset$Sal_Mor_Sex <- paste(VTA_pain_dataset$Sal_Mor, VTA_pain_dataset$Sex, sep = "_")
VTA_pain_dataset$Veh_CFA_Sal_Mor_Sex <- paste(VTA_pain_dataset$Veh_CFA_Sal_Mor, VTA_pain_dataset$Sex, sep = "_")

# set order level for Sex:
VTA_pain_dataset$Sex <- factor(VTA_pain_dataset$Sex, 
                            levels=c("female", "male"))

## remove un-needed variables:
VTA_pain_dataset$sex <- NULL
VTA_pain_dataset$treatment_A <- NULL
VTA_pain_dataset$treatment_B <- NULL
VTA_pain_dataset$treatment <- NULL
VTA_pain_dataset$cellType_original <- NULL
VTA_pain_dataset$cellType_reversed <- NULL

### NAc_TMP ###
# NOTE: object already V5, no need to convert beyond this point, just slim down the obj

#TODO: change var name to final variable name after name decision
NAc_TMP <- readRDS(file = "./analysis_objects/5_NAc0069_integrated_clean_NAcOnly.rds")

## Modify a few metadata variables for consistency and add new ones
NAc_TMP$Stim <- NAc_TMP$treatment
NAc_TMP$CellType <- NAc_TMP$cellType
NAc_TMP$Sex <- NAc_TMP$sex
NAc_TMP$Stim_Sex <- paste(NAc_TMP$Stim, NAc_TMP$Sex, sep = "_")

# set order level for Sex:
NAc_TMP$Sex <- factor(NAc_TMP$Sex, 
                      levels=c("female", "male"))

## remove un-needed variables:
NAc_TMP$sex <- NULL
NAc_TMP$treatment <- NULL
NAc_TMP$cellType <- NULL
NAc_TMP$scDblFinder_doublets <- NULL

## rename reduction to match others:

NAc_TMP@reductions[["umap"]] <- NAc_TMP@reductions[["rnaharmony.umap"]]
NAc_TMP@reductions[["rnaharmony.umap"]] <- NULL

# Set a new key for updated name
Key(NAc_TMP@reductions[["umap"]]) <- "UMAP_"

# Make the embedding column names match the new key
colnames(NAc_TMP@reductions[["umap"]]@cell.embeddings) <- paste0("UMAP_", 1:ncol(NAc_TMP@reductions[["umap"]]@cell.embeddings))

##### DietSeurat #####
# Retain the required data only

### VTA_pain_dataset ###
VTA_pain_dataset <- DietSeurat(VTA_pain_dataset,
                   layers = c("data"),
                   assays = c("RNA"), #harmony, thus no integrated assay
                   dimreducs = c("pca","umap","harmony_rna"))

# save as new RDS objects
saveRDS(VTA_pain_dataset, file = "./lean_datasets/VTA_pain_dataset.rds", compress = FALSE)

# fetch cluster names and save for UI
## VTA_pain: current Idents are `cellType`
cluster_names_VTA_pain <- sort(unique(as.character(Idents(VTA_pain_dataset))))
saveRDS(cluster_names_VTA_pain, file = "./lean_datasets/cluster_names_VTA_pain.rds", compress = FALSE)

### NAc_TMP ###
NAc_TMP <- DietSeurat(NAc_TMP,
                      layers = c("data"),
                      assays = c("RNA"), #harmony, thus no integrated assay
                      dimreducs = c("pca","umap","harmony_rna"))

# save as new RDS objects
saveRDS(NAc_TMP, file = "./lean_datasets/NAc_TMP.rds", compress = FALSE)

# fetch cluster names and save for UI
## NAc_TMP: current Idents are `CellType`
cluster_names_NAc_TMP<- sort(unique(as.character(Idents(NAc_TMP))))
saveRDS(cluster_names_NAc_TMP, file = "./lean_datasets/cluster_names_NAc_TMP.rds", compress = FALSE)


