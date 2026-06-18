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

### NAc_2026 ###
# NOTE: object already V5, no need to convert beyond this point, just slim down the obj

NAc_2026 <- readRDS(file = "./analysis_objects/5_NAc0069_integrated_clean_NAcOnly.rds")

## Modify a few metadata variables for consistency and add new ones
NAc_2026$CellType <- NAc_2026$cellType
NAc_2026$Sex <- NAc_2026$sex

# set order level for Sex:
NAc_2026$Sex <- factor(NAc_2026$Sex, 
                       levels=c("female", "male"))

## remove un-needed variables:
NAc_2026$sex <- NULL
NAc_2026$treatment <- NULL
NAc_2026$cellType <- NULL
NAc_2026$scDblFinder_doublets <- NULL

## rename reduction to match others:

NAc_2026@reductions[["umap"]] <- NAc_2026@reductions[["rnaharmony.umap"]]
NAc_2026@reductions[["rnaharmony.umap"]] <- NULL

# Set a new key for updated name
Key(NAc_2026@reductions[["umap"]]) <- "UMAP_"

# Make the embedding column names match the new key
colnames(NAc_2026@reductions[["umap"]]@cell.embeddings) <- paste0("UMAP_", 1:ncol(NAc_2026@reductions[["umap"]]@cell.embeddings))

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

### NAc_2026 ###
NAc_2026 <- DietSeurat(NAc_2026,
                      layers = c("data"),
                      assays = c("RNA"), #harmony, thus no integrated assay
                      dimreducs = c("pca","umap","harmony_rna"))

# save as new RDS objects
saveRDS(NAc_2026, file = "./lean_datasets/NAc_2026.rds", compress = FALSE)

# fetch cluster names and save for UI
## NAc_2026: current Idents are `CellType`
cluster_names_NAc_2026<- sort(unique(as.character(Idents(NAc_2026))))
saveRDS(cluster_names_NAc_2026, file = "./lean_datasets/cluster_names_NAc_2026.rds", compress = FALSE)


