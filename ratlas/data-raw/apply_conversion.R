# apply convert_assay_v5 across all existing Ratlas datasets

library(BPCells)
library(Seurat)
library(SeuratObject)

options(Seurat.object.assay.version = "v5")
source("./data-raw/convert_assay_v5.R")

# load all objects

All_Groups_log_rn6_rn7 <- readRDS(file = "./lean_datasets/All_Groups_log_rn7.rds")
Culture_log_rn6_rn7 <- readRDS(file = "./lean_datasets/Culture_log_rn7.rds")
VTA_dataset_rn6_rn7 <- readRDS(file = "./lean_datasets/VTA_dataset_rn7.rds")
MCN_dataset <- readRDS(file = "./lean_datasets/MCN_dataset.rds") # rn7 only

# run convert_assay_v5 to save updated Seurat objects.
# this is done once per assay to be kept (thus in some cases twice - RNA and RNArn7)
# NOTE: for cases where we do 2 conversions the on-disk data directory is copied
# into the latest run (only need to retain the latest run)

#------------------------------ convert --------------------------------------
# calling is separately - each run has a different name  etc... and a one timer, so keeping it simple

### All_Groups_log_rn6_rn7 ####

# RNA
convert_assay_v5(dataset = All_Groups_log_rn6_rn7, assay = "RNA",
                 tmp_output_path = "./lean_datasets/TMPS/All_Groups_log_rn6_rn7_RNA",
                 seurat_output_path = "./lean_datasets/TMPS/All_Groups_log_rn6_rn7_TMP_RNA", # RNA assay, tmp
                 updated_dataset_name = "All_Groups_log_rn6_rn7_assay5_TMP.rds")

# remove current obj from the environment, 
# and load the one created above for re-run on rn7 assay

rm(All_Groups_log_rn6_rn7)
All_Groups_log_rn6_rn7 <- readRDS("./lean_datasets/TMPS/All_Groups_log_rn6_rn7_TMP_RNA/All_Groups_log_rn6_rn7_assay5_TMP.rds")

# RNArn7
convert_assay_v5(dataset = All_Groups_log_rn6_rn7, assay = "RNArn7",
                 tmp_output_path = "./lean_datasets/TMPS/All_Groups_log_rn6_rn7_RNArn7",
                 seurat_output_path = "./lean_datasets/TMPS/All_Groups_log_rn6_rn7_TMP_RNArn7", # RNArn7 assay, tmp
                 updated_dataset_name = "All_Groups_log_rn6_rn7_assay5_TMP.rds")

rm(All_Groups_log_rn6_rn7)
All_Groups_log_rn6_rn7 <- readRDS("./lean_datasets/TMPS/All_Groups_log_rn6_rn7_TMP_RNArn7/All_Groups_log_rn6_rn7_assay5_TMP.rds")

# integrated
convert_assay_v5(dataset = All_Groups_log_rn6_rn7, assay = "integrated",
                 tmp_output_path = "./lean_datasets/TMPS/All_Groups_log_rn6_rn7_integrated",
                 seurat_output_path = "./lean_datasets/on_disk/All_Groups_log_rn6_rn7", # RNA + RNArn7 + integrated assay, final
                 updated_dataset_name = "All_Groups_log_rn6_rn7_assay5.rds")

### Culture_log_rn6_rn7 ####

# RNA
convert_assay_v5(dataset = Culture_log_rn6_rn7, assay = "RNA",
                 tmp_output_path = "./lean_datasets/TMPS/Culture_log_rn6_rn7_RNA",
                 seurat_output_path = "./lean_datasets/TMPS/Culture_log_rn6_rn7_TMP_RNA", # RNA assay, tmp
                 updated_dataset_name = "Culture_log_rn6_rn7_assay5_TMP.rds")

# remove current obj from the environment, 
# and load the one created above for re-run on rn7 assay

rm(Culture_log_rn6_rn7)
Culture_log_rn6_rn7 <- readRDS("./lean_datasets/TMPS/Culture_log_rn6_rn7_TMP_RNA/Culture_log_rn6_rn7_assay5_TMP.rds")

# RNArn7
convert_assay_v5(dataset = Culture_log_rn6_rn7, assay = "RNArn7",
                 tmp_output_path = "./lean_datasets/TMPS/Culture_log_rn6_rn7_RNArn7",
                 seurat_output_path = "./lean_datasets/TMPS/Culture_log_rn6_rn7_TMP_RNArn7", # RNArn7 assay, tmp
                 updated_dataset_name = "Culture_log_rn6_rn7_assay5_TMP.rds")

rm(Culture_log_rn6_rn7)
Culture_log_rn6_rn7 <- readRDS("./lean_datasets/TMPS/Culture_log_rn6_rn7_TMP_RNArn7/Culture_log_rn6_rn7_assay5_TMP.rds")

# integrated
convert_assay_v5(dataset = Culture_log_rn6_rn7, assay = "integrated",
                 tmp_output_path = "./lean_datasets/TMPS/Culture_log_rn6_rn7_integrated",
                 seurat_output_path = "./lean_datasets/on_disk/Culture_log_rn6_rn7", # RNA + RNArn7 + integrated assay, final
                 updated_dataset_name = "Culture_log_rn6_rn7_assay5.rds")

### VTA_dataset_rn6_rn7 ####

# RNA
convert_assay_v5(dataset = VTA_dataset_rn6_rn7, assay = "RNA",
                 tmp_output_path = "./lean_datasets/TMPS/VTA_dataset_rn6_rn7_RNA",
                 seurat_output_path = "./lean_datasets/TMPS/VTA_dataset_rn6_rn7_TMP_RNA", # RNA assay, tmp
                 updated_dataset_name = "VTA_dataset_rn6_rn7_assay5_TMP.rds")

# remove current obj from the environment, 
# and load the one created above for re-run on rn7 assay

rm(VTA_dataset_rn6_rn7)
VTA_dataset_rn6_rn7 <- readRDS("./lean_datasets/TMPS/VTA_dataset_rn6_rn7_TMP_RNA/VTA_dataset_rn6_rn7_assay5_TMP.rds")

# RNArn7
convert_assay_v5(dataset = VTA_dataset_rn6_rn7, assay = "RNArn7",
                 tmp_output_path = "./lean_datasets/TMPS/VTA_dataset_rn6_rn7_RNArn7",
                 seurat_output_path = "./lean_datasets/TMPS/VTA_dataset_rn6_rn7_TMP_RNArn7", # RNArn7 assay, tmp
                 updated_dataset_name = "VTA_dataset_rn6_rn7_assay5_TMP.rds")

rm(VTA_dataset_rn6_rn7)
VTA_dataset_rn6_rn7 <- readRDS("./lean_datasets/TMPS/VTA_dataset_rn6_rn7_TMP_RNArn7/VTA_dataset_rn6_rn7_assay5_TMP.rds")

# integrated
convert_assay_v5(dataset = VTA_dataset_rn6_rn7, assay = "integrated",
                 tmp_output_path = "./lean_datasets/TMPS/VTA_dataset_rn6_rn7_integrated",
                 seurat_output_path = "./lean_datasets/on_disk/VTA_dataset_rn6_rn7", # RNA + RNArn7 + integrated assay, final
                 updated_dataset_name = "VTA_dataset_rn6_rn7_assay5.rds")
### MCN_dataset ####

# RNA (rn7 only)
convert_assay_v5(dataset = MCN_dataset, assay = "RNA",
                 tmp_output_path = "./lean_datasets/TMPS/MCN_dataset_RNA",
                 seurat_output_path = "./lean_datasets/TMPS/MCN_dataset_TMP_RNA", # RNA only
                 updated_dataset_name = "MCN_dataset_assay5_TMP.rds")

rm(MCN_dataset)
MCN_dataset <- readRDS("./lean_datasets/TMPS/MCN_dataset_TMP_RNA/MCN_dataset_assay5_TMP.rds")

# integrated
convert_assay_v5(dataset = MCN_dataset, assay = "integrated",
                 tmp_output_path = "./lean_datasets/TMPS/MCN_dataset_integrated",
                 seurat_output_path = "./lean_datasets/on_disk/MCN_dataset", # RNA + integrated, final
                 updated_dataset_name = "MCN_dataset_assay5.rds")


