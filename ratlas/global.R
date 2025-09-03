#-------------------------------packages------------------------------------------------------
library(Seurat)
##library(BPCells)
library(shiny)
library(shinyjs)
library(ggplot2)
library(markdown)
library(shinyhelper)
library(dplyr)
library(bslib)

# Seurat V5 assay
options(Seurat.object.assay.version = "v5")

#--------------------------custom functions---------------------------------------------------

lapply(list.files("./R"), FUN = function(x) source(paste0("./R/", x)))

#--------------------------global objects/variables-------------------------------------------
main_panel_style <- "overflow-y:scroll; max-height: 1800px; border-top: solid; border-bottom: solid; border-color: #e8e8e8"

#TODO: consider renaming the first three with project num. as more datasets will be added derived from similar tissue etc.
All_Groups_log_rn6_rn7_path <- "./lean_datasets/All_Groups_log_rn7.rds"
MCN_dataset_path <- "./lean_datasets/MCN_dataset.rds" # rn7 only
Culture_log_rn6_rn7_path <- "./lean_datasets/Culture_log_rn7.rds"
VTA_dataset_rn6_rn7_path <- "./lean_datasets/VTA_dataset_rn7.rds"

# global variables for datasets

All_Groups_log_rn6_rn7 <- NULL
MCN_dataset <- NULL
Culture_log_rn6_rn7 <- NULL
VTA_dataset_rn6_rn7 <- NULL

#--------------------------- ADULT Seurat object ordering ----------------------
# ordering ident for violin
cluster_names_adult <- readRDS(file = "./lean_datasets/cluster_names_adult.rds")

#Idents(object = All_Groups_log_rn6_rn7) <- factor(Idents(All_Groups_log_rn6_rn7),levels = cluster_names_adult)

#--------------------------- MCN Seurat object ordering ----------------------
# ordering ident for violin
cluster_names_MCN <- readRDS(file = "./lean_datasets/cluster_names_MCN.rds")
#Idents(object = MCN_dataset) <- factor(Idents(MCN_dataset),levels = cluster_names_MCN)

#--------------------------- CULTURE Seurat object ordering ----------------------
# ordering ident for violin
cluster_names_cult <- readRDS(file = "./lean_datasets/cluster_names_cult.rds")
#Idents(object = Culture_log_rn6_rn7) <- factor(Idents(Culture_log_rn6_rn7), levels = cluster_names_cult)

#--------------------------- VTA adult Seurat object ordering ----------------------
# ordering ident for violin
cluster_names_VTA <- readRDS(file = "./lean_datasets/cluster_names_VTA.rds")
#Idents(object = VTA_dataset_rn6_rn7) <- factor(Idents(VTA_dataset_rn6_rn7),levels = cluster_names_VTA)

