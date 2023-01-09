#-------------------------------packages------------------------------------------------------
library(Seurat)
library(cowplot)
library(shiny)
library(shinyjs)
library(ggplot2)
library(markdown)
library(shinyhelper)
library(dplyr)

#--------------------------custom functions---------------------------------------------------

lapply(list.files("./R"), FUN = function(x) source(paste0("./R/", x)))

#--------------------------global objects/variables-------------------------------------------
#TODO: increase width
main_panel_style <- "overflow-y:scroll; max-height: 1250px; max-width: 1100px; border-top: solid; border-bottom: solid; border-color: #e8e8e8"

#TODO: consider renaming the first three with project num. as more datasets will be added derived from similar tissue etc.
All_Groups_log_rn6_rn7 <- readRDS(file = "./lean_datasets/All_Groups_log_rn7.rds")
MCN_dataset <- readRDS(file = "./lean_datasets/MCN_dataset.rds") # rn7 only
Culture_log_rn6_rn7 <- readRDS(file = "./lean_datasets/Culture_log_rn7.rds")
VTA_dataset_rn6_rn7 <- readRDS(file = "./lean_datasets/VTA_dataset_rn7.rds")

#--------------------------- ADULT Seurat object ordering ----------------------
# ordering ident for violin
cluster_names_adult <- sort(as.character(unique(All_Groups_log_rn6_rn7@meta.data$CellType)))
Idents(object = All_Groups_log_rn6_rn7) <- factor(Idents(All_Groups_log_rn6_rn7),levels = cluster_names_adult)

#--------------------------- MCN Seurat object ordering ----------------------
# ordering ident for violin
cluster_names_MCN <- sort(as.character(unique(MCN_dataset@meta.data$CellType)))
Idents(object = MCN_dataset) <- factor(Idents(MCN_dataset),levels = cluster_names_MCN)

#--------------------------- CULTURE Seurat object ordering ----------------------
# ordering ident for violin
cluster_names_cult <- sort(as.character(unique(Culture_log_rn6_rn7@meta.data$CellType)))
Idents(object = Culture_log_rn6_rn7) <- factor(Idents(Culture_log_rn6_rn7), levels = cluster_names_cult)

#--------------------------- VTA adult Seurat object ordering ----------------------
# ordering ident for violin
cluster_names_VTA <- sort(as.character(unique(VTA_dataset_rn6_rn7@meta.data$CellType)))
Idents(object = VTA_dataset_rn6_rn7) <- factor(Idents(VTA_dataset_rn6_rn7),levels = cluster_names_VTA)
