# A Daylab RShiny application of single-nuclei datasets
# Author: Lara Ianov | U-BDS
#------------ Configuration of the data --------------------
library(Seurat)
library(Signac)
library(cowplot)
library(shiny)
library(shinyjs)
library(ggplot2)
library(markdown)
library(shinyhelper)
library(dplyr)

# global non-interactive functions, info and datasets to be shared among all 
# connections per worker/process (https://shiny.rstudio.com/articles/scoping.html)
source("./home_description.R", local = TRUE)
source("./helper_functions.R", local = TRUE)
source("./app_modules.R", local = TRUE)

main_panel_style <- "overflow-y:scroll; max-height: 1250px; max-width: 1100px; border-top: solid; border-bottom: solid; border-color: #e8e8e8"

#TODO: consider renaming the first three with project num. as more datasets will be added derived from similar tissue etc.
All_Groups_log_rn6_rn7 <- readRDS(file = "./lean_datasets/All_Groups_log_rn7.rds")
Culture_log_rn6_rn7 <- readRDS(file = "./lean_datasets/Culture_log_rn7.rds")
VTA_dataset_rn6_rn7 <- readRDS(file = "./lean_datasets/VTA_dataset_rn7.rds")
NAc_snATAC_0045 <- readRDS(file = "./lean_datasets/NAc_snATAC_2020_0045.rds")

#--------------------------- ADULT Seurat object ordering ----------------------
# ordering ident for violin
cluster_names_adult <- sort(as.character(unique(All_Groups_log_rn6_rn7@meta.data$CellType)))
Idents(object = All_Groups_log_rn6_rn7) <- factor(Idents(All_Groups_log_rn6_rn7),levels = cluster_names_adult)

#--------------------------- CULTURE Seurat object ordering ----------------------
# ordering ident for violin
cluster_names_cult <- sort(as.character(unique(Culture_log_rn6_rn7@meta.data$CellType)))
Idents(object = Culture_log_rn6_rn7) <- factor(Idents(Culture_log_rn6_rn7), levels = cluster_names_cult)

#--------------------------- VTA adult Seurat object ordering ----------------------
# ordering ident for violin
cluster_names_VTA <- sort(as.character(unique(VTA_dataset_rn6_rn7@meta.data$CellType)))
Idents(object = VTA_dataset_rn6_rn7) <- factor(Idents(VTA_dataset_rn6_rn7),levels = cluster_names_VTA)

#--------------------------- NAc_snATAC_0045 Seurat object ordering ----------------------
cluster_names_NAc_snATAC_0045 <- sort(as.character(unique(NAc_snATAC_0045@meta.data$CellType)))
Idents(object = NAc_snATAC_0045) <- factor(Idents(NAc_snATAC_0045),levels = cluster_names_NAc_snATAC_0045)

#----------------------- app -------------------------------
ui <- function(){
  
  bootstrapPage("",
                useShinyjs(),
                navbarPage(title = "Ratlas",
                           inverse = TRUE,
                           home_description,
                           tabPanel(title = "Adult NAc",
                                    tabsetPanel(type = "tabs",
                                                tabPanel(title = "Adult NAc - rn6",          
                                                         sh_layout_UI(id = "adult",
                                                                      group_choices = adult_groups,
                                                                      plot_choices = all_plots,
                                                                      cluster_names = cluster_names_adult,
                                                                      correlation_label = contains_EES
                                                         )
                                                ),
                                                tabPanel(title = "Adult NAc - rn7",
                                                         sh_layout_UI(id = "adult_rn7",
                                                                      group_choices = adult_groups,
                                                                      plot_choices = all_plots,
                                                                      cluster_names = cluster_names_adult,
                                                                      correlation_label = contains_EES
                                                         )
                                                )
                                    )
                           ),
                           tabPanel(title = "Primary striatal culture",
                                    tabsetPanel(type = "tabs",
                                                tabPanel(title = "Primary striatal culture - rn6",
                                                         sh_layout_UI(id = "culture",
                                                                      group_choices = all_stim_groups,
                                                                      plot_choices = subset_plots,
                                                                      cluster_names = cluster_names_cult,
                                                                      correlation_label = no_EES
                                                         )
                                                ),
                                                tabPanel(title = "Primary striatal culture - rn7",
                                                         sh_layout_UI(id = "culture_rn7",
                                                                      group_choices = all_stim_groups,
                                                                      plot_choices = subset_plots,
                                                                      cluster_names = cluster_names_cult,
                                                                      correlation_label = no_EES
                                                         )
                                                )
                                    )
                           ),
                           tabPanel(title = "Adult VTA",
                                    tabsetPanel(type = "tabs",
                                                tabPanel(title = "Adult VTA - rn6",
                                                         sh_layout_UI(id = "vta",
                                                                      group_choices = all_VTA_groups,
                                                                      plot_choices = subset_plots,
                                                                      cluster_names = cluster_names_VTA,
                                                                      correlation_label = no_EES
                                                         )
                                                ),
                                                tabPanel(title = "Adult VTA - rn7",
                                                         sh_layout_UI(id = "vta_rn7",
                                                                      group_choices = all_VTA_groups,
                                                                      plot_choices = subset_plots,
                                                                      cluster_names = cluster_names_VTA,
                                                                      correlation_label = no_EES
                                                         )
                                                )
                                    )
                           ),

                           tabPanel(title = "Adult NAc - snATAC", #TODO: placeholder name; will need to think of more specific name as there may be other similar datasets in the future
                                    sh_layout_UI(id = "adult_nac_ATAC",
                                                 group_choices = adult_groups_browser,
                                                 plot_choices = atac_plots,
                                                 cluster_names = cluster_names_NAc_snATAC_0045,
                                                 correlation_label = no_EES
                                    )
                           )
                ),
                tags$head(
                  tags$style(HTML(".shiny-output-error-validation {
                color: black;
                                }")))
  )
}

# Reminder: objects inside server function are instantiated per session...
server <- function(input, output) {
  
  shinyhelper::observe_helpers(help_dir = "helpfiles", withMathJax = FALSE)
  
  callModule(sh_layout, id = "adult", 
             dataset = All_Groups_log_rn6_rn7, 
             UMAP_label = "The Rat rn6 NAc")
  
  callModule(sh_layout, id = "adult_rn7", 
             dataset = All_Groups_log_rn6_rn7, 
             UMAP_label = "The Rat rn7 NAc",
             assay = "RNArn7")
  
  callModule(sh_layout, id = "culture", 
             dataset = Culture_log_rn6_rn7, 
             UMAP_label = "Primary striatal neuron culture - rn6",
             EES_absent = "yes")
  
  callModule(sh_layout, id = "culture_rn7", 
             dataset = Culture_log_rn6_rn7, 
             UMAP_label = "Primary striatal neuron culture - rn7",
             EES_absent = "yes",
             assay = "RNArn7")
  
  callModule(sh_layout, id = "vta", 
             dataset = VTA_dataset_rn6_rn7, 
             UMAP_label = "The Rat VTA",
             EES_absent = "yes")
  
  callModule(sh_layout, id = "vta_rn7", 
             dataset = VTA_dataset_rn6_rn7, 
             UMAP_label = "The Rat VTA",
             EES_absent = "yes",
             assay = "RNArn7")
  
  callModule(sh_layout, id = "adult_nac_ATAC", 
             dataset = NAc_snATAC_0045, 
             UMAP_label = "The Rat NAc snATAC (7 consecutive day Veh/cocaine treatment)", #TODO: check with PI this label is good
             EES_absent = "yes",
             assay = "ATAC")
}

shinyApp(ui = ui, server = server)