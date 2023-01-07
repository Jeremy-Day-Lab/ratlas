# A Daylab RShiny application of single-nuclei datasets
# Author: Lara Ianov | U-BDS

source("./global.R")

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
                                                                      group_choices = adult_groups, #TODO: rename this
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
                           tabPanel(title = "Adult acute and repeated NAc", #TODO: double check name
                                    tabsetPanel(type = "tabs",
                                                tabPanel(title = "Adult acute and repeated NAc - rn7", # keeping tabs for consistency for now
                                                         sh_layout_UI(id = "adult_mcn",
                                                                      group_choices = adult_groups,
                                                                      plot_choices = subset_plots,
                                                                      cluster_names = cluster_names_MCN,
                                                                      correlation_label = no_EES
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
  
  callModule(sh_layout, id = "adult_mcn", 
             dataset = MCN_dataset, 
             UMAP_label = "The Rat rn7 acute and repeated NAc dataset",
             EES_absent = "yes")
  
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
}

shinyApp(ui = ui, server = server)