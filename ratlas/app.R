# A Daylab RShiny application of single-nuclei datasets
# Author: Lara Ianov | U-BDS

source("./global.R")

#----------------------- app -------------------------------
ui <- function(){
  
  bootstrapPage("",
                useShinyjs(),
                navbarPage(title = "Ratlas",
                           theme = bslib::bs_theme(version = 5, bootswatch = "cosmo", primary = "#232a30"),
                           home_description,
                           tabPanel(title = "Adult acute NAc",
                                    tabsetPanel(type = "tabs",
                                                tabPanel(title = "Adult NAc - rn6",          
                                                         sh_layout_UI(id = "adult",
                                                                      group_choices = adult_acute_groups,
                                                                      plot_choices = all_plots_EES,
                                                                      cluster_names = cluster_names_adult,
                                                                      correlation_label = contains_EES
                                                         )
                                                ),
                                                tabPanel(title = "Adult NAc - rn7",
                                                         sh_layout_UI(id = "adult_rn7",
                                                                      group_choices = adult_acute_groups,
                                                                      plot_choices = all_plots_EES,
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
                                                                      group_choices = adult_acute_repeated_groups,
                                                                      plot_choices = all_plots,
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
                                                                      plot_choices = all_plots,
                                                                      cluster_names = cluster_names_cult,
                                                                      correlation_label = no_EES
                                                         )
                                                ),
                                                tabPanel(title = "Primary striatal culture - rn7",
                                                         sh_layout_UI(id = "culture_rn7",
                                                                      group_choices = all_stim_groups,
                                                                      plot_choices = all_plots,
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
                                                                      plot_choices = all_plots,
                                                                      cluster_names = cluster_names_VTA,
                                                                      correlation_label = no_EES
                                                         )
                                                ),
                                                tabPanel(title = "Adult VTA - rn7",
                                                         sh_layout_UI(id = "vta_rn7",
                                                                      group_choices = all_VTA_groups,
                                                                      plot_choices = all_plots,
                                                                      cluster_names = cluster_names_VTA,
                                                                      correlation_label = no_EES
                                                         )
                                                )
                                    )
                           )
                ),
                tags$style(HTML(".irs--shiny .irs-bar {
                                background: #232a30;
                                border-top: 1px solid #232a30;
                                border-bottom: 1px solid #232a30;
                                }
                                .irs--shiny .irs-to, .irs--shiny .irs-from {
                                background-color: #232a30;
                                }
                                .irs--shiny .irs-single {
                                background: #232a30;
                                }")),
                tags$head(
                  tags$style(HTML(".shiny-output-error-validation {
                                  color: black;}"))))
}

# Reminder: objects inside server function are instantiated per session...
server <- function(input, output) {
  
  shinyhelper::observe_helpers(help_dir = "helpfiles", withMathJax = FALSE)
  
  sh_layout_server(id = "adult", 
                   dataset = All_Groups_log_rn6_rn7, 
                   UMAP_label = "The Rat acute NAc dataset - rn6")
  
  sh_layout_server(id = "adult_rn7", 
                   dataset = All_Groups_log_rn6_rn7, 
                   UMAP_label = "The Rat acute NAc dataset - rn7",
                   assay = "RNArn7")
  
  sh_layout_server(id = "adult_mcn", 
                   dataset = MCN_dataset, 
                   UMAP_label = "The Rat acute and repeated NAc dataset - rn7",
                   EES_absent = TRUE)
  
  sh_layout_server(id = "culture", 
                   dataset = Culture_log_rn6_rn7, 
                   UMAP_label = "Primary striatal neuron culture - rn6",
                   EES_absent = TRUE)
  
  sh_layout_server(id = "culture_rn7", 
                   dataset = Culture_log_rn6_rn7, 
                   UMAP_label = "Primary striatal neuron culture - rn7",
                   EES_absent = TRUE,
                   assay = "RNArn7")
  
  sh_layout_server(id = "vta", 
                   dataset = VTA_dataset_rn6_rn7, 
                   UMAP_label = "The Rat VTA dataset - rn6",
                   EES_absent = TRUE)
  
  sh_layout_server(id = "vta_rn7", 
                   dataset = VTA_dataset_rn6_rn7, 
                   UMAP_label = "The Rat VTA dataset - rn7",
                   EES_absent = TRUE,
                   assay = "RNArn7")
}

shinyApp(ui = ui, server = server)