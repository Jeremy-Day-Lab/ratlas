# A Daylab RShiny application of single-nuclei datasets
# Author: Lara Ianov | U-BDS

source("./global.R")

#----------------------- app -------------------------------
ui <- function(){

  bootstrapPage("",
                useShinyjs(),
                navbarPage(title = "Ratlas",
                           id = "main_nav",  # top-level navigation input
                           theme = bslib::bs_theme(version = 5, bootswatch = "cosmo", primary = "#232a30"),
                           home_description,
                           tabPanel(title = "Adult acute NAc", value = "nav_adult",
                                    tabsetPanel(id = "dataset_tabs_adult", type = "tabs",
                                                tabPanel(title = "Adult NAc - rn6", value = "adult_rn6_tab",
                                                         sh_layout_UI(id = "adult",
                                                                      group_choices = adult_acute_groups,
                                                                      plot_choices = all_plots_EES,
                                                                      cluster_names = cluster_names_adult,
                                                                      correlation_label = contains_EES
                                                         )
                                                ),
                                                tabPanel(title = "Adult NAc - rn7", value = "adult_rn7_tab",
                                                         sh_layout_UI(id = "adult_rn7",
                                                                      group_choices = adult_acute_groups,
                                                                      plot_choices = all_plots_EES,
                                                                      cluster_names = cluster_names_adult,
                                                                      correlation_label = contains_EES
                                                         )
                                                )
                                    )
                           ),
                           tabPanel(title = "Adult acute and repeated NAc", value = "nav_mcn",
                                    tabsetPanel(id = "dataset_tabs_mcn", type = "tabs",
                                                tabPanel(title = "Adult acute and repeated NAc - rn7", value = "adult_mcn_tab", # keeping tabs for consistency for now
                                                         sh_layout_UI(id = "adult_mcn",
                                                                      group_choices = adult_acute_repeated_groups,
                                                                      plot_choices = all_plots,
                                                                      cluster_names = cluster_names_MCN,
                                                                      correlation_label = no_EES
                                                         )
                                                )
                                    )
                           ),
                           tabPanel(title = "Primary striatal culture", value = "nav_culture",
                                    tabsetPanel(id = "dataset_tabs_culture", type = "tabs",
                                                tabPanel(title = "Primary striatal culture - rn6", value = "culture_rn6_tab",
                                                         sh_layout_UI(id = "culture",
                                                                      group_choices = all_stim_groups,
                                                                      plot_choices = all_plots,
                                                                      cluster_names = cluster_names_cult,
                                                                      correlation_label = no_EES
                                                         )
                                                ),
                                                tabPanel(title = "Primary striatal culture - rn7", value = "culture_rn7_tab",
                                                         sh_layout_UI(id = "culture_rn7",
                                                                      group_choices = all_stim_groups,
                                                                      plot_choices = all_plots,
                                                                      cluster_names = cluster_names_cult,
                                                                      correlation_label = no_EES
                                                         )
                                                )
                                    )
                           ),
                           tabPanel(title = "Adult VTA", value = "nav_vta",
                                    tabsetPanel(id = "dataset_tabs_VTA", type = "tabs",
                                                tabPanel(title = "Adult VTA - rn6", value = "vta_rn6_tab",
                                                         sh_layout_UI(id = "vta",
                                                                      group_choices = all_VTA_groups,
                                                                      plot_choices = all_plots,
                                                                      cluster_names = cluster_names_VTA,
                                                                      correlation_label = no_EES
                                                         )
                                                ),
                                                tabPanel(title = "Adult VTA - rn7", value = "vta_rn7_tab",
                                                         sh_layout_UI(id = "vta_rn7",
                                                                      group_choices = all_VTA_groups,
                                                                      plot_choices = all_plots,
                                                                      cluster_names = cluster_names_VTA,
                                                                      correlation_label = no_EES
                                                         )
                                                )
                                    )
                           ),
                           tabPanel(title = "VTA Pain", value = "nav_vta_pain",
                                    tabsetPanel(id = "dataset_tabs_VTA_pain", type = "tabs",
                                                tabPanel(title = "VTA pain - rn7", value = "vta_pain_tab",
                                                         sh_layout_UI(id = "vta_pain",
                                                                      group_choices = all_VTA_pain_groups,
                                                                      plot_choices = all_plots,
                                                                      cluster_names = cluster_names_VTA_pain,
                                                                      correlation_label = no_EES
                                                         )
                                                )
                                    )
                           ),
                           tabPanel(title = "NAc_2026", value = "nav_nac_2026",
                                    tabsetPanel(id = "dataset_tabs_NAc_2026", type = "tabs",
                                                tabPanel(title = "NAc_2026 - rn7", value = "nac_2026_tab",
                                                         sh_layout_UI(id = "nac_2026",
                                                                      group_choices = NAc_2026_groups,
                                                                      plot_choices = all_plots,
                                                                      cluster_names = cluster_names_NAc_2026,
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

  # env created fresh inside server().
  wired <- new.env(parent = emptyenv())

  observeEvent(input$main_nav, {

    nav <- input$main_nav

    # ---- Adult acute NAc ----
    if (nav == "nav_adult") {

      if (is.null(All_Groups_log_rn6_rn7)) {
        All_Groups_log_rn6_rn7 <<- readRDS(file = All_Groups_log_rn6_rn7_path)
        Idents(object = All_Groups_log_rn6_rn7) <<- factor(Idents(All_Groups_log_rn6_rn7), levels = cluster_names_adult)
      }

      if (is.null(wired[["nav_adult"]])) {
        wired[["nav_adult"]] <- TRUE

        sh_layout_server(id = "adult",
                         dataset = All_Groups_log_rn6_rn7,
                         UMAP_label = "The Rat acute NAc dataset - rn6",
                         cluster_names = cluster_names_adult)

        sh_layout_server(id = "adult_rn7",
                         dataset = All_Groups_log_rn6_rn7,
                         UMAP_label = "The Rat acute NAc dataset - rn7",
                         cluster_names = cluster_names_adult,
                         assay = "RNArn7")
      }
    }

    # ---- Adult acute and repeated NAc ----
    if (nav == "nav_mcn") {

      if (is.null(MCN_dataset)) {
        MCN_dataset <<- readRDS(file = MCN_dataset_path)
        Idents(object = MCN_dataset) <<- factor(Idents(MCN_dataset), levels = cluster_names_MCN)
      }

      if (is.null(wired[["nav_mcn"]])) {
        wired[["nav_mcn"]] <- TRUE

        sh_layout_server(id = "adult_mcn",
                         dataset = MCN_dataset,
                         UMAP_label = "The Rat acute and repeated NAc dataset - rn7",
                         cluster_names = cluster_names_MCN,
                         EES_absent = TRUE)
      }
    }

    # ---- Primary striatal culture ----
    if (nav == "nav_culture") {

      if (is.null(Culture_log_rn6_rn7)) {
        Culture_log_rn6_rn7 <<- readRDS(file = Culture_log_rn6_rn7_path)
        Idents(object = Culture_log_rn6_rn7) <<- factor(Idents(Culture_log_rn6_rn7), levels = cluster_names_cult)
      }
      
      if (is.null(wired[["nav_culture"]])) {
        wired[["nav_culture"]] <- TRUE

        sh_layout_server(id = "culture",
                         dataset = Culture_log_rn6_rn7,
                         UMAP_label = "Primary striatal neuron culture - rn6",
                         cluster_names = cluster_names_cult,
                         EES_absent = TRUE)

        sh_layout_server(id = "culture_rn7",
                         dataset = Culture_log_rn6_rn7,
                         UMAP_label = "Primary striatal neuron culture - rn7",
                         cluster_names = cluster_names_cult,
                         EES_absent = TRUE,
                         assay = "RNArn7")
      }
    }

    # ---- Adult VTA ----
    if (nav == "nav_vta") {

      if (is.null(VTA_dataset_rn6_rn7)) {
        VTA_dataset_rn6_rn7 <<- readRDS(file = VTA_dataset_rn6_rn7_path)
        Idents(object = VTA_dataset_rn6_rn7) <<- factor(Idents(VTA_dataset_rn6_rn7), levels = cluster_names_VTA)
      }

      if (is.null(wired[["nav_vta"]])) {
        wired[["nav_vta"]] <- TRUE

        sh_layout_server(id = "vta",
                         dataset = VTA_dataset_rn6_rn7,
                         UMAP_label = "The Rat VTA dataset - rn6",
                         cluster_names = cluster_names_VTA,
                         EES_absent = TRUE)

        sh_layout_server(id = "vta_rn7",
                         dataset = VTA_dataset_rn6_rn7,
                         UMAP_label = "The Rat VTA dataset - rn7",
                         cluster_names = cluster_names_VTA,
                         EES_absent = TRUE,
                         assay = "RNArn7")
      }
    }

    # ---- VTA Pain ----
    if (nav == "nav_vta_pain") {

      if (is.null(VTA_pain_dataset)) {
        VTA_pain_dataset <<- readRDS(file = VTA_pain_dataset_path)
        Idents(object = VTA_pain_dataset) <<- factor(Idents(VTA_pain_dataset), levels = cluster_names_VTA_pain)
      }

      if (is.null(wired[["nav_vta_pain"]])) {
        wired[["nav_vta_pain"]] <- TRUE

        sh_layout_server(id = "vta_pain",
                         dataset = VTA_pain_dataset,
                         UMAP_label = "The Rat VTA Pain dataset - rn7",
                         cluster_names = cluster_names_VTA_pain,
                         EES_absent = TRUE)
      }
    }

    # ---- NAc_2026 ----
    if (nav == "nav_nac_2026") {

      if (is.null(NAc_2026_dataset)) {
        NAc_2026_dataset <<- readRDS(file = NAc_2026_dataset_path)
        Idents(object = NAc_2026_dataset) <<- factor(Idents(NAc_2026_dataset), levels = cluster_names_NAc_2026)
      }

      if (is.null(wired[["nav_nac_2026"]])) {
        wired[["nav_nac_2026"]] <- TRUE

        sh_layout_server(id = "nac_2026",
                         dataset = NAc_2026_dataset,
                         UMAP_label = "The Rat NAc_2026 dataset - rn7",
                         cluster_names = cluster_names_NAc_2026,
                         EES_absent = TRUE)
      }
    }

  }, ignoreInit = TRUE)
}

shinyApp(ui = ui, server = server)
