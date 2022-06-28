# app modules
#-----------------------------------------global UI options------------------------------------------
options(spinner.type=1,spinner.color="#232a30", spinner.size=2)

#-------------------------------------------UI-----------------------------------------------------

##### group choices

# adult group choices
adult_groups <- c("All", "Stim","Sex","Stim_Sex")
adult_groups_browser <- c("All", "CellType.Stim", "CellType.Sex")

# choices for culture (and drd-1 subclustering if we add this)
all_stim_groups <- c("All", "Stim")

# choices for adult VTA
all_VTA_groups <- c("All", "Sex")

# correlation plot UI message based EES data being present
contains_EES <- "Type a gene or EES to correlate to gene name typed above"
no_EES <- "Type a gene to correlate to gene name selected above"

##### plots available

all_plots <- c("UMAP","FeaturePlot","Violin", "EES_FeaturePlot","EES_Violin", "Correlation_plot")
subset_plots <- c("UMAP","FeaturePlot","Violin","Correlation_plot")
atac_plots <- c("UMAP", "CoveragePlot", "CoveragePlot_genome")

sh_layout_UI <- function(id, group_choices, plot_choices, cluster_names, correlation_label) {
    ns <- NS(id)
  
    sidebarLayout(
        sidebarPanel(width = 3,
                     textInput(inputId = ns("gene"),
                               label = "Choose a gene",
                               placeholder = "Gad1") %>%
                       shinyhelper::helper(icon = "info-circle",
                                           colour ="#232a30",
                                           type = "markdown",
                                           content = "Gene_name_help"
                       ),
                     
                     actionButton(inputId = ns("go"),
                                  label = "Update gene!",
                                  icon("dna"),
                                  style="color: #ededed; background-color: #232a30"),
                     
                     selectInput(inputId = ns("group"),
                                 label = "Choose all or split by group",
                                 choices = group_choices,
                                 selected = "All",
                                 multiple = FALSE),
                     
                     checkboxGroupInput(inputId = ns("plots"),
                                        label = "Choose which plot(s) to display",
                                        choices = plot_choices,
                                        selected = c("UMAP"),
                                        inline = TRUE),
                     hr(),
            
                    # conditional panels based on user choice of plots
                    conditionalPanel(
                        condition = "input.plots.indexOf('FeaturePlot') > -1",
                        numericInput(inputId = ns("expression"),
                                     label = "FeaturePlot option: choose a maximum cutoff for expression scale",
                                     value = NA,
                                     min = 0) %>%
                        shinyhelper::helper(icon = "exclamation-circle",
                                            colour ="#232a30",
                                            type = "markdown",
                                            content = "Featureplot_cutoff_help"
                        ),

                        actionButton(inputId = ns("reset"),
                                     label = "Reset scale",
                                     icon("redo"),
                                     style="color: #ededed; background-color: #232a30"),
                        hr(),
                        ns = ns),
                    
                    conditionalPanel(
                        condition = "input.plots.indexOf('Violin') > -1 || input.plots.indexOf('EES_Violin') > -1",
                        checkboxGroupInput(inputId = ns("cluster"),
                                           label = "Violin plot option: choose all or show specific clusters",
                                           choices = cluster_names,
                                           selected = cluster_names),
                        
                        actionButton(inputId = ns("cluster_selection"),
                                     label = "Plot selected clusters",
                                     icon("check"),
                                     style="color: #ededed; background-color: #232a30") %>%
                          shinyhelper::helper(icon = "info-circle",
                                              colour ="#232a30",
                                              type = "markdown",
                                              content = "Plot_selected_clusters"
                          ),

                        actionButton(inputId = ns("reset_clusters"),
                                     label = "Clear all cluster choices",
                                     icon("redo"),
                                     style="color: #ededed; background-color: #232a30"),

                        actionButton(inputId = ns("select_all_clusters"),
                                     label = "Check all cluster choices",
                                     icon("check-double"),
                                     style="color: #ededed; background-color: #232a30"),
                        
                        hr(),
                        ns = ns),
        
                    conditionalPanel(
                      condition = "input.plots.indexOf('Correlation_plot') > -1",
                      selectInput(inputId = ns("cluster_corr"),
                                  label = "Choose cell type for correlation between features",
                                  choices = cluster_names,
                                  selected = "Drd1-MSN",
                                  multiple = FALSE),
                      
                      textInput(inputId = ns("feature_corr"),
                                label = correlation_label,
                                placeholder = "Sik3") %>%
                      shinyhelper::helper(icon = "info-circle",
                                          colour ="#232a30",
                                          type = "markdown",
                                          content = "Gene_name_help"
                      ),
                      
                      actionButton(inputId = ns("go_corr"),
                                   label = "Update correlation feature",
                                   icon("chart-line"),
                                   style="color: #ededed; background-color: #232a30"),
                      ns = ns),
                    conditionalPanel(
                      condition = "input.plots.indexOf('CoveragePlot') > -1 || input.plots.indexOf('CoveragePlot_genome') > -1",

                      checkboxGroupInput(inputId = ns("cluster_browser"),
                                         label = "CoveragePlot plot option: choose all or show specific clusters",
                                         choices = cluster_names,
                                         selected = cluster_names) %>%
                        shinyhelper::helper(icon = "info-circle",
                                            colour ="#232a30",
                                            type = "markdown",
                                            content = "CoveragePlot_help"
                        ),
                      
                      actionButton(inputId = ns("cluster_selection_browser"),
                                   label = "Plot selected clusters",
                                   icon("check"),
                                   style="color: #ededed; background-color: #232a30") %>%
                        shinyhelper::helper(icon = "info-circle",
                                            colour ="#232a30",
                                            type = "markdown",
                                            content = "Plot_selected_clusters"
                        ),
                      
                      actionButton(inputId = ns("reset_browser_clusters"),
                                   label = "Clear all cluster choices",
                                   icon("redo"),
                                   style="color: #ededed; background-color: #232a30"),
                      
                      actionButton(inputId = ns("select_all_browser_clusters"),
                                   label = "Check all cluster choices",
                                   icon("check-double"),
                                   style="color: #ededed; background-color: #232a30"),
                      
                      selectInput(inputId = ns("show_bulk"),
                                  label = "Plot pseudo-bulk track",
                                  choices = c(TRUE, FALSE),
                                  selected = FALSE,
                                  multiple = FALSE),
                      hr(),
                      ns = ns),
                    conditionalPanel(
                      condition = "input.plots.indexOf('CoveragePlot') > -1",
                      
                      selectInput(inputId = ns("extend_window"),
                                  label = "Choose number of bases to extend the region down/upstream of gene. Max: 5000",
                                  choices = c(0,1000,2000,3000,4000,5000),
                                  selected = 1000,
                                  multiple = FALSE),
                      
                      hr(),
                      ns = ns),
                    conditionalPanel(
                      condition = "input.plots.indexOf('CoveragePlot_genome') > -1",
                      textInput(inputId = ns("genome_coords"),
                                label = "Choose a genomic coordinate instead or in addition to a gene name",
                                placeholder = "1-172934-175664"),
                      
                      actionButton(inputId = ns("go_cords"),
                                   label = "Update coordinates!",
                                   icon("dna"),
                                   style="color: #ededed; background-color: #232a30"),
                      hr(),
                      ns = ns),
                      ),
        
        mainPanel(width = 9, style = main_panel_style,
                  conditionalPanel(
                    condition = "input.plots.indexOf('UMAP') > -1",
                    plotOutput(ns("UMAP")) %>% 
                      shinycssloaders::withSpinner(),
                    ns = ns),
                  conditionalPanel(
                    condition = "input.plots.indexOf('EES_FeaturePlot') > -1",
                    plotOutput(ns("EES_FeaturePlot")) %>% 
                      shinycssloaders::withSpinner(),
                    ns = ns),
                  conditionalPanel(
                    condition = "input.plots.indexOf('EES_Violin') > -1",
                    plotOutput(ns("EES_Violin")) %>% 
                      shinycssloaders::withSpinner(),
                    ns = ns),
                  conditionalPanel(
                    condition = "input.plots.indexOf('FeaturePlot') > -1",
                    plotOutput(ns("FeaturePlot")) %>% 
                      shinycssloaders::withSpinner(),
                    ns = ns),
                  conditionalPanel(
                    condition = "input.plots.indexOf('Violin') > -1",
                    plotOutput(ns("Violin")) %>% 
                      shinycssloaders::withSpinner(),
                    ns = ns),
                  conditionalPanel(
                    condition = "input.plots.indexOf('Correlation_plot') > -1",
                    plotOutput(ns("Correlation_plot")) %>% 
                      shinycssloaders::withSpinner(),
                    ns = ns),
                  conditionalPanel(
                    condition = "input.plots.indexOf('CoveragePlot') > -1",
                    plotOutput(ns("CoveragePlot"), height = "1000px") %>% 
                      shinycssloaders::withSpinner(),
                    ns = ns),
                  conditionalPanel(
                    condition = "input.plots.indexOf('CoveragePlot_genome') > -1",
                    plotOutput(ns("CoveragePlot_genome"), height = "1000px") %>% 
                      shinycssloaders::withSpinner(),
                    ns = ns))
    )
}


#-------------------------------------------SERVER-----------------------------------------------------
# not all datasets have EES in metadata, thus indicated by argument EES_absent

sh_layout <- function(input, output, session, dataset, UMAP_label, EES_absent = "no", assay = "RNA") {
    
    output$UMAP <- renderPlot({
      DimPlot(object = dataset, reduction = "umap", label = TRUE,
              label.size = 5) + NoLegend() + ggtitle(label = UMAP_label)
    })
    
    output$EES_FeaturePlot <- renderPlot({
      EES_FeaturePlot(Seurat_object = dataset, split_type = input$group)
    })
    
    #----------------------------------cluster selection outside of ATAC browser-----------------------------------------
    
    # #NOTE ignoreNULL = F to ensure that when there is no selection at launch, the user only needs to click on Update gene
    update_cluster <- eventReactive(input$cluster_selection, {input$cluster},
                                    ignoreNULL = FALSE)

    # option to clear or select all checkboxes from clusters

    observeEvent(input$reset_clusters, {
      updateCheckboxGroupInput(session, "cluster",
                               choices = sort(as.character(unique(dataset@meta.data$CellType))),
                               selected = NULL)
    })

    observeEvent(input$select_all_clusters, {
      updateCheckboxGroupInput(session, "cluster",
                               choices = sort(as.character(unique(dataset@meta.data$CellType))),
                               selected = sort(as.character(unique(dataset@meta.data$CellType))))
    })
    #-----------------------------------------------------------------------------------------------------------------------

    output$EES_Violin <- renderPlot({
      VlnPlot_single_dataset(Seurat_object = dataset, split_type = input$group, features = "EES", idents = update_cluster(),
                             assay = assay)
    })

    update_gene <- eventReactive(input$go, {rat_nomenclature(input$gene,dataset,assay)})

    output$FeaturePlot <- renderPlot({
      
      #change assay as needed for input:

      if (input$group == "All") {
        FeaturePlot(object = change_assay(dataset = dataset, assay = assay),
                    cols = c("lightgrey", "#0072B2"),
                    features = update_gene(),
                    max.cutoff = input$expression,
                    split.by = NULL)
        } else {
          FeaturePlot(object = change_assay(dataset = dataset, assay = assay),
                      cols = c("lightgrey", "#0072B2"),
                      features = update_gene(),
                      max.cutoff = input$expression,
                      split.by = input$group)
        }
    })
    
    observeEvent(input$reset, {shinyjs::reset("expression")})

    output$Violin <- renderPlot({
      VlnPlot_single_dataset(Seurat_object = dataset, split_type = input$group, features = update_gene(), idents = update_cluster(),
                             assay = assay)
    })
    
    update_feature_corr <- eventReactive(input$go_corr, {feature2_eval(input$feature_corr, dataset, EES_absent, assay = assay)})
    
    output$Correlation_plot <- renderPlot({
      Scatter_feature(Seurat_object = dataset, split_type = input$group, 
                      features = update_gene(), features2 = update_feature_corr(), idents = input$cluster_corr,
                      assay = assay)
    })
    
    #----------------------------------cluster selection __for__ ATAC browser-----------------------------------------
    
    #NOTE ignoreNULL = F to ensure that when there is no selection at launch, the user only needs to click on Update gene
    update_cluster_browser <- eventReactive(input$cluster_selection_browser, {input$cluster_browser},
                                            ignoreNULL = FALSE)
    
    # option to clear or select all checkboxes in ATAC browser
    
    observeEvent(input$reset_browser_clusters, {
      updateCheckboxGroupInput(session, "cluster_browser", 
                               choices = sort(as.character(unique(dataset@meta.data$CellType))), 
                               selected = NULL)
    })
    
    observeEvent(input$select_all_browser_clusters, {
      updateCheckboxGroupInput(session, "cluster_browser", 
                               choices = sort(as.character(unique(dataset@meta.data$CellType))), 
                               selected = sort(as.character(unique(dataset@meta.data$CellType))))
    })
    
    #-----------------------------------------------------------------------------------------------------------------------

    output$CoveragePlot <- renderPlot({
      CoveragePlot_ratlas(Seurat_object = dataset, split_type = input$group,
                          region = update_gene(), extend_window = as.numeric(input$extend_window),
                          show_bulk = input$show_bulk, idents = update_cluster_browser(),
                          assay = assay)
    })

    update_cord <- eventReactive(input$go_cords, {validate_cords(input$genome_coords,dataset, assay = assay)})
    
    # optionally provide a second plot for genomic coordinates
    output$CoveragePlot_genome <- renderPlot({
      CoveragePlot_ratlas(Seurat_object = dataset, split_type = input$group,
                          region = update_cord(), extend_window = 0,
                          show_bulk = input$show_bulk, idents = update_cluster_browser(),
                          assay = assay)
    })
    
}