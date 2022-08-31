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

# footer of saved figs
caption_label <- "Source: Ratlas | Day Lab"

##### plots available

all_plots <- c("UMAP","FeaturePlot","Violin", "EES_FeaturePlot","EES_Violin", "Correlation_plot")
subset_plots <- c("UMAP","FeaturePlot","Violin","Correlation_plot")

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
                                    selected = c("UMAP","FeaturePlot","Violin"),
                                    inline = TRUE),
                 
                 sliderInput(inputId = ns("width"),
                             label = "Width (pixels) for downloading plot",
                             value = 800,
                             min = 400,
                             max = 2000,
                             step = 100,
                             round = TRUE),
                 
                 sliderInput(inputId = ns("height"),
                             label = "Height (pixels) for downloading plot",
                             value = 800,
                             min = 400,
                             max = 2000,
                             step = 100,
                             round = TRUE),
                 
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
                downloadButton(ns("EES_FeaturePlot_downl"), label = "Download EES FeaturePlot"),
                ns = ns),
              conditionalPanel(
                condition = "input.plots.indexOf('EES_Violin') > -1",
                plotOutput(ns("EES_Violin")) %>% 
                  shinycssloaders::withSpinner(),
                downloadButton(ns("EES_Violin_downl"), label = "Download EES Violin"),
                ns = ns),
              conditionalPanel(
                condition = "input.plots.indexOf('FeaturePlot') > -1",
                plotOutput(ns("FeaturePlot")) %>% 
                  shinycssloaders::withSpinner(),
                ns = ns),
              conditionalPanel(
                condition = "output.FeaturePlot",
                downloadButton(ns("FeaturePlot_downl"), label = "Download FeaturePlot"),
                ns = ns), # downloadButton have their own conditional dependent upon output (for genes, not EES)
              conditionalPanel(
                condition = "input.plots.indexOf('Violin') > -1",
                plotOutput(ns("Violin")) %>% 
                  shinycssloaders::withSpinner(),
                ns = ns),
              conditionalPanel(
                condition = "output.Violin",
                downloadButton(ns("Violin_downl"), label = "Download Violin"),
                ns = ns),
              conditionalPanel(
                condition = "input.plots.indexOf('Correlation_plot') > -1",
                plotOutput(ns("Correlation_plot")) %>% 
                  shinycssloaders::withSpinner(),
                ns = ns),
              conditionalPanel(
                condition = "output.Correlation_plot",
                downloadButton(ns("Correlation_downl"), label = "Download correlation plot"),
                ns = ns)
    )
  )
}


#-------------------------------------------SERVER-----------------------------------------------------
# not all datasets have EES in metadata, thus indicated by argument EES_absent

sh_layout <- function(input, output, session, dataset, UMAP_label, EES_absent = "no", assay = "RNA") {
  
  output$UMAP <- renderPlot({
    DimPlot(object = dataset, reduction = "umap", label = TRUE,
            label.size = 5) + NoLegend() + ggtitle(label = UMAP_label)
  })
  
  # make plots in reactive code repetition when saving figs
  ees_featureplot <- reactive({
    EES_FeaturePlot(Seurat_object = dataset, split_type = input$group)
  })
  
  output$EES_FeaturePlot <- renderPlot({
    ees_featureplot()
  })
  
  output$EES_FeaturePlot_downl <- downloadHandler(
    filename = function() { "EES_FeaturePlot_Ratlas.png" },
    
    content = function(file) {
      png(file, height = input$height, width = input$width)
      
      # adding footer to plot_grid requires a different process:
      
      if (input$group == "All") {
        plot(ees_featureplot() + labs(caption = caption_label) + theme(plot.caption = element_text(size=18, face="bold")))
      } else {
        plot(ees_featureplot() + draw_label(caption_label, x = 1, y = 0, hjust = 1.5, vjust = 1, size = 18, fontface = "bold"))
      }
      
      dev.off()
    }
  )
  #----------------------------------cluster selection outside for non-ATAC visualization -----------------------------------------
  
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
  
  ees_violin <- reactive({
    VlnPlot_single_dataset(Seurat_object = dataset, split_type = input$group, features = "EES", idents = update_cluster(),
                           assay = assay)
  })
  
  output$EES_Violin <- renderPlot({
    ees_violin()
  })
  
  output$EES_Violin_downl <- downloadHandler(
    filename = function() { "EES_Violin_Ratlas.png" },
    
    content = function(file) {
      png(file, height = input$height, width = input$width)
      
      plot(ees_violin() + labs(caption = caption_label) + theme(plot.caption = element_text(size=18, face="bold")))
      
      dev.off()
    }
  )
  
  update_gene <- eventReactive(input$go, {rat_nomenclature(input$gene,dataset,assay)})
  
  featureplot <- reactive({
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
  
  output$FeaturePlot <- renderPlot({
    featureplot()
  })
  
  output$FeaturePlot_downl <- downloadHandler(
    filename = function() { paste0(rat_nomenclature(input$gene,dataset,assay),"_FeaturePlot_Ratlas.png") },
    
    content = function(file) {
      png(file, height = input$height, width = input$width)
      
      plot(featureplot() + labs(caption = caption_label) + theme(plot.caption = element_text(size=18, face="bold")))
      
      dev.off()
    }
  )
  
  observeEvent(input$reset, {shinyjs::reset("expression")})
  
  violin <- reactive({
    VlnPlot_single_dataset(Seurat_object = dataset, split_type = input$group, features = update_gene(), idents = update_cluster(),
                           assay = assay)
  })
  
  
  output$Violin <- renderPlot({
    violin()
  })
  
  
  output$Violin_downl <- downloadHandler(
    filename = function() { paste0(rat_nomenclature(input$gene,dataset,assay),"_Violin_Ratlas.png") },
    
    content = function(file) {
      png(file, height = input$height, width = input$width)
      
      plot(violin() + labs(caption = caption_label) + theme(plot.caption = element_text(size=18, face="bold")))
      
      dev.off()
    }
  )
  
  
  update_feature_corr <- eventReactive(input$go_corr, {feature2_eval(input$feature_corr, dataset, EES_absent, assay = assay)})
  
  corr_plot <- reactive({
    Scatter_feature(Seurat_object = dataset, split_type = input$group, 
                    features = update_gene(), features2 = update_feature_corr(), idents = input$cluster_corr,
                    assay = assay)
  })
  
  output$Correlation_plot <- renderPlot({
    corr_plot()
  })
  
  output$Correlation_downl <- downloadHandler(
    filename = function() { paste0(rat_nomenclature(input$gene,dataset,assay), "_",
                                   feature2_eval(input$feature_corr, dataset, EES_absent, assay = assay),
                                   "_Correlation_Ratlas.png") },
    
    content = function(file) {
      png(file, height = input$height, width = input$width)
      
      plot(corr_plot() + labs(caption = caption_label) + theme(plot.caption = element_text(size=18, face="bold")))
      
      dev.off()
    }
  )
}