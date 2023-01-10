# app modules
#-----------------------------------------global UI options------------------------------------------------------------------
options(spinner.type=1,spinner.color="#232a30", spinner.size=2)

#-----------------------------------------other UI variables-----------------------------------------------------------------
##### group choices
# NOTE: may to global vars
# adult group choices
adult_acute_groups <- c("All", "Stim","Sex","Stim_Sex")

# adult acute+repeated groups
adult_acute_repeated_groups <- c("All", "Stim","Sex", "Dataset","Stim_Sex", "Dataset_Stim", "Dataset_Sex", "Dataset_Stim_Sex")

# choices for culture dataset
all_stim_groups <- c("All", "Stim")

# choices for adult VTA
all_VTA_groups <- c("All", "Sex")

# correlation plot UI message based EES data being present
contains_EES <- "Type a gene or EES to correlate to gene name typed above"
no_EES <- "Type a gene to correlate to gene name selected above"

# footer of saved figs
caption_label <- "Source: Ratlas | Day Lab"

##### plots available

all_plots_EES <- c("UMAP","FeaturePlot","Violin", "EES_FeaturePlot","EES_Violin", "Correlation_plot")
all_plots <- c("UMAP","FeaturePlot","Violin","Correlation_plot")

#----------------------------------------------------UI-----------------------------------------------------------------------
sh_layout_UI <- function(id, group_choices, plot_choices, cluster_names, correlation_label) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(width = 3,
                 checkboxGroupInput(inputId = ns("plots"),
                                    label = "Choose which plot(s) to display",
                                    choices = plot_choices,
                                    selected = c("UMAP","FeaturePlot","Violin"),
                                    inline = TRUE),
                 
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
                              icon("dna")
                              ),
                 
                 selectInput(inputId = ns("group"),
                             label = "Choose all or split by group",
                             choices = group_choices,
                             selected = "All",
                             multiple = FALSE),

                 
                 radioButtons(inputId = ns("plot_type"),
                                    label = "Output type for downloading plot",
                                    choices = c("png","pdf"),
                                    selected = "png",
                                    inline = TRUE),
                 
                 hr(),
                 
                 # conditional panels based on user choices
                 
                 conditionalPanel(
                   condition = "input.plot_type.indexOf('png') > -1",
                   sliderInput(inputId = ns("png_width"),
                               label = "Width (pixels) for downloading PNG plot",
                               value = 800,
                               min = 400,
                               max = 3000,
                               step = 100,
                               round = TRUE),
                   
                   sliderInput(inputId = ns("png_height"),
                               label = "Height (pixels) for downloading PNG plot",
                               value = 800,
                               min = 400,
                               max = 3000,
                               step = 100,
                               round = TRUE),
                   hr(),
                   ns = ns),
                 
                 conditionalPanel(
                   condition = "input.plot_type.indexOf('pdf') > -1",
                   sliderInput(inputId = ns("pdf_width"),
                               label = "Width (inches) for downloading PDF plot",
                               value = 8,
                               min = 4,
                               max = 30,
                               step = 1,
                               round = TRUE),
                   
                   sliderInput(inputId = ns("pdf_height"),
                               label = "Height (inches) for downloading PDF plot",
                               value = 8,
                               min = 4,
                               max = 30,
                               step = 1,
                               round = TRUE),
                   hr(),
                   ns = ns),
                 
                 conditionalPanel(
                   condition = "input.plots.indexOf('FeaturePlot') > -1",
                   numericInput(inputId = ns("expression"),
                                label = HTML("FeaturePlot option: <br/> set a max scale cutoff"),
                                value = NA,
                                min = 0) %>%
                     shinyhelper::helper(icon = "exclamation-circle",
                                         colour ="#232a30",
                                         type = "markdown",
                                         content = "Featureplot_cutoff_help"
                     ),
                   
                   actionButton(inputId = ns("reset"),
                                label = "Reset scale",
                                icon("redo")
                                ),
                   hr(),
                   ns = ns),
                 
                 conditionalPanel(
                   condition = "input.plots.indexOf('Violin') > -1 || input.plots.indexOf('EES_Violin') > -1",
                   
                   checkboxInput(inputId = ns("pt_size"), 
                                 label = "Violin plot option: display single-nuclei as points", 
                                 value = FALSE),
                   
                   checkboxGroupInput(inputId = ns("cluster"),
                                      label = "Violin plot option: choose all or show specific clusters",
                                      choices = cluster_names,
                                      selected = cluster_names),
                   
                   actionButton(inputId = ns("cluster_selection"),
                                label = "Plot selected clusters",
                                icon("check")
                                ) %>%
                     shinyhelper::helper(icon = "info-circle",
                                         colour ="#232a30",
                                         type = "markdown",
                                         content = "Plot_selected_clusters"
                     ),
                   
                   actionButton(inputId = ns("reset_clusters"),
                                label = "Clear all cluster choices",
                                icon("redo")
                                ),
                   
                   actionButton(inputId = ns("select_all_clusters"),
                                label = "Check all cluster choices",
                                icon("check-double")
                                ),
                   
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
                                icon("chart-line")
                                ),
                   ns = ns),
    ),
    # downloadButton have their own conditional dependent upon output as well (for genes, not EES)
    # to ensure the button doesn't show up when no gene is selected or if an error is shown
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
                condition = "output.FeaturePlot && input.plots.indexOf('FeaturePlot') > -1",
                downloadButton(ns("FeaturePlot_downl"), label = "Download FeaturePlot"),
                ns = ns),
              conditionalPanel(
                condition = "input.plots.indexOf('Violin') > -1",
                plotOutput(ns("Violin")) %>% 
                  shinycssloaders::withSpinner(),
                ns = ns),
              conditionalPanel(
                condition = "output.Violin && input.plots.indexOf('Violin') > -1",
                downloadButton(ns("Violin_downl"), label = "Download Violin"),
                ns = ns),
              conditionalPanel(
                condition = "input.plots.indexOf('Correlation_plot') > -1",
                plotOutput(ns("Correlation_plot")) %>% 
                  shinycssloaders::withSpinner(),
                ns = ns),
              conditionalPanel(
                condition = "output.Correlation_plot && input.plots.indexOf('Correlation_plot') > -1",
                downloadButton(ns("Correlation_downl"), label = "Download correlation plot"),
                ns = ns)
    )
  )
}


#-------------------------------------------SERVER-----------------------------------------------------
# not all datasets have EES in metadata, thus indicated by argument EES_absent

sh_layout_server <- function(id, dataset, UMAP_label, EES_absent = FALSE, assay = "RNA") {
  
  moduleServer(
    id,
    function(input, output, session) {
      
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
        filename = function() { paste0("EES_FeaturePlot_Ratlas.", input$plot_type) },
        
        content = function(file) {
          
          height <- ifelse(input$plot_type == "png", input$png_height, input$pdf_height)
          width <- ifelse(input$plot_type == "png", input$png_width, input$pdf_width)
          
          # adding footer to plot_grid requires a different process:
          
          if (input$group == "All") {
            plot_save <- ees_featureplot() + labs(caption = caption_label) + theme(plot.caption = element_text(size=18, face="bold"))
          } else {
            plot_save <- ees_featureplot() + draw_label(caption_label, x = 1, y = 0, hjust = 1.5, vjust = 1, size = 18, fontface = "bold")
          }
          
          plot_png_pdf(file_name = file, plot = plot_save, height = height, width = width, image_format = input$plot_type)
          
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
        
        #evaluate violin options
        if (input$group != "All") {
          split_group <- input$group
        } else {
          split_group <- NULL
        }
        
        if (input$pt_size == FALSE) {
          size <- 0
        } else {
          size <- NULL # default pt.size
        }
        
        VlnPlot(object = dataset, split.by = split_group,
                features = "EES", idents = update_cluster(),
                pt.size = size, assay = assay, log = FALSE)
      })
      
      output$EES_Violin <- renderPlot({
        ees_violin()
      })
      
      output$EES_Violin_downl <- downloadHandler(
        filename = function() { paste0("EES_Violin_Ratlas.", input$plot_type) },
        
        content = function(file) {
          
          height <- ifelse(input$plot_type == "png", input$png_height, input$pdf_height)
          width <- ifelse(input$plot_type == "png", input$png_width, input$pdf_width)
          
          plot_save <- ees_violin() + labs(caption = caption_label) + theme(plot.caption = element_text(size=18, face="bold"))
          
          plot_png_pdf(file_name = file, plot = plot_save, height = height, width = width, image_format = input$plot_type)
          
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
        filename = function() { paste0(rat_nomenclature(input$gene,dataset,assay),"_FeaturePlot_Ratlas.", input$plot_type) },
        
        content = function(file) {
          
          height <- ifelse(input$plot_type == "png", input$png_height, input$pdf_height)
          width <- ifelse(input$plot_type == "png", input$png_width, input$pdf_width)
          
          plot_save <- featureplot() + labs(caption = caption_label) + theme(plot.caption = element_text(size=18, face="bold"))
          
          plot_png_pdf(file_name = file, plot = plot_save, height = height, width = width, image_format = input$plot_type)
          
        }
      )
      
      observeEvent(input$reset, {shinyjs::reset("expression")})
      
      violin <- reactive({
        
        #evaluate violin options
        if (input$group != "All") {
          split_group <- input$group
        } else {
          split_group <- NULL
        }
        
        if (input$pt_size == FALSE) {
          size <- 0
        } else {
          size <- NULL # default pt.size
        }
        
        VlnPlot(object = dataset, split.by = split_group,
                features = update_gene(), idents = update_cluster(),
                pt.size = size, assay = assay, log = FALSE)
      })
      
      
      output$Violin <- renderPlot({
        violin()
      })
      
      
      output$Violin_downl <- downloadHandler(
        filename = function() { paste0(rat_nomenclature(input$gene,dataset,assay),"_Violin_Ratlas.", input$plot_type) },
        
        content = function(file) {
          
          height <- ifelse(input$plot_type == "png", input$png_height, input$pdf_height)
          width <- ifelse(input$plot_type == "png", input$png_width, input$pdf_width)
          
          plot_save <- violin() + labs(caption = caption_label) + theme(plot.caption = element_text(size=18, face="bold"))
          
          plot_png_pdf(file_name = file, plot = plot_save, height = height, width = width, image_format = input$plot_type)
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
                                       "_Correlation_Ratlas.", input$plot_type) },
        
        content = function(file) {
          
          height <- ifelse(input$plot_type == "png", input$png_height, input$pdf_height)
          width <- ifelse(input$plot_type == "png", input$png_width, input$pdf_width)
          
          plot_save <- corr_plot() + labs(caption = caption_label) + theme(plot.caption = element_text(size=18, face="bold"))
          
          plot_png_pdf(file_name = file, plot = plot_save, height = height, width = width, image_format = input$plot_type)
          
        }
      )
    }
  )
}