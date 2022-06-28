# Helper functions

change_assay <- function(dataset, assay = assay) {
  DefaultAssay(dataset) <- assay
  
  return(dataset)
}

rat_nomenclature <- function(user_gene, dataset, assay = assay) {
  # first ensure the correct assay is selected under each mode
  
  dataset <- change_assay(dataset = dataset, assay = assay)
  
  # remove any empty spaces if present (since gene names don't have that)
  user_gene <- gsub(" ","",user_gene, fixed = TRUE)
  
  # ensure that the input is not empty
  validate(
    need(nchar(user_gene) > 0,
         message = "You did not type a gene name. Please choose a gene (e.g.: Drd2)")
  )
  
  # if user types gene name inside quotes, app will reads as is (no nomenclature conversion)
  #NOTE: for multimodal will likely need to modify this further...
  
  # I use '"' but can also use "\"", with quote escaped ...
  if ((substr(user_gene,1,1) == '"' & substr(user_gene,nchar(user_gene),nchar(user_gene)) == '"')) {
    
    # if user used quotes, we then now need to remove it to search it within the object as is
    
    user_gene <- gsub("\"","",user_gene, fixed = TRUE)
    
    if (TRUE %in% (c("RNA", "RNArn7") %in% Assays(dataset))) {
      
      if (DefaultAssay(dataset) == "RNA") {
        
        validate(
          need(user_gene %in% rownames(dataset@assays$RNA@data),
               message = paste0("The gene name ", user_gene, " was not found in the dataset"))
        )
        
        return(user_gene)
        
      } else if (DefaultAssay(dataset) == "RNArn7") {

        validate(
          need(user_gene %in% rownames(dataset@assays$RNArn7@data),
               message = paste0("The gene name ", user_gene, " was not found in the dataset"))
        )
        
        return(user_gene)
        
      }
      
    } else if (TRUE %in% (c("ATAC", "peaks") %in% Assays(dataset))) {
      validate(
        need(user_gene %in% Annotation(dataset)$gene_name,
             message = paste0("The gene name ", user_gene, " was not found in the dataset"))
      )
      
      return(user_gene)
    }
    
  } else {
    
    # convert any input gene to rat nomeclature
    user_gene <- paste0(toupper(substr(user_gene, 1,1)),
                        tolower(substr(user_gene,2,nchar(user_gene))))
    
    if ("RNA" %in% Assays(dataset)) {
      
      validate(
        need(user_gene %in% rownames(dataset@assays$RNA@data),
             message = paste0("The gene name ", user_gene, " was not found in the dataset"))
      )
      
      return(user_gene)
      
    } else if (TRUE %in% (c("ATAC", "peaks") %in% Assays(dataset))) {
      validate(
        need(user_gene %in% Annotation(dataset)$gene_name,
             message = paste0("The gene name ", user_gene, " was not found in the dataset"))
      )
      
      return(user_gene)
    }
  }
}

# validate coordinates format for datasets
# NOTE: check cannot be against ranges in object to provide flexibility to plot any genomic coordinate
# not just peaks. May need to consider a check for valid coordinates against rn6 genome.

validate_cords <- function(user_cords, dataset, assay = assay) {
  
  #change assay as needed:
  dataset <- change_assay(dataset = dataset, assay = assay)
  
  # remove any empty spaces if present
  user_cords <- gsub(" ","",user_cords, fixed = TRUE)
  
  # remove commas and : if present and provide warning to users
  
  if (length(grep(",",user_cords)) > 0) {
    showNotification("FYI: Coordinates should not contain commas, removing them now. Format should follow: chromossome_number-start-end",
                     type = "warning",
                     duration = 8)
    user_cords <- gsub(",","",user_cords, fixed = TRUE)
  }


  if (length(grep(":",user_cords)) > 0) {
    showNotification("FYI: Changing : to -. Format should follow: chromossome_number-start-end",
                     type = "warning",
                     duration = 8)
    user_cords <- gsub(":","-",user_cords, fixed = TRUE)
  }

  # remove chr if user adds it to coordinate (datasets use Ensembl annotation)
  
  if (length(grep("chr", ignore.case = TRUE,user_cords)) > 0) {
    user_cords <- gsub("chr","", ignore.case = TRUE, user_cords)
  }
  
  # ensure that the input is not empty
  validate(
    need(nchar(user_cords) > 0,
         message = "You did not provide a genomic coordinate (e.g.: 18-84984575-84985121)")
  )

  return(user_cords)
}

# ensure EES is capitalized or if gene, use rat nomenclature.
feature2_eval <- function(user_feature2, dataset, EES_absent = "no", assay = assay) {
  
  if (EES_absent == "yes") {
    err_empty_message <- "You did not type a second gene name for correlation. Note: EES not available for current selected dataset"
  } else {
    err_empty_message <- "You did not type a second gene name or EES for correlation. Please type a gene name or EES for correlation"
  }
  
  # ensure that the input is not empty
  validate(
    need(nchar(user_feature2) > 0,
         message = err_empty_message)
  )
  
  if (toupper(user_feature2) == "EES") {
    
    # EES not present in Culture/VTA, thus display clear message if user tries to plot it
    validate(
      need(EES_absent == "no",
           message = "EES data only available for Adult NAc dataset. Please choose gene names only in current dataset")
    )
    
    return(toupper(user_feature2))
    
  } else {
    rat_nomenclature(user_feature2, dataset, assay)
  }
}

# EES plots

EES_FeaturePlot <- function(Seurat_object, split_type = "All") {
  if (split_type == "All") {
    FeaturePlot(object = Seurat_object,
                features = "EES",
                split.by = NULL)  + scale_color_gradientn(colours = c("gray90", "gray80", "red"))
  } else {
    EES_ggplot_list <- FeaturePlot(object = Seurat_object,
                                   features = "EES",
                                   split.by = split_type,
                                   combine = FALSE)
    
    # applies themes as Seurat removes x-axis labels from first plot when combine=F;
    # also applies limits to scale color for each plot as min/max.cutoff fails when combine=F
    EES_UMAP_l <- list()
    for (n in 1:length(EES_ggplot_list)) {
      EES_UMAP_l[[n]] <- EES_ggplot_list[[n]] +
        scale_color_gradientn(colours = c("gray90", "gray80", "red"),
                              limits=c(round(min(Seurat_object@meta.data$EES), digits = 1),
                                       round(max(Seurat_object@meta.data$EES), digits = 1))) +
        theme(axis.text.x = element_text(size=13), axis.ticks.x = element_line(),
              axis.line.x = element_line(), panel.border = element_blank(),
              axis.title.y.right = element_text(size = 16),
              plot.caption = element_text(size = 16, hjust = 0.5)) +
        labs(title = "", x = NULL, caption = "UMAP_1")
    }
    
    plot_grid(plotlist = EES_UMAP_l, align = TRUE, labels = "EES")
    
  }
}

# Violin plot for split.by (to prevent needing to subset data); features should be "EES" or gene from update_gene()

VlnPlot_single_dataset <- function(Seurat_object, split_type = "All", features = features, idents = idents, assay = assay) {
  
  #change assay as needed:
  Seurat_object <- change_assay(dataset = Seurat_object, assay = assay)
  
  if (split_type == "All") {
    VlnPlot(object = Seurat_object,
            features = features,
            split.by = NULL, pt.size = 0, idents = idents, log = FALSE) + NoLegend()
  } else {
    # here is the modified secion for when split_type is used
    # grab the relevant data from the VlnPlot and feed it into standard geom_violin
    vln_data <- VlnPlot(object = Seurat_object,
                        features = features,
                        split.by = split_type,
                        idents = idents,
                        log = FALSE, adjust = 1)
    
    vln_data <- vln_data$data
    names(vln_data)[1] <- "gene_of_interest" # renaming first col
    
    ggplot(vln_data,
           aes(x=ident,
               y=gene_of_interest,
               fill=split)) + 
      geom_violin(scale = "width",
                  adjust = 1) +
      ggtitle(features) +
      xlab("Identity") +
      ylab("Expression Level") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size=14),
            axis.text.y = element_text(size=13),
            axis.title =  element_text(size=16),
            axis.line = element_line(),
            plot.title = element_text(size = 18, face = "bold"),
            panel.background = element_blank(),
            legend.text = element_text(size=14), legend.title = element_text(size=15))
  }
}

# Scatter plot of correlation between features

Scatter_feature <- function(Seurat_object, split_type = "All", features = features, features2 = features2, idents = idents, assay = assay) {
  
  #change assay as needed:
  Seurat_object <- change_assay(dataset = Seurat_object, assay = assay)
  
  if (split_type == "All") {
    fetched_scatter_data <- base::subset(FetchData(Seurat_object,
                                                   vars = c("CellType",features,features2),
                                                   slot = "data"),
                                         CellType == idents)
    aes_mapping_colour <- NULL
    
  } else {
    fetched_scatter_data <- base::subset(FetchData(Seurat_object,
                                                   vars = c("CellType",features,features2, split_type),
                                                   slot = "data"),
                                         CellType == idents)
    
    aes_mapping_colour <- split_type
  }
  
  # due to gene names that may start with numbers (issue #2), for simplicity rename cols
  
  names(fetched_scatter_data)[2] <- "gene_name1" #x-axis
  names(fetched_scatter_data)[3] <- "gene_name2" #y-axis
  
  ggplot(fetched_scatter_data,
         aes(x=gene_name1,
             y=gene_name2)) +
    geom_point(size = 1, alpha = 0.7,  aes_string(colour=aes_mapping_colour), stroke = 1) +
    theme(axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),
          axis.title =  element_text(size=16),
          axis.line = element_line(),
          plot.title = element_text(size = 18, face = "bold"),
          panel.background = element_blank(),
          legend.text = element_text(size=14),
          legend.title = element_text(size=15),
          legend.key = element_blank()) +
    xlab(features) +
    ylab(features2) +
    ggtitle(
      paste0("Celltype: ", idents, " \nPearson Correlation: ",
             round(cor(fetched_scatter_data[,"gene_name1"],
                       fetched_scatter_data[,"gene_name2"],
                       method = "pearson"), digits = 2))
    )
}
