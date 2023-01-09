#' Scatter_feature
#' Generates a scatter plot of correlation between two features from a Seurat object
#' @param Seurat_object A Seurat object
#' @param split_type Factor to split the groups by. Default is show all / do not split.
#' @param features First feature from 2 features to draw correlation from
#' @param features2 Second feature from 2 features to draw correlation from
#' @param idents Which cell type / ident to draw the correlation from
#' @param assay Assay within the Seurat object to search gene/features
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' Scatter_feature(Seurat_object = input_obj, split_type = "All", 
#' features = "Gad1", features2 = "Gad2", idents = "Drd1-MSN", assay = "RNA")
#' }
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