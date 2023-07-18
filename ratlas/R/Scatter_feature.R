#' Scatter_feature
#' Generates a scatter plot of correlation between two features from a Seurat object
#' @param Seurat_object A Seurat object
#' @param split_type Factor to split the groups by. Default is show all / do not split.
#' @param feature First feature from 2 features to draw correlation from
#' @param feature2 Second feature from 2 features to draw correlation from
#' @param idents Which cell type / ident to draw the correlation from
#' @param assay Assay within the Seurat object to search gene/features. Defaults to "RNA"
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' Scatter_feature(Seurat_object = input_obj, split_type = "All", 
#' feature = "Gad1", feature2 = "Gad2", idents = "Drd1-MSN", assay = "RNA")
#' }
Scatter_feature <- function(Seurat_object, split_type = "All", feature, feature2, idents, assay = "RNA") {
  
  #change assay as needed:
  Seurat_object <- change_assay(dataset = Seurat_object, assay = assay)
  
  if (split_type == "All") {
    fetched_scatter_data <- base::subset(FetchData(Seurat_object,
                                                   vars = c("CellType",feature,feature2),
                                                   layer = "data"), 
                                         CellType == idents)
    aes_mapping_colour <- NULL
    
  } else {
    fetched_scatter_data <- base::subset(FetchData(Seurat_object,
                                                   vars = c("CellType",feature,feature2, split_type),
                                                   layer = "data"),
                                         CellType == idents)
    
    aes_mapping_colour <- split_type
  }
  
  # scatter plot with correlation value
  scatter_plot <- ggplot(fetched_scatter_data,
                         aes(x=.data[[feature]],
                             y=.data[[feature2]])) +
    geom_point(size = 1, alpha = 0.7, stroke = 1) +
    theme(axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),
          axis.title =  element_text(size=16),
          axis.line = element_line(),
          plot.title = element_text(size = 18, face = "bold"),
          panel.background = element_blank(),
          legend.text = element_text(size=14),
          legend.title = element_text(size=15),
          legend.key = element_blank()) +
    ggtitle(
      paste0("Celltype: ", idents, " \nPearson Correlation: ",
             round(cor(fetched_scatter_data[,feature],
                       fetched_scatter_data[,feature2],
                       method = "pearson"), digits = 2))
    )
  
  if (!is.null(aes_mapping_colour)) scatter_plot <- scatter_plot + aes(colour = .data[[aes_mapping_colour]])
  
  return(scatter_plot)
}