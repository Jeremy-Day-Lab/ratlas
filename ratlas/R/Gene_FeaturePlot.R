#' Gene_FeaturePlot
#' A `FeaturePlot` wrapper which arranges plots with more than 4 plots as separate rows
#' with a maximun of 4 plot per row. Returns a single image which can be fit in most
#' user-screens from the app.
#' @param Seurat_object A Seurat object
#' @param feature Feature to plot
#' @param split_type Factor to split the groups by. Default is show all / do not split.
#' @param assay Assay within the Seurat object to search gene/features. Defaults to "RNA"
#' @param max.cutoff Maximum cutoff values for feature
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' Gene_FeaturePlot(Seurat_object = input_obj, split_type = "Stim_Sex", feature = "Gad1", assay = "RNA")
#' }
Gene_FeaturePlot <- function(Seurat_object, feature, split_type = "All", assay = "RNA", max.cutoff = NA) {
  
  #change assay as needed:
  Seurat_object <- change_assay(dataset = Seurat_object, assay = assay)
  
  # vector of length of groups and choice of split
  if (split_type != "All") {
    length_of_split <- length(unique(Seurat_object@meta.data[[split_type]]))
    split_group <- split_type
  } else {
    length_of_split <- 0 # "All" length is irrelevant
    split_group <- NULL
  }
  
  # any data that has more than 4 group in selected variable gets plotted
  # as a separate # for easier visibility within the app (mainly width)
  # If less than 4 or All, default to standard FeaturePlot
  if (split_type == "All" | length_of_split <= 4) {
    FeaturePlot(object = Seurat_object,
                cols = c("#D3D3D3", "#0072B2"),
                features = feature,
                max.cutoff = max.cutoff,
                split.by = split_group)
  } else if (length_of_split > 4) {
    
    ggplot_list <- FeaturePlot(object = Seurat_object,
                               features = feature,
                               split.by = split_group,
                               combine = FALSE)
    
    # applies themes as Seurat removes x-axis labels from first plot when combine=F;
    # also applies limits to scale color for each plot as min/max.cutoff fails when combine=F
    # oob = out-of-bounds values
    UMAP_l <- list()
    suppressMessages(
      for (n in 1:length(ggplot_list)) {
        UMAP_l[[n]] <- ggplot_list[[n]] +
          scale_color_gradientn(colours = c("#D3D3D3", "#0072B2"),
                                limits=c(0,max.cutoff),
                                oob = scales::squish) +
          theme(axis.text.x = element_text(size=13), axis.ticks.x = element_line(),
                axis.line.x = element_line(), panel.border = element_blank(),
                axis.title.y.right = element_text(size = 16),
                plot.caption = element_text(size = 16, hjust = 0.5)) +
          labs(title = "", x = NULL, caption = "UMAP_1")
      }
    )

    ggpubr::ggarrange(plotlist = UMAP_l, labels = feature, ncol = 4, nrow = 2)
  }
}
