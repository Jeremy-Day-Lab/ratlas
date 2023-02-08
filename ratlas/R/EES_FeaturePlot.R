#' EES_FeaturePlot
#' A wrapper function of `Seurat::FeaturePlot` which addresses needs when combine == FALSE: 
#' 1. applies consistent scale limits on split.by plots
#' 2. fixes an issue with the x-axis labels
#' Also applies custom colors to EES FeaturePlots
#' @param Seurat_object A Seurat object
#' @param split_type Factor to split the groups by. Default is show all / do not split.
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' EES_FeaturePlot(Seurat_object = input_obj, split_type = "Stim")
#' }
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
    suppressMessages(
      for (n in 1:length(EES_ggplot_list)) {
        EES_UMAP_l[[n]] <- EES_ggplot_list[[n]] +
          scale_color_gradientn(colours = c("gray90", "gray80", "red"),
                                limits=c(round(min(Seurat_object@meta.data$EES), digits = 1),
                                         round(max(Seurat_object@meta.data$EES), digits = 1)),
                                oob = scales::squish) +
          theme(axis.text.x = element_text(size=13), axis.ticks.x = element_line(),
                axis.line.x = element_line(), panel.border = element_blank(),
                axis.title.y.right = element_text(size = 16),
                plot.caption = element_text(size = 16, hjust = 0.5)) +
          labs(title = "", x = NULL, caption = "UMAP_1")
      }
    )
    
    ggpubr::ggarrange(plotlist = EES_UMAP_l, labels = "EES")
    
  }
}