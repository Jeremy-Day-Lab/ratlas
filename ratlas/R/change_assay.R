#' change_assay
#' A helper function to swap assays when/if needed
#' @param dataset A Seurat object
#' @param assay Assay within the Seurat object to be made default. Defaults to "RNA"
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' change_assay(dataset = input_obj, assay = "RNA")
#' }
change_assay <- function(dataset, assay = "RNA") {
  DefaultAssay(dataset) <- assay
  
  return(dataset)
}