#' convert_assay_v5
#' Converts objects to Seurat V5 with data layer being held in desk via BPCells
#'
#' @param dataset A Seurat object
#' @param assay Assay within the Seurat object to be converted to Assay5
#' @param tmp_output_path output path of on disk data - can be deleted as this will also be present in `seurat_output_path`
#' @param seurat_output_path output path of updated seurat object
#' @param updated_dataset_name name of output file for updated seurat object (include extension .rds in file name)
#'
#' @return an rds file containing updated Seurat object
#' @export
#'
#' @examples
convert_assay_v5 <- function(dataset, assay = "RNA", tmp_output_path, seurat_output_path, updated_dataset_name) {
  options(Seurat.object.assay.version = "v5")
  
  # Write the data layer to a directory
  # NOTE: main "slot"/layer is data (no counts for the app)
  # These will also be present in the final object dir, and thus can be deleted by the user
  write_matrix_dir(mat = dataset[[assay]]$data, dir = tmp_output_path)
  
  # load the matrix that was saved on disk
  counts_mat <- open_matrix_dir(dir = tmp_output_path)
  
  # convert the object's assay to Assay5
  dataset[[assay]] <- as(object = dataset[[assay]], Class = "Assay5")
  
  # replace data layer with on disk data
  dataset[[assay]]$data <- counts_mat
  
  # save updated object with `destdir` to ensure on-disk matrices are saved in the same location
  saveRDS(
    object = dataset,
    file = updated_dataset_name,
    destdir = seurat_output_path
  )
}