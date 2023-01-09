#' rat_nomenclature
#' Checks that feature/gene input is valid for active dataset. By default, input is converted to
#' rat gene nomenclature, but users have the choice to overwrite default by
#' typing feature name inside double quotes to search for features which do not yet follow standard nomenclature. 
#' @param user_gene Gene/Feature name
#' @param dataset A Seurat object
#' @param assay Assay within the Seurat object to search gene/features
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' rat_nomenclature(user_gene = 'Gad1', dataset = input_obj, assay = 'RNA')
#' }
#' \dontrun{
#' rat_nomenclature(user_gene = '"ENSRNOG00000020272"', dataset = input_obj, assay = 'RNA') #does not convert to lowercase
#' }
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