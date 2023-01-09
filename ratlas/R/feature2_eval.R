#' feature2_eval
#' Evaluates second feature to be used in correlation. Capitalizes ESS when needed or calls `rat_nomenclature`
#' @param user_feature2 Gene/Feature name
#' @param dataset A Seurat object
#' @param EES_absent indicates if dataset contains EES
#' @param assay Assay within the Seurat object to search gene/features
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' feature2_eval(user_feature2 = "Gad1", dataset = input_obj, EES_absent = FALSE, assay = "RNA")
#' }
feature2_eval <- function(user_feature2, dataset, EES_absent = FALSE, assay = assay) {
  
  if (EES_absent == TRUE) {
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
    
    # EES only present in adult acute NAc, thus display clear message if user tries to plot it
    validate(
      need(EES_absent == FALSE,
           message = "EES data only available for Adult acute NAc dataset. Please choose gene names only in current dataset")
    )
    
    return(toupper(user_feature2))
    
  } else {
    rat_nomenclature(user_feature2, dataset, assay)
  }
}