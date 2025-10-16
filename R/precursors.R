#' Filter precursor data
#'
#' @param diann_data DiaNNData R6 object (not DiaNNData)
#' @param proteotypic logical. Filter for proteotypic peptides only (default: TRUE)
#' @param remove_contaminants logical. Remove contaminant proteins (default: FALSE)
#'
#' @return data.frame. Filtered precursor data
#' @export
#' @importFrom dplyr filter left_join
#' @importFrom tidyr replace_na
precursors <- function(diann_data, proteotypic = TRUE, remove_contaminants = FALSE) {
  # Validate input
  if (!inherits(diann_data, "DiaNNData")) {
    stop("diann_data must be a DiaNNData R6 object", call. = FALSE)
  }
  
  result <- diann_data$input_data
  
  # Filter for proteotypic peptides
  if (proteotypic && "Proteotypic" %in% names(result)) {
    result <- dplyr::filter(result, Proteotypic == 1)
  }
  
  # Remove contaminants if requested
  if (remove_contaminants && length(diann_data$contaminants) > 0) {
    # Check if contaminants are already annotated in the data
    if ("is_contaminant" %in% names(result)) {
      # Use existing annotation
      result <- dplyr::filter(result, !is_contaminant)
    } else {
      # Filter directly using Protein.Ids
      result <- dplyr::filter(result, !(Protein.Ids %in% diann_data$contaminants))
    }
  }
  
  # Join with metadata if available
  if (!is.null(diann_data$metadata) && nrow(diann_data$metadata) > 0) {
    result <- dplyr::left_join(result, diann_data$metadata, by = c("File.Name", "Run"))
  }
  
  return(result)
}