#' Filter protein group data
#'
#' @param diann_data DiaNNData R6 object containing DIA-NN output data
#' @param proteotypic logical. Filter for proteotypic peptides only (default: TRUE)
#' @param min_peptides integer. Minimum number of peptides per protein (default: 2)
#' @param remove_contaminants logical. Remove contaminant proteins (default: FALSE)
#' @param filter_q_value numeric. Protein.Q.Value threshold for protein filtering (default: 0.01). Set to NULL to disable filtering.
#'
#' @return data.frame. Filtered protein group data
#' @export
#' @importFrom dplyr filter group_by summarise n_distinct left_join across all_of mutate
#' @importFrom tidyr replace_na
proteins <- function(diann_data, proteotypic = TRUE, min_peptides = 2, remove_contaminants = FALSE, filter_q_value = 0.01) {
  # Validate input
  if (!inherits(diann_data, "DiaNNData")) {
    stop("diann_data must be a DiaNNData R6 object", call. = FALSE)
  }
  
  result <- diann_data$input_data
  
  # Filter for proteotypic peptides first
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

  # Filter by Q-value if requested
  if (!is.null(filter_q_value) && "Protein.Q.Value" %in% names(result)) {
    result <- dplyr::filter(result, Protein.Q.Value <= filter_q_value)
  }
  
  # Check required columns exist
  required_cols <- c("File.Name", "Run", "Protein.Ids", "Stripped.Sequence")
  missing_cols <- setdiff(required_cols, names(result))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  
  # Group and summarize at protein level
  group_cols <- intersect(
    c("File.Name", "Run", "Protein.Ids", "Protein.Names", "PG.Quantity", "PG.Normalised", "PG.MaxLFQ"),
    names(result)
  )
  
  result <- result |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      peptides = dplyr::n_distinct(Stripped.Sequence),
      proteins = dplyr::n_distinct(Protein.Ids),
      .groups = "drop"
    ) |>
    dplyr::filter(peptides >= min_peptides)
  
  # Join with metadata if available
  if (!is.null(diann_data$metadata) && nrow(diann_data$metadata) > 0) {
    result <- dplyr::left_join(result, diann_data$metadata, by = c("File.Name", "Run"))
  }
  
  return(result)
}