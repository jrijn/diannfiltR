#' DIA-NN Data Processor Class
#'
#' @description
#' An R6 class for processing and filtering DIA-NN data using method-based approach.
#' This class encapsulates filtering operations for precursors, proteins, and genes.
#' 
#' @details
#' The DiaNNProcessor class provides a structured approach to filter DIA-NN data
#' at multiple levels (precursors, proteins, genes) while handling contaminant
#' removal and quality control filtering.
#' 
#' Methods:
#' \describe{
#'   \item{\code{initialize(diann_data)}}{Create a new DiaNNProcessor object with a DiaNNData object}
#'   \item{\code{annotate_contaminants()}}{Annotate contaminants in the data}
#'   \item{\code{filter_precursors(proteotypic = TRUE, remove_contaminants = FALSE, digest_trypsin = FALSE)}}{
#'     Filter precursor-level data with options for proteotypic filtering, contaminant removal, and tryptic digest calculation}
#'   \item{\code{filter_proteins(proteotypic = TRUE, min_peptides = 2, remove_contaminants = FALSE, filter_q_value = 0.01)}}{
#'     Filter and summarize protein group data}
#'   \item{\code{filter_genes(proteotypic = TRUE, min_peptides = 2, remove_contaminants = FALSE, filter_q_value = 0.01)}}{
#'     Filter and summarize gene-level data}
#' }
#'
#' @field data data.frame. DIA-NN data to be processed
#' @field contaminants data.frame or character vector. Contaminant protein annotations
#' @field parameters list. Available data types detected in input
#' @field metadata data.frame. File and run metadata
#'
#' @examples
#' \dontrun{
#' # Create processor with DiaNNData object
#' diann_data <- DiaNNData$new(input_data = your_data)
#' processor <- DiaNNProcessor$new(diann_data)
#' 
#' # Filter precursors with tryptic digest
#' precursor_result <- processor$filter_precursors(
#'   proteotypic = TRUE, 
#'   remove_contaminants = TRUE,
#'   digest_trypsin = TRUE
#' )
#' 
#' # Filter proteins
#' protein_result <- processor$filter_proteins(
#'   min_peptides = 2,
#'   remove_contaminants = TRUE,
#'   filter_q_value = 0.01
#' )
#' 
#' # Filter genes
#' gene_result <- processor$filter_genes(
#'   min_peptides = 2, 
#'   remove_contaminants = TRUE
#' )
#' }
#'
#' @importFrom R6 R6Class
#' @importFrom dplyr filter group_by summarise n_distinct left_join across all_of mutate select
#' @importFrom tidyr replace_na
#' @export
DiaNNProcessor <- R6Class(
  classname = "DiaNNProcessor",
  
  public = list(
    # Public fields
    data = NULL,
    contaminants = NULL,
    parameters = NULL,
    metadata = NULL,
    
    #' @description
    #' Initialize DiaNNProcessor object
    #' @param diann_data DiaNNData R6 object containing DIA-NN output data
    #' @return A new DiaNNProcessor object
    initialize = function(diann_data) {
      if (!inherits(diann_data, "DiaNNData")) {
        stop("diann_data must be a DiaNNData R6 object", call. = FALSE)
      }
      
      # Store the input data for processing
      self$data <- diann_data$input_data
      
      # Store the contaminants for processing
      self$contaminants <- diann_data$contaminants
      
      # Store parameters and metadata
      self$parameters <- diann_data$parameters
      self$metadata <- diann_data$metadata
      
      # Annotate contaminants in the data
      self$annotate_contaminants()
    },

    #' Annotate contaminants in the data
    #' @return Invisible self
    annotate_contaminants = function() {
      if (is.null(self$contaminants) || length(self$contaminants) == 0) {
        warning("No contaminants data available to annotate", call. = FALSE)
        return(invisible(self))
      }
      
      # Annotate contaminants in the data
      if (is.data.frame(self$contaminants) && "Protein.Ids" %in% names(self$contaminants)) {
        self$data <- self$data |>
          dplyr::left_join(
            dplyr::select(self$contaminants, Protein.Ids, Source.of.Contamination),
            by = "Protein.Ids"
          ) |>
          dplyr::mutate(
            is_contaminant = !is.na(Source.of.Contamination)
          )
      }
      
      invisible(self)
    },

    #' @description
    #' Filter precursor-level data
    #' @param proteotypic logical. Filter for proteotypic peptides only (default: TRUE)
    #' @param remove_contaminants logical. Remove contaminant proteins (default: FALSE)
    #' @param digest_trypsin logical. Calculate missed cleavages for tryptic digest (default: FALSE)
    #' @return data.frame. Filtered precursor data
    filter_precursors = function(proteotypic = TRUE, remove_contaminants = FALSE, digest_trypsin = FALSE) {
      if (!self$parameters$precursors) {
        stop("Precursor data not available in input data", call. = FALSE)
      }
      
      result <- self$data
      
      # Filter for proteotypic peptides
      if (proteotypic && "Proteotypic" %in% names(result)) {
        result <- dplyr::filter(result, Proteotypic == 1)
      }
      
      # Remove contaminants if requested
      if (remove_contaminants && !is.null(self$contaminants)) {
        # Check if contaminants are already annotated in the data
        if ("is_contaminant" %in% names(result)) {
          # Use existing annotation
          result <- dplyr::filter(result, !is_contaminant)
        } else {
          # Filter directly using Protein.Ids from contaminants data.frame
          if (is.data.frame(self$contaminants) && "Protein.Ids" %in% names(self$contaminants)) {
            contaminant_ids <- self$contaminants$Protein.Ids
            result <- dplyr::filter(result, !(Protein.Ids %in% contaminant_ids))
          } else if (is.character(self$contaminants)) {
            # Handle character vector case
            result <- dplyr::filter(result, !(Protein.Ids %in% self$contaminants))
          }
        }
      }

      # Calculate missed cleavages if requested
      if (digest_trypsin) {
        if (!"Stripped.Sequence" %in% names(result)) {
          stop("Stripped.Sequence column not found in precursor data", call. = FALSE)
        }
        
        # Calculate missed cleavages
        result$missed_cleavages <- sapply(
          result$Stripped.Sequence, 
          function(x) {
            if (is.na(x) || x == "") return(0)
            # Count K and R not followed by P (tryptic cleavage sites)
            fragments <- strsplit(x, "(?!P)(?<=[RK])", perl = TRUE)[[1]]
            max(0, length(fragments) - 1)  # Ensure non-negative
          }
        )
      }
      
      # Join with metadata if available
      if (!is.null(self$metadata) && nrow(self$metadata) > 0) {
        result <- dplyr::left_join(result, self$metadata, by = c("File.Name", "Run"))
      }
      
      return(result)
    },

    #' @description
    #' Filter and summarize protein group data
    #' @param proteotypic logical. Filter for proteotypic peptides only (default: TRUE)
    #' @param min_peptides integer. Minimum number of peptides per protein (default: 2)
    #' @param remove_contaminants logical. Remove contaminant proteins (default: FALSE)
    #' @return data.frame. Filtered protein group data
    filter_proteins = function(proteotypic = TRUE, min_peptides = 2, remove_contaminants = FALSE) {
      if (!self$parameters$protein_group) {
        stop("Protein group data not available in input data", call. = FALSE)
      }
      
      result <- self$data
      
      # Filter for proteotypic peptides first
      if (proteotypic && "Proteotypic" %in% names(result)) {
        result <- dplyr::filter(result, Proteotypic == 1)
      }

      # Remove contaminants if requested
      if (remove_contaminants && !is.null(self$contaminants)) {
        if (is.data.frame(self$contaminants)) {
          # If contaminants is a data.frame (from Frankenfield dataset)
          contaminant_ids <- self$contaminants$Protein.Ids
        } else {
          # If contaminants is a character vector
          contaminant_ids <- self$contaminants
        }
        
        if (length(contaminant_ids) > 0) {
          # Check if contaminants are already annotated in the data
          if ("is_contaminant" %in% names(result)) {
            # Use existing annotation
            result <- dplyr::filter(result, !is_contaminant)
          } else {
            # Filter directly using Protein.Ids
            result <- dplyr::filter(result, !(Protein.Ids %in% contaminant_ids))
          }
        }
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
      if (!is.null(self$metadata) && nrow(self$metadata) > 0) {
        result <- dplyr::left_join(result, self$metadata, by = c("File.Name", "Run"))
      }
      
      return(result)
    },

    #' @description
    #' Filter and summarize gene-level data
    #' @param proteotypic logical. Filter for proteotypic peptides only (default: TRUE)
    #' @param min_peptides integer. Minimum number of peptides per gene (default: 2)
    #' @param remove_contaminants logical. Remove contaminant proteins (default: FALSE)
    #' @return data.frame. Filtered gene 
    filter_genes = function(proteotypic = TRUE, min_peptides = 2, remove_contaminants = FALSE) {
      if (!self$parameters$genes) {
        stop("Gene data not available in input data", call. = FALSE)
      }
      
      result <- self$data
      
      # Filter for proteotypic peptides first
      if (proteotypic && "Proteotypic" %in% names(result)) {
        result <- dplyr::filter(result, Proteotypic == 1)
      }
      
      # Remove contaminants if requested
      if (remove_contaminants && !is.null(self$contaminants)) {
        if (is.data.frame(self$contaminants)) {
          # If contaminants is a data.frame (from Frankenfield dataset)
          contaminant_ids <- self$contaminants$Protein.Ids
        } else {
          # If contaminants is a character vector
          contaminant_ids <- self$contaminants
        }
        
        if (length(contaminant_ids) > 0) {
          # Check if contaminants are already annotated in the data
          if ("is_contaminant" %in% names(result)) {
            # Use existing annotation
            result <- dplyr::filter(result, !is_contaminant)
          } else {
            # Filter directly using Protein.Ids
            result <- dplyr::filter(result, !(Protein.Ids %in% contaminant_ids))
          }
        }
      }
      
      # Check required columns exist
      required_cols <- c("File.Name", "Run", "Genes", "Stripped.Sequence")
      missing_cols <- setdiff(required_cols, names(result))
      if (length(missing_cols) > 0) {
        stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
      }
      
      # Group and summarize at gene level
      group_cols <- intersect(
        c("File.Name", "Run", "Genes", "Genes.Quantity", "Genes.Normalised", "Genes.MaxLFQ"),
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
      if (!is.null(self$metadata) && nrow(self$metadata) > 0) {
        result <- dplyr::left_join(result, self$metadata, by = c("File.Name", "Run"))
      }
      
      return(result)
    }
  )
)