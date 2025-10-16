#' DiaNNData R6 Class for DIA-NN Output Analysis
#'
#' @description
#' An R6 class for processing and analyzing DIA-NN (Data Independent Acquisition
#' Neural Networks) proteomics output data. Provides methods for data filtering,
#' contaminant removal, and proteomics-specific transformations.
#'
#' @details
#' This class provides a structured approach to handle DIA-NN output data with
#' methods for filtering at different levels (precursors, proteins, genes) and
#' handling contaminant proteins.
#'
#' @examples
#' \dontrun{
#' # Create a DiaNNData object
#' diann_data <- DiaNNData$new(input_data = your_diann_output)
#' 
#' # Load contaminants and filter data
#' diann_data$load_contaminants("contaminants.txt")
#' diann_data$filter_precursors(proteotypic = TRUE, remove_contaminants = TRUE)
#' 
#' # Print summary
#' print(diann_data)
#' }
#'
#' @importFrom R6 R6Class
#' @importFrom dplyr filter group_by left_join mutate n_distinct summarise
#' @importFrom readr read_csv read_table
#' @importFrom tidyr replace_na
#' @export
DiaNNData <- R6::R6Class(
  classname = "DiaNNData",
  
  public = list(
    #' @field input_data Original DIA-NN output data
    input_data = NULL,

    #' @field object_type Type of object
    object_type = "DiaNNData",
    
    #' @field metadata File and run metadata
    metadata = NULL,
    
    #' @field contaminants Vector of contaminant protein IDs
    contaminants = NULL,
    
    #' @field contaminant_annotations Detailed contaminant annotations
    contaminant_annotations = NULL,
    
    #' @field parameters List of available data types in input
    parameters = NULL,
    
    # #' @field precursors Filtered precursor data
    # precursors = NULL,
    
    # #' @field proteins Filtered protein data
    # proteins = NULL,
    
    # #' @field protein_group Filtered protein group data
    # protein_group = NULL,
    
    # #' @field genes Filtered gene data
    # genes = NULL,
    
    #' Initialize DiaNNData object
    #'
    #' @param input_data A data.frame containing DIA-NN output data
    #' @param contaminants_file Optional path to contaminants file
    #'
    #' @return A new DiaNNData object
    initialize = function(input_data, contaminants_file = NULL) {
      # Validate input
      if (!is.data.frame(input_data)) {
        stop("input_data must be a data.frame", call. = FALSE)
      }
      
      self$input_data <- input_data
      
      # Create metadata
      if (all(c("File.Name", "Run") %in% names(input_data))) {
        files <- unique(input_data$File.Name)
        runs <- unique(input_data$Run)
        self$metadata <- data.frame(
          File.Name = files,
          Run = runs
        )
      } else {
        warning("File.Name and/or Run columns not found in input data")
        self$metadata <- data.frame()
      }
      
      # Load contaminants if file is provided
      if (!is.null(contaminants_file)) {
        self$load_contaminants(contaminants_file)
      } else {
        self$contaminants <- self$load_contaminants()
      }
      
      # Define parameters based on column presence
      self$parameters <- list(
        precursors = any(grepl("^Precursor\\.", names(input_data))),
        proteins = any(grepl("^Protein\\.", names(input_data))),
        protein_group = any(grepl("^PG\\.", names(input_data))),
        genes = any(grepl("^Genes\\.", names(input_data)))
      )
      
    },
    
    #' Load contaminants from file
    #'
    #' @param file_path Path to contaminants file (CSV or text format)
    #'
    #' @return Self (invisibly) for method chaining
    load_contaminants = function(file_path = NULL) {
      # If no file path provided, use default contaminants
      if (is.null(file_path)) {
        self$contaminants <- Frankenfield_et_al_2022$Protein.Ids
        self$contaminant_annotations <- Frankenfield_et_al_2022
        return(invisible(self))

      # Check if file exists
      } else if (!file.exists(file_path)) {
        warning(paste("Contaminants file not found:", file_path))
        self$contaminants <- character(0)
        return(invisible(self))
      }
      
      # Check file extension to determine format
      file_ext <- tolower(sub(".*\\.", "", basename(file_path)))
      if (file_ext == "csv") {
        # Handle CSV format with Protein.Ids and contaminant columns
        contaminants_data <- readr::read_csv(file_path, show_col_types = FALSE)
        
        # Store both the IDs and the annotation data
        if ("Protein.Ids" %in% names(contaminants_data)) {
          self$contaminants <- contaminants_data$Protein.Ids
          self$contaminant_annotations <- contaminants_data
        } else {
          stop("CSV file must contain 'Protein.Ids' column", call. = FALSE)
        }
      } else if (file_ext %in% c("txt", "tsv")) {
        # Handle tab-separated text file with Protein.Ids and contaminant columns
        contaminants_data <- readr::read_tsv(file_path, col_names = TRUE, show_col_types = FALSE, skip = 1, name_repair = "universal")
        
        # Store both the IDs and the annotation data
        if ("Protein.Ids" %in% names(contaminants_data)) {
          self$contaminants <- contaminants_data$Protein.Ids
          self$contaminant_annotations <- contaminants_data
        } else {
          stop("Text file must contain 'Protein.Ids' column", call. = FALSE)
        }
      } else {
        # Handle single column text file format
        self$contaminants <- readr::read_table(
          file_path, 
          col_names = FALSE, 
          show_col_types = FALSE
        ) |> 
          dplyr::pull()
        self$contaminant_annotations <- NULL
      }
      
      invisible(self)
    },
    
    #' Print method for DiaNNData objects
    #'
    #' @return Self (invisibly)
    print = function() {
      cat("DiaNNData Object:\n")
      cat("Parameters:\n")
      cat(" - precursors:\t", self$parameters$precursors, "\n")
      cat(" - proteins:\t", self$parameters$proteins, "\n")
      cat(" - groups:\t", self$parameters$protein_group, "\n")
      cat(" - genes:\t", self$parameters$genes, "\n")
      
      cat("Input Data:\n")
      cat(" - ", ncol(self$input_data), " columns, ", nrow(self$input_data), " rows\n")
      
      # Safe access to columns that might not exist
      if ("Stripped.Sequence" %in% names(self$input_data)) {
        cat(" - ", length(unique(self$input_data$Stripped.Sequence)), " peptides\n")
      }
      if ("Protein.Names" %in% names(self$input_data)) {
        cat(" - ", length(unique(self$input_data$Protein.Names)), " proteins\n")
      }
      if ("Protein.Group" %in% names(self$input_data)) {
        cat(" - ", length(unique(self$input_data$Protein.Group)), " protein groups\n")
      }
      if ("Genes" %in% names(self$input_data)) {
        cat(" - ", length(unique(self$input_data$Genes)), " genes\n")
      }
      
      invisible(self)
    },
    
    #' Annotate contaminants in input data
    #'
    #' @return Self (invisibly) for method chaining
    annotate_contaminants = function() {
      if (length(self$contaminants) == 0) {
        stop("No contaminants loaded. Use load_contaminants() first.", call. = FALSE)
      }
      
      if (!"Protein.Ids" %in% names(self$input_data)) {
        stop("Protein.Ids column not found in input data", call. = FALSE)
      }
      
      # Create contaminant data frame for joining
      if (!is.null(self$contaminant_annotations)) {
        # Use the full annotation data with categories
        contaminant_df <- self$contaminant_annotations |>
          dplyr::mutate(is_contaminant = TRUE)
      } else {
        # Use simple binary annotation
        contaminant_df <- data.frame(
          Protein.Ids = self$contaminants,
          is_contaminant = TRUE
        )
      }
      
      # Join and mark contaminants
      self$input_data <- self$input_data |>
        dplyr::left_join(contaminant_df, by = "Protein.Ids") |>
        dplyr::mutate(is_contaminant = tidyr::replace_na(is_contaminant, FALSE))
      
      invisible(self)
    },
    
    #' Calculate missed cleavages for tryptic peptides
    #'
    #' @return Self (invisibly) for method chaining
    digest_trypsin = function() {
      if (nrow(self$precursors) == 0) {
        stop("No precursor data available. Run filter_precursors() first.", call. = FALSE)
      }
      
      if (!"Stripped.Sequence" %in% names(self$precursors)) {
        stop("Stripped.Sequence column not found in precursor data", call. = FALSE)
      }
      
      # Calculate missed cleavages
      self$precursors$missed_cleavages <- sapply(
        self$precursors$Stripped.Sequence, 
        function(x) {
          # Count K and R not followed by P
          fragments <- strsplit(x, "(?!P)(?<=[RK])", perl = TRUE)[[1]]
          length(fragments) - 1
        }
      )
      
      invisible(self)
    }
  ),
  
  private = list(
    # Remove contaminants from data
    remove_contaminants = function(data) {
      if (length(self$contaminants) == 0) {
        warning("No contaminants loaded")
        return(data)
      }
      
      if (!"Protein.Ids" %in% names(data)) {
        stop("Protein.Ids column not found in data", call. = FALSE)
      }
      
      # If contaminants haven't been annotated yet, do it temporarily
      if (!"is_contaminant" %in% names(data)) {
        contaminant_df <- data.frame(
          Protein.Ids = self$contaminants,
          is_contaminant = TRUE
        )
        
        data <- data |>
          dplyr::left_join(contaminant_df, by = "Protein.Ids") |>
          dplyr::mutate(is_contaminant = tidyr::replace_na(is_contaminant, FALSE))
      }
      
      dplyr::filter(data, !is_contaminant)
    }
  )
)