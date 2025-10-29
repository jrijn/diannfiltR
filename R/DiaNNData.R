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
#'
#' @examples
#' \dontrun{
#' # Create a DiaNNData object
#' diann_data <- DiaNNData$new(input_data = your_diann_output)
#' 
#' # Load contaminants and filter data
#' diann_data$load_contaminants("contaminants.txt")
#' diann_data$annotate_contaminants()
#' 
#' # Print summary
#' print(diann_data)
#' }
#'
#' @importFrom R6 R6Class
#' @importFrom dplyr filter group_by left_join mutate n_distinct summarise select
#' @importFrom readr read_csv read_table read_tsv
#' @importFrom tidyr replace_na
#' @export
DiaNNData <- R6Class(
  classname = "DiaNNData",
  
  public = list(
    # Public fields
    #' @field input_data data.frame. Original DIA-NN output data
    input_data = NULL,
    #' @field metadata data.frame. File and run metadata extracted from input data
    metadata = NULL,
    #' @field contaminants data.frame. Contaminant protein annotations (default: Frankenfield et al. 2022)
    contaminants = NULL,
    #' @field parameters list. Available data types detected in input (precursors, proteins, protein_group, genes)
    parameters = NULL,
    
    #' @description
    #' Initialize DiaNNData object
    #' @param input_data A data.frame containing DIA-NN output data
    #' @param contaminants_file Optional path to contaminants file
    #' @return A new DiaNNData object
    initialize = function(input_data, contaminants_file = NULL) {
      # Validate input
      if (!is.data.frame(input_data)) {
        stop("input_data must be a data.frame", call. = FALSE)
      }
      
      self$input_data <- input_data
      
      # Create metadata from File.Name and Run columns
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
      
      # Load contaminants - custom file takes priority, otherwise use default
      if (!is.null(contaminants_file)) {
        self$load_contaminants(contaminants_file)
      } else {
        # Load default Frankenfield et al. 2022 contaminants
        data("Frankenfield_et_al_2022", envir = environment())
        self$contaminants <- Frankenfield_et_al_2022
      }
      
      # Detect available data types based on column patterns
      self$parameters <- list(
        precursors = any(grepl("^Precursor\\.", names(input_data))),
        proteins = any(grepl("^Protein\\.", names(input_data))),
        protein_group = any(grepl("^PG\\.", names(input_data))),
        genes = any(grepl("^Genes\\.", names(input_data)))
      )
    },

    #' @description
    #' Print summary of DiaNNData object
    #' @return Invisibly returns the DiaNNData object
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
      
      cat("Metadata:\n")

      # Show contaminant info
      if (!is.null(self$contaminants)) {
        if (is.data.frame(self$contaminants)) {
          cat(" - ", nrow(self$contaminants), " contaminants loaded\n")
        } else {
          cat(" - ", length(self$contaminants), " contaminants loaded\n")
        }
      }

      # Metadata summary
      if (!is.null(self$metadata) && nrow(self$metadata) > 0) {
        cat(" - ", nrow(self$metadata), " metadata entries\n")
      } else {
        cat(" - No metadata available\n")
      }
      
      invisible(self)
    },

    #' @description
    #' Load contaminants from file
    #' @param file_path Optional path to contaminants file (CSV or tab-separated text)
    #' @return Invisibly returns the DiaNNData object
    load_contaminants = function(file_path = NULL) {
      # If no file path provided, use default contaminants
      if (is.null(file_path)) {
        data("Frankenfield_et_al_2022", envir = environment())
        self$contaminants <- Frankenfield_et_al_2022
        return(invisible(self))
      }

      # Check if file exists
      if (!file.exists(file_path)) {
        warning(paste("Contaminants file not found:", file_path))
        self$contaminants <- data.frame()
        return(invisible(self))
      }
      
      # Determine file format from extension
      file_ext <- tolower(sub(".*\\.", "", basename(file_path)))
      
      if (file_ext == "csv") {
        # Handle CSV format
        contaminants_data <- readr::read_csv(file_path, show_col_types = FALSE)
        
        # Validate required column
        if ("Protein.Ids" %in% names(contaminants_data)) {
          self$contaminants <- contaminants_data
        } else {
          stop("CSV file must contain 'Protein.Ids' column", call. = FALSE)
        }
        
      } else if (file_ext %in% c("txt", "tsv")) {
        # Handle tab-separated text file
        contaminants_data <- readr::read_tsv(
          file_path, 
          col_names = TRUE, 
          show_col_types = FALSE, 
          skip = 1, 
          name_repair = "universal"
        )
        
        # Validate required column
        if ("Protein.Ids" %in% names(contaminants_data)) {
          self$contaminants <- contaminants_data
        } else {
          stop("Text file must contain 'Protein.Ids' column", call. = FALSE)
        }
        
      } else {
        stop("Unsupported contaminants file format. Use CSV or tab-separated text file.", call. = FALSE)
      }
      
      invisible(self)
    }
  )
)
