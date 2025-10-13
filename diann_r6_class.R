library(tidyverse, quietly = TRUE)
library(R6)

DiaNNData <- R6Class(
  classname = "DiaNNData",
  
  public = list(
    # Public fields
    input_data = NULL,
    metadata = NULL,
    contaminants = NULL,
    contaminant_annotations = NULL,
    parameters = NULL,
    precursors = NULL,
    proteins = NULL,
    protein_group = NULL,
    genes = NULL,
    
    # Constructor
    initialize = function(input_data, contaminants_file = NULL) {
      # Validate input
      if (!is.data.frame(input_data)) {
        stop("input_data must be a data.frame")
      }
      
      self$input_data <- input_data
      
      # Create metadata
      files <- unique(input_data$File.Name)
      runs <- unique(input_data$Run)
      self$metadata <- data.frame(
        File.Name = files,
        Run = runs
      )
      
      # Load contaminants if file is provided
      if (!is.null(contaminants_file)) {
        self$load_contaminants(contaminants_file)
      } else {
        self$contaminants <- character(0)
      }
      
      # Define parameters based on column presence
      self$parameters <- list(
        precursors = any(grepl("^Precursor\\.", names(input_data))),
        proteins = any(grepl("^Protein\\.", names(input_data))),
        protein_group = any(grepl("^PG\\.", names(input_data))),
        genes = any(grepl("^Genes\\.", names(input_data)))
      )
      
      # Initialize result data frames
      self$precursors <- data.frame()
      self$proteins <- data.frame()
      self$protein_group <- data.frame()
      self$genes <- data.frame()
    },
    
    # Load contaminants from file
    load_contaminants = function(file_path) {
      if (!file.exists(file_path)) {
        warning(paste("Contaminants file not found:", file_path))
        self$contaminants <- character(0)
        return(invisible(self))
      }
      
      # Check file extension to determine format
      if (tools::file_ext(file_path) == "csv") {
        # Handle CSV format with Protein.Ids and contaminant columns
        contaminants_data <- read_csv(file_path, show_col_types = FALSE)
        
        # Store both the IDs and the annotation data
        self$contaminants <- contaminants_data$Protein.Ids
        self$contaminant_annotations <- contaminants_data
      } else {
        # Handle single column text file format
        self$contaminants <- read_table(
          file_path, 
          col_names = FALSE, 
          show_col_types = FALSE
        ) |> 
          pull()
        self$contaminant_annotations <- NULL
      }
      
      invisible(self)
    },
    
    # Show method
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
    
    # Filter precursors
    filter_precursors = function(proteotypic = TRUE, remove_contaminants = FALSE) {
      if (!self$parameters$precursors) {
        stop("Precursor data not available in input data")
      }
      
      result <- self$input_data
      
      # Filter for proteotypic peptides
      if (proteotypic && "Proteotypic" %in% names(result)) {
        result <- filter(result, Proteotypic == 1)
      }
      
      # Remove contaminants if requested
      if (remove_contaminants) {
        result <- private$remove_contaminants(result)
      }
      
      # Join with metadata
      self$precursors <- left_join(result, self$metadata, by = c("File.Name", "Run"))
      
      invisible(self)
    },
    
    # Filter proteins
    filter_proteins = function(proteotypic = TRUE, min_peptides = 2, remove_contaminants = FALSE) {
      if (!self$parameters$protein_group) {
        stop("Protein group data not available in input data")
      }
      
      result <- self$input_data
      
      # Remove contaminants if requested
      if (remove_contaminants) {
        result <- private$remove_contaminants(result)
      }
      
      # Filter for proteotypic peptides if requested
      if (proteotypic && "Proteotypic" %in% names(result)) {
        result <- filter(result, Proteotypic == 1)
      }
      
      # Group and summarize at protein level
      self$protein_group <- result |>
        group_by(
          File.Name, Run, Protein.Ids, Protein.Names,
          PG.Quantity, PG.Normalised, PG.MaxLFQ
        ) |>
        summarise(
          peptides = n_distinct(Stripped.Sequence),
          proteins = n_distinct(Protein.Ids),
          .groups = "keep"
        ) |>
        filter(peptides >= min_peptides) |>
        left_join(self$metadata, by = c("File.Name", "Run"))
      
      invisible(self)
    },
    
    # Filter genes
    filter_genes = function(proteotypic = TRUE, min_peptides = 2, remove_contaminants = FALSE) {
      if (!self$parameters$genes) {
        stop("Genes data not available in input data")
      }
      
      result <- self$input_data
      
      # Remove contaminants if requested
      if (remove_contaminants) {
        result <- private$remove_contaminants(result)
      }
      
      # Filter for proteotypic peptides if requested
      if (proteotypic && "Proteotypic" %in% names(result)) {
        result <- filter(result, Proteotypic == 1)
      }
      
      # Group and summarize at gene level
      self$genes <- result |>
        group_by(
          File.Name, Run, Genes,
          Genes.Quantity, Genes.Normalised, Genes.MaxLFQ
        ) |>
        summarise(
          peptides = n_distinct(Stripped.Sequence),
          proteins = n_distinct(Protein.Ids),
          .groups = "keep"
        ) |>
        filter(peptides >= min_peptides) |>
        left_join(self$metadata, by = c("File.Name", "Run"))
      
      invisible(self)
    },
    
    # Annotate contaminants
    annotate_contaminants = function() {
      if (length(self$contaminants) == 0) {
        stop("No contaminants loaded. Use load_contaminants() first.")
      }
      
      # Create contaminant data frame for joining
      if (!is.null(self$contaminant_annotations)) {
        # Use the full annotation data with categories
        contaminant_df <- self$contaminant_annotations |>
          mutate(is_contaminant = TRUE)
      } else {
        # Use simple binary annotation
        contaminant_df <- data.frame(
          Protein.Ids = self$contaminants,
          is_contaminant = TRUE
        )
      }
      
      # Join and mark contaminants
      self$input_data <- self$input_data |>
        left_join(contaminant_df, by = "Protein.Ids") |>
        mutate(is_contaminant = replace_na(is_contaminant, FALSE))
      
      invisible(self)
    },
    
    # Digest trypsin
    digest_trypsin = function() {
      if (nrow(self$precursors) == 0) {
        stop("No precursor data available. Run filter_precursors() first.")
      }
      
      if (!"Stripped.Sequence" %in% names(self$precursors)) {
        stop("Stripped.Sequence column not found in precursor data")
      }
      
      # Calculate missed cleavages
      self$precursors$missed_cleavages <- sapply(
        self$precursors$Stripped.Sequence, 
        function(x) {
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
      
      # If contaminants haven't been annotated yet, do it temporarily
      if (!"is_contaminant" %in% names(data)) {
        contaminant_df <- data.frame(
          Protein.Ids = self$contaminants,
          is_contaminant = TRUE
        )
        
        data <- data |>
          left_join(contaminant_df, by = "Protein.Ids") |>
          mutate(is_contaminant = replace_na(is_contaminant, FALSE))
      }
      
      filter(data, !is_contaminant)
    }
  )
)
