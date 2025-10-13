library(tidyverse, quietly = T)

contaminants <- read_table(file = "C:\\Users\\jrijn8\\OneDrive - Prinses Maxima Centrum\\Lab_data\\T-ALL surface proteomics\\231108_JR-052_patient samples surfaceome\\contaminant_protein_ids.txt", 
col_names = F, show_col_types = F) %>%
  pull()

# Define the reference class
DiaNNData <- setRefClass(
  Class = "DiaNNData",

  fields = list(
    inputData = "data.frame",  # Input data used for predictions
    metadata = "data.frame",   # Metadata about the file names
    contaminants = "data.frame",    # vector of contaminant Uniprot IDs
    parameters = "list",       # List of precursor, protein, protein.group, genes level availability
    precursors = "data.frame", # Precursor level output
    proteins = "data.frame",   # Protein level output
    protein.group = "data.frame", # Protein group level
    genes = "data.frame"       # Genes level output  
    ),

  methods = list(

    # Define the metadata
    initialize = function(inputData) {
      files <- unique(inputData$File.Name)
      runs <- unique(inputData$Run)
      metadata <<- data.frame(
        File.Name = files,
        Run = runs
      )

      # Define the contaminants data frame.
      contaminants <<- data.frame(
        Protein.Ids = c(),
        Genes = c()
      )

      # Define the parameters
      parameters <<- list(
        precursors = any(grepl("^Precursor.", names(inputData))),
        proteins = any(grepl("^Protein.", names(inputData))),
        protein.group = any(grepl("^PG.", names(inputData))),
        genes = any(grepl("^Genes.", names(inputData)))
      )
      
      # Initialize the fields
      inputData <<- inputData
    },

    show = function() {
      cat("DiaNNData Object:\n")
      cat("Input type:", metadata$modelType, "\n")
      cat("Parameters:\n",
          " - precursors:\t", parameters$precursors, "\n",
          " - proteins:\t", parameters$proteins, "\n",
          " - groups:\t", parameters$protein.group, "\n",
          " - genes:\t", parameters$genes, "\n")
      cat("Input Data:\n",
          " - ", ncol(inputData), "columns,", nrow(inputData), "rows\n",
          " - ", length(unique(inputData$Stripped.Sequence)), "peptides\n",
          " - ", length(unique(inputData$Protein.Names)), "proteins\n",
          " - ", length(unique(inputData$Protein.Group)), "protein groups\n",
          " - ", length(unique(inputData$Genes)), "genes\n")
    },

    filterPrecursors = function(proteotypic = TRUE, remove.contaminants = F) {
      
      if (parameters$precursors & proteotypic) {
        precursors <<- filter(inputData, Proteotypic == 1)
        precursors <<- left_join(precursors, metadata, by = c("File.Name", "Run"))

      } else if (parameters$protein.group & !proteotypic) {
        precursors <<- left_join(precursors, metadata, by = c("File.Name", "Run"))

      } else if (!parameters$protein.group) {
        stop("Genes not available.")
      }

    },

    filterProteins = function(proteotypic = TRUE, min.peptides = 2, remove.contaminants = F) {
      
      if (remove.contaminants & is.null(contaminants)) {
        stop("object$contaminants == Null")
      } else if (!is.logical(remove.contaminants)) {
        stop("remove.contaminants must be logical")
      }

      if (remove.contaminants) {
        protein.group <<- filter(inputData, is.na(contaminant))
      } else {
        protein.group <<- inputData
      }
      
      if (parameters$protein.group & proteotypic) {
        protein.group <<- inputData %>%
          filter(Proteotypic == 1) %>%
          group_by(File.Name, Run, Protein.Ids, Protein.Names,
                    PG.Quantity, PG.Normalised, PG.MaxLFQ, Proteotypic, contaminant) %>%
          summarise(
            peptides = n_distinct(Stripped.Sequence),
            proteins = n_distinct(Protein.Ids),
            .groups = "keep")
        protein.group <<- filter(protein.group, peptides >= min.peptides)
        protein.group <<- left_join(protein.group, metadata, by = c("File.Name", "Run"))

      } else if (parameters$protein.group & !proteotypic) {
        genes <<- inputData %>%
          group_by(File.Name, Run, Protein.Ids, Protein.Names,
                    PG.Quantity, PG.Normalised, PG.MaxLFQ, Proteotypic, contaminant) %>%
          summarise(
            peptides = n_distinct(Stripped.Sequence),
            proteins = n_distinct(Protein.Ids))
        protein.group <<- filter(protein.group, peptides >= min.peptides)
        protein.group <<- left_join(protein.group, metadata, by = c("File.Name", "Run"))

      } else if (!parameters$protein.group) {
        stop("Genes not available.")
      }

    },
                          
    filterGenes = function(proteotypic = TRUE, min.peptides = 2, remove.contaminants = F) {
        
      if (remove.contaminants & is.null(contaminants)) {
        stop("object$contaminants == Null")
      } else if (!is.logical(remove.contaminants)) {
        stop("remove.contaminants must be logical")
      }

      if (remove.contaminants) {
        genes <<- filter(inputData, is.na(contaminant))
      } else {
        protein.group <<- inputData
      }
      
      if (parameters$genes & proteotypic) {
        genes <<- inputData %>%
          filter(Proteotypic == 1) %>%
          group_by(File.Name, Run, Genes,
                    Genes.Quantity, Genes.Normalised, Genes.MaxLFQ, Proteotypic, contaminant) %>%
          summarise(
            peptides = n_distinct(Stripped.Sequence),
            proteins = n_distinct(Protein.Ids),
            .groups = "keep")
        genes <<- filter(genes, peptides >= min.peptides)
        genes <<- left_join(genes, metadata, by = c("File.Name", "Run"))

      } else if (parameters$genes & !proteotypic) {
        genes <<- inputData %>%
          group_by(File.Name, Run, Genes,
                    Genes.Quantity, Genes.Normalised, Genes.MaxLFQ, Proteotypic, contaminant) %>%
          summarise(
            peptides = n_distinct(Stripped.Sequence),
            proteins = n_distinct(Protein.Ids))
        genes <<- filter(genes, peptides >= min.peptides)
        genes <<- left_join(genes, metadata, by = c("File.Name", "Run"))

      } else if (!parameters$genes) {
        stop("Genes not available.")
      }

    },

    annotateContaminants = function() {
      if (!is.null(contaminants)) {
      inputData <<- inputData %>%
        left_join(., contaminants, by = c("Protein.Group" = "Protein.Ids"), relationship = "many-to-many")
      } else {
        stop("Contaminants are not defined")
      }
    },

    digestTrypsin = function() {
      y <- sapply(precursors$Stripped.Sequence, function(x) strsplit(x, "(?!P)(?<=[RK])", perl = T))
      y <- unlist(lapply(y, length))
      precursors$missed.cleavages <<- y-1
    }

  ))
