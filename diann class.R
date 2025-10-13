# Load the methods package
library(methods)

# Define the S4 class
setClass("DiaNNData",
         slots = list(
          inputData = "data.frame",  # Input data used for predictions
          metadata = "list",       # Metadata about the model or data
          parameters = "list", # list of precursor, protein, protein.group, genes level availability
          precursors = "data.frame",   # Precursor level output
          proteins = "data.frame",   # Protein level output
          protein.group = "data.frame",   # Protein group level 
          genes = "data.frame"   # Genes level output
         ))

# Define a constructor function
DiaNNData <- function(inputData) {

  # Define the metadata
  metadata <- list(
    modelType = "DiaNN"
  )

  # Define the parameters
  parameters <- list( 
    precursors = any(grepl("^Precursor.", names(inputData))), 
    proteins = any(grepl("^Protein.", names(inputData))), 
    protein.group = any(grepl("^PG.", names(inputData))), 
    genes = any(grepl("^Genes.", names(inputData))) )

  # Construct the object
  new("DiaNNData",
      inputData = inputData,
      metadata = metadata,
      parameters = parameters
      )
}

# Define accessor methods
setMethod("show", "DiaNNData", function(object) {
  cat("DiaNNData Object:\n")
  cat("Input type:", object@metadata$modelType, "\n")
  cat("Parameters:\n", 
  " - precursors:\t", object@parameters$precursors, "\n",
  " - proteins:\t", object@parameters$proteins, "\n",
  " - groups:\t", object@parameters$protein.group, "\n",
  " - genes:\t", object@parameters$genes, "\n")
  cat("Input Data:\n",
  " - ", ncol(object@inputData), "columns,", nrow(object@inputData), "rows\n",
  " - ", length(unique(object@inputData$Stripped.Sequence)), "peptides\n", 
  " - ", length(unique(object@inputData$Protein.Names)), "proteins\n",
  " - ", length(unique(object@inputData$Protein.Group)), "protein groups\n",
  " - ", length(unique(object@inputData$Genes)), "genes\n")
})

setMethod("precursors", "DiaNNData", function(object) {
  if (object@parameters$precursors) {
    return(object@inputData)
  } else {
    stop("Precursors not available.")
  }
})

setMethod("filterProteins", "DiaNNData", function(object, proteotypic = T, min.peptides = 2) {
 
  if (object@parameters$proteins & proteotypic) {

    out <- object@inputData %>%
      filter(Proteotypic == T) %>%
      group_by(File.Name, Run, Protein.Ids, Protein.Names, 
      PG.Quantity, PG.Normalised, PG.MaxLFQ, Proteotypic) %>%
      summarise(peptides = n_distinct(Stripped.Sequence))
    out <- filter(out, peptides >= min.peptides)
    return(out)

  } else if (object@parameters$proteins & !proteotypic) {

    out <- object@inputData %>%
      # filter(Proteotypic == T) %>%
      group_by(File.Name, Run, Protein.Ids, Protein.Names, 
      PG.Quantity, PG.Normalised, PG.MaxLFQ, Proteotypic) %>%
      summarise(peptides = n_distinct(Stripped.Sequence))
    out <- filter(out, peptides >= min.peptides)
    return(out)

  } else if (!object@parameters$proteins) {
    stop("Proteins not available.")
  }
})

setMethod("filterGenes", "DiaNNData", function(object, proteotypic = T, min.peptides = 2) {
 
  if (object@parameters$genes & proteotypic) {

    out <- object@inputData %>%
      filter(Proteotypic == T) %>%
      group_by(File.Name, Run, Genes,
      Genes.Quantity, Genes.Normalised, Genes.MaxLFQ, Proteotypic) %>%
      summarise(
        peptides = n_distinct(Stripped.Sequence),
        proteins = n_distinct(Protein.Ids))
    out <- filter(out, peptides >= min.peptides)
    return(out)

  } else if (object@parameters$genes & !proteotypic) {

    out <- object@inputData %>%
      # filter(Proteotypic == T) %>%
      group_by(File.Name, Run, Genes,
        Genes.Quantity, Genes.Normalised, Genes.MaxLFQ, Proteotypic) %>%
        summarise(
          peptides = n_distinct(Stripped.Sequence),
          proteins = n_distinct(Protein.Ids))
    out <- filter(out, peptides >= min.peptides)
    return(out)

  } else if (!object@parameters$genes) {
    stop("Genes not available.")
  }
})

# Example usage
# metadata <- list(modelName = "DiaNN", version = "1.0")
# parameters <- data.frame(param1 = c(0.1, 0.2), param2 = c(0.3, 0.4))
# predictions <- matrix(c(0.5, 0.6, 0.7, 0.8), nrow = 2)
# inputData <- data.frame(feature1 = c(1, 2), feature2 = c(3, 4))

# Create an instance of the class
# diaNNDataInstance <- DiaNNData(inputData)

# Show the instance
# show(diaNNDataInstance)
