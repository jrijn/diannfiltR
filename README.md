# diannfiltR

An R package for filtering and analyzing DIA-NN proteomics data using R6 class architecture.

## Installation

```r
# Install from GitHub
devtools::install_github("your-username/diannfiltR")
library(diannfiltR)
```

## Data Structure

Your DIA-NN data should contain these key columns:

```r
# Create sample dataset for testing
sample_data <- data.frame(
  File.Name = rep(c("sample1.raw", "sample2.raw"), each = 100),
  Run = rep(c("run1", "run2"), each = 100),
  Protein.Ids = paste0("P", rep(1:50, 4)),
  Protein.Names = paste0("Protein_", rep(1:50, 4)),
  Stripped.Sequence = paste0("PEPTIDE", rep(1:100, 2)),
  Proteotypic = sample(0:1, 200, replace = TRUE),
  PG.Quantity = rnorm(200, 1000, 200),
  PG.Normalised = rnorm(200, 1000, 200),
  PG.MaxLFQ = rnorm(200, 1000, 200),
  Genes = paste0("Gene_", rep(1:50, 4)),
  Genes.Quantity = rnorm(200, 1000, 200),
  Genes.Normalised = rnorm(200, 1000, 200),
  Genes.MaxLFQ = rnorm(200, 1000, 200)
)
```

## Creating DiaNNData Object

```r
# Initialize with your data
diann_obj <- DiaNNData$new(input_data = sample_data)

# View summary
print(diann_obj)
```

The object automatically creates metadata from `File.Name` and `Run` columns and detects available data types.

## Working with Contaminants

If not otherwise specified, the built-in dataset from Frankenfield et al. (2022) is used to annotate contaminants. 

```r
# Load built-in Frankenfield et al. 2022 contaminant database
data("Frankenfield_et_al_2022")

# Load contaminant IDs from a specified dataset
diann_obj$contaminants <- Frankenfield_et_al_2022$Protein.Ids

# Or load from custom file
# diann_obj$load_contaminants("path/to/contaminants.csv")

# Annotate contaminants in your data
diann_obj$annotate_contaminants()
```

## Filtering Data

### Filter Precursors

```r
precursor_data <- precursors(
  diann_data = diann_obj,
  proteotypic = TRUE,
  remove_contaminants = FALSE
)
```

### Filter Proteins

```r
protein_data <- proteins(
  diann_data = diann_obj,
  proteotypic = TRUE,
  min_peptides = 2,
  remove_contaminants = TRUE
)
```

### Filter Genes

```r
gene_data <- genes(
  diann_data = diann_obj,
  proteotypic = TRUE,
  min_peptides = 2,
  remove_contaminants = TRUE
)
```

## Example Workflow

```r
library(diannfiltR)

# Load your data
data <- read.csv("your_diann_output.csv")

# Initialize object and load contaminants
diann_obj <- DiaNNData$new(input_data = data)
data("Frankenfield_et_al_2022")
diann_obj$load_contaminants(Frankenfield_et_al_2022)

# Filter data at different levels
precursors_filtered <- precursors(diann_obj, proteotypic = TRUE)
proteins_filtered <- proteins(diann_obj, min_peptides = 2, remove_contaminants = TRUE)
genes_filtered <- genes(diann_obj, min_peptides = 2, remove_contaminants = TRUE)

# View results
cat("Precursors:", nrow(precursors_filtered), "rows\n")
cat("Proteins:", nrow(proteins_filtered), "rows\n")
cat("Genes:", nrow(genes_filtered), "rows\n")
```

# Load contaminants

diann$load_contaminants("contaminants.txt")

# Filter and process data

diann$
  filter_precursors(proteotypic = TRUE, remove_contaminants = TRUE)$
  digest_trypsin()

# View summary

diann

```

## Features

- **Data filtering**: Filter at precursor, protein, and gene levels
- **Contaminant handling**: Load and remove contaminant proteins
- **Proteomics calculations**: Calculate missed cleavages, peptide counts
- **Method chaining**: Fluent interface for data processing workflows

## Methods

### Core Methods

- `DiaNNData$new(input_data, contaminants_file)` - Create new object
- `load_contaminants(file_path)` - Load contaminant proteins
- `print()` - Display object summary

### Filtering Methods

- `filter_precursors(proteotypic, remove_contaminants)` - Filter precursor data
- `filter_proteins(proteotypic, min_peptides, remove_contaminants)` - Filter protein data
- `filter_genes(proteotypic, min_peptides, remove_contaminants)` - Filter gene data

### Analysis Methods

- `annotate_contaminants()` - Mark contaminant proteins in data
- `digest_trypsin()` - Calculate missed cleavages

## License

MIT License
===========

### Basic Workflow

1. **Initialize the analysis environment:**

```r
source("diann_r6_class.R")
data <- diann_load("report.tsv")
analysis <- DiaNNData$new(data)
```

2. **Configure metadata and contaminants:**

```r
# Metadata is automatically parsed from filenames
analysis$load_contaminants("contaminant_annotation.csv")
analysis$annotate_contaminants()
```

3. **Perform quality control:**

```r
analysis$filter_precursors(remove_contaminants = FALSE)
# Generate QC plots
```

4. **Filter and analyze data:**

```r
analysis$filter_proteins(remove_contaminants = TRUE)
analysis$filter_genes(remove_contaminants = TRUE)
analysis$digest_trypsin()
```

### Advanced Features

The R6 class architecture enables:

- **Method chaining** for streamlined workflows
- **Custom method addition** for specific analyses
- **State tracking** throughout the analysis pipeline
- **Consistent data handling** across different experiments

## Advantages of Class-Based Approach

### 1. **Encapsulation**

- All related data and methods contained in single object
- Prevents data inconsistencies and naming conflicts
- Simplifies complex analysis workflows

### 2. **Reproducibility**

- Standardized methods ensure consistent analysis
- Object state can be saved and restored
- Clear audit trail of applied transformations

### 3. **Scalability**

- Easy to extend with new methods and functionality
- Handles multiple datasets with consistent interface
- Supports batch processing workflows

### 4. **Maintainability**

- Centralized code organization
- Clear separation of concerns
- Easier debugging and testing

MIT License

## Contributing

Please report bugs and feature requests on the GitHub issues page.

## Citation

If you use `diannfiltR` in your research, please cite:

- The DIA-NN software: Demichev, V., Messner, C.B., Vernardis, S.I. et al. Nature Methods 17, 41–44 (2020).
- The contaminant database: Frankenfield et al. 2022 (if using built-in contaminants)

The algorithm generates comprehensive quality control plots:

- **Contaminant Analysis**: Sample-wise and type-specific contamination levels
- **Precursor Quality**: Charge state distributions and digestion efficiency
- **Protein Coverage**: Identification rates and intensity rankings
- **Sample Comparison**: Cross-condition protein detection rates

## Future Development

The R6 framework supports extension for:

- Differential expression analysis methods
- Statistical testing integration
- Pathway analysis capabilities
- Multi-dataset comparison tools
- Interactive visualization components

## License

This project uses open-source R packages and is intended for research purposes. Please cite appropriate packages and methods when using this code.
