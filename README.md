# diannr

An R package for analyzing DIA-NN (Data Independent Acquisition Neural Networks) proteomics output data.

## Installation

You can install the development version of diannr from GitHub:

```r
# install.packages("devtools")
devtools::install_github("yourusername/diannr")
```

## Usage

```r
library(diannr)

# Load DIA-NN output data
data <- read.csv("your_diann_output.csv")

# Create DiaNNData object
diann <- DiaNNData$new(input_data = data)

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