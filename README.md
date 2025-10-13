# DiaNNData R6 Class-Based Proteomics Analysis

A sophisticated R6 class-based approach for processing and analyzing DIA-NN (Data Independent Acquisition Neural Networks) mass spectrometry data with built-in quality control, contaminant filtering, and visualization capabilities.

## Overview

This project implements a class-based selection algorithm using a custom `DiaNNData` R6 class to streamline proteomics data analysis workflows. The approach provides an object-oriented framework for handling complex proteomics datasets with consistent methodology and reproducible results.

## Key Features

### R6 Class Architecture
- **Encapsulated Data Management**: All data components (precursors, proteins, genes, metadata) managed within a single object
- **Built-in Methods**: Integrated functions for filtering, annotation, and analysis
- **State Persistence**: Object maintains data state throughout the analysis pipeline
- **Extensible Design**: Easy to add new methods and functionality

### Core Capabilities
- DIA-NN report import and processing
- Automated metadata parsing and sample annotation
- Comprehensive contaminant detection and filtering
- Tryptic digestion analysis
- Quality control visualization suite
- Protein and gene-level data aggregation

## Class Methods

### Data Loading and Initialization
```r
# Initialize DiaNNData object
set1 <- DiaNNData$new(diann_data)

# Load contaminant databases
set1$load_contaminants("contaminant_annotation.csv")
set1$annotate_contaminants()
```

### Filtering Operations
```r
# Filter precursors (with/without contaminant removal)
set1$filter_precursors(remove_contaminants = FALSE)

# Filter proteins and genes
set1$filter_proteins(remove_contaminants = TRUE)
set1$filter_genes(remove_contaminants = TRUE)
```

### Analysis Methods
```r
# Perform tryptic digestion analysis
set1$digest_trypsin()

# Access processed data
precursor_data <- set1$precursors
protein_data <- set1$protein_group
metadata <- set1$metadata
```

## Sample Metadata Structure

The algorithm automatically parses sample information from DIA-NN filenames:

| Patient ID | Original ID | Tissue Type | Fractions |
|------------|-------------|-------------|-----------|
| PTN001     | RX847       | PB1, BM, PB | F1, F2, F3 |
| PTN002     | QM394       | NA          | F1, F2, F3 |
| PTN003     | KL672       | PB, BM      | F1, F2, F3 |
| PTN004     | NV138       | PB          | F1, F2, F3 |
| PTN005     | ZT925       | -           | F1, F2, F3 |
| Healthy    | BH561       | PB          | F1, F2, F3 |

**Tissue Types:**
- `PB` = Peripheral Blood
- `BM` = Bone Marrow
- `NA` = Not Available/Unspecified

## Quality Control Metrics

The class-based approach provides comprehensive QC analysis:

### 1. Contaminant Analysis
- Intensity-based contaminant quantification
- Type-specific contaminant distribution
- Sample-wise contamination assessment

### 2. Precursor Characteristics
- Charge state distribution analysis
- Missed cleavage assessment
- Tryptic digestion efficiency

### 3. Protein Identification
- Ranked intensity plots
- Sample-wise protein counts
- Cross-sample comparison metrics

## Dependencies

```r
# Core packages
library(tidyverse)      # Data manipulation and visualization
library(R6)             # Object-oriented programming

# Proteomics-specific
library(diann)          # DIA-NN data import
library(protti)         # Protein analysis toolkit
library(streamlineR)    # Analysis utilities

# Visualization and statistics
library(viridis)        # Color palettes
library(colorspace)     # Color utilities
library(rstatix)        # Statistical analysis
library(forcats)        # Factor manipulation
```

## File Structure

```
project/
├── JR-052_class based selection algorithm2.qmd
├── diann_r6_class.R               # Custom R6 class definition
├── report.tsv                     # DIA-NN output file
├── contaminant_annotation.csv     # Contaminant database
└── outputs/
    ├── contaminants_per_sample.png
    ├── contaminants_per_type.png
    ├── precursor_charge.png
    └── missed_cleavages.png
```

## Usage

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

## Output Visualizations

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
