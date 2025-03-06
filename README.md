---
title: "Making Bioinformatics Analyses Reproducible"
author: "Jamie Billington"
format: 
  revealjs:
    theme: default
    slide-number: true
    transition: none
---

## Workshop Overview

1. Why Reproducibility Matters
2. Project Organization
3. Documentation
4. Environment Management
5. Hands-on Project Cleanup

---

---
If you only take away three thing's from today:
- Use `here()` and relative paths in your R scripts (No more ~/Downloads/da1/etc)
- Use `renv()` to record your projects dependencies
- There's early you tidy up your projects (and put them in git) the better


---


# 1. Why Reproducibility Matters

> "An article about computational science **with a means for reporducing it out** in a scientific publication is **not** the scholarship itself, it is merely **advertising** of the scholarship..."
>
> -- Jonathan Buckheit and David Donoho, 1995 / Me, 2025

---
## Benefits of reproducible code

* **Scientific integrity**: Results can be verified
* **Collaboration**: Others can build on your work
* **Efficiency**: Less time spent recreating past work
* **Impact**: More citations for usable code

---

## Evolving Journal Requirements

* **2010s**: Code "available upon request"
* **2015-2020**: Code in supplementary materials or repository
* **2020+**: Executable code, containers, and code capsules
* **Today**: Some journals running automated testing of submitted code

---

## Journal Requirements Examples

* **Nature**: Code must be "available in a persistent repository"
* **Science**: Code availability statement, encouraged to use containers
* **PLOS**: Code must be in a repository with a DOI
* **eLife**: Docker containers encouraged, Code Execution Editor
* **Genome Biology**: Code must be under version control, containers preferred

---
## The current state of things
CodeOcean: [Code Ocean](https://codeocean.com/capsule/9835100/tree/v1)

---

## Reproducing Your Own Work

### Common scenarios:

* Reviewer asks for changes 6 months after submission
* Need to apply existing analysis to new data
* Collaborator wants to run your analysis
* Your future self needs to revisit analysis


# 2. Project Organization

## Project Structure Principles

1. Intuitive for a stranger
2. "Separation of concerns"
3. Consistent naming
4. Self-documentation

---

## Directory Structure

```
project/
├── README.md
├── LICENSE
├── data/
│   ├── raw/
│   └── processed/
├── src/
│   ├── preprocessing/
│   ├── analysis/
│   └── visualization/
├── results/
│   ├── figures/
│   └── tables/
└── environment/
```

---
Sequencing 
```
01_do_x.sh
02_do_y.sh
03_do_z.sh
```
---

## What to Include & Exclude

**Include**:

* Code
* Small metadata files
* Configuration
* Documentation (At the very least a README)
* Environment files

**Exclude** (use `.gitignore`):

* Raw genomic data
* Large intermediate files
* Sensitive information
* System files

---

## README.md Best Practices

A good README should include:

1. Project title and description
2. Installation instructions
3. Usage instructions
4. Input and output descriptions
5. Examples
6. Dependencies
7. Contact information

---

## Example README.md

```markdown
# Differential Expression Analysis of COVID-19 Lung Samples

## Overview
This repository contains code for analyzing RNA-seq data from COVID-19 patients.

## Requirements
- R 4.1.0 or higher
- Bioconductor 3.13
- DESeq2, edgeR, limma packages

## Installation
1. Clone this repository
2. Install R dependencies: `Rscript install_deps.R`

## Usage
1. Place fastq files in `data/raw/`
2. Run preprocessing: `Rscript src/preprocessing/quality_control.R`
3. Run differential expression: `Rscript src/analysis/diff_exp.R`

## Results
Results will be generated in the `results/` directory.
```

---

## Naming Conventions

* **Files**: lowercase, no spaces, descriptive
  * Good: `preprocess_counts.R`, `differential_expression_analysis.R`
  * Bad: `script1.R`, `DE stuff.R`, `final_FINAL_v2_REALLY_FINAL.R`

* **Directories**: logical grouping, consistent pattern

---

## Separating Code, Data, and Results

* **Data**: Raw data should be read-only
* **Code**: Scripts should be modular and focused
* **Results**: Generated programmatically, not manually

::: {.callout-tip}
Use relative paths instead of absolute paths:
```r
# Good
data_path <- file.path("data", "raw", "counts.csv")

# Bad
data_path <- "/home/user/projects/analysis/data/raw/counts.csv"
```
:::

---

## Exercise: Organizing a Project

**Scenario**: You have a directory with:
- Fastq files
- R scripts
- Python scripts
- Excel files with metadata
- PDFs of papers
- Figures in various formats
- Notes in text files

**Task**: Create a proper structure for these files, explaining your decisions.

---
## What Should You Document?

* **Data**: Sources, cleaning steps, versions
* **Methods**: Algorithms, parameters, rationale
* **Results**: Interpretation, validation, limitations
* **Workflow**: Order of operations, dependencies
* **Environment**: Software versions, configurations

---

# 4. Environment Management

## The Problem: Dependency Hell

:::: {.columns}
::: {.column width="45%"}
**Common issues**:
- "Works on my machine"
- Package version conflicts
- System dependencies
- Unspecified requirements
:::

::: {.column width="45%"}
![](https://imgs.xkcd.com/comics/python_environment.png){width=80%}
*Source: XKCD #1987*
:::
::::

---

## Capturing Environments in R

**Using renv**:

```r
# Initialize project
renv::init()

# Install packages
install.packages("tidyverse")
install.packages("DESeq2")

# Snapshot dependencies
renv::snapshot()
```

This creates `renv.lock` which records exact package versions.

---

## Capturing Environment in Python

**Using pip and requirements.txt**:

```bash
# Generate requirements.txt
pip freeze > requirements.txt
```

**Using virtual environments**:
```bash
# Create virtual environment
python -m venv myenv

# Activate
source myenv/bin/activate  # Linux/Mac
myenv\Scripts\activate     # Windows

# Install packages
pip install pandas scikit-learn
```

---

## Using Conda for Environment Management

```bash
# Create environment
conda create -n myproject python=3.9

# Activate
conda activate myproject

# Install packages
conda install pandas scikit-learn

# Export environment
conda env export > environment.yml
```

---

## Introduction to Containers

**What are containers?**
* Isolated, portable environments
* Include OS, dependencies, code
* Consistent across computing environments
* "Write once, run anywhere"

---

## Docker vs. Singularity

:::: {.columns}
::: {.column width="45%"}
**Docker**:
- Industry standard
- Not always allowed on HPC
- Root privileges
- Dockerfile format
:::

::: {.column width="45%"}
**Singularity**:
- HPC-friendly
- No root privileges required
- Can convert Docker containers
- Commonly used in bioinformatics
:::
::::

---

## Anatomy of a Dockerfile

```dockerfile
FROM python:3.9-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    libz-dev

# Install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy code
COPY src/ /app/src/

# Default command
CMD ["python", "src/main.py"]
```

---

## Exercise: Creating Environment Files

**Task**: For a project using:
- R 4.1 with DESeq2, ggplot2, dplyr
- Python 3.9 with pandas, matplotlib, scikit-learn
- samtools, bedtools, and STAR aligner

Create:
1. An renv.lock (for R)
2. A requirements.txt (for Python)
3. A basic Dockerfile

---

# 5. Hands-on Project Reproducibility

## Project Cleanup Checklist

1. **Organization**: Directory structure, file naming
2. **Documentation**: README, in-code comments, notebooks
3. **Environment**: Record the Dependency specification
4. **Code quality**: Consistent style, error handling
5. **Data management**: Raw vs. processed, paths

---

## Common Issues to Address

* Hard-coded paths
* Missing documentation
* Unspecified parameters
* Insufficient error handling
* Incomplete dependency list
* Undocumented data transformations

---

## Example R Script - Before

```r
# DE analysis
setwd("/home/user/projects/covid_analysis/")
data <- read.csv("data.csv")
meta <- read.csv("metadata.csv")
library(DESeq2)
library(ggplot2)

for (i in 1:ncol(data)) {
  if (i > 1) {
    data[,i] = as.numeric(data[,i])
  }
}

dds <- DESeqDataSetFromMatrix(countData=data[,2:ncol(data)], 
                              colData=meta, 
                              design=~condition)
dds <- DESeq(dds)
res <- results(dds)

pdf("plot.pdf")
plotMA(res)
dev.off()

write.csv(res, "results.csv")
```

---

## Example R Script - After

```r
#!/usr/bin/env Rscript
#
# Title: Differential Expression Analysis for COVID-19 Study
# Author: Your Name
# Date: 2023-03-15
# Description: Performs DE analysis on RNA-seq data from COVID-19 patients

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(here)  # For project-relative paths
})

# Define paths using relative references
project_dir <- here::here()
data_path <- file.path(project_dir, "data", "raw", "counts.csv")
meta_path <- file.path(project_dir, "data", "raw", "metadata.csv")
results_dir <- file.path(project_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Read and validate input data
message("Reading input data...")
counts <- read.csv(data_path, row.names = 1)
metadata <- read.csv(meta_path, row.names = 1)

# Ensure input data is valid
if (ncol(counts) != nrow(metadata)) {
  stop("Number of samples in count data does not match metadata")
}

# Convert count columns to numeric
counts[] <- lapply(counts, function(x) as.numeric(as.character(x)))

# Create DESeq2 dataset
message("Creating DESeq2 dataset...")
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~condition
)

# Run DE analysis
message("Running differential expression analysis...")
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)  # Set significance threshold

# Create MA plot
message("Generating visualizations...")
ma_plot_path <- file.path(results_dir, "MA_plot.pdf")
pdf(ma_plot_path, width = 8, height = 6)
plotMA(res, main = "MA Plot of Differential Expression")
dev.off()

# Save results
message("Saving results...")
results_path <- file.path(results_dir, "differential_expression_results.csv")
write.csv(res, results_path)

message("Analysis complete. Results saved to: ", results_dir)
```

---

## Bring Your Own Project!

* Apply what we've learned to your own projects
* Identify areas for improvement
* Work on specific issues
* Get help with implementation

---
