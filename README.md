# Inferential Seamless Designs for Treatment Selection during a Phase II/III Study in Oncology

Student: Manuel Pfister

University of Zurich and Bristol Myers Squibb Ltd.

Year: 2024

## Overview

This repository contains the research and code associated with the STA495 Master Thesis by Manuel Pfister. The thesis explores efficient statistical designs for multi-arm multi-stage (MAMS) clinical trials, focusing on the "Drop-the-Losers" approach to improve drug development in oncology.

The repository includes:

- R scripts for generating data, running models, and simulations.
- The final thesis report compiled using knitr and LaTeX.
- Necessary datasets and model results for reproducibility.
- Use HS24_Masterthesis_BMS.Rproj Project in the report folder for compilation and running code.
- renv folder for managing the R environment to ensure reproducibility.


**NOTE:**  
The simulation data files included in this release are too large to be hosted directly on GitHub. Therefore, you can download them from the following [[Google Drive link](https://drive.google.com/drive/folders/1DX-D1qHVUcbXPk3v5d8naI3KoSB7UR7c?usp=sharing)].  

### Instructions for Use:
1. **For Report Compilation:**  
   - Copy the entire `data` folder into the `report` folder. This is required to compile the report successfully.

2. **For Data Analysis:**  
   - Alternatively, copy the `data` folder into the main directory if you wish to use the files for general data analysis.

Please note that the data files are large `.RData` files generated from extensive simulations.  
- Ensure sufficient system resources (e.g., memory) when loading and processing these datasets in R.  
- **Important:** Make sure to unzip the data file before use!


### Project Setup with renv
renv helps manage R package dependencies for reproducible projects.

1. Install renv:
```r
install.packages("renv")
```
2. Initialize the project:
```r
renv::init()
```
3. Save dependencies:
```r
renv::snapshot()
```

4. Restore environment on a new system:
```r
renv::restore()
```
5. Add renv.lock to Git and ignore renv/library in .gitignore.

For more, see the official documentation: https://rstudio.github.io/renv/

### Install Packages manually
```r
required_packages <- c("mvtnorm", "rpact", "doParallel", "foreach", 
                       "survival", "renv", "knitr", "dplyr", 
                       "ggplot2", "xtable")

# Install all packages
install.packages(required_packages)
```

## Repository Structure

- **`HS24_Masterthesis_BMS.Rproj`**: The R project file for organizing the working environment.

  - **`admin/`**: Contains proposals and organizational documents.
  - **`code/`**: Contains R scripts for data generation, simulations, and calculations.
  - **`data/`**: Contains datasets used in the analysis.
  - **`report/`**: Contains the LaTeX and `.Rnw` files for the thesis report.
    - **`report_master.Rnw`**: The main file for compiling the thesis.
    - **`figures/`**: Folder containing figures for the report.
  - **`literature/`**: Contains relevant references.
  - **`presentation/`**: Contains presentation material.


## Contact

For any questions or collaboration inquiries, please contact Manuel Pfister at [manuel.pfister@uzh.ch](mailto:manuel.pfister@uzh.ch) or [manuel.pfister@bms.com](mailto:manuel.pfister@bms.com).
