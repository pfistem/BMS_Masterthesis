# Inferential Seamless Designs for Treatment Selection during a Phase II/III Study in Oncology

## Overview

This repository contains the research and code associated with the STA495 Master Thesis by Manuel Pfister. The thesis explores efficient statistical designs for multi-arm multi-stage (MAMS) clinical trials, focusing on the "drop-the-losers" approach to improve drug development in oncology.


## Introduction to `renv`

`renv` is an R package that helps manage project-specific R package libraries. It allows you to create isolated environments for your R projects, ensuring that the packages used in one project do not interfere with those used in another. This is particularly useful for maintaining reproducibility and consistency in data analysis workflows.

## Getting Started with `renv`

### Installation

You can install `renv` from CRAN using the following command:

```r
install.packages("renv")
```
### Initializing a Project
To start using renv in a project, navigate to your project's directory and initialize renv:

```r
renv::init()
```

This will set up a new project-specific library and snapshot the current state of your R packages. It will also create an renv folder in your project directory to store the environment configuration.

### Adding and Managing Packages
To install a new package in your project-specific library, use renv::install():

```r
renv::install("dplyr")
```

This will install 'dplyr' and add it to your project's renv.lock file, which records the state of your project's library.

### Snapshotting the Project
After installing or updating packages, you should snapshot the project's library to update the renv.lock file:
```r
renv::snapshot()
```
This ensures that the renv.lock file reflects the current state of the project library.

### Restoring the Project Library
If you clone a project that uses renv, you can restore the project library as specified in the renv.lock file by running:

```r
renv::restore()
```
This will install the required package versions as listed in the lockfile, ensuring reproducibility.

### Managing Dependencies
To check for any missing dependencies or update the renv.lock file, you can use:
```r
renv::status()
```
This provides an overview of the current state of your project's dependencies.


## Contact

For any questions or collaboration inquiries, please contact Manuel Pfister at [manuel.pfister@uzh.ch](mailto:manuel.pfister@uzh.ch) or [manuel.pfister@bms.com](mailto:manuel.pfister@bms.com).
