
<!-- README.md is generated from README.Rmd. Please edit that file -->

# NanoStringNCTools

## Overview

The NanoStringNCTools is a package that contains tools for analyzing
data from NanoString nCounter Technology. It provides functions to read,
perform quality control (QC) and normalization on Nanostring RCC files
generated from the NanoString nCounter platform.

It contains the definition of the NanoStringRCCSet which inherits from
Biobase’s ExpressionSet class and complimentary functions that will help
in starting analysis.

## Installation

You can download the the package from bioconductor from this link
<https://bioconductor.org/packages/devel/bioc/html/NanoStringNCTools.html>

``` r
# Install from Bioconductor
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

BiocManager::install(version="devel")

BiocManager::install("NanoStringNCTools")

# Or the development version from GitHub
# install.packages("devtools")
devtools::install_github("Nanostring-Biostats/NanoStringNCTools", 
                         build_vignettes = TRUE)
```

## Documentation

To learn how to start using NanoStringNCTools, view documentation for
the version of this package installed in your system, start R and enter:

``` r
browseVignettes("NanoStringNCTools")
```
