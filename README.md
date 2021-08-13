
# NanoStringNCTools

## Overview

The NanoStringNCTools package contains tools for analyzing data from 
NanoString nCounter Technology. It provides functions to read,
quality control (QC) and normalize starting from Nanostring 
RCC and RLF files generated from the NanoString nCounter.

It contains the definition of the NanoStringRCCSet which inherits from
Biobase’s ExpressionSet class and complimentary functions that will help
in starting analysis.

## Installation

### Download the release version from Bioconductor
<https://bioconductor.org/packages/release/bioc/html/NanoStringNCTools.html>

### Install the release version from Bioconductor
``` r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

BiocManager::install(version="release")

BiocManager::install("NanoStringNCTools")
```

### Install the development version from GitHub
``` r
devtools::install_github("Nanostring-Biostats/NanoStringNCTools", 
                         build_vignettes = TRUE, ref = "dev")
```

## Documentation

To learn how to start using NanoStringNCTools, view documentation for
the version of this package installed in your system, start R and enter:

``` r
browseVignettes("NanoStringNCTools")
```

## Branches
The release version on Bioconductor is the stable version.
<https://bioconductor.org/packages/release/bioc/html/NanoStringNCTools.html>

The devel version on Bioconductor is upstream of master on GitHub.
It is under active development and no guarantee is made on usability
at any given time.

The dev branch on GitHub is under active development and no guarantee 
is made on usability at any given time.
