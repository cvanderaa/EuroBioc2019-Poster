# EuroBioc2019-Poster

The poster for the EuroBioc2019 conference about the 'scpdata' package. The poster was presented on the December 9, 2019. It is available as a PDF file (`Poster.pdf`) in the main directory of this repository, it can be downloaded from [here](http://doi.org/10.5281/zenodo.3566018) or the poster can be generated from the latex file.

## Installation

In order to reproduce the poster you will need:

* A working `R` installation with the required packages also installed. The required packages are listed at the begining of the R script `/R/utils-1.0.0.R`. Note that the package `scpdata` might not (yet) be available on CRAN or Bioconductor. The easiest way is to install it using `devtools::install_github("cvanderaa/scpdata")`.
* A working `TeX` distribution

## Compiling poster

The repository contains the required files to generate the poster. The poster is generated by compiling the `Poster.tex` file. The figures used in the poster are generated from the R script `R/utils-1.0.0.R`. 

## Reference

You can cite the poster as:

> Vanderaa, Christophe, & Gatto, Laurent. (2019). scpdata: a data package for single-cell proteomics. Zenodo. http://doi.org/10.5281/zenodo.3566018

