# omicRexposome

## Summary

`omicRexposome` is an R package for extending `rexposome` capabilities and inclide exposome-omic data analysis and integration. It depends in a series of third party R packages to provide:

  1. [UNDER DEVELOPMENT] A series of pipelines to test exposome-omic and diseasome-omic associations.
    * Basic GWAS pipeline based on `snpStats`
    * Gene Expression Association Analysis based on `limma`
    * Methylation Levels Association Analysis based on `MEAL`
    * Proteine Level Association Analysis based on `limma`
  2. [UNDER DEVELOPMENT] Two different approaches to integrate exposome with omic data.

## Installation

`omicRexposome` requires R version equal or newer than 3.3.0. The following script allows to install `rexposome` dependencies:

```r
source( "http://bioconductor.org/biocLite.R" )

packages = c('Biobase', 'MultiDataSet', 'BiocInstaller', 'gridExtra', 'stringr', 
    'ggplot2', 'reshape2', 'snpStats', 'MEAL', 'limma', 'sva', 'glmnet', 
    'omicade4', 'ggrepel', 'PMA', 'qqman'
)
for( pkg in packages ) {
  if( !pkg %in% rownames( installed.packages() ) ) {
    message( "Installing ", pkg )
    biocLite( pkg )
  }
}
```

The package can be installed using the R package `devtools`. `devtools` can be installed win the following code:

```r
install.packages("devtools")
```

Once `devtools` and the dependences are installed, the following code installs `omicRexposome` and the basi dependence `rexposome`:

```r
devtools::install_github("isglobal-brge/rexposome")
devtools::install_github("isglobal-brge/omicRexposome")
```

### Details

## Basic Guide

### Loading Exposome

### Loading Omics

### Exposome-Omic Association

  * Function `assocSNP` allows to perform a basic GWAS.
  * Function `assocGE` allows to perform a DE analysis in base of the exposures.
  * Function `assocME` allows to perform an EWAS for each exposure.
  * Function `assocPRT` allows to test the association of the protein quantification with the exposures levels.
  * `plotAssociation` allows to plot the result of all _assoc*_ functions.

### Exposome-Omic Integration

  * Function `crossomics` allows to perform a multi-omic integration join exposome by selecting one of the available methods (`"mcia"` or `"mcca"`).
  * `plotIntegration` allows to plot the results of `crossomics`.
