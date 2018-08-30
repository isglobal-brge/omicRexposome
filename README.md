# omicRexposome

## Summary

`omicRexposome` is an R package for extending `rexposome` capabilities and include exposome-omic data analysis and integration. It depends in a series of third party R packages to provide:

  1. A series of pipelines to test exposome-omic and diseasome-omic associations.
    * [UNDER DEVELOPMENT] Basic GWAS pipeline based on `snpStats`
    * Methylome, Transcriptome and Proteome Association Analysis based on `limma`
  2. Two different approaches to integrate exposome with omic data are implemented using *multiple co-inertia analysis* from `omicade4` and *multi canonical correlation analysis* from `PMA`

## Installation

`omicRexposome` requires R version equal or newer than 3.3.0. The following script allows to install `rexposome` dependencies:

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

packages = c('Biobase', 'methods', 'snpStats', 'limma', 'sva', 'ggplot2',
    'ggrepel', 'PMA', 'omicade4', 'ggplot2', 'qqman', 'gridExtra'
)
for( pkg in packages ) {
  if( !pkg %in% rownames( installed.packages() ) ) {
    message( "Installing ", pkg )
    BiocManager::install( pkg )
  }
}
```

The package can be installed using the R package `devtools`. `devtools` can be installed win the following code:

```r
install.packages("devtools")
```

Once `devtools` and the dependences are installed, the following code installs `omicRexposome` and the basic dependence `rexposome`:

```r
devtools::install_github("isglobal-brge/rexposome")
devtools::install_github("isglobal-brge/omicRexposome")
```

## Basic Guide

Exposome-Omic Association is done using the function `assocES`. This function requires an argument `x` being an `ExposomeSet` and an argument `y` being an `ExpressionSet` with the correct omic data (gene expression for transcriptome, betas or Ms for methylome, and protein level for proteome).

  * `plotAssociation` allows to plot the result of all _assoc*_ functions having an argument `type` that can takes:
    * `"manhattan"` to draw a typical Manhattan plot
    * `"protein"` to draw an adapted version of a Manhattan plot for protein data
    * `"volcano"` to draw a volcano plot, having the option to fill the arguments `tPV` (significant P-Value) and `tFC` (significant fold change)
    * `"qq"` to draw a standard QQ plot

Function `crossomics` allows to perform a multi-omic integration join exposome by selecting one of the available methods (`"mcia"` or `"mcca"`). The main argument, called `list`, must be filled with a list of `ExpressionSet`s (plus `ExposomeSet`s).

  * `plotIntegration` allows to plot the results of `crossomics`, having a proper visualization for each method.
