#' omicRexposome: Package for exposome and omic data associatin and integration
#'
#' @section exposome-omic data association study:
#' The packages offers the function \code{\link{assocES}} that allows to perform
#' an association study using transcriptome, methylome, etc. as dependent
#' variable and exposome data as independent variable. The function relies on
#' \code{limma} pipeline and generates an object of class \code{\link{ResultSet}},
#' that can be ploted using \code{\link{plotAssociation}}.
#'
#' @section exposome-omic data integration study:
#' The packages offers the function \code{\link{crossomics}} that allows to perform
#' two types of integration study: Multi Canonical Correlation Analysis and
#' Multi Co-Inertia Analysis. The function allos to use any type and number of
#' datasets (aka. exposome transcriptome, methylome, etc.). The function generates an
#' object of class \code{\link{ResultSet}}, that can be ploted using
#' \code{\link{plotIntegration}}.
#'
#' @docType package
#' @name omicRexposome
#'
#' @import utils
#' @import methods
#'
# @importClassesFrom MultiDataSet MultiDataSet
#' @importClassesFrom rexposome ExposomeSet ExposomeClust
#'
# @importFrom MultiDataSet commonSamples betas getMs
#' @importFrom Biobase sampleNames pData fData
#' @importFrom rexposome exposureNames phenotypeNames
#' @importFrom sva num.sv sva
#' @importFrom limma lmFit eBayes
#' @importFrom ggplot2 theme theme_bw
#' @importFrom ggplot2 aes geom_polygon geom_point geom_hline geom_vline
#' @importFrom ggplot2 scale_colour_manual
#' @importFrom ggplot2 xlab ylab
#' @importFrom ggrepel geom_text_repel
#' @importFrom omicade4 mcia
#' @importFrom PMA MultiCCA MultiCCA.permute
#' @importFrom grDevices rainbow
#' @importFrom graphics plot
#' @importFrom methods as new
#' @importFrom stats as.formula model.matrix qbeta qchisq qnorm
#' @importFrom gridExtra grid.arrange
#' @importFrom qqman manhattan
#' @importFrom grDevices rainbow
#' @importFrom methods as new
#' @importFrom stats as.formula model.matrix qbeta qchisq qnorm
NULL
