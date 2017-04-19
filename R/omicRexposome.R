#' omicRexposome: Package for
#'
#' @section
#'
#' @docType package
#' @name rexposome
#'
#' @import utils
#'
#' @importClassesFrom MultiDataSet MultiDataSet
#' @importClassFrom rexposome ExposomeSet ExposomeClust
#'
#' @importFrom MultiDataSet commonSamples betas getMs
#' @importFrom Biobase sampleNames phenotypeNames pData fData
#' @importFrom rexposome exposureNames
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
