#' @describeIn ResultSet Draws a plot depending on the integration analysis
#' stored in the \code{Resultset}.
#' @param cmpX Component to be drawn on X-axis
#' @param cmpY Component to be drawn on Y-axis
#' @param lb.th Threshold to considere to place labels or not
#' @param legend.show Indicates if legend should be placed or not
#' @param colors Labeled vector for colors of the features
#' @param ... Other arguments passed to internal methods
setMethod(
    f = "plotIntegration",
    signature = "ResultSet",
    definition = function(object, cmpX=1, cmpY=2, lb.th=0.20, legend.show=TRUE,
                          colors, ...) {
        if(missing(colors)) {
            colors <- rainbow(length(Biobase::fData(object)))
            names(colors) <- names(object)
        }

        if(object@options$method == "MultiCCA") {
            .plot_integration_mcca(object, tcolors=colors, lb.th=lb.th,
                                   legend.show=legend.show, ...)
        } else if(object@options$method == "mcia") {
            .plot_integration_mcia(object,
                                   cmpX=cmpX, cmpY=cmpY, tcolors=colors)
        }
})
