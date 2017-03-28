#' @describeIn ResultSet Draws a plot depending on the integration analysis
#' stored in the \code{Resultset}.
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
            .plot_integration_mcia(object, cmpX=cmpX, cmpY=cmpY, ...)
        }
})
