#' @describeIn ResultSet Draws a plot depending on the integration analysis
#' stored in the \code{Resultset}.
setMethod(
    f = "plotIntegration",
    signature = "ResultSet",
    definition = function(object, cmpX=1, cmpY=2, colors, lb.th=0.20, legend.show=TRUE,
                          ...) {
        if(missing(colors)) {
            tcolors <- rainbow(length(featureData(object)))
            names(colors) <- names(object)
        }

        if(object@class_origin == "<m:mcca>") {
            .plot_integration_mcca(object, tcolors=colors, lb.th=lb.th,
                                   legend.show=legend.show, ...)
        } else if(object@class_origin == "<m:mcia>") {
            .plot_integration_mcia(object, cmpX=cmpX, cmpY=cmpY, ...)
        }
})
