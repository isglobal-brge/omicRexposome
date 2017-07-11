setMethod(
    f = "plotIntegration",
    signature = "ResultSet",
    definition = function(object, cmpX=1, cmpY=2, lb.th=0.20, legend.show=TRUE,
                          colors, ...) {
        if(missing(colors)) {
            colors <- rainbow(length(Biobase::fData(object)))
            names(colors) <- names(object)
        } else {
            if(sum(names(object@fData) %in% names(colors)) != length(colors)) {
                stop("Invalid colors given. Colors must be a vector with names",
                     " equal to original list of data-sets.")
            } else {
                colors <- colors[names(object@fData)]
            }
        }

        if(object@options$method == "MultiCCA") {
            .plot_integration_mcca(object, tcolors=colors, lb.th=lb.th,
                                   legend.show=legend.show, ...)
        } else if(object@options$method == "mcia") {
            .plot_integration_mcia(object,
                                   cmpX=cmpX, cmpY=cmpY, tcolors=colors)
        }
})
