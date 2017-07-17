#' @aliases getIntegration
#' @rdname getIntegration-methods
setMethod(
    f = "getIntegration",
    signature = "ResultSet",
    definition = function(object, ...) {
        object@results[[1]]$result
    }
)
