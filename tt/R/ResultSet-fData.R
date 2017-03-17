#' @describeIn ResultSet Returns \code{AnnotatedDaatFrame} with feature's data.
setMethod(
    f = "featureData",
    signature = "ResultSet",
    definition = function(object) {
        return(object@fData)
    }
)

#' @describeIn ResultSet Returns \code{data.frame} with feature's data.
setMethod(
    f = "fData",
    signature = "ResultSet",
    definition = function(object) {
        return(lapply(object@fData, function(x) as(x, "data.frame")))
    }
)
