#' @describeIn ResultSet Returns \code{data.frame} with feature's data.
#' @param object A \code{ResultSet} object.
setMethod(
    f = "fData",
    signature = "ResultSet",
    definition = function(object) {
        return(lapply(object@fData, function(x) as(x, "data.frame")))
    }
)
