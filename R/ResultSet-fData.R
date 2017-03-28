#' @describeIn ResultSet Returns \code{data.frame} with feature's data.
setMethod(
    f = "fData",
    signature = "ResultSet",
    definition = function(object) {
        return(lapply(object@fData, function(x) as(x, "data.frame")))
    }
)
