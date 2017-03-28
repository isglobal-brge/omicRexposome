#' @describeIn ResultSet Returns a list with the options used to create the
#' \code{ResultSet}
setMethod(
    f = "options",
    signature = "ResultSet",
    definition = function(object) {
        c(fun_origin=object@fun_origin, object@options)
    }
)
