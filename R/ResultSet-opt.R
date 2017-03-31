#' @describeIn ResultSet Returns a list with the options used to create the
#' \code{ResultSet}
setMethod(
    f = "opt",
    signature = "ResultSet",
    definition = function(object) {
        c(fun_origin=object@fun_origin, object@options)
    }
)
