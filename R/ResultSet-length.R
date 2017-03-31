#' @describeIn ResultSet Returns the amoung of analyses stored in the
#' \code{ResultSet}.
#' @param x A \code{ResultSet} object.
setMethod(
    f = "length",
    signature="ResultSet",
    definition = function(x) {
        return(length(x@results))
    }
)
