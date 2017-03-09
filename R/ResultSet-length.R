#' @describeIn ResultSet Returns the amoung of analyses stored in the
#' \code{ResultSet}.
setMethod(
    f = "length",
    signature="ResultSet",
    definition = function(x) {
        return(length(x@results))
    }
)
