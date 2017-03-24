#' @describeIn ResultSet Returns the names of the omics data used to create
#' the \code{ResultSet}.
setMethod(
    f = "names",
    signature="ResultSet",
    definition = function(x) {
        return(names(x@fData))
    }
)
