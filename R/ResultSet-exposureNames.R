#' @describeIn ResultSet Getter to obtain the exposures's names.
#' @importMethodsFrom rexposome exposureNames
setMethod(
    f = "exposureNames",
    signature="ResultSet",
    definition = function(object) {
        return(names(object@results))
    }
)
