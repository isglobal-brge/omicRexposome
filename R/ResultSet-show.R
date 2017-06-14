setMethod(
    f = "show",
    signature="ResultSet",
    definition = function(object) {
        cat("Object of class 'ResultSet'\n", sep="")
        cat(" . created with:", object@fun_origin, "\n")
        if(object@fun_origin == "association") {
            cat("    . via:", paste0(rst_2@names, collapse=" and "), "\n")
        }
        cat(" . sva: ", ifelse(sum(is.na(sapply(object@results, "[[", "sva.num")) != 0 ), "no", "yes"), "\n")
        if(object@fun_origin == "crossomics") {
            cat("    . method:", object@options$method, " (", object@options$package, ")\n")
        }
        cat(" . #results:", length(object@results), "( error:", sum(!is.na(sapply(object@results, "[[", "error"))), ")\n")
        cat(" . featureData: ", length(object@fData), "\n")
        for(nm in names(object@fData)) {
            nr <- nrow(object@fData[[nm]])
            nc <- ncol(object@fData[[nm]])
            cat("    . ", nm, ": ", nr, "x", nc, "\n", sep="")
        }
    }
)
