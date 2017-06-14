#' @describeIn ResultSet Getter to obtain the raw \code{data.frame} from
#' association and integration analysis.
#' @param sort Indicates of the result should be ordered by P-Value
setMethod(
    f = "topTable",
    signature = "ResultSet",
    definition = function(object, rid, coef=2, sort=TRUE) {
        if(object@fun_origin == "association") {
            ff <- ifelse(object@options$eBayes, limma::topTable, limma::toptable)
            if(missing(rid)) {
                res <- lapply(names(object@results), function(nme) {
                    if(class(object@results[[nme]]$result) == "MArrayLM") {
                        tt <- ff(object@results[[nme]]$result, coef=coef, n=Inf)
                        tt$exposure <- nme
                        return(tt)
                    } else {
                        res <- data.frame(logFC=0, t=0, P.Value=0, adj.P.Val=0, B=0)
                        res[-1, ]
                    }
                })
                res <- do.call(rbind, res)
            } else {
                if(class(object@results[[rid]]$result) == "MArrayLM") {
                    res <- ff(object@results[[rid]]$result, coef=coef, n=Inf)
                } else {
                    res <- data.frame(logFC=0, t=0, P.Value=0, adj.P.Val=0, B=0)
                    res[-1, ]
                }
            }
        } else {
            stop("Invalid 'object'. Value for attribue 'fun_origin' (",
                 object@fun_origin, ") not recognized.")
        }

        return(res)
})
