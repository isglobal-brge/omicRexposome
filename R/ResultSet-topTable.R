#' @describeIn ResultSet Getter to obtain the raw \code{data.frame} from
#' association and integration analysis.
#' @param sort Indicates of the result should be ordered by P-Value
setMethod(
    f = "topTable",
    signature = "ResultSet",
    definition = function(object, rid, coef=2, contrast=1, sort=TRUE) {
        if(object@fun_origin == "assocES") {
            ff <- ifelse(object@options$eBayes, limma::topTable, limma::toptable)
            if(missing(rid)) {
                res <- lapply(names(object@results), function(nme) {
                    if(class(object@results[[nme]]$result) == "MArrayLM") {
                        message(nme)
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
                    if(class(res) == "list") {
                        res <- res[contrast]
                    }
                } else {
                    res <- data.frame(logFC=0, t=0, P.Value=0, adj.P.Val=0, B=0)
                    res[-1, ]
                }
            }
            #res <- cbind(res, object@fData[[2]][rownames(res), ])
        } else if(object@fun_origin == "assocSNP") {
            if(!missing(rid)) {
                warning("Given 'rid'. Invalid argument for assocSNP result.")
            }
            res <- object@results[[1]]$result
            res <- res[order(res$PValHWE), ]
        } else if(object@fun_origin == "crossomics") {
            res <- data.frame(do.call(rbind, lapply(1:length(names(object)), function(ii) {
                tbl <- cbind(object@results[[1]][[1]]$ws[[ii]], names(object)[ii])
                rownames(tbl) <- rownames(object@fData[[ii]])
                colnames(tbl) <- c("x", "y", "feature")
                tbl[ tbl[ , 1] != 0 | tbl[ , 2] != 0, ]
            })), stringsAsFactors = FALSE)
            res$x <- as.numeric(res$x)
            res$y <- as.numeric(res$y)
        } else if(object@fun_origin == "assocPRT") {
            if(missing(rid)) {
                res <- lapply(names(object@results), function(nme) {
                    tt <- object@results[[nme]]$result
                    tt$exposure <- nme
                    tt
                })
                res <- do.call(rbind, res)
            } else {
                res <- object@results[[rid]]$result
            }
        } else {
            stop("Invalid 'object'. Value for attribue 'fun_origin' (",
                 object@fun_origin, ") not recognized.")
        }

        return(res)
})
