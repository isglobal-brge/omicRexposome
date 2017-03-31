#' @describeIn ResultSet Returns a table with the lambda value of each
#' analyses stored in the \code{ResultSet}.
setMethod(
    f = "tableLambda",
    signature = "ResultSet",
    definition = function(object, trim=0.5) {
        data.frame(
            "exposure"=unique(rid(object)),
            "lambda"=sapply(rid(object), function(expo) {
                lambdaClayton(topTable(object, rid=expo)$P.Value, trim=trim)
                #qchisq(median(extract(object, rid=expo)$P.Value), df=2, lower.tail=FALSE)
            })
        )
    }
)
