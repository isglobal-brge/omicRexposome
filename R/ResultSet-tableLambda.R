#' @aliases tableLambda
#' @rdname tableLambda-methods
setMethod(
    f = "tableLambda",
    signature = "ResultSet",
    definition = function(object, trim=0.5) {
        data.frame(
            "exposure"=unique(names(object)),
            "lambda"=sapply(names(object), function(expo) {
                MultiDataSet::lambdaClayton(
                    MultiDataSet::getAssociation(object, rid=expo, fNames=NULL)$P.Value, trim=trim)
                #qchisq(median(extract(object, rid=expo)$P.Value), df=2, lower.tail=FALSE)
            }),
            stringsAsFactors = FALSE
        )
    }
)
