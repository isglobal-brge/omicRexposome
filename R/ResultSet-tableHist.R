#' @describeIn ResultSet Returns a table with the number of hits in each
#' result of the \code{ResultSet}.
setMethod(
    f = "tableHits",
    signature = "ResultSet",
    definition = function(object, th=0.05) {
        data.frame(
            exposure=rid(object),
            hits=sapply(rid(object), function(expo) {
                tt <- MultiDataSet::topTable(object, rid=expo)
                sum(tt$P.Value < th)
            })
        )
    }
)
