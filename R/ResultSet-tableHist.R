
setMethod(
    f = "tableHits",
    signature = "ResultSet",
    definition = function(object, th=0.05) {
        data.frame(
            exposure=rid(object),
            hits=sapply(rid(object), function(expo) {
                tt <- topTable(object, rid=expo)
                sum(tt$P.Value < th)
            })
        )
    }
)
