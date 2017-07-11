setMethod(
    f = "tableHits",
    signature = "ResultSet",
    definition = function(object, th=0.05) {
        data.frame(
            exposure=names(object),
            hits=sapply(names(object), function(expo) {
                tt <- MultiDataSet::getAssociation(object, rid=expo, fNames=NULL)
                sum(tt$P.Value < th)
            }),
            stringsAsFactors = FALSE
        )
    }
)
