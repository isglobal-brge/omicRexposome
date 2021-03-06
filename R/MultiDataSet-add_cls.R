#' @aliases add_cls
#' @rdname add_cls-methods
setMethod(
    f = "add_cls",
    signature = c("MultiDataSet", "ExposomeClust"),
    definition = function(object, clsSet, ...) {
        if("cluster" %in% colnames(pData(clsSet))) {
            object <- MultiDataSet::add_eset(object, clsSet,
                                             dataset.type = "cluster",
                                             GRanges = NA, ...)
        } else {
            stop("Column 'cluster' must be in cls'Set's phenotype data.")
        }
        return(object)
    }
)
