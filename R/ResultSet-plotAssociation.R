#' @describeIn ResultSet Allows to plot a series of plots (QQ plot, Manhattan
#' plot and Volcano plot) depending on the results stored in the
#' \code{ResultSet}.
setMethod(
    f = "plotAssociation",
    signature = "ResultSet",
    definition = function(object, rid = 1, coef = 2, contrast = 1, type, tPV, tFC) {
        ## plot.new()
        ## type = c("manhattan", "qq", "feature")

        if(missing(tPV)) tPV <- NULL
        if(missing(tFC)) tFC <- NULL

        #if(object@fun_origin %in% c("assocGE", "assocME")) {
        if(object@fun_origin == "assocES") {
            # if(sum(object@class_origin %in%  c("ExposomeSet", "ExpressionSet", "ExposomeClust")) != 2) {
            #     stop("Invalid object 'ResultSet'. Expected an object ",
            #          "obtained from 'ExposomeSet' and 'ExpressionSet'. ",
            #          "Given one from '", paste(object@class_origin,
            #                                    collapse="', '"), "'")
            # }
            return(.plot_assoc(object, rid, coef, contrast, type, tPV, tFC))
        # } else if(object@fun_origin == "assocME") {
        #     if(sum(object@class_origin %in%  c("ExposomeSet", "MethylationSet", "ExposomeClust")) != 2) {
        #         stop("Invalid object 'ResultSet'. Expected an object ",
        #              "obtained from 'ExposomeSet' and 'MethylationSet'. ",
        #              "Given one from '", paste(object@class_origin,
        #                                        collapse="', '"), "'")
        #     }
        #     .plot_assoc(object, rid, type, tPV, tFC)
        } else if(object@fun_origin == "assocSNP") {
            if(sum(object@class_origin %in%  c("ExposomeSet", "SnpSet", "ExposomeClust")) != 2) {
                stop("Invalid object 'ResultSet'. Expected an object ",
                     "obtained from 'ExposomeSet' and 'SnpSet'. ",
                     "Given one from '", paste(object@class_origin,
                                               collapse="', '"), "'")
            }
            .plot_assoc_snps(object, type, ...)
        } else if(object@fun_origin == "assocPRT") {
            if(sum(object@class_origin %in%  c("ExposomeSet", "ExpressionSet", "ExposomeClust")) != 2) {
                stop("Invalid object 'ResultSet'. Expected an object ",
                     "obtained from 'ExposomeSet' and 'ExpressionSet'. ",
                     "Given one from '", paste(object@class_origin,
                                               collapse="', '"), "'")
            }
            .plot_assoc_prot(object, rid, type, ...)
        } else {
            stop("Invalid 'object'. Value for attribue 'fun_origin' (",
                 object@fun_origin, ") not recognized.")
        }
    })
