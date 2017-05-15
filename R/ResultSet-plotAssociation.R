#' @describeIn ResultSet Allows to plot a series of plots (QQ plot, Manhattan
#' plot and Volcano plot) depending on the results stored in the
#' \code{ResultSet}.
#' @param rid Name or index of the internal result to be used
#' @param coef Coefficient to be returne, usually 2
#' @param contrast If coefficient to be used was multicategorical, number
#' of the contrast to be returned.
#' @param type Type of plot to be drawn
#' @param tPV Threshold for P-Value
#' @param tFC Threshold for log FC of effect
setMethod(
    f = "plotAssociation",
    signature = "ResultSet",
    definition = function(object, rid = 1, coef = 2, contrast = 1, type, tPV, tFC, show.effect=FALSE) {
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
            return(.plot_assoc(object, rid, coef, contrast, type, tPV, tFC, show.effect))
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
