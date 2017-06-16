setMethod(
    f = "crossomics",
    signature = "MultiDataSet",
    definition =  function(object, method="mcca", ncomponents=2, ..., na.rm=FALSE,
                       permute = c(100, 3), verbose=FALSE, warnings=TRUE) {
        list <- as_list_mds(object); rm( object )
        ## --------------------------------------------------------------------- ##
        ## GENERAL CHECKS
        method <- match.arg(method, choices = c("mcca", "mcia"))
        if(length(list) < 2) {
            stop("At last two different datasets are required for integration processes.")
        }
        if(is.null(names(list))) {
            names(list) <- paste("set", 1:length(list))
        }
        ## --------------------------------------------------------------------- ##

        ## --------------------------------------------------------------------- ##
        ## REDUCE DATASETS TO COMMON SAMPLES
        if(warnings | verbose) {
            warning("Sets in list will be reduced to common samples.")
        }

        sc <- Reduce(intersect, lapply(list, Biobase::sampleNames))
        list <- lapply(list, function(it) it[ , sc])
        ## --------------------------------------------------------------------- ##

        if(method == "mcca") {
            .crossomics_mcca_list(list, ncomponents=ncomponents, na.rm=na.rm,
                                  permute=permute, verbose=verbose, warnings=warnings, ...)
        } else if(method == "mcia") {
            .crossomics_mcia_list(list, verbose=verbose, warnings=warnings, ...)
        } else {
            stop("Invalid method (", method, ") was given.")
        }
    }
)
