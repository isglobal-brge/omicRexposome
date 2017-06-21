setMethod(
    f = "crossomics",
    signature = "MultiDataSet",
    definition =  function(object, method="mcca", ncomponents=2, ..., na.rm=FALSE,
                       permute = c(100, 3), verbose=FALSE, warnings=TRUE) {
        ## --------------------------------------------------------------------- ##
        ## GENERAL CHECKS
        method <- match.arg(method, choices = c("mcca", "mcia"))
        if(length(object) < 2) {
            stop("At last two different datasets are required for integration processes.")
        }
        ## --------------------------------------------------------------------- ##

        ## --------------------------------------------------------------------- ##
        ## REDUCE DATASETS TO COMMON SAMPLES
        if(warnings | verbose) {
            warning("Sets from 'MultiDataSet' will be reduced to common samples")
        }

        l1 <- vapply(Biobase::sampleNames(object), length, FUN.VALUE = numeric(1))
        object <- MultiDataSet::commonSamples(object)
        l2 <- sapply(Biobase::sampleNames(object), length)
        l3 <- mapply('-', l1, l2, SIMPLIFY = FALSE)

        if(verbose) {
            message(paste(unlist(l3), names(l3),
                          sep = " samples were reduced from ", collapse = ", "))
        }
        ## --------------------------------------------------------------------- ##

        if(method == "mcca") {
            .crossomics_mcca_list(mds, ncomponents=ncomponents, na.rm=na.rm,
                permute=permute, verbose=verbose, warnings=warnings, ...)
        } else if(method == "mcia") {
            .crossomics_mcia_list(mds, ncomponents=ncomponents,
                verbose=verbose, warnings=warnings, ...)
        } else {
            stop("Invalid method (", method, ") was given.")
        }
    }
)
