.crossomics_mcca <- function(mds, ncomponents = 2, ..., na.rm = FALSE,
                             permute=NULL, verbose = FALSE, warnings = TRUE) {

    ## --------------------------------------------------------------------- ##
    ## CREATE LIST OF TABLES IN BASE OF THE TYPE OF DATA
    dta_list <- in_mds_for_crosomics(mds, na.rm = na.rm, verbose = verbose, warnings = warnings)
    ## --------------------------------------------------------------------- ##

    ## --------------------------------------------------------------------- ##
    ## PERFORM THE INTEGRATION WITH MCCA
    if(verbose) {
        message("Performing crossomics")
    }

    if(is.null(permute)) {
        int <- PMA::MultiCCA(dta_list, ncomponents = ncomponents, ...)
    } else {
        perm.out <- PMA::MultiCCA.permute(dta_list, nperms = permute[1],
            niter = permute[2])
        int <- PMA::MultiCCA(dta_list, ncomponents = ncomponents,
            penalty = perm.out$bestpenalties, ws=perm.out$ws.init, ...)
    }
    ## --------------------------------------------------------------------- ##

    names(fdt_list) <- names(list)
    options=list(
        N = nrow(dta_list[[1]]),
        S = length(list),
        names = names(list),
        ncomponents = ncomponents,
        na.rm = na.rm,
        permute = permute,
        method="MultiCCA",
        package="PMA"
    )

    MultiDataSet::create_resultset(
        fOrigin = "crossomics",
        lresults = list("crossomics"=list("result" = int)),
        fData = fdt_list,
        lOptions = options
    )
}
