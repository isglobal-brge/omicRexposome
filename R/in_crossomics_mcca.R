.crossomics_mcca <- function(object, ncomponents = 2, ..., na.rm = FALSE,
                             permute=NULL, verbose = FALSE, warnings = TRUE) {

    ## --------------------------------------------------------------------- ##
    ## CREATE LIST OF TABLES IN BASE OF THE TYPE OF DATA
    dta_list <- in_mds_for_crosomics(object,
        na.rm = na.rm, verbose = verbose, warnings = warnings)
    fdt_list <- dta_list[["fdata"]]
    dta_list <- dta_list[["adata"]]
    ## --------------------------------------------------------------------- ##

    ## --------------------------------------------------------------------- ##
    ## PERFORM THE INTEGRATION WITH MCCA
    if(verbose) {
        message("Performing crossomics (MCCA)")
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
        lResults = list("crossomics"=list("result" = int, error=NA)),
        fData = fdt_list,
        lOptions = options
    )
}
