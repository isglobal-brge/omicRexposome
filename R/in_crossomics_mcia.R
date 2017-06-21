.crossomics_mcia <- function(mds, ncomponents = 2, ...,
                             verbose = FALSE, warnings = TRUE) {

    ## --------------------------------------------------------------------- ##
    ## CREATE LIST OF TABLES IN BASE OF THE TYPE OF DATA
    dta_list <- in_mds_for_crosomics(mds, na.rm = FALSE, verbose = verbose, warnings = warnings)
    fdt_list <- dta_list[["fdata"]]
    dta_list <- dta_list[["adata"]]

    dta_list <- lapply(dta_list, t)
    ## --------------------------------------------------------------------- ##

    ## --------------------------------------------------------------------- ##
    ## PERFORM THE INTEGRATION WITH MCIA
    mres <- omicade4::mcia(dta_list, cia.nf = ncomponents, ...)
    ## --------------------------------------------------------------------- ##

    names(fdt_list) <- names(list)
    options = list(
        N = ncol(dta_list[[1]]),
        S = length(list),
        names = names(list),
        ncomponents = ncomponents,
        method = "mcia",
        package = "omicade4"
    )

    MultiDataSet::create_resultset(
        fOrigin = "crossomics",
        lresults = list("crossomics"=list("result" = mres)),
        fData = fdt_list,
        lOptions = options
    )
}
