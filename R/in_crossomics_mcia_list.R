.crossomics_mcia_list <- function(list, ncomponents = 2, ...,
                             verbose = FALSE, warnings = TRUE) {

    list2 <<- list
    mres <- omicade4::mcia(list, cia.nf = ncomponents, ...)
    new("ResultSet",
        fun_origin = "crossomics",
        results = list("crossomics"=list("result" = mres)),
        fData = lapply(list, Biobase::fData),
        options = list(
            ncomponents = ncomponents,
            method = "mcia",
            package = "omicade4"
        )
    )
}
