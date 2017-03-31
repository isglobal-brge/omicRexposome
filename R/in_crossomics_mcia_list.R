.crossomics_mcia_list <- function(list, ncomponents = 2, ...,
                             verbose = FALSE, warnings = TRUE) {

    ## --------------------------------------------------------------------- ##
    ## CREATE LIST OF TABLES IN BASE OF THE TYPE OF DATA
    dta_list <- list()
    fdt_list <- list()
    na.rm = FALSE

    ii <- 1
    for(set in list) {
        if(class(set) == "SnpSet")  {
            if(verbose) {
                message("Creating table for SnpSet.")
            }

            dta <- t(snpToContinuous(set, verbose))
            rownames(dta) <- Biobase::sampleNames(set)
            colnames(dta) <- Biobase::featureData(set)
            fdt <- fData(set)[colnames(dta), ]

            if(sum(is.na(dta)) != 0 & !na.rm) {
                 #stop("Table SnpSet contains NA values.")
            } else if(sum(is.na(dta)) != 0 & na.rm) {
                if(verbose) {
                    message("Removing SNPs with missing data.")
                }
                onc <- ncol(dta)
                dta <- dta[ , colSums(dta) == 0]
                if(verbose | warnings) {
                    r <- onc - ncol(dta)
                    warning("Removed ", r, " (", round(r/onc * 100, 2), "%) SNPs from original SnpSet.")
                }
                rm(onc)

                fdt <- fData(set)[colnames(dta), ]
            }

            dta_list[[ii]] <- dta
            fdt_list[[ii]] <- fdt

        } else if (class(set) == "ExpressionSet") {
            if(verbose) {
                message("Creating table for ExpressionSet.")
            }
            dta <- t(Biobase::exprs(set))
            fdt <- fData(set)[colnames(dta), ]

            if(sum(is.na(dta)) != 0 & !na.rm) {
                #stop("Table ExpressionSet contains NA values.")
            } else if(sum(is.na(dta)) != 0 & na.rm) {
                if(verbose) {
                    message("Removing probes with missing data.")
                }
                onc <- ncol(dta)
                dta <- dta[ , colSums(dta) == 0]
                if(verbose | warnings) {
                    r <- onc - ncol(dta)
                    warning("Removed ", r, " (", round(r/onc * 100, 2), "%) probes from original ExpressionSet.")
                }
                rm(onc)

                fdt <- fData(set)[colnames(dta), ]
            }

            dta_list[[ii]] <- dta
            fdt_list[[ii]] <- fdt

        } else if (class(set) == "ExposomeSet") {
            if(verbose) {
                message("Creating table for ExposomeSet")
            }
            dta <- rexposome::expos(set)
            if(verbose | warnings) {
                warning("Factor exposures will be discarded.")
            }
            sel <- exposureNames(set)[fData(set)[ , '_type', drop=FALSE] == "numeric"]
            dta <- dta[ , sel, drop=FALSE]
            fdt <- fData(set)[colnames(dta), ]

            if(sum(is.na(dta)) != 0 & !na.rm) {
                #stop("Table ExposomeSet contains NA values.")
            } else if(sum(is.na(dta)) != 0 & na.rm) {
                if(verbose) {
                    message("Removing exposures with missing data.")
                }
                onc <- ncol(dta)
                dta <- dta[ , colSums(dta) == 0]
                if(verbose | warnings) {
                    r <- onc - ncol(dta)
                    warning("Removed ", r, " (", round(r/onc * 100, 2), "%) exposures from original ExpressionSet.")
                }
                rm(onc)

                fdt <- fData(set)[colnames(dta), ]
            }

            dta_list[[ii]] <- dta
            fdt_list[[ii]] <- fdt

        }

        ii <- ii + 1
    }
    ## --------------------------------------------------------------------- ##

    dta_list <- lapply(dta_list, t)

    mres <- omicade4::mcia(dta_list, cia.nf = ncomponents, ...)
    new("ResultSet",
        fun_origin = "crossomics",
        results = list("crossomics"=list("result" = mres)),
        fData = fdt_list,
        options = list(
            N = ncol(dta_list[[1]]),
            S = length(list),
            names = names(list),
            ncomponents = ncomponents,
            method = "mcia",
            package = "omicade4"
        )
    )
}
