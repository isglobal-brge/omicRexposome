.crossomics_mcca_list <- function(list, ncomponents = 2, ..., na.rm = FALSE,
                             permute=NULL, verbose = FALSE, warnings = TRUE) {

    ## --------------------------------------------------------------------- ##
    ## CREATE LIST OF TABLES IN BASE OF THE TYPE OF DATA
    dta_list <- list()
    fdt_list <- list()

    ii <- 1
    for(set in list) {
        if(class(set) == "SnpSet")  {
            if(verbose) {
                message("Creating table for SnpSet.")
            }

            dta <- t(snpToContinuous(set, verbose))
            rownames(dta) <- Biobase::sampleNames(set)
            colnames(dta) <- Biobase::featureData(set)
            fdt <- Biobase::fData(set)[colnames(dta), ]

            if(sum(is.na(dta)) != 0 & !na.rm) {
                stop("Table SnpSet contains NA values.")
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

                fdt <- Biobase::fData(set)[colnames(dta), ]
            }

            dta_list[[ii]] <- dta
            fdt_list[[ii]] <- fdt

        } else if (class(set) == "ExpressionSet") {
            if(verbose) {
                message("Creating table for ExpressionSet.")
            }
            dta <- t(Biobase::exprs(set))
            fdt <- Biobase::fData(set)[colnames(dta), ]

            if(sum(is.na(dta)) != 0 & !na.rm) {
                stop("Table ExpressionSet contains NA values.")
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

                fdt <- Biobase::fData(set)[colnames(dta), ]
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
            sel <- rexposome::exposureNames(set)[Biobase::fData(set)[ , '_type', drop=FALSE] == "numeric"]
            dta <- dta[ , sel, drop=FALSE]
            fdt <- Biobase::fData(set)[colnames(dta), ]

            if(sum(is.na(dta)) != 0 & !na.rm) {
                stop("Table ExposomeSet contains NA values.")
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

                fdt <- Biobase::fData(set)[colnames(dta), ]
            }

            dta_list[[ii]] <- dta
            fdt_list[[ii]] <- fdt

        }

        ii <- ii + 1
    }
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
    new("ResultSet",
        fun_origin = "crossomics",
        results = list("crossomics"=list("result" = int)),
        fData = fdt_list,
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
    )
}
