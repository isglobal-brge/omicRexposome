setMethod(
    f = "association",
    signature = "MultiDataSet",
    definition = function(object, formula, expset, omicset, set = "exposures",
            method = "ls", ..., sva = FALSE, ebayes = TRUE, verbose = FALSE,
            warnings = TRUE) {
        ## CHEKS
        ## --------------------------------------------------------------------
        if(missing(expset) | missing(omicset)) {
            stop("Arguments 'expset' and 'omicset' must be fileld selecting ",
                " an 'ExposomeSet' and an omic-set.")
        }
        if(length(expset) > 1) {
            warning("Given 'expset' contains more than one element. Only ",
                    "first element will be used.")
            expset <- expset[1]
        }
        if(length(omicset) > 1) {
            warning("Given 'omicset' contains more than one element. Only ",
                    "first element will be used.")
            omicset <- omicset[1]
        }
        if(!expset %in% names(object)) {
            stop("Given 'expset' not in MultiDayaSet.")
        }
        if(!omicset %in% names(object)) {
            stop("Given 'omicset' not in MultiDayaSet.")
        }
        set <- match.arg(set, choices = c("exposures", "phenotypes"))
        if(length(object) < 2)  {
            stop("Given 'MultiDataSet' must contains an 'ExposomeSet' and ",
                 "another omic data-set.")
        }
        ## --------------------------------------------------------------------

        ## REDUCTION
        ## --------------------------------------------------------------------
        if(warnings | verbose) {
            warning("Sets from 'MultiDataSet' will be reduced to common samples")
        }

        l1 <- sapply(Biobase::sampleNames(object)[c(omicset, expset)], length)
        object <- MultiDataSet::commonSamples(object)
        l2 <- sapply(Biobase::sampleNames(object)[c(omicset, expset)], length)
        l3 <- mapply('-', l1, l2, SIMPLIFY = FALSE)

        if(verbose) {
            message(paste(unlist(l3), names(l3),
                          sep = " samples were reduced from ", collapse = ", "))
        }
        ## --------------------------------------------------------------------


        ## CONVERT FORMULA
        ## --------------------------------------------------------------------
        formula <- as.character(as.formula(formula))
        es <- object[[expset]]
        exp.dt <- data.frame(pData(es), rexposome::expos(es),
            pData(object[[omicset]]))
        rm(es)

        if(length(grep("cluster", names(mds))) == 1) {
            ## EXPOSOME CLUSTER ANALYSIS
            ## ----------------------------------------------------------------
            if(warnings | verbose) {
                warning("Association test will be performed on clustering result.")
            }
            select <- "cluster"
        } else {
            ## EXPOSOME SET ANALYSIS
            ## ----------------------------------------------------------------
            if(set == "exposures") {
                select <- rexposome::exposureNames(object[[expset]])
            } else { ## set == "phenotypes"
                select <- Biobase::phenotypeNames(object[[expset]])
            }
        }

        omic <- as_list_mds(object[ , omicset])[[1]]
        ## --------------------------------------------------------------------


        ## EXTRACT SETS AND PERFORM ANALYSIS
        ## --------------------------------------------------------------------
        results <- lapply(select, function(ex) {
            design <- as.formula(paste0("~", ex, "+", formula[2]))
            if(verbose) {
                message("Evaluating model '", as.character(design), "'.")
            }

            na.loc <- rowSums(apply(exp.dt, 2, is.na))
            na.loc <- which(na.loc != 0)
            if(length(na.loc) != 0) {
                if(warnings | verbose) {
                    warning("There are missing values. ", length(na.loc),
                            " samples will be removed.")
                }
                exp.dt <- exp.dt[-na.loc, , drop=FALSE]
            }
            omic <- omic[ , rownames(exp.dt), drop = FALSE]

            tbl <- sapply(apply(exp.dt[ , all.vars(design)], 2, table), length)
            if(sum(!sapply(tbl, ">", 1)) != 0) {
                warning("When testing for '", ex, "', at last one covariate ",
                        "is constant (",
                        paste(paste(names(tbl), tbl, sep=": "), collapse=", "),
                        ")")
                list(
                    N=NA,
                    sva.num=NA,
                    ebayes=ebayes,
                    design=NA,
                    result=NA,
                    error=paste0("Covariate in model '", as.character(design) ,"' is constant")
                )
            } else {
                # Design model
                design.mm <- model.matrix(formula(design), data = exp.dt)
                # If required, apply SVA
                n.sv <- NA
                if(sva) {
                    ## Determine number of surrogate variables
                    Y.r <- t(stats::resid(stats::lm(t(omic) ~ exp.dt[ , ex])))
                    #,data=exp.dt)))
                    n.sv <- isva::EstDimRMT(Y.r, FALSE)$dim + 1

                    if (n.sv > 0) {
                        sv.obj <- SmartSVA::smartsva.cpp(omic,
                            design.mm, mod0 = NULL, n.sv = n.sv)
                            # mod0=design.mm[ , -1, drop=FALSE]
                        design.mm <- cbind(design.mm, sv.obj$sv)
                    }
                }

                # Fit the model
                if (verbose){
                    message("Fitting the model.")
                }

                fit <- limma::lmFit(omic, design.mm, method=method, ...)
                if(ebayes) {
                    fit <- limma::eBayes(fit)
                }

                list(
                    N=nrow(exp.dt),
                    sva.num=n.sv,
                    ebayes=ebayes,
                    error=NA,
                    design=design,
                    result=fit
                )
            }
        })
        names(results) <- select

        class_origin = c(
            ifelse(select[1] == "cluster", "ExposomeClust", "ExposomeSet"),
            "ExpressionSet")
        new("ResultSet",
            fun_origin = "association",
            results = results,
            fData = fData(object)[c(expset, omicset)],
            options = list(
                method="assocES",
                sva=sva,
                eBayes=ebayes,
                class_origin=class_origin,
                names = c(expset, omicset)
            )
        )
        ## --------------------------------------------------------------------
    }
)


as_list_mds <- function(x) {
    ll <- lapply(names(x), function(dtype) {
        elm <- assayDataElementNames(assayData(x)[[dtype]])[1]
        assayDataElement(assayData(x)[[dtype]], elm)
    })
    names(ll) <- names(x)
    return(ll)
}
