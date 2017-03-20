setMethod(
    f = "assocME",
    signature = "MultiDataSet",
    definition = function(object, formula, select, set="exposures",
                          area.test=FALSE, method="ls", betas=TRUE, ...,
                          sva=FALSE, vfilter=NULL, ncores=1, verbose=FALSE,
                          warnings=TRUE) {
        ## ----------------------------------------------------------------- ##
        ## CHEKS
        ## ----------------------------------------------------------------- ##
        set <- match.arg(set, choices=c("exposures", "phenotypes"))

        tomic <- grep("methylation", names(object))
        texp <- grep("exposures", names(object))
        tcls <- grep("cluster", names(object))
        tomic <- names(object)[tomic]

        if(length(texp) != 0 & length(tcls) != 0) {
            stop("Given 'MultiDataSet' contains both 'ExposomeSet' and 'ExposomeClust'.")
        }

        tann <- ifelse(length(texp) == 0, "cluster", "exposure")
        texp <- ifelse(length(texp) == 0, names(object)[tcls], names(object)[texp])
        rm(tcls)

        if(length(tomic) != 1) {
            stop("One of the tables must exists in 'MultiDataSet' as ",
                "'methylation'")
        }
        if(length(texp) != 1) {
            stop("One of the tables must exists in 'MultiDataSet' as ",
                "'exposures'")
        }
        if(verbose) {
            message("The following tables will be used in the association ",
                "process: ", paste0("'", paste(c(tomic, texp), collapse="', '"),
            "'"))
        }

        if(warnings | verbose) {
            warning("Sets from 'MultiDataSet' will be reduced to common samples")
        }

        l1 <- sapply(Biobase::sampleNames(object)[c(tomic, texp)], length)
        object <- MultiDataSet::commonSamples(object)
        l2 <- sapply(Biobase::sampleNames(object)[c(tomic, texp)], length)
        l3 <- mapply('-', l1, l2, SIMPLIFY = FALSE)

        if(verbose) {
            message(paste(unlist(l3), names(l3),
                sep = " samples were reduced from ", collapse = ", "))
        }
        ## ----------------------------------------------------------------- ##

        ## ----------------------------------------------------------------- ##
        ## CONVERT FORMULA
        ## ----------------------------------------------------------------- ##
        formula <- as.character(as.formula(formula))
        es <- object[[texp]]
        exp.dt <- data.frame(pData(es), expos(es))
        fd <- object@featureData[[tomic]]@data
        rm(es)

        if(tann == "cluster") {
            ## ------------------------------------------------------------- ##
            ## EXPOSOME CLUSTER ANALYSIS
            ## ------------------------------------------------------------- ##
            if(warnings | verbose) {
                warning("Given 'MultiDataSet' contains an 'ExposomeClust'. Association test will be performed on clustering result.")
            }
            if(!missing(select) & (warnings | verbose)) {
                warning("Given values in 'select' argument. They will be dropped.")
            }
            select <- "cluster"
        } else {
            ## ------------------------------------------------------------- ##
            ## EXPOSOME SET ANALYSIS
            ## ------------------------------------------------------------- ##
            if(missing(select)) {
                if(set == "exposures") {
                    warning("No given 'select'. association will be computed for all exposures")
                    select <- rexposome::exposureNames(object[[texp]])
                } else { ## set == "phenotypes"
                    warning("No given 'select'. association will be computed for all phenotypes")
                    select <- Biobase::phenotypeNames(object[[texp]])
                }
            }
        }
        ## ----------------------------------------------------------------- ##


        ## ----------------------------------------------------------------- ##
        ## EXTRACT SETS AND PERFORM ANALYSIS
        ## ----------------------------------------------------------------- ##
        results <- lapply(select, function(ex) {
            design <- as.formula(paste0("~", ex, "+", formula[2]))
            if(verbose) {
                message("Evaluating model '", as.character(design), "'.")
            }
            pheno <- .create_p(
                expo.dt = exp.dt,
                omic.p = pData(object[[tomic]]),
                select = all.vars(design)
            )

            if(class(pheno) == "character") {
                stop("Invalid value '", pheno, "' in 'exposures' or 'covariates'")
            }

            na.loc <- rowSums(apply(pheno, 2, is.na))
            na.loc <- which(na.loc != 0)
            if(length(na.loc) != 0){
                warning("There are missing values. ", length(na.loc),
                        " samples will be removed.")
                pheno <- pheno[-na.loc, , drop=FALSE]
            }


            if(sum(!sapply(sapply(apply(pheno, 2, table), length), ">", 1))) {
                warning("When testing for '", ex, "', at last one covariate ",
                    "is constant")
                return(list(
                    N=NA,
                    sva.num=NA,
                    design=NA,
                    result=NULL,
                    error=paste0("Covariate in model '", as.character(design) ,"' is constant")
                ))
            } else {
                tryCatch({
                    # Design model for association to classification
                    design.mm <- model.matrix(formula(design), data = pheno)

                    # Extract betas or Ms
                    if(betas) {
                        methy <- MultiDataSet::betas(object[[tomic]][ , rownames(pheno)])
                    } else {
                        methy <- MultiDataSet::getMs(object[[tomic]][ , rownames(pheno)])
                    }

                    # If required, apply SVA
                    if(sva) {
                        if (verbose | warnings){
                            if(is.null(vfilter)) {
                                message("Computing SVA. This step can be very time consuming.",
                                        "Try using argument 'vfilter'.")
                            } else {
                                message("Computing SVA. This step can be very time consuming.")
                            }
                        }
                        n.sv <- sva::num.sv(methy, design.mm, vfilter=vfilter)
                        if (n.sv > 0){
                            svobj <- sva::sva(methy, design.mm,
                                              design.mm[ , -1, drop=FALSE],
                                              n.sv=n.sv, vfilter=vfilter)
                            design.mm <- cbind(design.mm, svobj$sv)
                        }
                        rm(svobj)
                        suppressMessages(gc())
                    }

                    # Fit the model
                    if(area.test) {
                        result <- MEAL::DARegion(set= methy, model=design.mm,
                            methods=c("blockFinder", "bumphunter", "DMRcate"),
                            ...)
                    } else {
                        result <- MEAL::DAProbe(set=methy, model=design.mm,
                            method=method, ...)
                    }
                    # -----------------------------------------------------
                    suppressMessages(gc())


                    result <- cbind(result, fd[rownames(result), ])
                    list(
                        N=nrow(pheno),
                        sva.num=n.sv,
                        design=design,
                        result=result,
                        error=NA
                    )
                }, error=function(e){
                    return(list(
                        N=NA,
                        sva.num=NA,
                        design=NA,
                        result=NULL,
                        error=e
                    ))
                })
            }
        })
        names(results) <- select

        class_origin = c(
            ifelse(tann == "cluster", "ExposomeClust", "ExposomeSet"),
            "MethylationSet")

        new("ResultSet",
            fun_origin = "assocME",
            class_origin = class_origin,
            names = c(texp, tomic),
            results = results,
            fData = fData(object)[c(texp, tomic)],
            sva=ifelse(sva, 1, 0)
        )
})
