#' Function to perform a Transcriptome-Wide Association Study
#'
#' This function allows to perform a Transcriptome-Wide Association Study
#' by using an \code{ExposomeSet} and an \code{ExpressionSet}. It
#' allows to perform an adjustment using Surrogate Variable Analysis (from
#' R package \code{sva}).
#'
#' @param x An \code{\link{ExposomeSet}} object or an
#' \code{\link{ExposomeClust}} object.
#' @param y An \code{\link{ExpressionSet}} object.
#' @param select (optional) Character vector of exposures or fenotypes (see
#' argument \code{set}). If given, only this exposusures or phenotypes will be
#' tested, otherwise all the exposures or phenotypes in \code{x}  will be
#' tested.
#' @param set (default \code{"exposures"}) Set to be used for the association
#' analysis.
#' @param sva (default \code{FALSE}) If set to true, Surrogate Variable
#' Analysis will be done on transcriptome level in each test and the model
#' will be corrected by those surrogate variables.
#' @param vfilter (defaul \code{NULL}) Number of probes from transcriptome
#' to perform the Surrogate Variable Analysis. If \code{NULL}, all the probes
#' are used.
#' @param eBayes (default \code{TRUE}) If \code{TRUE} applies \code{eBayes}
#' after \code{lmFit}.
#' @param verbose (default \code{FALSE}) If set to \code{TRUE}, a series of
#' messages descriving the process are shown.
#' @param warnings (default \code{TRUE}) If set to \code{TRUE}, a series of
#' warnings are shown when required user atention.
#' @return An object of class \code{\link{ResultSet}}.
#' @examples
#' data(prot_r)
#' data(exp_r)
#' assocES(exp_r, prot_r, eBayes=FALSE)
#' @export
assocES <- function(x, y, formula, select, set="exposures", sva=FALSE,
    vfilter=NULL, eBayes=TRUE, verbose=FALSE, warnings=TRUE, ...) {

    ## ----------------------------------------------------------------- ##
    # CHECK FOR INPUT CLASSES
    if(class(y) %in% c("ExposomeSet", "ExposomeClust") &
        class(x) == "ExpressionSet") {
        t <- x
        x <- y
        y <- x
        rm(t)
    }
    if(!(class(x) %in% c("ExposomeSet", "ExposomeClust") & class(y) == "ExpressionSet")) {
        stop("Required object of class 'ExposomeSet'/'ExposomeClust' and ",
             "'ExpressionSet' for arguments 'x, and 'y'.")
    }
    ## ----------------------------------------------------------------- ##

    ## ----------------------------------------------------------------- ##
    # TAKE THE COMMON SAMPLES
    if(warnings | verbose) {
        warning("Both input sets will be reduced to common samples.")
    }

    samples_selection <- intersect(
        Biobase::sampleNames(x),
        Biobase::sampleNames(y)
    )

    l1 <- sapply(list(x, y), function(n) length(Biobase::sampleNames(n)))
    x <- x[ , samples_selection]
    y <- y[ , samples_selection]
    l2 <- sapply(list(x, y), function(n) length(Biobase::sampleNames(n)))
    l3 <- mapply('-', l1, l2, SIMPLIFY = FALSE)
    names(l3) <- c(class(x), class(y))

    if(verbose) {
        message(paste(unlist(l3), names(l3),
                      sep = " samples were reduced from ", collapse = ", "))
    }

    if(length(samples_selection) == 0) {
        stop("No samles in common between both datasets.")
    }
    ## ----------------------------------------------------------------- ##

    ## ----------------------------------------------------------------- ##
    ## CONVERT FORMULA
    formula <- as.character(as.formula(formula))
    exp.dt <- data.frame(pData(x), expos(x))

    if(class(x) == "ExposomeClust") {
        ## ------------------------------------------------------------- ##
        ## EXPOSOME CLUSTER ANALYSIS
        if(warnings | verbose) {
            warning("Given an object of class 'ExposomeClust'. Association ",
                "test will be performed on clustering result.")
        }
        if(!missing(select) & (warnings | verbose)) {
            warning("Given values in 'select' argument. They will be dropped.")
        }
        select <- "cluster"
    } else {
        ## ------------------------------------------------------------- ##
        ## EXPOSOME SET ANALYSIS
        if(missing(select)) {
            if(set == "exposures") {
                warning("No given 'select'. Association tests will be ",
                    "computed for all exposures")
                select <- rexposome::exposureNames(x)
            } else { ## set == "phenotypes"
                warning("No given 'select'. Association test will be ",
                    "computed for all phenotypes")
                select <- Biobase::phenotypeNames(x)
            }
        }
    }
    ## ----------------------------------------------------------------- ##

    ## ----------------------------------------------------------------- ##
    ## EXTRACT SETS AND PERFORM ANALYSIS
    results <- lapply(select, function(ex) {
        design <- as.formula(paste0("~", ex, "+", formula[2]))
        if(verbose) {
            message("Evaluating model '", as.character(design), "'.")
        }
        pheno <- .create_p(
            expo.dt = exp.dt,
            omic.p = Biobase::pData(y),
            select = all.vars(design)
        )

        if(class(pheno) == "character") {
            stop("Invalid value '", pheno, "' in 'exposures' or 'covariates'")
        }

        na.loc <- rowSums(apply(pheno, 2, is.na))
        na.loc <- which(na.loc != 0)
        if(length(na.loc) != 0) {
            if(warnings | verbose) {
                warning("There are missing values. ", length(na.loc),
                        " samples will be removed.")
            }
            pheno <- pheno[-na.loc, , drop=FALSE]
        }


        if(sum(!sapply(sapply(apply(pheno, 2, table), length), ">", 1)) != 0) {
            warning("When testing for '", ex, "', at last one covariate ",
                    "is constant")
            list(
                N=NA,
                sva.num=NA,
                design=NA,
                result=NA,
                error=paste0("Covariate in model '", as.character(design) ,"' is constant")
            )
        } else {
            tryCatch({
                # Design model
                design.mm <- model.matrix(formula(design), data = pheno)
                y <- y[ , rownames(pheno), drop=FALSE]
                # If required, apply SVA
                n.sv <- NA
                if(sva) {
                    if (verbose | warnings){
                        if(is.null(vfilter)) {
                            message("Computing SVA. This step can be very time consuming.",
                                    "Try using argument 'vfilter'.")
                        } else {
                            message("Computing SVA. This step can be very time consuming.")
                        }
                    }

                    n.sv <- sva::num.sv(Biobase::exprs(y),
                                        design.mm, vfilter=vfilter)
                    if (n.sv > 0){
                        svobj <- sva::sva(Biobase::exprs(y), design.mm,
                                          design.mm[ , -1, drop=FALSE],
                                          n.sv=n.sv, vfilter=vfilter)
                        design.mm <- cbind(design.mm, svobj$sv)
                    }
                    rm(svobj)
                    suppressMessages(gc())
                }

                # Fit the model
                if (verbose){
                    message("Fitting the model.")
                }

                fit <- limma::lmFit(y, design.mm, ...)
                if(eBayes) {
                    fit <- limma::eBayes(fit)
                }

                list(
                    N=nrow(pheno),
                    sva.num=n.sv,
                    error=NA,
                    design=design,
                    result=fit
                )
            }, error=function(e) {
                if(warnings) {
                    warning(e)
                }
                return(list(
                    N=NA,
                    sva.num=NA,
                    error=e,
                    design=NA,
                    result=NULL
                ))
            })
        }
    })
    names(results) <- select

    new("ResultSet",
        fun_origin = "assocES",
        results = results,
        fData = list(fData(x), fData(y)),
        options = list(
            method="assocES",
            sva=sva,
            eBayes=eBayes
        )
    )
    ## ------------------------------------------------------------- ##
}
