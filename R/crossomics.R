#' Function to perform a Transcriptome-Wide Association Study
#'
#' This function allows to perform a Transcriptome-Wide Association Study
#' by using an \code{ExposmeSet} and an \code{ExpressionSet}. It
#' allows to perform an adjustment using Surrogate Variable Analysis (from
#' R package \code{sva}).
#'
#' @param list A list containing at last two \code{eSet} based objects
#' like \code{ExposomeSet}, \code{ExpressionSet} or \code{MethylationSet}.
#' @param method (default \code{"mcca"}) It can takes values \code{"mcca"} for
#' Multiple Canonical Correlation Analysis or \code{"mcia"} for Multiple
#' Co-Inertia Analysis.
#' @param ncomponents (default \code{2}) Number of components to be estimated.
#' @param na.rm (default \code{FALSE}) If \code{method} was set to
#' \code{"mcca"} and \code{na.rm} was set to \code{TRUE}, features containing
#' missing values are removed.
#' @param permute (default \code{c(100, 3)}). If \code{method="mcca"} and this
#' agument is set to \code{NULL} no permutation test to tune-up the parameters
#' for \code{MultiCCA}. When filles, \code{permute[1]} corresponds to
#' the number permutations (default in \code{MultiCCa.permute} is \code{25})
#' and \code{permute[2]} the number of iterations
#' (default in \code{MultiCCA.permute} is 3).
#' @param ... Other arguments given to \code{mcia} (from \code{omicade4}) or
#' to \code{MultiCCA} (from \code{PMA}).
#' @param verbose (default \code{FALSE}) If set to \code{TRUE}, a series of
#' messages descriving the process are shown.
#' @param warnings (default \code{TRUE}) If set to \code{TRUE}, a series of
#' warnings are shown when required user atention.
#' @return An object of class \code{\link{ResultSet}}.
#' @examples
#' dta(methy_r)
#' data(gexp_r)
#' rst <- crossomics(list(methy=methy_r, gexp=gexp_r), permute=NULL)
#' rst
#' @export
crossomics <- function(list, method="mcca", ncomponents=2, ..., na.rm=FALSE,
        permute = c(100, 3), verbose=FALSE, warnings=TRUE) {
    ## --------------------------------------------------------------------- ##
    ## GENERAL CHECKS
    method <- match.arg(method, choices = c("mcca", "mcia"))
    if(length(list) < 2) {
        stop("At last two different datasets are required for integration processes.")
    }
    if(is.null(names(list))) {
        names(list) <- paste("set", 1:length(list))
    }
    ## --------------------------------------------------------------------- ##

    ## --------------------------------------------------------------------- ##
    ## REDUCE DATASETS TO COMMON SAMPLES
    if(warnings | verbose) {
        warning("Sets in list will be reduced to common samples.")
    }

    sc <- Reduce(intersect, lapply(list, Biobase::sampleNames))
    list <- lapply(list, function(it) it[ , sc])
    ## --------------------------------------------------------------------- ##

    if(method == "mcca") {
        .crossomics_mcca_list(list, ncomponents=ncomponents, na.rm=na.rm,
            permute=permute, verbose=verbose, warnings=warnings, ...)
    } else if(method == "mcia") {
        .crossomics_mcia_list(list, verbose=verbose, warnings=warnings, ...)
    } else {
        stop("Invalid method (", method, ") was given.")
    }
}
