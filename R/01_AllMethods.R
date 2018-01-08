#' Function to draw de result of an association study
#'
#' This function draws two type of plots for the ResultSet from association
#' functions
#'
#' @name plotAssociation
#' @rdname plotAssociation-methods
#' @aliases plotAssociation
#' @param object An object of class \link{ResultSet} obtained from assoc_*
#' functions.
#' @param rid (default \code{1}) Index or name of the test to be plotted.
#  Not used it the \code{ResultSet} comes from \link{assocSNP}.
#' @param coef (default \code{2}) Index of the coefficient to be extracted.
#' @param contrast (default \code{1}) When \code{code} corresponds to a
#' multicategorical variable, contasr selects the comparison.
#' @param type Can take \code{"volcano"}, \code{"qq"},  \code{"manhattan"} and
#' \code{"protein"}. \code{"protein"} lot is a type of Manhattan plot designed
#' for protein association analysis.
#' @param tPV (optional) Threshold for P.Value when \code{type="volcano"}.
#' @param tFC (optional) Threshold for Fold Change or Effect when
#' \code{type="volcano"}.
#' @param show.effect (default \code{FALSE}) If set to \code{TRUE}, when
#' \code{type="volcano"} the X-axis will show \code{2^logFC} instead of
#' \code{logFC}.
#' @examples
#' data("asr", package = "omicRexposome")
#' plotAssociation(asr, type = "qq")
#' plotAssociation(asr, type = "volcano")
#' @return A ggplot2 object
#' @export plotAssociation
#' @seealso \code{\link{plotIntegration}} for plotting integration results.
#' \code{\link{association}} to create a \code{ResultSet} to be passed to
#' this function.
setGeneric("plotAssociation", function(object,  rid = 1, coef = 2, contrast = 1,
        type = c("manhattan", "qq", "volcano"), tPV = NULL, tFC = NULL,
        show.effect=FALSE)
    standardGeneric("plotAssociation")
)

#' Function to draw de result of an integration study
#'
#' This function draws a plots for the ResultSet from integration function
#'
#' @name plotIntegration
#' @rdname plotIntegration-methods
#' @aliases plotIntegration
#' @param object An object of class \link{ResultSet} obtained from
#' \link{crossomics}.
#' @param cmpX (default \code{1}) Value of the X-axis when ploting rsults
#' from \link{mcia}.
#' @param cmpY (default \code{2}) Value of the Y-axis when ploting rsults
#' from \link{mcia}.
#' @param lb.th (default \code{0.20}) Threshold to place labels on radar chart
#' drawn when ploting results from \link{MultiCCA}.
#' @param legend.show (default \code{TRUE}) If set to FALSE, right legend
#' of radar plot is hidden when ploting results from \link{MultiCCA}.
#' @param colors (optional) Names vector with the colors sued to draw
#' each dataset. Used when ploting results from \link{MultiCCA}. If missing,
#' random colores are chosen.
#' @param ... Optional arguments are given to \code{plot} from \link{omicade4}
#' pacage (argument \code{axes} is filled with values from \code{cmpX} and
#' \code{cmpY}).
#' @examples
#' data("crs", package = "omicRexposome")
#' plotIntegration(crs)
#' @return A ggplot2 object
#' @export plotIntegration
#' @seealso \code{\link{plotAssociation}} for plotting association results.
#' \code{\link{crossomics}} to create a \code{ResultSet} to be passed to
#' this function.
setGeneric("plotIntegration", function(object, cmpX=1, cmpY=2, lb.th=0.20,
                                       legend.show=TRUE, colors, ...)
    standardGeneric("plotIntegration")
)

# -----------------------------------------------------------------------------

#' Method to extrat integration-feature result from a ResultSet
#'
#' Homologous methods from \code{MultiDataSet} (\code{getAssociation}) but
#' for \code{ResultsSet} created by \code{\link{crossomics}}. It Resturns a
#' \code{data.frame} with the result from \code{mcia} (\code{omicade4}) or
#' from \code{MultiCCA} (\code{PMA}).
#'
#' @name getIntegration
#' @rdname getIntegration-methods
#' @aliases getIntegration
#' @param object An object of class \link{ResultSet} obtained from
#' @param ... NOT USED
#' @return A \code{data.frame}
#' @examples
#' data("crs", package = "omicRexposome")
#' class(getIntegration(crs))
#' @export getIntegration
setGeneric("getIntegration", function(object, ...)
    standardGeneric("getIntegration")
)


# -----------------------------------------------------------------------------

#' Compute a lambda score on the results stored in a ResultSet
#'
#' Compute lambda score on each result in the given \link{ResultSet} by using
#' \code{lambdaClayton}.
#'
#' @name tableLambda
#' @rdname tableLambda-methods
#' @aliases tableLambda
#' @param object An object of class \link{ResultSet}
#' @return Returns a \code{data.frame} having the exposures and the computed
#' lambda score.
#' @param trim (default \code{0.5}) percentage of right omited values for
#' \link{lambdaClayton}.
#' @return A labeled numeric vector with the lambda score for each exposure.
#' @examples
#' data("asr", package = "omicRexposome")
#' tableLambda(asr)
#' @export tableLambda
#' @seealso \code{\link{tableHits}} for the number of hits per analysys,
#' \code{\link{plotHits}} for a graphical representation of the hists
#' per analysys, \code{\link{plotLambda}} for a graphical representation of
#' the lambda score per analysys
setGeneric("tableLambda", function(object, trim=0.5)
    standardGeneric("tableLambda")
)

#' Counts the number of hits on the results stored in a ResultSet
#'
#' Given a threshold it counts the number of hits in each result in the
#' given \link{ResultSet}.
#'
#' @name tableHits
#' @rdname tableHits-methods
#' @aliases tableHits
#' @param object An object of class \link{ResultSet}
#' @param th (default \code{0.05}) Threshold (p-value) to considere a result
#' as a hit.
#' @return A labeled numeric vector with the exposures and the number of hits.
#' @examples
#' data("asr", package = "omicRexposome")
#' tableHits(asr)
#' @export tableHits
#' @seealso \code{\link{tableLambda}} for the lambda score per analysys,
#' \code{\link{plotLambda}} for a graphical representation of
#' the lambda score per analysys, \code{\link{plotHits}} for a graphical
#' representation of the hists per analysys
setGeneric("tableHits", function(object, th=0.05)
    standardGeneric("tableHits")
)

#' Plot lambda score for all results in a ResultSet
#'
#' This method draws a baplor with the lambda score of each result in the
#' given \link{ResultSet}.
#'
#' @name plotLambda
#' @rdname plotLambda-methods
#' @aliases plotLambda
#' @param object An object of class \link{ResultSet}
#' @param width (default \code{0.70}) width of the bar
#' @examples
#' data("asr", package = "omicRexposome")
#' plotLambda(asr)
#' @return A ggplot2 object
#' @seealso \code{\link{plotHits}} for a graphical representation of
#' the hits per analysys, \code{\link{tableLambda}} for the lambda
#' score per analysys, \code{\link{tableHits}} for the hists per analysys
#' @export plotLambda
setGeneric("plotLambda", function(object, width=0.75)
    standardGeneric("plotLambda")
)

#' Plot number of hits per result in ResultSet
#'
#' This method draws a barplot with the number of hits in each result
#' stored in the given \link{ResultSet}.
#'
#' @name plotHits
#' @rdname plotHits-methods
#' @aliases plotHits
#' @param object An object of class \link{ResultSet}
#' @param th (default \code{0.05}) Threshold (p-value) to considere a result
#' as a hit.
#' @param width (default \code{0.70}) width of the bar
#' @return A ggplot2 object
#' @examples
#' data(asr, package = "omicRexposome")
#' plotHits(asr)
#' @seealso \code{\link{plotLambda}} for a graphical representation of
#' the lambda score per analysys, \code{\link{tableLambda}} for the lambda
#' score per analysys, \code{\link{tableHits}} for the hists per analysys
#' @export plotHits
setGeneric("plotHits", function(object, th=0.05, width=0.75)
    standardGeneric("plotHits")
)


# -----------------------------------------------------------------------------

#' Method to add an ExposomeSet to a MultiDataSet
#'
#' This method allows to insert an object of class \link{ExposomeSet} as an
#' independent dataset into an object of class \link{MultiDataSet}.
#'
#' @name add_exp
#' @rdname add_exp-methods
#' @aliases add_exp
#' @param object An object of class \link{MultiDataSet}.
#' @param expoSet An object of class \link{ExposomeSet}.
#' @param warnings (default \code{TRUE}) If set to \code{FALSE} warnings will
#' not be displayed.
#' @param ... Arguments given to \link{add_eset} from \link{MultiDataSet}.
#' @return A \link{MultiDataSet} with the \link{ExpressionSet} added as an
#' independent dataset.
#' @examples
#' data("exposome", package = "rexposome")
#' library(MultiDataSet)
#' md <- new("MultiDataSet")
#' names(md)
#' md <- add_exp(md, expo)
#' names(md)
#' @export add_exp
setGeneric("add_exp", function(object, expoSet, warnings = TRUE, ...)
    standardGeneric("add_exp")
)

#' Method to add an ExposomeClust to a MultiDataSet
#'
#' This method allows to insert an object of class \link{ExposomeClust} as an
#' independent dataset into an object of class \link{MultiDataSet}.
#'
#' @name add_cls
#' @rdname add_cls-methods
#'
#' @param object An object of class \link{MultiDataSet}.
#' @param clsSet An object of class \link{ExposomeClust}.
#' @param ... Arguments given to \link{add_eset} from \link{MultiDataSet}.
#' @return A \link{MultiDataSet} with the \link{ExpressionSet} added as an
#' independent dataset.
#' @examples
#' data("eclust", package = "rexposome")
#' library(MultiDataSet)
#' md <- new("MultiDataSet")
#' names(md)
#' md <- add_cls(md, expo_c)
#' names(md)
#' @export add_cls
setGeneric("add_cls", function(object, clsSet, ...)
    standardGeneric("add_cls")
)

#' # -----------------------------------------------------------------------------

#' Method to perform an association study between transcriptome and exposom
#'
#' This function allows to perform an association study between gene
#' expression from microarray and the exposome. An \code{ExpresionSet} is
#' the object storing the gene expresion and an \code{ExposomeSet} the one
#' storing the exposome. Both of them needs to be encapsulated in a
#' \code{MultiDataSet}. The association study is perform through standard
#' \code{limma} pipeline. The function allows to perform multiple tests using
#' the argument \code{exposures}.
#'
#' @name association
#' @rdname association-methods
#' @aliases association
#' @param object A \code{MultiDataSet} object containing at last one omic
#' data-sets like \code{ExpressionSet}, \code{MethylationSet}... and, at last,
#' one \code{ExposomeSet}.
#' @param formula formula to be evaluated by each exposure (or phenotype, see
#' \code{set} argument). It should not contain any exposures (or phenotype),
#' it will be added automatically when evaluated.
#' @param expset Name of the \code{ExposomeSet} in \code{object}.
#' @param omicset Name of the omic data-set in \code{object}.
#' @param set (default \code{"exposures"}) Can take value \code{"exposures"}
#' to test the association of the exposures in the \code{ExposomeSet} vs.
#' the features in the omic data-set. If takes \code{"phenotypes"} all
#' phenotypes in \code{ExposomeSet} are tested.
#' @param method (default \code{"lm"}) Check \code{limma} help pages.
#' @param ... Arguments passed to \code{limma}'s \code{lmFit}.
#' @param baselevels (optional) If set, must be a labeled vector with the
#' default base level for categorical exposures.
#' @param sva (default \code{"none"}). This argument can take value
#' \code{"none"} to do not apply SVA. Value \code{"fast"} will run SVA
#' using \code{isva} and \code{SmartSVA}. Value \code{"slow"}
#' will run SVA using \code{sva}.
#' @param vfilter (default \code{NULL}). Only used when \code{sva = "slow"}.
#' Numeric number of probes used in \link{sva}. Recomended ~10\% of real probes.
#' @param verbose (default \code{FALSE}) If set to \code{TRUE}, a series of
#' messages descriving the process are shown.
#' @param warnings (default \code{TRUE}) If set to \code{TRUE}, a series of
#' warnings are shown when required user atention.
#' @return An object of class \code{\link{ResultSet}}.
#' @examples
#' library(MultiDataSet)
#' data(brge_prot, package = "brgedata")
#' data(brge_expo, package = "brgedata")
#' mds <- createMultiDataSet()
#' mds <- add_eset(mds, brge_prot, dataset.type = "proteines")
#' mds <- add_eset(mds, brge_expo, dataset.type = "exposures", GRanges = NA)
#'
#' asr <- association(mds, formula = Asthma ~ Sex + Age,
#'   expset = "exposures", omicset = "proteines")
#' asr
#' @export association
setGeneric("association", function(object, formula, expset, omicset,
        set = "exposures", method = "ls", ..., baselevels, sva = "none",
        vfilter = NULL, verbose = FALSE, warnings = TRUE)
    standardGeneric("association")
)


#' Function to perform a Transcriptome-Wide Association Study
#'
#' This function allows to perform a Transcriptome-Wide Association Study
#' by using an \code{ExposmeSet} and an \code{ExpressionSet}. It
#' allows to perform an adjustment using Surrogate Variable Analysis (from
#' R package \code{sva}).
#'
#' @name crossomics
#' @rdname crossomics-methods
#' @aliases crossomics
#' @param object A \code{MultiDataSet} object containing at last two data-sets
#' like \code{ExposomeSet}, \code{ExpressionSet}, \code{MethylationSet}...
#' @param method (default \code{"mcca"}) It can takes values \code{"mcca"} for
#' Multiple Canonical Correlation Analysis or \code{"mcia"} for Multiple
#' Co-Inertia Analysis.
#' @param ncomponents (default \code{2}) Number of components to be estimated.
#' @param ... Other arguments given to \code{mcia} (from \code{omicade4}) or
#' to \code{MultiCCA} (from \code{PMA}).
#' @param na.rm (default \code{FALSE}) If \code{method} was set to
#' \code{"mcca"} and \code{na.rm} was set to \code{TRUE}, features containing
#' missing values are removed.
#' @param permute (default \code{c(100, 3)}). If \code{method="mcca"} and this
#' agument is set to \code{NULL} no permutation test to tune-up the parameters
#' for \code{MultiCCA}. When filles, \code{permute[1]} corresponds to
#' the number permutations (default in \code{MultiCCa.permute} is \code{25})
#' and \code{permute[2]} the number of iterations
#' (default in \code{MultiCCA.permute} is 3).
#' @param verbose (default \code{FALSE}) If set to \code{TRUE}, a series of
#' messages descriving the process are shown.
#' @param warnings (default \code{TRUE}) If set to \code{TRUE}, a series of
#' warnings are shown when required user atention.
#' @return An object of class \code{\link{ResultSet}}.
#' @examples
#' library(MultiDataSet)
#' library(rexposome)
#' data(brge_prot, package = "brgedata")
#' data(brge_expo, package = "brgedata")
#' mds <- createMultiDataSet()
#' mds <- add_eset(mds, brge_prot, dataset.type = "proteines")
#' mds <- add_eset(mds, imputation(brge_expo),
#'     dataset.type = "exposures", GRanges = NA)
#'
#' crs <- crossomics(mds, method = "mcia")
#' crs
#' @export crossomics
setGeneric("crossomics", function(object, method="mcca", ncomponents=2, ...,
        na.rm=FALSE, permute = c(100, 3), verbose=FALSE, warnings=TRUE)
    standardGeneric("crossomics")
)


