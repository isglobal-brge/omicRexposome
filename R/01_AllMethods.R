#' Function to draw de result of an association study
#'
#' This function draws two type of plots for the ResultSet from association
#' functions
#'
#' @name plotAssociation
#' @rdname plotAssociation-methods
#' @aliases plotAssociation,ResultSet-methods
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
#' data(prot_r)
#' data(gexp_r)
#' data(exp_r)
#'
#' # Manhattan like plot
#' rst <- assocES(exp_r, prot_r, formula=~sex+age, eBayes=FALSE)
#' plotAssociation(rst, type="protein")
#'
#' # Volcano plot
#' rst <- assocES(exp_r, gexp_r, formula=~sex+age)
#' plotAssociation(rst, rid="Cotinine", type="qq")
#' @return An association plot
#' @export plotAssociation
#' @seealso \link{plotIntegration} for plotting integration results
setGeneric("plotAssociation", function(object,  rid = 1, coef = 2, contrast = 1,
                                       type = c("manhattan", "qq"), tPV, tFC,
                                       show.effect=FALSE)
    standardGeneric("plotAssociation")
)

#' Function to draw de result of an integration study
#'
#' This function draws a plots for the ResultSet from integration function
#'
#' @name plotIntegration
#' @rdname plotIntegration-methods
#' @aliases plotIntegration,ResultSet-methods
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
#' data(methy_r)
#' data(gexp_r)
#' rst <- crossomics(list(methy=methy_r, gexp=gexp_r), permute=NULL)
#' plotIntegration(rst)
#' @return An integration plot
#' @export plotIntegration
#' @seealso \link{plotAssociation} for plotting association results
setGeneric("plotIntegration", function(object, cmpX=1, cmpY=2, lb.th=0.20,
                                       legend.show=TRUE, colors, ...)
    standardGeneric("plotIntegration")
)

# -----------------------------------------------------------------------------

#' Obtain the "rid"s from a ResultSet
#'
#' This method resunts as character with the \code{rid}
#' in a given \link{ResultSet}.
#'
#' @name rid
#' @rdname rid-methods
#' @aliases rid
#' @param object An object of class \link{ResultSet}
#' @return A character vector of \code{rid}s.
#' @examples
#' data(prot_r)
#' data(exp_r)
#' rst <- assocES(exp_r, prot_r, formula=~sex+age, eBayes=FALSE)
#' rid(rst)
#' @export rid
setGeneric("rid", function(object)
    standardGeneric("rid")
)

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
#' data(gexp_r)
#' data(exp_r)
#' rst <- assocES(exp_r, gexp_r, formula=~sex+age)
#' tableLambda(rst)
#' @export tableLambda
#' @seealso \code{\link{tableHist}} for the number of hits per analysys,
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
#' @name tableHist
#' @rdname tableHist-methods
#' @aliases tableHits
#' @param object An object of class \link{ResultSet}
#' @param th (default \code{0.05}) Threshold (p-value) to considere a result
#' as a hit.
#' @return A labeled numeric vector with the exposures and the number of hits.
#' @examples
#' data(gexp_r)
#' data(exp_r)
#' rst <- assocES(exp_r, gexp_r, formula=~sex+age)
#' tableHits(rst, th=0.05)
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
#' data(gexp_r)
#' data(exp_r)
#' rst <- assocES(exp_r, gexp_r, formula=~sex+age)
#' plotLambda(rst)
#' @return An lambda plot
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
#' @examples
#' data(gexp_r)
#' data(exp_r)
#' rst <- assocES(exp_r, gexp_r, formula=~sex+age)
#' plotHits(rst)
#' @return An hits plot
#' @export plotHits
setGeneric("plotHits", function(object, th=0.05, width=0.75)
    standardGeneric("plotHits")
)

#' Method to extrat feature result from a ResultSet
#'
#' Homologous methods from \code{limma}, \code{topTable} resturns a
#' \code{data.frame} with the \code{logFC} and  \code{PValue} per
#' featrue for the selcted \code{coef} and for given result (\code{rid}).
#'
#' @name topTable
#' @rdname topTable-methods
#' @aliases topTable
#' @param object A \code{\link{ResultSet}} object.
#' @param rid The name or index of the result to be extracted.
#' @param coef (default \code{2}) Index of the coefficient to be extracted.
#' @param contrast (default \code{1}) When \code{code} corresponds to a
#' multicategorical variable, contasr selects the comparison.
#' @param sort (default \code{TRUE}) If \code{TRUE}, results are ordered
#' by P-Value.
#' @return A \code{data.frame} with the result of the association study,
#' including P-Value and Fold Change.
#' @examples
#' data(gexp_r)
#' data(exp_r)
#' rst <- assocES(exp_r, gexp_r, formula=~sex+age)
#' topTable(rst, rid=1)
#' @export topTable
setGeneric("topTable", function(object, rid, coef=2, contrast=1, sort = TRUE)
    standardGeneric("topTable")
)

#' Method to get the options sued to create the ResultSet
#'
#' Method that returns a list with the options used to create the
#' \code{ResultSet}.
#'
#' @name opt
#' @rdname opt-methods
#' @aliases opt
#' @param object A \code{\link{ResultSet}} object.
#' @return A list with the options used to create the \code{ResultSet}.
#' @examples
#' data(gexp_r)
#' data(exp_r)
#' rst <- assocES(exp_r, gexp_r, formula=~sex+age)
#' opt(rst)
#' @export opt
setGeneric("opt", function(object)
    standardGeneric("opt")
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
#' data("exposome")
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
#' @aliases add_cls
#' @param object An object of class \link{MultiDataSet}.
#' @param clsSet An object of class \link{ExposomeClust}.
#' @param ... Arguments given to \link{add_eset} from \link{MultiDataSet}.
#' @return A \link{MultiDataSet} with the \link{ExpressionSet} added as an
#' independent dataset.
#' @examples
#' data("eclust")
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
#' expression from microarray and the exposome. An \link{ExpresionSet} is
#' the object storing the gene expresion and an \link{ExposomeSet} the one
#' storing the exposome. Both of them needs to be encapsulated in a
#' \link{MultiDataSet}. The association study is perform through standard
#' \link{limma} pipeline. The function allows to perform multiple tests using
#' the argument \code{exposures}.
#'
#' @name association
#' @rdname association-methods
#' @export association
setGeneric("association", function(object, formula, expset, omicset, set = "exposures",
        method = "ls", ..., sva = FALSE, ebayes = TRUE, verbose = FALSE, warnings = TRUE)
    standardGeneric("association")
)


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
#' data(methy_r)
#' data(gexp_r)
#' rst <- crossomics(list(methy=methy_r, gexp=gexp_r), permute=NULL)
#' rst
#' @name crossomics
#' @rdname crossomics-methods
#' @export crossomics
setGeneric("crossomics", function(object, method="mcca", ncomponents=2, ...,
        na.rm=FALSE, permute = c(100, 3), verbose=FALSE, warnings=TRUE)
    standardGeneric("crossomics")
)


