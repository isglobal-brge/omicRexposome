#' Method to add an ExpressionSet to a MultiDataSet tagged as
#' protein data container.
#'
#' This method allows to insert an object of class \link{ExpressionSet} as an
#' independent dataset into an object of class \link{MultiDataSet}. The
#' \link{ExpressionSet} will be tagged, in the \link{MultiDataSet} as
#' \code{"protein"}.
#'
#' @param object An object of class \link{MultiDataSet}.
#' @param expoSet An object of class \link{ExpressionSet} with protein data.
#' @param warnings (default \code{TRUE}) If set to \code{FALSE} warnings will
#' not be displayed.
#' @param ... Arguments given to \link{add_eset} from \link{MultiDataSet}.
#' @return A \link{MultiDataSet} with the \link{ExpressionSet} added as an
#' independent dataset.
#' @export add_prot
setGeneric("add_prot", function(object, expoSet, ...)
    standardGeneric("add_prot")
)

# -----------------------------------------------------------------------------

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
#' @name assocGE
#' @rdname assocGE-methods
#' @aliases assocGE,MultiDataSet-rexposome
#' @param object An object of class \link{MultiDataSet} containing an
#' \code{ExposomeSet} and an \link{ExposomeSet}.
#' @param formula Model (without exposure/s) to be tested.
#' @param select One or multiple exposures existing into the
#' \link{ExposomeSet} (encapsulated into the \link{MultiDataSet}).
#' @param set (default \code{"exposures"}) Can take values \code{"exposures"}
#' of \code{"phenotype"}.
#' @param ... Arguments passed to \code{\link{lmFit}}.
#' @param sva (default \code{FALSE}) If \code{TRUE} SVAnalysis is done.
#' @param vfilter (default \code{NULL}) Numeric number of probes used
#' in \link{sva} if argument \code{sva} is set to \code{TRUE}.de
#' @param ncores (default \code{1}) Number of test performed at a time.
#' @param verbose (default \code{FALSE}) If set to \code{TRUE} information
#' about the process is show.
#' @param warnings (default \code{TRUE}) If set to \code{FALSE} warnings will
#' not be displayed.
#' @return An object of class \code{ResultSet} containing the results of the
#' test and the associated feature-data.
#' @export assocGE
#' @seealso \link{assocME} to test the association between epigenome and
#' exposome, \link{assocSNP} to test the association between genome and
#' exposome, \link{assocPRT} to test the association between proteome and
#' exposome, package \link{limma} for more details, \link{plotAssociation} to
#' plot the results
setGeneric("assocGE", function(object, formula, select, set = "exposures",
                               ..., sva = FALSE, vfilter = NULL, ncores = 1,
                               verbose = FALSE, warnings = TRUE)
    standardGeneric("assocGE")
)

#' Method to perform an association study between proteome and exposom
#'
#' This function allows to perform an association study between protein
#' expression from microarray and the exposome. An \link{ExpresionSet} is
#' the object storing the protein expresion and an \link{ExposomeSet} the one
#' storing the exposome. Both of them needs to be encapsulated in a
#' \link{MultiDataSet}. The association study is perform through standard
#' \link{limma} pipeline. The function allows to perform multiple tests using
#' the argument \code{exposures}.
#'
#' @name assocPRT
#' @rdname assocPRT-methods
#' @aliases assocPRT,MultiDataSet-rexposome
#' @param object An object of class \link{MultiDataSet} containing an
#' \code{ExposomeSet} and an \link{ExposomeSet}.
#' @param formula Model (without exposure/s) to be tested.
#' @param select One or multiple exposures existing into the
#' \link{ExposomeSet} (encapsulated into the \link{MultiDataSet}).
#' @param set (default \code{"exposures"}) if set to \code{"phenotypes"} the
#' variables used in the association test are the phenotypes instead of the
#' exposures.
#' @param ... ...
#' @param ebayes (default \code{TRUE}) If \code{TRUE}, fitted limma modes is
#' given to \link{eBayes} before extrating results.
#' @param sva (default \code{FALSE}) if set to \code{TRUE} SVa is applied on
#' the model and the surrogate variables are incldue in the association test.
#' @param vfilter (default \code{NULL}) Numeric number of probes used
#' in \link{sva} if argument \code{sva} is set to \code{TRUE}.
#' @param ncores (default \code{1}) Number of test performed at a time.
#' @param verbose (default \code{FALSE}) If set to \code{TRUE} information
#' about the process is show.
#' @param warnings (default \code{TRUE}) If set to \code{FALSE} warnings will
#' not be displayed.
#' @return An object of class \code{ResultSet} containing the results of the
#' test and the associated feature-data.
#' @export assocPRT
#' @seealso \link{assocME} to test the association between epigenome and
#' exposome, \link{assocSNP} to test the association between genome and
#' exposone, \link{assocGE} to test the association between gene expression
#' and exposome, package \link{limma} for more details,
#' \link{plotAssociation} to plot the results
setGeneric("assocPRT", function(object, formula, exposures, ..., sva = TRUE,
                                vfilter = NULL, ncores = 1,
                                verbose = FALSE, warnings = TRUE)
    standardGeneric("assocPRT")
)

#' Method to perform an association study between epigenome and exposome
#'
#' This function allows to perform epigenetic wide association studies
#' between methylation and expsome. A \link{MethylationSet} is the object
#' storing the CpGs and an \link{ExposomeSet} the one storing the exposures.
#' The association study relies on package \link{MEAL}. The function allows to
#' perform multiple tests using the argument \code{exposures}.
#'
#' @name assocME
#' @rdname assocME-methods
#' @aliases assocME,MultiDataSet-rexposome
#' @param object An object of class \link{MultiDataSet} containing an
#' \code{ExposomeSet} and a \link{MethylationSet}.
#' #' @param formula Model (without exposure/s) to be tested.
#' @param select One or multiple exposures existing into the
#' \link{ExposomeSet} (encapsulated into the \link{MultiDataSet}).
#' @param set (default \code{"exposures"}) Can take values \code{"exposures"}
#' of \code{"phenotype"}.
#' @param area.test (default \code{FALSE}) If set to \code{TRUE} it tests
#' area instead of CpG.
#' @param method (default \code{"ls"}) Can takes \code{"ls"} or \code{"robust"}
#' and it will be the method used by \link{MEAL}/\code{limma} to test
#' the association.
#' @param betas (default \code{FALSE}) By default M-values are used. If set
#' to \code{TRUE} beta-vaues are used instead.
#' @param ... Other arguments passed to \link{DAProbe} or \link{DARegion}.
#' @param sva (default \code{FALSE}) If set to \code{TRUE}, surrogate
#' variable analysis is done before calling \link{DAProbe}, or \link{DARegion},
#' and obtained surrogate variables are added.
#' @param vfilter (default \code{NULL}) Numeric number of probes used
#' in \link{sva} if argument \code{sva} is set to \code{TRUE}.
#' @param ncores (default \code{1}) nomber of analysis that will be done
#' at the same time (see \link{mclapply}).
#' @param verbose (default \code{FALSE}) If set to \code{TRUE} information
#' about the process is show.
#' @param warnings (default \code{TRUE}) If set to \code{FALSE} warnings will
#' not be displayed.
#' @return An object of class \code{ResultSet} containing the results of the
#' test and the associated feature-data.
#' @export assocME
#' @importClassesFrom MultiDataSet MethylationSet
#' @seealso \link{assocGE} to test the association between transcriptome
#' and exposome, \link{assocSNP} to test the association between genome and
#' exposome, package \link{MEAL} for more details, \link{plotAssociation} to
#' plot the results
setGeneric("assocME", function(object, formula, select, set="exposures",
                               area.test=FALSE, method="ls", betas=FALSE,
                               ..., sva=FALSE, vfilter = NULL,
                               ncores=1, verbose=FALSE,
                               warnings=TRUE)
    standardGeneric("assocME")
)

#' Method to perform an association study between genome and diseasome
#' (exposome's phenotypes)
#'
#' This function allows to perform genome wide association studies between
#' SNPs and phenotypes. A \link{SnpSet} is the object storing the genotypes and
#' an \link{ExposomeSet} the phenotypes. The association study is performed
#' through \link{snp.rhs.tests} from \link{snpStats}.
#'
#' @section Warning:
#' A difference between \code{assocSNP} and both \code{assocGE}
#' and \code{assocME}, is that \code{assocSNP} cannot perform multiple
#' test in a single run.
#'
#' @name assocSNP
#' @rdname assocSNP-methods
#' @aliases assocSNP,MultiDataSet-rexposome
#' @param object An object of class \link{MultiDataSet} containing an
#' \code{ExposomeSet} and a \link{SnpSet}.
#' @param design Formula to be tested at \link{snp.rhs.tests}.
#' @param family (default \code{"binomial"}) A string defining the generalized
#' linear model family. This can take \code{"binomial"}, \code{"Poisson"},
#' \code{"Gaussian"} or \code{"gamma"}. It is given to \link{snp.rhs.tests}.
#' @param filter.callrate (default \code{0.98}) Filter for SNPs call-rate. Ones
#' under the threshold are removed.
#' @param filter.maf (default \code{0.05}) Filter for minimum allele frequency
#' (MAF). Ones under the threshold are removed.
#' @param filter.pvalHEW (default \code{0.001}) Filter for Hardyâ€“Weinberg. Ones
#' under the threshold are removed.
#' @param ... Argument given to \link{snp.rhs.tests} from \link{snpStats}.
#' @export assocSNP
#' @seealso \link{assocGE} to test the association between transcriptome
#' and exposome, \link{assocME} to test the association between epigenome
#' and exposome, package \link{snpStats} for more details,
#' \link{plotAssociation} to plot the results
setGeneric("assocSNP", function(object, design, family = "binomial",
                                 filter.callrate = 0.98, filter.maf = 0.05, filter.pvalHWE = 0.001, ...)
    standardGeneric("assocSNP")
)

#' #' Method to perform a crossomics integration analysis.
#' #'
#' #' This functions performs a crossomics analysis using all the datasets stored
#' #' into a \link{MultiDataSet} object. The process of integration performs a
#' #' sparse multiple canonical correlation analysis with all the available
#' #' tables. The genome (SNPs) data is transformed to a continuous variable.
#' #'
#' #' @section Warning:
#' #' Current accepted datasets are: \link{ExposomeSer} for exposures data.
#' #' \link{ExpressionSet} for micro/gene expression data, \link{MethylationSet}
#' #' for methylome (CpGs) data, \link{SnpSet} for genome (SNPs) data. Proteome
#' #' data is allowed as an \link{eSet} while tagged as \code{"proteome"} into the
#' #' \link{MultiDataSet} object.
#' #'
#' #' @name crossomics
#' #' @rdname crossomics-methods
#' #' @aliases crossomics,MultiDataSet-rexposome
#' #' @param object An object of class \link{MultiDataSet} containing the multiple
#' #' sets to be integrated
#' #' @param ncomponents (default \code{2})
#' #' @param na.rm (default \code{FALSE})
#' #' @param verbose (default \code{FALSE}) If set to \code{TRUE} information
#' #' about the process is show.
#' #' @param warnings (default \code{TRUE}) If set to \code{FALSE} warnings will
#' #' not be displayed.
#' #' @references
#' #' Witten, DM and Tibshirani, R and Hastie, T. A penalized matrix
#' #' decomposition, with applications to sparse principal components and canonical
#' #' correlation analysis. Biostatistics. 2009. 10.1093/biostatistics/kxp008
#' #'
#' #' Abraham, G and Inouye M. Fast Principal Component Analysis of Large-Scale
#' #' Genome-Wide Data. PLOS ONE. http://dx.doi.org/10.1371/journal.pone.0093766
#' #' @export crossomics
#' #' @seealso the used functions \link{MultiCCA} and \link{MultiCCA.permut} or
#' #' the package \link{PMA} for details, \link{plotIntegration} to plot the
#' #' results
#' setGeneric("crossomics", function(object, ncomponents = 2, ..., na.rm = FALSE,
#'                                   verbose = FALSE, warnings = TRUE)
#'     standardGeneric("crossomics")
#' )

# -----------------------------------------------------------------------------

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
#' @param rid (default \code{1}) Index or name of the test to be plotted. Not
#' used it the \code{ResultSet} comes from \link{assocSNP}.
#' @param type Can take \code{"manhattan"} or \code{"qq"}.
#' @param ... Other arguments used in the different plots.
#' @export plotAssociation
#' @seealso \link{plotIntegration} for plotting integration results
setGeneric("plotAssociation", function(object,  rid = 1,
                                       type = c("manhattan", "qq"), ...)
    standardGeneric("plotAssociation")
)

#' #' Function to draw de result of an integration study
#' #'
#' #' This function draws a plots for the ResultSet from integration function
#' #'
#' #' @name plotIntegration
#' #' @rdname plotIntegration-methods
#' #' @aliases plotIntegration,ResultSet-methods
#' #' @param object An object of class \link{ResultSet} obtained from
#' #' \link{crossomics}.
#' #' @param cmpX (default \code{1}) Value of the X-axis when ploting rsults
#' #' from \link{mcia}.
#' #' @param cmpY (default \code{2}) Value of the Y-axis when ploting rsults
#' #' from \link{mcia}.
#' #' @param tcolors (optional) Names vector with the colors sued to draw
#' #' each dataset. Used when ploting results from \link{MultiCCA}. If missing,
#' #' random colores are chosen.
#' #' @param lb.th (default \code{0.20}) Threshold to place labels on radar chart
#' #' drawn when ploting results from \link{MultiCCA}.
#' #' @param legend.show (default \code{TRUE}) If set to FALSE, right legend
#' #' of radar plot is hidden when ploting results from \link{MultiCCA}.
#' #' @param ... Optional arguments are given to \code{plot} from \link{omicade4}
#' #' pacage (argument \code{axes} is filled with values from \code{cmpX} and
#' #' \code{cmpY}).
#' #' @export plotIntegration
#' #' @seealso \link{plotAssociation} for plotting association results
#' setGeneric("plotIntegration", function(object, cmpX=1, cmpY=2, tcolors,
#'         lb.th=0.20, legend.show=TRUE, ...)
#'     standardGeneric("plotIntegration")
#' )

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
#' @export tableLambda
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
#' @export tableHits
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
#' @param trim (default \code{0.5}) percentage of right omited values for
#' \link{lambdaClayton}.
#' @param width (default \code{0.70}) width of the bar
setGeneric("plotLambda", function(object, trim=0.5, width=0.75)
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
#' @param sort (default \code{TRUE}) If \code{TRUE}, results are ordered
#' by P-Value.
#' @export
setGeneric("topTable", function(object, rid, coef=2, sort = TRUE)
    standardGeneric("topTable")
)
