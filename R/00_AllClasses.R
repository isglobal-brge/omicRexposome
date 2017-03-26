#' Class ResultSet
#'
#' Class \code{ResultSet} encapsulates the results from the association and
#' integration processes developed in \link{\code{rexposome}}. This include the
#' functions \code{\link{assocSNP}}, \code{\link{assocGE}},
#' \code{\link{assocME}}, \code{\link{assocPRT}} and \code{\link{crossomics}}.
#'
#' @name ResultSet
#' @aliases ResultSet-class
#' @rdname ResultSet-class
#' @exportClass ResultSet
#' @slot fun_origin Character containing the function that creates the object.
# @slot class_origin Character containing the class of the objects used to
# create the object.
# @slot names Character with the names of the datasets in the
# \link{MultiDataSet} used to create the object
#' @slot results List containing the results of the association/integration.
#' @slot fData List containing the feature-data of the original objects.
#' @slot options list of options used to create the \code{ResultSet}.
#' @return An object of class \code{ResultSet}
setClass(
  Class = "ResultSet",
  representation = representation(
    fun_origin = "character",
    results = "list",
    fData = "list",
    options = "list"
  )
)
