#' @describeIn ResultSet Draws a bar plot with the lambda value of each
#' analyses stored in the \code{ResultSet}.
setMethod(
    f = "plotLambda",
    signature = "ResultSet",
    definition = function(object, width=0.75) {
        tt <- tableLambda(object)

        ggplot2::ggplot(tt, ggplot2::aes(x=exposure, y=lambda, width=width)) +
            ggplot2::geom_bar(stat="identity") +
            ggplot2::ylab("lambda score") + ggplot2::xlab("") +
            ggplot2::theme(
                axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)
            ) + ggplot2::theme_bw()
    }
)
