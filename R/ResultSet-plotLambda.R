#' @aliases plotLambda
#' @rdname plotLambda-methods
setMethod(
    f = "plotLambda",
    signature = "ResultSet",
    definition = function(object, width=0.75) {
        tt <- tableLambda(object)

        ggplot2::ggplot(tt, ggplot2::aes_string(x="exposure", y="lambda")) +
            ggplot2::geom_bar(stat="identity", width=width) +
            ggplot2::ylab("lambda score") + ggplot2::xlab("") +
            ggplot2::theme(
                axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)
            ) + ggplot2::theme_bw()
    }
)
