#' @aliases plotHits
#' @rdname plotHits-methods
setMethod(
    f = "plotHits",
    signature = "ResultSet",
    definition = function(object, th=0.05, width=0.75) {
        tt <- tableHits(object, th)
        ggplot2::ggplot(tt, ggplot2::aes_string(x="exposure", y="hits")) +
            ggplot2::geom_bar(stat="identity",  width=width) +
            ggplot2::ylab("#hits") + ggplot2::xlab("") +
            ggplot2::theme(
                axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)
            ) + ggplot2::theme_bw()
    }
)
