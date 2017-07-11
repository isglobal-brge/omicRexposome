setMethod(
    f = "plotAssociation",
    signature = "ResultSet",
    definition = function(object, rid = 1, coef = 2, contrast = NULL, type,
            tPV, tFC, show.effect=FALSE) {
        if(type == "protein") {
            dta <- MultiDataSet::getAssociation(object,
                rid=rid, coef=coef, contrast=contrast)
            if(nrow(dta) == 1) {
                stop("Invalid data obtained from 'topTable.'")
            }
            dta$PV <- -log10(dta$P.Value)
            dta$protein <- rownames(dta)
            ggplot2::ggplot(dta,
                    ggplot2::aes_string(x="protein", y="PV", size="PV",
                    fill="PV", color="PV")) +
                ggplot2::theme_bw() +
                ggplot2::geom_point(alpha=0.7) +
                ggplot2::theme(
                    axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                    legend.position = "none"
                ) +
                ggplot2::xlab("") +
                ggplot2::ylab(expression(-log[10](P-Value))) +
                ggplot2::scale_colour_gradientn(
                    colours = c("darkgray", "darkblue"))

        } else {
            MultiDataSet::plot(object, rid = rid, coef = coef,
                contrast = contrast, type = type, tPV = tPV, tFC = tFC,
                show.effect = show.effect)
        }
    }
)
