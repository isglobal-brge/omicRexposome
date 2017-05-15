.plot_assoc <- function(object, rid, coef, contrast, type, tPV, tFC, effect) {
    dta <- topTable(object, rid=rid, coef=coef, contrast=contrast)

    if(type == "qq") {
        qq_plot(topTable(object, rid=rid, coef=coef)$P.Value)
    } else if(type == "manhattan") {
        if("MethylationSet" %in% object@class_origin) {
            dta <- dta[ , c("P.Value", "chromosome", "position")]
        } else {
            dta <- dta[ , c("P.Value", "chromosome", "start")]
        }
        dta$SNP <- rownames(object@results[[rid]]$result)
        colnames(dta) <- c("P", "CHR", "BP", "SNP")
        dta$CHR <- gsub("chr", "", sapply(strsplit(dta$CHR, "_"), "[[", 1))
        dta$CHR <- gsub("X", "23", gsub("Y", "24", dta$CHR))
        dta <- dta[dta$CHR %in% as.character(1:24), ]
        dta$CHR <- as.numeric(dta$CHR)
        dta$BP <- as.numeric(dta$BP)

        qqman::manhattan(dta, ylab="-log10(P.Value)", ...)
    } else if(type == "volcano") {
        volcano_plot(
            pval=dta$P.Value,
            fc=dta$logFC,
            names=rownames(dta),
            tFC=tFC,
            tPV=tPV,
            show.effect=effect
        )
    } else if(type == "protein") {
        dta$PV <- -log10(dta$P.Value)
        dta$protein <- rownames(dta)
        if(nrow(dta) == 1) {
            stop("Invalid data obtained from 'topTable.'")
        }
        ggplot2::ggplot(dta,
            ggplot2::aes_string(x="protein", y="PV", size="PV", fill="PV", color="PV")) +

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
        stop("Invalid type of plot ('", type, "').")
    }
}

