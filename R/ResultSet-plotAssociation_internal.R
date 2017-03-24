.plot_assoc <- function(object, rid, coef, contrast, type, tPV, tFC, ...) {
    dta <- topTable(object, rid=rid, coef=coef, contrast=contrast)

    if(type == "qq") {
        qq_plot(object@results[[rid]]$result$P.Value)
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
            tPV=tPV
        )
    } else {
        stop("Invalid type of plot ('", type, "').")
    }
}

