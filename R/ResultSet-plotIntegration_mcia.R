.plot_integration_mcia <- function(object, cmpX, cmpY, tcolors, ...) {
    if(object@results[[1]]$result$coa[[1]]$nf < 2) {
        stop("Input ResultSet was created for less than 2 components.")
    }

    if(sum(c(cmpX, cmpY) <= object@results[[1]]$result$coa[[1]]$nf) != 2) {
        stop("Given component X or component Y (cmpX, cmpY) higher than number of axis in ResultSet.")
    }


    plot_mcia(object@results[[1]]$result, cmpX, cmpY, colors=tcolors)
}

#colors <- c("red", "blue", "yellow")
#names(colors) <- c("set1", "set2", "set3")


splot_variables <- function (x, axis1 = 1, axis2 = 2, colors) {
    co <- x$mcoa$Tco[, c(axis1, axis2)]
    c <- as.numeric(x$mcoa$TC$T)
    c <- names(colors)[c]
    co <- cbind(co, color=c, feature=rownames(co))
    colnames(co)[1:2] <- c("axis1", "axis2")

    ggplot2::ggplot(co, ggplot2::aes_string(x="axis1", y="axis2")) +
        ggplot2::geom_point(shape=19, alpha=0.5, ggplot2::aes_string(color="color")) +
        ggplot2::scale_color_manual(values=colors, name="Data Set") +
        ggplot2::theme(legend.position="right") +
        ggplot2::xlab("") +
        ggplot2::ylab("") +
        ggplot2::ggtitle("Feaure Space") +
        ggplot2::theme_bw()
}

splot_samples <- function (x, axis1 = 1, axis2 = 2, colors) {
    dfxy2 <- x$mcoa$Tl1
    dfxy2$ID <- sapply(strsplit(rownames(dfxy2), "\\."), "[[", 1)
    dfxy2$label <- dfxy2$ID
    dfxy2[duplicated(dfxy2$label), "label"] <- ""
    dfxy2$color <- names(colors)[as.numeric(gsub("df", "", sapply(strsplit(rownames(dfxy2), "\\."), "[[", 2)))]
    colnames(dfxy2)[1:2] <- c("axis1", "axis2")

    ggplot2::ggplot(dfxy2, ggplot2::aes_string(x="axis1", y="axis2", group="ID")) +
        ggplot2::geom_point(shape=19, size=2, ggplot2::aes_string(color="color")) +
        ggplot2::scale_color_manual(values=colors) +
        ggplot2::geom_line(color="DarkGray") +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position="none") +
        ggplot2::xlab("") +
        ggplot2::ylab("") +
        ggplot2::ggtitle("Sample Space") +
        ggplot2::geom_text(ggplot2::aes_string(label="label"))
}



plot_mcia <- function (mcoin, cmpX=1, cmpY=2, colors) {
    eig <- mcoin$mcoa$pseudoeig
    eig <- data.frame(
        Component=paste("Comp ", 1:length(eig)),
        Explained=eig
    )

    eigen_plot <- ggplot2::ggplot(eig, ggplot2::aes_string(x="Component", y="Explained")) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::xlab("Eigen Value") +
        ggplot2::ylab("") +
        ggplot2::theme_bw()

    cov2 <- as.data.frame(mcoin$mcoa$cov2[, c(cmpX, cmpY)])
    colnames(cov2) <- paste0("ax", 1:2)
    cov2$label <- names(colors)

    set_plot <- ggplot2::ggplot(cov2, ggplot2::aes_string(x="ax1", y="ax2")) +
        ggplot2::geom_point(shape=19, size=5, ggplot2::aes_string(color="label")) +
        #ggplot2::geom_text(ggplot2::aes_string(label="label")) +
        ggplot2::xlab(paste("pseudoeig", cmpX)) +
        ggplot2::ylab(paste("pseudoeig", cmpY)) +
        ggplot2::scale_color_manual(values=colors) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position="none")

    gridExtra::grid.arrange(
        splot_samples(mcoin, axis1 = cmpX, axis2 = cmpY, colors = colors),
        splot_variables(mcoin, axis1 = cmpX, axis2 = cmpY, colors = colors),
        eigen_plot,
        set_plot,
        ncol = 2
    )
}
