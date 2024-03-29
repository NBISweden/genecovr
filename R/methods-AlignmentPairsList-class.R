##' AlignmentPairsList
##'
##' @export
##' @rdname AlignmentPairsList-class
##'
##' @importFrom methods new
##'
setMethod("AlignmentPairsList", "list",
          function(obj) new("AlignmentPairsList", listData = obj))


##' as.data.frame
##'
##' @description Convert AlignmentPairsList to data.frame.
##'
##' @param x AlignmentPairsList object
##' @param ... additional parameters to lapply
##' @param .id name of id column added by dplyr::bind_rows
##'
##' @return data.frame
##'
##' @importFrom dplyr bind_rows
##'
##' @export
##'
setMethod("as.data.frame", signature = "AlignmentPairsList",
          function(x, ..., .id="id") {
    data <- bind_rows(lapply(x, as.data.frame, ...), .id = .id)
    if (!is.null(names(x)))
        data[[.id]] <- factor(data[[.id]], levels=names(x))
    data
})



##' insertionSummary
##'
##' @description Summarize insertions
##'
##' @details Summarize number of insertions
##'
##' @param x AlignmentPairsList object
##' @param reduce reduce data to transcripts
##' @param bpparam biocparallel parameters
##' @param ... additional parameters
##'
##' @return data.frame
##'
##'
##' @export
##' @rdname insertionSummary
##'
##' @examples
##' fn <- system.file("extdata", "transcripts2polished.psl",
##'       package="genecovr")
##' x <- AlignmentPairsList(list(readPsl(fn)))
##' insertionSummary(x)
##'
insertionSummary <- function(x, reduce=TRUE, bpparam=NULL, ...) {
    stopifnot(inherits(x, "AlignmentPairsList"))
    if (reduce) {
        x <- geneBodyCoverage(x, bpparam=bpparam, ...)
    } else {
        x <- as.data.frame(x)
        colnames(x) <- gsub("query.", "", colnames(x))
    }
    x$cuts <- cut(x$NumInsert, c(-1:3, Inf), include.lowest=FALSE)
    levels(x$cuts) <- c(0:3, ">3")
    x$cuts <- factor(x$cuts, levels=c(">3", 3:0))
    x
}

##' geneBodyCoverage
##' @rdname geneBodyCoverage
##'
##' @param bpparam BiocParallel parameter object
##'
##' @importFrom BiocParallel bplapply
##' @importFrom dplyr bind_rows
##'
setMethod("geneBodyCoverage", signature = c("AlignmentPairsList"),
          definition = function(obj, min.match=0.9, bpparam=NULL) {
    if (is.null(bpparam)) {
        gbc <- lapply(obj, function(x) {geneBodyCoverage(x, min.match=min.match)})
    } else {
        gbc <- BiocParallel::bplapply(obj, function(x) {geneBodyCoverage(x, min.match=min.match)},
                                      BPPARAM=bpparam)
    }
    x <- bind_rows(lapply(gbc, data.frame), .id="id")
    if (!is.null(names(obj)))
        x[["id"]] <- factor(x[["id"]], levels=names(obj))
    x
})


##' autoplot.AlignmentPairsList
##'
##' @importFrom ggplot2 autoplot ggplot geom_point geom_boxplot
##'     geom_violin geom_density facet_wrap geom_bar
##' @importFrom rlang enquo quo_set_expr as_label get_expr
##' @importFrom tidyr pivot_longer
##'
##' @param object AlignmentPairsList
##' @param aes aes mapping
##' @param vars variable mapping to facet plots
##' @param ... additional parameters to ggplot function
##' @param which which plot to make. 'grid' option makes a scatter
##'     plot with marginal densities
##'
##' @export
##'
autoplot.AlignmentPairsList <- function(object, aes, vars, ..., which="point") {
    data <- as.data.frame(object)
    n <- 1
    if (grepl("enquo", as_label(aes$y)))
        n <- length(eval(get_expr(aes$y)))
    if (n > 1) {
        data <- pivot_longer(data, cols = !!aes$y)
        aes$y <- quo_set_expr(aes$y, as.symbol("value"))
    }
    p <- ggplot(data, {{ aes }})
    which <- match.arg(which, c("point", "boxplot", "violin", "density"))
    if (which == "point")
        p <- p + geom_point(...)
    else if (which == "boxplot")
        p <- p + geom_boxplot(...)
    else if (which == "violin")
        p <- p + geom_violin(...)
    else if (which == "density")
        p <- p + geom_density(...)
    else if (which == "bar")
        p <- p + geom_bar(...)
    if (!missing(vars))
        p <- p + facet_wrap({{ vars }})
    p
}


##' plot
##'
##' @description plot an AlignmentPairsList
##'
##' @param x object to plot
##' @param ... additional arguments for autoplot
##'
##' @export
##' @importFrom graphics plot
##'
plot.AlignmentPairsList <- function(x, ...) {
    print(autoplot(x, ...))
}
