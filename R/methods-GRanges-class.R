##' @rdname insertions-methods
##' @aliases insertions,GRanges,ANY
##'
##' @param as.bases Get insertions as base pairs
##'
setMethod("insertions", "GRanges",
          function(x, as.bases=FALSE, ...) {
    vnames <- names(values(x))
    retval <- NULL
    if (as.bases) {
        if ("BaseInsert" %in% vnames)
            retval <- x$BaseInsert
    } else {
        if ("NumInsert" %in% vnames)
            retval <- x$NumInsert
    }
    retval
})
