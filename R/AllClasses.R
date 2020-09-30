setClassUnion("characterOrMissing", c("character", "missing", "logical"))
setClassUnion("factorOrcharacter", c("factor", "character"))
setClassUnion("integerOrMissing", c("integer", "missing", "logical"))


.valid.AlignmentPairs <- function(object) {
    if (length(object)) {
        if (!("query" %in% colnames(elementMetadata(object))))
            return(paste0("AlignmentPairs must have a 'query' column"))
        if (!("subject" %in% colnames(elementMetadata(object))))
            return(paste0("AlignmentPairs must have a 'subject' column"))
    }
}

##' Representation of an alignment pair
##'
##' @description Pairs subclass
##'
##'
##' @details The AlignmentPairs class extends the
##'     \code{\link[S4Vectors]{Pairs-class}}
##'
##' @export
##' @rdname AlignmentPairs-class
##'
##' @import S4Vectors
##'
##' @seealso \code{\link[S4Vectors]{Pairs-class}}
##'
setClass("AlignmentPairs",
         contains = c("Pairs"),
         validity = .valid.AlignmentPairs)


##' @importFrom methods validObject callNextMethod
setMethod("initialize", "AlignmentPairs", function(.Object, ...) {
    .Object <- callNextMethod()
    validObject(.Object)
    .Object
})


##' List of AlignmentPairs instances
##'
##' @description Subclass of S4Vectors SimpleList, where each entry is
##'     an AlignmentPairs object
##'
##' @export
##' @rdname AlignmentPairsList-class
##'
##' @import S4Vectors
##'
setClass("AlignmentPairsList",
         contains = "SimpleList",
         prototype = prototype(elementType = "AlignmentPairs")
         )
