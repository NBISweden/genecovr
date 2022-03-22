##' @section Constructor:
##'
##' \code{AlignmentPairs(query, subject, ...)}: Constructs an AlignmentPairs
##' object by aligning the vectors \code{query} and \code{subject}.
##'
##' @param query vectorized object representing query
##' @param subject vectorized object representing subject
##' @param ... Arguments to pass to constructor
##'
##' @return AlignmentPairs object
##'
##' @export
##' @rdname AlignmentPairs-class
##'
setGeneric("AlignmentPairs", signature = c("query", "subject"),
           function(query, subject, ...)
    standardGeneric("AlignmentPairs"))


##' Retrieve score from object
##'
##' @description For any object with a score slot retrieve the score
##'     object.
##'
##' @param x object representing an alignment
##' @param ... Additional arguments
##'
##' @return numeric
##'
##' @export
##' @rdname score-methods
##'
##'
setGeneric("score", function(x, ...) standardGeneric("score"))

##' Retrieve query from object
##'
##' @description For any object with a query slot retrieve the query
##'     object.
##'
##' @param x object representing an alignment
##' @param ... Additional arguments
##'
##' @return query object
##'
##' @export
##' @rdname query-methods
##'
setGeneric("query", function(x, ...) standardGeneric("query"))

##' Set the query of an object
##'
##' @description For any object with a query slot set the
##'     query.
##'
##' @param x object representing an alignment
##' @param value object representing a query
##'
##' @return object representing an alignment with an updated query slot
##'
##' @export
##' @rdname query-methods
##'
setGeneric("query<-", function(x, value) standardGeneric("query<-"))

##' Retrieve subject from object
##'
##' @description For any object with a subject slot retrieve the
##'     subject object.
##'
##' @param x object representing an alignment
##' @param ... Additional arguments
##'
##' @return object representing a subject
##'
##' @export
##' @rdname sbjct-methods
##'
setGeneric("sbjct", function(x, ...) standardGeneric("sbjct"))

##' Set the subject of an object
##'
##' @description For any object with a subject slot set the
##'     subject.
##'
##' @param x object representing an alignment
##' @param value object representing a subject
##'
##' @return object representing an alignment with an updated subject
##'     slot
##'
##' @export
##' @rdname sbjct-methods
##'
setGeneric("sbjct<-", function(x, value) standardGeneric("sbjct<-"))

##' Get divergence from an alignment
##'
##' @description Get divergence from an alignment.
##'
##' @param x Object representation of a sequence alignment
##' @param ... Additional argument list
##'
##' @return numeric
##'
##' @export
##' @rdname divergence-methods
##'
##'
setGeneric("divergence", function(x, ...) standardGeneric("divergence"))

##' Get deletions from an alignment.
##'
##' @description Get deletions from an alignment.
##'
##' @param x Object representation of a sequence alignment
##' @param ... Additional argument list
##'
##' @return numeric
##'
##' @export
##' @docType methods
##' @rdname deletions-methods
##'
setGeneric("deletions", function(x, ...) standardGeneric("deletions"))

##' Get insertions from an alignment
##'
##' @description Get insertions from an alignment.
##'
##' @param x Object representation of a sequence alignment
##' @param ... Additional argument list
##'
##' @return numeric
##'
##' @export
##' @docType methods
##' @rdname insertions-methods
##'
setGeneric("insertions", function(x, ...) standardGeneric("insertions"))

##' matches
##'
##' @description Get matches from an alignment.
##'
##' @param x Object representation of a sequence alignment
##' @param ... Additional argument list
##'
##' @return numeric
##'
##' @export
##' @rdname matches-methods
##'
setGeneric("matches", function(x, ...) standardGeneric("matches"))

##' Get mismatches from an alignment
##'
##' @description Get mismatches from an alignment.
##'
##' @param x Object representation of a sequence alignment
##' @param ... Additional argument list
##'
##' @return numeric
##' @export
##' @rdname matches-methods
##'
setGeneric("mismatches", function(x, ...) standardGeneric("mismatches"))

##' Get matches that are parts of repeats from an alignment
##'
##' @description Get matches that are parts of repeats from an alignment.
##'
##' @param x Object representation of a sequence alignment
##' @param ... Additional argument list
##'
##' @return numeric
##'
##' @export
##' @rdname matches-methods
##'
setGeneric("repmatches", function(x, ...) standardGeneric("repmatches"))


##' @section Constructor:
##'
##' \code{AlignmentPairsList(obj, ...)}: Constructs an
##' AlignmentPairsList object from the input
##'
##' @param obj list
##' @param ... ellipsis
##'
##' @export
##' @rdname AlignmentPairsList-class
##'
##' @return AlignmentPairsList object
##'
setGeneric("AlignmentPairsList", signature = c("obj"),
           function(obj, ...)
    standardGeneric("AlignmentPairsList"))

##' Get gene body coverage
##'
##' \code{geneBodyCoverage(obj, ...)}: calculate gene body coverage from input
##'
##' @param obj alignment object
##' @param ... additional parameters
##'
##'
##' @export
##' @rdname geneBodyCoverage
##'
setGeneric("geneBodyCoverage", signature = c("obj"),
           function(obj, ...)
    standardGeneric("geneBodyCoverage"))
