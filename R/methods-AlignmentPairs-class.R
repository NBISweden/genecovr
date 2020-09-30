##' AlignmentPairs
##'
##' @description Constructor for AlignmentPairs
##'
##' @export
##' @rdname AlignmentPairs-class
##'
##' @importFrom GenomicRanges GRanges
##' @importFrom S4Vectors DataFrame
##'
setMethod("AlignmentPairs", signature = c("GRanges", "GRanges"),
          definition = function(query, subject, ...) {
    if (!missing(...)) {
        elementMetadata <- DataFrame(...)
    } else {
        elementMetadata <- new("DataFrame", nrows = length(query))
    }
    elementMetadata$query <- query
    elementMetadata$subject <- subject
    new("AlignmentPairs",
        first = seq_along(query),
        second = seq_along(subject),
        elementMetadata = elementMetadata)
})

##############################
## Getters
##############################

##'
##' @export
##' @rdname query-methods
##'
setMethod("query", "AlignmentPairs", function(x) mcols(x)$query)

##'
##' @export
##' @rdname query-methods
##'
##' @importFrom methods validObject
##'
setMethod("query<-",
          signature = c("AlignmentPairs", "GRanges"),
          function(x, value) {
    mcols(x)$query <- value
    validObject(x)
    x
})

##'
##' @export
##' @rdname sbjct-methods
##'
setMethod("sbjct", "AlignmentPairs", function(x) mcols(x)$subject)

##'
##' @export
##' @rdname sbjct-methods
##'
##' @importFrom methods validObject
##'
setMethod("sbjct<-",
          signature = c("AlignmentPairs", "GRanges"),
          function(x, value) {
    mcols(x)$subject <- value
    validObject(x)
    x
})


##'
##' @export
##' @rdname matches-methods
##'
setMethod("matches", "AlignmentPairs",
          function(x) {
    retval <- NULL
    if ("matches" %in% names(values(x)))
        retval <- values(x)$matches
    retval
})

##'
##' @export
##' @rdname matches-methods
##'
setMethod("mismatches", "AlignmentPairs",
          function(x) {
    retval <- NULL
    if ("mismatches" %in% names(values(x)))
        retval <- values(x)$mismatches
    else if ("misMatches" %in% names(values(x)))
        retval <- values(x)$misMatches
    retval
})

##'
##' @export
##' @rdname matches-methods
##'
setMethod("repmatches", "AlignmentPairs",
          function(x) {
    retval <- NULL
    if ("repmatches" %in% names(values(x)))
        retval <- values(x)$repmatches
    else if ("repMatches" %in% names(values(x)))
        retval <- values(x)$repMatches
    retval
})



##' reduceHitCoverage
##'
##' @description Reduce hit coverage for a seqname
##'
##' @details Given an AlignmentPairs object, reduce overlapping query
##'     regions. For instance, if one hit spans coordinates 100-400
##'     and another 200-500, we merge to 100-500. The function
##'     calculates how many overlaps that are merged, the rationale
##'     being that if there are >1 overlaps for a seqname, the
##'     sequence is associated with multiple distinct regions in the
##'     subject sequence. Thus, the unit of interest here is the
##'     query, and the function seeks to identify reduced disjoint
##'     regions of a query sequence that map to a subject.
##'
##'     The function returns a GRanges object with reduced regions,
##'     i.e. overlaps are merged. Information about matches and
##'     mismatches is currently dropped as it usually is not possible
##'     to infer where mismatches occur. Instead, the data column
##'     `coverage` simply holds the ratio of the width of the region
##'     to the transcript length. Summing up coverages from disjoint
##'     regions then gives total coverage of the transcript.
##'
##'     The cutoff is used to filter regions based on the ratio of
##'     matches to the width of the region.
##'
##'     Since the association between query and subject regions is
##'     removed, the return value is a GRanges object consisting of
##'     the reduced query ranges with a revmap and coverage attribute.
##'
##' @rdname reduceHitCoverage
##' @export
##'
##' @importFrom IRanges CharacterList reduce width
##' @importFrom S4Vectors subset
##'
##' @param x AlignmentPairs object
##' @param min.match filter out hits with fraction of matching bases
##'     less than min.match
##'
##' @return reduced and filtered GRanges object
##'
##' @examples
##' ranges <- IRanges::IRanges(
##'           start=c(100, 200, 700),
##'           end=c(400, 500, 1000)
##' )
##' qry <- GenomicRanges::GRanges(
##'           ranges=ranges,
##'           seqnames=c("t1"),
##'           seqinfo=GenomeInfoDb::Seqinfo(seqnames=c("t1"),
##'                                         seqlengths=c(1050))
##' )
##' ranges <- IRanges::IRanges(
##'           start=c(1000, 1000, 1000),
##'           end=c(1300, 1300, 1300)
##' )
##' sbj <- GenomicRanges::GRanges(ranges=ranges, seqnames=c("c1", "c2", "c1"))
##' x <- AlignmentPairs(query=qry, subject=sbj, matches=c(300, 290, 280))
##' gr <- reduceHitCoverage(x, 0.1)
##'
reduceHitCoverage <- function(x, min.match=0.9) {
    .Deprecated()
    stopifnot(inherits(x, "AlignmentPairs"))
    message("Reducing hits for ", summary(x), ", min.match ", min.match)
    y <- subset(x, matches(x) / width(query(x)) > min.match)
    q <- reduce(query(y), ignore.strand = TRUE,
                with.revmap = TRUE)
    q
}

##' breadthOfCoverage
##'
##' @description Calculate breadth of coverage by sequence id
##'
##' @param x GRangesList object
##'
##' @rdname breadthOfCoverage
##' @export
##'
##' @importFrom GenomicRanges GRangesList
##' @importFrom IRanges width reduce
##'
##' @return integer representing the sum of range widths
##'
## see https://stackoverflow.com/questions/39401376/aggregate-bins-in-large-granges-efficiently
breadthOfCoverage <- function(x) {
    stopifnot(inherits(x, "GRangesList"))
    sum(width(reduce(x)))
}

##' depthOfCoverage
##'
##' @description Calculate depth of coverage by sequence id
##'
##' @param x GRangesList object
##'
##' @rdname depthOfCoverage
##' @export
##'
##' @importFrom GenomicRanges GRangesList
##' @importFrom IRanges width
##'
##' @return integer representing the sum of range widths
##'
## see https://stackoverflow.com/questions/39401376/aggregate-bins-in-large-granges-efficiently
depthOfCoverage <- function(x) {
    stopifnot(inherits(x, "GRangesList"))
    sum(width(x))
}


##' geneBodyCoverage
##'
##' @description Calculate geneBodyCoverage of query sequences
##'
##' @details Given an AlignmentPairs object, calculate gene body
##'     coverages of query sequences.
##'     regions. For instance, if one hit spans coordinates 100-400
##'     and another 200-500, we merge to 100-500. The function
##'     calculates how many overlaps that are merged, the rationale
##'     being that if there are >1 overlaps for a seqname, the
##'     sequence is associated with multiple distinct regions in the
##'     subject sequence. Thus, the unit of interest here is the
##'     query, and the function seeks to identify reduced disjoint
##'     regions of a query sequence that map to a subject.
##'
##'     The function returns a GRanges object with reduced regions,
##'     i.e. overlaps are merged. Information about matches and
##'     mismatches is currently dropped as it usually is not possible
##'     to infer where mismatches occur. Instead, the data column
##'     `coverage` simply holds the ratio of the width of the region
##'     to the transcript length. Summing up coverages from disjoint
##'     regions then gives total coverage of the transcript.
##'
##'     The cutoff is used to filter regions based on the ratio of
##'     matches to the width of the region.
##'
##'     Since the association between query and subject regions is
##'     removed, the return value is a GRanges object consisting of
##'     the reduced query ranges with a revmap and coverage attribute.
##'
##'
##' @rdname geneBodyCoverage
##' @export
##'
##' @param x AlignmentPairs object
##' @param min.match filter out hits with fraction matching bases less
##'     than min.match
##'
##' @importFrom IRanges IntegerList
##' @importFrom S4Vectors DataFrame
##'
##' @return reduced and filtered DataFrame object where each row is a
##'     transcript. See @details for information on seqnames and
##'     seqlengths columns are collected from the corresponding
##'     seqinfo object. breadthOfCoverage is the total width of all
##'     reduced regions of the transcript and corresponds to how much
##'     of a transcript is covered. Dividing breadthOfCoverage by the
##'     total transcript length gives the coverage fraction.
##'     depthOfCoverage sums all ranges and is an estimate of how many
##'     times each range is present. Dividing the depthOfCoverage by
##'     breadthOfCoverage indicates the multiplicity of the query. For
##'     queries that map multiple subjects, a value close to 1
##'     indicates the query has been split in several subjects,
##'     whereas a higher value indicates sequence duplication at the
##'     subject level.
##'
##'     The revmap column maps the output ranges to the input ranges
##'     as lists of numerical ids. These ids can be used to retrieve
##'     the corresponding AlignmentPairs ranges providing a link to
##'     the subjects. The hitCoverage lists the width of each reduced
##'     hit, and hitStart and hitEnd provide the transcript
##'     coordinates of these hits.
##'
geneBodyCoverage <- function(x, min.match=0.9) {
    stopifnot(inherits(x, "AlignmentPairs"))
    message("Calculating gene body coverage for ",
            summary(x), ", min.match ", min.match)
    x <- subset(x, matches(x) / width(query(x)) > min.match)
    y <- reduce(query(x), with.revmap=TRUE, ignore.strand=TRUE)
    grl <- split(query(x), seqnames(query(x)))
    revmapToList <- function(x) unique(unlist(x$revmap))
    data <- DataFrame(
        seqnames = seqnames(seqinfo(query(x))),
        seqlengths = seqlengths(seqinfo(query(x))),
        breadthOfCoverage = breadthOfCoverage(grl),
        depthOfCoverage = depthOfCoverage(grl),
        revmap = IRanges::IntegerList(lapply(split(y, seqnames(y)), revmapToList)),
        hitCoverage = width(reduce(grl)),
        hitStart = start(grl),
        hitEnd = end(grl)
    )
    data$coverage <- data$breadthOfCoverage / data$seqlengths
    data$revmap.count <- unlist(lapply(data$revmap, length))
    data$n.subjects <- unlist(lapply(data$revmap, function(j) {length(seqnames(sbjct(x[j])))}))
    data
}


##' summarizeGeneBodyCoverage
##'
##' @rdname summarizeGeneBodyCoverage
##' @export
##'
##' @importFrom S4Vectors DataFrame
##'
##' @param x AlignmentPairs object
##' @param min.coverage coverage cutoffs to apply to transcripts
##' @param min.match filter out hits with fraction matching bases less
##'     than min.match
##'
summarizeGeneBodyCoverage <- function(x, min.coverage=seq(0, 1, 0.05),
                                      min.match=0.9) {
    stopifnot(inherits(x, "AlignmentPairs"))
    message("Summarizing gene body coverage for ",
            summary(x), ", ", length(min.coverage), " coverage cutoffs")
    data <- geneBodyCoverage(x, min.match)
    DataFrame(
        count = rev(cumsum(rev(apply(
            table(data$coverage, cut(data$coverage, min.coverage)), 2, sum)))),
        total = nrow(data),
        min.coverage = min.coverage[2:length(min.coverage)],
        min.match = min.match
    )
}


## FIXME: the method needs a rewrite!
##' countSubjectsByCoverage
##'
##' @rdname countSubjectsByCoverage
##' @export
##'
##' @importFrom S4Vectors DataFrame
##'
##' @param x AlignmentPairs object
##' @param min.coverage coverage cutoffs to apply to transcripts
##' @param min.match filter out hits with fraction matching bases less
##'     than min.match
##' @param nmax maximum number of subjects to count
##'
countSubjectsByCoverage <- function(x, min.coverage=seq(0.25, 1, 0.25),
                                    min.match=0.9, nmax=5) {
    stopifnot(inherits(x, "AlignmentPairs"))
    message("Counting number of subjects by coverage for ",
            summary(x), ", ", length(min.coverage), " coverage cutoffs")
    data <- geneBodyCoverage(x, min.match = min.match)
    tab <- table(factor(
        data$n.subjects,
        levels = sort(unique(c(0, data$n.subjects)))),
        cut(data$coverage, min.coverage))[, (length(min.coverage) - 1):1]
    z <- t(as.table(apply(tab, 1, cumsum)))
    z[1, ] <- nrow(data) - apply(z[2:nrow(z), ], 2, sum)
    z.df <- DataFrame(z)
    colnames(z.df) <- c("n.subjects", "min.coverage", "Freq")
    z.df$min.match <- min.match
    if (nmax < max(as.numeric(z.df$n.subjects)))
        nmax <- max(as.numeric(z.df$n.subjects))
    z.df$n.subjects <- factor(z.df$n.subjects, levels = c(0, nmax:1))
    levels(z.df$n.subjects)[1] <- "filtered"
    levels(z.df$min.coverage) <- rev(min.coverage[1:(length(min.coverage) - 1)])
    z.df$min.coverage <- factor(
        z.df$min.coverage,
        levels = min.coverage[1:(length(min.coverage) - 1)])
    z.df
}
