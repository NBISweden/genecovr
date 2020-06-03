##' methods-psl
##'
##' @description
##' Methods for dealing with psl data
##'
##' @details
##'
##' @name psl methods
##'
NULL

.header <- c(
    "matches", # Number of bases that match that aren't repeats
    "misMatches", # Number of bases that don't match
    "repMatches", # Number of bases that match but are part of repeats
    "nCount", # Number of "N" bases
    "qNumInsert", # Number of inserts in query
    "qBaseInsert", # Number of bases inserted in query
    "tNumInsert", # Number of inserts in target
    "tBaseInsert", # Number of bases inserted in target
    "strand", # "+" or "-" for query strand. For translated alignments, second "+"or "-" is for target genomic strand.
    "qName", # Query sequence name
    "qSize", # Query sequence size.
    "qStart", # Alignment start position in query
    "qEnd", # Alignment end position in query
    "tName", # Target sequence name
    "tSize", # Target sequence size
    "tStart", # Alignment start position in target
    "tEnd", # Alignment end position in target
    "blockCount", # Number of blocks in the alignment (a block contains no gaps)
    "blockSizes", # Comma-separated list of sizes of each block. If the query is a protein and the target the genome, blockSizes are in amino acids. See below for more information on protein query PSLs.
    "qStarts", # Comma-separated list of starting positions of each block in query
    "tStarts" # Comma-separated list of starting positions of each block in target
)


.create_seqinfo <- function(names, lengths, which="query") {
    message("manually inferring seqinfo for ", which)
    sinfo <- unique(data.frame(names = as.character(names),
                               lengths = lengths,
                               stringsAsFactors = FALSE))
    if (any(duplicated(sinfo$names)))
        sinfo <- sinfo[!duplicated(sinfo$names),]
    GenomeInfoDb::Seqinfo(sinfo$names, sinfo$lengths)
}

##' readPsl
##'
##' @description Read psl output
##'
##' @details
##'
##' @param filename input filename
##' @param seqinfo.query Seqinfo object for query (transcripts)
##' @param seqinfo.sbjct Seqinfo object for subject (reference)
##' @param metadata metadata for AlignmentPairs result
##' @param ... additional parameters
##'
##' @return AlignmentPairs object
##'
##' @importFrom GenomicRanges GRanges
##' @import S4Vectors
##' @import IRanges
##' @importFrom utils read.table
##'
##' @export
##' @rdname readPsl
##'
##' @examples
##' fn <- system.file("extdata", "transcripts2polished.psl",
##'       package="genecovr")
##' ap <- readPsl(fn)
##'
readPsl <- function(filename, seqinfo.query=NULL, seqinfo.sbjct=NULL, metadata=list(), ...) {
    message("reading file ", filename)
    start_time <- Sys.time()
    con <- file(filename, "r")
    on.exit(close(con))
    data <- read.table(con, header=FALSE)
    colnames(data) <- .header
    ## Create query and subject
    if (is.null(seqinfo.query))
        seqinfo.query <- .create_seqinfo(data$qName, data$qSize)
    if (is.null(seqinfo.sbjct))
        seqinfo.sbjct <- .create_seqinfo(data$tName, data$tSize, "subject")
    query <- GRanges(
        seqnames = S4Vectors::Rle(data$qName),
        ## NB: psl is 0-based
        ranges = IRanges::IRanges(start = data$qStart + 1,
                                  end = data$qEnd),
        strand = "+",
        seqinfo = seqinfo.query,
        NumInsert = data$qNumInsert,
        BaseInsert = data$qBaseInsert,
        blockCount = data$blockCount,
        blockSizes = data$blockSizes,
        blockStart = data$qStarts)
    sbjct <- GRanges(
        seqnames = S4Vectors::Rle(data$tName),
        ranges = IRanges::IRanges(start = data$tStart + 1,
                                  end = data$tEnd),
        strand = S4Vectors::Rle(data$strand),
        seqinfo = seqinfo.sbjct,
        NumInsert = data$tNumInsert,
        BaseInsert = data$tBaseInsert,
        blockCount = data$blockCount,
        blockSizes = data$blockSizes,
        blockStart = data$tStarts)
    ap <- AlignmentPairs(query, sbjct)
    ## Add data columns
    i.values <- match(c("matches", "misMatches", "repMatches", "nCount"), .header)
    values(ap)[, .header[i.values]] <- data[, .header[i.values]]
    message("Processed ", nrow(data), " lines in ",
            format(Sys.time() - start_time, digits = 2))
    ap
}
