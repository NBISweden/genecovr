##' methods-sequence-data
##'
##' @description Methods for dealing with sequence data
##'
##' @details
##'
##' @name sequence data methods
##'
NULL


##' readFastaIndex
##'
##' @description Read a fasta index file
##'
##' @param fai fasta index file
##' @param ... additional parameters to pass to GenomeInfoDb::Seqinfo
##'
##' @return Seqinfo object
##'
##' @import GenomeInfoDb
##'
##' @export
##' @rdname readFastaIndex
##'
##' @examples
##' fai <- system.file("extdata", "polished.fai",
##'        package="genecovr")
##' sinfo <- readFastaIndex(fai)
##'
readFastaIndex <- function(fai, ...) {
    data <- read.table(fai, header = FALSE,
                       col.names = c("NAME", "LENGTH", "OFFSET",
                                     "LINEBASES", "LINEWIDTH"),
                       as.is = TRUE)
    GenomeInfoDb::Seqinfo(data$NAME, data$LENGTH, ...)
}
