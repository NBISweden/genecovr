##################################################
## Small queries and subjects
##
## Cover as many corner cases as possible with simple examples
##################################################

##############################
## Multiple queries to one subject
##############################
ap1 <- AlignmentPairs(
    query = GenomicRanges::GRanges(
                               ranges = IRanges::IRanges(
                                                     start = c(1, 1, 90),
                                                     end = c(95, 98, 98)),
                               seqnames = c("t1.1", "t1.2", "t1.1"),
                               seqinfo = GenomeInfoDb::Seqinfo(seqnames = c("t1.1", "t1.2"),
                                                               seqlengths = c(100, 110))
                           ),
    subject = GenomicRanges::GRanges(
                                 ranges = IRanges::IRanges(
                                                       start = c(300, 300, 100),
                                                       end = c(400, 400, 200)),
                                 seqnames = c("ctg1.1", "ctg1.1", "ctg1.2"),
                                 seqinfo = GenomeInfoDb::Seqinfo(seqnames = c("ctg1.1", "ctg1.2"),
                                                                 seqlengths = c(500, 250))),
    matches = c(90, 90, 4))


##############################
## Multiple subjects to one query
##
## depthOfCoverage / breadthOfCoverage ~ 2 (i.e. duplication)
##############################
ap2 <- AlignmentPairs(
    query = GenomicRanges::GRanges(
                               ranges = IRanges::IRanges(
                                                     start = c(1, 3),
                                                     end = c(95, 98)),
                               seqnames = c("t2.1", "t2.1"),
                               seqinfo = GenomeInfoDb::Seqinfo(seqnames = c("t2.1"),
                                                               seqlengths = c(100))
                           ),
    subject = GenomicRanges::GRanges(
                                 ranges = IRanges::IRanges(
                                                       start = c(300, 300),
                                                       end = c(400, 400)),
                                 seqnames = c("ctg2.1", "ctg2.2"),
                                 seqinfo = GenomeInfoDb::Seqinfo(seqnames = c("ctg2.1", "ctg2.2"),
                                                                 seqlengths = c(500, 600))),
    matches = c(80, 90))


##############################
## Query split on different subjects
##
## depthOfCoverage / breadthOfCoverage ~ 1 (i.e. split)
##############################
ap3 <- AlignmentPairs(
    query = GenomicRanges::GRanges(
                               ranges = IRanges::IRanges(
                                                     start = c(1, 51),
                                                     end = c(54, 98)),
                               seqnames = c("t3.1", "t3.1"),
                               seqinfo = GenomeInfoDb::Seqinfo(seqnames = c("t3.1"),
                                                               seqlengths = c(100))
                           ),
    subject = GenomicRanges::GRanges(
                                 ranges = IRanges::IRanges(
                                                       start = c(100, 100),
                                                       end = c(160, 160)),
                                 seqnames = c("ctg3.1", "ctg3.2"),
                                 seqinfo = GenomeInfoDb::Seqinfo(seqnames = c("ctg3.1", "ctg3.2"),
                                                                 seqlengths = c(200, 200))),
    matches = c(50, 45))



##############################
## Slightly larger AlignmentPairs
##############################
ap4 <- AlignmentPairs(
    query = GenomicRanges::GRanges(
                               ranges = IRanges::IRanges(
                                                     start = c(1, 1, 90, 1, 115, 1, 190),
                                                     end = c(95, 98, 98, 120, 190, 170, 200)),
                               seqnames = c("t4.1", "t4.2", "t4.1", "t4.3", "t4.3", "t4.3", "t4.3"),
                               seqinfo = GenomeInfoDb::Seqinfo(seqnames = c("t4.1", "t4.2", "t4.3"),
                                                               seqlengths = c(100, 110, 200))
                           ),
    subject = GenomicRanges::GRanges(
                                 ranges = IRanges::IRanges(
                                                       start = c(300, 300, 100, 100, 220, 200, 300),
                                                       end = c(400, 400, 200, 225, 310, 380, 310)),
                                 seqnames = c("ctg4.1", "ctg4.1", "ctg4.2", "ctg4.3", "ctg4.3", "ctg4.4", "ctg4.5"),
                                 seqinfo = GenomeInfoDb::Seqinfo(seqnames = c("ctg4.1", "ctg4.2", "ctg4.3", "ctg4.4", "ctg4.5"),
                                                                 seqlengths = c(500, 250, 400, 500, 400)),
                                 NumInsert = c(),
                                 BaseInsert = c(),
                             ),
    matches = c(90, 90, 4, 100, 60, 160, 10))


##############################
## AlignmentPairs with numinserts
##############################
ap5 <- AlignmentPairs(
    query = GenomicRanges::GRanges(
                               ranges = IRanges::IRanges(
                                                     start = c(1, 81, 10, 90, 1, 50, 150),
                                                     end = c(80, 100, 100, 110, 100, 170, 200)),
                               seqnames = c("t5.1", "t5.1", "t5.2", "t5.2", "t5.3", "t5.3", "t5.3"),
                               seqinfo = GenomeInfoDb::Seqinfo(seqnames = c("t5.1", "t5.2", "t5.3"),
                                                               seqlengths = c(100, 110, 200)),
                               NumInsert = c(2, 0, 1, 4, 2, 3, 0),
                               BaseInsert = c(2, 0, 2, 6, 3, 3, 0),
                           ),
    subject = GenomicRanges::GRanges(
                                 ranges = IRanges::IRanges(
                                                       start = c(300, 401, 100, 200, 221, 200, 300),
                                                       end = c(380, 420, 190, 220, 320, 320, 350)),
                                 seqnames = c("ctg5.1", "ctg5.1", "ctg5.2", "ctg5.3", "ctg5.3", "ctg5.4", "ctg5.5"),
                                 seqinfo = GenomeInfoDb::Seqinfo(seqnames = c("ctg5.1", "ctg5.2", "ctg5.3", "ctg5.4", "ctg5.5"),
                                                                 seqlengths = c(500, 250, 400, 500, 400)),
                                 NumInsert = c(5, 1, 1, 5, 4, 3, 0),
                                 BaseInsert = c(6, 2, 1, 6, 4, 5, 0),
                             ),
    matches = c(78, 20, 90, 5, 97, 116, 50)
)




##############################
## AligmnentPairsList
##############################
apl <- AlignmentPairsList(list(ap1, ap2, ap3))

##################################################
## System files
##################################################
assembly_fai_fn <- list(
    nonpol = system.file("extdata", "nonpolished.fai", package = "genecovr"),
    pol = system.file("extdata", "polished.fai", package = "genecovr")
)
transcripts_fai_fn <- list(
    nonpol = system.file("extdata", "transcripts.fai", package = "genecovr"),
    pol = system.file("extdata", "transcripts.fai", package = "genecovr")
)

psl_fn <- list(
    nonpol = system.file("extdata", "transcripts2nonpolished.psl",
                         package = "genecovr"),
    pol = system.file("extdata", "transcripts2polished.psl",
                      package = "genecovr")
)


withr::defer(rm(assembly_fai_fn, transcripts_fai_fn, psl_fn))
