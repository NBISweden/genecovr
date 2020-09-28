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
                                                                 seqlengths = c(500, 250, 400, 500, 400))),
    matches = c(90, 90, 4, 100, 60, 160, 10))




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


withr::defer(rm(c(assembly_fai_fn, transcripts_fai_fn, psl_fn)))
