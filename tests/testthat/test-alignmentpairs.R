test_that("query returns query object correctly", {
    expect_equal(as.character(seqnames(query(ap1))), c("t1.1", "t1.2", "t1.1"))
})

test_that("matches returns correct matches", {
    expect_equal(matches(ap1), c(90, 90, 4))
})

## TODO: add function for APL
test_that("breadthOfCoverage calculates breadth of coverage by sequence identifier", {
    grl <- GRangesList(lapply(apl, query))
    ## NB: Will  fail when we add reduce option
    expect_equal(breadthOfCoverage(grl), c(202, 191, 102))
})

test_that("revmapList returns unique list of revmap entries", {
    y <- reduceHitCoverage(ap1, 0.1)
    grl <- split(y, seqnames(y))
    expect_equal(as.vector(unlist(lapply(grl, revmapList))), c(1, 3, 2))
})

test_that("reduceHitCoverage redoces hit coverage where applicable", {
    expect_equal(width(reduceHitCoverage(ap1)), c(95, 98))
    expect_equal(width(reduceHitCoverage(ap1, 0.1)), c(98, 98))
    expect_equal(length(reduceHitCoverage(ap1, 0.99)), 0)
})

test_that("geneBodyCoverage transcripts are equal to query seqinfo", {
    expect_equal(geneBodyCoverage(ap1)$seqnames, seqnames(seqinfo(query(ap1))))
})

test_that("geneBodyCoverage calculates correct coverages", {
    y <- geneBodyCoverage(ap1)
    expect_equal(as.vector(y$coverage), c(0.95, 0.89090909))
    expect_equal(as.vector(y$revmap.count), c(1, 1))
    y <- geneBodyCoverage(ap1, 0.1)
    expect_equal(as.vector(y$coverage), c(0.98, 0.89090909))
    expect_equal(as.vector(y$revmap.count), c(2, 1))
})

test_that("summarizeGeneBodyCoverage calculates correct summaries", {
    y <- summarizeGeneBodyCoverage(ap1, min.match = 0.1)
    expect_equal(as.vector(tail(y, 4)$count), c(2, 2, 1, 1))
    y <- summarizeGeneBodyCoverage(ap1)
    expect_equal(as.vector(tail(y, 4)$count), c(2, 2, 1, 0))
    y <- summarizeGeneBodyCoverage(ap1, min.match = 0.94)
    expect_equal(as.vector(tail(y, 4)$count), c(1, 1, 1, 0))
})

test_that("countSubjectsByCoverage generates correct table", {
    y <- countSubjectsByCoverage(ap4)
    expect_equal(levels(y$n.subjects), as.character(c("filtered", 5:1)))
    expect_equal(subset(y, n.subjects == 1)$Freq, c(2, 2, 2))
    expect_equal(subset(y, n.subjects == "filtered")$Freq, c(0, 0, 0))
    expect_equal(subset(y, n.subjects == 2)$Freq, c(1, 1, 1))
    expect_equal(subset(y, n.subjects == 3)$Freq, integer(0))
    y <- countSubjectsByCoverage(ap4, min.coverage=c(0.9, 0.95, 1.0))
    expect_equal(subset(y, n.subjects == "filtered")$Freq, c(3, 1))
})
