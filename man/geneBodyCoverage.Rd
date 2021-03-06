% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods-AlignmentPairs-class.R
\name{geneBodyCoverage}
\alias{geneBodyCoverage}
\title{geneBodyCoverage}
\usage{
geneBodyCoverage(x, min.match = 0.9)
}
\arguments{
\item{x}{AlignmentPairs object}

\item{min.match}{filter out hits with fraction matching bases less
than min.match}
}
\value{
reduced and filtered DataFrame object where each row is a
transcript. See @details for information on seqnames and
seqlengths columns are collected from the corresponding
seqinfo object. breadthOfCoverage is the total width of all
reduced regions of the transcript and corresponds to how much
of a transcript is covered. Dividing breadthOfCoverage by the
total transcript length gives the coverage fraction.
depthOfCoverage sums all ranges and is an estimate of how many
times each range is present. Dividing the depthOfCoverage by
breadthOfCoverage indicates the multiplicity of the query. For
queries that map multiple subjects, a value close to 1
indicates the query has been split in several subjects,
whereas a higher value indicates sequence duplication at the
subject level.\preformatted{The revmap column maps the output ranges to the input ranges
as lists of numerical ids. These ids can be used to retrieve
the corresponding AlignmentPairs ranges providing a link to
the subjects. The hitCoverage lists the width of each reduced
hit, and hitStart and hitEnd provide the transcript
coordinates of these hits.
}
}
\description{
Calculate geneBodyCoverage of query sequences
}
\details{
Given an AlignmentPairs object, calculate gene body
coverages of query sequences.
regions. For instance, if one hit spans coordinates 100-400
and another 200-500, we merge to 100-500. The function
calculates how many overlaps that are merged, the rationale
being that if there are >1 overlaps for a seqname, the
sequence is associated with multiple distinct regions in the
subject sequence. Thus, the unit of interest here is the
query, and the function seeks to identify reduced disjoint
regions of a query sequence that map to a subject.\preformatted{The function returns a GRanges object with reduced regions,
i.e. overlaps are merged. Information about matches and
mismatches is currently dropped as it usually is not possible
to infer where mismatches occur. Instead, the data column
`coverage` simply holds the ratio of the width of the region
to the transcript length. Summing up coverages from disjoint
regions then gives total coverage of the transcript.

The cutoff is used to filter regions based on the ratio of
matches to the width of the region.

Since the association between query and subject regions is
removed, the return value is a GRanges object consisting of
the reduced query ranges with a revmap and coverage attribute.
}
}
