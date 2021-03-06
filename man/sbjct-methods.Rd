% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/methods-AlignmentPairs-class.R
\name{sbjct}
\alias{sbjct}
\alias{sbjct<-}
\alias{sbjct,AlignmentPairs-method}
\alias{sbjct<-,AlignmentPairs,GRanges-method}
\title{Retrieve subject from object}
\usage{
sbjct(x, ...)

sbjct(x) <- value

\S4method{sbjct}{AlignmentPairs}(x)

\S4method{sbjct}{AlignmentPairs,GRanges}(x) <- value
}
\arguments{
\item{x}{object representing an alignment}

\item{...}{Additional arguments}

\item{value}{object representing a subject}
}
\value{
object representing a subject

object representing an alignment with an updated subject
slot
}
\description{
For any object with a subject slot retrieve the
subject object.

For any object with a subject slot set the
subject.
}
