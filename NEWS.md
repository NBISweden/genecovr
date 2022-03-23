# Release 0.0.0.9013

- fix factor level ordering for geneBodyCoverage plot
- save geneBodyCoverage as tsv

# Release 0.0.0.9012

- adjust factor levels for number of inserts (#4)
- summarize number of inserts by transcript (#5)

# Release 0.0.0.9011

- fix order of factors

# Release 0.0.0.9010

- remove duplicate entries in psl input

# Release 0.0.0.9009

- add plot of transcript length distributions conditioned on number of
  mapped contigs
- fix y axis for histogram plots
- increase point size of some plots
- summarizeGeneBodyCoverage and countSubjectsByCoverage now accept
  DataFrame inputs, obviating the need to rerun geneBodyCoverage
  multiple times in genecovr script


# Release 0.0.0.9008

- Remove characters trailing first space in fasta headers

# Release 0.0.0.9007

- Fix conversion of DNAStringSet to Seqinfo
- Make sure geneBodyCoverage table has nmax levels


# Release 0.0.0.9006

- add depthOfCoverage function and analysis to vignette and script
- reduceHitCoverage is deprecated
- improve some docs
- add wrapper for saving plots
- add tests mainly for alignmentpairs and test setup
