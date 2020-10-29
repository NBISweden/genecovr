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
