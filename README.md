
<!-- README.md is generated from README.Rmd with devtools::build_readme(). Please edit that file -->

# genecovr

<!-- badges: start -->

[![R build
status](https://github.com/NBISweden/genecovr/workflows/R-CMD-check/badge.svg)](https://github.com/NBISweden/genecovr/actions)
<!-- badges: end -->

Perform gene body coverage analyses in R to evaluate genome assembly
quality.

## Installation

You can install the released version of genecovr from [NBIS
GitHub](https://github.com/nbis) with:

``` r
# install.packages("devtools")
devtools::install_github("NBISweden/genecovr")
```

## Quick usage

There is a helper script for generating basic plots located in
PACKAGE\_DIR/bin/genecovr. Create a data input csv-delimited file with
columns

1.  data label
2.  mapping file (supported formats: psl)
3.  assembly file (fasta or fasta index)
4.  transcript file (fasta or fasta index)

Columns 3 and 4 can be set to missing value (NA) in which case sequence
sizes will be inferred from the alignment files. Then run the script to
generate plots:

``` shell
PACKAGE_DIR/bin/genecovr indata.csv
```

## Vignette

Alternatively, import the library as usual in an R script and use the
package functions. See the vignette for a minimum working example.
