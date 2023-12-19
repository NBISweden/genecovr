
<!-- README.md is generated from README.Rmd with devtools::build_readme(). Please edit that file -->

# genecovr

<!-- badges: start -->

[![R build
status](https://github.com/NBISweden/genecovr/workflows/R-CMD-check/badge.svg)](https://github.com/NBISweden/genecovr/actions)
<!-- badges: end -->

`genecovr` is an `R` package that provides plotting functions that
summarize gene transcript to genome alignments. The main purpose is to
assess the effect of polishing and scaffolding operations has on the
quality of a genome assembly. The gene transcript set is a large
sequence set consisting of assembled transcripts from RNA-seq data
generated in relation to a genome assembly project. Therefore,
`genecovr` serves as a complement to software such as
[BUSCO](https://busco.ezlab.org/), which evaluates genome assembly
quality using a smaller set of well-defined single-copy orthologs.

## Installation

You can install the released version of genecovr from [NBIS
GitHub](https://github.com/nbis) with:

``` r
# If necessary, uncomment to install devtools
# install.packages("devtools")
devtools::install_github("NBISweden/genecovr")
```

## Usage

### genecovr script quick start

There is a helper script for generating basic plots located in
PACKAGE_DIR/bin/genecovr. Create a data input csv-delimited file with
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

#### Example

There are example files located in PACKAGE_DIR/inst/extdata consisting
of two psl alignment files containing gmap alignments and fasta indices
for the transcript sequences and two for different assembly versions:

- nonpolished.fai - fasta index for raw assembly
- polished.fai - fasta index for polished assembly
- transcripts.fai - fasta index for transcript sequences
- transcripts2nonpolished.psl - gmap alignments, transcripts to raw
  assembly
- transcripts2polished.psl - gmap alignments, transcripts to polished
  assembly

Using these files and the labels `non` and `pol` for the different
assemblies, a `genecovr` input file (called e.g., `assemblies.csv`)
would look as follows:

    nonpol,transcripts2nonpolished.psl,nonpolished.fai,transcripts.fai
    pol,transcripts2polished.psl,polished.fai,transcripts.fai

and the command to run would be:

    genecovr assemblies.csv

#### genecovr options

To list genecovr script options, type ’genecovr -h\`:

    usage: genecovr [-h] [-v] [-p number]
                                 [-d OUTPUT_DIRECTORY] [--height HEIGHT]
                                 [--width WIDTH]
                                 csvfile

    positional arguments:
      csvfile               csv-delimited file with columns
                                1. data label
                                2. mapping file (supported formats: psl)
                                3. assembly file (fasta or fasta index)
                                4. transcript file (fasta or fasta index)

    optional arguments:
      -h, --help            show this help message and exit
      -v, --verbose         print extra output
      -p number, --cpus number
                            number of cpus [default 1]
      -d OUTPUT_DIRECTORY, --output-directory OUTPUT_DIRECTORY
                            output directory
      --height HEIGHT       figure height in inches [default 6.0]
      --width WIDTH         figure width in inches [default 6.0]

### R package vignette

Alternatively, import the library in an R script and use the package
functions. See [Get started](articles/genecovr.html) or run
`vignette("genecovr")` for a minimum working example.
