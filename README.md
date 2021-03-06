# tasselr

[![Build Status](https://travis-ci.org/vsbuffalo/tasselr.svg?branch=master)](https://travis-ci.org/vsbuffalo/tasselr)

[Tassel](http://www.maizegenetics.net/index.php?option=com_content&task=view&id=89&Itemid=119)
outputs HDF5 files full of GBS genotyping data. This package is an R interface
to a subset of the information in these files, so users can quickly load and
work with these data.

## Required Packages

From CRAN:

 - methods
 - Rcpp

From [Bioconductor](http://bioconductor.org):

 - GenomicRanges
 - Biobase
 - rhdf5

Install these packages, then using
[devtools](https://github.com/hadley/devtools) (which you can install via
CRAN), do:

    install_github("vsbuffalo/tasselr")

## Loading Tassel HDF5 GBS into R

First, we initialize the HDF5 file with `initTasselHDF5`. By default,
`initTasselHDF5()` assumes the HDF5 file *is in Tassel5's schema*; set the
`version='4'` to change this. `initTasselHDF5()` loads in the loci positions
as a `GRanges` object, and stores reference and alternate alleles (which you
can access with the `ref` and `alt` accessor functions, respectively):

    > gbs <- initTasselHDF5("path/to/mygbs.h5")
    > gbs
		Tassel HDF5 object at 'path/to/mygbs.h5'
		955690 loci x 2060 samples
		Number of chromosomes: 11
		Object size: 127.022 Mb

    > head(alt(gbs), 20)
      3   4   5   6   8  11  13  14  15  17  18
     "T" "A" "C" "G" "T" "C" "T" "T" "G" "G" "T"   [...]

    > head(ref(gbs))
    [1] "C" "C" "C" "C" "C" "A"

## Loading Genotypes

Genotypes are loaded and decoded into the number of alternate alleles they have
(0, 1, 2) for biallelic loci by the method `loadBiallelicGenotypes()`:

    > gbs <- loadBiallelicGenotypes(gbs)

The conversion methods are written in C++ with [Rcpp](http://www.rcpp.org/) so
they're fast-ish.

The accessor function `geno()` can be used to extract this genotype matrix. Note
that the number of loci, and the reference and alternate alleles will change,
as only biallelic loci are kept. The object will always have internal
consistency, and you should always use accessor functions to access data in
slots.

Accessor functions:

 - `granges()`: for loci positions.
 - `alt()`: for alternate alleles.
 - `ref()`: for reference alleles.
 - `geno()`: for genotype matrix.

## Warnings

 - This is tied to Tassel's HDF5 formats (Tassel 4 and 5 only). If these
   change, there are no guarantees the data loaded in will be correct (sorry).
   Always run a sanity checks!

## Todo

 - Lazier loading, e.g. of specific samples or chromosomes at a time through
   reading in specific sample datasets or subsets.

 - Interface to the `/depth` datasets, which contain a lot of information.

 - Possibly re-looking into storing genotyping data using R's raw atomic
   vector, which will use 1/4 of the memory.

 - Apply functions, which load parts from HDF5 files only as needed and apply
   function on them, e.g. functions for taking coverage depth statistics.

## Dev Notes

 - Initially I used R to find and convert missing values 0xFF values into NAs.
This took 261.093 seconds (261.093 user, 129.891 system, 454.579 elapsed). Now
this is done in C++ and it takes about 189.709 seconds (189.709 user, 59.941
system, 290.936 elapsed).

 - Initially I used Rcpp to encode genotypes from tasselr. However, Rcpp's
   `IntegerMatrix` doesn't support matrices that are larger than 2^31 - 1
elements. To remedy this, I wrote everything in C and used a `do.call()`.
