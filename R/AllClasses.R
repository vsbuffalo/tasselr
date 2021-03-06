# AllClasses.R
# Copyright (C) 2014 Vince Buffalo <vsbuffaloAAAAA@gmail.com>
# Distributed under terms of the BSD license.

#' An S4 class that stores links to Tassel HDF5 objects, loci positions, and genotypes
#'
#' @slot filename path to HDF5 file
#' @slot ranges a \code{GenomicRanges} object of loci
#' @slot ref a Tassel-encoded integer vector of reference alleles
#' @slot alt a Tassel-encoded \code{IntegerList} of alternate alleles
#' @slot genotypes a matrix of bialleic genotypes
#' @slot allele_states character vector of allele states
#' @slot samples sample names
#' @slot version string Tassel version
#'
#' @exportClass TasselHDF5
setClass("TasselHDF5",
         slots=list(filename="character",
                    #seqnames="Rle",
                    #positions="Rle",
                    ranges="GRanges",
                    #ref="character",
                    #alt="CharacterList",
                    ref="integer",
                    alt="IntegerList",
                    allele_states="character",
                    genotypes="matrix",
                    samples="character",
                    version="character"))
