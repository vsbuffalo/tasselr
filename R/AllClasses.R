# AllClasses.R
# Copyright (C) 2014 Vince Buffalo <vsbuffalo@gmail.com>
# Distributed under terms of the BSD license.

#' An S4 class that stores links to Tassel HDF5 objects, loci positions, and genotypes
#'
#' @slot filename path to HDF5 file
#' @slot ranges a \code{GenomicRanges} object of loci
#' @slot ref a Tassel-encoded integer vector of reference alleles
#' @slot alt a Tassel-encoded \code{IntegerList} of alternate alleles
#' @slot genotypes a matrix of bialleic genotypes
#' @slot samples sample names
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
										genotypes="matrix",
										samples="character"))
