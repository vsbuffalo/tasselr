# AllGenerics.R -- 
# Copyright (C) 2014 Vince Buffalo <vsbuffalo@gmail.com>
# Distributed under terms of the BSD license.

#' Extract genotype information from a TasselHDF5 object
#'
#' @export
setGeneric("geno", function(x) {
					 standardGeneric("geno")
						 })

#' Load and decode biallelic genotypes from HDF5 file
#'
#' @param x a TasselHDF5 object
#' @param verbose a logical describing whether to be verbose during loading
#' @export
setGeneric("loadBiallelicGenotypes", function(x, verbose=TRUE) {
					 standardGeneric("loadBiallelicGenotypes")
						 })


