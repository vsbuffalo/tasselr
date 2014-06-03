# AllGenerics.R --
# Copyright (C) 2014 Vince Buffalo <vsbuffaloAAAAA@gmail.com>
# Distributed under terms of the BSD license.


#' Load biallelic genotype information
#'
setGeneric("loadBiallelicGenotypes", function(x, verbose=TRUE) {
           standardGeneric("loadBiallelicGenotypes")
             })

#' Accessor for genotype matrix
#'
setGeneric("geno", function(x) {
           standardGeneric("geno")
             })

#' Accessor for alternate alleles
#'
setGeneric("alt", function(x, as_char=TRUE) {
           standardGeneric("alt")
             })

#' Accessor for reference alleles
#'
setGeneric("ref", function(x, as_char=TRUE) {
           standardGeneric("ref")
             })

