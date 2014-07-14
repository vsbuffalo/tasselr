# AllGenerics.R --
# Copyright (C) 2014 Vince Buffalo <vsbuffaloAAAAA@gmail.com>
# Distributed under terms of the BSD license.


#' Load biallelic genotype information
#'
#' @param x a TasselHDF5 object
#' @param verbose logical indicating whether to print status mesages during loading
setGeneric("loadBiallelicGenotypes", function(x, verbose=TRUE) {
           standardGeneric("loadBiallelicGenotypes")
             })

#' Accessor for genotype matrix
#'
#' @param x a TasselHDF5 object
#' @export
setGeneric("geno", function(x) {
           standardGeneric("geno")
             })

#' Accessor for alternate alleles
#' @param x a TasselHDF5 object
#' @param as_char logical indicating whether to return as \code{CharacterList}
#' or as Tassel encodes them, as an \code{IntegerList}
#' @export
#'
setGeneric("alt", function(x, as_char=TRUE) {
           standardGeneric("alt")
             })

#' Accessor for reference alleles
#' @param x a TasselHDF5 object
#' @param as_char logical indicating whether to return as character
#' or as Tassel encodes them, as integer
#' @export
setGeneric("ref", function(x, as_char=TRUE) {
           standardGeneric("ref")
             })

