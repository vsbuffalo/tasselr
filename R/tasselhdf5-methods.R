# tasselhdf5-methods.R --
# Copyright (C) 2014 Vince Buffalo <vsbuffaloAAAAAA@gmail.com>
# Distributed under terms of the BSD license.

#' Pretty-print a TasselHDF5 object
#'
#' @param object a TasselHDF5 object
#'
#' @export
setMethod("show",
          c(object="TasselHDF5"),
          function(object) {
            cat(sprintf("Tassel HDF5 object at '%s'\n%d loci x %d samples\n",
                        object@filename,
                        length(object@ranges), length(object@samples)))
            cat(sprintf("Number of chromosomes: %d\nObject size: %s Mb\n",
                        length(seqlevels(object@ranges)),
                        round(object.size(object)/1024^2, 3)))
          })

#' Accessor for genotype information from a TasselHDF5 object
#'
#' @param x a TasselHDF5 object
#' @export
setMethod("geno",
          c(x="TasselHDF5"),
          function(x) {
            return(x@genotypes)
          })

#' Accessor for genomic ranges from a TasselHDF5 object
#'
#' @param x a TasselHDF5 object
#' @export
setMethod("granges",
          c(x="TasselHDF5"),
          function(x, use.mcols=FALSE, ...) {
            return(x@ranges)
          })

#' Accessor for alterate alleles from a TasselHDF5 object
#'
#' @param x a TasselHDF5 object
#' @param as_char logical indicating whether to return as \code{CharacterList}
#' or as Tassel encodes them, as an \code{IntegerList}
#' @export
setMethod("alt",
          c(x="TasselHDF5"),
          function(x, as_char=TRUE) {
            if (!as_char)
              return(x@alt)
            out <- sapply(x@alt, function(x) TASSELL_ALLELES[x+1L])
            out
          })

#' Accessor for reference alleles from a TasselHDF5 object
#'
#' @param x a TasselHDF5 object
#' @param as_char logical indicating whether to return as character
#' or as Tassel encodes them, as integer
#' @export
setMethod("ref",
          c(x="TasselHDF5"),
          function(x, as_char=TRUE) {
            if (!as_char)
              return(x@ref)
            out <- sapply(x@ref, function(x) TASSELL_ALLELES[x+1L])
            return(out)
          })


#' Accessor for samples from a TasselHDF5 object
#'
#' @param object a TasselHDF5 object
#'
#' @importMethodsFrom Biobase samples
#' @export
setMethod("samples",
          c(object="TasselHDF5"),
          function(object) {
              return(object@samples)
          })













