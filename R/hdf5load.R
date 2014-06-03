# hdf5load.R --
# Copyright (C) 2014 Vince Buffalo <vsbuffaloAAAAA@gmail.com>
# Distributed under terms of the BSD license.

getRefAlleles <- function(x) {
  # x is a AlleleFreqOrder matrix
  return(ifelse(x[, 1] == -1L, NA, x[, 1]))
}

getAltAlleles <- function(x) {
  # is is an AlleleFreqOrder matrix
  lapply(split(x[, -1], seq_len(nrow(x))), function(alts) {
         alts[alts != 0xFL]
         })
}

#' @export
initTasselHDF5 <- function(file) {
  seqlevels <- h5read(file, "/SeqRegion")
  tmp <- h5read(file, "/SeqRegionIndices")
  snpnames <- as.character(h5read(file, "/SnpIds"))
  seqnames <- factor(seqlevels[tmp+1L], levels=seqlevels)
  positions <- h5read(file, "/Positions")
  samples <- names(h5read(file, "/Genotypes", count=1))
  allele_mat <- extractAlleleFreqOrder(file)
  ref <- getRefAlleles(allele_mat)
  alt <- as(getAltAlleles(allele_mat), "IntegerList")
  ranges <- setNames(GRanges(seqnames, IRanges(start=positions, width=1)),
                     snpnames)
  obj <- new("TasselHDF5", filename=file,
             #seqnames=seqnames, positions=positions,
             ranges=ranges,
             ref=ref, alt=alt,
             samples=samples)
  return(obj)
}

#' @export
setMethod("loadBiallelicGenotypes",
          c(x="TasselHDF5"),
          function(x, verbose=TRUE) {
            vmessage <- function(x) {
              if (verbose)
                message(x, appendLF=FALSE)
            }
            vmessage(sprintf("loading in genotypes from HDF5 file '%s'... ",
                             basename(x@filename)))
            glist <- h5read(x@filename, "/Genotypes")
            vmessage("done.\n")
            vmessage("coercing to matrix... ")
            gmat <- do.call(cbind, glist)
            vmessage("done.\n")
            vmessage("coercing 0xFF to NA_integer_... ")
            gmat[gmat == -1L] <- NA_integer_
            vmessage("done.\n")
            # filter, keeping biallelic only
            i <- which(elementLengths(x@alt) == 1)
            x@ranges <- x@ranges[i]
            x@ref <- x@ref[i]
            x@alt <- x@alt[i]
            x@genotypes <- encodeNumAltAlleles(gmat[i, ], x@ref, unlist(x@alt))
            rownames(x@genotypes) <- names(x@ranges)
            colnames(x@genotypes) <- x@samples
            return(x)
          })

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

#' @export
setMethod("geno",
          c(x="TasselHDF5"),
          function(x) {
            return(x@genotypes)
          })

#' @export
setMethod("granges",
          c(x="TasselHDF5"),
          function(x, use.mcols=FALSE, ...) {
            return(x@ranges)
          })

