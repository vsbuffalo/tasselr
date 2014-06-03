# hdf5load.R -- Functions for loading data in and out of Tassel HDF5 files
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

extractAlleleFreqOrder <- function(file) {
	# Tassel's HDF5 encodes ref/alt alleles as a ncoci x 6 matrix. The first
	# column are all ref alleles, and the other columns are all alternate alleles
	# (0xF if not present). This works because there can be at most 6 alleles.
	h5read(file,"/SiteDesc/AlleleFreqOrder")
}

#' Initialize a TasselHDF5 object, connected to an HDF5 file
#'
#' @param file a path to an HDF5 file
#'
#' @export
initTasselHDF5 <- function(file) {
  file <- path.expand(file)
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


#' Load and decode biallelic genotypes from HDF5 file
#'
#' @param x a TasselHDF5 object
#' @param verbose a logical describing whether to be verbose during loading
#'
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
						# note: replacing -1L with NA now done in C++
            #vmessage("coercing 0xFF to NA_integer_... ")
            #gmat[gmat == -1L] <- NA_integer_
            #vmessage("done.\n")
            # filter, keeping biallelic only
            vmessage("filtering biallelic loci... ")
            i <- which(elementLengths(x@alt) == 1)
            x@ranges <- x@ranges[i]
            x@ref <- x@ref[i]
            x@alt <- x@alt[i]
            stopifnot(length(x@ref) == length(unlist(x@alt)))
            stopifnot(length(x@ranges) == length(x@ref))
            stopifnot(length(i) == length(x@alt))
            vmessage("done.\n")
            vmessage("encoding genotypes... ")
						alt <- as(x@alt, "integer")
						stopifnot(is(alt, "integer"))
            x@genotypes <- encodeNumAltAlleles(gmat[i, ], x@ref, alt)
            vmessage("done.\n")
            rownames(x@genotypes) <- names(x@ranges)
            colnames(x@genotypes) <- x@samples
            return(x)
          })

