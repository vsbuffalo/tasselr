# hdf5load.R -- Functions for loading data in and out of Tassel HDF5 files
# Copyright (C) 2014 Vince Buffalo <vsbuffaloAAAAA@gmail.com>
# Distributed under terms of the BSD license.

schema <- list("4"=list(seqnames="/SeqRegion", indices="/SeqRegionIndices",
                        snpnames="/SnpIds", positions="/Positions",
                        genotypes="/Genotypes"),
               "5"=list(seqnames="/Positions/Chromosomes",
                        indices="/Positions/ChromosomeIndices",
                        snpnames="/Positions/SnpIds",
                        positions="/Positions/Positions",
                        genotypes="/Genotypes",
                        taxa_order="/Taxa/TaxaOrder",
                        taxa="/Taxa"))


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
	h5read(file,"/Genotypes/_Descriptors/AlleleFreqOrder")
}

#' Initialize a TasselHDF5 object, connected to an HDF5 file
#'
#' @param file a path to an HDF5 file
#'
#' @export
initTasselHDF5 <- function(file, version="5") {
  schm <- schema[[version]]
  file <- path.expand(file)
  seqlevels <- h5read(file, schm$seqnames)
  tmp <- h5read(file, schm$indices)
  snpnames <- as.character(h5read(file, schm$snpnames))
  seqnames <- factor(seqlevels[tmp+1L], levels=seqlevels)
  positions <- h5read(file, schm$positions)
  if (version == "5")
      samples <- as.character(h5read(file, schm$taxa_order))
  if (version == "4")
      samples <- names(h5read(file, schm$genotypes, count=1))
  allele_mat <- extractAlleleFreqOrder(file)
  ref <- getRefAlleles(allele_mat)
  alt <- as(getAltAlleles(allele_mat), "IntegerList")
  ranges <- setNames(GRanges(seqnames, IRanges(start=positions, width=1)),
                     snpnames)
  obj <- new("TasselHDF5", filename=file,
             #seqnames=seqnames, positions=positions,
             ranges=ranges,
             ref=ref, alt=alt,
             samples=samples, version=version)
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
            schm <- schema[[x@version]]
            vmessage(sprintf("loading in genotypes from HDF5 file '%s'... ",
                             basename(x@filename)))
            glist <- lapply(h5read(x@filename, schm$genotypes),
                            function(x) as.integer(x[[1]]))
            glist <- glist[x@samples]
            vmessage("done.\n")
            vmessage("binding samples together into matrix... ")
            #gmat <- do.call(data.frame, glist) # dataframe avoids R mem issues
            gmat <- bindGenotypeList(glist)
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



