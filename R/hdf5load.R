# hdf5load.R -- Functions for loading data in and out of Tassel HDF5 files
# Copyright (C) 2014 Vince Buffalo <vsbuffaloAAAAA@gmail.com>
# Distributed under terms of the BSD license.

schema <- list("4"=list(seqnames="/SeqRegion", indices="/SeqRegionIndices",
                 snpnames="/SnpIds", positions="/Positions",
                 allele_freq_order="/SiteDesc/AlleleFreqOrder", 
                 genotypes="/Genotypes",
                 allele_states="/AlleleStates"),
               "5"=list(seqnames="/Positions/Chromosomes",
                 indices="/Positions/ChromosomeIndices",
                 snpnames="/Positions/SnpIds",
                 positions="/Positions/Positions",
                 genotypes="/Genotypes",
                 taxa_order="/Taxa/TaxaOrder",
                 allele_freq_order="/Genotypes/_Descriptors/AlleleFreqOrder",
                 # Note: Tassel's HDF5 encodes ref/alt alleles as a ncoci x 6
                 # matrix. The first column are all ref alleles, and the other
                 # columns are all alternate alleles (0xF if not
                 # present). This works because there can be at most 6
                 # alleles.
                 taxa="/Taxa",
                 allele_states="/Genotypes/AlleleStates"))

getLoadGenotypeFun <- function(version) {
  schm <- schema[[version]]
  loadGenotypesFuns <- setNames(list( # separate sep names step since number names
  v4=function(x) {
      h5read(x@filename, schm$genotypes)
  },
  v5=function(x) {
      lapply(x@samples, function(sn) {
          # it's best to load in each sample individually in Tassel5's layout
          dt <- paste(schm$genotypes, sn, "calls", sep="/")
          h5read(x@filename, dt)
      })
  }), c('4', '5'))
  loadGenotypesFuns[[version]]
}


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

#' Initialize a TasselHDF5 object, connected to an HDF5 file
#'
#' @param file a path to an HDF5 file
#'
#' @export
initTasselHDF5 <- function(file, version="5") {
  schm <- schema[[version]]
  if (is.null(schm)) stop(sprintf("Did not find schema for version '%s'", version))
  file <- path.expand(file)
  msg <- "could not read /Positions/Chromosomes; did you mean to use version='4'?"
  tryCatch(seqnames <- h5read(file, schm$seqnames), error=function(e) stop(paste(msg, e, sep="\n")))
  tmp <- h5read(file, schm$indices)
  snpnames <- as.character(h5read(file, schm$snpnames))
  seqnames <- factor(seqnames[tmp+1L], levels=seqnames)
  positions <- h5read(file, schm$positions)
  if (version == "5")
      samples <- as.character(h5read(file, schm$taxa_order))
  if (version == "4")
      samples <- names(h5read(file, schm$genotypes, count=1))
  allele_mat <- h5read(file, schm$allele_freq_order)
  ref <- getRefAlleles(allele_mat)
  alt <- as(getAltAlleles(allele_mat), "IntegerList")
  allele_states <- h5read(teo@filename, schm$allele_states)[, 1]
  ranges <- setNames(GRanges(seqnames, IRanges(start=positions, width=1)),
                     snpnames)
  obj <- new("TasselHDF5", filename=file,
             #seqnames=seqnames, positions=positions,
             ranges=ranges,
             ref=ref, alt=alt,
             allele_states=allele_states,
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
            loadFun <- getLoadGenotypeFun(x@version) # get this version's load func
            glist <- loadFun(x)
            vmessage("done.\n")
            vmessage("binding samples together into matrix... ")
            #gmat <- do.call(data.frame, glist) # dataframe avoids R mem issues
            gmat <- do.call(cbind, glist)
            rm(glist); gc() # these are big objects...
            #gmat <- bindGenotypeList(glist)
            #vmessage("done.\n")
            # note: replacing -1L with NA now done in C++
            #vmessage("coercing 0xFF to NA_integer_... ")
            #gmat[gmat == -1L] <- NA_integer_
            #vmessage("done.\n")
            # filter, keeping biallelic only
            vmessage("done.\n")
            vmessage("filtering biallelic loci... ")
            elens <- elementLengths(x@alt) 
            i <- which(elens == 1)
            nremoved <- sum(elens != 1)
            if (nremoved > 0)
                warning(sprintf("Removed %d loci non-biallelic.", nremoved))
            gmat <- gmat[i, ]
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
            #x@genotypes <- encodeNumAltAlleles(gmat, x@ref, alt)
            x@genotypes <- .Call("encodeNumAltAlleles2", nrow(gmat), ncol(gmat),
                                 gmat, x@ref, alt)
            vmessage("done.\n")
            rownames(x@genotypes) <- names(x@ranges)
            colnames(x@genotypes) <- x@samples
            gc() # collecting garbage frees mem (a lot since these objs are large)
            return(x)
          })



