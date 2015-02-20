
#' A function which binds together a list of genotypes from Tassel HDF5
#'
bindGenotypeList <- function(x) {
 if (length(unique(sapply(x, length))) > 1)
   stop("All elements of x must have the same length.")
 if (unique(sapply(x, typeof)) != "integer")
   stop("All elements of x must be integer vectors.")
 .bindCols <- cfunction(c(x="LISTSXP"), "
  SEXP out;
  int nind = length(x);
  int nloci = length(VECTOR_ELT(x, 0));
  out = PROTECT(allocMatrix(INTSXP, nloci, nind));

  for (int i=0; i < nind; i++) {
    for (int l=0; l < nloci; l++) {
      INTEGER(out)[l + nloci*i] = INTEGER(VECTOR_ELT(x, i))[l];
    }
  }
  UNPROTECT(1);
  return out;
")
 .bindCols(x)
 
}
