#include <stdint.h>
#include <R.h>
#include <Rinternals.h>

SEXP encodeNumAltAlleles2(SEXP nrow, SEXP ncol, SEXP x, SEXP ref, SEXP alt) {
  unsigned int nrows = INTEGER(nrow)[0], ncols = INTEGER(ncol)[0];
  SEXP nalt = PROTECT(allocMatrix(INTSXP, nrows, ncols));
  int8_t a1, a2, geno;
  for (unsigned int i = 0; i < nrows; i++) {
    for (unsigned int j = 0; j < ncols; j++) {
      //printf("loop start: i = %d; j = %d; nrows = %d; i + nrows*j = %d, nalt: %d x %d\n", i, j, nrows, i + nrows*j, nrows, ncols);
      if (INTEGER(x)[i + nrows*j] == -1) {
        // 0xFF (-1 once converted to R's 32-bit integers) is missing genotype
	INTEGER(nalt)[i + nrows*j] = NA_INTEGER;
      } else {
	geno = INTEGER(x)[i + nrows*j];
	a1 = geno >> 4;
	a2 = geno & 0x0F;
	if ((a1 != INTEGER(ref)[i] && a1 != INTEGER(alt)[i]) ||
	    (a2 != INTEGER(ref)[i] && a2 != INTEGER(alt)[i])) {
	  error("genotype includes allele that's not ref or alt");
	}
	INTEGER(nalt)[i + nrows*j] = (a1 == INTEGER(alt)[i]) +
	  (a2 == INTEGER(alt)[i]);
      }
    }
  }
  UNPROTECT(1);
  return nalt;
}









