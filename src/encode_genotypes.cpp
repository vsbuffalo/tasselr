#include <Rcpp.h>
#include <cstdint>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix encodeNumAltAlleles(IntegerMatrix x, IntegerVector ref, IntegerVector alt) {
	int nrow = x.nrow();
	int ncol = x.ncol();
	IntegerMatrix nalt(nrow, ncol);
	int8_t a1, a2, geno;

	for (int i = 0; i < nrow; i++) {
		for (int j = 0; j < ncol; j++) {
			if (IntegerVector::is_na(x(i, j))) {
				nalt(i, j) = NA_INTEGER;
			} else {
				geno = x(i, j);
				a1 = geno >> 4;
				a2 = geno & 0xF;
				if ((a1 != ref[i] && a1 != alt[i]) || (a2 != ref[i] && a2 != alt[i])) {
					// TODO error
					std::cout << "geno: " << static_cast<int16_t>(geno) << std::endl <<
						"(a1,a2): (" << static_cast<int16_t>(a1) << "," << static_cast<int16_t>(a2) << ")" << std::endl <<
						"ref: " << static_cast<int16_t>(ref[i]) << " alt: " << static_cast<int16_t>(alt[i]) << std::endl <<
						"(i,j): (" << i << "," << j << ")" << std::endl;
					throw Rcpp::exception("genotype includes allele that's not ref or alt");
				}
				nalt(i, j) = (a1 == alt[i]) + (a2 == alt[i]);
			}
		}
	}
	return nalt;
}

