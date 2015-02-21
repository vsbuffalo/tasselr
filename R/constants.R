# constants.R -- constants to handle conversions from Tassel
# Copyright (C) 2014 Vince Buffalo <vsbuffaloAAAAA@gmail.com>
# Distributed under terms of the BSD license.

# Tassel alleles take up four bits:
# A = 0x0
# C = 0x1
# G = 0x2
# T = 0x3
# Ins(+) = 0x4
# Del(-) = 0x5
# N = 0xF
# Each of these has a position in a vector with length 16. Since R is
# 1-indexed, the appoprtiate allele is the bit value in decimal - 1.
TASSELL_ALLELES <- c("A", "C", "G", "T", "+", "-", NA, NA,
                     NA, NA, NA, NA, NA, NA, NA, "N")


