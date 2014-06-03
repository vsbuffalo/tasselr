import sys
import gzip
from collections import Counter

with gzip.open(sys.argv[1]) as hapmap_file:
    rs = sys.argv[2]
    for line in hapmap_file:
        parts = line.strip().split("\t")
        if parts[0] == rs:
            cc = Counter(parts[11:])
            print parts[1], cc, sum([x for k, x in cc.items()])
            break
