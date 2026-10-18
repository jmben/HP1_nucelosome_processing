#!/usr/bin/env python3
"""
invert_bigwig.py  —  Multiply all values in a bigwig by -1
=============================================================
Usage:
    python invert_bigwig.py input.bw output.bw
=============================================================
"""
import sys
import pyBigWig
import numpy as np

def invert_bigwig(in_path: str, out_path: str):
    bw_in  = pyBigWig.open(in_path)
    bw_out = pyBigWig.open(out_path, "w")
    chroms = bw_in.chroms()
    bw_out.addHeader(list(chroms.items()))
    for chrom, chrom_len in chroms.items():
        print(f"  Processing {chrom} ({chrom_len:,} bp) ...")
        intervals = bw_in.intervals(chrom)
        if intervals is None or len(intervals) == 0:
            continue
        starts = [iv[0] for iv in intervals]
        ends   = [iv[1] for iv in intervals]
        values = [-iv[2] for iv in intervals]   # invert sign
        bw_out.addEntries([chrom] * len(starts), starts, ends=ends, values=values)
    bw_in.close()
    bw_out.close()
    print(f"✓  Saved inverted bigwig → {out_path}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python invert_bigwig.py input.bw output.bw")
        sys.exit(1)
    invert_bigwig(sys.argv[1], sys.argv[2])
