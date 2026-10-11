#!/bin/sh
# Purpose:   data for Conjecture 2 of Cheon et al. (2019) at k = 6: run build/bell on all 2^31
#            io-decomposable Riordan graphs of the Bell type G_64(g, tg) ([t^j]g = [t^(2j+1)]g),
#            recording the diameter distribution, Wiener-index range and every graph of diameter 6.
#            Runs `bell all 64 --threads 1 --no-alpha --part i 32` for i = 0..31.
# Usage:     scripts/bell_all64.sh
# Arguments: none.
# Inputs:    build/bell (must already be built: `make build/bell`).
# Outputs:   data/bell/all64/part<i>.txt, i = 0..31, each the `bell all` output for one slice of
#            2^26 graphs ("n=64 graphs=... diam: ... wiener:[...]" and the diameter-6 list);
#            data/bell/all64/tmp<i> while slice i is running. The summary over all
#            slices is data/bell/all_64.txt (see docs/COMPUTATIONS.md, section 4.4).
# Cost:      hours in total (32 single-threaded slices, run one after another).
# Notes:     idempotent and resumable: a slice with a non-empty part<i>.txt is skipped; an
#            interrupted slice leaves only tmp<i> and is redone. Do not run two copies at once
#            (both would work on the same slice and tmp file); to use several cores, start
#            slices by hand with `build/bell all 64 --threads 1 --no-alpha --part i 32`.
#            Works from any directory (changes to the repository root).
cd "$(dirname "$0")/.."; D=data/bell/all64
mkdir -p $D
for i in $(seq 0 31); do
  [ -s $D/part$i.txt ] || { ./build/bell all 64 --threads 1 --no-alpha --part $i 32 > $D/tmp$i && mv $D/tmp$i $D/part$i.txt; }
done
