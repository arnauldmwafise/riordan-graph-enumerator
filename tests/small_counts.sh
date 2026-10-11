#!/bin/sh
# Purpose:   regression test: the three independent methods for u(n) (number of unlabeled
#            Riordan graphs) must reproduce the known values for n = 1..8:
#              A  build/rgcount     (enumeration of the L(n) = (4^(n-1)+2)/3 labeled graphs G_n(g,f)
#                                    + nauty canonical forms)
#              B  build/tree_count  (vertex-by-vertex generation tree, different canonical hash)
#              C  build/geng | build/recognize_all  (all graphs of order n, Riordan-labeling search)
# Usage:     sh tests/small_counts.sh          (run by `make check`, from the repository root)
# Arguments: none. Environment: B = build directory (default: build).
# Inputs:    the binaries $B/rgcount, $B/tree_count, $B/geng, $B/recognize_all (must be built).
# Outputs:   stdout only, one line per n: "n=N: u(n)=U (methods A, B, C agree)". No files written.
#            On a mismatch: "n=N: expected U, got A=a B=b C=c" and exit status 1.
# Cost:      a few seconds.
# Notes:     read-only and idempotent. Expected values: docs/COMPUTATIONS.md, section 2.
set -e
B=${B:-build}
expected="1 2 4 11 30 109 406 1864"
n=1
for u in $expected; do
  a=$($B/rgcount $n --threads 1 2>/dev/null | awk '{print $3}')
  b=$($B/tree_count $n --threads 1 2>/dev/null | awk '{print $3}')
  c=$($B/geng -q $n | $B/recognize_all | awk '{print $3}')
  if [ "$a" != "$u" ] || [ "$b" != "$u" ] || [ "$c" != "$u" ]; then
    echo "n=$n: expected $u, got A=$a B=$b C=$c"; exit 1
  fi
  echo "n=$n: u(n)=$u (methods A, B, C agree)"
  n=$((n+1))
done
