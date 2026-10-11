#!/bin/sh
# Purpose:   method C for u(n): generate every graph of order n (nauty geng, one per isomorphism
#            class) and search each for a Riordan labeling, i.e. a vertex order in which it equals
#            some labeled Riordan graph G_n(g,f) (build/recognize_all). Independent of the
#            enumeration of the L(n) = (4^(n-1)+2)/3 labeled graphs and of canonical forms.
# Usage:     scripts/method_c.sh [NMAX]
# Arguments: NMAX  orders n = 1..NMAX are processed (default 9)
# Inputs:    build/geng, build/recognize_all (built with `make` if missing).
# Outputs:   data/lists/R_<n>.g6, NR_<n>.g6  Riordan / non-Riordan graphs of order n (graph6);
#            data/runs/C<n>.txt  summary line "n graphs riordan nonriordan" (also printed);
#            data/runs/C<n>.err  search-node statistics (--stats).
#            The lists are the input of build/structure (analysis/structure.cpp).
# Cost:      n <= 9 in about 2 minutes; n = 10 (1.2·10^7 graphs) takes a few hours on one core and
#            NR_10.g6 is large (stored compressed as data/lists/NR_10.g6.gz).
# Notes:     idempotent and resumable per order: n is skipped if data/runs/C<n>.txt is non-empty;
#            an interrupted order is redone from scratch. Exits on the first failing command
#            (set -e). Works from any directory (changes to the repository root).
set -e
cd "$(dirname "$0")/.."
NMAX=${1:-9}
mkdir -p data/lists data/runs
[ -x build/geng ] && [ -x build/recognize_all ] || make
for n in $(seq 1 $NMAX); do
  [ -s data/runs/C$n.txt ] && continue
  ./build/geng -q $n | ./build/recognize_all --stats --riordan data/lists/R_$n.g6 --nonriordan data/lists/NR_$n.g6 > data/runs/C$n.txt 2> data/runs/C$n.err
  cat data/runs/C$n.txt
done
