#!/bin/sh
# Purpose:   compute u(n), the number of unlabeled Riordan graphs (isomorphism classes of the
#            L(n) = (4^(n-1)+2)/3 labeled Riordan graphs G_n(g,f)), for n = NMIN..NMAX with the two
#            independent exhaustive methods (A: build/rgcount, B: build/tree_count), compare them,
#            and write data/counts.csv.
# Usage:     scripts/compute_counts.sh [NMIN] [NMAX] [THREADS]
# Arguments: NMIN     smallest order (default 1)
#            NMAX     largest order (default 14)
#            THREADS  worker threads per run (default: all online cores)
#            Passes (memory split) are chosen automatically: 1 for n <= 15, 4 for n = 16,
#            12 for n = 17, 48 for n >= 18.
# Inputs:    build/rgcount, build/tree_count (built with `make` if missing).
# Outputs:   data/runs/A<n>.txt, B<n>.txt  final line "n L(n) u(n)" of each method;
#            data/runs/A<n>.ckpt, B<n>.ckpt per-pass checkpoints; data/runs/A<n>.err, B<n>.err
#            stderr (appended);
#            data/counts.csv  rebuilt from all available runs (n = 1..30) with header
#            "n,labeled,unlabeled_method_A,unlabeled_method_B,agree" (agree = yes|no), also
#            printed to stdout.
# Cost:      approximately (one core) n=14 1.5 min, n=15 8 min, n=16 45 min, n=17 4 h per
#            method; memory about 16 bytes per class and pass. See docs/COMPUTATIONS.md, section 2.
# Notes:     idempotent and resumable: an order whose result file is non-empty is skipped, and
#            interrupted runs resume from their checkpoints. Exits on the first failing command
#            (set -e). Works from any directory (changes to the repository root).
set -e
cd "$(dirname "$0")/.."
NMIN=${1:-1}; NMAX=${2:-14}; T=${3:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)}
R=data/runs; mkdir -p $R
[ -x build/rgcount ] && [ -x build/tree_count ] || make
for n in $(seq $NMIN $NMAX); do
  P=1; [ $n -ge 16 ] && P=4; [ $n -ge 17 ] && P=12; [ $n -ge 18 ] && P=48
  [ -s $R/A$n.txt ] || ./build/rgcount $n --threads $T --passes $P --checkpoint $R/A$n.ckpt > $R/A$n.txt 2>> $R/A$n.err
  [ -s $R/B$n.txt ] || ./build/tree_count $n --threads $T --passes $P --checkpoint $R/B$n.ckpt > $R/B$n.txt 2>> $R/B$n.err
done
echo "n,labeled,unlabeled_method_A,unlabeled_method_B,agree" > data/counts.csv
for n in $(seq 1 30); do
  [ -s $R/A$n.txt ] && [ -s $R/B$n.txt ] || continue
  set -- $(cat $R/A$n.txt); L=$2; a=$3
  set -- $(cat $R/B$n.txt); b=$3
  ok=no; [ "$a" = "$b" ] && ok=yes
  echo "$n,$L,$a,$b,$ok" >> data/counts.csv
done
cat data/counts.csv
