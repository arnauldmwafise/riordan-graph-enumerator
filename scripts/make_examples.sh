#!/bin/sh
# Purpose:   regenerate the browsable listings in examples/ with build/rglist: Markdown tables
#            with SVG drawings of labeled Riordan graphs G_n(g,f) and of unlabeled ones, and CSV.
# Usage:     scripts/make_examples.sh
# Arguments: none.
# Inputs:    build/rglist (built with `make build/rglist` if missing).
# Outputs:   (the directories below are deleted and recreated)
#            examples/unlabeled/order_<n>.md + order_<n>/*.svg            n = 1..6
#            examples/labeled/order_<n>.md + SVGs                          n = 2..4
#            examples/complement/order_<n>.md, labeled_order_<n>.md + SVGs  n = 4..6
#                (graphs whose complement is Riordan, unlabeled and labeled)
#            examples/self_complementary/order_<n>.md + SVGs               n = 4, 5, 8, 9, 12
#            examples/csv/unlabeled_order_<n>.csv (n = 1..9), labeled_order_<n>.csv (n = 1..7)
#            and the line "examples regenerated" on stdout.
# Cost:      a few seconds.
# Notes:     idempotent: the output is fully regenerated on each run (deterministic). Exits on the
#            first failing command (set -e). Works from any directory.
set -e
cd "$(dirname "$0")/.."
[ -x build/rglist ] || make build/rglist
X=examples
rm -rf $X/unlabeled $X/labeled $X/complement $X/self_complementary $X/csv
mkdir -p $X/unlabeled $X/labeled $X/complement $X/self_complementary $X/csv
for n in 1 2 3 4 5 6; do
  ./build/rglist $n --unlabeled --format markdown --svg $X/unlabeled/order_$n > $X/unlabeled/order_$n.md
done
for n in 2 3 4; do
  ./build/rglist $n --labeled --format markdown --svg $X/labeled/order_$n > $X/labeled/order_$n.md
done
for n in 4 5 6; do
  ./build/rglist $n --unlabeled --filter complement --format markdown --svg $X/complement/order_$n > $X/complement/order_$n.md
  ./build/rglist $n --labeled --filter complement --format markdown --svg $X/complement/labeled_order_$n > $X/complement/labeled_order_$n.md
done
for n in 4 5 8 9 12; do
  ./build/rglist $n --unlabeled --filter selfcomplementary --format markdown --svg $X/self_complementary/order_$n > $X/self_complementary/order_$n.md
done
for n in 1 2 3 4 5 6 7 8 9; do
  ./build/rglist $n --unlabeled --format csv > $X/csv/unlabeled_order_$n.csv
done
for n in 1 2 3 4 5 6 7; do
  ./build/rglist $n --labeled --format csv > $X/csv/labeled_order_$n.csv
done
echo "examples regenerated"
