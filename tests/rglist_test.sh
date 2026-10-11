#!/bin/sh
# Purpose:   regression test of build/rglist (listing of labeled/unlabeled Riordan graphs G_n(g,f)):
#            it must reproduce, for n = 1..8, the known counts of
#              --unlabeled                        u(n)
#              --labeled                          L(n) = (4^(n-1)+2)/3
#              --labeled --filter complement      labeled graphs whose complement is labeled Riordan
#              --filter complement                unlabeled graphs whose complement is Riordan
#              --filter selfcomplementary         self-complementary Riordan graphs
#            check that the class sizes in the CSV listing (column 8) sum to L(n) for n = 4, 6, 8,
#            and that --show prints the expected edge list of G_4(1+t^2, t).
# Usage:     sh tests/rglist_test.sh           (run by `make check`, from the repository root)
# Arguments: none. Environment: B = build directory (default: build).
# Inputs:    the binary $B/rglist.
# Outputs:   stdout only: "rglist: all checks passed" on success; otherwise a message naming the
#            failed check and exit status 1. No files written.
# Cost:      a few seconds.
# Notes:     read-only and idempotent. Expected values: docs/COMPUTATIONS.md, sections 2 and 4.2.
set -e
B=${B:-build}
check() { got=$($B/rglist $1 $2 --count); [ "$got" = "$3" ] || { echo "rglist $1 $2: expected $3, got $got"; exit 1; }; }
n=1; for u in 1 2 4 11 30 109 406 1864; do check $n "--unlabeled" $u; n=$((n+1)); done
n=1; for L in 1 2 6 22 86 342 1366 5462; do check $n "--labeled" $L; n=$((n+1)); done
n=1; for c in 1 2 4 12 28 92 220 732; do check $n "--labeled --filter complement" $c; n=$((n+1)); done
n=1; for c in 1 2 4 11 26 90 230 707; do check $n "--filter complement" $c; n=$((n+1)); done
n=1; for s in 1 0 0 1 2 0 0 1; do check $n "--filter selfcomplementary" $s; n=$((n+1)); done
# the class sizes of the unlabeled listing add up to L(n)
for n in 4 6 8; do
  s=$($B/rglist $n --format csv | tail -n +2 | awk -F, '{s+=$8} END{print s}')
  L=$($B/rglist $n --labeled --count)
  [ "$s" = "$L" ] || { echo "class sizes at n=$n sum to $s, expected $L"; exit 1; }
done
# --show reproduces a known graph
$B/rglist 4 --show "1+t^2" t | grep -q "edges (4): 12 14 23 34" || { echo "--show failed"; exit 1; }
echo "rglist: all checks passed"
