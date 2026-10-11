#!/bin/sh
# tests/rgsample_test.sh -- statistical and functional self-tests of the sampler (about 15 s).
#
# Input : environment variable B = build directory (default: build).
# Output: one line per check; exits non-zero on the first failure.
#
# Checks:
#  1. Unlabeled model: the acceptance rate estimates u(8)/L(8) = 1864/5462 = 0.3413; with 6000
#     accepted samples its standard error is about 0.004, and the test allows +-0.02.
#  2. Estimator: mean 1/|C(H)| over uniform labeled graphs estimates u(9)/L(9) = 0.3739 (+-0.03).
#  3. Structure filters: every graph accepted with --avoids K3 has 0 triangles, and every graph
#     accepted with --contains C4 --connected is connected.
#  4. Reproducibility: the same seed gives the same output.
set -e
B=${B:-build}
rate=$($B/rgsample 8 --samples 6000 --model unlabeled --summary --seed 1 | sed -n 's/.*(rate \([0-9.]*\)).*/\1/p')
awk -v r="$rate" 'BEGIN{ if (r < 0.3213 || r > 0.3613) { print "unlabeled acceptance rate " r " not near 0.3413"; exit 1 } }'
echo "unlabeled sampler: acceptance rate $rate (expected 0.3413)"
est=$($B/rgsample 9 --samples 3000 --class-size --summary --seed 2 | sed -n 's/.*mean 1\/|C(H)| = \([0-9.]*\),.*/\1/p')
awk -v e="$est" 'BEGIN{ if (e < 0.3439 || e > 0.4039) { print "estimate " e " not near 0.3739"; exit 1 } }'
echo "estimator: u(9)/L(9) estimated as $est (exact 0.3739)"
bad=$($B/rgsample 12 --samples 300 --avoids K3 --format csv --seed 3 2>/dev/null | awk -F, 'NR>1 && $8 != 0' | wc -l)
[ "$bad" -eq 0 ] || { echo "--avoids K3 accepted graphs with triangles"; exit 1; }
bad=$($B/rgsample 12 --samples 300 --contains C4 --connected --format csv --seed 4 2>/dev/null | awk -F, 'NR>1 && $9 != 1' | wc -l)
[ "$bad" -eq 0 ] || { echo "--connected accepted a disconnected graph"; exit 1; }
echo "structure filters: ok"
a=$($B/rgsample 20 --samples 20 --seed 9 --format csv 2>/dev/null | md5sum)
b=$($B/rgsample 20 --samples 20 --seed 9 --format csv 2>/dev/null | md5sum)
[ "$a" = "$b" ] || { echo "same seed gave different output"; exit 1; }
echo "rgsample: all checks passed"
