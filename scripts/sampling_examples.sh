#!/bin/sh
# scripts/sampling_examples.sh -- regenerate the sampler showcase in data/sampling/.
#
# Purpose : illustrate rgsample at orders beyond the exhaustive range (n > 17) and validate it
#           where exact values are known (n <= 17).
# Usage   : scripts/sampling_examples.sh
# Inputs  : build/rgsample (built by make if missing)
# Outputs : data/sampling/*.txt, one summary per experiment (overwritten on every run)
# Cost    : about 2 minutes on one core.
# Notes   : seeds are fixed, so a run is reproducible on the same platform. The C++ standard
#           library may implement random distributions differently, so other platforms can give
#           different (but statistically equivalent) samples.
set -e
cd "$(dirname "$0")/.."
[ -x build/rgsample ] || make build/rgsample
D=data/sampling; mkdir -p $D
# 1. validation: estimate u(12) and u(14) (exact: 671656 and 11481042)
build/rgsample 12 --samples 2000 --class-size --summary --seed 1 > $D/estimate_u12.txt
build/rgsample 14 --samples 1000 --class-size --summary --seed 1 > $D/estimate_u14.txt
# 2. typical labeled Riordan graphs of larger order
for n in 20 32 48 64; do
  build/rgsample $n --samples 5000 --summary --seed 1 --report K4 --report C4 --report claw --report cochair > $D/labeled_n$n.txt
done
# 3. exactly uniform unlabeled Riordan graphs (rejection with exact class sizes)
build/rgsample 12 --samples 1000 --model unlabeled --summary --seed 1 --report K3 --report C4 > $D/unlabeled_n12.txt
# 4. structure-conditioned samples: K4-free connected graphs of order 30, and graphs of order 30
#    containing an induced 5-cycle but no K4
build/rgsample 30 --samples 300 --avoids K4 --connected --summary --seed 1 > $D/K4free_connected_n30.txt
build/rgsample 30 --samples 300 --contains-induced C5 --avoids K4 --summary --seed 1 > $D/inducedC5_K4free_n30.txt
# 5. graphs with f'(0) = 1 of order 22: class sizes (reversal pairs)
build/rgsample 22 --samples 100 --fprime 1 --class-size --summary --seed 1 > $D/fprime1_classes_n22.txt
ls -l $D
