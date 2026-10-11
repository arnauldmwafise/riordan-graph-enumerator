# Architecture and API reference

This document describes how the code is organised. It covers the header-only library
`include/rgraph/`, the executables built on it, and the flow of data from the programs into `data/`
and `examples/`. Every source file has a header comment covering purpose, method, usage, inputs and
outputs, and every function is documented in place with Doxygen-style comments. With
[Doxygen](https://www.doxygen.nl) installed, `doxygen Doxyfile` builds an HTML reference in
`docs/html/`.

## 1. Conventions used throughout

| item | convention |
|---|---|
| polynomial over Z₂ | a 64-bit word `u64`; bit k is the coefficient of tᵏ |
| graph | `rgraph::Graph`: `n ≤ 64` and `adj[v]`, the 64-bit neighbourhood mask of vertex v |
| vertex numbers | 0-based in code; vertex v in code is vertex v+1 in the paper and in all printed output |
| labeled Riordan graph | `G_n(g,f)`: for i > j (1-based), ij is an edge iff `[t^(i-2)] g f^(j-1) = 1` |
| canonical parameters | `g mod t^(n-1)` and, if `g ≠ 0` has lowest term tʳ, `f mod t^(n-1-r)`. Every labeled Riordan graph has exactly one such pair |
| Riordan matrix | `B[i][j] = [tⁱ] g fʲ`, 0 ≤ j ≤ i ≤ n−2. The adjacency lower triangle is `A[i][j] = B[i−2][j−1]` |
| class C(H) | the set of labeled Riordan graphs isomorphic to H; `u(n) = Σ_H 1/\|C(H)\|` |

## 2. Library `include/rgraph/` (header-only)

| header | responsibility | main API |
|---|---|---|
| `riordan_graph.hpp` | definition and labeled enumeration | `Graph`; `riordan_graph(n,g,f)` builds G_n(g,f) column by column; `mul_trunc`, `low_mask`; `labeled_count(n) = (4^(n-1)+2)/3`; `Enumeration` numbers all labeled Riordan graphs `0..L(n)-1` (`pair(idx,g,f)`, `graph(idx)`) |
| `canon.hpp` | isomorphism classes (nauty) | `canonical_graph(G)` (canonical form with degree partition); `class_hash(G)`, a 128-bit fingerprint of the isomorphism class; `degree_invariant(G)` (work splitting); `complement(G)` |
| `recognize.hpp` | labeled recognition | `recognize_labeled(G)` reads g from column 1 and g·f from column 2, solves for f, rebuilds and compares (O(n²) word operations); `is_riordan_labeled`; `reversed(G)` (i ↦ n+1−i) |
| `labelings.hpp` | unlabeled recognition and class sizes | `riordan_labelings(G, node_limit, first_only)` enumerates all Riordan labelings of an unlabeled graph by backtracking: at most 4 admissible rows per position, twin pruning, capacity pruning. It returns the set C(G) and a completeness flag |
| `properties.hpp` | invariants and pattern search | `edge_count`, `component_count`, `is_connected`, `triangle_count`, `distance_stats` (diameter, Wiener index), `degree_sequence`, `contains_subgraph(H,P,induced)`, `pattern_graph(name)` (named small graphs or graph6) |
| `io.hpp` | formatting | `poly_str`/`parse_poly` (text ↔ polynomial), `canonical_params`, `riordan_matrix`, `edge_list`, `graph6`/`parse_graph6`, `write_svg` |

Dependencies: `riordan_graph.hpp` ← `canon.hpp` (plus nauty), `recognize.hpp`, `labelings.hpp`,
`io.hpp` ← `properties.hpp`. Only `canon.hpp` needs the nauty object files at link time.

## 3. Executables

| program | source | role | reads | writes |
|---|---|---|---|---|
| `rgcount` | `tools/rgcount.cpp` | **Method A** for u(n) | arguments | `N L(N) u(N) (time)`; checkpoint file |
| `tree_count` | `verify/tree_count.cpp` | **Method B** for u(n), independent code | arguments | `N L(N) u(N)`; checkpoint file |
| `recognize_all` | `verify/recognize_all.cpp` | **Method C**: tests every graph for a Riordan labeling | graph6 on stdin (from `geng`) | `N graphs riordan nonriordan`; optional graph6 lists |
| `geng` | `third_party/nauty/geng.c` | nauty's generator of all graphs of order N | arguments | graph6 on stdout |
| `rglist` | `tools/rglist.cpp` | **listing**: all labeled or unlabeled Riordan graphs of order N, or one graph | arguments | text / Markdown / CSV / JSON, SVG drawings |
| `rgsample` | `tools/rgsample.cpp` | **sampling**: random Riordan graphs (labeled, unlabeled, proper, Appell, Bell), structure filters, statistics, estimate of u(N) | arguments | text / CSV / JSON, summary, SVG drawings |
| `complement_labeled` | `analysis/` | counts labeled Riordan graphs whose complement is labeled Riordan | arguments | one line per n |
| `complement_unlabeled` | `analysis/` | unlabeled complement counts, self-complementary graphs | arguments | one line per n |
| `structure` | `analysis/` | minimal non-Riordan graphs, failures of heredity, complements | `data/lists/R_n.g6`, `NR_n.g6` | one line per n; `minimal_n.g6` |
| `labelstats`, `labelstats2` | `analysis/` | distribution of class sizes \|C(H)\| and the reversal statistics | arguments | one block per n |
| `bell` | `analysis/` | io-decomposable Bell-type graphs: diameter, Wiener index, independence number | arguments | text tables |
| `quick_check`, `recognize_test` | `tests/` | unit tests of enumeration and recognition | – | pass/fail |

## 4. Data flow

```
                      ┌──────────── rgcount (A) ────────┐
riordan_graph.hpp ────┤                                 ├──► data/runs/A*,B*  ──► scripts/compute_counts.sh ──► data/counts.csv
 (Enumeration)        └──────────── tree_count (B) ─────┘
geng ──► recognize_all (C) ──► data/runs/C*, data/lists/R_n, NR_n ──► structure ──► minimal_n.g6
rglist   ──► examples/  (scripts/make_examples.sh)        docs/figures/ (README drawings)
rgsample ──► data/sampling/ (scripts/sampling_examples.sh), data/rigidity/
bell     ──► data/bell/ (scripts/bell_all64.sh)
```

## 5. Tests (`make check`)

| test | checks |
|---|---|
| `tests/quick_check.cpp` | the enumeration yields L(n) pairwise distinct labeled graphs (n ≤ 11) and prints u(n) |
| `tests/recognize_test.cpp` | labeled recognition agrees with the enumeration on all 2^(n(n−1)/2) labeled graphs for n ≤ 7, plus random graphs and perturbations for n ≤ 20 |
| `tests/small_counts.sh` | u(1..8) by methods A, B and C |
| `tests/rglist_test.sh` | listing counts (labeled, unlabeled, complement, self-complementary), class sizes summing to L(n), `--show` |
| `tests/rgsample_test.sh` | unlabeled acceptance rate ≈ u(8)/L(8); the estimator ≈ u(9)/L(9); filters respected; reproducible seeds |

The CI workflow (`.github/workflows/ci.yml`) runs the build and `make check` on Linux and macOS,
and checks u(11) by methods A and B.
