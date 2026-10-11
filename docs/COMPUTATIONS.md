# Computations: methods, verification and data

This document describes how every number in this repository was computed, how it was checked, and
where the raw data are stored. It deliberately contains no mathematical proofs. It is a record of
the computations.

Notation: `G_n(g,f)` is the labeled Riordan graph on vertices 1..n (see the README). For i > j,
`ij` is an edge iff `[t^(i-2)] g f^(j-1) = 1` over Z₂. `L(n) = (4^(n-1)+2)/3` is the number of labeled
Riordan graphs, and `u(n)` the number of unlabeled ones (isomorphism classes).

---

## 1. Enumeration of labeled Riordan graphs

A labeled Riordan graph depends only on `g mod t^(n-1)` and, if `g ≠ 0` with lowest term `t^r`, on
`f mod t^(n-1-r)`. Every labeled Riordan graph therefore arises exactly once from:

* the null graph (`g = 0`), and
* for `r = 0..n-2` and `k = n-2-r`: `g = t^r + g_(r+1) t^(r+1) + … + g_(n-2) t^(n-2)` and
  `f = f_1 t + … + f_k t^k`, with free bits (`4^k` choices).

`include/rgraph/riordan_graph.hpp` (`struct Enumeration`) numbers these graphs `0 .. L(n)-1`.
`tests/quick_check` confirms for n ≤ 11 that the L(n) generated graphs are pairwise distinct.
`tests/recognize_test` confirms for n ≤ 7 that they are exactly the graphs accepted by the
independent recognizer, out of all `2^(n(n-1)/2)` labeled graphs.

## 2. Three methods for u(n)

| | Method A — `build/rgcount` | Method B — `build/tree_count` | Method C — `build/geng` + `build/recognize_all` |
|---|---|---|---|
| Source | `tools/rgcount.cpp`, `include/rgraph/` | `verify/tree_count.cpp` (shares no code with A) | `verify/recognize_all.cpp` |
| Generates | the L(n) labeled Riordan graphs, from the enumeration above | the L(n) labeled Riordan graphs, vertex by vertex: each new vertex fixes one coefficient of g and at most one of f, so the search tree has ≤ 4 children per node | all graphs of order n (nauty `geng`, one per isomorphism class) |
| Decides isomorphism by | nauty canonical form, initial partition by degree; 128-bit hash | nauty canonical form, unit partition; different 128-bit hash | no canonical forms: an exact search for a Riordan labeling of each graph |
| Work split | passes by hash of the sorted degree sequence | passes by sorted (degree, triangles) pairs; subtrees claimed by threads | – |
| Built-in check | number of generated labeled graphs = L(n) | same | number of graphs = OEIS A000088 |
| Used for | n ≤ 17 | n ≤ 17 | n ≤ 10 |

Methods A and B are separately written programs. Both use nauty for canonical forms, but with
different initial partitions and different hash functions. Method C uses neither generation nor
canonical forms, so it is independent of both. A collision of two isomorphism classes under a
128-bit fingerprint would make A and B disagree, unless the same two classes also collided under
the other hash function. With fewer than 10⁹ classes per order, the probability of an undetected
collision is negligible.

### Results

| n | L(n) | u(n) | A | B | C |
|---:|---:|---:|:-:|:-:|:-:|
| 1 | 1 | 1 | ✓ | ✓ | ✓ |
| 2 | 2 | 2 | ✓ | ✓ | ✓ |
| 3 | 6 | 4 | ✓ | ✓ | ✓ |
| 4 | 22 | 11 | ✓ | ✓ | ✓ |
| 5 | 86 | 30 | ✓ | ✓ | ✓ |
| 6 | 342 | 109 | ✓ | ✓ | ✓ |
| 7 | 1 366 | 406 | ✓ | ✓ | ✓ |
| 8 | 5 462 | 1 864 | ✓ | ✓ | ✓ |
| 9 | 21 846 | 8 169 | ✓ | ✓ | ✓ |
| 10 | 87 382 | 36 897 | ✓ | ✓ | ✓ |
| 11 | 349 526 | 158 351 | ✓ | ✓ | |
| 12 | 1 398 102 | 671 656 | ✓ | ✓ | |
| 13 | 5 592 406 | 2 657 424 | ✓ | ✓ | |
| 14 | 22 369 622 | 11 481 042 | ✓ | ✓ | |
| 15 | 89 478 486 | 46 683 483 | ✓ | ✓ | |
| 16 | 357 913 942 | 190 658 218 | ✓ | ✓ | |
| 17 | 1 431 655 766 | 789 399 022 | ✓ | ✓ | |

The machine-readable summary is `data/counts.csv` and the run logs are `data/runs/`:

* `A<n>.txt`, `B<n>.txt`: final line `n L(n) u(n)`;
* `*.ckpt`: per-pass checkpoints, `n passes pass classes labeled`. For n = 17 both methods used 12
  passes, and the per-pass class counts sum to the total;
* `C<n>.txt`: `n graphs riordan nonriordan`. `C<n>.err` gives search-node statistics.

### Cost (measured on a 2-core cloud VM, 8 GB RAM)

| n | L(n) | Method A | Method B | memory per pass |
|---:|---:|---|---|---|
| 14 | 2.2·10⁷ | 1.5 min (1 thread) | 1.5 min | 0.2 GB |
| 15 | 8.9·10⁷ | ≈ 9 min (1 thread) | ≈ 10 min | 0.8 GB |
| 16 | 3.6·10⁸ | ≈ 41 min (1 thread, 4 passes) | ≈ 50 min | 0.8 GB |
| 17 | 1.4·10⁹ | ≈ 3.8 h (2 threads, 12 passes) | ≈ 3 h | 1–2 GB |
| 18 | 5.7·10⁹ | estimated about 18 h on 2 cores per method; ≈ 3·10⁹ classes; use `--passes 48` or more RAM | | |

Throughput is about 3–4 µs per labeled graph per core, dominated by nauty. Every pass regenerates
all labeled graphs, which costs about 0.4 µs each, so fewer passes on a machine with more memory
is faster.

## 3. Method C details and recognition

`recognize_all` decides for each input graph whether some ordering of its vertices is a labeled
Riordan graph. It fills positions 1, 2, 3, … in turn. At position m, the coefficients found so far
leave at most four admissible adjacency patterns (one new bit of g, at most one new bit of f), so
only matching vertices are tried. Vertices that are twins (equal neighbourhoods apart from each
other) are tried only once.

Search nodes per graph, over all graphs of order n:

| n | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|
| mean | 15.9 | 64.5 | 305.7 | 1 028 | 2 727 | 6 013 |
| max | 59 | 336 | 1 069 | 4 019 | 16 886 | 81 514 |

Labeled recognition (`include/rgraph/recognize.hpp`) works as follows. Read g from column 1 and
g·f from column 2 of the adjacency matrix, solve for f, rebuild `G_n(g,f)` and compare it with the
input. The `rglist` tool uses this to decide whether a complement is labeled Riordan.

## 4. Further sequences

All of the following are produced by the programs in `analysis/` and `tools/rglist`. Every entry
was cross-checked by at least two programs; the checks are listed after the tables.

### 4.1 Riordan graphs among all graphs (method C, `build/structure`)

| n | graphs | Riordan u(n) | non-Riordan | minimal non-Riordan¹ | connected Riordan | Riordan trees / all trees | Riordan forests |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 1 | 0 | 0 | 1 | 1/1 | 1 |
| 2 | 2 | 2 | 0 | 0 | 1 | 1/1 | 2 |
| 3 | 4 | 4 | 0 | 0 | 2 | 1/1 | 3 |
| 4 | 11 | 11 | 0 | 0 | 6 | 2/2 | 6 |
| 5 | 34 | 30 | 4 | 4 | 19 | 3/3 | 10 |
| 6 | 156 | 109 | 47 | 17 | 78 | 6/6 | 20 |
| 7 | 1 044 | 406 | 638 | 111 | 316 | 10/11 | 36 |
| 8 | 12 346 | 1 864 | 10 482 | 123 | 1 502 | 19/23 | 70 |
| 9 | 274 668 | 8 169 | 266 499 | 232 | 6 969 | 40/47 | 140 |
| 10 | 12 005 168 | 36 897 | 11 968 271 | 303 | 32 382 | 60/106 | 253 |

¹ A non-Riordan graph is *minimal* if each of its vertex-deleted subgraphs is Riordan. The lists are
in `data/lists/minimal_<n>.g6`.

The four non-Riordan graphs of order 5 are C₄ ∪ K₁, K₄ ∪ K₁, and the complements of the chair and of
K₁,₃ ∪ K₁. The smallest non-Riordan tree is the spider S(2,2,2) on 7 vertices.

The class is not closed under vertex deletion. The number of Riordan graphs having a vertex whose
deletion gives a non-Riordan graph is 0, 23, 198, 1510, 7529, 35949 for n = 5..10. The smallest
example is G₆(1+t, t²+t³), with edges 12 13 24 26 36: deleting vertex 4 leaves C₄ ∪ K₁. The second is
G₆(1, t+t²), with edges 12 23 24 34 36 45 46 56: deleting vertex 5 leaves the complement of the chair.

### 4.2 Complements

| n | labeled: G and its complement both labeled Riordan, c(n) | unlabeled: complement also Riordan, u_c(n) | self-complementary Riordan graphs |
|---:|---:|---:|---:|
| 1 | 1 | 1 | 1 |
| 2 | 2 | 2 | 0 |
| 3 | 4 | 4 | 0 |
| 4 | 12 | 11 | 1 |
| 5 | 28 | 26 | 2 |
| 6 | 92 | 90 | 0 |
| 7 | 220 | 230 | 0 |
| 8 | 732 | 707 | 1 |
| 9 | 1 756 | 1 591 | 1 |
| 10 | 5 852 | 3 804 | 0 |
| 11 | 14 044 | 9 438 | 0 |
| 12 | 46 812 | 26 326 | 2 |
| 13 | 112 348 | 64 544 | 2 |
| 14 | 374 492 | 201 696 | 0 |

The tools are `build/complement_labeled 14`, `build/complement_unlabeled 14`, and
`build/rglist n --filter complement|selfcomplementary --count`. The unlabeled counts for n ≤ 10
agree with method C lists (`build/structure`). The labeled counts agree with an independent brute
force for n ≤ 11.

### 4.3 How many labelings does an unlabeled Riordan graph have?

For a labeled Riordan graph H, let C(H) be the set of labeled Riordan graphs isomorphic to H. Then
u(n) = Σ 1/|C(H)| over all labeled H. `build/labelstats` and `build/labelstats2` report the
distribution of |C(H)|.

The split used below is by `f'(0)`. Graphs with `f'(0) = 1` are reversible: renumbering the vertices
`i ↦ n+1−i` gives another labeled Riordan graph. They make up about half of all labeled graphs.

| n | 10 | 11 | 12 | 13 | 14 | 15 |
|---|---|---|---|---|---|---|
| u(n)/L(n) | 0.422 | 0.453 | 0.480 | 0.475 | 0.513 | 0.522 |
| f'(0)=1: fraction in classes of size exactly 2 | 0.84 | 0.94 | 0.94 | 0.98 | 0.98 | 0.99 |
| f'(0)=0: fraction alone in their class | 0.22 | 0.25 | 0.35 | 0.33 | 0.42 | 0.43 |
| graphs equal to their own reversal | 1 942 | 5 526 | 9 622 | 32 150 | 64 918 | 187 798 |

Random samples with exact class sizes, for graphs with `f'(0)=1` only. The data are in
`data/rigidity/`; see `data/rigidity/README.md` for provenance. The same experiment is run with
`build/rgsample n --samples S --fprime 1 --class-size --node-limit 20000000`. Searches over
2·10⁷ nodes were stopped.

| n | completed samples | \|C(H)\| ≤ 2 |
|---:|---:|---:|
| 18 | 398 | 398 |
| 20 | 397 | 395 |
| 24 | 292 | 291 |
| 28 | 180 | 180 |

The ratio `u(n)/4^(n−2)` is 0.634, 0.684, 0.696, 0.710, 0.735 for n = 13..17.

### 4.4 io-decomposable Riordan graphs of the Bell type (`build/bell`)

Problem 6 and Conjecture 2 of Cheon et al. concern the io-decomposable graphs `G_n(g, tg)`, for which
`[t^j]g = [t^(2j+1)]g`. There are `2^(⌈n/2⌉−1)` of them. Data (`data/bell/`):

* diameter distribution and Wiener-index range over all of them for n ≤ 34, n = 40 and n = 64
  (2³¹ graphs, `scripts/bell_all64.sh`, 32 resumable slices);
* the Catalan graph `CG_n` has diameter k at n = 2^k for every k ≤ 13. Among all io-decomposable
  graphs of order 2^k it is the only one with diameter k for k = 4, 5, 6, but not for k = 2
  (PG₄ also has diameter 2) or k = 3 (G₈(1+t+t³+t⁶, tg) also has diameter 3);
* Wiener index and diameter of the Pascal and Catalan graphs for n ≤ 64 and n = 2^k ≤ 8192;
* independence number ⌊n/2⌋ for every io-decomposable graph with 2 ≤ n ≤ 34.

## 5. Sampling beyond the exhaustive range (`build/rgsample`)

Exact enumeration stops at n = 17. `rgsample` draws Riordan graphs at random. Each model is exactly
uniform on its set:

* `labeled`: uniform over all labeled Riordan graphs. Canonical parameters are drawn by index for
  n ≤ 32, and by r ∝ 4^(n−2−r) for n > 32.
* `proper`, `appell`, `bell`: uniform over those families.
* `unlabeled`: uniform over isomorphism classes. A labeled graph H is accepted with probability
  1/|C(H)|.

Filters are applied by rejection: `--fprime`, `--val-g`, `--contains`, `--avoids`,
`--contains-induced`, `--avoids-induced`, `--connected`, `--edges`. With `--class-size` the tool
computes `|C(H)|` exactly by enumerating all Riordan labelings (`include/rgraph/labelings.hpp`), and
estimates `u(n) = L(n)·E[1/|C(H)|]`.

### Validation against exact values

| check | result |
|---|---|
| labeled model, n = 5: all 86 labeled graphs, 43 000 samples | χ² = 98.9 with 85 degrees of freedom (consistent with uniform) |
| unlabeled model, n = 6: all 109 classes, 21 800 samples | χ² = 79.2 with 108 degrees of freedom (consistent with uniform) |
| unlabeled acceptance rate, n = 8 | 0.3430 (20 000 samples); exact u(8)/L(8) = 0.3413 |
| estimator, n = 10 (2 × 5000 samples) | 0.4175 and 0.4156 ± 0.008; exact u(10)/L(10) = 0.4222 |
| estimator, n = 12 (2000 samples) | u(12) ∈ [657 019, 694 253]; exact 671 656 |
| estimator, n = 14 (1000 samples) | u(14) ∈ [1.128·10⁷, 1.212·10⁷]; exact 1.148·10⁷ |
| class sizes from the labeling search, all labeled graphs n ≤ 10 | identical to the fingerprint classes; Σ 1/\|C(H)\| = u(n) exactly |
| subgraph / induced-subgraph filters, 5 232 accepted graphs (n = 8, 10 patterns) | 0 disagreements with networkx |

### Typical labeled Riordan graphs of larger order (5000 samples each, `data/sampling/`)

| n | edges (mean) | max degree (mean) | triangles (mean) | connected | diameter (mean, connected) | contains K₄ |
|---:|---:|---:|---:|---:|---:|---:|
| 20 | 70.7 | 12.2 | 94 | 0.803 | 2.73 | 0.711 |
| 32 | 181.1 | 19.4 | 384 | 0.865 | 2.62 | 0.857 |
| 48 | 406.2 | 28.7 | 1 292 | 0.911 | 2.61 | 0.927 |
| 64 | 711.7 | 37.7 | 2 961 | 0.936 | 2.62 | 0.963 |

The mean number of edges is roughly n²/6 (slightly above). Typical graphs become connected and have diameter 2–3.
For 100 random graphs with `f'(0) = 1` at n = 22 (`data/sampling/fprime1_classes_n22.txt`), every
class was exactly `{H, rev H}`.

### Limits

* Rejection filters cost 1/(probability of the filter) draws per accepted sample. For example,
  triangle-free graphs of order 30 have probability about 10⁻⁵.
* The labeling search is fast for graphs with `f'(0) = 1`. Graphs with `f'(0) = 0` and `val f ≥ 3`
  can need many millions of search nodes from n ≈ 15 on.
* Searches over `--node-limit` are reported. The u(n) estimate is then given as a conservative
  interval, and the `unlabeled` model skips those draws, which the summary flags as a possible bias.
* For n ≥ 16, use `--fprime 1`, or accept wider intervals.

## 6. Cross-checks performed

| Quantity | Checked by |
|---|---|
| u(n), n ≤ 17 | methods A and B (independent code, different hashes and invariants) |
| u(n), n ≤ 10 | method C (no generation, no canonical forms) |
| u(n), n ≤ 7 | Python reference implementation `proto/riordan_graphs_proto.py` (brute-force isomorphism) |
| u(n), n ≤ 14 | `rglist --count` and `complement_unlabeled` (re-derived from the same fingerprints as method A; a consistency check of the tools, not an independent method) |
| labeled enumeration | `tests/quick_check` (distinctness, n ≤ 11), `tests/recognize_test` (all labeled graphs, n ≤ 7) |
| complement counts | `complement_labeled` vs an independent brute force (n ≤ 11); `complement_unlabeled` vs method C lists (n ≤ 10) |
| structure table (4.1) | C++ `structure` vs a separate networkx/labelg implementation (n ≤ 8) |
| examples in the README | regenerated by `scripts/make_examples.sh` and checked against an independent networkx script |
