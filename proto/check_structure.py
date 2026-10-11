"""Prototype checks of the enumeration of labeled Riordan graphs used by the C++ code.

A labeled Riordan graph G_n(g, f) has vertices 1..n; for i > j, ij is an edge iff
[t^(i-2)] g f^(j-1) = 1 over GF(2), with f(0) = 0. There are L(n) = (4^(n-1)+2)/3 of them.

Two checks, against the brute-force set of proto/riordan_graphs_proto.py:

1. n = 1..8: the bijective enumeration (null graph; or g = t^r + free higher bits up to t^(n-2),
   f = free bits f_1..f_(n-2-r)), which is the one implemented by rgraph::Enumeration in
   include/rgraph/riordan_graph.hpp, yields L(n) pairwise distinct graphs, equal to the
   brute-force set.
2. m = 3..8: for each labeled Riordan graph on m-1 vertices (the "prefix"), the possible adjacency
   rows of vertex m in labeled Riordan graphs on m vertices extending it. Prints the maximum number
   of such extensions per prefix and whether they differ only in the adjacency to vertices 1 and 2.

Usage:
    python3 proto/check_structure.py      (imports riordan_graphs_proto from its own
                                           directory; takes under a second)

Inputs: none. Outputs (stdout):
    "n=<n>: bijection enumeration gives <count> graphs, all distinct and equal to brute-force set: <True|False>"
    "m=<m>: max extensions per prefix = <k>; extensions differ only in adjacency to vertices 1,2: <True|False>"
"""
from collections import defaultdict
from riordan_graphs_proto import riordan_graph, all_labeled_riordan_graphs

def bijection_graphs(n):
    """Bijective enumeration of the labeled Riordan graphs on n vertices (as in Prop. 2.8).

    The null graph, and for r = 0..n-2 with k = n-2-r: g = t^r + free bits g_(r+1)..g_(n-2),
    f = free bits f_1..f_k (4^k choices for each r).

    Args:
        n: number of vertices (n >= 1).

    Returns:
        list of edge frozensets, one per parameter choice (duplicates kept, so that
        distinctness can be tested); its length is L(n) = (4^(n-1)+2)/3.
    """
    out = [frozenset()]
    for r in range(0, n - 1):
        k = n - 2 - r
        for x in range(4 ** k):
            gb, fb = x & ((1 << k) - 1), x >> k
            g = [0] * r + [1] + [(gb >> i) & 1 for i in range(k)]
            f = [0] + [(fb >> i) & 1 for i in range(k)] + [0] * (r + 1)
            out.append(riordan_graph(g, f, n))
    return out

for n in range(1, 9):
    lst = bijection_graphs(n)
    ok = len(lst) == len(set(lst)) == (4 ** (n - 1) + 2) // 3 and set(lst) == all_labeled_riordan_graphs(n)
    print(f"n={n}: bijection enumeration gives {len(lst)} graphs, all distinct and equal to brute-force set: {ok}")

# Two free entries: extensions of each labeled Riordan graph on m-1 vertices to m vertices
for m in range(3, 9):
    ext = defaultdict(set)
    for e in all_labeled_riordan_graphs(m):
        prefix = frozenset((a, b) for a, b in e if b < m)
        row = frozenset(a for a, b in e if b == m)
        ext[prefix].add(row)
    worst = max(len(v) for v in ext.values())
    differ_only_12 = all(len({r - {1, 2} for r in v}) == 1 for v in ext.values())
    print(f"m={m}: max extensions per prefix = {worst}; extensions differ only in adjacency to vertices 1,2: {differ_only_12}")
