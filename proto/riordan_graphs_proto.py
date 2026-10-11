"""Reference prototype for Riordan graphs (Cheon, Jung, Kitaev, Mojallal, LAA 579 (2019)).

Written directly from Definition 2.4, as simply as possible, to serve as an independent
check of the fast C++ code. Not optimised.

Definition 2.4: vertices 1..n; for i > j >= 1,  ij is an edge  iff  [t^(i-2)] g f^(j-1) = 1 (mod 2),
with g, f power series over GF(2) and f(0) = 0. There are L(n) = (4^(n-1)+2)/3 labeled Riordan
graphs; u(n) is the number of isomorphism classes (unlabeled Riordan graphs).

Method: all pairs (g mod t^(n-1), f mod t^(n-1)) are enumerated by brute force (2^(2n-3) pairs),
the distinct edge sets give the labeled graphs, and the unlabeled count uses a brute-force
canonical form (minimum over all n! relabelings). No external packages are needed.

Usage:
    python3 proto/riordan_graphs_proto.py      # prints n = 1..7 (about 20 s; n! growth)

    or, as a module (used by proto/check_structure.py; proto/ must be on sys.path):
    from riordan_graphs_proto import riordan_graph, all_labeled_riordan_graphs

Inputs: none. Outputs (stdout, as a script), one line per n = 1..7:
    "<n> labeled: <number of distinct labeled graphs> formula: <L(n)> unlabeled: <u(n)>"
Expected u(n): 1 2 4 11 30 109 406, as in docs/COMPUTATIONS.md (section 5 lists this check).

Graphs are represented as frozensets of edges (j, i) with 1 <= j < i <= n.
"""
from itertools import permutations, product


def poly_mul(a, b, n):
    """Product of two GF(2) power series, truncated to n terms.

    Args:
        a: coefficient list of the first factor (a[k] = [t^k] a, values 0/1).
        b: coefficient list of the second factor.
        n: number of coefficients to keep.

    Returns:
        List c of length n with c[k] = [t^k] (a * b) mod 2.
    """
    c = [0] * n
    for i, x in enumerate(a[:n]):
        if x:
            for j, y in enumerate(b[: n - i]):
                c[i + j] ^= y
    return c


def riordan_graph(g, f, n):
    """Edge set of the labeled Riordan graph G_n(g, f), directly from Definition 2.4.

    Args:
        g: coefficient list of g over GF(2) (only the first n-1 entries are used; shorter lists
           are padded with zeros).
        f: coefficient list of f over GF(2); f[0] must be 0 (asserted).
        n: number of vertices.

    Returns:
        frozenset of edges (j, i), 1 <= j < i <= n, with ij an edge iff [t^(i-2)] g f^(j-1) = 1.
    """
    assert f[0] == 0
    m = max(n - 1, 0)
    edges = set()
    col = list(g[:m]) + [0] * (m - len(g[:m]))  # g f^(j-1), starting at j = 1
    for j in range(1, n + 1):
        for i in range(j + 1, n + 1):
            if i - 2 < m and col[i - 2]:
                edges.add((j, i))
        col = poly_mul(col, f, m)
    return frozenset(edges)


def all_labeled_riordan_graphs(n):
    """Every labeled Riordan graph on n vertices, by brute force over all (g, f) mod t^(n-1).

    Args:
        n: number of vertices (n >= 1).

    Returns:
        set of edge frozensets (as from riordan_graph); its size should be L(n) = (4^(n-1)+2)/3.
    """
    m = max(n - 1, 1)
    seen = set()
    for gbits in product((0, 1), repeat=m):
        for fbits in product((0, 1), repeat=m - 1):
            seen.add(riordan_graph(list(gbits), [0] + list(fbits), n))
    return seen


def canonical(edges, n):
    """Canonical form of a graph: the lexicographically smallest sorted edge list over all relabelings.

    Brute force over all n! permutations, so only for small n (n <= 8).

    Args:
        edges: iterable of edges (a, b) with vertices in 1..n.
        n: number of vertices.

    Returns:
        tuple of sorted edge pairs; equal for two graphs iff they are isomorphic.
    """
    best = None
    for p in permutations(range(1, n + 1)):
        e = tuple(sorted(tuple(sorted((p[a - 1], p[b - 1]))) for a, b in edges))
        if best is None or e < best:
            best = e
    return best


def unlabeled_count(n):
    """Number u(n) of unlabeled Riordan graphs on n vertices (brute force).

    Args:
        n: number of vertices.

    Returns:
        Number of distinct canonical forms among all labeled Riordan graphs on n vertices.
    """
    return len({canonical(e, n) for e in all_labeled_riordan_graphs(n)})


if __name__ == "__main__":
    for n in range(1, 8):
        lab = all_labeled_riordan_graphs(n)
        print(n, "labeled:", len(lab), "formula:", (4 ** (n - 1) + 2) // 3, "unlabeled:", unlabeled_count(n), flush=True)
