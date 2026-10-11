/**
 * @file   rgraph/riordan_graph.hpp
 * @brief  Labeled Riordan graphs G_n(g, f) for n <= 64: polynomial arithmetic over Z_2, the graph
 *         type, construction of G_n(g, f), and the enumeration of all labeled Riordan graphs.
 *
 * Definition (Cheon, Jung, Kitaev, Mojallal, "Riordan graphs I", LAA 579 (2019), Def. 2.4)
 * -----------------------------------------------------------------------------------------
 * G_n(g, f) has vertices 1..n; for i > j >= 1, ij is an edge iff [t^(i-2)] g f^(j-1) = 1 (mod 2),
 * where g, f are power series over Z_2 = GF(2) and f(0) = 0 (f'(0) may be 0 or 1, g may be 0).
 * Only g mod t^(n-1) and f mod t^(n-1) influence the graph. There are
 * L(n) = (4^(n-1) + 2) / 3 distinct labeled Riordan graphs of order n (Prop. 2.8).
 *
 * Conventions
 * -----------
 * - A polynomial over Z_2 is a 64-bit word `u64`; bit k is the coefficient of t^k.
 * - A graph is an array of n adjacency rows; row v is the bitmask of the neighbours of vertex v.
 * - Vertices are 0-based in code: vertex i of the paper (1..n) is index i-1.
 *
 * Everything is header-only. Functions are pure, except that riordan_graph() and
 * Enumeration::pair() throw on invalid arguments.
 */
#pragma once
#include <bit>
#include <cstdint>
#include <stdexcept>

namespace rgraph {

/// 64-bit word used both for Z_2 polynomials (bit k = coefficient of t^k) and vertex bitmasks.
using u64 = std::uint64_t;

/**
 * @brief Mask of the lowest `bits` bits, i.e. the truncation "mod t^bits".
 * @param bits  number of bits (0 gives 0; values >= 64 give all ones)
 * @return      (1 << bits) - 1, saturated at 64 bits
 */
inline u64 low_mask(int bits) { return bits >= 64 ? ~u64{0} : ((u64{1} << bits) - 1); }

/**
 * @brief Product of two Z_2 polynomials truncated to the bits in `mask`.
 *
 * Bits of the untruncated product beyond t^63 are lost, so the result is exact whenever `mask`
 * only selects bits below 64 (always true here).
 * @param a, b  factors (bit k = coefficient of t^k)
 * @param mask  bits to keep, typically low_mask(m) for "mod t^m"
 * @return      (a * b) & mask over Z_2
 */
inline u64 mul_trunc(u64 a, u64 b, u64 mask) {
    u64 r = 0;
    for (; b; b &= b - 1) r ^= a << std::countr_zero(b);
    return r & mask;
}

/**
 * @brief Simple undirected graph on n <= 64 vertices (0-based), stored as adjacency bitmasks.
 *
 * Invariants maintained by all producers in this library: adj is symmetric, has no loops
 * (bit v of adj[v] is 0), and rows v >= n as well as bits >= n are zero.
 */
struct Graph {
    int n = 0;          ///< number of vertices (0 <= n <= 64)
    u64 adj[64] = {};   ///< adj[v]: bit w set iff vertices v and w are adjacent
    /** @brief True iff u and v are adjacent (0-based, u, v < n). */
    bool has_edge(int u, int v) const { return (adj[u] >> v) & 1; }
    /** @brief Degree of vertex v (0-based). */
    int degree(int v) const { return std::popcount(adj[v]); }
    /** @brief Equality as labeled graphs: same order and identical rows 0..n-1. */
    bool operator==(const Graph& o) const {
        if (n != o.n) return false;
        for (int v = 0; v < n; ++v)
            if (adj[v] != o.adj[v]) return false;
        return true;
    }
};

/**
 * @brief The labeled Riordan graph G_n(g, f). Only g mod t^(n-1) and f mod t^(n-1) matter.
 *
 * Column j (paper numbering) of the strictly lower adjacency triangle is g f^(j-1) mod t^(n-1);
 * the columns are produced by repeated truncated multiplication by f, stopping early once a column
 * vanishes. Cost: at most n-1 truncated products, i.e. O(n^2) word operations.
 * @param n  order, 1 <= n <= 64
 * @param g  first generating function (any value, including 0)
 * @param f  second generating function; must satisfy f(0) = 0 (bit 0 clear)
 * @return   the graph, with vertex i of the paper at index i-1
 * @throws std::invalid_argument if n is out of range or f(0) = 1
 */
inline Graph riordan_graph(int n, u64 g, u64 f) {
    if (n < 1 || n > 64) throw std::invalid_argument("riordan_graph: need 1 <= n <= 64");
    if (f & 1) throw std::invalid_argument("riordan_graph: need f(0) = 0");
    Graph G;
    G.n = n;
    const u64 mask = low_mask(n - 1);  // coefficients t^0 .. t^(n-2)
    u64 col = g & mask;                // g f^(j-1) for the current vertex j (0-based j0 = j-1)
    for (int j0 = 0; j0 + 1 < n && col; ++j0) {
        // bit k of col  <->  edge between vertex j0 and vertex k+1  (paper: i = k+2, j = j0+1);
        // val(col) >= j0 because val(f) >= 1, so all k+1 > j0 automatically.
        const u64 lower = col << 1;
        G.adj[j0] |= lower;
        for (u64 x = lower; x; x &= x - 1) G.adj[std::countr_zero(x)] |= u64{1} << j0;
        col = mul_trunc(col, f, mask);
    }
    return G;
}

/**
 * @brief Number of labeled Riordan graphs of order n: L(n) = (4^(n-1) + 2) / 3 (Prop. 2.8).
 * @param n  order; exact for 1 <= n <= 32 (the value must fit in 64 bits)
 * @return   L(n), or 0 for n < 1
 */
inline u64 labeled_count(int n) { return n < 1 ? 0 : ((u64{1} << (2 * (n - 1))) + 2) / 3; }

/**
 * @brief Enumeration of every labeled Riordan graph of order n exactly once, following the proof
 *        of Prop. 2.8.
 *
 * Order: the null graph (g = 0), then for r = 0..n-2 (the lowest term t^r of g) the
 * 4^(n-2-r) choices of the free bits g_(r+1..n-2) and f_(1..n-2-r). The graphs are numbered
 * 0 .. labeled_count(n)-1; index 0 is the null graph. The (g, f) produced are the canonical
 * parameters (g mod t^(n-1), f mod t^(n-1-r)), so distinct indices give distinct graphs.
 * Used by method A (tools/rgcount.cpp). Requires 1 <= n <= 32.
 */
struct Enumeration {
    int n;  ///< order of the graphs
    /** @brief Enumeration of the labeled Riordan graphs of order n_ (1 <= n_ <= 32). */
    explicit Enumeration(int n_) : n(n_) {}
    /** @brief Number of graphs, L(n). */
    u64 size() const { return labeled_count(n); }
    /**
     * @brief (g, f) of graph number idx.
     *
     * Within the block of a given r, the low k = n-2-r bits of the local index are
     * g_(r+1..n-2) and the next k bits are f_(1..k).
     * @param idx  index, 0 <= idx < size()
     * @param g, f out: canonical parameters of the graph
     * @throws std::out_of_range if idx >= size()
     */
    void pair(u64 idx, u64& g, u64& f) const {
        if (idx == 0) { g = 0; f = 0; return; }
        --idx;
        for (int r = 0; r <= n - 2; ++r) {
            const int k = n - 2 - r;
            const u64 block = u64{1} << (2 * k);
            if (idx < block) {
                const u64 gb = idx & low_mask(k), fb = idx >> k;
                g = (u64{1} << r) | (gb << (r + 1));
                f = fb << 1;
                return;
            }
            idx -= block;
        }
        throw std::out_of_range("Enumeration::pair");
    }
    /**
     * @brief Graph number idx, i.e. riordan_graph(n, g, f) for (g, f) = pair(idx).
     * @param idx  index, 0 <= idx < size()
     */
    Graph graph(u64 idx) const {
        u64 g, f;
        pair(idx, g, f);
        return riordan_graph(n, g, f);
    }
};

}  // namespace rgraph
