/**
 * @file   rgraph/recognize.hpp
 * @brief  Recognition of LABELED Riordan graphs (Problem 4, labeled case), in O(n^2) word
 *         operations.
 *
 * Given a labeled graph G on vertices 1..n (0-based in code), decide whether G = G_n(g, f) for
 * some g, f (G_n(g, f): for i > j, ij is an edge iff [t^(i-2)] g f^(j-1) = 1 over Z_2, f(0) = 0),
 * and if so return the canonical parameters (g mod t^(n-1), f mod t^(n-1-r)).
 *
 * Algorithm (columns 1 and 2 of the adjacency matrix are g and g*f mod t^(n-1)):
 *  1. g_k = adj(k+2, 1) for k = 0..n-2 (paper numbering), i.e. column 1 of the matrix is g.
 *  2. If g = 0, G is Riordan iff G is the null graph (all columns g f^(j-1) vanish).
 *  3. Otherwise let r = val(g). Column 2 is h = g f (mod t^(n-1)); f is determined uniquely
 *     modulo t^(n-1-r) by solving h_(r+k) = sum_a g_a f_(r+k-a) for f_k, k = 1, 2, ...
 *     (g_r = 1). Bits f_k with k >= n-1-r do not influence G_n(g, f).
 *  4. G is Riordan iff G_n(g, f) == G.
 * Cost: O(n^2) bit operations for steps 1-3 and n-1 truncated polynomial products for step 4,
 * i.e. O(n^3) bit operations (O(n^2) word operations for n <= 64).
 *
 * This is the labeled problem only; deciding whether some relabeling of G is Riordan (the
 * unlabeled problem) is done by rgraph/labelings.hpp and verify/recognize_all.cpp.
 * Conventions: bit k of a u64 polynomial = coefficient of t^k. All functions are pure.
 */
#pragma once
#include "riordan_graph.hpp"

#include <optional>
#include <utility>

namespace rgraph {

/** @brief Parameters (g, f) of a labeled Riordan graph G_n(g, f). */
struct RiordanPair {
    u64 g = 0, f = 0;  ///< generating functions (bit k = coefficient of t^k; f(0) = 0)
};

/**
 * @brief Decide whether G (as labeled) equals some G_n(g, f), and find (g, f).
 * @param G  labeled graph, n <= 64, vertex index v = paper vertex v+1
 * @return   the canonical parameters (g mod t^(n-1), f mod t^(n-1-r), r = val(g); g = f = 0 for
 *           the null graph and for n <= 1), or std::nullopt if G is not a labeled Riordan graph
 */
inline std::optional<RiordanPair> recognize_labeled(const Graph& G) {
    const int n = G.n;
    if (n <= 1) return RiordanPair{};
    u64 g = 0;
    for (int k = 0; k <= n - 2; ++k) g |= u64((G.adj[k + 1] >> 0) & 1) << k;  // vertex k+2 vs vertex 1
    if (g == 0) {
        for (int v = 0; v < n; ++v)
            if (G.adj[v]) return std::nullopt;
        return RiordanPair{};
    }
    const int r = std::countr_zero(g);
    u64 h = 0;  // column 2: h_k = adj(vertex k+2, vertex 2), k >= 1
    if (n >= 2)
        for (int k = 1; k <= n - 2; ++k) h |= u64((G.adj[k + 1] >> 1) & 1) << k;
    u64 f = 0;
    for (int k = 1; r + k <= n - 2; ++k) {  // solve for f_k
        int s = 0;
        for (int a = r + 1; a <= r + k - 1; ++a) s ^= int((g >> a) & 1) & int((f >> (r + k - a)) & 1);
        const int want = int((h >> (r + k)) & 1);
        if (want ^ s) f |= u64{1} << k;
    }
    if (!(riordan_graph(n, g, f) == G)) return std::nullopt;
    return RiordanPair{g, f};
}

/**
 * @brief True iff G, with its given labeling, is a labeled Riordan graph.
 * @param G  labeled graph, n <= 64
 */
inline bool is_riordan_labeled(const Graph& G) { return recognize_labeled(G).has_value(); }

/**
 * @brief The same graph with the labels reversed (paper vertex i -> n+1-i; index v -> n-1-v).
 * @param G  graph, n <= 64
 * @return   the relabeled graph (isomorphic to G)
 */
inline Graph reversed(const Graph& G) {
    Graph R;
    R.n = G.n;
    for (int v = 0; v < G.n; ++v)
        for (u64 x = G.adj[v]; x; x &= x - 1) R.adj[G.n - 1 - v] |= u64{1} << (G.n - 1 - std::countr_zero(x));
    return R;
}

}  // namespace rgraph
