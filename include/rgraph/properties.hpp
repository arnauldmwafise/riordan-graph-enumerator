/**
 * @file   rgraph/properties.hpp
 * @brief  Elementary graph invariants and (induced) subgraph search for graphs with n <= 64.
 *
 * Graphs are rgraph::Graph values: adjacency rows as 64-bit masks, 0-based vertices. Every function
 * is header-only and pure. Bit-parallel breadth-first search makes the distance computations
 * O(n^2) word operations per graph.
 *
 * Subgraph search (contains_subgraph) asks whether a small pattern P (k <= 10 vertices) occurs in a
 * host graph H. It looks for an injective map phi: V(P) -> V(H) such that
 *  - (subgraph)          every edge uv of P is mapped to an edge phi(u)phi(v) of H;
 *  - (induced subgraph)  in addition, every non-edge of P is mapped to a non-edge of H.
 * It uses backtracking over pattern vertices in BFS order, keeping candidate sets as bitmasks. The
 * worst case is exponential in k, which is fine for the small patterns used here.
 */
#pragma once
#include "io.hpp"
#include "riordan_graph.hpp"

#include <algorithm>
#include <bit>
#include <map>
#include <string>
#include <vector>

namespace rgraph {

/** @brief Number of edges of G. */
inline int edge_count(const Graph& G) {
    int m = 0;
    for (int v = 0; v < G.n; ++v) m += std::popcount(G.adj[v]);
    return m / 2;
}

/** @brief Number of connected components (isolated vertices count as components). */
inline int component_count(const Graph& G) {
    u64 left = low_mask(G.n);
    int c = 0;
    while (left) {
        u64 seen = left & -left, fr = seen;
        while (fr) {
            u64 nx = 0;
            for (u64 x = fr; x; x &= x - 1) nx |= G.adj[std::countr_zero(x)];
            fr = nx & ~seen;
            seen |= nx;
        }
        left &= ~seen;
        ++c;
    }
    return c;
}

/** @brief True iff G is connected (the graph on 0 or 1 vertices counts as connected). */
inline bool is_connected(const Graph& G) { return G.n <= 1 || component_count(G) == 1; }

/** @brief Number of triangles: the sum over vertices of the edges among their neighbours, divided by 3. */
inline long long triangle_count(const Graph& G) {
    long long t = 0;
    for (int v = 0; v < G.n; ++v)
        for (u64 x = G.adj[v] & ~low_mask(v + 1); x; x &= x - 1) {
            const int w = std::countr_zero(x);
            t += std::popcount(G.adj[v] & G.adj[w] & ~low_mask(w + 1));
        }
    return t;
}

/**
 * @brief Diameter and Wiener index by bit-parallel BFS from every vertex.
 * @param G         graph
 * @param diameter  out: maximum distance, or -1 if G is disconnected
 * @param wiener    out: sum of distances over unordered pairs (connected graphs only, else 0)
 */
inline void distance_stats(const Graph& G, int& diameter, long long& wiener) {
    diameter = 0;
    wiener = 0;
    const u64 all = low_mask(G.n);
    for (int s = 0; s < G.n; ++s) {
        u64 seen = u64{1} << s, fr = seen;
        int d = 0;
        while (seen != all) {
            u64 nx = 0;
            for (u64 x = fr; x; x &= x - 1) nx |= G.adj[std::countr_zero(x)];
            nx &= ~seen;
            if (!nx) { diameter = -1; wiener = 0; return; }
            ++d;
            wiener += (long long)d * std::popcount(nx);
            seen |= nx;
            fr = nx;
        }
        diameter = std::max(diameter, d);
    }
    wiener /= 2;
}

/** @brief Degree sequence in non-increasing order. */
inline std::vector<int> degree_sequence(const Graph& G) {
    std::vector<int> d(G.n);
    for (int v = 0; v < G.n; ++v) d[v] = G.degree(v);
    std::sort(d.rbegin(), d.rend());
    return d;
}

/**
 * @brief True iff the pattern P occurs in H as a subgraph (induced == false) or as an induced
 *        subgraph (induced == true).
 * @param H        host graph (n <= 64)
 * @param P        pattern graph (k <= 10 vertices recommended)
 * @param induced  require non-edges of P to map to non-edges of H
 */
inline bool contains_subgraph(const Graph& H, const Graph& P, bool induced) {
    const int k = P.n, n = H.n;
    if (k == 0) return true;
    if (k > n) return false;
    // order the pattern vertices so that each one (after the first of its component) has an
    // already-placed neighbour; more constrained vertices go first
    std::vector<int> order;
    u64 placed = 0;
    while ((int)order.size() < k) {
        int best = -1, bestc = -1;
        for (int v = 0; v < k; ++v) {
            if (placed >> v & 1) continue;
            const int c = std::popcount(P.adj[v] & placed) * 64 + P.degree(v);
            if (c > bestc) { bestc = c; best = v; }
        }
        order.push_back(best);
        placed |= u64{1} << best;
    }
    std::vector<int> phi(k, -1);
    const u64 allH = low_mask(n);
    // recursive lambda via explicit stack of candidate masks
    std::vector<u64> cand(k + 1, 0);
    auto candidates = [&](int depth) {
        const int v = order[depth];
        u64 c = allH;
        u64 used = 0;
        for (int d = 0; d < depth; ++d) used |= u64{1} << phi[order[d]];
        c &= ~used;
        for (int d = 0; d < depth; ++d) {
            const int u = order[d];
            if (P.has_edge(u, v)) c &= H.adj[phi[u]];
            else if (induced) c &= ~H.adj[phi[u]];
        }
        // degree pruning: an image must have at least the pattern degree
        u64 ok = 0;
        for (u64 x = c; x; x &= x - 1) {
            const int w = std::countr_zero(x);
            if (H.degree(w) >= P.degree(v) && (!induced || n - 1 - H.degree(w) >= k - 1 - P.degree(v))) ok |= u64{1} << w;
        }
        return ok;
    };
    int depth = 0;
    cand[0] = candidates(0);
    while (depth >= 0) {
        if (!cand[depth]) {
            --depth;
            continue;
        }
        const int w = std::countr_zero(cand[depth]);
        cand[depth] &= cand[depth] - 1;
        phi[order[depth]] = w;
        if (depth + 1 == k) return true;
        ++depth;
        cand[depth] = candidates(depth);
    }
    return false;
}

/**
 * @brief Named small patterns for the sampler, or a graph6 string.
 *
 * Names: K2, K3 (triangle), K4, K5, P3, P4, P5, C4, C5, C6, claw (K1,3), star4 (K1,4), paw, diamond
 * (K4-e), 2K2, bull, house, K23 (K2,3), chair (fork), cochair, C4+K1, K4+K1.
 * @param name  pattern name (case-sensitive) or graph6 string
 * @param P     out: pattern graph
 * @return      false if the name is unknown and not valid graph6
 */
inline bool pattern_graph(const std::string& name, Graph& P) {
    static const std::map<std::string, std::pair<int, std::vector<std::pair<int, int>>>> named = {
        {"K2", {2, {{0, 1}}}},
        {"K3", {3, {{0, 1}, {1, 2}, {0, 2}}}},
        {"triangle", {3, {{0, 1}, {1, 2}, {0, 2}}}},
        {"K4", {4, {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}}},
        {"K5", {5, {{0, 1}, {0, 2}, {0, 3}, {0, 4}, {1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4}}}},
        {"P3", {3, {{0, 1}, {1, 2}}}},
        {"P4", {4, {{0, 1}, {1, 2}, {2, 3}}}},
        {"P5", {5, {{0, 1}, {1, 2}, {2, 3}, {3, 4}}}},
        {"C4", {4, {{0, 1}, {1, 2}, {2, 3}, {3, 0}}}},
        {"C5", {5, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}}}},
        {"C6", {6, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}}}},
        {"claw", {4, {{0, 1}, {0, 2}, {0, 3}}}},
        {"star4", {5, {{0, 1}, {0, 2}, {0, 3}, {0, 4}}}},
        {"paw", {4, {{0, 1}, {1, 2}, {0, 2}, {2, 3}}}},
        {"diamond", {4, {{0, 1}, {0, 2}, {1, 2}, {1, 3}, {2, 3}}}},
        {"2K2", {4, {{0, 1}, {2, 3}}}},
        {"bull", {5, {{0, 1}, {1, 2}, {0, 2}, {0, 3}, {1, 4}}}},
        {"house", {5, {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}, {1, 4}}}},
        {"K23", {5, {{0, 2}, {0, 3}, {0, 4}, {1, 2}, {1, 3}, {1, 4}}}},
        {"chair", {5, {{0, 1}, {1, 2}, {1, 3}, {3, 4}}}},
        {"cochair", {5, {{0, 2}, {0, 3}, {0, 4}, {1, 4}, {2, 3}, {2, 4}}}},
        {"C4+K1", {5, {{0, 1}, {1, 2}, {2, 3}, {3, 0}}}},
        {"K4+K1", {5, {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}}},
    };
    auto it = named.find(name);
    if (it != named.end()) {
        P = Graph{};
        P.n = it->second.first;
        for (auto [a, b] : it->second.second) { P.adj[a] |= u64{1} << b; P.adj[b] |= u64{1} << a; }
        return true;
    }
    return parse_graph6(name, P);
}

}  // namespace rgraph
