/**
 * @file   rgraph/labelings.hpp
 * @brief  All Riordan labelings of a given (unlabeled) graph, by exact backtracking search.
 *
 * For a graph G on n <= 64 vertices, a *Riordan labeling* is an ordering v_1, ..., v_n of its
 * vertices such that the relabeled graph equals some G_n(g, f). The search fills positions
 * 0, 1, 2, ... in order. Placing position m >= 1 fixes the coefficient g_(m-1), which is the
 * adjacency to position 0. Once g != 0 with lowest term t^r, it also fixes at most one new
 * coefficient f_(m-1-r). Every other adjacency of the new vertex to earlier positions is then
 * forced: it is [t^(m-1)] g f^j to position j. So at most four adjacency patterns are admissible
 * at each step, and only vertices with a matching pattern are tried. Twin vertices (identical
 * neighbourhoods apart from each other) give the same labeled graph, so only one per twin class is
 * tried at each position.
 *
 * Output: the set C(G) of distinct labeled Riordan graphs isomorphic to G, each given by its
 * canonical parameters (g mod t^(n-1), f mod t^(n-1-r)). In particular:
 *  - G is Riordan  <=>  C(G) is non-empty;
 *  - u(n) = sum over all labeled Riordan graphs H of 1/|C(H)|.
 *
 * Cost: exponential in the worst case; a node limit stops the search (the result is then flagged
 * as incomplete). Typical random Riordan graphs with n <= 24 finish in well under a second.
 */
#pragma once
#include "riordan_graph.hpp"

#include <bit>
#include <cstring>
#include <set>
#include <utility>

namespace rgraph {

/** Result of a labeling search. */
struct LabelingResult {
    std::set<std::pair<u64, u64>> labelings;  ///< canonical (g, f) of every Riordan labeling found
    unsigned long long nodes = 0;             ///< search nodes visited
    bool complete = true;                     ///< false if the node limit stopped the search
};

/**
 * @brief Enumerate all Riordan labelings of G (up to graph automorphisms and twin swaps).
 * @param G           input graph (vertex numbering irrelevant), 1 <= n <= 64
 * @param node_limit  maximum number of search nodes before giving up
 * @param first_only  stop at the first labeling found (recognition only)
 * @return            the distinct labeled Riordan graphs isomorphic to G, with search statistics
 */
inline LabelingResult riordan_labelings(const Graph& G, unsigned long long node_limit = 100000000ULL,
                                        bool first_only = false) {
    struct Search {
        int n;
        u64 adjm[64];
        int twin[64];
        LabelingResult res;
        unsigned long long limit;
        bool first_only, stop = false;
        u64 posmask[64];  // vertex -> positions of its placed neighbours
        u64 used = 0;

        static u64 pmul(u64 a, u64 b, int bits) {
            u64 r = 0;
            for (; b; b &= b - 1) r ^= a << std::countr_zero(b);
            return bits >= 64 ? r : (r & ((u64{1} << bits) - 1));
        }
        // forced adjacency of position m to positions 1..m-1, for known g, f
        static u64 predicted(u64 g, u64 f, int m) {
            u64 pat = 0, col = g;
            for (int j = 1; j < m; ++j) { col = pmul(col, f, m); pat |= ((col >> (m - 1)) & 1) << j; }
            return pat;
        }
        int pos_of[64];  // position -> vertex (valid for placed positions)

        // Capacity pruning. Column j of the Riordan matrix (the later neighbours of position j)
        // has its first non-zero coefficient at index >= r + j*v, where r = val g and v = val f.
        // Once g != 0, r is known exactly, and v >= v_lb, with v_lb = val of the known part of f,
        // or kf+1 if f_1..f_kf are all zero. Position j can therefore have later neighbours only at
        // positions p with p-1 >= max(r + j*v_lb, m-1), p <= n-1. Check this bound for the new
        // vertex w at position m and for every placed position j < m.
        bool later_ok(int m, int w, u64 g, u64 f, int kf) const {
            if (g == 0) return true;
            const int r = std::countr_zero(g);
            const u64 fk = f & low_mask(kf + 1);
            const int vlb = fk ? std::countr_zero(fk) : kf + 1;
            const int N = n - 1;
            const u64 placed_after = used | (u64{1} << w);
            auto cap = [&](int j, int first_free_pos) {  // max later neighbours of position j at positions >= first_free_pos
                const long long lo = std::max<long long>((long long)r + (long long)j * vlb, first_free_pos - 1);
                return (int)std::max<long long>(0, N - lo);
            };
            // the new vertex: its unplaced neighbours must fit after position m
            if (std::popcount(adjm[w] & ~placed_after) > cap(m, m + 1)) return false;
            // the placed vertices: their unplaced neighbours must fit after position m
            for (int j = 0; j < m; ++j) {
                const int u = pos_of[j];
                if (std::popcount(adjm[u] & ~placed_after) > cap(j, m + 1)) return false;
            }
            return true;
        }
        void put(int m, int w) { pos_of[m] = w; put_(m, w); }
        void put_(int m, int w) { used |= u64{1} << w; for (u64 x = adjm[w]; x; x &= x - 1) posmask[std::countr_zero(x)] |= u64{1} << m; }
        void unput(int m, int w) { used &= ~(u64{1} << w); for (u64 x = adjm[w]; x; x &= x - 1) posmask[std::countr_zero(x)] &= ~(u64{1} << m); }

        void place(int m, u64 g, u64 f, int kf) {
            if (stop) return;
            if (++res.nodes > limit) { res.complete = false; stop = true; return; }
            if (m == n) {
                u64 gg = g & low_mask(n - 1), ff = 0;
                if (gg) ff = f & low_mask(n - 1 - std::countr_zero(gg)) & ~u64{1};
                res.labelings.insert({gg, ff});
                if (first_only) stop = true;
                return;
            }
            const u64 mid = (m >= 2) ? (low_mask(m) & ~u64{1}) : 0;
            struct Opt { u64 g, f; int kf; u64 pat; } opts[4];
            int nopt = 0;
            for (int b = 0; b < 2; ++b) {
                const u64 g2 = g | (u64(b) << (m - 1));
                if (g2 == 0) { opts[nopt++] = {g2, f, kf, 0}; continue; }
                const int need = m - 1 - std::countr_zero(g2);
                if (need <= kf) { opts[nopt++] = {g2, f, kf, predicted(g2, f, m)}; continue; }
                for (int fb = 0; fb < 2; ++fb) {
                    const u64 f2 = f | (u64(fb) << need);
                    opts[nopt++] = {g2, f2, need, predicted(g2, f2, m)};
                }
            }
            u64 tried = 0;
            for (int w = 0; w < n && !stop; ++w) {
                if ((used >> w & 1) || (tried >> twin[w] & 1)) continue;
                const u64 b = posmask[w] & 1, row = posmask[w] & mid;
                bool any = false;
                for (int k = 0; k < nopt && !stop; ++k) {
                    if (((opts[k].g >> (m - 1)) & 1) != b || opts[k].pat != row) continue;
                    any = true;
                    if (!later_ok(m, w, opts[k].g, opts[k].f, opts[k].kf)) continue;
                    put(m, w);
                    place(m + 1, opts[k].g, opts[k].f, opts[k].kf);
                    unput(m, w);
                }
                if (any) tried |= u64{1} << twin[w];
            }
        }
    } S;
    S.n = G.n;
    S.limit = node_limit;
    S.first_only = first_only;
    std::memset(S.posmask, 0, sizeof S.posmask);
    for (int v = 0; v < G.n; ++v) S.adjm[v] = G.adj[v];
    for (int u = 0; u < G.n; ++u) {
        S.twin[u] = u;
        for (int w = 0; w < u; ++w) {
            const u64 mask = ~((u64{1} << u) | (u64{1} << w));
            if ((S.adjm[u] & mask) == (S.adjm[w] & mask)) { S.twin[u] = S.twin[w]; break; }
        }
    }
    if (G.n == 1) { S.res.labelings.insert({0, 0}); return S.res; }
    u64 tried = 0;
    for (int v0 = 0; v0 < G.n && !S.stop; ++v0) {
        if (tried >> S.twin[v0] & 1) continue;
        tried |= u64{1} << S.twin[v0];
        S.put(0, v0);
        S.place(1, 0, 0, 0);
        S.unput(0, v0);
    }
    return S.res;
}

}  // namespace rgraph
