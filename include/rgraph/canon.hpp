/**
 * @file   rgraph/canon.hpp
 * @brief  Isomorphism-class fingerprints of graphs (n <= 64) via nauty canonical forms, used by
 *         method A (tools/rgcount.cpp).
 *
 * Method
 * ------
 * The canonical form is computed by nauty's densenauty() relative to the ordered partition of the
 * vertices by degree (cells in increasing degree). Since that partition is isomorphism-invariant,
 * isomorphic graphs get the same canonical graph. The class fingerprint is a 128-bit hash of the
 * order n and the canonical adjacency rows (which also determine the degree sequence). Distinct
 * classes could in principle collide; method B (verify/tree_count.cpp) uses a different initial
 * partition and a different hash, so a collision would show up as a disagreement between A and B.
 *
 * Build requirements: nauty must be compiled with MAXN=WORDSIZE=64 (one setword per row) and with
 * thread-local storage (USE_TLS), so that canonical_graph() may be called from several threads.
 *
 * Conventions: rgraph::Graph rows with bit w = vertex w (0-based); nauty numbers bits from the
 * most significant end (vertex v <-> bit[v] = bit 63 - v); the conversion is done here.
 */
#pragma once
#include "riordan_graph.hpp"
#include "nauty.h"

#include <algorithm>

namespace rgraph {

/** @brief 128-bit hash value; totally ordered so that vectors of hashes can be sorted and deduplicated. */
struct Hash128 {
    u64 a = 0, b = 0;  ///< the two 64-bit halves
    /** @brief Lexicographic order on (a, b). */
    bool operator<(const Hash128& o) const { return a != o.a ? a < o.a : b < o.b; }
    /** @brief Equality of both halves. */
    bool operator==(const Hash128& o) const { return a == o.a && b == o.b; }
};

/**
 * @brief splitmix64 finaliser: a bijective 64-bit mixing function.
 * @param x  input word
 * @return   mixed word
 */
inline u64 mix64(u64 x) {  // splitmix64 finaliser
    x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27; x *= 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

/**
 * @brief Two independently seeded chains -> 128-bit hash of a word sequence.
 *
 * Feed words with add() in order, then read the hash with get(). The position counter i makes the
 * hash order-sensitive and length-sensitive.
 */
struct Hasher {
    /// a, b: states of the two chains; i: number of words added so far
    u64 a = 0x243F6A8885A308D3ULL, b = 0x13198A2E03707344ULL, i = 0;
    /** @brief Append the word w to the hashed sequence. */
    void add(u64 w) {
        ++i;
        a = mix64(a ^ mix64(w + i * 0x9E3779B97F4A7C15ULL));
        b = mix64(b + mix64(w ^ (i * 0xD6E8FEB86659FD93ULL + 0x632BE59BD9B4E019ULL)));
    }
    /** @brief The 128-bit hash of the words added so far (does not modify the state). */
    Hash128 get() const { return {mix64(a ^ i), mix64(b + i)}; }
};

/**
 * @brief Hash of n and the sorted degree sequence: an isomorphism invariant (used to split work
 *        into passes in method A).
 * @param G  graph (n <= 64)
 * @return   64-bit hash; isomorphic graphs give equal values
 */
inline u64 degree_invariant(const Graph& G) {
    int d[64];
    for (int v = 0; v < G.n; ++v) d[v] = G.degree(v);
    std::sort(d, d + G.n);
    u64 h = 0x452821E638D01377ULL ^ static_cast<u64>(G.n);
    for (int v = 0; v < G.n; ++v) h = mix64(h ^ (static_cast<u64>(d[v]) + 0x9E3779B97F4A7C15ULL * (v + 1)));
    return h;
}

/**
 * @brief Canonical relabeling of G: the returned graph is the same for all graphs isomorphic to G.
 *
 * Calls nauty (densenauty, getcanon) with the initial ordered partition of the vertices by
 * increasing degree. Thread-safe given a USE_TLS build of nauty (the option block is thread-local).
 * @param G  graph, 0 <= n <= 64
 * @return   canonical graph C (C.n = G.n), in rgraph bit convention
 */
inline Graph canonical_graph(const Graph& G) {
    static thread_local DEFAULTOPTIONS_GRAPH(opts);
    const int n = G.n;
    graph g[MAXN], cg[MAXN];
    int lab[MAXN], ptn[MAXN], orbits[MAXN];
    int deg[MAXN];
    for (int v = 0; v < n; ++v) { deg[v] = G.degree(v); lab[v] = v; }
    std::sort(lab, lab + n, [&](int x, int y) { return deg[x] != deg[y] ? deg[x] < deg[y] : x < y; });
    for (int k = 0; k < n; ++k) ptn[k] = (k + 1 < n && deg[lab[k + 1]] == deg[lab[k]]) ? 1 : 0;
    // nauty numbers bits from the most significant end: vertex v <-> bit (63 - v).
    for (int v = 0; v < n; ++v) {
        setword row = 0;
        for (u64 x = G.adj[v]; x; x &= x - 1) row |= bit[std::countr_zero(x)];
        g[v] = row;
    }
    optionblk o = opts;
    o.getcanon = TRUE;
    o.defaultptn = FALSE;
    statsblk stats;
    densenauty(g, lab, ptn, orbits, &o, &stats, 1, n, cg);
    Graph C;
    C.n = n;
    for (int v = 0; v < n; ++v) {
        u64 row = 0;
        for (setword x = cg[v]; x; ) {
            const int w = FIRSTBITNZ(x);
            x ^= bit[w];
            row |= u64{1} << w;
        }
        C.adj[v] = row;
    }
    return C;
}

/**
 * @brief 128-bit fingerprint of the isomorphism class of G: Hasher over (n, canonical rows).
 * @param G  graph, 0 <= n <= 64
 * @return   equal for isomorphic graphs; distinct classes differ except for hash collisions
 */
inline Hash128 class_hash(const Graph& G) {
    const Graph C = canonical_graph(G);
    Hasher h;
    h.add(static_cast<u64>(G.n));
    for (int v = 0; v < C.n; ++v) h.add(C.adj[v]);
    return h.get();
}

/**
 * @brief Complement of G (no loops): same vertices, uv an edge iff it is not an edge of G.
 * @param G  graph, n <= 64
 * @return   the complement graph
 */
inline Graph complement(const Graph& G) {
    Graph H;
    H.n = G.n;
    const u64 all = low_mask(G.n);
    for (int v = 0; v < G.n; ++v) H.adj[v] = ~G.adj[v] & all & ~(u64{1} << v);
    return H;
}

}  // namespace rgraph
