/**
 * @file   tests/recognize_test.cpp
 * @brief  recognize_test — exhaustive test of the labeled recognizer rgraph::recognize_labeled
 *         against the enumeration rgraph::Enumeration.
 *
 * Context
 * -------
 * A labeled Riordan graph G_n(g, f) has vertices 1..n; for i > j, ij is an edge iff
 * [t^(i-2)] g f^(j-1) = 1 over Z_2, with f(0) = 0. There are L(n) = (4^(n-1)+2)/3 of them.
 * The recognizer (include/rgraph/recognize.hpp) decides whether a labeled graph equals some
 * G_n(g, f); the enumeration generates all of them. The two are written independently.
 *
 * Method
 * ------
 *  1. n = 1..7: every generated graph must be recognized, and G_n(g, f) rebuilt from the returned
 *     (g, f) must equal it. Then all 2^(n(n-1)/2) labeled graphs on n vertices are tested: the
 *     recognizer must accept exactly the generated ones, so the accepted count equals L(n).
 *  2. n = 8..20: 200000 pseudo-random generated graphs (64-bit LCG, seed 12345 + n) must be
 *     recognized. Each is then modified by flipping one random vertex pair; if the result is still
 *     accepted, rebuilding from the returned (g, f) must reproduce it.
 *
 * Usage:  recognize_test          (no arguments; run by `make check`)
 *
 * Output (stdout):
 *   "n=N accepted A of 2^(n(n-1)/2), L=L(N)"   per n <= 7 (A must equal L(N))
 *   "n=N ok (F flips still Riordan)"            per n = 8..20 (F = accepted flipped graphs)
 * On the first failure one of "FAIL gen n=N", "FAIL n=N bits=B a=A b=B", "FAIL n=N" or "FAIL2"
 * is printed and the exit status is 1; otherwise 0.
 *
 * Cost: about 3 s (dominated by the 2^21 labeled graphs of order 7).
 *
 * See docs/COMPUTATIONS.md, sections 1, 3 and 5.
 */
#include "rgraph/recognize.hpp"
#include <cstdio>
#include <set>
#include <vector>
using namespace rgraph;

/**
 * @brief Run both test stages; see the file header.
 * @return 0 if all checks pass, 1 at the first failure
 */
int main() {
    // exhaustive: n <= 7, every labeled graph
    for (int n = 1; n <= 7; ++n) {
        std::set<std::vector<u64>> gen;
        Enumeration E(n);
        for (u64 i = 0; i < E.size(); ++i) {
            Graph G = E.graph(i);
            gen.insert(std::vector<u64>(G.adj, G.adj + n));
            auto p = recognize_labeled(G);
            if (!p || !(riordan_graph(n, p->g, p->f) == G)) {
                printf("FAIL gen n=%d\n", n);
                return 1;
            }
        }
        int m = n * (n - 1) / 2;
        u64 acc = 0;
        for (u64 bits = 0; bits < (u64{1} << m); ++bits) {
            Graph G;
            G.n = n;
            int k = 0;
            for (int j = 1; j < n; ++j)
                for (int i = 0; i < j; ++i, ++k)
                    if (bits >> k & 1) {
                        G.adj[i] |= u64{1} << j;
                        G.adj[j] |= u64{1} << i;
                    }
            bool a = is_riordan_labeled(G), b = gen.count(std::vector<u64>(G.adj, G.adj + n));
            if (a != b) {
                printf("FAIL n=%d bits=%llu a=%d b=%d\n", n, (unsigned long long)bits, a, b);
                return 1;
            }
            acc += a;
        }
        printf("n=%d accepted %llu of %llu, L=%llu\n", n, (unsigned long long)acc, (unsigned long long)(u64{1} << m),
               (unsigned long long)labeled_count(n));
    }
    // larger n: all generated recognized, and random single-edge flips rejected or re-recognized consistently
    for (int n = 8; n <= 20; ++n) {
        Enumeration E(n);
        u64 s = 12345 + n, bad = 0;
        for (int t = 0; t < 200000; ++t) {
            s = s * 6364136223846793005ULL + 1442695040888963407ULL;
            u64 i = (s >> 11) % E.size();
            Graph G = E.graph(i);
            if (!is_riordan_labeled(G)) {
                printf("FAIL n=%d\n", n);
                return 1;
            }
            int a = (s >> 3) % n, b = (s >> 20) % n;
            if (a == b) continue;
            G.adj[a] ^= u64{1} << b;
            G.adj[b] ^= u64{1} << a;
            auto p = recognize_labeled(G);
            if (p) {
                if (!(riordan_graph(n, p->g, p->f) == G)) {
                    printf("FAIL2\n");
                    return 1;
                }
                ++bad;
            }
        }
        printf("n=%d ok (%llu flips still Riordan)\n", n, (unsigned long long)bad);
    }
}
