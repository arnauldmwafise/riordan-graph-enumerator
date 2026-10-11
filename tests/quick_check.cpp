/**
 * @file   tests/quick_check.cpp
 * @brief  quick_check — self-test of the labeled enumeration: the L(n) generated graphs are
 *         pairwise distinct, and the number of isomorphism classes u(n) for n <= 11.
 *
 * Context
 * -------
 * A labeled Riordan graph G_n(g, f) has vertices 1..n; for i > j, ij is an edge iff
 * [t^(i-2)] g f^(j-1) = 1 over Z_2, with f(0) = 0. rgraph::Enumeration numbers the
 * L(n) = (4^(n-1)+2)/3 labeled Riordan graphs by their canonical parameters 0..L(n)-1.
 *
 * Method
 * ------
 * For n = 1..11, every enumerated graph is inserted into a set of adjacency matrices (labeled
 * distinctness) and its class fingerprint rgraph::class_hash into a set of classes (u(n)).
 * The check passes when "labeled distinct" equals the formula value in every line; this is
 * confirmed by reading the output (the program itself always exits with status 0).
 *
 * Usage:  quick_check            (no arguments; run by `make check`)
 *
 * Output (stdout), one line per n:
 *   "n=NN  labeled distinct=D (formula L)  unlabeled=U"      (n with width 2)
 * Expected unlabeled values: 1 2 4 11 30 109 406 1864 8169 36897 158351.
 *
 * Cost: about 2 s (n = 11 dominates, 3.5·10^5 canonical forms).
 *
 * See docs/COMPUTATIONS.md, section 1 (labeled enumeration) and section 5 (cross-checks).
 */
#include "rgraph/canon.hpp"
#include <cstdio>
#include <set>
#include <vector>
using namespace rgraph;

/**
 * @brief Run the check for n = 1..11 and print one line per order.
 * @return 0
 */
int main() {
    for (int n = 1; n <= 11; ++n) {
        Enumeration E(n);
        std::set<std::vector<u64>> labeled;
        std::set<Hash128> classes;
        for (u64 i = 0; i < E.size(); ++i) {
            Graph G = E.graph(i);
            labeled.insert(std::vector<u64>(G.adj, G.adj + n));
            classes.insert(class_hash(G));
        }
        std::printf("n=%2d  labeled distinct=%llu (formula %llu)  unlabeled=%zu\n", n,
                    (unsigned long long)labeled.size(), (unsigned long long)E.size(), classes.size());
    }
}
