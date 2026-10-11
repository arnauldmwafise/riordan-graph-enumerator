/**
 * @file   analysis/complement_unlabeled.cpp
 * @brief  complement_unlabeled — count unlabeled Riordan graphs whose complement is also Riordan
 *         (in some labeling), and self-complementary Riordan graphs.
 *
 * Context
 * -------
 * A labeled Riordan graph G_n(g, f) has vertices 1..n; for i > j, ij is an edge iff
 * [t^(i-2)] g f^(j-1) = 1 over Z_2, with f(0) = 0. There are L(n) = (4^(n-1)+2)/3 labeled ones
 * and u(n) unlabeled ones (isomorphism classes). This is the unlabeled version of the complement
 * question ("Problem 3"); the labeled version is analysis/complement_labeled.cpp.
 *
 * Method
 * ------
 * For each n, all L(n) labeled Riordan graphs G are generated (rgraph::Enumeration). For each, the
 * 128-bit isomorphism-class fingerprint (rgraph::class_hash: nauty canonical form, then hash) of G
 * and of its complement are recorded. With S = set of class fingerprints of Riordan graphs and
 * C = set of fingerprints of their complements:
 *   riordan             = |S|        (= u(n))
 *   complement_riordan  = |S ∩ C|    (classes whose complement is also a Riordan class)
 *   self_complementary  = number of distinct classes with hash(G) = hash(complement G)
 * The fingerprint vectors are sorted and deduplicated whenever S exceeds 2^24 entries, to bound
 * memory.
 *
 * Usage
 * -----
 *   complement_unlabeled NMAX      (n = 1..NMAX; required; NMAX <= 64, in practice <= 15)
 *
 * Output (stdout), one line per n:
 *   "n=N riordan=U complement_riordan=UC self_complementary=SC"
 *
 * Cost: two nauty canonical forms per labeled graph, about 7 µs per graph per core (n = 14: a few
 * minutes); memory about 32 bytes per labeled graph before deduplication.
 *
 * Results (n <= 14) are recorded in docs/COMPUTATIONS.md, section 4.2; for n <= 10 they agree with
 * the method C lists processed by analysis/structure.cpp.
 */
#include "rgraph/canon.hpp"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <vector>
using namespace rgraph;

/**
 * @brief Sort a vector of fingerprints and remove duplicates, in place.
 * @param v  fingerprints
 */
static void uniq(std::vector<Hash128>& v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
}

/**
 * @brief Entry point; see the file header for arguments and output.
 * @return 0
 */
int main(int argc, char** argv) {
    const int NMAX = std::atoi(argv[1]);
    for (int n = 1; n <= NMAX; ++n) {
        Enumeration E(n);
        std::vector<Hash128> S, C, SC;
        for (u64 i = 0; i < E.size(); ++i) {
            const Graph G = E.graph(i);
            const Hash128 a = class_hash(G), b = class_hash(complement(G));
            S.push_back(a);
            C.push_back(b);
            if (a == b) SC.push_back(a);
            if (S.size() > (1u << 24)) {
                uniq(S);
                uniq(C);
            }
        }
        uniq(S);
        uniq(C);
        uniq(SC);
        std::vector<Hash128> I;
        std::set_intersection(S.begin(), S.end(), C.begin(), C.end(), std::back_inserter(I));
        std::printf("n=%d riordan=%zu complement_riordan=%zu self_complementary=%zu\n", n, S.size(), I.size(), SC.size());
        std::fflush(stdout);
    }
}
