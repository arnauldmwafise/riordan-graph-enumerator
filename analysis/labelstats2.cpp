/**
 * @file   analysis/labelstats2.cpp
 * @brief  labelstats2 — the role of the reversal involution in u(n): how the labeled Riordan
 *         graphs with f_1 = 1 and the others contribute to the number of isomorphism classes.
 *
 * Context
 * -------
 * A labeled Riordan graph G_n(g, f) has vertices 1..n; for i > j, ij is an edge iff
 * [t^(i-2)] g f^(j-1) = 1 over Z_2, with f(0) = 0. There are L(n) = (4^(n-1)+2)/3 labeled ones
 * and u(n) unlabeled ones; u(n) = sum over labeled H of 1/|C(H)|, where C(H) is the set of
 * labeled Riordan graphs isomorphic to H. Graphs with f'(0) = f_1 = 1 are reversible: renumbering
 * i -> n+1-i gives another labeled Riordan graph (Cheon et al. 2019, Thm 3.14). This program
 * splits the labeled graphs into
 *   R1: g != 0 and f_1 = 1 (reversal acts on R1),     R0: all others (null graph included).
 *
 * Method
 * ------
 * For each n, all labeled Riordan graphs are generated (rgraph::Enumeration); for those in R1 it is
 * checked whether the reversal rgraph::reversed(G) equals G. Each graph is stored as
 * (class fingerprint rgraph::class_hash, in R1?) and the list is sorted to group classes.
 *
 * Usage
 * -----
 *   labelstats2 NMIN NMAX      (n = NMIN..NMAX; both required; in practice NMAX <= 15)
 *
 * Output (stdout), one line per n:
 *   "n=NN L=L u=U u/L=X | R1=A fixed=B R1-in-2-classes=Y | R0=C alone=Z | u-share from R0=P*L R1=Q*L"
 *   n            order, printed with width 2
 *   L, u         L(n) and u(n) (number of distinct fingerprints)
 *   R1, fixed    |R1| and the number of G in R1 with reversed(G) = G
 *   R1-in-2-classes   fraction of R1 lying in classes of size exactly 2 made of two R1 graphs
 *   R0, alone    |R0| and the fraction of R0 graphs that are alone in their class (|C(H)| = 1)
 *   u-share      sum over R0 (resp. R1) of 1/|C(H)|, divided by L; the two add up to u/L
 * All fractions have 4 decimals.
 *
 * Cost: one nauty canonical form per labeled graph (about 3.5 µs; n = 14 about 1.5 min) and about
 * 24 bytes of memory per labeled graph.
 *
 * Results (n = 10..15) are tabulated in docs/COMPUTATIONS.md, section 4.3.
 */
#include "rgraph/canon.hpp"
#include "rgraph/recognize.hpp"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <vector>
using namespace rgraph;

/** @brief A labeled Riordan graph: class fingerprint and t = 1 if in R1, ordered by (h, t). */
struct E {
    Hash128 h;
    unsigned char t;
    /** @brief Lexicographic order on (h, t), so that equal classes are adjacent after sorting. */
    bool operator<(const E& o) const { return h < o.h || (h == o.h && t < o.t); }
};

/**
 * @brief Entry point; see the file header for arguments and output.
 * @return 0
 */
int main(int argc, char** argv) {
    const int A = std::atoi(argv[1]), B = std::atoi(argv[2]);
    for (int n = A; n <= B; ++n) {
        Enumeration En(n);
        std::vector<E> v;
        v.reserve(En.size());
        u64 r1 = 0, fixed = 0;
        for (u64 i = 0; i < En.size(); ++i) {
            u64 g, f;
            En.pair(i, g, f);
            const Graph G = En.graph(i);
            const bool inR1 = (g != 0) && ((f >> 1) & 1);
            if (inR1) {
                ++r1;
                if (reversed(G) == G) ++fixed;
            }
            v.push_back({class_hash(G), (unsigned char)inR1});
        }
        std::sort(v.begin(), v.end());
        // contribution of each labeled graph to u(n): 1/(class size); sum by type
        double contrib[2] = {0, 0};
        u64 alone0 = 0, n0 = 0, pair1 = 0, classes = 0, pure1pairs = 0;
        for (size_t i = 0; i < v.size();) {
            size_t j = i;
            int c[2] = {0, 0};
            while (j < v.size() && v[j].h == v[i].h) {
                c[v[j].t]++;
                ++j;
            }
            const int m = int(j - i);
            ++classes;
            contrib[0] += double(c[0]) / m;
            contrib[1] += double(c[1]) / m;
            n0 += c[0];
            if (m == 1 && c[0]) ++alone0;
            if (m == 2 && c[1] == 2) {
                pair1 += 2;
                ++pure1pairs;
            }
            i = j;
        }
        const u64 L = En.size();
        std::printf("n=%2d L=%llu u=%llu u/L=%.4f | R1=%llu fixed=%llu R1-in-2-classes=%.4f | R0=%llu alone=%.4f | u-share from R0=%.4f*L R1=%.4f*L\n",
                    n, (unsigned long long)L, (unsigned long long)classes, double(classes) / L, (unsigned long long)r1, (unsigned long long)fixed,
                    double(pair1) / r1, (unsigned long long)n0, double(alone0) / n0, contrib[0] / L, contrib[1] / L);
        std::fflush(stdout);
    }
}
