/**
 * @file   analysis/labelstats.cpp
 * @brief  labelstats — distribution of the number of Riordan labelings per isomorphism class,
 *         split by the type of the parameters (g, f).
 *
 * Context
 * -------
 * A labeled Riordan graph G_n(g, f) has vertices 1..n; for i > j, ij is an edge iff
 * [t^(i-2)] g f^(j-1) = 1 over Z_2, with f(0) = 0. There are L(n) = (4^(n-1)+2)/3 labeled ones,
 * and u(n) isomorphism classes. For a labeled H let C(H) be the set of labeled Riordan graphs
 * isomorphic to H (its "multiplicity" |C(H)|); then u(n) = sum over H of 1/|C(H)|. This program
 * describes how the L(n) labeled graphs are spread over the u(n) classes.
 *
 * Each labeled graph, with canonical parameters (g, f) and r = lowest degree of g, gets a type:
 *   T0: null graph (g = 0);   T1: r >= 1;   T2: r = 0 and f_1 = 0;   T3: r = 0 and f_1 = 1 (proper).
 *
 * Method
 * ------
 * For each n, all labeled Riordan graphs are generated (rgraph::Enumeration) and stored as pairs
 * (class fingerprint rgraph::class_hash, type). Sorting groups each isomorphism class together;
 * one pass over the groups collects the statistics below.
 *
 * Usage
 * -----
 *   labelstats NMIN NMAX      (n = NMIN..NMAX; both required; in practice NMAX <= 15)
 *
 * Output (stdout), per n, six lines (the label texts are exactly as printed):
 *   "n=N L=L(N) u=U L/u=RATIO"
 *   "  multiplicity: m:k m:k ..."        k classes have exactly m labelings (only m <= 8 shown)
 *   "  labeled per type: T0=x(alone a) T1=... T2=... T3=..."
 *                                         x labeled graphs of each type; a = fraction of them that
 *                                         are alone in their class (|C(H)| = 1), 3 decimals
 *   "  in classes of size<=2: T1 p T2 p T3 p"   fraction of each type with |C(H)| <= 2
 *   "  classes containing only type t: T0 k T1 k T2 k T3 k"
 *   "  top signatures (T1,T2,T3 counts -> classes): (a,b,c):k ..."
 *                                         the 10 most frequent per-class type counts
 *                                         (numbers of T1, T2, T3 members), by decreasing number
 *                                         of classes (ties: larger signature first)
 *
 * Cost: one nauty canonical form per labeled graph (about 3.5 µs per graph; n = 14 about 1.5 min)
 * and about 24 bytes of memory per labeled graph.
 *
 * Results are summarised in docs/COMPUTATIONS.md, section 4.3 (together with labelstats2).
 */
#include "rgraph/canon.hpp"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <vector>
using namespace rgraph;

/** @brief A labeled Riordan graph: class fingerprint and type (0..3), ordered by (h, t). */
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
        for (u64 i = 0; i < En.size(); ++i) {
            u64 g, f;
            En.pair(i, g, f);
            int t = g == 0 ? 0 : (g & 1) == 0 ? 1 : (f & 2) ? 3 : 2;
            v.push_back({class_hash(En.graph(i)), (unsigned char)t});
        }
        std::sort(v.begin(), v.end());
        std::map<int, u64> mult;              // multiplicity -> #classes
        std::map<std::vector<int>, u64> sig;  // type signature (counts per type) -> #classes
        u64 classes = 0, lab_in_single[4] = {0}, lab[4] = {0}, lab_le2[4] = {0}, pure[4] = {0};
        for (size_t i = 0; i < v.size();) {
            size_t j = i;
            std::vector<int> c(4, 0);
            while (j < v.size() && v[j].h == v[i].h) {
                c[v[j].t]++;
                ++j;
            }
            int m = int(j - i);
            ++classes;
            ++mult[m];
            ++sig[c];
            int ntypes = 0;
            for (int t = 0; t < 4; ++t) ntypes += c[t] > 0;
            for (int t = 0; t < 4; ++t) {
                lab[t] += c[t];
                if (m == 1) lab_in_single[t] += c[t];
                if (m <= 2) lab_le2[t] += c[t];
                if (ntypes == 1 && c[t]) pure[t]++;
            }
            i = j;
        }
        std::printf("n=%d L=%llu u=%llu L/u=%.4f\n  multiplicity:", n, (unsigned long long)En.size(), (unsigned long long)classes, double(En.size()) / classes);
        for (auto& [m, k] : mult)
            if (m <= 8) std::printf(" %d:%llu", m, (unsigned long long)k);
        std::printf("\n  labeled per type:");
        for (int t = 0; t < 4; ++t)
            std::printf(" T%d=%llu(alone %.3f)", t, (unsigned long long)lab[t], lab[t] ? double(lab_in_single[t]) / lab[t] : 0.0);
        std::printf("\n  in classes of size<=2:");
        for (int t = 1; t < 4; ++t) std::printf(" T%d %.3f", t, double(lab_le2[t]) / lab[t]);
        std::printf("\n  classes containing only type t:");
        for (int t = 0; t < 4; ++t) std::printf(" T%d %llu", t, (unsigned long long)pure[t]);
        std::printf("\n  top signatures (T1,T2,T3 counts -> classes):");
        std::vector<std::pair<u64, std::vector<int>>> s;
        for (auto& [c, k] : sig) s.push_back({k, c});
        std::sort(s.rbegin(), s.rend());
        for (size_t i = 0; i < s.size() && i < 10; ++i)
            std::printf(" (%d,%d,%d):%llu", s[i].second[1], s[i].second[2], s[i].second[3], (unsigned long long)s[i].first);
        std::printf("\n");
        std::fflush(stdout);
    }
}
