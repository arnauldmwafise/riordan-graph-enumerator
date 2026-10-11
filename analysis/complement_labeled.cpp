/**
 * @file   analysis/complement_labeled.cpp
 * @brief  complement_labeled — count the labeled Riordan graphs whose complement, with the same
 *         vertex labeling, is again a labeled Riordan graph (the sequence c(n)).
 *
 * Context
 * -------
 * A labeled Riordan graph G_n(g, f) has vertices 1..n; for i > j, ij is an edge iff
 * [t^(i-2)] g f^(j-1) = 1 over Z_2, with f(0) = 0. There are L(n) = (4^(n-1)+2)/3 of them
 * (u(n) up to isomorphism). This program counts those G whose complement (same labels) equals some
 * G_n(g', f'); this is the labeled version of the complement question ("Problem 3"). The unlabeled
 * version is analysis/complement_unlabeled.cpp.
 *
 * Method
 * ------
 * For each n = 1..NMAX, every labeled Riordan graph is generated once from its canonical
 * parameters (rgraph::Enumeration), complemented, and tested with the labeled recognizer
 * rgraph::recognize_labeled (read g and g f from the first two columns, solve for f, rebuild and
 * compare).
 *
 * Usage
 * -----
 *   complement_labeled NMAX [SHOW]
 *
 *   NMAX  largest order n (n = 1..NMAX are processed; NMAX <= 64, in practice <= 15). Required.
 *   SHOW  optional order n at which every complement-Riordan pair is printed (default 0: none).
 *
 * Output (stdout)
 * ---------------
 *   For n = SHOW, before the summary line of that n, one line per graph:
 *     " g=BITS f=BITS  -> g'=BITS f'=BITS"
 *   where each BITS is the n-1 coefficients [t^0], [t^1], ..., [t^(n-2)] as '0'/'1' characters
 *   (lowest degree first), (g, f) are the canonical parameters of G and (g', f') those found for
 *   its complement.
 *   For every n, one summary line:
 *     "n=N labeled L(N) complement-Riordan C"
 *
 * Cost: L(n) generations and recognitions of O(n^2) word operations each, about 0.5 µs per graph
 * on one core; n = 14 (2.2·10^7 graphs) takes about 10 s. Memory is constant.
 *
 * Results (c(n) for n <= 14) are recorded in docs/COMPUTATIONS.md, section 4.2, and checked there
 * against `rglist --labeled --filter complement --count` and an independent brute force.
 */
#include "rgraph/recognize.hpp"
#include "rgraph/canon.hpp"
#include <cstdio>
#include <cstdlib>
using namespace rgraph;

/**
 * @brief Print the low @p bits coefficients of a Z_2 polynomial as '0'/'1', lowest degree first.
 * @param p     coefficient bits (bit k = [t^k])
 * @param bits  number of coefficients to print
 */
static void pr(u64 p, int bits) {
    for (int k = 0; k < bits; ++k) putchar(p >> k & 1 ? '1' : '0');
}

/**
 * @brief Entry point; see the file header for arguments and output.
 * @return 0
 */
int main(int argc, char** argv) {
    int NMAX = atoi(argv[1]), SHOW = argc > 2 ? atoi(argv[2]) : 0;
    for (int n = 1; n <= NMAX; ++n) {
        Enumeration E(n);
        u64 c = 0;
        for (u64 i = 0; i < E.size(); ++i) {
            u64 g, f;
            E.pair(i, g, f);
            Graph G = E.graph(i);
            auto p = recognize_labeled(complement(G));
            if (!p) continue;
            ++c;
            if (n == SHOW) {
                printf(" g=");
                pr(g, n - 1);
                printf(" f=");
                pr(f, n - 1);
                printf("  -> g'=");
                pr(p->g, n - 1);
                printf(" f'=");
                pr(p->f, n - 1);
                puts("");
            }
        }
        printf("n=%d labeled %llu complement-Riordan %llu\n", n, (unsigned long long)E.size(), (unsigned long long)c);
        fflush(stdout);
    }
}
