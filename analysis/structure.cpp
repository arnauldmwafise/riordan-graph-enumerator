/**
 * @file   analysis/structure.cpp
 * @brief  structure — minimal non-Riordan graphs, failures of heredity and complements, from the
 *         complete lists of Riordan / non-Riordan graphs produced by method C.
 *
 * Context
 * -------
 * A labeled Riordan graph G_n(g, f) has vertices 1..n; for i > j, ij is an edge iff
 * [t^(i-2)] g f^(j-1) = 1 over Z_2, with f(0) = 0. There are L(n) = (4^(n-1)+2)/3 labeled ones;
 * a graph is (unlabeled) Riordan if some ordering of its vertices is a labeled Riordan graph, and
 * u(n) counts these up to isomorphism. Method C (scripts/method_c.sh) tests every graph of order n
 * (nauty geng) and splits them into two graph6 lists:
 *
 *   geng -q N | recognize_all --riordan DIR/R_N.g6 --nonriordan DIR/NR_N.g6
 *
 * For each n this program reports ("Problems 1, 3 and 5"):
 *  - minimal non-Riordan graphs: non-Riordan G such that every G - v is Riordan;
 *  - Riordan graphs G having a vertex v with G - v non-Riordan (failures of heredity);
 *  - Riordan graphs whose complement is Riordan, and self-complementary Riordan graphs.
 *
 * Method
 * ------
 * Isomorphism classes are compared by the 128-bit nauty fingerprints rgraph::class_hash of
 * include/rgraph/canon.hpp. The Riordan list of order n is loaded into a fingerprint set; the
 * non-Riordan list is streamed. A vertex-deleted subgraph G - v is Riordan iff its fingerprint is in
 * the set of order n-1 (kept from the previous iteration), so the orders must be processed from 1.
 *
 * Usage
 * -----
 *   structure DIR NMAX
 *
 *   DIR   directory with R_n.g6 and NR_n.g6 for n = 1..NMAX (uncompressed; a missing file is
 *         treated as empty, so e.g. data/lists/NR_10.g6.gz must be decompressed first)
 *   NMAX  largest order (n <= 62, graph6 short form)
 * Fewer than two arguments: "usage: structure DIR NMAX" on stderr, exit status 2.
 *
 * Inputs: DIR/R_n.g6, DIR/NR_n.g6 (one graph6 string per line).
 * Outputs:
 *  - file DIR/minimal_n.g6 for every n (overwritten; the minimal non-Riordan graphs, graph6);
 *  - stdout, per n:
 *      "n=N riordan=R nonriordan=NR minimal_nonriordan=M (connected MC) riordan_with_nonriordan_deletion=H complement_riordan=CR self_complementary=SC"
 *    and, if H > 0, the first such graph:
 *      "  example: G6 minus vertex V = G6"     (V 0-based, G6 = graph6 strings)
 *
 * Cost: at most n canonical forms per graph in NR_n.g6 (vertex deletions, stopping at the first
 * non-Riordan one) and up to n + 3 per graph in R_n.g6. The non-Riordan list dominates: at most
 * 1.2·10^8 canonical forms for n = 10 (minutes); n <= 9 takes about a second.
 * Memory: the Riordan lists of orders n and n-1.
 *
 * Results (n <= 10, run on data/lists) are in docs/COMPUTATIONS.md, sections 4.1 and 4.2; the
 * minimal graphs are in data/lists/minimal_<n>.g6.
 */
#include "rgraph/canon.hpp"

#include <cstdio>
#include <set>
#include <string>
#include <vector>

using namespace rgraph;

/**
 * @brief Decode a short-form graph6 string (n <= 62).
 * @param s  graph6 text without newline
 * @param G  output graph
 * @return   false if the first character is not a valid order byte (the body is not validated)
 */
static bool parse_g6(const std::string& s, Graph& G) {
    if (s.empty() || s[0] < 63 || s[0] > 126) return false;
    G = Graph{};
    G.n = s[0] - 63;
    std::size_t p = 1;
    int k = 6, x = 0;
    for (int j = 1; j < G.n; ++j)
        for (int i = 0; i < j; ++i) {
            if (k == 6) { x = s[p++] - 63; k = 0; }
            if ((x >> (5 - k)) & 1) { G.adj[i] |= u64{1} << j; G.adj[j] |= u64{1} << i; }
            ++k;
        }
    return true;
}

/**
 * @brief graph6 encoding of a graph (short form, n <= 62).
 * @param G  graph
 * @return   graph6 string without newline
 */
static std::string to_g6(const Graph& G) {
    std::string s(1, char(63 + G.n));
    int k = 0, x = 0;
    for (int j = 1; j < G.n; ++j)
        for (int i = 0; i < j; ++i) {
            x = (x << 1) | int(G.has_edge(i, j));
            if (++k == 6) { s += char(63 + x); k = 0; x = 0; }
        }
    if (k) s += char(63 + (x << (6 - k)));
    return s;
}

/**
 * @brief Read all graphs of a graph6 file (one per line; CR/LF stripped; invalid lines skipped).
 * @param path  file name
 * @return      the graphs; empty if the file cannot be opened
 */
static std::vector<Graph> load(const std::string& path) {
    std::vector<Graph> v;
    FILE* f = std::fopen(path.c_str(), "r");
    if (!f) return v;
    char buf[256];
    while (std::fgets(buf, sizeof buf, f)) {
        std::string s(buf);
        while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
        Graph G;
        if (parse_g6(s, G)) v.push_back(G);
    }
    std::fclose(f);
    return v;
}

/**
 * @brief The vertex-deleted subgraph G - v (remaining vertices keep their relative order).
 * @param G  graph
 * @param v  vertex to delete (0-based)
 * @return   graph on G.n - 1 vertices
 */
static Graph delete_vertex(const Graph& G, int v) {
    Graph H;
    H.n = G.n - 1;
    for (int a = 0, i = 0; a < G.n; ++a) {
        if (a == v) continue;
        u64 row = 0;
        for (int b = 0, j = 0; b < G.n; ++b) {
            if (b == v) continue;
            if (G.has_edge(a, b)) row |= u64{1} << j;
            ++j;
        }
        H.adj[i++] = row;
    }
    return H;
}

/**
 * @brief Entry point; see the file header for arguments, inputs and outputs.
 * @return 0 on success, 2 on a usage error
 */
int main(int argc, char** argv) {
    if (argc < 3) { std::fprintf(stderr, "usage: structure DIR NMAX\n"); return 2; }
    const std::string dir = argv[1];
    const int NMAX = std::atoi(argv[2]);
    std::set<Hash128> riordan_prev;
    for (int n = 1; n <= NMAX; ++n) {
        const auto R = load(dir + "/R_" + std::to_string(n) + ".g6");
        // the non-Riordan list can be large: it is streamed, not loaded
        FILE* nrf = std::fopen((dir + "/NR_" + std::to_string(n) + ".g6").c_str(), "r");
        u64 nr_count = 0;
        std::set<Hash128> riordan;
        for (const auto& G : R) riordan.insert(class_hash(G));
        u64 minimal = 0, nonhered = 0, comp = 0, selfc = 0, minimal_conn = 0;
        std::string example;
        FILE* out = std::fopen((dir + "/minimal_" + std::to_string(n) + ".g6").c_str(), "w");
        if (n >= 2) {
            char line[256];
            while (nrf && std::fgets(line, sizeof line, nrf)) {
                std::string sline(line);
                while (!sline.empty() && (sline.back() == '\n' || sline.back() == '\r')) sline.pop_back();
                Graph G;
                if (!parse_g6(sline, G)) continue;
                ++nr_count;
                bool all_r = true;
                for (int v = 0; v < n && all_r; ++v) all_r = riordan_prev.count(class_hash(delete_vertex(G, v)));
                if (all_r) {
                    ++minimal;
                    if (out) std::fprintf(out, "%s\n", to_g6(G).c_str());
                    // connected?
                    u64 seen = 1, fr = 1;
                    while (fr) {
                        u64 nx = 0;
                        for (u64 x = fr; x; x &= x - 1) nx |= G.adj[std::countr_zero(x)];
                        fr = nx & ~seen;
                        seen |= nx;
                    }
                    if (std::popcount(seen) == n) ++minimal_conn;
                }
            }
            for (const auto& G : R)
                for (int v = 0; v < n; ++v)
                    if (!riordan_prev.count(class_hash(delete_vertex(G, v)))) {
                        if (!nonhered++) example = to_g6(G) + " minus vertex " + std::to_string(v) + " = " + to_g6(delete_vertex(G, v));
                        break;
                    }
        }
        if (out) std::fclose(out);
        if (nrf) std::fclose(nrf);
        for (const auto& G : R) {
            const Hash128 hc = class_hash(complement(G));
            if (riordan.count(hc)) ++comp;
            if (hc == class_hash(G)) ++selfc;
        }
        std::printf("n=%d riordan=%zu nonriordan=%llu minimal_nonriordan=%llu (connected %llu) "
                    "riordan_with_nonriordan_deletion=%llu complement_riordan=%llu self_complementary=%llu\n",
                    n, R.size(), (unsigned long long)nr_count, (unsigned long long)minimal, (unsigned long long)minimal_conn,
                    (unsigned long long)nonhered, (unsigned long long)comp, (unsigned long long)selfc);
        if (!example.empty()) std::printf("  example: %s\n", example.c_str());
        std::fflush(stdout);
        riordan_prev.swap(riordan);
    }
    return 0;
}
