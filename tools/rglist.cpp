/**
 * @file   tools/rglist.cpp
 * @brief  rglist — list every labeled or unlabeled Riordan graph of a given order, with its pair of
 *         generating functions, Riordan matrix, adjacency matrix, edge list, graph6 code and drawing.
 *
 * Purpose
 * -------
 * Makes the objects counted by u(n) explicit. For small orders, anyone can inspect each unlabeled
 * Riordan graph, see a Riordan labeling of it, and check the counts by hand.
 *
 * Method
 * ------
 * All L(n) = (4^(n-1)+2)/3 labeled Riordan graphs are generated from their canonical parameters
 * (rgraph::Enumeration). In unlabeled mode they are grouped by the 128-bit nauty fingerprint of
 * their isomorphism class (rgraph/canon.hpp), the same fingerprint as counting method A. The
 * representative of a class is the member G_n(g,f) with the fewest non-zero coefficients in g and f
 * (ties: smallest g, then smallest f). Whether the complement is Riordan is decided by looking up
 * the fingerprint of the complement (unlabeled mode), or by labeled recognition of the complement
 * (labeled mode, rgraph/recognize.hpp).
 *
 * Usage
 * -----
 *   rglist N [--unlabeled | --labeled] [--filter all|complement|selfcomplementary|connected]
 *            [--format text|markdown|csv|json] [--svg DIR] [--count]
 *   rglist N --show G F [--format ...] [--svg DIR]
 *
 *   --unlabeled (default)  one entry per isomorphism class (1 <= N <= 14), with the number of
 *                          labeled Riordan graphs in the class and whether the complement is Riordan.
 *   --labeled              every labeled Riordan graph once, in enumeration order (1 <= N <= 16).
 *   --filter complement    only graphs whose complement is Riordan (unlabeled: in some labeling;
 *                          labeled: with the same labeling).
 *   --filter selfcomplementary   only self-complementary graphs (unlabeled mode).
 *   --filter connected     only connected graphs.
 *   --format               text (default), markdown, csv (one row per graph) or json.
 *   --svg DIR              write a drawing DIR/rg_<N>_<k>.svg for each listed graph.
 *   --count                print only the number of graphs that pass the filter.
 *   --show G F             a single graph G_N(g, f); polynomials are written like "1+t+t^3".
 *
 * Outputs (stdout)
 * ----------------
 *   text/markdown: per graph: (g, f) in canonical form; edges; degree sequence; graph6; class
 *                  size and complement information; the Riordan matrix B[i][j] = [t^i] g f^j
 *                  (rows 0..N-2) next to the lower triangle of the adjacency matrix (rows 2..N).
 *   csv:           k,n,g,f,edges,graph6,degree_sequence, then either
 *                  riordan_labelings,complement_riordan,self_complementary (unlabeled) or
 *                  complement_riordan_same_labeling (labeled).
 *   json:          an array of objects with the same fields plus "riordan_matrix".
 *   Exit code 0 on success, 2 on a usage error.
 *
 * Cost: generates L(N) graphs and computes L(N) canonical forms (about 3 µs each): N = 10 takes
 * well under a second, N = 14 about 1–2 minutes and 2 GB of memory.
 */
#include "rgraph/canon.hpp"
#include "rgraph/io.hpp"
#include "rgraph/properties.hpp"
#include "rgraph/recognize.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <tuple>

using namespace rgraph;

namespace {

/** One listed graph. */
struct Item {
    u64 g = 0, f = 0;          ///< canonical parameters of the listed Riordan labeling
    Graph G;                   ///< the labeled graph G_n(g, f)
    u64 labelings = 1;         ///< unlabeled mode: number of labeled Riordan graphs in the class
    int comp_riordan = -1;     ///< 1 if the complement is Riordan, 0 if not, -1 unknown
    bool selfcomp = false;     ///< unlabeled mode: the graph is isomorphic to its complement
};

/**
 * @brief Write the SVG drawing of one listed graph.
 * @param path  output file
 * @param it    the graph and its parameters
 * @param k     position in the listing (1-based), shown in the caption; 0 for --show (no number)
 */
void draw(const std::string& path, const Item& it, int k) {
    const std::string note = (k ? "#" + std::to_string(k) + ", " : std::string()) + std::to_string(edge_list(it.G).size()) + " edges";
    write_svg(path, it.G, it.g, it.f, note);
}

}  // namespace

/**
 * @brief Entry point: parse options, build the list, apply the filter, print it.
 * @return 0 on success, 2 on a usage error
 */
int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: rglist N [--unlabeled|--labeled] [--filter all|complement|selfcomplementary|connected]\n"
                             "              [--format text|markdown|csv|json] [--svg DIR] [--count]\n");
        return 2;
    }
    const int n = std::atoi(argv[1]);
    bool unlabeled = true, count_only = false;
    std::string filter = "all", format = "text", svgdir;
    bool show = false;
    u64 show_g = 0, show_f = 0;
    for (int i = 2; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--unlabeled") unlabeled = true;
        else if (a == "--labeled") unlabeled = false;
        else if (a == "--filter" && i + 1 < argc) filter = argv[++i];
        else if (a == "--format" && i + 1 < argc) format = argv[++i];
        else if (a == "--svg" && i + 1 < argc) svgdir = argv[++i];
        else if (a == "--count") count_only = true;
        else if (a == "--show" && i + 2 < argc) {
            try { show_g = parse_poly(argv[i + 1]); show_f = parse_poly(argv[i + 2]); }
            catch (const std::exception& e) { std::fprintf(stderr, "rglist: %s\n", e.what()); return 2; }
            if (show_f & 1) { std::fprintf(stderr, "rglist: f(0) must be 0\n"); return 2; }
            show = true; unlabeled = false; i += 2;
        }
        else { std::fprintf(stderr, "unknown option %s\n", a.c_str()); return 2; }
    }
    if (n < 1 || n > (unlabeled ? 14 : 16)) {
        std::fprintf(stderr, "rglist: n must be between 1 and %d for this mode\n", unlabeled ? 14 : 16);
        return 2;
    }
    if (!unlabeled && filter == "selfcomplementary") { std::fprintf(stderr, "selfcomplementary needs --unlabeled\n"); return 2; }

    const Enumeration E(n);
    std::vector<Item> items;
    if (show) {
        Item x;
        x.g = show_g; x.f = show_f;
        canonical_params(n, x.g, x.f);
        x.G = riordan_graph(n, x.g, x.f);
        x.comp_riordan = is_riordan_labeled(complement(x.G)) ? 1 : 0;
        items.push_back(x);
    } else if (unlabeled) {
        std::map<Hash128, Item> cls;  // class fingerprint -> best representative
        for (u64 idx = 0; idx < E.size(); ++idx) {
            u64 g, f;
            E.pair(idx, g, f);
            canonical_params(n, g, f);
            const Graph G = riordan_graph(n, g, f);
            const Hash128 h = class_hash(G);
            auto it = cls.find(h);
            if (it == cls.end()) { Item x; x.g = g; x.f = f; x.G = G; cls.emplace(h, x); continue; }
            Item& b = it->second;
            ++b.labelings;
            auto key = [](u64 gg, u64 ff) { return std::make_tuple(std::popcount(gg) + std::popcount(ff), gg, ff); };
            if (key(g, f) < key(b.g, b.f)) { b.g = g; b.f = f; b.G = G; }
        }
        for (auto& [h, x] : cls) {
            const Hash128 hc = class_hash(complement(x.G));
            x.comp_riordan = cls.count(hc) ? 1 : 0;
            x.selfcomp = (hc == h);
            items.push_back(x);
        }
        std::sort(items.begin(), items.end(), [](const Item& a, const Item& b) {
            const auto ea = edge_list(a.G).size(), eb = edge_list(b.G).size();
            if (ea != eb) return ea < eb;
            return std::make_tuple(std::popcount(a.g) + std::popcount(a.f), a.g, a.f) <
                   std::make_tuple(std::popcount(b.g) + std::popcount(b.f), b.g, b.f);
        });
    } else {
        for (u64 idx = 0; idx < E.size(); ++idx) {
            Item x;
            E.pair(idx, x.g, x.f);
            canonical_params(n, x.g, x.f);
            x.G = riordan_graph(n, x.g, x.f);
            x.comp_riordan = is_riordan_labeled(complement(x.G)) ? 1 : 0;
            items.push_back(x);
        }
    }
    std::vector<Item> out;
    for (auto& x : items) {
        if (filter == "complement" && x.comp_riordan != 1) continue;
        if (filter == "selfcomplementary" && !x.selfcomp) continue;
        if (filter == "connected" && !is_connected(x.G)) continue;
        if (filter != "all" && filter != "complement" && filter != "selfcomplementary" && filter != "connected") {
            std::fprintf(stderr, "unknown filter %s\n", filter.c_str()); return 2;
        }
        out.push_back(x);
    }
    if (count_only) { std::printf("%zu\n", out.size()); return 0; }
    if (!svgdir.empty()) std::filesystem::create_directories(svgdir);

    if (format == "csv") std::printf("k,n,g,f,edges,graph6,degree_sequence%s\n", unlabeled ? ",riordan_labelings,complement_riordan,self_complementary" : ",complement_riordan_same_labeling");
    if (format == "json") std::printf("[\n");
    if (format == "markdown") std::printf("# %s Riordan graphs of order %d (filter: %s): %zu\n\n", unlabeled ? "Unlabeled" : "Labeled", n, filter.c_str(), out.size());
    if (format == "text") {
        if (show) std::printf("Labeled Riordan graph of order %d\n", n);
        else std::printf("%s Riordan graphs of order %d (filter: %s): %zu\n", unlabeled ? "Unlabeled" : "Labeled", n, filter.c_str(), out.size());
    }

    for (size_t k = 0; k < out.size(); ++k) {
        const Item& x = out[k];
        const auto B = riordan_matrix(n, x.g, x.f);
        const auto e = edge_list(x.G);
        std::vector<int> deg(n);
        for (int v = 0; v < n; ++v) deg[v] = x.G.degree(v);
        std::vector<int> ds = deg; std::sort(ds.rbegin(), ds.rend());
        std::ostringstream es, dss;
        for (size_t i = 0; i < e.size(); ++i) es << (i ? " " : "") << e[i].first << e[i].second;
        for (size_t i = 0; i < ds.size(); ++i) dss << (i ? " " : "") << ds[i];
        const std::string name = "rg_" + std::to_string(n) + "_" + std::to_string(k + 1);
        if (!svgdir.empty()) draw(svgdir + "/" + name + ".svg", x, show ? 0 : int(k + 1));

        if (format == "csv") {
            std::printf("%zu,%d,%s,%s,%s,%s,%s", k + 1, n, poly_str(x.g).c_str(), poly_str(x.f).c_str(), es.str().c_str(), graph6(x.G).c_str(), dss.str().c_str());
            if (unlabeled) std::printf(",%llu,%d,%d\n", (unsigned long long)x.labelings, x.comp_riordan, int(x.selfcomp));
            else std::printf(",%d\n", x.comp_riordan);
            continue;
        }
        if (format == "json") {
            std::printf("  {\"k\": %zu, \"n\": %d, \"g\": \"%s\", \"f\": \"%s\", \"edges\": [", k + 1, n, poly_str(x.g).c_str(), poly_str(x.f).c_str());
            for (size_t i = 0; i < e.size(); ++i) std::printf("%s[%d,%d]", i ? "," : "", e[i].first, e[i].second);
            std::printf("], \"graph6\": \"%s\", \"riordan_matrix\": [", graph6(x.G).c_str());
            for (int i = 0; i < n - 1; ++i) {
                std::printf("%s[", i ? "," : "");
                for (int j = 0; j <= i; ++j) std::printf("%s%d", j ? "," : "", B[i][j]);
                std::printf("]");
            }
            std::printf("]");
            if (unlabeled) std::printf(", \"riordan_labelings\": %llu, \"complement_riordan\": %s, \"self_complementary\": %s",
                                       (unsigned long long)x.labelings, x.comp_riordan ? "true" : "false", x.selfcomp ? "true" : "false");
            else std::printf(", \"complement_riordan_same_labeling\": %s", x.comp_riordan ? "true" : "false");
            std::printf("}%s\n", k + 1 < out.size() ? "," : "");
            continue;
        }
        const bool md = format == "markdown";
        if (md) std::printf("## %zu. G_%d(%s, %s)\n\n", k + 1, n, poly_str(x.g).c_str(), poly_str(x.f).c_str());
        else std::printf("\n[%zu] G_%d(g, f) with g = %s,  f = %s\n", k + 1, n, poly_str(x.g).c_str(), poly_str(x.f).c_str());
        if (md && !svgdir.empty()) std::printf("![%s](%s/%s.svg)\n\n", name.c_str(), std::filesystem::path(svgdir).filename().c_str(), name.c_str());
        if (md) std::printf("- edges (%zu): %s\n- degree sequence: %s\n- graph6: `%s`\n", e.size(), es.str().c_str(), dss.str().c_str(), graph6(x.G).c_str());
        else std::printf("    edges (%zu): %s\n    degree sequence: %s\n    graph6: %s\n", e.size(), es.str().c_str(), dss.str().c_str(), graph6(x.G).c_str());
        if (unlabeled) {
            if (md) std::printf("- labeled Riordan graphs in this class: %llu\n- complement is Riordan: %s%s\n", (unsigned long long)x.labelings,
                                x.comp_riordan ? "yes" : "no", x.selfcomp ? " (self-complementary)" : "");
            else std::printf("    labeled Riordan graphs in this class: %llu;  complement Riordan: %s%s\n", (unsigned long long)x.labelings,
                             x.comp_riordan ? "yes" : "no", x.selfcomp ? " (self-complementary)" : "");
        } else {
            std::printf(md ? "- complement is a labeled Riordan graph (same labeling): %s\n" : "    complement labeled Riordan (same labeling): %s\n", x.comp_riordan ? "yes" : "no");
        }
        if (n >= 2) {
            std::printf(md ? "\nRiordan matrix B = [ [t^i] g f^j ] (rows i = 0..%d)      adjacency matrix A (lower triangle, rows 2..%d)\n\n```\n"
                           : "    Riordan matrix B (rows i = 0..%d)      adjacency lower triangle (rows 2..%d)\n", n - 2, n);
            for (int i = 0; i < n - 1; ++i) {
                std::string l = md ? "" : "    ";
                for (int j = 0; j <= i; ++j) l += (j ? " " : "") + std::to_string(B[i][j]);
                l += std::string(2 * (n - 1 - i), ' ') + "        ";
                l += std::to_string(i + 2) + ": ";
                for (int j = 0; j <= i; ++j) l += (j ? " " : "") + std::to_string(int(x.G.has_edge(i + 1, j)));
                std::printf("%s\n", l.c_str());
            }
            if (md) std::printf("```\n\n");
        }
    }
    if (format == "json") std::printf("]\n");
    return 0;
}
