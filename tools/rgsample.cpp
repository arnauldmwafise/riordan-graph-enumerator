/**
 * @file   tools/rgsample.cpp
 * @brief  rgsample — random sampling of labeled and unlabeled Riordan graphs, optionally
 *         conditioned on subgraph structure, with summary statistics.
 *
 * Purpose
 * -------
 * Exhaustive enumeration of unlabeled Riordan graphs currently stops at n = 17 (u(17) = 789399022).
 * Beyond that, random sampling is the practical way to study typical Riordan graphs. rgsample draws
 * graphs from well-defined distributions and reports their invariants. It also gives an unbiased
 * estimate of u(n) for orders where exact counting is out of reach.
 *
 * Sampling models (--model)
 * -------------------------
 *  labeled    Uniform over all L(n) = (4^(n-1)+2)/3 labeled Riordan graphs G_n(g,f). Each graph has
 *             unique canonical parameters (g mod t^(n-1), f mod t^(n-1-r)). The lowest degree r of g
 *             is drawn with probability proportional to 4^(n-2-r), and the remaining coefficient
 *             bits uniformly. For n <= 32 this is exactly uniform, the null graph included. For
 *             n > 32 the null graph (probability < 4^-31) is omitted.
 *  proper     Uniform over proper graphs: g(0) = 1, f'(0) = 1, other bits free (2^(2n-5) graphs).
 *  appell     Uniform over Appell-type graphs G_n(g, t) (Toeplitz graphs), g mod t^(n-1) free.
 *  bell       Uniform over Bell-type graphs G_n(g, t g) with g(0) = 1.
 *  unlabeled  Uniform over isomorphism classes (unlabeled Riordan graphs). A graph H is drawn from
 *             the labeled model, the size |C(H)| of its class is computed exactly with the
 *             labeling search of rgraph/labelings.hpp, and H is accepted with probability 1/|C(H)|.
 *             Each class is then accepted with probability |class|/L(n) * 1/|class| = 1/L(n), i.e.
 *             uniformly. The acceptance rate is u(n)/L(n), about 0.5. Feasible up to n ~ 25.
 *
 * Restrictions and structure filters (rejection sampling; uniform on the restricted set)
 * --------------------------------------------------------------------------------------
 *  --fprime 0|1            only graphs with f'(0) = 0 (resp. 1)
 *  --val-g R               only graphs whose g has lowest term t^R (vertex 1 adjacent first to R+2)
 *  --contains P            must contain pattern P as a subgraph       (repeatable)
 *  --contains-induced P    must contain P as an induced subgraph      (repeatable)
 *  --avoids P              must not contain P as a subgraph           (repeatable)
 *  --avoids-induced P      must not contain P as an induced subgraph  (repeatable)
 *  --connected             must be connected
 *  --edges A:B             number of edges between A and B
 *  P is a name (K3, K4, C4, C5, P4, claw, paw, diamond, 2K2, bull, house, K23, chair, cochair, …,
 *  see rgraph/properties.hpp) or a graph6 string. Filters are isomorphism-invariant, so with
 *  --model unlabeled the result is uniform over the unlabeled Riordan graphs satisfying them.
 *
 * Outputs
 * -------
 *  Per sample (text, csv or json): (g, f), edges, minimum/maximum degree, triangles, components,
 *  diameter (-1 if disconnected), Wiener index, graph6, and with --class-size the number |C(H)|
 *  of labeled Riordan graphs isomorphic to H, a completeness flag, and whether C(H) = {H, rev H}
 *  (the class consists of H and its reversal, the generic situation when f'(0) = 1).
 *  Summary: acceptance rate and mean / standard deviation / range of each invariant; the fraction
 *  of samples containing each --report pattern; and, for --model labeled --class-size without
 *  filters, the unbiased estimate u(n) ≈ L(n)·mean(1/|C(H)|) with a 95% confidence interval.
 *
 * Notes
 * -----
 *  - Filters work by rejection: the cost per accepted sample is 1/(probability of the filter). For
 *    rare structures (e.g. triangle-free graphs of order 30, probability about 10^-5) raise
 *    --max-attempts or relax the condition; the summary reports the observed acceptance rate.
 *  - Class sizes (--class-size, --model unlabeled) need a search over Riordan labelings. It is fast
 *    for graphs with f'(0) = 1 (a few ms at n = 20). Graphs with f'(0) = 0 and val f >= 3 can need
 *    many millions of nodes from n ~ 15 on; searches over --node-limit are reported, and the u(n)
 *    estimate then becomes a conservative interval.
 *
 * Usage
 * -----
 *   rgsample N [--samples S] [--seed X] [--model labeled|proper|appell|bell|unlabeled]
 *              [--fprime 0|1] [--val-g R] [--contains P] [--contains-induced P] [--avoids P]
 *              [--avoids-induced P] [--connected] [--edges A:B] [--report P]
 *              [--class-size] [--node-limit L] [--max-attempts M]
 *              [--format text|csv|json] [--summary] [--svg DIR]
 */
#include "rgraph/io.hpp"
#include "rgraph/labelings.hpp"
#include "rgraph/properties.hpp"
#include "rgraph/recognize.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <random>
#include <string>
#include <vector>

using namespace rgraph;

namespace {

/** Command-line configuration. */
struct Config {
    int n = 0;
    long long samples = 10;
    unsigned long long seed = 1;
    std::string model = "labeled";
    int fprime = -1;     ///< -1: no restriction
    int val_g = -1;      ///< -1: no restriction
    std::vector<std::pair<std::string, Graph>> contains, contains_ind, avoids, avoids_ind, report;
    bool connected = false;
    int emin = 0, emax = 1 << 30;
    bool class_size = false;
    unsigned long long node_limit = 50000000ULL;
    long long max_attempts = -1;  ///< -1: 1000 * samples + 10^6
    std::string format = "text", svgdir;
    bool summary_only = false;
};

/** Running mean / variance / range (Welford). */
struct Stat {
    long long k = 0;
    double mean = 0, m2 = 0, lo = 1e300, hi = -1e300;
    void add(double x) {
        ++k;
        const double d = x - mean;
        mean += d / k;
        m2 += d * (x - mean);
        lo = std::min(lo, x);
        hi = std::max(hi, x);
    }
    double sd() const { return k > 1 ? std::sqrt(m2 / (k - 1)) : 0; }
};

/**
 * @brief Draw canonical parameters (g, f) for the chosen model.
 * @param c    configuration (model, n, fprime, val_g)
 * @param rng  random generator
 * @param g, f out: parameters; the graph is G_n(g, f)
 */
void draw_params(const Config& c, std::mt19937_64& rng, u64& g, u64& f) {
    const int n = c.n, N = n - 1;
    auto bits = [&](int k) -> u64 { return k <= 0 ? 0 : (rng() & low_mask(k)); };
    if (c.model == "appell") { g = bits(N); f = 2; canonical_params(n, g, f); return; }
    if (c.model == "bell") { g = 1 | (bits(N - 1) << 1); f = (g << 1) & low_mask(N); canonical_params(n, g, f); return; }
    if (c.model == "proper") { g = 1 | (bits(N - 1) << 1); f = 2 | (bits(N - 2) << 2); canonical_params(n, g, f); return; }
    // labeled / unlabeled: uniform over canonical parameters
    int r;
    if (c.val_g >= 0) {
        r = c.val_g;
    } else if (n <= 32) {
        const u64 L = labeled_count(n);
        const u64 idx = std::uniform_int_distribution<u64>(0, L - 1)(rng);
        Enumeration(n).pair(idx, g, f);
        return;
    } else {
        // P(k) ∝ 4^k for k = 0..N-1 (k = N-1-r): from the top, stop at k with prob 3/4
        int k = N - 1;
        while (k > 0 && (rng() & 3) == 0) --k;
        r = N - 1 - k;
    }
    const int k = N - 1 - r;
    g = (u64{1} << r) | (bits(k) << (r + 1));
    f = bits(k) << 1;
}

/** @brief True iff (g, f) satisfies the parameter restrictions --fprime and --val-g. */
bool params_ok(const Config& c, u64 g, u64 f) {
    if (c.val_g >= 0 && (g == 0 || std::countr_zero(g) != c.val_g)) return false;
    if (c.fprime >= 0) {
        const bool has_f1 = g != 0 && (c.n - 2 - std::countr_zero(g)) >= 1;  // f_1 is a free coefficient
        const int f1 = has_f1 ? int(f >> 1 & 1) : 0;
        if (c.fprime == 1 && !(has_f1 && f1 == 1)) return false;
        if (c.fprime == 0 && f1 == 1) return false;
    }
    return true;
}

/** @brief True iff the graph passes every structure filter. */
bool structure_ok(const Config& c, const Graph& G) {
    if (c.connected && !is_connected(G)) return false;
    const int m = edge_count(G);
    if (m < c.emin || m > c.emax) return false;
    for (auto& [nm, P] : c.contains) if (!contains_subgraph(G, P, false)) return false;
    for (auto& [nm, P] : c.contains_ind) if (!contains_subgraph(G, P, true)) return false;
    for (auto& [nm, P] : c.avoids) if (contains_subgraph(G, P, false)) return false;
    for (auto& [nm, P] : c.avoids_ind) if (contains_subgraph(G, P, true)) return false;
    return true;
}

int usage() {
    std::fprintf(stderr,
                 "usage: rgsample N [--samples S] [--seed X] [--model labeled|proper|appell|bell|unlabeled]\n"
                 "                  [--fprime 0|1] [--val-g R] [--contains P] [--contains-induced P] [--avoids P]\n"
                 "                  [--avoids-induced P] [--connected] [--edges A:B] [--report P]\n"
                 "                  [--class-size] [--node-limit L] [--max-attempts M]\n"
                 "                  [--format text|csv|json] [--summary] [--svg DIR]\n"
                 "patterns P: K2 K3 K4 K5 P3 P4 P5 C4 C5 C6 claw star4 paw diamond 2K2 bull house K23 chair cochair\n"
                 "            C4+K1 K4+K1, or a graph6 string\n");
    return 2;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) return usage();
    Config c;
    c.n = std::atoi(argv[1]);
    auto pattern = [&](const char* s, std::vector<std::pair<std::string, Graph>>& v) {
        Graph P;
        if (!pattern_graph(s, P)) { std::fprintf(stderr, "rgsample: unknown pattern '%s'\n", s); std::exit(2); }
        v.push_back({s, P});
    };
    for (int i = 2; i < argc; ++i) {
        const std::string a = argv[i];
        const bool has = i + 1 < argc;
        if (a == "--samples" && has) c.samples = std::atoll(argv[++i]);
        else if (a == "--seed" && has) c.seed = std::strtoull(argv[++i], nullptr, 10);
        else if (a == "--model" && has) c.model = argv[++i];
        else if (a == "--fprime" && has) c.fprime = std::atoi(argv[++i]);
        else if (a == "--val-g" && has) c.val_g = std::atoi(argv[++i]);
        else if (a == "--contains" && has) pattern(argv[++i], c.contains);
        else if (a == "--contains-induced" && has) pattern(argv[++i], c.contains_ind);
        else if (a == "--avoids" && has) pattern(argv[++i], c.avoids);
        else if (a == "--avoids-induced" && has) pattern(argv[++i], c.avoids_ind);
        else if (a == "--report" && has) pattern(argv[++i], c.report);
        else if (a == "--connected") c.connected = true;
        else if (a == "--edges" && has) {
            const std::string r = argv[++i];
            const auto p = r.find(':');
            if (p == std::string::npos) return usage();
            c.emin = std::atoi(r.substr(0, p).c_str());
            c.emax = std::atoi(r.substr(p + 1).c_str());
        }
        else if (a == "--class-size") c.class_size = true;
        else if (a == "--node-limit" && has) c.node_limit = std::strtoull(argv[++i], nullptr, 10);
        else if (a == "--max-attempts" && has) c.max_attempts = std::atoll(argv[++i]);
        else if (a == "--format" && has) c.format = argv[++i];
        else if (a == "--summary") c.summary_only = true;
        else if (a == "--svg" && has) c.svgdir = argv[++i];
        else { std::fprintf(stderr, "rgsample: unknown option %s\n", a.c_str()); return usage(); }
    }
    if (c.n < 2 || c.n > 64) { std::fprintf(stderr, "rgsample: need 2 <= N <= 64\n"); return 2; }
    if (c.model != "labeled" && c.model != "proper" && c.model != "appell" && c.model != "bell" && c.model != "unlabeled")
        return usage();
    if (c.val_g >= c.n - 1) { std::fprintf(stderr, "rgsample: --val-g must be < N-1\n"); return 2; }
    if (c.max_attempts < 0) c.max_attempts = 1000 * c.samples + 1000000;
    const bool need_class = c.class_size || c.model == "unlabeled";
    if (!c.svgdir.empty()) std::filesystem::create_directories(c.svgdir);

    std::mt19937_64 rng(c.seed * 0x9E3779B97F4A7C15ULL + c.n);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    Stat st_edges, st_maxdeg, st_mindeg, st_tri, st_comp, st_diam, st_wiener, st_cls, st_inv_lo, st_inv_hi;
    long long accepted = 0, attempts = 0, incomplete = 0, connected_count = 0, revpairs = 0, class_known = 0;
    std::vector<long long> rep_sub(c.report.size(), 0), rep_ind(c.report.size(), 0);

    FILE* out = stdout;
    const bool text = c.format == "text";
    if (!c.summary_only) {
        if (c.format == "csv")
            std::fprintf(out, "k,n,g,f,edges,min_degree,max_degree,triangles,components,diameter,wiener,graph6%s\n",
                         need_class ? ",class_size,class_complete,reversal_pair" : "");
        if (c.format == "json") std::fprintf(out, "[\n");
    }
    while (accepted < c.samples && attempts < c.max_attempts) {
        ++attempts;
        u64 g, f;
        draw_params(c, rng, g, f);
        if (!params_ok(c, g, f)) continue;
        const Graph G = riordan_graph(c.n, g, f);
        if (!structure_ok(c, G)) continue;
        long long cls = -1;
        bool complete = true;
        int revpair = -1;  // 1 if C(H) = {H, rev H} with rev H != H (the reversal symmetry), 0 otherwise
        if (need_class) {
            const LabelingResult R = riordan_labelings(G, c.node_limit);
            complete = R.complete;
            cls = (long long)R.labelings.size();
            if (complete) {
                const Graph Rv = reversed(G);
                const auto pr = recognize_labeled(Rv);
                revpair = (cls == 2 && pr && !(Rv == G) && R.labelings.count({pr->g, pr->f})) ? 1 : 0;
            }
            if (!complete) {
                ++incomplete;
                if (c.model == "unlabeled") continue;  // cannot decide acceptance; reported in summary
            }
            if (c.model == "unlabeled" && U(rng) >= 1.0 / double(cls)) continue;
        }
        ++accepted;
        const int m = edge_count(G);
        const auto ds = degree_sequence(G);
        const long long tri = triangle_count(G);
        const int comp = component_count(G);
        int diam;
        long long wien;
        distance_stats(G, diam, wien);
        st_edges.add(m); st_maxdeg.add(ds.front()); st_mindeg.add(ds.back()); st_tri.add(double(tri)); st_comp.add(comp);
        if (diam >= 0) { st_diam.add(diam); st_wiener.add(double(wien)); ++connected_count; }
        if (need_class && complete) { ++class_known; revpairs += (revpair == 1); }
        if (need_class && complete) { st_cls.add(double(cls)); st_inv_lo.add(1.0 / double(cls)); st_inv_hi.add(1.0 / double(cls)); }
        if (need_class && !complete) {  // |C(H)| >= labelings found so far: 1/|C(H)| lies in (0, 1/found]
            st_inv_lo.add(0.0);
            st_inv_hi.add(1.0 / double(std::max<long long>(cls, 1)));
        }
        for (std::size_t i = 0; i < c.report.size(); ++i) {
            rep_sub[i] += contains_subgraph(G, c.report[i].second, false);
            rep_ind[i] += contains_subgraph(G, c.report[i].second, true);
        }
        if (!c.svgdir.empty())
            write_svg(c.svgdir + "/sample_" + std::to_string(c.n) + "_" + std::to_string(accepted) + ".svg", G, g, f,
                      "sample " + std::to_string(accepted) + ", " + std::to_string(m) + " edges");
        if (c.summary_only) continue;
        const std::string g6 = graph6(G);
        if (c.format == "csv") {
            std::fprintf(out, "%lld,%d,%s,%s,%d,%d,%d,%lld,%d,%d,%lld,%s", accepted, c.n, poly_str(g).c_str(), poly_str(f).c_str(), m,
                         ds.back(), ds.front(), tri, comp, diam, wien, g6.c_str());
            if (need_class) std::fprintf(out, ",%lld,%d,%d", cls, int(complete), revpair);
            std::fprintf(out, "\n");
        } else if (c.format == "json") {
            std::fprintf(out, "%s  {\"k\": %lld, \"n\": %d, \"g\": \"%s\", \"f\": \"%s\", \"edges\": %d, \"min_degree\": %d, "
                              "\"max_degree\": %d, \"triangles\": %lld, \"components\": %d, \"diameter\": %d, \"wiener\": %lld, "
                              "\"graph6\": \"%s\"",
                         accepted > 1 ? ",\n" : "", accepted, c.n, poly_str(g).c_str(), poly_str(f).c_str(), m, ds.back(), ds.front(), tri,
                         comp, diam, wien, g6.c_str());
            if (need_class) std::fprintf(out, ", \"class_size\": %lld, \"class_complete\": %s, \"reversal_pair\": %d", cls, complete ? "true" : "false", revpair);
            std::fprintf(out, "}");
        } else {
            std::fprintf(out, "[%lld] G_%d(%s, %s)\n      edges %d, degrees %d..%d, triangles %lld, components %d, diameter %d, graph6 %s",
                         accepted, c.n, poly_str(g).c_str(), poly_str(f).c_str(), m, ds.back(), ds.front(), tri, comp, diam, g6.c_str());
            if (need_class) std::fprintf(out, ", class size %lld%s%s", cls, complete ? "" : " (search incomplete)", revpair == 1 ? " = {H, rev H}" : "");
            std::fprintf(out, "\n");
        }
    }
    if (!c.summary_only && c.format == "json") std::fprintf(out, "\n]\n");

    // ---- summary (stdout for text/--summary, stderr for csv/json so the data stay parseable)
    FILE* s = (text || c.summary_only) ? stdout : stderr;
    std::fprintf(s, "%s# rgsample n=%d model=%s seed=%llu: %lld samples accepted out of %lld draws (rate %.4f)\n",
                 text && !c.summary_only ? "\n" : "", c.n, c.model.c_str(), c.seed, accepted, attempts,
                 attempts ? double(accepted) / attempts : 0.0);
    if (accepted < c.samples) std::fprintf(s, "# WARNING: stopped at --max-attempts before reaching %lld samples\n", c.samples);
    if (incomplete)
        std::fprintf(s, "# WARNING: %lld labeling searches hit --node-limit%s\n", incomplete,
                     c.model == "unlabeled" ? " and were skipped (slight bias against large classes)" : "");
    auto row = [&](const char* name, const Stat& x) {
        if (x.k) std::fprintf(s, "#   %-22s mean %10.4f   sd %9.4f   min %8.0f   max %8.0f   (%lld values)\n", name, x.mean, x.sd(), x.lo, x.hi, x.k);
    };
    row("edges", st_edges);
    row("minimum degree", st_mindeg);
    row("maximum degree", st_maxdeg);
    row("triangles", st_tri);
    row("components", st_comp);
    row("diameter (connected)", st_diam);
    row("Wiener index (conn.)", st_wiener);
    row("class size |C(H)|", st_cls);
    if (accepted) std::fprintf(s, "#   connected fraction     %.4f\n", double(connected_count) / accepted);
    if (class_known)
        std::fprintf(s, "#   fraction with C(H) = {H, rev H}: %.4f (of %lld complete class computations)\n", double(revpairs) / class_known, class_known);
    for (std::size_t i = 0; i < c.report.size(); ++i)
        std::fprintf(s, "#   pattern %-10s contained as subgraph %.4f, as induced subgraph %.4f\n", c.report[i].first.c_str(),
                     accepted ? double(rep_sub[i]) / accepted : 0.0, accepted ? double(rep_ind[i]) / accepted : 0.0);
    const bool unfiltered = c.fprime < 0 && c.val_g < 0 && c.contains.empty() && c.contains_ind.empty() && c.avoids.empty() &&
                            c.avoids_ind.empty() && !c.connected && c.emin == 0 && c.emax == (1 << 30);
    if (c.model == "labeled" && c.class_size && unfiltered && st_inv_lo.k > 1) {
        // u(n) = L(n) * E[1/|C(H)|] for H uniform over labeled Riordan graphs (exact identity).
        const long double L = (std::pow(4.0L, c.n - 1) + 2.0L) / 3.0L;
        const double lo = st_inv_lo.mean - 1.96 * st_inv_lo.sd() / std::sqrt(double(st_inv_lo.k));
        const double hi = st_inv_hi.mean + 1.96 * st_inv_hi.sd() / std::sqrt(double(st_inv_hi.k));
        if (incomplete == 0) {
            std::fprintf(s, "#   estimate of u(n)/L(n) = mean 1/|C(H)| = %.5f, 95%% interval [%.5f, %.5f]\n", st_inv_lo.mean, lo, hi);
            std::fprintf(s, "#   estimate of u(%d) = %.6Lg, 95%% interval [%.6Lg, %.6Lg]\n", c.n, L * st_inv_lo.mean, L * lo, L * hi);
        } else {
            std::fprintf(s, "#   incomplete searches only bound 1/|C(H)| in (0, 1/found]; conservative 95%% interval:\n");
            std::fprintf(s, "#   u(n)/L(n) in [%.5f, %.5f],  u(%d) in [%.6Lg, %.6Lg]\n", lo, hi, c.n, L * lo, L * hi);
            std::fprintf(s, "#   (raise --node-limit to tighten; see the documentation for the cost)\n");
        }
    }
    return 0;
}
