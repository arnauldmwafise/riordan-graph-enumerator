/**
 * @file   analysis/bell.cpp
 * @brief  bell — diameter, Wiener index and independence number of io-decomposable Riordan graphs
 *         of the Bell type G_n(g, tg), including the Catalan and Pascal graphs.
 *
 * Context
 * -------
 * A labeled Riordan graph G_n(g, f) has vertices 1..n; for i > j, ij is an edge iff
 * [t^(i-2)] g f^(j-1) = 1 over Z_2 (g, f power series over Z_2, f(0) = 0). There are
 * L(n) = (4^(n-1)+2)/3 labeled Riordan graphs and u(n) unlabeled ones. This program studies the
 * Bell type f = t g with g(0) = 1. Such a graph is io-decomposable when
 * [t^j] g = [t^(2j+1)] g for every j >= 0 (Cheon et al. 2019, Lemma 4.4). Then the odd
 * coefficients of g are determined by the even ones: g_j = g_(q-1), where q is the odd part of
 * j + 1. With g_0 = 1, the free bits are g_2, g_4, ..., g_(2m) with 2m <= n-2, so there are
 * 2^(ceil(n/2)-1) io-decomposable Bell graphs of order n. Two of them are named:
 *  - Catalan graph CG_n: g_k = 1 iff k+1 is a power of 2 (g = sum of t^(2^i - 1));
 *  - Pascal graph PG_n:  g = 1/(1+t), i.e. every coefficient is 1.
 * Problem 6 and Conjecture 2 of Cheon et al. (2019) concern these graphs; Conjecture 2 is about
 * the diameter of CG_n at n = 2^k (diameter k). The "all" mode lists, for n = 2^k, every
 * io-decomposable graph with diameter k, which is the data used for that conjecture.
 *
 * Method
 * ------
 * The graph is built column by column: column j (0-based) of the Riordan matrix is t^j g^(j+1)
 * mod t^(n-1), obtained by repeated multiplication by tg; bit k of column j is the edge between
 * vertices k+1 and j (0-based). For n <= 64 every adjacency row is one 64-bit word; for larger n
 * (only modes named/pow2) rows are bitsets of ceil(n/64) words. Distances come from bit-parallel
 * BFS from every vertex. The independence number (n <= 64 only) is computed exactly by branch and
 * bound: branch on a vertex of minimum degree among the candidates, bound by a greedy partition of
 * the candidates into cliques.
 *
 * Usage
 * -----
 *   bell named KIND NMAX
 *   bell pow2  KIND KMAX
 *   bell all   N [--threads T] [--no-alpha] [--part I P]
 *
 * Arguments and options
 * ---------------------
 *   KIND          catalan | pascal (anything else: "unknown kind" on stderr, exit 2)
 *   named ... NMAX   orders n = 1..NMAX
 *   pow2  ... KMAX   orders n = 2^k, k = 1..KMAX (any size; n > 64 uses bitset graphs)
 *   all N         all 2^(ceil(N/2)-1) io-decomposable Bell graphs of order N (1 <= N <= 64)
 *   --threads T   worker threads for "all" (default: hardware concurrency)
 *   --no-alpha    "all": skip the independence number (it dominates the cost for large N)
 *   --part I P    "all": process only slice I (0-based) of P equal slices of the range of free
 *                 bits, so a long run can be split into resumable pieces
 *   BELL_BIG      environment variable: if set, named/pow2 use the bitset path even for n <= 64
 *                 (no independence number), as a cross-check of the two implementations
 * With fewer than two arguments the program prints "see header" on stderr and exits with 2; an
 * unknown mode exits with 2 silently.
 *
 * Output (stdout)
 * ---------------
 *  named/pow2: a header line "# KIND graphs: n edges diam wiener [alpha]", then one line per n:
 *      "n edges diam wiener alpha"   (n <= 64, word path)
 *      "n edges diam wiener"         (bitset path)
 *    diam is the largest finite distance (the graphs here are connected for n >= 2).
 *  all: one line
 *      "n=N graphs=COUNT diam: d:c d:c ... wiener:[MIN,MAX] alpha: a:c a:c ..."
 *    (the alpha part is omitted with --no-alpha; COUNT is the size of the processed slice), and,
 *    when N = 2^K, a second line
 *      "  graphs with diam = log2 n = K: C g=HEX g=HEX(Catalan) ..."
 *    listing in increasing order every g (bits of g mod t^(N-1), hexadecimal) of diameter K; the
 *    Catalan graph is marked.
 *
 * Cost
 * ----
 *  "all" is linear in the number of graphs 2^(ceil(N/2)-1), about N^2 word operations each, plus the
 *  (worst-case exponential) independence number. N = 34 (2^16 graphs) takes under a second on
 *  one core; N = 64 (2^31 graphs, --no-alpha) takes hours and is run as 32 slices by
 *  scripts/bell_all64.sh. named/pow2 use BFS from every vertex, O(n^3/64) word operations on
 *  bitsets; "pow2 catalan 13" (up to n = 8192) takes a few seconds.
 *
 * Results are recorded in data/bell/ (all_small.txt, all_pow2.txt, all_64.txt and all64/part*.txt,
 * catalan_named.txt, catalan_pow2.txt, pascal_named.txt, pascal_pow2.txt) and summarised in
 * docs/COMPUTATIONS.md, section 4.4.
 */
#include <algorithm>
#include <atomic>
#include <bit>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

using u64 = std::uint64_t;

// ---------------------------------------------------------------- general n (bitset graphs)

/** @brief Undirected graph of arbitrary order with bitset adjacency rows (n rows of W words). */
struct BigGraph {
    int n, W;
    std::vector<u64> a;  // n rows of W words
    /** @brief Edgeless graph on @p n_ vertices. */
    BigGraph(int n_) : n(n_), W((n_ + 63) / 64), a(std::size_t(n_) * W, 0) {}
    /** @brief Adjacency row of vertex @p v (W words; bit u = edge uv). */
    u64* row(int v) { return &a[std::size_t(v) * W]; }
    const u64* row(int v) const { return &a[std::size_t(v) * W]; }
    /** @brief Add the edge uv (0-based vertices). */
    void add(int u, int v) { row(u)[v >> 6] |= u64{1} << (v & 63); row(v)[u >> 6] |= u64{1} << (u & 63); }
};

using Poly = std::vector<u64>;  // bits 0..L-1

/** @brief Coefficient of t^k in the Z_2 polynomial @p p. @return true iff the bit is 1 */
static bool getb(const Poly& p, int k) { return (p[k >> 6] >> (k & 63)) & 1; }
/** @brief Set the coefficient of t^k in @p p to 1. */
static void setb(Poly& p, int k) { p[k >> 6] |= u64{1} << (k & 63); }

/**
 * @brief Product of two Z_2 polynomials truncated to t^L: r = (a * b) mod t^L.
 * @param a, b  multi-word polynomials of equal word count
 * @param L     truncation length (number of kept coefficients)
 * @return      product with the same word count as @p a, bits >= L cleared
 */
static Poly mul(const Poly& a, const Poly& b, int L) {
    const int W = int(a.size());
    Poly r(W, 0);
    for (int s = 0; s < L; ++s) {
        if (!getb(b, s)) continue;
        const int ws = s >> 6, bs = s & 63;
        for (int i = W - 1; i >= ws; --i) {
            u64 x = a[i - ws] << bs;
            if (bs && i - ws - 1 >= 0) x |= a[i - ws - 1] >> (64 - bs);
            r[i] ^= x;
        }
    }
    const int rem = L & 63;
    if (rem) r[W - 1] &= (u64{1} << rem) - 1;
    for (int i = (L + 63) / 64; i < W; ++i) r[i] = 0;
    return r;
}

/**
 * @brief g for a named Bell graph, as a multi-word polynomial with L coefficients.
 * @param kind  "pascal" (all coefficients 1) or "catalan" (g_k = 1 iff k+1 is a power of 2);
 *              anything else prints "unknown kind" and exits with status 2
 * @param L     number of coefficients (n-1 for order n)
 * @return      g, with one spare word
 */
static Poly make_g(const std::string& kind, int L) {
    Poly g((L + 63) / 64 + 1, 0);
    for (int k = 0; k < L; ++k) {
        if (kind == "pascal") setb(g, k);
        else if (kind == "catalan") { if (std::has_single_bit(unsigned(k + 1))) setb(g, k); }
        else { std::fprintf(stderr, "unknown kind\n"); std::exit(2); }
    }
    return g;
}

/**
 * @brief Build the Bell-type Riordan graph G_n(g, tg) of any order on bitsets.
 *
 * Column j0 (0-based) of the Riordan matrix is t^j0 g^(j0+1) mod t^(n-1); bit k of it is the edge
 * between 0-based vertices k+1 and j0.
 * @param n  order of the graph
 * @param g  first generating function (at least n-1 coefficients)
 * @return   the graph
 */
static BigGraph bell_graph(int n, const Poly& g) {
    BigGraph G(n);
    const int L = n - 1;
    Poly tg(g.size(), 0);
    for (int k = 0; k + 1 < L; ++k) if (getb(g, k)) setb(tg, k + 1);
    Poly col = g;
    for (int i = (L + 63) / 64; i < int(col.size()); ++i) col[i] = 0;
    if (L & 63) col[(L - 1) >> 6] &= (u64{1} << (L & 63)) - 1;
    for (int j0 = 0; j0 + 1 < n; ++j0) {
        for (int k = j0; k < L; ++k) if (getb(col, k)) G.add(k + 1, j0);
        col = mul(col, tg, L);
    }
    return G;
}

/** @brief Distance statistics of a graph. */
struct DistStats { int diam = 0; u64 wiener = 0; bool connected = true; };

/**
 * @brief Diameter and Wiener index of a bitset graph by BFS from every vertex.
 * @param G  graph
 * @return   diam = maximum eccentricity over reachable vertices; wiener = sum of distances over
 *           unordered pairs of mutually reachable vertices; connected = false if some BFS did
 *           not reach all vertices
 */
static DistStats distances(const BigGraph& G) {
    const int n = G.n, W = G.W;
    DistStats st;
    std::vector<u64> vis(W), fr(W), nx(W);
    for (int s = 0; s < n; ++s) {
        std::fill(vis.begin(), vis.end(), 0); std::fill(fr.begin(), fr.end(), 0);
        vis[s >> 6] |= u64{1} << (s & 63); fr = vis;
        int d = 0, seen = 1;
        for (;;) {
            std::fill(nx.begin(), nx.end(), 0);
            for (int w = 0; w < W; ++w)
                for (u64 x = fr[w]; x; x &= x - 1) {
                    const u64* r = G.row(w * 64 + std::countr_zero(x));
                    for (int i = 0; i < W; ++i) nx[i] |= r[i];
                }
            int c = 0;
            for (int i = 0; i < W; ++i) { nx[i] &= ~vis[i]; vis[i] |= nx[i]; c += std::popcount(nx[i]); }
            if (!c) break;
            ++d; seen += c; st.wiener += u64(d) * c;
            fr = nx;
        }
        if (seen < n) st.connected = false;
        st.diam = std::max(st.diam, d);
    }
    st.wiener /= 2;
    return st;
}

// ---------------------------------------------------------------- n <= 64 (one word per row)
/** @brief Graph with n <= 64 vertices; adj[v] bit u = edge uv. */
struct Small { int n; u64 adj[64]; };

/**
 * @brief Product of two single-word Z_2 polynomials, truncated by @p mask.
 * @param a, b  polynomials (bit k = coefficient of t^k)
 * @param mask  low-bit mask of the coefficients to keep
 * @return      (a * b) & mask
 */
static u64 mul64(u64 a, u64 b, u64 mask) { u64 r = 0; for (; b; b &= b - 1) r ^= a << std::countr_zero(b); return r & mask; }

/**
 * @brief Build G_n(g, tg) for n <= 64 with one word per adjacency row.
 * @param n  order (1..64)
 * @param g  first generating function (bit k = g_k); only g mod t^(n-1) is used
 * @return   the graph
 */
static Small small_bell(int n, u64 g) {
    Small G; G.n = n; std::memset(G.adj, 0, sizeof G.adj);
    const u64 mask = n - 1 >= 64 ? ~u64{0} : ((u64{1} << (n - 1)) - 1);
    const u64 f = (g << 1) & mask;
    u64 col = g & mask;
    for (int j0 = 0; j0 + 1 < n && col; ++j0) {
        const u64 lower = col << 1;
        G.adj[j0] |= lower;
        for (u64 x = lower; x; x &= x - 1) G.adj[std::countr_zero(x)] |= u64{1} << j0;
        col = mul64(col, f, mask);
    }
    return G;
}

/**
 * @brief Diameter and Wiener index of a graph with n <= 64 by bit-parallel BFS from every vertex.
 * @param G  graph
 * @return   as distances(); a BFS that stops early marks the graph disconnected
 */
static DistStats small_dist(const Small& G) {
    DistStats st;
    const u64 all = G.n == 64 ? ~u64{0} : ((u64{1} << G.n) - 1);
    for (int s = 0; s < G.n; ++s) {
        u64 vis = u64{1} << s, fr = vis;
        int d = 0;
        while (vis != all) {
            u64 nx = 0;
            for (u64 x = fr; x; x &= x - 1) nx |= G.adj[std::countr_zero(x)];
            nx &= ~vis;
            if (!nx) { st.connected = false; break; }
            ++d; st.wiener += u64(d) * std::popcount(nx);
            vis |= nx; fr = nx;
        }
        st.diam = std::max(st.diam, d);
    }
    st.wiener /= 2;
    return st;
}

/**
 * @brief Branch-and-bound step of the maximum independent set search.
 *
 * Bound: cur plus the number of cliques in a greedy clique partition of @p cand (a clique cover
 * bounds any independent set from above). Branch on a candidate of minimum degree within @p cand:
 * take it (remove its neighbours) or skip it; if it has no neighbours in @p cand only "take" is
 * explored.
 * @param G     graph
 * @param cand  vertices still allowed in the independent set
 * @param cur   size of the independent set chosen so far
 * @param best  best size found so far
 * @return      max(best, size of the largest independent set extending the current one)
 */
static int alpha_rec(const Small& G, u64 cand, int cur, int best) {
    if (!cand) return std::max(cur, best);
    if (cur + std::popcount(cand) <= best) return best;
    // bound: clique cover of cand in the graph = colouring of the complement; use a cheap bound:
    // greedily partition cand into cliques of G
    int cover = 0;
    for (u64 c = cand; c; ++cover) {
        u64 q = c & -c, rest = c & ~q;
        u64 common = G.adj[std::countr_zero(q)] & rest;
        while (common) { const u64 b = common & -common; q |= b; common &= G.adj[std::countr_zero(b)]; }
        c &= ~q;
    }
    if (cur + cover <= best) return best;
    // branch on a vertex of minimum degree inside cand
    int v = -1, md = 1 << 30;
    for (u64 x = cand; x; x &= x - 1) {
        const int w = std::countr_zero(x), d = std::popcount(G.adj[w] & cand);
        if (d < md) { md = d; v = w; }
    }
    best = alpha_rec(G, cand & ~G.adj[v] & ~(u64{1} << v), cur + 1, best);  // take v
    if (md == 0) return best;
    best = alpha_rec(G, cand & ~(u64{1} << v), cur, best);  // skip v
    return best;
}
/**
 * @brief Independence number (exact) of a graph with n <= 64.
 * @param G  graph
 * @return   size of a maximum independent set
 */
static int alpha(const Small& G) {
    const u64 all = G.n == 64 ? ~u64{0} : ((u64{1} << G.n) - 1);
    return alpha_rec(G, all, 0, 0);
}

/**
 * @brief g mod t^(n-1) of a named Bell graph as a single word (n <= 64).
 * @param kind  "pascal" sets every coefficient; any other value gives the Catalan g
 *              (g_k = 1 iff k+1 is a power of 2)
 * @param n     order of the graph
 * @return      coefficient bits g_0..g_(n-2)
 */
static u64 small_g(const std::string& kind, int n) {
    u64 g = 0;
    for (int k = 0; k <= n - 2; ++k)
        if (kind == "pascal" || std::has_single_bit(unsigned(k + 1))) g |= u64{1} << k;
    return g;
}

/**
 * @brief The io-decomposable g mod t^(n-1) with given free even coefficients.
 *
 * g_0 = 1; g_(2m) = bit m-1 of @p e for m >= 1; for odd j, g_j = g_(q-1) where q is the odd part
 * of j+1, which realises [t^j] g = [t^(2j+1)] g.
 * @param n  order of the graph (n <= 64)
 * @param e  free bits: bit m-1 is g_(2m), for 1 <= m <= (n-2)/2
 * @return   coefficient bits g_0..g_(n-2)
 */
static u64 iodec_g(int n, u64 e) {
    u64 g = 1 | 0;
    auto even = [&](int m) { return m == 0 ? 1 : int((e >> (m - 1)) & 1); };  // value of g_(2m)
    for (int j = 1; j <= n - 2; ++j) {
        int b;
        if (j % 2 == 0) b = even(j / 2);
        else { int q = j + 1; while (q % 2 == 0) q /= 2; b = even((q - 1) / 2); }
        if (b) g |= u64{1} << j;
    }
    return g;
}

/**
 * @brief Entry point: dispatch on the mode (named, pow2, all); see the file header.
 * @return 0 on success, 2 on a usage error
 */
int main(int argc, char** argv) {
    if (argc < 3) { std::fprintf(stderr, "see header\n"); return 2; }
    const std::string mode = argv[1];
    if (mode == "named" || mode == "pow2") {
        const std::string kind = argv[2];
        if (kind != "catalan" && kind != "pascal") { std::fprintf(stderr, "unknown kind\n"); return 2; }
        if (argc < 4) { std::fprintf(stderr, "see header\n"); return 2; }
        const int M = std::atoi(argv[3]);
        std::vector<int> ns;
        if (mode == "named") for (int n = 1; n <= M; ++n) ns.push_back(n);
        else for (int k = 1; k <= M; ++k) ns.push_back(1 << k);
        std::printf("# %s graphs: n edges diam wiener [alpha]\n", kind.c_str());
        for (int n : ns) {
            if (n <= 64 && !getenv("BELL_BIG")) {
                const Small G = small_bell(n, small_g(kind, n));
                u64 m = 0; for (int v = 0; v < n; ++v) m += std::popcount(G.adj[v]);
                const DistStats d = small_dist(G);
                std::printf("%d %llu %d %llu %d\n", n, (unsigned long long)(m / 2), d.diam, (unsigned long long)d.wiener, alpha(G));
            } else {
                const BigGraph G = bell_graph(n, make_g(kind, n - 1));
                u64 m = 0; for (u64 w : G.a) m += std::popcount(w);
                const DistStats d = distances(G);
                std::printf("%d %llu %d %llu\n", n, (unsigned long long)(m / 2), d.diam, (unsigned long long)d.wiener);
            }
            std::fflush(stdout);
        }
        return 0;
    }
    if (mode == "all") {
        const int n = std::atoi(argv[2]);
        int T = int(std::max(1u, std::thread::hardware_concurrency()));
        bool no_alpha = false;
        int part = 0, parts = 1;  // --part I P: only the I-th of P equal slices of the g range
        for (int i = 3; i < argc; ++i) {
            if (!std::strcmp(argv[i], "--threads") && i + 1 < argc) T = std::atoi(argv[++i]);
            else if (!std::strcmp(argv[i], "--no-alpha")) no_alpha = true;
            else if (!std::strcmp(argv[i], "--part") && i + 2 < argc) { part = std::atoi(argv[i + 1]); parts = std::atoi(argv[i + 2]); i += 2; }
        }
        const int free_bits = (n - 2) / 2;  // g_2, g_4, ..., g_(2*free_bits) with 2m <= n-2
        const u64 all_g = u64{1} << std::max(0, free_bits);
        const u64 first = all_g / parts * part, total = (part + 1 == parts) ? all_g : all_g / parts * (part + 1);
        const int K = std::has_single_bit(unsigned(n)) ? std::countr_zero(unsigned(n)) : -1;
        std::mutex mu;
        std::map<int, u64> diam_hist, alpha_hist;
        u64 wmin = ~u64{0}, wmax = 0;
        std::vector<u64> diamK;
        std::atomic<u64> next{first};
        const u64 catalan_g = small_g("catalan", n);
        auto work = [&]() {
            std::map<int, u64> dh, ah;
            u64 lo_w = ~u64{0}, hi_w = 0;
            std::vector<u64> dk;
            for (;;) {
                const u64 lo = next.fetch_add(4096);
                if (lo >= total) break;
                const u64 hi = std::min(total, lo + 4096);
                for (u64 e = lo; e < hi; ++e) {
                    const u64 g = iodec_g(n, e);
                    const Small G = small_bell(n, g);
                    const DistStats d = small_dist(G);
                    ++dh[d.diam];
                    lo_w = std::min(lo_w, d.wiener); hi_w = std::max(hi_w, d.wiener);
                    if (!no_alpha) ++ah[alpha(G)];
                    if (K >= 0 && d.diam == K) dk.push_back(g);
                }
            }
            std::lock_guard<std::mutex> lk(mu);
            for (auto& [k, v] : dh) diam_hist[k] += v;
            for (auto& [k, v] : ah) alpha_hist[k] += v;
            wmin = std::min(wmin, lo_w); wmax = std::max(wmax, hi_w);
            diamK.insert(diamK.end(), dk.begin(), dk.end());
        };
        std::vector<std::thread> pool;
        for (int t = 1; t < T; ++t) pool.emplace_back(work);
        work();
        for (auto& th : pool) th.join();
        std::printf("n=%d graphs=%llu diam:", n, (unsigned long long)(total - first));
        for (auto& [k, v] : diam_hist) std::printf(" %d:%llu", k, (unsigned long long)v);
        std::printf(" wiener:[%llu,%llu]", (unsigned long long)wmin, (unsigned long long)wmax);
        if (!no_alpha) { std::printf(" alpha:"); for (auto& [k, v] : alpha_hist) std::printf(" %d:%llu", k, (unsigned long long)v); }
        std::printf("\n");
        if (K >= 0) {
            std::sort(diamK.begin(), diamK.end());
            std::printf("  graphs with diam = log2 n = %d: %zu", K, diamK.size());
            for (u64 g : diamK) std::printf(" g=%llx%s", (unsigned long long)g, g == catalan_g ? "(Catalan)" : "");
            std::printf("\n");
        }
        return 0;
    }
    return 2;
}
