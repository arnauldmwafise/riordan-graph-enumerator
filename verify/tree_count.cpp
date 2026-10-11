/**
 * @file   verify/tree_count.cpp
 * @brief  tree_count -- method B: u(n), the number of unlabeled Riordan graphs, by a generator
 *         written independently of include/rgraph (method A, tools/rgcount.cpp).
 *
 * Definition used: G_n(g, f) has vertices 1..n; for i > j, ij is an edge iff
 * [t^(i-2)] g f^(j-1) = 1 over Z_2, with f(0) = 0. There are L(n) = (4^(n-1)+2)/3 labeled
 * Riordan graphs of order n.
 *
 * Generation: labeled Riordan graphs are built vertex by vertex. Placing vertex m (0-based)
 * fixes the coefficient g_(m-1) (= adjacency to vertex 0) and, once g != 0 with lowest term
 * t^r, the coefficient f_(m-1-r); all other entries of the new row are forced:
 * adjacency(m, j) = [t^(m-1)] g f^j for 1 <= j < m. Every choice gives a different row, so the
 * leaves of the tree are the labeled Riordan graphs, each exactly once. Each node has at most 4
 * children.
 *
 * Isomorphism classes: nauty canonical form with the unit partition (no initial colouring),
 * FNV/murmur-style 128-bit hash (different from method A's partition and hash). Work split:
 * invariant = sorted multiset of (degree, triangles at the vertex); pass p keeps only leaves whose
 * invariant hash is p mod P. Threads claim subtrees at depth --split with atomic flags
 * (deterministic DFS order): every thread walks the tree down to depth D and descends only into
 * the subtrees it claims first.
 *
 * Usage
 * -----
 *   tree_count N [--threads T] [--passes P] [--checkpoint FILE] [--split D]
 *
 *  N                  order, 1 <= N <= 31 (otherwise exit code 2, no message)
 *  --threads T        worker threads (default: std::thread::hardware_concurrency(), at least 1)
 *  --passes P         number of passes (default 1)
 *  --checkpoint FILE  after each pass the line "N P pass classes labeled" (pass 0-based) is
 *                     appended; a rerun with the same N and P skips recorded passes
 *  --split D          depth at which subtrees are distributed to threads (default 10, clamped
 *                     to 1..N-1)
 * Options are read in (name, value) pairs; unknown options and a trailing option without value
 * are ignored.
 *
 * Output
 * ------
 *  stdout  one line "N L u": order, number of labeled graphs generated (leaves), number of
 *          classes, e.g. "10 87382 36897". For N = 1 it prints "1 1 1" directly.
 *  stderr  " pass p/P: c classes" per pass run; error messages.
 *
 * Built-in check: every pass generates all labeled graphs; the tool stops (exit code 1) if the
 * number of leaves differs between passes or from a checkpoint ("passes disagree on labeled
 * count", "checkpoint disagrees"). The printed L is compared with (4^(N-1)+2)/3 by the caller.
 *
 * Exit codes: 0 success; 1 inconsistent labeled counts; 2 usage error.
 *
 * Cost: about L(N) leaves per pass (generation is repeated in every pass) and L(N) nauty
 * canonicalisations in total; memory about 16 bytes per stored class fingerprint.
 *
 * Role in the verification strategy: method B, independent of method A (shares no code, grows a
 * search tree instead of enumerating parameters, uses a different nauty partition, hash and
 * work-splitting invariant). Its u(N) must agree with method A, and for N <= 10 with method C
 * (geng | verify/recognize_all); see tests/small_counts.sh.
 */
#include "nauty.h"

#include <algorithm>
#include <atomic>
#include <bit>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <functional>
#include <memory>
#include <string>
#include <thread>
#include <vector>

/// Z_2 polynomial (bit k = coefficient of t^k) or vertex bitmask (bit v = vertex v, 0-based).
using u64 = std::uint64_t;

namespace {

/** @brief 128-bit class fingerprint, ordered for sort + unique. */
struct H {
    u64 a, b;  ///< the two halves
    bool operator<(const H& o) const { return a != o.a ? a < o.a : b < o.b; }
    bool operator==(const H& o) const { return a == o.a && b == o.b; }
};
/** @brief murmur3 64-bit finaliser (bijective mixing). */
inline u64 fm(u64 k) {
    k ^= k >> 33; k *= 0xff51afd7ed558ccdULL;
    k ^= k >> 33; k *= 0xc4ceb9fe1a85ec53ULL;
    return k ^ (k >> 33);
}

int N, P = 1, SPLIT = 10;  ///< order; number of passes; subtree split depth (--split)
std::size_t T = 1;         ///< number of worker threads

/**
 * @brief Product of two Z_2 polynomials modulo t^bits.
 * @param a, b  factors (bit k = coefficient of t^k)
 * @param bits  number of low coefficients kept (>= 64 keeps all)
 */
u64 polymul(u64 a, u64 b, int bits) {
    u64 r = 0;
    for (; b; b &= b - 1) r ^= a << std::countr_zero(b);
    return bits >= 64 ? r : (r & ((u64{1} << bits) - 1));
}

/** @brief State of a node of the generation tree: the rows placed so far and the known coefficients. */
struct Node {
    u64 row[64];  ///< row[m] = adjacency of vertex m to vertices < m
    u64 g = 0, f = 0;  ///< known coefficients g_0..g_(m-1) and f_1..f_kf
    int kf = 0;        ///< highest known f index (0: none known yet)
};

/** @brief Per-thread results of a pass. */
struct Worker {
    std::vector<H> mem;  ///< fingerprints of the leaves in this pass (deduplicated when full)
    u64 leaves = 0;      ///< all leaves (labeled graphs) visited, whatever their pass
    /** @brief Append h; when at capacity, sort + unique first and grow by max(2^20, size/4). */
    void push(const H& h) {
        if (mem.size() == mem.capacity()) {
            std::sort(mem.begin(), mem.end());
            mem.erase(std::unique(mem.begin(), mem.end()), mem.end());
            std::vector<H> nx;
            nx.reserve(mem.size() + std::max<std::size_t>(std::size_t{1} << 20, mem.size() / 4));
            nx.assign(mem.begin(), mem.end());
            mem.swap(nx);
        }
        mem.push_back(h);
    }
};

/**
 * @brief Isomorphism invariant used to split classes into passes: hash of the sorted multiset of
 *        (degree, number of triangles at the vertex).
 * @param adj  symmetric adjacency rows of a graph on N vertices
 */
u64 split_invariant(const u64* adj) {
    u64 key[64];
    for (int v = 0; v < N; ++v) {
        int tri = 0;
        for (u64 x = adj[v]; x; x &= x - 1) tri += std::popcount(adj[v] & adj[std::countr_zero(x)]);
        key[v] = (u64(std::popcount(adj[v])) << 32) | u64(tri / 2);
    }
    std::sort(key, key + N);
    u64 h = 0x6a09e667f3bcc908ULL;
    for (int v = 0; v < N; ++v) h = fm(h + key[v] * 0x9e3779b97f4a7c15ULL + v);
    return h;
}

/**
 * @brief 128-bit fingerprint of the isomorphism class: nauty canonical form (unit partition),
 *        hashed together with N. Thread-safe with a USE_TLS build of nauty.
 * @param adj  symmetric adjacency rows of a graph on N vertices
 */
H leaf_hash(const u64* adj) {
    static thread_local DEFAULTOPTIONS_GRAPH(opts);
    graph g[MAXN], cg[MAXN];
    int lab[MAXN], ptn[MAXN], orb[MAXN];
    for (int v = 0; v < N; ++v) {
        setword s = 0;
        for (u64 x = adj[v]; x; x &= x - 1) s |= bit[std::countr_zero(x)];
        g[v] = s;
    }
    optionblk o = opts;
    o.getcanon = TRUE;  // defaultptn stays TRUE: unit partition
    statsblk st;
    densenauty(g, lab, ptn, orb, &o, &st, 1, N, cg);
    u64 a = 0xcbf29ce484222325ULL, b = 0x84222325cbf29ce4ULL;
    auto eat = [&](u64 w) {
        for (int s = 0; s < 64; s += 8) a = (a ^ ((w >> s) & 0xff)) * 0x100000001b3ULL;
        b = fm(b ^ w) + 0x9e3779b97f4a7c15ULL;
    };
    eat(u64(N));
    for (int v = 0; v < N; ++v) eat(u64(cg[v]));
    return {fm(a), fm(b ^ a)};
}

/**
 * @brief Process a leaf (a complete labeled Riordan graph): count it, and if it belongs to this
 *        pass, store its fingerprint.
 * @param nd    node with all N rows placed
 * @param w     worker receiving the result
 * @param pass  0-based pass index
 */
void leaf(Node& nd, Worker& w, int pass) {
    ++w.leaves;
    u64 adj[64] = {};
    for (int m = 0; m < N; ++m)
        for (u64 x = nd.row[m]; x; x &= x - 1) {
            const int j = std::countr_zero(x);
            adj[m] |= u64{1} << j;
            adj[j] |= u64{1} << m;
        }
    if (P > 1 && split_invariant(adj) % u64(P) != u64(pass)) return;
    w.push(leaf_hash(adj));
}

/**
 * @brief All admissible rows for vertex m; calls visit(row, g, f, kf) for each (at most 4).
 *
 * row bit 0 is the new coefficient g_(m-1); bit j (1 <= j < m) is [t^(m-1)] g f^j. If the new
 * row needs the not yet known coefficient f_(m-1-r), both values are tried.
 * @param nd     current node (rows 0..m-1 placed)
 * @param m      vertex to place, 1 <= m < N
 * @param visit  callback receiving the row and the updated g, f, kf
 */
// all admissible rows for vertex m; calls visit(row, g, f, kf)
template <class F>
void children(const Node& nd, int m, F&& visit) {
    for (int b = 0; b < 2; ++b) {
        const u64 g = nd.g | (u64(b) << (m - 1));
        auto row_for = [&](u64 f) {
            u64 row = u64(b), col = g;
            for (int j = 1; j < m; ++j) {
                col = polymul(col, f, m);
                row |= ((col >> (m - 1)) & 1) << j;
            }
            return row;
        };
        if (g == 0) { visit(u64(0), g, nd.f, nd.kf); continue; }
        const int need = m - 1 - std::countr_zero(g);
        if (need <= nd.kf) { visit(row_for(nd.f), g, nd.f, nd.kf); continue; }
        for (int fb = 0; fb < 2; ++fb) {
            const u64 f = nd.f | (u64(fb) << need);
            visit(row_for(f), g, f, need);
        }
    }
}

/** @brief Shared state of one pass. */
struct Run {
    int pass;                                      ///< 0-based pass index
    std::unique_ptr<std::atomic<bool>[]> claimed;  ///< claimed[k]: subtree k at depth SPLIT taken by a thread
    std::size_t nclaim;                            ///< number of subtrees at depth SPLIT
};

/**
 * @brief Depth-first traversal from vertex m. At depth SPLIT an unowned walk claims the next
 *        subtree (numbered in DFS order via idx) or skips it if another thread claimed it.
 * @param nd     current node (rows overwritten during the walk; g, f, kf restored afterwards)
 * @param m      next vertex to place
 * @param w      worker of this thread
 * @param run    shared pass state
 * @param idx    running index of depth-SPLIT nodes seen by this thread
 * @param owned  true inside a claimed subtree; leaves are processed only when owned
 */
void dfs(Node& nd, int m, Worker& w, Run& run, std::size_t& idx, bool owned) {
    if (m == N) { if (owned) leaf(nd, w, run.pass); return; }
    if (m == SPLIT && !owned) {
        const std::size_t k = idx++;
        if (run.claimed[k].exchange(true)) return;  // another thread does this subtree
        owned = true;
    }
    children(nd, m, [&](u64 row, u64 g, u64 f, int kf) {
        const u64 g0 = nd.g, f0 = nd.f;
        const int kf0 = nd.kf;
        nd.row[m] = row; nd.g = g; nd.f = f; nd.kf = kf;
        dfs(nd, m + 1, w, run, idx, owned);
        nd.g = g0; nd.f = f0; nd.kf = kf0;
    });
}

/** @brief Number of nodes of the generation tree at depth SPLIT (= number of claimable subtrees). */
std::size_t count_split_nodes() {  // number of nodes at depth SPLIT
    std::size_t k = 0;
    std::function<void(Node&, int)> go = [&](Node& nd, int m) {
        if (m == SPLIT) { ++k; return; }
        children(nd, m, [&](u64 row, u64 g, u64 f, int kf) {
            Node c = nd;
            c.row[m] = row; c.g = g; c.f = f; c.kf = kf;
            go(c, m + 1);
        });
    };
    Node root;
    std::memset(root.row, 0, sizeof root.row);
    go(root, 1);
    return k;
}

}  // namespace

/**
 * @brief Parse options, resume from the checkpoint, run the remaining passes and print "N L u".
 * @return exit code (see file header)
 */
int main(int argc, char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: tree_count N [--threads T] [--passes P] [--checkpoint F] [--split D]\n"); return 2; }
    N = std::atoi(argv[1]);
    T = std::max(1u, std::thread::hardware_concurrency());
    std::string ckpt;
    for (int i = 2; i + 1 < argc; i += 2) {
        const std::string a = argv[i];
        if (a == "--threads") T = std::max(1, std::atoi(argv[i + 1]));
        else if (a == "--passes") P = std::max(1, std::atoi(argv[i + 1]));
        else if (a == "--checkpoint") ckpt = argv[i + 1];
        else if (a == "--split") SPLIT = std::max(1, std::atoi(argv[i + 1]));
    }
    if (N < 1 || N > 31) return 2;
    nauty_check(WORDSIZE, 1, N, NAUTYVERSIONID);
    if (N == 1) { std::printf("1 1 1\n"); return 0; }
    SPLIT = std::max(1, std::min(SPLIT, N - 1));
    std::vector<char> done(P, 0);
    u64 classes = 0, labeled = 0;
    if (!ckpt.empty() && std::filesystem::exists(ckpt)) {
        FILE* f = std::fopen(ckpt.c_str(), "r");
        unsigned long long cn, cp, cq, cr, cl;
        while (f && std::fscanf(f, "%llu %llu %llu %llu %llu", &cn, &cp, &cq, &cr, &cl) == 5)
            if (int(cn) == N && int(cp) == P && int(cq) < P && !done[cq]) {
                if (labeled && labeled != cl) { std::fprintf(stderr, "checkpoint disagrees\n"); return 1; }
                done[cq] = 1; classes += cr; labeled = cl;
            }
        if (f) std::fclose(f);
    }
    const std::size_t K = count_split_nodes();
    for (int pass = 0; pass < P; ++pass) {
        if (done[pass]) continue;
        Run run{pass, std::unique_ptr<std::atomic<bool>[]>(new std::atomic<bool>[K]), K};
        for (std::size_t k = 0; k < K; ++k) run.claimed[k] = false;
        std::vector<Worker> ws(T);
        auto work = [&](std::size_t t) {
            Node root;
            std::memset(root.row, 0, sizeof root.row);
            std::size_t idx = 0;
            dfs(root, 1, ws[t], run, idx, false);
        };
        std::vector<std::thread> pool;
        for (std::size_t t = 1; t < T; ++t) pool.emplace_back(work, t);
        work(0);
        for (auto& th : pool) th.join();
        u64 lab = 0;
        std::vector<H> all;
        for (auto& w : ws) { lab += w.leaves; all.insert(all.end(), w.mem.begin(), w.mem.end()); std::vector<H>().swap(w.mem); }
        std::sort(all.begin(), all.end());
        const u64 rp = u64(std::unique(all.begin(), all.end()) - all.begin());
        if (labeled && labeled != lab) { std::fprintf(stderr, "passes disagree on labeled count\n"); return 1; }
        labeled = lab; classes += rp;
        std::fprintf(stderr, " pass %d/%d: %llu classes\n", pass + 1, P, (unsigned long long)rp);
        if (!ckpt.empty()) {
            FILE* f = std::fopen(ckpt.c_str(), "a");
            std::fprintf(f, "%d %d %d %llu %llu\n", N, P, pass, (unsigned long long)rp, (unsigned long long)lab);
            std::fclose(f);
        }
    }
    std::printf("%d %llu %llu\n", N, (unsigned long long)labeled, (unsigned long long)classes);
    return 0;
}
