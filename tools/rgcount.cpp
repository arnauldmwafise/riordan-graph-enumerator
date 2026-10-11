/**
 * @file   tools/rgcount.cpp
 * @brief  rgcount -- method A: the number u(n) of non-isomorphic (unlabeled) Riordan graphs of
 *         order N, by generating every labeled Riordan graph and counting isomorphism classes.
 *
 * Method
 * ------
 * A labeled Riordan graph G_n(g, f) has vertices 1..n; for i > j, ij is an edge iff
 * [t^(i-2)] g f^(j-1) = 1 over Z_2, with f(0) = 0. Every labeled Riordan graph is generated
 * exactly once by rgraph::Enumeration (the enumeration of the proof of Prop. 2.8 of
 * Cheon-Jung-Kitaev-Mojallal 2019); there are L(N) = (4^(N-1)+2)/3 of them. Each is reduced to a
 * 128-bit fingerprint of its isomorphism class (rgraph::class_hash: nauty canonical form with the
 * initial partition by degree), and distinct fingerprints are counted (sort + unique).
 *
 * Pass mode: pass p handles only graphs whose degree-sequence hash (rgraph::degree_invariant) is
 * p mod P; isomorphic graphs have the same degree sequence, so every class is counted in exactly
 * one pass and every graph is canonicalised exactly once overall (only the cheap generation is
 * repeated). Memory per pass is about 16 bytes x 2.5 x u(N)/P.
 *
 * Threads take chunks of C consecutive enumeration indices from a shared atomic counter; each
 * thread keeps its own hash vector, compacted (sorted, deduplicated) whenever it is full.
 *
 * Usage
 * -----
 *   rgcount N [--threads T] [--passes P] [--checkpoint FILE] [--chunk C]
 *
 *  N                  order, 1 <= N <= 31
 *  --threads T        worker threads (default: std::thread::hardware_concurrency(), at least 1)
 *  --passes P         number of passes (default 1); see "Pass mode"
 *  --checkpoint FILE  after each pass the line "N P pass classes labeled" is appended (pass is
 *                     0-based); a rerun with the same N and P skips the recorded passes
 *  --chunk C          enumeration indices claimed per step by a thread (default 16384)
 * Options are read in (name, value) pairs; a trailing option without a value is ignored, and
 * non-numeric values abort with an uncaught std::invalid_argument.
 *
 * Output
 * ------
 *  stdout  one line "N L(N) u(N)  (S s)", e.g. "10 87382 36897  (0.1 s)" (S = elapsed seconds
 *          for the passes run now, one decimal; u(N) includes classes from resumed passes).
 *  stderr  a header "rgcount N=...: L labeled graphs, T thread(s), P pass(es)", one line
 *          " pass p/P: c classes (s s)" per pass, and " resuming: ..." for checkpointed passes.
 *
 * Built-in check: every pass counts all labeled graphs it generates; the tool stops if that
 * number differs from L(N) (hence also if it differs between passes), and rejects checkpoint
 * lines whose labeled count is not L(N).
 *
 * Exit codes: 0 success; 1 labeled-count mismatch, bad checkpoint or unwritable checkpoint file;
 * 2 usage error (no arguments, unknown option, N out of range).
 *
 * Cost: L(N) graph generations per pass plus L(N) nauty canonicalisations in total, about
 * 3-4 microseconds per labeled graph per core (see docs/COMPUTATIONS.md).
 *
 * Role in the verification strategy: method A. Its value of u(N) is compared with method B
 * (verify/tree_count.cpp, written independently and sharing no code) and, for N <= 10, with
 * method C (geng | verify/recognize_all.cpp); see tests/small_counts.sh.
 */
#include "rgraph/canon.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

using namespace rgraph;

namespace {

/** @brief Command-line configuration. */
struct Options {
    int N = 0;                            ///< order of the graphs
    std::size_t threads = 1, passes = 1;  ///< worker threads; number of passes P
    u64 chunk = 1 << 14;                  ///< enumeration indices claimed per step
    std::string checkpoint;               ///< checkpoint file, empty = none
};

/**
 * @brief Sort v and remove duplicates in place.
 * @param v  vector of fingerprints
 */
void compact(std::vector<Hash128>& v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
}

/**
 * @brief Append h to v; when v is at capacity, first compact it and then grow the capacity by
 *        max(2^20, size/4) elements, which bounds memory to roughly the number of distinct values.
 * @param v  per-thread fingerprint vector (may contain duplicates between compactions)
 * @param h  fingerprint to add
 */
void push_bounded(std::vector<Hash128>& v, const Hash128& h) {
    if (v.size() == v.capacity()) {  // compact, then grow by a bounded step
        compact(v);
        std::vector<Hash128> next;
        next.reserve(v.size() + std::max<std::size_t>(std::size_t{1} << 20, v.size() / 4));
        next.assign(v.begin(), v.end());
        v.swap(next);
    }
    v.push_back(h);
}

/** @brief Outcome of one pass. */
struct PassResult {
    u64 classes = 0, labeled = 0;  ///< distinct fingerprints in this pass; labeled graphs generated
};

/**
 * @brief Run pass number `pass`: generate all L(N) labeled graphs with o.threads threads and count
 *        the distinct class fingerprints of those whose degree_invariant() is pass mod o.passes.
 * @param o     configuration
 * @param pass  0-based pass index, < o.passes
 * @return      number of classes found in this pass and number of labeled graphs generated
 */
PassResult run_pass(const Options& o, std::size_t pass) {
    const Enumeration E(o.N);
    const u64 total = E.size();
    std::atomic<u64> next{0};
    std::vector<std::vector<Hash128>> mem(o.threads);
    std::vector<u64> labeled(o.threads, 0);
    auto work = [&](std::size_t t) {
        auto& v = mem[t];
        for (;;) {
            const u64 lo = next.fetch_add(o.chunk);
            if (lo >= total) break;
            const u64 hi = std::min(total, lo + o.chunk);
            for (u64 i = lo; i < hi; ++i) {
                const Graph G = E.graph(i);
                ++labeled[t];
                if (o.passes > 1 && degree_invariant(G) % o.passes != pass) continue;
                push_bounded(v, class_hash(G));
            }
        }
    };
    std::vector<std::thread> pool;
    for (std::size_t t = 1; t < o.threads; ++t) pool.emplace_back(work, t);
    work(0);
    for (auto& th : pool) th.join();
    std::vector<Hash128> all;
    std::size_t sz = 0;
    for (auto& v : mem) { compact(v); sz += v.size(); }
    all.reserve(sz);
    for (auto& v : mem) { all.insert(all.end(), v.begin(), v.end()); std::vector<Hash128>().swap(v); }
    compact(all);
    PassResult r;
    r.classes = all.size();
    for (u64 x : labeled) r.labeled += x;
    return r;
}

}  // namespace

/**
 * @brief Parse options, resume from the checkpoint, run the remaining passes and print the result.
 * @return exit code (see file header)
 */
int main(int argc, char** argv) {
    std::vector<std::string> a(argv + 1, argv + argc);
    if (a.empty()) {
        std::cerr << "usage: rgcount N [--threads T] [--passes P] [--checkpoint FILE] [--chunk C]\n";
        return 2;
    }
    Options o;
    o.N = std::stoi(a[0]);
    o.threads = std::max(1u, std::thread::hardware_concurrency());
    for (std::size_t i = 1; i + 1 < a.size(); i += 2) {
        if (a[i] == "--threads") o.threads = std::max<std::size_t>(1, std::stoul(a[i + 1]));
        else if (a[i] == "--passes") o.passes = std::max<std::size_t>(1, std::stoul(a[i + 1]));
        else if (a[i] == "--checkpoint") o.checkpoint = a[i + 1];
        else if (a[i] == "--chunk") o.chunk = std::max<u64>(1, std::stoull(a[i + 1]));
        else { std::cerr << "unknown option " << a[i] << "\n"; return 2; }
    }
    if (o.N < 1 || o.N > 31) { std::cerr << "rgcount: need 1 <= N <= 31\n"; return 2; }
    nauty_check(WORDSIZE, 1, o.N, NAUTYVERSIONID);

    const u64 expected = labeled_count(o.N);
    std::vector<char> done(o.passes, 0);
    u64 classes = 0;
    if (!o.checkpoint.empty() && std::filesystem::exists(o.checkpoint)) {
        FILE* f = std::fopen(o.checkpoint.c_str(), "r");
        unsigned long long cn, cp, cq, cr, cl;
        while (f && std::fscanf(f, "%llu %llu %llu %llu %llu", &cn, &cp, &cq, &cr, &cl) == 5)
            if (cn == static_cast<unsigned long long>(o.N) && cp == o.passes && cq < o.passes && !done[cq]) {
                if (cl != expected) { std::cerr << "checkpoint: wrong labeled count\n"; return 1; }
                done[cq] = 1;
                classes += cr;
                std::cerr << " resuming: pass " << cq + 1 << " already done (" << cr << " classes)\n";
            }
        if (f) std::fclose(f);
    }
    const auto t0 = std::chrono::steady_clock::now();
    std::cerr << "rgcount N=" << o.N << ": " << expected << " labeled graphs, " << o.threads << " thread(s), "
              << o.passes << " pass(es)\n";
    for (std::size_t p = 0; p < o.passes; ++p) {
        if (done[p]) continue;
        const PassResult r = run_pass(o, p);
        if (r.labeled != expected) {
            std::cerr << "rgcount: generated " << r.labeled << " labeled graphs, expected " << expected << "\n";
            return 1;
        }
        classes += r.classes;
        const double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::cerr << " pass " << p + 1 << "/" << o.passes << ": " << r.classes << " classes (" << s << " s)\n";
        if (!o.checkpoint.empty()) {
            FILE* f = std::fopen(o.checkpoint.c_str(), "a");
            if (!f) { std::cerr << "cannot write checkpoint\n"; return 1; }
            std::fprintf(f, "%d %zu %zu %llu %llu\n", o.N, o.passes, p, (unsigned long long)r.classes,
                         (unsigned long long)r.labeled);
            std::fclose(f);
        }
    }
    const double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("%d %llu %llu  (%.1f s)\n", o.N, (unsigned long long)expected, (unsigned long long)classes, s);
    return 0;
}
