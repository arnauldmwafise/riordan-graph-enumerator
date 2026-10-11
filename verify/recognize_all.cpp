/**
 * @file   verify/recognize_all.cpp
 * @brief  recognize_all -- method C: decide, for every graph read from stdin (graph6, e.g. from
 *         nauty's geng), whether it is an (unlabeled) Riordan graph, i.e. whether some labeling of
 *         its vertices is a labeled Riordan graph G_n(g, f).
 *
 * Independent of the generators of methods A and B: no Riordan graphs are generated and no
 * canonical forms are used. Self-contained (does not include rgraph/ or nauty).
 *
 * Labeled Riordan graphs (Cheon et al., LAA 579 (2019), Def. 2.4): for positions i > j
 * (1-based) the pair is adjacent iff [t^(i-2)] g f^(j-1) = 1 mod 2, f(0) = 0.
 *
 * Search: vertices are placed one position at a time. With 0-based positions, placing position
 * m >= 1 fixes g_(m-1) = adj(w, v_0); for j >= 1 the entry to position j must equal
 * [t^(m-1)] g f^j. These coefficients use only g_0..g_(m-1) and f_1..f_(m-1-r) (r = lowest
 * index with g_r = 1); at most one f bit, f_(m-1-r), is new at this step and is tried both ways.
 * If g_0..g_(m-1) are all 0, every entry to positions >= 1 must be 0.
 * A complete placement is a labeling whose adjacency matrix is exactly G_n(g, f).
 * Twin vertices (equal neighbourhoods apart from each other) are tried only once per position,
 * and only one vertex per twin class is tried at position 0. The search stops at the first
 * Riordan labeling found. Cost: exponential in the worst case; the mean number of search nodes
 * per graph is about 6000 for n = 10 (see --stats and docs/COMPUTATIONS.md).
 *
 * Usage
 * -----
 *   geng -q N | recognize_all [--nonriordan FILE] [--riordan FILE] [--stats]
 *
 *  stdin               graph6 lines (short form, n <= 62, lines < 4096 bytes); lines starting
 *                      with '>' (headers) and lines that are not graph6 are skipped. Input lines
 *                      are assumed to be well-formed graph6 (no length check).
 *  --nonriordan FILE   write the input lines of the non-Riordan graphs to FILE (verbatim)
 *  --riordan FILE      write the input lines of the Riordan graphs to FILE (verbatim)
 *  --stats             print search-node statistics to stderr
 * Unknown arguments are ignored; if a FILE cannot be opened, that list is silently not written.
 *
 * Output
 * ------
 *  stdout  one line "n total riordan nonriordan", where n is the order of the last graph read
 *          (-1 if none), e.g. "7 1044 406 638" for geng -q 7.
 *  stderr  with --stats: "search nodes: mean X, max Y" (X with one decimal).
 * Exit code: always 0.
 *
 * Role in the verification strategy: method C. Run on all graphs of order n from geng (one per
 * isomorphism class), "riordan" is u(n) and must agree with methods A (rgcount) and B
 * (tree_count); "total" must equal the number of graphs of order n (OEIS A000088). Used for
 * n <= 10 (scripts/method_c.sh, tests/small_counts.sh).
 */
#include <bit>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

/// Z_2 polynomial (bit k = coefficient of t^k) or vertex/position bitmask.
using u64 = std::uint64_t;

namespace {

int n_g;       ///< order of the current input graph (n <= 62)
u64 adjm[64];  ///< adjacency rows of the current input graph (bit w of adjm[v]: vw is an edge)
u64 nodes;  ///< search nodes visited (statistics)

// GF(2) polynomial helpers on 64-bit words (bit k = coefficient of t^k)
/**
 * @brief Product of two Z_2 polynomials modulo t^bits.
 * @param a, b  factors
 * @param bits  number of low coefficients kept (>= 64 keeps all)
 */
u64 pmul(u64 a, u64 b, int bits) {
    u64 r = 0;
    for (; b; b &= b - 1) r ^= a << std::countr_zero(b);
    return bits >= 64 ? r : (r & ((u64{1} << bits) - 1));
}

// Twin classes: u, w are twins if their neighbourhoods agree outside {u, w}; swapping them is
// an automorphism fixing all other vertices, so only one representative per class is tried.
int twin_rep[64];  ///< twin_rep[u]: representative (smallest index found) of the twin class of u

/** @brief Search state: the partial labeling and the coefficients of g and f known so far. */
struct State {
    int pos[64];     ///< position -> vertex
    u64 posmask[64]; ///< for each vertex: bit j set iff it is adjacent to the vertex at position j
    u64 used = 0;    ///< vertices placed
    u64 g = 0;       ///< known coefficients g_0..g_(m-1)
    u64 f = 0;       ///< known coefficients f_1..f_kf
    int kf = 0;      ///< highest known f index
};

bool place(State& s, int m);  // forward declaration

/**
 * @brief Place vertex w at position m and update posmask of its neighbours.
 * @param s  state (modified)
 * @param m  position
 * @param w  unused vertex
 */
void put(State& s, int m, int w) {
    s.pos[m] = w;
    s.used |= u64{1} << w;
    for (u64 x = adjm[w]; x; x &= x - 1) s.posmask[std::countr_zero(x)] |= u64{1} << m;
}

/** @brief Undo put(s, m, w) (pos[m] is left stale; it is overwritten before it is read again). */
void unput(State& s, int m, int w) {
    s.used &= ~(u64{1} << w);
    for (u64 x = adjm[w]; x; x &= x - 1) s.posmask[std::countr_zero(x)] &= ~(u64{1} << m);
}

/**
 * @brief Predicted adjacency of the new position m to positions 1..m-1.
 * @param g  coefficients g_0..g_(m-1)
 * @param f  known coefficients of f (f(0) = 0)
 * @param m  new position, m >= 1
 * @return   position bitmask, bit j set iff [t^(m-1)] g f^j = 1 (1 <= j < m)
 */
// Predicted adjacency of the new position m to positions 1..m-1 (as a position bitmask),
// for given g (bits 0..m-1) and f.
u64 predicted(u64 g, u64 f, int m) {
    u64 pat = 0, col = g;
    for (int j = 1; j < m; ++j) {
        col = pmul(col, f, m);
        pat |= ((col >> (m - 1)) & 1) << j;
    }
    return pat;
}

/**
 * @brief Recursive backtracking: try every unused vertex (one per twin class) at position m
 *        whose adjacency to positions 0..m-1 matches one of the at most 4 admissible patterns.
 * @param s  state with positions 0..m-1 filled (restored on failure)
 * @param m  position to fill
 * @return   true iff the placement can be completed to a Riordan labeling; on success the state
 *           is left at the complete labeling. Increments the global counter `nodes`.
 */
bool place(State& s, int m) {
    ++nodes;
    if (m == n_g) return true;
    const u64 posmask_mid = (m >= 2) ? (((u64{1} << m) - 1) & ~u64{1}) : 0;  // positions 1..m-1
    // the up to 4 admissible rows: choice of the new g bit b, and of the new f bit if one is needed
    struct Opt { u64 g, f; int kf; u64 pat; } opts[4];
    int nopt = 0;
    for (int b = 0; b < 2; ++b) {
        const u64 g = s.g | (u64(b) << (m - 1));
        if (g == 0) { opts[nopt++] = {g, s.f, s.kf, 0}; continue; }
        const int need = m - 1 - std::countr_zero(g);
        if (need <= s.kf) { opts[nopt++] = {g, s.f, s.kf, predicted(g, s.f, m)}; continue; }
        for (int fb = 0; fb < 2; ++fb) {
            const u64 f = s.f | (u64(fb) << need);
            opts[nopt++] = {g, f, need, predicted(g, f, m)};
        }
    }
    u64 tried = 0;  // twin classes already tried at this position
    for (int w = 0; w < n_g; ++w) {
        if (s.used >> w & 1) continue;
        if (tried >> twin_rep[w] & 1) continue;
        const u64 b = s.posmask[w] & 1;  // adjacency to position 0 = new g bit
        const u64 row = s.posmask[w] & posmask_mid;
        bool any = false;
        for (int k = 0; k < nopt; ++k) {
            if (((opts[k].g >> (m - 1)) & 1) != b || opts[k].pat != row) continue;
            any = true;
            const u64 g0 = s.g, f0 = s.f;
            const int kf0 = s.kf;
            put(s, m, w);
            s.g = opts[k].g; s.f = opts[k].f; s.kf = opts[k].kf;
            if (place(s, m + 1)) return true;
            unput(s, m, w);  // undo in place instead of copying the state
            s.g = g0; s.f = f0; s.kf = kf0;
        }
        if (any) tried |= u64{1} << twin_rep[w];
    }
    return false;
}

/** @brief Fill twin_rep[] for the current graph (O(n^2) word operations). */
void compute_twins() {
    for (int u = 0; u < n_g; ++u) {
        twin_rep[u] = u;
        for (int w = 0; w < u; ++w) {
            const u64 mask = ~((u64{1} << u) | (u64{1} << w));
            if ((adjm[u] & mask) == (adjm[w] & mask)) { twin_rep[u] = twin_rep[w]; break; }
        }
    }
}

/**
 * @brief True iff the current graph (n_g, adjm) has a Riordan labeling. Tries each twin-class
 *        representative at position 0 and calls place() for the rest.
 */
bool is_riordan_unlabeled() {
    if (n_g <= 1) return true;
    compute_twins();
    u64 tried = 0;
    for (int v0 = 0; v0 < n_g; ++v0) {
        if (tried >> twin_rep[v0] & 1) continue;
        tried |= u64{1} << twin_rep[v0];
        State s;
        std::memset(s.posmask, 0, sizeof s.posmask);
        put(s, 0, v0);
        if (place(s, 1)) return true;
    }
    return false;
}

/**
 * @brief Decode a short-form graph6 line into n_g and adjm.
 * @param line  NUL-terminated text; must be valid graph6 (data bytes are not range-checked)
 * @return      false if the first character is not a valid order byte (63..126)
 */
bool read_graph6(const char* line) {
    const unsigned char* p = reinterpret_cast<const unsigned char*>(line);
    if (*p < 63 || *p > 126) return false;
    n_g = *p++ - 63;  // n <= 62 in the short form
    std::memset(adjm, 0, sizeof adjm);
    int k = 6, x = 0;
    for (int j = 1; j < n_g; ++j)
        for (int i = 0; i < j; ++i) {
            if (k == 6) { x = *p++ - 63; k = 0; }
            if ((x >> (5 - k)) & 1) { adjm[i] |= u64{1} << j; adjm[j] |= u64{1} << i; }
            ++k;
        }
    return true;
}

}  // namespace

/**
 * @brief Read graph6 lines from stdin, classify each one and print the summary line.
 * @return 0
 */
int main(int argc, char** argv) {
    FILE* out = nullptr; FILE* rout = nullptr;
    bool stats = false;
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--nonriordan") && i + 1 < argc) out = std::fopen(argv[++i], "w");
        else if (!std::strcmp(argv[i], "--riordan") && i + 1 < argc) rout = std::fopen(argv[++i], "w");
        else if (!std::strcmp(argv[i], "--stats")) stats = true;
    }
    char buf[4096];
    u64 total = 0, riordan = 0, max_nodes = 0, sum_nodes = 0;
    int n_seen = -1;
    while (std::fgets(buf, sizeof buf, stdin)) {
        if (buf[0] == '>' || !read_graph6(buf)) continue;
        n_seen = n_g;
        ++total;
        nodes = 0;
        const bool ok = is_riordan_unlabeled();
        sum_nodes += nodes;
        if (nodes > max_nodes) max_nodes = nodes;
        if (ok) { ++riordan; if (rout) std::fputs(buf, rout); }
        else if (out) std::fputs(buf, out);
    }
    if (out) std::fclose(out);
    if (rout) std::fclose(rout);
    std::printf("%d %llu %llu %llu\n", n_seen, (unsigned long long)total, (unsigned long long)riordan,
                (unsigned long long)(total - riordan));
    if (stats)
        std::fprintf(stderr, "search nodes: mean %.1f, max %llu\n", total ? double(sum_nodes) / total : 0.0,
                     (unsigned long long)max_nodes);
    return 0;
}
