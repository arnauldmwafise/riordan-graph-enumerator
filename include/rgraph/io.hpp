/**
 * @file   rgraph/io.hpp
 * @brief  Input/output helpers shared by the command-line tools: polynomials over Z_2, graph6,
 *         Riordan matrices, edge lists and SVG drawings.
 *
 * Conventions
 * -----------
 * - A polynomial over Z_2 is a 64-bit word `u64`; bit k is the coefficient of t^k.
 * - A labeled Riordan graph G_n(g, f) has vertices 1..n (0-based in code). For i > j, the pair ij
 *   is an edge iff [t^(i-2)] g f^(j-1) = 1 (Cheon, Jung, Kitaev, Mojallal, LAA 579 (2019),
 *   Def. 2.4).
 * - The Riordan matrix printed by the tools is the (n-1)x(n-1) lower-triangular matrix
 *   B[i][j] = [t^i] g f^j (0 <= j <= i <= n-2). The strictly lower triangle of the adjacency
 *   matrix is B shifted down by one row: A[i][j] = B[i-2][j-1] for i > j (1-based).
 *
 * All functions are header-only and have no side effects, except write_svg(), which writes a file.
 */
#pragma once
#include "riordan_graph.hpp"

#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace rgraph {

/**
 * @brief Human-readable form of a Z_2 polynomial, e.g. 0b1011 -> "1 + t + t^3".
 * @param p  coefficient bits (bit k = coefficient of t^k)
 * @return   "0" for the zero polynomial, otherwise the sum of its monomials in increasing degree
 */
inline std::string poly_str(u64 p) {
    if (!p) return "0";
    std::string s;
    for (int k = 0; k < 64; ++k)
        if (p >> k & 1) {
            if (!s.empty()) s += " + ";
            s += k == 0 ? "1" : k == 1 ? "t" : "t^" + std::to_string(k);
        }
    return s;
}

/**
 * @brief Same as poly_str() but with Unicode superscript exponents ("1 + t³"), for drawings.
 * @param p  coefficient bits
 * @return   UTF-8 string
 */
inline std::string poly_unicode(u64 p) {
    static const char* sup[] = {"⁰", "¹", "²", "³", "⁴",
                                "⁵", "⁶", "⁷", "⁸", "⁹"};
    if (!p) return "0";
    std::string s;
    for (int k = 0; k < 64; ++k)
        if (p >> k & 1) {
            if (!s.empty()) s += " + ";
            if (k == 0) { s += "1"; continue; }
            s += "t";
            if (k > 1) for (char c : std::to_string(k)) s += sup[c - '0'];
        }
    return s;
}

/**
 * @brief Unicode subscript digits for an integer (used for "G₅").
 * @param n  non-negative integer
 * @return   UTF-8 string of subscript digits
 */
inline std::string subscript_unicode(int n) {
    static const char* sub[] = {"₀", "₁", "₂", "₃", "₄",
                                "₅", "₆", "₇", "₈", "₉"};
    std::string s;
    for (char c : std::to_string(n)) s += sub[c - '0'];
    return s;
}

/**
 * @brief Parse a polynomial written like "1 + t + t^3" (spaces optional; "0" is zero).
 *
 * Repeated monomials cancel, since coefficients are taken mod 2.
 * @param in  text to parse; monomials are "1", "t" and "t^k" with 0 <= k <= 62
 * @return    coefficient bits
 * @throws std::invalid_argument on malformed input
 */
inline u64 parse_poly(const std::string& in) {
    u64 p = 0;
    std::string s;
    for (char c : in) if (c != ' ') s += c;
    if (s == "0" || s.empty()) return 0;
    std::size_t i = 0;
    while (i < s.size()) {
        std::size_t j = s.find('+', i);
        if (j == std::string::npos) j = s.size();
        const std::string term = s.substr(i, j - i);
        int k;
        if (term == "1") k = 0;
        else if (term == "t") k = 1;
        else if (term.rfind("t^", 0) == 0) k = std::stoi(term.substr(2));
        else throw std::invalid_argument("cannot parse term '" + term + "'");
        if (k < 0 || k > 62) throw std::invalid_argument("exponent out of range");
        p ^= u64{1} << k;
        i = j + 1;
    }
    return p;
}

/**
 * @brief Reduce (g, f) to the canonical parameters that determine G_n(g, f).
 *
 * Only g mod t^(n-1) matters, and, if g != 0 has lowest term t^r, only f mod t^(n-1-r).
 * After this call, every labeled Riordan graph has exactly one (g, f).
 * @param n  order of the graph
 * @param g  in/out: first generating function
 * @param f  in/out: second generating function (f(0) = 0)
 */
inline void canonical_params(int n, u64& g, u64& f) {
    g &= low_mask(n - 1);
    if (!g) { f = 0; return; }
    f &= low_mask(n - 1 - std::countr_zero(g));
}

/**
 * @brief The (n-1)x(n-1) lower-triangular Riordan matrix B[i][j] = [t^i] g f^j over Z_2.
 * @param n  order of the graph (n >= 1)
 * @param g, f  generating functions
 * @return   B as a vector of rows; entries above the diagonal are 0
 */
inline std::vector<std::vector<int>> riordan_matrix(int n, u64 g, u64 f) {
    const int m = n - 1;
    std::vector<std::vector<int>> B(m, std::vector<int>(m, 0));
    const u64 mask = low_mask(m);
    u64 col = g & mask;
    for (int j = 0; j < m; ++j) {
        for (int i = j; i < m; ++i) B[i][j] = int(col >> i & 1);
        col = mul_trunc(col, f, mask);
    }
    return B;
}

/**
 * @brief Edge list of a graph with 1-based vertex numbers, as pairs (a, b) with a < b, sorted by b
 *        within a and then by a.
 */
inline std::vector<std::pair<int, int>> edge_list(const Graph& G) {
    std::vector<std::pair<int, int>> e;
    for (int j = 0; j < G.n; ++j)
        for (int i = j + 1; i < G.n; ++i)
            if (G.has_edge(i, j)) e.push_back({j + 1, i + 1});
    return e;
}

/**
 * @brief graph6 encoding (B. McKay's format; n <= 62), readable by nauty, networkx, SageMath.
 * @param G  graph
 * @return   graph6 string without header or newline
 */
inline std::string graph6(const Graph& G) {
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
 * @brief Decode a graph6 string (short form, n <= 62).
 * @param s  graph6 text (trailing newline allowed)
 * @param G  output graph
 * @return   false if the string is not valid short-form graph6
 */
inline bool parse_graph6(const std::string& s, Graph& G) {
    if (s.empty() || s[0] < 63 || s[0] > 125) return false;
    G = Graph{};
    G.n = s[0] - 63;
    std::size_t p = 1;
    int k = 6, x = 0;
    for (int j = 1; j < G.n; ++j)
        for (int i = 0; i < j; ++i) {
            if (k == 6) {
                if (p >= s.size() || s[p] < 63 || s[p] > 126) return false;
                x = s[p++] - 63; k = 0;
            }
            if ((x >> (5 - k)) & 1) { G.adj[i] |= u64{1} << j; G.adj[j] |= u64{1} << i; }
            ++k;
        }
    return true;
}

/**
 * @brief Write an SVG drawing of a labeled Riordan graph: vertices clockwise in Riordan order on a
 *        circle (vertex 1 at the top), edges as straight lines, caption "G_n(g, f)".
 * @param path   output file
 * @param G      the graph G_n(g, f)
 * @param g, f   its parameters (used for the caption)
 * @param note   second caption line, e.g. "#7, 5 edges" (may be empty)
 */
inline void write_svg(const std::string& path, const Graph& G, u64 g, u64 f, const std::string& note) {
    const int n = G.n;
    const double W = 260, H = 290, cx = 130, cy = 130, R = n == 1 ? 0 : 95, PI = 3.14159265358979323846;
    std::ofstream o(path);
    o << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << W << "\" height=\"" << H << "\" viewBox=\"0 0 " << W
      << " " << H << "\" font-family=\"Helvetica,Arial,sans-serif\">\n<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
    std::vector<double> X(n), Y(n);
    for (int v = 0; v < n; ++v) {
        const double a = -PI / 2 + 2 * PI * v / (n > 0 ? n : 1);
        X[v] = cx + R * std::cos(a);
        Y[v] = cy + R * std::sin(a);
    }
    for (auto [a, b] : edge_list(G))
        o << "<line x1=\"" << X[a - 1] << "\" y1=\"" << Y[a - 1] << "\" x2=\"" << X[b - 1] << "\" y2=\"" << Y[b - 1]
          << "\" stroke=\"#1f3b73\" stroke-width=\"2.2\"/>\n";
    const double r = n > 24 ? 9 : 13;
    for (int v = 0; v < n; ++v)
        o << "<circle cx=\"" << X[v] << "\" cy=\"" << Y[v] << "\" r=\"" << r
          << "\" fill=\"#f5f7fb\" stroke=\"#1f3b73\" stroke-width=\"1.8\"/>\n<text x=\"" << X[v] << "\" y=\"" << Y[v] + 4
          << "\" font-size=\"" << (n > 24 ? 9 : 14) << "\" text-anchor=\"middle\" fill=\"#111\">" << v + 1 << "</text>\n";
    o << "<text x=\"" << cx << "\" y=\"262\" font-size=\"13\" text-anchor=\"middle\" fill=\"#111\">G" << subscript_unicode(n)
      << "(" << poly_unicode(g) << ", " << poly_unicode(f) << ")</text>\n";
    if (!note.empty())
        o << "<text x=\"" << cx << "\" y=\"280\" font-size=\"11\" text-anchor=\"middle\" fill=\"#555\">" << note << "</text>\n";
    o << "</svg>\n";
}

}  // namespace rgraph
